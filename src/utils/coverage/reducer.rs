use anyhow::{Context, Result};
use fxhash::{FxBuildHasher, FxHashMap};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::io::{BufRead, BufReader, Write};

use crate::utils::coverage::tiled_run::{format_number_simplify, parse_tile_index, round_to};
use crate::utils::coverage::window_results::CoverageWindowAction;

/// One item from a per-tile partial stream
#[derive(Clone)]
struct Row {
    idx: u64,
    sum: f64,
    allowed: u64,
    blacklisted: u64,
}

/// Stream wrapper that yields `Row` from a compressed or plain TSV file
struct PartStream {
    reader: BufReader<Box<dyn std::io::Read + Send>>,
    buffer: String,
    /// Tile index for diagnostics
    tile_index: u32,
}

impl PartStream {
    fn open(path: &std::path::Path, tile_index: u32) -> Result<Self> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("Opening partials file {}", path.display()))?;
        // Detect .zst by extension. Support plain TSV for tests.
        let boxed: Box<dyn std::io::Read + Send> = if path
            .extension()
            .and_then(|s| s.to_str())
            .map(|e| e.eq_ignore_ascii_case("zst"))
            .unwrap_or(false)
        {
            Box::new(zstd::Decoder::new(file).context("Opening zstd decoder")?)
        } else {
            Box::new(file)
        };
        Ok(Self {
            reader: BufReader::new(boxed),
            buffer: String::new(),
            tile_index,
        })
    }

    /// Read next row, returning None on EOF
    fn next(&mut self) -> Result<Option<Row>> {
        self.buffer.clear();
        let n = self.reader.read_line(&mut self.buffer)?;
        if n == 0 {
            return Ok(None);
        }
        let raw = self.buffer.trim_end_matches('\n');
        if raw.is_empty() {
            return self.next(); // Skip empty lines
        }
        let mut it = raw.split('\t');
        let idx: u64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing idx column"))?
            .parse()
            .with_context(|| format!("Invalid idx in tile {}", self.tile_index))?;
        let sum: f64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing sum column"))?
            .parse()
            .with_context(|| format!("Invalid sum in tile {}", self.tile_index))?;
        let allowed: u64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing allowed column"))?
            .parse()
            .with_context(|| format!("Invalid allowed in tile {}", self.tile_index))?;
        let blacklisted: u64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing blacklisted column"))?
            .parse()
            .with_context(|| format!("Invalid blacklisted in tile {}", self.tile_index))?;
        Ok(Some(Row {
            idx,
            sum,
            allowed,
            blacklisted,
        }))
    }
}

/// Accumulator for a window aggregated across all contributing tiles
#[derive(Default)]
struct Acc {
    sum: f64,
    allowed: u64,
    blacklisted: u64,
    seen: u32,
}

/// Min-heap item for the K-way merge
#[derive(Clone)]
struct HeapItem {
    idx: u64,
    stream_id: usize,
    row: Row,
}

// Implement Ord for min-heap via Reverse wrapper at push time

/// Reduce `--by-bed` aggregates using sidecar `.cross` files to compute expected counts,
/// K-way merging the per-tile partials by `orig_idx` to guarantee output ordering.
///
/// Inputs:
/// - `partials_prefix`: Prefix used when writing `{prefix}.{chr}.{tile_index}.tsv.zst`
/// - `windows_chr`: Start-sorted `(start, end, orig_idx)` for this chromosome. `orig_idx` must be the
///                  start-sorted rank (0..n-1) so that per-tile partials are sorted by `orig_idx`.
///
/// Output columns (already written header outside):
/// `chromosome  start  end  value  blacklisted_positions`
pub fn reduce_bed_with_sidecars_for_chr(
    chr: &str,
    temp_dir: &std::path::Path,
    partials_prefix: &str,
    windows_chr: &[(u64, u64, u64)],
    masked: bool,
    mode: CoverageWindowAction, // Average | Total
    decimals: i32,
    out: &mut std::io::BufWriter<std::fs::File>,
) -> Result<()> {
    anyhow::ensure!(
        matches!(
            mode,
            CoverageWindowAction::Average | CoverageWindowAction::Total
        ),
        "This reducer supports only 'average' or 'total'"
    );

    // Build a quick lookup from orig_idx -> (start, end)
    // Assumes orig_idx == position in this array; we still guard against OOB.
    let n = windows_chr.len();
    let mut coords: Vec<(u64, u64)> = vec![(0, 0); n];
    for &(_s, _e, idx) in windows_chr {
        let i = idx as usize;
        anyhow::ensure!(i < n, "orig_idx {} out of bounds for {}", idx, chr);
    }
    for &(s, e, idx) in windows_chr {
        coords[idx as usize] = (s, e);
    }

    // Discover per-tile partial and sidecar files for this chromosome
    #[derive(Default, Clone)]
    struct Paths {
        part: Option<std::path::PathBuf>,
        cross: Option<std::path::PathBuf>,
    }
    let mut by_tile: FxHashMap<u32, Paths> = FxHashMap::with_hasher(FxBuildHasher::default());

    for entry in std::fs::read_dir(temp_dir)
        .with_context(|| format!("Listing temp dir {}", temp_dir.display()))?
    {
        let p = entry?.path();
        if !p.is_file() {
            continue;
        }
        let fname = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
        if !fname.starts_with(partials_prefix) || !fname.contains(&format!(".{chr}.")) {
            continue;
        }
        if let Some(idx) = parse_tile_index(fname) {
            if fname.ends_with(".tsv") || fname.ends_with(".tsv.zst") {
                by_tile.entry(idx).or_default().part = Some(p);
            } else if fname.ends_with(".cross") || fname.ends_with(".cross.zst") {
                by_tile.entry(idx).or_default().cross = Some(p);
            }
        }
    }

    // Compute expected contribution counts per orig_idx from sidecars
    let mut expected: FxHashMap<u64, u32> = FxHashMap::with_hasher(FxBuildHasher::default()); // only indices that cross need >1
    for (_t, paths) in by_tile.iter() {
        if let Some(cross_path) = &paths.cross {
            let f = std::fs::File::open(cross_path)
                .with_context(|| format!("Opening sidecar {}", cross_path.display()))?;
            let reader: Box<dyn std::io::Read + Send> =
                if cross_path.extension().and_then(|s| s.to_str()) == Some("zst") {
                    Box::new(zstd::Decoder::new(f)?)
                } else {
                    Box::new(f)
                };
            let mut r = BufReader::new(reader);
            let mut line = String::new();
            loop {
                line.clear();
                if r.read_line(&mut line)? == 0 {
                    break;
                }
                let raw = line.trim_end_matches('\n');
                if raw.is_empty() {
                    continue;
                }
                let idx: u64 = raw.parse().with_context(|| {
                    format!("Invalid orig_idx in sidecar {}", cross_path.display())
                })?;
                *expected.entry(idx).or_insert(0) += 1;
            }
        }
    }

    // Add 1 to all that appeared in sidecars, since sidecar counts tiles where the
    // window is not fully contained; fully-contained windows have implicit expected=1
    for (_idx, exp) in expected.iter_mut() {
        *exp = (*exp).max(1);
    }

    // Prepare K-way merge across all per-tile partial streams
    let mut streams: Vec<PartStream> = Vec::new();
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new(); // (idx, stream_id)
    for (tile_idx, paths) in by_tile.iter() {
        let Some(part_path) = &paths.part else {
            continue;
        };
        let mut ps = PartStream::open(part_path, *tile_idx)?;
        if let Some(row) = ps.next()? {
            let sid = streams.len();
            heap.push(Reverse((row.idx, sid)));
            // Store the row temporarily by pushing it back into the stream wrapper? Better approach: keep a side vector.
            // We keep a slot for "current row" per stream.
            streams.push(ps);
            // We also need a parallel vec of "current rows"
        } else {
            // Empty stream, skip
            continue;
        }
    }

    // Parallel storage for the current row of each stream
    let mut curr: Vec<Option<Row>> = vec![None; streams.len()];
    // Initialize current rows for streams we pushed
    // We need to re-open because we consumed one row above. Instead, take advantage:
    // Rewind: we already read the first row; let’s store it properly.
    // To do that, we must re-read it when we popped from heap. Adjust logic:
    //
    // Fix: Rebuild streams/heap with a helper that pushes the first row into `curr` and heap.

    // Reinitialize properly
    streams.clear();
    heap.clear();
    curr.clear();

    for (tile_idx, paths) in by_tile.iter() {
        let Some(part_path) = &paths.part else {
            continue;
        };
        let mut ps = PartStream::open(part_path, *tile_idx)?;
        if let Some(row) = ps.next()? {
            let sid = streams.len();
            streams.push(ps);
            curr.push(Some(row));
            heap.push(Reverse((curr[sid].as_ref().unwrap().idx, sid)));
        }
    }

    // Accumulator for indices currently in flight
    let mut acc: FxHashMap<u64, Acc> = FxHashMap::with_hasher(FxBuildHasher::default());

    // Emit helper
    let mut emit_idx = |idx: u64, acc_row: Acc| -> Result<()> {
        let (s, e) = coords[idx as usize];
        let value = match mode {
            CoverageWindowAction::Average => {
                if masked {
                    if acc_row.allowed == 0 {
                        0.0
                    } else {
                        acc_row.sum / acc_row.allowed as f64
                    }
                } else {
                    let span = (e - s) as f64;
                    if span == 0.0 { 0.0 } else { acc_row.sum / span }
                }
            }
            CoverageWindowAction::Total => acc_row.sum,
            _ => unreachable!(),
        };
        let value = round_to(value, decimals);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            chr,
            s,
            e,
            format_number_simplify(value, decimals),
            acc_row.blacklisted
        )?;
        Ok(())
    };

    // Merge loop
    while let Some(Reverse((_, sid))) = heap.pop() {
        // Take current row from that stream
        let row = curr[sid]
            .take()
            .expect("Heap and curr out of sync: missing current row");

        // Accumulate
        let ent = acc.entry(row.idx).or_insert_with(Acc::default);
        ent.sum += row.sum;
        ent.allowed += row.allowed;
        ent.blacklisted += row.blacklisted;
        ent.seen += 1;

        // Expected contributions for this idx
        let need = *expected.get(&row.idx).unwrap_or(&1);

        // If complete, emit and clear
        if ent.seen == need {
            // Safety: ensure coordinates exist
            anyhow::ensure!(
                (row.idx as usize) < coords.len(),
                "orig_idx {} out of bounds for {}",
                row.idx,
                chr
            );
            let done = acc.remove(&row.idx).unwrap();
            emit_idx(row.idx, done)?;
        }

        // Advance stream sid and push next row if any
        if let Some(next_row) = streams[sid].next()? {
            curr[sid] = Some(next_row);
            heap.push(Reverse((curr[sid].as_ref().unwrap().idx, sid)));
        }
    }

    // Ensure nothing remains partially accumulated
    anyhow::ensure!(
        acc.is_empty(),
        "Incomplete windows remain for {}: {}",
        chr,
        acc.len()
    );

    Ok(())
}
