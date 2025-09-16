use crate::utils::coverage::tiled_run::{finalize_value, round_to};
use crate::utils::coverage::window_results::CoverageWindowAction;
use crate::utils::coverage::writer::write_final_row;
use anyhow::{Context, Result};
use fxhash::{FxBuildHasher, FxHashMap};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::io::{BufRead, BufReader, Write};

/// Row parsed from a per-tile partials file
struct PartialsRow {
    orig_idx: u64,
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
}

/// Lightweight reader over a compressed or plain per-tile partials file
struct PartialsStream {
    reader: BufReader<Box<dyn std::io::Read + Send>>,
    line_buf: String,
    tile_index: u32, // For diagnostics
}

impl PartialsStream {
    fn open(path: &std::path::Path, tile_index: u32) -> Result<Self> {
        let file = std::fs::File::open(path)
            .with_context(|| format!("Opening partials file {}", path.display()))?;
        // Detect zstd by extension; allow plain TSV for tests
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
            line_buf: String::new(),
            tile_index,
        })
    }

    /// Read next row, or Ok(None) on EOF
    fn next_row(&mut self) -> Result<Option<PartialsRow>> {
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            return Ok(None);
        }
        let raw = self.line_buf.trim_end_matches('\n');
        if raw.is_empty() {
            return self.next_row(); // Skip blank lines
        }
        let mut cols = raw.split('\t');

        // Expected columns in per-tile partials
        //   orig_idx   sum   allowed   blacklisted
        let orig_idx: u64 = cols
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing orig_idx in partials"))?
            .parse()
            .with_context(|| format!("Invalid orig_idx in tile {}", self.tile_index))?;
        let sum: f64 = cols
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing sum in partials"))?
            .parse()
            .with_context(|| format!("Invalid sum in tile {}", self.tile_index))?;
        let allowed_positions: u64 = cols
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing allowed in partials"))?
            .parse()
            .with_context(|| format!("Invalid allowed in tile {}", self.tile_index))?;
        let blacklisted_positions: u64 = cols
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing blacklisted in partials"))?
            .parse()
            .with_context(|| format!("Invalid blacklisted in tile {}", self.tile_index))?;

        Ok(Some(PartialsRow {
            orig_idx,
            sum,
            allowed_positions,
            blacklisted_positions,
        }))
    }
}

/// Accumulator for a window across all contributing tiles
#[derive(Default)]
struct WindowAccum {
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    seen_contributions: u32,
}

/// Reduce aggregates for one chromosome using:
///  * Per-tile **partials** files: rows are sorted by `orig_idx`
///  * Per-tile **cross-index** files: list `orig_idx` that are NOT fully contained in that tile core
///
/// Goal
///  * Emit final rows strictly in **increasing `orig_idx`** order without
///    buffering all windows or scanning every window
///
/// How ordering is guaranteed
///  * Each tile’s partials are sorted by `orig_idx`
///  * We open all tile streams and perform a **K-way merge** on `orig_idx`
///
/// What is a K-way merge and why use `BinaryHeap<Reverse<...>>`
///  * A K-way merge picks the smallest next key across K sorted inputs
///  * A binary heap is a priority queue. In Rust, `BinaryHeap` is a max-heap by default
///  * Wrapping keys in `Reverse(..)` flips the order so the heap behaves like a **min-heap**
///    (smallest key at the top), ideal for always extracting the lowest `orig_idx` next
///
/// Cross-index logic
///  * For windows fully contained in a single tile core: they appear in exactly one partials file
///    and are absent from all cross-index files -> expected contributions = 1
///  * For windows that cross tile core boundaries: the window appears in each overlapped tile’s
///    partials file and is listed in each of those tiles’ cross-index files -> expected contributions
///    equals the total number of tiles it overlaps
///
/// Requirements
///  * `windows_chr` must be start-sorted and must have `orig_idx == start-sorted rank`
///    so that per-tile partials are naturally sorted by `orig_idx` with no extra work
///
/// Output columns already have their header written by the caller:
/// `chromosome  start  end  value  blacklisted_positions`
pub fn reduce_bed_with_cross_index_for_chr<W: Write>(
    chr: &str,
    temp_dir: &std::path::Path,
    partials_prefix: &str,
    windows_chr: &[(u64, u64, u64)], // (start, end, orig_idx) reindexed to 0..n-1
    masked: bool,
    mode: CoverageWindowAction, // Average | Total
    decimals: i32,
    final_writer: &mut W,
) -> Result<()> {
    anyhow::ensure!(
        matches!(
            mode,
            CoverageWindowAction::Average | CoverageWindowAction::Total
        ),
        "Reducer supports only 'average' or 'total'"
    );

    // Map orig_idx -> (start, end) so we can compute averages without storing starts/ends in partials
    // Assumes orig_idx is the local rank 0..n-1; we still assert bounds for safety
    let n_windows = windows_chr.len();
    let mut coords_by_idx: Vec<(u64, u64)> = vec![(0, 0); n_windows];
    for &(start, end, orig_idx) in windows_chr {
        let i = orig_idx as usize;
        anyhow::ensure!(
            i < n_windows,
            "orig_idx {} out of bounds for {}",
            orig_idx,
            chr
        );
        coords_by_idx[i] = (start, end);
    }

    // Extract files from temp dir
    let files_by_tile = discover_tile_files_for_chr(temp_dir, chr, partials_prefix)?;

    // Compute expected contribution counts per orig_idx from cross-index files
    // Windows not present in any cross-index file implicitly expect exactly 1 contribution
    let mut expected_contribs: FxHashMap<u64, u32> =
        FxHashMap::with_hasher(FxBuildHasher::default());
    for (_tile_idx, tfs) in files_by_tile.iter() {
        if let Some(cross_path) = &tfs.cross_index_path {
            let f = std::fs::File::open(cross_path)
                .with_context(|| format!("Opening cross-index {}", cross_path.display()))?;
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
                let idx: u64 = raw
                    .parse()
                    .with_context(|| format!("Invalid orig_idx in {}", cross_path.display()))?;
                *expected_contribs.entry(idx).or_insert(0) += 1;
            }
        }
    }
    // Note: Missing => 1, Present => exact number of overlapping tiles
    let expected_for = |idx: u64| -> u32 { *expected_contribs.get(&idx).unwrap_or(&1) };

    // Prepare K-way merge across all per-tile partials streams
    //
    // Data structure choice
    //  * `BinaryHeap` is a priority queue; by wrapping the key in `Reverse((key, id))` we
    //    make it behave as a min-heap, so `pop()` always returns the smallest `orig_idx`
    //
    // Invariant
    //  * Each stream yields rows in ascending `orig_idx` order
    let mut streams: Vec<PartialsStream> = Vec::new();
    let mut current_row: Vec<Option<PartialsRow>> = Vec::new();
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new(); // Reverse to get min-heap behavior

    // Open streams and push their first row into the heap
    for (tile_idx, tfs) in files_by_tile.iter() {
        let Some(partials_path) = &tfs.partials_path else {
            continue;
        };
        let mut ps = PartialsStream::open(partials_path, *tile_idx)?;
        if let Some(row) = ps.next_row()? {
            let stream_id = streams.len();
            streams.push(ps);
            current_row.push(Some(row));
            let key = current_row[stream_id].as_ref().unwrap().orig_idx;
            heap.push(Reverse((key, stream_id)));
        }
    }

    // Accumulators for indices currently “in flight”
    let mut accum_by_idx: FxHashMap<u64, WindowAccum> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    // Emit helper
    let mut emit_idx = |orig_idx: u64, acc: WindowAccum| -> Result<()> {
        let (start, end) = coords_by_idx[orig_idx as usize];
        let unmasked_span_bp = (end - start) as u64;
        let value = finalize_value(
            acc.sum,
            acc.allowed_positions,
            unmasked_span_bp,
            masked,
            &mode,
        );
        let value = round_to(value, decimals);
        write_final_row(
            final_writer,
            &chr,
            start,
            end,
            value,
            acc.blacklisted_positions,
            decimals,
        )?;
        Ok(())
    };

    // Merge loop: always take the smallest available orig_idx across streams
    while let Some(Reverse((_, stream_id))) = heap.pop() {
        let row = current_row[stream_id]
            .take()
            .expect("Heap and current_row out of sync");

        // Accumulate this contribution
        let entry = accum_by_idx.entry(row.orig_idx).or_default();
        entry.sum += row.sum;
        entry.allowed_positions += row.allowed_positions;
        entry.blacklisted_positions += row.blacklisted_positions;
        entry.seen_contributions += 1;

        // If we have collected all expected contributions for this window, emit immediately
        if entry.seen_contributions == expected_for(row.orig_idx) {
            let done = accum_by_idx.remove(&row.orig_idx).unwrap();
            emit_idx(row.orig_idx, done)?;
        }

        // Advance this stream and re-insert into the heap
        if let Some(next_row) = streams[stream_id].next_row()? {
            current_row[stream_id] = Some(next_row);
            let next_key = current_row[stream_id].as_ref().unwrap().orig_idx;
            heap.push(Reverse((next_key, stream_id)));
        }
    }

    // Safety check
    anyhow::ensure!(
        accum_by_idx.is_empty(),
        "Incomplete windows remain for {}: {}",
        chr,
        accum_by_idx.len()
    );

    Ok(())
}

/* By-size reducer (when windows don't align) */

/// One row from a size-based partials file.
struct SizePartialsRow {
    start: u64,
    end: u64,
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
}

/// Lightweight streaming reader for partials (compressed or plain).
struct SizePartialsStream {
    reader: BufReader<Box<dyn std::io::Read + Send>>,
    line_buf: String,
    tile_index: u32,
}

impl SizePartialsStream {
    fn open(path: &std::path::Path, tile_index: u32) -> Result<Self> {
        let f = std::fs::File::open(path)
            .with_context(|| format!("Opening partials file {}", path.display()))?;
        let boxed: Box<dyn std::io::Read + Send> = if path
            .extension()
            .and_then(|s| s.to_str())
            .map(|e| e.eq_ignore_ascii_case("zst"))
            .unwrap_or(false)
        {
            Box::new(zstd::Decoder::new(f).context("Opening zstd decoder")?)
        } else {
            Box::new(f)
        };
        Ok(Self {
            reader: BufReader::new(boxed),
            line_buf: String::new(),
            tile_index: tile_index,
        })
    }

    fn next_row(&mut self) -> Result<Option<SizePartialsRow>> {
        self.line_buf.clear();
        if self.reader.read_line(&mut self.line_buf)? == 0 {
            return Ok(None);
        }
        let raw = self.line_buf.trim_end_matches('\n');
        if raw.is_empty() {
            return self.next_row();
        }
        let mut it = raw.split('\t');
        let start: u64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing start"))?
            .parse()
            .with_context(|| format!("Invalid start in tile {}", self.tile_index))?;
        let end: u64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing end"))?
            .parse()
            .with_context(|| format!("Invalid end in tile {}", self.tile_index))?;
        let sum: f64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing sum"))?
            .parse()
            .with_context(|| format!("Invalid sum in tile {}", self.tile_index))?;
        let allowed_positions: u64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing allowed"))?
            .parse()
            .with_context(|| format!("Invalid allowed in tile {}", self.tile_index))?;
        let blacklisted_positions: u64 = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("Missing blacklisted"))?
            .parse()
            .with_context(|| format!("Invalid blacklisted in tile {}", self.tile_index))?;
        Ok(Some(SizePartialsRow {
            start,
            end,
            sum,
            allowed_positions,
            blacklisted_positions,
        }))
    }
}

/// Accumulator per fixed-size bin.
#[derive(Default)]
struct BinAccum {
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    seen_contributions: u32,
}

/// Reduce `--by-size` partials for one chromosome in strictly ascending `start` order.
///
/// Ordering is guaranteed by a K-way merge across sorted per-tile partials.
/// A priority queue (`BinaryHeap`) is used as a min-heap via `Reverse((start, stream_id))`:
/// the smallest start is popped first. This keeps peak memory low while preserving order.
///
/// The cross-index counts how many tiles contribute to each bin start:
/// - If a bin is not listed in any cross-index file, we expect exactly 1 contribution.
/// - If it appears N times, we expect N contributions before emitting that bin.
pub fn reduce_aggregates_by_size_with_cross_index_for_chr<W: Write>(
    chr: &str,
    temp_dir: &std::path::Path,
    partials_prefix: &str,
    masked: bool,
    mode: CoverageWindowAction, // Average | Total
    decimals: i32,
    out: &mut W,
) -> Result<()> {
    anyhow::ensure!(
        matches!(
            mode,
            CoverageWindowAction::Average | CoverageWindowAction::Total
        ),
        "Reducer supports only 'average' or 'total'"
    );

    // Extract files from temp dir
    let files_by_tile = discover_tile_files_for_chr(temp_dir, chr, partials_prefix)?;

    // Build expected contribution counts per bin start from cross-index files
    let mut expected_contribs: FxHashMap<u64, u32> =
        FxHashMap::with_hasher(FxBuildHasher::default());
    for (_idx, tfs) in files_by_tile.iter() {
        if let Some(cross_path) = &tfs.cross_index_path {
            let f = std::fs::File::open(cross_path)
                .with_context(|| format!("Opening cross-index {}", cross_path.display()))?;
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
                let start: u64 = raw
                    .parse()
                    .with_context(|| format!("Invalid start in {}", cross_path.display()))?;
                *expected_contribs.entry(start).or_insert(0) += 1;
            }
        }
    }
    let expected_for = |start: u64| -> u32 { *expected_contribs.get(&start).unwrap_or(&1) };

    // Prepare the merge structures
    let mut streams: Vec<SizePartialsStream> = Vec::new();
    let mut current_row: Vec<Option<SizePartialsRow>> = Vec::new();
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new(); // Min-heap by start

    // Open each partials stream and seed the heap with its first row
    for (tile_idx, tfs) in files_by_tile.iter() {
        let Some(partials_path) = &tfs.partials_path else {
            continue;
        };
        let mut ps = SizePartialsStream::open(partials_path, *tile_idx)?;
        if let Some(row) = ps.next_row()? {
            let sid = streams.len();
            streams.push(ps);
            current_row.push(Some(row));
            let start_key = current_row[sid].as_ref().unwrap().start;
            heap.push(Reverse((start_key, sid)));
        }
    }

    // Accumulate contributions per start bin until the expected count is reached
    let mut accum_by_start: FxHashMap<u64, BinAccum> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    // Emit helper for one completed bin
    let mut emit_bin = |start: u64, acc: BinAccum, end: u64| -> Result<()> {
        let unmasked_span_bp = (end - start) as u64;
        let value = finalize_value(
            acc.sum,
            acc.allowed_positions,
            unmasked_span_bp,
            masked,
            &mode,
        );
        let value = round_to(value, decimals);
        write_final_row(
            out,
            &chr,
            start,
            end,
            value,
            acc.blacklisted_positions,
            decimals,
        )?;
        Ok(())
    };

    // K-way merge loop
    while let Some(Reverse((_, sid))) = heap.pop() {
        let row = current_row[sid]
            .take()
            .expect("heap and current_row out of sync");

        let entry = accum_by_start.entry(row.start).or_default();
        entry.sum += row.sum;
        entry.allowed_positions += row.allowed_positions;
        entry.blacklisted_positions += row.blacklisted_positions;
        entry.seen_contributions += 1;

        // Emit when we have all expected contributions for this bin start
        if entry.seen_contributions == expected_for(row.start) {
            let done = accum_by_start.remove(&row.start).unwrap();
            // Use the end from the last seen row (all rows for a bin share the same [start,end) by construction)
            emit_bin(row.start, done, row.end)?;
        }

        // Advance this stream and push next row if present
        if let Some(next_row) = streams[sid].next_row()? {
            current_row[sid] = Some(next_row);
            let next_key = current_row[sid].as_ref().unwrap().start;
            heap.push(Reverse((next_key, sid)));
        }
    }

    anyhow::ensure!(
        accum_by_start.is_empty(),
        "Incomplete size bins remain for {}: {}",
        chr,
        accum_by_start.len()
    );

    Ok(())
}

/// Sidecar type information for one tile
#[derive(Default, Clone)]
struct TileFiles {
    pub partials_path: Option<std::path::PathBuf>,
    pub cross_index_path: Option<std::path::PathBuf>,
}

/// Find per-tile files for a chromosome
///
/// Definitions
/// - **Partials**: the per-tile contributions for windows/bins (tsv or tsv.zst)
/// - **Cross-index**: a side list that marks which windows/bins cross tile core boundaries
///   The reducer uses it to know how many tile contributions to expect for each window/bin
///
/// Filenames
/// - Must start with `per_tile_prefix` and contain `.{chr}.`
/// - We detect `.cross.` in the name to classify the sidecar
fn discover_tile_files_for_chr(
    temp_dir: &std::path::Path,
    chr: &str,
    per_tile_prefix: &str,
) -> anyhow::Result<FxHashMap<u32, TileFiles>> {
    let mut files_by_tile: FxHashMap<u32, TileFiles> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    for entry in std::fs::read_dir(temp_dir)? {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        let fname = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
        if !fname.starts_with(per_tile_prefix) || !fname.contains(&format!(".{chr}.")) {
            continue;
        }
        if let Some(tile_idx) = crate::utils::coverage::tiled_run::parse_tile_index(fname) {
            // Recognize cross files by the marker in the name (simple and robust)
            if fname.contains(".cross.") {
                files_by_tile.entry(tile_idx).or_default().cross_index_path = Some(path);
            } else if fname.ends_with(".tsv") || fname.ends_with(".tsv.zst") {
                files_by_tile.entry(tile_idx).or_default().partials_path = Some(path);
            }
        }
    }
    Ok(files_by_tile)
}
