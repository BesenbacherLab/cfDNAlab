use anyhow::{Context, Result};
use fxhash::{FxBuildHasher, FxHashMap};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::io::{BufRead, BufReader, Write};

use crate::commands::fcoverage::tiling::finalize_value;
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use crate::commands::fcoverage::writers::write_final_row;
use crate::shared::formatters::round_to;
use crate::shared::interval::{IndexedInterval, Interval};

fn open_partials_reader(path: &std::path::Path) -> Result<BufReader<Box<dyn std::io::Read + Send>>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("Opening partials file {}", path.display()))?;
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
    Ok(BufReader::new(boxed))
}

/// Row parsed from a non-summary BED partials file.
///
/// These rows carry only the fields needed to finish `average` and `total` outputs.
struct BasicPartialsRow {
    orig_idx: u64,
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
}

/// Row parsed from a summary-stats BED partials file.
///
/// Summary-stats reducers need the extra raw moments to derive variance, SD, CV, and covered
/// fraction without revisiting per-base coverage.
struct SummaryPartialsRow {
    orig_idx: u64,
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
}

/// Lightweight reader over a non-summary BED partials file.
struct BasicPartialsStream {
    reader: BufReader<Box<dyn std::io::Read + Send>>,
    line_buf: String,
    chr: String,
    line_number: u64,
    tile_index: u32, // For diagnostics
}

impl BasicPartialsStream {
    fn open(path: &std::path::Path, chr: &str, tile_index: u32) -> Result<Self> {
        Ok(Self {
            reader: open_partials_reader(path)?,
            line_buf: String::new(),
            chr: chr.to_string(),
            line_number: 0,
            tile_index,
        })
    }

    /// Read next row, or Ok(None) on EOF
    fn next_row(&mut self) -> Result<Option<BasicPartialsRow>> {
        self.line_buf.clear();
        let next_line_number = self.line_number + 1;
        let bytes_read = self.reader.read_line(&mut self.line_buf).with_context(|| {
            format!(
                "Reading partials for chromosome '{}' tile {} line {}",
                self.chr, self.tile_index, next_line_number
            )
        })?;
        if bytes_read == 0 {
            return Ok(None);
        }
        self.line_number = next_line_number;
        let raw = self.line_buf.trim_end_matches('\n');
        if raw.is_empty() {
            return self.next_row(); // Skip blank lines
        }
        let mut cols = raw.split('\t');

        // Expected columns in non-summary per-tile partials
        //   orig_idx   sum   allowed   blacklisted
        let orig_idx: u64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing orig_idx in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid orig_idx in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let sum: f64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing sum in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid sum in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let allowed_positions: u64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing allowed in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid allowed in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let blacklisted_positions: u64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing blacklisted in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid blacklisted in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;

        Ok(Some(BasicPartialsRow {
            orig_idx,
            sum,
            allowed_positions,
            blacklisted_positions,
        }))
    }
}

/// Lightweight reader over a summary-stats BED partials file.
struct SummaryPartialsStream {
    reader: BufReader<Box<dyn std::io::Read + Send>>,
    line_buf: String,
    chr: String,
    line_number: u64,
    tile_index: u32, // For diagnostics
}

impl SummaryPartialsStream {
    fn open(path: &std::path::Path, chr: &str, tile_index: u32) -> Result<Self> {
        Ok(Self {
            reader: open_partials_reader(path)?,
            line_buf: String::new(),
            chr: chr.to_string(),
            line_number: 0,
            tile_index,
        })
    }

    /// Read next row, or Ok(None) on EOF
    fn next_row(&mut self) -> Result<Option<SummaryPartialsRow>> {
        self.line_buf.clear();
        let next_line_number = self.line_number + 1;
        let bytes_read = self.reader.read_line(&mut self.line_buf).with_context(|| {
            format!(
                "Reading partials for chromosome '{}' tile {} line {}",
                self.chr, self.tile_index, next_line_number
            )
        })?;
        if bytes_read == 0 {
            return Ok(None);
        }
        self.line_number = next_line_number;
        let raw = self.line_buf.trim_end_matches('\n');
        if raw.is_empty() {
            return self.next_row(); // Skip blank lines
        }
        let mut cols = raw.split('\t');

        // Expected columns in summary-stats per-tile partials
        //   orig_idx   sum   allowed   blacklisted   nonzero_positions   coverage_sum_of_squares
        let orig_idx: u64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing orig_idx in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid orig_idx in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let sum: f64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing sum in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid sum in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let allowed_positions: u64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing allowed in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid allowed in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let blacklisted_positions: u64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing blacklisted in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid blacklisted in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let nonzero_positions: u64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing nonzero_positions in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid nonzero_positions in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let coverage_sum_of_squares: f64 = cols
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing coverage_sum_of_squares in partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid coverage_sum_of_squares in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;

        Ok(Some(SummaryPartialsRow {
            orig_idx,
            sum,
            allowed_positions,
            blacklisted_positions,
            nonzero_positions,
            coverage_sum_of_squares,
        }))
    }
}

/// Accumulator for a non-summary BED window across all contributing tiles.
#[derive(Default)]
struct BasicWindowAccum {
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    seen_contributions: u32,
}

/// Accumulator for a summary-stats BED window across all contributing tiles.
#[derive(Default)]
struct SummaryWindowAccum {
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
    seen_contributions: u32,
}

/// Final raw aggregate row after cross-tile reduction.
#[derive(Debug, Clone, Copy)]
pub struct ReducedAggregateRow {
    pub idx: u64,
    pub interval: Interval<u64>,
    pub coverage_sum: f64,
    pub eligible_positions: u64,
    pub blacklisted_positions: u64,
    pub nonzero_positions: u64,
    pub coverage_sum_of_squares: f64,
}

/// Reduce non-summary BED partials for one chromosome into complete raw window rows.
///
/// This follows the same cross-tile bookkeeping as the summary-stats BED reducer, but it keeps
/// only the fields that `average` and `total` still need after tile counting. The callback sees
/// one finalized raw row per original BED window, with the summary-only fields set to zero on
/// purpose.
///
/// The important rule is the same as for the summary reducer: accumulate by original BED index
/// until the window has received every tile contribution that the cross-index sidecars say to
/// expect. Windows that never appear in any cross-index file are assumed to have stayed inside one
/// tile core and therefore expect exactly one contribution.
///
/// Parameters
/// ----------
/// - `chr`:
///     Chromosome being reduced.
/// - `temp_dir`:
///     Directory that holds this chromosome's partial and cross-index files.
/// - `partials_prefix`:
///     Shared filename prefix used to discover the partial files for this reducer pass.
/// - `windows_chr`:
///     Original BED windows for this chromosome. Their `idx` values must match the `orig_idx`
///     that tile counting wrote into the partial rows.
/// - `on_row`:
///     Callback invoked once for each fully reduced raw row.
///
/// Returns
/// -------
/// - `()`:
///     Completes once every discovered partial row has been accounted for and every finished
///     window row has been handed to `on_row`.
pub(crate) fn reduce_bed_basic_with_cross_index_for_chr_rows(
    chr: &str,
    temp_dir: &std::path::Path,
    partials_prefix: &str,
    windows_chr: &[IndexedInterval<u64>],
    mut on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()> {
    // Keep a direct lookup from the stable BED row id to genomic coordinates.
    // The reducer only sees `orig_idx` in the partial rows, so it needs this side map to rebuild
    // the final interval once all contributions for one window have been merged.
    let mut coords_by_idx: FxHashMap<u64, Interval<u64>> =
        FxHashMap::with_capacity_and_hasher(windows_chr.len(), FxBuildHasher::default());
    for window in windows_chr {
        anyhow::ensure!(
            coords_by_idx
                .insert(window.idx(), window.interval)
                .is_none(),
            "duplicate orig_idx {} for {}",
            window.idx(),
            chr
        );
    }

    let files_by_tile = discover_tile_files_for_chr(temp_dir, chr, partials_prefix)?;

    // Count how many tile partials we should see for each `orig_idx`.
    //
    // Cross-index sidecars list only the windows that crossed tile-core boundaries. Everything
    // else defaults to one expected partial row.
    let mut expected_contribs: FxHashMap<u64, u32> =
        FxHashMap::with_hasher(FxBuildHasher::default());
    for (_tile_idx, tfs) in &files_by_tile {
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
    // If a window never appears in any cross-index file, it fit entirely inside one tile core
    // and therefore contributes exactly one partial row
    let expected_for = |idx: u64| -> u32 { *expected_contribs.get(&idx).unwrap_or(&1) };

    // Open one stream per partials file and keep only that stream's current row in memory.
    // The heap chooses the smallest visible `orig_idx` across those rows.
    //
    // `BinaryHeap` is a max-heap, so `Reverse((orig_idx, stream_id))` turns it into a min-heap.
    //
    // The heap does not prove that a window is complete. It only picks the next row to read.
    // Final correctness comes from `accum_by_idx` plus `expected_for(...)`, so this reducer does
    // not depend on the partials files arriving in perfectly grouped `orig_idx` blocks.
    let mut streams: Vec<BasicPartialsStream> = Vec::new();
    let mut current_row: Vec<Option<BasicPartialsRow>> = Vec::new();
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new(); // Min-heap via `Reverse`

    for (tile_idx, tfs) in &files_by_tile {
        let Some(partials_path) = &tfs.partials_path else {
            continue;
        };
        let mut ps = BasicPartialsStream::open(partials_path, chr, *tile_idx)?;
        if let Some(row) = ps.next_row()? {
            let stream_id = streams.len();
            streams.push(ps);
            current_row.push(Some(row));
            let key = current_row[stream_id].as_ref().unwrap().orig_idx;
            heap.push(Reverse((key, stream_id)));
        }
    }

    // One accumulator per BED window that has started receiving tile contributions but is not yet
    // complete. The basic reducer keeps only sum, eligible positions, and blacklist counts.
    let mut accum_by_idx: FxHashMap<u64, BasicWindowAccum> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    while let Some(Reverse((_, stream_id))) = heap.pop() {
        // Pull the current row for the stream that currently shows the smallest visible `orig_idx`.
        let row = current_row[stream_id].take().ok_or_else(|| {
            anyhow::anyhow!(
                "Reducer heap and current_row fell out of sync for chromosome '{}' stream {}",
                chr,
                stream_id
            )
        })?;

        // Add this tile's contribution into the running totals for the window.
        let entry = accum_by_idx.entry(row.orig_idx).or_default();
        entry.sum += row.sum;
        entry.allowed_positions += row.allowed_positions;
        entry.blacklisted_positions += row.blacklisted_positions;
        entry.seen_contributions += 1;

        // Write the row only once the cross-index says this window is complete.
        // This is what turns several tile-local overlaps back into one final BED row.
        if entry.seen_contributions == expected_for(row.orig_idx) {
            let done = accum_by_idx.remove(&row.orig_idx).ok_or_else(|| {
                anyhow::anyhow!(
                    "Reducer lost accumulated window state for chromosome '{}' orig_idx {}",
                    chr,
                    row.orig_idx
                )
            })?;
            let interval = *coords_by_idx.get(&row.orig_idx).ok_or_else(|| {
                anyhow::anyhow!(
                    "Reducer is missing interval coordinates for chromosome '{}' orig_idx {}",
                    chr,
                    row.orig_idx
                )
            })?;
            on_row(ReducedAggregateRow {
                idx: row.orig_idx,
                interval,
                coverage_sum: done.sum,
                eligible_positions: done.allowed_positions,
                blacklisted_positions: done.blacklisted_positions,
                // The shared callback row type also serves summary reducers. In the basic path
                // these summary-only fields were never written, so keep them explicitly zeroed.
                nonzero_positions: 0,
                coverage_sum_of_squares: 0.0,
            })?;
        }

        // Advance only the stream we just consumed from, then push its next visible row back into
        // the heap. Memory stays bounded by the number of open tile files.
        if let Some(next_row) = streams[stream_id].next_row()? {
            let next_key = next_row.orig_idx;
            current_row[stream_id] = Some(next_row);
            heap.push(Reverse((next_key, stream_id)));
        }
    }

    anyhow::ensure!(
        accum_by_idx.is_empty(),
        "Incomplete windows remain for {}: {}",
        chr,
        accum_by_idx.len()
    );

    Ok(())
}

/// Reduce summary-stats BED partials for one chromosome into complete raw window rows.
///
/// This is the BED summary-stats counterpart to the lighter aggregate reducers. Each tile writes
/// one partial row per overlapping window segment, and this function stitches those tile-local
/// contributions back into one exact row per original BED window. The callback receives raw
/// finalized moments, not already-derived variance or SD, so downstream writers can decide how to
/// format or further combine them.
///
/// The key idea is "accumulate by original BED index until we have seen every tile that was
/// supposed to contribute to that window". Cross-index sidecars tell us how many tile partials to
/// expect for boundary-crossing windows. Windows that stayed fully inside one tile core do not
/// appear in the cross-index and therefore default to one expected contribution.
///
/// Parameters
/// ----------
/// - `chr`:
///     Chromosome being reduced.
/// - `temp_dir`:
///     Directory that holds this chromosome's partial and cross-index files.
/// - `partials_prefix`:
///     Shared filename prefix used to discover the partial files for this reducer pass.
/// - `windows_chr`:
///     Original BED windows for this chromosome. Their `idx` values must match the `orig_idx`
///     that tile counting wrote into the partial rows.
/// - `on_row`:
///     Callback invoked once for each fully reduced raw row.
///
/// Returns
/// -------
/// - `()`:
///     Completes once every discovered partial row has been accounted for and every finished
///     window row has been handed to `on_row`.
pub(crate) fn reduce_bed_with_cross_index_for_chr_rows(
    chr: &str,
    temp_dir: &std::path::Path,
    partials_prefix: &str,
    windows_chr: &[IndexedInterval<u64>],
    mut on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()> {
    //
    // Keep a direct lookup from the stable BED row id to genomic coordinates.
    // The reducer only sees `orig_idx` in the partial rows, so it needs this side map to rebuild
    // the final interval once all contributions for one window have been merged.
    let mut coords_by_idx: FxHashMap<u64, Interval<u64>> =
        FxHashMap::with_capacity_and_hasher(windows_chr.len(), FxBuildHasher::default());
    for window in windows_chr {
        anyhow::ensure!(
            coords_by_idx
                .insert(window.idx(), window.interval)
                .is_none(),
            "duplicate orig_idx {} for {}",
            window.idx(),
            chr
        );
    }

    let files_by_tile = discover_tile_files_for_chr(temp_dir, chr, partials_prefix)?;

    // Count how many tile partials we should see for each `orig_idx`.
    //
    // Windows that crossed tile-core boundaries are listed once per contributing tile in the
    // cross-index sidecars. Ordinary windows that stayed inside one core never show up there, so
    // they implicitly expect exactly one partial row.
    let mut expected_contribs: FxHashMap<u64, u32> =
        FxHashMap::with_hasher(FxBuildHasher::default());
    for (_tile_idx, tfs) in &files_by_tile {
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
    // If a window never appears in any cross-index file, it fit entirely inside one tile core
    // and therefore contributes exactly one partial row
    let expected_for = |idx: u64| -> u32 { *expected_contribs.get(&idx).unwrap_or(&1) };

    // Open one stream per partials file and keep only that stream's current row in memory.
    // The heap chooses the smallest visible `orig_idx` across those rows.
    //
    // `BinaryHeap` is a max-heap, so `Reverse((orig_idx, stream_id))` turns it into a min-heap.
    //
    // The heap does not prove that a window is complete. It only picks the next row to read.
    // Final correctness comes from `accum_by_idx` plus `expected_for(...)`, so this reducer does
    // not depend on the partials files arriving in perfectly grouped `orig_idx` blocks.
    let mut streams: Vec<SummaryPartialsStream> = Vec::new();
    let mut current_row: Vec<Option<SummaryPartialsRow>> = Vec::new();
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new(); // Min-heap via `Reverse`

    for (tile_idx, tfs) in &files_by_tile {
        let Some(partials_path) = &tfs.partials_path else {
            continue;
        };
        let mut ps = SummaryPartialsStream::open(partials_path, chr, *tile_idx)?;
        if let Some(row) = ps.next_row()? {
            let stream_id = streams.len();
            streams.push(ps);
            current_row.push(Some(row));
            let key = current_row[stream_id].as_ref().unwrap().orig_idx;
            heap.push(Reverse((key, stream_id)));
        }
    }

    // One accumulator per BED window that has started receiving contributions but is not yet
    // complete. Summary-stats needs the full raw moments, so we sum every additive field here.
    let mut accum_by_idx: FxHashMap<u64, SummaryWindowAccum> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    while let Some(Reverse((_, stream_id))) = heap.pop() {
        // Pull the current row for the stream that currently shows the smallest visible `orig_idx`.
        let row = current_row[stream_id].take().ok_or_else(|| {
            anyhow::anyhow!(
                "Reducer heap and current_row fell out of sync for chromosome '{}' stream {}",
                chr,
                stream_id
            )
        })?;

        // Add this one tile contribution into the running totals for the window.
        // Every field here is additive across tiles because tile writing already clipped each row
        // to the tile core and preserved the original window identity in `orig_idx`.
        let entry = accum_by_idx.entry(row.orig_idx).or_default();
        entry.sum += row.sum;
        entry.allowed_positions += row.allowed_positions;
        entry.blacklisted_positions += row.blacklisted_positions;
        entry.nonzero_positions += row.nonzero_positions;
        entry.coverage_sum_of_squares += row.coverage_sum_of_squares;
        entry.seen_contributions += 1;

        // Write the row only after every expected tile contribution has arrived.
        // This is the key guard that lets cross-tile windows be reduced exactly once instead of
        // leaking partial intermediate rows downstream.
        if entry.seen_contributions == expected_for(row.orig_idx) {
            let done = accum_by_idx.remove(&row.orig_idx).ok_or_else(|| {
                anyhow::anyhow!(
                    "Reducer lost accumulated window state for chromosome '{}' orig_idx {}",
                    chr,
                    row.orig_idx
                )
            })?;
            let interval = *coords_by_idx.get(&row.orig_idx).ok_or_else(|| {
                anyhow::anyhow!(
                    "Reducer is missing interval coordinates for chromosome '{}' orig_idx {}",
                    chr,
                    row.orig_idx
                )
            })?;
            on_row(ReducedAggregateRow {
                idx: row.orig_idx,
                interval,
                coverage_sum: done.sum,
                eligible_positions: done.allowed_positions,
                blacklisted_positions: done.blacklisted_positions,
                nonzero_positions: done.nonzero_positions,
                coverage_sum_of_squares: done.coverage_sum_of_squares,
            })?;
        }

        // Advance only the stream we just consumed from, then push its next visible row back into
        // the heap. This is the standard K-way merge rhythm and keeps memory bounded by the number
        // of open tile files rather than the number of total rows.
        if let Some(next_row) = streams[stream_id].next_row()? {
            let next_key = next_row.orig_idx;
            current_row[stream_id] = Some(next_row);
            heap.push(Reverse((next_key, stream_id)));
        }
    }

    anyhow::ensure!(
        accum_by_idx.is_empty(),
        "Incomplete windows remain for {}: {}",
        chr,
        accum_by_idx.len()
    );

    Ok(())
}

/// Reduce non-summary BED aggregates for one chromosome using:
///  * Per-tile **partials** files
///  * Per-tile **cross-index** files: list `orig_idx` that are NOT fully contained in that tile core
///
/// Goal
///  * Merge cross-tile contributions back into full windows without buffering all windows or
///    scanning every possible index
///
/// About ordering
///  * The heap key is the current row's `orig_idx` from each open stream
///  * Correctness does not depend on every partials stream being globally sorted by `orig_idx`
///  * Rows are accumulated by `orig_idx` and written only after the expected contribution count
///    has been reached for that window
///  * When window indices were reindexed into coordinate order upstream, that also gives an
///    increasing-`orig_idx` output order
///  * Ordinary BED windows keep their original file indices, so callers must not rely on the
///    final output being sorted by `orig_idx`
///
/// Cross-index logic
///  * For windows fully contained in a single tile core: they appear in exactly one partials file
///    and are absent from all cross-index files -> expected contributions = 1
///  * For windows that cross tile core boundaries: the window appears in each overlapped tile's
///    partials file and is listed in each of those tiles' cross-index files -> expected contributions
///    equals the total number of tiles it overlaps
///
/// Requirements
///  * `windows_chr` must describe the same window identities that were written into the partials
///    files for this chromosome
///
/// Output columns already have their header written by the caller:
/// `chromosome  start  end  value  blacklisted_positions`
pub fn reduce_bed_with_cross_index_for_chr<W: Write>(
    chr: &str,
    temp_dir: &std::path::Path,
    partials_prefix: &str,
    windows_chr: &[IndexedInterval<u64>], // reindexed to 0..n-1
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

    reduce_bed_basic_with_cross_index_for_chr_rows(
        chr,
        temp_dir,
        partials_prefix,
        windows_chr,
        |row| {
            let interval = row.interval;
            let acc_sum = row.coverage_sum;
            let acc_allowed_positions = row.eligible_positions;
            let acc_blacklisted_positions = row.blacklisted_positions;
            let unmasked_span_bp = interval.len();
            let value = finalize_value(
                acc_sum,
                acc_allowed_positions,
                unmasked_span_bp,
                masked,
                &mode,
            );
            let value = round_to(value, decimals);
            write_final_row(
                final_writer,
                chr,
                interval,
                value,
                acc_blacklisted_positions,
                decimals,
            )?;
            Ok(())
        },
    )
}

/* By-size reducer (when windows don't align) */

/// One row from a non-summary size partials file.
struct BasicSizePartialsRow {
    interval: Interval<u64>,
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
}

/// One row from a summary-stats size partials file.
struct SummarySizePartialsRow {
    interval: Interval<u64>,
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
}

/// Lightweight streaming reader for non-summary size partials.
struct BasicSizePartialsStream {
    reader: BufReader<Box<dyn std::io::Read + Send>>,
    line_buf: String,
    chr: String,
    line_number: u64,
    tile_index: u32,
}

impl BasicSizePartialsStream {
    fn open(path: &std::path::Path, chr: &str, tile_index: u32) -> Result<Self> {
        Ok(Self {
            reader: open_partials_reader(path)?,
            line_buf: String::new(),
            chr: chr.to_string(),
            line_number: 0,
            tile_index,
        })
    }

    fn next_row(&mut self) -> Result<Option<BasicSizePartialsRow>> {
        self.line_buf.clear();
        let next_line_number = self.line_number + 1;
        let bytes_read = self.reader.read_line(&mut self.line_buf).with_context(|| {
            format!(
                "Reading size partials for chromosome '{}' tile {} line {}",
                self.chr, self.tile_index, next_line_number
            )
        })?;
        if bytes_read == 0 {
            return Ok(None);
        }
        self.line_number = next_line_number;
        let raw = self.line_buf.trim_end_matches('\n');
        if raw.is_empty() {
            return self.next_row();
        }
        let mut it = raw.split('\t');
        let start: u64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing start in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid start in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let end: u64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing end in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid end in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let sum: f64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing sum in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid sum in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let allowed_positions: u64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing allowed in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid allowed in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let blacklisted_positions: u64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing blacklisted in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid blacklisted in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        Ok(Some(BasicSizePartialsRow {
            interval: Interval::new(start, end)?,
            sum,
            allowed_positions,
            blacklisted_positions,
        }))
    }
}

/// Lightweight streaming reader for summary-stats size partials.
struct SummarySizePartialsStream {
    reader: BufReader<Box<dyn std::io::Read + Send>>,
    line_buf: String,
    chr: String,
    line_number: u64,
    tile_index: u32,
}

impl SummarySizePartialsStream {
    fn open(path: &std::path::Path, chr: &str, tile_index: u32) -> Result<Self> {
        Ok(Self {
            reader: open_partials_reader(path)?,
            line_buf: String::new(),
            chr: chr.to_string(),
            line_number: 0,
            tile_index,
        })
    }

    fn next_row(&mut self) -> Result<Option<SummarySizePartialsRow>> {
        self.line_buf.clear();
        let next_line_number = self.line_number + 1;
        let bytes_read = self.reader.read_line(&mut self.line_buf).with_context(|| {
            format!(
                "Reading size partials for chromosome '{}' tile {} line {}",
                self.chr, self.tile_index, next_line_number
            )
        })?;
        if bytes_read == 0 {
            return Ok(None);
        }
        self.line_number = next_line_number;
        let raw = self.line_buf.trim_end_matches('\n');
        if raw.is_empty() {
            return self.next_row();
        }
        let mut it = raw.split('\t');
        let start: u64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing start in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid start in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let end: u64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing end in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid end in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let sum: f64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing sum in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid sum in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let allowed_positions: u64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing allowed in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid allowed in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let blacklisted_positions: u64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing blacklisted in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid blacklisted in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let nonzero_positions: u64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing nonzero_positions in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid nonzero_positions in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        let coverage_sum_of_squares: f64 = it
            .next()
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "Missing coverage_sum_of_squares in size partials for chromosome '{}' tile {} line {}",
                    self.chr,
                    self.tile_index,
                    self.line_number
                )
            })?
            .parse()
            .with_context(|| {
                format!(
                    "Invalid coverage_sum_of_squares in chromosome '{}' tile {} line {}",
                    self.chr, self.tile_index, self.line_number
                )
            })?;
        Ok(Some(SummarySizePartialsRow {
            interval: Interval::new(start, end)?,
            sum,
            allowed_positions,
            blacklisted_positions,
            nonzero_positions,
            coverage_sum_of_squares,
        }))
    }
}

/// Accumulator per non-summary fixed-size bin.
#[derive(Default)]
struct BasicBinAccum {
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    seen_contributions: u32,
}

/// Accumulator per summary-stats fixed-size bin.
#[derive(Default)]
struct SummaryBinAccum {
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
    seen_contributions: u32,
}

/// Reduce summary-stats `--by-size` partials for one chromosome into complete raw bin rows.
///
/// Tile counting writes one partial row per fixed-size bin that overlaps a tile core. This
/// reducer merges those partial rows back into one raw summary-stats row per bin start. The
/// callback receives additive raw moments, not already-derived mean or variance, so later code can
/// derive the final statistics once the full bin is known.
///
/// The reducer keys bins by their full `start` coordinate. Cross-index sidecars tell us how
/// many tiles contributed to each bin start. When the last bin would extend past the chromosome
/// end, the reducer clips that final interval before handing it to the callback.
///
/// Parameters
/// ----------
/// - `chr`:
///     Chromosome being reduced.
/// - `temp_dir`:
///     Directory that holds this chromosome's partial and cross-index files.
/// - `partials_prefix`:
///     Shared filename prefix used to discover the partial files for this reducer pass.
/// - `chrom_len`:
///     True chromosome length, used to clip the final bin when the fixed bin size overruns the end.
/// - `on_row`:
///     Callback invoked once for each fully reduced raw bin row.
///
/// Returns
/// -------
/// - `()`:
///     Completes once every discovered partial row has been accounted for and every finished bin
///     row has been handed to `on_row`.
pub(crate) fn reduce_aggregates_by_size_with_cross_index_for_chr_rows(
    chr: &str,
    temp_dir: &std::path::Path,
    partials_prefix: &str,
    chrom_len: u64,
    mut on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()> {
    let files_by_tile = discover_tile_files_for_chr(temp_dir, chr, partials_prefix)?;

    // Count how many tile partials we should see for each full bin start.
    // Bins not listed in any cross-index file are assumed to have stayed inside one tile core.
    let mut expected_contribs: FxHashMap<u64, u32> =
        FxHashMap::with_hasher(FxBuildHasher::default());
    for (_idx, tfs) in &files_by_tile {
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
    // If a bin never appears in any cross-index file, it fit entirely inside one tile core and
    // therefore contributes exactly one partial row
    let expected_for = |start: u64| -> u32 { *expected_contribs.get(&start).unwrap_or(&1) };

    // Open one stream per partials file and keep only that stream's current row in memory.
    // The heap chooses the smallest visible full bin start across those rows.
    //
    // `BinaryHeap` is a max-heap, so `Reverse((bin_start, stream_id))` turns it into a min-heap.
    //
    // The heap key must be the full bin start, not this tile's clipped overlap. Different tile
    // contributions for the same bin need to meet in the same accumulator entry. Final correctness
    // comes from `accum_by_start` plus `expected_for(...)`, not from assuming one stream already
    // contains every row for that bin.
    let mut streams: Vec<SummarySizePartialsStream> = Vec::new();
    let mut current_row: Vec<Option<SummarySizePartialsRow>> = Vec::new();
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new(); // Min-heap via `Reverse`

    for (tile_idx, tfs) in &files_by_tile {
        let Some(partials_path) = &tfs.partials_path else {
            continue;
        };
        let mut ps = SummarySizePartialsStream::open(partials_path, chr, *tile_idx)?;
        if let Some(row) = ps.next_row()? {
            let sid = streams.len();
            streams.push(ps);
            current_row.push(Some(row));
            let start_key = current_row[sid].as_ref().unwrap().interval.start();
            heap.push(Reverse((start_key, sid)));
        }
    }

    // One accumulator per full bin start that has begun receiving tile contributions but is not
    // complete yet. Summary-stats needs every additive field.
    let mut accum_by_start: FxHashMap<u64, SummaryBinAccum> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    while let Some(Reverse((_, sid))) = heap.pop() {
        // Pull the current row for the stream that currently shows the smallest visible bin start.
        let row = current_row[sid].take().ok_or_else(|| {
            anyhow::anyhow!(
                "Reducer heap and current_row fell out of sync for chromosome '{}' stream {}",
                chr,
                sid
            )
        })?;

        // Add this tile contribution into the raw running totals for the bin keyed by `start`.
        let entry = accum_by_start.entry(row.interval.start()).or_default();
        entry.sum += row.sum;
        entry.allowed_positions += row.allowed_positions;
        entry.blacklisted_positions += row.blacklisted_positions;
        entry.nonzero_positions += row.nonzero_positions;
        entry.coverage_sum_of_squares += row.coverage_sum_of_squares;
        entry.seen_contributions += 1;

        // Write the row only once the bin has received all expected tile contributions.
        if entry.seen_contributions == expected_for(row.interval.start()) {
            let done = accum_by_start
                .remove(&row.interval.start())
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "Reducer lost accumulated size-bin state for chromosome '{}' start {}",
                        chr,
                        row.interval.start()
                    )
                })?;
            // The final fixed-size bin can overrun the chromosome end. Clip it here so downstream
            // writers see the true genomic interval length for the last row.
            let clipped_interval =
                Interval::new(row.interval.start(), row.interval.end().min(chrom_len))?;
            on_row(ReducedAggregateRow {
                idx: row.interval.start(),
                interval: clipped_interval,
                coverage_sum: done.sum,
                eligible_positions: done.allowed_positions,
                blacklisted_positions: done.blacklisted_positions,
                nonzero_positions: done.nonzero_positions,
                coverage_sum_of_squares: done.coverage_sum_of_squares,
            })?;
        }

        // Advance only the stream we just consumed from, then push its next visible row back into
        // the heap.
        if let Some(next_row) = streams[sid].next_row()? {
            let next_key = next_row.interval.start();
            current_row[sid] = Some(next_row);
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

/// Reduce non-summary `--by-size` partials for one chromosome into complete raw bin rows.
///
/// This is the lighter fixed-bin counterpart to the summary-stats size reducer. It merges tile
/// partials back into one final row per bin start while keeping only the fields that `average`
/// and `total` still need. The callback receives raw totals, with the summary-only fields set to
/// zero on purpose.
///
/// Bins are keyed by their full `start` coordinate, not by the overlap clipped to one tile.
/// Cross-index sidecars still define how many tile contributions each bin start is expected to
/// receive. As in the summary reducer, the final bin is clipped to the true chromosome end before
/// it is handed to the callback.
///
/// Parameters
/// ----------
/// - `chr`:
///     Chromosome being reduced.
/// - `temp_dir`:
///     Directory that holds this chromosome's partial and cross-index files.
/// - `partials_prefix`:
///     Shared filename prefix used to discover the partial files for this reducer pass.
/// - `chrom_len`:
///     True chromosome length, used to clip the final bin when the fixed bin size overruns the end.
/// - `on_row`:
///     Callback invoked once for each fully reduced raw bin row.
///
/// Returns
/// -------
/// - `()`:
///     Completes once every discovered partial row has been accounted for and every finished bin
///     row has been handed to `on_row`.
pub(crate) fn reduce_aggregates_by_size_basic_with_cross_index_for_chr_rows(
    chr: &str,
    temp_dir: &std::path::Path,
    partials_prefix: &str,
    chrom_len: u64,
    mut on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()> {
    let files_by_tile = discover_tile_files_for_chr(temp_dir, chr, partials_prefix)?;

    // Count how many tile partials we should see for each full bin start.
    let mut expected_contribs: FxHashMap<u64, u32> =
        FxHashMap::with_hasher(FxBuildHasher::default());
    for (_idx, tfs) in &files_by_tile {
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
    // If a bin never appears in any cross-index file, it fit entirely inside one tile core and
    // therefore contributes exactly one partial row
    let expected_for = |start: u64| -> u32 { *expected_contribs.get(&start).unwrap_or(&1) };

    // Open one stream per partials file and keep only that stream's current row in memory.
    // The heap chooses the smallest visible full bin start across those rows.
    //
    // `BinaryHeap` is a max-heap, so `Reverse((bin_start, stream_id))` turns it into a min-heap.
    //
    // The heap key must be the full bin start, not this tile's clipped overlap. Different tile
    // contributions for the same bin need to meet in the same accumulator entry. Final correctness
    // comes from `accum_by_start` plus `expected_for(...)`, not from assuming one stream already
    // contains every row for that bin.
    let mut streams: Vec<BasicSizePartialsStream> = Vec::new();
    let mut current_row: Vec<Option<BasicSizePartialsRow>> = Vec::new();
    let mut heap: BinaryHeap<Reverse<(u64, usize)>> = BinaryHeap::new(); // Min-heap via `Reverse`

    for (tile_idx, tfs) in &files_by_tile {
        let Some(partials_path) = &tfs.partials_path else {
            continue;
        };
        let mut ps = BasicSizePartialsStream::open(partials_path, chr, *tile_idx)?;
        if let Some(row) = ps.next_row()? {
            let sid = streams.len();
            streams.push(ps);
            current_row.push(Some(row));
            let start_key = current_row[sid].as_ref().unwrap().interval.start();
            heap.push(Reverse((start_key, sid)));
        }
    }

    // One accumulator per full bin start that has begun receiving tile contributions but is not
    // complete yet.
    let mut accum_by_start: FxHashMap<u64, BasicBinAccum> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    while let Some(Reverse((_, sid))) = heap.pop() {
        // Pull the current row for the stream that currently shows the smallest visible bin start.
        let row = current_row[sid].take().ok_or_else(|| {
            anyhow::anyhow!(
                "Reducer heap and current_row fell out of sync for chromosome '{}' stream {}",
                chr,
                sid
            )
        })?;

        // Add this tile contribution into the running totals for the bin keyed by `start`.
        let entry = accum_by_start.entry(row.interval.start()).or_default();
        entry.sum += row.sum;
        entry.allowed_positions += row.allowed_positions;
        entry.blacklisted_positions += row.blacklisted_positions;
        entry.seen_contributions += 1;

        // Write the row only once the cross-index says this bin has received every contributing
        // tile row.
        if entry.seen_contributions == expected_for(row.interval.start()) {
            let done = accum_by_start
                .remove(&row.interval.start())
                .ok_or_else(|| {
                    anyhow::anyhow!(
                        "Reducer lost accumulated size-bin state for chromosome '{}' start {}",
                        chr,
                        row.interval.start()
                    )
                })?;
            // The final fixed-size bin can overrun the chromosome end. Clip it here so downstream
            // writers see the true genomic interval length for the last row.
            let clipped_interval =
                Interval::new(row.interval.start(), row.interval.end().min(chrom_len))?;
            on_row(ReducedAggregateRow {
                idx: row.interval.start(),
                interval: clipped_interval,
                coverage_sum: done.sum,
                eligible_positions: done.allowed_positions,
                blacklisted_positions: done.blacklisted_positions,
                // The shared callback row type also serves summary reducers. In the basic path
                // these summary-only fields were never written, so keep them explicitly zeroed.
                nonzero_positions: 0,
                coverage_sum_of_squares: 0.0,
            })?;
        }

        // Advance only the stream we just consumed from, then push its next visible row back into
        // the heap.
        if let Some(next_row) = streams[sid].next_row()? {
            let next_key = next_row.interval.start();
            current_row[sid] = Some(next_row);
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

/// Reduce non-summary `--by-size` partials for one chromosome in strictly ascending `start` order.
///
/// Ordering is guaranteed by a K-way merge across sorted per-tile partials.
/// A priority queue (`BinaryHeap`) is used as a min-heap via `Reverse((start, stream_id))`:
/// the smallest start is popped first. This keeps peak memory low while preserving order.
///
/// The cross-index counts how many tiles contribute to each full bin start:
/// - If a bin is not listed in any cross-index file, we expect exactly 1 contribution.
/// - If it appears N times, we wait for N contributions before writing that bin.
///
/// Important
/// - The partial rows must carry the logical `bin_start` and `bin_end`, not the clipped
///   tile-local overlap for that bin. The reducer keys on `start`, so changing the
///   partial row bounds to clipped pieces will silently break cross-tile merging.
///
/// The final bin is truncated to the chromosome end; it may be shorter than window_bp.
pub fn reduce_aggregates_by_size_with_cross_index_for_chr<W: Write>(
    chr: &str,
    temp_dir: &std::path::Path,
    partials_prefix: &str,
    masked: bool,
    mode: CoverageWindowAction, // Average | Total
    chrom_len: u64,
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

    reduce_aggregates_by_size_basic_with_cross_index_for_chr_rows(
        chr,
        temp_dir,
        partials_prefix,
        chrom_len,
        |row| {
            let interval = row.interval;
            let unmasked_span_bp = interval.len();
            debug_assert!(unmasked_span_bp >= 1);
            let value = finalize_value(
                row.coverage_sum,
                row.eligible_positions,
                unmasked_span_bp,
                masked,
                &mode,
            );
            let value = round_to(value, decimals);
            write_final_row(
                out,
                chr,
                interval,
                value,
                row.blacklisted_positions,
                decimals,
            )?;
            Ok(())
        },
    )
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
        if let Some(tile_idx) = crate::shared::tiled_run::parse_tile_index(fname) {
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
