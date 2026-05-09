use anyhow::{Context, Result};
use fxhash::{FxBuildHasher, FxHashMap};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fmt::Display;
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use crate::shared::interval::{IndexedInterval, Interval};
use crate::shared::io::open_text_reader;

type StreamHeap = BinaryHeap<Reverse<(u64, usize)>>;

/// Returned aggregate temp paths for one tile.
///
/// `partials_path` points at tile-local contribution rows, not BED records. The row identity
/// depends on the reducer family. BED rows are keyed by stable `orig_idx`, and fixed-size rows are
/// keyed by the logical full-bin start. `cross_index_path` is optional because aligned single-tile
/// contributions use the reducer rule that missing cross-index entries mean one contribution.
#[derive(Debug, Clone)]
pub(crate) struct TileAggregateTempFiles {
    pub tile_index: u32,
    pub partials_path: PathBuf,
    pub cross_index_path: Option<PathBuf>,
}

/// Parse one tab-delimited column with consistent missing/invalid diagnostics.
///
/// The field names still stay explicit at each call site, but the repetitive parse scaffolding no
/// longer has to be copied into every schema branch. This is the main low-risk parser cleanup from
/// the refactor plan.
fn parse_col<T: FromStr>(
    cols: &mut std::str::Split<'_, char>,
    field_name: &str,
    partials_label: &str,
    chr: &str,
    tile_index: u32,
    line_number: u64,
) -> Result<T>
where
    T::Err: Display,
{
    cols.next()
        .ok_or_else(|| {
            anyhow::anyhow!(
                "Missing {} in {} for chromosome '{}' tile {} line {}",
                field_name,
                partials_label,
                chr,
                tile_index,
                line_number
            )
        })?
        .parse()
        .map_err(|parse_error| {
            anyhow::anyhow!(
                "Invalid {} in chromosome '{}' tile {} line {}: {}",
                field_name,
                chr,
                tile_index,
                line_number,
                parse_error
            )
        })
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum PartialsSchema {
    BedBasic,
    BedSummary,
    SizeBasic,
    SizeSummary,
}

impl PartialsSchema {
    fn partials_label(self) -> &'static str {
        match self {
            Self::BedBasic | Self::BedSummary => "partials",
            Self::SizeBasic | Self::SizeSummary => "size partials",
        }
    }

    /// Return whether this schema carries the summary-only raw moment columns.
    ///
    /// This affects only the on-disk column layout. BED vs size row identity stays separate and is
    /// still handled explicitly in `parse_row`.
    fn is_summary(self) -> bool {
        matches!(self, Self::BedSummary | Self::SizeSummary)
    }
}

/// Parsed row from one partials line, normalized into the reducer's shared in-memory shape.
///
/// BED rows keep `interval: None` because BED partial files only persist the stable `orig_idx`.
/// The BED merge engine still has to recover the final interval from `coords_by_idx`.
///
/// Size rows carry `interval: Some(full_bin_interval)` because the fixed-size partial files
/// already persist the full bin coordinates, even when one tile wrote only a clipped overlap.
///
/// Basic temp files stay narrow on disk. Their summary-only fields are filled with zeroes only
/// after parsing so the in-memory reducer field set matches the summary reducer and grouped fold.
#[derive(Debug, Clone, Copy)]
struct ParsedPartialRow {
    key: u64,
    interval: Option<Interval<u64>>,
    coverage_sum: f64,
    eligible_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
}

/// Shared streaming reader for all four partial schemas.
///
/// The IO mechanics are identical across BED/size and basic/summary partials:
/// - open plain or zstd-compressed files
/// - skip blank lines
/// - keep line numbers for diagnostics
/// - parse the current schema into one normalized row shape
struct PartialsStream {
    reader: Box<dyn BufRead>,
    line_buf: String,
    chr: String,
    line_number: u64,
    tile_index: u32,
    schema: PartialsSchema,
}

impl PartialsStream {
    fn open(path: &Path, chr: &str, tile_index: u32, schema: PartialsSchema) -> Result<Self> {
        Ok(Self {
            reader: open_text_reader(path)
                .with_context(|| format!("Opening partials file {}", path.display()))?,
            line_buf: String::new(),
            chr: chr.to_string(),
            line_number: 0,
            tile_index,
            schema,
        })
    }

    /// Read the next non-blank parsed row, or `Ok(None)` on EOF.
    fn next_row(&mut self) -> Result<Option<ParsedPartialRow>> {
        loop {
            self.line_buf.clear();
            let next_line_number = self.line_number + 1;
            let bytes_read = self.reader.read_line(&mut self.line_buf).with_context(|| {
                format!(
                    "Reading {} for chromosome '{}' tile {} line {}",
                    self.schema.partials_label(),
                    self.chr,
                    self.tile_index,
                    next_line_number
                )
            })?;
            if bytes_read == 0 {
                return Ok(None);
            }

            self.line_number = next_line_number;
            let raw = self.line_buf.trim_end_matches('\n').trim_end_matches('\r');
            if raw.is_empty() {
                continue;
            }

            return self.parse_row(raw).map(Some);
        }
    }

    /// Parse one non-blank partials line into the shared reducer row shape.
    ///
    /// The function keeps the two real schema families visible:
    /// - BED partials, keyed by stable `orig_idx`
    /// - fixed-size partials, keyed by full bin `start` and carrying full bin coordinates
    ///
    /// Within each family, basic and summary schemas share the leading columns and differ only by
    /// the summary-only tail. That is why the code splits first by BED vs size, then gates the
    /// extra raw moment columns behind `summary`.
    fn parse_row(&self, raw: &str) -> Result<ParsedPartialRow> {
        let mut cols = raw.split('\t');
        let partials_label = self.schema.partials_label();
        let summary = self.schema.is_summary();

        match self.schema {
            PartialsSchema::BedBasic | PartialsSchema::BedSummary => {
                // Expected columns in BED partials:
                //   orig_idx   sum   allowed   blacklisted
                //   [summary-only] nonzero_positions   coverage_sum_of_squares
                let key: u64 = parse_col(
                    &mut cols,
                    "orig_idx",
                    partials_label,
                    &self.chr,
                    self.tile_index,
                    self.line_number,
                )?;
                let coverage_sum: f64 = parse_col(
                    &mut cols,
                    "sum",
                    partials_label,
                    &self.chr,
                    self.tile_index,
                    self.line_number,
                )?;
                let eligible_positions: u64 = parse_col(
                    &mut cols,
                    "allowed",
                    partials_label,
                    &self.chr,
                    self.tile_index,
                    self.line_number,
                )?;
                let blacklisted_positions: u64 = parse_col(
                    &mut cols,
                    "blacklisted",
                    partials_label,
                    &self.chr,
                    self.tile_index,
                    self.line_number,
                )?;
                let (nonzero_positions, coverage_sum_of_squares) = if summary {
                    (
                        parse_col(
                            &mut cols,
                            "nonzero_positions",
                            partials_label,
                            &self.chr,
                            self.tile_index,
                            self.line_number,
                        )?,
                        parse_col(
                            &mut cols,
                            "coverage_sum_of_squares",
                            partials_label,
                            &self.chr,
                            self.tile_index,
                            self.line_number,
                        )?,
                    )
                } else {
                    // Basic temp files deliberately do not carry summary-only columns.
                    (0, 0.0)
                };

                Ok(ParsedPartialRow {
                    key,
                    interval: None,
                    coverage_sum,
                    eligible_positions,
                    blacklisted_positions,
                    nonzero_positions,
                    coverage_sum_of_squares,
                })
            }
            PartialsSchema::SizeBasic | PartialsSchema::SizeSummary => {
                // Expected columns in size partials:
                //   start   end   sum   allowed   blacklisted
                //   [summary-only] nonzero_positions   coverage_sum_of_squares
                let start: u64 = parse_col(
                    &mut cols,
                    "start",
                    partials_label,
                    &self.chr,
                    self.tile_index,
                    self.line_number,
                )?;
                let end: u64 = parse_col(
                    &mut cols,
                    "end",
                    partials_label,
                    &self.chr,
                    self.tile_index,
                    self.line_number,
                )?;
                let coverage_sum: f64 = parse_col(
                    &mut cols,
                    "sum",
                    partials_label,
                    &self.chr,
                    self.tile_index,
                    self.line_number,
                )?;
                let eligible_positions: u64 = parse_col(
                    &mut cols,
                    "allowed",
                    partials_label,
                    &self.chr,
                    self.tile_index,
                    self.line_number,
                )?;
                let blacklisted_positions: u64 = parse_col(
                    &mut cols,
                    "blacklisted",
                    partials_label,
                    &self.chr,
                    self.tile_index,
                    self.line_number,
                )?;
                let (nonzero_positions, coverage_sum_of_squares) = if summary {
                    (
                        parse_col(
                            &mut cols,
                            "nonzero_positions",
                            partials_label,
                            &self.chr,
                            self.tile_index,
                            self.line_number,
                        )?,
                        parse_col(
                            &mut cols,
                            "coverage_sum_of_squares",
                            partials_label,
                            &self.chr,
                            self.tile_index,
                            self.line_number,
                        )?,
                    )
                } else {
                    // Basic temp files deliberately do not carry summary-only columns.
                    (0, 0.0)
                };

                Ok(ParsedPartialRow {
                    key: start,
                    interval: Some(Interval::new(start, end)?),
                    coverage_sum,
                    eligible_positions,
                    blacklisted_positions,
                    nonzero_positions,
                    coverage_sum_of_squares,
                })
            }
        }
    }
}

/// Shared additive accumulator for one BED row keyed by `orig_idx` or one fixed-size bin keyed by
/// full bin `start`.
///
/// The summary-only fields are intentionally present even in basic-mode reduction. They remain
/// zero there, which keeps the reducer output shape aligned with grouped folding and summary
/// writers without widening the temp-file schema.
#[derive(Debug, Clone, Copy, Default)]
struct AggregateAccum {
    coverage_sum: f64,
    eligible_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
    seen_contributions: u32,
}

impl AggregateAccum {
    fn add_row(&mut self, row: ParsedPartialRow) {
        self.coverage_sum += row.coverage_sum;
        self.eligible_positions += row.eligible_positions;
        self.blacklisted_positions += row.blacklisted_positions;
        self.nonzero_positions += row.nonzero_positions;
        self.coverage_sum_of_squares += row.coverage_sum_of_squares;
        self.seen_contributions += 1;
    }

    fn into_reduced_row(self, idx: u64, interval: Interval<u64>) -> ReducedAggregateRow {
        ReducedAggregateRow {
            idx,
            interval,
            coverage_sum: self.coverage_sum,
            eligible_positions: self.eligible_positions,
            blacklisted_positions: self.blacklisted_positions,
            nonzero_positions: self.nonzero_positions,
            coverage_sum_of_squares: self.coverage_sum_of_squares,
        }
    }
}

/// Final raw aggregate row after cross-tile reduction.
///
/// This is the contract between reducer code and writer code:
/// - reducers stop at exact additive raw values
/// - writers derive `average`, `total`, variance, SD, CV, and other presentation-layer fields
///
/// BED and fixed-size reducers both produce this same row shape, even though they recover the
/// interval differently.
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

/// Build the BED interval lookup keyed by stable original row index.
///
/// BED partial rows intentionally keep only `orig_idx` on disk, so interval recovery has to stay
/// explicit in the BED reducer rather than being hidden inside the parsed row shape.
fn build_bed_coords_by_idx(
    chr: &str,
    windows_chr: &[IndexedInterval<u64>],
) -> Result<FxHashMap<u64, Interval<u64>>> {
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
    Ok(coords_by_idx)
}

/// Count the expected number of tile contributions for each reducer row key.
///
/// Cross-index sidecars list only boundary-crossing rows. Any key missing from every sidecar still
/// expects exactly one contribution, but this helper deliberately returns only the explicit counts.
/// The reducer applies the `default = 1` policy at the read point so that rule stays obvious.
///
/// Key meaning depends on reducer family:
/// - BED reducers use stable `orig_idx`
/// - fixed-size reducers use full bin `start`
fn load_expected_contributions(
    files_by_tile: &FxHashMap<u32, TileFiles>,
) -> Result<FxHashMap<u64, u32>> {
    let mut expected_contributions: FxHashMap<u64, u32> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    for tile_files in files_by_tile.values() {
        let Some(cross_path) = &tile_files.cross_index_path else {
            continue;
        };

        let mut reader = open_text_reader(cross_path)
            .with_context(|| format!("Opening cross-index {}", cross_path.display()))?;
        let mut line = String::new();
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }

            let raw = line.trim_end_matches('\n').trim_end_matches('\r');
            if raw.is_empty() {
                continue;
            }

            let key: u64 = raw
                .parse()
                .with_context(|| format!("Invalid cross-index key in {}", cross_path.display()))?;
            *expected_contributions.entry(key).or_insert(0) += 1;
        }
    }

    Ok(expected_contributions)
}

#[inline]
/// Return how many tile rows must be seen before one reduced row is complete.
///
/// Keys missing from all cross-index sidecars default to one contribution because they stayed
/// inside one tile core and therefore produced exactly one partial row.
fn expected_contribution_count(expected_contributions: &FxHashMap<u64, u32>, key: u64) -> u32 {
    *expected_contributions.get(&key).unwrap_or(&1)
}

/// Open one partials stream per tile and seed the merge heap with the current visible row.
///
/// The returned heap is keyed only by the next visible row key. It is not a completeness proof.
/// Final correctness still comes from the stable row identity plus the expected
/// contribution counts loaded from the cross-index sidecars.
fn open_partials_streams(
    chr: &str,
    files_by_tile: &FxHashMap<u32, TileFiles>,
    schema: PartialsSchema,
) -> Result<(
    Vec<PartialsStream>,
    Vec<Option<ParsedPartialRow>>,
    StreamHeap,
)> {
    // Open one stream per partials file and keep only that stream's current row in memory.
    // The heap chooses the smallest visible row key across those rows.
    //
    // `BinaryHeap` is a max-heap, so `Reverse((key, stream_id))` turns it into a min-heap.
    //
    // The heap does not prove that a BED row or fixed-size bin is complete. It only picks the
    // next stream to read from. Final correctness still comes from the stable row identity plus
    // the expected contribution count loaded from the cross-index sidecars.
    let mut streams: Vec<PartialsStream> = Vec::new();
    let mut current_rows: Vec<Option<ParsedPartialRow>> = Vec::new();
    let mut heap: StreamHeap = BinaryHeap::new();

    for (tile_index, tile_files) in files_by_tile {
        let Some(partials_path) = &tile_files.partials_path else {
            continue;
        };

        let mut stream = PartialsStream::open(partials_path, chr, *tile_index, schema)?;
        if let Some(row) = stream.next_row()? {
            let stream_id = streams.len();
            streams.push(stream);
            current_rows.push(Some(row));
            heap.push(Reverse((row.key, stream_id)));
        }
    }

    Ok((streams, current_rows, heap))
}

fn files_by_tile_from_outputs(
    tile_outputs: &[TileAggregateTempFiles],
) -> Result<FxHashMap<u32, TileFiles>> {
    let mut files_by_tile: FxHashMap<u32, TileFiles> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    for output in tile_outputs {
        let previous = files_by_tile.insert(
            output.tile_index,
            TileFiles {
                partials_path: Some(output.partials_path.clone()),
                cross_index_path: output.cross_index_path.clone(),
            },
        );
        anyhow::ensure!(
            previous.is_none(),
            "duplicate tile index {} in explicit aggregate tile outputs",
            output.tile_index
        );
    }

    Ok(files_by_tile)
}

/// Advance one stream after its current row was consumed and push the next visible row into the heap.
///
/// Keeping this as a tiny helper makes the merge loops read as "consume current row, maybe emit,
/// then advance that same stream" instead of repeating the stream bookkeeping in every reducer.
fn push_next_row_for_stream(
    streams: &mut [PartialsStream],
    current_rows: &mut [Option<ParsedPartialRow>],
    heap: &mut StreamHeap,
    stream_id: usize,
) -> Result<()> {
    if let Some(next_row) = streams[stream_id].next_row()? {
        current_rows[stream_id] = Some(next_row);
        heap.push(Reverse((next_row.key, stream_id)));
    }
    Ok(())
}

/// Reduce BED partial rows for one chromosome into complete raw aggregate rows.
///
/// The reducer consumes the exact tile paths returned by tile processing. It does not discover
/// files from a temp directory, so stale or decoy files cannot enter the reduction.
///
/// The row identity rule is unchanged by explicit paths. BED reduction always groups by stable
/// `orig_idx`, then recovers the final interval from `windows_chr` after every expected tile
/// contribution has arrived.
///
/// About ordering:
/// - the merge heap reads the smallest visible `orig_idx` from the open streams
/// - final emission waits for the expected number of contributions for that `orig_idx`
/// - correctness does not depend on every partials stream being globally sorted
/// - ordinary BED windows keep their original file indices, so callers must not assume final
///   output is coordinate-sorted unless the windows were explicitly reindexed upstream
///
/// Cross-index logic:
/// - sidecars list only boundary-crossing rows
/// - windows fully contained in one tile are absent from all sidecars, so the reducer expects one
///   contribution
/// - boundary windows appear in each crossed tile's sidecar, and that sidecar count is the
///   expected number of partial rows
///
/// This engine intentionally stays separate from the fixed-size engine. Both engines share the
/// same high-level merge pattern, but BED reduction has one BED-specific responsibility that size
/// reduction does not: recover the final interval from `orig_idx` after reduction.
///
/// Parameters
/// ----------
/// - `chr`:
///     Chromosome label used in diagnostics.
/// - `tile_outputs`:
///     Returned partial and optional cross-index paths for this chromosome.
/// - `windows_chr`:
///     BED windows with the same `orig_idx` identities written in the partial rows.
/// - `summary`:
///     Selects the on-disk partial schema only. It does not change row identity.
/// - `on_row`:
///     Callback receiving exact additive raw rows.
pub(crate) fn reduce_bed_rows(
    chr: &str,
    tile_outputs: &[TileAggregateTempFiles],
    windows_chr: &[IndexedInterval<u64>],
    summary: bool,
    mut on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()> {
    let coords_by_idx = build_bed_coords_by_idx(chr, windows_chr)?;
    let files_by_tile = files_by_tile_from_outputs(tile_outputs)?;
    let expected_contributions = load_expected_contributions(&files_by_tile)?;
    let schema = if summary {
        PartialsSchema::BedSummary
    } else {
        PartialsSchema::BedBasic
    };
    let (mut streams, mut current_rows, mut heap) =
        open_partials_streams(chr, &files_by_tile, schema)?;

    // One accumulator per BED row that has started receiving tile contributions but is not yet
    // complete.
    let mut accum_by_idx: FxHashMap<u64, AggregateAccum> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    while let Some(Reverse((_, stream_id))) = heap.pop() {
        let row = current_rows[stream_id].take().ok_or_else(|| {
            anyhow::anyhow!(
                "Reducer heap and current_row fell out of sync for chromosome '{}' stream {}",
                chr,
                stream_id
            )
        })?;
        let key = row.key;

        // Add this tile's contribution into the running totals for the BED row keyed by stable
        // `orig_idx`.
        let entry = accum_by_idx.entry(key).or_default();
        entry.add_row(row);

        // Write the row only once the cross-index says every tile contribution has arrived.
        // Windows that never appear in any cross-index file default to one expected contribution.
        if entry.seen_contributions == expected_contribution_count(&expected_contributions, key) {
            let done = accum_by_idx.remove(&key).ok_or_else(|| {
                anyhow::anyhow!(
                    "Reducer lost accumulated window state for chromosome '{}' orig_idx {}",
                    chr,
                    key
                )
            })?;
            let interval = *coords_by_idx.get(&key).ok_or_else(|| {
                anyhow::anyhow!(
                    "Reducer is missing interval coordinates for chromosome '{}' orig_idx {}",
                    chr,
                    key
                )
            })?;
            on_row(done.into_reduced_row(key, interval))?;
        }

        // Advance only the stream we just consumed from, then push its next visible row back into
        // the heap. Memory stays bounded by the number of open tile files.
        push_next_row_for_stream(&mut streams, &mut current_rows, &mut heap, stream_id)?;
    }

    anyhow::ensure!(
        accum_by_idx.is_empty(),
        "Incomplete windows remain for {}: {}",
        chr,
        accum_by_idx.len()
    );

    Ok(())
}

/* By-size reducer (when windows don't align) */

/// Reduce fixed-size partial rows for one chromosome into complete raw bin rows.
///
/// The reducer consumes the exact tile paths returned by tile processing. It does not discover
/// files from a temp directory.
///
/// The row identity rule stays fixed. Size reduction groups by the logical full-bin `start`, not
/// by any clipped tile-local overlap. Partial rows must therefore carry the full `bin_start` and
/// `bin_end`. Changing those bounds to clipped pieces would break cross-tile merging.
///
/// Cross-index sidecars count how many tiles contribute to each full bin start. If a bin is not
/// listed in any cross-index file, the reducer expects exactly one contribution. Missing
/// cross-index files are valid for aligned single-contribution partials.
///
/// Ordering comes from the same streaming heap as the BED reducer. `BinaryHeap` is a max-heap, so
/// the reader heap stores `Reverse((start, stream_id))` to read the smallest visible logical bin
/// start first.
///
/// The final fixed-size bin can extend past chromosome end. The reducer clips that interval after
/// all contributions have been combined so downstream writers see the true genomic span.
///
/// This engine intentionally stays separate from the BED engine. Both engines share the same
/// merge rhythm, but fixed-size reduction has one size-specific responsibility that BED reduction
/// does not: clip the final bin after reduction using the true chromosome end.
///
/// Parameters
/// ----------
/// - `chr`:
///     Chromosome label used in diagnostics.
/// - `tile_outputs`:
///     Returned partial and optional cross-index paths for this chromosome.
/// - `chrom_len`:
///     True chromosome length used to clip the final bin.
/// - `summary`:
///     Selects the on-disk partial schema only. It does not change row identity.
/// - `on_row`:
///     Callback receiving exact additive raw rows.
pub(crate) fn reduce_size_rows(
    chr: &str,
    tile_outputs: &[TileAggregateTempFiles],
    chrom_len: u64,
    summary: bool,
    mut on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()> {
    let files_by_tile = files_by_tile_from_outputs(tile_outputs)?;
    let expected_contributions = load_expected_contributions(&files_by_tile)?;
    let schema = if summary {
        PartialsSchema::SizeSummary
    } else {
        PartialsSchema::SizeBasic
    };
    let (mut streams, mut current_rows, mut heap) =
        open_partials_streams(chr, &files_by_tile, schema)?;

    // One accumulator per full bin start that has begun receiving tile contributions but is not
    // complete yet.
    let mut accum_by_start: FxHashMap<u64, AggregateAccum> =
        FxHashMap::with_hasher(FxBuildHasher::default());

    while let Some(Reverse((_, stream_id))) = heap.pop() {
        let row = current_rows[stream_id].take().ok_or_else(|| {
            anyhow::anyhow!(
                "Reducer heap and current_row fell out of sync for chromosome '{}' stream {}",
                chr,
                stream_id
            )
        })?;
        let key = row.key;
        let full_interval = row.interval.ok_or_else(|| {
            anyhow::anyhow!(
                "Size reducer row for chromosome '{}' start {} is missing its full bin interval",
                chr,
                key
            )
        })?;
        debug_assert_eq!(full_interval.start(), key);

        // Add this tile contribution into the running totals for the fixed-size bin keyed by full
        // bin `start`.
        let entry = accum_by_start.entry(key).or_default();
        entry.add_row(row);

        if entry.seen_contributions == expected_contribution_count(&expected_contributions, key) {
            let done = accum_by_start.remove(&key).ok_or_else(|| {
                anyhow::anyhow!(
                    "Reducer lost accumulated size-bin state for chromosome '{}' start {}",
                    chr,
                    key
                )
            })?;
            // The final fixed-size bin can overrun the chromosome end. Clip it here so downstream
            // writers see the true genomic interval length for the last row.
            let clipped_interval = full_interval.clip_upper(chrom_len).ok_or_else(|| {
                anyhow::anyhow!(
                    "Size reducer produced an empty clipped interval for chromosome '{}' start {} with chrom_len {}",
                    chr,
                    key,
                    chrom_len
                )
            })?;
            on_row(done.into_reduced_row(key, clipped_interval))?;
        }

        push_next_row_for_stream(&mut streams, &mut current_rows, &mut heap, stream_id)?;
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
    pub partials_path: Option<PathBuf>,
    pub cross_index_path: Option<PathBuf>,
}

#[cfg(test)]
mod tests {
    include!("reducer_tests.rs");
}
