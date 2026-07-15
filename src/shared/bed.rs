use crate::shared::interval::IndexedInterval;
#[cfg(feature = "cmd_fcoverage")]
use crate::shared::interval::Interval;
use crate::shared::interval::{ScoredInterval, Span};
#[cfg(feature = "cmd_fcoverage")]
use crate::shared::interval::{TouchingMergePolicy, push_merged_interval};
use crate::shared::io::{open_text_reader, open_text_reader_in_background};
use anyhow::{Context, Result, bail, ensure};
use fxhash::{FxHashMap, FxHashSet};
#[cfg(feature = "cmd_fcoverage")]
use rayon::prelude::*;
#[cfg(feature = "cmd_fcoverage")]
use std::fs::File;
#[cfg(feature = "cmd_fcoverage")]
use std::io::{BufWriter, Write};
use std::{io::BufRead, path::Path};

// Small bounded prefix used to reject obvious non-BED input before whitelist skipping
const BED_FORMAT_SNIFF_ROWS: usize = 16;

/// Load windows from a BED file into a per-chromosome map.
///
/// The original window index is added. Any valid window increases the index,
/// even if they are filtered from the returned windows.
///
/// Parameters
/// ----------
///  - bed: Path to BED file.
///  - chromosomes: Names of chromosomes to include in output,
///    even when not present in the BED file.
///  - filter_fn: Function for deciding whether to include
///    an interval. Should take in the `chr,start,end` values
///    and return `true` (keep) or `false` (discard).
///  - exp_num_windows: Optional number of expected windows
///    in the BED file before filtering. Returns an error
///    if the incorrect number of windows are observed.
///  - read_in_background: Run file reading and decompression on a background thread while the
///    calling thread parses rows.
///
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, original window index)'.
pub(crate) fn load_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    filter_fn: Option<&dyn Fn(&str, u64, u64) -> bool>,
    exp_num_windows: Option<u64>,
    read_in_background: bool,
) -> Result<FxHashMap<String, Windows>> {
    let mut reader = if read_in_background {
        open_text_reader_in_background(bed.as_ref())?
    } else {
        open_text_reader(bed.as_ref())?
    };

    // Optional whitelist of chromosomes
    let mut windows_by_chromosome: FxHashMap<String, Vec<(u64, u64, u64)>> = FxHashMap::default();
    let allowed_chromosomes: Option<FxHashSet<&str>> = chromosomes.map(|chr_list| {
        let mut allowed = FxHashSet::with_capacity_and_hasher(chr_list.len(), Default::default());
        for chr in chr_list {
            allowed.insert(chr.as_str());
            windows_by_chromosome.entry(chr.clone()).or_default();
        }
        allowed
    });

    // Reuse a single buffer for all lines
    let mut buf = String::new();
    let mut lineno: usize = 0;
    let mut orig_win_idx: u64 = 0; // Counter for all *valid* windows whether filtered out or not
    let mut sniffed_rows = 0usize;

    loop {
        buf.clear();
        let bytes_read = reader.read_line(&mut buf)?;
        if bytes_read == 0 {
            break;
        }
        lineno += 1;

        // Fast skips
        if buf.as_bytes().first().is_some_and(|b| *b == b'#') {
            continue;
        }
        let line = buf.trim_end_matches(['\n', '\r']);

        if line.is_empty() {
            continue;
        }

        // Skip UCSC header directives in BED files
        let trimmed_line_start = line.trim_start();
        if trimmed_line_start.starts_with("track") || trimmed_line_start.starts_with("browser") {
            continue;
        }

        if sniffed_rows < BED_FORMAT_SNIFF_ROWS {
            sniff_bed3_line(line, lineno)?;
            sniffed_rows += 1;
        }

        // Strict parse of first 3 BED columns without allocating a Vec
        let mut fields = line.split_ascii_whitespace();

        // Get the three coordinate fields to ensure at least three values are present
        // We need to know if the line represents a valid window (whether we want it or not)

        let chr = match fields.next() {
            Some(s) => s,
            None => continue,
        };

        // NOTE: "!Allowed" != "!Valid"
        // A chromosome name can be valid with the line representing a genomic window
        // but disallowed due to not being in the chromosome whitelist
        // but any valid window is considered for the original index incrementing
        let is_allowed_chrom = if let Some(allowed_chroms) = &allowed_chromosomes {
            allowed_chroms.contains(chr)
        } else {
            // If we don't have a whitelist, we must assume it's an allowed chromosome name (whether valid or not)
            true
        };

        let start_str = match fields.next() {
            Some(s) => s,
            None => {
                if is_allowed_chrom {
                    // If initial value is known to be a valid chromosome (when it's in the whitelist)
                    // or we don't know whether it's valid (we have no whitelist)
                    // we expect the start field to exist
                    bail!(
                        "BED parse error at line {}: missing start for chromosome: '{}'",
                        lineno,
                        chr
                    );
                }
                continue;
            }
        };

        let end_str = match fields.next() {
            Some(s) => s,
            None => {
                if is_allowed_chrom {
                    // If initial value is known to be a valid chromosome (when it's in the whitelist)
                    // or we don't know whether it's valid (we have no whitelist)
                    // we expect the end field to exist
                    bail!(
                        "BED parse error at line {}: missing end for chromosome: '{}'",
                        lineno,
                        chr
                    );
                }
                continue;
            }
        };

        // Skip if not an allowed chromosome
        // Increment iff window is assumed valid
        if !is_allowed_chrom {
            // Only increment original window index coordinates make a valid window
            if start_end_are_valid_coordinates(start_str, end_str).is_some() {
                orig_win_idx += 1;
            }
            continue;
        }

        // Parse coordinates
        let start: u64 = start_str.parse().with_context(|| {
            format!(
                "BED parse error at line {}: invalid start '{}'",
                lineno, start_str
            )
        })?;
        let end: u64 = end_str.parse().with_context(|| {
            format!(
                "BED parse error at line {}: invalid end '{}'",
                lineno, end_str
            )
        })?;

        ensure!(
            end > start,
            "BED parse error at line {}: end ({}) must be greater than start ({})",
            lineno,
            end,
            start
        );

        // At this point, we assume the window to be a valid window
        // and increment the counter for the original index
        let current_orig_win_idx = orig_win_idx;
        orig_win_idx += 1;

        if let Some(filterer) = filter_fn
            && !filterer(chr, start, end)
        {
            continue;
        }

        windows_by_chromosome
            .entry(chr.to_string())
            .or_default()
            .push((start, end, current_orig_win_idx));
    }

    if let Some(expected_num_windows) = exp_num_windows {
        ensure!(
            expected_num_windows == orig_win_idx,
            "the BED file did not contain the correct number of windows: obs: {} != exp: {}",
            orig_win_idx,
            expected_num_windows
        );
    }

    let windows_mapping: FxHashMap<String, Windows> = windows_by_chromosome
        .into_iter()
        .map(|(chr, windows)| Ok((chr, Windows::from_tuples(&windows)?)))
        .collect::<crate::Result<_>>()?;

    Ok(windows_mapping)
}

/// Returns `Option` indicating whether start and end strings are parseable coordinates
/// making up a valid window (only checks end > start).
fn start_end_are_valid_coordinates(start_str: &str, end_str: &str) -> Option<()> {
    let start = start_str.parse::<u64>().ok()?;
    let end = end_str.parse::<u64>().ok()?;
    if end <= start {
        return None;
    }
    Some(())
}

fn sniff_bed3_line(line: &str, lineno: usize) -> Result<()> {
    ensure!(
        !line.as_bytes().contains(&0),
        "BED parse error at line {}: input appears to be binary, not BED text",
        lineno
    );

    let mut fields = line.split_ascii_whitespace();
    let chr = fields
        .next()
        .with_context(|| format!("BED parse error at line {}: missing chromosome", lineno))?;
    let start_str = fields.next().with_context(|| {
        format!(
            "BED parse error at line {}: missing start for chromosome: '{}'",
            lineno, chr
        )
    })?;
    let end_str = fields.next().with_context(|| {
        format!(
            "BED parse error at line {}: missing end for chromosome: '{}'",
            lineno, chr
        )
    })?;

    let start: u64 = start_str.parse().with_context(|| {
        format!(
            "BED parse error at line {}: invalid start '{}'",
            lineno, start_str
        )
    })?;
    let end: u64 = end_str.parse().with_context(|| {
        format!(
            "BED parse error at line {}: invalid end '{}'",
            lineno, end_str
        )
    })?;

    ensure!(
        end > start,
        "BED parse error at line {}: end ({}) must be greater than start ({})",
        lineno,
        end,
        start
    );
    Ok(())
}

/// Owned collection of half-open windows with a cached genomic span.
///
/// Invariants
/// ----------
/// - `windows` should be sorted by start (ascending order).
/// - Coordinates are half-open: `[start, end)`.
#[derive(Debug, Clone)]
pub(crate) struct Windows {
    pub(crate) windows: Vec<IndexedInterval<u64>>,
    /// Cached outer envelope across all windows.
    #[allow(dead_code)]
    span: Span<i64>,
}

impl Windows {
    /// Construct from any window list (may be unsorted/overlapping).
    /// Ensures start- and end-sorted order (does not retain initial order)
    /// and computes span as `min(start)` .. `max(end)`.
    pub(crate) fn new(mut windows: Vec<IndexedInterval<u64>>) -> Self {
        windows.sort_unstable_by_key(|window| (window.start(), window.end()));
        Windows::from_sorted(windows)
    }

    /// Construct from raw `(start, end, idx)` tuples.
    ///
    /// Use this when a loader or tiny fixture still naturally produces tuples.
    /// Prefer `new` and `from_sorted` when the windows are already checked.
    pub(crate) fn from_tuples(windows: &[(u64, u64, u64)]) -> crate::Result<Self> {
        Ok(Self::new(IndexedInterval::from_tuples(windows)?))
    }

    /// Construct from a list you guarantee is already sorted by start (non-decreasing).
    /// Computes span as `min(start)` .. `max(end)` (robust to irregular ends).
    pub(crate) fn from_sorted(indexed_windows: Vec<IndexedInterval<u64>>) -> Self {
        debug_assert!(
            is_sorted_by_start_indexed(&indexed_windows),
            "windows must be start-sorted"
        );
        let span = if indexed_windows.is_empty() {
            Span::from_ordered(0, 0)
        } else {
            let min_start = indexed_windows[0].start() as i64;
            let max_end = indexed_windows
                .iter()
                .map(|window| window.end())
                .max()
                .unwrap() as i64;
            Span::from_ordered(min_start, max_end)
        };
        Self {
            windows: indexed_windows,
            span,
        }
    }

    /// Number of windows.
    #[allow(
        dead_code,
        reason = "feature-limited builds compile shared BED window containers without every accessor being used"
    )]
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.windows.len()
    }

    /// True if there are no windows.
    #[allow(
        dead_code,
        reason = "feature-limited builds compile shared BED window containers without every accessor being used"
    )]
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.windows.is_empty()
    }

    /// Borrow the underlying windows.
    #[allow(
        dead_code,
        reason = "feature-limited builds compile shared BED window containers without every accessor being used"
    )]
    #[inline]
    pub(crate) fn as_slice(&self) -> &[IndexedInterval<u64>] {
        &self.windows
    }

    /// Consume and return the inner vector.
    #[inline]
    pub(crate) fn into_inner(self) -> Vec<IndexedInterval<u64>> {
        self.windows
    }

    /// Span start (inclusive).
    /// This is the most-left coordinate covered by any of the windows.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn span_start(&self) -> i64 {
        self.span.start()
    }

    /// Span end (exclusive).
    /// This is the most-right coordinate covered by any of the windows.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn span_end(&self) -> i64 {
        self.span.end()
    }

    /// Return the collection span.
    ///
    /// This is the outer envelope of the stored windows: the smallest start and
    /// largest end seen in the collection. It returns `Span<i64>` rather than
    /// `Interval<i64>` because empty collections use the empty span `[0, 0)`.
    ///
    /// There are no guarantees that every position inside this span is covered
    /// by a window.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn span(&self) -> Span<i64> {
        self.span
    }

    /// Merge/flatten touching or overlapping windows and reindex them sequentially **in-place**.
    ///
    /// Summary
    /// -------
    /// Consumes `self`, merges `[start, end)` windows that overlap or touch, and reassigns
    /// new indices starting at `start_idx`. Reuses the original allocation to avoid
    /// peak-memory spikes when window lists are huge.
    ///
    /// Parameters
    /// ----------
    /// - start_idx:
    ///   First index to assign to the merged interval series.
    ///
    /// Returns
    /// -------
    /// - (merged, next_start_idx):
    ///     - `merged`: New `Windows` with merged, start-sorted `(start, end, new_idx)` tuples.
    ///     - `next_start_idx`: `start_idx + merged.len()`; pass this to the next chromosome.
    #[cfg(flattens_bed_windows)]
    pub(crate) fn into_flattened_reindexed(self, start_idx: u64) -> (Windows, u64) {
        let mut merged_windows = self.windows; // Take ownership; reuse allocation
        if merged_windows.is_empty() {
            return (
                Windows {
                    windows: merged_windows,
                    span: Span::from_ordered(0, 0),
                },
                start_idx,
            );
        }

        debug_assert!(
            is_sorted_by_start_indexed(&merged_windows),
            "windows must be start-sorted"
        );

        // In-place compaction with two indices: read cursor and write cursor
        let mut write_index: usize = 0;
        let mut current_start = merged_windows[0].start();
        let mut current_end = merged_windows[0].end();

        for read_index in 1..merged_windows.len() {
            let next_start = merged_windows[read_index].start();
            let next_end = merged_windows[read_index].end();
            if next_start <= current_end {
                if next_end > current_end {
                    current_end = next_end;
                }
            } else {
                // Write merged block at the current write position with a new index
                merged_windows[write_index] = IndexedInterval::new(
                    current_start,
                    current_end,
                    start_idx + write_index as u64,
                )
                .expect("merged windows must remain valid non-empty intervals");
                write_index += 1;
                current_start = next_start;
                current_end = next_end;
            }
        }
        // Write the final block
        merged_windows[write_index] =
            IndexedInterval::new(current_start, current_end, start_idx + write_index as u64)
                .expect("merged windows must remain valid non-empty intervals");
        write_index += 1;

        // Shrink to the number of merged intervals
        merged_windows.truncate(write_index);

        // Since the windows are start-sorted and merged, first and last bound the span
        let span_start = merged_windows
            .first()
            .map(|window| window.start() as i64)
            .unwrap_or(0);
        let span_end = merged_windows
            .last()
            .map(|window| window.end() as i64)
            .unwrap_or(0);

        let next_idx = start_idx + write_index as u64;

        (
            Windows {
                windows: merged_windows,
                span: Span::from_ordered(span_start, span_end),
            },
            next_idx,
        )
    }
}

impl AsRef<[IndexedInterval<u64>]> for Windows {
    fn as_ref(&self) -> &[IndexedInterval<u64>] {
        self.as_slice()
    }
}

fn is_sorted_by_start_indexed(ws: &[IndexedInterval<u64>]) -> bool {
    ws.windows(2)
        .all(|window_pair| window_pair[0].start() <= window_pair[1].start())
}

#[allow(dead_code)]
#[inline]
fn is_sorted_by_start_with_scores(ws: &[ScoredInterval<u64>]) -> bool {
    ws.windows(2)
        .all(|window_pair| window_pair[0].start() <= window_pair[1].start())
}

/* GROUPED bed files */

#[cfg(loads_grouped_bed)]
const GROUPED_BED_STRAND_SAMPLE_ROWS: usize = 200;

/// Load *grouped* windows from a BED file into a per-chromosome map.
///
/// Parameters
/// ----------
///  - bed: Path to BED file with group names in the fourth column.
///  - chromosomes: Names of chromosomes to include in output,
///    even when not present in the BED file.
///  - read_strands: Detect and read interval strandedness from the BED fields. For files with six
///    or more columns, only column 6 is read as strand. Column 5 is accepted only for five-column
///    grouped files.
///  - filter_fn: Function for deciding whether to include
///    an interval. Should take in the `chr,start,end` values
///    and return `true` (keep) or `false` (discard).
///  - exp_num_windows: Optional number of expected windows
///    in the BED file before filtering. Returns an error
///    if the incorrect number of windows are observed.
///  - read_in_background: Run file reading and decompression on a background thread while the
///    calling thread parses rows.
///
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, group index)'.
///
///  - Mapping of 'group index -> group name'.
///
///  - Optional strand-detection metadata when `read_strands` is enabled.
#[cfg(loads_grouped_bed)]
pub(crate) fn load_grouped_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    read_strands: bool,
    filter_fn: Option<&dyn Fn(&str, u64, u64) -> bool>,
    exp_num_windows: Option<u64>,
    read_in_background: bool,
) -> Result<(
    FxHashMap<String, GroupedWindows>,
    FxHashMap<u64, String>,
    Option<GroupedBedStrandDetection>,
)> {
    let bed_path = bed.as_ref();
    let strand_detection = if read_strands {
        Some(detect_grouped_bed_strand_column(bed_path)?)
    } else {
        None
    };
    let detected_strand_column = strand_detection
        .as_ref()
        .and_then(|detection| detection.column);
    let mut reader = if read_in_background {
        open_text_reader_in_background(bed_path)?
    } else {
        open_text_reader(bed_path)?
    };

    // Initialize maps
    let mut grouped_windows_by_chromosome: FxHashMap<String, Vec<IndexedInterval<u64>>> =
        FxHashMap::default();
    let mut strand_mapping: Option<FxHashMap<String, Vec<Strand>>> =
        detected_strand_column.is_some().then(FxHashMap::default);

    // Optional whitelist of chromosomes
    let allowed_chromosomes: Option<FxHashSet<&str>> = chromosomes.map(|chr_list| {
        let mut allowed = FxHashSet::with_capacity_and_hasher(chr_list.len(), Default::default());
        for chr in chr_list {
            allowed.insert(chr.as_str());
            grouped_windows_by_chromosome
                .entry(chr.clone())
                .or_default();
            if let Some(strands_by_chr) = strand_mapping.as_mut() {
                strands_by_chr.entry(chr.clone()).or_default();
            }
        }
        allowed
    });

    // Enumeration of group names
    let mut group_name_to_idx: FxHashMap<String, u64> = FxHashMap::default();
    let mut next_group_idx: u64 = 0;

    // Reuse a single buffer for all lines
    let mut buf = String::new();
    let mut lineno: usize = 0;
    let mut orig_win_idx: u64 = 0; // Counter for all *valid* windows whether filtered out or not
    let mut sniffed_rows = 0usize;

    loop {
        buf.clear();
        let bytes_read = reader.read_line(&mut buf)?;
        if bytes_read == 0 {
            break;
        }
        lineno += 1;

        let line = buf.trim_end_matches(['\n', '\r']);

        // Fast skips
        if should_skip_bed_line(line) {
            continue;
        }

        if sniffed_rows < BED_FORMAT_SNIFF_ROWS {
            sniff_bed3_line(line, lineno)?;
            sniffed_rows += 1;
        }

        // Strict parse of first 3 BED columns without allocating a Vec
        let mut fields = line.split_ascii_whitespace();

        let chr = match fields.next() {
            Some(s) => s,
            None => continue,
        };

        // NOTE: "!Allowed" != "!Valid"
        // A chromosome name can be valid with the line representing a genomic window
        // but disallowed due to not being in the chromosome whitelist
        // but any valid window is considered for the original index incrementing
        let is_allowed_chrom = if let Some(allowed_chroms) = &allowed_chromosomes {
            allowed_chroms.contains(chr)
        } else {
            // If we don't have a whitelist, we must assume it's an allowed chromosome name (whether valid or not)
            true
        };

        let start_str = match fields.next() {
            Some(s) => s,
            None => {
                if is_allowed_chrom {
                    // If initial value is known to be a valid chromosome (when it's in the whitelist)
                    // or we don't know whether it's valid (we have no whitelist)
                    // we expect the start field to exist
                    bail!(
                        "BED parse error at line {}: missing start for chromosome: '{}'",
                        lineno,
                        chr
                    );
                }
                continue;
            }
        };

        let end_str = match fields.next() {
            Some(s) => s,
            None => {
                if is_allowed_chrom {
                    // If initial value is known to be a valid chromosome (when it's in the whitelist)
                    // or we don't know whether it's valid (we have no whitelist)
                    // we expect the end field to exist
                    bail!(
                        "BED parse error at line {}: missing end for chromosome: '{}'",
                        lineno,
                        chr
                    );
                }
                continue;
            }
        };

        // Skip if not an allowed chromosome
        // Increment iff window is assumed valid
        if !is_allowed_chrom {
            // Only increment original window index coordinates make a valid window
            if start_end_are_valid_coordinates(start_str, end_str).is_some() {
                orig_win_idx += 1;
            }
            continue;
        }

        // Get group ID from fourth column
        let group = fields
            .next()
            .with_context(|| format!("BED parse error at line {}: missing group name", lineno))?;

        // Get group idx (enumerate and insert if first occurence)
        // We use this if/else approach only allocate a String once per unique group name
        let group_idx = if let Some(&existing_group_idx) = group_name_to_idx.get(group) {
            existing_group_idx
        } else {
            let id = next_group_idx;
            next_group_idx += 1;
            group_name_to_idx.insert(group.to_owned(), id); // Only allocate here
            id
        };

        let start: u64 = start_str.parse().with_context(|| {
            format!(
                "BED parse error at line {}: invalid start '{}'",
                lineno, start_str
            )
        })?;
        let end: u64 = end_str.parse().with_context(|| {
            format!(
                "BED parse error at line {}: invalid end '{}'",
                lineno, end_str
            )
        })?;

        ensure!(
            end > start,
            "BED parse error at line {}: end ({}) must be greater than start ({})",
            lineno,
            end,
            start
        );

        // At this point, we assume the window to be a valid window
        // and increment the counter for the original index
        orig_win_idx += 1;

        // Apply passed filtering function
        if let Some(filterer) = filter_fn
            && !filterer(chr, start, end)
        {
            continue;
        }

        let strand = match detected_strand_column {
            None => Strand::Unstranded,
            Some(GroupedBedStrandColumn::Column5) => {
                let column5 = fields.next();
                let column6 = fields.next();
                ensure!(
                    column6.is_none(),
                    "BED parse error at line {}: inconsistent grouped BED column count. Strands were detected in column 5 from 5-column rows, but this row has 6 or more columns",
                    lineno
                );
                parse_grouped_bed_strand_value(column5, lineno, 5)?
            }
            Some(GroupedBedStrandColumn::Column6) => {
                let _column5 = fields.next();
                let column6 = fields.next();
                parse_grouped_bed_strand_value(column6, lineno, 6)?
            }
        };

        let checked_window = IndexedInterval::new(start, end, group_idx).with_context(|| {
            format!(
                "BED parse error at line {}: invalid interval [{start},{end})",
                lineno
            )
        })?;

        grouped_windows_by_chromosome
            .entry(chr.to_string())
            .or_default()
            .push(checked_window);
        if let Some(strands_by_chr) = strand_mapping.as_mut() {
            strands_by_chr
                .entry(chr.to_string())
                .or_default()
                .push(strand);
        }
    }

    if let Some(expected_num_windows) = exp_num_windows {
        ensure!(
            expected_num_windows == orig_win_idx,
            "the BED file did not contain the correct number of windows: obs: {} != exp: {}",
            orig_win_idx,
            expected_num_windows
        );
    }
    ensure!(
        orig_win_idx > 0,
        "grouped BED file did not contain any interval rows. Empty BED files are not valid input"
    );

    // Convert parsed checked intervals into grouped windows.
    // GroupedWindows::new sorts internally and caches the chromosome span.
    let windows_mapping: FxHashMap<String, GroupedWindows> = grouped_windows_by_chromosome
        .into_iter()
        .map(|(chr, windows)| {
            let strands = strand_mapping
                .as_mut()
                .and_then(|strands_by_chr| strands_by_chr.remove(&chr));
            (chr, GroupedWindows::new(windows, strands))
        })
        .collect();

    // Invert the group mapping to allow getting the group name from the group index
    let group_idx_to_name: FxHashMap<u64, String> = group_name_to_idx
        .iter()
        .map(|(name, &idx)| (idx, name.clone()))
        .collect();

    Ok((windows_mapping, group_idx_to_name, strand_detection))
}

/// Site orientation read from a BED file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg(loads_grouped_bed)]
pub(crate) enum Strand {
    /// Site has no directional interpretation.
    Unstranded,
    /// Site is forward-oriented.
    Forward,
    /// Site is reverse-oriented.
    Reverse,
}

/// Check whether a sampled BED field is one of the supported UCSC strand tokens.
///
/// This is a classifier, not strict parsing. The detector uses `None` to mean "this sampled
/// column is not strand-like". Once a column has been selected, `parse_bed_strand_token` turns
/// invalid values into user-facing errors.
#[cfg(loads_grouped_bed)]
fn classify_bed_strand_token(token: &str) -> Option<Strand> {
    match token {
        "." => Some(Strand::Unstranded),
        "+" => Some(Strand::Forward),
        "-" => Some(Strand::Reverse),
        _ => None,
    }
}

/// Parse a selected grouped BED strand field.
///
/// This function is strict. It is only used after the loader has decided which column is the
/// strand column, so any value other than `+`, `-`, or `.` is a malformed BED row and returns an
/// error with the line and column number.
#[cfg(loads_grouped_bed)]
fn parse_bed_strand_token(value: &str, lineno: usize, column_number: usize) -> Result<Strand> {
    classify_bed_strand_token(value).with_context(|| {
        format!(
            "BED parse error at line {}: invalid strand '{}' in column {}. Expected '+', '-', or '.'",
            lineno, value, column_number
        )
    })
}

/// Strand column found during grouped BED sampling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg(loads_grouped_bed)]
pub(crate) enum GroupedBedStrandColumn {
    Column5,
    Column6,
}

/// What the grouped BED sampling pass found about interval strands.
///
/// `sampled_rows` is the number of non-header BED rows used for detection. `column` is `Some`
/// only when the sampled rows identified a usable strand column. `saw_column6` lets commands warn
/// when a wide BED-like file was treated as unstranded because column 6 did not contain `+`, `-`,
/// or `.`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg(loads_grouped_bed)]
pub(crate) struct GroupedBedStrandDetection {
    pub(crate) sampled_rows: usize,
    pub(crate) column: Option<GroupedBedStrandColumn>,
    pub(crate) saw_column6: bool,
}

/// Return whether a BED line should be ignored before field parsing.
///
/// Blank lines, comments, and UCSC `track` or `browser` directives do not represent intervals and
/// should not count toward strand-column detection or parsed window totals.
#[cfg(loads_grouped_bed)]
fn should_skip_bed_line(line: &str) -> bool {
    if line.as_bytes().first().is_some_and(|value| *value == b'#') {
        return true;
    }
    if line.is_empty() {
        return true;
    }

    let trimmed_line_start = line.trim_start();
    trimmed_line_start.starts_with("track") || trimmed_line_start.starts_with("browser")
}

/// Detect whether grouped BED rows contain a strand column.
///
/// The grouped BED loader samples a bounded prefix of data rows instead of scanning the full file
/// up front. Column 6 is preferred because it matches standard BED6 layout. Column 5 is accepted
/// only when no sampled row has a column 6. If column 5 looks stranded but column 6 exists and is
/// not a strand column, the file is treated as ambiguous and rejected.
#[cfg(loads_grouped_bed)]
fn detect_grouped_bed_strand_column(bed: &Path) -> Result<GroupedBedStrandDetection> {
    let mut reader = open_text_reader(bed)?;
    let mut buffer = String::new();
    let mut sampled_rows = 0usize;
    let mut any_column5 = false;
    let mut any_column6 = false;
    let mut all_sampled_column5_are_strands = true;
    let mut all_sampled_column6_are_strands = true;

    while sampled_rows < GROUPED_BED_STRAND_SAMPLE_ROWS {
        buffer.clear();
        if reader.read_line(&mut buffer)? == 0 {
            break;
        }

        let line = buffer.trim_end_matches(['\n', '\r']);
        if should_skip_bed_line(line) {
            continue;
        }

        let columns: Vec<&str> = line.split_ascii_whitespace().collect();
        if columns.is_empty() {
            continue;
        }

        sampled_rows += 1;
        match columns.get(4) {
            Some(value) => {
                any_column5 = true;
                if classify_bed_strand_token(value).is_none() {
                    all_sampled_column5_are_strands = false;
                }
            }
            None => all_sampled_column5_are_strands = false,
        }

        match columns.get(5) {
            Some(value) => {
                any_column6 = true;
                if classify_bed_strand_token(value).is_none() {
                    all_sampled_column6_are_strands = false;
                }
            }
            None => all_sampled_column6_are_strands = false,
        }
    }

    if sampled_rows == 0 {
        bail!(
            "grouped BED file did not contain any interval rows. Empty BED files are not valid input"
        );
    }

    if any_column6 && all_sampled_column6_are_strands {
        return Ok(GroupedBedStrandDetection {
            sampled_rows,
            column: Some(GroupedBedStrandColumn::Column6),
            saw_column6: any_column6,
        });
    }

    if any_column6 && any_column5 && all_sampled_column5_are_strands {
        bail!(
            "grouped BED strand column is ambiguous: column 5 contains only '+', '-', or '.', but column 6 exists and is not a strand column. When 6 or more BED columns are supplied, put strands in column 6"
        );
    }

    if !any_column6 && any_column5 && all_sampled_column5_are_strands {
        return Ok(GroupedBedStrandDetection {
            sampled_rows,
            column: Some(GroupedBedStrandColumn::Column5),
            saw_column6: any_column6,
        });
    }
    Ok(GroupedBedStrandDetection {
        sampled_rows,
        column: None,
        saw_column6: any_column6,
    })
}

/// Parse one already-selected grouped BED strand field.
///
/// The loader decides which field to pass based on the detected file layout. This helper only
/// checks that the selected field exists and contains a supported strand token.
#[cfg(loads_grouped_bed)]
fn parse_grouped_bed_strand_value(
    value: Option<&str>,
    lineno: usize,
    column_number: usize,
) -> Result<Strand> {
    let value = value.with_context(|| {
        format!(
            "BED parse error at line {}: missing strand in column {}",
            lineno, column_number
        )
    })?;
    parse_bed_strand_token(value, lineno, column_number)
}

/// Owned collection of half-open windows with a cached genomic span.
///
/// Invariants
/// ----------
/// - `windows` should be sorted by start (ascending order).
/// - Coordinates are half-open: `[start, end)`.
/// - `strands`, when present, is one-to-one with `windows` and uses the same order.
#[derive(Debug, Clone)]
#[cfg(loads_grouped_bed)]
pub(crate) struct GroupedWindows {
    pub(crate) windows: Vec<IndexedInterval<u64>>, // (start, end, group idx)
    pub(crate) strands: Option<Vec<Strand>>,       // strands in the same order
    /// Cached outer envelope across all windows.
    #[allow(dead_code)]
    span: Span<i64>,
}

#[cfg(loads_grouped_bed)]
impl GroupedWindows {
    /// Construct from any window list (may be unsorted/overlapping).
    /// Ensures start- and end-sorted order (does not retain initial order)
    /// and computes span as `min(start)` .. `max(end)`.
    pub(crate) fn new(
        mut windows: Vec<IndexedInterval<u64>>,
        strands: Option<Vec<Strand>>,
    ) -> Self {
        match strands {
            Some(strands) => {
                assert_eq!(
                    windows.len(),
                    strands.len(),
                    "grouped window strands must match grouped window count"
                );
                let mut windows_and_strands: Vec<(IndexedInterval<u64>, Strand)> =
                    windows.drain(..).zip(strands).collect();
                windows_and_strands
                    .sort_unstable_by_key(|(window, _)| (window.start(), window.end()));
                let (sorted_windows, sorted_strands) = windows_and_strands.into_iter().unzip();
                GroupedWindows::from_sorted(sorted_windows, Some(sorted_strands))
            }
            None => {
                windows.sort_unstable_by_key(|window| (window.start(), window.end()));
                GroupedWindows::from_sorted(windows, None)
            }
        }
    }

    /// Construct from raw `(start, end, group_idx)` tuples.
    ///
    /// Use this when grouped BED parsing still works in tuple space.
    /// Prefer `new` and `from_sorted` when the windows are already checked.
    #[allow(dead_code)]
    pub(crate) fn from_tuples(
        windows: &[(u64, u64, u64)],
        strands: Option<Vec<Strand>>,
    ) -> crate::Result<Self> {
        Ok(Self::new(IndexedInterval::from_tuples(windows)?, strands))
    }

    #[cfg(feature = "cmd_midpoints")]
    pub(crate) fn mut_windows_and_strands(
        &mut self,
    ) -> (&mut Vec<IndexedInterval<u64>>, &mut Option<Vec<Strand>>) {
        (&mut self.windows, &mut self.strands)
    }

    /// Construct from a list you guarantee is already sorted by start (non-decreasing).
    /// Computes span as `min(start)` .. `max(end)` (robust to irregular ends).
    pub(crate) fn from_sorted(
        grouped_windows: Vec<IndexedInterval<u64>>,
        strands: Option<Vec<Strand>>,
    ) -> Self {
        if let Some(strands) = strands.as_ref() {
            assert_eq!(
                grouped_windows.len(),
                strands.len(),
                "grouped window strands must match grouped window count"
            );
        }
        debug_assert!(
            is_sorted_by_start_indexed(&grouped_windows),
            "windows must be start-sorted"
        );
        let span = if grouped_windows.is_empty() {
            Span::from_ordered(0, 0)
        } else {
            let min_start = grouped_windows[0].start() as i64;
            let max_end = grouped_windows
                .iter()
                .map(|window| window.end())
                .max()
                .unwrap() as i64;
            Span::from_ordered(min_start, max_end)
        };
        Self {
            windows: grouped_windows,
            strands,
            span,
        }
    }

    /// Number of windows.
    #[inline]
    #[cfg(feature = "cmd_midpoints")]
    pub(crate) fn len(&self) -> usize {
        self.windows.len()
    }

    /// True if there are no windows.
    #[inline]
    #[cfg(feature = "cmd_fcoverage")]
    pub(crate) fn is_empty(&self) -> bool {
        self.windows.is_empty()
    }

    /// Borrow the underlying windows.
    #[inline]
    pub(crate) fn windows_as_slice(&self) -> &[IndexedInterval<u64>] {
        &self.windows
    }

    /// Consume and return the inner vector.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn into_inner(self) -> Vec<IndexedInterval<u64>> {
        self.windows
    }

    /// Span start (inclusive).
    /// This is the most-left coordinate covered by any of the windows.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn span_start(&self) -> i64 {
        self.span.start()
    }

    /// Span end (exclusive).
    /// This is the most-right coordinate covered by any of the windows.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn span_end(&self) -> i64 {
        self.span.end()
    }

    /// Return the collection span.
    ///
    /// This is the outer envelope of the stored windows: the smallest start and
    /// largest end seen in the collection. It returns `Span<i64>` rather than
    /// `Interval<i64>` because empty collections use the empty span `[0, 0)`.
    ///
    /// There are no guarantees that every position inside this span is covered
    /// by a window.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn span(&self) -> Span<i64> {
        self.span
    }

    /// The intervals have strand information.
    #[allow(dead_code)]
    #[inline]
    pub(crate) fn has_strands(&self) -> bool {
        self.strands.is_some()
    }
}

#[cfg(loads_grouped_bed)]
impl AsRef<[IndexedInterval<u64>]> for GroupedWindows {
    fn as_ref(&self) -> &[IndexedInterval<u64>] {
        self.windows_as_slice()
    }
}

#[cfg(loads_grouped_bed)]
impl Default for GroupedWindows {
    fn default() -> Self {
        Self {
            windows: Vec::new(),
            strands: None,
            span: Span::from_ordered(0, 0),
        }
    }
}

/// Segment layout used when grouped BED rows need unique internal identities for reduction.
///
/// Each stored segment carries a stable `segment_idx` in the `IndexedInterval`, while the
/// same position in `group_idx_by_segment_idx` preserves the grouped row identity.
#[derive(Debug, Clone)]
#[cfg(feature = "cmd_fcoverage")]
pub(crate) struct GroupedCoverageLayout {
    pub(crate) segments_by_chr: FxHashMap<String, Windows>,
    pub(crate) group_idx_by_segment_idx: Vec<u64>,
    pub(crate) group_span_positions: FxHashMap<u64, u64>,
    pub(crate) group_idx_to_name: FxHashMap<u64, String>,
}

/// Chromosome-local grouped segments prepared before global segment indices are assigned.
#[derive(Debug)]
#[cfg(feature = "cmd_fcoverage")]
struct PreparedGroupedCoverageChromosome {
    chromosome: String,
    windows: Windows,
    group_indices: Vec<u64>,
    group_span_positions: FxHashMap<u64, u64>,
    segment_idx_offset: u64,
}

/// Build a grouped coverage layout for grouped `fcoverage` outputs.
///
/// Plain grouped actions keep every loaded interval as its own segment, while unique-base grouped
/// actions merge same-group overlaps and touches before assigning new internal segment indices.
///
/// Parameters
/// ----------
/// - `grouped_windows_by_chr`:
///   Grouped BED coordinates keyed by chromosome.
/// - `group_idx_to_name`:
///   Stable `group_idx -> group_name` mapping from the grouped BED loader.
/// - `chromosomes`:
///   Chromosome order used to assign deterministic segment indices.
/// - `unique_bases`:
///   Merge same-group overlaps and touches before indexing.
///
/// Returns
/// -------
/// - `layout`:
///   Per-chromosome segments plus stable maps back to the grouped row identity.
#[cfg(feature = "cmd_fcoverage")]
pub(crate) fn build_grouped_coverage_layout(
    grouped_windows_by_chr: &FxHashMap<String, GroupedWindows>,
    group_idx_to_name: &FxHashMap<u64, String>,
    chromosomes: &[String],
    unique_bases: bool,
) -> Result<GroupedCoverageLayout> {
    // Output shape:
    // - `segments_by_chr` is what downstream coverage code iterates over
    // - `group_idx_by_segment_idx` lets reducers recover the original grouped BED row identity
    // - `group_span_positions` tracks the total number of positions represented by each group

    // Chromosomes are independent until global segment indices are assigned
    let mut prepared_chromosomes: Vec<Option<PreparedGroupedCoverageChromosome>> = chromosomes
        .par_iter()
        .map(|chromosome| {
            let Some(grouped_windows) = grouped_windows_by_chr.get(chromosome) else {
                return None;
            };

            Some(prepare_grouped_coverage_chromosome(
                chromosome,
                grouped_windows,
                unique_bases,
            ))
        })
        .collect();

    // Follow the requested chromosome order when assigning each chromosome's global index range
    let mut next_segment_idx = 0_u64;
    for prepared in prepared_chromosomes.iter_mut().flatten() {
        prepared.segment_idx_offset = next_segment_idx;
        let chromosome_segment_count = u64::try_from(prepared.windows.len())
            .context("grouped coverage has more segments than can be indexed by u64")?;
        next_segment_idx = next_segment_idx
            .checked_add(chromosome_segment_count)
            .context("grouped coverage segment index overflow")?;
    }

    // Apply each chromosome's global segment index offset in parallel
    prepared_chromosomes
        .par_iter_mut()
        .try_for_each(|prepared| -> Result<()> {
            let Some(prepared) = prepared else {
                return Ok(());
            };
            ensure!(
                prepared.windows.len() == prepared.group_indices.len(),
                "grouped coverage prepared {} segments but {} group indices for chromosome '{}'",
                prepared.windows.len(),
                prepared.group_indices.len(),
                prepared.chromosome
            );

            for (local_segment_idx, segment) in prepared.windows.windows.iter_mut().enumerate() {
                let local_segment_idx = u64::try_from(local_segment_idx).context(
                    "grouped coverage chromosome has more segments than can be indexed by u64",
                )?;
                segment.idx = prepared
                    .segment_idx_offset
                    .checked_add(local_segment_idx)
                    .context("grouped coverage segment index overflow")?;
            }
            Ok(())
        })?;

    // Combine chromosome-local group spans in parallel
    let group_span_positions = prepared_chromosomes
        .par_iter_mut()
        .filter_map(|prepared| prepared.as_mut())
        .map(|prepared| std::mem::take(&mut prepared.group_span_positions))
        .reduce(
            FxHashMap::default,
            |mut combined_spans, chromosome_spans| {
                for (group_idx, span_positions) in chromosome_spans {
                    *combined_spans.entry(group_idx).or_insert(0) += span_positions;
                }
                combined_spans
            },
        );

    let total_segments = usize::try_from(next_segment_idx)
        .context("grouped coverage has more segments than fit in memory")?;
    let mut segments_by_chr: FxHashMap<String, Windows> =
        FxHashMap::with_capacity_and_hasher(grouped_windows_by_chr.len(), Default::default());
    let mut group_idx_by_segment_idx: Vec<u64> = Vec::with_capacity(total_segments);

    // Assemble chromosome results in the same order used to assign their global index ranges
    for mut prepared in prepared_chromosomes.into_iter().flatten() {
        group_idx_by_segment_idx.append(&mut prepared.group_indices);
        segments_by_chr.insert(prepared.chromosome, prepared.windows);
    }

    Ok(GroupedCoverageLayout {
        segments_by_chr,
        group_idx_by_segment_idx,
        group_span_positions,
        group_idx_to_name: group_idx_to_name.clone(),
    })
}

/// Build and sort the coverage segments for one chromosome.
#[cfg(feature = "cmd_fcoverage")]
fn prepare_grouped_coverage_chromosome(
    chromosome: &str,
    grouped_windows: &GroupedWindows,
    unique_bases: bool,
) -> PreparedGroupedCoverageChromosome {
    // Build the chromosome-local segments that coverage will be summed over:
    // - raw mode keeps every original grouped BED interval as its own segment
    // - unique-base mode first merges overlaps and touching intervals within each group
    let mut indexed_segments: Vec<IndexedInterval<u64>> = if unique_bases {
        merged_group_segments(grouped_windows.windows_as_slice())
    } else {
        grouped_windows.windows_as_slice().to_vec()
    };

    // Sort before assigning fresh local indices so the final global identifiers stay stable
    indexed_segments
        .sort_unstable_by_key(|segment| (segment.start(), segment.end(), segment.idx()));

    let mut group_indices = Vec::with_capacity(indexed_segments.len());
    let mut group_span_positions: FxHashMap<u64, u64> = FxHashMap::default();
    for segment in &indexed_segments {
        let group_idx = segment.idx();
        group_indices.push(group_idx);
        *group_span_positions.entry(group_idx).or_insert(0) += segment.len();
    }

    PreparedGroupedCoverageChromosome {
        chromosome: chromosome.to_string(),
        windows: Windows::from_sorted(indexed_segments),
        group_indices,
        group_span_positions,
        segment_idx_offset: 0,
    }
}

#[cfg(feature = "cmd_fcoverage")]
fn merged_group_segments(grouped_windows: &[IndexedInterval<u64>]) -> Vec<IndexedInterval<u64>> {
    // First regroup the chromosome-local windows by their original grouped BED row id.
    // We only merge within a group. Different groups must stay separate even if they overlap.
    //
    // After this regrouping step, each hashmap entry can be processed independently as:
    //   "all intervals that belong to one grouped BED row on this chromosome".
    let mut intervals_by_group: FxHashMap<u64, Vec<Interval<u64>>> = FxHashMap::default();

    for window in grouped_windows {
        intervals_by_group
            .entry(window.idx())
            .or_default()
            .push(window.interval);
    }

    let mut merged_segments: Vec<IndexedInterval<u64>> = Vec::new();
    for (group_idx, mut intervals) in intervals_by_group {
        // `push_merged_interval` assumes start-sorted input, so sort once per group before
        // collapsing intervals.
        intervals.sort_unstable_by_key(|interval| (interval.start(), interval.end()));

        // This mutable accumulator is the current merged representation for one group.
        // We feed intervals into it from left to right. `push_merged_interval` either:
        // - extends the last merged interval in-place when the new interval overlaps or touches it
        // - or appends a brand-new merged interval when there is a gap
        //
        // That is why the argument is `&mut merged_for_group`: the helper updates the current
        // merged state directly instead of returning a new vector at every step.
        let mut merged_for_group: Vec<Interval<u64>> = Vec::with_capacity(intervals.len());
        for interval in intervals {
            // Unique-base grouped coverage treats touching intervals as one continuous covered span.
            // Example: [10, 20) and [20, 30) become [10, 30).
            push_merged_interval(
                &mut merged_for_group,
                interval,
                TouchingMergePolicy::MergeTouching,
            );
        }

        // Preserve the originating group id next to each merged segment so the caller can later
        // assign fresh segment ids while still knowing which grouped BED row each segment belongs to.
        // At this point, `merged_for_group` contains the minimal non-touching interval list for
        // this one group on this one chromosome.
        merged_segments.extend(
            merged_for_group
                .into_iter()
                .map(|interval| IndexedInterval::from_interval(interval, group_idx)),
        );
    }

    merged_segments
}

/// Write a TSV mapping from `group_idx` -> `group_name`.
///
/// - Output has a header: `group_idx\tgroup_name`
/// - Rows are sorted by `group_idx` ascending for determinism.
/// - Creates the parent directory if needed.
#[cfg(feature = "cmd_fcoverage")]
pub(crate) fn write_group_idx_to_name_tsv<P: AsRef<Path>>(
    output_path: P,
    group_idx_to_name: &FxHashMap<u64, String>,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path).with_context(|| format!("Creating TSV file {:?}", path))?;
    let mut writer = BufWriter::new(file);

    // Header
    writeln!(writer, "group_idx\tgroup_name")
        .with_context(|| format!("Writing header to {:?}", path))?;

    // Collect and sort by index for stable output
    let mut entries: Vec<(u64, &str)> = group_idx_to_name
        .iter()
        .map(|(idx, name)| (*idx, name.as_str()))
        .collect();
    entries.sort_unstable_by_key(|(idx, _)| *idx);

    // Write rows
    for (idx, name) in entries {
        // Sanitize tabs/newlines to keep TSV well-formed (should not be needed but may reduce errors)
        let name = name.replace('\t', "    ").replace('\n', " ");
        writeln!(writer, "{idx}\t{name}")
            .with_context(|| format!("Writing row for group_idx {idx} to {:?}", path))?;
    }

    Ok(())
}

/* Scored windows */

/// Load *scored* windows from a BED file into a per-chromosome map.
///
/// The original window index is added. Any valid window increases the index,
/// even if they are filtered from the returned windows.
///
/// Parameters
/// ----------
///  - bed: Path to BED file with scores (float64) in the fourth column.
///  - chromosomes: Names of chromosomes to include in output,
///    even when not present in the BED file.
///  - filter_fn: Function for deciding whether to include
///    an interval. Should take in the `chr,start,end,score` values
///    and return `true` (keep) or `false` (discard).
///  - exp_num_windows: Optional number of expected windows
///    in the BED file before filtering. Returns an error
///    if the incorrect number of windows are observed.
///    
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, score, original index)'.
#[allow(dead_code)]
pub(crate) fn load_scored_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    filter_fn: Option<&dyn Fn(&str, u64, u64, f64) -> bool>,
    exp_num_windows: Option<u64>,
) -> Result<FxHashMap<String, ScoredWindows>> {
    let mut reader = open_text_reader(bed.as_ref())?; // Works with &Path, PathBuf, &str

    // Optional whitelist of chromosomes
    let mut scored_windows_by_chromosome: FxHashMap<String, Vec<(u64, u64, u64, f64)>> =
        FxHashMap::default();
    let allowed_chromosomes: Option<FxHashSet<&str>> = chromosomes.map(|chr_list| {
        let mut allowed = FxHashSet::with_capacity_and_hasher(chr_list.len(), Default::default());
        for chr in chr_list {
            allowed.insert(chr.as_str());
            scored_windows_by_chromosome.entry(chr.clone()).or_default();
        }
        allowed
    });

    // Reuse a single buffer for all lines
    let mut buf = String::new();
    let mut lineno: usize = 0;
    let mut orig_win_idx: u64 = 0; // Counter for all *valid* windows whether filtered out or not
    let mut sniffed_rows = 0usize;

    loop {
        buf.clear();
        let bytes_read = reader.read_line(&mut buf)?;
        if bytes_read == 0 {
            break;
        }
        lineno += 1;

        // Fast skips
        if buf.as_bytes().first().is_some_and(|b| *b == b'#') {
            continue;
        }
        let line = buf.trim_end_matches(['\n', '\r']);

        if line.is_empty() {
            continue;
        }

        // Skip UCSC header directives in BED files
        let trimmed_line_start = line.trim_start();
        if trimmed_line_start.starts_with("track") || trimmed_line_start.starts_with("browser") {
            continue;
        }

        if sniffed_rows < BED_FORMAT_SNIFF_ROWS {
            sniff_bed3_line(line, lineno)?;
            sniffed_rows += 1;
        }

        // Strict parse of first 3 BED columns without allocating a Vec
        let mut fields = line.split_ascii_whitespace();

        // Get the three coordinate fields to ensure at least three values are present
        // We need to know if the line represents a valid window (whether we want it or not)

        let chr = match fields.next() {
            Some(s) => s,
            None => continue,
        };

        // NOTE: "!Allowed" != "!Valid"
        // A chromosome name can be valid with the line representing a genomic window
        // but disallowed due to not being in the chromosome whitelist
        // but any valid window is considered for the original index incrementing
        let is_allowed_chrom = if let Some(allowed_chroms) = &allowed_chromosomes {
            allowed_chroms.contains(chr)
        } else {
            // If we don't have a whitelist, we must assume it's an allowed chromosome name (whether valid or not)
            true
        };

        let start_str = match fields.next() {
            Some(s) => s,
            None => {
                if is_allowed_chrom {
                    // If initial value is known to be a valid chromosome (when it's in the whitelist)
                    // or we don't know whether it's valid (we have no whitelist)
                    // we expect the start field to exist
                    bail!(
                        "BED parse error at line {}: missing start for chromosome: '{}'",
                        lineno,
                        chr
                    );
                }
                continue;
            }
        };

        let end_str = match fields.next() {
            Some(s) => s,
            None => {
                if is_allowed_chrom {
                    // If initial value is known to be a valid chromosome (when it's in the whitelist)
                    // or we don't know whether it's valid (we have no whitelist)
                    // we expect the end field to exist
                    bail!(
                        "BED parse error at line {}: missing end for chromosome: '{}'",
                        lineno,
                        chr
                    );
                }
                continue;
            }
        };

        // Skip if not an allowed chromosome
        // Increment iff window is assumed valid
        if !is_allowed_chrom {
            // Only increment original window index coordinates make a valid window
            if start_end_are_valid_coordinates(start_str, end_str).is_some() {
                orig_win_idx += 1;
            }
            continue;
        }

        // Parse coordinates
        let start: u64 = start_str.parse().with_context(|| {
            format!(
                "BED parse error at line {}: invalid start '{}'",
                lineno, start_str
            )
        })?;
        let end: u64 = end_str.parse().with_context(|| {
            format!(
                "BED parse error at line {}: invalid end '{}'",
                lineno, end_str
            )
        })?;

        ensure!(
            end > start,
            "BED parse error at line {}: end ({}) must be greater than start ({})",
            lineno,
            end,
            start
        );

        // Extract score field
        let score_str = fields
            .next()
            .with_context(|| format!("BED parse error at line {}: missing score", lineno))?;

        let score: f64 = score_str.parse().with_context(|| {
            format!(
                "BED parse error at line {}: invalid score '{}'",
                lineno, score_str
            )
        })?;

        // At this point, we assume the window to be a valid window
        // and increment the counter for the original index
        let current_orig_win_idx = orig_win_idx;
        orig_win_idx += 1;

        // Apply passed filtering function
        if let Some(filterer) = filter_fn
            && !filterer(chr, start, end, score)
        {
            continue;
        }

        scored_windows_by_chromosome
            .entry(chr.to_string())
            .or_default()
            .push((start, end, current_orig_win_idx, score));
    }

    if let Some(expected_num_windows) = exp_num_windows {
        ensure!(
            expected_num_windows == orig_win_idx,
            "the BED file did not contain the correct number of windows: obs: {} != exp: {}",
            orig_win_idx,
            expected_num_windows
        );
    }

    // Convert parsed tuples into typed scored windows.
    // ScoredWindows::from_tuples delegates to ScoredWindows::new, which sorts internally.
    let windows_mapping: FxHashMap<String, ScoredWindows> = scored_windows_by_chromosome
        .into_iter()
        .map(|(chr, windows)| Ok((chr.to_string(), ScoredWindows::from_tuples(&windows)?)))
        .collect::<crate::Result<_>>()?;

    Ok(windows_mapping)
}

/// Owned collection of half-open windows with a cached genomic span.
///
/// Invariants
/// ----------
/// - `windows` should be sorted by start (ascending order).
/// - Coordinates are half-open: `[start, end)`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub(crate) struct ScoredWindows {
    pub(crate) windows: Vec<ScoredInterval<u64>>, // (start, end, original_idx, score)
    /// Cached outer envelope across all windows.
    span: Span<i64>,
}

#[allow(dead_code)]
impl ScoredWindows {
    /// Construct from any window list (may be unsorted/overlapping).
    /// Ensures start- and end-sorted order (does not retain initial order)
    /// and computes span as `min(start)` .. `max(end)`.
    pub(crate) fn new(mut windows: Vec<ScoredInterval<u64>>) -> Self {
        windows.sort_unstable_by_key(|window| (window.start(), window.end()));
        ScoredWindows::from_sorted(windows)
    }

    /// Construct from raw `(start, end, idx, score)` tuples.
    pub(crate) fn from_tuples(windows: &[(u64, u64, u64, f64)]) -> crate::Result<Self> {
        Ok(Self::new(ScoredInterval::from_tuples(windows)?))
    }

    /// Convert to Windows collection by dropping the score.
    pub(crate) fn to_windows(&self) -> Windows {
        Windows::from_sorted(self.windows.iter().map(|window| window.window).collect())
    }

    /// Construct from a list you guarantee is already sorted by start (non-decreasing).
    /// Computes span as `min(start)` .. `max(end)` (robust to irregular ends).
    pub(crate) fn from_sorted(windows: Vec<ScoredInterval<u64>>) -> Self {
        debug_assert!(
            is_sorted_by_start_with_scores(&windows),
            "windows must be start-sorted"
        );
        let span = if windows.is_empty() {
            Span::from_ordered(0, 0)
        } else {
            let min_start = windows[0].start() as i64;
            let max_end = windows.iter().map(|window| window.end()).max().unwrap() as i64;
            Span::from_ordered(min_start, max_end)
        };
        Self { windows, span }
    }

    /// Number of windows.
    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.windows.len()
    }

    /// True if there are no windows.
    #[inline]
    pub(crate) fn is_empty(&self) -> bool {
        self.windows.is_empty()
    }

    /// Borrow the underlying windows.
    #[inline]
    pub(crate) fn as_slice(&self) -> &[ScoredInterval<u64>] {
        &self.windows
    }

    /// Consume and return the inner vector.
    #[inline]
    pub(crate) fn into_inner(self) -> Vec<ScoredInterval<u64>> {
        self.windows
    }

    /// Span start (inclusive).
    /// This is the most-left coordinate covered by any of the windows.
    #[inline]
    pub(crate) fn span_start(&self) -> i64 {
        self.span.start()
    }

    /// Span end (exclusive).
    /// This is the most-right coordinate covered by any of the windows.
    #[inline]
    pub(crate) fn span_end(&self) -> i64 {
        self.span.end()
    }

    /// Return the collection span.
    ///
    /// This is the outer envelope of the stored windows: the smallest start and
    /// largest end seen in the collection. It returns `Span<i64>` rather than
    /// `Interval<i64>` because empty collections use the empty span `[0, 0)`.
    ///
    /// There are no guarantees that every position inside this span is covered
    /// by a window.
    #[inline]
    pub(crate) fn span(&self) -> Span<i64> {
        self.span
    }
}

/* Other utilities */

/// Check whether line looks like a header or an observation
#[cfg(feature = "cmd_prepare_windows")]
pub(crate) fn line_looks_like_header(line: &str, separator: char) -> bool {
    let trimmed = line.trim_start();
    if trimmed.starts_with('#') {
        return true;
    }
    let fields: Vec<&str> = line
        .trim_end_matches(['\n', '\r'])
        .split(separator)
        .collect();
    if fields.len() < 3 {
        return true;
    }
    let start_ok = fields[1].trim().parse::<u64>().is_ok();
    let end_ok = fields[2].trim().parse::<u64>().is_ok();
    !(start_ok && end_ok)
}

// TODO: Generalize and test
/// Detect whether a file appears to have a header by peeking the first non-comment line.
///
/// Parameters
/// ----------
/// - path:
///     Path to file.
/// - separator:
///     Field separator.
///
/// Returns
/// -------
/// - has_header:
///     True if a header is likely present.
#[cfg(feature = "cmd_prepare_windows")]
pub(crate) fn detect_header(path: &Path, separator: char) -> Result<bool> {
    let mut reader = open_text_reader(path)?;
    let mut line = String::new();
    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            return Ok(false);
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        return Ok(line_looks_like_header(&line, separator));
    }
}

#[cfg(test)]
mod tests {
    include!("bed_tests.rs");
}
