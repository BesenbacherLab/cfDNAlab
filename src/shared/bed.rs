use crate::shared::interval::{IndexedInterval, ScoredInterval, Span};
use crate::shared::io::open_text_reader;
use anyhow::{Context, Result, bail, ensure};
use fxhash::{FxHashMap, FxHashSet};
use std::{
    fs::File,
    io::{BufRead, BufWriter, Write},
    path::Path,
};

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
///
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, original window index)'.
pub fn load_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    filter_fn: Option<&dyn Fn(&str, u64, u64) -> bool>,
    exp_num_windows: Option<u64>,
) -> Result<FxHashMap<String, Windows>> {
    let mut reader = open_text_reader(bed.as_ref())?;

    // Optional whitelist of chromosomes
    let mut vec_mapping: FxHashMap<String, Vec<(u64, u64, u64)>> = FxHashMap::default();
    let allowed_chromosomes: Option<FxHashSet<&str>> = chromosomes.map(|chr_list| {
        let mut allowed = FxHashSet::with_capacity_and_hasher(chr_list.len(), Default::default());
        for chr in chr_list {
            allowed.insert(chr.as_str());
            vec_mapping.entry(chr.clone()).or_default();
        }
        allowed
    });

    // Reuse a single buffer for all lines
    let mut buf = String::new();
    let mut lineno: usize = 0;
    let mut orig_win_idx: u64 = 0; // Counter for all *valid* windows whether filtered out or not

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

        vec_mapping
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

    let windows_mapping: FxHashMap<String, Windows> = vec_mapping
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

/// Owned collection of half-open windows with a cached genomic span.
///
/// Invariants
/// ----------
/// - `windows` should be sorted by start (ascending order).
/// - Coordinates are half-open: `[start, end)`.
#[derive(Debug, Clone)]
pub struct Windows {
    pub windows: Vec<IndexedInterval<u64>>,
    /// Cached outer envelope across all windows.
    span: Span<i64>,
}

impl Windows {
    /// Construct from any window list (may be unsorted/overlapping).
    /// Ensures start- and end-sorted order (does not retain initial order)
    /// and computes span as `min(start)` .. `max(end)`.
    pub fn new(mut windows: Vec<IndexedInterval<u64>>) -> Self {
        windows.sort_unstable_by_key(|window| (window.start(), window.end()));
        Windows::from_sorted(windows)
    }

    /// Construct from raw `(start, end, idx)` tuples.
    ///
    /// Use this when a loader or tiny fixture still naturally produces tuples.
    /// Prefer `new` and `from_sorted` when the windows are already checked.
    pub fn from_tuples(windows: &[(u64, u64, u64)]) -> crate::Result<Self> {
        Ok(Self::new(IndexedInterval::from_tuples(windows)?))
    }

    /// Construct from a list you guarantee is already sorted by start (non-decreasing).
    /// Computes span as `min(start)` .. `max(end)` (robust to irregular ends).
    pub fn from_sorted(indexed_windows: Vec<IndexedInterval<u64>>) -> Self {
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
    #[inline]
    pub fn len(&self) -> usize {
        self.windows.len()
    }

    /// True if there are no windows.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.windows.is_empty()
    }

    /// Borrow the underlying windows.
    #[inline]
    pub fn as_slice(&self) -> &[IndexedInterval<u64>] {
        &self.windows
    }

    /// Consume and return the inner vector.
    #[inline]
    pub fn into_inner(self) -> Vec<IndexedInterval<u64>> {
        self.windows
    }

    /// Span start (inclusive).
    /// This is the most-left coordinate covered by any of the windows.
    #[inline]
    pub fn span_start(&self) -> i64 {
        self.span.start()
    }

    /// Span end (exclusive).
    /// This is the most-right coordinate covered by any of the windows.
    #[inline]
    pub fn span_end(&self) -> i64 {
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
    pub fn span(&self) -> Span<i64> {
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
    ///     First index to assign to the merged interval series.
    ///
    /// Returns
    /// -------
    /// - (merged, next_start_idx):
    ///     - `merged`: New `Windows` with merged, start-sorted `(start, end, new_idx)` tuples.
    ///     - `next_start_idx`: `start_idx + merged.len()`; pass this to the next chromosome.
    pub fn into_flattened_reindexed(self, start_idx: u64) -> (Windows, u64) {
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

    /// Borrowing variant: leaves `self` intact and returns a flattened copy.
    ///
    /// Parameters
    /// ----------
    /// - start_idx:
    ///     Starting index for the first merged interval.
    ///
    /// Returns
    /// -------
    /// - (merged, next_start_idx):
    ///     See `into_flattened_reindexed`.
    pub fn flattened_reindexed(&self, start_idx: u64) -> (Windows, u64) {
        // Clone once, then consume in the main routine to avoid duplicating logic.
        self.clone().into_flattened_reindexed(start_idx)
    }
}

fn is_sorted_by_start_indexed(ws: &[IndexedInterval<u64>]) -> bool {
    ws.windows(2)
        .all(|window_pair| window_pair[0].start() <= window_pair[1].start())
}

#[inline]
fn is_sorted_by_start_with_scores(ws: &[ScoredInterval<u64>]) -> bool {
    ws.windows(2)
        .all(|window_pair| window_pair[0].start() <= window_pair[1].start())
}

/* GROUPED bed files */

/// Load *grouped* windows from a BED file into a per-chromosome map.
///
/// Parameters
/// ----------
///  - bed: Path to BED file with group names in the fourth column.
///  - chromosomes: Names of chromosomes to include in output,
///    even when not present in the BED file.
///  - filter_fn: Function for deciding whether to include
///    an interval. Should take in the `chr,start,end` values
///    and return `true` (keep) or `false` (discard).
///  - exp_num_windows: Optional number of expected windows
///    in the BED file before filtering. Returns an error
///    if the incorrect number of windows are observed.
///
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, group index)'.
///
///  - Mapping of 'group index -> group name'.
pub fn load_grouped_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    filter_fn: Option<&dyn Fn(&str, u64, u64) -> bool>,
    exp_num_windows: Option<u64>,
) -> Result<(FxHashMap<String, GroupedWindows>, FxHashMap<u64, String>)> {
    let mut reader = open_text_reader(bed.as_ref())?; // Works with &Path, PathBuf, &str

    // Optional whitelist of chromosomes
    let mut vec_mapping: FxHashMap<String, Vec<(u64, u64, u64)>> = FxHashMap::default();
    let allowed_chromosomes: Option<FxHashSet<&str>> = chromosomes.map(|chr_list| {
        let mut allowed = FxHashSet::with_capacity_and_hasher(chr_list.len(), Default::default());
        for chr in chr_list {
            allowed.insert(chr.as_str());
            vec_mapping.entry(chr.clone()).or_default();
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

        vec_mapping
            .entry(chr.to_string())
            .or_default()
            .push((start, end, group_idx));
    }

    if let Some(expected_num_windows) = exp_num_windows {
        ensure!(
            expected_num_windows == orig_win_idx,
            "the BED file did not contain the correct number of windows: obs: {} != exp: {}",
            orig_win_idx,
            expected_num_windows
        );
    }

    // Convert parsed tuples into typed grouped windows.
    // GroupedWindows::from_tuples delegates to GroupedWindows::new, which sorts internally.
    let windows_mapping: FxHashMap<String, GroupedWindows> = vec_mapping
        .into_iter()
        .map(|(chr, windows)| Ok((chr.to_string(), GroupedWindows::from_tuples(&windows)?)))
        .collect::<crate::Result<_>>()?;

    // Invert the group mapping to allow getting the group name from the group index
    let group_idx_to_name: FxHashMap<u64, String> = group_name_to_idx
        .iter()
        .map(|(name, &idx)| (idx, name.clone()))
        .collect();

    Ok((windows_mapping, group_idx_to_name))
}

/// Owned collection of half-open windows with a cached genomic span.
///
/// Invariants
/// ----------
/// - `windows` should be sorted by start (ascending order).
/// - Coordinates are half-open: `[start, end)`.
#[derive(Debug, Clone)]
pub struct GroupedWindows {
    pub windows: Vec<IndexedInterval<u64>>, // (start, end, group idx)
    /// Cached outer envelope across all windows.
    span: Span<i64>,
}

impl GroupedWindows {
    /// Construct from any window list (may be unsorted/overlapping).
    /// Ensures start- and end-sorted order (does not retain initial order)
    /// and computes span as `min(start)` .. `max(end)`.
    pub fn new(mut windows: Vec<IndexedInterval<u64>>) -> Self {
        windows.sort_unstable_by_key(|window| (window.start(), window.end()));
        GroupedWindows::from_sorted(windows)
    }

    /// Construct from raw `(start, end, group_idx)` tuples.
    ///
    /// Use this when grouped BED parsing still works in tuple space.
    /// Prefer `new` and `from_sorted` when the windows are already checked.
    pub fn from_tuples(windows: &[(u64, u64, u64)]) -> crate::Result<Self> {
        Ok(Self::new(IndexedInterval::from_tuples(windows)?))
    }

    /// Construct from a list you guarantee is already sorted by start (non-decreasing).
    /// Computes span as `min(start)` .. `max(end)` (robust to irregular ends).
    pub fn from_sorted(grouped_windows: Vec<IndexedInterval<u64>>) -> Self {
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
            span,
        }
    }

    /// Number of windows.
    #[inline]
    pub fn len(&self) -> usize {
        self.windows.len()
    }

    /// True if there are no windows.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.windows.is_empty()
    }

    /// Borrow the underlying windows.
    #[inline]
    pub fn as_slice(&self) -> &[IndexedInterval<u64>] {
        &self.windows
    }

    /// Consume and return the inner vector.
    #[inline]
    pub fn into_inner(self) -> Vec<IndexedInterval<u64>> {
        self.windows
    }

    /// Span start (inclusive).
    /// This is the most-left coordinate covered by any of the windows.
    #[inline]
    pub fn span_start(&self) -> i64 {
        self.span.start()
    }

    /// Span end (exclusive).
    /// This is the most-right coordinate covered by any of the windows.
    #[inline]
    pub fn span_end(&self) -> i64 {
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
    pub fn span(&self) -> Span<i64> {
        self.span
    }
}

/// Write a TSV mapping from `group_idx` -> `group_name`.
///
/// - Output has a header: `group_idx\tgroup_name`
/// - Rows are sorted by `group_idx` ascending for determinism.
/// - Creates the parent directory if needed.
pub fn write_group_idx_to_name_tsv<P: AsRef<Path>>(
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
pub fn load_scored_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    filter_fn: Option<&dyn Fn(&str, u64, u64, f64) -> bool>,
    exp_num_windows: Option<u64>,
) -> Result<FxHashMap<String, ScoredWindows>> {
    let mut reader = open_text_reader(bed.as_ref())?; // Works with &Path, PathBuf, &str

    // Optional whitelist of chromosomes
    let mut vec_mapping: FxHashMap<String, Vec<(u64, u64, u64, f64)>> = FxHashMap::default();
    let allowed_chromosomes: Option<FxHashSet<&str>> = chromosomes.map(|chr_list| {
        let mut allowed = FxHashSet::with_capacity_and_hasher(chr_list.len(), Default::default());
        for chr in chr_list {
            allowed.insert(chr.as_str());
            vec_mapping.entry(chr.clone()).or_default();
        }
        allowed
    });

    // Reuse a single buffer for all lines
    let mut buf = String::new();
    let mut lineno: usize = 0;
    let mut orig_win_idx: u64 = 0; // Counter for all *valid* windows whether filtered out or not

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

        vec_mapping.entry(chr.to_string()).or_default().push((
            start,
            end,
            current_orig_win_idx,
            score,
        ));
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
    let windows_mapping: FxHashMap<String, ScoredWindows> = vec_mapping
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
pub struct ScoredWindows {
    pub windows: Vec<ScoredInterval<u64>>, // (start, end, original_idx, score)
    /// Cached outer envelope across all windows.
    span: Span<i64>,
}

impl ScoredWindows {
    /// Construct from any window list (may be unsorted/overlapping).
    /// Ensures start- and end-sorted order (does not retain initial order)
    /// and computes span as `min(start)` .. `max(end)`.
    pub fn new(mut windows: Vec<ScoredInterval<u64>>) -> Self {
        windows.sort_unstable_by_key(|window| (window.start(), window.end()));
        ScoredWindows::from_sorted(windows)
    }

    /// Construct from raw `(start, end, idx, score)` tuples.
    pub fn from_tuples(windows: &[(u64, u64, u64, f64)]) -> crate::Result<Self> {
        Ok(Self::new(ScoredInterval::from_tuples(windows)?))
    }

    /// Convert to Windows collection by dropping the score.
    pub fn to_windows(&self) -> Windows {
        Windows::from_sorted(self.windows.iter().map(|window| window.window).collect())
    }

    /// Construct from a list you guarantee is already sorted by start (non-decreasing).
    /// Computes span as `min(start)` .. `max(end)` (robust to irregular ends).
    pub fn from_sorted(windows: Vec<ScoredInterval<u64>>) -> Self {
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
    pub fn len(&self) -> usize {
        self.windows.len()
    }

    /// True if there are no windows.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.windows.is_empty()
    }

    /// Borrow the underlying windows.
    #[inline]
    pub fn as_slice(&self) -> &[ScoredInterval<u64>] {
        &self.windows
    }

    /// Consume and return the inner vector.
    #[inline]
    pub fn into_inner(self) -> Vec<ScoredInterval<u64>> {
        self.windows
    }

    /// Span start (inclusive).
    /// This is the most-left coordinate covered by any of the windows.
    #[inline]
    pub fn span_start(&self) -> i64 {
        self.span.start()
    }

    /// Span end (exclusive).
    /// This is the most-right coordinate covered by any of the windows.
    #[inline]
    pub fn span_end(&self) -> i64 {
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
    pub fn span(&self) -> Span<i64> {
        self.span
    }
}

/* Other utilities */

/// Check whether line looks like a header or an observation
pub fn line_looks_like_header(line: &str, separator: char) -> bool {
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
pub fn detect_header(path: &Path, separator: char) -> Result<bool> {
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
