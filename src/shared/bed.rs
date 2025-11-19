use crate::shared::io::open_text_reader;
use anyhow::{Context, Result, ensure};
use fxhash::{FxHashMap, FxHashSet};
use std::{
    fs::File,
    io::{BufRead, BufWriter, Write},
    path::Path,
};

/// Load windows from a BED file into a per-chromosome map.
///
/// Parameters
/// ----------
///  - bed: Path to BED file.
///  - chromosomes: Names of chromosomes to include in output,
///    even when not present in the BED file.
///  - filter_fn: Function for deciding whether to include
///    an interval. Should take in the `chr,start,end` values
///    and return `true` (keep) or `false` (discard).
///
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, original window index)'.
pub fn load_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    filter_fn: Option<&dyn Fn(&str, u64, u64) -> bool>,
) -> Result<FxHashMap<String, Windows>> {
    let mut reader = open_text_reader(bed.as_ref())?;

    // Optional whitelist of chromosomes
    let mut vec_mapping: FxHashMap<String, Vec<(u64, u64, u64)>> = FxHashMap::default();
    let allowed: Option<FxHashSet<&str>> = chromosomes.map(|chr_list| {
        let mut allowed = FxHashSet::with_capacity_and_hasher(chr_list.len(), Default::default());
        for chr in chr_list {
            allowed.insert(chr.as_str());
            vec_mapping.entry(chr.clone()).or_default();
        }
        allowed
    });

    // Reuse a single buffer for all lines
    let mut buf = String::new();
    let mut win_idx: u64 = 0;
    let mut lineno: usize = 0;

    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
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
        let ls = line.trim_start();
        if ls.starts_with("track") || ls.starts_with("browser") {
            continue;
        }

        // Strict parse of first 3 BED columns without allocating a Vec
        let mut it = line.split_ascii_whitespace();

        let chr = match it.next() {
            Some(s) => s,
            None => continue,
        };
        if let Some(allowed_chroms) = &allowed {
            if !allowed_chroms.contains(chr) {
                continue;
            }
        }

        let start_str = it
            .next()
            .with_context(|| format!("BED parse error at line {}: missing start", lineno))?;
        let end_str = it
            .next()
            .with_context(|| format!("BED parse error at line {}: missing end", lineno))?;

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

        if let Some(filterer) = filter_fn {
            if !filterer(chr, start, end) {
                continue;
            }
        }

        vec_mapping
            .entry(chr.to_string())
            .or_default()
            .push((start, end, win_idx));
        win_idx += 1;
    }

    let windows_mapping: FxHashMap<String, Windows> = vec_mapping
        .into_iter()
        .map(|(chr, v)| (chr, Windows::new(v)))
        .collect();

    Ok(windows_mapping)
}

/// Owned collection of half-open windows with a cached genomic span.
///
/// Invariants
/// ----------
/// - `windows` should be sorted by start (ascending order).
/// - Coordinates are half-open: `[start, end)`.
#[derive(Debug, Clone)]
pub struct Windows {
    pub windows: Vec<(u64, u64, u64)>, // (start, end, original_idx)
    /// Span start (inclusive) across all windows, as `i64`.
    /// This is the most-left coordinate covered by any of the windows.
    span_start: i64,
    /// Span end (exclusive) across all windows, as `i64`.
    /// This is the most-right coordinate covered by any of the windows.
    span_end: i64,
}

impl Windows {
    /// Construct from any window list (may be unsorted/overlapping).
    /// Ensures start- and end-sorted order (does not retain initial order)
    /// and computes span as `min(start)` .. `max(end)`.
    pub fn new(mut windows: Vec<(u64, u64, u64)>) -> Self {
        windows.sort_unstable_by_key(|w| (w.0, w.1));
        Windows::from_sorted(windows)
    }

    /// Construct from a list you guarantee is already sorted by start (non-decreasing).
    /// Computes span as `min(start)` .. `max(end)` (robust to irregular ends).
    pub fn from_sorted(windows: Vec<(u64, u64, u64)>) -> Self {
        debug_assert!(is_sorted_by_start(&windows), "windows must be start-sorted");
        let (span_start, span_end) = if windows.is_empty() {
            (0, 0)
        } else {
            let min_start = windows[0].0 as i64;
            let max_end = windows.iter().map(|w| w.1).max().unwrap() as i64;
            (min_start, max_end)
        };
        Self {
            windows,
            span_start,
            span_end,
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
    pub fn as_slice(&self) -> &[(u64, u64, u64)] {
        &self.windows
    }

    /// Consume and return the inner vector.
    #[inline]
    pub fn into_inner(self) -> Vec<(u64, u64, u64)> {
        self.windows
    }

    /// Span start (inclusive).
    /// This is the most-left coordinate covered by any of the windows.
    #[inline]
    pub fn span_start(&self) -> i64 {
        self.span_start
    }

    /// Span end (exclusive).
    /// This is the most-right coordinate covered by any of the windows.
    #[inline]
    pub fn span_end(&self) -> i64 {
        self.span_end
    }

    /// Span tuple `(start, end)`.
    /// These are the most-left and most-right coordinates covered by any of the windows.
    ///
    /// There are no guarantees that all positions between these two coordinates
    /// are covered by the windows.
    #[inline]
    pub fn span(&self) -> (i64, i64) {
        (self.span_start, self.span_end)
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
        let mut v = self.windows; // Take ownership; reuse allocation
        if v.is_empty() {
            return (
                Windows {
                    windows: v,
                    span_start: 0,
                    span_end: 0,
                },
                start_idx,
            );
        }

        debug_assert!(is_sorted_by_start(&v), "windows must be start-sorted");

        // In-place compaction with two indices: read cursor (i) and write cursor (w)
        let mut w: usize = 0;
        let mut cur_s = v[0].0;
        let mut cur_e = v[0].1;

        for i in 1..v.len() {
            let (s, e, _) = v[i];
            if s <= cur_e {
                if e > cur_e {
                    cur_e = e;
                }
            } else {
                // Write merged block at position w with new index
                v[w] = (cur_s, cur_e, start_idx + w as u64);
                w += 1;
                cur_s = s;
                cur_e = e;
            }
        }
        // Write the final block
        v[w] = (cur_s, cur_e, start_idx + w as u64);
        w += 1;

        // Shrink to the number of merged intervals
        v.truncate(w);

        // Since v is start-sorted and merged, first/last bound the span
        let span_start = v.first().map(|t| t.0 as i64).unwrap_or(0);
        let span_end = v.last().map(|t| t.1 as i64).unwrap_or(0);

        let next_idx = start_idx + w as u64;

        (
            Windows {
                windows: v,
                span_start,
                span_end,
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

#[inline]
fn is_sorted_by_start(ws: &[(u64, u64, u64)]) -> bool {
    ws.windows(2).all(|w| w[0].0 <= w[1].0)
}

#[inline]
fn is_sorted_by_start_with_scores(ws: &[(u64, u64, u64, f64)]) -> bool {
    ws.windows(2).all(|w| w[0].0 <= w[1].0)
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
///    
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, group index)'.
///
///  - Mapping of 'group index -> group name'.
pub fn load_grouped_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: &Vec<String>,
    filter_fn: Option<&dyn Fn(&str, u64, u64) -> bool>,
) -> Result<(FxHashMap<String, GroupedWindows>, FxHashMap<u64, String>)> {
    let mut reader = open_text_reader(bed.as_ref())?; // Works with &Path, PathBuf, &str

    // Pre-seed output map with requested chromosomes
    let mut vec_mapping: FxHashMap<&str, Vec<(u64, u64, u64)>> =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());
    for chr in chromosomes {
        vec_mapping.entry(chr.as_str()).or_default();
    }

    // Quick-hashing set of chromosomes to include
    let mut allowed_chromosomes: FxHashSet<&str> =
        FxHashSet::with_capacity_and_hasher(chromosomes.len(), Default::default());
    for chr in chromosomes {
        allowed_chromosomes.insert(chr.as_str());
    }

    // Enumeration of group names
    let mut group_name_to_idx: FxHashMap<String, u64> = FxHashMap::default();
    let mut next_group_idx: u64 = 0;

    // Reuse a single buffer for all lines
    let mut buf = String::new();
    let mut lineno: usize = 0;

    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
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
        let ls = line.trim_start();
        if ls.starts_with("track") || ls.starts_with("browser") {
            continue;
        }

        // Strict parse of first 3 BED columns without allocating a Vec
        let mut it = line.split_ascii_whitespace();

        let chr = match it.next() {
            Some(s) => s,
            None => continue, // or bail; here we skip blank/whitespace-only lines
        };
        if !allowed_chromosomes.contains(chr) {
            continue;
        }

        let start_str = it
            .next()
            .with_context(|| format!("BED parse error at line {}: missing start", lineno))?;
        let end_str = it
            .next()
            .with_context(|| format!("BED parse error at line {}: missing end", lineno))?;
        let group = it
            .next()
            .with_context(|| format!("BED parse error at line {}: missing group name", lineno))?;

        // Get group idx (enumerate and insert if first occurence)
        // We use this if/else approach only allocate a String once per unique group name
        let group_idx = if let Some(&i) = group_name_to_idx.get(group) {
            i
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

        // Apply passed filtering function
        if let Some(filterer) = filter_fn {
            if !filterer(chr, start, end) {
                continue;
            }
        }

        vec_mapping
            .get_mut(chr)
            .unwrap()
            .push((start, end, group_idx));
    }

    // Convert to Windows collections (Windows::new sorts internally)
    let windows_mapping: FxHashMap<String, GroupedWindows> = vec_mapping
        .into_iter()
        .map(|(chr, v)| (chr.to_string(), GroupedWindows::new(v)))
        .collect();

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
    pub windows: Vec<(u64, u64, u64)>, // (start, end, original_idx)
    /// Span start (inclusive) across all windows, as `i64`.
    /// This is the most-left coordinate covered by any of the windows.
    span_start: i64,
    /// Span end (exclusive) across all windows, as `i64`.
    /// This is the most-right coordinate covered by any of the windows.
    span_end: i64,
}

impl GroupedWindows {
    /// Construct from any window list (may be unsorted/overlapping).
    /// Ensures start- and end-sorted order (does not retain initial order)
    /// and computes span as `min(start)` .. `max(end)`.
    pub fn new(mut windows: Vec<(u64, u64, u64)>) -> Self {
        windows.sort_unstable_by_key(|w| (w.0, w.1));
        GroupedWindows::from_sorted(windows)
    }

    /// Construct from a list you guarantee is already sorted by start (non-decreasing).
    /// Computes span as `min(start)` .. `max(end)` (robust to irregular ends).
    pub fn from_sorted(windows: Vec<(u64, u64, u64)>) -> Self {
        debug_assert!(is_sorted_by_start(&windows), "windows must be start-sorted");
        let (span_start, span_end) = if windows.is_empty() {
            (0, 0)
        } else {
            let min_start = windows[0].0 as i64;
            let max_end = windows.iter().map(|w| w.1).max().unwrap() as i64;
            (min_start, max_end)
        };
        Self {
            windows,
            span_start,
            span_end,
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
    pub fn as_slice(&self) -> &[(u64, u64, u64)] {
        &self.windows
    }

    /// Consume and return the inner vector.
    #[inline]
    pub fn into_inner(self) -> Vec<(u64, u64, u64)> {
        self.windows
    }

    /// Span start (inclusive).
    /// This is the most-left coordinate covered by any of the windows.
    #[inline]
    pub fn span_start(&self) -> i64 {
        self.span_start
    }

    /// Span end (exclusive).
    /// This is the most-right coordinate covered by any of the windows.
    #[inline]
    pub fn span_end(&self) -> i64 {
        self.span_end
    }

    /// Span tuple `(start, end)`.
    /// These are the most-left and most-right coordinates covered by any of the windows.
    ///
    /// There are no guarantees that all positions between these two coordinates
    /// are covered by the windows.
    #[inline]
    pub fn span(&self) -> (i64, i64) {
        (self.span_start, self.span_end)
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
    let mut w = BufWriter::new(file);

    // Header
    writeln!(w, "group_idx\tgroup_name")
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
        writeln!(w, "{idx}\t{name}")
            .with_context(|| format!("Writing row for group_idx {idx} to {:?}", path))?;
    }

    Ok(())
}

/* Scored windows */

/// Load *scored* windows from a BED file into a per-chromosome map.
///
/// Parameters
/// ----------
///  - bed: Path to BED file with scores (float64) in the fourth column.
///  - chromosomes: Names of chromosomes to include in output,
///    even when not present in the BED file.
///  - filter_fn: Function for deciding whether to include
///    an interval. Should take in the `chr,start,end,score` values
///    and return `true` (keep) or `false` (discard).
///    
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, score, original index)'.
pub fn load_scored_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: Option<&[String]>,
    filter_fn: Option<&dyn Fn(&str, u64, u64, f64) -> bool>,
) -> Result<FxHashMap<String, ScoredWindows>> {
    let mut reader = open_text_reader(bed.as_ref())?; // Works with &Path, PathBuf, &str

    // Optional whitelist of chromosomes
    let mut vec_mapping: FxHashMap<String, Vec<(u64, u64, u64, f64)>> = FxHashMap::default();
    let allowed: Option<FxHashSet<&str>> = chromosomes.map(|chr_list| {
        let mut allowed = FxHashSet::with_capacity_and_hasher(chr_list.len(), Default::default());
        for chr in chr_list {
            allowed.insert(chr.as_str());
            vec_mapping.entry(chr.clone()).or_default();
        }
        allowed
    });

    // Reuse a single buffer for all lines
    let mut buf = String::new();
    let mut win_idx: u64 = 0;
    let mut lineno: usize = 0;

    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
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
        let ls = line.trim_start();
        if ls.starts_with("track") || ls.starts_with("browser") {
            continue;
        }

        // Strict parse of first 3 BED columns without allocating a Vec
        let mut it = line.split_ascii_whitespace();

        let chr = match it.next() {
            Some(s) => s,
            None => continue, // or bail; here we skip blank/whitespace-only lines
        };
        if let Some(allowed_chroms) = &allowed {
            if !allowed_chroms.contains(chr) {
                continue;
            }
        }

        let start_str = it
            .next()
            .with_context(|| format!("BED parse error at line {}: missing start", lineno))?;
        let end_str = it
            .next()
            .with_context(|| format!("BED parse error at line {}: missing end", lineno))?;
        let score_str = it
            .next()
            .with_context(|| format!("BED parse error at line {}: missing group name", lineno))?;

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

        let score: f64 = score_str.parse().with_context(|| {
            format!(
                "BED parse error at line {}: invalid score '{}'",
                lineno, score_str
            )
        })?;

        // Apply passed filtering function
        if let Some(filterer) = filter_fn {
            if !filterer(chr, start, end, score) {
                continue;
            }
        }

        vec_mapping
            .get_mut(chr)
            .unwrap()
            .push((start, end, win_idx, score));
        win_idx += 1;
    }

    // Convert to Windows collections (Windows::new sorts internally)
    let windows_mapping: FxHashMap<String, ScoredWindows> = vec_mapping
        .into_iter()
        .map(|(chr, v)| (chr.to_string(), ScoredWindows::new(v)))
        .collect();

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
    pub windows: Vec<(u64, u64, u64, f64)>, // (start, end, original_idx, score)
    /// Span start (inclusive) across all windows, as `i64`.
    /// This is the most-left coordinate covered by any of the windows.
    span_start: i64,
    /// Span end (exclusive) across all windows, as `i64`.
    /// This is the most-right coordinate covered by any of the windows.
    span_end: i64,
}

impl ScoredWindows {
    /// Construct from any window list (may be unsorted/overlapping).
    /// Ensures start- and end-sorted order (does not retain initial order)
    /// and computes span as `min(start)` .. `max(end)`.
    pub fn new(mut windows: Vec<(u64, u64, u64, f64)>) -> Self {
        windows.sort_unstable_by_key(|w| (w.0, w.1));
        ScoredWindows::from_sorted(windows)
    }

    /// Convert to Windows collection by dropping the score.
    pub fn to_windows(&self) -> Windows {
        Windows::from_sorted(
            self.windows
                .iter()
                .map(|(start, end, idx, _)| (*start, *end, *idx))
                .collect(),
        )
    }

    /// Construct from a list you guarantee is already sorted by start (non-decreasing).
    /// Computes span as `min(start)` .. `max(end)` (robust to irregular ends).
    pub fn from_sorted(windows: Vec<(u64, u64, u64, f64)>) -> Self {
        debug_assert!(
            is_sorted_by_start_with_scores(&windows),
            "windows must be start-sorted"
        );
        let (span_start, span_end) = if windows.is_empty() {
            (0, 0)
        } else {
            let min_start = windows[0].0 as i64;
            let max_end = windows.iter().map(|w| w.1).max().unwrap() as i64;
            (min_start, max_end)
        };
        Self {
            windows,
            span_start,
            span_end,
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
    pub fn as_slice(&self) -> &[(u64, u64, u64, f64)] {
        &self.windows
    }

    /// Consume and return the inner vector.
    #[inline]
    pub fn into_inner(self) -> Vec<(u64, u64, u64, f64)> {
        self.windows
    }

    /// Span start (inclusive).
    /// This is the most-left coordinate covered by any of the windows.
    #[inline]
    pub fn span_start(&self) -> i64 {
        self.span_start
    }

    /// Span end (exclusive).
    /// This is the most-right coordinate covered by any of the windows.
    #[inline]
    pub fn span_end(&self) -> i64 {
        self.span_end
    }

    /// Span tuple `(start, end)`.
    /// These are the most-left and most-right coordinates covered by any of the windows.
    ///
    /// There are no guarantees that all positions between these two coordinates
    /// are covered by the windows.
    #[inline]
    pub fn span(&self) -> (i64, i64) {
        (self.span_start, self.span_end)
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
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            return Ok(false);
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        return Ok(line_looks_like_header(&line, separator));
    }
}
