use anyhow::{Context, Result, ensure};
use fxhash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::{
    collections::HashMap,
    io::{BufRead, BufReader},
    path::Path,
};

/// Load windows from a BED file into a per-chromosome map.
///
/// Parameters
/// ----------
///  - bed: Path to BED file.
///  - chromosomes: Names of chromosomes to include in output,
///    even when not present in the BED file.
///
/// Returns
/// -------
///  - Mapping of 'chromosome -> sorted window coordinates (start, end, original window index)'.
pub fn load_windows_from_bed(
    bed: impl AsRef<Path>,
    chromosomes: &Vec<String>,
) -> Result<HashMap<String, Windows>> {
    let f = File::open(bed.as_ref()).context("Opening BED file with windows/intervals")?; // Works with &Path, PathBuf, &str
    let mut reader = BufReader::with_capacity(1 << 20, f);

    // Pre-seed output map with requested chromosomes.
    // let mut vec_mapping: HashMap<String, Vec<(u64, u64, u64)>> =
    //     HashMap::with_capacity(chromosomes.len());
    let mut vec_mapping: FxHashMap<&str, Vec<(u64, u64, u64)>> =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());

    for chr in chromosomes {
        vec_mapping.entry(chr.as_str()).or_default();
    }

    // O(1) membership checks without per-line allocation.
    let allowed_chromosomes: FxHashSet<&str> =
        FxHashSet::with_capacity_and_hasher(chromosomes.len(), Default::default());
    // let allowed_chromosomes: HashSet<&str> = chromosomes.iter().map(String::as_str).collect();

    // Reuse a single buffer for all lines.
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

        // Fast skips.
        if buf.as_bytes().first().is_some_and(|b| *b == b'#') {
            continue;
        }
        let line = buf.trim_end_matches(['\n', '\r']);

        if line.is_empty() {
            continue;
        }

        // Strict parse of first 3 BED columns without allocating a Vec.
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

        vec_mapping
            .get_mut(chr)
            .unwrap()
            .push((start, end, win_idx));
        win_idx += 1;
    }

    // Convert to Windows collections (Windows::new sorts internally).
    let windows_mapping: HashMap<String, Windows> = vec_mapping
        .into_iter()
        .map(|(chr, v)| (chr.to_string(), Windows::new(v)))
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
    windows: Vec<(u64, u64, u64)>, // (start, end, original_idx)
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
}

#[inline]
fn is_sorted_by_start(ws: &[(u64, u64, u64)]) -> bool {
    ws.windows(2).all(|w| w[0].0 <= w[1].0)
}
