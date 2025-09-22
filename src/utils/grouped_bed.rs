use anyhow::{Context, Result, ensure};
use fxhash::{FxHashMap, FxHashSet};
use std::{
    fs::{File, create_dir_all},
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

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
    let f = File::open(bed.as_ref()).context("Opening BED file with windows/intervals")?; // Works with &Path, PathBuf, &str
    let mut reader = BufReader::with_capacity(1 << 20, f);

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
    windows: Vec<(u64, u64, u64)>, // (start, end, original_idx)
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

#[inline]
fn is_sorted_by_start(ws: &[(u64, u64, u64)]) -> bool {
    ws.windows(2).all(|w| w[0].0 <= w[1].0)
}

/// Get window length and ensure it's the same for ALL windows.
pub fn ensure_uniform_window_len(
    windows_by_chr: &FxHashMap<String, GroupedWindows>,
) -> Result<usize> {
    let mut reference_len: Option<usize> = None;

    for (chr, gw) in windows_by_chr {
        for (start, end, _) in &gw.windows {
            let len = end.checked_sub(*start).with_context(|| {
                format!("Invalid window on {chr}: end ({end}) < start ({start})")
            })? as usize;

            match reference_len {
                None => reference_len = Some(len),
                Some(ref_len) if (len) != ref_len => {
                    anyhow::bail!(
                        "Non-uniform window length detected on {chr}: [{start},{end}) has len {}, expected {}",
                        len,
                        ref_len
                    );
                }
                _ => {}
            }
        }
    }

    reference_len.context("No windows found when checking uniform window length")
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
