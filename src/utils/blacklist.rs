use anyhow::{Context, Result};
use std::{collections::HashMap, path::PathBuf};

/// Load blacklist intervals into a `HashMap` keyed by chromosome name.
///
/// * Uses **only** the first three columns (`chrom`, `start`, `end`) and
///   ignores any additional BED fields.
/// * Lines that begin with `#`, `track`, `browser`, or are blank are skipped.
/// * `chromosomes` is usually the autosome whitelist (e.g. `["chr1", … "chr22"]`).
pub fn load_blacklist(
    bed: &PathBuf,
    min_size: u64,
    chromosomes: &Vec<String>,
) -> Result<HashMap<String, Vec<(u64, u64)>>> {
    // Create a map from chromosome name to its blacklist intervals
    let mut map: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    let content =
        std::fs::read_to_string(bed).context(format!("Error reading blacklist BED {:?}", bed))?;
    for line in content.lines().map(str::trim) {
        // Skip comments, headers, empty lines
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        // Take only the first three whitespace-separated fields
        let mut fields = line.split_whitespace();
        let chr = match fields.next() {
            Some(c) => c.to_string(),
            None => continue, // Malformed line
        };
        // Skip non-autosomes
        if !chromosomes.contains(&chr) {
            continue;
        }
        // Parse start and end; skip line if either fails
        let start: u64 = match fields.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue, // non-numeric or missing
        };
        let end: u64 = match fields.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue, // non-numeric or missing
        };
        // Keep interval if length ≥ min_size
        if end > start && (end - start) >= min_size {
            map.entry(chr.clone()).or_default().push((start, end));
        }
    }
    // Sort intervals for each chromosome
    for iv in map.values_mut() {
        iv.sort_unstable();
    }

    Ok(map)
}

// TODO: Test properly loaded, concatenated and sorted

/// Load blacklisted genomic intervals from *one or more* BED files.
///
/// Overlapping intervals are merged (flattened).
///
/// Returns sorted, flattened intervals per chromosome.
/// 
/// * Uses **only** the first three columns (`chrom`, `start`, `end`) and
///   ignores any additional BED fields.
/// * Lines that begin with `#`, `track`, `browser`, or are blank are skipped.
/// * `chromosomes` is usually the autosome whitelist (e.g. `["chr1", … "chr22"]`).
pub fn load_blacklists(
    beds: &[PathBuf],
    min_size: u64,
    chromosomes: &Vec<String>,
) -> Result<HashMap<String, Vec<(u64, u64)>>> {
    let mut merged: HashMap<String, Vec<(u64, u64)>> = HashMap::new();

    // Load each BED file
    for bed in beds {
        let single = load_blacklist(bed, min_size, chromosomes)?;
        for (chr, mut ivs) in single {
            merged.entry(chr).or_default().append(&mut ivs);
        }
    }

    // Sort and merge per chromosome
    for ivs in merged.values_mut() {
        ivs.sort_unstable();
        *ivs = merge_intervals(std::mem::take(ivs));
    }
    Ok(merged)
}

/// Check if the full fragment lies within an interval
pub fn is_full(intervals: &[(u64, u64)], start: u64, end: u64, ptr: &mut usize) -> bool {
    // Skip any intervals that end entirely before our fragment start
    while *ptr < intervals.len() && intervals[*ptr].1 <= start {
        *ptr += 1;
    }
    // If there's an interval here, check full overlap of [start,end)
    if let Some(&(s, e)) = intervals.get(*ptr) {
        s <= start && e >= end
    } else {
        false
    }
}

// TODO: Test, rename, and improve documentation!

/// Advance `ptr` to skip any intervals ending before `start`, then
/// sum up how many bases of [start,end) overlap the intervals.
/// `ptr` is left of the first interval that might overlap the next bin.
///
/// Intervals must be sorted by start position and be non‐overlapping.
pub fn compute_blacklist_overlap(
    intervals: &[(u64, u64)],
    start: u64,
    end: u64,
    ptr: &mut usize,
) -> f64 {
    // Skip intervals that end at or before the bin start
    while *ptr < intervals.len() && intervals[*ptr].1 <= start {
        *ptr += 1;
    }
    // Sum all overlap lengths for this bin
    let mut covered = 0;
    let mut i = *ptr;
    while i < intervals.len() && intervals[i].0 < end {
        let (s, e) = intervals[i];
        covered += e.min(end).saturating_sub(s.max(start));
        i += 1;
    }
    covered as f64 / (end - start) as f64
}

/// Merge intervals when they touch or overlap.
///
/// * ivs: Intervals sorted by start and end positions.
pub fn merge_intervals(ivs: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    if ivs.is_empty() {
        return ivs;
    }
    // Already sorted by caller (`sort_unstable`)
    let mut merged = Vec::with_capacity(ivs.len());
    let mut cur = ivs[0];
    // Find
    for (s, e) in ivs.into_iter().skip(1) {
        if s <= cur.1 {
            // overlap or touch
            cur.1 = cur.1.max(e); // extend the current block
        } else {
            merged.push(cur);
            cur = (s, e);
        }
    }
    merged.push(cur);
    merged
}

// -- Ref sequence position blacklisting --

/// Byte used for blacklisted bases in the reference sequence
pub const BLACKLIST_BYTE: u8 = b'X';

/// Mask every base that falls inside a blacklist interval with `BLACKLIST_BYTE`.
///
/// * `seq`: mutable byte slice of the reference chromosome  
/// * `intervals`: merged, **sorted**, non-overlapping `[start, end)` pairs  
///
/// Runs in **O(total interval length)** – no per-base scanning.
pub fn apply_blacklist_mask_to_seq(seq: &mut [u8], intervals: &[(u64, u64)]) {
    for &(start, end) in intervals {
        let s = start as usize;
        let e = end as usize;
        // Silent bounds-check: some BEDs can extend past chromosome end
        if s >= seq.len() {
            break;
        }
        let e = e.min(seq.len());
        seq[s..e].fill(BLACKLIST_BYTE);
    }
}
