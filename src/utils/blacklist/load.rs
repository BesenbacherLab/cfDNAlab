use crate::utils::bed::load_windows_from_bed;
use anyhow::Result;
use std::{collections::HashMap, path::Path};

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
pub fn load_blacklists<P: AsRef<Path>>(
    beds: &[P],
    min_size: u64,
    chromosomes: &Vec<String>,
) -> Result<HashMap<String, Vec<(u64, u64)>>> {
    let mut merged: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    let filter_fn = move |_chr: &str, s: u64, e: u64| (e - s) >= min_size;

    // Load each BED file
    for bed in beds {
        let single = load_windows_from_bed(bed, chromosomes, Some(&filter_fn))?;
        for (chr, ivs) in single {
            let mut v: Vec<(u64, u64)> = ivs
                .into_inner()
                .into_iter()
                .map(|(s, e, _)| (s, e))
                .collect();
            merged.entry(chr).or_default().append(&mut v);
        }
    }

    // Sort and merge per chromosome
    for ivs in merged.values_mut() {
        ivs.sort_unstable();
        *ivs = merge_intervals(std::mem::take(ivs));
    }
    Ok(merged)
}

// TODO: Test, rename, and improve documentation!

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
