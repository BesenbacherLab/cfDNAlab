use crate::shared::bed::{Windows, load_windows_from_bed};
use anyhow::Result;
use fxhash::FxHashMap;
use std::path::Path;

/// Load blacklisted genomic intervals from *one or more* BED files.
///
/// Intervals shorter than `min_size` are discarded. When `halo_bp > 0`, each
/// interval is expanded by that halo on both sides *before* merging, ensuring
/// the coalescing step accounts for the expanded span.
///
/// # Parameters
/// - `beds`: BED files to read.
/// - `min_size`: Minimum interval length (bp) to keep.
/// - `halo_bp`: Halo to expand on both sides before merging.
/// - `chromosomes`: Optional whitelist of chromosome names to retain.
///
/// TODO: plumb through a streaming reader so we do not materialise every
/// interval before merging when very large inputs are used.
pub fn load_blacklists<P: AsRef<Path>>(
    beds: &[P],
    min_size: u64,
    halo_bp: u64,
    chromosomes: Option<&[String]>,
) -> Result<FxHashMap<String, Vec<(u64, u64)>>> {
    if beds.is_empty() {
        return Ok(FxHashMap::default());
    }

    let mut merged: FxHashMap<String, Vec<(u64, u64)>> = FxHashMap::default();

    for bed in beds {
        let single = load_windows_from_bed(bed, chromosomes, None)?;
        accumulate_blacklist_windows(&mut merged, single, min_size, halo_bp);
    }

    for ivs in merged.values_mut() {
        ivs.sort_unstable();
        *ivs = merge_intervals(std::mem::take(ivs));
    }

    Ok(merged)
}

pub fn merge_intervals(ivs: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    if ivs.is_empty() {
        return ivs;
    }
    let mut merged = Vec::with_capacity(ivs.len());
    let mut cur = ivs[0];
    for (s, e) in ivs.into_iter().skip(1) {
        if s <= cur.1 {
            cur.1 = cur.1.max(e);
        } else {
            merged.push(cur);
            cur = (s, e);
        }
    }
    merged.push(cur);
    merged
}

fn accumulate_blacklist_windows(
    merged: &mut FxHashMap<String, Vec<(u64, u64)>>,
    windows_map: FxHashMap<String, Windows>,
    min_size: u64,
    halo_bp: u64,
) {
    for (chr, ivs) in windows_map {
        let mut out: Vec<(u64, u64)> = Vec::new();
        for (start, end, _) in ivs.into_inner() {
            if end <= start {
                continue;
            }
            if (end - start) < min_size {
                continue;
            }
            let halo_start = start.saturating_sub(halo_bp);
            let halo_end = end.saturating_add(halo_bp);
            out.push((halo_start, halo_end));
        }
        if !out.is_empty() {
            merged.entry(chr).or_default().extend(out);
        }
    }
}
