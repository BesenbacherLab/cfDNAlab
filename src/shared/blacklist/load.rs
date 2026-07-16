use crate::shared::bed::{Windows, load_windows_from_bed};
use crate::shared::interval::Interval;
use anyhow::Result as AnyResult;
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
/// - `read_in_background`: Whether file reading and decompression should overlap with BED parsing.
///
/// TODO: plumb through a streaming reader so we do not materialise every
/// interval before merging when very large inputs are used.
pub fn load_blacklists<P: AsRef<Path>>(
    beds: &[P],
    min_size: u64,
    halo_bp: u64,
    chromosomes: Option<&[String]>,
    read_in_background: bool,
) -> AnyResult<FxHashMap<String, Vec<Interval<u64>>>> {
    if beds.is_empty() {
        return Ok(FxHashMap::default());
    }

    let mut merged: FxHashMap<String, Vec<Interval<u64>>> = FxHashMap::default();

    for bed in beds {
        let single = load_windows_from_bed(bed, chromosomes, None, None, read_in_background)?;
        accumulate_blacklist_windows(&mut merged, single, min_size, halo_bp);
    }

    for ivs in merged.values_mut() {
        let intervals = std::mem::take(ivs);
        *ivs = merge_intervals(intervals)?;
    }

    Ok(merged)
}

/// Merge touching or overlapping half-open intervals.
///
/// The input may be unsorted. The returned intervals are sorted by start and
/// coalesced so that adjacent intervals with `end == next.start` become one
/// interval.
pub fn merge_intervals(intervals: Vec<Interval<u64>>) -> crate::Result<Vec<Interval<u64>>> {
    if intervals.is_empty() {
        return Ok(intervals);
    }

    let mut intervals = intervals;
    intervals.sort_unstable_by_key(|interval| (interval.start(), interval.end()));

    let mut merged = Vec::with_capacity(intervals.len());
    let mut current = intervals[0];

    for interval in intervals.into_iter().skip(1) {
        if interval.start() <= current.end() {
            current = Interval::new(current.start(), current.end().max(interval.end()))?;
        } else {
            merged.push(current);
            current = interval;
        }
    }
    merged.push(current);
    Ok(merged)
}

fn accumulate_blacklist_windows(
    merged: &mut FxHashMap<String, Vec<Interval<u64>>>,
    windows_map: FxHashMap<String, Windows>,
    min_size: u64,
    halo_bp: u64,
) {
    for (chr, ivs) in windows_map {
        let mut out: Vec<Interval<u64>> = Vec::new();
        for window in ivs.into_inner() {
            let start = window.start();
            let end = window.end();
            if end <= start {
                continue;
            }
            if (end - start) < min_size {
                continue;
            }
            let halo_start = start.saturating_sub(halo_bp);
            let halo_end = end.saturating_add(halo_bp);
            out.push(
                Interval::new(halo_start, halo_end)
                    .expect("blacklist halo expansion must preserve non-empty intervals"),
            );
        }
        if !out.is_empty() {
            merged.entry(chr).or_default().extend(out);
        }
    }
}

#[cfg(test)]
mod tests {
    include!("load_tests.rs");
}
