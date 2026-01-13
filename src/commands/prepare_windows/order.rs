use crate::commands::prepare_windows::{config::CoordinateSet, prepare_windows::Window};
use rayon::prelude::*;

// Parallel sorts only pay off on larger batches, so keep a threshold.
const PARALLEL_SORT_MIN_LEN: usize = 50_000;

/// Canonical sort orders used throughout the prepare_windows pipeline.
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum WindowSortOrder {
    /// Group-first ordering used by minimum-distance filtering/deduplication `(group_key, chrom, start, end)`.
    GroupChromStartEnd,
    /// Genomic ordering that uses `group_key` as a tie-breaker `(chrom, start, end, group_key)`.
    ChromStartEndGroup,
    /// Genomic ordering that ignores `group_key` `(chrom, start, end)`.
    ChromStartEnd,
}

/// Sort a slice of windows in place according to the requested ordering.
pub fn sort_windows_in_place(
    windows: &mut [Window],
    order: WindowSortOrder,
    coord_set: CoordinateSet,
) {
    let use_parallel = windows.len() >= PARALLEL_SORT_MIN_LEN;

    match order {
        WindowSortOrder::GroupChromStartEnd => {
            let cmp = |a: &Window, b: &Window| {
                a.group_key
                    .cmp(&b.group_key)
                    .then(a.chrom.cmp(&b.chrom))
                    .then(a.start_for(coord_set).cmp(&b.start_for(coord_set)))
                    .then(a.end_for(coord_set).cmp(&b.end_for(coord_set)))
            };
            if use_parallel {
                windows.par_sort_unstable_by(cmp);
            } else {
                windows.sort_unstable_by(cmp);
            }
        }
        WindowSortOrder::ChromStartEndGroup => {
            let cmp = |a: &Window, b: &Window| {
                a.chrom
                    .cmp(&b.chrom)
                    .then(a.start_for(coord_set).cmp(&b.start_for(coord_set)))
                    .then(a.end_for(coord_set).cmp(&b.end_for(coord_set)))
                    .then(a.group_key.cmp(&b.group_key))
            };
            if use_parallel {
                windows.par_sort_unstable_by(cmp);
            } else {
                windows.sort_unstable_by(cmp);
            }
        }
        WindowSortOrder::ChromStartEnd => {
            let cmp = |a: &Window, b: &Window| {
                a.chrom
                    .cmp(&b.chrom)
                    .then(a.start_for(coord_set).cmp(&b.start_for(coord_set)))
                    .then(a.end_for(coord_set).cmp(&b.end_for(coord_set)))
            };
            if use_parallel {
                windows.par_sort_unstable_by(cmp);
            } else {
                windows.sort_unstable_by(cmp);
            }
        }
    }
}

/// Convenience helper that consumes and returns a sorted vector.
#[inline]
pub fn sort_windows_vec(mut windows: Vec<Window>, order: WindowSortOrder) -> Vec<Window> {
    sort_windows_in_place(&mut windows, order, CoordinateSet::Resized);
    windows
}
