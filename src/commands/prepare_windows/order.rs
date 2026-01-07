use crate::commands::prepare_windows::{config::CoordinateSet, prepare_windows::FinalWindow};

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
    windows: &mut [FinalWindow],
    order: WindowSortOrder,
    coord_set: CoordinateSet,
) {
    match order {
        WindowSortOrder::GroupChromStartEnd => windows.sort_unstable_by(|a, b| {
            a.group_key
                .cmp(&b.group_key)
                .then(a.chrom.cmp(&b.chrom))
                .then(a.start_for(coord_set).cmp(&b.start_for(coord_set)))
                .then(a.end_for(coord_set).cmp(&b.end_for(coord_set)))
        }),
        WindowSortOrder::ChromStartEndGroup => windows.sort_unstable_by(|a, b| {
            a.chrom
                .cmp(&b.chrom)
                .then(a.start_for(coord_set).cmp(&b.start_for(coord_set)))
                .then(a.end_for(coord_set).cmp(&b.end_for(coord_set)))
                .then(a.group_key.cmp(&b.group_key))
        }),
        WindowSortOrder::ChromStartEnd => windows.sort_unstable_by(|a, b| {
            a.chrom
                .cmp(&b.chrom)
                .then(a.start_for(coord_set).cmp(&b.start_for(coord_set)))
                .then(a.end_for(coord_set).cmp(&b.end_for(coord_set)))
        }),
    }
}

/// Convenience helper that consumes and returns a sorted vector.
#[inline]
pub fn sort_windows_vec(mut windows: Vec<FinalWindow>, order: WindowSortOrder) -> Vec<FinalWindow> {
    sort_windows_in_place(&mut windows, order, CoordinateSet::Resized);
    windows
}
