use crate::commands::prepare_windows::prepare_windows::FinalWindow;

/// Canonical sort orders used throughout the prepare_windows pipeline.
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub enum WindowSortOrder {
    /// Group-first ordering used by spacing/deduplication `(group, chrom, start, end)`.
    GroupChromStartEnd,
    /// Genomic ordering used for output `(chrom, start, end, group)`.
    ChromStartEndGroup,
}

/// Sort a slice of windows in place according to the requested ordering.
pub fn sort_windows_in_place(windows: &mut [FinalWindow], order: WindowSortOrder) {
    match order {
        WindowSortOrder::GroupChromStartEnd => windows.sort_unstable_by(|a, b| {
            a.group
                .cmp(&b.group)
                .then(a.chrom.cmp(&b.chrom))
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
        }),
        WindowSortOrder::ChromStartEndGroup => windows.sort_unstable_by(|a, b| {
            a.chrom
                .cmp(&b.chrom)
                .then(a.start.cmp(&b.start))
                .then(a.end.cmp(&b.end))
                .then(a.group.cmp(&b.group))
        }),
    }
}

/// Convenience helper that consumes and returns a sorted vector.
#[inline]
pub fn sort_windows_vec(mut windows: Vec<FinalWindow>, order: WindowSortOrder) -> Vec<FinalWindow> {
    sort_windows_in_place(&mut windows, order);
    windows
}
