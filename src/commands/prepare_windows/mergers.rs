use crate::commands::prepare_windows::{
    config::{CoordinateSet, MergeLabel, MergeScope},
    labels::normalize_label_tuples,
    order::{WindowSortOrder, sort_windows_in_place},
    prepare_windows::Window,
};

/// Merge windows that are separated by at most `merge_gap_bp`.
///
/// When merging within a group, process windows grouped by `group_key`.
///
/// When merging across groups, process all windows together (merging labels as configured).
///
/// Parameters
/// ----------
/// - windows:
///     Windows to merge.
/// - merge_scope:
///     Whether to merge within groups or across groups.
/// - merge_gap_bp:
///     Merge gap threshold in base pairs (None to disable).
/// - merge_label:
///     Policy for merged label tuples.
/// - merge_on:
///     Coordinate set used for merge comparisons.
/// - presorted:
///     Whether the output is already sorted properly for the `merge_scope` (by coordinate or by group).
///     When `merge_scope=within`, the windows should be sorted by (group_key, chrom, start, end).
///     When `merge_scope=across`, the windows should be sorted by (chrom, start, end, group_key).
///
/// Returns
/// -------
/// - merged:
///     Merged windows with label tuples composed according to policy.
pub(crate) fn merge_windows(
    windows: Vec<Window>,
    merge_scope: MergeScope,
    merge_gap_bp: Option<u32>,
    merge_label: MergeLabel,
    merge_on: CoordinateSet,
    presorted: bool,
) -> Vec<Window> {
    if merge_gap_bp.is_none() || windows.is_empty() || matches!(merge_scope, MergeScope::None) {
        return windows;
    }
    let gap = merge_gap_bp.unwrap();

    match merge_scope {
        MergeScope::None => unreachable!(),
        MergeScope::Within => merge_within_groups(windows, gap, merge_label, merge_on, presorted),
        MergeScope::Across => merge_across_groups(windows, gap, merge_label, merge_on, presorted),
    }
}

pub(crate) fn merge_within_groups(
    mut windows: Vec<Window>,
    gap: u32,
    merge_label: MergeLabel,
    merge_on: CoordinateSet,
    presorted: bool,
) -> Vec<Window> {
    if !presorted {
        sort_windows_in_place(&mut windows, WindowSortOrder::GroupChromStartEnd, merge_on);
    }
    let mut result: Vec<Window> = Vec::with_capacity(windows.len());
    let mut i = 0usize;
    while i < windows.len() {
        let mut current = windows[i].clone();
        let mut merged = false;
        i += 1;
        while i < windows.len()
            && windows[i].group_key == current.group_key
            && windows[i].chrom == current.chrom
            && windows[i].start_for(merge_on) <= current.end_for(merge_on).saturating_add(gap)
        {
            merged = true;
            current.expand_intervals_to_include(&windows[i]);
            if let MergeLabel::Join = merge_label {
                current
                    .label_tuples
                    .extend(windows[i].label_tuples.iter().cloned());
                normalize_label_tuples(&mut current.label_tuples);
            }
            i += 1;
        }
        if merged {
            current.merged = true;
        }
        result.push(current);
    }
    result
}

pub(crate) fn merge_across_groups(
    mut windows: Vec<Window>,
    gap: u32,
    merge_label: MergeLabel,
    merge_on: CoordinateSet,
    presorted: bool,
) -> Vec<Window> {
    if !presorted {
        sort_windows_in_place(&mut windows, WindowSortOrder::ChromStartEndGroup, merge_on);
    }

    let mut result: Vec<Window> = Vec::with_capacity(windows.len());
    let mut i = 0usize;
    while i < windows.len() {
        let mut current = windows[i].clone();
        let mut merged = false;
        i += 1;
        while i < windows.len()
            && windows[i].chrom == current.chrom
            && windows[i].start_for(merge_on) <= current.end_for(merge_on).saturating_add(gap)
        {
            merged = true;
            current.expand_intervals_to_include(&windows[i]);
            match merge_label {
                MergeLabel::Join => {
                    current
                        .label_tuples
                        .extend(windows[i].label_tuples.iter().cloned());
                    normalize_label_tuples(&mut current.label_tuples);
                }
                MergeLabel::First => {
                    // Keep first group's label
                }
            }
            i += 1;
        }
        if merged {
            current.merged = true;
        }
        result.push(current);
    }
    result
}

#[cfg(test)]
mod tests {
    include!("mergers_tests.rs");
}
