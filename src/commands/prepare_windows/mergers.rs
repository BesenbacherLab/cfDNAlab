use crate::commands::prepare_windows::{
    config::{MergeLabel, MergeScope},
    order::{WindowSortOrder, sort_windows_in_place},
    prepare_windows::FinalWindow,
};

/// Merge windows that are separated by at most `merge_gap_bp`.
///
/// When merging within a group, process windows grouped by `group`.
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
///     Policy for merged group label.
/// - presorted:
///     Whether the output is already sorted properly for the `merge_scope` (by coordinate or by group).
///     When `merge_scope=within`, the windows should be sortedby (group, chrom, start, end).
///     When `merge_scope=across`, the windows should be sortedby (chrom, start, end, group).
///
/// Returns
/// -------
/// - merged:
///     Merged windows with group labels composed according to policy.
pub fn merge_windows(
    windows: Vec<FinalWindow>,
    merge_scope: MergeScope,
    merge_gap_bp: Option<u32>,
    merge_label: MergeLabel,
    presorted: bool,
) -> Vec<FinalWindow> {
    if merge_gap_bp.is_none() || windows.is_empty() || matches!(merge_scope, MergeScope::None) {
        return windows;
    }
    let gap = merge_gap_bp.unwrap();

    match merge_scope {
        MergeScope::None => unreachable!(),
        MergeScope::Within => merge_within_groups(windows, gap, merge_label, presorted),
        MergeScope::Across => merge_across_groups(windows, gap, merge_label, presorted),
    }
}

pub fn merge_within_groups(
    mut windows: Vec<FinalWindow>,
    gap: u32,
    merge_label: MergeLabel,
    presorted: bool,
) -> Vec<FinalWindow> {
    if !presorted {
        sort_windows_in_place(&mut windows, WindowSortOrder::GroupChromStartEnd);
    }
    let mut result: Vec<FinalWindow> = Vec::with_capacity(windows.len());
    let mut i = 0usize;
    while i < windows.len() {
        let mut current = windows[i].clone();
        let mut merged = false;
        i += 1;
        while i < windows.len()
            && windows[i].group == current.group
            && windows[i].chrom == current.chrom
            && windows[i].start <= current.end.saturating_add(gap)
        {
            merged = true;
            if windows[i].end > current.end {
                current.end = windows[i].end;
            }
            if let MergeLabel::Join = merge_label {
                if !windows[i].group.is_empty() && windows[i].group != current.group {
                    current.group = format!("{}__{}", current.group, windows[i].group);
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

pub fn merge_across_groups(
    mut windows: Vec<FinalWindow>,
    gap: u32,
    merge_label: MergeLabel,
    presorted: bool,
) -> Vec<FinalWindow> {
    if !presorted {
        sort_windows_in_place(&mut windows, WindowSortOrder::ChromStartEndGroup);
    }

    let mut result: Vec<FinalWindow> = Vec::with_capacity(windows.len());
    let mut i = 0usize;
    while i < windows.len() {
        let mut current = windows[i].clone();
        let mut merged = false;
        i += 1;
        while i < windows.len()
            && windows[i].chrom == current.chrom
            && windows[i].start <= current.end.saturating_add(gap)
        {
            merged = true;
            if windows[i].end > current.end {
                current.end = windows[i].end;
            }
            match merge_label {
                MergeLabel::Join => {
                    if current.group.is_empty() {
                        current.group = windows[i].group.clone();
                    } else if !windows[i].group.is_empty()
                        && !current.group.contains(&windows[i].group)
                    {
                        current.group = format!("{}__{}", current.group, windows[i].group);
                    }
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
