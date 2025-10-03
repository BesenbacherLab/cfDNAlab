use crate::commands::prepare_windows::{
    config::{MergeLabel, MergeScope},
    prepare_windows::FinalWindow,
};

/// Merge windows that are separated by at most `merge_gap_bp`.
///
/// When merging within a group, process windows grouped by `group`.
/// When merging across groups, process all windows together (merging labels as configured).
///
/// Parameters
/// ----------
/// - windows:
///     Windows sorted by (group, chrom, start, end) if `merge_scope=within`,
///     or by (chrom, start, end, group) if `merge_scope=across`.
/// - merge_scope:
///     Whether to merge within groups or across groups.
/// - merge_gap_bp:
///     Merge gap threshold in base pairs (None to disable).
/// - merge_label:
///     Policy for merged group label.
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
) -> Vec<FinalWindow> {
    if merge_gap_bp.is_none() || windows.is_empty() {
        return windows;
    }
    let gap = merge_gap_bp.unwrap();

    match merge_scope {
        MergeScope::None => windows,
        MergeScope::Within => merge_within_groups(windows, gap, merge_label),
        MergeScope::Across => merge_across_groups(windows, gap, merge_label),
    }
}

pub fn merge_within_groups(
    windows: Vec<FinalWindow>,
    gap: u32,
    merge_label: MergeLabel,
) -> Vec<FinalWindow> {
    // Windows are expected sorted by (group, chrom, start, end)
    let mut result: Vec<FinalWindow> = Vec::with_capacity(windows.len());
    let mut i = 0usize;
    while i < windows.len() {
        let mut current = windows[i].clone();
        i += 1;
        while i < windows.len()
            && windows[i].group == current.group
            && windows[i].chrom == current.chrom
            && windows[i].start <= current.end.saturating_add(gap)
        {
            // Merge by extending end and possibly adjusting label
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
        result.push(current);
    }
    result
}

pub fn merge_across_groups(
    mut windows: Vec<FinalWindow>,
    gap: u32,
    merge_label: MergeLabel,
) -> Vec<FinalWindow> {
    // Sort by (chrom, start, end, group) to merge across groups deterministically
    windows.sort_unstable_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.group.cmp(&b.group))
    });

    let mut result: Vec<FinalWindow> = Vec::with_capacity(windows.len());
    let mut i = 0usize;
    while i < windows.len() {
        let mut current = windows[i].clone();
        i += 1;
        while i < windows.len()
            && windows[i].chrom == current.chrom
            && windows[i].start <= current.end.saturating_add(gap)
        {
            // Merge interval span
            if windows[i].end > current.end {
                current.end = windows[i].end;
            }
            // Merge label if requested
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
        result.push(current);
    }
    result
}
