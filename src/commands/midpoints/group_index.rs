use crate::shared::{interval::IndexedInterval, io::create_text_writer};
use anyhow::{Context, Result};
use fxhash::FxHashMap;
use std::{io::Write, path::Path};

/// Count profile-eligible input intervals per group.
///
/// The midpoint profile tensor is summed over intervals. This metadata records how many intervals
/// remain in each group after chromosome filtering and interval-level blacklist prefiltering, even
/// when no fragment midpoint later overlaps an interval.
pub(crate) fn eligible_interval_counts_by_group(
    indexed_intervals_by_chromosome: &FxHashMap<String, Vec<IndexedInterval<u64>>>,
    group_idx_to_name: &FxHashMap<u64, String>,
) -> FxHashMap<u64, usize> {
    let mut counts: FxHashMap<u64, usize> = group_idx_to_name
        .keys()
        .map(|&group_idx| (group_idx, 0usize))
        .collect();
    for intervals in indexed_intervals_by_chromosome.values() {
        for interval in intervals {
            *counts.entry(interval.idx()).or_insert(0) += 1;
        }
    }
    counts
}

/// Write midpoint group metadata next to the profile array.
///
/// `eligible_intervals` is the number of output intervals retained in each group after
/// command-level interval filtering. Users can divide summed group profiles by this count when it
/// is nonzero and they want mean profiles over eligible intervals.
pub(crate) fn write_midpoint_group_index_tsv(
    output_path: &Path,
    group_idx_to_name: &FxHashMap<u64, String>,
    eligible_interval_counts: &FxHashMap<u64, usize>,
) -> Result<()> {
    let mut writer = create_text_writer(output_path)
        .with_context(|| format!("creating midpoint group index {}", output_path.display()))?;
    writeln!(writer, "group_idx\tgroup_name\teligible_intervals").with_context(|| {
        format!(
            "writing midpoint group index header to {}",
            output_path.display()
        )
    })?;

    let mut entries: Vec<(u64, &str)> = group_idx_to_name
        .iter()
        .map(|(idx, name)| (*idx, name.as_str()))
        .collect();
    entries.sort_unstable_by_key(|(idx, _)| *idx);

    for (group_idx, group_name) in entries {
        let group_name = group_name.replace('\t', "    ").replace('\n', " ");
        let eligible_intervals = eligible_interval_counts
            .get(&group_idx)
            .copied()
            .unwrap_or(0);
        writeln!(writer, "{group_idx}\t{group_name}\t{eligible_intervals}").with_context(|| {
            format!(
                "writing midpoint group index row for group_idx {} to {}",
                group_idx,
                output_path.display()
            )
        })?;
    }

    writer
        .finish()
        .with_context(|| format!("finalizing midpoint group index {}", output_path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("group_index_tests.rs");
}
