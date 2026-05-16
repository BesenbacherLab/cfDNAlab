use crate::shared::{bed::GroupedWindows, io::create_text_writer};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use std::{io::Write, path::Path};

/// Count profile-eligible input intervals per group.
///
/// The midpoint profile tensor is summed over intervals. This metadata records how many intervals
/// remain in each group after chromosome filtering and interval-level blacklist prefiltering, even
/// when no fragment midpoint later overlaps an interval.
pub(crate) fn eligible_interval_counts_by_group(
    grouped_windows_by_chromosome: &FxHashMap<String, GroupedWindows>,
    group_idx_to_name: &FxHashMap<u64, String>,
) -> FxHashMap<u64, usize> {
    let mut counts: FxHashMap<u64, usize> = group_idx_to_name
        .keys()
        .map(|&group_idx| (group_idx, 0usize))
        .collect();
    for grouped_windows in grouped_windows_by_chromosome.values() {
        for interval in grouped_windows.windows.iter() {
            *counts.entry(interval.idx()).or_insert(0) += 1;
        }
    }
    counts
}

/// Ordered midpoint group metadata shared by all midpoint outputs.
///
/// Midpoint counts use `group_idx` directly as the row index in the profile tensor. The
/// group-index TSV and Zarr group coordinate must therefore describe groups in exactly the same
/// order. This summary is the single source of that ordering.
#[derive(Debug)]
pub(crate) struct MidpointGroupSummary<'a> {
    pub(crate) group_idx: u64,
    pub(crate) group_name: &'a str,
    pub(crate) eligible_intervals: usize,
}

/// Return midpoint group metadata in count-row order.
///
/// The BED loader assigns `group_idx` values that are expected to be contiguous from zero. This
/// function walks the expected count rows directly, then rejects gaps so callers cannot silently
/// reinterpret count rows. Missing eligible interval counts default to zero because a group can
/// have no retained intervals after filtering.
pub(crate) fn ordered_midpoint_group_summaries<'a>(
    group_idx_to_name: &'a FxHashMap<u64, String>,
    eligible_interval_counts: &FxHashMap<u64, usize>,
) -> Result<Vec<MidpointGroupSummary<'a>>> {
    let mut summaries = Vec::with_capacity(group_idx_to_name.len());
    for row_index in 0..group_idx_to_name.len() {
        let group_idx = row_index as u64;
        let Some(group_name) = group_idx_to_name.get(&group_idx) else {
            let observed_indices = sorted_group_indices(group_idx_to_name);
            bail!(
                "midpoint group indices must match count rows 0..{} but observed {:?}",
                group_idx_to_name.len().saturating_sub(1),
                observed_indices
            );
        };
        summaries.push(MidpointGroupSummary {
            group_idx,
            group_name: group_name.as_str(),
            eligible_intervals: eligible_interval_counts
                .get(&group_idx)
                .copied()
                .unwrap_or(0),
        });
    }

    Ok(summaries)
}

fn sorted_group_indices(group_idx_to_name: &FxHashMap<u64, String>) -> Vec<u64> {
    let mut indices: Vec<u64> = group_idx_to_name.keys().copied().collect();
    indices.sort_unstable();
    indices
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

    let group_summaries =
        ordered_midpoint_group_summaries(group_idx_to_name, eligible_interval_counts)?;

    for group_summary in group_summaries {
        let group_idx = group_summary.group_idx;
        let eligible_intervals = group_summary.eligible_intervals;
        write!(writer, "{group_idx}\t").with_context(|| {
            format!(
                "writing midpoint group index row to {}",
                output_path.display()
            )
        })?;
        write_validated_tsv_field(
            &mut writer,
            output_path,
            group_summary.group_name,
            "group_name",
        )?;
        writeln!(writer, "\t{eligible_intervals}").with_context(|| {
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

/// Write a text field that must remain valid as one TSV cell.
///
/// Control characters are rejected rather than rewritten. Rewriting tabs or newlines would merge
/// distinct group names into the same public value.
fn write_validated_tsv_field(
    writer: &mut impl Write,
    output_path: &Path,
    value: &str,
    field_name: &str,
) -> Result<()> {
    ensure!(
        !value.chars().any(char::is_control),
        "{field_name} contains a control character and cannot be written to TSV without changing its value"
    );
    write!(writer, "{value}").with_context(|| format!("write {}", output_path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("group_index_tests.rs");
}
