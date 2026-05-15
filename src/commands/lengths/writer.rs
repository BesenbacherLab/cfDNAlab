use crate::{
    commands::{
        cli_common::{DistributionWindowSpec, WindowAssigner},
        gc_bias::correct::{GCLengthRange, MarginalizeLengthsWeightingScheme},
        lengths::{config::LengthsConfig, counting::LengthCounts},
    },
    shared::{
        bed::GroupedWindows,
        blacklist::compute_blacklist_overlap,
        clip_mode::ClipMode,
        formatters::CompactNumber,
        indel_mode::IndelMode,
        interval::Interval,
        io::create_text_writer,
        length_axis::{LengthAxis, LengthAxisSettings},
        windowing::WindowBinInfo,
    },
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use serde::Serialize;
use std::{io::Write, path::Path};

const BLACKLISTED_FRACTION_DECIMALS: i32 = 3;

/// Interpretation metadata for a fragment-length count table.
///
/// This sidecar records the information needed to understand the count columns
/// and row semantics without inspecting the command line that produced the
/// output. Ordinary filters stay out of the sidecar unless they change how the
/// table should be interpreted downstream.
#[derive(Serialize)]
struct LengthSettings<'a> {
    length_axis: LengthAxisSettings<'a>,
    aggregation_level: &'static str,
    window_mode: &'static str,
    indel_mode: &'static str,
    clip_mode: &'static str,
    max_soft_clips: u16,
    max_deletion_bases: u16,
    assign_by: String,
    decimals: u8,
    gc_length_weighting: &'static str,
    gc_length_range: &'static str,
    gc_length_trim_rare: f64,
    gc_correction_used: bool,
    scaling_factors_used: bool,
    blacklist_used: bool,
}

/// Write the JSON sidecar that describes a fragment-length count table.
///
/// The settings are intentionally focused on interpretation. They include the
/// length axis, row aggregation mode, assignment mode, and whether optional
/// weighting inputs affected the counts.
pub(super) fn write_length_settings_json(
    settings_path: &std::path::Path,
    opt: &LengthsConfig,
    window_opt: &DistributionWindowSpec,
    length_axis: &LengthAxis,
) -> Result<()> {
    let settings = LengthSettings {
        length_axis: length_axis.settings(),
        aggregation_level: aggregation_level_name(window_opt),
        window_mode: window_mode_name(window_opt),
        indel_mode: indel_mode_name(opt.indel_mode),
        clip_mode: clip_mode_name(opt.clip_mode),
        max_soft_clips: opt.max_soft_clips,
        max_deletion_bases: opt.max_deletion_bases,
        assign_by: window_assigner_name(opt.window_assignment.assign_by),
        decimals: opt.decimals,
        gc_length_weighting: gc_length_weighting_name(opt.gc_length_weighting),
        gc_length_range: gc_length_range_name(opt.gc_length_range),
        gc_length_trim_rare: opt.gc_length_trim_rare,
        gc_correction_used: opt.gc.gc_file.is_some(),
        scaling_factors_used: opt.scale_genome.scaling_factors.is_some(),
        blacklist_used: opt.blacklist.is_some(),
    };

    let mut settings_writer = create_text_writer(settings_path)
        .with_context(|| format!("create {}", settings_path.display()))?;
    serde_json::to_writer_pretty(&mut settings_writer, &settings)
        .with_context(|| format!("write {}", settings_path.display()))?;
    writeln!(settings_writer).with_context(|| format!("write {}", settings_path.display()))?;
    settings_writer
        .finish()
        .with_context(|| format!("finalize {}", settings_path.display()))?;
    Ok(())
}

/// Write the public fragment-length count table.
///
/// The table is intentionally wide: each output unit is one row, and each fragment-length bin is
/// one count column. Single-bp bins use `count_<length>`, and wider bins use
/// `count_<start>_<end>`. This keeps row metadata from being repeated for every length bin while
/// remaining directly readable as a TSV in R and Python. `decimals` controls only the written count
/// representation, not the internal aggregation precision. Blacklist fractions are always rounded
/// to three decimals.
pub(super) fn write_length_counts_tsv(
    output_path: &Path,
    counts: &[LengthCounts],
    length_axis: &LengthAxis,
    decimals: u8,
    row_metadata: LengthCountRowMetadata<'_>,
) -> Result<()> {
    let decimals = i32::from(decimals);
    let count_headers = length_count_column_headers(length_axis);
    ensure!(
        counts.iter().all(
            |length_counts| length_counts.axis.edges() == length_axis.edges()
                && length_counts.counts.len() == count_headers.len()
        ),
        "length count rows do not all match the output length axis and column count"
    );

    let mut writer = create_text_writer(output_path)
        .with_context(|| format!("create {}", output_path.display()))?;

    match row_metadata {
        LengthCountRowMetadata::Global => {
            ensure!(
                counts.len() == 1,
                "global length count output should have one row, found {}",
                counts.len()
            );
            write_count_header(&mut writer, output_path, &[], &count_headers)?;
            write_count_values(&mut writer, output_path, &counts[0].counts, decimals, false)?;
        }
        LengthCountRowMetadata::Windows {
            windows,
            include_blacklisted_fraction,
        } => {
            ensure!(
                counts.len() == windows.len(),
                "window metadata entries ({}) did not match length count rows ({})",
                windows.len(),
                counts.len()
            );
            let metadata_headers = if include_blacklisted_fraction {
                &["chrom", "start", "end", "blacklisted_fraction"][..]
            } else {
                &["chrom", "start", "end"][..]
            };
            write_count_header(&mut writer, output_path, metadata_headers, &count_headers)?;
            for (window, length_counts) in windows.iter().zip(counts.iter()) {
                write!(
                    writer,
                    "{}\t{}\t{}",
                    window.chromosome, window.start, window.end
                )
                .with_context(|| format!("write {}", output_path.display()))?;
                if include_blacklisted_fraction {
                    write!(
                        writer,
                        "\t{}",
                        CompactNumber {
                            v: window.blacklisted_fraction,
                            decimals: BLACKLISTED_FRACTION_DECIMALS,
                        }
                    )
                    .with_context(|| format!("write {}", output_path.display()))?;
                }
                write_count_values(
                    &mut writer,
                    output_path,
                    &length_counts.counts,
                    decimals,
                    true,
                )?;
            }
        }
        LengthCountRowMetadata::Groups {
            group_idx_to_name,
            chromosomes,
            grouped_windows_map,
            blacklist_map,
            include_blacklisted_fraction,
        } => {
            let group_summaries = grouped_length_count_row_summaries(
                group_idx_to_name,
                chromosomes,
                grouped_windows_map,
                blacklist_map,
                include_blacklisted_fraction,
            )?;
            ensure!(
                counts.len() == group_summaries.len(),
                "group metadata entries ({}) did not match length count rows ({})",
                group_summaries.len(),
                counts.len()
            );
            ensure_group_summaries_match_count_rows(&group_summaries)?;
            let metadata_headers = if include_blacklisted_fraction {
                &["group_name", "eligible_windows", "blacklisted_fraction"][..]
            } else {
                &["group_name", "eligible_windows"][..]
            };
            write_count_header(&mut writer, output_path, metadata_headers, &count_headers)?;
            for (summary, length_counts) in group_summaries.iter().zip(counts.iter()) {
                write_validated_tsv_field(
                    &mut writer,
                    output_path,
                    summary.group_name,
                    "group_name",
                )?;
                write!(writer, "\t{}", summary.eligible_windows)
                    .with_context(|| format!("write {}", output_path.display()))?;
                if include_blacklisted_fraction {
                    write!(
                        writer,
                        "\t{}",
                        CompactNumber {
                            v: summary.blacklisted_fraction,
                            decimals: BLACKLISTED_FRACTION_DECIMALS,
                        }
                    )
                    .with_context(|| format!("write {}", output_path.display()))?;
                }
                write_count_values(
                    &mut writer,
                    output_path,
                    &length_counts.counts,
                    decimals,
                    true,
                )?;
            }
        }
    }

    writer
        .finish()
        .with_context(|| format!("finalize {}", output_path.display()))?;
    Ok(())
}

/// Row metadata needed to write the public length-count table.
///
/// The count matrix has only numeric values. This enum supplies the row keys
/// that make those values usable without sidecar joins.
pub(super) enum LengthCountRowMetadata<'a> {
    /// One global row with only count columns.
    Global,
    /// One row per fixed-size or BED window.
    Windows {
        /// Window coordinates in the same order as the count rows.
        windows: &'a [WindowBinInfo],
        /// Whether the table should expose the per-window blacklist fraction.
        include_blacklisted_fraction: bool,
    },
    /// One row per grouped-BED group.
    Groups {
        /// Mapping from internal group index to public group name.
        group_idx_to_name: &'a FxHashMap<u64, String>,
        /// Chromosomes retained by the run, in processing order.
        chromosomes: &'a [String],
        /// Grouped windows retained after chromosome filtering.
        grouped_windows_map: &'a FxHashMap<String, GroupedWindows>,
        /// Blacklist intervals keyed by chromosome.
        blacklist_map: &'a FxHashMap<String, Vec<Interval<u64>>>,
        /// Whether the table should expose the per-group blacklist fraction.
        include_blacklisted_fraction: bool,
    },
}

/// Metadata for one grouped output row.
///
/// The `group_idx` is retained even though it is not written. It lets the
/// writer verify that the summary row order still matches the count row order.
struct GroupedLengthCountRowSummary<'a> {
    group_idx: u64,
    group_name: &'a str,
    eligible_windows: usize,
    blacklisted_fraction: f64,
}

/// Build the public `count_*` column names for a length axis.
///
/// Single-bp bins are named `count_<length>`. Wider half-open bins are named
/// `count_<start>_<end>`, where `end` is exclusive.
fn length_count_column_headers(length_axis: &LengthAxis) -> Vec<String> {
    length_axis
        .edges()
        .windows(2)
        .map(|edges| {
            let start = edges[0];
            let end = edges[1];
            if end == start + 1 {
                format!("count_{start}")
            } else {
                format!("count_{start}_{end}")
            }
        })
        .collect()
}

/// Write the header row for a length-count table.
///
/// Metadata columns come first, followed by count columns. Global outputs pass
/// no metadata headers, which makes the first count header the first column.
fn write_count_header(
    writer: &mut impl Write,
    output_path: &Path,
    metadata_headers: &[&str],
    count_headers: &[String],
) -> Result<()> {
    let mut is_first_column = true;
    for header in metadata_headers {
        if !is_first_column {
            write!(writer, "\t").with_context(|| format!("write {}", output_path.display()))?;
        }
        write!(writer, "{header}").with_context(|| format!("write {}", output_path.display()))?;
        is_first_column = false;
    }
    for header in count_headers {
        if !is_first_column {
            write!(writer, "\t").with_context(|| format!("write {}", output_path.display()))?;
        }
        write!(writer, "{header}").with_context(|| format!("write {}", output_path.display()))?;
        is_first_column = false;
    }
    writeln!(writer).with_context(|| format!("write {}", output_path.display()))?;
    Ok(())
}

/// Write the numeric count portion of one output row.
///
/// `has_leading_metadata` tells the helper whether a tab is needed before the
/// first count. This keeps global rows and metadata-keyed rows on the same
/// code path.
fn write_count_values(
    writer: &mut impl Write,
    output_path: &Path,
    counts: &[f64],
    decimals: i32,
    has_leading_metadata: bool,
) -> Result<()> {
    let mut is_first_count = !has_leading_metadata;
    for count in counts {
        if !is_first_count {
            write!(writer, "\t").with_context(|| format!("write {}", output_path.display()))?;
        }
        write!(
            writer,
            "{}",
            CompactNumber {
                v: *count,
                decimals
            }
        )
        .with_context(|| format!("write {}", output_path.display()))?;
        is_first_count = false;
    }
    writeln!(writer).with_context(|| format!("write {}", output_path.display()))?;
    Ok(())
}

/// Build grouped row summaries in count-row order.
///
/// Grouped counts are indexed by the internal `group_idx`. The returned vector
/// is sorted by that index, and the caller verifies that the sorted index is
/// exactly the row number before writing counts next to metadata.
fn grouped_length_count_row_summaries<'a>(
    group_idx_to_name: &'a FxHashMap<u64, String>,
    chromosomes: &[String],
    grouped_windows_map: &FxHashMap<String, GroupedWindows>,
    blacklist_map: &FxHashMap<String, Vec<Interval<u64>>>,
    include_blacklisted_fraction: bool,
) -> Result<Vec<GroupedLengthCountRowSummary<'a>>> {
    let mut eligible_windows_by_group: FxHashMap<u64, usize> = group_idx_to_name
        .keys()
        .map(|&group_idx| (group_idx, 0usize))
        .collect();
    let mut total_bp_by_group: FxHashMap<u64, u64> = FxHashMap::default();
    let mut blacklisted_bp_by_group: FxHashMap<u64, f64> = FxHashMap::default();

    for chromosome in chromosomes {
        let windows = grouped_windows_map
            .get(chromosome)
            .map(|windows| windows.windows_as_slice())
            .unwrap_or(&[]);
        let blacklist_intervals = blacklist_map
            .get(chromosome)
            .map(|intervals| intervals.as_slice())
            .unwrap_or(&[]);
        let mut blacklist_ptr = 0usize;
        for window in windows {
            let (start, end, group_idx) = window.as_tuple();
            ensure!(
                group_idx_to_name.contains_key(&group_idx),
                "grouped window references group_idx {group_idx}, but no group name was registered"
            );
            *eligible_windows_by_group
                .get_mut(&group_idx)
                .expect("validated group_idx should be pre-seeded") += 1;
            if include_blacklisted_fraction {
                let window_bp = end
                    .checked_sub(start)
                    .context("grouped window end must be >= start")?;
                let blacklisted_fraction = compute_blacklist_overlap(
                    blacklist_intervals,
                    Interval::new(start, end)?,
                    0,
                    &mut blacklist_ptr,
                );
                *total_bp_by_group.entry(group_idx).or_insert(0) += window_bp;
                *blacklisted_bp_by_group.entry(group_idx).or_insert(0.0) +=
                    blacklisted_fraction * window_bp as f64;
            }
        }
    }

    let mut entries: Vec<(u64, &str)> = group_idx_to_name
        .iter()
        .map(|(group_idx, group_name)| (*group_idx, group_name.as_str()))
        .collect();
    entries.sort_unstable_by_key(|(group_idx, _)| *group_idx);

    Ok(entries
        .into_iter()
        .map(|(group_idx, group_name)| {
            let eligible_windows = eligible_windows_by_group
                .get(&group_idx)
                .copied()
                .unwrap_or(0);
            let total_bp = total_bp_by_group.get(&group_idx).copied().unwrap_or(0);
            let blacklisted_fraction = if total_bp == 0 {
                0.0
            } else {
                blacklisted_bp_by_group
                    .get(&group_idx)
                    .copied()
                    .unwrap_or(0.0)
                    / total_bp as f64
            };
            GroupedLengthCountRowSummary {
                group_idx,
                group_name,
                eligible_windows,
                blacklisted_fraction,
            }
        })
        .collect())
}

/// Verify that grouped metadata can be zipped with grouped count rows.
///
/// The sorted set of internal group indices must be exactly `0..n_rows`.
/// Otherwise zipping grouped summaries with count rows would attach counts to
/// the wrong public group name.
fn ensure_group_summaries_match_count_rows(
    group_summaries: &[GroupedLengthCountRowSummary<'_>],
) -> Result<()> {
    for (row_index, summary) in group_summaries.iter().enumerate() {
        ensure!(
            summary.group_idx == row_index as u64,
            "grouped length count row {row_index} corresponds to group_idx {}, expected group_idx {row_index}",
            summary.group_idx
        );
    }
    Ok(())
}

/// Write a text field that must remain valid as one TSV cell.
///
/// Control characters are rejected rather than rewritten. Rewriting tabs or
/// newlines would merge distinct group names into the same public value.
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

/// Return the row aggregation label stored in the settings sidecar.
fn aggregation_level_name(window_opt: &DistributionWindowSpec) -> &'static str {
    match window_opt {
        DistributionWindowSpec::Global => "global",
        DistributionWindowSpec::GroupedBed(_) => "groups",
        DistributionWindowSpec::Size(_) | DistributionWindowSpec::Bed(_) => "windows",
    }
}

/// Return the windowing mode label stored in the settings sidecar.
fn window_mode_name(window_opt: &DistributionWindowSpec) -> &'static str {
    match window_opt {
        DistributionWindowSpec::Global => "global",
        DistributionWindowSpec::Size(_) => "by-size",
        DistributionWindowSpec::Bed(_) => "by-bed",
        DistributionWindowSpec::GroupedBed(_) => "by-grouped-bed",
    }
}

/// Return the indel handling label stored in the settings sidecar.
fn indel_mode_name(indel_mode: IndelMode) -> &'static str {
    match indel_mode {
        IndelMode::Ignore => "ignore",
        IndelMode::Adjust => "adjust",
        IndelMode::Skip => "skip",
    }
}

/// Return the soft-clip handling label stored in the settings sidecar.
fn clip_mode_name(clip_mode: ClipMode) -> &'static str {
    match clip_mode {
        ClipMode::Aligned => "aligned",
        ClipMode::Adjust => "adjust",
        ClipMode::Skip => "skip",
    }
}

/// Return the window assignment label stored in the settings sidecar.
fn window_assigner_name(assigner: WindowAssigner) -> String {
    match assigner {
        WindowAssigner::CountOverlap => "count-overlap".to_string(),
        WindowAssigner::Any => "any".to_string(),
        WindowAssigner::All => "all".to_string(),
        WindowAssigner::Midpoint => "midpoint".to_string(),
        WindowAssigner::Proportion(threshold) => format!("proportion={threshold}"),
    }
}

/// Return the GC length-weighting label stored in the settings sidecar.
fn gc_length_weighting_name(weighting: MarginalizeLengthsWeightingScheme) -> &'static str {
    match weighting {
        MarginalizeLengthsWeightingScheme::Equal => "equal",
        MarginalizeLengthsWeightingScheme::Frequency => "frequency",
        MarginalizeLengthsWeightingScheme::MaxFrequency => "max-frequency",
    }
}

/// Return the GC package range label stored in the settings sidecar.
fn gc_length_range_name(range: GCLengthRange) -> &'static str {
    match range {
        GCLengthRange::Requested => "requested",
        GCLengthRange::Package => "package",
    }
}

#[cfg(test)]
mod tests {
    include!("writer_tests.rs");
}
