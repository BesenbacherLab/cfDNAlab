use crate::{
    commands::{
        cli_common::{DistributionWindowSpec, WindowAssigner},
        gc_bias::correct::{GCLengthRange, MarginalizeLengthsWeightingScheme},
        lengths::{config::LengthsConfig, counting::LengthAxis},
    },
    shared::{clip_mode::ClipMode, indel_mode::IndelMode, io::create_text_writer},
};
use anyhow::{Context, Result};
use serde::Serialize;
use std::io::Write;

#[derive(Serialize)]
struct LengthAxisSettings<'a> {
    column_intervals: &'static str,
    min_fragment_length: u32,
    max_fragment_length: u32,
    n_bins: usize,
    single_bp_bins: bool,
    bin_definition: LengthBinDefinition<'a>,
}

#[derive(Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum LengthBinDefinition<'a> {
    SteppedRange { start: u32, end: u32, step: u32 },
    ExplicitEdges { edges: &'a [u32] },
}

#[derive(Serialize)]
struct FragmentLengthSettings<'a> {
    length_axis: LengthAxisSettings<'a>,
    aggregation_level: &'static str,
    window_mode: &'static str,
    indel_mode: &'static str,
    clip_mode: &'static str,
    max_soft_clips: u16,
    max_deletion_bases: u16,
    assign_by: String,
    gc_length_weighting: &'static str,
    gc_length_range: &'static str,
    gc_length_trim_rare: f64,
    gc_correction_used: bool,
    scaling_factors_used: bool,
}

pub(super) fn write_fragment_length_settings_json(
    settings_path: &std::path::Path,
    opt: &LengthsConfig,
    window_opt: &DistributionWindowSpec,
    length_axis: &LengthAxis,
) -> Result<()> {
    let settings = FragmentLengthSettings {
        length_axis: LengthAxisSettings {
            column_intervals: "half_open",
            min_fragment_length: length_axis.min_fragment_length(),
            max_fragment_length: length_axis.max_fragment_length(),
            n_bins: length_axis.num_bins(),
            single_bp_bins: length_axis.is_single_bp_bins(),
            bin_definition: length_bin_definition(length_axis.edges()),
        },
        aggregation_level: aggregation_level_name(window_opt),
        window_mode: window_mode_name(window_opt),
        indel_mode: indel_mode_name(opt.indel_mode),
        clip_mode: clip_mode_name(opt.clip_mode),
        max_soft_clips: opt.max_soft_clips,
        max_deletion_bases: opt.max_deletion_bases,
        assign_by: window_assigner_name(opt.window_assignment.assign_by),
        gc_length_weighting: gc_length_weighting_name(opt.gc_length_weighting),
        gc_length_range: gc_length_range_name(opt.gc_length_range),
        gc_length_trim_rare: opt.gc_length_trim_rare,
        gc_correction_used: opt.gc.gc_file.is_some(),
        scaling_factors_used: opt.scale_genome.scaling_factors.is_some(),
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

fn length_bin_definition(edges: &[u32]) -> LengthBinDefinition<'_> {
    if let Some(step) = stepped_range_step(edges) {
        LengthBinDefinition::SteppedRange {
            start: edges[0],
            end: *edges.last().expect("validated length axis has edges"),
            step,
        }
    } else {
        LengthBinDefinition::ExplicitEdges { edges }
    }
}

fn stepped_range_step(edges: &[u32]) -> Option<u32> {
    let step = edges.get(1)?.checked_sub(*edges.first()?)?;
    if step == 0 {
        return None;
    }

    let last_bin_index = edges.len().checked_sub(2)?;
    for (bin_index, edge_pair) in edges.windows(2).enumerate() {
        let width = edge_pair[1].checked_sub(edge_pair[0])?;
        if bin_index == last_bin_index {
            if width == 0 || width > step {
                return None;
            }
        } else if width != step {
            return None;
        }
    }

    Some(step)
}

fn aggregation_level_name(window_opt: &DistributionWindowSpec) -> &'static str {
    match window_opt {
        DistributionWindowSpec::Global => "global",
        DistributionWindowSpec::GroupedBed(_) => "groups",
        DistributionWindowSpec::Size(_) | DistributionWindowSpec::Bed(_) => "windows",
    }
}

fn window_mode_name(window_opt: &DistributionWindowSpec) -> &'static str {
    match window_opt {
        DistributionWindowSpec::Global => "global",
        DistributionWindowSpec::Size(_) => "by-size",
        DistributionWindowSpec::Bed(_) => "by-bed",
        DistributionWindowSpec::GroupedBed(_) => "by-grouped-bed",
    }
}

fn indel_mode_name(indel_mode: IndelMode) -> &'static str {
    match indel_mode {
        IndelMode::Ignore => "ignore",
        IndelMode::Adjust => "adjust",
        IndelMode::Skip => "skip",
    }
}

fn clip_mode_name(clip_mode: ClipMode) -> &'static str {
    match clip_mode {
        ClipMode::Aligned => "aligned",
        ClipMode::Adjust => "adjust",
        ClipMode::Skip => "skip",
    }
}

fn window_assigner_name(assigner: WindowAssigner) -> String {
    match assigner {
        WindowAssigner::CountOverlap => "count-overlap".to_string(),
        WindowAssigner::Any => "any".to_string(),
        WindowAssigner::All => "all".to_string(),
        WindowAssigner::Midpoint => "midpoint".to_string(),
        WindowAssigner::Proportion(threshold) => format!("proportion={threshold}"),
    }
}

fn gc_length_weighting_name(weighting: MarginalizeLengthsWeightingScheme) -> &'static str {
    match weighting {
        MarginalizeLengthsWeightingScheme::Equal => "equal",
        MarginalizeLengthsWeightingScheme::Frequency => "frequency",
        MarginalizeLengthsWeightingScheme::MaxFrequency => "max-frequency",
    }
}

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
