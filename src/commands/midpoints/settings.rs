use crate::{
    commands::midpoints::{
        config::MidpointsConfig, postprocess::ProfileLayout, smoothing::SAVGOL_POLYNOMIAL_ORDER,
    },
    shared::{
        io::create_text_writer,
        length_axis::{LengthAxis, LengthAxisSettings},
    },
};
use anyhow::{Context, Result};
use serde::Serialize;
use std::{io::Write, path::Path};

#[derive(Serialize)]
struct PositionAxisSettings {
    coordinate_frame: &'static str,
    column_intervals: &'static str,
    output_interval_length_bp: usize,
    counted_interval_length_bp: usize,
    n_bins: usize,
    bin_size_bp: u32,
    bin_aggregation: &'static str,
    last_bin_width_bp: u32,
}

#[derive(Serialize)]
struct SmoothingSettings {
    method: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    polynomial_order: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    window_bp: Option<u32>,
    computation_flank_bp: u32,
    applied_before_binning: bool,
}

#[derive(Serialize)]
struct IntervalBlacklistPrefilterSettings {
    enabled: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    margin_bp: Option<u64>,
}

#[derive(Serialize)]
struct MidpointSettings<'a> {
    array_axes: [&'static str; 3],
    length_axis: LengthAxisSettings<'a>,
    position_axis: PositionAxisSettings,
    smoothing: SmoothingSettings,
    gc_correction_used: bool,
    scaling_factors_used: bool,
    fragment_blacklist_used: bool,
    interval_blacklist_prefilter: IntervalBlacklistPrefilterSettings,
}

/// Write the midpoint command settings sidecar.
///
/// The Zarr store carries the axis metadata needed for downstream loading. This JSON file records
/// command settings in a plain human-readable sidecar.
///
/// Parameters
/// ----------
/// - `settings_path`:
///     Destination path in the command's final-output temp directory.
/// - `opt`:
///     Runtime configuration used for filters and correction flags.
/// - `length_axis`:
///     Resolved fragment length bins used for axis 1 of the output array.
/// - `profile_layout`:
///     Resolved counted and written position dimensions.
/// - `interval_blacklist_margin`:
///     Margin used when interval-level blacklist prefiltering is enabled.
/// - `use_blacklist_prefilter`:
///     Whether interval-level blacklist prefiltering was active for this run.
pub(crate) fn write_midpoint_settings_json(
    settings_path: &Path,
    opt: &MidpointsConfig,
    length_axis: &LengthAxis,
    profile_layout: ProfileLayout,
    interval_blacklist_margin: u64,
    use_blacklist_prefilter: bool,
) -> Result<()> {
    let settings = MidpointSettings {
        array_axes: ["group", "length_bin", "position"],
        length_axis: length_axis.settings(),
        position_axis: PositionAxisSettings {
            coordinate_frame: "interval_relative_zero_based",
            column_intervals: "half_open",
            output_interval_length_bp: profile_layout.output_len,
            counted_interval_length_bp: profile_layout.flanked_length,
            n_bins: profile_layout.output_positions,
            bin_size_bp: profile_layout.bin_size,
            bin_aggregation: "mean",
            last_bin_width_bp: last_position_bin_width(profile_layout),
        },
        smoothing: smoothing_settings(profile_layout),
        gc_correction_used: opt.gc.gc_file.is_some() || opt.gc.gc_tag.is_some(),
        scaling_factors_used: opt.scale_genome.scaling_factors.is_some(),
        fragment_blacklist_used: opt.blacklist.is_some(),
        interval_blacklist_prefilter: IntervalBlacklistPrefilterSettings {
            enabled: use_blacklist_prefilter,
            margin_bp: use_blacklist_prefilter.then_some(interval_blacklist_margin),
        },
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

fn smoothing_settings(profile_layout: ProfileLayout) -> SmoothingSettings {
    match profile_layout.smoothing_window {
        None => SmoothingSettings {
            method: "none",
            polynomial_order: None,
            window_bp: None,
            computation_flank_bp: profile_layout.smoothing_flank,
            applied_before_binning: false,
        },
        Some(window_bp) => SmoothingSettings {
            method: "savitzky_golay",
            polynomial_order: Some(SAVGOL_POLYNOMIAL_ORDER),
            window_bp: Some(window_bp),
            computation_flank_bp: profile_layout.smoothing_flank,
            applied_before_binning: true,
        },
    }
}

fn last_position_bin_width(profile_layout: ProfileLayout) -> u32 {
    let remainder = profile_layout.output_len % profile_layout.bin_size as usize;
    if remainder == 0 {
        profile_layout.bin_size
    } else {
        remainder as u32
    }
}

#[cfg(test)]
mod tests {
    include!("settings_tests.rs");
}
