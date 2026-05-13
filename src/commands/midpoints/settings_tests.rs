use crate::{
    commands::{
        cli_common::{ChromosomeArgs, IOCArgs},
        midpoints::{
            config::MidpointsConfig, postprocess::ProfileLayout, smoothing::MidpointSmoothing,
        },
    },
    shared::length_axis::LengthAxis,
};
use serde_json::Value;
use std::path::PathBuf;
use tempfile::TempDir;

use super::{last_position_bin_width, smoothing_settings};

fn settings_test_config() -> MidpointsConfig {
    MidpointsConfig::new(
        IOCArgs {
            bam: PathBuf::from("input.bam"),
            output_dir: PathBuf::from("out"),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        PathBuf::from("sites.bed"),
    )
}

#[test]
fn smoothing_settings_describes_raw_and_savgol_layouts() {
    let raw_layout =
        ProfileLayout::resolve(11, 1, MidpointSmoothing::Raw).expect("raw layout should resolve");
    let savgol_layout = ProfileLayout::resolve(7, 3, MidpointSmoothing::SavGol { window_bp: 7 })
        .expect("savgol layout should resolve");

    let raw = smoothing_settings(raw_layout);
    let savgol = smoothing_settings(savgol_layout);

    assert_eq!(raw.method, "raw");
    assert_eq!(raw.polynomial_order, None);
    assert_eq!(raw.window_bp, None);
    assert_eq!(raw.computation_flank_bp, 0);
    assert!(!raw.applied_before_binning);

    assert_eq!(savgol.method, "savitzky_golay");
    assert_eq!(savgol.polynomial_order, Some(3));
    assert_eq!(savgol.window_bp, Some(7));
    assert_eq!(savgol.computation_flank_bp, 3);
    assert!(savgol.applied_before_binning);
}

#[test]
fn last_position_bin_width_uses_real_width_of_final_bin() {
    // Output length 7 with bin size 3 produces [0,3), [3,6), and [6,7), so the last bin is 1 bp.
    let partial_final_bin =
        ProfileLayout::resolve(7, 3, MidpointSmoothing::Raw).expect("layout should resolve");
    assert_eq!(last_position_bin_width(partial_final_bin), 1);

    // Output length 9 with bin size 3 has no partial final bin, so the last bin keeps the full width.
    let full_final_bin =
        ProfileLayout::resolve(9, 3, MidpointSmoothing::Raw).expect("layout should resolve");
    assert_eq!(last_position_bin_width(full_final_bin), 3);
}

#[test]
fn midpoint_settings_json_records_savgol_and_final_position_binning() {
    // Arrange:
    // A 7 bp output interval with a 7 bp SavGol window has a 3 bp smoothing flank on each side,
    // so counting spans 13 bp. Final 3 bp binning yields three position bins, where the last bin
    // averages one retained base because 7 = 3 + 3 + 1.
    let temp = TempDir::new().expect("temp dir should be created");
    let settings_path = temp.path().join("sites.midpoint_profile_settings.json");
    let mut config = settings_test_config();
    config.set_smoothing(MidpointSmoothing::SavGol { window_bp: 7 });
    config.set_bin_size(3);
    config.blacklist = Some(vec![PathBuf::from("blacklist.bed")]);
    let length_axis =
        LengthAxis::new(vec![30, 40, 50]).expect("valid stepped length axis should resolve");
    let profile_layout = ProfileLayout::resolve(7, 3, config.smooth)
        .expect("profile layout should resolve for 7 bp SavGol over a 7 bp interval");

    // Act
    super::write_midpoint_profile_settings_json(
        &settings_path,
        &config,
        &length_axis,
        profile_layout,
        503,
        true,
    )
    .expect("settings json should write");

    // Assert
    let settings_text =
        std::fs::read_to_string(settings_path).expect("settings json should be readable");
    let settings: Value =
        serde_json::from_str(&settings_text).expect("settings json should parse");
    assert_eq!(settings["array_axes"], serde_json::json!(["group", "length_bin", "position"]));
    assert_eq!(settings["length_axis"]["column_intervals"], "half_open");
    assert_eq!(settings["length_axis"]["min_fragment_length"], 30);
    assert_eq!(settings["length_axis"]["max_fragment_length"], 49);
    assert_eq!(settings["length_axis"]["n_bins"], 2);
    assert_eq!(settings["length_axis"]["single_bp_bins"], false);
    assert_eq!(settings["length_axis"]["bin_definition"]["kind"], "stepped_range");
    assert_eq!(settings["length_axis"]["bin_definition"]["start"], 30);
    assert_eq!(settings["length_axis"]["bin_definition"]["end"], 50);
    assert_eq!(settings["length_axis"]["bin_definition"]["step"], 10);
    assert_eq!(
        settings["position_axis"]["coordinate_frame"],
        "interval_relative_zero_based"
    );
    assert_eq!(settings["position_axis"]["column_intervals"], "half_open");
    assert_eq!(settings["position_axis"]["output_interval_length_bp"], 7);
    assert_eq!(settings["position_axis"]["counted_interval_length_bp"], 13);
    assert_eq!(settings["position_axis"]["n_bins"], 3);
    assert_eq!(settings["position_axis"]["bin_size_bp"], 3);
    assert_eq!(settings["position_axis"]["bin_aggregation"], "mean");
    assert_eq!(settings["position_axis"]["last_bin_width_bp"], 1);
    assert_eq!(settings["smoothing"]["method"], "savitzky_golay");
    assert_eq!(settings["smoothing"]["polynomial_order"], 3);
    assert_eq!(settings["smoothing"]["window_bp"], 7);
    assert_eq!(settings["smoothing"]["computation_flank_bp"], 3);
    assert_eq!(settings["smoothing"]["applied_before_binning"], true);
    assert_eq!(settings["fragment_blacklist_used"], true);
    assert_eq!(settings["interval_blacklist_prefilter"]["enabled"], true);
    assert_eq!(settings["interval_blacklist_prefilter"]["margin_bp"], 503);
}

#[test]
fn midpoint_settings_json_records_raw_explicit_length_edges_without_savgol_fields() {
    // Arrange:
    // Raw full-resolution output has no smoothing window and the counted interval is identical to
    // the public interval. Explicit non-uniform length edges must be written back as explicit
    // edges so downstream readers do not infer a fake step size.
    let temp = TempDir::new().expect("temp dir should be created");
    let settings_path = temp.path().join("raw.midpoint_profile_settings.json");
    let config = settings_test_config();
    let length_axis =
        LengthAxis::new(vec![30, 80, 150]).expect("valid explicit length axis should resolve");
    let profile_layout = ProfileLayout::resolve(11, 1, config.smooth)
        .expect("raw profile layout should resolve");

    // Act
    super::write_midpoint_profile_settings_json(
        &settings_path,
        &config,
        &length_axis,
        profile_layout,
        0,
        false,
    )
    .expect("settings json should write");

    // Assert
    let settings_text =
        std::fs::read_to_string(settings_path).expect("settings json should be readable");
    let settings: Value =
        serde_json::from_str(&settings_text).expect("settings json should parse");
    assert_eq!(settings["length_axis"]["bin_definition"]["kind"], "explicit_edges");
    assert_eq!(
        settings["length_axis"]["bin_definition"]["edges"],
        serde_json::json!([30, 80, 150])
    );
    assert_eq!(settings["position_axis"]["output_interval_length_bp"], 11);
    assert_eq!(settings["position_axis"]["counted_interval_length_bp"], 11);
    assert_eq!(settings["position_axis"]["n_bins"], 11);
    assert_eq!(settings["position_axis"]["bin_size_bp"], 1);
    assert_eq!(settings["position_axis"]["last_bin_width_bp"], 1);
    assert_eq!(settings["smoothing"]["method"], "raw");
    assert!(settings["smoothing"].get("polynomial_order").is_none());
    assert!(settings["smoothing"].get("window_bp").is_none());
    assert_eq!(settings["smoothing"]["computation_flank_bp"], 0);
    assert_eq!(settings["smoothing"]["applied_before_binning"], false);
    assert_eq!(settings["fragment_blacklist_used"], false);
    assert_eq!(settings["interval_blacklist_prefilter"]["enabled"], false);
    assert!(
        settings["interval_blacklist_prefilter"]
            .get("margin_bp")
            .is_none()
    );
}
