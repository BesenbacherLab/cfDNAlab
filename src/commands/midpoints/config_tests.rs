use super::MidpointsConfig;
use crate::{
    commands::cli_common::{ChromosomeArgs, IOCArgs},
    shared::constants::MAX_SUPPORTED_FRAGMENT_LENGTH,
};
use std::path::PathBuf;

fn config_for_length_bin_resolution() -> MidpointsConfig {
    MidpointsConfig::new(
        IOCArgs {
            bam: PathBuf::from("input.bam"),
            output_dir: PathBuf::from("out"),
            n_threads: 1,
        },
        ChromosomeArgs::default(),
        PathBuf::from("intervals.tsv"),
    )
}

#[test]
fn resolve_length_bins_keeps_midpoints_default_as_one_broad_bin() {
    let config = config_for_length_bin_resolution();

    let edges = config
        .resolve_length_bins()
        .expect("default midpoint length bins should resolve");

    assert_eq!(edges, vec![30, 1001]);
}

#[test]
fn resolve_length_bins_accepts_explicit_per_bp_range_spec() {
    let mut config = config_for_length_bin_resolution();
    config.set_length_bins_spec("30:33:1");

    let edges = config
        .resolve_length_bins()
        .expect("explicit dense midpoint length bins should resolve");

    assert_eq!(edges, vec![30, 31, 32, 33]);
}

#[test]
fn resolve_length_bins_accepts_edge_for_max_supported_fragment_length() {
    let mut config = config_for_length_bin_resolution();
    config.set_length_bins(vec![10, MAX_SUPPORTED_FRAGMENT_LENGTH + 1]);

    let edges = config
        .resolve_length_bins()
        .expect("exclusive edge at max supported length + 1 should be valid");

    assert_eq!(edges, vec![10, MAX_SUPPORTED_FRAGMENT_LENGTH + 1]);
}

#[test]
fn resolve_length_bins_rejects_edge_past_max_supported_fragment_length() {
    let invalid_edge = MAX_SUPPORTED_FRAGMENT_LENGTH + 2;
    let mut config = config_for_length_bin_resolution();
    config.set_length_bins(vec![10, invalid_edge]);

    let error = config
        .resolve_length_bins()
        .expect_err("edges beyond the supported fragment length cap should fail");
    let message = error.to_string();

    assert!(
        message.contains(&format!(
            "length bin edges must be <= {}",
            MAX_SUPPORTED_FRAGMENT_LENGTH + 1
        )),
        "unexpected error: {message}"
    );
}

#[test]
fn resolve_length_bins_rejects_range_past_max_supported_fragment_length() {
    let invalid_end = MAX_SUPPORTED_FRAGMENT_LENGTH + 2;
    let mut config = config_for_length_bin_resolution();
    config.set_length_bins_spec(format!("10:{invalid_end}:1"));

    let error = config
        .resolve_length_bins()
        .expect_err("range specs beyond the supported fragment length cap should fail");
    let message = error.to_string();

    assert!(
        message.contains(&format!(
            "length-bins end ({invalid_end}) must be <= max fragment length + 1 ({})",
            MAX_SUPPORTED_FRAGMENT_LENGTH + 1
        )),
        "unexpected error: {message}"
    );
}
