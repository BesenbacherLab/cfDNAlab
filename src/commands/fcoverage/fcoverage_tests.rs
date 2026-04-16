use super::{
    internal_residual_coverage_floor, minimum_positive_base_weight,
    minimum_positive_gc_weight, minimum_positive_pre_scaling_support,
};
use crate::commands::cli_common::{ApplyGCArgs, ChromosomeArgs, IOCArgs};
use crate::commands::fcoverage::config::FCoverageConfig;
use crate::shared::gc_tag::MIN_REASONABLE_GC_WEIGHT;
use std::path::PathBuf;

fn base_config() -> FCoverageConfig {
    FCoverageConfig::new(
        IOCArgs {
            bam: PathBuf::from("input.bam"),
            output_dir: PathBuf::from("out"),
            n_threads: 1,
        },
        ChromosomeArgs::default(),
    )
}

#[test]
fn minimum_positive_support_is_one_without_gc_or_length_normalization() {
    let opt = base_config();

    assert_eq!(minimum_positive_base_weight(&opt), 1.0);
    assert_eq!(minimum_positive_gc_weight(&opt), 1.0);
    assert_eq!(minimum_positive_pre_scaling_support(&opt), 1.0);
    assert_eq!(internal_residual_coverage_floor(&opt), 0.5);
}

#[test]
fn minimum_positive_support_uses_max_fragment_length_when_length_normalized() {
    let mut opt = base_config();
    opt.set_normalize_by_length(true);
    opt.fragment_lengths_mut().max_fragment_length = 500;

    let expected_min_support = 1.0 / 500.0;

    assert_eq!(minimum_positive_base_weight(&opt), expected_min_support);
    assert_eq!(minimum_positive_gc_weight(&opt), 1.0);
    assert_eq!(
        minimum_positive_pre_scaling_support(&opt),
        expected_min_support
    );
    assert_eq!(
        internal_residual_coverage_floor(&opt),
        expected_min_support / 2.0
    );
}

#[test]
fn minimum_positive_support_uses_gc_lower_bound_for_gc_file_runs() {
    let mut opt = base_config();
    opt.set_gc(ApplyGCArgs {
        gc_file: Some(PathBuf::from("gc_bias_correction.npz")),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });

    assert_eq!(minimum_positive_base_weight(&opt), 1.0);
    assert_eq!(minimum_positive_gc_weight(&opt), MIN_REASONABLE_GC_WEIGHT);
    assert_eq!(
        minimum_positive_pre_scaling_support(&opt),
        MIN_REASONABLE_GC_WEIGHT
    );
    assert_eq!(
        internal_residual_coverage_floor(&opt),
        MIN_REASONABLE_GC_WEIGHT / 2.0
    );
}

#[test]
fn minimum_positive_support_uses_gc_lower_bound_for_gc_tag_runs() {
    let mut opt = base_config();
    opt.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    });

    assert_eq!(minimum_positive_base_weight(&opt), 1.0);
    assert_eq!(minimum_positive_gc_weight(&opt), MIN_REASONABLE_GC_WEIGHT);
    assert_eq!(
        minimum_positive_pre_scaling_support(&opt),
        MIN_REASONABLE_GC_WEIGHT
    );
}

#[test]
fn internal_cleanup_floor_stays_below_theoretical_minimum_with_gc_and_length_normalization() {
    let mut opt = base_config();
    opt.set_normalize_by_length(true);
    opt.fragment_lengths_mut().max_fragment_length = 1000;
    opt.set_gc(ApplyGCArgs {
        gc_file: Some(PathBuf::from("gc_bias_correction.npz")),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });

    // With --normalize-by-length, the smallest real positive per-base mass comes from the
    // longest allowed fragment. GC correction can lower that further down to the minimum
    // supported positive GC weight.
    let min_support = (1.0 / 1000.0) * MIN_REASONABLE_GC_WEIGHT;
    let cleanup_floor = internal_residual_coverage_floor(&opt);

    assert_eq!(minimum_positive_pre_scaling_support(&opt), min_support);
    assert_eq!(cleanup_floor, min_support / 2.0);
    assert!(cleanup_floor > 0.0);
    assert!(cleanup_floor < min_support);
}
