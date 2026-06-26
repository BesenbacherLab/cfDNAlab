use super::{
    ApplyGCArgFileOnly, ApplyGCArgs, ChromosomeArgs, ContigSource, FragmentLengthArgs,
    WindowAssigner, min_overlap_fraction_for_window_assignment, parse_length_bins,
    parse_output_prefix, parse_sam_aux_tag_name, resolve_length_bin_edges,
    validate_max_soft_clips, validate_output_prefix,
};
use crate::shared::constants::{MAX_MAX_SOFT_CLIPS, MAX_SUPPORTED_FRAGMENT_LENGTH};
use std::io::Write;
use std::path::{Path, PathBuf};

fn assert_close(actual: f64, expected: f64) {
    let delta = (actual - expected).abs();
    assert!(
        delta <= 1e-12,
        "expected {expected}, got {actual}, delta {delta}"
    );
}

#[test]
fn fragment_length_args_rejects_inverted_range() {
    let args = FragmentLengthArgs {
        min_fragment_length: 500,
        max_fragment_length: 100,
    };

    let error = args
        .validate()
        .expect_err("inverted fragment length range should fail");
    let message = error.to_string();

    assert!(
        message.contains("--min-fragment-length (500) must be <= --max-fragment-length (100)"),
        "unexpected error: {message}"
    );
}

#[test]
fn fragment_length_args_accepts_inclusive_single_length_range() {
    let args = FragmentLengthArgs {
        min_fragment_length: 10,
        max_fragment_length: 10,
    };

    args.validate()
        .expect("single-length fragment range should be valid");
}

#[test]
fn fragment_length_args_accepts_max_supported_fragment_length() {
    let args = FragmentLengthArgs {
        min_fragment_length: MAX_SUPPORTED_FRAGMENT_LENGTH,
        max_fragment_length: MAX_SUPPORTED_FRAGMENT_LENGTH,
    };

    args.validate()
        .expect("the configured maximum supported fragment length should be valid");
}

#[test]
fn validate_max_soft_clips_accepts_configured_limit() {
    validate_max_soft_clips(MAX_MAX_SOFT_CLIPS)
        .expect("the configured maximum soft-clip limit should be valid");
}

#[test]
fn validate_max_soft_clips_rejects_values_above_configured_limit() {
    let too_large = MAX_MAX_SOFT_CLIPS + 1;

    let error = validate_max_soft_clips(too_large)
        .expect_err("soft-clip limits above the configured maximum should fail");
    let message = error.to_string();

    assert!(
        message.contains(&format!(
            "--max-soft-clips ({too_large}) must be <= {MAX_MAX_SOFT_CLIPS}"
        )),
        "unexpected error: {message}"
    );
}

#[test]
fn min_overlap_fraction_for_window_assignment_uses_any_overlap_thresholds() {
    let span = 4;

    let count_overlap =
        min_overlap_fraction_for_window_assignment(WindowAssigner::CountOverlap, span);
    let any = min_overlap_fraction_for_window_assignment(WindowAssigner::Any, span);

    // One base of overlap across a 4 bp assignment interval is 1/4. A threshold of 1/5 is lower
    // than the smallest nonzero overlap fraction, so both modes keep every nonzero overlap.
    assert_close(count_overlap, 1.0 / 5.0);
    assert_close(any, 1.0 / 5.0);
}

#[test]
fn min_overlap_fraction_for_window_assignment_uses_full_overlap_thresholds() {
    let span = 4;

    let all = min_overlap_fraction_for_window_assignment(WindowAssigner::All, span);
    let midpoint = min_overlap_fraction_for_window_assignment(WindowAssigner::Midpoint, span);

    // For full-span assignment, three of four bases is 3/4 and full overlap is 1. A threshold of
    // 4/5 rejects every partial 4 bp overlap and keeps full containment. Midpoint callers project
    // to a 1 bp interval separately, where a covered midpoint has overlap fraction 1.
    assert_close(all, 4.0 / 5.0);
    assert_close(midpoint, 4.0 / 5.0);
}

#[test]
fn min_overlap_fraction_for_window_assignment_preserves_proportion_threshold() {
    let threshold = 0.375;

    let actual =
        min_overlap_fraction_for_window_assignment(WindowAssigner::Proportion(threshold), 4);

    assert_close(actual, threshold);
}

#[test]
#[should_panic(expected = "window assignment span must be positive")]
fn min_overlap_fraction_for_window_assignment_rejects_zero_span() {
    min_overlap_fraction_for_window_assignment(WindowAssigner::Any, 0);
}

#[test]
fn parse_output_prefix_trims_valid_filename_stem_prefix() {
    let prefix = parse_output_prefix(" sample.v1_alpha-2 ")
        .expect("simple filename-stem prefix should parse");

    assert_eq!(prefix, "sample.v1_alpha-2");
}

#[test]
fn validate_output_prefix_accepts_empty_and_single_separator_prefixes() {
    validate_output_prefix("").expect("empty prefix should keep unprefixed output names valid");
    validate_output_prefix(".").expect("single dot should be allowed");
    validate_output_prefix("_").expect("single underscore should be allowed");
    validate_output_prefix("-").expect("single dash should be allowed");
}

#[test]
fn validate_output_prefix_rejects_parent_directory_sequence() {
    let error = validate_output_prefix("sample..v2")
        .expect_err("double dots should not be accepted in output prefixes");
    let message = error.to_string();

    assert!(
        message.contains("--output-prefix cannot contain '..'"),
        "unexpected error: {message}"
    );
}

#[test]
fn validate_output_prefix_rejects_path_separators() {
    for prefix in ["nested/sample", r"nested\sample"] {
        let error = validate_output_prefix(prefix)
            .expect_err("path separators should not be accepted in output prefixes");
        let message = error.to_string();

        assert!(
            message.contains("--output-prefix cannot contain path separators"),
            "unexpected error for {prefix}: {message}"
        );
    }
}

#[test]
fn validate_output_prefix_rejects_non_filename_stem_symbols() {
    for prefix in ["sample name", "sample:v1", "sample*"] {
        let error = validate_output_prefix(prefix)
            .expect_err("non-token symbols should not be accepted in output prefixes");
        let message = error.to_string();

        assert!(
            message.contains("--output-prefix contains invalid character"),
            "unexpected error for {prefix}: {message}"
        );
    }
}

#[test]
fn fragment_length_args_rejects_above_max_supported_fragment_length() {
    let too_large = MAX_SUPPORTED_FRAGMENT_LENGTH + 1;
    let args = FragmentLengthArgs {
        min_fragment_length: 10,
        max_fragment_length: too_large,
    };

    let error = args
        .validate()
        .expect_err("fragment lengths above the supported cap should fail");
    let message = error.to_string();

    assert!(
        message.contains(&format!(
            "--max-fragment-length ({too_large}) must be <= {MAX_SUPPORTED_FRAGMENT_LENGTH}"
        )),
        "unexpected error: {message}"
    );
}

#[test]
fn parse_length_bins_accepts_bin_ending_after_max_supported_fragment_length() {
    let spec = format!(
        "{}:{}:1",
        MAX_SUPPORTED_FRAGMENT_LENGTH,
        MAX_SUPPORTED_FRAGMENT_LENGTH + 1
    );

    let bins = parse_length_bins(Some(&spec), 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect("the final exclusive edge may be max supported length + 1");

    assert_eq!(
        bins.to_edges(),
        vec![
            MAX_SUPPORTED_FRAGMENT_LENGTH,
            MAX_SUPPORTED_FRAGMENT_LENGTH + 1
        ]
    );
}

#[test]
fn parse_length_bins_rejects_range_past_max_supported_fragment_length() {
    let invalid_end = MAX_SUPPORTED_FRAGMENT_LENGTH + 2;
    let spec = format!("10:{invalid_end}:1");

    let error = parse_length_bins(Some(&spec), 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("length-bin ranges beyond the supported cap should fail before expansion");
    let message = error.to_string();

    assert!(
        message.contains(&format!(
            "length-bins end ({invalid_end}) must be <= max fragment length + 1 ({})",
            MAX_SUPPORTED_FRAGMENT_LENGTH + 1
        )),
        "unexpected error: {message}"
    );
}

#[test]
fn parse_length_bins_rejects_step_larger_than_max_fragment_length() {
    let spec = format!("30:1001:{}", u32::MAX);

    let error = parse_length_bins(Some(&spec), 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("steps above the maximum fragment length should fail");
    let message = error.to_string();

    assert!(
        message.contains(&format!(
            "length-bins step ({}) must be <= max fragment length ({})",
            u32::MAX,
            MAX_SUPPORTED_FRAGMENT_LENGTH
        )),
        "unexpected error: {message}"
    );
}

#[test]
fn parse_length_bins_rejects_non_numeric_step_with_plain_message() {
    let spec = "30:1001:1.5";

    let error = parse_length_bins(Some(spec), 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("non-integer steps should fail with a field-specific message");
    let message = format!("{error:#}");

    assert!(
        message.contains("length-bins step must be a positive whole number"),
        "unexpected error: {message}"
    );
}

#[test]
fn parse_length_bins_rejects_negative_step_with_plain_message() {
    let spec = "30:1001:-1";

    let error = parse_length_bins(Some(spec), 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("negative steps should fail with a field-specific message");
    let message = format!("{error:#}");

    assert!(
        message.contains("length-bins step must be a positive whole number"),
        "unexpected error: {message}"
    );
}

#[test]
fn parse_length_bins_rejects_unparseably_large_step_with_plain_message() {
    let spec = "30:1001:999999999999999999999999999999999999999999999999999";

    let error = parse_length_bins(Some(spec), 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("unparseably large steps should fail with a field-specific message");
    let message = format!("{error:#}");

    assert!(
        message.contains("length-bins step is too large"),
        "unexpected error: {message}"
    );
}

#[test]
fn resolve_length_bin_edges_accepts_compact_range_spec() {
    let raw_values = vec!["30:101:10".to_string()];

    let edges = resolve_length_bin_edges(&raw_values, 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect("compact length-bin range should resolve");

    assert_eq!(edges, vec![30, 40, 50, 60, 70, 80, 90, 100, 101]);
}

#[test]
fn resolve_length_bin_edges_accepts_explicit_edge_values() {
    let raw_values = vec!["10".to_string(), "151".to_string(), "221".to_string()];

    let edges = resolve_length_bin_edges(&raw_values, 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect("explicit length-bin edges should resolve");

    assert_eq!(edges, vec![10, 151, 221]);
}

#[test]
fn resolve_length_bin_edges_accepts_final_exclusive_edge_after_max_supported_length() {
    let raw_values = vec![
        MAX_SUPPORTED_FRAGMENT_LENGTH.to_string(),
        (MAX_SUPPORTED_FRAGMENT_LENGTH + 1).to_string(),
    ];

    let edges = resolve_length_bin_edges(&raw_values, 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect("explicit final edge may be max supported length + 1");

    assert_eq!(
        edges,
        vec![
            MAX_SUPPORTED_FRAGMENT_LENGTH,
            MAX_SUPPORTED_FRAGMENT_LENGTH + 1
        ]
    );
}

#[test]
fn resolve_length_bin_edges_rejects_explicit_edge_after_final_exclusive_cap() {
    let too_large = MAX_SUPPORTED_FRAGMENT_LENGTH + 2;
    let raw_values = vec!["10".to_string(), too_large.to_string()];

    let error = resolve_length_bin_edges(&raw_values, 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("explicit edges above the final exclusive cap should fail");
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
fn resolve_length_bin_edges_rejects_non_increasing_edges() {
    let raw_values = vec!["10".to_string(), "151".to_string(), "151".to_string()];

    let error = resolve_length_bin_edges(&raw_values, 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("non-increasing edges should fail");
    let message = error.to_string();

    assert!(
        message.contains("length bin edges must be strictly increasing"),
        "unexpected error: {message}"
    );
}

#[test]
fn resolve_length_bin_edges_rejects_edge_below_minimum() {
    let raw_values = vec!["9".to_string(), "151".to_string()];

    let error = resolve_length_bin_edges(&raw_values, 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("edges below the minimum allowed length should fail");
    let message = error.to_string();

    assert!(
        message.contains("length bin edges must be >= 10"),
        "unexpected error: {message}"
    );
}

#[test]
fn resolve_length_bin_edges_rejects_non_numeric_explicit_edge_with_plain_message() {
    let raw_values = vec!["10".to_string(), "151.5".to_string()];

    let error = resolve_length_bin_edges(&raw_values, 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("non-integer explicit edges should fail with a field-specific message");
    let message = format!("{error:#}");

    assert!(
        message.contains("length-bins edge must be a positive whole number"),
        "unexpected error: {message}"
    );
}

#[test]
fn resolve_length_bin_edges_rejects_negative_explicit_edge_with_plain_message() {
    let raw_values = vec!["10".to_string(), "-151".to_string()];

    let error = resolve_length_bin_edges(&raw_values, 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("negative explicit edges should fail with a field-specific message");
    let message = format!("{error:#}");

    assert!(
        message.contains("length-bins edge must be a positive whole number"),
        "unexpected error: {message}"
    );
}

#[test]
fn resolve_length_bin_edges_rejects_start_end_list_format() {
    let raw_values = vec!["30-80,80-150".to_string()];

    let error = resolve_length_bin_edges(&raw_values, 10, MAX_SUPPORTED_FRAGMENT_LENGTH)
        .expect_err("start-end lists are intentionally unsupported");
    let message = error.to_string();

    assert!(
        message.contains("explicit start-end lists are not supported"),
        "unexpected error: {message}"
    );
}

#[test]
fn resolve_chromosomes_all_uses_contig_source_order() {
    let mut chrom_sizes = tempfile::NamedTempFile::new().unwrap();
    writeln!(chrom_sizes, "chrB\t20").unwrap();
    writeln!(chrom_sizes, "chrA\t10").unwrap();

    let args = ChromosomeArgs {
        chromosomes: Some(vec!["all".to_string()]),
        chromosomes_file: None,
    };

    let chromosomes = args
        .resolve_chromosomes(Some(ContigSource::chrom_sizes(chrom_sizes.path())))
        .expect("all should resolve from the provided contig source");

    assert_eq!(chromosomes, vec!["chrB", "chrA"]);
}

#[test]
fn resolve_chromosomes_all_without_contig_source_fails() {
    let args = ChromosomeArgs {
        chromosomes: Some(vec!["all".to_string()]),
        chromosomes_file: None,
    };

    let error = args
        .resolve_chromosomes(None)
        .expect_err("all should require a contig source");
    let message = error.to_string();

    assert!(
        message.contains("`--chromosomes all` requires a contig source"),
        "unexpected error: {message}"
    );
}

#[test]
fn apply_gc_args_rejects_gc_file_without_ref_2bit() {
    let args = ApplyGCArgs {
        gc_file: Some(PathBuf::from("gc_bias_correction.zarr")),
        gc_tag: None,
        neutralize_invalid_gc: false,
    };

    let error = args
        .validate(None)
        .expect_err("GC file should require a reference genome");
    let message = error.to_string();

    assert!(
        message.contains("--gc-file requires --ref-2bit"),
        "unexpected error: {message}"
    );
}

#[test]
fn apply_gc_args_allows_gc_tag_without_ref_2bit() {
    let args = ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    };

    args.validate(None)
        .expect("GC-tag correction should not require a reference genome");
}

#[test]
fn apply_gc_args_allows_lowercase_local_gc_tag() {
    let args = ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("gc".to_string()),
        neutralize_invalid_gc: false,
    };

    args.validate(None)
        .expect("lowercase local AUX tags should be valid SAM/BAM tags");
}

#[test]
fn apply_gc_args_rejects_overlong_gc_tag() {
    let args = ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GCP".to_string()),
        neutralize_invalid_gc: false,
    };

    let error = args
        .validate(None)
        .expect_err("overlong AUX tags should fail before BAM lookup");
    let message = error.to_string();

    assert!(
        message.contains("exactly two ASCII bytes"),
        "unexpected error: {message}"
    );
}

#[test]
fn parse_sam_aux_tag_name_rejects_invalid_shapes() {
    // `GC` is the common external GC-weight tag. `cw` is a local/private tag shape.
    assert_eq!(parse_sam_aux_tag_name("GC").unwrap(), "GC");
    assert_eq!(parse_sam_aux_tag_name("cw").unwrap(), "cw");

    for invalid_tag in ["", "G", "GCP", "1G", "G_", "åC"] {
        let error = parse_sam_aux_tag_name(invalid_tag)
            .expect_err("invalid SAM/BAM AUX tag shapes should be rejected");
        let message = error.to_string();
        assert!(
            message.contains("SAM/BAM AUX tag"),
            "unexpected error for {invalid_tag:?}: {message}"
        );
    }
}

#[test]
fn apply_gc_file_only_args_rejects_gc_file_without_ref_2bit() {
    let args = ApplyGCArgFileOnly {
        gc_file: Some(PathBuf::from("gc_bias_correction.zarr")),
        neutralize_invalid_gc: false,
    };

    let error = args
        .validate(None)
        .expect_err("GC file should require a reference genome");
    let message = error.to_string();

    assert!(
        message.contains("--gc-file requires --ref-2bit"),
        "unexpected error: {message}"
    );
}

#[test]
fn apply_gc_file_only_args_allows_gc_file_with_ref_2bit() {
    let args = ApplyGCArgFileOnly {
        gc_file: Some(PathBuf::from("gc_bias_correction.zarr")),
        neutralize_invalid_gc: false,
    };

    args.validate(Some(Path::new("reference.2bit")))
        .expect("GC file should be valid when a reference genome is configured");
}
