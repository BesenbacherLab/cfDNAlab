use super::{
    ApplyGCArgFileOnly, ApplyGCArgs, FragmentLengthArgs, MAX_SUPPORTED_FRAGMENT_LENGTH,
    parse_length_bins,
};
use std::path::{Path, PathBuf};

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

    assert_eq!(bins.to_edges(), vec![MAX_SUPPORTED_FRAGMENT_LENGTH, MAX_SUPPORTED_FRAGMENT_LENGTH + 1]);
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
fn apply_gc_args_rejects_gc_file_without_ref_2bit() {
    let args = ApplyGCArgs {
        gc_file: Some(PathBuf::from("gc_bias_correction.npz")),
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
fn apply_gc_file_only_args_rejects_gc_file_without_ref_2bit() {
    let args = ApplyGCArgFileOnly {
        gc_file: Some(PathBuf::from("gc_bias_correction.npz")),
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
        gc_file: Some(PathBuf::from("gc_bias_correction.npz")),
        neutralize_invalid_gc: false,
    };

    args.validate(Some(Path::new("reference.2bit")))
        .expect("GC file should be valid when a reference genome is configured");
}
