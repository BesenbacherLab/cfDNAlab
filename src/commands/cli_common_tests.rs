use super::{ApplyGCArgFileOnly, ApplyGCArgs, FragmentLengthArgs};
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
