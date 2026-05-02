use super::{
    WindowMetadataEntry, configured_max_fragment_reach_bp, reorder_bed_outputs_by_original_index,
};
use crate::commands::cli_common::{ChromosomeArgs, IOCArgs};
use crate::commands::lengths::counting::{LengthAxis, LengthCounts};
use crate::commands::lengths::config::LengthsConfig;
use crate::shared::{clip_mode::ClipMode, indel_mode::IndelMode};
use std::sync::Arc;

fn counts_with_value(value: f64) -> LengthCounts {
    let axis = Arc::new(LengthAxis::new(vec![10, 11]).expect("test axis should be valid"));
    let mut counts = LengthCounts::new(axis);
    counts.counts[0] = value;
    counts
}

fn bin_info_entry(original_index: u64) -> WindowMetadataEntry {
    ("chr1".to_string(), 0, 10, original_index, 0.0)
}

fn test_config() -> LengthsConfig {
    LengthsConfig::new(
        IOCArgs {
            bam: "input.bam".into(),
            output_dir: ".".into(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    )
}

#[test]
fn configured_max_fragment_reach_uses_larger_active_cap() {
    // Arrange:
    // The length axis allows adjusted lengths through 100. Deletion adjustment can make the
    // aligned reference span longer, while clip adjustment can move the assignment interval
    // outside the aligned span. The reach uses the larger active cap, not their sum.
    let length_axis = LengthAxis::new(vec![30, 101]).expect("test axis should be valid");
    let mut config = test_config();
    config.set_indel_mode(IndelMode::Adjust);
    config.clip_mode = ClipMode::Adjust;
    config.max_soft_clips = 7;
    config.max_deletion_bases = 11;

    // Act
    let reach_bp = configured_max_fragment_reach_bp(&config, &length_axis);

    // Assert
    assert_eq!(reach_bp, 111);
}

#[test]
fn configured_max_fragment_reach_ignores_inactive_caps() {
    // Arrange:
    // In the default modes, neither cap changes the reference-coordinate reach because lengths
    // are already based on aligned reference spans.
    let length_axis = LengthAxis::new(vec![30, 101]).expect("test axis should be valid");
    let mut config = test_config();
    config.max_soft_clips = 7;
    config.max_deletion_bases = 11;

    // Act
    let reach_bp = configured_max_fragment_reach_bp(&config, &length_axis);

    // Assert
    assert_eq!(reach_bp, 100);
}

#[test]
fn bed_output_reorder_rejects_metadata_count_length_mismatch() {
    // Arrange
    let mut bin_info = vec![bin_info_entry(0)];
    let mut all_bins = vec![counts_with_value(1.0), counts_with_value(2.0)];

    // Act
    let error = reorder_bed_outputs_by_original_index(&mut bin_info, &mut all_bins)
        .expect_err("mismatched BED metadata and count vectors should fail before zipping");
    let message = error.to_string();

    // Assert
    assert!(
        message.contains("BED metadata entries (1) did not match length count vectors (2)"),
        "unexpected error: {message}"
    );
}

#[test]
fn bed_output_reorder_keeps_counts_paired_with_original_window_index() {
    // Arrange
    let mut bin_info = vec![bin_info_entry(7), bin_info_entry(3)];
    let mut all_bins = vec![counts_with_value(7.0), counts_with_value(3.0)];

    // Act
    reorder_bed_outputs_by_original_index(&mut bin_info, &mut all_bins)
        .expect("matching BED metadata and count vectors should reorder together");

    // Assert
    let original_indices: Vec<u64> = bin_info.iter().map(|entry| entry.3).collect();
    let first_column_counts: Vec<f64> = all_bins.iter().map(|counts| counts.counts[0]).collect();
    assert_eq!(original_indices, vec![3, 7]);
    assert_eq!(first_column_counts, vec![3.0, 7.0]);
}
