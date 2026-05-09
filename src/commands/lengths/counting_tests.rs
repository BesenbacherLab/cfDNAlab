use super::{LengthAxis, LengthCounts, stack_length_counts};
use crate::shared::constants::MAX_SUPPORTED_FRAGMENT_LENGTH;
use std::sync::Arc;

fn axis(edges: Vec<u32>) -> Arc<LengthAxis> {
    Arc::new(LengthAxis::new(edges).expect("test length axis should be valid"))
}

#[test]
fn length_axis_maps_lengths_to_half_open_bins() {
    // Arrange: two wider bins, [30,40) and [40,50).
    let length_axis = LengthAxis::new(vec![30, 40, 50]).expect("axis should be valid");

    // Assert: starts and ends follow half-open interval semantics.
    assert_eq!(length_axis.bin_index(29), None);
    assert_eq!(length_axis.bin_index(30), Some(0));
    assert_eq!(length_axis.bin_index(39), Some(0));
    assert_eq!(length_axis.bin_index(40), Some(1));
    assert_eq!(length_axis.bin_index(49), Some(1));
    assert_eq!(length_axis.bin_index(50), None);
    assert_eq!(length_axis.min_fragment_length(), 30);
    assert_eq!(length_axis.max_fragment_length(), 49);
    assert!(!length_axis.is_single_bp_bins());
}

#[test]
fn length_axis_detects_single_bp_bins() {
    let length_axis = LengthAxis::new(vec![30, 31, 32]).expect("axis should be valid");

    assert!(length_axis.is_single_bp_bins());
}

#[test]
fn length_axis_rejects_edges_below_minimum_fragment_length() {
    let error =
        LengthAxis::new(vec![9, 10]).expect_err("length edges below 10 bp should fail");

    assert!(
        error.to_string().contains("length bin edges must be >= 10"),
        "unexpected error: {error}"
    );
}

#[test]
fn length_axis_rejects_non_increasing_edges() {
    let error =
        LengthAxis::new(vec![10, 20, 20]).expect_err("non-increasing edges should fail");

    assert!(
        error
            .to_string()
            .contains("length bin edges must be strictly increasing"),
        "unexpected error: {error}"
    );
}

#[test]
fn length_axis_rejects_final_edge_past_supported_exclusive_cap() {
    let error = LengthAxis::new(vec![10, MAX_SUPPORTED_FRAGMENT_LENGTH + 2])
        .expect_err("final edge past supported cap should fail");

    assert!(
        error.to_string().contains(&format!(
            "length bin edges must be <= {}",
            MAX_SUPPORTED_FRAGMENT_LENGTH + 1
        )),
        "unexpected error: {error}"
    );
}

#[test]
fn length_counts_increment_wider_bin_by_absolute_length() {
    // Arrange: lengths 35 and 39 both belong to [30,40).
    let mut counts = LengthCounts::new(axis(vec![30, 40, 50]));

    // Act
    counts
        .incr_weighted(35, 1.5)
        .expect("length 35 should be in range");
    counts
        .incr_weighted(39, 2.5)
        .expect("length 39 should be in range");
    counts
        .incr_weighted(40, 3.0)
        .expect("length 40 should be in range");

    // Assert
    assert_eq!(counts.counts, vec![4.0, 3.0]);
    assert_eq!(counts.get(35), Some(4.0));
    assert_eq!(counts.get(39), Some(4.0));
    assert_eq!(counts.get(40), Some(3.0));
}

#[test]
fn length_counts_errors_for_length_outside_axis() {
    let mut counts = LengthCounts::new(axis(vec![10, 20]));

    let error = counts
        .incr_weighted(20, 1.0)
        .expect_err("final exclusive edge should be out of range");

    assert!(
        error
            .to_string()
            .contains("fragment length 20 did not map to any configured length bin"),
        "unexpected error: {error}"
    );
}

#[test]
fn length_counts_merge_rejects_different_axes_with_same_width() {
    let mut left = LengthCounts::new(axis(vec![10, 20]));
    let right = LengthCounts::new(axis(vec![20, 30]));

    let error = left
        .merge_from(&right)
        .expect_err("different axes with the same width should not merge");

    assert!(
        error.to_string().contains("incompatible LengthCounts"),
        "unexpected error: {error}"
    );
}

#[test]
fn length_counts_merge_accepts_same_edges_from_different_axes() {
    let mut left = LengthCounts::new(axis(vec![10, 20]));
    let mut right = LengthCounts::new(axis(vec![10, 20]));
    right
        .incr_weighted(10, 2.0)
        .expect("length 10 should be in range");

    left.merge_from(&right)
        .expect("same edges should be compatible even with distinct Arc allocations");

    assert_eq!(left.counts, vec![2.0]);
}

#[test]
fn stack_length_counts_rejects_empty_input() {
    let error = stack_length_counts(&[]).expect_err("empty counters should fail");

    assert!(
        error
            .to_string()
            .contains("stack_length_counts requires at least one counter"),
        "unexpected error: {error}"
    );
}

#[test]
fn stack_length_counts_rejects_different_axes_with_same_width() {
    let left = LengthCounts::new(axis(vec![10, 20]));
    let right = LengthCounts::new(axis(vec![20, 30]));

    let error = stack_length_counts(&[left, right])
        .expect_err("same-width counters with different axes should fail");

    assert!(
        error
            .to_string()
            .contains("length count entry 1 has incompatible length axis"),
        "unexpected error: {error}"
    );
}
