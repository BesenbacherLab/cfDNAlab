use super::*;
use crate::commands::ends::config_structs::BaseQualityComparisonOp;
use crate::shared::gc_tag::GcTagValue;
use rust_htslib::bam::record::{Cigar, CigarString, Record};

fn base_read_info(qualities: &[u8]) -> EndReadInfo {
    base_read_info_with_clips(qualities, 0, 0)
}

fn base_read_info_with_clips(
    qualities: &[u8],
    left_soft_clip_bp: u32,
    right_soft_clip_bp: u32,
) -> EndReadInfo {
    EndReadInfo {
        tid: 0,
        interval: Interval::new(10, 24).expect("valid interval"),
        is_reverse: false,
        left_soft_clip_bp,
        right_soft_clip_bp,
        left_motif_has_indels: false,
        right_motif_has_indels: false,
        has_hard_clip: false,
        seq: vec![b'A'; qualities.len()],
        qualities: Some(qualities.to_vec()),
        gc_tag: GcTagValue::default(),
    }
}

fn base_end(boundary_pos: u32) -> ResolvedFragmentEnd {
    ResolvedFragmentEnd {
        boundary_pos,
        inside_bases: vec![b'A', b'C', b'G'],
        inside_reference_validation_bp: 3,
    }
}

#[test]
fn apply_base_quality_filters_distinguishes_min_mean_and_max_for_end_filters() {
    // Arrange: one end with qualities [35, 20, 40].
    //
    // Mental derivation:
    // - min = 20, so `min in end >= 30` fails and should drop the end
    // - mean = (35 + 20 + 40) / 3 = 31.666..., so `mean in end >= 30` passes
    // - max = 40, so `max in end < 30` fails and should drop the end
    let read = base_read_info(&[35, 20, 40]);
    let left_end = Some(base_end(10));

    let min_filter = [BaseQualityFilter {
        aggregation: BaseQualityAggregation::Min,
        scope: BaseQualityFilterScope::End,
        op: BaseQualityComparisonOp::Ge,
        threshold: 30.0,
    }];
    let mean_filter = [BaseQualityFilter {
        aggregation: BaseQualityAggregation::Mean,
        scope: BaseQualityFilterScope::End,
        op: BaseQualityComparisonOp::Ge,
        threshold: 30.0,
    }];
    let max_filter = [BaseQualityFilter {
        aggregation: BaseQualityAggregation::Max,
        scope: BaseQualityFilterScope::End,
        op: BaseQualityComparisonOp::Lt,
        threshold: 30.0,
    }];

    // Act / Assert
    assert!(
        apply_base_quality_filters(
            &read,
            &read,
            left_end.clone(),
            None,
            ClipStrategy::Aligned,
            3,
            &min_filter
        )
        .is_none()
    );
    let mean_result = apply_base_quality_filters(
        &read,
        &read,
        left_end.clone(),
        None,
        ClipStrategy::Aligned,
        3,
        &mean_filter,
    );
    assert!(mean_result.is_some());
    let (kept_left_end, kept_right_end) = mean_result.expect("mean filter should keep the end");
    assert!(kept_right_end.is_none());
    let kept_left_end = kept_left_end.expect("left end should stay present");
    let original_left_end = left_end
        .clone()
        .expect("fixture should contain a left end");
    assert_eq!(kept_left_end.boundary_pos, original_left_end.boundary_pos);
    assert_eq!(kept_left_end.inside_bases, original_left_end.inside_bases);
    assert_eq!(
        kept_left_end.inside_reference_validation_bp,
        original_left_end.inside_reference_validation_bp
    );
    assert!(
        apply_base_quality_filters(
            &read,
            &read,
            left_end,
            None,
            ClipStrategy::Aligned,
            3,
            &max_filter
        )
        .is_none()
    );
}

#[test]
fn apply_base_quality_filters_is_a_no_op_when_no_filters_are_present() {
    // Arrange
    let left_read = base_read_info(&[35, 20, 40]);
    let right_read = base_read_info(&[10, 20, 30]);
    let left_end = Some(base_end(10));
    let right_end = Some(base_end(20));

    // Act
    let result = apply_base_quality_filters(
        &left_read,
        &right_read,
        left_end.clone(),
        right_end.clone(),
        ClipStrategy::Aligned,
        3,
        &[],
    );

    // Assert
    assert!(result.is_some());
    let (kept_left_end, kept_right_end) = result.expect("no filters should preserve both ends");
    let kept_left_end = kept_left_end.expect("left end should stay present");
    let kept_right_end = kept_right_end.expect("right end should stay present");
    let original_left_end = left_end.expect("fixture should contain a left end");
    let original_right_end = right_end.expect("fixture should contain a right end");
    assert_eq!(kept_left_end.boundary_pos, original_left_end.boundary_pos);
    assert_eq!(kept_left_end.inside_bases, original_left_end.inside_bases);
    assert_eq!(
        kept_left_end.inside_reference_validation_bp,
        original_left_end.inside_reference_validation_bp
    );
    assert_eq!(kept_right_end.boundary_pos, original_right_end.boundary_pos);
    assert_eq!(kept_right_end.inside_bases, original_right_end.inside_bases);
    assert_eq!(
        kept_right_end.inside_reference_validation_bp,
        original_right_end.inside_reference_validation_bp
    );
}

#[test]
fn fragment_scope_filters_run_before_end_scope_filters_remove_failing_ends() {
    // Arrange: left qualities [40], right qualities [10].
    //
    // Mental derivation:
    // - fragment mean = (40 + 10) / 2 = 25, so `mean in fragment >= 30` fails
    // - left end still passes `min in end >= 30`
    // - because fragment filters are applied before end filters, the fragment must be dropped
    let left_read = base_read_info(&[40]);
    let right_read = base_read_info(&[10]);
    let filters = [
        BaseQualityFilter {
            aggregation: BaseQualityAggregation::Min,
            scope: BaseQualityFilterScope::End,
            op: BaseQualityComparisonOp::Ge,
            threshold: 30.0,
        },
        BaseQualityFilter {
            aggregation: BaseQualityAggregation::Mean,
            scope: BaseQualityFilterScope::Fragment,
            op: BaseQualityComparisonOp::Ge,
            threshold: 30.0,
        },
    ];

    // Act / Assert
    assert!(
        apply_base_quality_filters(
            &left_read,
            &right_read,
            Some(base_end(10)),
            Some(base_end(20)),
            ClipStrategy::Aligned,
            1,
            &filters
        )
        .is_none()
    );
}

#[test]
fn apply_base_quality_filters_drop_the_fragment_when_fragment_scope_passes_but_both_end_filters_fail()
{
    // Arrange: both ends have quality 20.
    //
    // Mental derivation:
    // - fragment mean = (20 + 20) / 2 = 20, so `mean in fragment >= 20` passes
    // - `min in end >= 30` fails for both ends
    // - once both ends are dropped, the helper must return `None`
    let left_read = base_read_info(&[20]);
    let right_read = base_read_info(&[20]);
    let filters = [
        BaseQualityFilter {
            aggregation: BaseQualityAggregation::Mean,
            scope: BaseQualityFilterScope::Fragment,
            op: BaseQualityComparisonOp::Ge,
            threshold: 20.0,
        },
        BaseQualityFilter {
            aggregation: BaseQualityAggregation::Min,
            scope: BaseQualityFilterScope::End,
            op: BaseQualityComparisonOp::Ge,
            threshold: 30.0,
        },
    ];

    // Act / Assert
    assert!(
        apply_base_quality_filters(
            &left_read,
            &right_read,
            Some(base_end(10)),
            Some(base_end(20)),
            ClipStrategy::Aligned,
            1,
            &filters
        )
        .is_none()
    );
}

#[test]
fn fragment_scope_filters_distinguish_min_mean_and_max_for_k_inside_gt_one() {
    // Arrange: the fragment has left qualities [40, 35, 30] and right qualities [20, 20, 20].
    //
    // Mental derivation across the raw candidate fragment:
    // - min = 20, so `min in fragment > 20` fails
    // - mean = (40 + 35 + 30 + 20 + 20 + 20) / 6 = 27.5, so `mean in fragment >= 27.5` passes
    // - max = 40, so `max in fragment < 35` fails
    let left_read = base_read_info(&[40, 35, 30]);
    let right_read = base_read_info(&[20, 20, 20]);
    let left_end = Some(base_end(10));
    let right_end = Some(base_end(20));

    let min_filter = [BaseQualityFilter {
        aggregation: BaseQualityAggregation::Min,
        scope: BaseQualityFilterScope::Fragment,
        op: BaseQualityComparisonOp::Gt,
        threshold: 20.0,
    }];
    let mean_filter = [BaseQualityFilter {
        aggregation: BaseQualityAggregation::Mean,
        scope: BaseQualityFilterScope::Fragment,
        op: BaseQualityComparisonOp::Ge,
        threshold: 27.5,
    }];
    let max_filter = [BaseQualityFilter {
        aggregation: BaseQualityAggregation::Max,
        scope: BaseQualityFilterScope::Fragment,
        op: BaseQualityComparisonOp::Lt,
        threshold: 35.0,
    }];

    // Act / Assert
    assert!(
        apply_base_quality_filters(
            &left_read,
            &right_read,
            left_end.clone(),
            right_end.clone(),
            ClipStrategy::Aligned,
            3,
            &min_filter
        )
        .is_none()
    );
    let mean_result = apply_base_quality_filters(
        &left_read,
        &right_read,
        left_end.clone(),
        right_end.clone(),
        ClipStrategy::Aligned,
        3,
        &mean_filter,
    );
    assert!(mean_result.is_some());
    let (kept_left_end, kept_right_end) =
        mean_result.expect("mean fragment filter should keep both ends");
    let kept_left_end = kept_left_end.expect("left end should stay present");
    let kept_right_end = kept_right_end.expect("right end should stay present");
    let original_left_end = left_end
        .clone()
        .expect("fixture should contain a left end");
    let original_right_end = right_end
        .clone()
        .expect("fixture should contain a right end");
    assert_eq!(kept_left_end.boundary_pos, original_left_end.boundary_pos);
    assert_eq!(kept_left_end.inside_bases, original_left_end.inside_bases);
    assert_eq!(
        kept_left_end.inside_reference_validation_bp,
        original_left_end.inside_reference_validation_bp
    );
    assert_eq!(kept_right_end.boundary_pos, original_right_end.boundary_pos);
    assert_eq!(kept_right_end.inside_bases, original_right_end.inside_bases);
    assert_eq!(
        kept_right_end.inside_reference_validation_bp,
        original_right_end.inside_reference_validation_bp
    );
    assert!(
        apply_base_quality_filters(
            &left_read,
            &right_read,
            left_end,
            right_end,
            ClipStrategy::Aligned,
            3,
            &max_filter
        )
        .is_none()
    );
}

#[test]
fn extract_inside_qualities_respects_clip_strategy_for_k_inside_gt_one() {
    // Arrange: 2S10M2S with per-base qualities 10, 20, ..., 140.
    //
    // Mental derivation for k_inside=3:
    // - aligned left skips the left clips -> [30, 40, 50]
    // - aligned right skips the right clips -> [100, 110, 120]
    // - raw left starts at the raw read edge -> [10, 20, 30]
    // - raw right ends at the raw read edge -> [120, 130, 140]
    let read = base_read_info_with_clips(
        &[10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140],
        2,
        2,
    );

    // Act / Assert
    assert_eq!(
        extract_inside_qualities(&read, FragmentEndSide::Left, ClipStrategy::Aligned, 3),
        Some(vec![30, 40, 50])
    );
    assert_eq!(
        extract_inside_qualities(&read, FragmentEndSide::Right, ClipStrategy::Aligned, 3),
        Some(vec![100, 110, 120])
    );
    assert_eq!(
        extract_inside_qualities(
            &read,
            FragmentEndSide::Left,
            ClipStrategy::RawAlignedBoundary,
            3
        ),
        Some(vec![10, 20, 30])
    );
    assert_eq!(
        extract_inside_qualities(
            &read,
            FragmentEndSide::Right,
            ClipStrategy::RawAlignedBoundary,
            3
        ),
        Some(vec![120, 130, 140])
    );
    assert_eq!(
        extract_inside_qualities(
            &read,
            FragmentEndSide::Left,
            ClipStrategy::RawShiftedBoundary,
            3
        ),
        Some(vec![10, 20, 30])
    );
    assert_eq!(
        extract_inside_qualities(
            &read,
            FragmentEndSide::Right,
            ClipStrategy::RawShiftedBoundary,
            3
        ),
        Some(vec![120, 130, 140])
    );
}

#[test]
fn inside_slice_bounds_return_the_expected_indices_for_each_strategy_and_side() {
    // Arrange: len=14 with 2 soft-clipped bases on the left and 3 on the right.
    //
    // Mental derivation for k_inside=3:
    // - aligned/skip left starts after the left clips -> [2, 5)
    // - aligned/skip right ends before the right clips -> [8, 11)
    // - raw-aligned/raw-shifted left starts at the raw read edge -> [0, 3)
    // - raw-aligned/raw-shifted right ends at the raw read edge -> [11, 14)
    let len = 14;
    let left_soft_clip_bp = 2;
    let right_soft_clip_bp = 3;
    let k_inside = 3;

    // Act / Assert
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Left,
            ClipStrategy::Aligned,
            k_inside
        ),
        Some((2, 5))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Right,
            ClipStrategy::Aligned,
            k_inside
        ),
        Some((8, 11))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Left,
            ClipStrategy::Skip,
            k_inside
        ),
        Some((2, 5))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Right,
            ClipStrategy::Skip,
            k_inside
        ),
        Some((8, 11))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Left,
            ClipStrategy::RawAlignedBoundary,
            k_inside
        ),
        Some((0, 3))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Right,
            ClipStrategy::RawAlignedBoundary,
            k_inside
        ),
        Some((11, 14))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Left,
            ClipStrategy::RawShiftedBoundary,
            k_inside
        ),
        Some((0, 3))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Right,
            ClipStrategy::RawShiftedBoundary,
            k_inside
        ),
        Some((11, 14))
    );
}

#[test]
fn inside_slice_bounds_return_zero_width_slices_when_k_inside_is_zero() {
    // Arrange: zero inside bases should produce an empty half-open slice at the side-specific
    // start/end position rather than failing.
    let len = 14;
    let left_soft_clip_bp = 2;
    let right_soft_clip_bp = 3;

    // Act / Assert
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Left,
            ClipStrategy::Aligned,
            0
        ),
        Some((2, 2))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Right,
            ClipStrategy::Aligned,
            0
        ),
        Some((11, 11))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Left,
            ClipStrategy::RawShiftedBoundary,
            0
        ),
        Some((0, 0))
    );
    assert_eq!(
        inside_slice_bounds(
            len,
            left_soft_clip_bp,
            right_soft_clip_bp,
            FragmentEndSide::Right,
            ClipStrategy::RawShiftedBoundary,
            0
        ),
        Some((14, 14))
    );
}

#[test]
fn inside_slice_bounds_return_none_when_the_requested_inside_span_does_not_fit() {
    // Arrange / Assert
    //
    // Mental derivation:
    // - len=4, left clip=2, aligned left start is 2, so k_inside=3 would need [2, 5) and fail
    // - len=4, right clip=2, aligned right end is 2, so k_inside=3 would need [-1, 2) and fail
    // - len=2 with raw boundaries cannot serve k_inside=3 from either end
    assert_eq!(
        inside_slice_bounds(4, 2, 0, FragmentEndSide::Left, ClipStrategy::Aligned, 3),
        None
    );
    assert_eq!(
        inside_slice_bounds(4, 0, 2, FragmentEndSide::Right, ClipStrategy::Aligned, 3),
        None
    );
    assert_eq!(
        inside_slice_bounds(
            2,
            0,
            0,
            FragmentEndSide::Left,
            ClipStrategy::RawAlignedBoundary,
            3
        ),
        None
    );
    assert_eq!(
        inside_slice_bounds(
            2,
            0,
            0,
            FragmentEndSide::Right,
            ClipStrategy::RawShiftedBoundary,
            3
        ),
        None
    );
}

#[test]
fn from_record_with_gc_tag_skips_loading_qualities_when_not_requested() {
    // Arrange: qualities of 255 denote missing QVs in BAM, but this should not matter when no
    // base-quality filter is active and the hot path intentionally avoids loading them.
    let mut record = Record::new();
    record.set_tid(0);
    record.set_pos(10);
    record.set(
        b"missing_quals_allowed_without_filter",
        Some(&CigarString(vec![Cigar::Match(4)])),
        b"ACGT",
        &[255, 255, 255, 255],
    );

    // Act
    let read_info =
        EndReadInfo::from_record_with_gc_tag(&record, None, ClipStrategy::Aligned, 1, false)
            .expect("qualities should stay unloaded when BQ filters are absent");

    // Assert
    assert_eq!(read_info.qualities, None);
}

#[test]
fn from_record_with_gc_tag_errors_on_missing_base_qualities_when_requested() {
    // Arrange: BAM encodes missing qualities as 255 placeholders.
    let mut record = Record::new();
    record.set_tid(0);
    record.set_pos(10);
    record.set(
        b"missing_quals",
        Some(&CigarString(vec![Cigar::Match(4)])),
        b"ACGT",
        &[255, 255, 255, 255],
    );

    // Act
    let error = EndReadInfo::from_record_with_gc_tag(&record, None, ClipStrategy::Aligned, 1, true)
        .expect_err("missing BAM qualities should fail when BQ filters are active");

    // Assert
    assert!(error.to_string().contains("missing base qualities"));
    assert!(error.to_string().contains("--bq-filter"));
}
