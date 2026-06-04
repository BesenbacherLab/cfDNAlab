use cfdnalab::interval::{
    IndexedInterval, Interval, Span, TouchingMergePolicy, merge_intervals, merge_sorted_intervals,
    push_merged_interval,
};

#[test]
fn creates_valid_half_open_interval_and_reports_bounds() -> anyhow::Result<()> {
    // Arrange
    let expected_start = 10_u32;
    let expected_end = 25_u32;

    // Act
    let interval = Interval::new(expected_start, expected_end)?;

    // Assert
    assert_eq!(interval.start(), expected_start);
    assert_eq!(interval.end(), expected_end);
    assert_eq!(interval.as_tuple(), (expected_start, expected_end));
    assert_eq!(interval.into_inner(), (expected_start, expected_end));
    Ok(())
}

#[test]
fn creates_valid_half_open_span_and_reports_bounds() -> anyhow::Result<()> {
    // Arrange: spans may be empty, so start == end should be accepted here.
    let expected_start = 25_i64;
    let expected_end = 25_i64;

    // Act
    let span = Span::new(expected_start, expected_end)?;

    // Assert
    assert_eq!(span.start(), expected_start);
    assert_eq!(span.end(), expected_end);
    assert_eq!(span.as_tuple(), (expected_start, expected_end));
    assert_eq!(span.into_inner(), (expected_start, expected_end));
    assert!(span.is_empty());
    Ok(())
}

#[test]
fn rejects_inverted_span() {
    // Arrange
    let start = 30_u32;
    let end = 20_u32;

    // Act
    let error = Span::new(start, end).expect_err("expected inverted span to fail");

    // Assert
    assert_eq!(
        error.to_string(),
        "span end (20) must be greater than or equal to start (30)"
    );
}

#[test]
fn converts_tuple_bounds_with_try_from() -> anyhow::Result<()> {
    // Arrange
    let bounds = (100_u32, 125_u32);

    // Act
    let interval = Interval::try_from(bounds)?;

    // Assert
    assert_eq!(interval.start(), bounds.0);
    assert_eq!(interval.end(), bounds.1);
    Ok(())
}

#[test]
fn reports_interval_length_from_bounds() -> anyhow::Result<()> {
    // Arrange
    let start = 100_u32;
    let end = 125_u32;

    // Act
    let interval = Interval::new(start, end)?;

    // Assert
    assert_eq!(interval.len(), 25);
    Ok(())
}

#[test]
fn reports_indexed_interval_length_and_index() -> anyhow::Result<()> {
    // Arrange
    let start = 200_u64;
    let end = 260_u64;
    let original_index = 17_u64;

    // Act
    let indexed_interval = IndexedInterval::new(start, end, original_index)?;

    // Assert
    assert_eq!(indexed_interval.len(), 60);
    assert_eq!(indexed_interval.idx(), original_index);
    assert_eq!(indexed_interval.as_tuple(), (start, end, original_index));
    assert_eq!(indexed_interval.into_tuple(), (start, end, original_index));
    Ok(())
}

#[test]
fn rejects_empty_interval_when_start_equals_end() {
    // Arrange
    let coordinate = 42_u32;

    // Act
    let error = Interval::new(coordinate, coordinate).expect_err("expected empty interval to fail");

    // Assert
    assert_eq!(
        error.to_string(),
        "interval end (42) must be greater than start (42)"
    );
}

#[test]
fn rejects_interval_when_end_is_before_start() {
    // Arrange
    let start = 25_u32;
    let end = 10_u32;

    // Act
    let error = Interval::new(start, end).expect_err("expected inverted interval to fail");

    // Assert
    assert_eq!(
        error.to_string(),
        "interval end (10) must be greater than start (25)"
    );
}

#[test]
fn reports_set_relationships_for_intervals() -> anyhow::Result<()> {
    // Arrange:
    // - outer = [10,30), inner = [12,18), touching = [30,35), overlapping = [20,40).
    // - Half-open semantics mean 30 is excluded from outer, so touching must not intersect.
    // - The overlap with [20,40) is exactly [20,30).
    let outer = Interval::new(10_u32, 30_u32)?;
    let inner = Interval::new(12_u32, 18_u32)?;
    let touching = Interval::new(30_u32, 35_u32)?;
    let overlapping = Interval::new(20_u32, 40_u32)?;

    assert!(outer.contains_point(10));
    assert!(outer.contains_point(29));
    assert!(!outer.contains_point(30));
    assert!(outer.contains_interval(inner));
    assert!(!outer.contains_interval(overlapping));
    assert!(outer.intersects(overlapping));
    assert!(!outer.intersects(touching));
    assert_eq!(
        outer.intersection(overlapping),
        Some(Interval::new(20, 30)?)
    );
    assert_eq!(outer.intersection(touching), None);
    assert_eq!(outer.clip_to(overlapping), Some(Interval::new(20, 30)?));
    assert_eq!(outer.clip_to(touching), None);
    assert_eq!(outer.clip_lower(20), Some(Interval::new(20, 30)?));
    assert_eq!(outer.clip_lower(30), None);
    assert_eq!(outer.clip_upper(20), Some(Interval::new(10, 20)?));
    assert_eq!(outer.clip_upper(10), None);
    Ok(())
}

#[test]
fn expands_interval_to_cover_other_interval() -> anyhow::Result<()> {
    // Arrange
    let left = Interval::new(12_u32, 20_u32)?;
    let right = Interval::new(5_u32, 18_u32)?;

    // Act
    let expanded = left.expand_to_include(right);

    // Assert: the combined covered span is from the earlier start 5 to the later end 20.
    assert_eq!(expanded, Interval::new(5, 20)?);
    Ok(())
}

#[test]
fn expands_interval_on_both_sides() -> anyhow::Result<()> {
    let interval = Interval::new(12_u64, 20_u64)?;

    let expanded = interval.expand(3)?;

    assert_eq!(expanded, Interval::new(9_u64, 23_u64)?);
    Ok(())
}

#[test]
fn contracts_interval_on_both_sides() -> anyhow::Result<()> {
    let interval = Interval::new(12_u64, 20_u64)?;

    let contracted = interval.contract(3);

    assert_eq!(contracted, Some(Interval::new(15_u64, 17_u64)?));
    Ok(())
}

#[test]
fn returns_none_when_contracted_interval_would_be_empty() -> anyhow::Result<()> {
    let interval = Interval::new(12_u64, 20_u64)?;

    let contracted = interval.contract(4);

    assert_eq!(contracted, None);
    Ok(())
}

#[test]
fn shifts_unsigned_interval_left() -> anyhow::Result<()> {
    // Arrange: shifting [20,40) left by 5 should preserve the length 20 and move both bounds
    // together to [15,35).
    let interval = Interval::new(20_u64, 40_u64)?;

    // Act
    let shifted = interval.shift_left(5)?;

    // Assert
    assert_eq!(shifted, Interval::new(15_u64, 35_u64)?);
    Ok(())
}

#[test]
fn shifts_unsigned_interval_right() -> anyhow::Result<()> {
    // Arrange: shifting [20,40) right by 5 should preserve the length 20 and move both bounds
    // together to [25,45).
    let interval = Interval::new(20_u64, 40_u64)?;

    // Act
    let shifted = interval.shift_right(5)?;

    // Assert
    assert_eq!(shifted, Interval::new(25_u64, 45_u64)?);
    assert_eq!(shifted.len(), interval.len());
    Ok(())
}

#[test]
fn rejects_unsigned_left_shift_that_would_underflow() -> anyhow::Result<()> {
    // Arrange: shifting [3,8) left by 5 would try to move the start below zero, so the checked
    // shift must fail instead of wrapping the coordinates.
    let interval = Interval::new(3_u64, 8_u64)?;

    // Act
    let error = interval
        .shift_left(5)
        .expect_err("expected underflowing interval shift to fail");

    // Assert
    assert_eq!(
        error.to_string(),
        "offsetting interval [3, 8) by 5 would go out of bounds"
    );
    Ok(())
}

#[test]
fn offsets_signed_interval_by_signed_delta() -> anyhow::Result<()> {
    // Arrange: signed intervals keep the generic `offset(...)` API, so shifting [20,40) by -5
    // should move it to [15,35).
    let interval = Interval::new(20_i64, 40_i64)?;

    // Act
    let shifted = interval.offset(-5)?;

    // Assert
    assert_eq!(shifted, Interval::new(15_i64, 35_i64)?);
    Ok(())
}

#[test]
fn converts_tuple_slice_into_checked_intervals() -> anyhow::Result<()> {
    // Arrange: all tuples are already valid half-open intervals, so conversion should preserve the
    // original order and bounds exactly.
    let bounds = [(5_u64, 10_u64), (10_u64, 20_u64), (25_u64, 40_u64)];

    // Act
    let intervals = Interval::from_tuples(&bounds)?;

    // Assert
    assert_eq!(intervals.len(), 3);
    assert_eq!(intervals[0].as_tuple(), bounds[0]);
    assert_eq!(intervals[1].as_tuple(), bounds[1]);
    assert_eq!(intervals[2].as_tuple(), bounds[2]);
    Ok(())
}

#[test]
fn rejects_invalid_tuple_slice_when_building_intervals() {
    // Arrange
    let bounds = [(5_u64, 10_u64), (20_u64, 20_u64)];

    // Act
    let error = Interval::from_tuples(&bounds).expect_err("expected invalid interval tuple");

    // Assert
    assert_eq!(
        error.to_string(),
        "interval end (20) must be greater than start (20)"
    );
}

#[test]
fn converts_unsigned_interval_to_signed_interval() -> anyhow::Result<()> {
    // Arrange: both bounds fit in i64, so only the numeric type should change.
    let interval = Interval::new(25_u64, 40_u64)?;

    // Act
    let signed_interval = interval.try_to_i64()?;

    // Assert
    assert_eq!(signed_interval, Interval::new(25_i64, 40_i64)?);
    Ok(())
}

#[test]
fn converts_u32_interval_to_u64_interval() -> anyhow::Result<()> {
    // Arrange: widening the coordinate type should preserve both bounds exactly.
    let interval = Interval::new(12_u32, 34_u32)?;

    // Act
    let widened_interval = interval.try_to_u64()?;

    // Assert
    assert_eq!(widened_interval, Interval::new(12_u64, 34_u64)?);
    Ok(())
}

#[test]
fn converts_u64_interval_to_usize_interval() -> anyhow::Result<()> {
    // Arrange: both bounds fit in usize, so the checked conversion should preserve them exactly.
    let interval = Interval::new(25_u64, 40_u64)?;

    // Act
    let usize_interval = interval.try_to_usize()?;

    // Assert
    assert_eq!(usize_interval, Interval::new(25_usize, 40_usize)?);
    Ok(())
}

#[test]
fn converts_i64_interval_to_usize_interval() -> anyhow::Result<()> {
    // Arrange: a non-negative signed interval should convert cleanly to usize.
    let interval = Interval::new(7_i64, 19_i64)?;

    // Act
    let usize_interval = interval.try_to_usize()?;

    // Assert
    assert_eq!(usize_interval, Interval::new(7_usize, 19_usize)?);
    Ok(())
}

#[test]
fn rejects_signed_interval_that_does_not_fit_in_usize() -> anyhow::Result<()> {
    // Arrange: negative bounds are invalid for usize conversion.
    let interval = Interval::new(-2_i64, -1_i64)?;

    // Act
    let error = interval
        .try_to_usize()
        .expect_err("expected negative interval conversion to fail");

    // Assert
    assert_eq!(
        error.to_string(),
        "converting interval [-2, -1) to usize failed"
    );
    Ok(())
}

#[test]
fn rejects_unsigned_interval_that_does_not_fit_in_i64() -> anyhow::Result<()> {
    // Arrange: both bounds exceed i64::MAX, so the checked conversion must fail instead of
    // silently wrapping or truncating.
    let interval = Interval::new((i64::MAX as u64) + 1, (i64::MAX as u64) + 2)?;

    // Act
    let error = interval
        .try_to_i64()
        .expect_err("expected out-of-range interval conversion to fail");

    // Assert
    assert_eq!(
        error.to_string(),
        format!(
            "converting interval [{}, {}) to i64 failed",
            (i64::MAX as u64) + 1,
            (i64::MAX as u64) + 2
        )
    );
    Ok(())
}

#[test]
fn pushes_sorted_intervals_with_explicit_touching_policy() -> anyhow::Result<()> {
    // Arrange:
    // - Start with [10,20).
    // - [15,25) overlaps and should always merge into [10,25).
    // - [25,30) only touches that merged interval at one boundary.
    let first = Interval::new(10_u32, 20_u32)?;
    let overlapping = Interval::new(15_u32, 25_u32)?;
    let touching = Interval::new(25_u32, 30_u32)?;

    // Act
    let mut keep_touching_separate = vec![first];
    push_merged_interval(
        &mut keep_touching_separate,
        overlapping,
        TouchingMergePolicy::KeepTouchingSeparate,
    );
    push_merged_interval(
        &mut keep_touching_separate,
        touching,
        TouchingMergePolicy::KeepTouchingSeparate,
    );

    let mut merge_touching = vec![first];
    push_merged_interval(
        &mut merge_touching,
        overlapping,
        TouchingMergePolicy::MergeTouching,
    );
    push_merged_interval(
        &mut merge_touching,
        touching,
        TouchingMergePolicy::MergeTouching,
    );

    // Assert:
    // - overlap merges in both cases, so both paths first become [10,25)
    // - the touching interval stays separate in the first path
    // - the touching interval collapses into [10,30) in the second path
    assert_eq!(
        keep_touching_separate,
        vec![
            Interval::new(10_u32, 25_u32)?,
            Interval::new(25_u32, 30_u32)?
        ]
    );
    assert_eq!(merge_touching, vec![Interval::new(10_u32, 30_u32)?]);
    Ok(())
}

#[test]
fn merges_sorted_intervals_without_resorting() -> anyhow::Result<()> {
    // Arrange:
    // - The input is already start-sorted.
    // - [5,10) and [10,12) touch, so they merge only under the touching-merge policy.
    // - [20,25) stays separate either way.
    let sorted_intervals = vec![
        Interval::new(5_u32, 10_u32)?,
        Interval::new(10_u32, 12_u32)?,
        Interval::new(20_u32, 25_u32)?,
    ];

    // Act
    let keep_touching_separate = merge_sorted_intervals(
        sorted_intervals.clone(),
        TouchingMergePolicy::KeepTouchingSeparate,
    );
    let merge_touching =
        merge_sorted_intervals(sorted_intervals, TouchingMergePolicy::MergeTouching);

    // Assert
    assert_eq!(
        keep_touching_separate,
        vec![
            Interval::new(5_u32, 10_u32)?,
            Interval::new(10_u32, 12_u32)?,
            Interval::new(20_u32, 25_u32)?,
        ]
    );
    assert_eq!(
        merge_touching,
        vec![
            Interval::new(5_u32, 12_u32)?,
            Interval::new(20_u32, 25_u32)?
        ]
    );
    Ok(())
}

#[test]
fn sorts_then_merges_unsorted_intervals() -> anyhow::Result<()> {
    // Arrange:
    // - The intervals arrive unsorted as [20,30), [5,10), [8,12).
    // - Sorting by (start, end) gives [5,10), [8,12), [20,30).
    // - The first two overlap and therefore merge into [5,12).
    let unsorted_intervals = vec![
        Interval::new(20_u32, 30_u32)?,
        Interval::new(5_u32, 10_u32)?,
        Interval::new(8_u32, 12_u32)?,
    ];

    // Act
    let merged = merge_intervals(
        unsorted_intervals,
        TouchingMergePolicy::KeepTouchingSeparate,
    );

    // Assert
    assert_eq!(
        merged,
        vec![
            Interval::new(5_u32, 12_u32)?,
            Interval::new(20_u32, 30_u32)?
        ]
    );
    Ok(())
}
