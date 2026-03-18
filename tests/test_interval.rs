use cfdnalab::shared::interval::{IndexedInterval, Interval, Span};

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
            "interval end ({}) must be greater than start ({})",
            (i64::MAX as u64) + 2,
            (i64::MAX as u64) + 1
        )
    );
    Ok(())
}
