use cfdnalab::shared::interval::{IndexedInterval, Interval};

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
    assert_eq!(interval.into_inner(), (expected_start, expected_end));
    Ok(())
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
