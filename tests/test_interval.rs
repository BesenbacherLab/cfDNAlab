use cfdnalab::shared::interval::Interval;

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
