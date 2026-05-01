use super::{
    SummaryStatsRow, coverage_cv_overflow_label, derive_coefficient_of_variation_coverage,
    derive_nonnegative_variance_coverage, derive_summary_stats, write_bedgraph_runs,
    write_summary_stats_row, write_windowed_runs,
};
use crate::shared::base::{ZEROISH_F32_TOLERANCE, ZEROISH_F64_TOLERANCE};
use crate::shared::interval::Interval;
use std::io::{Error, ErrorKind, Result as IoResult, Write};

struct FailingWriter;

impl Write for FailingWriter {
    fn write(&mut self, _buf: &[u8]) -> IoResult<usize> {
        Err(Error::new(
            ErrorKind::BrokenPipe,
            "synthetic positional write failure",
        ))
    }

    fn flush(&mut self) -> IoResult<()> {
        Ok(())
    }
}

#[test]
fn derive_summary_stats_returns_nan_fields_when_no_positions_are_eligible() {
    // Arrange
    // With zero eligible positions, mean, variance, SD, CV, and covered fraction are all
    // undefined. This is not a numeric failure, it is just missing support for the statistic
    let stats = derive_summary_stats(20, 20, 0, 0, 0.0, 0.0).expect("zero-eligible row should succeed");

    // Assert
    assert_eq!(stats.span_positions, 20);
    assert_eq!(stats.blacklisted_positions, 20);
    assert_eq!(stats.eligible_positions, 0);
    assert_eq!(stats.nonzero_positions, 0);
    assert_eq!(stats.coverage_sum, 0.0);
    assert_eq!(stats.coverage_sum_of_squares, 0.0);
    assert_eq!(stats.total_coverage, 0.0);
    assert!(stats.average_coverage.is_nan());
    assert!(stats.variance_coverage.is_nan());
    assert!(stats.sd_coverage.is_nan());
    assert!(stats.coefficient_of_variation_coverage.is_nan());
    assert!(stats.covered_fraction.is_nan());
}

#[test]
fn derive_nonnegative_variance_coverage_snaps_tiny_negative_cancellation_to_zero() {
    // Arrange
    // Use a mean of 1.0 and choose `coverage_sum_of_squares` so the raw variance becomes
    // `-ZEROISH_F32_TOLERANCE / 2`. That value is mathematically impossible, but small enough
    // that we intentionally classify it as a cancellation residue from `E[x^2] - E[x]^2`
    let average_coverage = 1.0_f64;
    let tiny_negative_variance = -(ZEROISH_F32_TOLERANCE as f64) / 2.0;
    let coverage_sum_of_squares = 1.0 + tiny_negative_variance;

    // Act
    let variance = derive_nonnegative_variance_coverage(1, coverage_sum_of_squares, average_coverage)
        .expect("tiny negative variance should be repaired");

    // Assert
    assert_eq!(variance, 0.0);
}

#[test]
fn derive_nonnegative_variance_coverage_errors_on_material_negative_values() {
    // Arrange
    // Push the raw variance well beyond the allowed tolerance so this is treated as a real
    // invariant violation rather than recoverable floating-point residue
    let average_coverage = 1.0_f64;
    let material_negative_variance = -10.0 * ZEROISH_F32_TOLERANCE as f64;
    let coverage_sum_of_squares = 1.0 + material_negative_variance;

    // Act
    let err = derive_nonnegative_variance_coverage(1, coverage_sum_of_squares, average_coverage)
        .expect_err("materially negative variance should fail");

    // Assert
    let err_text = err.to_string();
    assert!(err_text.contains("negative variance_coverage"));
    assert!(err_text.contains("eligible_positions=1"));
}

#[test]
fn derive_summary_stats_keeps_tiny_positive_variance_and_finite_sd() {
    // Arrange
    // This aggregate is consistent with two coverage values symmetric around 1.0:
    //   x1 = 1 + a
    //   x2 = 1 - a
    // Their mean is exactly 1.0 and their variance is exactly `a^2`.
    // Choosing `a^2 = 1e-9` gives a genuinely tiny positive variance, which should remain
    // positive and produce a finite SD instead of being forced to zero
    let expected_variance = 1.0e-9_f64;
    let variance_tolerance = 2.0 * f64::EPSILON;
    let sd_squared_tolerance = 16.0 * f64::EPSILON * expected_variance;
    let stats = derive_summary_stats(2, 0, 2, 2, 2.0, 2.0 + 2.0 * expected_variance)
        .expect("tiny positive variance should remain valid");

    // Assert
    // The input aggregate is hand-derived, but the literal `1e-9` and the derived `2.0 + 2.0 * a^2`
    // are still represented as `f64`. That means we should assert the intended value with a tight
    // tolerance instead of requiring exact decimal equality from binary floating-point
    //
    // For SD, the strongest behavior-level check is not "does it equal the idealized decimal
    // sqrt(1e-9) bit-for-bit", but "is it the square root of the variance that this function
    // actually produced"
    assert_eq!(stats.average_coverage, 1.0);
    assert!(
        (stats.variance_coverage - expected_variance).abs() <= variance_tolerance,
        "expected variance within {variance_tolerance}, got {} vs {}",
        stats.variance_coverage,
        expected_variance
    );
    let variance_recovered_from_sd = stats.sd_coverage * stats.sd_coverage;
    assert!(
        (variance_recovered_from_sd - stats.variance_coverage).abs() <= sd_squared_tolerance,
        "expected SD^2 within {sd_squared_tolerance} of the derived variance, got {} vs {}",
        variance_recovered_from_sd,
        stats.variance_coverage
    );
    assert!(stats.sd_coverage.is_finite());
}

#[test]
fn derive_nonnegative_variance_coverage_keeps_legitimate_positive_values() {
    // Arrange
    // Choose an exact positive variance:
    //   eligible_positions = 2
    //   mean = 1
    //   E[x^2] = 1.25
    // so variance = 1.25 - 1^2 = 0.25
    let variance = derive_nonnegative_variance_coverage(2, 2.5, 1.0)
        .expect("positive variance should pass through unchanged");

    // Assert
    assert_eq!(variance, 0.25);
}

#[test]
fn derive_coefficient_of_variation_coverage_returns_nan_for_zeroish_mean() {
    // Arrange
    // A mean inside the shared `f64` zeroish tolerance is not a meaningful coverage signal.
    // Treating it as zero avoids manufacturing a gigantic ratio from denominator noise
    let zeroish_mean = ZEROISH_F64_TOLERANCE / 2.0;
    let finite_sd = 1.0_f64;

    // Act
    let coefficient_of_variation =
        derive_coefficient_of_variation_coverage(zeroish_mean, finite_sd);

    // Assert
    assert!(coefficient_of_variation.is_nan());
}

#[test]
fn write_summary_stats_row_marks_extreme_finite_cv_as_greater_than_one_e6() {
    // Arrange
    // Keep every other field ordinary so the test isolates the CV presentation rule.
    // The exact CV value does not matter once it crosses the readability threshold.
    // This row is intentionally synthetic and only meant to exercise formatting. Its derived
    // fields do not need to be mathematically consistent with the raw fields here
    let stats = SummaryStatsRow {
        span_positions: 10,
        blacklisted_positions: 0,
        eligible_positions: 10,
        nonzero_positions: 1,
        coverage_sum: 1.0,
        coverage_sum_of_squares: 2.0,
        average_coverage: 0.1,
        total_coverage: 1.0,
        variance_coverage: 0.01,
        sd_coverage: 0.1,
        coefficient_of_variation_coverage: 1.0e6 + 1.0,
        covered_fraction: 0.1,
    };
    let mut out = Vec::new();

    // Act
    write_summary_stats_row(
        &mut out,
        "chr1",
        Interval::new(0_u64, 10_u64).expect("test interval should be valid"),
        stats,
        3,
    )
    .expect("summary row should format");
    let rendered = String::from_utf8(out).expect("summary row should stay valid UTF-8");

    // Assert
    assert!(
        rendered.contains("\t>1e6\t"),
        "expected CV field to render as >1e6, got: {rendered}"
    );
}

#[test]
fn coverage_cv_overflow_label_is_derived_from_the_threshold_constant() {
    // Arrange / Act / Assert
    // Keep this explicit so a future threshold change has to update the expected label too
    assert_eq!(coverage_cv_overflow_label(), ">1e6");
}

#[test]
fn write_bedgraph_runs_returns_immediate_writer_errors() {
    // Arrange: a single non-zero run forces one bedGraph row write. The writer fails immediately
    // on `write`, while `flush` would succeed, so this catches swallowed callback errors.
    let coverage = [1.0_f32, 1.0_f32];
    let mut writer = FailingWriter;

    // Act
    let err = write_bedgraph_runs("chr1", &coverage, None, 0, coverage.len(), 100, 3, false, &mut writer)
        .expect_err("bedGraph writer should return the row write error");

    // Assert
    assert!(
        err.to_string()
            .contains("synthetic positional write failure"),
        "expected immediate write error, got {err:#}"
    );
}

#[test]
fn write_windowed_runs_returns_immediate_writer_errors() {
    // Arrange: the indexed windowed writer path has a distinct row format from bedGraph, so it
    // should independently propagate the same immediate `Write` failure.
    let coverage = [2.0_f32, 2.0_f32];
    let mut writer = FailingWriter;

    // Act
    let err = write_windowed_runs(
        "chr1",
        &coverage,
        None,
        0,
        coverage.len(),
        200,
        Some(7),
        3,
        false,
        &mut writer,
    )
    .expect_err("windowed writer should return the row write error");

    // Assert
    assert!(
        err.to_string()
            .contains("synthetic positional write failure"),
        "expected immediate write error, got {err:#}"
    );
}
