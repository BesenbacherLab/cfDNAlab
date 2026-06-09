use super::{
    SummaryStatsRow, coverage_cv_overflow_label, derive_coefficient_of_variation_coverage,
    derive_nonnegative_variance_coverage, derive_summary_stats, write_bedgraph_runs,
    write_bed_aggregate_output, write_grouped_bed_aggregate_output,
    write_summary_stats_row, write_windowed_runs,
};
use crate::commands::fcoverage::{
    config::COVERAGE_SIGNAL_LABEL, reducer::TileAggregateTempFiles,
    window_results::CoverageWindowAction,
};
use crate::shared::base::{ZEROISH_F32_TOLERANCE, ZEROISH_F64_TOLERANCE};
use crate::shared::{
    bed::{GroupedCoverageLayout, Windows},
    interval::{IndexedInterval, Interval},
    io::open_text_reader,
};
use anyhow::Result;
use fxhash::FxHashMap;
use std::io::{Error, ErrorKind, Read, Result as IoResult, Write};
use std::path::Path;
use tempfile::TempDir;

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

fn write_text(path: &Path, text: &str) -> Result<()> {
    std::fs::write(path, text)?;
    Ok(())
}

fn read_text(path: &Path) -> Result<String> {
    let mut reader = open_text_reader(path)?;
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    Ok(text)
}

fn grouped_layout_for_writer_tests() -> Result<GroupedCoverageLayout> {
    let mut segments_by_chr = FxHashMap::default();
    segments_by_chr.insert(
        "chr1".to_string(),
        Windows::new(vec![
            IndexedInterval::new(0_u64, 10_u64, 0_u64)?,
            IndexedInterval::new(20_u64, 30_u64, 1_u64)?,
        ]),
    );

    let mut segment_idx_to_group_idx = FxHashMap::default();
    segment_idx_to_group_idx.insert(0, 5);
    segment_idx_to_group_idx.insert(1, 5);

    let mut group_span_positions = FxHashMap::default();
    group_span_positions.insert(5, 20);

    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(5, "alpha".to_string());

    Ok(GroupedCoverageLayout {
        segments_by_chr,
        segment_idx_to_group_idx,
        group_span_positions,
        group_idx_to_name,
    })
}

#[test]
fn grouped_writer_folds_explicit_segment_tile_outputs_into_group_rows() -> Result<()> {
    // Arrange
    // Grouped output has two stages: explicit tile segment rows are reduced, then segment rows are
    // folded into stable group rows. This test pins the second stage at the final writer boundary.
    //
    // Group 5 has two segments with total represented span 20. The tile rows contribute:
    //   segment 0: coverage_sum 4, eligible 10
    //   segment 1: coverage_sum 6, eligible 10
    // So group 5 total_coverage = 10 and eligible_positions = 20.
    let temp_dir = TempDir::new()?;
    let out_dir = TempDir::new()?;
    let partials_path = temp_dir.path().join("group_segment_rows");
    write_text(&partials_path, "0\t4\t10\t0\n1\t6\t10\t0\n")?;
    write_text(
        &temp_dir.path().join("run.part.chrom-000000.9.tsv"),
        "0\t999\t999\t0\n",
    )?;

    let mut tile_outputs_by_chr = FxHashMap::default();
    tile_outputs_by_chr.insert(
        "chr1".to_string(),
        vec![TileAggregateTempFiles {
            tile_index: 0,
            partials_path,
            cross_index_path: None,
        }],
    );
    let grouped_layout = grouped_layout_for_writer_tests()?;
    let final_path = out_dir.path().join("grouped_total.tsv.zst");

    // Act
    write_grouped_bed_aggregate_output(
        &final_path,
        &tile_outputs_by_chr,
        &grouped_layout,
        &["chr1".to_string()],
        CoverageWindowAction::Total,
        0,
        1,
        COVERAGE_SIGNAL_LABEL,
        None,
    )?;
    let text = read_text(&final_path)?;

    // Assert
    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\ttotal_coverage",
            "5\t20\t0\t20\t10",
        ]
    );

    Ok(())
}

#[test]
fn grouped_summary_writer_folds_explicit_raw_moments_before_deriving_stats() -> Result<()> {
    // Arrange
    // Summary grouped output must carry raw moments through segment reduction and group folding
    // before deriving final statistics. The two segment rows represent coverage values:
    //   segment 0: ten eligible positions, coverage_sum=4, sum_of_squares=4
    //   segment 1: ten eligible positions, coverage_sum=6, sum_of_squares=10
    // Group totals:
    //   span_positions = eligible_positions = 20
    //   nonzero_positions = 8
    //   covered_fraction = 8 / 20 = 0.4
    //   total_coverage = 10
    //   total_squared_coverage = 14
    // Derived:
    //   average = 10 / 20 = 0.5
    //   variance = 14 / 20 - 0.5^2 = 0.45
    let temp_dir = TempDir::new()?;
    let out_dir = TempDir::new()?;
    let partials_path = temp_dir.path().join("group_segment_summary_rows");
    write_text(&partials_path, "0\t4\t10\t0\t3\t4\n1\t6\t10\t0\t5\t10\n")?;

    let mut tile_outputs_by_chr = FxHashMap::default();
    tile_outputs_by_chr.insert(
        "chr1".to_string(),
        vec![TileAggregateTempFiles {
            tile_index: 0,
            partials_path,
            cross_index_path: None,
        }],
    );
    let grouped_layout = grouped_layout_for_writer_tests()?;
    let final_path = out_dir.path().join("grouped_summary.tsv.zst");

    // Act
    write_grouped_bed_aggregate_output(
        &final_path,
        &tile_outputs_by_chr,
        &grouped_layout,
        &["chr1".to_string()],
        CoverageWindowAction::SummaryStats,
        6,
        1,
        COVERAGE_SIGNAL_LABEL,
        None,
    )?;
    let text = read_text(&final_path)?;
    let rows: Vec<Vec<_>> = text
        .lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect();

    // Assert
    assert_eq!(
        rows[0],
        vec![
            "group_idx",
            "span_positions",
            "blacklisted_positions",
            "eligible_positions",
            "nonzero_positions",
            "covered_fraction",
            "total_coverage",
            "total_squared_coverage",
            "average_coverage",
            "variance_coverage",
            "sd_coverage",
            "coefficient_of_variation_coverage",
        ]
    );
    assert_eq!(rows[1][0..5], ["5", "20", "0", "20", "8"]);
    assert_eq!(rows[1][5], "0.4");
    assert_eq!(rows[1][6], "10");
    assert_eq!(rows[1][7], "14");
    assert_eq!(rows[1][8], "0.5");
    assert_eq!(rows[1][9], "0.45");

    Ok(())
}

#[test]
fn bed_writer_errors_when_windows_have_no_returned_tile_outputs() -> Result<()> {
    // Arrange
    // With explicit returned paths, a chromosome that has BED windows should also have aggregate
    // tile outputs. Missing outputs are an upstream contract error, not a valid zero-coverage row.
    let out_dir = TempDir::new()?;
    let mut windows_by_chr = FxHashMap::default();
    windows_by_chr.insert(
        "chr1".to_string(),
        Windows::new(vec![IndexedInterval::new(0_u64, 10_u64, 0_u64)?]),
    );
    let tile_outputs_by_chr = FxHashMap::default();

    // Act
    let err = write_bed_aggregate_output(
        &out_dir.path().join("missing_outputs.tsv.zst"),
        &tile_outputs_by_chr,
        &windows_by_chr,
        &["chr1".to_string()],
        true,
        CoverageWindowAction::Total,
        0,
        1,
        COVERAGE_SIGNAL_LABEL,
        None,
    )
    .expect_err("BED windows without returned aggregate tile outputs should fail");

    // Assert
    let message = err.to_string();
    assert!(message.contains("No returned BED aggregate tile outputs"));
    assert!(message.contains("chr1"));

    Ok(())
}

#[test]
fn grouped_writer_errors_when_segments_have_no_returned_tile_outputs() -> Result<()> {
    // Arrange
    // Group span is seeded from the layout before reduction. If segment outputs are missing, the
    // writer must fail rather than write the seeded span with zero coverage-derived fields.
    let out_dir = TempDir::new()?;
    let tile_outputs_by_chr = FxHashMap::default();
    let grouped_layout = grouped_layout_for_writer_tests()?;

    // Act
    let err = write_grouped_bed_aggregate_output(
        &out_dir.path().join("missing_grouped_outputs.tsv.zst"),
        &tile_outputs_by_chr,
        &grouped_layout,
        &["chr1".to_string()],
        CoverageWindowAction::Total,
        0,
        1,
        COVERAGE_SIGNAL_LABEL,
        None,
    )
    .expect_err("grouped segments without returned aggregate tile outputs should fail");

    // Assert
    let message = err.to_string();
    assert!(message.contains("No returned grouped BED aggregate tile outputs"));
    assert!(message.contains("chr1"));

    Ok(())
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
    assert!(stats.covered_fraction.is_nan());
    assert_eq!(stats.total_coverage, 0.0);
    assert_eq!(stats.total_squared_coverage, 0.0);
    assert!(stats.average_coverage.is_nan());
    assert!(stats.variance_coverage.is_nan());
    assert!(stats.sd_coverage.is_nan());
    assert!(stats.coefficient_of_variation_coverage.is_nan());
}

#[test]
fn derive_nonnegative_variance_coverage_snaps_tiny_negative_cancellation_to_zero() {
    // Arrange
    // Use a mean of 1.0 and choose `total_squared_coverage` so the raw variance becomes
    // `-ZEROISH_F32_TOLERANCE / 2`. That value is mathematically impossible, but small enough
    // that we intentionally classify it as a cancellation residue from `E[x^2] - E[x]^2`
    let average_coverage = 1.0_f64;
    let tiny_negative_variance = -(ZEROISH_F32_TOLERANCE as f64) / 2.0;
    let total_squared_coverage = 1.0 + tiny_negative_variance;

    // Act
    let variance = derive_nonnegative_variance_coverage(1, total_squared_coverage, average_coverage)
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
    let total_squared_coverage = 1.0 + material_negative_variance;

    // Act
    let err = derive_nonnegative_variance_coverage(1, total_squared_coverage, average_coverage)
        .expect_err("materially negative variance should fail");

    // Assert
    let err_text = err.to_string();
    assert!(err_text.contains("negative variance_coverage"));
    assert!(err_text.contains("eligible_positions=1"));
    assert!(err_text.contains("total_squared_coverage="));
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
        covered_fraction: 0.1,
        total_coverage: 1.0,
        total_squared_coverage: 2.0,
        average_coverage: 0.1,
        variance_coverage: 0.01,
        sd_coverage: 0.1,
        coefficient_of_variation_coverage: 1.0e6 + 1.0,
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
fn write_bedgraph_runs_encodes_whole_positional_coverage_as_nonzero_runs() -> Result<()> {
    // Arrange
    // Tile-local coverage:
    //   positions 100..102 => 0, omitted because keep_zero_runs=false
    //   positions 102..104 => 1, one run
    //   position  104..105 => 0, omitted
    //   positions 105..108 => 2, one run
    let coverage = [0.0_f32, 0.0, 1.0, 1.0, 0.0, 2.0, 2.0, 2.0];
    let mut output = Vec::new();

    // Act
    write_bedgraph_runs(
        "chr1",
        &coverage,
        None,
        0,
        coverage.len(),
        100,
        3,
        false,
        &mut output,
    )?;

    // Assert
    assert_eq!(
        String::from_utf8(output)?,
        "chr1\t102\t104\t1\nchr1\t105\t108\t2\n"
    );
    Ok(())
}

#[test]
fn write_windowed_runs_skips_masked_positions_and_preserves_window_index() -> Result<()> {
    // Arrange
    // Window [8,13) has coverage at every base, but positions 9,10,11 are masked.
    // Current positional output omits masked bases, so the unmasked bases become two runs:
    // [8,9) and [12,13), both carrying the original window index.
    let coverage = [1.0_f32; 13];
    let mut mask = [0_u8; 13];
    mask[9] = 1;
    mask[10] = 1;
    mask[11] = 1;
    let mut output = Vec::new();

    // Act
    write_windowed_runs(
        "chr1",
        &coverage,
        Some(&mask),
        8,
        13,
        0,
        Some(7),
        3,
        false,
        &mut output,
    )?;

    // Assert
    assert_eq!(
        String::from_utf8(output)?,
        "chr1\t8\t9\t1\t7\nchr1\t12\t13\t1\t7\n"
    );
    Ok(())
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
