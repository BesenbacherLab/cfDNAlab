use super::*;
use anyhow::Result;
use crate::commands::fcoverage::config::LengthNormalizationMode;
use std::fs::File;
use std::io::Write;
use tempfile::TempDir;

fn make_run_result(path: PathBuf) -> FCoverageRunResult {
    FCoverageRunResult {
        counters: crate::commands::counters::FCoverageCounters::default(),
        mean_normalization_length: None,
        final_out_path: path,
    }
}

fn write_plain_tsv(tempdir: &TempDir, file_name: &str, lines: &[&str]) -> Result<PathBuf> {
    let path = tempdir.path().join(file_name);
    let mut file = File::create(&path)?;
    for line in lines {
        writeln!(file, "{line}")?;
    }
    Ok(path)
}

fn write_zstd_tsv(tempdir: &TempDir, file_name: &str, lines: &[&str]) -> Result<PathBuf> {
    let path = tempdir.path().join(file_name);
    let file = File::create(&path)?;
    let mut encoder = zstd::Encoder::new(file, 3)?;
    for line in lines {
        writeln!(encoder, "{line}")?;
    }
    encoder.finish()?;
    Ok(path)
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_reads_contiguous_bins_from_zstd() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_zstd_tsv(
        &tempdir,
        "coverage.average.tsv.zst",
        &[
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.25\t0",
            "chr1\t10\t20\t2.5\t3",
            "chr2\t0\t5\t0\t0",
        ],
    )?;

    // Act
    let bins_by_chr = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string(), "chr2".to_string()],
        "average_coverage",
    )?;

    // Assert
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert_eq!(chr1_bins.len(), 2);
    assert_eq!(chr1_bins[0].start(), 0);
    assert_eq!(chr1_bins[0].end(), 10);
    assert_eq!(chr1_bins[0].average_coverage, 1.25);
    assert_eq!(chr1_bins[1].start(), 10);
    assert_eq!(chr1_bins[1].end(), 20);
    assert_eq!(chr1_bins[1].average_coverage, 2.5);
    assert_eq!(chr1_bins[1].average_overlap_coverage, 0.0);
    assert_eq!(chr1_bins[1].scaling_factor, 0.0);

    let chr2_bins = bins_by_chr.get("chr2").expect("chr2 bins should exist");
    assert_eq!(chr2_bins.len(), 1);
    assert_eq!(chr2_bins[0].start(), 0);
    assert_eq!(chr2_bins[0].end(), 5);
    assert_eq!(chr2_bins[0].average_coverage, 0.0);

    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_empty_file() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(&tempdir, "coverage.average.tsv", &[])?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "average_coverage",
    )
    .expect_err("empty file should fail");

    // Assert
    assert!(err.to_string().contains("empty file; header required"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_unexpected_header() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.average.tsv",
        &[
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "average_coverage",
    )
    .expect_err("wrong header should fail");

    // Assert
    assert!(err.to_string().contains("unexpected fcoverage header"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_reads_requested_value_column() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.total.tsv",
        &[
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t10\t12.5\t0",
            "chr1\t10\t20\t25\t0",
        ],
    )?;

    // Act
    let bins_by_chr = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "total_coverage",
    )?;

    // Assert
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert_eq!(chr1_bins.len(), 2);
    assert_eq!(chr1_bins[0].average_coverage, 12.5);
    assert_eq!(chr1_bins[1].average_coverage, 25.0);
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_short_rows() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.average.tsv",
        &[
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.25",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "average_coverage",
    )
    .expect_err("short row should fail");

    // Assert
    assert!(err.to_string().contains("expected 5 tab-separated columns"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_invalid_start() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.average.tsv",
        &[
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\tzero\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "average_coverage",
    )
    .expect_err("invalid start should fail");

    // Assert
    assert!(err.to_string().contains("invalid start"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_invalid_interval() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.average.tsv",
        &[
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t10\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "average_coverage",
    )
    .expect_err("zero-length interval should fail");

    // Assert
    assert!(
        err.to_string().contains("must be <")
            || err.to_string().contains("must be strictly less")
            || err.to_string().contains("start")
    );
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_missing_requested_chromosome() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.average.tsv",
        &[
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string(), "chr2".to_string()],
        "average_coverage",
    )
    .expect_err("missing requested chromosome should fail");

    // Assert
    assert!(err.to_string().contains("missing chromosome 'chr2'"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_first_bin_not_starting_at_zero() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.average.tsv",
        &[
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t5\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "average_coverage",
    )
    .expect_err("first bin starting after zero should fail");

    // Assert
    assert!(err.to_string().contains("did not start at 0"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_gap_between_bins() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.average.tsv",
        &[
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.0\t0",
            "chr1\t12\t20\t2.0\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "average_coverage",
    )
    .expect_err("gap between bins should fail");

    // Assert
    assert!(err.to_string().contains("non-contiguous stride bins"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_overlap_between_bins() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.average.tsv",
        &[
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.0\t0",
            "chr1\t8\t20\t2.0\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "average_coverage",
    )
    .expect_err("overlap between bins should fail");

    // Assert
    assert!(err.to_string().contains("non-contiguous stride bins"));
    Ok(())
}

#[test]
fn normalize_average_overlap_by_global_mean_ignores_bins_below_support_floor() -> Result<()> {
    // Arrange
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                average_coverage: 1.0,
                average_overlap_coverage: 0.5,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 20)?,
                average_coverage: 0.0,
                average_overlap_coverage: 5e-11,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_average_overlap_by_global_mean(&mut bins_by_chr, false, true)?;

    // Assert
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert!((mean - 0.5).abs() <= 1e-10, "expected mean 0.5, got {mean}");
    assert!(
        (chr1_bins[0].scaling_factor - 1.0).abs() <= 1e-6,
        "supported bin should normalize to 1.0, got {}",
        chr1_bins[0].scaling_factor
    );
    assert_eq!(
        chr1_bins[1].scaling_factor, 0.0,
        "below-floor support should be treated as zero"
    );

    Ok(())
}

#[test]
fn normalize_average_overlap_by_global_mean_explains_all_zero_smoothed_mass() -> Result<()> {
    // Arrange
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                average_coverage: 0.0,
                average_overlap_coverage: 0.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 20)?,
                average_coverage: 0.0,
                average_overlap_coverage: 5e-11,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let err = normalize_average_overlap_by_global_mean(&mut bins_by_chr, true, true)
        .expect_err("all-zero smoothed mass should fail");

    // Assert
    let message = err.to_string();
    assert!(
        message.contains("no finite non-zero smoothed fragment mass after filtering"),
        "unexpected error message: {message}"
    );
    assert!(
        message.contains("2 stride bins"),
        "unexpected error message: {message}"
    );
    assert!(
        message.contains("--chromosomes") && message.contains("--min-mapq"),
        "unexpected error message: {message}"
    );

    Ok(())
}

#[test]
fn normalize_average_overlap_by_global_mean_keeps_bins_above_support_floor() -> Result<()> {
    // Arrange
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                average_coverage: 1.0,
                average_overlap_coverage: 1.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 20)?,
                average_coverage: 0.0,
                average_overlap_coverage: 2e-9,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_average_overlap_by_global_mean(&mut bins_by_chr, false, true)?;

    // Assert
    let expected_mean = ((1.0 + 2e-9) / 2.0) as f32;
    let expected_small_scaling = expected_mean / 2e-9;
    let expected_large_scaling = expected_mean / 1.0;
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert!(
        (mean - expected_mean).abs() <= 1e-12_f32,
        "expected mean {expected_mean}, got {mean}"
    );
    assert!(
        (chr1_bins[0].scaling_factor - expected_large_scaling).abs() <= 1e-6_f32,
        "expected large-bin scaling close to {expected_large_scaling}, got {}",
        chr1_bins[0].scaling_factor
    );
    assert!(
        chr1_bins[1].scaling_factor > 0.0,
        "above-floor support should remain non-zero"
    );
    assert!(
        (chr1_bins[1].scaling_factor - expected_small_scaling).abs() <= 1.0_f32,
        "expected small-bin scaling close to {expected_small_scaling}, got {}",
        chr1_bins[1].scaling_factor
    );

    Ok(())
}

#[test]
fn build_fcoverage_stride_config_uses_unit_mass_and_total_for_fragment_count_weights() {
    let tempdir = TempDir::new().expect("tempdir should exist");
    let args = ScalingWeightsArgs::new(
        crate::commands::cli_common::IOCArgs {
            bam: PathBuf::from("input.bam"),
            output_dir: tempdir.path().to_path_buf(),
            n_threads: 1,
        },
        crate::commands::cli_common::ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    );

    let cfg = build_fcoverage_stride_config(
        &args,
        tempdir.path(),
        true,
        ScalingWeightsCommand::FragmentCount,
        false,
    );

    assert_eq!(
        cfg.normalize_by_length_mode,
        LengthNormalizationMode::UnitMass
    );
    assert_eq!(cfg.per_window, CoverageWindowAction::Total);
}
