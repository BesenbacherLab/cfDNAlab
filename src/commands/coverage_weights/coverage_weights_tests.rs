use super::*;
use anyhow::Result;
use std::fs::File;
use std::io::Write;
use tempfile::TempDir;

fn make_run_result(path: PathBuf) -> FCoverageRunResult {
    FCoverageRunResult {
        counters: crate::commands::counters::FCoverageCounters::default(),
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
fn load_stride_bins_from_fcoverage_average_tsv_reads_contiguous_bins_from_zstd() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_zstd_tsv(
        &tempdir,
        "coverage.avg.tsv.zst",
        &[
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.25\t0",
            "chr1\t10\t20\t2.5\t3",
            "chr2\t0\t5\t0\t0",
        ],
    )?;

    // Act
    let bins_by_chr = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string(), "chr2".to_string()],
    )?;

    // Assert
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert_eq!(chr1_bins.len(), 2);
    assert_eq!(chr1_bins[0].start(), 0);
    assert_eq!(chr1_bins[0].end(), 10);
    assert_eq!(chr1_bins[0].avg_coverage, 1.25);
    assert_eq!(chr1_bins[1].start(), 10);
    assert_eq!(chr1_bins[1].end(), 20);
    assert_eq!(chr1_bins[1].avg_coverage, 2.5);
    assert_eq!(chr1_bins[1].avg_overlap_coverage, 0.0);
    assert_eq!(chr1_bins[1].scaling_factor, 0.0);

    let chr2_bins = bins_by_chr.get("chr2").expect("chr2 bins should exist");
    assert_eq!(chr2_bins.len(), 1);
    assert_eq!(chr2_bins[0].start(), 0);
    assert_eq!(chr2_bins[0].end(), 5);
    assert_eq!(chr2_bins[0].avg_coverage, 0.0);

    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_average_tsv_rejects_empty_file() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(&tempdir, "coverage.avg.tsv", &[])?;

    // Act
    let err = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
    )
    .expect_err("empty file should fail");

    // Assert
    assert!(err.to_string().contains("empty file; header required"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_average_tsv_rejects_unexpected_header() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.avg.tsv",
        &[
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
    )
    .expect_err("wrong header should fail");

    // Assert
    assert!(err.to_string().contains("unexpected fcoverage header"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_average_tsv_rejects_short_rows() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.avg.tsv",
        &[
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.25",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
    )
    .expect_err("short row should fail");

    // Assert
    assert!(err.to_string().contains("expected 5 tab-separated columns"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_average_tsv_rejects_invalid_start() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.avg.tsv",
        &[
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
            "chr1\tzero\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
    )
    .expect_err("invalid start should fail");

    // Assert
    assert!(err.to_string().contains("invalid start"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_average_tsv_rejects_invalid_interval() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.avg.tsv",
        &[
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
            "chr1\t10\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
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
fn load_stride_bins_from_fcoverage_average_tsv_rejects_missing_requested_chromosome() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.avg.tsv",
        &[
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string(), "chr2".to_string()],
    )
    .expect_err("missing requested chromosome should fail");

    // Assert
    assert!(err.to_string().contains("missing chromosome 'chr2'"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_average_tsv_rejects_first_bin_not_starting_at_zero() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.avg.tsv",
        &[
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
            "chr1\t5\t10\t1.25\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
    )
    .expect_err("first bin starting after zero should fail");

    // Assert
    assert!(err.to_string().contains("did not start at 0"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_average_tsv_rejects_gap_between_bins() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.avg.tsv",
        &[
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.0\t0",
            "chr1\t12\t20\t2.0\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
    )
    .expect_err("gap between bins should fail");

    // Assert
    assert!(err.to_string().contains("non-contiguous stride bins"));
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_average_tsv_rejects_overlap_between_bins() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.avg.tsv",
        &[
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.0\t0",
            "chr1\t8\t20\t2.0\t0",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_average_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
    )
    .expect_err("overlap between bins should fail");

    // Assert
    assert!(err.to_string().contains("non-contiguous stride bins"));
    Ok(())
}
