use super::*;
use crate::commands::fcoverage::config::LengthNormalizationMode;
use anyhow::Result;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
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
        10,
    )?;

    // Assert
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert_eq!(chr1_bins.len(), 2);
    assert_eq!(chr1_bins[0].start(), 0);
    assert_eq!(chr1_bins[0].end(), 10);
    assert_eq!(chr1_bins[0].stride_value, 1.25);
    assert_eq!(chr1_bins[1].start(), 10);
    assert_eq!(chr1_bins[1].end(), 20);
    assert_eq!(chr1_bins[1].eligible_positions, 7);
    assert!((chr1_bins[1].support_ratio - 0.7).abs() <= 1e-12);
    assert_eq!(chr1_bins[1].stride_value, 2.5);
    assert_eq!(chr1_bins[1].smoothed_value, 0.0);
    assert_eq!(chr1_bins[1].scaling_factor, 0.0);

    let chr2_bins = bins_by_chr.get("chr2").expect("chr2 bins should exist");
    assert_eq!(chr2_bins.len(), 1);
    assert_eq!(chr2_bins[0].start(), 0);
    assert_eq!(chr2_bins[0].end(), 5);
    assert!((chr2_bins[0].support_ratio - 0.5).abs() <= 1e-12);
    assert_eq!(chr2_bins[0].stride_value, 0.0);

    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_marks_fully_blacklisted_rows_as_nan() -> Result<()> {
    // Arrange:
    // Total outputs from fcoverage use raw sums, so a fully blacklisted stride can appear as
    // finite zero even though it has no eligible denominator. The loader must convert that row
    // to missing support before smoothing.
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.total.tsv",
        &[
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t10\t0\t10",
            "chr1\t10\t20\t5\t0",
        ],
    )?;

    // Act
    let bins_by_chr = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "total_coverage",
        10,
    )?;

    // Assert
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert_eq!(chr1_bins[0].eligible_positions, 0);
    assert_eq!(chr1_bins[0].support_ratio, 0.0);
    assert!(
        chr1_bins[0].stride_value.is_nan(),
        "fully blacklisted rows must be missing support"
    );
    assert_eq!(chr1_bins[1].eligible_positions, 10);
    assert_eq!(chr1_bins[1].support_ratio, 1.0);
    assert_eq!(chr1_bins[1].stride_value, 5.0);
    Ok(())
}

#[test]
fn load_stride_bins_from_fcoverage_tsv_rejects_blacklisted_positions_above_span() -> Result<()> {
    // Arrange
    let tempdir = TempDir::new()?;
    let path = write_plain_tsv(
        &tempdir,
        "coverage.average.tsv",
        &[
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t10\t1.25\t11",
        ],
    )?;

    // Act
    let err = load_stride_bins_from_fcoverage_tsv(
        &make_run_result(path),
        &["chr1".to_string()],
        "average_coverage",
        10,
    )
    .expect_err("blacklisted count above row span should fail");

    // Assert
    assert!(
        err.to_string()
            .contains("blacklisted_positions 11 exceeds row span"),
        "unexpected error message: {err}"
    );
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
        10,
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
        10,
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
        10,
    )?;

    // Assert
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert_eq!(chr1_bins.len(), 2);
    assert_eq!(chr1_bins[0].stride_value, 12.5);
    assert_eq!(chr1_bins[1].stride_value, 25.0);
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
        10,
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
        10,
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
        10,
    )
    .expect_err("zero-length interval should fail");

    // Assert
    let error_chain = err
        .chain()
        .map(|source| source.to_string())
        .collect::<Vec<_>>()
        .join("\n");
    assert!(
        error_chain.contains("invalid stride-bin interval 10..10")
            && error_chain.contains("interval end (10) must be greater than start (10)"),
        "unexpected error chain: {error_chain}"
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
        10,
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
        10,
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
        10,
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
        10,
    )
    .expect_err("overlap between bins should fail");

    // Assert
    assert!(err.to_string().contains("non-contiguous stride bins"));
    Ok(())
}

#[test]
fn fill_triangular_overlap_preserves_constant_support_with_short_final_bin() -> Result<()> {
    // Arrange:
    // Four stride bins cover a 35 bp chromosome with stride 10:
    // - three full bins: [0,10), [10,20), [20,30)
    // - one short final bin: [30,35)
    //
    // Every base has the same support, so every bin has stride value 4.0.
    // Smoothing should not change a constant signal: the weighted numerator and
    // denominator should represent the same in-chromosome bases, including the
    // half-length final bin.
    let mut bins = vec![
        StrideBin {
            interval: Interval::new(0, 10)?,
            eligible_positions: 10,
            support_ratio: 1.0,
            stride_value: 4.0,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
        StrideBin {
            interval: Interval::new(10, 20)?,
            eligible_positions: 10,
            support_ratio: 1.0,
            stride_value: 4.0,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
        StrideBin {
            interval: Interval::new(20, 30)?,
            eligible_positions: 10,
            support_ratio: 1.0,
            stride_value: 4.0,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
        StrideBin {
            interval: Interval::new(30, 35)?,
            eligible_positions: 5,
            support_ratio: 0.5,
            stride_value: 4.0,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
    ];

    // Act
    fill_triangular_overlap(&mut bins, 30, 10);

    // Assert
    for bin in &bins {
        assert!(
            (bin.smoothed_value - 4.0).abs() <= 1e-6,
            "constant support should remain 4.0 after smoothing for {:?}, got {}",
            bin.interval,
            bin.smoothed_value
        );
    }

    Ok(())
}

#[test]
fn fill_triangular_overlap_skips_nan_stride_averages_when_smoothing() -> Result<()> {
    // Arrange:
    // A fully blacklisted stride from fcoverage is `NaN`, not zero.
    // That bin has no denominator, so it should not contribute either value
    // or weight to neighboring smoothed bins.
    let mut bins = vec![
        StrideBin {
            interval: Interval::new(0, 10)?,
            eligible_positions: 10,
            support_ratio: 1.0,
            stride_value: 0.0,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
        StrideBin {
            interval: Interval::new(10, 20)?,
            eligible_positions: 10,
            support_ratio: 1.0,
            stride_value: f32::NAN,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
        StrideBin {
            interval: Interval::new(20, 30)?,
            eligible_positions: 10,
            support_ratio: 1.0,
            stride_value: 1.0,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
    ];

    // Act
    fill_triangular_overlap(&mut bins, 20, 10);

    // Assert:
    // With bin-size=20 and stride=10, the triangular kernel is [1,2,1].
    // - row 0: truncated [2,1] over [0,NaN], skipping NaN -> 0 / 2 = 0
    // - row 1: [1,2,1] over [0,NaN,1], skipping NaN -> 1 / (1+1) = 1/2
    // - row 2: truncated [1,2] over [NaN,1], skipping NaN -> 2 / 2 = 1
    let expected_smoothed_values = [0.0_f32, 0.5, 1.0];
    for (bin_index, expected) in expected_smoothed_values.iter().enumerate() {
        let actual = bins[bin_index].smoothed_value;
        assert!(
            (actual - expected).abs() <= 1e-6,
            "expected smoothed value {expected} at bin {bin_index}, got {actual}"
        );
    }

    Ok(())
}

#[test]
fn fill_triangular_overlap_weights_stride_averages_by_eligible_positions() -> Result<()> {
    // Arrange:
    // With bin-size=20 and stride=10, the center row uses triangular weights [1,2,1].
    // The center stride has only one eligible base, so its support is 1/10 of a full stride:
    //   weighted sum = 0*1 + 10*(2*0.1) + 0*1 = 2
    //   weight sum   = 1   +    (2*0.1) + 1   = 2.2
    //   smoothed center = 2 / 2.2 = 10/11
    let mut bins = vec![
        StrideBin {
            interval: Interval::new(0, 10)?,
            eligible_positions: 10,
            support_ratio: 1.0,
            stride_value: 0.0,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
        StrideBin {
            interval: Interval::new(10, 20)?,
            eligible_positions: 1,
            support_ratio: 0.1,
            stride_value: 10.0,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
        StrideBin {
            interval: Interval::new(20, 30)?,
            eligible_positions: 10,
            support_ratio: 1.0,
            stride_value: 0.0,
            smoothed_value: 0.0,
            scaling_factor: 0.0,
        },
    ];

    // Act
    fill_triangular_overlap(&mut bins, 20, 10);

    // Assert
    assert!(
        (bins[1].smoothed_value - (10.0 / 11.0)).abs() <= 1e-6,
        "partly blacklisted center stride should be downweighted by eligible support, got {}",
        bins[1].smoothed_value
    );

    Ok(())
}

#[test]
fn normalize_weighted_average_overlap_by_global_mean_ignores_bins_below_support_floor() -> Result<()> {
    // Arrange
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: 1.0,
                smoothed_value: 0.5,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 20)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: 0.0,
                smoothed_value: 5e-11,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_weighted_average_overlap_by_global_mean(&mut bins_by_chr, false, true)?;

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
fn normalize_weighted_average_overlap_by_global_mean_ignores_nan_raw_or_smoothed_bins() -> Result<()> {
    // Arrange:
    // A NaN raw stride value means the stride had no eligible denominator.
    // A NaN smoothed value means no finite neighboring support was available.
    // Both cases stay visible in the diagnostic columns, but neither should participate in
    // the global mean or receive a usable scaling factor.
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: f32::NAN,
                smoothed_value: 4.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 20)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: 1.0,
                smoothed_value: f32::NAN,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(20, 30)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: 0.0,
                smoothed_value: 0.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(30, 40)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: 2.0,
                smoothed_value: 2.0,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_weighted_average_overlap_by_global_mean(&mut bins_by_chr, false, true)?;

    // Assert
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert!((mean - 2.0).abs() <= 1e-6, "expected mean 2.0, got {mean}");
    assert_eq!(
        chr1_bins[0].scaling_factor, 0.0,
        "NaN raw support should get scaling 0 even when smoothed support is finite"
    );
    assert_eq!(
        chr1_bins[1].scaling_factor, 0.0,
        "NaN smoothed support should get scaling 0"
    );
    assert_eq!(
        chr1_bins[2].scaling_factor, 0.0,
        "zero smoothed support should get scaling 0"
    );
    // Since only this bin is included in the mean, the normalization leads
    // to a scaling factor of 2/2=1
    assert!(
        (chr1_bins[3].scaling_factor - 1.0).abs() <= 1e-6,
        "only finite non-zero bin should normalize to 1.0, got {}",
        chr1_bins[3].scaling_factor
    );

    Ok(())
}

#[test]
fn normalize_weighted_average_overlap_by_global_mean_weights_by_eligible_positions() -> Result<()> {
    // Arrange:
    // Both bins have the same physical length, but different eligible support after masking:
    //   mean = (1*1 + 3*9) / (1 + 9) = 28/10 = 2.8
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                eligible_positions: 1,
                support_ratio: 0.1,
                stride_value: 1.0,
                smoothed_value: 1.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 20)?,
                eligible_positions: 9,
                support_ratio: 0.9,
                stride_value: 3.0,
                smoothed_value: 3.0,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_weighted_average_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

    // Assert
    assert!((mean - 2.8).abs() <= 1e-6, "expected mean 2.8, got {mean}");
    let chr1_bins = bins_by_chr.get("chr1").expect("chr1 bins should exist");
    assert!(
        (chr1_bins[0].scaling_factor - 2.8).abs() <= 1e-6,
        "expected low-support bin scaling 2.8, got {}",
        chr1_bins[0].scaling_factor
    );
    assert!(
        (chr1_bins[1].scaling_factor - (2.8 / 3.0)).abs() <= 1e-6,
        "expected high-support bin scaling 2.8/3, got {}",
        chr1_bins[1].scaling_factor
    );

    Ok(())
}

#[test]
fn normalize_weighted_average_overlap_by_global_mean_explains_all_zero_smoothed_mass() -> Result<()> {
    // Arrange
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: 0.0,
                smoothed_value: 0.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 20)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: 0.0,
                smoothed_value: 5e-11,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let err = normalize_weighted_average_overlap_by_global_mean(&mut bins_by_chr, true, true)
        .expect_err("all-zero smoothed mass should fail");

    // Assert
    let message = err.to_string();
    assert!(
        message.contains("no usable finite non-zero smoothed fragment mass after filtering"),
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
fn normalize_weighted_average_overlap_by_global_mean_keeps_bins_above_support_floor() -> Result<()> {
    // Arrange
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: 1.0,
                smoothed_value: 1.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 20)?,
                eligible_positions: 10,
                support_ratio: 1.0,
                stride_value: 0.0,
                smoothed_value: 2e-9,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_weighted_average_overlap_by_global_mean(&mut bins_by_chr, false, true)?;

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
