#![cfg(feature = "cmd_coverage_weights")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs};
use cfdnalab::commands::coverage_weights::config::CoverageWeightsConfig;
use cfdnalab::commands::coverage_weights::coverage_weights::run;
use fixtures::{ReadSpec, bam_from_specs, simple_inward_bam};
use tempfile::TempDir;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
        chromosomes_file: None,
    }
}

#[derive(Debug)]
struct ScalingRow {
    chromosome: String,
    start: u64,
    end: u64,
    avg_pos_cov: f64,
    avg_overlapping_pos_cov: f64,
    scaling_factor: f64,
}

fn parse_scaling_rows(tsv_path: &std::path::Path) -> Result<Vec<ScalingRow>> {
    let content = std::fs::read_to_string(tsv_path)?;
    let mut lines = content.lines();
    let header = lines.next().unwrap_or("");
    assert_eq!(
        header,
        "chromosome\tstart\tend\tavg_pos_cov\tavg_overlapping_pos_cov\tscaling_factor"
    );

    let mut rows = Vec::new();
    for line in lines {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 6, "Unexpected column count for line: {line}");
        rows.push(ScalingRow {
            chromosome: parts[0].to_string(),
            start: parts[1].parse()?,
            end: parts[2].parse()?,
            avg_pos_cov: parts[3].parse()?,
            avg_overlapping_pos_cov: parts[4].parse()?,
            scaling_factor: parts[5].parse()?,
        });
    }
    Ok(rows)
}

fn assert_approx(actual: f64, expected: f64, tolerance: f64, label: &str) {
    let difference = (actual - expected).abs();
    assert!(
        difference <= tolerance,
        "{label}: expected {expected}, got {actual} (difference {difference}, tolerance {tolerance})"
    );
}

fn make_simple_coverage_weights_config(
    out_dir: &std::path::Path,
    bam: &std::path::Path,
) -> CoverageWeightsConfig {
    let mut cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.to_path_buf(),
            output_dir: out_dir.to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_bin_size(40);
    cfg.set_stride(20);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_output_prefix("coverage".to_string());
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }
    cfg
}

fn single_read_fragment_bam(name: &str) -> Result<fixtures::BamFixture> {
    bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: 20,
            cigar: vec![('M', 60)],
            seq: vec![b'A'; 60],
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
        name,
    )
}

#[test]
fn coverage_scaling_written_with_expected_ranges() -> Result<()> {
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let cfg = make_simple_coverage_weights_config(out_dir.path(), &bam.bam);

    run(&cfg)?;

    let tsv_path = out_dir.path().join("coverage.scaling_factors.tsv");
    assert!(tsv_path.exists());
    let rows = parse_scaling_rows(&tsv_path)?;
    let mut saw_zero = false;
    let mut saw_non_zero = false;
    for row in rows {
        if row.scaling_factor == 0.0 {
            saw_zero = true;
        }
        if row.start >= 20 && row.start < 80 && row.scaling_factor > 0.0 {
            saw_non_zero = true;
        }
    }
    assert!(saw_zero, "expected uncovered stride bin with scaling 0");
    assert!(
        saw_non_zero,
        "expected covered stride bin with positive scaling"
    );

    Ok(())
}

#[test]
fn given_simple_fragment_when_coverage_weights_run_then_output_bins_cover_chromosome_without_gaps()
-> Result<()> {
    // Arrange
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let cfg = make_simple_coverage_weights_config(out_dir.path(), &bam.bam);

    // Act
    run(&cfg)?;
    let tsv_path = out_dir.path().join("coverage.scaling_factors.tsv");
    let rows = parse_scaling_rows(&tsv_path)?;

    // Assert
    assert_eq!(
        rows.len(),
        10,
        "Expected 10 stride bins for chr1 length 200 and stride 20"
    );
    assert_eq!(rows[0].chromosome, "chr1");
    assert_eq!(rows[0].start, 0);
    assert_eq!(rows[0].end, 20);

    for (row_index, row) in rows.iter().enumerate() {
        assert_eq!(row.chromosome, "chr1");
        assert_eq!(
            row.end - row.start,
            20,
            "Expected fixed stride-sized bins in this fixture"
        );
        if row_index > 0 {
            let previous = &rows[row_index - 1];
            assert_eq!(
                row.start, previous.end,
                "Expected bins to be contiguous without gaps or overlaps"
            );
        }
    }
    assert_eq!(rows.last().unwrap().end, 200);

    Ok(())
}

#[test]
fn given_simple_fragment_when_coverage_weights_run_then_scaling_values_match_hand_derivation()
-> Result<()> {
    // Arrange
    //
    // Fixture:
    // - Chromosome length = 200
    // - One fragment spans [20,80), so average stride-bin coverage (stride=20) is:
    //   [0, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    // - With bin-size=40 and stride=20, half-window=(40/20)-1=1
    //   triangular weights are [1,2,1]
    // - Hand-derived avg-overlapping-position-coverage per bin:
    //   [1/3, 3/4, 1, 3/4, 1/4, 0, 0, 0, 0, 0]
    // - Mean over non-zero bins:
    //   (1/3 + 3/4 + 1 + 3/4 + 1/4) / 5 = 37/60
    // - Scaling factors are inverted normalized overlap:
    //   scaling = 1 / (overlap / (37/60))
    //   -> [37/20, 37/45, 37/60, 37/45, 37/15, 0, 0, 0, 0, 0]
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let cfg = make_simple_coverage_weights_config(out_dir.path(), &bam.bam);

    // Act
    run(&cfg)?;
    let tsv_path = out_dir.path().join("coverage.scaling_factors.tsv");
    let rows = parse_scaling_rows(&tsv_path)?;

    // Assert
    let expected_avg_pos_cov = [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let expected_avg_overlap = [
        1.0 / 3.0,
        3.0 / 4.0,
        1.0,
        3.0 / 4.0,
        1.0 / 4.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    let expected_scaling = [
        37.0 / 20.0,
        37.0 / 45.0,
        37.0 / 60.0,
        37.0 / 45.0,
        37.0 / 15.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];

    assert_eq!(rows.len(), 10);
    for row_index in 0..rows.len() {
        let row = &rows[row_index];
        assert_approx(
            row.avg_pos_cov,
            expected_avg_pos_cov[row_index],
            1e-6,
            &format!("avg_pos_cov at row {row_index}"),
        );
        assert_approx(
            row.avg_overlapping_pos_cov,
            expected_avg_overlap[row_index],
            1e-6,
            &format!("avg_overlapping_pos_cov at row {row_index}"),
        );
        assert_approx(
            row.scaling_factor,
            expected_scaling[row_index],
            1e-6,
            &format!("scaling_factor at row {row_index}"),
        );
    }

    Ok(())
}

#[test]
fn given_unpaired_read_fragment_when_coverage_weights_run_then_scaling_matches_same_fragment_span()
-> Result<()> {
    // Arrange
    //
    // Hand derivation:
    // - The unpaired fixture has one aligned read covering [20,80) on chr1
    // - In `reads_are_fragments` mode, the command counts that read as exactly one fragment
    // - That produces the same stride-bin coverage as the paired simple fixture:
    //   [0, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    // - The triangular smoothing and scaling factors are therefore identical:
    //   overlap  = [1/3, 3/4, 1, 3/4, 1/4, 0, 0, 0, 0, 0]
    //   scaling  = [37/20, 37/45, 37/60, 37/45, 37/15, 0, 0, 0, 0, 0]
    let bam = single_read_fragment_bam("coverage_weights_unpaired_single_read")?;
    let out_dir = TempDir::new()?;
    let mut cfg = make_simple_coverage_weights_config(out_dir.path(), &bam.bam);
    cfg.unpaired.reads_are_fragments = true;

    // Act
    run(&cfg)?;
    let tsv_path = out_dir.path().join("coverage.scaling_factors.tsv");
    let rows = parse_scaling_rows(&tsv_path)?;

    // Assert
    let expected_avg_pos_cov = [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let expected_avg_overlap = [
        1.0 / 3.0,
        3.0 / 4.0,
        1.0,
        3.0 / 4.0,
        1.0 / 4.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    let expected_scaling = [
        37.0 / 20.0,
        37.0 / 45.0,
        37.0 / 60.0,
        37.0 / 45.0,
        37.0 / 15.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];

    assert_eq!(rows.len(), 10);
    for row_index in 0..rows.len() {
        let row = &rows[row_index];
        assert_approx(
            row.avg_pos_cov,
            expected_avg_pos_cov[row_index],
            1e-6,
            &format!("unpaired avg_pos_cov at row {row_index}"),
        );
        assert_approx(
            row.avg_overlapping_pos_cov,
            expected_avg_overlap[row_index],
            1e-6,
            &format!("unpaired avg_overlapping_pos_cov at row {row_index}"),
        );
        assert_approx(
            row.scaling_factor,
            expected_scaling[row_index],
            1e-6,
            &format!("unpaired scaling_factor at row {row_index}"),
        );
    }

    Ok(())
}

#[test]
fn check_bin_sizes_rejects_invalid_stride() {
    let mut cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: std::path::PathBuf::new(),
            output_dir: std::path::PathBuf::new(),
            n_threads: 1,
        },
        ChromosomeArgs::default(),
    );
    cfg.set_bin_size(30);
    cfg.set_stride(40);
    let err = cfg.check_bin_sizes().unwrap_err();
    assert!(format!("{err}").contains("stride"));
}
