#![cfg(feature = "cmd_coverage_weights")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs};
use cfdnalab::commands::coverage_weights::config::CoverageWeightsConfig;
use cfdnalab::commands::coverage_weights::coverage_weights::run;
use cfdnalab::commands::coverage_weights::striding::{
    StrideBin, normalize_average_overlap_by_global_mean,
};
#[cfg(feature = "cmd_fragment_count_weights")]
use cfdnalab::commands::fragment_count_weights::{
    config::FragmentCountWeightsConfig, fragment_count_weights::run as run_fragment_count_weights,
};
use cfdnalab::shared::interval::Interval;
use fixtures::{
    FragmentSpec, ReadSpec, bam_from_specs, bam_from_specs_strict_identity, paired_fragment,
    simple_inward_bam,
};
use fxhash::FxHashMap;
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
    average_pos_coverage: f64,
    average_overlapping_pos_coverage: f64,
    scaling_factor: f64,
}

fn parse_scaling_rows(tsv_path: &std::path::Path) -> Result<Vec<ScalingRow>> {
    let content = std::fs::read_to_string(tsv_path)?;
    let mut lines = content.lines();
    let first_line = lines.next().unwrap_or("");
    let header = if let Some(gc_mode) = first_line.strip_prefix("# gc_mode=") {
        assert!(
            !gc_mode.is_empty(),
            "expected non-empty gc_mode metadata in {}",
            tsv_path.display()
        );
        lines.next().unwrap_or("")
    } else {
        first_line
    };
    assert_eq!(
        header,
        "chromosome\tstart\tend\taverage_pos_coverage\taverage_overlapping_pos_coverage\tscaling_factor"
    );

    let mut rows = Vec::new();
    for line in lines {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 6, "Unexpected column count for line: {line}");
        rows.push(ScalingRow {
            chromosome: parts[0].to_string(),
            start: parts[1].parse()?,
            end: parts[2].parse()?,
            average_pos_coverage: parts[3].parse()?,
            average_overlapping_pos_coverage: parts[4].parse()?,
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

fn scaling_row_chromosomes(rows: &[ScalingRow]) -> Vec<String> {
    rows.iter().map(|row| row.chromosome.clone()).collect()
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

#[cfg(feature = "cmd_fragment_count_weights")]
fn make_simple_fragment_count_weights_config(
    out_dir: &std::path::Path,
    bam: &std::path::Path,
) -> FragmentCountWeightsConfig {
    let mut cfg = FragmentCountWeightsConfig::new(
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
    cfg.set_output_prefix("counts".to_string());
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

fn multi_chrom_order_bam(name: &str) -> Result<fixtures::BamFixture> {
    // BAM header order is intentionally non-lexicographic:
    //   chr2, chr10, chr1
    //
    // We add one 15 bp fragment on each chromosome using 5 bp reads, so the paired-read fixture
    // is physically valid (`read_len <= fragment_len`) while still giving one simple fragment per
    // chromosome. With `bin_size = stride = 20`, each
    // chromosome contributes exactly one TSV row, so the written chromosome sequence directly
    // exposes the command's row-order contract.
    bam_from_specs(
        vec![
            ("chr2".to_string(), 20),
            ("chr10".to_string(), 20),
            ("chr1".to_string(), 20),
        ],
        vec![
            fragment_on_tid(paired_fragment(0, 15, 5), 0),
            fragment_on_tid(paired_fragment(0, 15, 5), 1),
            fragment_on_tid(paired_fragment(0, 15, 5), 2),
        ],
        Vec::new(),
        name,
    )
}

fn fragment_on_tid(mut fragment: FragmentSpec, tid: usize) -> FragmentSpec {
    fragment.forward.tid = tid;
    fragment.reverse.tid = tid;
    fragment.forward.mate_tid = Some(tid);
    fragment.reverse.mate_tid = Some(tid);
    fragment
}

#[cfg(feature = "cmd_fragment_count_weights")]
fn mixed_length_fragment_bam(name: &str) -> Result<fixtures::BamFixture> {
    // Two fragments on chr1:
    // - one short fragment [0,20)
    // - one long fragment [20,80)
    //
    // This makes the two weighting schemes diverge cleanly when bin_size == stride == 20.
    bam_from_specs(
        vec![("chr1".to_string(), 100)],
        vec![paired_fragment(0, 20, 10), paired_fragment(20, 60, 10)],
        Vec::new(),
        name,
    )
}

#[test]
fn coverage_scaling_written_with_expected_ranges() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let cfg = make_simple_coverage_weights_config(out_dir.path(), &bam.bam);

    run(&cfg)?;

    let tsv_path = out_dir.path().join("coverage.coverage.scaling_factors.tsv");
    assert!(tsv_path.exists());
    let rows = parse_scaling_rows(&tsv_path)?;
    let expected_non_zero_rows = [
        true, true, true, true, true, false, false, false, false, false,
    ];

    assert_eq!(
        rows.len(),
        10,
        "expected one stride bin per 20 bp across chr1"
    );
    for (row_index, row) in rows.iter().enumerate() {
        assert_eq!(row.chromosome, "chr1");
        assert_eq!(row.start, (row_index as u64) * 20);
        assert_eq!(row.end, row.start + 20);

        let should_be_non_zero = expected_non_zero_rows[row_index];
        if should_be_non_zero {
            assert!(
                row.average_overlapping_pos_coverage > 0.0,
                "row {row_index} should overlap the smoothed coverage support"
            );
            assert!(
                row.scaling_factor > 0.0,
                "row {row_index} should have a positive scaling factor"
            );
        } else {
            assert_eq!(
                row.average_overlapping_pos_coverage, 0.0,
                "row {row_index} should be outside the smoothed coverage support"
            );
            assert_eq!(
                row.scaling_factor, 0.0,
                "row {row_index} should have zero scaling outside coverage support"
            );
        }
    }

    Ok(())
}

#[cfg(feature = "cmd_fragment_count_weights")]
#[test]
fn fragment_count_scaling_written_with_expected_ranges() -> Result<()> {
    // Human verification status: unverified
    // Arrange: simple_inward_bam has chr1 length 200 and one fragment spanning [20,80).
    // With stride 20 this yields exactly 10 stride bins in the written scaling TSV.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let cfg = make_simple_fragment_count_weights_config(out_dir.path(), &bam.bam);

    // Act
    run_fragment_count_weights(&cfg)?;

    // Assert
    let tsv_path = out_dir
        .path()
        .join("counts.fragment_counts.scaling_factors.tsv");
    assert!(tsv_path.exists());
    let rows = parse_scaling_rows(&tsv_path)?;

    assert_eq!(
        rows.len(),
        10,
        "expected one stride bin per 20 bp across chr1"
    );
    assert_eq!(rows[0].chromosome, "chr1");
    assert_eq!(rows[0].start, 0);
    assert_eq!(rows.last().unwrap().end, 200);

    Ok(())
}

#[test]
fn given_simple_fragment_when_coverage_weights_run_then_output_bins_cover_chromosome_without_gaps()
-> Result<()> {
    // Human verification status: unverified
    // Arrange
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let cfg = make_simple_coverage_weights_config(out_dir.path(), &bam.bam);

    // Act
    run(&cfg)?;
    let tsv_path = out_dir.path().join("coverage.coverage.scaling_factors.tsv");
    let rows = parse_scaling_rows(&tsv_path)?;

    // Assert
    assert_eq!(
        rows.len(),
        10,
        "Expected 10 stride bins for chr1 length 200 and stride 20"
    );
    for (row_index, row) in rows.iter().enumerate() {
        assert_eq!(row.chromosome, "chr1");
        assert_eq!(
            row.start,
            (row_index as u64) * 20,
            "Expected the row starts to enumerate the full chromosome in stride steps"
        );
        assert_eq!(
            row.end,
            ((row_index as u64) + 1) * 20,
            "Expected the row ends to enumerate the full chromosome in stride steps"
        );
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

#[cfg(feature = "cmd_fragment_count_weights")]
#[test]
fn fragment_count_weights_differs_from_coverage_weights_for_mixed_fragment_lengths() -> Result<()> {
    // Human verification status: unverified
    // Arrange
    //
    // We choose bin_size == stride == 20 so the smoothed value equals the raw stride-bin average.
    // The BAM has two fragments on chr1:
    // - short: [0,20), length 20
    // - long:  [20,80), length 60
    //
    // Coverage weights:
    // - each covered base gets weight 1.0
    // - bins [0,20), [20,40), [40,60), [60,80) each get average_pos_coverage = 1.0
    // - global mean over non-zero bins = 1.0
    // - scaling_factor for each covered bin = 1.0
    //
    // Fragment-count weights:
    // - short fragment contributes 1/20 per base in [0,20)       -> average_pos_coverage = 0.05
    // - long fragment contributes 1/60 per base in [20,80)       -> average_pos_coverage = 1/60
    // - global mean over the four non-zero bins is:
    //     (0.05 + 3*(1/60)) / 4 = (0.05 + 0.05) / 4 = 0.025
    // - scaling factors become:
    //     short bin: 0.025 / 0.05    = 0.5
    //     long bins: 0.025 / (1/60)  = 1.5
    let bam = mixed_length_fragment_bam("weights_mixed_lengths")?;
    let coverage_out = TempDir::new()?;
    let counts_out = TempDir::new()?;

    let mut coverage_cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: coverage_out.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1"]),
    );
    coverage_cfg.set_bin_size(20);
    coverage_cfg.set_stride(20);
    coverage_cfg.set_min_mapq(0);
    coverage_cfg.set_require_proper_pair(false);
    coverage_cfg.set_output_prefix("coverage".to_string());
    {
        let frag = coverage_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }

    let mut counts_cfg = FragmentCountWeightsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: counts_out.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1"]),
    );
    counts_cfg.set_bin_size(20);
    counts_cfg.set_stride(20);
    counts_cfg.set_min_mapq(0);
    counts_cfg.set_require_proper_pair(false);
    counts_cfg.set_output_prefix("counts".to_string());
    {
        let frag = counts_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }

    // Act
    run(&coverage_cfg)?;
    run_fragment_count_weights(&counts_cfg)?;

    let coverage_rows = parse_scaling_rows(
        &coverage_out
            .path()
            .join("coverage.coverage.scaling_factors.tsv"),
    )?;
    let counts_rows = parse_scaling_rows(
        &counts_out
            .path()
            .join("counts.fragment_counts.scaling_factors.tsv"),
    )?;

    // Assert
    assert_eq!(coverage_rows.len(), 5);
    assert_eq!(counts_rows.len(), 5);

    // Covered bins stay at scaling factor 1.0 in coverage mode because all four
    // non-zero bins have the same average support.
    for row_index in 0..4 {
        assert_approx(
            coverage_rows[row_index].average_pos_coverage,
            1.0,
            1e-6,
            "coverage average_pos_coverage",
        );
        assert_approx(
            coverage_rows[row_index].scaling_factor,
            1.0,
            1e-6,
            "coverage scaling_factor",
        );
    }
    assert_eq!(coverage_rows[4].average_pos_coverage, 0.0);
    assert_eq!(coverage_rows[4].scaling_factor, 0.0);

    // In fragment-count mode the short fragment dominates its one bin less than a long
    // fragment would in coverage mode, while the long fragment is spread across three bins.
    assert_approx(
        counts_rows[0].average_pos_coverage,
        0.05,
        1e-6,
        "fragment-count short-bin average_pos_coverage",
    );
    for row_index in 1..4 {
        assert_approx(
            counts_rows[row_index].average_pos_coverage,
            1.0 / 60.0,
            1e-6,
            "fragment-count long-bin average_pos_coverage",
        );
    }
    assert_approx(
        counts_rows[0].scaling_factor,
        0.5,
        1e-6,
        "fragment-count short-bin scaling_factor",
    );
    for row_index in 1..4 {
        assert_approx(
            counts_rows[row_index].scaling_factor,
            1.5,
            1e-6,
            "fragment-count long-bin scaling_factor",
        );
    }
    assert_eq!(counts_rows[4].average_pos_coverage, 0.0);
    assert_eq!(counts_rows[4].scaling_factor, 0.0);

    Ok(())
}

#[test]
fn given_simple_fragment_when_coverage_weights_run_then_scaling_values_match_hand_derivation()
-> Result<()> {
    // Human verification status: unverified
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
    let tsv_path = out_dir.path().join("coverage.coverage.scaling_factors.tsv");
    let rows = parse_scaling_rows(&tsv_path)?;

    // Assert
    let expected_average_pos_coverage = [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let expected_average_overlap = [
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
            row.average_pos_coverage,
            expected_average_pos_coverage[row_index],
            1e-6,
            &format!("average_pos_coverage at row {row_index}"),
        );
        assert_approx(
            row.average_overlapping_pos_coverage,
            expected_average_overlap[row_index],
            1e-6,
            &format!("average_overlapping_pos_coverage at row {row_index}"),
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
    // Human verification status: unverified
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
    let tsv_path = out_dir.path().join("coverage.coverage.scaling_factors.tsv");
    let rows = parse_scaling_rows(&tsv_path)?;

    // Assert
    let expected_average_pos_coverage = [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let expected_average_overlap = [
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
            row.average_pos_coverage,
            expected_average_pos_coverage[row_index],
            1e-6,
            &format!("unpaired average_pos_coverage at row {row_index}"),
        );
        assert_approx(
            row.average_overlapping_pos_coverage,
            expected_average_overlap[row_index],
            1e-6,
            &format!("unpaired average_overlapping_pos_coverage at row {row_index}"),
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
    // Human verification status: unverified
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

#[test]
fn check_bin_sizes_accepts_valid_stride_values() {
    // Human verification status: unverified
    let mut divisible_cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: std::path::PathBuf::new(),
            output_dir: std::path::PathBuf::new(),
            n_threads: 1,
        },
        ChromosomeArgs::default(),
    );
    divisible_cfg.set_bin_size(40);
    divisible_cfg.set_stride(20);
    divisible_cfg
        .check_bin_sizes()
        .expect("stride dividing bin_size should be accepted");

    let mut equal_cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: std::path::PathBuf::new(),
            output_dir: std::path::PathBuf::new(),
            n_threads: 1,
        },
        ChromosomeArgs::default(),
    );
    equal_cfg.set_bin_size(40);
    equal_cfg.set_stride(40);
    equal_cfg
        .check_bin_sizes()
        .expect("stride equal to bin_size should be accepted");
}

#[test]
fn normalize_average_overlap_keeps_sparse_non_zero_scaling_finite() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Use three stride bins with unequal lengths on one chromosome:
    // - a very sparse but non-zero bin with avg-overlap coverage 0.0001
    // - a typical covered bin with avg-overlap coverage 1.0 and 4x more genomic span
    // - an uncovered bin with avg-overlap coverage 0.0
    //
    // With length-weighted normalization, the zero bin is ignored and the mean is:
    //   mean = (10 * 0.0001 + 40 * 1.0) / (10 + 40)
    //        = (0.001 + 40) / 50
    //        = 40.001 / 50
    //        = 0.80002
    //
    // With inversion enabled, scaling becomes:
    //   sparse bin   = 1 / (0.0001 / 0.80002) = 8000.2
    //   covered bin  = 1 / (1.0    / 0.80002) = 0.80002
    //   zero bin     = 0 by explicit zero-preserving logic
    //
    // This distinguishes length-weighted from simple averaging. A wrong unweighted mean would be
    // 0.50005 instead of 0.80002.
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                average_coverage: 0.0001,
                average_overlap_coverage: 0.0001,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 50)?,
                average_coverage: 1.0,
                average_overlap_coverage: 1.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(50, 60)?,
                average_coverage: 0.0,
                average_overlap_coverage: 0.0,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_average_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

    // Assert
    assert_approx(mean as f64, 0.80002, 1e-7, "global mean before inversion");
    let bins = bins_by_chr
        .get("chr1")
        .expect("chr1 bins should remain present");
    assert!(
        bins[0].scaling_factor.is_finite(),
        "sparse non-zero bin should produce a large but finite scaling factor"
    );
    assert_approx(
        bins[0].scaling_factor as f64,
        8000.2,
        1e-3,
        "sparse-bin scaling factor",
    );
    assert_approx(
        bins[1].scaling_factor as f64,
        0.80002,
        1e-6,
        "covered-bin scaling factor",
    );
    assert_eq!(
        bins[2].scaling_factor, 0.0,
        "zero-overlap bins should remain zero after normalization"
    );

    Ok(())
}

#[test]
fn normalize_average_overlap_support_floor_excludes_tiny_bin_from_mean_and_inversion() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Use three bins with unequal lengths on one chromosome:
    // - one extremely tiny bin with avg-overlap coverage 1e-40
    // - one ordinary covered bin with avg-overlap coverage 1.0 and 4x more genomic span
    // - one uncovered zero bin
    //
    // The support floor now treats the tiny bin as zero-support before both global-mean
    // accumulation and inversion. That means both the tiny bin and the explicit zero bin are
    // excluded from the denominator, leaving only the ordinary covered bin:
    //   mean = (40 * 1.0) / 40 = 1.0
    //
    // With inversion enabled, scaling becomes:
    // - tiny bin     -> 0.0 because it is below the support floor
    // - covered bin  -> 1 / (1.0 / 1.0) = 1.0
    // - zero bin     -> 0.0 by explicit zero-preserving logic
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 10)?,
                average_coverage: 1.0e-40_f32,
                average_overlap_coverage: 1.0e-40_f32,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(10, 50)?,
                average_coverage: 1.0,
                average_overlap_coverage: 1.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(50, 60)?,
                average_coverage: 0.0,
                average_overlap_coverage: 0.0,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_average_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

    // Assert
    assert!(
        mean.is_finite() && mean > 0.0,
        "global mean should remain finite and positive, got {mean}"
    );
    assert_approx(
        mean as f64,
        1.0,
        1e-6,
        "length-weighted global mean with tiny bin below support floor",
    );
    let bins = bins_by_chr
        .get("chr1")
        .expect("chr1 bins should remain present");
    assert_eq!(
        bins[0].scaling_factor, 0.0,
        "extremely tiny bins below the support floor should not be inverted"
    );
    assert_approx(
        bins[1].scaling_factor as f64,
        1.0,
        1e-6,
        "ordinary covered-bin scaling factor with tiny bin excluded",
    );
    assert_eq!(
        bins[2].scaling_factor, 0.0,
        "zero-overlap bins should remain zero with support-floor exclusion"
    );

    Ok(())
}

#[test]
fn normalize_average_overlap_weights_short_final_bin_in_global_mean() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Use two bins on one chromosome, where the final bin is half as long:
    // - [0, 20): avg-overlap coverage 1.0
    // - [20, 30): avg-overlap coverage 3.0
    //
    // With length weighting enabled, the global mean is base-weighted:
    //   mean = (1.0 * 20 + 3.0 * 10) / (20 + 10) = 50 / 30 = 5/3
    //
    // With inversion enabled, scaling becomes:
    //   first bin = 1 / (1 / (5/3)) = 5/3
    //   short bin = 1 / (3 / (5/3)) = 5/9
    //
    // This proves the short final bin does not get counted as a full-length bin in the
    // global denominator, which would incorrectly yield mean = (1 + 3) / 2 = 2.
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 20)?,
                average_coverage: 1.0,
                average_overlap_coverage: 1.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(20, 30)?,
                average_coverage: 3.0,
                average_overlap_coverage: 3.0,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_average_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

    // Assert
    assert_approx(mean as f64, 5.0 / 3.0, 1e-6, "length-weighted global mean");
    let bins = bins_by_chr
        .get("chr1")
        .expect("chr1 bins should remain present");
    assert_approx(
        bins[0].scaling_factor as f64,
        5.0 / 3.0,
        1e-6,
        "full-length bin scaling factor",
    );
    assert_approx(
        bins[1].scaling_factor as f64,
        5.0 / 9.0,
        1e-6,
        "short final bin scaling factor",
    );

    Ok(())
}

#[test]
fn multi_chromosome_scaling_uses_one_shared_global_mean() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Build two chromosomes with different local overlap profiles under the same stride/bin settings.
    //
    // chr1 has one fragment over [20, 80), so the hand-derived avg-overlap coverage is:
    //   [1/3, 3/4, 1, 3/4, 1/4, 0, ...]
    //
    // chr2 has one shorter fragment over [20, 40), so the hand-derived avg-overlap coverage is:
    //   [1/3, 1/2, 1/4, 0, ...]
    //
    // With bin_size=40 and stride=20, all stride bins have equal length, so the global mean is the
    // plain mean over all non-zero avg-overlap bins from both chromosomes:
    //
    //   chr1 non-zero sum = 1/3 + 3/4 + 1 + 3/4 + 1/4 = 37/12
    //   chr2 non-zero sum = 1/3 + 1/2 + 1/4 = 13/12
    //   total non-zero sum = 50/12 = 25/6
    //   number of non-zero bins = 8
    //   shared global mean = (25/6) / 8 = 25/48
    //
    // Inverted scaling factors are therefore:
    //   chr1 at overlap 1     -> (25/48) / 1     = 25/48
    //   chr2 at overlap 1/2   -> (25/48) / (1/2) = 25/24
    //   chr2 at overlap 1/3   -> (25/48) / (1/3) = 25/16
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![
            paired_fragment(20, 60, 20),
            fragment_on_tid(paired_fragment(20, 20, 10), 1),
        ],
        Vec::new(),
        "coverage_weights_multi_chr_shared_mean",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1", "chr2"]),
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

    // Act
    run(&cfg)?;
    let rows = parse_scaling_rows(&out_dir.path().join("coverage.coverage.scaling_factors.tsv"))?;

    // Assert
    let chr1_rows: Vec<_> = rows.iter().filter(|row| row.chromosome == "chr1").collect();
    let chr2_rows: Vec<_> = rows.iter().filter(|row| row.chromosome == "chr2").collect();

    assert_eq!(chr1_rows.len(), 10, "chr1 should contribute 10 stride bins");
    assert_eq!(chr2_rows.len(), 10, "chr2 should contribute 10 stride bins");

    let expected_chr1_overlap = [
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
    let expected_chr1_scaling = [
        25.0 / 16.0,
        25.0 / 36.0,
        25.0 / 48.0,
        25.0 / 36.0,
        25.0 / 12.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    let expected_chr2_overlap = [
        1.0 / 3.0,
        1.0 / 2.0,
        1.0 / 4.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    let expected_chr2_scaling = [
        25.0 / 16.0,
        25.0 / 24.0,
        25.0 / 12.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];

    for (row_index, row) in chr1_rows.iter().enumerate() {
        assert_eq!(row.start, (row_index as u64) * 20);
        assert_eq!(row.end, row.start + 20);
        assert_approx(
            row.average_overlapping_pos_coverage,
            expected_chr1_overlap[row_index],
            1e-6,
            &format!("chr1 avg overlap row {row_index}"),
        );
        assert_approx(
            row.scaling_factor,
            expected_chr1_scaling[row_index],
            1e-6,
            &format!("chr1 scaling row {row_index}"),
        );
    }

    for (row_index, row) in chr2_rows.iter().enumerate() {
        assert_eq!(row.start, (row_index as u64) * 20);
        assert_eq!(row.end, row.start + 20);
        assert_approx(
            row.average_overlapping_pos_coverage,
            expected_chr2_overlap[row_index],
            1e-6,
            &format!("chr2 avg overlap row {row_index}"),
        );
        assert_approx(
            row.scaling_factor,
            expected_chr2_scaling[row_index],
            1e-6,
            &format!("chr2 scaling row {row_index}"),
        );
    }

    Ok(())
}

#[test]
fn coverage_weights_default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero()
-> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Build two otherwise identical inward fragments on the same span [20,80):
    // - one with MAPQ 60
    // - one with MAPQ 20
    //
    // `coverage-weights` intentionally defaults to `min_mapq = 30`, so:
    // - default config must keep only the MAPQ-60 fragment
    // - explicit `min_mapq = 30` must match exactly
    // - explicit `min_mapq = 0` must keep both fragments
    //
    // With `bin_size = 40` and `stride = 20`, keeping exactly one fragment would produce the same
    // hand-derived overlap profile as the standard simple fixture:
    //   avg-pos-cov = [0, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    //   avg-overlap = [1/3, 3/4, 1, 3/4, 1/4, 0, 0, 0, 0, 0]
    //   scaling     = [37/20, 37/45, 37/60, 37/45, 37/15, 0, 0, 0, 0, 0]
    //
    // Keeping both identical fragments doubles `average_pos_coverage` and `average_overlapping_pos_coverage`
    // element-wise, but leaves `scaling_factor` unchanged because both the numerator and the
    // global mean double.
    let high_mapq = paired_fragment(20, 60, 20);
    let mut low_mapq = paired_fragment(20, 60, 20);
    low_mapq.forward.mapq = 20;
    low_mapq.reverse.mapq = 20;
    let bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 200)],
        vec![high_mapq, low_mapq],
        Vec::new(),
        "coverage_weights_default_min_mapq",
    )?;

    let default_out = TempDir::new()?;
    let explicit_thirty_out = TempDir::new()?;
    let explicit_zero_out = TempDir::new()?;

    let mut default_cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: default_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    default_cfg.set_bin_size(40);
    default_cfg.set_stride(20);
    default_cfg.set_output_prefix("default".to_string());
    default_cfg.set_require_proper_pair(false);
    {
        let frag = default_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    let mut explicit_thirty_cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: explicit_thirty_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    explicit_thirty_cfg.set_bin_size(40);
    explicit_thirty_cfg.set_stride(20);
    explicit_thirty_cfg.set_output_prefix("thirty".to_string());
    explicit_thirty_cfg.set_require_proper_pair(false);
    explicit_thirty_cfg.set_min_mapq(30);
    {
        let frag = explicit_thirty_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    let mut explicit_zero_cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: explicit_zero_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    explicit_zero_cfg.set_bin_size(40);
    explicit_zero_cfg.set_stride(20);
    explicit_zero_cfg.set_output_prefix("zero".to_string());
    explicit_zero_cfg.set_require_proper_pair(false);
    explicit_zero_cfg.set_min_mapq(0);
    {
        let frag = explicit_zero_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    // Act
    run(&default_cfg)?;
    run(&explicit_thirty_cfg)?;
    run(&explicit_zero_cfg)?;

    // Assert
    let default_rows = parse_scaling_rows(
        &default_out
            .path()
            .join("default.coverage.scaling_factors.tsv"),
    )?;
    let explicit_thirty_rows = parse_scaling_rows(
        &explicit_thirty_out
            .path()
            .join("thirty.coverage.scaling_factors.tsv"),
    )?;
    let explicit_zero_rows = parse_scaling_rows(
        &explicit_zero_out
            .path()
            .join("zero.coverage.scaling_factors.tsv"),
    )?;

    assert_eq!(default_rows.len(), 10);
    assert_eq!(default_rows.len(), explicit_thirty_rows.len());
    assert_eq!(default_rows.len(), explicit_zero_rows.len());

    for (default_row, explicit_row) in default_rows.iter().zip(explicit_thirty_rows.iter()) {
        assert_eq!(default_row.chromosome, explicit_row.chromosome);
        assert_eq!(default_row.start, explicit_row.start);
        assert_eq!(default_row.end, explicit_row.end);
        assert_approx(
            default_row.average_pos_coverage,
            explicit_row.average_pos_coverage,
            1e-12,
            "default vs explicit-30 average_pos_coverage",
        );
        assert_approx(
            default_row.average_overlapping_pos_coverage,
            explicit_row.average_overlapping_pos_coverage,
            1e-12,
            "default vs explicit-30 average_overlap",
        );
        assert_approx(
            default_row.scaling_factor,
            explicit_row.scaling_factor,
            1e-12,
            "default vs explicit-30 scaling",
        );
    }

    let expected_default_average_pos_coverage = [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let expected_default_average_overlap = [
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
    let expected_non_zero_scaling = [
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

    for (row_index, row) in default_rows.iter().enumerate() {
        assert_approx(
            row.average_pos_coverage,
            expected_default_average_pos_coverage[row_index],
            1e-6,
            &format!("default average_pos_coverage at row {row_index}"),
        );
        assert_approx(
            row.average_overlapping_pos_coverage,
            expected_default_average_overlap[row_index],
            1e-6,
            &format!("default average_overlap at row {row_index}"),
        );
        assert_approx(
            row.scaling_factor,
            expected_non_zero_scaling[row_index],
            1e-6,
            &format!("default scaling factor at row {row_index}"),
        );
    }

    for (row_index, row) in explicit_zero_rows.iter().enumerate() {
        assert_approx(
            row.average_pos_coverage,
            expected_default_average_pos_coverage[row_index] * 2.0,
            1e-6,
            &format!("explicit-zero average_pos_coverage at row {row_index}"),
        );
        assert_approx(
            row.average_overlapping_pos_coverage,
            expected_default_average_overlap[row_index] * 2.0,
            1e-6,
            &format!("explicit-zero average_overlap at row {row_index}"),
        );
        assert_approx(
            row.scaling_factor,
            expected_non_zero_scaling[row_index],
            1e-6,
            &format!("explicit-zero scaling factor at row {row_index}"),
        );
    }

    Ok(())
}

#[test]
fn explicit_chromosome_order_controls_scaling_tsv_row_order() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // The command writes rows by iterating the resolved chromosome list.
    // For explicit `--chromosomes`, `ChromosomeArgs::resolve_chromosomes` returns the user-supplied
    // order unchanged.
    //
    // We therefore build a BAM with header order [chr2, chr10, chr1] but request:
    //   --chromosomes chr10,chr2
    //
    // With `bin_size = stride = 20` and one fragment per chromosome, each selected chromosome
    // yields exactly one TSV row. The output chromosome sequence must therefore be exactly:
    //   [chr10, chr2]
    let bam = multi_chrom_order_bam("coverage_weights_explicit_order_bam")?;
    let out_dir = TempDir::new()?;

    let mut cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 2,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr10".to_string(), "chr2".to_string()]),
            chromosomes_file: None,
        },
    );
    cfg.set_output_prefix("coverage".to_string());
    cfg.set_bin_size(20);
    cfg.set_stride(20);
    cfg.set_min_mapq(0);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    // Act
    run(&cfg)?;
    let rows = parse_scaling_rows(&out_dir.path().join("coverage.coverage.scaling_factors.tsv"))?;

    // Assert
    assert_eq!(
        scaling_row_chromosomes(&rows),
        vec!["chr10".to_string(), "chr2".to_string()]
    );

    Ok(())
}

#[test]
fn chromosomes_all_uses_bam_header_order_for_scaling_tsv_rows() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // `ChromosomeArgs::resolve_chromosomes` documents that `--chromosomes all` uses BAM header
    // order for BAM-backed commands.
    //
    // The BAM header order here is intentionally:
    //   [chr2, chr10, chr1]
    //
    // With one stride/bin row per chromosome, the scaling TSV must preserve that exact order.
    let bam = multi_chrom_order_bam("coverage_weights_all_order_bam")?;
    let out_dir = TempDir::new()?;

    let mut cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 2,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["all".to_string()]),
            chromosomes_file: None,
        },
    );
    cfg.set_output_prefix("coverage".to_string());
    cfg.set_bin_size(20);
    cfg.set_stride(20);
    cfg.set_min_mapq(0);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    // Act
    run(&cfg)?;
    let rows = parse_scaling_rows(&out_dir.path().join("coverage.coverage.scaling_factors.tsv"))?;

    // Assert
    assert_eq!(
        scaling_row_chromosomes(&rows),
        vec!["chr2".to_string(), "chr10".to_string(), "chr1".to_string()]
    );

    Ok(())
}

#[test]
fn blacklist_masking_changes_scaling_profile_and_excludes_zeroed_bins_from_global_mean()
-> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Start from the same simple fixture with one fragment [20,80), then blacklist [20,40).
    //
    // Per-stride average coverage with stride=20 and blacklist exclusion becomes:
    // - [0,20):   0
    // - [20,40):  0   because every covered base in this stride bin is blacklisted
    // - [40,60):  1   fully covered and unmasked
    // - [60,80):  1   fully covered and unmasked
    // - later bins: 0
    //
    // With bin-size=40 and stride=20, the triangular kernel is [1,2,1].
    // Hand-derived avg-overlapping-position-coverage:
    // - row 0: truncated [2,1] over [0,0]          -> 0
    // - row 1: [1,2,1] over [0,0,1]                -> 1/4
    // - row 2: [1,2,1] over [0,1,1]                -> 3/4
    // - row 3: [1,2,1] over [1,1,0]                -> 3/4
    // - row 4: [1,2,1] over [1,0,0]                -> 1/4
    // - later rows: 0
    //
    // The global mean ignores zeros, so:
    //   mean = (1/4 + 3/4 + 3/4 + 1/4) / 4 = 2 / 4 = 1/2
    //
    // Inverted scaling factors therefore become:
    // - 0       -> 0
    // - 1/4     -> (1/2) / (1/4) = 2
    // - 3/4     -> (1/2) / (3/4) = 2/3
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let blacklist_path = out_dir.path().join("blacklist.bed");
    std::fs::write(&blacklist_path, "chr1\t20\t40\n")?;

    let mut cfg = make_simple_coverage_weights_config(out_dir.path(), &bam.bam);
    cfg.blacklist = Some(vec![blacklist_path]);

    // Act
    run(&cfg)?;
    let rows = parse_scaling_rows(&out_dir.path().join("coverage.coverage.scaling_factors.tsv"))?;

    // Assert
    let expected_average_pos_coverage = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let expected_average_overlap = [
        0.0,
        1.0 / 4.0,
        3.0 / 4.0,
        3.0 / 4.0,
        1.0 / 4.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    let expected_scaling = [0.0, 2.0, 2.0 / 3.0, 2.0 / 3.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    assert_eq!(rows.len(), 10);
    for row_index in 0..rows.len() {
        let row = &rows[row_index];
        assert_approx(
            row.average_pos_coverage,
            expected_average_pos_coverage[row_index],
            1e-6,
            &format!("blacklist average_pos_coverage at row {row_index}"),
        );
        assert_approx(
            row.average_overlapping_pos_coverage,
            expected_average_overlap[row_index],
            1e-6,
            &format!("blacklist average_overlapping_pos_coverage at row {row_index}"),
        );
        assert_approx(
            row.scaling_factor,
            expected_scaling[row_index],
            1e-6,
            &format!("blacklist scaling_factor at row {row_index}"),
        );
    }

    Ok(())
}
