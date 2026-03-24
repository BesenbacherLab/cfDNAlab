#![cfg(feature = "cmd_coverage_weights")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs};
use cfdnalab::commands::coverage_weights::config::CoverageWeightsConfig;
use cfdnalab::commands::coverage_weights::coverage_weights::run;
use cfdnalab::commands::coverage_weights::striding::{
    StrideBin, normalize_avg_overlap_by_global_mean,
};
use cfdnalab::shared::interval::Interval;
use fxhash::FxHashMap;
use fixtures::{FragmentSpec, ReadSpec, bam_from_specs, paired_fragment, simple_inward_bam};
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
    // We add one 10 bp fragment on each chromosome. With `bin_size = stride = 20`, each
    // chromosome contributes exactly one TSV row, so the written chromosome sequence directly
    // exposes the command's row-order contract.
    bam_from_specs(
        vec![
            ("chr2".to_string(), 20),
            ("chr10".to_string(), 20),
            ("chr1".to_string(), 20),
        ],
        vec![
            fragment_on_tid(paired_fragment(0, 10, 60), 0),
            fragment_on_tid(paired_fragment(0, 10, 60), 1),
            fragment_on_tid(paired_fragment(0, 10, 60), 2),
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

#[test]
fn normalize_avg_overlap_keeps_sparse_non_zero_scaling_finite() -> Result<()> {
    // Arrange:
    // Use three equal-length stride bins on one chromosome:
    // - a very sparse but non-zero bin with avg-overlap coverage 0.0001
    // - a typical covered bin with avg-overlap coverage 1.0
    // - an uncovered bin with avg-overlap coverage 0.0
    //
    // The normalization ignores the zero bin when computing the global mean:
    //   mean = (0.0001 + 1.0) / 2 = 1.0001 / 2 = 0.50005
    //
    // With inversion enabled, scaling becomes:
    //   sparse bin   = 1 / (0.0001 / 0.50005) = 5000.5
    //   covered bin  = 1 / (1.0    / 0.50005) = 0.50005
    //   zero bin     = 0 by explicit zero-preserving logic
    //
    // This pins the important near-zero regime: the factor is huge but still finite.
    let mut bins_by_chr = FxHashMap::default();
    bins_by_chr.insert(
        "chr1".to_string(),
        vec![
            StrideBin {
                interval: Interval::new(0, 20)?,
                avg_coverage: 0.0001,
                avg_overlap_coverage: 0.0001,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(20, 40)?,
                avg_coverage: 1.0,
                avg_overlap_coverage: 1.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(40, 60)?,
                avg_coverage: 0.0,
                avg_overlap_coverage: 0.0,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_avg_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

    // Assert
    assert_approx(mean as f64, 0.50005, 1e-7, "global mean before inversion");
    let bins = bins_by_chr.get("chr1").expect("chr1 bins should remain present");
    assert!(
        bins[0].scaling_factor.is_finite(),
        "sparse non-zero bin should produce a large but finite scaling factor"
    );
    assert_approx(
        bins[0].scaling_factor as f64,
        5000.5,
        1e-3,
        "sparse-bin scaling factor",
    );
    assert_approx(
        bins[1].scaling_factor as f64,
        0.50005,
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
fn normalize_avg_overlap_weights_short_final_bin_in_global_mean() -> Result<()> {
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
                avg_coverage: 1.0,
                avg_overlap_coverage: 1.0,
                scaling_factor: 0.0,
            },
            StrideBin {
                interval: Interval::new(20, 30)?,
                avg_coverage: 3.0,
                avg_overlap_coverage: 3.0,
                scaling_factor: 0.0,
            },
        ],
    );

    // Act
    let mean = normalize_avg_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

    // Assert
    assert_approx(mean as f64, 5.0 / 3.0, 1e-6, "length-weighted global mean");
    let bins = bins_by_chr.get("chr1").expect("chr1 bins should remain present");
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
    let rows = parse_scaling_rows(&out_dir.path().join("coverage.scaling_factors.tsv"))?;

    // Assert
    let chr1_row_40 = rows
        .iter()
        .find(|row| row.chromosome == "chr1" && row.start == 40)
        .expect("missing chr1 row at start 40");
    assert_approx(
        chr1_row_40.avg_overlapping_pos_cov,
        1.0,
        1e-6,
        "chr1 central overlap coverage",
    );
    assert_approx(
        chr1_row_40.scaling_factor,
        25.0 / 48.0,
        1e-6,
        "chr1 central scaling factor from shared global mean",
    );

    let chr2_row_0 = rows
        .iter()
        .find(|row| row.chromosome == "chr2" && row.start == 0)
        .expect("missing chr2 row at start 0");
    assert_approx(
        chr2_row_0.avg_overlapping_pos_cov,
        1.0 / 3.0,
        1e-6,
        "chr2 edge overlap coverage",
    );
    assert_approx(
        chr2_row_0.scaling_factor,
        25.0 / 16.0,
        1e-6,
        "chr2 edge scaling factor from shared global mean",
    );

    let chr2_row_20 = rows
        .iter()
        .find(|row| row.chromosome == "chr2" && row.start == 20)
        .expect("missing chr2 row at start 20");
    assert_approx(
        chr2_row_20.avg_overlapping_pos_cov,
        1.0 / 2.0,
        1e-6,
        "chr2 central overlap coverage",
    );
    assert_approx(
        chr2_row_20.scaling_factor,
        25.0 / 24.0,
        1e-6,
        "chr2 central scaling factor from shared global mean",
    );

    Ok(())
}

#[test]
fn default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero() -> Result<()> {
    // Arrange:
    // Build one inward fragment whose mates both have MAPQ 20.
    //
    // `coverage-weights` intentionally defaults to `min_mapq = 30`, so:
    // - default config must drop the fragment completely
    // - explicit `min_mapq = 30` must match exactly
    // - explicit `min_mapq = 0` must keep the fragment
    //
    // With `bin_size = 40` and `stride = 20`, keeping the fragment would produce the same
    // hand-derived overlap profile as the standard simple fixture:
    //   avg-overlap = [1/3, 3/4, 1, 3/4, 1/4, 0, 0, 0, 0, 0]
    //   scaling     = [37/20, 37/45, 37/60, 37/45, 37/15, 0, 0, 0, 0, 0]
    //
    // Dropping the fragment yields zero coverage everywhere, so every row in the TSV should be 0.
    let mut low_mapq = paired_fragment(20, 60, 20);
    low_mapq.forward.mapq = 20;
    low_mapq.reverse.mapq = 20;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![low_mapq],
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
    let default_rows = parse_scaling_rows(&default_out.path().join("default.scaling_factors.tsv"))?;
    let explicit_thirty_rows =
        parse_scaling_rows(&explicit_thirty_out.path().join("thirty.scaling_factors.tsv"))?;
    let explicit_zero_rows =
        parse_scaling_rows(&explicit_zero_out.path().join("zero.scaling_factors.tsv"))?;

    assert_eq!(default_rows.len(), 10);
    assert_eq!(default_rows.len(), explicit_thirty_rows.len());
    assert_eq!(default_rows.len(), explicit_zero_rows.len());

    for (default_row, explicit_row) in default_rows.iter().zip(explicit_thirty_rows.iter()) {
        assert_eq!(default_row.chromosome, explicit_row.chromosome);
        assert_eq!(default_row.start, explicit_row.start);
        assert_eq!(default_row.end, explicit_row.end);
        assert_approx(default_row.avg_pos_cov, explicit_row.avg_pos_cov, 1e-12, "default vs explicit-30 avg_pos_cov");
        assert_approx(
            default_row.avg_overlapping_pos_cov,
            explicit_row.avg_overlapping_pos_cov,
            1e-12,
            "default vs explicit-30 avg_overlap",
        );
        assert_approx(
            default_row.scaling_factor,
            explicit_row.scaling_factor,
            1e-12,
            "default vs explicit-30 scaling",
        );
        assert_eq!(default_row.avg_pos_cov, 0.0);
        assert_eq!(default_row.avg_overlapping_pos_cov, 0.0);
        assert_eq!(default_row.scaling_factor, 0.0);
    }

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
    for (row_index, row) in explicit_zero_rows.iter().enumerate() {
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
        frag.min_fragment_length = 1;
        frag.max_fragment_length = 200;
    }

    // Act
    run(&cfg)?;
    let rows = parse_scaling_rows(&out_dir.path().join("coverage.scaling_factors.tsv"))?;

    // Assert
    assert_eq!(
        scaling_row_chromosomes(&rows),
        vec!["chr10".to_string(), "chr2".to_string()]
    );

    Ok(())
}

#[test]
fn chromosomes_all_uses_bam_header_order_for_scaling_tsv_rows() -> Result<()> {
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
        frag.min_fragment_length = 1;
        frag.max_fragment_length = 200;
    }

    // Act
    run(&cfg)?;
    let rows = parse_scaling_rows(&out_dir.path().join("coverage.scaling_factors.tsv"))?;

    // Assert
    assert_eq!(
        scaling_row_chromosomes(&rows),
        vec![
            "chr2".to_string(),
            "chr10".to_string(),
            "chr1".to_string()
        ]
    );

    Ok(())
}

#[test]
fn blacklist_masking_changes_scaling_profile_and_excludes_zeroed_bins_from_global_mean()
-> Result<()> {
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
    let rows = parse_scaling_rows(&out_dir.path().join("coverage.scaling_factors.tsv"))?;

    // Assert
    let expected_avg_pos_cov = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let expected_avg_overlap = [
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
    let expected_scaling = [
        0.0,
        2.0,
        2.0 / 3.0,
        2.0 / 3.0,
        2.0,
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
            &format!("blacklist avg_pos_cov at row {row_index}"),
        );
        assert_approx(
            row.avg_overlapping_pos_cov,
            expected_avg_overlap[row_index],
            1e-6,
            &format!("blacklist avg_overlapping_pos_cov at row {row_index}"),
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
