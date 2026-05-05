#![cfg(feature = "cmd_fcoverage")]

mod fixtures;

use anyhow::Result;
#[cfg(feature = "cmd_bam_to_bam")]
use cfdnalab::commands::bam_to_bam::{
    bam_to_bam::run_inner as run_bam_to_bam, config::BamToBamConfig,
};
use cfdnalab::commands::cli_common::{
    ApplyGCArgs, AssignToWindowArgs, ChromosomeArgs, DistributionWindowsArgs, IOCArgs,
    ScaleGenomeArgs,
};
#[cfg(feature = "cmd_coverage_weights")]
use cfdnalab::commands::coverage_weights::{
    config::CoverageWeightsConfig, coverage_weights::run as run_coverage_weights,
};
use cfdnalab::commands::fcoverage::config::{FCoverageConfig, LengthNormalizationMode};
use cfdnalab::commands::fcoverage::fcoverage::{run, run_inner};
use cfdnalab::commands::fcoverage::window_results::CoverageWindowAction;
use cfdnalab::commands::gc_bias::package::GCCorrectionPackage;
use cfdnalab::commands::lengths::{config::LengthsConfig, lengths::run as run_lengths};
use cfdnalab::shared::fragment::minimal_fragment::{MinimalReadInfo, collect_fragment};
use cfdnalab::shared::indel_mode::IndelMode;
use cfdnalab::shared::io::dot_join;
use cfdnalab::shared::read::default_include_read_paired_end;
use cfdnalab::shared::{
    constants::GC_CORRECTION_SCHEMA_VERSION,
    reference::{ContigFootprintEntry, twobit_contig_footprint},
};
use fixtures::{
    BamFixture, FragmentSpec, LONG_FRAGMENT_LENGTH, LONG_FRAGMENT_STARTS, ReadSpec, bam_from_specs,
    bam_from_specs_strict_identity, build_real_neutral_gc_package,
    build_real_neutral_gc_package_for_range, build_real_non_neutral_gc_package,
    late_origin_gc_reference_sequence, long_fragment_bam, paired_fragment, read_zst_to_string,
    simple_inward_bam, simple_reference_twobit, twobit_from_sequences, write_bed,
    write_scaling_factors, write_two_bin_gc_package,
};
use ndarray::Array2;
use ndarray::array;
use ndarray_npy::read_npy;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read, Reader};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

fn collect_fragment_from_records_for_test(
    a: &rust_htslib::bam::Record,
    b: &rust_htslib::bam::Record,
) -> Option<cfdnalab::shared::fragment::minimal_fragment::Fragment> {
    collect_fragment(
        &MinimalReadInfo::from_record_with_gc_tag(a, None).ok()?,
        &MinimalReadInfo::from_record_with_gc_tag(b, None).ok()?,
    )
}

#[derive(Debug)]
struct TaggedBamFixture {
    _tempdir: TempDir,
    bam: PathBuf,
}

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
        chromosomes_file: None,
    }
}

fn base_config(bam_path: &Path, output_dir: &Path) -> FCoverageConfig {
    let mut cfg = FCoverageConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_output_prefix("testcov");
    cfg.set_tile_size(1_000);
    cfg.set_ignore_gap(false);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }
    cfg
}

#[test]
fn fcoverage_rejects_inverted_fragment_length_range_before_reading_inputs() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let mut cfg = base_config(Path::new("missing.bam"), temp_dir.path());
    {
        let fragment_lengths = cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 500;
        fragment_lengths.max_fragment_length = 100;
    }

    let error = match run_inner(&cfg) {
        Ok(_) => panic!("inverted fragment length range should fail"),
        Err(error) => error,
    };
    let message = error.to_string();

    assert!(
        message.contains("--min-fragment-length (500) must be <= --max-fragment-length (100)"),
        "unexpected error: {message}"
    );

    Ok(())
}

#[test]
fn fcoverage_rejects_gc_file_without_ref_2bit_before_output_setup() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let output_dir = temp_dir.path().join("not_created");
    let mut cfg = base_config(Path::new("missing.bam"), &output_dir);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(temp_dir.path().join("missing_gc_package.npz")),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });

    let error = match run_inner(&cfg) {
        Ok(_) => panic!("GC file without --ref-2bit should fail"),
        Err(error) => error,
    };
    let message = error.to_string();

    assert!(
        message.contains("--gc-file requires --ref-2bit"),
        "unexpected error: {message}"
    );
    assert!(
        !output_dir.exists(),
        "configuration validation should fail before creating the output directory"
    );

    Ok(())
}

fn set_restore_mean_length_normalization(cfg: &mut FCoverageConfig) {
    cfg.set_normalize_by_length_mode(LengthNormalizationMode::RestoreMean);
}

fn mixed_length_fragment_bam() -> Result<BamFixture> {
    bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment(20, 40, 20), paired_fragment(100, 80, 20)],
        Vec::new(),
        "fcoverage_restore_mean_mixed_lengths",
    )
}

fn mixed_length_three_chromosome_bam() -> Result<BamFixture> {
    bam_from_specs(
        vec![
            ("chr1".to_string(), 200),
            ("chr2".to_string(), 200),
            ("chr3".to_string(), 200),
        ],
        vec![
            paired_fragment_on_tid(0, 20, 40, 20),
            paired_fragment_on_tid(1, 10, 60, 20),
            paired_fragment_on_tid(2, 40, 80, 20),
        ],
        Vec::new(),
        "fcoverage_restore_mean_three_chr",
    )
}

fn empty_bam_fixture(name: &str) -> Result<BamFixture> {
    bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        Vec::new(),
        name,
    )
}

fn dense_bedgraph_for_chromosome(
    text: &str,
    chromosome: &str,
    chromosome_length: usize,
) -> Vec<f64> {
    let mut coverage = vec![0.0; chromosome_length];
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let cols: Vec<_> = line.split('\t').collect();
        assert!(
            cols.len() >= 4,
            "expected at least four bedGraph columns, got '{}'",
            line
        );
        if cols[0] != chromosome {
            continue;
        }
        let start: usize = cols[1].parse().expect("bedGraph start must parse");
        let end: usize = cols[2].parse().expect("bedGraph end must parse");
        let value: f64 = cols[3].parse().expect("bedGraph value must parse");
        for position in start..end {
            coverage[position] = value;
        }
    }
    coverage
}

fn parse_tsv(text: &str) -> Vec<Vec<&str>> {
    text.lines()
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect()
}

fn assert_close(observed: f64, expected: f64, tolerance: f64) {
    assert!(
        (observed - expected).abs() <= tolerance,
        "expected {expected}, got {observed}"
    );
}

fn assert_close_to_written_precision(observed: f64, expected: f64, decimals: i32) {
    let tolerance = 0.5 * 10_f64.powi(-decimals) + 2e-6;
    assert_close(observed, expected, tolerance);
}

fn pearson_r_from_vectors(coverage: &[f64], mask: &[f64]) -> f64 {
    assert_eq!(
        coverage.len(),
        mask.len(),
        "coverage and mask must have the same length"
    );
    let n = coverage.len();
    assert!(n > 0, "cannot compute Pearson R on empty vectors");

    let average_coverage = coverage.iter().sum::<f64>() / n as f64;
    let mean_mask = mask.iter().sum::<f64>() / n as f64;

    let mut covariance_numerator = 0.0;
    let mut coverage_variance_numerator = 0.0;
    let mut mask_variance_numerator = 0.0;

    for (&coverage_value, &mask_value) in coverage.iter().zip(mask.iter()) {
        let centered_coverage = coverage_value - average_coverage;
        let centered_mask = mask_value - mean_mask;
        covariance_numerator += centered_coverage * centered_mask;
        coverage_variance_numerator += centered_coverage * centered_coverage;
        mask_variance_numerator += centered_mask * centered_mask;
    }

    let denominator = (coverage_variance_numerator * mask_variance_numerator).sqrt();
    if denominator == 0.0 {
        f64::NAN
    } else {
        covariance_numerator / denominator
    }
}

fn pearson_r_from_summary_stats_rows(global_row: &[String], group_row: &[String]) -> Result<f64> {
    let global_eligible_positions = global_row[3].parse::<f64>()?;
    let global_coverage_sum = global_row[5].parse::<f64>()?;
    let global_coverage_sum_of_squares = global_row[6].parse::<f64>()?;
    let group_eligible_positions = group_row[3].parse::<f64>()?;
    let group_coverage_sum = group_row[5].parse::<f64>()?;

    let coverage_term =
        global_eligible_positions * global_coverage_sum_of_squares - global_coverage_sum.powi(2);
    let mask_term =
        global_eligible_positions * group_eligible_positions - group_eligible_positions.powi(2);
    let denominator = (coverage_term * mask_term).sqrt();

    if denominator == 0.0 {
        Ok(f64::NAN)
    } else {
        Ok((global_eligible_positions * group_coverage_sum
            - global_coverage_sum * group_eligible_positions)
            / denominator)
    }
}

fn grouped_rows_by_name(
    output_text: &str,
    sidecar_text: &str,
) -> std::collections::BTreeMap<String, Vec<String>> {
    let output_rows = parse_tsv(output_text);
    let sidecar_rows = parse_tsv(sidecar_text);

    let mut idx_to_name = std::collections::BTreeMap::new();
    for row in sidecar_rows.into_iter().skip(1) {
        idx_to_name.insert(row[0].to_string(), row[1].to_string());
    }

    let mut rows_by_name = std::collections::BTreeMap::new();
    for row in output_rows.into_iter().skip(1) {
        let group_name = idx_to_name
            .get(row[0])
            .expect("group_idx in output must exist in sidecar")
            .clone();
        rows_by_name.insert(
            group_name,
            row.into_iter().map(|value| value.to_string()).collect(),
        );
    }

    rows_by_name
}

fn overlapping_fragment_bam() -> Result<fixtures::BamFixture> {
    bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment(20, 60, 20), paired_fragment(30, 40, 20)],
        Vec::new(),
        "overlapping_fcoverage",
    )
}

fn single_read_fragment_bam(name: &str) -> Result<BamFixture> {
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

fn single_read_fragment_bam_at(
    name: &str,
    fragment_start: i64,
    fragment_len: u32,
) -> Result<BamFixture> {
    bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: fragment_start,
            cigar: vec![('M', fragment_len)],
            seq: vec![b'A'; fragment_len as usize],
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

fn build_bai_for_test_bam(bam_path: &Path) -> Result<()> {
    let bai_path = bam_path.with_extension("bam.bai");
    bam::index::build(bam_path, None, bam::index::Type::Bai, 1)?;
    let target = bam_path.with_extension("bai");
    if bai_path.exists() {
        std::fs::rename(&bai_path, &target)?;
    }
    Ok(())
}

fn bam_with_gc_tags(base_bam: &Path, name: &str, tags: &[Option<f32>]) -> Result<TaggedBamFixture> {
    let tempdir = TempDir::new()?;
    let bam_path = tempdir.path().join(format!("{name}.bam"));

    let mut reader = Reader::from_path(base_bam)?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam)?;

    for (record_index, record_result) in reader.records().enumerate() {
        let mut record = record_result?;
        if let Some(Some(tag_value)) = tags.get(record_index) {
            record.push_aux(b"GC", Aux::Float(*tag_value))?;
        }
        writer.write(&record)?;
    }

    drop(writer);
    build_bai_for_test_bam(&bam_path)?;

    Ok(TaggedBamFixture {
        _tempdir: tempdir,
        bam: bam_path,
    })
}

fn build_gc_package(
    path: &Path,
    end_offset: u64,
    reference_contig_footprint: Vec<ContigFootprintEntry>,
) -> Result<()> {
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset,
        length_edges: vec![10, 60, 200],
        gc_edges: vec![0, 50, 101],
        length_bin_frequencies: array![1.0_f64, 3.0_f64],
        reference_contig_footprint,
        correction_matrix: array![[1.0_f64, 1.0_f64], [2.0_f64, 10.0_f64]],
    };
    package.write_npz(path)?;
    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
fn make_simple_coverage_weights_config(
    out_dir: &std::path::Path,
    bam: &std::path::Path,
) -> CoverageWeightsConfig {
    let mut cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam.to_path_buf(),
            output_dir: out_dir.to_path_buf(),
            n_threads: 1,
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

fn paired_fragment_on_tid(
    tid: usize,
    start: i64,
    fragment_len: i64,
    read_len: i64,
) -> FragmentSpec {
    const FLAG_FIRST_MATE: u16 = 0x40;
    const FLAG_SECOND_MATE: u16 = 0x80;
    const FLAG_PROPER_PAIR: u16 = 0x2;
    const FLAG_MATE_REVERSE: u16 = 0x20;

    let reverse_start = start + fragment_len - read_len;
    let insert_size = fragment_len;
    FragmentSpec {
        forward: ReadSpec {
            tid,
            pos: start,
            cigar: vec![('M', read_len as u32)],
            seq: vec![b'A'; read_len as usize],
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
            mate_tid: Some(tid),
            mate_pos: Some(reverse_start),
            insert_size,
        },
        reverse: ReadSpec {
            tid,
            pos: reverse_start,
            cigar: vec![('M', read_len as u32)],
            seq: vec![b'T'; read_len as usize],
            qual: 40,
            is_reverse: true,
            mapq: 60,
            flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
            mate_tid: Some(tid),
            mate_pos: Some(start),
            insert_size: -insert_size,
        },
    }
}

#[test]
fn per_position_outputs_basic_fragment() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);

    let mut reader = Reader::from_path(&bam.bam)?;
    let mut accepted = (0, 0);
    let mut pair_store = Vec::new();
    for (idx, result) in reader.records().enumerate() {
        let rec = result?;
        if default_include_read_paired_end(&rec, cfg.require_proper_pair, cfg.min_mapq) {
            if rec.is_reverse() {
                accepted.1 += 1;
            } else {
                accepted.0 += 1;
            }
            pair_store.push(rec);
        } else if idx == 0 {
            eprintln!(
                "forward read filtered: flags={:#x}, is_reverse={}, mate_reverse={}, mapq={}",
                rec.flags(),
                rec.is_reverse(),
                rec.is_mate_reverse(),
                rec.mapq()
            );
        }
    }
    assert_eq!(
        accepted,
        (1, 1),
        "expected both mates accepted, got {accepted:?}"
    );
    assert_eq!(pair_store.len(), 2);
    let frag = collect_fragment_from_records_for_test(&pair_store[0], &pair_store[1]);
    assert!(frag.is_some(), "expected fragment collection to succeed");

    run(&cfg)?;

    let bedgraph = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    assert!(bedgraph.exists(), "expected positional bedgraph output");
    let text = read_zst_to_string(&bedgraph)?;
    assert_eq!(text, "chr1\t20\t80\t1\n");

    Ok(())
}

#[test]
fn normalize_by_length_keeps_fractional_positional_output_without_other_weights() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(4);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);

    // Manual expectations:
    // - The single fragment covers [20, 80), length 60.
    // - With `--normalize-by-length`, each covered base gets weight 1 / 60.
    // - No scaling or GC correction is applied, so the whole covered run stays constant.
    // - Rounded to 4 decimals:
    //     1 / 60 = 0.016666... -> 0.0167
    run(&cfg)?;

    let bedgraph = out_dir.path().join(dot_join(&[
        "testcov",
        "length_normalized",
        "fcoverage.per_position.bedgraph.zst",
    ]));
    assert!(bedgraph.exists(), "expected positional bedgraph output");
    let text = read_zst_to_string(&bedgraph)?;
    assert_eq!(text, "chr1\t20\t80\t0.0167\n");

    Ok(())
}

#[test]
fn normalize_by_length_by_size_total_counts_each_fragment_as_one() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_size = Some(200);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(5);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(windows);
    cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);

    // Manual expectations:
    // - The fragment spans 60 bases and contributes 1 / 60 to each covered base.
    // - The only size window is [0, 200), which contains the full fragment.
    // - Total coverage in that window is therefore:
    //     60 * (1 / 60) = 1
    // - This is the intended fragment-count interpretation of the total output.
    run(&cfg)?;

    let output_path = out_dir.path().join(dot_join(&[
        "testcov",
        "length_normalized",
        "fcoverage.total.tsv.zst",
    ]));
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t200\t1\t0",
        ]
    );

    Ok(())
}

#[test]
fn normalize_by_length_and_gc_file_weights_multiply_per_position() -> Result<()> {
    // Human verification status: unverified
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment_on_tid(0, 20, 61, 20)],
        Vec::new(),
        "fcoverage_normalize_by_length_gc_file",
    )?;
    let ref_twobit = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("constant_gc_pkg.npz");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![61, 62],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
        correction_matrix: array![[3.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 61;
        frag.max_fragment_length = 61;
    }

    // Manual expectations:
    // - The fragment spans [20, 81), length 61.
    // - `--normalize-by-length` gives each covered base weight 1 / 61.
    // - The constant GC package multiplies every accepted fragment by 3.0.
    // - Final per-base coverage is therefore 3 / 61 across the full span.
    run(&cfg)?;

    let output_path = out_dir.path().join(dot_join(&[
        "testcov",
        "length_normalized",
        "fcoverage.per_position.bedgraph.zst",
    ]));
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines.len(), 1, "expected one constant bedGraph run");

    let parts: Vec<_> = lines[0].split('\t').collect();
    assert_eq!(parts.len(), 4, "unexpected bedGraph row: {}", lines[0]);
    assert_eq!(parts[0], "chr1");
    assert_eq!(parts[1].parse::<u64>()?, 20);
    assert_eq!(parts[2].parse::<u64>()?, 81);
    let value = parts[3].parse::<f64>()?;
    let expected = 3.0_f64 / 61.0_f64;
    assert!(
        (value - expected).abs() <= 1e-6,
        "expected value {expected} for row {}, got {value}",
        lines[0]
    );

    Ok(())
}

#[test]
fn gc_file_windowed_late_tile_uses_reference_coordinates_after_fetch_narrowing() -> Result<()> {
    // Arrange:
    // - The tile core starts at 0, but the only BED window is [930,941).
    // - The reference is shorter than the BAM chromosome, but long enough for the narrowed
    //   window-derived fetch span. Reading the full tile reference would overrun the reference.
    // - The accepted fragment interval [900,961) is all C, so it lands in the high-GC correction
    //   bin with weight 7.0. Using prefix-local origin 0 would see A-only sequence instead.
    // - In unique-position mode, only the 11 covered bases inside [930,941) are written, each at
    //   coverage 7.0.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 1_500)],
        vec![paired_fragment_on_tid(0, 900, 61, 20)],
        Vec::new(),
        "fcoverage_late_tile_gc_origin",
    )?;
    let reference = twobit_from_sequences(
        "fcoverage_late_tile_gc_origin_ref",
        vec![("chr1".to_string(), late_origin_gc_reference_sequence())],
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("late_window.bed");
    let gc_path = out_dir.path().join("two_bin_gc_package.npz");
    write_bed(&bed_path, &[("chr1", 930, 941, "late")])?;
    write_two_bin_gc_package(
        &gc_path,
        61,
        2.0,
        7.0,
        twobit_contig_footprint(&reference.path)?,
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_keep_zero_runs(false);
    cfg.set_per_window(CoverageWindowAction::OnlyIncludeThesePositionsUnique);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 61;
        frag.max_fragment_length = 61;
    }

    // Act
    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    // Assert
    assert_eq!(text, "chr1\t930\t941\t7\n");
    Ok(())
}

#[test]
fn normalize_by_length_uses_counted_segment_length_for_gapped_fragments() -> Result<()> {
    // Human verification status: unverified
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: 20,
            cigar: vec![('M', 20), ('D', 10), ('M', 20)],
            seq: vec![b'A'; 40],
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
        "fcoverage_normalize_by_length_gapped_unpaired",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(4);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.unpaired.reads_are_fragments = true;
    cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);

    // Manual expectations:
    // - The read spans [20, 70) on the reference with a 10 bp deletion in the middle.
    // - `fcoverage` respects reference-supported segments, so the counted spans are:
    //     [20, 40) and [50, 70)
    // - Their total counted length is therefore 40 bp, not the outer span length 50 bp.
    // - With `--normalize-by-length`, each counted base must get weight:
    //     1 / 40 = 0.025
    // - If the implementation incorrectly divided by the outer span length, this would instead be
    //   1 / 50 = 0.02.
    run(&cfg)?;

    let output_path = out_dir.path().join(dot_join(&[
        "testcov",
        "length_normalized",
        "fcoverage.per_position.bedgraph.zst",
    ]));
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t40\t0.025", "chr1\t50\t70\t0.025"]);

    Ok(())
}

#[test]
fn normalize_by_length_ignore_gap_renormalizes_over_remaining_counted_span() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;

    let run_with_ignore_gap = |ignore_gap: bool| -> Result<String> {
        let out_dir = TempDir::new()?;
        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(4);
        cfg.set_per_window(CoverageWindowAction::Average);
        cfg.set_keep_zero_runs(false);
        cfg.set_ignore_gap(ignore_gap);
        cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);

        run(&cfg)?;

        let output_path = out_dir.path().join(dot_join(&[
            "testcov",
            "length_normalized",
            "fcoverage.per_position.bedgraph.zst",
        ]));
        read_zst_to_string(&output_path)
    };

    // Manual expectations:
    // - The paired fragment has outer span [20, 80), length 60.
    // - Without `--ignore-gap`, that whole span is counted, so each base gets:
    //     1 / 60 = 0.016666... -> 0.0167
    // - With `--ignore-gap`, only the read-covered segments remain:
    //     [20, 40) and [60, 80)
    //   Their total counted length is 40, so each counted base gets:
    //     1 / 40 = 0.025
    let without_ignore_gap = run_with_ignore_gap(false)?;
    let with_ignore_gap = run_with_ignore_gap(true)?;

    assert_eq!(without_ignore_gap, "chr1\t20\t80\t0.0167\n");
    assert_eq!(
        with_ignore_gap,
        "chr1\t20\t40\t0.025\nchr1\t60\t80\t0.025\n"
    );

    Ok(())
}

#[test]
fn normalize_by_length_matches_between_paired_and_unpaired_for_same_span() -> Result<()> {
    // Human verification status: unverified
    let paired_bam = simple_inward_bam()?;
    let unpaired_bam = single_read_fragment_bam("fcoverage_normalize_by_length_unpaired_parity")?;
    let paired_out = TempDir::new()?;
    let unpaired_out = TempDir::new()?;

    let mut paired_cfg = base_config(&paired_bam.bam, paired_out.path());
    paired_cfg.set_decimals(4);
    paired_cfg.set_per_window(CoverageWindowAction::Average);
    paired_cfg.set_keep_zero_runs(false);
    paired_cfg.set_output_prefix("paired");
    paired_cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);

    let mut unpaired_cfg = base_config(&unpaired_bam.bam, unpaired_out.path());
    unpaired_cfg.set_decimals(4);
    unpaired_cfg.set_per_window(CoverageWindowAction::Average);
    unpaired_cfg.set_keep_zero_runs(false);
    unpaired_cfg.set_output_prefix("unpaired");
    unpaired_cfg.unpaired.reads_are_fragments = true;
    unpaired_cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);

    // Manual expectations:
    // - Both inputs represent the same counted span [20, 80), length 60.
    // - With `--normalize-by-length`, each counted base gets:
    //     1 / 60 = 0.016666... -> 0.0167
    // - Because there is no mate-gap exclusion or segment logic here, paired and unpaired modes
    //   must therefore write identical normalized coverage over the same interval.
    run(&paired_cfg)?;
    run(&unpaired_cfg)?;

    let paired_text = read_zst_to_string(&paired_out.path().join(dot_join(&[
        "paired",
        "length_normalized",
        "fcoverage.per_position.bedgraph.zst",
    ])))?;
    let unpaired_text = read_zst_to_string(&unpaired_out.path().join(dot_join(&[
        "unpaired",
        "length_normalized",
        "fcoverage.per_position.bedgraph.zst",
    ])))?;
    let expected = "chr1\t20\t80\t0.0167\n";
    assert_eq!(paired_text, expected);
    assert_eq!(unpaired_text, expected);

    Ok(())
}

#[test]
fn normalize_by_length_uses_paired_counted_segments_when_ignore_gap_is_on() -> Result<()> {
    // Human verification status: unverified
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![FragmentSpec {
            forward: ReadSpec {
                tid: 0,
                pos: 20,
                cigar: vec![('M', 20), ('D', 10), ('M', 20)],
                seq: vec![b'A'; 40],
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0x40 | 0x20 | 0x2,
                mate_tid: Some(0),
                mate_pos: Some(60),
                insert_size: 60,
            },
            reverse: ReadSpec {
                tid: 0,
                pos: 60,
                cigar: vec![('M', 20)],
                seq: vec![b'T'; 20],
                qual: 40,
                is_reverse: true,
                mapq: 60,
                flags: 0x80 | 0x2,
                mate_tid: Some(0),
                mate_pos: Some(20),
                insert_size: -60,
            },
        }],
        Vec::new(),
        "fcoverage_normalize_by_length_paired_gapped_ignore_gap",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(4);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.set_ignore_gap(true);
    cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);

    // Manual expectations:
    // - The fragment outer span is [20, 80), length 60.
    // - The forward read has a 10 bp deletion, so its counted reference-supported segments are:
    //     [20, 40) and [50, 70)
    // - The reverse read contributes [60, 80).
    // - With `--ignore-gap`, the inter-mate gap [40, 60) is not added.
    // - Merging overlapping counted spans gives:
    //     [20, 40) and [50, 80)
    //   with total counted length 20 + 30 = 50.
    // - With `--normalize-by-length`, each counted base must therefore get:
    //     1 / 50 = 0.02
    run(&cfg)?;

    let output_path = out_dir.path().join(dot_join(&[
        "testcov",
        "length_normalized",
        "fcoverage.per_position.bedgraph.zst",
    ]));
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t40\t0.02", "chr1\t50\t80\t0.02"]);

    Ok(())
}

#[test]
fn normalize_by_length_uses_counted_segment_length_for_refskip_fragments() -> Result<()> {
    // Human verification status: unverified
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: 20,
            cigar: vec![('M', 20), ('N', 10), ('M', 20)],
            seq: vec![b'A'; 40],
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
        "fcoverage_normalize_by_length_refskip_unpaired",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(4);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.unpaired.reads_are_fragments = true;
    cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);

    // Manual expectations:
    // - The read spans [20, 70) on the reference with a 10 bp ref-skip in the middle.
    // - `fcoverage` treats that like counted reference segments:
    //     [20, 40) and [50, 70)
    // - Their total counted length is therefore 40 bp.
    // - With `--normalize-by-length`, each counted base gets:
    //     1 / 40 = 0.025
    run(&cfg)?;

    let output_path = out_dir.path().join(dot_join(&[
        "testcov",
        "length_normalized",
        "fcoverage.per_position.bedgraph.zst",
    ]));
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t40\t0.025", "chr1\t50\t70\t0.025"]);

    Ok(())
}

#[test]
fn normalize_by_length_segmented_fragment_still_multiplies_gc_and_scaling() -> Result<()> {
    // Human verification status: unverified
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: 20,
            cigar: vec![('M', 20), ('D', 10), ('M', 20)],
            seq: vec![b'A'; 40],
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
        "fcoverage_normalize_by_length_gapped_with_gc_and_scaling",
    )?;
    let ref_twobit = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("constant_gc_pkg.npz");
    let scaling_path = out_dir.path().join("scaling.tsv");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![50, 51],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
        correction_matrix: array![[2.0_f64]],
    };
    package.write_npz(&gc_path)?;
    write_scaling_factors(&scaling_path, &[("chr1", 0, 200, 5.0)])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(4);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.unpaired.reads_are_fragments = true;
    cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let mut scale_genome = ScaleGenomeArgs::default();
        scale_genome.scaling_factors = Some(scaling_path);
        cfg.set_scale_genome(scale_genome);
    }
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 50;
        frag.max_fragment_length = 50;
    }

    // Manual expectations:
    // - The read has counted segments [20, 40) and [50, 70), total counted length 40.
    // - Intrinsic normalized base weight is therefore:
    //     1 / 40 = 0.025
    // - The GC package contributes a fragment weight of 2.0.
    // - The scaling TSV contributes a per-base factor of 5.0 everywhere.
    // - Final per-base coverage is therefore:
    //     (1 / 40) * 2 * 5 = 0.25
    run(&cfg)?;

    let output_path = out_dir.path().join(dot_join(&[
        "testcov",
        "length_normalized",
        "fcoverage.per_position.bedgraph.zst",
    ]));
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t40\t0.25", "chr1\t50\t70\t0.25"]);

    Ok(())
}

#[test]
fn restore_mean_positional_output_uses_mean_normalization_length_for_mixed_fragment_lengths()
-> Result<()> {
    // Future-facing test for `--normalize-by-length=restore-mean`.
    //
    // Manual expectations:
    // - Two fragments are counted:
    //     [20, 60) with normalization length 40
    //     [100, 180) with normalization length 80
    // - Mean normalization length is therefore:
    //     (40 + 80) / 2 = 60
    // - Unit-mass per-base weights would be:
    //     first fragment:  1 / 40
    //     second fragment: 1 / 80
    // - `restore-mean` multiplies those by 60, giving:
    //     first fragment:  60 / 40 = 1.5
    //     second fragment: 60 / 80 = 0.75
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    assert!(
        result
            .final_out_path
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.contains("length_normalized.restored_mean")),
        "restore-mean output path should include the restore-mean marker, got {}",
        result.final_out_path.display()
    );
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(text, "chr1\t20\t60\t1.5\nchr1\t100\t180\t0.75\n");

    Ok(())
}

#[test]
fn restore_mean_positional_output_is_basewise_tile_size_invariant_across_tile_boundaries()
-> Result<()> {
    // Future-facing test for the positional scaled-merge path.
    //
    // Manual expectations:
    // - Same mixed-length fixture as above:
    //     [20, 60)  -> 1.5
    //     [100, 180) -> 0.75
    // - The final bedGraph segmentation may differ across tile sizes, but the basewise coverage
    //   must stay identical.
    let bam = mixed_length_fragment_bam()?;
    let mut observed_coverages = Vec::new();

    for tile_size in [33_u32, 1_000_u32] {
        let out_dir = TempDir::new()?;
        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(2);
        cfg.set_per_window(CoverageWindowAction::Average);
        cfg.set_keep_zero_runs(true);
        cfg.set_tile_size(tile_size);
        set_restore_mean_length_normalization(&mut cfg);

        let result = run_inner(&cfg)?;
        let text = read_zst_to_string(&result.final_out_path)?;
        observed_coverages.push(dense_bedgraph_for_chromosome(&text, "chr1", 200));
    }

    let mut expected = vec![0.0_f64; 200];
    for position in 20..60 {
        expected[position] = 1.5;
    }
    for position in 100..180 {
        expected[position] = 0.75;
    }

    assert_eq!(observed_coverages[0], expected);
    assert_eq!(observed_coverages[1], expected);

    Ok(())
}

#[test]
fn restore_mean_unique_positions_windowed_output_scales_selected_runs() -> Result<()> {
    // Future-facing test for `OnlyIncludeThesePositionsUnique`.
    //
    // Manual expectations:
    // - Mean normalization length is still 60.
    // - BED windows:
    //     [15, 50)  intersects the first fragment as [20, 50) with value 1.5
    //     [120, 190) intersects the second fragment as [120, 180) with value 0.75
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("restore_mean_unique_windows.bed");
    write_bed(
        &bed_path,
        &[("chr1", 15, 50, "first"), ("chr1", 120, 190, "second")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_keep_zero_runs(false);
    cfg.set_per_window(CoverageWindowAction::OnlyIncludeThesePositionsUnique);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });
    cfg.set_decimals(2);
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(text, "chr1\t20\t50\t1.5\nchr1\t120\t180\t0.75\n");

    Ok(())
}

#[test]
fn restore_mean_indexed_positions_keep_window_indices_and_scaled_values() -> Result<()> {
    // Future-facing test for `OnlyIncludeThesePositionsIndexed`.
    //
    // Manual expectations:
    // - Mean normalization length is 60.
    // - Window 0 [15, 50) contributes:
    //     [20, 50) -> 1.5
    // - Window 1 [40, 130) contributes:
    //     [40, 60) -> 1.5
    //     [100, 130) -> 0.75
    // - Window 2 [120, 190) contributes:
    //     [120, 180) -> 0.75
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("restore_mean_indexed_windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 15, 50, "window_a"),
            ("chr1", 40, 130, "window_b"),
            ("chr1", 120, 190, "window_c"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_keep_zero_runs(false);
    cfg.set_per_window(CoverageWindowAction::OnlyIncludeThesePositionsIndexed);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });
    cfg.set_decimals(2);
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chr1\t20\t50\t1.5\t0",
            "chr1\t40\t60\t1.5\t1",
            "chr1\t100\t130\t0.75\t1",
            "chr1\t120\t180\t0.75\t2",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_by_size_total_counts_mean_normalization_length_per_fragment() -> Result<()> {
    // Future-facing test for by-size total aggregation.
    //
    // Manual expectations:
    // - Mean normalization length = 60.
    // - Each fragment contributes total mass 60 after `restore-mean`.
    // - The single 200 bp window therefore gets:
    //     60 + 60 = 120
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: Some(200),
        by_bed: None,
        by_grouped_bed: None,
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t200\t120\t0",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_by_size_average_restores_global_mean_level_for_single_full_window() -> Result<()> {
    // Future-facing test for by-size average aggregation.
    //
    // Manual expectations:
    // - Restored total coverage across the chromosome is 120.
    // - The only window spans 200 bp.
    // - Average coverage is therefore:
    //     120 / 200 = 0.6
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: Some(200),
        by_bed: None,
        by_grouped_bed: None,
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t200\t0.6\t0",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_by_size_total_aligned_fast_path_matches_general_path() -> Result<()> {
    // Future-facing fast-path/general-path parity test.
    //
    // Manual expectations for 40 bp windows:
    // - [0, 40):   20 covered bases from the 40 bp fragment at value 1.5 -> 30
    // - [40, 80):  20 covered bases from the 40 bp fragment at value 1.5 -> 30
    // - [80, 120): 20 covered bases from the 80 bp fragment at value 0.75 -> 15
    // - [120,160): 40 covered bases from the 80 bp fragment at value 0.75 -> 30
    // - [160,200): 20 covered bases from the 80 bp fragment at value 0.75 -> 15
    let bam = mixed_length_fragment_bam()?;
    let mut outputs = Vec::new();

    for tile_size in [40_u32, 55_u32] {
        let out_dir = TempDir::new()?;
        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(2);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_tile_size(tile_size);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: Some(40),
            by_bed: None,
            by_grouped_bed: None,
        });
        set_restore_mean_length_normalization(&mut cfg);

        let result = run_inner(&cfg)?;
        outputs.push(read_zst_to_string(&result.final_out_path)?);
    }

    let expected = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t0\t40\t30\t0\n",
        "chr1\t40\t80\t30\t0\n",
        "chr1\t80\t120\t15\t0\n",
        "chr1\t120\t160\t30\t0\n",
        "chr1\t160\t200\t15\t0\n",
    );

    assert_eq!(outputs[0], expected);
    assert_eq!(outputs[1], expected);

    Ok(())
}

#[test]
fn restore_mean_by_size_summary_stats_derives_scaled_raw_and_derived_values() -> Result<()> {
    // Future-facing test for by-size summary stats.
    //
    // Manual expectations for the full 200 bp chromosome window:
    // - Coverage values:
    //     40 positions at 1.5
    //     80 positions at 0.75
    //     80 positions at 0
    // - coverage_sum = 40*1.5 + 80*0.75 = 120
    // - coverage_sum_of_squares = 40*(1.5^2) + 80*(0.75^2) = 90 + 45 = 135
    // - average_coverage = 120 / 200 = 0.6
    // - variance_coverage = (135 / 200) - 0.6^2 = 0.675 - 0.36 = 0.315
    // - sd_coverage = sqrt(0.315)
    // - coefficient_of_variation_coverage = sd / mean
    // - covered_fraction = 120 / 200 = 0.6
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(5);
    cfg.set_per_window(CoverageWindowAction::SummaryStats);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: Some(200),
        by_bed: None,
        by_grouped_bed: None,
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;
    let rows = parse_tsv(&text);

    assert_eq!(
        rows[0],
        vec![
            "chromosome",
            "start",
            "end",
            "span_positions",
            "blacklisted_positions",
            "eligible_positions",
            "nonzero_positions",
            "coverage_sum",
            "coverage_sum_of_squares",
            "average_coverage",
            "total_coverage",
            "variance_coverage",
            "sd_coverage",
            "coefficient_of_variation_coverage",
            "covered_fraction",
        ]
    );
    assert_eq!(
        rows[1][0..7],
        ["chr1", "0", "200", "200", "0", "200", "120"]
    );
    assert_close_to_written_precision(rows[1][7].parse::<f64>()?, 120.0, 5);
    assert_close_to_written_precision(rows[1][8].parse::<f64>()?, 135.0, 5);
    assert_close_to_written_precision(rows[1][9].parse::<f64>()?, 0.6, 5);
    assert_close_to_written_precision(rows[1][10].parse::<f64>()?, 120.0, 5);
    assert_close_to_written_precision(rows[1][11].parse::<f64>()?, 0.315, 5);
    assert_close_to_written_precision(rows[1][12].parse::<f64>()?, 0.315_f64.sqrt(), 5);
    assert_close(
        rows[1][13].parse::<f64>()?,
        0.315_f64.sqrt() / 0.6_f64,
        1e-5,
    );
    assert_close_to_written_precision(rows[1][14].parse::<f64>()?, 0.6, 5);

    Ok(())
}

#[test]
fn restore_mean_per_position_handles_three_chromosomes_in_global_mode() -> Result<()> {
    // Future-facing three-chromosome positional test.
    //
    // Manual expectations:
    // - chr1 fragment length 40 -> 60 / 40 = 1.5 on [20, 60)
    // - chr2 fragment length 60 -> 60 / 60 = 1.0 on [10, 70)
    // - chr3 fragment length 80 -> 60 / 80 = 0.75 on [40, 120)
    let bam = mixed_length_three_chromosome_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.chromosomes = base_chromosomes(&["chr1", "chr2", "chr3"]);
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chr1\t20\t60\t1.5",
            "chr2\t10\t70\t1",
            "chr3\t40\t120\t0.75",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_by_bed_average_handles_three_chromosomes_with_global_window_indices() -> Result<()>
{
    // Future-facing three-chromosome BED-average test.
    //
    // Manual expectations:
    // - Mean normalization length across fragment lengths 40, 60, 80 is 60.
    // - chr1 window [20, 60): full 40 bp fragment at 1.5 -> average 1.5
    // - chr2 window [0, 40): 30 covered bases at 1.0 -> average 30 / 40 = 0.75
    // - chr3 window [50, 100): 50 covered bases at 0.75 -> average 37.5 / 50 = 0.75
    let bam = mixed_length_three_chromosome_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("restore_mean_three_chr_windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 20, 60, "chr1_window"),
            ("chr2", 0, 40, "chr2_window"),
            ("chr3", 50, 100, "chr3_window"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.chromosomes = base_chromosomes(&["chr1", "chr2", "chr3"]);
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t20\t60\t1.5\t0",
            "chr2\t0\t40\t0.75\t0",
            "chr3\t50\t100\t0.75\t0",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_by_bed_average_skips_chromosomes_without_windows_and_keeps_later_chromosomes()
-> Result<()> {
    // Future-facing BED-average chromosome-skipping test.
    //
    // Manual expectations:
    // - We keep chromosomes chr1 and chr2.
    // - Only chr2 has a BED window: [0, 40).
    // - The chr2 fragment spans [10, 70) at restored value 1.0, so overlap inside [0,40) is 30.
    // - Average coverage is therefore 30 / 40 = 0.75.
    let bam = mixed_length_three_chromosome_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("restore_mean_chr2_only.bed");
    write_bed(&bed_path, &[("chr2", 0, 40, "chr2_window")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.chromosomes = base_chromosomes(&["chr1", "chr2"]);
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr2\t0\t40\t0.75\t0",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_by_size_total_handles_three_chromosomes() -> Result<()> {
    // Future-facing by-size multi-chromosome test.
    //
    // Manual expectations:
    // - Fragment normalization lengths are 40, 60, and 80, so mean_normalization_length = 60.
    // - Each fragment therefore contributes total mass 60 to its chromosome-wide 200 bp window.
    let bam = mixed_length_three_chromosome_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.chromosomes = base_chromosomes(&["chr1", "chr2", "chr3"]);
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: Some(200),
        by_bed: None,
        by_grouped_bed: None,
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t200\t60\t0",
            "chr2\t0\t200\t60\t0",
            "chr3\t0\t200\t60\t0",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_grouped_total_on_unique_bases_merges_same_group_overlaps_after_scaling()
-> Result<()> {
    // Future-facing grouped-total-on-unique-bases test.
    //
    // Manual expectations:
    // - Mean normalization length = 60.
    // - Group `beta` intervals:
    //     [20, 50) and [40, 60)   -> merge to [20, 60) with 40 bp at 1.5
    //     [100,150) and [140,180) -> merge to [100,180) with 80 bp at 0.75
    //   Unique union span = 120, total_coverage = 40*1.5 + 80*0.75 = 120
    // - Group `gamma` stays a zero row on [0, 10)
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("restore_mean_grouped_unique_total.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 40, 60, "beta"),
            ("chr1", 100, 150, "beta"),
            ("chr1", 140, 180, "beta"),
            ("chr1", 0, 10, "gamma"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(5);
    cfg.set_per_window(CoverageWindowAction::TotalOnUniqueBases);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;
    let rows = parse_tsv(&text);

    assert_eq!(
        rows,
        vec![
            vec![
                "group_idx",
                "span_positions",
                "blacklisted_positions",
                "eligible_positions",
                "total_coverage",
            ],
            vec!["0", "120", "0", "120", "120"],
            vec!["1", "10", "0", "10", "0"],
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_grouped_summary_stats_on_unique_bases_writes_scaled_rows() -> Result<()> {
    // Future-facing grouped summary-stat test.
    //
    // Manual expectations for `beta` after merging same-group overlaps:
    // - eligible = nonzero = span = 120
    // - coverage_sum = 120
    // - coverage_sum_of_squares = 40*(1.5^2) + 80*(0.75^2) = 135
    // - average_coverage = 120 / 120 = 1
    // - variance_coverage = (135 / 120) - 1^2 = 0.125
    // - sd_coverage = sqrt(0.125)
    // - coefficient_of_variation_coverage = sqrt(0.125)
    // - covered_fraction = 1
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir
        .path()
        .join("restore_mean_grouped_unique_summary.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 40, 60, "beta"),
            ("chr1", 100, 150, "beta"),
            ("chr1", 140, 180, "beta"),
            ("chr1", 0, 10, "gamma"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(5);
    cfg.set_per_window(CoverageWindowAction::SummaryStatsOnUniqueBases);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;
    let rows = parse_tsv(&text);

    assert_eq!(rows[1][0..5], ["0", "120", "0", "120", "120"]);
    assert_close_to_written_precision(rows[1][5].parse::<f64>()?, 120.0, 5);
    assert_close_to_written_precision(rows[1][6].parse::<f64>()?, 135.0, 5);
    assert_close_to_written_precision(rows[1][7].parse::<f64>()?, 1.0, 5);
    assert_close_to_written_precision(rows[1][8].parse::<f64>()?, 120.0, 5);
    assert_close_to_written_precision(rows[1][9].parse::<f64>()?, 0.125, 5);
    assert_close_to_written_precision(rows[1][10].parse::<f64>()?, 0.125_f64.sqrt(), 5);
    assert_close_to_written_precision(rows[1][11].parse::<f64>()?, 0.125_f64.sqrt(), 5);
    assert_close_to_written_precision(rows[1][12].parse::<f64>()?, 1.0, 5);

    assert_eq!(rows[2][0..5], ["1", "10", "0", "10", "0"]);
    assert_close_to_written_precision(rows[2][5].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows[2][6].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows[2][7].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows[2][8].parse::<f64>()?, 0.0, 5);
    assert!(rows[2][11].eq_ignore_ascii_case("NaN"));
    assert_close_to_written_precision(rows[2][12].parse::<f64>()?, 0.0, 5);

    Ok(())
}

#[test]
fn restore_mean_grouped_plain_summary_stats_writes_scaled_rows() -> Result<()> {
    // Future-facing plain grouped summary-stat test.
    //
    // Manual expectations for `beta`:
    // - Two separate same-group intervals are kept separate in plain grouped mode:
    //     [20, 50)  -> 30 bases at 1.5
    //     [100,160) -> 60 bases at 0.75
    // - span = eligible = nonzero = 90
    // - coverage_sum = 45 + 45 = 90
    // - coverage_sum_of_squares = 30*(1.5^2) + 60*(0.75^2) = 67.5 + 33.75 = 101.25
    // - average = 90 / 90 = 1
    // - variance = 101.25 / 90 - 1 = 0.125
    // - sd = sqrt(0.125)
    // - cv = sqrt(0.125)
    // - covered_fraction = 1
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir
        .path()
        .join("restore_mean_grouped_plain_summary.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 100, 160, "beta"),
            ("chr1", 0, 10, "gamma"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(5);
    cfg.set_per_window(CoverageWindowAction::SummaryStats);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let output = read_zst_to_string(&result.final_out_path)?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(
        rows_by_name["beta"][0..5]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "90", "0", "90", "90"]
    );
    assert_close_to_written_precision(rows_by_name["beta"][5].parse::<f64>()?, 90.0, 5);
    assert_close_to_written_precision(rows_by_name["beta"][6].parse::<f64>()?, 101.25, 5);
    assert_close_to_written_precision(rows_by_name["beta"][7].parse::<f64>()?, 1.0, 5);
    assert_close_to_written_precision(rows_by_name["beta"][8].parse::<f64>()?, 90.0, 5);
    assert_close_to_written_precision(rows_by_name["beta"][9].parse::<f64>()?, 0.125, 5);
    assert_close(
        rows_by_name["beta"][10].parse::<f64>()?,
        0.125_f64.sqrt(),
        1e-5,
    );
    assert_close(
        rows_by_name["beta"][11].parse::<f64>()?,
        0.125_f64.sqrt(),
        1e-5,
    );
    assert_close_to_written_precision(rows_by_name["beta"][12].parse::<f64>()?, 1.0, 5);

    assert_eq!(
        rows_by_name["gamma"][0..5]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["1", "10", "0", "10", "0"]
    );
    assert_close_to_written_precision(rows_by_name["gamma"][5].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows_by_name["gamma"][6].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows_by_name["gamma"][7].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows_by_name["gamma"][8].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows_by_name["gamma"][9].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows_by_name["gamma"][10].parse::<f64>()?, 0.0, 5);
    assert!(rows_by_name["gamma"][11].eq_ignore_ascii_case("NaN"));
    assert_close_to_written_precision(rows_by_name["gamma"][12].parse::<f64>()?, 0.0, 5);

    Ok(())
}

#[test]
fn restore_mean_grouped_plain_summary_stats_is_invariant_when_segments_cross_tiles() -> Result<()> {
    // Future-facing grouped summary tile-size invariance test.
    let bam = mixed_length_fragment_bam()?;
    let mut outputs = Vec::new();
    let mut sidecars = Vec::new();

    for tile_size in [33_u32, 1_000_u32] {
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir.path().join(format!(
            "restore_mean_grouped_plain_summary_{tile_size}.bed"
        ));
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 20, 50, "beta"),
                ("chr1", 100, 160, "beta"),
                ("chr1", 0, 10, "gamma"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_tile_size(tile_size);
        cfg.set_decimals(6);
        cfg.set_per_window(CoverageWindowAction::SummaryStats);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        set_restore_mean_length_normalization(&mut cfg);

        let result = run_inner(&cfg)?;
        outputs.push(read_zst_to_string(&result.final_out_path)?);
        sidecars.push(std::fs::read_to_string(
            out_dir.path().join("testcov.group_index.tsv"),
        )?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(sidecars[0], sidecars[1]);

    Ok(())
}

#[test]
fn restore_mean_segmented_fragment_still_multiplies_gc_and_scaling() -> Result<()> {
    // Future-facing restore-mean version of the segmented GC+scaling test.
    //
    // Manual expectations:
    // - The counted segments are [20, 40) and [50, 70), total normalization length 40.
    // - Because there is only one counted fragment, `mean_normalization_length = 40`.
    // - Unit-mass base weight is 1 / 40.
    // - GC weight = 2 and scaling factor = 5 everywhere, so unit-mass positional value is:
    //     (1 / 40) * 2 * 5 = 0.25
    // - `restore-mean` multiplies by 40, giving:
    //     0.25 * 40 = 10
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: 20,
            cigar: vec![('M', 20), ('D', 10), ('M', 20)],
            seq: vec![b'A'; 40],
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
        "fcoverage_restore_mean_gapped_with_gc_and_scaling",
    )?;
    let ref_twobit = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("restore_mean_gc_pkg.npz");
    let scaling_path = out_dir.path().join("restore_mean_scaling.tsv");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![50, 51],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
        correction_matrix: array![[2.0_f64]],
    };
    package.write_npz(&gc_path)?;
    write_scaling_factors(&scaling_path, &[("chr1", 0, 200, 5.0)])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.unpaired.reads_are_fragments = true;
    set_restore_mean_length_normalization(&mut cfg);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let mut scale_genome = ScaleGenomeArgs::default();
        scale_genome.scaling_factors = Some(scaling_path);
        cfg.set_scale_genome(scale_genome);
    }
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 50;
        frag.max_fragment_length = 50;
    }

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec!["chr1\t20\t40\t10", "chr1\t50\t70\t10"]
    );

    Ok(())
}

#[test]
fn restore_mean_keep_zero_runs_writes_zero_flanks_and_zero_gaps() -> Result<()> {
    // Future-facing test for `keep_zero_runs` under restore-mean.
    //
    // Manual expectations:
    // - Covered restore-mean runs are:
    //     [20, 60)   -> 1.5
    //     [100, 180) -> 0.75
    // - With `keep_zero_runs=true`, the uncovered flanks and middle gap remain as explicit zeros.
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(true);
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chr1\t0\t20\t0",
            "chr1\t20\t60\t1.5",
            "chr1\t60\t100\t0",
            "chr1\t100\t180\t0.75",
            "chr1\t180\t200\t0",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_blacklist_masks_positions_in_positional_output() -> Result<()> {
    // Future-facing blacklist test for positional restore-mean output.
    //
    // Manual expectations:
    // - Base restore-mean runs are [20, 60) -> 1.5 and [100, 180) -> 0.75.
    // - Blacklisting [30, 35) and [130, 150) removes those positions entirely.
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let blacklist_path = out_dir.path().join("restore_mean_blacklist.bed");
    write_bed(
        &blacklist_path,
        &[("chr1", 30, 35, "masked_a"), ("chr1", 130, 150, "masked_b")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.set_blacklist(Some(vec![blacklist_path]));
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chr1\t20\t30\t1.5",
            "chr1\t35\t60\t1.5",
            "chr1\t100\t130\t0.75",
            "chr1\t150\t180\t0.75",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_with_no_counted_fragments_keeps_zero_signal_and_returns_no_mean() -> Result<()> {
    // Future-facing zero-count bookkeeping test.
    //
    // Manual expectations:
    // - No fragments contribute, so all positional coverage is zero.
    // - `keep_zero_runs=true` therefore yields one full-chromosome zero run.
    // - `mean_normalization_length` is undefined, so the internal result should return `None`.
    let bam = empty_bam_fixture("fcoverage_restore_mean_empty")?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(true);
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(text, "chr1\t0\t200\t0\n");
    assert_eq!(result.mean_normalization_length, None);

    Ok(())
}

#[test]
fn restore_mean_by_bed_total_is_invariant_when_windows_cross_tile_boundaries() -> Result<()> {
    // Future-facing BED total tile-size invariance test.
    //
    // Manual expectations:
    // - Window [0, 40):   20 bases at 1.5 -> 30
    // - Window [20, 80):  40 bases at 1.5 -> 60
    // - Window [120,170): 50 bases at 0.75 -> 37.5
    let bam = mixed_length_fragment_bam()?;
    let mut outputs = Vec::new();

    for tile_size in [33_u32, 1_000_u32] {
        let out_dir = TempDir::new()?;
        let bed_path = out_dir
            .path()
            .join(format!("restore_mean_bed_total_{tile_size}.bed"));
        write_bed(
            &bed_path,
            &[
                ("chr1", 0, 40, "window_a"),
                ("chr1", 20, 80, "window_b"),
                ("chr1", 120, 170, "window_c"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(2);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_tile_size(tile_size);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        set_restore_mean_length_normalization(&mut cfg);

        let result = run_inner(&cfg)?;
        outputs.push(read_zst_to_string(&result.final_out_path)?);
    }

    let expected = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t0\t40\t30\t0\n",
        "chr1\t20\t80\t60\t0\n",
        "chr1\t120\t170\t37.5\t0\n",
    );
    assert_eq!(outputs[0], expected);
    assert_eq!(outputs[1], expected);

    Ok(())
}

#[test]
fn restore_mean_by_bed_summary_stats_is_invariant_when_windows_cross_tiles() -> Result<()> {
    // Future-facing BED summary-stats tile-size invariance test.
    //
    // Manual expectations:
    // - Window [0, 80):
    //     40 bases at 1.5 and 40 bases at 0
    //     coverage_sum = 60
    //     coverage_sum_of_squares = 90
    //     average = 0.75
    //     variance = 0.5625
    //     sd = 0.75
    //     cv = 1
    //     covered_fraction = 0.5
    // - Window [100, 180):
    //     80 bases at 0.75
    //     coverage_sum = 60
    //     coverage_sum_of_squares = 45
    //     average = 0.75
    //     variance = sd = cv = 0
    //     covered_fraction = 1
    let bam = mixed_length_fragment_bam()?;
    let mut outputs = Vec::new();

    for tile_size in [33_u32, 1_000_u32] {
        let out_dir = TempDir::new()?;
        let bed_path = out_dir
            .path()
            .join(format!("restore_mean_bed_summary_{tile_size}.bed"));
        write_bed(
            &bed_path,
            &[("chr1", 0, 80, "window_a"), ("chr1", 100, 180, "window_b")],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(5);
        cfg.set_per_window(CoverageWindowAction::SummaryStats);
        cfg.set_tile_size(tile_size);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        set_restore_mean_length_normalization(&mut cfg);

        let result = run_inner(&cfg)?;
        outputs.push(read_zst_to_string(&result.final_out_path)?);
    }

    assert_eq!(outputs[0], outputs[1]);
    let rows = parse_tsv(&outputs[0]);
    assert_eq!(rows[1][0..7], ["chr1", "0", "80", "80", "0", "80", "40"]);
    assert_close_to_written_precision(rows[1][7].parse::<f64>()?, 60.0, 5);
    assert_close_to_written_precision(rows[1][8].parse::<f64>()?, 90.0, 5);
    assert_close_to_written_precision(rows[1][9].parse::<f64>()?, 0.75, 5);
    assert_close_to_written_precision(rows[1][10].parse::<f64>()?, 60.0, 5);
    assert_close_to_written_precision(rows[1][11].parse::<f64>()?, 0.5625, 5);
    assert_close_to_written_precision(rows[1][12].parse::<f64>()?, 0.75, 5);
    assert_close_to_written_precision(rows[1][13].parse::<f64>()?, 1.0, 5);
    assert_close_to_written_precision(rows[1][14].parse::<f64>()?, 0.5, 5);

    assert_eq!(rows[2][0..7], ["chr1", "100", "180", "80", "0", "80", "80"]);
    assert_close_to_written_precision(rows[2][7].parse::<f64>()?, 60.0, 5);
    assert_close_to_written_precision(rows[2][8].parse::<f64>()?, 45.0, 5);
    assert_close_to_written_precision(rows[2][9].parse::<f64>()?, 0.75, 5);
    assert_close_to_written_precision(rows[2][10].parse::<f64>()?, 60.0, 5);
    assert_close_to_written_precision(rows[2][11].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows[2][12].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows[2][13].parse::<f64>()?, 0.0, 5);
    assert_close_to_written_precision(rows[2][14].parse::<f64>()?, 1.0, 5);

    Ok(())
}

#[test]
fn restore_mean_by_bed_total_keeps_coordinate_sorted_output_when_same_start_windows_cross_tiles()
-> Result<()> {
    // Future-facing coordinate-order regression test.
    //
    // Manual expectations:
    // - Window [0, 40) has restored total 30.
    // - Window [0, 100) has restored total 60.
    // - Final output should stay coordinate-sorted, so [0,40) precedes [0,100).
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("restore_mean_same_start_windows.bed");
    write_bed(
        &bed_path,
        &[("chr1", 0, 100, "wide"), ("chr1", 0, 40, "narrow")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_tile_size(33);
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t40\t30\t0",
            "chr1\t0\t100\t60\t0",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_by_bed_total_halo_only_window_is_not_double_counted_across_tiles() -> Result<()> {
    // Future-facing CoreOverlap regression test under restore-mean.
    //
    // Manual expectations:
    // - One 20 bp fragment spans [5, 25), so with restore-mean based on that one fragment the
    //   restored value is still 1 across the covered span.
    // - BED windows [5,10), [15,20), [20,25) each overlap 5 covered bases and must therefore
    //   each report total_coverage = 5 exactly once.
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 40)],
        vec![paired_fragment(5, 20, 10)],
        Vec::new(),
        "fcoverage_restore_mean_halo_only",
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("restore_mean_halo_only_windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 5, 10, "a"),
            ("chr1", 15, 20, "b"),
            ("chr1", 20, 25, "c"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_tile_size(10);
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let text = read_zst_to_string(&result.final_out_path)?;

    assert_eq!(
        text.lines().collect::<Vec<_>>(),
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t5\t10\t5\t0",
            "chr1\t15\t20\t5\t0",
            "chr1\t20\t25\t5\t0",
        ]
    );

    Ok(())
}

#[test]
fn restore_mean_grouped_bed_total_uses_site_weighted_group_semantics() -> Result<()> {
    // Future-facing grouped plain-total semantics test.
    //
    // Manual expectations:
    // - Group `beta` keeps its two loaded intervals separate:
    //     [20, 50)  -> 30 * 1.5 = 45
    //     [100,160) -> 60 * 0.75 = 45
    //   span_positions = eligible_positions = 90
    //   total_coverage = 90
    // - Group `alpha` = [20, 40) -> 20 * 1.5 = 30
    // - Group `gamma` = [0, 10) -> 0
    let bam = mixed_length_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("restore_mean_grouped_plain_total.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 20, 40, "alpha"),
            ("chr1", 100, 160, "beta"),
            ("chr1", 0, 10, "gamma"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(5);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let output = read_zst_to_string(&result.final_out_path)?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(
        rows_by_name["beta"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "90", "0", "90", "90"]
    );
    assert_eq!(
        rows_by_name["alpha"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["1", "20", "0", "20", "30"]
    );
    assert_eq!(
        rows_by_name["gamma"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["2", "10", "0", "10", "0"]
    );

    Ok(())
}

#[test]
fn restore_mean_grouped_bed_total_is_invariant_when_plain_group_segments_cross_tiles() -> Result<()>
{
    // Future-facing grouped plain-total tile-size invariance test.
    let bam = mixed_length_fragment_bam()?;
    let mut outputs = Vec::new();
    let mut sidecars = Vec::new();

    for tile_size in [33_u32, 1_000_u32] {
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir
            .path()
            .join(format!("restore_mean_grouped_plain_{tile_size}.bed"));
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 20, 50, "beta"),
                ("chr1", 20, 40, "alpha"),
                ("chr1", 100, 160, "beta"),
                ("chr1", 0, 10, "gamma"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_tile_size(tile_size);
        cfg.set_decimals(6);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });
        set_restore_mean_length_normalization(&mut cfg);

        let result = run_inner(&cfg)?;
        outputs.push(read_zst_to_string(&result.final_out_path)?);
        sidecars.push(std::fs::read_to_string(
            out_dir.path().join("testcov.group_index.tsv"),
        )?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(sidecars[0], sidecars[1]);

    Ok(())
}

#[test]
fn restore_mean_grouped_bed_ignores_groups_on_filtered_out_chromosomes() -> Result<()> {
    // Future-facing filtered-chromosome grouped test.
    //
    // Manual expectations:
    // - We count only `chr1`.
    // - Only the 40 bp chr1 fragment is counted, so the observed mean normalization length is 40.
    //   `alpha` on chr1 covers [0,100) and therefore gets total_coverage = 40.
    // - `beta` on chr2 must not appear because chr2 is filtered out.
    let bam = mixed_length_three_chromosome_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("restore_mean_grouped_filtered_chr.bed");
    write_bed(
        &grouped_bed,
        &[("chr1", 0, 100, "alpha"), ("chr2", 0, 120, "beta")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.chromosomes = base_chromosomes(&["chr1"]);
    cfg.set_decimals(5);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });
    set_restore_mean_length_normalization(&mut cfg);

    let result = run_inner(&cfg)?;
    let output = read_zst_to_string(&result.final_out_path)?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(rows_by_name.len(), 1);
    assert_eq!(
        rows_by_name["alpha"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "100", "0", "100", "40"]
    );

    Ok(())
}

#[test]
fn per_position_keep_zero_runs_toggles_zero_segments() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;

    let out_dir_without_zeros = TempDir::new()?;
    let mut cfg_without_zeros = base_config(&bam.bam, out_dir_without_zeros.path());
    cfg_without_zeros.set_decimals(0);
    cfg_without_zeros.set_per_window(CoverageWindowAction::Average);
    cfg_without_zeros.set_keep_zero_runs(false);

    // Manual expectations:
    // - The single fragment covers [20, 80) with coverage 1.
    // - The rest of chr1, with length 200, has coverage 0.
    // - When keep_zero_runs=false, only the covered run is written.
    run(&cfg_without_zeros)?;

    let without_zeros_path = out_dir_without_zeros
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let without_zeros_text = read_zst_to_string(&without_zeros_path)?;
    let without_zeros_lines: Vec<_> = without_zeros_text.lines().collect();
    assert_eq!(without_zeros_lines, vec!["chr1\t20\t80\t1"]);

    let out_dir_with_zeros = TempDir::new()?;
    let mut cfg_with_zeros = base_config(&bam.bam, out_dir_with_zeros.path());
    cfg_with_zeros.set_decimals(0);
    cfg_with_zeros.set_per_window(CoverageWindowAction::Average);
    cfg_with_zeros.set_keep_zero_runs(true);

    // With keep_zero_runs=true, the zero-coverage flanks are kept as separate runs.
    run(&cfg_with_zeros)?;

    let with_zeros_path = out_dir_with_zeros
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let with_zeros_text = read_zst_to_string(&with_zeros_path)?;
    let with_zeros_lines: Vec<_> = with_zeros_text.lines().collect();
    assert_eq!(
        with_zeros_lines,
        vec!["chr1\t0\t20\t0", "chr1\t20\t80\t1", "chr1\t80\t200\t0"]
    );

    Ok(())
}

#[test]
fn ignore_gap_removes_inter_mate_gap_from_positional_output() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.set_ignore_gap(true);

    // Manual expectations:
    // - The paired reads cover [20, 40) and [60, 80).
    // - Their fragment span is [20, 80), but the inter-mate gap is [40, 60).
    // - With ignore_gap=true, only the read-covered segments remain.
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t40\t1", "chr1\t60\t80\t1"]);

    Ok(())
}

#[test]
fn ignore_gap_keeps_scaling_and_blacklist_on_genomic_coordinates() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let blacklist_path = out_dir.path().join("ignore_gap_blacklist.bed");
    let scaling_path = out_dir.path().join("ignore_gap_scaling.tsv");
    write_bed(&blacklist_path, &[("chr1", 25, 35, "masked")])?;
    write_scaling_factors(
        &scaling_path,
        &[
            ("chr1", 0, 30, 2.0_f32),
            ("chr1", 30, 70, 3.0_f32),
            ("chr1", 70, 200, 5.0_f32),
        ],
    )?;

    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_keep_zero_runs(false);
    cfg.set_ignore_gap(true);
    cfg.set_scale_genome(scale_genome);
    cfg.blacklist = Some(vec![blacklist_path]);

    // Manual expectations:
    // - With ignore_gap=true, only the read-covered segments remain: [20, 40) and [60, 80).
    // - Scaling still applies on genomic coordinates:
    //   [20, 30): factor 2, [30, 70): factor 3, [70, 80): factor 5.
    // - Blacklisting [25, 35) still masks those reference positions after coverage/scaling.
    // - The final positional runs are therefore:
    //   [20, 25): 2
    //   [35, 40): 3
    //   [60, 70): 3
    //   [70, 80): 5
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chr1\t20\t25\t2",
            "chr1\t35\t40\t3",
            "chr1\t60\t70\t3",
            "chr1\t70\t80\t5",
        ]
    );

    Ok(())
}

#[test]
fn blacklist_inside_inter_mate_gap_only_matters_without_ignore_gap() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let blacklist_rows = [("chr1", 45, 55, "gap_mask")];

    let run_with_ignore_gap = |ignore_gap: bool| -> Result<Vec<String>> {
        let out_dir = TempDir::new()?;
        let blacklist_path = out_dir
            .path()
            .join(format!("gap_blacklist_{ignore_gap}.bed"));
        write_bed(&blacklist_path, &blacklist_rows)?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.set_per_window(CoverageWindowAction::Average);
        cfg.set_keep_zero_runs(false);
        cfg.set_ignore_gap(ignore_gap);
        cfg.blacklist = Some(vec![blacklist_path]);

        run(&cfg)?;

        let output_path = out_dir
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst");
        let text = read_zst_to_string(&output_path)?;
        Ok(text.lines().map(|line| line.to_string()).collect())
    };

    // Manual expectations:
    // - The paired reads cover [20, 40) and [60, 80), with an inter-mate gap [40, 60).
    // - The blacklist [45, 55) lies entirely inside that gap.
    // - Without ignore-gap, coverage is over the full fragment span [20, 80), so the blacklist
    //   removes part of the covered interval and splits the run.
    // - With ignore-gap, the gap is already absent from coverage, so masking inside the gap
    //   must not change the surviving read-covered segments.
    let without_ignore_gap = run_with_ignore_gap(false)?;
    let with_ignore_gap = run_with_ignore_gap(true)?;

    assert_eq!(
        without_ignore_gap,
        vec!["chr1\t20\t45\t1", "chr1\t55\t80\t1"]
    );
    assert_eq!(with_ignore_gap, vec!["chr1\t20\t40\t1", "chr1\t60\t80\t1"]);

    Ok(())
}

#[test]
fn unpaired_single_read_matches_fragment_span_output() -> Result<()> {
    // Human verification status: unverified
    let bam = single_read_fragment_bam("fcoverage_unpaired_single_read")?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.unpaired.reads_are_fragments = true;

    // Manual expectations:
    // - In unpaired mode, each read is treated as one fragment spanning its aligned
    //   reference interval.
    // - The single read spans [20, 80), so output should match a fragment over that span.
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t80\t1"]);

    Ok(())
}

#[test]
fn unpaired_single_read_matches_paired_fragment_output_for_same_span() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Compare two representations of the same physical fragment span [20, 80):
    // - paired-end fixture `simple_inward_bam()`
    // - one unpaired read with aligned span [20, 80)
    //
    // In paired mode the fragment span is [forward.pos, reverse.end).
    // In unpaired `reads_are_fragments` mode the fragment span is [read.pos, read.end).
    // So both inputs should yield the same positional coverage bedGraph:
    //   chr1  20  80  1
    let paired_bam = simple_inward_bam()?;
    let unpaired_bam = single_read_fragment_bam("fcoverage_unpaired_parity")?;
    let paired_out = TempDir::new()?;
    let unpaired_out = TempDir::new()?;

    let mut paired_cfg = base_config(&paired_bam.bam, paired_out.path());
    paired_cfg.set_decimals(0);
    paired_cfg.set_output_prefix("paired");

    let mut unpaired_cfg = base_config(&unpaired_bam.bam, unpaired_out.path());
    unpaired_cfg.set_decimals(0);
    unpaired_cfg.set_output_prefix("unpaired");
    unpaired_cfg.unpaired.reads_are_fragments = true;

    // Act
    run(&paired_cfg)?;
    run(&unpaired_cfg)?;

    // Assert
    let paired_text = read_zst_to_string(
        &paired_out
            .path()
            .join("paired.fcoverage.per_position.bedgraph.zst"),
    )?;
    let unpaired_text = read_zst_to_string(
        &unpaired_out
            .path()
            .join("unpaired.fcoverage.per_position.bedgraph.zst"),
    )?;
    let expected = "chr1\t20\t80\t1\n";
    assert_eq!(paired_text, expected);
    assert_eq!(unpaired_text, expected);

    Ok(())
}

#[test]
fn unpaired_mode_rejects_ignore_gap() -> Result<()> {
    // Human verification status: unverified
    let bam = single_read_fragment_bam("fcoverage_unpaired_ignore_gap")?;
    let baseline_out_dir = TempDir::new()?;
    let mut baseline_cfg = base_config(&bam.bam, baseline_out_dir.path());
    baseline_cfg.unpaired.reads_are_fragments = true;
    run(&baseline_cfg)?;
    let baseline_text = read_zst_to_string(
        &baseline_out_dir
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst"),
    )?;
    assert_eq!(baseline_text, "chr1\t20\t80\t1\n");

    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.unpaired.reads_are_fragments = true;
    cfg.set_ignore_gap(true);

    let err = run(&cfg).expect_err("unpaired mode should reject ignore_gap");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("--ignore-gap cannot be used with --reads-are-fragments"),
        "unexpected error: {msg}"
    );

    Ok(())
}

#[test]
fn unpaired_mode_rejects_require_proper_pair() -> Result<()> {
    // Human verification status: unverified
    let bam = single_read_fragment_bam("fcoverage_unpaired_require_pp")?;
    let baseline_out_dir = TempDir::new()?;
    let mut baseline_cfg = base_config(&bam.bam, baseline_out_dir.path());
    baseline_cfg.unpaired.reads_are_fragments = true;
    run(&baseline_cfg)?;
    let baseline_text = read_zst_to_string(
        &baseline_out_dir
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst"),
    )?;
    assert_eq!(baseline_text, "chr1\t20\t80\t1\n");

    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.unpaired.reads_are_fragments = true;
    cfg.set_require_proper_pair(true);

    let err = run(&cfg).expect_err("unpaired mode should reject require_proper_pair");
    let msg = format!("{err:#}");
    assert!(
        msg.contains("--require-proper-pair cannot be used with --reads-are-fragments"),
        "unexpected error: {msg}"
    );

    Ok(())
}

#[test]
fn blacklist_masks_positions_in_positional_output() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let blacklist_path = out_dir.path().join("blacklist.bed");
    write_bed(&blacklist_path, &[("chr1", 30, 35, "masked")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.blacklist = Some(vec![blacklist_path]);

    // Manual expectations:
    // - Base coverage without masking is 1 on [20, 80).
    // - The blacklist removes [30, 35) from the output entirely, so the run splits.
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t30\t1", "chr1\t35\t80\t1"]);

    Ok(())
}

#[test]
fn blacklist_masks_positions_in_positional_output_across_tile_boundary() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let blacklist_path = out_dir.path().join("blacklist_cross_tile_positional.bed");
    write_bed(&blacklist_path, &[("chr1", 30, 35, "masked")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_tile_size(33);
    cfg.blacklist = Some(vec![blacklist_path]);

    // Manual expectations:
    // - Base coverage without masking is 1 on [20, 80).
    // - The blacklist removes [30, 35) entirely.
    // - With tile_size=33, the blacklist itself crosses a tile boundary and the surviving
    //   covered span [35, 80) is also split at the later tile boundary 66.
    // - The positional output therefore keeps the same masked coordinates as the large-tile case,
    //   but the second covered run is split into [35, 66) and [66, 80).
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec!["chr1\t20\t30\t1", "chr1\t35\t66\t1", "chr1\t66\t80\t1"]
    );

    Ok(())
}

#[test]
fn blacklist_reduces_by_size_totals_and_reports_blacklisted_positions() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let blacklist_path = out_dir.path().join("blacklist.bed");
    write_bed(&blacklist_path, &[("chr1", 30, 35, "masked")])?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_size = Some(40);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(windows);
    cfg.blacklist = Some(vec![blacklist_path]);

    // Manual expectations:
    // - Window [0, 40): covered bases are [20, 30) and [35, 40), so total = 15.
    //   The masked interval contributes 5 blacklisted positions.
    // - Window [40, 80): fully covered, total = 40 and blacklisted_positions = 0.
    // - Remaining windows on chr1 have zero coverage and no blacklisted positions.
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t40\t15\t5",
            "chr1\t40\t80\t40\t0",
            "chr1\t80\t120\t0\t0",
            "chr1\t120\t160\t0\t0",
            "chr1\t160\t200\t0\t0",
        ]
    );

    Ok(())
}

#[test]
fn fcoverage_default_min_mapq_matches_explicit_thirty_and_differs_from_explicit_zero() -> Result<()>
{
    // Human verification status: unverified
    // Arrange:
    // Build three 60 bp fragments on a single 200 bp chromosome:
    // - fragment A: [20, 80), min MAPQ 60
    // - fragment B: [40, 100), min MAPQ 0
    // - fragment C: [120, 180), min MAPQ 30
    //
    // We reduce to one by-size window [0, 200) and request total coverage, so the
    // expected window totals are just the sum of per-fragment covered bases:
    // - default `min_mapq = 30`: A + C = 60 + 60 = 120
    // - explicit `min_mapq = 30`: same as default
    // - explicit `min_mapq = 0`: A + B + C = 180
    let fragment_with_mapq = |start: i64, mapq: u8| -> FragmentSpec {
        let mut fragment = fixtures::paired_fragment(start, 60, 20);
        fragment.forward.mapq = mapq;
        fragment.reverse.mapq = mapq;
        fragment
    };
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![
            fragment_with_mapq(20, 60),
            fragment_with_mapq(40, 0),
            fragment_with_mapq(120, 30),
        ],
        Vec::new(),
        "fcoverage_default_min_mapq",
    )?;
    let out_default = TempDir::new()?;
    let out_thirty = TempDir::new()?;
    let out_zero = TempDir::new()?;

    let make_cfg = |out_dir: &Path, prefix: &str| {
        let mut cfg = FCoverageConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
        );
        cfg.set_output_prefix(prefix);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: Some(200),
            by_bed: None,
            by_grouped_bed: None,
        });
        cfg.set_decimals(0);
        cfg.set_require_proper_pair(false);
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 10;
            frag.max_fragment_length = 200;
        }
        cfg
    };

    let default_cfg = make_cfg(out_default.path(), "default");
    let mut explicit_thirty_cfg = make_cfg(out_thirty.path(), "explicit_thirty");
    explicit_thirty_cfg.set_min_mapq(30);
    let mut explicit_zero_cfg = make_cfg(out_zero.path(), "explicit_zero");
    explicit_zero_cfg.set_min_mapq(0);

    // Act
    run(&default_cfg)?;
    run(&explicit_thirty_cfg)?;
    run(&explicit_zero_cfg)?;

    // Assert
    let default_text =
        read_zst_to_string(&out_default.path().join("default.fcoverage.total.tsv.zst"))?;
    let explicit_thirty_text = read_zst_to_string(
        &out_thirty
            .path()
            .join("explicit_thirty.fcoverage.total.tsv.zst"),
    )?;
    let explicit_zero_text = read_zst_to_string(
        &out_zero
            .path()
            .join("explicit_zero.fcoverage.total.tsv.zst"),
    )?;

    let expected_filtered = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t0\t200\t120\t0\n",
    );
    let expected_unfiltered = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t0\t200\t180\t0\n",
    );

    assert_eq!(default_text, expected_filtered);
    assert_eq!(explicit_thirty_text, expected_filtered);
    assert_eq!(explicit_zero_text, expected_unfiltered);

    Ok(())
}

#[test]
fn fcoverage_and_lengths_agree_on_the_single_fragment_that_survives_mapq_filtering() -> Result<()> {
    // Arrange:
    // Build two 60 bp fragments on one chromosome:
    // - fragment A: [20,80), MAPQ 60 -> kept by both commands
    // - fragment B: [100,160), MAPQ 0 -> removed by both commands at `min_mapq = 30`
    //
    // These commands report different quantities, so the parity check is about *which fragment
    // survives*, not about identical output formats:
    // - `lengths` should count exactly one fragment in the 60 bp bin
    // - `fcoverage` should write exactly one coverage run [20,80) with value 1
    let fragment_with_mapq = |start: i64, mapq: u8| -> FragmentSpec {
        let mut fragment = fixtures::paired_fragment(start, 60, 20);
        fragment.forward.mapq = mapq;
        fragment.reverse.mapq = mapq;
        fragment
    };
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![fragment_with_mapq(20, 60), fragment_with_mapq(100, 0)],
        Vec::new(),
        "fcoverage_lengths_mapq_parity",
    )?;

    let lengths_out = TempDir::new()?;
    let mut lengths_cfg = LengthsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: lengths_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    lengths_cfg.set_indel_mode(IndelMode::Ignore);
    lengths_cfg.set_windows(DistributionWindowsArgs::default());
    lengths_cfg.set_window_assignment(AssignToWindowArgs::default());
    lengths_cfg.set_min_mapq(30);
    lengths_cfg.set_require_proper_pair(false);
    lengths_cfg.set_per_bp_length_bins(10, 200);

    let fcoverage_out = TempDir::new()?;
    let mut fcoverage_cfg = base_config(&bam.bam, fcoverage_out.path());
    fcoverage_cfg.set_decimals(0);
    fcoverage_cfg.set_min_mapq(30);

    // Act
    run_lengths(&lengths_cfg)?;
    run(&fcoverage_cfg)?;

    // Assert
    let lengths_path = lengths_out.path().join(dot_join(&[
        lengths_cfg.output_prefix.trim(),
        "length_counts.npy",
    ]));
    let lengths_arr: Array2<f64> = read_npy(&lengths_path)?;
    assert_eq!(lengths_arr.shape(), &[1, 191]);
    let len60_idx = 60 - 10;
    assert!((lengths_arr[(0, len60_idx)] - 1.0).abs() < 1e-6);
    assert!((lengths_arr.sum() - 1.0).abs() < 1e-6);

    let fcoverage_text = read_zst_to_string(
        &fcoverage_out
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst"),
    )?;
    let fcoverage_lines: Vec<_> = fcoverage_text.lines().collect();
    assert_eq!(fcoverage_lines, vec!["chr1\t20\t80\t1"]);

    Ok(())
}

#[test]
fn blacklist_crossing_tile_boundary_keeps_same_by_size_output() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let tile_sizes = [33_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let blacklist_path = out_dir.path().join("blacklist_cross_tile.bed");
        write_bed(&blacklist_path, &[("chr1", 30, 35, "masked")])?;

        let mut windows = DistributionWindowsArgs::default();
        windows.by_size = Some(40);

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_windows(windows);
        cfg.set_tile_size(tile_size);
        cfg.blacklist = Some(vec![blacklist_path]);

        // Manual expectations:
        // - The single fragment covers [20, 80) with coverage 1.
        // - The blacklist masks [30, 35), and this interval crosses the 33 bp tile boundary.
        // - Window [0, 40): covered bases are [20, 30) and [35, 40), so total = 15.
        //   The masked interval contributes 5 blacklisted positions.
        // - Window [40, 80): fully covered, total = 40 and blacklisted_positions = 0.
        // - Remaining windows on chr1 have zero coverage and no blacklisted positions.
        // - Because this is a reduced by-size output, changing tile size must not change
        //   the final rows, even when the blacklist itself crosses a tile boundary.
        run(&cfg)?;

        let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
        outputs.push(read_zst_to_string(&output_path)?);
    }

    let expected = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t0\t40\t15\t5\n",
        "chr1\t40\t80\t40\t0\n",
        "chr1\t80\t120\t0\t0\n",
        "chr1\t120\t160\t0\t0\n",
        "chr1\t160\t200\t0\t0\n",
    );
    assert_eq!(outputs, vec![expected.to_string(), expected.to_string()]);

    Ok(())
}

#[test]
fn blacklist_average_uses_only_unmasked_positions_in_denominator() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let blacklist_path = out_dir.path().join("blacklist_average.bed");
    write_bed(&blacklist_path, &[("chr1", 30, 35, "masked")])?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_size = Some(40);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(3);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(windows);
    cfg.blacklist = Some(vec![blacklist_path]);

    // Manual expectations:
    // - The single fragment covers [20, 80) with coverage 1.
    // - The blacklist masks [30, 35), so window [0, 40) has:
    //   covered allowed bases = 15, blacklisted bases = 5, allowed span = 35.
    // - The average must therefore use 35 as denominator, giving 15 / 35 = 0.428571... -> 0.429.
    // - Window [40, 80): 40 / 40 = 1
    // - Remaining windows: 0
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.average.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t40\t0.429\t5",
            "chr1\t40\t80\t1\t0",
            "chr1\t80\t120\t0\t0",
            "chr1\t120\t160\t0\t0",
            "chr1\t160\t200\t0\t0",
        ]
    );

    Ok(())
}

#[test]
fn blacklist_average_writes_nan_when_window_has_no_eligible_positions() -> Result<()> {
    // Human verification status: verified
    // Manual expectations:
    // - Window [0,40) is fully blacklisted, so it has no eligible positions.
    // - No denominator exists for average coverage, so the scalar average must be undefined
    //   rather than true zero.
    // - Window [40,80) is unmasked and covered by the single fragment, so average = 1.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let blacklist_path = out_dir.path().join("blacklist_average_nan.bed");
    write_bed(&blacklist_path, &[("chr1", 0, 40, "fully_masked")])?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_size = Some(40);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(3);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(windows);
    cfg.blacklist = Some(vec![blacklist_path]);

    // Act
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.average.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();

    // Assert
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t40\tNaN\t40",
            "chr1\t40\t80\t1\t0",
            "chr1\t80\t120\t0\t0",
            "chr1\t120\t160\t0\t0",
            "chr1\t160\t200\t0\t0",
        ]
    );

    Ok(())
}

#[test]
fn per_position_and_by_size_totals_conserve_total_covered_bases() -> Result<()> {
    // Human verification status: unverified
    let bam = long_fragment_bam("fcoverage_conservation_fixture")?;
    let expected_total_coverage =
        (LONG_FRAGMENT_STARTS.len() as u64) * (LONG_FRAGMENT_LENGTH as u64);
    let tile_sizes = [1_000_u32, 1_100, 1_500, 1_700, 2_300];
    let mut positional_totals = Vec::new();
    let mut by_size_totals = Vec::new();

    // Manual expectations:
    // - The fixture contains 10 fragments.
    // - Each fragment spans 600 bp.
    // - Total covered bases are therefore 10 * 600 = 6000, regardless of overlaps.
    // - Changing tile size must not change that total for either positional output
    //   or by-size total windows.
    for tile_size in tile_sizes {
        let positional_out_dir = TempDir::new()?;
        let mut positional_cfg = base_config(&bam.bam, positional_out_dir.path());
        positional_cfg.set_decimals(0);
        positional_cfg.set_per_window(CoverageWindowAction::Average);
        positional_cfg.set_keep_zero_runs(false);
        positional_cfg.set_tile_size(tile_size);
        {
            let frag = positional_cfg.fragment_lengths_mut();
            frag.min_fragment_length = 100;
            frag.max_fragment_length = 700;
        }

        run(&positional_cfg)?;

        let positional_path = positional_out_dir
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst");
        let positional_text = read_zst_to_string(&positional_path)?;
        let positional_total: u64 = positional_text
            .lines()
            .map(|line| {
                let cols: Vec<_> = line.split('\t').collect();
                let start = cols[1].parse::<u64>().expect("bedgraph start should parse");
                let end = cols[2].parse::<u64>().expect("bedgraph end should parse");
                let value = cols[3]
                    .parse::<u64>()
                    .expect("bedgraph value should parse as integer");
                (end - start) * value
            })
            .sum();
        assert_eq!(
            positional_total, expected_total_coverage,
            "unexpected positional total for tile_size={tile_size}"
        );
        positional_totals.push(positional_total);

        let by_size_out_dir = TempDir::new()?;
        let mut by_size_windows = DistributionWindowsArgs::default();
        by_size_windows.by_size = Some(500);

        let mut by_size_cfg = base_config(&bam.bam, by_size_out_dir.path());
        by_size_cfg.set_decimals(0);
        by_size_cfg.set_per_window(CoverageWindowAction::Total);
        by_size_cfg.set_keep_zero_runs(true);
        by_size_cfg.set_tile_size(tile_size);
        by_size_cfg.set_windows(by_size_windows);
        {
            let frag = by_size_cfg.fragment_lengths_mut();
            frag.min_fragment_length = 100;
            frag.max_fragment_length = 700;
        }

        run(&by_size_cfg)?;

        let totals_path = by_size_out_dir
            .path()
            .join("testcov.fcoverage.total.tsv.zst");
        let totals_text = read_zst_to_string(&totals_path)?;
        let by_size_total: u64 = totals_text
            .lines()
            .skip(1)
            .map(|line| {
                let cols: Vec<_> = line.split('\t').collect();
                cols[3]
                    .parse::<u64>()
                    .expect("total_coverage should parse as integer")
            })
            .sum();
        assert_eq!(
            by_size_total, expected_total_coverage,
            "unexpected by-size total for tile_size={tile_size}"
        );
        by_size_totals.push(by_size_total);
    }

    assert_eq!(
        positional_totals,
        vec![expected_total_coverage; tile_sizes.len()]
    );
    assert_eq!(
        by_size_totals,
        vec![expected_total_coverage; tile_sizes.len()]
    );

    Ok(())
}

#[test]
fn global_positional_output_matches_single_full_chromosome_window_totals() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // `simple_inward_bam()` contains one fragment spanning [20, 80) on a 200 bp chromosome.
    //
    // Compare three logically equivalent ways of describing "the whole chromosome":
    // - global positional output
    // - one by-size window [0, 200)
    // - one BED window [0, 200)
    //
    // Hand-derived expectation:
    // - Global positional output is one bedGraph segment:
    //     chr1 20 80 1
    //   so the total covered mass is:
    //     (80 - 20) * 1 = 60
    // - The single full-chromosome by-size window must therefore report total_coverage = 60
    // - The single full-chromosome BED window must report the same total_coverage = 60
    let bam = simple_inward_bam()?;
    let global_out = TempDir::new()?;
    let by_size_out = TempDir::new()?;
    let by_bed_out = TempDir::new()?;
    let bed_path = by_bed_out.path().join("whole_chr_window.bed");
    write_bed(&bed_path, &[("chr1", 0, 200, "whole_chr")])?;

    let mut global_cfg = base_config(&bam.bam, global_out.path());
    global_cfg.set_decimals(0);

    let mut by_size_cfg = base_config(&bam.bam, by_size_out.path());
    by_size_cfg.set_decimals(0);
    by_size_cfg.set_per_window(CoverageWindowAction::Total);
    by_size_cfg.set_windows(DistributionWindowsArgs {
        by_size: Some(200),
        by_bed: None,
        by_grouped_bed: None,
    });

    let mut by_bed_cfg = base_config(&bam.bam, by_bed_out.path());
    by_bed_cfg.set_decimals(0);
    by_bed_cfg.set_per_window(CoverageWindowAction::Total);
    by_bed_cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });

    // Act
    run(&global_cfg)?;
    run(&by_size_cfg)?;
    run(&by_bed_cfg)?;

    // Assert
    let global_text = read_zst_to_string(
        &global_out
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst"),
    )?;
    let global_lines: Vec<_> = global_text.lines().collect();
    assert_eq!(global_lines, vec!["chr1\t20\t80\t1"]);

    let global_total: u64 = global_lines
        .iter()
        .map(|line| {
            let cols: Vec<_> = line.split('\t').collect();
            let start = cols[1].parse::<u64>().unwrap();
            let end = cols[2].parse::<u64>().unwrap();
            let value = cols[3].parse::<u64>().unwrap();
            (end - start) * value
        })
        .sum();
    assert_eq!(global_total, 60);

    let by_size_text =
        read_zst_to_string(&by_size_out.path().join("testcov.fcoverage.total.tsv.zst"))?;
    let by_bed_text =
        read_zst_to_string(&by_bed_out.path().join("testcov.fcoverage.total.tsv.zst"))?;
    let expected_window_total = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t0\t200\t60\t0\n",
    );

    assert_eq!(by_size_text, expected_window_total);
    assert_eq!(by_bed_text, expected_window_total);

    Ok(())
}

#[test]
fn by_bed_total_is_invariant_when_windows_cross_tile_boundaries() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let tile_sizes = [33_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let bed_path = out_dir
            .path()
            .join(format!("aggregate_windows_{tile_size}.bed"));
        write_bed(
            &bed_path,
            &[
                ("chr1", 0, 40, "window_a"),
                ("chr1", 20, 80, "window_b"),
                ("chr1", 70, 90, "window_c"),
            ],
        )?;

        let mut windows = DistributionWindowsArgs::default();
        windows.by_bed = Some(bed_path);

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_windows(windows);
        cfg.set_tile_size(tile_size);

        // Manual expectations:
        // - The single fragment covers [20, 80) with coverage 1.
        // - Window [0, 40): total covered bases = 20
        // - Window [20, 80): total covered bases = 60
        // - Window [70, 90): total covered bases = 10
        // - With tile_size=33, windows [0, 40) and [20, 80) cross one or more tile boundaries.
        //   The by-bed reducer must still merge the partials back to the same final rows.
        run(&cfg)?;

        let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
        outputs.push(read_zst_to_string(&output_path)?);
    }

    let expected = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t0\t40\t20\t0\n",
        "chr1\t20\t80\t60\t0\n",
        "chr1\t70\t90\t10\t0\n",
    );
    assert_eq!(outputs, vec![expected.to_string(), expected.to_string()]);

    Ok(())
}

#[test]
fn by_bed_total_mixed_core_and_downstream_windows_is_tile_size_invariant() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // - one unpaired fragment span [19,29), coverage 1 across those ten bases
    // - BED window [10,11) is upstream and should stay zero
    // - BED window [22,23) is downstream and should total 1
    // - with tile_size=10, the two windows fall into different tiles; with tile_size=1000 they do
    //   not. The final output must still be identical.
    let bam = single_read_fragment_bam_at("fcoverage_mixed_core_and_downstream_bed", 19, 10)?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let bed_path = out_dir
            .path()
            .join(format!("mixed_windows_{tile_size}.bed"));
        write_bed(
            &bed_path,
            &[
                ("chr1", 10, 11, "core_row"),
                ("chr1", 22, 23, "downstream_row"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.unpaired.reads_are_fragments = true;
        cfg.set_decimals(0);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        cfg.set_tile_size(tile_size);
        {
            let fragment_lengths = cfg.fragment_lengths_mut();
            fragment_lengths.min_fragment_length = 10;
            fragment_lengths.max_fragment_length = 10;
        }

        run(&cfg)?;

        let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
        outputs.push(read_zst_to_string(&output_path)?);
    }

    let expected = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t10\t11\t0\t0\n",
        "chr1\t22\t23\t1\t0\n",
    );
    assert_eq!(outputs, vec![expected.to_string(), expected.to_string()]);

    Ok(())
}

#[test]
fn by_bed_total_handles_window_spanning_more_than_two_tiles() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let tile_sizes = [30_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let bed_path = out_dir
            .path()
            .join(format!("aggregate_large_window_{tile_size}.bed"));
        write_bed(&bed_path, &[("chr1", 0, 100, "wide_window")])?;

        let mut windows = DistributionWindowsArgs::default();
        windows.by_bed = Some(bed_path);

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_windows(windows);
        cfg.set_tile_size(tile_size);

        // Manual expectations:
        // - The single fragment covers [20, 80), so window [0, 100) should sum to 60.
        // - With tile_size=30, the BED window spans four tiles:
        //   [0,30), [30,60), [60,90), and [90,120).
        // - The reducer must therefore merge 4 partial contributions back into one final row,
        //   exactly matching the large-tile case.
        run(&cfg)?;

        let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
        outputs.push(read_zst_to_string(&output_path)?);
    }

    let expected = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t0\t100\t60\t0\n",
    );
    assert_eq!(outputs, vec![expected.to_string(), expected.to_string()]);

    Ok(())
}

#[test]
fn by_size_total_counts_covered_bases_per_window() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_size = Some(40);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_keep_zero_runs(true);
    cfg.set_windows(windows);

    run(&cfg)?;

    let totals = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
    assert!(totals.exists(), "expected per-window totals output");
    let text = read_zst_to_string(&totals)?;
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap_or(""),
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions"
    );
    assert_eq!(lines.next().unwrap_or(""), "chr1\t0\t40\t20\t0");
    assert_eq!(lines.next().unwrap_or(""), "chr1\t40\t80\t40\t0");

    Ok(())
}

#[test]
fn by_size_average_reduces_across_non_aligned_tiles() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_size = Some(40);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(windows);
    cfg.set_tile_size(55);

    // Manual expectations:
    // - The single fragment covers [20, 80) with coverage 1.
    // - Window [0, 40): covered length 20 of 40 -> average 0.5
    // - Window [40, 80): covered length 40 of 40 -> average 1
    // - Remaining 40 bp windows have zero covered length -> average 0
    // - tile_size=55 forces the non-aligned reducer path, but should not change values.
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.average.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t40\t0.5\t0",
            "chr1\t40\t80\t1\t0",
            "chr1\t80\t120\t0\t0",
            "chr1\t120\t160\t0\t0",
            "chr1\t160\t200\t0\t0",
        ]
    );

    Ok(())
}

#[test]
fn by_size_total_keeps_full_bin_coordinates_when_bins_cross_tile_boundaries() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - One fragment covers [20, 80) with coverage 1.
    // - We request 40 bp bins on a 200 bp chromosome:
    //     [0, 40), [40, 80), [80, 120), [120, 160), [160, 200)
    // - With tile_size=33, both covered bins cross tile boundaries:
    //     [0, 40) overlaps tiles [0,33) and [33,66)
    //     [40, 80) overlaps tiles [33,66) and [66,99)
    // - The partial rows must still preserve the full bin coordinates, because the reducer merges
    //   by `bin_start`. If tile-local clipped coordinates leaked into the partial rows, the final
    //   output would split these into extra rows such as [0,33), [33,40), [40,66), [66,80).
    // - The correct final output therefore still has exactly one row per 40 bp bin:
    //     [0, 40) -> total 20
    //     [40, 80) -> total 40
    //     remaining bins -> total 0
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_size = Some(40);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(windows);
    cfg.set_tile_size(33);

    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t40\t20\t0",
            "chr1\t40\t80\t40\t0",
            "chr1\t80\t120\t0\t0",
            "chr1\t120\t160\t0\t0",
            "chr1\t160\t200\t0\t0",
        ]
    );

    Ok(())
}

#[test]
fn by_size_total_aligned_fast_path_matches_general_path_with_blacklist_scaling_and_gc() -> Result<()>
{
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let tile_sizes = [40_u32, 55_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let blacklist_path = out_dir
            .path()
            .join(format!("fast_path_blacklist_{tile_size}.bed"));
        let scaling_path = out_dir
            .path()
            .join(format!("fast_path_scaling_{tile_size}.tsv"));
        let gc_path = out_dir.path().join(format!("fast_path_gc_{tile_size}.npz"));
        write_bed(&blacklist_path, &[("chr1", 30, 35, "masked")])?;
        write_scaling_factors(
            &scaling_path,
            &[("chr1", 0, 50, 2.0_f32), ("chr1", 50, 200, 3.0_f32)],
        )?;
        build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;

        let mut scale_genome = ScaleGenomeArgs::default();
        scale_genome.scaling_factors = Some(scaling_path);

        let mut windows = DistributionWindowsArgs::default();
        windows.by_size = Some(40);

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_windows(windows);
        cfg.set_tile_size(tile_size);
        cfg.set_scale_genome(scale_genome);
        cfg.blacklist = Some(vec![blacklist_path]);
        cfg.set_gc(ApplyGCArgs {
            gc_file: Some(gc_path),
            gc_tag: None,
            neutralize_invalid_gc: false,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

        // Manual expectations:
        // - The fragment spans [20, 80), so the test GC package gives weight 10 across the fragment.
        // - Scaling then applies per genomic position:
        //   [20, 50): 10 * 2 = 20 coverage units per base
        //   [50, 80): 10 * 3 = 30 coverage units per base
        // - The blacklist removes [30, 35) entirely.
        // - By-size totals for 40 bp windows are therefore:
        //   [0, 40): [20, 30) contributes 10 * 20 = 200
        //            [35, 40) contributes  5 * 20 = 100
        //            total = 300, blacklisted_positions = 5
        //   [40, 80): [40, 50) contributes 10 * 20 = 200
        //            [50, 80) contributes 30 * 30 = 900
        //            total = 1100, blacklisted_positions = 0
        //   Remaining windows have zero coverage and zero blacklisted positions.
        // - tile_size=40 aligns exactly with by-size windows and should exercise the fast path.
        // - tile_size=55 forces the general reducer path.
        run(&cfg)?;

        let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
        outputs.push(read_zst_to_string(&output_path)?);
    }

    let expected = concat!(
        "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n",
        "chr1\t0\t40\t300\t5\n",
        "chr1\t40\t80\t1100\t0\n",
        "chr1\t80\t120\t0\t0\n",
        "chr1\t120\t160\t0\t0\n",
        "chr1\t160\t200\t0\t0\n",
    );
    assert_eq!(outputs, vec![expected.to_string(), expected.to_string()]);

    Ok(())
}

#[test]
fn by_bed_average_matches_manual_window_means() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("aggregate_windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 0, 40, "window_a"),
            ("chr1", 20, 80, "window_b"),
            ("chr1", 70, 90, "window_c"),
        ],
    )?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_bed = Some(bed_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(windows);

    // Manual expectations:
    // - The single fragment covers [20, 80) with coverage 1.
    // - Window [0, 40): covered length 20 of 40 -> average 0.5
    // - Window [20, 80): covered length 60 of 60 -> average 1
    // - Window [70, 90): covered length 10 of 20 -> average 0.5
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.average.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t40\t0.5\t0",
            "chr1\t20\t80\t1\t0",
            "chr1\t70\t90\t0.5\t0",
        ]
    );

    Ok(())
}

#[test]
fn by_bed_average_is_invariant_when_overlapping_windows_cross_tiles() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The single fragment covers [20, 80) with coverage 1.
    // - We use overlapping BED windows that exercise several average cases:
    //     [0, 40)   -> 20 covered bases out of 40 total, average = 20 / 40 = 1 / 2
    //     [10, 95)  -> 60 covered bases out of 85 total, average = 60 / 85 = 12 / 17
    //     [20, 80)  -> 60 covered bases out of 60 total, average = 1
    //     [70, 90)  -> 10 covered bases out of 20 total, average = 10 / 20 = 1 / 2
    // - `tile_size=33` forces the first three windows to cross one or more tile boundaries,
    //   while `tile_size=1000` keeps them in one tile. The final BED-average output must be
    //   identical because averages are properties of the full windows, not the temporary tile
    //   decomposition.
    let bam = simple_inward_bam()?;
    let tile_sizes = [33_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let bed_path = out_dir
            .path()
            .join(format!("average_windows_{tile_size}.bed"));
        write_bed(
            &bed_path,
            &[
                ("chr1", 0, 40, "window_a"),
                ("chr1", 10, 95, "window_b"),
                ("chr1", 20, 80, "window_c"),
                ("chr1", 70, 90, "window_d"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_tile_size(tile_size);
        cfg.set_decimals(6);
        cfg.set_per_window(CoverageWindowAction::Average);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });

        run(&cfg)?;

        outputs.push(read_zst_to_string(
            &out_dir.path().join("testcov.fcoverage.average.tsv.zst"),
        )?);
    }

    assert_eq!(
        outputs[0], outputs[1],
        "BED average output should not depend on tile size when windows cross tiles"
    );

    let rows = parse_tsv(&outputs[0]);
    assert_eq!(rows.len(), 5);
    assert_eq!(
        rows[0],
        vec![
            "chromosome",
            "start",
            "end",
            "average_coverage",
            "blacklisted_positions",
        ]
    );
    assert_eq!(rows[1][0..3], ["chr1", "0", "40"]);
    assert_close(rows[1][3].parse::<f64>()?, 1.0 / 2.0, 1e-9);
    assert_eq!(rows[1][4], "0");

    assert_eq!(rows[2][0..3], ["chr1", "10", "95"]);
    assert_close(rows[2][3].parse::<f64>()?, 12.0 / 17.0, 1e-6);
    assert_eq!(rows[2][4], "0");

    assert_eq!(rows[3][0..3], ["chr1", "20", "80"]);
    assert_close(rows[3][3].parse::<f64>()?, 1.0, 1e-9);
    assert_eq!(rows[3][4], "0");

    assert_eq!(rows[4][0..3], ["chr1", "70", "90"]);
    assert_close(rows[4][3].parse::<f64>()?, 1.0 / 2.0, 1e-9);
    assert_eq!(rows[4][4], "0");

    Ok(())
}

#[test]
fn by_bed_average_handles_three_chromosomes_with_global_window_indices() -> Result<()> {
    // Human verification status: unverified
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 200),
            ("chr2".to_string(), 200),
            ("chr3".to_string(), 200),
        ],
        vec![
            paired_fragment_on_tid(0, 20, 60, 20),
            paired_fragment_on_tid(1, 10, 40, 20),
            paired_fragment_on_tid(2, 40, 50, 20),
        ],
        Vec::new(),
        "fcoverage_three_chr_bed_avg",
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("aggregate_windows_three_chr.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 0, 40, "chr1_window"),
            ("chr2", 0, 40, "chr2_window"),
            ("chr3", 50, 100, "chr3_window"),
        ],
    )?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_bed = Some(bed_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.chromosomes = base_chromosomes(&["chr1", "chr2", "chr3"]);
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(windows);

    // Manual expectations:
    // - BED loading assigns one original window index across the whole file:
    //   chr1 -> 0, chr2 -> 1, chr3 -> 2.
    // - This test checks that aggregate BED reduction still works for later
    //   chromosomes instead of assuming chromosome-local dense indices.
    // - chr1 fragment spans [20, 80), so window [0, 40) has 20 covered bases -> 20 / 40 = 0.5.
    // - chr2 fragment spans [10, 50), so window [0, 40) has 30 covered bases -> 30 / 40 = 0.75.
    // - chr3 fragment spans [40, 90), so window [50, 100) has 40 covered bases -> 40 / 50 = 0.8.
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.average.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t0\t40\t0.5\t0",
            "chr2\t0\t40\t0.75\t0",
            "chr3\t50\t100\t0.8\t0",
        ]
    );

    Ok(())
}

#[test]
fn by_bed_average_skips_chromosomes_without_windows_and_keeps_later_chromosomes() -> Result<()> {
    // Human verification status: unverified
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![
            paired_fragment_on_tid(0, 20, 60, 20),
            paired_fragment_on_tid(1, 10, 40, 20),
        ],
        Vec::new(),
        "fcoverage_bed_average_skip_empty_chr",
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("aggregate_windows_chr2_only.bed");
    write_bed(&bed_path, &[("chr2", 0, 40, "chr2_window")])?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_bed = Some(bed_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.chromosomes = base_chromosomes(&["chr1", "chr2"]);
    cfg.set_decimals(2);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_windows(windows);

    // Manual expectations:
    // - chr1 has a fragment but no BED windows, so BED mode should skip that chromosome entirely.
    // - chr2 has one fragment on [10, 50) and one BED window [0, 40), so the covered overlap is
    //   30 bases and the mean coverage is 30 / 40 = 0.75.
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.average.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr2\t0\t40\t0.75\t0",
        ]
    );

    Ok(())
}

#[test]
fn per_position_handles_three_chromosomes_in_global_mode() -> Result<()> {
    // Human verification status: unverified
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 200),
            ("chr2".to_string(), 200),
            ("chr3".to_string(), 200),
        ],
        vec![
            paired_fragment_on_tid(0, 20, 60, 20),
            paired_fragment_on_tid(1, 10, 40, 20),
            paired_fragment_on_tid(2, 40, 50, 20),
        ],
        Vec::new(),
        "fcoverage_three_chr_global",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.chromosomes = base_chromosomes(&["chr1", "chr2", "chr3"]);
    cfg.set_decimals(0);
    cfg.set_keep_zero_runs(false);

    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec!["chr1\t20\t80\t1", "chr2\t10\t50\t1", "chr3\t40\t90\t1"]
    );

    Ok(())
}

#[test]
fn by_size_total_handles_three_chromosomes() -> Result<()> {
    // Human verification status: unverified
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 200),
            ("chr2".to_string(), 200),
            ("chr3".to_string(), 200),
        ],
        vec![
            paired_fragment_on_tid(0, 20, 60, 20),
            paired_fragment_on_tid(1, 10, 40, 20),
            paired_fragment_on_tid(2, 40, 50, 20),
        ],
        Vec::new(),
        "fcoverage_three_chr_size",
    )?;
    let out_dir = TempDir::new()?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_size = Some(200);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.chromosomes = base_chromosomes(&["chr1", "chr2", "chr3"]);
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(windows);

    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t200\t60\t0",
            "chr2\t0\t200\t40\t0",
            "chr3\t0\t200\t50\t0",
        ]
    );

    Ok(())
}

#[test]
fn by_bed_total_matches_manual_window_sums() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("aggregate_windows_total.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 0, 40, "window_a"),
            ("chr1", 20, 80, "window_b"),
            ("chr1", 70, 90, "window_c"),
        ],
    )?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_bed = Some(bed_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(windows);

    // Manual expectations:
    // - The single fragment covers [20, 80) with coverage 1.
    // - Window [0, 40): total covered bases = 20
    // - Window [20, 80): total covered bases = 60
    // - Window [70, 90): total covered bases = 10
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t40\t20\t0",
            "chr1\t20\t80\t60\t0",
            "chr1\t70\t90\t10\t0",
        ]
    );

    Ok(())
}

#[test]
fn by_bed_unique_positions_merge_overlapping_windows() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 15, 30, "window_a"),
            ("chr1", 25, 40, "window_b"),
            ("chr1", 70, 90, "window_c"),
        ],
    )?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_bed = Some(bed_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_per_window(CoverageWindowAction::OnlyIncludeThesePositionsUnique);
    cfg.set_keep_zero_runs(false);
    cfg.set_windows(windows);

    // Manual expectations:
    // - The single fragment covers [20, 80) with coverage 1 everywhere.
    // - unique-positions flattens overlapping BED windows before counting:
    //   [15, 30) and [25, 40) -> [15, 40), while [70, 90) stays separate.
    // - Intersecting those flattened windows with the covered span gives
    //   [20, 40) and [70, 80).
    // - Zero runs outside those intersections are omitted.
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t40\t1", "chr1\t70\t80\t1"]);

    Ok(())
}

#[test]
fn by_bed_indexed_positions_keep_window_indices_and_overlap_duplicates() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 15, 30, "window_a"),
            ("chr1", 25, 40, "window_b"),
            ("chr1", 70, 90, "window_c"),
        ],
    )?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_bed = Some(bed_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_per_window(CoverageWindowAction::OnlyIncludeThesePositionsIndexed);
    cfg.set_keep_zero_runs(false);
    cfg.set_windows(windows);

    // Manual expectations:
    // - The fragment still covers [20, 80) with value 1.
    // - indexed-positions keeps each original BED window separate and appends
    //   its original 0-based window index.
    // - Window 0 contributes [20, 30), window 1 contributes [25, 40), and
    //   window 2 contributes [70, 80).
    // - The overlap [25, 30) is intentionally duplicated because the windows
    //   are reported independently.
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position_per_window.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chr1\t20\t30\t1\t0",
            "chr1\t25\t40\t1\t1",
            "chr1\t70\t80\t1\t2",
        ]
    );

    Ok(())
}

#[test]
fn by_size_rejects_positional_per_window_modes() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    for action in [
        CoverageWindowAction::OnlyIncludeThesePositionsUnique,
        CoverageWindowAction::OnlyIncludeThesePositionsIndexed,
    ] {
        let mut windows = DistributionWindowsArgs::default();
        windows.by_size = Some(40);

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_per_window(action);
        cfg.set_windows(windows);

        let err = run(&cfg).expect_err("by-size positional window mode should fail");
        let err_text = err.to_string();
        assert!(
            err_text.contains(
                "in --by-size mode, --per-window can only be 'average', 'total', or 'summary-stats'"
            ),
            "unexpected error for {action:?}: {err_text}"
        );
    }

    Ok(())
}

#[test]
fn global_mode_rejects_total_per_window_choice() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - Without `--by-size`, `--by-bed`, or `--by-grouped-bed`, `fcoverage` writes positional
    //   coverage, not reduced per-window aggregates.
    // - `--per-window total` therefore reads like a meaningful aggregation choice but would be
    //   ignored in global mode.
    // - The command should fail fast instead of accepting a silent no-op.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_per_window(CoverageWindowAction::Total);

    let err = run(&cfg).expect_err("global mode should reject --per-window total");
    let err_text = err.to_string();
    assert!(
        err_text.contains(
            "without windowing, --per-window total is not supported because fcoverage writes positional coverage"
        ),
        "unexpected error: {err_text}"
    );

    Ok(())
}

#[test]
fn scaling_keeps_fractional_outputs_and_applies_rounding() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let scaling_path = out_dir.path().join("scaling.tsv");
    write_scaling_factors(
        &scaling_path,
        &[("chr1", 0, 50, 1.24), ("chr1", 50, 200, 1.26)],
    )?;

    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(1);
    cfg.set_keep_zero_runs(false);
    cfg.set_scale_genome(scale_genome);

    // Manual expectations:
    // - Unscaled positional coverage is 1 on [20, 80).
    // - Scaling bins split that covered region at 50bp:
    //   [20, 50) gets 1 * 1.24 = 1.24
    //   [50, 80) gets 1 * 1.26 = 1.26
    // - With decimals=1, those become 1.2 and 1.3 in the bedGraph output.
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t50\t1.2", "chr1\t50\t80\t1.3"]);

    Ok(())
}

#[test]
fn scaling_tsv_must_cover_requested_chromosome_end_in_fcoverage() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // `simple_inward_bam()` uses chr1 length 200.
    // A valid scaling TSV must cover the chromosome contiguously from 0 up to that exact length.
    //
    // Here we intentionally stop at 100:
    //   [0,100) factor 2.0
    // so the artifact is malformed for this chromosome even though the covered fragment itself
    // lies inside the provided span. The command should fail while loading the scaling artifact,
    // before any counting starts.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let scaling_path = out_dir.path().join("truncated_scaling.tsv");
    write_scaling_factors(&scaling_path, &[("chr1", 0, 100, 2.0)])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);
    cfg.set_scale_genome(scale_genome);

    // Act
    let err = run(&cfg).expect_err("truncated scaling TSV should fail");

    // Assert:
    // The shared scaling loader validates exact chromosome coverage, so the error should pin the
    // missing tail explicitly rather than allowing partial scaling.
    let err_text = format!("{err:#}");
    assert!(
        err_text.contains("scaling TSV: bins on 'chr1' must end at chrom_len=200 (got end=100)"),
        "unexpected error message: {err_text}"
    );

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn real_coverage_weights_tsv_changes_fcoverage_per_base_not_by_fragment_average() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let weights_out_dir = out_dir.path().join("weights_out");
    std::fs::create_dir_all(&weights_out_dir)?;
    let weights_cfg = make_simple_coverage_weights_config(&weights_out_dir, &bam.bam);

    // Manual expectations:
    // - `coverage-weights` on the simple fixture yields stride-bin scaling factors:
    //   [0,20): 37/20
    //   [20,40): 37/45
    //   [40,60): 37/60
    //   [60,80): 37/45
    //   [80,100): 37/15
    //   remaining bins: 0
    // - `fcoverage` applies scaling per covered base in place, not as one fragment-average
    //   multiplier.
    // - The only covered positions are [20,80), so the final positional output must be:
    //     [20,40): 37/45
    //     [40,60): 37/60
    //     [60,80): 37/45
    // - With `decimals = 6`, those become:
    //     [20,40): 0.822222
    //     [40,60): 0.616667
    //     [60,80): 0.822222
    // - If `fcoverage` incorrectly used full-fragment averaging like `midpoints` or the
    //   converters, this would instead collapse to one constant run with value 407/540.
    run_coverage_weights(&weights_cfg)?;
    let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");

    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_keep_zero_runs(false);
    cfg.set_scale_genome(scale_genome);

    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines.len(), 3, "expected three scaled bedGraph runs");

    let expected = [
        (20_u64, 40_u64, 37.0_f64 / 45.0_f64),
        (40_u64, 60_u64, 37.0_f64 / 60.0_f64),
        (60_u64, 80_u64, 37.0_f64 / 45.0_f64),
    ];

    for (line, (expected_start, expected_end, expected_value)) in lines.iter().zip(expected) {
        let parts: Vec<_> = line.split('\t').collect();
        assert_eq!(parts.len(), 4, "unexpected bedGraph row: {line}");
        assert_eq!(parts[0], "chr1");
        assert_eq!(parts[1].parse::<u64>()?, expected_start);
        assert_eq!(parts[2].parse::<u64>()?, expected_end);
        let value = parts[3].parse::<f64>()?;
        assert!(
            (value - expected_value).abs() <= 1e-6,
            "expected value {expected_value} for row {line}, got {value}"
        );
    }

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn real_ref_gc_bias_gc_bias_and_coverage_weights_chain_is_coherent_in_fcoverage() -> Result<()> {
    // Arrange:
    // Build the smallest real artifact chain that reaches a released downstream consumer:
    //   ref-gc-bias -> gc-bias -> coverage-weights -> fcoverage
    //
    // Use the standard repeated-ACGT reference and the standard single fragment [20,80):
    // - the real GC producer chain for this fixture is neutral, because every 60 bp fragment has
    //   GC%=50 in both the reference package and the sample BAM
    // - the real coverage-weights TSV is non-trivial and already hand-derived elsewhere:
    //     [20,40): 37/45
    //     [40,60): 37/60
    //     [60,80): 37/45
    //
    // Because the GC package is neutral, the final `fcoverage` output must match the
    // scaling-only expectation exactly. This makes the chain a strong release-spine check:
    // if any producer writes incompatible semantics, the downstream bedGraph changes.
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let weights_out_dir = out_dir.path().join("weights_out");
    std::fs::create_dir_all(&weights_out_dir)?;
    let mut weights_cfg = make_simple_coverage_weights_config(&weights_out_dir, &bam.bam);
    let weights_gc_path = build_real_neutral_gc_package_for_range(
        &bam.bam,
        &ref_twobit.path,
        out_dir.path(),
        10,
        200,
    )?;
    let gc_path = build_real_neutral_gc_package(&bam.bam, &ref_twobit.path, out_dir.path(), 60)?;

    weights_cfg.set_gc(ApplyGCArgs {
        gc_file: Some(weights_gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    weights_cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

    run_coverage_weights(&weights_cfg)?;
    let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");

    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_keep_zero_runs(false);
    cfg.set_scale_genome(scale_genome);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let fragment_lengths = cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 60;
        fragment_lengths.max_fragment_length = 60;
    }

    // Act
    run(&cfg)?;

    // Assert
    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    let expected = [
        (20_u64, 40_u64, 37.0_f64 / 45.0_f64),
        (40_u64, 60_u64, 37.0_f64 / 60.0_f64),
        (60_u64, 80_u64, 37.0_f64 / 45.0_f64),
    ];

    assert_eq!(lines.len(), expected.len());
    for (line, (expected_start, expected_end, expected_value)) in lines.iter().zip(expected) {
        let parts: Vec<_> = line.split('\t').collect();
        assert_eq!(parts.len(), 4, "unexpected bedGraph row: {line}");
        assert_eq!(parts[0], "chr1");
        assert_eq!(parts[1].parse::<u64>()?, expected_start);
        assert_eq!(parts[2].parse::<u64>()?, expected_end);
        let value = parts[3].parse::<f64>()?;
        assert!(
            (value - expected_value).abs() <= 1e-6,
            "expected value {expected_value} for row {line}, got {value}"
        );
    }

    Ok(())
}

#[test]
fn near_zero_scaling_tsv_stays_finite_and_correct_in_fcoverage() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // `normalize_average_overlap_keeps_sparse_non_zero_scaling_finite` already pins the helper-level
    // normalization result for the near-zero regime:
    // - sparse bin factor   = 5000.5
    // - ordinary covered bin factor = 0.50005
    //
    // This test carries that artifact contract through a real downstream consumer.
    //
    // Use the standard single fragment [20, 80), so only positions inside that span are covered.
    // Feed `fcoverage` a contiguous scaling TSV with:
    //   [0,20):   0
    //   [20,40):  5000.5
    //   [40,60):  0.50005
    //   [60,80):  0.50005
    //   [80,200): 0
    //
    // `fcoverage` applies scaling per covered base, so the exact bedGraph must be:
    //   chr1 20 40 5000.5
    //   chr1 40 80 0.50005
    //
    // The important downstream contract is that the huge factor remains finite and is emitted as a
    // normal numeric run rather than collapsing to zero, NaN, or inf during counting/output.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let scaling_path = out_dir.path().join("near_zero_scaling.tsv");
    write_scaling_factors(
        &scaling_path,
        &[
            ("chr1", 0, 20, 0.0_f32),
            ("chr1", 20, 40, 5000.5_f32),
            ("chr1", 40, 60, 0.50005_f32),
            ("chr1", 60, 80, 0.50005_f32),
            ("chr1", 80, 200, 0.0_f32),
        ],
    )?;

    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_keep_zero_runs(false);
    cfg.set_scale_genome(scale_genome);

    // Act
    run(&cfg)?;

    // Assert
    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines.len(), 2, "expected two finite scaled bedGraph runs");

    let expected = [(20_u64, 40_u64, 5000.5_f64), (40_u64, 80_u64, 0.50005_f64)];
    for (line, (expected_start, expected_end, expected_value)) in lines.iter().zip(expected) {
        let parts: Vec<_> = line.split('\t').collect();
        assert_eq!(parts.len(), 4, "unexpected bedGraph row: {line}");
        assert_eq!(parts[0], "chr1");
        assert_eq!(parts[1].parse::<u64>()?, expected_start);
        assert_eq!(parts[2].parse::<u64>()?, expected_end);
        let value = parts[3].parse::<f64>()?;
        assert!(
            value.is_finite(),
            "expected finite scaled coverage, got {value}"
        );
        assert!(
            (value - expected_value).abs() <= 1e-6,
            "expected value {expected_value} for row {line}, got {value}"
        );
    }

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn real_multi_chromosome_coverage_weights_tsv_is_applied_per_chromosome_in_fcoverage() -> Result<()>
{
    // Human verification status: unverified
    // Arrange:
    // Build a real multi-chromosome scaling TSV and consume it through per-position coverage.
    //
    // Producer BAM:
    // - chr1 has one 61 bp fragment [20, 81)
    // - chr2 has two identical 61 bp fragments [20, 81)
    //
    // With `bin_size = stride = 20`, each scaling bin is just the average positional coverage in
    // that 20 bp interval.
    //
    // Non-zero producer bins are:
    // - chr1:
    //     [20,40): 1
    //     [40,60): 1
    //     [60,80): 1
    //     [80,100): 1/20
    // - chr2:
    //     [20,40): 2
    //     [40,60): 2
    //     [60,80): 2
    //     [80,100): 1/10
    //
    // Shared global mean over the 8 non-zero bins:
    //   ((3 * 1) + 1/20 + (3 * 2) + 1/10) / 8
    // = (3 + 1/20 + 6 + 1/10) / 8
    // = 183/160.
    //
    // Written scaling factors = mean / average_pos_coverage:
    // - chr1 [20,80):  (183/160) / 1    = 183/160
    // - chr1 [80,100): (183/160) / 1/20 = 183/8
    // - chr2 [20,80):  (183/160) / 2    = 183/320
    // - chr2 [80,100): (183/160) / 1/10 = 183/16
    //
    // Consumer BAM:
    // - one 61 bp fragment [20, 81) on chr1
    // - one 61 bp fragment [20, 81) on chr2
    //
    // `fcoverage` applies scaling per covered base, so the expected bedGraph is:
    // - chr1  [20,80): 183/160
    // - chr1  [80,81): 183/8
    // - chr2  [20,80): 183/320
    // - chr2  [80,81): 183/16
    // chr2 deliberately stacks two identical fragments at one start. Use strict identity so the
    // producer really contains three molecules and the derived scaling TSV is correct.
    let producer_bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![
            paired_fragment_on_tid(0, 20, 61, 20),
            paired_fragment_on_tid(1, 20, 61, 20),
            paired_fragment_on_tid(1, 20, 61, 20),
        ],
        Vec::new(),
        "fcoverage_multichrom_scaling_producer",
    )?;
    let consumer_bam = bam_from_specs(
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![
            paired_fragment_on_tid(0, 20, 61, 20),
            paired_fragment_on_tid(1, 20, 61, 20),
        ],
        Vec::new(),
        "fcoverage_multichrom_scaling_consumer",
    )?;
    let temp = TempDir::new()?;
    let weights_out_dir = temp.path().join("coverage_weights");
    std::fs::create_dir_all(&weights_out_dir)?;

    let mut scaling_cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: producer_bam.bam.clone(),
            output_dir: weights_out_dir.clone(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1", "chr2"]),
    );
    scaling_cfg.set_bin_size(20);
    scaling_cfg.set_stride(20);
    scaling_cfg.set_min_mapq(0);
    scaling_cfg.set_require_proper_pair(false);
    scaling_cfg.set_output_prefix("coverage".to_string());
    {
        let frag = scaling_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    // Act
    run_coverage_weights(&scaling_cfg)?;
    let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");

    let mut cfg = FCoverageConfig::new(
        IOCArgs {
            bam: consumer_bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1", "chr2"]),
    );
    cfg.set_output_prefix("testcov");
    cfg.set_tile_size(1_000);
    cfg.set_ignore_gap(false);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_decimals(6);
    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);
    cfg.set_scale_genome(scale_genome);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }
    run(&cfg)?;

    // Assert
    let output_path = temp
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    let expected = [
        ("chr1", 20_u64, 80_u64, 183.0_f64 / 160.0_f64),
        ("chr1", 80_u64, 81_u64, 183.0_f64 / 8.0_f64),
        ("chr2", 20_u64, 80_u64, 183.0_f64 / 320.0_f64),
        ("chr2", 80_u64, 81_u64, 183.0_f64 / 16.0_f64),
    ];
    assert_eq!(lines.len(), expected.len());
    for (line, (expected_chr, expected_start, expected_end, expected_value)) in
        lines.iter().zip(expected.iter())
    {
        let parts: Vec<_> = line.split('\t').collect();
        assert_eq!(parts.len(), 4, "unexpected bedGraph row: {line}");
        assert_eq!(parts[0], *expected_chr);
        assert_eq!(parts[1].parse::<u64>()?, *expected_start);
        assert_eq!(parts[2].parse::<u64>()?, *expected_end);
        let actual = parts[3].parse::<f64>()?;
        assert!(
            (actual - expected_value).abs() <= 1e-6,
            "expected value {expected_value} for row {line}, got {actual}"
        );
    }

    Ok(())
}

#[test]
fn gc_tag_weights_unpaired_positional_output() -> Result<()> {
    // Human verification status: unverified
    let base_bam = single_read_fragment_bam("fcoverage_gc_tag_base")?;
    let tagged_bam = bam_with_gc_tags(&base_bam.bam, "fcoverage_gc_tag_valid", &[Some(2.5)])?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&tagged_bam.bam, out_dir.path());
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_decimals(1);
    cfg.unpaired.reads_are_fragments = true;
    cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    });

    // Manual expectations:
    // - The single read spans [20, 80).
    // - In unpaired mode, its GC tag weight is applied directly to the fragment.
    // - Coverage is therefore 2.5 across the whole span.
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t80\t2.5"]);

    Ok(())
}

#[test]
fn normalize_by_length_and_gc_tag_weights_multiply_per_position() -> Result<()> {
    // Human verification status: unverified
    let base_bam = single_read_fragment_bam("fcoverage_normalize_by_length_gc_tag_base")?;
    let tagged_bam = bam_with_gc_tags(
        &base_bam.bam,
        "fcoverage_normalize_by_length_gc_tag_valid",
        &[Some(2.5)],
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&tagged_bam.bam, out_dir.path());
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_decimals(4);
    cfg.set_keep_zero_runs(false);
    cfg.unpaired.reads_are_fragments = true;
    cfg.set_normalize_by_length_mode(LengthNormalizationMode::UnitMass);
    cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    });

    // Manual expectations:
    // - The single read spans [20, 80), so the counted length is 60.
    // - `--normalize-by-length` gives each counted base weight 1 / 60.
    // - In unpaired mode, the GC tag weight 2.5 is applied directly to the fragment.
    // - Final per-base coverage is therefore:
    //     2.5 / 60 = 0.041666... -> 0.0417
    run(&cfg)?;

    let output_path = out_dir.path().join(dot_join(&[
        "testcov",
        "length_normalized",
        "fcoverage.per_position.bedgraph.zst",
    ]));
    let text = read_zst_to_string(&output_path)?;
    assert_eq!(text, "chr1\t20\t80\t0.0417\n");

    Ok(())
}

#[test]
fn gc_tag_averages_valid_mate_weights_in_paired_mode() -> Result<()> {
    // Human verification status: unverified
    let base_bam = simple_inward_bam()?;
    let tagged_bam = bam_with_gc_tags(
        &base_bam.bam,
        "fcoverage_gc_tag_paired_avg",
        &[Some(2.0), Some(4.0)],
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&tagged_bam.bam, out_dir.path());
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_decimals(0);
    cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    });

    // Manual expectations:
    // - The paired fragment spans [20, 80).
    // - Mate GC tags are 2.0 and 4.0.
    // - Paired GC-tag mode combines them as their average, so the fragment weight is 3.0.
    // - Coverage is therefore 3 across the full fragment span.
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t80\t3"]);

    Ok(())
}

#[cfg(feature = "cmd_bam_to_bam")]
#[test]
fn bam_to_bam_gc_file_output_drives_fcoverage_gc_tag_same_as_original_gc_file() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // `simple_inward_bam()` contains one paired fragment spanning [20, 80), length 60.
    // We build the smallest GC package that assigns a constant weight 3.0 to every 60 bp
    // fragment:
    // - length_edges = [60, 61]
    // - gc_edges     = [0, 101]
    // - correction_matrix = [[3.0]]
    //
    // Then we compare two logically equivalent released workflows:
    // 1. original BAM -> `fcoverage --gc-file <pkg>`
    // 2. original BAM -> `bam-to-bam --gc-file <pkg>` -> `fcoverage --gc-tag GC`
    //
    // The only fragment covers [20, 80), so both workflows must yield the same positional output:
    //   chr1  20  80  3
    let bam = simple_inward_bam()?;
    let reference = simple_reference_twobit()?;
    let temp = TempDir::new()?;
    let tagged_bam_path = temp.path().join("tagged_gc.bam");
    let gc_path = temp.path().join("constant_gc_pkg.npz");

    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![60, 61],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&reference.path)?,
        correction_matrix: array![[3.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut bam_to_bam_cfg = BamToBamConfig::new(
        bam.bam.clone(),
        tagged_bam_path.clone(),
        base_chromosomes(&["chr1"]),
    );
    bam_to_bam_cfg.skip_chromosome_sort = true;
    bam_to_bam_cfg.min_mapq = 0;
    bam_to_bam_cfg.set_gc(cfdnalab::commands::cli_common::ApplyGCArgFileOnly {
        gc_file: Some(gc_path.clone()),
        neutralize_invalid_gc: false,
    });
    bam_to_bam_cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let frag = bam_to_bam_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    let original_out = TempDir::new()?;
    let tagged_out = TempDir::new()?;
    let mut original_cfg = base_config(&bam.bam, original_out.path());
    original_cfg.set_per_window(CoverageWindowAction::Average);
    original_cfg.set_decimals(0);
    original_cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path.clone()),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    original_cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let frag = original_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    // Act 1: write the tagged BAM from the real producer and index it for downstream fetching.
    run_bam_to_bam(&bam_to_bam_cfg)?;
    build_bai_for_test_bam(&tagged_bam_path)?;

    // Act 2: compare original `--gc-file` consumption with downstream `--gc-tag GC`.
    run(&original_cfg)?;

    let mut tagged_cfg = base_config(&tagged_bam_path, tagged_out.path());
    tagged_cfg.set_per_window(CoverageWindowAction::Average);
    tagged_cfg.set_decimals(0);
    tagged_cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    });
    {
        let frag = tagged_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }
    run(&tagged_cfg)?;

    // Assert
    let original_text = read_zst_to_string(
        &original_out
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst"),
    )?;
    let tagged_text = read_zst_to_string(
        &tagged_out
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst"),
    )?;
    assert_eq!(original_text, tagged_text);
    assert_eq!(original_text.trim(), "chr1\t20\t80\t3");

    Ok(())
}

#[test]
fn gc_tag_paired_edge_cases_follow_fragment_combination_rules() -> Result<()> {
    // Human verification status: unverified
    let scenarios = [
        (
            "invalid_mate_neutralized",
            &[Some(2.0), Some(2_000.0)][..],
            true,
            false,
            vec!["chr1\t20\t80\t1"],
            1_u64,
            1_u64,
            1_u64,
        ),
        (
            "zero_mate_counts_with_zero_weight",
            &[Some(0.0), Some(4.0)][..],
            false,
            true,
            vec!["chr1\t0\t200\t0"],
            0_u64,
            0_u64,
            1_u64,
        ),
    ];

    for (
        name,
        tags,
        neutralize_invalid_gc,
        keep_zero_runs,
        expected_lines,
        expected_gc_failed_fragments,
        expected_gc_out_of_range_tags,
        expected_counted_fragments,
    ) in scenarios
    {
        let base_bam = simple_inward_bam()?;
        let tagged_bam =
            bam_with_gc_tags(&base_bam.bam, &format!("fcoverage_gc_tag_{name}"), tags)?;
        let out_dir = TempDir::new()?;

        let mut cfg = base_config(&tagged_bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.set_keep_zero_runs(keep_zero_runs);
        cfg.set_gc(ApplyGCArgs {
            gc_file: None,
            gc_tag: Some("GC".to_string()),
            neutralize_invalid_gc,
        });

        // Manual expectations:
        // - Scenario invalid_mate_neutralized:
        //   mate tags 2.0 and 2000.0 average to 1001.0 at the fragment level, which is above the
        //   supported GC-tag range.
        //   With neutralize_invalid_gc=true, fcoverage keeps the fragment with neutral weight 1.0
        //   -> [20, 80) with value 1.
        //   Because the fragment-level GC tag is invalid and out of range:
        //     gc_failed_fragments = 1
        //     gc_out_of_range_tags = 1
        //     counted_fragments = 1 after neutralization
        // - Scenario zero_mate_counts_with_zero_weight:
        //   an explicit zero on either mate takes precedence in `combine_gc_tag_values`, so the
        //   fragment-level GC tag is `Usable(0.0)`, not missing or invalid.
        //   The fragment still overlaps the tile core and is therefore counted as a fragment, but
        //   it adds zero coverage everywhere:
        //     gc_failed_fragments = 0
        //     gc_out_of_range_tags = 0
        //     counted_fragments = 1
        let run_result = run_inner(&cfg)?;

        let output_path = out_dir
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst");
        let text = read_zst_to_string(&output_path)?;
        let lines: Vec<_> = text.lines().collect();
        assert_eq!(
            lines, expected_lines,
            "unexpected lines for scenario {name}"
        );
        assert_eq!(
            run_result.counters.gc_failed_fragments, expected_gc_failed_fragments,
            "unexpected gc_failed_fragments for scenario {name}"
        );
        assert_eq!(
            run_result.counters.gc_out_of_range_tags, expected_gc_out_of_range_tags,
            "unexpected gc_out_of_range_tags for scenario {name}"
        );
        assert_eq!(
            run_result.counters.base.counted_fragments, expected_counted_fragments,
            "unexpected counted_fragments for scenario {name}"
        );
    }

    Ok(())
}

#[test]
fn gc_tag_missing_or_invalid_values_skip_by_default_or_neutralize() -> Result<()> {
    // Human verification status: unverified
    let scenarios = [
        ("missing_skipped", None, false, Vec::<&str>::new()),
        ("missing_neutralized", None, true, vec!["chr1\t20\t80\t1"]),
        (
            "out_of_range_skipped",
            Some(2_000.0),
            false,
            Vec::<&str>::new(),
        ),
        (
            "out_of_range_neutralized",
            Some(2_000.0),
            true,
            vec!["chr1\t20\t80\t1"],
        ),
    ];

    for (name, tag_value, neutralize_invalid_gc, expected_lines) in scenarios {
        let base_bam = single_read_fragment_bam(&format!("fcoverage_gc_tag_{name}_base"))?;
        let tagged_bam = bam_with_gc_tags(
            &base_bam.bam,
            &format!("fcoverage_gc_tag_{name}"),
            &[tag_value],
        )?;
        let out_dir = TempDir::new()?;

        let mut cfg = base_config(&tagged_bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.unpaired.reads_are_fragments = true;
        cfg.set_gc(ApplyGCArgs {
            gc_file: None,
            gc_tag: Some("GC".to_string()),
            neutralize_invalid_gc,
        });

        // Manual expectations:
        // - Missing tags and out-of-range tags both produce no usable GC weight.
        // - With neutralize_invalid_gc=false, the fragment is skipped by default.
        // - With neutralize_invalid_gc=true, fcoverage keeps it with neutral weight 1.0.
        run(&cfg)?;

        let output_path = out_dir
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst");
        let text = read_zst_to_string(&output_path)?;
        let lines: Vec<_> = text.lines().collect();
        assert_eq!(
            lines, expected_lines,
            "unexpected lines for scenario {name}"
        );
    }

    Ok(())
}

#[test]
fn gc_file_requires_ref_2bit() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("gc_pkg.npz");
    let ref_twobit = simple_reference_twobit()?;
    build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });

    let err = run(&cfg).expect_err("GC correction should require --ref-2bit");
    let msg = format!("{err:#}");
    assert!(msg.contains("--ref-2bit"), "unexpected error: {msg}");

    Ok(())
}

#[test]
fn gc_file_weights_positional_output_from_reference_package() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("gc_pkg.npz");
    build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

    // Manual expectations:
    // - The fragment length is 60, which lands in the [60, 200) length bin.
    // - The simple reference has 50% GC over [20, 80), which lands in the [50, 101) GC bin.
    // - The test GC package assigns weight 10.0 to that bin, so coverage is 10 on [20, 80).
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t80\t10"]);

    Ok(())
}

#[test]
fn gc_file_rejects_package_when_fragment_length_range_is_outside_supported_range() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // `simple_inward_bam()` contains one fragment of length 60.
    // Give `fcoverage` a GC package that only covers fragment lengths 10..=59:
    //   length_edges = [10, 59]
    // With command fragment length bounds set to exactly [60, 60], the consumer must reject the
    // package before any coverage counting starts.
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("gc_pkg_short.npz");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![10, 59],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
        correction_matrix: array![[1.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    // Act
    let err = run(&cfg).expect_err("out-of-range GC package should fail");

    // Assert
    let msg = err.to_string();
    assert!(
        msg.contains("fragment length range [60-60] is outside the range covered by the correction package [10-59]"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn gc_file_rejects_package_with_schema_version_mismatch() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Build the smallest possible syntactically valid GC correction package, but make the schema
    // version intentionally incompatible. The command should fail while loading the package, before
    // any reference lookup or coverage counting begins.
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("gc_pkg_bad_version.npz");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION + 1,
        end_offset: 0,
        length_edges: vec![10, 200],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
        correction_matrix: array![[1.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

    // Act
    let err = run(&cfg).expect_err("schema version mismatch should fail");

    // Assert
    let msg = err.to_string();
    assert!(
        msg.contains("GC correction package schema version mismatch"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn real_ref_gc_bias_then_gc_bias_package_is_neutral_in_single_bin_case_for_fcoverage() -> Result<()>
{
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let gc_path = build_real_neutral_gc_package(&bam.bam, &ref_twobit.path, out_dir.path(), 60)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let fragment_lengths = cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 60;
        fragment_lengths.max_fragment_length = 60;
    }

    // Manual expectations:
    // - `ref-gc-bias` is run for exactly one fragment length: 60 bp.
    // - On the simple reference ("ACGT" repeated), every 60 bp fragment has exactly 30 G/C bases,
    //   so every reference fragment lands in the same GC%=50 cell.
    // - `gc-bias` is then run on `simple_inward_bam`, which contains exactly one 60 bp fragment
    //   over the same repeated reference, so the cfDNA counts also land in that same single cell.
    // - With one populated cfDNA cell and one populated reference cell:
    //   - mean normalization of each 1x1 matrix gives 1.0
    //   - their ratio is 1.0
    //   - inversion keeps the correction at 1.0
    // - Therefore the real produced GC package must be neutral for downstream coverage:
    //   positional coverage remains 1 on [20, 80).
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t80\t1"]);

    Ok(())
}

#[test]
fn real_ref_gc_bias_then_gc_bias_package_changes_fcoverage_in_expected_direction() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Use the same real two-bin producer logic as the non-neutral `gc-bias` command test:
    //
    // Reference genome:
    // - chr1[0,100)   = all A
    // - chr1[100,200) = all C
    //
    // `ref-gc-bias` is run for length 10 with:
    // - all 191 valid starts sampled
    // - BED windows [0,91) and [100,191)
    // - under the `ref-gc-bias` fit rule they count starts 0..=81 and 100..=181, yielding
    //   balanced reference support: 82 pure-A starts and 82 pure-C starts
    //
    // `gc-bias` is then run on a sample BAM with:
    // - one A-only fragment [10,20)   -> GC%=0
    // - nine C-only fragments in [110,200) -> GC%=100
    //
    // The resulting real correction package is hand-derived as:
    // - GC bin edges [0, 1, 100]
    // - correction weights [5.0, 5/9]
    //
    // Downstream `fcoverage` should therefore apply:
    // - weight 5.0 to the A-only fragment [10,20)
    // - weight 5/9 to each C-only fragment
    //
    // The nine C fragments are disjoint 10 bp runs, so the final bedGraph should be:
    // - chr1 10 20 5
    // - chr1 110 200 5/9
    let reference = fixtures::twobit_from_sequences(
        "fcoverage_real_non_neutral_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;
    let starts = [10_i64, 110, 120, 130, 140, 150, 160, 170, 180, 190];
    let fragments = starts
        .into_iter()
        .map(|start| paired_fragment(start, 10, 5))
        .collect();
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        fragments,
        Vec::new(),
        "fcoverage_real_non_neutral_bam",
    )?;

    let out_dir = TempDir::new()?;
    // Reference is 200 bp and fragment length is 10, so there are exactly:
    //   200 - 10 + 1 = 191 valid starts.
    // Sampling all valid starts plus BED windows `[0,91)` and `[100,191)` still makes the
    // reference-side masses exactly balanced between the pure-A and pure-C bins under the
    // `ref-gc-bias` fit rule.
    let gc_path = build_real_non_neutral_gc_package(
        &bam.bam,
        &reference.path,
        out_dir.path(),
        10,
        "chr1\t0\t91\nchr1\t100\t191\n",
        191,
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);

    // Act
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 10;
    }
    run(&cfg)?;

    // Assert
    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines.len(), 2);

    let expected = [
        ("chr1", 10_u64, 20_u64, 5.0_f64),
        ("chr1", 110_u64, 200_u64, 5.0_f64 / 9.0_f64),
    ];
    for (line, (expected_chr, expected_start, expected_end, expected_value)) in
        lines.iter().zip(expected.iter())
    {
        let parts: Vec<_> = line.split('\t').collect();
        assert_eq!(parts.len(), 4, "unexpected bedGraph row: {line}");
        assert_eq!(parts[0], *expected_chr);
        assert_eq!(parts[1].parse::<u64>()?, *expected_start);
        assert_eq!(parts[2].parse::<u64>()?, *expected_end);
        let actual = parts[3].parse::<f64>()?;
        assert!(
            (actual - expected_value).abs() <= 1e-6,
            "expected value {expected_value} for row {line}, got {actual}"
        );
    }

    Ok(())
}

#[test]
fn gc_file_invalid_weights_skip_by_default_or_neutralize() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let scenarios = [
        ("skipped_by_default", false, Vec::<&str>::new()),
        ("neutralized", true, vec!["chr1\t20\t80\t1"]),
    ];

    for (name, neutralize_invalid_gc, expected_lines) in scenarios {
        let out_dir = TempDir::new()?;
        let gc_path = out_dir.path().join(format!("gc_pkg_{name}.npz"));
        build_gc_package(&gc_path, 26, twobit_contig_footprint(&ref_twobit.path)?)?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.set_gc(ApplyGCArgs {
            gc_file: Some(gc_path),
            gc_tag: None,
            neutralize_invalid_gc,
        });
        cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
        {
            let frag = cfg.fragment_lengths_mut();
            frag.min_fragment_length = 53;
        }

        // Manual expectations:
        // - The GC package uses end_offset 26, so the command validator requires
        //   min_fragment_length > 52. We set it to 53 so the run reaches GC weighting.
        // - Length 60 with end_offset 26 leaves only 8 bp for GC counting.
        // - The corrector requires at least 10 A/C/G/T bases, so it returns no weight.
        // - With neutralize_invalid_gc=false, the fragment is skipped by default.
        // - With neutralize_invalid_gc=true, fcoverage keeps it with neutral weight 1.0.
        run(&cfg)?;

        let output_path = out_dir
            .path()
            .join("testcov.fcoverage.per_position.bedgraph.zst");
        let text = read_zst_to_string(&output_path)?;
        let lines: Vec<_> = text.lines().collect();
        assert_eq!(
            lines, expected_lines,
            "unexpected lines for scenario {name}"
        );
    }

    Ok(())
}

#[test]
fn unique_positions_split_one_merged_window_into_multiple_runs() -> Result<()> {
    // Human verification status: unverified
    let bam = overlapping_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("windows_multi_run.bed");
    write_bed(
        &bed_path,
        &[("chr1", 15, 45, "window_a"), ("chr1", 45, 85, "window_b")],
    )?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_bed = Some(bed_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_per_window(CoverageWindowAction::OnlyIncludeThesePositionsUnique);
    cfg.set_keep_zero_runs(false);
    cfg.set_windows(windows);

    // Manual expectations:
    // - Fragment 1 covers [20, 80)
    // - Fragment 2 covers [30, 70)
    // - Combined positional coverage is:
    //   [20, 30) -> 1
    //   [30, 70) -> 2
    //   [70, 80) -> 1
    // - unique-positions first merges the touching BED windows
    //   [15, 45) and [45, 85) -> [15, 85)
    // - Intersecting [15, 85) with the covered span keeps the same three runs
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec!["chr1\t20\t30\t1", "chr1\t30\t70\t2", "chr1\t70\t80\t1",]
    );

    Ok(())
}

#[test]
fn indexed_positions_repeat_window_index_for_each_run_and_duplicate_overlap() -> Result<()> {
    // Human verification status: unverified
    let bam = overlapping_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("windows_multi_run_indexed.bed");
    write_bed(
        &bed_path,
        &[("chr1", 15, 75, "window_a"), ("chr1", 25, 85, "window_b")],
    )?;

    let mut windows = DistributionWindowsArgs::default();
    windows.by_bed = Some(bed_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_per_window(CoverageWindowAction::OnlyIncludeThesePositionsIndexed);
    cfg.set_keep_zero_runs(false);
    cfg.set_windows(windows);

    // Manual expectations:
    // - The two fragments still create coverage runs
    //   [20, 30) -> 1, [30, 70) -> 2, [70, 80) -> 1
    // - indexed-positions keeps each BED window separate:
    //   Window 0 = [15, 75) intersects those runs as
    //     [20, 30) -> 1, [30, 70) -> 2, [70, 75) -> 1
    //   Window 1 = [25, 85) intersects those runs as
    //     [25, 30) -> 1, [30, 70) -> 2, [70, 80) -> 1
    // - Each emitted run carries the original 0-based window index
    // - The overlap between the two BED windows is intentionally duplicated
    run(&cfg)?;

    let output_path = out_dir
        .path()
        .join("testcov.fcoverage.per_position_per_window.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chr1\t20\t30\t1\t0",
            "chr1\t30\t70\t2\t0",
            "chr1\t70\t75\t1\t0",
            "chr1\t25\t30\t1\t1",
            "chr1\t30\t70\t2\t1",
            "chr1\t70\t80\t1\t1",
        ]
    );

    Ok(())
}

// ─── CoreOverlap double-counting guard tests ─────────────────────────────────
//
// fcoverage uses the CoreOverlap model: a tile only emits coverage for BED
// windows that overlap its core region. If someone accidentally switched to a
// fragment-reach model (like `lengths` uses), a fragment spanning a tile
// boundary could cause a halo-only window to be counted by the wrong tile,
// inflating coverage. These tests pin the correct behavior.

#[test]
fn by_bed_total_halo_only_window_is_not_double_counted_across_tiles() -> Result<()> {
    // Scenario:
    //
    //   Chromosome:  0         10        20        30
    //                |---------|---------|---------|
    //   Tile cores:  [  tile0  )[  tile1  )[  tile2  )   (tile_size = 10)
    //
    //   Fragment:         [=========)                     unpaired read [5, 25)
    //
    //   BED windows:      [-----)                        window A: [5, 10)
    //                              [-----)               window B: [15, 20)
    //                                    [-----)         window C: [20, 25)
    //
    // Fragment [5, 25) physically covers all three windows.
    //
    // Tile 0 core [0, 10): overlaps window A [5, 10). Fetches the fragment
    //   (fragment starts at 5 which is in the core). Counts coverage in A only.
    //   Window B [15, 20) is outside tile 0's core — must NOT be counted here.
    //
    // Tile 1 core [10, 20): overlaps window B [15, 20). The fragment is fetched
    //   via the halo (it starts at 5, before the core). Counts coverage in B only.
    //   Window C [20, 25) is outside tile 1's core — must NOT be counted here.
    //
    // Tile 2 core [20, 30): overlaps window C [20, 25). The fragment is fetched
    //   via the halo. Counts coverage in C only.
    //
    // Correct output: each window gets total = 5 (5 positions × depth 1,
    //   contributed by exactly one tile).
    // Regression (fragment-reach model): windows B and/or C could get total = 10
    //   because tile 0 or tile 1 would also count them.

    let bam = single_read_fragment_bam_at("fcoverage_core_overlap_guard", 5, 20)?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("core_overlap_guard.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 5, 10, "window_a"),
            ("chr1", 15, 20, "window_b"),
            ("chr1", 20, 25, "window_c"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.unpaired.reads_are_fragments = true;
    cfg.set_decimals(0);
    cfg.set_tile_size(10);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });
    {
        let fl = cfg.fragment_lengths_mut();
        fl.min_fragment_length = 10;
        fl.max_fragment_length = 30;
    }

    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.fcoverage.total.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t5\t10\t5\t0",  // 5 positions covered, each with coverage 1
            "chr1\t15\t20\t5\t0", // same: 5 positions, coverage 1 each
            "chr1\t20\t25\t5\t0", // same
        ],
        "Each window must be counted exactly once. Double-counting from a \
         fragment-reach model would produce totals > 5."
    );

    Ok(())
}

#[test]
fn by_bed_average_halo_only_window_is_counted_once_regardless_of_tile_size() -> Result<()> {
    // Same fragment/window setup as the total test above, but we verify
    // tile-size invariance of the average output. This catches the
    // double-counting regression from a different angle: if a halo-only
    // window is counted by two tiles, the reduced average would change
    // when tile_size changes (because the duplicate contribution only
    // appears in the small-tile run).

    let bam = single_read_fragment_bam_at("fcoverage_core_overlap_tsi", 5, 20)?;
    let tile_sizes = [10_u32, 100_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let bed_path = out_dir
            .path()
            .join(format!("core_overlap_tsi_{tile_size}.bed"));
        write_bed(
            &bed_path,
            &[
                ("chr1", 5, 10, "window_a"),
                ("chr1", 15, 20, "window_b"),
                ("chr1", 20, 25, "window_c"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.unpaired.reads_are_fragments = true;
        cfg.set_decimals(2);
        cfg.set_tile_size(tile_size);
        cfg.set_per_window(CoverageWindowAction::Average);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });
        {
            let fl = cfg.fragment_lengths_mut();
            fl.min_fragment_length = 10;
            fl.max_fragment_length = 30;
        }

        run(&cfg)?;

        let output_path = out_dir.path().join("testcov.fcoverage.average.tsv.zst");
        outputs.push(read_zst_to_string(&output_path)?);
    }

    assert_eq!(
        outputs[0], outputs[1],
        "BED average output must be identical for tile_size={} and tile_size={}. \
         A difference indicates that a halo-only window was counted by an \
         additional tile in the small-tile run (double-counting regression).",
        tile_sizes[0], tile_sizes[1]
    );

    // Also verify the values are correct: each window has 5 positions,
    // all with coverage 1, so average = 1.00
    let lines: Vec<_> = outputs[0].lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\taverage_coverage\tblacklisted_positions",
            "chr1\t5\t10\t1\t0",
            "chr1\t15\t20\t1\t0",
            "chr1\t20\t25\t1\t0",
        ]
    );

    Ok(())
}

#[test]
fn grouped_bed_total_uses_site_weighted_group_semantics() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The single fragment covers [20, 80) with per-base coverage 1.
    // - Group `beta` has two same-group intervals:
    //     [20, 50) -> 30 covered bases
    //     [40, 80) -> 40 covered bases
    //   Plain grouped totals keep loaded intervals separate, so:
    //     span_positions = 30 + 40 = 70
    //     total_coverage = 30 + 40 = 70
    // - Group `alpha` has [20, 40), so total_coverage = 20.
    // - Group `gamma` has [0, 10), outside the fragment, so total_coverage = 0.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("grouped_site_weighted.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 20, 40, "alpha"),
            ("chr1", 40, 80, "beta"),
            ("chr1", 0, 10, "gamma"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(&out_dir.path().join("testcov.fcoverage.total.tsv.zst"))?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let output_rows = parse_tsv(&output);
    let sidecar_rows = parse_tsv(&sidecar);

    assert_eq!(
        output_rows,
        vec![
            vec![
                "group_idx",
                "span_positions",
                "blacklisted_positions",
                "eligible_positions",
                "total_coverage",
            ],
            vec!["0", "70", "0", "70", "70"],
            vec!["1", "20", "0", "20", "20"],
            vec!["2", "10", "0", "10", "0"],
        ]
    );
    assert_eq!(
        sidecar_rows,
        vec![
            vec!["group_idx", "group_name"],
            vec!["0", "beta"],
            vec!["1", "alpha"],
            vec!["2", "gamma"],
        ]
    );

    Ok(())
}

#[test]
fn grouped_bed_total_is_invariant_when_plain_group_segments_cross_tiles() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The single fragment covers [20, 80) with per-base coverage 1.
    // - Plain grouped totals keep loaded intervals separate, even when multiple intervals share
    //   the same group.
    // - Group `beta` therefore sums two interval contributions:
    //     [20, 50) -> 30 covered bases
    //     [40, 80) -> 40 covered bases
    //   So:
    //     span_positions = eligible_positions = total_coverage = 70
    // - Group `alpha` is [20, 40), so:
    //     span_positions = eligible_positions = total_coverage = 20
    // - Group `delta` is [0, 60), so only [20, 60) is covered:
    //     span_positions = eligible_positions = 60
    //     total_coverage = 40
    // - Group `gamma` is [0, 10), outside the fragment:
    //     span_positions = eligible_positions = 10
    //     total_coverage = 0
    // - `tile_size=33` forces the segment rows for `beta`, `alpha`, and `delta` through the
    //   cross-tile BED-basic reducer, while `tile_size=1000` keeps them in one tile. The final
    //   grouped totals and stable group-index sidecar must still be identical.
    let bam = simple_inward_bam()?;
    let tile_sizes = [33_u32, 1_000_u32];
    let mut outputs = Vec::new();
    let mut sidecars = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir
            .path()
            .join(format!("grouped_plain_total_{tile_size}.bed"));
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 20, 50, "beta"),
                ("chr1", 20, 40, "alpha"),
                ("chr1", 40, 80, "beta"),
                ("chr1", 0, 60, "delta"),
                ("chr1", 0, 10, "gamma"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_tile_size(tile_size);
        cfg.set_decimals(0);
        cfg.set_per_window(CoverageWindowAction::Total);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });

        run(&cfg)?;

        outputs.push(read_zst_to_string(
            &out_dir.path().join("testcov.fcoverage.total.tsv.zst"),
        )?);
        sidecars.push(std::fs::read_to_string(
            out_dir.path().join("testcov.group_index.tsv"),
        )?);
    }

    assert_eq!(
        outputs[0], outputs[1],
        "plain grouped total output should not depend on tile size"
    );
    assert_eq!(
        sidecars[0], sidecars[1],
        "group sidecar should not depend on tile size"
    );

    let rows_by_name = grouped_rows_by_name(&outputs[0], &sidecars[0]);

    assert_eq!(
        rows_by_name["beta"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "70", "0", "70", "70"]
    );
    assert_eq!(
        rows_by_name["alpha"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["1", "20", "0", "20", "20"]
    );
    assert_eq!(
        rows_by_name["delta"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["2", "60", "0", "60", "40"]
    );
    assert_eq!(
        rows_by_name["gamma"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["3", "10", "0", "10", "0"]
    );

    Ok(())
}

#[test]
fn grouped_bed_errors_when_chromosome_filter_excludes_every_group() -> Result<()> {
    // Human verification status: verified
    // Manual expectations:
    // - `base_config` selects only `chr1`.
    // - The grouped BED contains one valid grouped row on `chr2`.
    // - Chromosome filtering removes that row before grouping, leaving no selected grouped
    //   windows. The command should fail directly instead of writing header-only outputs.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("grouped_all_rows_filtered.bed");
    write_bed(&grouped_bed, &[("chr2", 0, 100, "filtered")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    let error = match run_inner(&cfg) {
        Ok(_) => panic!("grouped BED with no selected rows should fail"),
        Err(error) => error,
    };
    let message = error.to_string();

    assert!(
        message.contains(
            "grouped BED file did not contain any valid windows on the selected chromosomes"
        ),
        "unexpected error: {message}"
    );

    Ok(())
}

#[test]
fn grouped_bed_ignores_groups_on_filtered_out_chromosomes() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The BAM contains one 60 bp fragment on each of `chr1` and `chr2`.
    // - The command explicitly selects only `chr1`.
    // - The grouped BED contains:
    //     `alpha` on `chr1` with span [0, 100)
    //     `beta` on `chr2` with span [0, 120)
    // - Only the `chr1` chromosome should participate in grouped loading and counting.
    //   Therefore:
    //     sidecar contains only `alpha`
    //     output contains only one grouped row
    //     `alpha` has span_positions = eligible_positions = 100 and total_coverage = 60
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![
            paired_fragment_on_tid(0, 20, 60, 20),
            paired_fragment_on_tid(1, 40, 60, 20),
        ],
        Vec::new(),
        "grouped_filtered_chromosomes",
    )?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("grouped_filtered_chromosomes.bed");
    write_bed(
        &grouped_bed,
        &[("chr1", 0, 100, "alpha"), ("chr2", 0, 120, "beta")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(&out_dir.path().join("testcov.fcoverage.total.tsv.zst"))?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;

    assert_eq!(
        parse_tsv(&output),
        vec![
            vec![
                "group_idx",
                "span_positions",
                "blacklisted_positions",
                "eligible_positions",
                "total_coverage",
            ],
            vec!["0", "100", "0", "100", "60"],
        ]
    );
    assert_eq!(
        parse_tsv(&sidecar),
        vec![vec!["group_idx", "group_name"], vec!["0", "alpha"]]
    );

    Ok(())
}

#[test]
fn grouped_bed_total_on_unique_bases_merges_same_group_overlaps() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The fragment again covers [20, 80).
    // - Group `beta` has [20, 50) and [40, 80). Under `total-on-unique-bases`,
    //   same-group overlaps collapse to the union [20, 80), so:
    //     span_positions = 60
    //     total_coverage = 60
    // - `gamma` stays a zero row.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("grouped_unique_bases.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 40, 80, "beta"),
            ("chr1", 0, 10, "gamma"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::TotalOnUniqueBases);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.total_on_unique_bases.tsv.zst"),
    )?;
    let rows = parse_tsv(&output);

    assert_eq!(
        rows,
        vec![
            vec![
                "group_idx",
                "span_positions",
                "blacklisted_positions",
                "eligible_positions",
                "total_coverage",
            ],
            vec!["0", "60", "0", "60", "60"],
            vec!["1", "10", "0", "10", "0"],
        ]
    );

    Ok(())
}

#[test]
fn grouped_bed_average_on_unique_bases_merges_same_group_overlaps() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The fragment covers [20, 80) with per-base coverage 1.
    // - Group `beta` has [20, 50) and [40, 80). Under `average-on-unique-bases`,
    //   same-group overlaps collapse to the union [20, 80), so:
    //     span_positions = eligible_positions = 60
    //     average_coverage = 60 / 60 = 1
    // - `gamma` covers [0, 20) with no fragment support, so:
    //     span_positions = eligible_positions = 20
    //     average_coverage = 0
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("grouped_unique_bases_avg.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 40, 80, "beta"),
            ("chr1", 0, 20, "gamma"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::AverageOnUniqueBases);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.average_on_unique_bases.tsv.zst"),
    )?;
    let rows = parse_tsv(&output);

    assert_eq!(
        rows,
        vec![
            vec![
                "group_idx",
                "span_positions",
                "blacklisted_positions",
                "eligible_positions",
                "average_coverage",
            ],
            vec!["0", "60", "0", "60", "1"],
            vec!["1", "20", "0", "20", "0"],
        ]
    );

    Ok(())
}

#[test]
fn grouped_bed_total_with_blacklist_uses_site_weighted_group_semantics() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The fragment covers [20, 80) with per-base coverage 1.
    // - The blacklist masks [45, 55), so masked positions must not contribute to coverage.
    // - Group `beta` keeps its two loaded intervals separate in plain grouped mode:
    //     [20, 50) contributes 25 eligible covered bases and 5 blacklisted bases
    //     [40, 80) contributes 30 eligible covered bases and 10 blacklisted bases
    //   Therefore:
    //     span_positions = 30 + 40 = 70
    //     blacklisted_positions = 5 + 10 = 15
    //     eligible_positions = 25 + 30 = 55
    //     total_coverage = 25 + 30 = 55
    // - Group `alpha` is [20, 40), fully outside the blacklist:
    //     span_positions = eligible_positions = total_coverage = 20
    // - Group `gamma` is [0, 10), so it stays a zero row with no blacklisted positions.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("grouped_site_weighted_blacklist.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 20, 40, "alpha"),
            ("chr1", 40, 80, "beta"),
            ("chr1", 0, 10, "gamma"),
        ],
    )?;
    let blacklist_bed = out_dir
        .path()
        .join("grouped_site_weighted_blacklist_mask.bed");
    write_bed(&blacklist_bed, &[("chr1", 45, 55, "masked")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_blacklist(Some(vec![blacklist_bed]));
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(&out_dir.path().join("testcov.fcoverage.total.tsv.zst"))?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(
        rows_by_name["beta"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "70", "15", "55", "55"]
    );
    assert_eq!(
        rows_by_name["alpha"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["1", "20", "0", "20", "20"]
    );
    assert_eq!(
        rows_by_name["gamma"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["2", "10", "0", "10", "0"]
    );

    Ok(())
}

#[test]
fn grouped_bed_total_on_unique_bases_with_blacklist_merges_same_group_overlap_once() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The fragment covers [20, 80).
    // - Group `beta` merges [20, 50) and [40, 80) to the union [20, 80).
    // - The blacklist masks [45, 55) once inside that union, so:
    //     span_positions = 60
    //     blacklisted_positions = 10
    //     eligible_positions = 50
    //     total_coverage = 50
    // - Group `gamma` is [0, 10), so it remains a zero row.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir
        .path()
        .join("grouped_unique_bases_blacklist_total.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 40, 80, "beta"),
            ("chr1", 0, 10, "gamma"),
        ],
    )?;
    let blacklist_bed = out_dir
        .path()
        .join("grouped_unique_bases_blacklist_total_mask.bed");
    write_bed(&blacklist_bed, &[("chr1", 45, 55, "masked")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::TotalOnUniqueBases);
    cfg.set_blacklist(Some(vec![blacklist_bed]));
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.total_on_unique_bases.tsv.zst"),
    )?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(
        rows_by_name["beta"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "60", "10", "50", "50"]
    );
    assert_eq!(
        rows_by_name["gamma"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["1", "10", "0", "10", "0"]
    );

    Ok(())
}

#[test]
fn grouped_bed_average_on_unique_bases_with_blacklist_uses_only_eligible_positions() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The fragment covers [20, 80) with per-base coverage 1.
    // - The blacklist masks [30, 35).
    // - Group `beta` merges to [20, 80):
    //     span_positions = 60
    //     blacklisted_positions = 5
    //     eligible_positions = 55
    //     total eligible covered bases = 55
    //     average_coverage = 55 / 55 = 1
    // - Group `delta` is [0, 40):
    //     span_positions = 40
    //     blacklisted_positions = 5
    //     eligible_positions = 35
    //     covered eligible bases are [20, 30) and [35, 40), so total_coverage = 15
    //     average_coverage = 15 / 35 = 3 / 7
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir
        .path()
        .join("grouped_unique_bases_blacklist_avg.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 40, 80, "beta"),
            ("chr1", 0, 40, "delta"),
        ],
    )?;
    let blacklist_bed = out_dir
        .path()
        .join("grouped_unique_bases_blacklist_average_mask.bed");
    write_bed(&blacklist_bed, &[("chr1", 30, 35, "masked")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::AverageOnUniqueBases);
    cfg.set_blacklist(Some(vec![blacklist_bed]));
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.average_on_unique_bases.tsv.zst"),
    )?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(
        rows_by_name["beta"][0..4]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "60", "5", "55"]
    );
    assert_close(rows_by_name["beta"][4].parse::<f64>()?, 1.0, 1e-9);

    assert_eq!(
        rows_by_name["delta"][0..4]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["1", "40", "5", "35"]
    );
    assert_close(rows_by_name["delta"][4].parse::<f64>()?, 3.0 / 7.0, 1e-6);

    Ok(())
}

#[test]
fn grouped_bed_average_writes_nan_when_group_has_no_eligible_positions() -> Result<()> {
    // Human verification status: verified
    // Manual expectations:
    // - The single fragment covers [20, 80) with per-base coverage 1.
    // - Group `masked` covers [20, 50) and is fully blacklisted:
    //     span_positions = blacklisted_positions = 30
    //     eligible_positions = 0
    //     average_coverage = NaN because no denominator exists
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("grouped_average_fully_masked.bed");
    write_bed(&grouped_bed, &[("chr1", 20, 50, "masked")])?;
    let blacklist_bed = out_dir
        .path()
        .join("grouped_average_fully_masked_blacklist.bed");
    write_bed(&blacklist_bed, &[("chr1", 20, 50, "fully_masked")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(3);
    cfg.set_per_window(CoverageWindowAction::Average);
    cfg.set_blacklist(Some(vec![blacklist_bed]));
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(&out_dir.path().join("testcov.fcoverage.average.tsv.zst"))?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(
        rows_by_name["masked"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "30", "30", "0", "NaN"]
    );

    Ok(())
}

#[test]
fn grouped_bed_average_on_unique_bases_writes_nan_when_group_has_no_eligible_positions()
-> Result<()> {
    // Human verification status: verified
    // Manual expectations:
    // - The single fragment covers [20, 80) with per-base coverage 1.
    // - Group `masked` has [20, 50) and [40, 80). Under `average-on-unique-bases`,
    //   same-group overlaps collapse to the union [20, 80).
    // - The blacklist fully masks [20, 80), so:
    //     span_positions = blacklisted_positions = 60
    //     eligible_positions = 0
    //     average_coverage = NaN because no denominator exists
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir
        .path()
        .join("grouped_unique_bases_average_fully_masked.bed");
    write_bed(
        &grouped_bed,
        &[("chr1", 20, 50, "masked"), ("chr1", 40, 80, "masked")],
    )?;
    let blacklist_bed = out_dir
        .path()
        .join("grouped_unique_bases_average_fully_masked_blacklist.bed");
    write_bed(&blacklist_bed, &[("chr1", 20, 80, "fully_masked")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(3);
    cfg.set_per_window(CoverageWindowAction::AverageOnUniqueBases);
    cfg.set_blacklist(Some(vec![blacklist_bed]));
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.average_on_unique_bases.tsv.zst"),
    )?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(
        rows_by_name["masked"]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "60", "60", "0", "NaN"]
    );

    Ok(())
}

#[test]
fn grouped_summary_stats_on_unique_bases_with_blacklist_excludes_masked_positions_from_raw_and_derived_stats()
-> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - The fragment covers [20, 80) with per-base coverage 1.
    // - The blacklist masks [30, 35).
    // - Group `beta` merges to [20, 80):
    //     span_positions = 60
    //     blacklisted_positions = 5
    //     eligible_positions = nonzero_positions = coverage_sum = coverage_sum_of_squares = 55
    //     mean = 55 / 55 = 1
    //     variance = sd = cv = 0
    //     covered_fraction = 1
    // - Group `delta` is [0, 40):
    //     span_positions = 40
    //     blacklisted_positions = 5
    //     eligible_positions = 35
    //     nonzero_positions = 15
    //     coverage_sum = coverage_sum_of_squares = 15
    //     mean = 15 / 35 = 3 / 7
    //     variance = 15 / 35 - (3 / 7)^2 = 12 / 49
    //     sd = sqrt(12 / 49) = sqrt(12) / 7
    //     cv = (sqrt(12) / 7) / (3 / 7) = sqrt(12) / 3
    //     covered_fraction = 15 / 35 = 3 / 7
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir
        .path()
        .join("grouped_unique_bases_blacklist_summary.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 40, 80, "beta"),
            ("chr1", 0, 40, "delta"),
        ],
    )?;
    let blacklist_bed = out_dir
        .path()
        .join("grouped_unique_bases_blacklist_summary_mask.bed");
    write_bed(&blacklist_bed, &[("chr1", 30, 35, "masked")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::SummaryStatsOnUniqueBases);
    cfg.set_blacklist(Some(vec![blacklist_bed]));
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.summary_stats_on_unique_bases.tsv.zst"),
    )?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(
        rows_by_name["beta"][0..5]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "60", "5", "55", "55"]
    );
    assert_close(rows_by_name["beta"][5].parse::<f64>()?, 55.0, 1e-9);
    assert_close(rows_by_name["beta"][6].parse::<f64>()?, 55.0, 1e-9);
    assert_close(rows_by_name["beta"][7].parse::<f64>()?, 1.0, 1e-9);
    assert_close(rows_by_name["beta"][8].parse::<f64>()?, 55.0, 1e-9);
    assert_close(rows_by_name["beta"][9].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["beta"][10].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["beta"][11].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["beta"][12].parse::<f64>()?, 1.0, 1e-9);

    assert_eq!(
        rows_by_name["delta"][0..5]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["1", "40", "5", "35", "15"]
    );
    assert_close(rows_by_name["delta"][5].parse::<f64>()?, 15.0, 1e-9);
    assert_close(rows_by_name["delta"][6].parse::<f64>()?, 15.0, 1e-9);
    assert_close(rows_by_name["delta"][7].parse::<f64>()?, 3.0 / 7.0, 1e-6);
    assert_close(rows_by_name["delta"][8].parse::<f64>()?, 15.0, 1e-9);
    assert_close(rows_by_name["delta"][9].parse::<f64>()?, 12.0 / 49.0, 1e-6);
    assert_close(
        rows_by_name["delta"][10].parse::<f64>()?,
        (12.0_f64).sqrt() / 7.0,
        1e-6,
    );
    assert_close(
        rows_by_name["delta"][11].parse::<f64>()?,
        (12.0_f64).sqrt() / 3.0,
        1e-6,
    );
    assert_close(rows_by_name["delta"][12].parse::<f64>()?, 3.0 / 7.0, 1e-6);

    Ok(())
}

#[test]
fn grouped_summary_stats_with_blacklist_is_invariant_when_plain_group_segments_cross_tiles()
-> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - One fragment covers [20, 80) with per-base coverage 1.
    // - The blacklist masks [45, 55).
    // - Group `beta` keeps its two loaded intervals separate in plain grouped mode:
    //     [20, 50) contributes 25 eligible covered bases and 5 blacklisted bases
    //     [40, 80) contributes 30 eligible covered bases and 10 blacklisted bases
    //   Therefore:
    //     span_positions = 30 + 40 = 70
    //     blacklisted_positions = 5 + 10 = 15
    //     eligible_positions = nonzero_positions = coverage_sum = coverage_sum_of_squares = 55
    //     mean = 55 / 55 = 1
    //     variance = sd = cv = 0
    //     covered_fraction = 1
    // - Group `delta` is a single interval [0, 60):
    //     span_positions = 60
    //     blacklisted_positions = 10
    //     eligible_positions = 50
    //     nonzero_positions = coverage_sum = coverage_sum_of_squares = 30
    //     mean = 30 / 50 = 3 / 5
    //     variance = 30 / 50 - (3 / 5)^2 = 6 / 25
    //     sd = sqrt(6 / 25) = sqrt(6) / 5
    //     cv = (sqrt(6) / 5) / (3 / 5) = sqrt(6) / 3
    //     covered_fraction = 30 / 50 = 3 / 5
    // - Group `gamma` is [0, 10), so:
    //     span_positions = eligible_positions = 10
    //     blacklisted_positions = nonzero_positions = coverage_sum = coverage_sum_of_squares = 0
    //     mean = variance = sd = covered_fraction = 0
    //     cv = NaN because the mean is zero
    // - `tile_size=33` forces the group segments above to cross tile boundaries, while
    //   `tile_size=1000` keeps them inside one tile. The final grouped summary-stats output
    //   must still match exactly
    let bam = simple_inward_bam()?;
    let tile_sizes = [33_u32, 1_000_u32];
    let mut outputs = Vec::new();
    let mut sidecars = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let grouped_bed = out_dir
            .path()
            .join(format!("grouped_plain_summary_blacklist_{tile_size}.bed"));
        write_bed(
            &grouped_bed,
            &[
                ("chr1", 20, 50, "beta"),
                ("chr1", 40, 80, "beta"),
                ("chr1", 0, 60, "delta"),
                ("chr1", 0, 10, "gamma"),
            ],
        )?;
        let blacklist_bed = out_dir.path().join(format!(
            "grouped_plain_summary_blacklist_mask_{tile_size}.bed"
        ));
        write_bed(&blacklist_bed, &[("chr1", 45, 55, "masked")])?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_tile_size(tile_size);
        cfg.set_decimals(6);
        cfg.set_per_window(CoverageWindowAction::SummaryStats);
        cfg.set_blacklist(Some(vec![blacklist_bed]));
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        });

        run(&cfg)?;

        outputs.push(read_zst_to_string(
            &out_dir
                .path()
                .join("testcov.fcoverage.summary_stats.tsv.zst"),
        )?);
        sidecars.push(std::fs::read_to_string(
            out_dir.path().join("testcov.group_index.tsv"),
        )?);
    }

    assert_eq!(
        outputs[0], outputs[1],
        "plain grouped summary-stats output should not depend on tile size"
    );
    assert_eq!(
        sidecars[0], sidecars[1],
        "group sidecar should not depend on tile size"
    );

    let rows_by_name = grouped_rows_by_name(&outputs[0], &sidecars[0]);

    assert_eq!(
        rows_by_name["beta"][0..5]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "70", "15", "55", "55"]
    );
    assert_close(rows_by_name["beta"][5].parse::<f64>()?, 55.0, 1e-9);
    assert_close(rows_by_name["beta"][6].parse::<f64>()?, 55.0, 1e-9);
    assert_close(rows_by_name["beta"][7].parse::<f64>()?, 1.0, 1e-9);
    assert_close(rows_by_name["beta"][8].parse::<f64>()?, 55.0, 1e-9);
    assert_close(rows_by_name["beta"][9].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["beta"][10].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["beta"][11].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["beta"][12].parse::<f64>()?, 1.0, 1e-9);

    assert_eq!(
        rows_by_name["delta"][0..5]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["1", "60", "10", "50", "30"]
    );
    assert_close(rows_by_name["delta"][5].parse::<f64>()?, 30.0, 1e-9);
    assert_close(rows_by_name["delta"][6].parse::<f64>()?, 30.0, 1e-9);
    assert_close(rows_by_name["delta"][7].parse::<f64>()?, 3.0 / 5.0, 1e-6);
    assert_close(rows_by_name["delta"][8].parse::<f64>()?, 30.0, 1e-9);
    assert_close(rows_by_name["delta"][9].parse::<f64>()?, 6.0 / 25.0, 1e-6);
    assert_close(
        rows_by_name["delta"][10].parse::<f64>()?,
        6.0_f64.sqrt() / 5.0,
        1e-6,
    );
    assert_close(
        rows_by_name["delta"][11].parse::<f64>()?,
        6.0_f64.sqrt() / 3.0,
        1e-6,
    );
    assert_close(rows_by_name["delta"][12].parse::<f64>()?, 3.0 / 5.0, 1e-6);

    assert_eq!(
        rows_by_name["gamma"][0..5]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["2", "10", "0", "10", "0"]
    );
    assert_close(rows_by_name["gamma"][5].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["gamma"][6].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["gamma"][7].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["gamma"][8].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["gamma"][9].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["gamma"][10].parse::<f64>()?, 0.0, 1e-9);
    assert!(rows_by_name["gamma"][11].parse::<f64>()?.is_nan());
    assert_close(rows_by_name["gamma"][12].parse::<f64>()?, 0.0, 1e-9);

    Ok(())
}

#[test]
fn by_size_summary_stats_writes_expected_raw_and_derived_values() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - One fragment covers [20, 80) inside the only size window [0, 200).
    // - Raw stats:
    //     span_positions = eligible_positions = 200
    //     nonzero_positions = 60
    //     coverage_sum = 60
    //     coverage_sum_of_squares = 60
    // - Derived stats:
    //     mean = 60 / 200 = 0.3
    //     variance = 60 / 200 - 0.3^2 = 0.21
    //     sd = sqrt(0.21)
    //     cv = sd / 0.3
    //     covered_fraction = 60 / 200 = 0.3
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::SummaryStats);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: Some(200),
        by_bed: None,
        by_grouped_bed: None,
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.summary_stats.tsv.zst"),
    )?;
    let rows = parse_tsv(&output);

    assert_eq!(
        rows[0],
        vec![
            "chromosome",
            "start",
            "end",
            "span_positions",
            "blacklisted_positions",
            "eligible_positions",
            "nonzero_positions",
            "coverage_sum",
            "coverage_sum_of_squares",
            "average_coverage",
            "total_coverage",
            "variance_coverage",
            "sd_coverage",
            "coefficient_of_variation_coverage",
            "covered_fraction",
        ]
    );
    assert_eq!(rows[1][0..7], ["chr1", "0", "200", "200", "0", "200", "60"]);
    assert_close(rows[1][7].parse::<f64>()?, 60.0, 1e-9);
    assert_close(rows[1][8].parse::<f64>()?, 60.0, 1e-9);
    assert_close(rows[1][9].parse::<f64>()?, 0.3, 1e-6);
    assert_close(rows[1][10].parse::<f64>()?, 60.0, 1e-9);
    assert_close(rows[1][11].parse::<f64>()?, 0.21, 1e-6);
    assert_close(rows[1][12].parse::<f64>()?, 0.21_f64.sqrt(), 1e-6);
    assert_close(rows[1][13].parse::<f64>()?, 0.21_f64.sqrt() / 0.3, 1e-6);
    assert_close(rows[1][14].parse::<f64>()?, 0.3, 1e-6);

    Ok(())
}

#[test]
fn by_size_summary_stats_is_invariant_when_windows_cross_tile_boundaries() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - One fragment covers [20, 80) with per-base coverage 1.
    // - We request 40 bp windows on a 200 bp chromosome, giving:
    //   [0, 40), [40, 80), [80, 120), [120, 160), [160, 200)
    // - Window [0, 40):
    //     20 zeros + 20 ones
    //     span_positions = eligible_positions = 40
    //     nonzero_positions = 20
    //     coverage_sum = coverage_sum_of_squares = 20
    //     mean = 20 / 40 = 0.5
    //     variance = 20 / 40 - 0.5^2 = 0.25
    //     sd = 0.5
    //     cv = 1
    //     covered_fraction = 20 / 40 = 0.5
    // - Window [40, 80):
    //     40 ones
    //     span_positions = eligible_positions = 40
    //     nonzero_positions = 40
    //     coverage_sum = coverage_sum_of_squares = 40
    //     mean = total = 1
    //     variance = sd = cv = 0
    //     covered_fraction = 1
    // - Remaining windows contain only zeros:
    //     nonzero_positions = 0
    //     coverage_sum = coverage_sum_of_squares = 0
    //     mean = total = variance = sd = covered_fraction = 0
    //     cv = NaN because the mean is exactly zero
    // - `tile_size=33` forces cross-tile reduction for several windows, while `tile_size=1000`
    //   keeps the whole chromosome in one tile. The final summary-stats output must still match
    let bam = simple_inward_bam()?;
    let tile_sizes = [33_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_tile_size(tile_size);
        cfg.set_decimals(6);
        cfg.set_per_window(CoverageWindowAction::SummaryStats);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: Some(40),
            by_bed: None,
            by_grouped_bed: None,
        });

        run(&cfg)?;

        outputs.push(read_zst_to_string(
            &out_dir
                .path()
                .join("testcov.fcoverage.summary_stats.tsv.zst"),
        )?);
    }

    assert_eq!(
        outputs[0], outputs[1],
        "summary-stats output should not depend on tile size"
    );

    let rows = parse_tsv(&outputs[0]);
    assert_eq!(
        rows[0],
        vec![
            "chromosome",
            "start",
            "end",
            "span_positions",
            "blacklisted_positions",
            "eligible_positions",
            "nonzero_positions",
            "coverage_sum",
            "coverage_sum_of_squares",
            "average_coverage",
            "total_coverage",
            "variance_coverage",
            "sd_coverage",
            "coefficient_of_variation_coverage",
            "covered_fraction",
        ]
    );

    assert_eq!(rows[1][0..7], ["chr1", "0", "40", "40", "0", "40", "20"]);
    assert_close(rows[1][7].parse::<f64>()?, 20.0, 1e-9);
    assert_close(rows[1][8].parse::<f64>()?, 20.0, 1e-9);
    assert_close(rows[1][9].parse::<f64>()?, 0.5, 1e-9);
    assert_close(rows[1][10].parse::<f64>()?, 20.0, 1e-9);
    assert_close(rows[1][11].parse::<f64>()?, 0.25, 1e-9);
    assert_close(rows[1][12].parse::<f64>()?, 0.5, 1e-9);
    assert_close(rows[1][13].parse::<f64>()?, 1.0, 1e-9);
    assert_close(rows[1][14].parse::<f64>()?, 0.5, 1e-9);

    assert_eq!(rows[2][0..7], ["chr1", "40", "80", "40", "0", "40", "40"]);
    assert_close(rows[2][7].parse::<f64>()?, 40.0, 1e-9);
    assert_close(rows[2][8].parse::<f64>()?, 40.0, 1e-9);
    assert_close(rows[2][9].parse::<f64>()?, 1.0, 1e-9);
    assert_close(rows[2][10].parse::<f64>()?, 40.0, 1e-9);
    assert_close(rows[2][11].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[2][12].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[2][13].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[2][14].parse::<f64>()?, 1.0, 1e-9);

    let expected_zero_windows = [("80", "120"), ("120", "160"), ("160", "200")];
    for (row, (start, end)) in rows[3..].iter().zip(expected_zero_windows) {
        assert_eq!(row[0], "chr1");
        assert_eq!(row[1], start);
        assert_eq!(row[2], end);
        assert_eq!(row[3], "40");
        assert_eq!(row[4], "0");
        assert_eq!(row[5], "40");
        assert_eq!(row[6], "0");
        assert_close(row[7].parse::<f64>()?, 0.0, 1e-9);
        assert_close(row[8].parse::<f64>()?, 0.0, 1e-9);
        assert_close(row[9].parse::<f64>()?, 0.0, 1e-9);
        assert_close(row[10].parse::<f64>()?, 0.0, 1e-9);
        assert_close(row[11].parse::<f64>()?, 0.0, 1e-9);
        assert_close(row[12].parse::<f64>()?, 0.0, 1e-9);
        assert!(row[13].parse::<f64>()?.is_nan());
        assert_close(row[14].parse::<f64>()?, 0.0, 1e-9);
    }

    Ok(())
}

#[test]
fn by_bed_summary_stats_derives_variance_from_overlapping_fragments() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - Fragment 1 covers [20, 80) with coverage 1.
    // - Fragment 2 covers [30, 70) with coverage 1 on top of fragment 1.
    // - Window [0, 100) therefore sees:
    //     [0, 20)   -> 0 for 20 bp
    //     [20, 30)  -> 1 for 10 bp
    //     [30, 70)  -> 2 for 40 bp
    //     [70, 80)  -> 1 for 10 bp
    //     [80, 100) -> 0 for 20 bp
    // - Raw stats:
    //     span_positions = eligible_positions = 100
    //     nonzero_positions = 60
    //     coverage_sum = 10*1 + 40*2 + 10*1 = 100
    //     coverage_sum_of_squares = 10*1^2 + 40*2^2 + 10*1^2 = 180
    // - Derived stats:
    //     mean = 100 / 100 = 1
    //     variance = 180 / 100 - 1^2 = 0.8
    //     sd = sqrt(0.8)
    //     cv = sqrt(0.8) / 1 = sqrt(0.8)
    //     covered_fraction = 60 / 100 = 0.6
    let bam = overlapping_fragment_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("overlapping_summary_window.bed");
    write_bed(&bed_path, &[("chr1", 0, 100, "variance_window")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::SummaryStats);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.summary_stats.tsv.zst"),
    )?;
    let rows = parse_tsv(&output);

    assert_eq!(rows.len(), 2);
    assert_eq!(rows[1][0..7], ["chr1", "0", "100", "100", "0", "100", "60"]);
    assert_close(rows[1][7].parse::<f64>()?, 100.0, 1e-9);
    assert_close(rows[1][8].parse::<f64>()?, 180.0, 1e-9);
    assert_close(rows[1][9].parse::<f64>()?, 1.0, 1e-9);
    assert_close(rows[1][10].parse::<f64>()?, 100.0, 1e-9);
    assert_close(rows[1][11].parse::<f64>()?, 0.8, 1e-9);
    assert_close(rows[1][12].parse::<f64>()?, 0.8_f64.sqrt(), 1e-6);
    assert_close(rows[1][13].parse::<f64>()?, 0.8_f64.sqrt(), 1e-6);
    assert_close(rows[1][14].parse::<f64>()?, 0.6, 1e-9);

    Ok(())
}

#[test]
fn by_bed_summary_stats_is_invariant_when_overlapping_windows_cross_tiles() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - Fragment 1 covers [20, 80) with coverage 1.
    // - Fragment 2 covers [30, 70) with another coverage 1.
    // - The resulting per-base coverage is:
    //     [0, 20)   -> 0
    //     [20, 30)  -> 1
    //     [30, 70)  -> 2
    //     [70, 80)  -> 1
    //     [80, 200) -> 0
    // - Window [0, 100):
    //     span_positions = eligible_positions = 100
    //     nonzero_positions = 60
    //     coverage_sum = 10*1 + 40*2 + 10*1 = 100
    //     coverage_sum_of_squares = 10*1^2 + 40*2^2 + 10*1^2 = 180
    //     mean = 100 / 100 = 1
    //     variance = 180 / 100 - 1^2 = 4 / 5
    //     sd = sqrt(4 / 5)
    //     cv = sqrt(4 / 5)
    //     covered_fraction = 60 / 100 = 3 / 5
    // - Window [20, 80):
    //     span_positions = eligible_positions = 60
    //     nonzero_positions = 60
    //     coverage_sum = 100
    //     coverage_sum_of_squares = 180
    //     mean = 100 / 60 = 5 / 3
    //     variance = 180 / 60 - (5 / 3)^2 = 2 / 9
    //     sd = sqrt(2) / 3
    //     cv = sqrt(2) / 5
    //     covered_fraction = 1
    // - Window [30, 70):
    //     span_positions = eligible_positions = nonzero_positions = 40
    //     coverage_sum = 80
    //     coverage_sum_of_squares = 160
    //     mean = 2
    //     variance = sd = cv = 0
    //     covered_fraction = 1
    // - `tile_size=33` forces all three windows through cross-tile BED summary reduction, while
    //   `tile_size=1000` keeps them in one tile. The final summary-stats rows must still match
    //   exactly because the raw moments are additive across tiles.
    let bam = overlapping_fragment_bam()?;
    let tile_sizes = [33_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let bed_path = out_dir
            .path()
            .join(format!("overlap_summary_{tile_size}.bed"));
        write_bed(
            &bed_path,
            &[
                ("chr1", 0, 100, "window_a"),
                ("chr1", 20, 80, "window_b"),
                ("chr1", 30, 70, "window_c"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_tile_size(tile_size);
        cfg.set_decimals(6);
        cfg.set_per_window(CoverageWindowAction::SummaryStats);
        cfg.set_windows(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
            by_grouped_bed: None,
        });

        run(&cfg)?;

        outputs.push(read_zst_to_string(
            &out_dir
                .path()
                .join("testcov.fcoverage.summary_stats.tsv.zst"),
        )?);
    }

    assert_eq!(
        outputs[0], outputs[1],
        "BED summary-stats output should not depend on tile size when windows cross tiles"
    );

    let rows = parse_tsv(&outputs[0]);
    assert_eq!(rows.len(), 4);
    assert_eq!(
        rows[0],
        vec![
            "chromosome",
            "start",
            "end",
            "span_positions",
            "blacklisted_positions",
            "eligible_positions",
            "nonzero_positions",
            "coverage_sum",
            "coverage_sum_of_squares",
            "average_coverage",
            "total_coverage",
            "variance_coverage",
            "sd_coverage",
            "coefficient_of_variation_coverage",
            "covered_fraction",
        ]
    );

    assert_eq!(rows[1][0..7], ["chr1", "0", "100", "100", "0", "100", "60"]);
    assert_close(rows[1][7].parse::<f64>()?, 100.0, 1e-9);
    assert_close(rows[1][8].parse::<f64>()?, 180.0, 1e-9);
    assert_close(rows[1][9].parse::<f64>()?, 1.0, 1e-9);
    assert_close(rows[1][10].parse::<f64>()?, 100.0, 1e-9);
    assert_close(rows[1][11].parse::<f64>()?, 4.0 / 5.0, 1e-9);
    assert_close(rows[1][12].parse::<f64>()?, (4.0_f64 / 5.0).sqrt(), 1e-6);
    assert_close(rows[1][13].parse::<f64>()?, (4.0_f64 / 5.0).sqrt(), 1e-6);
    assert_close(rows[1][14].parse::<f64>()?, 3.0 / 5.0, 1e-9);

    assert_eq!(rows[2][0..7], ["chr1", "20", "80", "60", "0", "60", "60"]);
    assert_close(rows[2][7].parse::<f64>()?, 100.0, 1e-9);
    assert_close(rows[2][8].parse::<f64>()?, 180.0, 1e-9);
    assert_close(rows[2][9].parse::<f64>()?, 5.0 / 3.0, 1e-6);
    assert_close(rows[2][10].parse::<f64>()?, 100.0, 1e-9);
    assert_close(rows[2][11].parse::<f64>()?, 2.0 / 9.0, 1e-6);
    assert_close(rows[2][12].parse::<f64>()?, 2.0_f64.sqrt() / 3.0, 1e-6);
    assert_close(rows[2][13].parse::<f64>()?, 2.0_f64.sqrt() / 5.0, 1e-6);
    assert_close(rows[2][14].parse::<f64>()?, 1.0, 1e-9);

    assert_eq!(rows[3][0..7], ["chr1", "30", "70", "40", "0", "40", "40"]);
    assert_close(rows[3][7].parse::<f64>()?, 80.0, 1e-9);
    assert_close(rows[3][8].parse::<f64>()?, 160.0, 1e-9);
    assert_close(rows[3][9].parse::<f64>()?, 2.0, 1e-9);
    assert_close(rows[3][10].parse::<f64>()?, 80.0, 1e-9);
    assert_close(rows[3][11].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[3][12].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[3][13].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[3][14].parse::<f64>()?, 1.0, 1e-9);

    Ok(())
}

#[test]
fn by_bed_total_keeps_coordinate_sorted_output_when_same_start_windows_cross_tiles() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - One fragment covers [20, 80) with coverage 1.
    // - The BED input deliberately lists the wider window before the narrower one:
    //     line 0 -> [0, 100)
    //     line 1 -> [0, 40)
    // - `load_windows_from_bed` preserves those original indices, but `Windows::new` then sorts
    //   by `(start, end)`, so tile processing sees:
    //     [0, 40)  with orig_idx 1
    //     [0, 100) with orig_idx 0
    // - Raw totals are:
    //     [0, 40)  -> 20
    //     [0, 100) -> 60
    // - With `tile_size=33`, both windows cross tile boundaries. This test pins that reduction
    //   still succeeds and that final output stays in coordinate order rather than original BED
    //   line order
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("same_start_reverse_bed_order.bed");
    write_bed(
        &bed_path,
        &[("chr1", 0, 100, "wide"), ("chr1", 0, 40, "narrow")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_tile_size(33);
    cfg.set_decimals(0);
    cfg.set_per_window(CoverageWindowAction::Total);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
        by_grouped_bed: None,
    });

    run(&cfg)?;

    let output = read_zst_to_string(&out_dir.path().join("testcov.fcoverage.total.tsv.zst"))?;
    let lines: Vec<_> = output.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t40\t20\t0",
            "chr1\t0\t100\t60\t0",
        ]
    );

    Ok(())
}

#[test]
fn grouped_summary_stats_on_unique_bases_writes_expected_rows() -> Result<()> {
    // Human verification status: unverified
    // Manual expectations:
    // - `beta` merges to [20, 80), so:
    //     span = eligible = 60
    //     nonzero = 60
    //     coverage_sum = coverage_sum_of_squares = 60
    //     mean = total / eligible = 1
    //     variance = sd = 0
    //     cv = 0
    //     covered_fraction = 1
    // - `gamma` is [0, 10) with no coverage, so:
    //     span = eligible = 10
    //     nonzero = 0
    //     coverage_sum = coverage_sum_of_squares = 0
    //     mean = variance = sd = 0
    //     cv = NaN because mean is zero
    //     covered_fraction = 0
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("grouped_summary_unique.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 40, 80, "beta"),
            ("chr1", 0, 10, "gamma"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::SummaryStatsOnUniqueBases);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.summary_stats_on_unique_bases.tsv.zst"),
    )?;
    let rows = parse_tsv(&output);

    assert_eq!(
        rows[0],
        vec![
            "group_idx",
            "span_positions",
            "blacklisted_positions",
            "eligible_positions",
            "nonzero_positions",
            "coverage_sum",
            "coverage_sum_of_squares",
            "average_coverage",
            "total_coverage",
            "variance_coverage",
            "sd_coverage",
            "coefficient_of_variation_coverage",
            "covered_fraction",
        ]
    );
    assert_eq!(rows[1][0..5], ["0", "60", "0", "60", "60"]);
    assert_close(rows[1][5].parse::<f64>()?, 60.0, 1e-9);
    assert_close(rows[1][6].parse::<f64>()?, 60.0, 1e-9);
    assert_close(rows[1][7].parse::<f64>()?, 1.0, 1e-9);
    assert_close(rows[1][8].parse::<f64>()?, 60.0, 1e-9);
    assert_close(rows[1][9].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[1][10].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[1][11].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[1][12].parse::<f64>()?, 1.0, 1e-9);

    assert_eq!(rows[2][0..5], ["1", "10", "0", "10", "0"]);
    assert_close(rows[2][5].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[2][6].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[2][7].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[2][8].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[2][9].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows[2][10].parse::<f64>()?, 0.0, 1e-9);
    assert!(rows[2][11].parse::<f64>()?.is_nan());
    assert_close(rows[2][12].parse::<f64>()?, 0.0, 1e-9);

    Ok(())
}

#[test]
fn grouped_summary_stats_with_global_row_treats_global_as_ordinary_site_weighted_row() -> Result<()>
{
    // Human verification status: unverified
    // Manual expectations:
    // - The fragment covers [20, 80) with per-base coverage 1.
    // - Group `beta` is loaded twice as [20, 50) and [40, 80). Plain grouped summary stats keep
    //   those intervals separate, so:
    //     span_positions = eligible_positions = 30 + 40 = 70
    //     nonzero_positions = 30 + 40 = 70
    //     coverage_sum = coverage_sum_of_squares = 70
    //     mean = total / eligible = 1
    // - Group `global` covers [0, 200):
    //     span_positions = eligible_positions = 200
    //     nonzero_positions = 60
    //     coverage_sum = coverage_sum_of_squares = 60
    //     mean = 60 / 200 = 0.3
    // - `global` is just another site-weighted grouped row here. It does not trigger any extra
    //   `fcoverage`-level correlation output.
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let grouped_bed = out_dir.path().join("grouped_summary_sites_with_global.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "beta"),
            ("chr1", 40, 80, "beta"),
            ("chr1", 0, 200, "global"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::SummaryStats);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.summary_stats.tsv.zst"),
    )?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    assert_eq!(
        rows_by_name
            .get("beta")
            .expect("beta row")
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>()[0..5],
        ["0", "70", "0", "70", "70"]
    );
    assert_close(rows_by_name["beta"][5].parse::<f64>()?, 70.0, 1e-9);
    assert_close(rows_by_name["beta"][6].parse::<f64>()?, 70.0, 1e-9);
    assert_close(rows_by_name["beta"][7].parse::<f64>()?, 1.0, 1e-9);
    assert_close(rows_by_name["beta"][8].parse::<f64>()?, 70.0, 1e-9);
    assert_close(rows_by_name["beta"][9].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["beta"][10].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["beta"][11].parse::<f64>()?, 0.0, 1e-9);
    assert_close(rows_by_name["beta"][12].parse::<f64>()?, 1.0, 1e-9);

    assert_eq!(
        rows_by_name
            .get("global")
            .expect("global row")
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>()[0..5],
        ["1", "200", "0", "200", "60"]
    );
    assert_close(rows_by_name["global"][5].parse::<f64>()?, 60.0, 1e-9);
    assert_close(rows_by_name["global"][6].parse::<f64>()?, 60.0, 1e-9);
    assert_close(rows_by_name["global"][7].parse::<f64>()?, 0.3, 1e-6);
    assert_close(rows_by_name["global"][8].parse::<f64>()?, 60.0, 1e-9);
    assert_close(rows_by_name["global"][9].parse::<f64>()?, 0.21, 1e-6);
    assert_close(
        rows_by_name["global"][10].parse::<f64>()?,
        0.21_f64.sqrt(),
        1e-6,
    );
    assert_close(
        rows_by_name["global"][11].parse::<f64>()?,
        0.21_f64.sqrt() / 0.3,
        1e-6,
    );
    assert_close(rows_by_name["global"][12].parse::<f64>()?, 0.3, 1e-6);

    Ok(())
}

#[test]
fn grouped_summary_stats_on_unique_bases_supports_downstream_pearson_with_gc_scaling_and_three_chromosomes()
-> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Build three 200 bp chromosomes with the repeating reference pattern "ACGT". Every 60 bp
    // fragment on that reference has exactly 50% GC, so the synthetic GC package below applies the
    // same correction weight 10.0 on all three chromosomes.
    //
    // One 60 bp fragment is placed at [20, 80) on each chromosome, and genomic scaling factors are
    // set to:
    // - chr1: 1.0
    // - chr2: 2.0
    // - chr3: 3.0
    //
    // Therefore the final per-base coverage is:
    // - chr1 [20,80): 10
    // - chr2 [20,80): 20
    // - chr3 [20,80): 30
    //
    // Group rows:
    // - `open` is [20, 50) on each chromosome, so under unique-base semantics:
    //     eligible_positions = 3 * 30 = 90
    //     coverage_sum = 30 * (10 + 20 + 30) = 1800
    //     coverage_sum_of_squares = 30 * (10^2 + 20^2 + 30^2) = 42000
    // - `global` is [0, 200) on each chromosome:
    //     eligible_positions = 3 * 200 = 600
    //     nonzero_positions = 3 * 60 = 180
    //     coverage_sum = 60 * (10 + 20 + 30) = 3600
    //     coverage_sum_of_squares = 60 * (10^2 + 20^2 + 30^2) = 84000
    //
    // `fcoverage` itself now stops at the summary statistics above. The test then derives Pearson R
    // downstream from the written `global` and `open` rows and checks that this matches the
    // ordinary direct positional formula over an explicit 600-position coverage vector and a 0/1
    // mask vector for the `open` row.
    let reference = fixtures::twobit_from_sequences(
        "fcoverage_grouped_summary_three_chr_reference",
        vec![
            ("chr1".to_string(), "ACGT".repeat(50)),
            ("chr2".to_string(), "ACGT".repeat(50)),
            ("chr3".to_string(), "ACGT".repeat(50)),
        ],
    )?;
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 200),
            ("chr2".to_string(), 200),
            ("chr3".to_string(), 200),
        ],
        vec![
            paired_fragment_on_tid(0, 20, 60, 20),
            paired_fragment_on_tid(1, 20, 60, 20),
            paired_fragment_on_tid(2, 20, 60, 20),
        ],
        Vec::new(),
        "fcoverage_grouped_summary_three_chr",
    )?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("three_chr_gc_pkg.npz");
    build_gc_package(&gc_path, 0, twobit_contig_footprint(&reference.path)?)?;
    let scaling_path = out_dir.path().join("three_chr_scaling.tsv");
    write_scaling_factors(
        &scaling_path,
        &[
            ("chr1", 0, 200, 1.0),
            ("chr2", 0, 200, 2.0),
            ("chr3", 0, 200, 3.0),
        ],
    )?;
    let grouped_bed = out_dir.path().join("grouped_summary_three_chr.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 20, 50, "open"),
            ("chr2", 20, 50, "open"),
            ("chr3", 20, 50, "open"),
            ("chr1", 0, 200, "global"),
            ("chr2", 0, 200, "global"),
            ("chr3", 0, 200, "global"),
        ],
    )?;

    let mut cfg = FCoverageConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1", "chr2", "chr3"]),
    );
    cfg.set_output_prefix("testcov");
    cfg.set_tile_size(1_000);
    cfg.set_ignore_gap(false);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_decimals(6);
    cfg.set_per_window(CoverageWindowAction::SummaryStatsOnUniqueBases);
    cfg.set_windows(DistributionWindowsArgs {
        by_size: None,
        by_bed: None,
        by_grouped_bed: Some(grouped_bed),
    });
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.set_scale_genome(ScaleGenomeArgs {
        scaling_factors: Some(scaling_path),
    });
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    run(&cfg)?;

    let output = read_zst_to_string(
        &out_dir
            .path()
            .join("testcov.fcoverage.summary_stats_on_unique_bases.tsv.zst"),
    )?;
    let sidecar = std::fs::read_to_string(out_dir.path().join("testcov.group_index.tsv"))?;
    let rows_by_name = grouped_rows_by_name(&output, &sidecar);

    let mut direct_coverage = vec![0.0_f64; 600];
    let mut direct_mask = vec![0.0_f64; 600];
    for (chromosome_offset, scaled_coverage) in [(0_usize, 10.0_f64), (200, 20.0), (400, 30.0)] {
        for position in 20..80 {
            direct_coverage[chromosome_offset + position] = scaled_coverage;
        }
        for position in 20..50 {
            direct_mask[chromosome_offset + position] = 1.0;
        }
    }
    let expected_pearson = pearson_r_from_vectors(&direct_coverage, &direct_mask);

    assert_eq!(
        rows_by_name["open"][0..5]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["0", "90", "0", "90", "90"]
    );
    assert_close(rows_by_name["open"][5].parse::<f64>()?, 1800.0, 1e-9);
    assert_close(rows_by_name["open"][6].parse::<f64>()?, 42000.0, 1e-6);
    assert_close(rows_by_name["open"][7].parse::<f64>()?, 20.0, 1e-9);
    assert_close(rows_by_name["open"][8].parse::<f64>()?, 1800.0, 1e-9);
    assert_close(rows_by_name["open"][9].parse::<f64>()?, 200.0 / 3.0, 1e-6);
    assert_close(
        rows_by_name["open"][10].parse::<f64>()?,
        (200.0_f64 / 3.0_f64).sqrt(),
        1e-6,
    );
    assert_close(
        rows_by_name["open"][11].parse::<f64>()?,
        (200.0_f64 / 3.0_f64).sqrt() / 20.0,
        1e-6,
    );
    assert_close(rows_by_name["open"][12].parse::<f64>()?, 1.0, 1e-9);

    assert_eq!(
        rows_by_name["global"][0..5]
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        vec!["1", "600", "0", "600", "180"]
    );
    assert_close(rows_by_name["global"][5].parse::<f64>()?, 3600.0, 1e-9);
    assert_close(rows_by_name["global"][6].parse::<f64>()?, 84000.0, 1e-6);
    assert_close(rows_by_name["global"][7].parse::<f64>()?, 6.0, 1e-9);
    assert_close(rows_by_name["global"][8].parse::<f64>()?, 3600.0, 1e-9);
    assert_close(rows_by_name["global"][9].parse::<f64>()?, 104.0, 1e-9);
    assert_close(
        rows_by_name["global"][10].parse::<f64>()?,
        104.0_f64.sqrt(),
        1e-6,
    );
    assert_close(
        rows_by_name["global"][11].parse::<f64>()?,
        104.0_f64.sqrt() / 6.0,
        1e-6,
    );
    assert_close(rows_by_name["global"][12].parse::<f64>()?, 0.3, 1e-9);

    let derived_pearson =
        pearson_r_from_summary_stats_rows(&rows_by_name["global"], &rows_by_name["open"])?;
    assert_close(derived_pearson, expected_pearson, 1e-6);

    Ok(())
}
