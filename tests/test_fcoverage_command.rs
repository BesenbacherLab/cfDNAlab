#![cfg(feature = "cmd_fcoverage")]

mod fixtures;

use anyhow::Result;
#[cfg(feature = "cmd_bam_to_bam")]
use cfdnalab::commands::bam_to_bam::{
    bam_to_bam::run_inner as run_bam_to_bam, config::BamToBamConfig,
};
use cfdnalab::commands::cli_common::{
    ApplyGCArgs, AssignToWindowArgs, ChromosomeArgs, IOCArgs, ScaleGenomeArgs, WindowsArgs,
};
#[cfg(feature = "cmd_coverage_weights")]
use cfdnalab::commands::coverage_weights::{
    config::CoverageWeightsConfig, coverage_weights::run as run_coverage_weights,
};
use cfdnalab::commands::fcoverage::config::FCoverageConfig;
use cfdnalab::commands::fcoverage::fcoverage::run;
use cfdnalab::commands::fcoverage::window_results::CoverageWindowAction;
use cfdnalab::commands::gc_bias::{GC_CORRECTION_SCHEMA_VERSION, package::GCCorrectionPackage};
use cfdnalab::commands::lengths::{config::LengthsConfig, lengths::run as run_lengths};
use cfdnalab::shared::fragment::minimal_fragment::{MinimalReadInfo, collect_fragment};
use cfdnalab::shared::indel_mode::IndelMode;
use cfdnalab::shared::io::dot_join;
use cfdnalab::shared::read::default_include_read_paired_end;
use fixtures::{
    BamFixture, FragmentSpec, LONG_FRAGMENT_LENGTH, LONG_FRAGMENT_STARTS, ReadSpec, bam_from_specs,
    bam_from_specs_strict_identity, build_real_neutral_gc_package,
    build_real_non_neutral_gc_package, long_fragment_bam, paired_fragment, read_zst_to_string,
    simple_inward_bam, simple_reference_twobit, write_bed, write_scaling_factors,
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

fn build_gc_package(path: &Path, end_offset: u64) -> Result<()> {
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset,
        length_edges: vec![10, 60, 200],
        gc_edges: vec![0, 50, 101],
        length_bin_frequencies: array![1.0_f64, 3.0_f64],
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

    let mut windows = WindowsArgs::default();
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
        cfg.set_windows(WindowsArgs {
            by_size: Some(200),
            by_bed: None,
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
    lengths_cfg.set_windows(WindowsArgs::default());
    lengths_cfg.set_window_assignment(AssignToWindowArgs::default());
    lengths_cfg.set_min_mapq(30);
    lengths_cfg.set_require_proper_pair(false);
    {
        let frag = lengths_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

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

        let mut windows = WindowsArgs::default();
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

    let mut windows = WindowsArgs::default();
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

    let output_path = out_dir.path().join("testcov.fcoverage.avg.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
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
        let mut by_size_windows = WindowsArgs::default();
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
    by_size_cfg.set_windows(WindowsArgs {
        by_size: Some(200),
        by_bed: None,
    });

    let mut by_bed_cfg = base_config(&bam.bam, by_bed_out.path());
    by_bed_cfg.set_decimals(0);
    by_bed_cfg.set_per_window(CoverageWindowAction::Total);
    by_bed_cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
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

        let mut windows = WindowsArgs::default();
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
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(bed_path),
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

        let mut windows = WindowsArgs::default();
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

    let mut windows = WindowsArgs::default();
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
    let _header = lines.next().unwrap_or("");
    let first = lines.next().unwrap_or("");
    assert!(first.starts_with("chr1\t0\t40"));
    assert!(
        first.ends_with("20\t0"),
        "expected total coverage 20, got: {first}"
    );
    let second = lines.next().unwrap_or("");
    assert!(second.starts_with("chr1\t40\t80"));
    assert!(
        second.ends_with("40\t0"),
        "expected total coverage 40, got: {second}"
    );

    Ok(())
}

#[test]
fn by_size_average_reduces_across_non_aligned_tiles() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;

    let mut windows = WindowsArgs::default();
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

    let output_path = out_dir.path().join("testcov.fcoverage.avg.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
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
        build_gc_package(&gc_path, 0)?;

        let mut scale_genome = ScaleGenomeArgs::default();
        scale_genome.scaling_factors = Some(scaling_path);

        let mut windows = WindowsArgs::default();
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
            skip_invalid_gc: false,
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

    let mut windows = WindowsArgs::default();
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

    let output_path = out_dir.path().join("testcov.fcoverage.avg.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
            "chr1\t0\t40\t0.5\t0",
            "chr1\t20\t80\t1\t0",
            "chr1\t70\t90\t0.5\t0",
        ]
    );

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

    let mut windows = WindowsArgs::default();
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

    let output_path = out_dir.path().join("testcov.fcoverage.avg.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
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
        "fcoverage_bed_avg_skip_empty_chr",
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("aggregate_windows_chr2_only.bed");
    write_bed(&bed_path, &[("chr2", 0, 40, "chr2_window")])?;

    let mut windows = WindowsArgs::default();
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

    let output_path = out_dir.path().join("testcov.fcoverage.avg.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\tavg_coverage\tblacklisted_positions",
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

    let mut windows = WindowsArgs::default();
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

    let mut windows = WindowsArgs::default();
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

    let mut windows = WindowsArgs::default();
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

    let mut windows = WindowsArgs::default();
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
        let mut windows = WindowsArgs::default();
        windows.by_size = Some(40);

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_per_window(action);
        cfg.set_windows(windows);

        let err = run(&cfg).expect_err("by-size positional window mode should fail");
        let err_text = err.to_string();
        assert!(
            err_text.contains("in --by-size mode, --per-window can only be 'average' or 'total'"),
            "unexpected error for {action:?}: {err_text}"
        );
    }

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
    let scaling_path = weights_out_dir.join("coverage.scaling_factors.tsv");

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
    let weights_cfg = make_simple_coverage_weights_config(&weights_out_dir, &bam.bam);
    let gc_path = build_real_neutral_gc_package(&bam.bam, &ref_twobit.path, out_dir.path(), 60)?;

    run_coverage_weights(&weights_cfg)?;
    let scaling_path = weights_out_dir.join("coverage.scaling_factors.tsv");

    let mut scale_genome = ScaleGenomeArgs::default();
    scale_genome.scaling_factors = Some(scaling_path);

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(6);
    cfg.set_keep_zero_runs(false);
    cfg.set_scale_genome(scale_genome);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        skip_invalid_gc: false,
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
    // `normalize_avg_overlap_keeps_sparse_non_zero_scaling_finite` already pins the helper-level
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
    // Written scaling factors = mean / avg_pos_cov:
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
    let scaling_path = weights_out_dir.join("coverage.scaling_factors.tsv");

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
        skip_invalid_gc: false,
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
        skip_invalid_gc: false,
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
        skip_invalid_gc: false,
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
        skip_invalid_gc: false,
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
        skip_invalid_gc: false,
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
            "invalid_mate_falls_back",
            &[Some(2.0), Some(2_000.0)][..],
            false,
            vec!["chr1\t20\t80\t1"],
        ),
        (
            "zero_mate_forces_zero_weight",
            &[Some(0.0), Some(4.0)][..],
            true,
            vec!["chr1\t0\t200\t0"],
        ),
    ];

    for (name, tags, keep_zero_runs, expected_lines) in scenarios {
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
            skip_invalid_gc: false,
        });

        // Manual expectations:
        // - Scenario invalid_mate_falls_back:
        //   mate tags 2.0 and 2000.0 make the fragment GC tag invalid.
        //   With skip_invalid_gc=false, fcoverage falls back to weight 1.0 -> [20, 80) with value 1.
        // - Scenario zero_mate_forces_zero_weight:
        //   mate tags 0.0 and 4.0 combine to fragment weight 0.0.
        //   With keep_zero_runs=true, the whole chromosome is a single zero-coverage segment.
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
fn gc_tag_missing_or_invalid_values_fall_back_or_drop() -> Result<()> {
    // Human verification status: unverified
    let scenarios = [
        ("missing", None, false, vec!["chr1\t20\t80\t1"]),
        ("missing_drop", None, true, Vec::<&str>::new()),
        (
            "out_of_range",
            Some(2_000.0),
            false,
            vec!["chr1\t20\t80\t1"],
        ),
        ("out_of_range_drop", Some(2_000.0), true, Vec::<&str>::new()),
    ];

    for (name, tag_value, skip_invalid_gc, expected_lines) in scenarios {
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
            skip_invalid_gc,
        });

        // Manual expectations:
        // - Missing tags and out-of-range tags both produce no usable GC weight.
        // - With skip_invalid_gc=false, fcoverage falls back to weight 1.0.
        // - With skip_invalid_gc=true, the fragment is skipped entirely.
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
    build_gc_package(&gc_path, 0)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        skip_invalid_gc: false,
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
    build_gc_package(&gc_path, 0)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_decimals(0);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        skip_invalid_gc: false,
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
    // With command fragment-length bounds set to exactly [60, 60], the consumer must reject the
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
        correction_matrix: array![[1.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        skip_invalid_gc: false,
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
        correction_matrix: array![[1.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        skip_invalid_gc: false,
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
        skip_invalid_gc: false,
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
        skip_invalid_gc: false,
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
fn gc_file_drop_invalid_controls_short_effective_length_fragments() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let scenarios = [
        ("fallback", false, vec!["chr1\t20\t80\t1"]),
        ("drop", true, Vec::<&str>::new()),
    ];

    for (name, skip_invalid_gc, expected_lines) in scenarios {
        let out_dir = TempDir::new()?;
        let gc_path = out_dir.path().join(format!("gc_pkg_{name}.npz"));
        build_gc_package(&gc_path, 26)?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.set_gc(ApplyGCArgs {
            gc_file: Some(gc_path),
            gc_tag: None,
            skip_invalid_gc,
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
        // - With skip_invalid_gc=false, the fragment falls back to weight 1.0.
        // - With skip_invalid_gc=true, the fragment is skipped.
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

    let mut windows = WindowsArgs::default();
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

    let mut windows = WindowsArgs::default();
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
