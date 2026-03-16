#![cfg(feature = "cmd_fcoverage")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{
    ApplyGCArgs, ChromosomeArgs, IOCArgs, ScaleGenomeArgs, WindowsArgs,
};
use cfdnalab::commands::fcoverage::config::FCoverageConfig;
use cfdnalab::commands::fcoverage::fcoverage::run;
use cfdnalab::commands::fcoverage::window_results::CoverageWindowAction;
use cfdnalab::commands::gc_bias::{GC_CORRECTION_SCHEMA_VERSION, package::GCCorrectionPackage};
use cfdnalab::shared::fragment::minimal_fragment::collect_fragment_from_records;
use cfdnalab::shared::read::default_include_read_paired_end;
use fixtures::{
    BamFixture, FragmentSpec, LONG_FRAGMENT_LENGTH, LONG_FRAGMENT_STARTS, ReadSpec, bam_from_specs,
    long_fragment_bam, paired_fragment, read_zst_to_string, simple_inward_bam,
    simple_reference_twobit, write_bed, write_scaling_factors,
};
use ndarray::array;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read, Reader};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

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
    let frag = collect_fragment_from_records(&pair_store[0], &pair_store[1]);
    assert!(frag.is_some(), "expected fragment collection to succeed");

    run(&cfg)?;

    let bedgraph = out_dir.path().join("testcov.per_position.bedgraph.zst");
    assert!(bedgraph.exists(), "expected positional bedgraph output");
    let text = read_zst_to_string(&bedgraph)?;
    assert!(
        text.contains("chr1\t20\t80\t1"),
        "expected contiguous coverage run, got: {text}"
    );

    Ok(())
}

#[test]
fn per_position_keep_zero_runs_toggles_zero_segments() -> Result<()> {
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
        .join("testcov.per_position.bedgraph.zst");
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
        .join("testcov.per_position.bedgraph.zst");
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

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t40\t1", "chr1\t60\t80\t1"]);

    Ok(())
}

#[test]
fn unpaired_single_read_matches_fragment_span_output() -> Result<()> {
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

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t80\t1"]);

    Ok(())
}

#[test]
fn unpaired_mode_rejects_ignore_gap() -> Result<()> {
    let bam = single_read_fragment_bam("fcoverage_unpaired_ignore_gap")?;
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
    let bam = single_read_fragment_bam("fcoverage_unpaired_require_pp")?;
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

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t30\t1", "chr1\t35\t80\t1"]);

    Ok(())
}

#[test]
fn blacklist_masks_positions_in_positional_output_across_tile_boundary() -> Result<()> {
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

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
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

    let output_path = out_dir.path().join("testcov.total.tsv.zst");
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
fn blacklist_crossing_tile_boundary_keeps_same_by_size_output() -> Result<()> {
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

        let output_path = out_dir.path().join("testcov.total.tsv.zst");
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

    let output_path = out_dir.path().join("testcov.avg.tsv.zst");
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
            .join("testcov.per_position.bedgraph.zst");
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

        let totals_path = by_size_out_dir.path().join("testcov.total.tsv.zst");
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
fn by_bed_total_is_invariant_when_windows_cross_tile_boundaries() -> Result<()> {
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

        let output_path = out_dir.path().join("testcov.total.tsv.zst");
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
fn by_size_total_counts_covered_bases_per_window() -> Result<()> {
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

    let totals = out_dir.path().join("testcov.total.tsv.zst");
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

    let output_path = out_dir.path().join("testcov.avg.tsv.zst");
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
fn by_bed_average_matches_manual_window_means() -> Result<()> {
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

    let output_path = out_dir.path().join("testcov.avg.tsv.zst");
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

    let output_path = out_dir.path().join("testcov.avg.tsv.zst");
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
fn by_bed_total_matches_manual_window_sums() -> Result<()> {
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

    let output_path = out_dir.path().join("testcov.total.tsv.zst");
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

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t40\t1", "chr1\t70\t80\t1"]);

    Ok(())
}

#[test]
fn by_bed_indexed_positions_keep_window_indices_and_overlap_duplicates() -> Result<()> {
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
        .join("testcov.per_position_per_window.tsv.zst");
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

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t50\t1.2", "chr1\t50\t80\t1.3"]);

    Ok(())
}

#[test]
fn gc_tag_weights_unpaired_positional_output() -> Result<()> {
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
        drop_invalid_gc: false,
    });

    // Manual expectations:
    // - The single read spans [20, 80).
    // - In unpaired mode, its GC tag weight is applied directly to the fragment.
    // - Coverage is therefore 2.5 across the whole span.
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t80\t2.5"]);

    Ok(())
}

#[test]
fn gc_tag_averages_valid_mate_weights_in_paired_mode() -> Result<()> {
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
        drop_invalid_gc: false,
    });

    // Manual expectations:
    // - The paired fragment spans [20, 80).
    // - Mate GC tags are 2.0 and 4.0.
    // - Paired GC-tag mode combines them as their average, so the fragment weight is 3.0.
    // - Coverage is therefore 3 across the full fragment span.
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t80\t3"]);

    Ok(())
}

#[test]
fn gc_tag_paired_edge_cases_follow_fragment_combination_rules() -> Result<()> {
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
            drop_invalid_gc: false,
        });

        // Manual expectations:
        // - Scenario invalid_mate_falls_back:
        //   mate tags 2.0 and 2000.0 make the fragment GC tag invalid.
        //   With drop_invalid_gc=false, fcoverage falls back to weight 1.0 -> [20, 80) with value 1.
        // - Scenario zero_mate_forces_zero_weight:
        //   mate tags 0.0 and 4.0 combine to fragment weight 0.0.
        //   With keep_zero_runs=true, the whole chromosome is a single zero-coverage segment.
        run(&cfg)?;

        let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
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

    for (name, tag_value, drop_invalid_gc, expected_lines) in scenarios {
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
            drop_invalid_gc,
        });

        // Manual expectations:
        // - Missing tags and out-of-range tags both produce no usable GC weight.
        // - With drop_invalid_gc=false, fcoverage falls back to weight 1.0.
        // - With drop_invalid_gc=true, the fragment is skipped entirely.
        run(&cfg)?;

        let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
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
    let bam = simple_inward_bam()?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("gc_pkg.npz");
    build_gc_package(&gc_path, 0)?;

    let mut cfg = base_config(&bam.bam, out_dir.path());
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        drop_invalid_gc: false,
    });

    let err = run(&cfg).expect_err("GC correction should require --ref-2bit");
    let msg = format!("{err:#}");
    assert!(msg.contains("--ref-2bit"), "unexpected error: {msg}");

    Ok(())
}

#[test]
fn gc_file_weights_positional_output_from_reference_package() -> Result<()> {
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
        drop_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

    // Manual expectations:
    // - The fragment length is 60, which lands in the [60, 200) length bin.
    // - The simple reference has 50% GC over [20, 80), which lands in the [50, 101) GC bin.
    // - The test GC package assigns weight 10.0 to that bin, so coverage is 10 on [20, 80).
    run(&cfg)?;

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines, vec!["chr1\t20\t80\t10"]);

    Ok(())
}

#[test]
fn gc_file_drop_invalid_controls_short_effective_length_fragments() -> Result<()> {
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let scenarios = [
        ("fallback", false, vec!["chr1\t20\t80\t1"]),
        ("drop", true, Vec::<&str>::new()),
    ];

    for (name, drop_invalid_gc, expected_lines) in scenarios {
        let out_dir = TempDir::new()?;
        let gc_path = out_dir.path().join(format!("gc_pkg_{name}.npz"));
        build_gc_package(&gc_path, 26)?;

        let mut cfg = base_config(&bam.bam, out_dir.path());
        cfg.set_decimals(0);
        cfg.set_gc(ApplyGCArgs {
            gc_file: Some(gc_path),
            gc_tag: None,
            drop_invalid_gc,
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
        // - With drop_invalid_gc=false, the fragment falls back to weight 1.0.
        // - With drop_invalid_gc=true, the fragment is skipped.
        run(&cfg)?;

        let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
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

    let output_path = out_dir.path().join("testcov.per_position.bedgraph.zst");
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
        .join("testcov.per_position_per_window.tsv.zst");
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
