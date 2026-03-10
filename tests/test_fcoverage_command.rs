#![cfg(feature = "cmd_fcoverage")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs, ScaleGenomeArgs, WindowsArgs};
use cfdnalab::commands::fcoverage::config::FCoverageConfig;
use cfdnalab::commands::fcoverage::fcoverage::run;
use cfdnalab::commands::fcoverage::window_results::CoverageWindowAction;
use cfdnalab::shared::fragment::minimal_fragment::collect_fragment_from_records;
use cfdnalab::shared::read::default_include_read_paired_end;
use fixtures::{
    LONG_FRAGMENT_LENGTH, LONG_FRAGMENT_STARTS, bam_from_specs, long_fragment_bam, paired_fragment,
    read_zst_to_string, simple_inward_bam, write_bed, write_scaling_factors,
};
use rust_htslib::bam::{Read, Reader};
use std::path::Path;
use tempfile::TempDir;

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
fn by_size_total_and_average_outputs() -> Result<()> {
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
