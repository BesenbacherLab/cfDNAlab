#![cfg(feature = "cmd_wps")]

mod fixtures;

use anyhow::{Context, Result, ensure};
use cfdnalab::RunOptions;
use cfdnalab::reference::twobit_contig_footprint;
use cfdnalab::run_like_cli::common::{ApplyGCArgs, ChromosomeArgs, IOCArgs, WindowsArgs};
use cfdnalab::run_like_cli::fcoverage::CoverageWindowAction;
use cfdnalab::run_like_cli::wps::{WPSConfig, run_wps as run_fn};
use cfdnalab::testing::{
    Bed4Row, Cigar, FragmentSpec, ReadSpec, TempBam, bam_from_fragments,
    long_inward_fragment_series_bam, read_zst_to_string, twobit_from_sequences, write_bed4,
    write_two_bin_gc_correction_package,
};
use fixtures::late_origin_gc_reference_sequence;
use std::cmp::max;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use tempfile::TempDir;
use zstd::stream::read::Decoder as ZstdDecoder;

const EPSILON: f32 = 1e-6;
const FLAG_FIRST_MATE: u16 = 0x40;
const FLAG_SECOND_MATE: u16 = 0x80;
const FLAG_PROPER_PAIR: u16 = 0x2;
const FLAG_MATE_REVERSE: u16 = 0x20;

fn run_wps_quiet(cfg: &WPSConfig) -> Result<()> {
    run_fn(cfg, RunOptions::new_quiet()).map(|_| ())
}
const WPS_WINDOW_SIZE_BP: u32 = 120;

#[derive(Debug, Clone, PartialEq)]
struct WpsRun {
    chromosome: String,
    start: u32,
    end: u32,
    value: f32,
}

fn wps_run(chr: &str, start: u32, end: u32, value: f32) -> WpsRun {
    WpsRun {
        chromosome: chr.to_string(),
        start,
        end,
        value,
    }
}

fn fragment_spec(start: u32, end: u32) -> FragmentSpec {
    let length = end - start;
    let read_len = max(length / 2, 2);
    let forward_pos = start as i64;
    let reverse_pos = end as i64 - read_len as i64;
    let fragment_span = end as i64 - start as i64;

    let forward = ReadSpec {
        tid: 0,
        pos: forward_pos,
        cigar: vec![Cigar::Match(read_len)],
        seq: vec![b'A'; read_len as usize],
        base_quality: 40,
        is_reverse: false,
        mapq: 60,
        flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
        mate_tid: Some(0),
        mate_pos: Some(reverse_pos),
        insert_size: fragment_span,
    };

    let reverse = ReadSpec {
        tid: 0,
        pos: reverse_pos,
        cigar: vec![Cigar::Match(read_len)],
        seq: vec![b'T'; read_len as usize],
        base_quality: 40,
        is_reverse: true,
        mapq: 60,
        flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
        mate_tid: Some(0),
        mate_pos: Some(forward_pos),
        insert_size: -fragment_span,
    };

    FragmentSpec { forward, reverse }
}

fn fragment_spec_on_tid(tid: usize, start: u32, end: u32) -> FragmentSpec {
    let mut fragment = fragment_spec(start, end);
    fragment.forward.tid = tid;
    fragment.reverse.tid = tid;
    fragment.forward.mate_tid = Some(tid);
    fragment.reverse.mate_tid = Some(tid);
    fragment
}

fn make_fixture(name: &str, fragments: &[(u32, u32)]) -> Result<TempBam> {
    let chrom_len = fragments
        .iter()
        .map(|(_, end)| end + 100)
        .max()
        .unwrap_or(500);
    let specs: Vec<FragmentSpec> = fragments
        .iter()
        .map(|(start, end)| fragment_spec(*start, *end))
        .collect();
    bam_from_fragments(
        name,
        vec![("chr1".to_string(), chrom_len)],
        specs,
        Vec::new(),
    )
}

fn make_three_chrom_fixture(name: &str, fragments: &[(u32, u32)]) -> Result<TempBam> {
    let chrom_len = 100u32;
    let specs: Vec<FragmentSpec> = fragments
        .iter()
        .enumerate()
        .map(|(tid, (start, end))| fragment_spec_on_tid(tid, *start, *end))
        .collect();
    bam_from_fragments(
        name,
        vec![
            ("chr1".to_string(), chrom_len),
            ("chr2".to_string(), chrom_len),
            ("chr3".to_string(), chrom_len),
        ],
        specs,
        Vec::new(),
    )
}

fn make_config(
    window_size: u32,
    keep_zero_runs: bool,
    bam_path: &Path,
    output_dir: &Path,
    prefix: &str,
) -> WPSConfig {
    let ioc = IOCArgs {
        bam: bam_path.to_path_buf(),
        output_dir: output_dir.to_path_buf(),
        n_threads: 1,
    };
    let chroms = ChromosomeArgs {
        chromosomes: Some(vec!["chr1".to_string()]),
        chromosomes_file: None,
    };

    let mut cfg = WPSConfig::new(
        ioc,
        chroms,
        Some(CoverageWindowAction::OnlyIncludeThesePositionsUnique), // No genomic windowing so doesn't currently matter
    );
    cfg.set_output_prefix(prefix.to_string());
    cfg.set_window_size(window_size);
    cfg.set_keep_zero_runs(keep_zero_runs);
    cfg.set_min_fragment_length(window_size);
    cfg.set_max_fragment_length(window_size + 200);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_decimals(0);
    cfg.set_tile_size(1_000);
    cfg
}

fn run_wps(cfg: &WPSConfig) -> Result<Vec<WpsRun>> {
    run_wps_with_chrom(cfg)
}

fn run_wps_with_chrom(cfg: &WPSConfig) -> Result<Vec<WpsRun>> {
    run_wps_quiet(cfg)?;
    let prefix = cfg.shared_args.output_prefix.trim();
    let bedgraph_path = cfg
        .shared_args
        .ioc
        .output_dir
        .join(format!("{prefix}.wps.per_position.bedgraph.zst"));

    let file = File::open(&bedgraph_path)
        .with_context(|| format!("opening WPS output {}", bedgraph_path.display()))?;
    let decoder =
        ZstdDecoder::new(file).with_context(|| format!("decoding {}", bedgraph_path.display()))?;
    let reader = BufReader::new(decoder);

    let mut runs = Vec::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("reading {}", bedgraph_path.display()))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut cols = line.split('\t');
        let chromosome = cols
            .next()
            .context("missing chromosome column in WPS output")?;
        let start_str = cols.next().context("missing start column in WPS output")?;
        let end_str = cols.next().context("missing end column in WPS output")?;
        let value_str = cols.next().context("missing value column in WPS output")?;

        let start = start_str
            .parse::<u32>()
            .with_context(|| format!("invalid start value '{start_str}'"))?;
        let end = end_str
            .parse::<u32>()
            .with_context(|| format!("invalid end value '{end_str}'"))?;
        let value = value_str
            .parse::<f32>()
            .with_context(|| format!("invalid value '{value_str}'"))?;

        runs.push(WpsRun {
            chromosome: chromosome.to_string(),
            start,
            end,
            value,
        });
    }

    Ok(runs)
}

#[test]
fn gc_file_late_tile_window_uses_reference_coordinates_after_fetch_narrowing() -> Result<()> {
    // Arrange:
    // - The fragment spans [900,961), and the reported WPS centers are restricted to [925,936).
    // - The reference is shorter than the BAM chromosome, but long enough for the narrowed
    //   window-derived fetch span. Reading the full tile reference would overrun the reference.
    // - The fragment interval [900,961) is all C, so it lands in the high-GC correction bin with
    //   weight 7.0. Using prefix-local origin 0 would see A-only sequence instead.
    // - A 10 bp WPS window around every center in [925,936) lies fully inside the fragment and
    //   contains neither endpoint, so the unweighted score is +1 throughout.
    // - The GC package therefore makes the full run +7.
    let bam = make_fixture("wps_late_tile_gc_origin", &[(900, 961)])?;
    let reference = twobit_from_sequences(
        "wps_late_tile_gc_origin_ref",
        vec![("chr1".to_string(), late_origin_gc_reference_sequence())],
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("late_window.bed");
    let gc_path = out_dir.path().join("two_bin_gc_package.zarr");
    write_bed4(&bed_path, &[Bed4Row::new("chr1", 925, 936, "late")])?;
    write_two_bin_gc_correction_package(
        &gc_path,
        61,
        2.0,
        7.0,
        twobit_contig_footprint(&reference.path)?,
    )?;

    let mut cfg = make_config(10, false, &bam.bam, out_dir.path(), "late_gc");
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
    });
    cfg.set_min_fragment_length(61);
    cfg.set_max_fragment_length(61);
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act
    let actual = run_wps(&cfg)?;

    // Assert
    assert_runs_equal(&actual, &[wps_run("chr1", 925, 936, 7.0)]);
    Ok(())
}

#[test]
fn single_fragment_produces_central_plateau() -> Result<()> {
    let fixture = make_fixture("wps_single_fragment", &[(10, 22)])?;
    let out_dir = TempDir::new()?;
    let cfg = make_config(4, false, &fixture.bam, out_dir.path(), "single_fragment");

    // Manual expectations:
    // - Window size 4 gives left_span = right_span = 2. A fragment counts as fully covering a
    //   center when the window [c-2, c+2) stays within [10, 22).
    //   This happens for c = 12..=20, yielding the +1 plateau [12, 21).
    // - Endpoints only subtract when they fall strictly inside the window:
    //   * The left endpoint at 10 affects centers 9, 10, 11 -> run [9, 12) at -1.
    //   * The right endpoint at 22 affects centers 21, 22 -> run [21, 23) at -1.
    // - All remaining centers stay at zero and are omitted because keep_zero_runs=false.
    let expected = vec![
        wps_run("chr1", 9, 12, -1.0),
        wps_run("chr1", 12, 21, 1.0),
        wps_run("chr1", 21, 23, -1.0),
    ];

    let actual = run_wps(&cfg)?;

    assert_runs_equal(&actual, &expected);
    Ok(())
}

#[test]
fn overlapping_fragments_stack_scores() -> Result<()> {
    let fixture = make_fixture("wps_overlapping", &[(0, 20), (4, 12)])?;
    let out_dir = TempDir::new()?;
    let cfg = make_config(4, false, &fixture.bam, out_dir.path(), "overlapping");

    // Manual expectations for two fragments:
    // - Fragment F1: [0, 20), fragment F2: [4, 12); window size 4 keeps left_span = right_span = 2.
    // - Full-span contributions:
    //     * F1 covers c = 2..=18.
    //     * F2 covers c = 6..=10.
    // - Endpoint penalties:
    //     * F1 endpoints reduce centers c = 19, 20.
    //     * F2 endpoints reduce centers c = 3, 4, 5 on the left and c = 11, 12 on the right.
    // - Combining both fragments yields the visible runs:
    //     * [2, 3) at +1 from the long fragment.
    //     * [6, 11) at +2 where both fragments span the window.
    //     * [13, 19) at +1 once only the long fragment remains.
    //     * [19, 21) at -1 from the long fragment's right endpoint.
    let expected = vec![
        wps_run("chr1", 2, 3, 1.0),
        wps_run("chr1", 6, 11, 2.0),
        wps_run("chr1", 13, 19, 1.0),
        wps_run("chr1", 19, 21, -1.0),
    ];

    let actual = run_wps(&cfg)?;

    assert_runs_equal(&actual, &expected);
    Ok(())
}

#[test]
fn keep_zero_runs_emits_flat_segments() -> Result<()> {
    let fixture = make_fixture("wps_keep_zero", &[(10, 22)])?;
    let out_dir = TempDir::new()?;
    let cfg = make_config(4, true, &fixture.bam, out_dir.path(), "keep_zero");

    // Same geometry as the single-fragment test, but keep_zero_runs=true means we retain zero
    // plateaus between non-zero segments:
    // - Leading zeros from the first valid center (2) up to the first penalty at 9.
    // - Trailing zeros that stretch beyond the region of interest; we assert up to c = 30.
    let expected = vec![
        wps_run("chr1", 2, 9, 0.0),
        wps_run("chr1", 9, 12, -1.0),
        wps_run("chr1", 12, 21, 1.0),
        wps_run("chr1", 21, 23, -1.0),
        wps_run("chr1", 23, 30, 0.0),
    ];

    let actual = run_wps(&cfg)?;

    let clipped = clip_runs(&actual, 30);
    assert_runs_equal(&clipped, &expected);
    Ok(())
}

#[test]
fn fragment_equal_to_window_removes_central_signal() -> Result<()> {
    let fixture = make_fixture("wps_equal_window", &[(10, 14)])?;
    let out_dir = TempDir::new()?;
    let cfg = make_config(4, false, &fixture.bam, out_dir.path(), "equal_window");

    // Fragment length exactly matches the window (4 bp):
    // - Center c = 12 is fully covered (edge-aligned) so contributes +1.
    // - Endpoints reduce windows that contain them strictly:
    //     * Left endpoint at 10 subtracts for c = 9, 10, 11.
    //     * Right endpoint at 14 subtracts for c = 13, 14, 15.
    // - Net result: shoulders at -1 on either side with the midpoint staying at +1.
    let expected = vec![
        wps_run("chr1", 9, 12, -1.0),
        wps_run("chr1", 12, 13, 1.0),
        wps_run("chr1", 13, 15, -1.0),
    ];

    let actual = run_wps(&cfg)?;

    assert_runs_equal(&actual, &expected);
    Ok(())
}

#[test]
fn fragment_equal_to_window_with_zero_runs_emits_shoulders() -> Result<()> {
    let fixture = make_fixture("wps_equal_window_zero_runs", &[(10, 14)])?;
    let out_dir = TempDir::new()?;
    let cfg = make_config(
        4,
        true,
        &fixture.bam,
        out_dir.path(),
        "equal_window_zero_runs",
    );

    // Fragment length equals window size:
    // - Full coverage contributes +1 at center 12.
    // - Endpoint windows carry -1 shoulders on both sides.
    // - Remaining centers stay at 0 and are kept because keep_zero_runs=true.
    let expected = vec![
        wps_run("chr1", 2, 9, 0.0),
        wps_run("chr1", 9, 12, -1.0),
        wps_run("chr1", 12, 13, 1.0),
        wps_run("chr1", 13, 15, -1.0),
        wps_run("chr1", 15, 30, 0.0),
    ];

    let actual = run_wps(&cfg)?;

    let clipped = clip_runs(&actual, 30);
    assert_runs_equal(&clipped, &expected);
    Ok(())
}

#[test]
fn empty_bam_emits_single_zero_run_per_chromosome() -> Result<()> {
    // Chromosomes long enough to admit two tiles each.
    let chrom_defs = vec![("chr1".to_string(), 400u32), ("chr2".to_string(), 400u32)];
    let tile_bp = 200u32;
    let fixture = bam_from_fragments("wps_empty", chrom_defs.clone(), Vec::new(), Vec::new())?;
    let out_dir = TempDir::new()?;

    let mut cfg = make_config(4, true, &fixture.bam, out_dir.path(), "empty_two_chr");
    cfg.shared_args.chromosomes.chromosomes = Some(vec!["chr1".to_string(), "chr2".to_string()]);
    cfg.set_tile_size(tile_bp);

    let runs = run_wps_with_chrom(&cfg)?;

    // Each chromosome spans two tiles; we intentionally expose the uncrossed tile boundaries to
    // keep merge_positional_tiles fast (simple stream copy).
    // Valid centers start 2 bp in and stop at 399, so every zero run begins at 2 and ends at 399 (exclusive)
    ensure!(
        runs.len() == 4,
        "expected exactly 4 runs (two per chromosome), got {runs:?}"
    );

    // Chromosomes are 400 bp and the 4 bp window means valid centers start at 2 and stop before 399 (exclusive)
    // Each 200 bp tile is emitted separately so every chromosome contributes two zero runs
    let expected = vec![
        ("chr1".to_string(), 2, 200, 0.0f32),
        ("chr1".to_string(), 200, 399, 0.0f32),
        ("chr2".to_string(), 2, 200, 0.0f32),
        ("chr2".to_string(), 200, 399, 0.0f32),
    ];

    for (run, exp) in runs.iter().zip(expected.into_iter()) {
        assert_eq!(
            (&run.chromosome, run.start, run.end, run.value),
            (&exp.0, exp.1, exp.2, exp.3),
            "unexpected run"
        );
    }

    Ok(())
}

#[test]
fn empty_bam_without_keep_zero_runs_outputs_nothing() -> Result<()> {
    let chrom_defs = vec![("chr1".to_string(), 400u32), ("chr2".to_string(), 400u32)];
    let fixture = bam_from_fragments(
        "wps_empty_nozeros",
        chrom_defs.clone(),
        Vec::new(),
        Vec::new(),
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = make_config(
        4,
        false,
        &fixture.bam,
        out_dir.path(),
        "empty_two_chr_nozeros",
    );
    cfg.shared_args.chromosomes.chromosomes = Some(vec!["chr1".to_string(), "chr2".to_string()]);
    cfg.set_tile_size(200);

    let runs = run_wps_with_chrom(&cfg)?;

    ensure!(
        runs.is_empty(),
        "expected no runs when keep_zero_runs=false, got {runs:?}"
    );

    Ok(())
}

#[test]
fn long_fragment_fixture_produces_expected_wps_runs() -> Result<()> {
    let fixture = long_inward_fragment_series_bam("wps_long_fragment_fixture")?;
    let out_dir = TempDir::new()?;
    let mut cfg = make_config(
        WPS_WINDOW_SIZE_BP,
        true,
        &fixture.bam,
        out_dir.path(),
        "long_fragment_wps",
    );
    // Allow the 1000bp inserts from the shared fixture.
    cfg.set_max_fragment_length(1_000);

    // Manual expectations (window size 120 -> left_span = right_span = 60):
    // - A fragment contributes +1 wherever the 120bp window fits entirely inside it,
    //   meaning centers [start+60, end-60). The first fragment therefore yields
    //   [60, 341) at +1.
    // - When two fragments are 400bp apart, their full-coverage bands overlap for
    //   81 bases, so the combined run reaches +2 (e.g., [460, 541)).
    // - Between fragments we see 0 plateaus once the left fragment ends and before
    //   the next one begins.
    // - We assert only the leading repeats because the pattern continues across
    //   the contig with the same spacings.
    let runs = run_wps(&cfg)?;
    let expected = vec![
        wps_run("chr1", 60, 341, 1.0),
        wps_run("chr1", 341, 460, 0.0),
        wps_run("chr1", 460, 541, 2.0),
        wps_run("chr1", 541, 659, 0.0),
        wps_run("chr1", 659, 741, 1.0),
        wps_run("chr1", 741, 860, 0.0),
        wps_run("chr1", 860, 941, 2.0),
        wps_run("chr1", 941, 1000, 0.0),
    ];

    assert!(
        runs.len() >= expected.len(),
        "expected at least {} runs, got {}",
        expected.len(),
        runs.len()
    );
    let prefix: Vec<WpsRun> = runs.iter().take(expected.len()).cloned().collect();
    assert_runs_equal(&prefix, &expected);
    Ok(())
}

#[test]
fn global_mode_handles_three_chromosomes() -> Result<()> {
    let fixture =
        make_three_chrom_fixture("wps_three_chr_global", &[(10, 22), (10, 22), (10, 22)])?;
    let out_dir = TempDir::new()?;
    let mut cfg = make_config(4, false, &fixture.bam, out_dir.path(), "three_chr_global");
    cfg.shared_args.chromosomes.chromosomes = Some(vec![
        "chr1".to_string(),
        "chr2".to_string(),
        "chr3".to_string(),
    ]);

    // Manual expectations per chromosome for one fragment [10, 22) and window size 4:
    // - [9, 12) at -1 from the left endpoint
    // - [12, 21) at +1 where the window fits fully inside the fragment
    // - [21, 23) at -1 from the right endpoint
    let expected = vec![
        wps_run("chr1", 9, 12, -1.0),
        wps_run("chr1", 12, 21, 1.0),
        wps_run("chr1", 21, 23, -1.0),
        wps_run("chr2", 9, 12, -1.0),
        wps_run("chr2", 12, 21, 1.0),
        wps_run("chr2", 21, 23, -1.0),
        wps_run("chr3", 9, 12, -1.0),
        wps_run("chr3", 12, 21, 1.0),
        wps_run("chr3", 21, 23, -1.0),
    ];

    let actual = run_wps_with_chrom(&cfg)?;
    assert_runs_equal(&actual, &expected);
    Ok(())
}

#[test]
fn by_size_total_handles_three_chromosomes() -> Result<()> {
    let fixture =
        make_three_chrom_fixture("wps_three_chr_by_size", &[(10, 22), (10, 22), (10, 22)])?;
    let out_dir = TempDir::new()?;
    let mut cfg = make_config(4, false, &fixture.bam, out_dir.path(), "three_chr_by_size");
    cfg.shared_args.chromosomes.chromosomes = Some(vec![
        "chr1".to_string(),
        "chr2".to_string(),
        "chr3".to_string(),
    ]);
    cfg.shared_args.windows = WindowsArgs {
        by_size: Some(100),
        by_bed: None,
    };
    cfg.per_window = Some(CoverageWindowAction::Total);

    // The per-chromosome WPS runs from the global test sum to:
    // - [9, 12)  -> -3
    // - [12, 21) -> +9
    // - [21, 23) -> -2
    // Total WPS over the chromosome window [0, 100) is therefore 4.
    run_wps_quiet(&cfg)?;

    let output_path = out_dir.path().join("three_chr_by_size.wps.total.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t100\t4\t0",
            "chr2\t0\t100\t4\t0",
            "chr3\t0\t100\t4\t0",
        ]
    );

    Ok(())
}

#[test]
fn by_size_total_non_aligned_tiles_reduce_crossing_bins_by_logical_start() -> Result<()> {
    let fixture = bam_from_fragments(
        "wps_non_aligned_by_size_reduce",
        vec![("chr1".to_string(), 300u32)],
        Vec::new(),
        Vec::new(),
    )?;
    let out_dir = TempDir::new()?;
    let mut cfg = make_config(
        4,
        false,
        &fixture.bam,
        out_dir.path(),
        "non_aligned_by_size",
    );
    cfg.shared_args.chromosomes.chromosomes = Some(vec!["chr1".to_string()]);
    cfg.set_tile_size(200);
    cfg.shared_args.windows = WindowsArgs {
        by_size: Some(150),
        by_bed: None,
    };
    cfg.per_window = Some(CoverageWindowAction::Total);

    // Manual expectations:
    // - tile_size=200 and by_size=150 do not align, so the reducer must combine cross-tile
    //   partials instead of concatenating aligned tile finals.
    // - The logical bins are [0,150) and [150,300). The second bin crosses the tile boundary at 200,
    //   so both tiles must contribute under the same logical start=150 key.
    // - Empty BAM means total WPS is 0 in both bins.
    // - With window size 4, invalid centers are 0, 1, and 299 because the WPS window would extend
    //   past chromosome bounds. That gives blacklisted_positions 2 for [0,150) and 1 for [150,300).
    run_wps_quiet(&cfg)?;

    let output_path = out_dir.path().join("non_aligned_by_size.wps.total.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t150\t0\t2",
            "chr1\t150\t300\t0\t1",
        ]
    );

    Ok(())
}

#[test]
fn by_bed_total_handles_three_chromosomes() -> Result<()> {
    let fixture =
        make_three_chrom_fixture("wps_three_chr_by_bed", &[(10, 22), (10, 22), (10, 22)])?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("three_chr_windows.bed");
    write_bed4(
        &bed_path,
        &[
            Bed4Row::new("chr1", 0, 100, "chr1_window"),
            Bed4Row::new("chr2", 0, 100, "chr2_window"),
            Bed4Row::new("chr3", 0, 100, "chr3_window"),
        ],
    )?;

    let mut cfg = make_config(4, false, &fixture.bam, out_dir.path(), "three_chr_by_bed");
    cfg.shared_args.chromosomes.chromosomes = Some(vec![
        "chr1".to_string(),
        "chr2".to_string(),
        "chr3".to_string(),
    ]);
    cfg.shared_args.windows = WindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
    };
    cfg.per_window = Some(CoverageWindowAction::Total);

    run_wps_quiet(&cfg)?;

    let output_path = out_dir.path().join("three_chr_by_bed.wps.total.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr1\t0\t100\t4\t0",
            "chr2\t0\t100\t4\t0",
            "chr3\t0\t100\t4\t0",
        ]
    );

    Ok(())
}

#[test]
fn by_bed_total_skips_chromosomes_without_windows_and_keeps_later_chromosomes() -> Result<()> {
    let fixture = bam_from_fragments(
        "wps_bed_skip_empty_chr",
        vec![("chr1".to_string(), 100), ("chr2".to_string(), 100)],
        vec![
            fragment_spec_on_tid(0, 10, 22),
            fragment_spec_on_tid(1, 10, 22),
        ],
        Vec::new(),
    )?;
    let out_dir = TempDir::new()?;
    let bed_path = out_dir.path().join("chr2_only_windows.bed");
    write_bed4(&bed_path, &[Bed4Row::new("chr2", 0, 100, "chr2_window")])?;

    let mut cfg = make_config(4, false, &fixture.bam, out_dir.path(), "chr2_only_by_bed");
    cfg.shared_args.chromosomes.chromosomes = Some(vec!["chr1".to_string(), "chr2".to_string()]);
    cfg.shared_args.windows = WindowsArgs {
        by_size: None,
        by_bed: Some(bed_path),
    };
    cfg.per_window = Some(CoverageWindowAction::Total);

    // Manual expectations:
    // - chr1 has a fragment but no BED windows, so BED mode should skip that chromosome entirely.
    // - chr2 has one fragment [10, 22) and one BED window [0, 100). The single-fragment
    //   derivation used elsewhere in this file gives a total WPS sum of 4 over that whole window.
    run_wps_quiet(&cfg)?;

    let output_path = out_dir.path().join("chr2_only_by_bed.wps.total.tsv.zst");
    let text = read_zst_to_string(&output_path)?;
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(
        lines,
        vec![
            "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions",
            "chr2\t0\t100\t4\t0",
        ]
    );

    Ok(())
}

fn assert_runs_equal(actual: &[WpsRun], expected: &[WpsRun]) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "expected {expected:?}, got {actual:?}"
    );
    for (idx, (act, exp)) in actual.iter().zip(expected.iter()).enumerate() {
        assert_eq!(
            act.chromosome, exp.chromosome,
            "run {idx} chromosome mismatch: expected {exp:?}, got {act:?}"
        );
        assert_eq!(
            act.start, exp.start,
            "run {idx} start mismatch: expected {exp:?}, got {act:?}"
        );
        assert_eq!(
            act.end, exp.end,
            "run {idx} end mismatch: expected {exp:?}, got {act:?}"
        );
        assert!(
            (act.value - exp.value).abs() < EPSILON,
            "run {idx} value mismatch: expected {exp:?}, got {act:?}"
        );
    }
}

fn clip_runs(runs: &[WpsRun], max_end: u32) -> Vec<WpsRun> {
    let mut out = Vec::new();
    for run in runs {
        if run.start >= max_end {
            break;
        }
        let end = run.end.min(max_end);
        out.push(WpsRun {
            chromosome: run.chromosome.clone(),
            start: run.start,
            end,
            value: run.value,
        });
        if run.end >= max_end {
            break;
        }
    }
    out
}
