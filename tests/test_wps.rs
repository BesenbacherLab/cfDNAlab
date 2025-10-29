mod fixtures;

use anyhow::{Context, Result};
use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs};
use cfdnalab::commands::fcoverage::window_results::CoverageWindowAction;
use cfdnalab::commands::wps::config::WPSConfig;
use cfdnalab::commands::wps::wps::run as run_fn;
use fixtures::{BamFixture, FragmentSpec, ReadSpec, bam_from_specs};
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

#[derive(Debug, Clone, PartialEq)]
struct WpsRun {
    start: u32,
    end: u32,
    value: f32,
}

fn wps_run(start: u32, end: u32, value: f32) -> WpsRun {
    WpsRun { start, end, value }
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
        cigar: vec![('M', read_len)],
        seq: vec![b'A'; read_len as usize],
        qual: 40,
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
        cigar: vec![('M', read_len)],
        seq: vec![b'T'; read_len as usize],
        qual: 40,
        is_reverse: true,
        mapq: 60,
        flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
        mate_tid: Some(0),
        mate_pos: Some(forward_pos),
        insert_size: -fragment_span,
    };

    FragmentSpec { forward, reverse }
}

fn make_fixture(name: &str, fragments: &[(u32, u32)]) -> Result<BamFixture> {
    let chrom_len = fragments
        .iter()
        .map(|(_, end)| end + 100)
        .max()
        .unwrap_or(500);
    let specs: Vec<FragmentSpec> = fragments
        .iter()
        .map(|(start, end)| fragment_spec(*start, *end))
        .collect();
    bam_from_specs(
        vec![("chr1".to_string(), chrom_len)],
        specs,
        Vec::new(),
        name,
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
    cfg.set_output_prefix(prefix);
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
    run_fn(cfg)?;
    let prefix = cfg.output_prefix.trim();
    let bedgraph_path = cfg
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
        let _chromosome = cols
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

        runs.push(WpsRun { start, end, value });
    }

    Ok(runs)
}

fn assert_runs_equal(actual: &[WpsRun], expected: &[WpsRun]) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "expected {expected:?}, got {actual:?}"
    );
    for (idx, (act, exp)) in actual.iter().zip(expected.iter()).enumerate() {
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

#[ignore = "Awaiting implementation of wps::run_inner"]
#[test]
fn single_fragment_produces_central_plateau() -> Result<()> {
    let fixture = make_fixture("wps_single_fragment", &[(10, 22)])?;
    let out_dir = TempDir::new()?;
    let cfg = make_config(4, false, &fixture.bam, out_dir.path(), "single_fragment");

    // Manual expectations:
    // - Window size 4 gives left_span = right_span = 2, so each WPS sample is centred at c and
    //   covers [c - 2, c + 2). A fragment satisfies the geometric test for full coverage when
    //   c is between 12 and 20 inclusive, but centres 12 and 20 still see an endpoint, leaving
    //   only 13..=19 with a lasting +1.
    // - Each fragment end subtracts 1 wherever the window still contains that endpoint.
    //   * Left end at 10 affects centres 9..=12 (windows [7,11), [8,12), [9,13), [10,14)), so
    //     centre 12 receives both +1 (full span) and -1 (left end) -> net 0.
    //   * Right end at 21 affects centres 20..=23 (windows [18,22), [19,23), [20,24), [21,25)),
    //     and centre 20 again cancels to 0.
    // - Putting it together: 9-12 and 20-23 carry -1 (only endpoint contribution), 13-19 carries
    //   +1 (full coverage without endpoints), and all other centres remain 0. In BED form the
    //   positive plateau therefore emits [13,20) (centres 13..=19) and the right-side dip starts
    //   at [20,24), because centre 20 cancels to 0. Centre 8 does NOT include base 10 because its
    //   window is [6,10) -- the upper bound is exclusive -- so the left-end penalty begins exactly
    //   at 9.
    let expected = vec![
        wps_run(9, 13, -1.0),
        wps_run(13, 20, 1.0),
        wps_run(20, 24, -1.0),
    ];

    let actual = run_wps(&cfg)?;

    assert_runs_equal(&actual, &expected);
    Ok(())
}

#[ignore = "Awaiting implementation of wps::run_inner"]
#[test]
fn overlapping_fragments_stack_scores() -> Result<()> {
    let fixture = make_fixture("wps_overlapping", &[(0, 20), (4, 12)])?;
    let out_dir = TempDir::new()?;
    let cfg = make_config(4, false, &fixture.bam, out_dir.path(), "overlapping");

    // Manual expectations for two fragments:
    // - Fragment F1: [0, 20), fragment F2: [4, 12); window size 4 keeps left_span = right_span = 2.
    // - Strict full-span rule means:
    //     * F1 contributes +1 only for centres 3..=17 (window strictly inside [0,20)).
    //     * F2 contributes +1 only for centres 7..=9 (window strictly inside [4,12)).
    // - Endpoint penalties:
    //     * F1 left end at 0 hits centres -1..=2; right end at 19 hits 18..=21.
    //     * F2 left end at 4 hits centres 3..=6; right end at 11 hits 10..=13.
    // - Net behaviour across 0..24:
    //     * 0-3: only left-end penalties -> -1.
    //     * 7-10: both fragments span -> +2.
    //     * 14-18: only F1 spans -> +1.
    //     * 18-22: only right-end penalties remain -> -1.
    let expected = vec![
        wps_run(0, 3, -1.0),
        wps_run(7, 10, 2.0),
        wps_run(14, 18, 1.0),
        wps_run(18, 22, -1.0),
    ];

    let actual = run_wps(&cfg)?;

    assert_runs_equal(&actual, &expected);
    Ok(())
}

#[ignore = "Awaiting implementation of wps::run_inner"]
#[test]
fn keep_zero_runs_emits_flat_segments() -> Result<()> {
    let fixture = make_fixture("wps_keep_zero", &[(10, 22)])?;
    let out_dir = TempDir::new()?;
    let cfg = make_config(4, true, &fixture.bam, out_dir.path(), "keep_zero");

    // Same geometry as the single-fragment test, but keep_zero_runs=true means we retain zero
    // plateaus between non-zero segments:
    // - Leading zeros before the first penalty (0-9).
    // - Trailing zeros (24-30).
    // Non-zero spans follow the strict-span expectations from the first test.
    let expected = vec![
        wps_run(0, 9, 0.0),
        wps_run(9, 13, -1.0),
        wps_run(13, 20, 1.0),
        wps_run(20, 24, -1.0),
        wps_run(24, 30, 0.0),
    ];

    let actual = run_wps(&cfg)?;

    assert_runs_equal(&actual, &expected);
    Ok(())
}

#[ignore = "Awaiting implementation of wps::run_inner"]
#[test]
fn fragment_equal_to_window_removes_central_signal() -> Result<()> {
    let fixture = make_fixture("wps_equal_window", &[(10, 14)])?;
    let out_dir = TempDir::new()?;
    let cfg = make_config(4, false, &fixture.bam, out_dir.path(), "equal_window");

    // Fragment length exactly matches the window (4 bp):
    // - No centre achieves full coverage, because the window span cannot sit strictly inside the
    //   fragment once both edges must remain within [10,14).
    // - Both fragment ends still subtract wherever the window contains them. With left_span/right_span=2,
    //   the affected centres are 9..=11 for the left end at 10 and 13..=15 for the right end at 13.
    // - The overlapping penalties merge into the single continuous dip 9-16 at -1, with zeros elsewhere.
    let expected = vec![wps_run(9, 16, -1.0)];

    let actual = run_wps(&cfg)?;

    assert_runs_equal(&actual, &expected);
    Ok(())
}
