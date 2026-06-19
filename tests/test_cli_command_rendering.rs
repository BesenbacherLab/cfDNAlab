#![cfg(feature = "cli")]

//! Smoke tests for the exported CLI rendering API.
//!
//! These tests exercise the public `run_like_cli` config types and the public
//! `ToCliCommand` trait from outside the crate. The scope is intentionally
//! narrower than the crate-local roundtrip tests in `src/cli_app_roundtrip_tests.rs`.
//! This file checks that each exported config can render a full `cfdna <command> ...`
//! argv, that the rendered command is accepted by the public Clap command tree,
//! and that representative default arguments are present in the output.
//!
//! These tests do not compare the parsed config back to the original config.
//! Exact config equality is covered by the crate-local roundtrip tests, which
//! can inspect the private command enum without widening the public API.

use std::{ffi::OsString, path::PathBuf};

use cfdnalab::{ToCliCommand, build_docs_command};

use cfdnalab::run_like_cli::common::{ChromosomeArgs, IOCArgs};

/// Assert that rendered argv has the expected command shape and parses through Clap.
///
/// This helper deliberately uses `build_docs_command()` because this test module
/// is outside the crate and should exercise the same public command tree that
/// downstream Rust callers can rely on. It returns the rendered strings so each
/// test can also assert representative argument/value pairs without repeating
/// lossy `OsString` conversion code.
fn assert_cli_accepts(args: Vec<OsString>, expected_subcommand: &str) -> Vec<String> {
    let rendered: Vec<String> = args
        .iter()
        .map(|value| value.to_string_lossy().into_owned())
        .collect();
    assert_eq!(rendered[0], "cfdna");
    assert_eq!(rendered[1], expected_subcommand);

    build_docs_command()
        .try_get_matches_from(args)
        .unwrap_or_else(|error| {
            panic!(
                "rendered CLI for {expected_subcommand} did not parse:\n{error}\nargv: {rendered:?}"
            )
        });
    rendered
}

/// Assert that an argv vector contains a flag immediately followed by its value.
///
/// The CLI renderer writes flags and values as separate arguments. Checking
/// adjacent pairs catches both missing defaults and wrong value formatting while
/// avoiding a brittle assertion on the full argument order.
fn assert_contains_pair(args: &[String], flag: &str, value: &str) {
    assert!(
        args.windows(2)
            .any(|window| window[0] == flag && window[1] == value),
        "expected rendered CLI to contain `{flag} {value}`, got {args:?}"
    );
}

fn ioc() -> IOCArgs {
    IOCArgs {
        bam: PathBuf::from("input.bam"),
        output_dir: PathBuf::from("out"),
        n_threads: 2,
    }
}

fn chromosomes() -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
        chromosomes_file: None,
    }
}

#[cfg(feature = "cmd_bam_to_bam")]
#[test]
fn bam_to_bam_config_renders_cli_call() {
    use cfdnalab::run_like_cli::bam_to_bam::BamToBamConfig;

    let mut config = BamToBamConfig::new(
        PathBuf::from("input.bam"),
        PathBuf::from("output.bam"),
        chromosomes(),
    );
    config.set_by_bed(Some(PathBuf::from("windows.bed")));

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "bam-to-bam");
    assert_contains_pair(&args, "--min-mapq", "0");
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn bam_to_frag_config_renders_cli_call() {
    use cfdnalab::run_like_cli::bam_to_frag::BamToFragConfig;

    let mut config = BamToFragConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample");
    config.set_by_bed(Some(PathBuf::from("windows.bed")));

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "bam-to-frag");
    assert_contains_pair(&args, "--min-mapq", "0");
}

#[cfg(feature = "cmd_frag_to_bam")]
#[test]
fn frag_to_bam_config_renders_cli_call() {
    use cfdnalab::run_like_cli::frag_to_bam::FragToBamConfig;

    let mut config = FragToBamConfig::new(
        PathBuf::from("sample.frag.tsv"),
        PathBuf::from("out"),
        chromosomes(),
        PathBuf::from("genome.chrom.sizes"),
    );
    config.set_output_prefix("sample");
    config.set_frag_header(Some(PathBuf::from("sample.frag.header.tsv")));

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "frag-to-bam");
    assert_contains_pair(&args, "--min-mapq", "0");
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn coverage_weights_config_renders_cli_call() {
    use cfdnalab::run_like_cli::coverage_weights::CoverageWeightsConfig;

    let mut config = CoverageWeightsConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample".to_string());
    config.set_ignore_gap(true);

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "coverage-weights");
    assert_contains_pair(&args, "--stride", "500000");
}

#[cfg(feature = "cmd_fragment_count_weights")]
#[test]
fn fragment_count_weights_config_renders_cli_call() {
    use cfdnalab::run_like_cli::fragment_count_weights::FragmentCountWeightsConfig;

    let mut config = FragmentCountWeightsConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample".to_string());

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "fragment-count-weights");
    assert_contains_pair(&args, "--stride", "500000");
}

#[cfg(feature = "cmd_fcoverage")]
#[test]
fn fcoverage_config_renders_cli_call() {
    use cfdnalab::run_like_cli::fcoverage::{CoverageWindowAction, FCoverageConfig};

    let mut config = FCoverageConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample");
    config.set_per_window(CoverageWindowAction::Total);
    config.set_windows(cfdnalab::run_like_cli::common::DistributionWindowsArgs {
        by_size: Some(1_000_000),
        by_bed: None,
        by_grouped_bed: None,
    });

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "fcoverage");
    assert_contains_pair(&args, "--decimals", "2");
}

#[cfg(feature = "cmd_gc_bias")]
#[test]
fn gc_bias_config_renders_cli_call() {
    use cfdnalab::run_like_cli::gc_bias::GCConfig;

    let mut config = GCConfig::new(
        ioc(),
        PathBuf::from("ref.2bit"),
        PathBuf::from("ref_gc.zarr"),
        chromosomes(),
    );
    config.set_output_prefix("sample".to_string());

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "gc-bias");
    assert_contains_pair(&args, "--outlier-method", "iqr");
}

#[cfg(feature = "cmd_gc_bias")]
#[test]
fn ref_gc_bias_config_renders_cli_call() {
    use cfdnalab::run_like_cli::common::{FragmentLengthArgs, Ref2BitRequiredArgs};
    use cfdnalab::run_like_cli::ref_gc_bias::{RefGCBiasConfig, RefGCWindowsArgs};

    let config = RefGCBiasConfig {
        ref_genome: Ref2BitRequiredArgs {
            ref_2bit: PathBuf::from("ref.2bit"),
        },
        output_dir: PathBuf::from("out"),
        output_prefix: "hg38".to_string(),
        n_threads: 2,
        n_positions: 10_000,
        seed: Some(7),
        windows: RefGCWindowsArgs {
            by_bed: Some(PathBuf::from("regions.bed")),
        },
        chromosomes: chromosomes(),
        blacklist: Some(vec![PathBuf::from("blacklist.bed")]),
        fragment_lengths: FragmentLengthArgs::default(),
        end_offset: 10,
        skip_interpolation: false,
        smoothing_sigma: 0.8,
        smoothing_radius: 2,
        skip_smoothing: false,
        tile_size: 10_000_000,
        logging: Default::default(),
    };

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "ref-gc-bias");
    assert_contains_pair(&args, "--end-offset", "10");
}

#[cfg(feature = "cmd_lengths")]
#[test]
fn lengths_config_renders_cli_call() {
    use cfdnalab::run_like_cli::lengths::LengthsConfig;

    let mut config = LengthsConfig::new(ioc(), chromosomes());
    config.output_prefix = "sample".to_string();
    config.set_windows(cfdnalab::run_like_cli::common::DistributionWindowsArgs {
        by_size: Some(1_000_000),
        by_bed: None,
        by_grouped_bed: None,
    });
    config.set_length_bins_spec("30:101:10");

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "lengths");
    assert_contains_pair(&args, "--decimals", "6");
}

#[cfg(feature = "cmd_midpoints")]
#[test]
fn midpoints_config_renders_cli_call() {
    use cfdnalab::run_like_cli::midpoints::MidpointsConfig;

    let mut config = MidpointsConfig::new(ioc(), chromosomes(), PathBuf::from("sites.bed"));
    config.set_output_prefix("sample");
    config.set_length_bins(vec![30, 80, 151]);

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "midpoints");
    assert_contains_pair(&args, "--bin-size", "1");
}

#[cfg(feature = "cmd_ends")]
#[test]
fn ends_config_renders_cli_call() {
    use cfdnalab::run_like_cli::ends::EndsConfig;

    let mut config = EndsConfig::new(ioc(), chromosomes(), 2, 2);
    config.set_min_mapq(20);
    config.set_tile_size(1_000_000);

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "ends");
    assert_contains_pair(&args, "--source-inside", "read");
}

#[cfg(feature = "cmd_fragment_kmers")]
#[test]
fn fragment_kmers_config_renders_cli_call() {
    use cfdnalab::run_like_cli::common::Ref2BitRequiredArgs;
    use cfdnalab::run_like_cli::fragment_kmers::FragmentKmersConfig;

    let mut config = FragmentKmersConfig::new(
        ioc(),
        Ref2BitRequiredArgs {
            ref_2bit: PathBuf::from("ref.2bit"),
        },
        chromosomes(),
    );
    config.set_output_prefix("sample".to_string());
    config.set_kmer_sizes(vec![3, 5]);

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "fragment-kmers");
    assert_contains_pair(&args, "--frame", "left");
}

#[cfg(feature = "cmd_transitions")]
#[test]
fn transitions_config_renders_cli_call() {
    use cfdnalab::run_like_cli::common::Ref2BitRequiredArgs;
    use cfdnalab::run_like_cli::transitions::TransitionsConfig;

    let mut config = TransitionsConfig::new(
        ioc(),
        Ref2BitRequiredArgs {
            ref_2bit: PathBuf::from("ref.2bit"),
        },
        chromosomes(),
    );
    config.set_output_prefix("sample".to_string());
    config.set_orders(vec![1, 3]);

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "transitions");
    assert_contains_pair(&args, "--indel-mode", "ignore");
}

#[cfg(feature = "cmd_wps")]
#[test]
fn wps_config_renders_cli_call() {
    use cfdnalab::run_like_cli::wps::WPSConfig;

    let mut config = WPSConfig::new(ioc(), chromosomes(), None);
    config.set_output_prefix("sample".to_string());
    config.set_window_size(121);

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "wps");
    assert_contains_pair(&args, "--decimals", "2");
}

#[cfg(feature = "cmd_wps_peaks")]
#[test]
fn wps_peaks_config_renders_cli_call() {
    use cfdnalab::run_like_cli::wps_peaks::{PeaksWindowAction, WPSPeaksConfig};

    let mut config = WPSPeaksConfig::new(ioc(), chromosomes(), Some(PeaksWindowAction::Stats));
    config.set_window_size(121);
    config.set_min_peak_height(6.5);

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "wps-peaks");
    assert_contains_pair(&args, "--normalize-bp", "1000");
}

#[cfg(feature = "cmd_prepare_windows")]
#[test]
fn prepare_windows_config_renders_cli_call() {
    use cfdnalab::run_like_cli::prepare_windows::PrepareConfig;

    let config = PrepareConfig {
        input: PathBuf::from("windows.tsv"),
        output: PathBuf::from("prepared.bed"),
        group_cols: vec!["3".to_string()],
        ..PrepareConfig::default()
    };

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "prep-windows");
    assert_contains_pair(&args, "--header", "auto");
}

#[cfg(feature = "cmd_visualize_positions")]
#[test]
fn visualize_positions_config_renders_cli_call() {
    use cfdnalab::run_like_cli::common::{BaseSelectionArgs, FragmentPositionSelectionArgs};
    use cfdnalab::run_like_cli::visualize_positions::{Style, VisualizePositionsConfig};
    use cfdnalab::{positioning::BasesFrom, positioning::MismatchBasesFrom};

    let config = VisualizePositionsConfig {
        position_selection: FragmentPositionSelectionArgs::default(),
        base_selection: BaseSelectionArgs {
            bases_from: BasesFrom::Reference,
            mismatch_bases_from: MismatchBasesFrom::NearestRead,
        },
        work_dir: PathBuf::from("work"),
        lengths: Some(vec![100, 120]),
        length_range: None,
        kmer_sizes: Some(vec![3]),
        style: Style::Svg,
        width: Some(80),
        height: Some(240),
        output: Some(PathBuf::from("figure.svg")),
        label: Some("demo".to_string()),
        hide_index: true,
        show_half: true,
        hide_mid: false,
    };

    let args = assert_cli_accepts(config.to_cli_args().unwrap(), "visualize-positions");
    assert_contains_pair(&args, "--bases-from", "reference");
}
