#![cfg(not(feature = "cli"))]

//! No-CLI smoke tests for the exported command rendering API.
//!
//! Downstream Rust packages should be able to construct exported `run_like_cli`
//! configs and render equivalent command strings without enabling the `cli`
//! feature or pulling in Clap. These tests compile and run only when `cli` is
//! disabled, which catches accidental dependencies from `ToCliCommand` or the
//! config renderers back to CLI parsing code.
//!
//! True parse roundtrips are not possible in this feature set because the Clap
//! command tree is intentionally absent. The crate-local roundtrip tests cover
//! parse equality when `cli` is enabled. This file checks the no-CLI equivalent:
//! each command renders a full `cfdna <command> ...` argv and representative
//! all-argument/default values are present.

use std::{ffi::OsString, path::PathBuf};

use cfdnalab::{
    ToCliCommand,
    run_like_cli::common::{ChromosomeArgs, IOCArgs},
};

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

fn rendered_strings(args: &[OsString]) -> Vec<String> {
    args.iter()
        .map(|value| value.to_string_lossy().into_owned())
        .collect()
}

/// Render a config and assert the public no-CLI command shape.
///
/// This helper avoids Clap entirely. It checks the exported `ToCliCommand`
/// implementation by inspecting the returned argv and the shell-quoted display
/// string that downstream packages can use in logs or provenance records.
fn render_no_cli_command(
    config: &impl ToCliCommand,
    expected_subcommand: &str,
) -> (Vec<String>, String) {
    let args = config
        .to_cli_args()
        .expect("config should render without the cli feature");
    let rendered = rendered_strings(&args);

    assert!(
        rendered.len() >= 2,
        "rendered CLI should include program and subcommand, got {rendered:?}"
    );
    assert_eq!(rendered[0], "cfdna");
    assert_eq!(rendered[1], expected_subcommand);

    let command = config
        .to_cli_string()
        .expect("config should render to a display command without the cli feature");
    assert!(
        command.starts_with(&format!("cfdna {expected_subcommand}")),
        "display command should start with `cfdna {expected_subcommand}`, got {command:?}"
    );

    (rendered, command)
}

/// Assert that an argv vector contains a flag immediately followed by its value.
///
/// The renderer returns tokenized argv rather than one raw shell string, so
/// adjacent-pair checks verify values without relying on the total argument
/// ordering.
fn assert_contains_pair(args: &[String], flag: &str, value: &str) {
    assert!(
        args.windows(2)
            .any(|window| window[0] == flag && window[1] == value),
        "expected rendered CLI to contain `{flag} {value}`, got {args:?}"
    );
}

#[cfg(feature = "cmd_bam_to_bam")]
#[test]
fn bam_to_bam_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::bam_to_bam::BamToBamConfig;

    let mut config = BamToBamConfig::new(
        PathBuf::from("input.bam"),
        PathBuf::from("output.bam"),
        chromosomes(),
    );
    config.set_by_bed(Some(PathBuf::from("windows.bed")));

    let (args, _command) = render_no_cli_command(&config, "bam-to-bam");
    assert_contains_pair(&args, "--min-mapq", "0");
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn bam_to_frag_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::bam_to_frag::BamToFragConfig;

    let mut config = BamToFragConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample");
    config.set_by_bed(Some(PathBuf::from("windows.bed")));

    let (args, _command) = render_no_cli_command(&config, "bam-to-frag");
    assert_contains_pair(&args, "--min-mapq", "0");
}

#[cfg(feature = "cmd_frag_to_bam")]
#[test]
fn frag_to_bam_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::frag_to_bam::FragToBamConfig;

    let mut config = FragToBamConfig::new(
        PathBuf::from("sample.frag.tsv"),
        PathBuf::from("out"),
        chromosomes(),
        PathBuf::from("genome.chrom.sizes"),
    );
    config.set_output_prefix("sample");
    config.set_frag_header(Some(PathBuf::from("sample.frag.header.tsv")));

    let (args, _command) = render_no_cli_command(&config, "frag-to-bam");
    assert_contains_pair(&args, "--min-mapq", "0");
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn coverage_weights_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::coverage_weights::CoverageWeightsConfig;

    let mut config = CoverageWeightsConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample".to_string());
    config.set_ignore_gap(true);

    let (args, _command) = render_no_cli_command(&config, "coverage-weights");
    assert_contains_pair(&args, "--stride", "500000");
}

#[cfg(feature = "cmd_fragment_count_weights")]
#[test]
fn fragment_count_weights_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::fragment_count_weights::FragmentCountWeightsConfig;

    let mut config = FragmentCountWeightsConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample".to_string());

    let (args, _command) = render_no_cli_command(&config, "fragment-count-weights");
    assert_contains_pair(&args, "--stride", "500000");
}

#[cfg(feature = "cmd_fcoverage")]
#[test]
fn fcoverage_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::{
        common::DistributionWindowsArgs,
        fcoverage::{CoverageWindowAction, FCoverageConfig},
    };

    let mut config = FCoverageConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample");
    config.set_per_window(CoverageWindowAction::Total);
    config.set_windows(DistributionWindowsArgs {
        by_size: Some(1_000_000),
        by_bed: None,
        by_grouped_bed: None,
    });

    let (args, _command) = render_no_cli_command(&config, "fcoverage");
    assert_contains_pair(&args, "--decimals", "2");
}

#[cfg(feature = "cmd_gc_bias")]
#[test]
fn gc_bias_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::gc_bias::GCConfig;

    let mut config = GCConfig::new(
        ioc(),
        PathBuf::from("ref.2bit"),
        PathBuf::from("ref_gc.zarr"),
        chromosomes(),
    );
    config.set_output_prefix("sample".to_string());

    let (args, _command) = render_no_cli_command(&config, "gc-bias");
    assert_contains_pair(&args, "--outlier-method", "iqr");
}

#[cfg(feature = "cmd_ref_gc_bias")]
#[test]
fn ref_gc_bias_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::{
        common::{FragmentLengthArgs, Ref2BitRequiredArgs},
        ref_gc_bias::{RefGCBiasConfig, RefGCWindowsArgs},
    };

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

    let (args, _command) = render_no_cli_command(&config, "ref-gc-bias");
    assert_contains_pair(&args, "--end-offset", "10");
}

#[cfg(feature = "cmd_lengths")]
#[test]
fn lengths_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::{common::DistributionWindowsArgs, lengths::LengthsConfig};

    let mut config = LengthsConfig::new(ioc(), chromosomes());
    config.output_prefix = "sample".to_string();
    config.set_windows(DistributionWindowsArgs {
        by_size: Some(1_000_000),
        by_bed: None,
        by_grouped_bed: None,
    });
    config.set_length_bins_spec("30:101:10");

    let (args, command) = render_no_cli_command(&config, "lengths");
    assert_contains_pair(&args, "--decimals", "6");
    assert_contains_pair(&args, "--length-bins", "30:101:10");
    assert!(
        command.contains("30:101:10"),
        "display command should preserve length-bin range specs, got {command:?}"
    );
}

#[cfg(feature = "cmd_midpoints")]
#[test]
fn midpoints_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::midpoints::MidpointsConfig;

    let mut config = MidpointsConfig::new(ioc(), chromosomes(), PathBuf::from("sites.bed"));
    config.set_output_prefix("sample");
    config.set_length_bins_spec("30:151:10");

    let (args, command) = render_no_cli_command(&config, "midpoints");
    assert_contains_pair(&args, "--bin-size", "1");
    assert_contains_pair(&args, "--length-bins", "30:151:10");
    assert!(
        command.contains("30:151:10"),
        "display command should preserve length-bin range specs, got {command:?}"
    );
}

#[cfg(feature = "cmd_ends")]
#[test]
fn ends_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::ends::EndsConfig;

    let mut config = EndsConfig::new(ioc(), chromosomes(), 2, 2);
    config.set_min_mapq(20);
    config.set_tile_size(1_000_000);

    let (args, _command) = render_no_cli_command(&config, "ends");
    assert_contains_pair(&args, "--source-inside", "read");
}

#[cfg(feature = "cmd_fragment_kmers")]
#[test]
fn fragment_kmers_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::{
        common::Ref2BitRequiredArgs, fragment_kmers::FragmentKmersConfig,
    };

    let mut config = FragmentKmersConfig::new(
        ioc(),
        Ref2BitRequiredArgs {
            ref_2bit: PathBuf::from("ref.2bit"),
        },
        chromosomes(),
    );
    config.set_output_prefix("sample".to_string());
    config.set_kmer_sizes(vec![3, 5]);

    let (args, _command) = render_no_cli_command(&config, "fragment-kmers");
    assert_contains_pair(&args, "--frame", "left");
}

#[cfg(feature = "cmd_transitions")]
#[test]
fn transitions_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::{common::Ref2BitRequiredArgs, transitions::TransitionsConfig};

    let mut config = TransitionsConfig::new(
        ioc(),
        Ref2BitRequiredArgs {
            ref_2bit: PathBuf::from("ref.2bit"),
        },
        chromosomes(),
    );
    config.set_output_prefix("sample".to_string());
    config.set_orders(vec![1, 3]);

    let (args, _command) = render_no_cli_command(&config, "transitions");
    assert_contains_pair(&args, "--indel-mode", "ignore");
}

#[cfg(feature = "cmd_wps")]
#[test]
fn wps_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::wps::WPSConfig;

    let mut config = WPSConfig::new(ioc(), chromosomes(), None);
    config.set_output_prefix("sample".to_string());
    config.set_window_size(121);

    let (args, _command) = render_no_cli_command(&config, "wps");
    assert_contains_pair(&args, "--decimals", "2");
}

#[cfg(feature = "cmd_wps_peaks")]
#[test]
fn wps_peaks_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::wps_peaks::{PeaksWindowAction, WPSPeaksConfig};

    let mut config = WPSPeaksConfig::new(ioc(), chromosomes(), Some(PeaksWindowAction::Stats));
    config.set_window_size(121);
    config.set_min_peak_height(6.5);

    let (args, _command) = render_no_cli_command(&config, "wps-peaks");
    assert_contains_pair(&args, "--normalize-bp", "1000");
}

#[cfg(feature = "cmd_prepare_windows")]
#[test]
fn prepare_windows_renders_without_cli_feature() {
    use cfdnalab::run_like_cli::prepare_windows::PrepareConfig;

    let config = PrepareConfig {
        input: PathBuf::from("windows.tsv"),
        output: PathBuf::from("prepared.bed"),
        group_cols: vec!["3".to_string()],
        ..PrepareConfig::default()
    };

    let (args, _command) = render_no_cli_command(&config, "prep-windows");
    assert_contains_pair(&args, "--header", "auto");
}

#[cfg(feature = "cmd_visualize_positions")]
#[test]
fn visualize_positions_renders_without_cli_feature() {
    use cfdnalab::{
        positioning::{BasesFrom, MismatchBasesFrom},
        run_like_cli::{
            common::{BaseSelectionArgs, FragmentPositionSelectionArgs},
            visualize_positions::{Style, VisualizePositionsConfig},
        },
    };

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

    let (args, _command) = render_no_cli_command(&config, "visualize-positions");
    assert_contains_pair(&args, "--bases-from", "reference");
}
