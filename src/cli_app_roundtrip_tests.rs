//! Tests that rendered CLI commands are exact config roundtrips.
//!
//! These tests protect the `ToCliCommand` implementations from drifting away
//! from Clap parsing. Each fixture starts with a concrete config, renders the
//! full `cfdna <command> ...` argv, parses that argv through the real CLI
//! command, and requires the parsed config to equal the original config.
//! This is stronger than checking that the rendered argv merely parses.

use std::{ffi::OsString, fmt::Debug, path::PathBuf};

use clap::FromArgMatches;

use crate::{
    ToCliCommand,
    commands::cli_common::{ChromosomeArgs, IOCArgs},
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

/// Assert that a config renders to a CLI call that Clap parses back to the same config.
///
/// The macro is intentionally thin. `config_roundtrip_result()` does the real work and
/// returns structured failure messages that can be tested directly. This wrapper keeps the
/// per-command roundtrip tests concise while still failing as a normal Rust test assertion.
macro_rules! assert_config_roundtrips {
    ($config:expr, $variant:ident, $subcommand:literal) => {{
        config_roundtrip_result($config, $subcommand, |cmd| match cmd {
            super::Cmd::$variant(parsed_config) => Some(parsed_config),
            _ => None,
        })
        .unwrap_or_else(|message| panic!("{message}"));
    }};
}

fn config_roundtrip_result<T>(
    config: T,
    expected_subcommand: &str,
    extract_config: impl FnOnce(super::Cmd) -> Option<T>,
) -> Result<(), String>
where
    T: ToCliCommand + PartialEq + Debug,
{
    // Render the config through the public CLI-command renderer
    let args = config
        .to_cli_args()
        .map_err(|error| format!("config did not render to CLI arguments: {error}"))?;

    // Parse the rendered argv through the same Clap command tree used by the binary
    let parsed_command = parse_rendered_cli_args(args, expected_subcommand)?;

    // Extract the expected command variant from the parsed CLI enum
    let parsed_config = extract_config(parsed_command).ok_or_else(|| {
        format!("rendered CLI for {expected_subcommand} parsed as a different command")
    })?;

    // Require the parsed config to match the original exactly
    if parsed_config == config {
        Ok(())
    } else {
        Err(format!(
            "rendered CLI for {expected_subcommand} parsed to a different config\nexpected: {config:#?}\nparsed: {parsed_config:#?}"
        ))
    }
}

fn parse_rendered_cli_args(
    args: Vec<OsString>,
    expected_subcommand: &str,
) -> Result<super::Cmd, String> {
    // Keep a string copy for readable failure messages and basic command-shape checks
    let rendered = rendered_strings(&args);
    let program = rendered
        .first()
        .ok_or_else(|| format!("rendered CLI for {expected_subcommand} was empty"))?;
    if program != "cfdna" {
        return Err(format!(
            "rendered CLI for {expected_subcommand} must start with `cfdna`, got `{program}`\nargv: {rendered:?}"
        ));
    }

    // Check the subcommand before Clap parsing so command-routing failures are explicit
    let rendered_subcommand = rendered.get(1).ok_or_else(|| {
        format!("rendered CLI for {expected_subcommand} did not include a subcommand\nargv: {rendered:?}")
    })?;
    if rendered_subcommand != expected_subcommand {
        return Err(format!(
            "rendered CLI for {expected_subcommand} used subcommand `{rendered_subcommand}`\nargv: {rendered:?}"
        ));
    }

    // Parse through the real terminal command definition
    let matches = super::build_terminal_command()
        .try_get_matches_from(args)
        .map_err(|error| {
            format!(
                "rendered CLI for {expected_subcommand} did not parse:\n{error}\nargv: {rendered:?}"
            )
        })?;

    // Build the internal CLI enum from the parsed Clap matches
    let parsed_cli = super::Cli::from_arg_matches(&matches).map_err(|error| {
        format!("parsed matches for {expected_subcommand} did not build Cli: {error}")
    })?;

    Ok(parsed_cli.cmd)
}

#[cfg(feature = "cmd_bam_to_frag")]
fn bam_to_frag_roundtrip_config() -> crate::commands::bam_to_frag::config::BamToFragConfig {
    let mut config =
        crate::commands::bam_to_frag::config::BamToFragConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample");
    config.set_by_bed(Some(PathBuf::from("windows.bed")));
    config
}

#[cfg(feature = "cmd_bam_to_frag")]
fn bam_to_frag_roundtrip_args() -> Vec<OsString> {
    bam_to_frag_roundtrip_config()
        .to_cli_args()
        .expect("fixture config should render to CLI arguments")
}

fn expect_roundtrip_error<T>(result: Result<T, String>, context: &str) -> String {
    match result {
        Ok(_) => panic!("{context}"),
        Err(error) => error,
    }
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn roundtrip_helper_reports_wrong_program_name() {
    // Arrange: Simulate a renderer that does not start with the binary name.
    let mut args = bam_to_frag_roundtrip_args();
    args[0] = OsString::from("not-cfdna");

    // Act: Ask the parse helper to validate the rendered command shape.
    let error = expect_roundtrip_error(
        parse_rendered_cli_args(args, "bam-to-frag"),
        "wrong binary name should be reported as a roundtrip failure",
    );

    // Assert: The failure names both the invariant and the offending value.
    assert!(error.contains("must start with `cfdna`"));
    assert!(error.contains("not-cfdna"));
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn roundtrip_helper_reports_wrong_subcommand_name() {
    // Arrange: Simulate a renderer that emits a valid binary name but wrong subcommand.
    let mut args = bam_to_frag_roundtrip_args();
    args[1] = OsString::from("lengths");

    // Act: Validate the rendered command before Clap gets to route it.
    let error = expect_roundtrip_error(
        parse_rendered_cli_args(args, "bam-to-frag"),
        "wrong subcommand should be reported as a roundtrip failure",
    );

    // Assert: The failure points at command routing, not at config equality.
    assert!(error.contains("used subcommand `lengths`"));
    assert!(error.contains("bam-to-frag"));
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn roundtrip_helper_reports_clap_parse_errors() {
    // Arrange: Simulate a renderer that emits an unknown CLI flag.
    let mut args = bam_to_frag_roundtrip_args();
    args.push(OsString::from("--definitely-not-real"));

    // Act: Parse through the real Clap command tree.
    let error = expect_roundtrip_error(
        parse_rendered_cli_args(args, "bam-to-frag"),
        "unknown argument should be reported as a parse failure",
    );

    // Assert: The failure keeps the Clap error and argv context.
    assert!(error.contains("did not parse"));
    assert!(error.contains("--definitely-not-real"));
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn roundtrip_helper_reports_command_variant_mismatches() {
    // Arrange: Use a valid rendered command but reject the parsed variant.
    let config = bam_to_frag_roundtrip_config();

    // Act: The extractor returns None, matching the macro path for a wrong variant.
    let error = expect_roundtrip_error(
        config_roundtrip_result(config, "bam-to-frag", |_| None),
        "variant mismatch should be reported as a roundtrip failure",
    );

    // Assert: The failure distinguishes command routing from config-value mismatch.
    assert!(error.contains("parsed as a different command"));
    assert!(error.contains("bam-to-frag"));
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn roundtrip_helper_reports_config_value_mismatches() {
    // Arrange: Parse a valid command, then substitute a different config in the extractor.
    let original_config = bam_to_frag_roundtrip_config();
    let mut different_config = original_config.clone();
    different_config.set_output_prefix("different");

    // Act: Compare the intentionally different parsed config with the original config.
    let error = expect_roundtrip_error(
        config_roundtrip_result(original_config, "bam-to-frag", |_| Some(different_config)),
        "config mismatch should be reported as a roundtrip failure",
    );

    // Assert: The failure includes the mismatch class and enough detail to diagnose it.
    assert!(error.contains("parsed to a different config"));
    assert!(error.contains("different"));
}

#[cfg(feature = "cmd_bam_to_bam")]
#[test]
fn bam_to_bam_config_roundtrips_through_rendered_cli() {
    let mut config = crate::commands::bam_to_bam::config::BamToBamConfig::new(
        PathBuf::from("input.bam"),
        PathBuf::from("output.bam"),
        chromosomes(),
    );
    config.set_by_bed(Some(PathBuf::from("windows.bed")));

    assert_config_roundtrips!(config, BamToBam, "bam-to-bam");
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn bam_to_frag_config_roundtrips_through_rendered_cli() {
    assert_config_roundtrips!(bam_to_frag_roundtrip_config(), BamToFrag, "bam-to-frag");
}

#[cfg(feature = "cmd_frag_to_bam")]
#[test]
fn frag_to_bam_config_roundtrips_through_rendered_cli() {
    let mut config = crate::commands::frag_to_bam::config::FragToBamConfig::new(
        PathBuf::from("sample.frag.tsv"),
        PathBuf::from("out"),
        chromosomes(),
        PathBuf::from("genome.chrom.sizes"),
    );
    config.set_output_prefix("sample");
    config.set_frag_header(Some(PathBuf::from("sample.frag.header.tsv")));

    assert_config_roundtrips!(config, FragToBam, "frag-to-bam");
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn coverage_weights_config_roundtrips_through_rendered_cli() {
    let mut config =
        crate::commands::coverage_weights::config::CoverageWeightsConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample".to_string());
    config.set_ignore_gap(true);

    assert_config_roundtrips!(config, CoverageWeights, "coverage-weights");
}

#[cfg(feature = "cmd_fragment_count_weights")]
#[test]
fn fragment_count_weights_config_roundtrips_through_rendered_cli() {
    let mut config =
        crate::commands::fragment_count_weights::config::FragmentCountWeightsConfig::new(
            ioc(),
            chromosomes(),
        );
    config.set_output_prefix("sample".to_string());

    assert_config_roundtrips!(config, FragmentCountWeights, "fragment-count-weights");
}

#[cfg(feature = "cmd_fcoverage")]
#[test]
fn fcoverage_config_roundtrips_through_rendered_cli() {
    use crate::{
        commands::cli_common::DistributionWindowsArgs,
        commands::fcoverage::{config::FCoverageConfig, window_results::CoverageWindowAction},
    };

    let mut config = FCoverageConfig::new(ioc(), chromosomes());
    config.set_output_prefix("sample");
    config.set_per_window(CoverageWindowAction::Total);
    config.set_windows(DistributionWindowsArgs {
        by_size: Some(1_000_000),
        by_bed: None,
        by_grouped_bed: None,
    });

    assert_config_roundtrips!(config, Fcoverage, "fcoverage");
}

#[cfg(feature = "cmd_gc_bias")]
#[test]
fn gc_bias_config_roundtrips_through_rendered_cli() {
    let mut config = crate::commands::gc_bias::config::GCConfig::new(
        ioc(),
        PathBuf::from("ref.2bit"),
        PathBuf::from("ref_gc.zarr"),
        chromosomes(),
    );
    config.set_output_prefix("sample".to_string());

    assert_config_roundtrips!(config, GCBias, "gc-bias");
}

#[cfg(feature = "cmd_ref_gc_bias")]
#[test]
fn ref_gc_bias_config_roundtrips_through_rendered_cli() {
    use crate::commands::{
        cli_common::{FragmentLengthArgs, Ref2BitRequiredArgs},
        ref_gc_bias::config::{RefGCBiasConfig, RefGCWindowsArgs},
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

    assert_config_roundtrips!(config, RefGcBias, "ref-gc-bias");
}

#[cfg(feature = "cmd_lengths")]
#[test]
fn lengths_config_roundtrips_through_rendered_cli() {
    use crate::{
        commands::cli_common::DistributionWindowsArgs, commands::lengths::config::LengthsConfig,
    };

    let mut config = LengthsConfig::new(ioc(), chromosomes());
    config.output_prefix = "sample".to_string();
    config.set_windows(DistributionWindowsArgs {
        by_size: Some(1_000_000),
        by_bed: None,
        by_grouped_bed: None,
    });
    config.set_length_bins_spec("30:101:10");

    assert_config_roundtrips!(config, Lengths, "lengths");
}

#[cfg(feature = "cmd_midpoints")]
#[test]
fn midpoints_config_roundtrips_through_rendered_cli() {
    let mut config = crate::commands::midpoints::config::MidpointsConfig::new(
        ioc(),
        chromosomes(),
        PathBuf::from("sites.bed"),
    );
    config.set_output_prefix("sample");
    config.set_length_bins(vec![30, 80, 151]);

    assert_config_roundtrips!(config, Midpoints, "midpoints");
}

#[cfg(feature = "cmd_ends")]
#[test]
fn ends_config_roundtrips_through_rendered_cli() {
    let mut config = crate::commands::ends::config::EndsConfig::new(ioc(), chromosomes(), 2, 2);
    config.set_min_mapq(20);
    config.set_tile_size(1_000_000);

    assert_config_roundtrips!(config, Ends, "ends");
}

#[cfg(feature = "cmd_fragment_kmers")]
#[test]
fn fragment_kmers_config_roundtrips_through_rendered_cli() {
    use crate::commands::{
        cli_common::Ref2BitRequiredArgs, fragment_kmers::config::FragmentKmersConfig,
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

    assert_config_roundtrips!(config, FragmentKmers, "fragment-kmers");
}

#[cfg(feature = "cmd_transitions")]
#[test]
fn transitions_config_roundtrips_through_rendered_cli() {
    use crate::commands::{
        cli_common::Ref2BitRequiredArgs, transitions::config::TransitionsConfig,
    };

    let mut config = TransitionsConfig::new(
        ioc(),
        Ref2BitRequiredArgs {
            ref_2bit: PathBuf::from("ref.2bit"),
        },
        chromosomes(),
    );
    config.set_output_prefix("sample".to_string());
    config.set_orders(vec![1, 3]);

    assert_config_roundtrips!(config, Transitions, "transitions");
}

#[cfg(feature = "cmd_wps")]
#[test]
fn wps_config_roundtrips_through_rendered_cli() {
    let mut config = crate::commands::wps::config::WPSConfig::new(ioc(), chromosomes(), None);
    config.set_output_prefix("sample".to_string());
    config.set_window_size(121);

    assert_config_roundtrips!(config, WPS, "wps");
}

#[cfg(feature = "cmd_wps_peaks")]
#[test]
fn wps_peaks_config_roundtrips_through_rendered_cli() {
    use crate::commands::wps_peaks::{
        config::WPSPeaksConfig, window_peak_results::PeaksWindowAction,
    };

    let mut config = WPSPeaksConfig::new(ioc(), chromosomes(), Some(PeaksWindowAction::Stats));
    config.set_window_size(121);
    config.set_min_peak_height(6.5);

    assert_config_roundtrips!(config, WPSPeaks, "wps-peaks");
}

#[cfg(feature = "cmd_prepare_windows")]
#[test]
fn prepare_windows_config_roundtrips_through_rendered_cli() {
    let config = crate::commands::prepare_windows::config::PrepareConfig {
        input: PathBuf::from("windows.tsv"),
        output: PathBuf::from("prepared.bed"),
        group_cols: vec!["3".to_string()],
        ..crate::commands::prepare_windows::config::PrepareConfig::default()
    };

    assert_config_roundtrips!(config, PrepWindows, "prep-windows");
}

#[cfg(feature = "cmd_visualize_positions")]
#[test]
fn visualize_positions_config_roundtrips_through_rendered_cli() {
    use crate::{
        commands::{
            cli_common::{BaseSelectionArgs, FragmentPositionSelectionArgs},
            visualize_positions::{config::VisualizePositionsConfig, model::Style},
        },
        shared::positioning::{BasesFrom, MismatchBasesFrom},
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

    assert_config_roundtrips!(config, VisualizePositions, "visualize-positions");
}
