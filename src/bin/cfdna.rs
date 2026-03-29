use cfdnalab::cli_app::{Cli, Cmd, build_terminal_command};

#[cfg(not(feature = "cli"))]
fn main() {
    eprintln!("This binary requires --features cli");
    std::process::exit(1);
}

#[cfg(feature = "cli")]
fn main() {
    use clap::FromArgMatches;

    let command = build_terminal_command();
    let matches = command.clone().get_matches();
    let cli = Cli::from_arg_matches(&matches).expect("parse");

    let result: anyhow::Result<()> = match cli.cmd {
        #[cfg(feature = "cmd_gc_bias")]
        Cmd::GCBias(config) => cfdnalab::commands::gc_bias::gc_bias::run(&config),
        #[cfg(feature = "cmd_ref_gc_bias")]
        Cmd::RefGcBias(config) => cfdnalab::commands::ref_gc_bias::ref_gc_bias::run(&config),
        #[cfg(feature = "cmd_transitions")]
        Cmd::Transitions(config) => cfdnalab::commands::transitions::transitions::run(&config),
        #[cfg(feature = "cmd_coverage_weights")]
        Cmd::CoverageWeights(config) => {
            cfdnalab::commands::coverage_weights::coverage_weights::run(&config)
        }
        #[cfg(feature = "cmd_ends")]
        Cmd::Ends(config) => cfdnalab::commands::ends::ends::run(&config),
        #[cfg(feature = "cmd_lengths")]
        Cmd::Lengths(config) => cfdnalab::commands::lengths::lengths::run(&config),
        #[cfg(feature = "cmd_fcoverage")]
        Cmd::Fcoverage(config) => cfdnalab::commands::fcoverage::fcoverage::run(&config),
        #[cfg(feature = "cmd_wps")]
        Cmd::WPS(config) => cfdnalab::commands::wps::wps::run(&config),
        #[cfg(feature = "cmd_wps_peaks")]
        Cmd::WPSPeaks(config) => cfdnalab::commands::wps_peaks::wps_peaks::run(&config),
        #[cfg(feature = "cmd_midpoints")]
        Cmd::Midpoints(config) => cfdnalab::commands::midpoints::midpoints::run(&config),
        #[cfg(feature = "cmd_fragment_kmers")]
        Cmd::FragmentKmers(config) => {
            cfdnalab::commands::fragment_kmers::fragment_kmers::run(&config)
        }
        #[cfg(feature = "cmd_prepare_windows")]
        Cmd::PrepWindows(config) => {
            cfdnalab::commands::prepare_windows::prepare_windows::run(&config)
        }
        #[cfg(feature = "cmd_visualize_positions")]
        Cmd::VisualizePositions(config) => {
            cfdnalab::commands::visualize_positions::visualize_positions::run(&config)
        }
        #[cfg(feature = "cmd_bam_to_bam")]
        Cmd::BamToBam(config) => cfdnalab::commands::bam_to_bam::bam_to_bam::run(&config),
        #[cfg(feature = "cmd_bam_to_frag")]
        Cmd::BamToFrag(config) => cfdnalab::commands::bam_to_frag::bam_to_frag::run(&config),
        #[cfg(feature = "cmd_frag_to_bam")]
        Cmd::FragToBam(config) => cfdnalab::commands::frag_to_bam::frag_to_bam::run(&config),
    };

    if let Err(error) = result {
        eprintln!("{:#}", error);
        std::process::exit(1);
    }
}
