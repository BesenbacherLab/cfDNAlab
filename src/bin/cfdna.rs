use cfdnalab::{
    cli_app::{Cli, Cmd, build_terminal_command},
    shared::{cli_output, logging},
};

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
    let command_name = matches.subcommand_name().unwrap_or("help").to_string();
    let cli = Cli::from_arg_matches(&matches).expect("parse");

    let (log_spec, default_output_dir) = match &cli.cmd {
        #[cfg(feature = "cmd_coverage_weights")]
        Cmd::CoverageWeights(config) => (
            config.shared.logging.log.clone(),
            Some(config.shared.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_fragment_count_weights")]
        Cmd::FragmentCountWeights(config) => (
            config.shared.logging.log.clone(),
            Some(config.shared.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_fcoverage")]
        Cmd::Fcoverage(config) => (
            config.logging.log.clone(),
            Some(config.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_gc_bias")]
        Cmd::GCBias(config) => (
            config.logging.log.clone(),
            Some(config.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_ref_gc_bias")]
        Cmd::RefGcBias(config) => (
            config.logging.log.clone(),
            Some(config.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_transitions")]
        Cmd::Transitions(config) => (
            config.shared_args.logging.log.clone(),
            Some(config.shared_args.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_ends")]
        Cmd::Ends(config) => (
            config.logging.log.clone(),
            Some(config.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_lengths")]
        Cmd::Lengths(config) => (
            config.logging.log.clone(),
            Some(config.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_wps")]
        Cmd::WPS(config) => (
            config.shared_args.logging.log.clone(),
            Some(config.shared_args.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_wps_peaks")]
        Cmd::WPSPeaks(config) => (
            config.shared_args.logging.log.clone(),
            Some(config.shared_args.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_midpoints")]
        Cmd::Midpoints(config) => (
            config.logging.log.clone(),
            Some(config.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_fragment_kmers")]
        Cmd::FragmentKmers(config) => (
            config.shared_args.logging.log.clone(),
            Some(config.shared_args.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_bam_to_bam")]
        Cmd::BamToBam(config) => (config.logging.log.clone(), config.out_bam.parent()),
        #[cfg(feature = "cmd_bam_to_frag")]
        Cmd::BamToFrag(config) => (
            config.logging.log.clone(),
            Some(config.ioc.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_frag_to_bam")]
        Cmd::FragToBam(config) => (
            config.logging.log.clone(),
            Some(config.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_prepare_windows")]
        Cmd::PrepWindows(_config) => (logging::LogSpec::Stdout, None),
        #[cfg(feature = "cmd_visualize_positions")]
        Cmd::VisualizePositions(_config) => (logging::LogSpec::Stdout, None),
    };

    if let Err(error) = logging::init_cli_logging(&command_name, &log_spec, default_output_dir) {
        eprintln!("{:#}", error);
        std::process::exit(1);
    }

    cli_output::print_command_banner(&command_name);
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
        #[cfg(feature = "cmd_fragment_count_weights")]
        Cmd::FragmentCountWeights(config) => {
            cfdnalab::commands::fragment_count_weights::fragment_count_weights::run(&config)
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
    cli_output::print_command_footer();

    if let Err(error) = result {
        let rendered_error = format!("{:#}", error);
        eprintln!("{}", rendered_error);
        logging::duplicate_stderr_line_to_file(&rendered_error);
        std::process::exit(1);
    }
}
