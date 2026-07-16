#[cfg(has_cli_commands)]
use crate::command_run::RunOptions;
#[cfg(feature = "cmd_bam_to_bam")]
use crate::commands::bam_to_bam::config::BamToBamConfig;
#[cfg(feature = "cmd_bam_to_frag")]
use crate::commands::bam_to_frag::config::BamToFragConfig;
#[cfg(feature = "cmd_coverage_weights")]
use crate::commands::coverage_weights::config::CoverageWeightsConfig;
#[cfg(feature = "cmd_ends")]
use crate::commands::ends::config::EndsConfig;
#[cfg(feature = "cmd_fcoverage")]
use crate::commands::fcoverage::config::FCoverageConfig;
#[cfg(feature = "cmd_frag_to_bam")]
use crate::commands::frag_to_bam::config::FragToBamConfig;
#[cfg(feature = "cmd_fragment_count_weights")]
use crate::commands::fragment_count_weights::config::FragmentCountWeightsConfig;
#[cfg(feature = "cmd_fragment_kmers")]
use crate::commands::fragment_kmers::config::FragmentKmersConfig;
#[cfg(feature = "cmd_gc_bias")]
use crate::commands::gc_bias::config::GCConfig;
#[cfg(feature = "cmd_lengths")]
use crate::commands::lengths::config::LengthsConfig;
#[cfg(feature = "cmd_midpoints")]
use crate::commands::midpoints::config::MidpointsConfig;
#[cfg(feature = "cmd_prepare_windows")]
use crate::commands::prepare_windows::config::PrepareConfig;
#[cfg(feature = "cmd_gc_bias")]
use crate::commands::ref_gc_bias::config::RefGCBiasConfig;
#[cfg(feature = "cmd_ref_kmers")]
use crate::commands::ref_kmers::config::RefKmersConfig;
#[cfg(feature = "cmd_transitions")]
use crate::commands::transitions::config::TransitionsConfig;
#[cfg(feature = "cmd_visualize_positions")]
use crate::commands::visualize_positions::config::VisualizePositionsConfig;
#[cfg(feature = "cmd_wps")]
use crate::commands::wps::config::WPSConfig;
#[cfg(feature = "cmd_wps_peaks")]
use crate::commands::wps_peaks::config::WPSPeaksConfig;
use clap::CommandFactory;
use clap::builder::styling::{AnsiColor, Style, Styles};

pub(crate) const CLI_SEPARATOR_WIDTH: usize = 48;

#[cfg(all(feature = "cli", not(has_cli_commands)))]
compile_error!("Building the CLI requires enabling at least one cmd_* feature.");

#[cfg_attr(feature = "cli", derive(clap::Parser))]
#[command(name = "cfdna", version, about = env!("CARGO_PKG_DESCRIPTION"))]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub(crate) cmd: Cmd,
}

#[cfg_attr(feature = "cli", derive(clap::Subcommand))]
pub(crate) enum Cmd {
    #[cfg(feature = "cmd_gc_bias")]
    GCBias(GCConfig),
    #[cfg(feature = "cmd_gc_bias")]
    RefGcBias(RefGCBiasConfig),
    #[cfg(feature = "cmd_ref_kmers")]
    RefKmers(RefKmersConfig),
    #[cfg(feature = "cmd_transitions")]
    Transitions(TransitionsConfig),
    #[cfg(feature = "cmd_coverage_weights")]
    CoverageWeights(CoverageWeightsConfig),
    #[cfg(feature = "cmd_fragment_count_weights")]
    FragmentCountWeights(FragmentCountWeightsConfig),
    #[cfg(feature = "cmd_ends")]
    Ends(EndsConfig),
    #[cfg(feature = "cmd_lengths")]
    Lengths(LengthsConfig),
    #[cfg(feature = "cmd_fcoverage")]
    Fcoverage(FCoverageConfig),
    #[cfg(feature = "cmd_wps")]
    WPS(WPSConfig),
    #[cfg(feature = "cmd_wps_peaks")]
    WPSPeaks(WPSPeaksConfig),
    #[cfg(feature = "cmd_midpoints")]
    Midpoints(MidpointsConfig),
    #[cfg(feature = "cmd_fragment_kmers")]
    FragmentKmers(FragmentKmersConfig),
    #[cfg(feature = "cmd_prepare_windows")]
    PrepWindows(PrepareConfig),
    #[cfg(feature = "cmd_visualize_positions")]
    VisualizePositions(VisualizePositionsConfig),
    #[cfg(feature = "cmd_bam_to_bam")]
    BamToBam(BamToBamConfig),
    #[cfg(feature = "cmd_bam_to_frag")]
    BamToFrag(BamToFragConfig),
    #[cfg(feature = "cmd_frag_to_bam")]
    FragToBam(FragToBamConfig),
}

#[cfg(all(feature = "cli", has_cli_commands))]
pub(crate) fn run_cli() {
    use clap::FromArgMatches;

    #[cfg(uses_temp_dirs)]
    if crate::shared::tiled_run::run_temp_dir_cleanup_helper_if_requested() {
        return;
    }

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
        #[cfg(feature = "cmd_gc_bias")]
        Cmd::RefGcBias(config) => (
            config.logging.log.clone(),
            Some(config.output_dir.as_path()),
        ),
        #[cfg(feature = "cmd_ref_kmers")]
        Cmd::RefKmers(config) => (
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
        Cmd::PrepWindows(_config) => (crate::shared::logging::LogSpec::Stdout, None),
        #[cfg(feature = "cmd_visualize_positions")]
        Cmd::VisualizePositions(_config) => (crate::shared::logging::LogSpec::Stdout, None),
    };

    if let Err(error) =
        crate::shared::logging::init_cli_logging(&command_name, &log_spec, default_output_dir)
    {
        eprintln!("{:#}", error);
        std::process::exit(1);
    }
    let run_options = if matches!(log_spec, crate::shared::logging::LogSpec::Quiet) {
        RunOptions::new_quiet()
    } else {
        RunOptions::new_cli()
    };

    crate::shared::cli_output::print_command_banner(&command_name);
    let result: anyhow::Result<()> = match cli.cmd {
        #[cfg(feature = "cmd_gc_bias")]
        Cmd::GCBias(config) => {
            crate::commands::gc_bias::gc_bias::run_gc_bias(&config, run_options).map(|_| ())
        }
        #[cfg(feature = "cmd_gc_bias")]
        Cmd::RefGcBias(config) => {
            crate::commands::ref_gc_bias::ref_gc_bias::run_ref_gc_bias(&config, run_options)
                .map(|_| ())
        }
        #[cfg(feature = "cmd_ref_kmers")]
        Cmd::RefKmers(config) => {
            crate::commands::ref_kmers::ref_kmers::run_ref_kmers(&config, run_options).map(|_| ())
        }
        #[cfg(feature = "cmd_transitions")]
        Cmd::Transitions(config) => {
            crate::commands::transitions::transitions::run_transitions(&config, run_options)
                .map(|_| ())
        }
        #[cfg(feature = "cmd_coverage_weights")]
        Cmd::CoverageWeights(config) => {
            crate::commands::coverage_weights::coverage_weights::run_coverage_weights(
                &config,
                run_options,
            )
            .map(|_| ())
        }
        #[cfg(feature = "cmd_fragment_count_weights")]
        Cmd::FragmentCountWeights(config) => {
            crate::commands::fragment_count_weights::fragment_count_weights::run_fragment_count_weights(
                &config,
                run_options,
            )
            .map(|_| ())
        }
        #[cfg(feature = "cmd_ends")]
        Cmd::Ends(config) => {
            crate::commands::ends::ends::run_ends(&config, run_options).map(|_| ())
        }
        #[cfg(feature = "cmd_lengths")]
        Cmd::Lengths(config) => {
            crate::commands::lengths::lengths::run_lengths(&config, run_options).map(|_| ())
        }
        #[cfg(feature = "cmd_fcoverage")]
        Cmd::Fcoverage(config) => {
            crate::commands::fcoverage::fcoverage::run_fcoverage(&config, run_options).map(|_| ())
        }
        #[cfg(feature = "cmd_wps")]
        Cmd::WPS(config) => {
            crate::commands::wps::wps::run_wps(&config, run_options).map(|_| ())
        }
        #[cfg(feature = "cmd_wps_peaks")]
        Cmd::WPSPeaks(config) => {
            crate::commands::wps_peaks::wps_peaks::run_wps_peaks(&config, run_options)
                .map(|_| ())
        }
        #[cfg(feature = "cmd_midpoints")]
        Cmd::Midpoints(config) => {
            crate::commands::midpoints::midpoints::run_midpoints(&config, run_options)
                .map(|_| ())
        }
        #[cfg(feature = "cmd_fragment_kmers")]
        Cmd::FragmentKmers(config) => {
            crate::commands::fragment_kmers::fragment_kmers::run_fragment_kmers(
                &config,
                run_options,
            )
            .map(|_| ())
        }
        #[cfg(feature = "cmd_prepare_windows")]
        Cmd::PrepWindows(config) => {
            crate::commands::prepare_windows::prepare_windows::run_prepare_windows(
                &config,
                run_options,
            )
            .map(|_| ())
        }
        #[cfg(feature = "cmd_visualize_positions")]
        Cmd::VisualizePositions(config) => {
            crate::commands::visualize_positions::visualize_positions::run_visualize_positions(
                &config,
                run_options,
            )
            .map(|_| ())
        }
        #[cfg(feature = "cmd_bam_to_bam")]
        Cmd::BamToBam(config) => {
            crate::commands::bam_to_bam::bam_to_bam::run_bam_to_bam(&config, run_options)
                .map(|_| ())
        }
        #[cfg(feature = "cmd_bam_to_frag")]
        Cmd::BamToFrag(config) => {
            crate::commands::bam_to_frag::bam_to_frag::run_bam_to_frag(&config, run_options)
                .map(|_| ())
        }
        #[cfg(feature = "cmd_frag_to_bam")]
        Cmd::FragToBam(config) => {
            crate::commands::frag_to_bam::frag_to_bam::run_frag_to_bam(&config, run_options)
                .map(|_| ())
        }
    };
    crate::shared::cli_output::print_command_footer();

    if let Err(error) = result {
        let rendered_error = format!("{:#}", error);
        eprintln!("{}", rendered_error);
        crate::shared::logging::duplicate_stderr_line_to_file(&rendered_error);
        std::process::exit(1);
    }
}

#[cfg(all(feature = "cli", not(has_cli_commands)))]
pub(crate) fn run_cli() {
    unreachable!("the CLI requires at least one cmd_* feature");
}

/// Build terminal-oriented clap command with sanitized docs and branded signature
pub(crate) fn build_terminal_command() -> clap::Command {
    let mut command = Cli::command();
    let styles = Styles::styled()
        .header(AnsiColor::Yellow.on_default().bold())
        .usage(AnsiColor::Green.on_default().bold())
        .literal(AnsiColor::Blue.on_default().bold())
        .placeholder(AnsiColor::Cyan.on_default());
    command = command
        .help_template(
            "{before-help}{name} {version}\n{about}\n\n{usage-heading} {usage}\n\n{all-args}\n",
        )
        .styles(styles);
    command = sanitize_command(command);
    let signature = terminal_signature();
    add_signature(command, &signature, "cfdna")
}

/// Build docs-oriented clap command with raw help text
pub(crate) fn build_docs_command() -> clap::Command {
    Cli::command()
        .help_template("{name} {version}\n{about}\n\n{usage-heading} {usage}\n\n{all-args}\n")
}

/// Sanitize markdown-ish help for terminal rendering
fn sanitize_cli_text(text: &str) -> String {
    let mut output = String::with_capacity(text.len());
    let mut in_block = false;

    let block = Style::new().dimmed();
    let block_on = format!("{block}");
    let block_off = format!("{block:#}");

    for line in text.lines() {
        let trimmed = line.trim_start();
        if trimmed.starts_with("```") {
            in_block = !in_block;
            continue;
        }

        let normalized_line = line
            .replace('→', "->")
            .replace('’', "'")
            .replace(['“', '”'], "\"");

        if in_block {
            output.push_str("  ");
            output.push_str(&block_on);
            output.push_str(&normalized_line);
            output.push_str(&block_off);
        } else if let Some(rendered_header) = try_render_header(&normalized_line) {
            output.push_str(&rendered_header);
        } else {
            output.push_str(&stylize_inline(&normalized_line));
        }
        output.push('\n');
    }

    if output.ends_with('\n') {
        output.pop();
    }
    output
}

/// Turn **bold** and `inline code` markers into styled ANSI text
fn stylize_inline(line: &str) -> String {
    let mut output = String::with_capacity(line.len());
    let bold = Style::new().bold();
    let code = Style::new().dimmed().underline();
    let bold_on = format!("{bold}");
    let bold_off = format!("{bold:#}");
    let code_on = format!("{code}");
    let code_off = format!("{code:#}");

    let bytes = line.as_bytes();
    let mut index = 0usize;
    let mut in_bold = false;
    let mut in_code = false;

    while index < bytes.len() {
        if !in_code && index + 1 < bytes.len() && bytes[index] == b'*' && bytes[index + 1] == b'*' {
            if in_bold {
                output.push_str(&bold_off);
            } else {
                output.push_str(&bold_on);
            }
            in_bold = !in_bold;
            index += 2;
            continue;
        }
        if bytes[index] == b'`' {
            if in_code {
                output.push_str(&code_off);
            } else {
                output.push_str(&code_on);
            }
            in_code = !in_code;
            index += 1;
            continue;
        }
        output.push(bytes[index] as char);
        index += 1;
    }

    if in_code {
        output.push_str(&code_off);
    }
    if in_bold {
        output.push_str(&bold_off);
    }
    output
}

/// Render markdown-like headers on non-fenced lines
fn try_render_header(line: &str) -> Option<String> {
    let trimmed = line.trim_start();
    let mut hashes = 0usize;
    for byte in trimmed.as_bytes() {
        if *byte == b'#' {
            hashes += 1;
        } else {
            break;
        }
    }
    if hashes == 0 {
        return None;
    }
    let rest = &trimmed[hashes..];
    if !rest.starts_with(' ') {
        return None;
    }
    let text = rest.trim_start();

    let (style, underline): (Style, bool) = match hashes {
        1 => (
            Style::new().fg_color(Some(AnsiColor::Yellow.into())).bold(),
            true,
        ),
        2 => (Style::new().bold(), false),
        3 => (Style::new().underline(), false),
        _ => (Style::new().dimmed(), false),
    };

    let header_text = if hashes == 2 {
        text.to_uppercase()
    } else {
        text.to_string()
    };

    let style_on = format!("{style}");
    let style_off = format!("{style:#}");
    let mut output = String::with_capacity(line.len() + 32);
    if hashes == 2 {
        output.push('\n');
    }
    output.push_str(&style_on);
    output.push_str(&header_text);
    output.push_str(&style_off);
    if hashes != 2 {
        output.push('\n');
    }

    if underline {
        let bar_len = header_text.chars().count().min(64);
        output.push_str(&"─".repeat(bar_len));
    }
    Some(output)
}

/// Sanitize help/long_help for command, args, and all subcommands
fn sanitize_command(mut command: clap::Command) -> clap::Command {
    if let Some(about) = command.get_about().map(|value| value.to_string()) {
        command = command.about(sanitize_cli_text(&about));
    }
    if let Some(long_about) = command.get_long_about().map(|value| value.to_string()) {
        command = command.long_about(sanitize_cli_text(&long_about));
    }

    let argument_infos: Vec<(clap::Id, Option<String>, Option<String>)> = command
        .get_arguments()
        .map(|argument| {
            let id = argument.get_id().clone();
            let help = argument.get_help().map(|value| value.to_string());
            let long_help = argument.get_long_help().map(|value| value.to_string());
            (id, help, long_help)
        })
        .collect();

    for (argument_id, help, long_help) in argument_infos {
        if let Some(help_text) = help {
            let cleaned_help = sanitize_cli_text(&help_text);
            command = command.mut_arg(&argument_id, |argument| argument.help(cleaned_help));
        }
        if let Some(long_help_text) = long_help {
            let cleaned_long_help = sanitize_cli_text(&long_help_text);
            command = command.mut_arg(&argument_id, |argument| {
                argument.long_help(cleaned_long_help)
            });
        }
    }

    let subcommand_names: Vec<String> = command
        .get_subcommands()
        .map(|subcommand| subcommand.get_name().to_string())
        .collect();

    for subcommand_name in subcommand_names {
        command = command.mut_subcommand(&subcommand_name, sanitize_command);
    }

    command
}

/// Build the branded terminal signature shown in CLI help and command banners.
pub(crate) fn terminal_signature() -> String {
    let accent = Style::new().bold();
    terminal_signature_with_bars(&format!("{accent}"), &format!("{accent:#}"))
}

/// Build the plain-text terminal signature for non-terminal sinks such as log files.
pub(crate) fn plain_terminal_signature() -> String {
    terminal_signature_with_bars("", "")
}

fn terminal_signature_with_bars(style_on: &str, style_off: &str) -> String {
    let title = "cfDNAlab";
    let bar1 = "_".repeat(CLI_SEPARATOR_WIDTH);
    let bar2 = "─".repeat(CLI_SEPARATOR_WIDTH);
    format!("\n{style_on}{bar1}\n\n  {title}\n\n{bar2}{style_off}\n")
}

/// Apply the branded help header to the command and all subcommands.
///
/// Each header includes the full command path so standalone help output identifies which command
/// it describes. Building the path during recursion also keeps nested subcommands correct without
/// requiring command-specific help configuration.
fn add_signature(mut command: clap::Command, signature: &str, command_path: &str) -> clap::Command {
    let help_header = format!(
        "{signature}Command: {command_path}\n{}\n",
        "─".repeat(CLI_SEPARATOR_WIDTH)
    );
    command = command
        .before_help(help_header.clone())
        .before_long_help(help_header);

    let subcommand_names: Vec<String> = command
        .get_subcommands()
        .map(|subcommand| subcommand.get_name().to_string())
        .collect();

    for subcommand_name in subcommand_names {
        let subcommand_path = format!("{command_path} {subcommand_name}");
        command = command.mut_subcommand(&subcommand_name, |subcommand| {
            add_signature(subcommand, signature, &subcommand_path)
        });
    }
    command
}

#[cfg(all(test, feature = "cli"))]
#[path = "cli_app_roundtrip_tests.rs"]
mod cli_app_roundtrip_tests;

#[cfg(all(test, feature = "cli", feature = "cmd_ref_kmers"))]
#[path = "cli_app_tests.rs"]
mod cli_app_tests;
