#[cfg(feature = "cmd_bam_to_bam")]
use crate::commands::bam_to_bam::config::BamToBamConfig;
#[cfg(feature = "cmd_bam_to_frag")]
use crate::commands::bam_to_frag::config::BamToFragConfig;
#[cfg(feature = "cmd_coverage_weights")]
use crate::commands::coverage_weights::config::CoverageWeightsConfig;
#[cfg(feature = "cmd_fcoverage")]
use crate::commands::fcoverage::config::FCoverageConfig;
#[cfg(feature = "cmd_frag_to_bam")]
use crate::commands::frag_to_bam::config::FragToBamConfig;
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
#[cfg(feature = "cmd_ref_gc_bias")]
use crate::commands::ref_gc_bias::config::RefGCBiasConfig;
#[cfg(feature = "cmd_visualize_positions")]
use crate::commands::visualize_positions::config::VisualizePositionsConfig;
#[cfg(feature = "cmd_wps")]
use crate::commands::wps::config::WPSConfig;
#[cfg(feature = "cmd_wps_peaks")]
use crate::commands::wps_peaks::config::WPSPeaksConfig;
use clap::CommandFactory;
use clap::builder::styling::{AnsiColor, Style, Styles};

#[cfg(all(
    feature = "cli",
    not(any(
        feature = "cmd_bam_to_bam",
        feature = "cmd_bam_to_frag",
        feature = "cmd_frag_to_bam",
        feature = "cmd_coverage_weights",
        feature = "cmd_fcoverage",
        feature = "cmd_fragment_kmers",
        feature = "cmd_gc_bias",
        feature = "cmd_lengths",
        feature = "cmd_prepare_windows",
        feature = "cmd_midpoints",
        feature = "cmd_ref_gc_bias",
        feature = "cmd_visualize_positions",
        feature = "cmd_wps",
        feature = "cmd_wps_peaks"
    ))
))]
compile_error!("Building the CLI requires enabling at least one cmd_* feature.");

#[cfg_attr(feature = "cli", derive(clap::Parser))]
#[command(name = "cfdna", version)]
pub struct Cli {
    #[command(subcommand)]
    pub cmd: Cmd,
}

#[cfg_attr(feature = "cli", derive(clap::Subcommand))]
pub enum Cmd {
    #[cfg(feature = "cmd_gc_bias")]
    GCBias(GCConfig),
    #[cfg(feature = "cmd_ref_gc_bias")]
    RefGcBias(RefGCBiasConfig),
    #[cfg(feature = "cmd_coverage_weights")]
    CoverageWeights(CoverageWeightsConfig),
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

/// Build terminal-oriented clap command with sanitized docs and branded signature
pub fn build_terminal_command() -> clap::Command {
    let mut command = Cli::command();
    let styles = Styles::styled()
        .header(AnsiColor::Yellow.on_default().bold())
        .usage(AnsiColor::Green.on_default().bold())
        .literal(AnsiColor::Blue.on_default().bold())
        .placeholder(AnsiColor::Cyan.on_default());
    command = command
        .help_template("{name} {version}\n{about}\n\n{usage-heading} {usage}\n\n{all-args}\n")
        .styles(styles);
    command = sanitize_command(command);
    let signature = make_signature();
    add_signature(command, &signature)
}

/// Build docs-oriented clap command with raw help text
pub fn build_docs_command() -> clap::Command {
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

/// Build a first-line terminal signature
fn make_signature() -> String {
    let accent = Style::new().bold();
    let title = "cfDNAlab";
    let bar1 = "_".repeat(48);
    let bar2 = "─".repeat(48);
    format!("\n{accent}{bar1}\n\n  {title}\n\n{bar2}{accent:#}\n")
}

/// Apply signature to command and all subcommands
fn add_signature(mut command: clap::Command, signature: &str) -> clap::Command {
    command = command
        .before_help(signature.to_string())
        .before_long_help(signature.to_string());

    let subcommand_names: Vec<String> = command
        .get_subcommands()
        .map(|subcommand| subcommand.get_name().to_string())
        .collect();

    for subcommand_name in subcommand_names {
        command = command.mut_subcommand(&subcommand_name, |subcommand| {
            add_signature(
                subcommand
                    .before_help(signature.to_string())
                    .before_long_help(signature.to_string()),
                signature,
            )
        });
    }
    command
}
