use cfdnalab::gc::GCConfig;
use cfdnalab::lengths::LengthsConfig;
use cfdnalab::normalize_genome::NormalizeGenomeConfig;
use cfdnalab::refgc::RefGCConfig;
use clap::builder::styling::{AnsiColor, Styles};
use clap::{CommandFactory, FromArgMatches, Parser, Subcommand};

#[derive(Parser)]
#[command(name = "cfdna", version)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    GC(GCConfig),
    RefGC(RefGCConfig), // Extract reference GC counts
    NormalizeGenome(NormalizeGenomeConfig),
    Lengths(LengthsConfig),
    // Ends(EndsConfig),
}

fn main() {
    // Build Command from derive
    let mut cmd0 = Cli::command();

    // Optionally set styles/template here on cmd0 before sanitizing

    let styles = Styles::styled()
        .header(AnsiColor::Yellow.on_default().bold())
        .usage(AnsiColor::Green.on_default().bold())
        .literal(AnsiColor::Blue.on_default().bold())
        .placeholder(AnsiColor::Cyan.on_default());

    cmd0 = cmd0
        .help_template("{name} {version}\n{about}\n\n{usage-heading} {usage}\n\n{all-args}\n")
        .styles(styles);

    // Sanitize help/long_help pulled from your doc comments
    let cmd = sanitize_command(cmd0);

    // Parse using the sanitized command
    let matches = cmd.clone().get_matches();
    let cli = Cli::from_arg_matches(&matches).expect("parse");

    // Run selected subcommand and capture its Result (no `?` in main).
    let res: anyhow::Result<()> = match cli.cmd {
        Cmd::GC(cfg) => cfdnalab::gc::run(cfg),
        Cmd::RefGC(cfg) => cfdnalab::refgc::run(cfg),
        Cmd::NormalizeGenome(cfg) => cfdnalab::normalize_genome::run(cfg),
        Cmd::Lengths(cfg) => cfdnalab::lengths::run(cfg),
        // Cmd::Ends(cfg) => cfdnalab::ends::run(cfg),
    };

    if let Err(e) = res {
        eprintln!("{:#}", e);
        std::process::exit(1);
    }

    std::process::exit(0);
}

/// Sanitize help/long_help for a Command, its args, and all subcommands.
/// NOTE: Takes and RETURNS ownership to avoid borrow/move errors with clap's builder API.
fn sanitize_command(mut cmd: clap::Command) -> clap::Command {
    // Sanitize about / long_about (extract first to break borrows)
    if let Some(a) = cmd.get_about().map(|s| s.to_string()) {
        cmd = cmd.about(sanitize_cli_text(&a));
    }
    if let Some(a) = cmd.get_long_about().map(|s| s.to_string()) {
        cmd = cmd.long_about(sanitize_cli_text(&a));
    }

    // Collect arg IDs and their help strings up front (read-only pass)
    let arg_infos: Vec<(clap::Id, Option<String>, Option<String>)> = cmd
        .get_arguments()
        .map(|a| {
            let id = a.get_id().clone();
            let h = a.get_help().map(|s| s.to_string());
            let lh = a.get_long_help().map(|s| s.to_string());
            (id, h, lh)
        })
        .collect();

    // Rebuild args using mut_arg (consumes and returns Command)
    for (id, h, lh) in arg_infos {
        if let Some(hs) = h {
            let cleaned = sanitize_cli_text(&hs);
            cmd = cmd.mut_arg(&id, |a| a.help(cleaned));
        }
        if let Some(lhs) = lh {
            let cleaned = sanitize_cli_text(&lhs);
            cmd = cmd.mut_arg(&id, |a| a.long_help(cleaned));
        }
    }

    // Recurse into subcommands using mut_subcommand (also consumes self)
    let sub_names: Vec<String> = cmd
        .get_subcommands()
        .map(|sc| sc.get_name().to_string())
        .collect();

    for name in sub_names {
        cmd = cmd.mut_subcommand(&name, |sub| sanitize_command(sub));
    }

    cmd
}

/// Turn **bold** and `inline code` markers into styled ANSI text (not Markdown)
fn stylize_inline(s: &str) -> String {
    let mut out = String::with_capacity(s.len());

    // choose your styles
    let bold = Style::new().bold();
    let code = Style::new().dimmed().underline();

    let bold_on = format!("{bold}");
    let bold_off = format!("{bold:#}");
    let code_on = format!("{code}");
    let code_off = format!("{code:#}");

    let bytes = s.as_bytes();
    let mut i = 0usize;
    let mut in_bold = false;
    let mut in_code = false;

    while i < bytes.len() {
        // **bold**
        if !in_code && i + 1 < bytes.len() && bytes[i] == b'*' && bytes[i + 1] == b'*' {
            if in_bold {
                out.push_str(&bold_off);
            } else {
                out.push_str(&bold_on);
            }
            in_bold = !in_bold;
            i += 2;
            continue;
        }
        // `code`
        if bytes[i] == b'`' {
            if in_code {
                out.push_str(&code_off);
            } else {
                out.push_str(&code_on);
            }
            in_code = !in_code;
            i += 1;
            continue;
        }
        out.push(bytes[i] as char);
        i += 1;
    }

    // close any unclosed spans
    if in_code {
        out.push_str(&code_off);
    }
    if in_bold {
        out.push_str(&bold_off);
    }
    out
}

/// Sanitize Markdown-ish help for terminals:
/// - Treat ``` fences as block code (no inline styling inside)
/// - Apply inline **bold** and `code` elsewhere
/// - Normalize arrows/quotes
pub fn sanitize_cli_text(md: &str) -> String {
    let mut out = String::with_capacity(md.len());
    let mut in_block = false;

    // style used for code blocks (distinct from inline code)
    let block = Style::new().dimmed();
    let block_on = format!("{block}");
    let block_off = format!("{block:#}");

    for line in md.lines() {
        let trimmed = line.trim_start();

        // toggle code-block mode on lines that start with ``` (any language tag)
        if trimmed.starts_with("```") {
            in_block = !in_block;
            continue; // drop the fence line itself
        }

        // normalize a few typography chars for terminals
        let line = line
            .replace('→', "->")
            .replace('’', "'")
            .replace('“', "\"")
            .replace('”', "\"");

        if in_block {
            // in a fenced block: don't parse inline markers; optionally style/indent
            out.push_str("  "); // simple indent
            out.push_str(&block_on);
            out.push_str(&line);
            out.push_str(&block_off);
        } else {
            // outside a block: apply inline styling (**bold**, `code`)
            out.push_str(&stylize_inline(&line));
        }
        out.push('\n');
    }

    if out.ends_with('\n') {
        out.pop();
    }
    out
}
