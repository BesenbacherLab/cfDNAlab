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

/// Minimal Markdown -> terminal cleanup for CLI help
fn sanitize_cli_text(md: &str) -> String {
    let mut out = String::with_capacity(md.len());
    let mut in_code = false;
    for line in md.lines() {
        let trimmed = line.trim_start();
        if trimmed.starts_with("```") {
            in_code = !in_code;
            continue;
        }
        let mut s = line.to_string();
        if !in_code {
            s = s.replace('`', "");
            s = s.replace('→', "->");
            s = s.replace('’', "'").replace('“', "\"").replace('”', "\"");
        }
        if in_code {
            out.push_str("  ");
        } // indent code lines a bit
        out.push_str(&s);
        out.push('\n');
    }
    if out.ends_with('\n') {
        out.pop();
    }
    out
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
