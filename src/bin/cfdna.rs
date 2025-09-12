use cfdnalab::gc::GCConfig;
use cfdnalab::lengths::LengthsConfig;
use cfdnalab::normalize_genome::NormalizeGenomeConfig;
use cfdnalab::refgc::RefGCConfig;
use clap::{CommandFactory, Parser, Subcommand};

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
    // Build command from your derive
    let mut cmd = Cli::command();

    // Sanitize help/long_help pulled from your doc comments
    sanitize_command_help(&mut cmd);

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

/// Minimal Markdown -> terminal text sanitizer for CLI help
fn sanitize_cli_text(md: &str) -> String {
    let mut out = String::with_capacity(md.len());
    let mut in_code = false;

    for line in md.lines() {
        let trimmed = line.trim_start();

        // Toggle code fences ```lang
        if trimmed.starts_with("```") {
            in_code = !in_code;
            continue;
        }

        let mut s = line.to_string();

        if !in_code {
            // Drop inline backticks
            s = s.replace('`', "");
            // Replace fancy arrows with ->
            s = s.replace('→', "->");
            // Normalize curly quotes
            s = s.replace('’', "'").replace('“', "\"").replace('”', "\"");
        }

        if in_code {
            out.push_str("  "); // Indent code lines a bit
        }
        out.push_str(&s);
        out.push('\n');
    }

    if out.ends_with('\n') {
        out.pop();
    }
    out
}

/// Recursively apply sanitizer to a clap::Command tree
fn sanitize_command_help(cmd: &mut clap::Command) {
    if let Some(a) = cmd.get_about().map(|s| s.to_string()) {
        cmd.about(sanitize_cli_text(&a));
    }
    if let Some(a) = cmd.get_long_about().map(|s| s.to_string()) {
        cmd.long_about(sanitize_cli_text(&a));
    }
    for arg in cmd.get_arguments_mut() {
        if let Some(h) = arg.get_help().map(|s| s.to_string()) {
            arg.help(sanitize_cli_text(&h));
        }
        if let Some(h) = arg.get_long_help().map(|s| s.to_string()) {
            arg.long_help(sanitize_cli_text(&h));
        }
    }
    for sub in cmd.get_subcommands_mut() {
        sanitize_command_help(sub);
    }
}
