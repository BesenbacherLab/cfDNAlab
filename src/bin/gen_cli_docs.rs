use anyhow::{Context, Result, bail};
use cfdnalab::cli_app::build_docs_command;
use clap::Parser;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const GENERATED_MARKER: &str = "<!-- AUTO-GENERATED FILE - DO NOT EDIT -->";
const GENERATED_SOURCE: &str = "<!-- Source: cfdna Clap config and command tree -->";

#[derive(Debug, Clone, clap::ValueEnum)]
enum Scope {
    Release,
    All,
}

#[derive(Debug, Parser)]
#[command(name = "gen_cli_docs")]
#[command(about = "Generate CLI markdown docs for the docs website")]
struct Cli {
    /// Output directory for generated CLI markdown pages
    #[arg(long, value_parser)]
    out_dir: PathBuf,

    /// Command scope to include
    #[arg(long, value_enum, default_value = "release")]
    scope: Scope,

    /// Run git-diff validation after generation
    #[arg(long)]
    fail_on_drift: bool,
}

#[derive(Debug, Clone)]
struct CommandDoc {
    name: String,
    title: String,
    help_text: String,
}

#[cfg(not(all(feature = "cli", feature = "docs_gen")))]
fn main() {
    eprintln!("This binary requires --features cli,docs_gen");
    std::process::exit(1);
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn main() -> Result<()> {
    let args = Cli::parse();
    run(&args)
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn run(args: &Cli) -> Result<()> {
    let root_command = build_docs_command();
    let command_names = match args.scope {
        Scope::Release => release_command_names(),
        Scope::All => all_command_names(&root_command),
    };

    fs::create_dir_all(&args.out_dir)
        .with_context(|| format!("creating output dir {}", args.out_dir.display()))?;

    cleanup_generated_pages(&args.out_dir)?;

    let mut docs = Vec::new();
    for command_name in command_names {
        let help_text = command_help_text(&root_command, &command_name)?;
        docs.push(CommandDoc {
            title: format!("cfdna {command_name}"),
            name: command_name,
            help_text,
        });
    }
    docs.sort_by(|left, right| left.name.cmp(&right.name));

    write_generated_notice(&args.out_dir)?;
    write_overview_page(&args.out_dir, &docs)?;
    for command_doc in &docs {
        write_command_page(&args.out_dir, command_doc)?;
    }

    if args.fail_on_drift {
        enforce_git_clean(&args.out_dir)?;
    }

    println!(
        "Generated {} command page(s) in {}",
        docs.len(),
        args.out_dir.display()
    );
    Ok(())
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn release_command_names() -> Vec<String> {
    let mut names = Vec::new();
    #[cfg(feature = "cmd_gc_bias")]
    names.push("gc-bias".to_string());
    #[cfg(feature = "cmd_ref_gc_bias")]
    names.push("ref-gc-bias".to_string());
    #[cfg(feature = "cmd_coverage_weights")]
    names.push("coverage-weights".to_string());
    #[cfg(feature = "cmd_lengths")]
    names.push("lengths".to_string());
    #[cfg(feature = "cmd_fcoverage")]
    names.push("fcoverage".to_string());
    #[cfg(feature = "cmd_midpoints")]
    names.push("midpoints".to_string());
    #[cfg(feature = "cmd_bam_to_bam")]
    names.push("bam-to-bam".to_string());
    #[cfg(feature = "cmd_bam_to_frag")]
    names.push("bam-to-frag".to_string());
    #[cfg(feature = "cmd_frag_to_bam")]
    names.push("frag-to-bam".to_string());
    names
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn all_command_names(root_command: &clap::Command) -> Vec<String> {
    let mut names: Vec<String> = root_command
        .get_subcommands()
        .map(|subcommand| subcommand.get_name().to_string())
        .collect();
    names.sort();
    names
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn cleanup_generated_pages(out_dir: &Path) -> Result<()> {
    for entry_result in fs::read_dir(out_dir)
        .with_context(|| format!("reading output dir {}", out_dir.display()))?
    {
        let entry = entry_result?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let file_name = path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("");
        let is_generated_page = file_name.ends_with(".md");
        if is_generated_page {
            fs::remove_file(&path)
                .with_context(|| format!("removing stale generated file {}", path.display()))?;
        }
    }
    Ok(())
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn command_help_text(root_command: &clap::Command, command_name: &str) -> Result<String> {
    let mut command = root_command
        .find_subcommand(command_name)
        .cloned()
        .ok_or_else(|| anyhow::anyhow!("could not find command '{}'", command_name))?;

    let mut bytes = Vec::new();
    command
        .write_long_help(&mut bytes)
        .with_context(|| format!("rendering long help for {}", command_name))?;
    String::from_utf8(bytes).context("help text is not valid UTF-8")
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn write_generated_notice(out_dir: &Path) -> Result<()> {
    let notice_text = "AUTO-GENERATED DIRECTORY - DO NOT EDIT\nSource: cfdna Clap config and command tree\n\nRegenerate with:\n\ncargo run --bin gen_cli_docs --features cli,docs_gen,cmd_bam_to_bam,cmd_bam_to_frag,cmd_frag_to_bam,cmd_coverage_weights,cmd_fcoverage,cmd_gc_bias,cmd_lengths,cmd_midpoints,cmd_ref_gc_bias -- --out-dir website/docs/generated/cli --scope release\n";
    fs::write(out_dir.join("GENERATED_NOTICE.txt"), notice_text).with_context(|| {
        format!(
            "writing {}",
            out_dir.join("GENERATED_NOTICE.txt").display()
        )
    })?;
    Ok(())
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn write_overview_page(out_dir: &Path, docs: &[CommandDoc]) -> Result<()> {
    let mut body = String::new();
    body.push_str(GENERATED_MARKER);
    body.push('\n');
    body.push_str(GENERATED_SOURCE);
    body.push_str("\n\n# CLI Reference Overview\n\n");
    body.push_str("Auto-generated command reference pages.\n\n");
    for command_doc in docs {
        body.push_str(&format!(
            "- [`cfdna {0}`](./{1}.md)\n",
            command_doc.name, command_doc.name
        ));
    }
    fs::write(out_dir.join("overview.md"), body)
        .with_context(|| format!("writing {}", out_dir.join("overview.md").display()))?;
    Ok(())
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn write_command_page(out_dir: &Path, doc: &CommandDoc) -> Result<()> {
    let mut body = String::new();
    body.push_str(GENERATED_MARKER);
    body.push('\n');
    body.push_str(GENERATED_SOURCE);
    body.push_str("\n\n");
    body.push_str(&format!("# {}\n\n", doc.title));
    let fence = code_fence_for_content(doc.help_text.trim_end());
    body.push_str(&format!("{fence}text\n"));
    body.push_str(doc.help_text.trim_end());
    body.push('\n');
    body.push_str(&fence);
    body.push('\n');
    fs::write(out_dir.join(format!("{}.md", doc.name)), body)
        .with_context(|| format!("writing command page for {}", doc.name))?;
    Ok(())
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn code_fence_for_content(content: &str) -> String {
    let mut longest_run = 0usize;
    let mut current_run = 0usize;
    for character in content.chars() {
        if character == '`' {
            current_run += 1;
            if current_run > longest_run {
                longest_run = current_run;
            }
        } else {
            current_run = 0;
        }
    }

    let fence_len = std::cmp::max(3, longest_run + 1);
    "`".repeat(fence_len)
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn enforce_git_clean(out_dir: &Path) -> Result<()> {
    let status = Command::new("git")
        .arg("diff")
        .arg("--exit-code")
        .arg("--")
        .arg(out_dir)
        .status()
        .with_context(|| format!("running git diff for {}", out_dir.display()))?;
    if !status.success() {
        bail!(
            "generated docs drift detected in {}. regenerate docs and commit changes",
            out_dir.display()
        );
    }
    Ok(())
}
