//! Maintainer-only CLI reference generator for the docs website.
//!
//! This lives under `examples/` instead of `src/bin/` on purpose. The generated
//! command reference must be buildable from the repository checkout for CI and
//! local website development, but it is not a public cfDNAlab command and should
//! not appear as an installed binary on crates.io. The package `include` list in
//! `Cargo.toml` excludes `examples/`, so this tool stays out of the published
//! crate tarball while `cargo run --example gen_cli_docs` remains available here.
//!
//! Run through `website/scripts/generate_cli_docs.sh`, which passes the release
//! command features and writes into `website/docs/generated/cli/`.

#[cfg(all(feature = "cli", feature = "docs_gen"))]
use anyhow::{Context, Result, bail};
#[cfg(all(feature = "cli", feature = "docs_gen"))]
use cfdnalab::build_docs_command;
#[cfg(all(feature = "cli", feature = "docs_gen"))]
use clap::Parser;
#[cfg(all(feature = "cli", feature = "docs_gen"))]
use std::fs;
#[cfg(all(feature = "cli", feature = "docs_gen"))]
use std::path::{Path, PathBuf};
#[cfg(all(feature = "cli", feature = "docs_gen"))]
use std::process::Command;

#[cfg(all(feature = "cli", feature = "docs_gen"))]
const GENERATED_MARKER: &str = "<!-- AUTO-GENERATED FILE - DO NOT EDIT -->";
#[cfg(all(feature = "cli", feature = "docs_gen"))]
const GENERATED_SOURCE: &str = "<!-- Source: cfdna Clap config and command tree -->";

#[cfg(all(feature = "cli", feature = "docs_gen"))]
#[derive(Debug, Clone, clap::ValueEnum)]
enum Scope {
    Release,
    All,
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
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

#[cfg(all(feature = "cli", feature = "docs_gen"))]
#[derive(Debug, Clone)]
struct CommandDoc {
    name: String,
    title: String,
    help_text: String,
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
#[derive(Debug, Clone)]
struct ParsedHelp {
    intro_markdown: String,
    usage: Option<String>,
    sections: Vec<ParsedSection>,
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
#[derive(Debug, Clone)]
struct ParsedSection {
    title: String,
    options: Vec<ParsedOption>,
    notes: Vec<String>,
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
#[derive(Debug, Clone)]
struct ParsedOption {
    signature: String,
    description_lines: Vec<String>,
}

#[cfg(not(all(feature = "cli", feature = "docs_gen")))]
fn main() {
    eprintln!("This example requires --features cli,docs_gen");
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
    #[cfg(feature = "cmd_fragment_count_weights")]
    names.push("fragment-count-weights".to_string());
    #[cfg(feature = "cmd_ends")]
    names.push("ends".to_string());
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
    let help_text = String::from_utf8(bytes).context("help text is not valid UTF-8")?;
    Ok(normalize_help_text_for_docs(&help_text))
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn normalize_help_text_for_docs(help_text: &str) -> String {
    // Clap expands `default_value_t` before rendering help text. For options like
    // `--n-threads`, that makes generated docs depend on the machine that ran the
    // generator, which breaks CI drift checks. Normalize those host-dependent
    // defaults to a stable docs-only label.
    let mut normalized_lines = Vec::with_capacity(help_text.lines().count());
    let mut in_auto_threads_option = false;

    for line in help_text.lines() {
        if is_option_signature_line(line) {
            in_auto_threads_option = line.contains("--n-threads");
            normalized_lines.push(line.to_string());
            continue;
        }

        let trimmed_line = line.trim_start();
        if in_auto_threads_option && is_default_value_line(trimmed_line) {
            let indentation = &line[..line.len() - trimmed_line.len()];
            normalized_lines.push(format!("{indentation}[default: auto]"));
            continue;
        }

        normalized_lines.push(line.to_string());
    }

    normalized_lines.join("\n")
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn is_default_value_line(line: &str) -> bool {
    line.starts_with("[default: ") && line.ends_with(']')
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn write_generated_notice(out_dir: &Path) -> Result<()> {
    let notice_text = "AUTO-GENERATED DIRECTORY - DO NOT EDIT\nSource: cfdna Clap config and command tree\n\nRegenerate with:\n\ncargo run --example gen_cli_docs --features cli,docs_gen,cmd_bam_to_bam,cmd_bam_to_frag,cmd_frag_to_bam,cmd_coverage_weights,cmd_fragment_count_weights,cmd_ends,cmd_fcoverage,cmd_gc_bias,cmd_lengths,cmd_midpoints,cmd_ref_gc_bias -- --out-dir website/docs/generated/cli --scope release\n";
    fs::write(out_dir.join("GENERATED_NOTICE.txt"), notice_text)
        .with_context(|| format!("writing {}", out_dir.join("GENERATED_NOTICE.txt").display()))?;
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
    let parsed_help = parse_help_text(&doc.help_text, &doc.name);

    let mut body = String::new();
    body.push_str(GENERATED_MARKER);
    body.push('\n');
    body.push_str(GENERATED_SOURCE);
    body.push_str("\n\n");
    body.push_str(&format!("# {}\n\n", doc.title));

    if !parsed_help.intro_markdown.is_empty() {
        body.push_str(&parsed_help.intro_markdown);
        body.push_str("\n\n");
    }

    if let Some(usage_line) = parsed_help.usage {
        body.push_str("\n<hr class=\"cli-usage-separator\" />\n\n");
        body.push_str("## Usage\n\n");
        body.push_str("`");
        body.push_str(&usage_line);
        body.push_str("`\n\n");
    }

    for section in parsed_help.sections {
        body.push_str(&format!("## {}\n\n", section.title));

        let has_notes = !section.notes.is_empty();
        let has_options = !section.options.is_empty();

        for note_line in &section.notes {
            body.push_str(&note_line);
            body.push('\n');
        }
        if has_notes && has_options {
            body.push('\n');
        }

        for option in section.options {
            body.push_str("- `");
            body.push_str(&option.signature);
            body.push_str("`\n");

            if !option.description_lines.is_empty() {
                body.push('\n');
                for description_line in option.description_lines {
                    if description_line.is_empty() {
                        body.push_str("  \n");
                    } else {
                        if is_nested_bullet_line(&description_line) {
                            body.push_str("    ");
                        } else {
                            body.push_str("  ");
                        }
                        body.push_str(&description_line);
                        body.push('\n');
                    }
                }
            }
            body.push('\n');
        }
    }

    fs::write(out_dir.join(format!("{}.md", doc.name)), body)
        .with_context(|| format!("writing command page for {}", doc.name))?;
    Ok(())
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn parse_help_text(help_text: &str, command_name: &str) -> ParsedHelp {
    let lines: Vec<&str> = help_text.lines().collect();
    let usage_line_index = lines
        .iter()
        .position(|line| line.trim_start().starts_with("Usage: "));

    let intro_markdown = usage_line_index
        .map(|index| lines[..index].join("\n"))
        .unwrap_or_else(|| help_text.to_string())
        .trim_end()
        .to_string();

    let mut usage = None;
    let mut start_index = 0usize;
    if let Some(index) = usage_line_index {
        let usage_suffix = lines[index]
            .trim_start()
            .strip_prefix("Usage:")
            .unwrap_or("")
            .trim();

        if !usage_suffix.is_empty() {
            usage = Some(format!("cfdna {}", usage_suffix));
        } else if !command_name.is_empty() {
            usage = Some(format!("cfdna {}", command_name));
        }
        start_index = index + 1;
    }

    let mut sections = Vec::new();
    let mut current_section: Option<ParsedSection> = None;
    let mut line_index = start_index;

    while line_index < lines.len() {
        let current_line = lines[line_index];
        let trimmed_line = current_line.trim();

        if trimmed_line.is_empty() {
            line_index += 1;
            continue;
        }

        if is_section_header_line(current_line) {
            if let Some(section) = current_section.take() {
                sections.push(section);
            }
            current_section = Some(ParsedSection {
                title: trimmed_line.trim_end_matches(':').to_string(),
                options: Vec::new(),
                notes: Vec::new(),
            });
            line_index += 1;
            continue;
        }

        if is_option_signature_line(current_line) {
            if current_section.is_none() {
                current_section = Some(ParsedSection {
                    title: "Options".to_string(),
                    options: Vec::new(),
                    notes: Vec::new(),
                });
            }

            let signature = trimmed_line.to_string();
            let mut description_lines = Vec::new();
            line_index += 1;

            while line_index < lines.len() {
                let next_line = lines[line_index];
                let trimmed_next = next_line.trim();

                if is_section_header_line(next_line) || is_option_signature_line(next_line) {
                    break;
                }

                if trimmed_next.is_empty() {
                    let has_content = description_lines
                        .last()
                        .map(|line: &String| !line.is_empty())
                        .unwrap_or(false);
                    if has_content {
                        description_lines.push(String::new());
                    }
                    line_index += 1;
                    continue;
                }

                description_lines.push(trimmed_next.to_string());
                line_index += 1;
            }

            if let Some(section) = current_section.as_mut() {
                section.options.push(ParsedOption {
                    signature,
                    description_lines,
                });
            }
            continue;
        }

        if current_section.is_none() {
            current_section = Some(ParsedSection {
                title: "Details".to_string(),
                options: Vec::new(),
                notes: Vec::new(),
            });
        }
        if let Some(section) = current_section.as_mut() {
            section.notes.push(trimmed_line.to_string());
        }
        line_index += 1;
    }

    if let Some(section) = current_section.take() {
        sections.push(section);
    }

    ParsedHelp {
        intro_markdown,
        usage,
        sections,
    }
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn is_section_header_line(line: &str) -> bool {
    let trimmed_line = line.trim();
    if trimmed_line.is_empty() || !trimmed_line.ends_with(':') {
        return false;
    }
    let has_no_leading_whitespace = line
        .chars()
        .next()
        .map(|character| !character.is_whitespace());
    has_no_leading_whitespace.unwrap_or(false)
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn is_option_signature_line(line: &str) -> bool {
    let trimmed_line = line.trim_start();
    if trimmed_line.is_empty() {
        return false;
    }
    let is_indented = line.len() > trimmed_line.len();
    if !is_indented || !trimmed_line.starts_with('-') {
        return false;
    }

    // Real clap option signature lines look like:
    //   -h, --help
    //   --some-option <VALUE>
    // while prose bullets look like:
    //   - something
    // Require a non-whitespace character immediately after the first dash.
    let mut characters = trimmed_line.chars();
    let _first_dash = characters.next();
    matches!(characters.next(), Some(character) if !character.is_whitespace())
}

#[cfg(all(feature = "cli", feature = "docs_gen"))]
fn is_nested_bullet_line(line: &str) -> bool {
    line.starts_with("- ") || line.starts_with("* ")
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
