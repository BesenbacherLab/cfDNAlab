use cfdnalab::commands::coverage_weights::config::CoverageWeightsConfig;
use cfdnalab::commands::fcoverage::config::FCoverageConfig;
use cfdnalab::commands::fragment_kmers::config::FragmentKmersConfig;
use cfdnalab::commands::gc_bias::config::GCConfig;
use cfdnalab::commands::lengths::config::LengthsConfig;
use cfdnalab::commands::prepare_windows::config::PrepareConfig;
use cfdnalab::commands::profile_groups::config::ProfileGroupsConfig;
use cfdnalab::commands::reference_gc::config::RefGCConfig;

#[cfg(feature = "cli")]
#[cfg_attr(feature = "cli", derive(clap::Parser))]
#[command(name = "cfdna", version)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[cfg(feature = "cli")]
#[cfg_attr(feature = "cli", derive(clap::Subcommand))]
enum Cmd {
    GCBias(GCConfig),
    ReferenceGC(RefGCConfig), // Extract reference GC counts
    CoverageWeights(CoverageWeightsConfig),
    Lengths(LengthsConfig),
    Fcoverage(FCoverageConfig),
    Profiles(ProfileGroupsConfig),
    FragmentKmers(FragmentKmersConfig),
    PrepWindows(PrepareConfig),
    // Ends(EndsConfig),
}

#[cfg(not(feature = "cli"))]
fn main() {
    // Library-only builds (no binary) — keep this minimal
    eprintln!("This binary requires --features cli");
    std::process::exit(1);
}

#[cfg(feature = "cli")]
fn main() {
    use clap::FromArgMatches;
    let cmd = pretty::build_cmd();

    // Parse using the sanitized command
    let matches = cmd.clone().get_matches();
    let cli = Cli::from_arg_matches(&matches).expect("parse");

    // Run selected subcommand and capture its Result (no `?` in main).
    let res: anyhow::Result<()> = match cli.cmd {
        Cmd::GCBias(cfg) => cfdnalab::commands::gc_bias::gc_bias::run(&cfg),
        Cmd::ReferenceGC(cfg) => cfdnalab::commands::reference_gc::reference_gc::run(&cfg),
        Cmd::CoverageWeights(cfg) => {
            cfdnalab::commands::coverage_weights::coverage_weights::run(&cfg)
        }
        Cmd::Lengths(cfg) => cfdnalab::commands::lengths::lengths::run(&cfg),
        Cmd::Fcoverage(cfg) => cfdnalab::commands::fcoverage::fcoverage::run(&cfg),
        Cmd::Profiles(cfg) => cfdnalab::commands::profile_groups::profile_groups::run(&cfg),
        Cmd::FragmentKmers(cfg) => cfdnalab::commands::fragment_kmers::fragment_kmers::run(&cfg),
        Cmd::PrepWindows(cfg) => {
            cfdnalab::commands::prepare_windows::prepare_windows::run(&cfg)
        } // Cmd::Ends(cfg) => cfdnalab::ends::run(cfg),
    };

    if let Err(e) = res {
        eprintln!("{:#}", e);
        std::process::exit(1);
    }

    std::process::exit(0);
}

#[cfg(all(feature = "cli"))]
mod pretty {
    use clap::CommandFactory;
    use clap::builder::styling::{AnsiColor, Style, Styles};

    /// Sanitize Markdown-ish help for terminals:
    /// - Treat ``` fences as block code (no inline styling inside)
    /// - Apply header styling (#, ##, ###, ####) outside fences
    /// - Apply inline **bold** and `code` elsewhere
    /// - Normalize arrows/quotes
    pub fn sanitize_cli_text(md: &str) -> String {
        let mut out = String::with_capacity(md.len());
        let mut in_block = false;

        // Style used for code blocks (distinct from inline code)
        let block = Style::new().dimmed();
        let block_on = format!("{block}");
        let block_off = format!("{block:#}");

        for line in md.lines() {
            let trimmed = line.trim_start();

            // Toggle code-block mode on lines that start with ``` (any language tag)
            if trimmed.starts_with("```") {
                in_block = !in_block;
                continue; // drop the fence line itself
            }

            // Normalize a few typography chars for terminals
            let line = line
                .replace('→', "->")
                .replace('’', "'")
                .replace('“', "\"")
                .replace('”', "\"");

            if in_block {
                // In a fenced block: don't parse inline markers; optionally style/indent
                out.push_str("  "); // simple indent
                out.push_str(&block_on);
                out.push_str(&line);
                out.push_str(&block_off);
            } else if let Some(hdr) = try_render_header(&line) {
                // Header line outside fences
                out.push_str(&hdr);
            } else {
                // Outside a block and not a header: apply inline styling (**bold**, `code`)
                out.push_str(&stylize_inline(&line));
            }
            out.push('\n');
        }

        if out.ends_with('\n') {
            out.pop();
        }
        out
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

    /// Render Markdown-like headers (#, ##, ###, ####) on non-fenced lines
    fn try_render_header(line: &str) -> Option<String> {
        let t = line.trim_start();
        let mut hashes = 0usize;
        for b in t.as_bytes() {
            if *b == b'#' {
                hashes += 1;
            } else {
                break;
            }
        }
        if hashes == 0 {
            return None;
        }
        let rest = &t[hashes..];
        if !rest.starts_with(' ') {
            return None;
        } // require a space after #'s
        let text = rest.trim_start();

        // Choose styles per header level
        let (sty, underline): (Style, bool) = match hashes {
            1 => (
                Style::new().fg_color(Some(AnsiColor::Yellow.into())).bold(),
                true,
            ), // H1
            2 => (Style::new().bold(), false),      // H2
            3 => (Style::new().underline(), false), // H3
            _ => (Style::new().dimmed(), false),    // H4+
        };

        let on = format!("{sty}");
        let off = format!("{sty:#}");

        let mut out = String::with_capacity(line.len() + 32);
        out.push_str(&on);
        out.push_str(text);
        out.push_str(&off);
        out.push('\n');

        if underline {
            let bar_len = text.chars().count().min(64);
            out.push_str(&"─".repeat(bar_len));
        }
        Some(out)
    }

    /// Sanitize help/long_help for a Command, its args, and all subcommands.
    /// **NOTE**: Takes and RETURNS ownership to avoid borrow/move errors with clap's builder API.
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
    pub fn build_cmd() -> clap::Command {
        let mut cmd = crate::Cli::command();
        let styles = Styles::styled()
            .header(AnsiColor::Yellow.on_default().bold())
            .usage(AnsiColor::Green.on_default().bold())
            .literal(AnsiColor::Blue.on_default().bold())
            .placeholder(AnsiColor::Cyan.on_default());
        cmd = cmd
            .help_template("{name} {version}\n{about}\n\n{usage-heading} {usage}\n\n{all-args}\n")
            .styles(styles);
        cmd = sanitize_command(cmd);

        // Prepend a signature line everywhere
        let sig = make_signature();
        cmd = add_signature(cmd, &sig);
        cmd
    }

    /// Build a styled first-line signature (logo or horizontal rule)
    fn make_signature() -> String {
        // Choose a style; italic isn’t universal, bold is safe
        let accent = Style::new().bold();

        // A simple horizontal rule + title
        let title = "cfDNAlab";
        let bar1 = "_".repeat(48); // or just "-".repeat(60) for pure ASCII
        let bar2 = "─".repeat(48); // or just "-".repeat(60) for pure ASCII

        // Start style, content, then reset
        format!("\n{accent}{bar1}\n\n  {title}\n\n{bar2}{accent:#}\n")
    }

    /// Apply the signature to a Command and all its subcommands.
    /// Uses before_help / before_long_help so it prints at the very top.
    fn add_signature(mut cmd: clap::Command, sig: &str) -> clap::Command {
        cmd = cmd
            .before_help(sig.to_string())
            .before_long_help(sig.to_string());

        // Recurse into subcommands
        let sub_names: Vec<String> = cmd
            .get_subcommands()
            .map(|sc| sc.get_name().to_string())
            .collect();

        for name in sub_names {
            cmd = cmd.mut_subcommand(&name, |sub| {
                add_signature(
                    sub.before_help(sig.to_string())
                        .before_long_help(sig.to_string()),
                    sig,
                )
            });
        }
        cmd
    }
}
