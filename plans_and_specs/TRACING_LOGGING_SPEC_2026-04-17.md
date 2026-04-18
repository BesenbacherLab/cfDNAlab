# Tracing Logging Spec

## Purpose

Define a shared logging design for all main `cfdna` commands.

The design should keep the CLI pleasant for normal users, while giving us a
real diagnostics layer for debugging and future maintenance.

## Current Situation

The codebase currently mixes three output styles:

- `println!` for normal command lifecycle messages and statistics
- `eprintln!` for warnings and some failures
- `indicatif` progress bars on `stderr`

That is workable for a small CLI, but it starts to fray in this package because
we now have many long-running commands, shared internal runners, temp-file
reducers, tiled processing, and command nesting such as `coverage-weights`
calling internal `fcoverage`.

The main problems are:

- There is no shared sink policy across commands
- Normal lifecycle output is mostly ad hoc `println!`
- Warnings and errors are not routed through one consistent layer
- Internal commands cannot easily log rich context without printing directly
- The output looks less deliberate than it should

## Goals

- Keep the normal CLI narrative readable and visible during successful runs
- Keep warnings, errors, and progress on `stderr`
- Make logging shared across commands rather than copy-pasted
- Add structured diagnostics without forcing normal users to read noisy logs
- Support logging to `stdout`, quiet mode, and an explicit log file path
- Avoid duplicate banners or duplicate log setup in nested internal runs
- Preserve clear final statistics blocks and other deliberate CLI presentation

## Non-Goals

- We do not want per-read or per-fragment logging at normal levels
- We do not want progress bars implemented through `tracing`
- We do not want each command to invent its own logging flags
- We do not want library code to initialize global logging on its own
- We do not want command output data mixed with warnings and errors on `stderr`

## Core Decision

Adopt `tracing` as the shared diagnostics layer, but do not force every piece
of user-facing CLI presentation through `tracing`.

We should split responsibilities like this:

- `tracing`
  Handles lifecycle logging, warnings, errors, debug diagnostics, and spans
- Explicit CLI output helpers
  Handle the branded header, command separators, and structured final summary
  blocks
- `indicatif`
  Continues to own progress bars and spinners on `stderr`

This gives us a professional logging system without flattening the whole CLI
into generic log lines.

## Why `tracing`

`tracing` is a better fit than raw `println!` or even `log` + `env_logger`
because this package has:

- long-running command phases
- nested command-like internal calls
- tiled and chromosome-wise processing
- optional deeper diagnostics that should stay off by default
- shared code paths that benefit from spans and structured fields

The important point is not async support. The important point is that `tracing`
lets us express:

- `info!` for normal lifecycle
- `warn!` for suspicious but recoverable situations
- `error!` for top-level failures
- `debug!` and `trace!` for deeper investigation
- spans for coarse phases and nested work

## CLI Surface

Use one shared `--log` argument.

The user explicitly asked for a single logging argument with inline file-path
selection rather than several separate logging flags. The spec should follow
that requirement.

### Proposed shared args

Add a shared `LoggingArgs` struct in `src/commands/cli_common.rs` and flatten it
into each main command config.

Proposed shape:

```rust
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct LoggingArgs {
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "stdout",
            value_parser = parse_log_spec,
            help_heading = "Logging"
        )
    )]
    pub log: LogSpec,
}
```

### Behavior

- `--log stdout`
  Default mode. Normal lifecycle logs go to `stdout`
- `--log quiet`
  Suppress normal lifecycle logs and disable progress bars. Warnings and errors
  still go to `stderr`
- `--log file`
  Send normal lifecycle logs to an auto-generated log file. Warnings and
  errors still go to `stderr`, and are also duplicated to the file
- `--log file=/path/to/run.log`
  Send normal lifecycle logs to the file. Warnings and errors still go to
  `stderr`, and are also duplicated to the file

### Proposed value model

Represent the parsed logging choice explicitly rather than leaving it as a raw
string.

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LogSpec {
    Stdout,
    Quiet,
    File(Option<PathBuf>),
}
```

This keeps the CLI surface compact while still giving the implementation one
clear resolved type.

The CLI surface is intentionally small. If we later need deeper diagnostics,
that should be added as a separate follow-up decision rather than bloating the
initial `--log` grammar.

### Constraints

- `--log quiet` disables progress bars
- `file=...` requires a non-empty path
- invalid `--log` values should produce a clear clap error that shows the
  accepted grammar
- `--help` behavior stays unchanged and continues to be driven by clap

### Accepted grammar

The parser should accept exactly:

- `stdout`
- `quiet`
- `file`
- `file=<path>`

This is intentionally compact and shell-friendly. It avoids an extra path flag
and keeps the logging contract small enough to actually remember.

### Auto-generated file path

When the user passes `--log file`, the program should create a log file path
automatically.

Use this rule:

- if the command has an output directory, use `<output_dir>/logs/`
- otherwise, if the command has a primary output file, use its parent directory
  plus `logs/`
- otherwise, fall back to the current working directory plus `logs/`

The generated file name should be:

- `<command_name>_<YYYY-MM-DD_HH-MM-SS>_<suffix>.log`

Where:

- `command_name` is the top-level CLI subcommand name
- `suffix` is 8 random alphanumeric characters

Example:

- `logs/fcoverage_2026-04-17_14-23-51_a1b2c3d4.log`

The implementation should create the `logs/` directory if needed.

## Stream Policy

The package should adopt one shared output contract.

### `stdout`

Used for the normal success-path run narrative when `--log stdout` is active.

This includes:

- command header and separators
- phase messages such as "Loading blacklists"
- completion messages such as "Saved output to ..."
- normal final statistics blocks

### `stderr`

Used for:

- warnings
- errors
- progress bars and spinners

This keeps `stderr` useful for "things that need attention", while still
allowing successful runs to show a readable command narrative elsewhere.

### File sink

When a file mode is active:

- normal lifecycle logs go to the file
- warnings and errors still appear on `stderr`
- warnings and errors are also duplicated to the file
- progress bars stay on `stderr`

This gives the user a complete run log on disk without losing immediate warning
visibility in the terminal.

## Formatting Rules

The logging output should look deliberate, not like default debug spew.

### Normal log format

The normal lifecycle stream should use a compact formatter with:

- no timestamps by default on `stdout`
- clear level labels only when they add value
- stable command/module labels
- no source file paths or line numbers in normal mode

Possible style:

```text
[cfDNAlab] [fcoverage] Loading blacklists
[cfDNAlab] [fcoverage] Counting per tile
[cfDNAlab] [fcoverage] Saved output to /path/out.tsv.gz
```

### Warning and error format

Warnings and errors should stay easy to scan on `stderr`.

Possible style:

```text
warning[gc_tag]: suppressing further GC tag weight warnings
error[bam_to_bam]: failed to open input BAM
```

### File format

File logs should be more complete than terminal logs.

Include:

- timestamp
- level
- command target
- message

Do not include source file and line numbers by default. Those are useful for
local debugging but too noisy for normal user logs.

## What Should Stay Explicit Instead Of `tracing`

Not everything should become a `tracing` event.

These outputs should remain owned by explicit helpers:

- the branded CLI header/logo
- horizontal separator lines around command execution
- structured statistics reports such as `print_fragment_run_statistics`

These are presentation blocks, not just log events.

The important change is that they should no longer write directly with
`println!`. They should write through a shared sink that follows the chosen
`LogSpec`.

## Proposed Shared Modules

### `src/shared/logging.rs`

Owns:

- `LoggingArgs` resolution into runtime policy
- `tracing_subscriber` setup
- writer routing for `stdout`, quiet, and file modes
- non-blocking file appender guard if file logging is active
- compact terminal formatter and file formatter

Possible public types:

- `LoggingConfig`
- `LoggingRuntime`
- `LogSpec`

### `src/shared/cli_output.rs`

Owns:

- writing the branded header
- writing command separators
- writing structured summary blocks to the primary sink

This keeps deliberate presentation out of the raw tracing formatter.

### `src/shared/progress.rs`

Remains the progress abstraction.

It should gain a straightforward hook so the progress factory can honor logging
policy, for example:

- disabled when `LogSpec::Quiet`
- otherwise still gated on `stderr` being a terminal

## Subscriber Setup

Tracing should be initialized once in `src/bin/cfdna.rs`, after clap parsing and
before command dispatch.

That means the binary needs access to the selected command's `LoggingArgs`.

One clean way is:

- add `logging: LoggingArgs` to all main command configs
- add `impl Cmd { fn logging_args(&self) -> &LoggingArgs }`
- in `main`, resolve logging policy from `cli.cmd.logging_args()`
- initialize tracing once
- print the branded header through `cli_output`
- dispatch the selected command

The command modules themselves must never initialize a subscriber.

## How Commands Should Use `tracing`

### Top-level command structure

Each top-level `run()` should use one command span and a small number of
coarse-grained lifecycle logs.

Example shape:

```rust
pub fn run(cfg: &FCoverageConfig) -> Result<()> {
    let command_span = tracing::info_span!(
        "command",
        command = "fcoverage",
        output_prefix = %cfg.output_prefix
    );
    let _entered = command_span.enter();

    tracing::info!("Loading blacklists");
    tracing::info!("Counting per tile");
    tracing::info!("Merging temporary tile files");

    Ok(())
}
```

### Good `info!` usage

Use `info!` for:

- phase transitions
- important file reads and writes
- high-level command decisions
- final success-path milestones

Do not use `info!` inside tight loops.

### Good `warn!` usage

Use `warn!` for:

- suspicious GC tag values that get neutralized
- retained temp directories
- dropped optional inputs
- degraded fallback paths (we likely shouldn't have these?)
- unusual but non-fatal conditions

### Good `debug!` usage

Use `debug!` for:

- tile counts
- window counts
- chosen chromosome subsets
- reducer merge counts
- extra derived config details

This is the right place for details that are useful during troubleshooting but
not in every successful run.

### Good `trace!` usage

Use `trace!` sparingly for:

- per-tile or per-chrom substeps
- detailed matching or filtering decisions
- unusually deep instrumentation during refactors

Trace should never become the default path for understanding a successful run.

## What Should Not Be Logged Through `tracing`

These should stay out of the tracing event stream:

- progress bar redraws
- dense final statistics tables
- large generated text payloads
- direct scientific output files

Those are better handled by dedicated renderers or file writers.

## Internal Command Reuse

Internal command reuse is where ad hoc printing gets ugly fast.

For example, `coverage-weights` should be able to call internal `fcoverage`
logic without:

- reinitializing logging
- reprinting the top-level command banner
- pretending it is a separate user-invoked command

The rule should be:

- only the binary prints the top-level header and command separator
- only the binary initializes tracing
- nested internal work uses spans and events within the same subscriber

That keeps one coherent run log.

## Migration Rules For Existing Output

We should migrate existing output by category.

### Replace with `info!`

Most `println!("Start: ...")` lines should become `info!`.

Examples include:

- loading blacklists
- loading windows
- counting per tile
- reducing temp files
- plotting outputs
- writing final files

### Replace with `warn!`

Most warning-like `eprintln!` calls should become `warn!`.

Examples include:

- temp directory retention notices
- GC tag warning suppression notices
- duplicate rule removal warnings
- "skipping because unavailable" style notices

### Keep as explicit structured output

These should stay as dedicated output helpers, but routed through the shared
primary sink:

- command header and separators in `src/bin/cfdna.rs`
- statistics blocks in `src/commands/run_statistics.rs`
- any future structured report blocks

### Keep as plain error propagation

Most failures should still return `anyhow::Result`.

At the top level, `main` remains responsible for the final user-visible error
line on `stderr`.

## Progress Bar Rules

Progress should remain separate from tracing.

Rules:

- progress bars always target `stderr`
- progress bars are disabled in `LogSpec::Quiet`
- progress bars are hidden when `stderr` is not a terminal
- warnings and errors must be printed in a progress-safe way so they do not
  shred the active progress display

This means logging setup and progress setup must share one runtime policy.

## Testing Plan

We should add integration tests for sink routing.

### Terminal routing tests

- default mode logs normal lifecycle to `stdout`
- default mode sends warnings to `stderr`
- default mode keeps errors on `stderr`
- quiet mode suppresses normal lifecycle output
- quiet mode still shows warnings and errors

### File logging tests

- `--log file=...` creates the file
- info logs go to the file, not `stdout`
- warnings still appear on `stderr`
- warnings are duplicated to the file

### Progress tests

- quiet mode disables progress
- non-TTY mode hides progress without affecting warnings and errors

## Suggested Rollout

### Phase 1

Add shared logging policy and tracing initialization.

Do not migrate every command immediately. First establish:

- `LoggingArgs`
- subscriber setup
- primary sink routing
- progress policy hook
- header and separator routing

### Phase 2

Migrate common lifecycle prints:

- `println!("Start: ...")` -> `info!`
- warning `eprintln!` -> `warn!`

Keep the existing summary blocks explicit.

### Phase 3

Add targeted `debug!` instrumentation in the commands that benefit most.

These deeper diagnostics do not need a public CLI flag in the first version.
They can remain available for development or future expansion without cluttering
the user-facing logging API.

Focus on:

- tiled counters
- reducers
- GC weighting paths
- blacklist loading and filtering setup

### Phase 4

Clean up any remaining direct `println!` and `eprintln!` calls that no longer
belong.

## Final Recommendation

Use `tracing` as the shared logging layer for command lifecycle and diagnostics.

Do not replace the whole CLI with generic log lines. Keep progress bars explicit
on `stderr`, and keep the header plus structured summary blocks as deliberate
CLI presentation routed through a shared sink.

That gives us:

- a more professional CLI
- shared output policy across commands
- clean separation between success-path narrative and warnings/errors
- a real path to deeper diagnostics when runs get weird
