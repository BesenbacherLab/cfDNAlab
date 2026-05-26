# Library API refactor plan

This plan sketches a Rust library API for running command workflows programmatically without exposing command internals like `run_inner`.

The goal is one real runner per command. The CLI and library API should call the same runner, with different run options.

## Goals

- Keep `run_inner` private or remove it where it only exists to separate execution from printing.
- Let downstream Rust packages run CLI-equivalent workflows through stable entry points.
- Return command-specific result objects instead of `Result<()>`.
- Avoid forcing all commands into one shared counter struct.
- Keep progress bars, status logging, and statistics printing as explicit runner options.
- Make public runner names command-specific, so call sites are readable.
- Avoid extra wrapper functions. The named runner is the command runner.

## Non-goals

- Do not expose low-level reducers, parsers, tiling helpers, temp-file orchestration, or private counters just because current tests import them.
- Do not expose the internal `crate::commands` module as a public API. Re-export only curated command config types, named runners, run options, and result types through the library API path.
- Do not preserve backwards compatibility with the current internal runner names.

## Shared run options

Use one shared options type for all command workflows. These controls are execution/reporting concerns, not command-specific semantics.

API shape:

```rust
pub struct RunOptions {
    pub report_statistics: bool,
    pub show_progress: bool,
    pub log_statuses: bool,
}

impl RunOptions {
    pub fn new_cli() -> Self {
        Self {
            report_statistics: true,
            show_progress: true,
            log_statuses: true,
        }
    }

    pub fn new_quiet() -> Self {
        Self {
            report_statistics: false,
            show_progress: false,
            log_statuses: false,
        }
    }
}
```

The field name is `log_statuses`, not `emit_statuses` or `emit_status_logs`.

The three booleans should stay independent because current command composition needs combinations that are not captured well by a binary reporting mode:

- CLI runs usually want all three enabled.
- Internal command composition usually wants statistics and progress disabled.
- Some future uses may want status logs without progress bars, or returned results without printed statistics.

## Runner naming

Do not expose every command as a generic `run()` function. Use command-specific runner names:

- `run_bam_to_bam`
- `run_bam_to_frag`
- `run_coverage_weights`
- `run_ends`
- `run_fcoverage`
- `run_frag_to_bam`
- `run_fragment_count_weights`
- `run_fragment_kmers`
- `run_gc_bias`
- `run_lengths`
- `run_midpoints`
- `run_prepare_windows`
- `run_ref_gc_bias`
- `run_transitions`
- `run_visualize_positions`
- `run_wps`
- `run_wps_peaks`

The CLI should call these functions with `RunOptions::new_cli()`.

Programmatic callers should call the same functions with `RunOptions::new_quiet()` or explicit options.

Example:

```rust
let result = cfdnalab::run_like_cli::midpoints::run_midpoints(
    &config,
    RunOptions::new_quiet(),
)?;
```

`run_like_cli` is still a reasonable namespace if these functions preserve CLI semantics around config interpretation and output files. The function names should make the command explicit.

## Result objects

Each command should return a command-specific result object:

```rust
pub struct MidpointsRunResult {
    pub counters: ProfileGroupsCounters,
    pub output_files: Vec<PathBuf>,
    pub primary_output: Option<PathBuf>,
}
```

The exact fields should be command-specific. Examples:

- `BamToBamRunResult`: output BAM path; counters.
- `BamToFragRunResult`: output fragment path; header path if written; counters.
- `FCoverageRunResult`: final output path; mean normalization length; counters.
- `CoverageWeightsRunResult`: scaling factors path; source `FCoverageRunResult` or selected source summary if needed.
- `FragmentKmersRunResult`: count matrix paths; motif path(s); counters.
- `TransitionsRunResult`: transition output paths; fragment-kmers source result or counters needed for statistics.

Do not force these into one shared result struct. The outputs are genuinely different.

## Shared result interface

Command-specific result objects should implement a shared trait so generic library code can inspect common parts without erasing command-specific details.

Proposed shape:

```rust
pub trait CommandRunResult {
    type Counters;

    fn counters(&self) -> &Self::Counters;
    fn output_files(&self) -> &[PathBuf];
    fn primary_output(&self) -> Option<&Path>;
}
```

Keep counters generic through the associated type. Do not introduce one public `RunCounters` struct unless there is a later, concrete need for it. Current command counters differ enough that a shared counter model would either be leaky or would freeze the wrong abstraction too early.

If generic statistics printing is useful internally, it can use an internal helper trait over the existing counter types. That does not need to be the public result interface.

## Refactor pattern

For each command, converge on this API shape:

```rust
pub fn run_midpoints(
    config: &MidpointsConfig,
    options: RunOptions,
) -> Result<MidpointsRunResult> {
    // validate config
    // run workflow
    // write outputs
    // collect result object
    // print statistics only if options.report_statistics
    // show progress only if options.show_progress
    // log status messages only if options.log_statuses
    // return result
}
```

CLI dispatch should call the named runner directly and ignore the returned result where the binary only needs success/failure:

```rust
run_midpoints(&config, RunOptions::new_cli()).map(|_| ())
```

The public re-export should expose the named runner, config, options, and result type, not `run_inner`.

## Current `run_inner` mapping

### Simple execution split cases

For `bam-to-bam` and `bam-to-frag`, `run_inner` currently does the real work and returns counters. The current command runner mostly measures elapsed time and prints final statistics around that call.

Refactor target:

- Move the real body into `run_bam_to_bam(config, options) -> BamToBamRunResult`.
- Move the real body into `run_bam_to_frag(config, options) -> BamToFragRunResult`.
- Print statistics inside the named runner only when `options.report_statistics` is true.
- Gate progress and status logging on `show_progress` and `log_statuses`.

Most integration tests that currently call `run_inner` for output files should migrate to the named runner with `RunOptions::new_quiet()`.

Tests that inspect counters directly should either:

- stay in integration tests if the result counters are intentionally exposed, or
- move module-local if they validate internal counting implementation rather than public behavior.

### Internal producer cases

`fcoverage` and `fragment-kmers` are also used as building blocks by other commands.

Refactor target:

- `coverage-weights` calls `run_fcoverage(&cfg, RunOptions::new_quiet())` and consumes the returned `FCoverageRunResult`.
- `fragment-count-weights` follows the same pattern through shared coverage-weights code.
- `transitions` calls `run_fragment_kmers(&cfg, RunOptions::new_quiet())` and uses the returned result/counters for its own outputs and statistics.
- `visualize-positions` calls `run_fragment_kmers(&cfg, RunOptions::new_quiet())`.

This removes `run_inner_silent` without losing the internal producer workflow.

### Private-only cases

`frag-to-bam` has a private `run_inner` used only by its own command runner.

Refactor target:

- Replace it with `run_frag_to_bam(config, options) -> FragToBamRunResult`.
- CLI dispatch calls the named runner with `RunOptions::new_cli()`.

## Public API path

The eventual public API should be curated. One possible shape:

```rust
pub mod run_like_cli {
    pub use crate::api::RunOptions;
    pub use crate::api::CommandRunResult;

    pub mod midpoints {
        pub use crate::commands::midpoints::config::MidpointsConfig;
        pub use crate::commands::midpoints::midpoints::{
            MidpointsRunResult,
            run_midpoints,
        };
    }
}
```

This keeps the command workflows separate from lower-level reusable modules like `interval`, `blacklist`, `gc_bias`, and `scale_genome`.

The exact public module name can still change. The important API invariant is:

- Config type
- Shared `RunOptions`
- Command-specific named runner
- Command-specific result object
- Shared result trait with an associated counter type

## Test migration impact

After this refactor, integration tests can be classified cleanly:

- Keep in `tests/` if they call public named runners, public lower-level APIs, or the binary CLI.
- Move module-local if they require private parsers, reducers, tiling helpers, state machines, or non-public counters.
- Rewrite current `run_inner` command-output tests to named runners before tagging them as keep-in-`tests/`.

This avoids pretending that a test is integration-ready while it still imports private internals.

## Migration order

1. Add `RunOptions`.
2. Add result structs and `CommandRunResult` for one command.
3. Refactor `bam-to-bam` first because its `run_inner` split is simple.
4. Refactor `bam-to-frag` next for the same reason.
5. Refactor `fcoverage`, then update `coverage-weights` and `fragment-count-weights` to call `run_fcoverage`.
6. Refactor `fragment-kmers`, then update `transitions` and `visualize-positions`.
7. Apply the pattern to the remaining commands.
8. Re-export the curated public API after the command runner shape is consistent.
9. Rewrite integration tests from `run_inner` to named runners and only then tag them as keep-in-`tests/`.
