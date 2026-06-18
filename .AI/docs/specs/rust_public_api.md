# Rust Public API

This spec describes the current Rust crate surface that downstream crates may use directly.

## Public Module Boundaries

- `crate::commands` and `crate::shared` remain internal implementation modules.
- Public Rust APIs are curated through re-export modules in `src/lib.rs`.
- Do not make command internals public just to satisfy tests. Move internal tests next to the owning module or add a deliberate public wrapper.
- Command-dependent public modules must be gated by the same cargo feature as the backing command module.

Stable public categories currently include:

- `cfdnalab::interval`
- `cfdnalab::blacklist`
- `cfdnalab::scale_genome`
- `cfdnalab::bam`
- `cfdnalab::fragment`
- `cfdnalab::reference`
- `cfdnalab::overlaps`
- `cfdnalab::positioning`
- `cfdnalab::indel_mode`
- `cfdnalab::clip_mode`
- `cfdnalab::constants`
- `cfdnalab::parsing`
- `cfdnalab::output_loaders`

`cfdnalab::gc_bias` is public only when `cmd_gc_bias` is enabled.

## Output Loaders

`cfdnalab::output_loaders` exposes Rust loaders for files written by cfDNAlab
commands. Each loader is gated by the cargo feature for the command that writes
that output:

- `cmd_lengths` exposes `load_lengths_output`.
- `cmd_fcoverage` exposes `load_fcoverage_output` and
  `load_fcoverage_output_with_group_index`.
- `cmd_ends` exposes `load_ends_output`.
- `cmd_midpoints` exposes `load_midpoints_output`.

Loader methods return `OutputLoaderResult<T>` with `OutputLoaderError`. The
error messages carry path, line, selector, and Zarr-array context. The public
error type is stable, but downstream code should not rely on exact message text
unless a behavior test in this repo pins it.

Loaders are strict readers for cfDNAlab outputs. They reject malformed schemas,
duplicate grouped labels or duplicate grouped fcoverage rows, impossible support
metadata, non-finite count values where the writer contract does not allow
them, and explicit row selectors on global outputs.

The detailed loader contract is in
[`output_loaders_spec.md`](output_loaders_spec.md).

## Command-Style Runners

Command workflows are exposed through `cfdnalab::run_like_cli`.

Each command module exports:

- The command config type.
- A command-specific named runner, such as `run_midpoints` or `run_ends`.
- A command-specific result type.

Runner names are command-specific. Do not expose generic `run()` functions for command workflows.

The CLI and the Rust API call the same command runners. CLI dispatch passes `RunOptions::new_cli()`. Programmatic callers should pass `RunOptions::new_quiet()` when they want the same files and computation without progress bars, status logging, or printed statistics.

`RunOptions` controls reporting side effects only:

- `report_statistics`
- `show_progress`
- `log_statuses`

These flags must not change scientific computation or output-file content.

## Command Results

Every command-style runner returns a command-specific result object.

Result objects implement `CommandRunResult`:

- `type Counters`
- `counters()`
- `output_files()`
- `primary_output()`

Keep result structs command-specific. Do not force all commands into one shared output or counter struct unless a concrete downstream use case proves that abstraction.

## Testing API

`cfdnalab::testing` is available under `#[cfg(any(feature = "testing", test))]`.

The `testing` cargo feature is opt-in and must not be enabled by default. Downstream crates opt in explicitly:

```toml
cfdnalab = { version = "...", features = ["testing"] }
```

The testing module is for tests and examples, not for production analyses. It provides small synthetic input builders and public-output readers.

Current testing submodules are:

- `testing::bam`
- `testing::bed`
- `testing::gc_packages`
- `testing::output_readers`
- `testing::reference`
- `testing::scaling`

Temporary fixture values own their temporary directories. Paths returned by those helpers remain valid while the owning value is alive and are removed when it is dropped.

Testing helpers should validate fixture inputs and fail clearly. In particular, fragment fixtures must respect the project-wide minimum possible fragment length of 10 bp.

Some testing helpers require command features in addition to `testing`. Keep these feature gates on the helper or submodule that needs them.

Integration tests that import `cfdnalab::testing` are external users of the crate. Compile or run them with `--features testing` or `--all-features`.

## Test Placement

- Keep tests in `tests/` when they use the public Rust API or the binary CLI.
- Move tests next to the owning module when they need private parsers, reducers, tiling helpers, state machines, or internal counters.
- Public testing fixtures may be used by both integration tests and module-local tests, but they must not leak production-private types only for assertion convenience.
