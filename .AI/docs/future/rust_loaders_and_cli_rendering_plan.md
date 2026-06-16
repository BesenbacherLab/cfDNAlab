# Rust Loaders And CLI Rendering Plan

Status: future idea, not an accepted implementation plan.

This note captures two related Rust public API extensions:

- exported Rust loaders for files produced by command runners
- a way to render a filled command config as the equivalent `cfdna ...` CLI call

The shared goal is that downstream Rust tools can call cfDNAlab runners, inspect
the produced outputs, and show a reproducible CLI command without reimplementing
file parsers or command argument formatting.

## Goals

- Expose supported output readers through curated public modules, not by making
  command internals public.
- Let Rust callers load outputs returned by `CommandRunResult::output_files()`
  and `primary_output()`.
- Keep loaders typed by output family so callers do not need to parse headers,
  Zarr metadata, or compressed TSV schemas themselves.
- Let every filled config struct render an equivalent CLI command.
- Make rendered CLI calls suitable for logs, reports, notebooks, workflow
  provenance, and user-facing diagnostics. This could include splitting each argument
  into its own line with "  --arg value \".

## Non-Goals

- Do not turn the Rust crate into a dataframe library.
- Do not eagerly load large positional tracks by default.
- Do not expose private reducer, tiling, or writer internals only to support
  loaders.
- Do not promise to reconstruct the exact CLI text originally typed by a user.
  The useful target is an equivalent command for the filled config object.
- Do not add backwards-compatibility clutter unless this becomes a stated
  requirement for a release.

## Public Module Shape

Add a public module such as:

```rust
pub mod outputs;
```

Possible layout:

```text
outputs/
  mod.rs
  compression.rs
  lengths.rs
  ends.rs
  midpoints.rs
  fcoverage.rs
  scaling_weights.rs
  gc_bias.rs
  command_outputs.rs
```

Command-dependent loaders should be gated by the same cargo feature as the
command that writes the output. For example, `outputs::lengths` should require
`cmd_lengths`, and `outputs::fcoverage` should require `cmd_fcoverage`.

Keep `run_like_cli` for running commands. Keep `outputs` for reading command
outputs. Avoid putting output loading methods directly on every command module
unless that proves clearer in practice.

## Loader API Sketch

Prefer explicit command-specific reader names:

```rust
use cfdnalab::outputs::lengths::read_lengths_counts;
use cfdnalab::outputs::midpoints::open_midpoint_profiles;

let result = cfdnalab::run_like_cli::lengths::run_lengths(&config, RunOptions::new_quiet())?;
let lengths = read_lengths_counts(result.primary_output().unwrap())?;

let midpoint_result =
    cfdnalab::run_like_cli::midpoints::run_midpoints(&midpoint_config, RunOptions::new_quiet())?;
let profiles = open_midpoint_profiles(midpoint_result.primary_output().unwrap())?;
```

Also consider an extension trait for command results:

```rust
trait LoadPrimaryOutput {
    type Loaded;

    fn load_primary_output(&self) -> Result<Self::Loaded>;
}
```

This is convenient, but it should not replace standalone functions. Standalone
functions are still needed when users receive a path from a workflow engine or a
previous run.

## Loader Scope By Output Type

Start with cfDNAlab-specific outputs where the crate already owns the schema.

High priority:

- `lengths`: read `.length_counts.tsv.zst` into typed length-bin metadata,
  row metadata, and count arrays.
- `midpoints`: open `.midpoint_profiles.zarr` with schema validation and typed
  accessors for groups, length bins, positions, and counts.
- `ends`: open `.end_motifs.zarr` with schema validation, motif metadata, row
  metadata, and dense or sparse count accessors.
- `fcoverage` aggregate outputs: read `.fcoverage.average.tsv.zst`,
  `.fcoverage.total.tsv.zst`, and `.fcoverage.summary_stats.tsv.zst` variants.
- `coverage-weights` and `fragment-count-weights`: read scaling-factor TSVs
  using the existing scaling-factor parser where possible.
- `gc-bias` correction packages: expose stable package readers that match the
  public correction APIs.

Lower priority or adapter-only:

- BAM outputs from `bam-to-bam` and `frag-to-bam`: return paths or small helper
  adapters around `rust-htslib`, not a custom BAM abstraction.
- Fragment TSV outputs from `bam-to-frag`: add a typed streaming reader if the
  format is used by Rust callers.
- Positional `fcoverage` bedGraph outputs: expose streaming iteration only.
  Do not provide a convenient eager full-load helper.
- Plot or image outputs from visualization commands: return paths and metadata,
  not image-processing loaders.

## Loader Behavior

Each loader should:

- validate that the path exists and has the expected file or directory shape
- validate schema metadata before exposing data
- fail with command-specific error context when required columns or arrays are
  missing
- expose row mode explicitly, such as global, windowed, grouped, or positional
- use Rust zero-based indexing for public Rust selectors
- use half-open genomic intervals for public Rust ranges
- preserve cfDNAlab vocabulary, including fragment length, pos, end, and
  reference_end where those concepts appear
- avoid silently materializing huge dense arrays

The loader API should not require callers to know whether a file is compressed
with zstd. Public functions should accept the command output path and handle the
compression layer internally.

## Relationship To Python And R Loaders

The Rust loaders should share schema decisions with the Python and R packages,
but they do not need to imitate those APIs exactly.

Rust should prefer:

- typed structs and enums over dynamic objects
- iterators and borrowed views where practical
- `Result<T>` errors with source context
- zero-based selectors
- explicit materialization methods for large arrays

Python and R can remain data-frame-oriented. Rust loaders should provide enough
structured access for downstream crates to build their own data frames or
domain-specific objects without reparsing files.

## CLI Rendering API

Add a trait implemented by each filled command config:

```rust
pub trait ToCliCommand {
    fn to_cli_args(&self) -> Result<Vec<std::ffi::OsString>>;
    fn to_cli_string(&self) -> Result<String>;
}
```

`to_cli_args()` is the canonical API because paths and values are argument
tokens, not shell text. `to_cli_string()` should be a display helper that shell
quotes tokens for logs and reports.

The rendered command should include the binary and subcommand by default:

```text
cfdna lengths --bam sample.bam --output-dir out --output-prefix sample ...
```

If users need only the argument tail, add an explicit helper:

```rust
fn to_cli_argv_tail(&self) -> Result<Vec<OsString>>;
```

Do not overload `Display` for configs unless the output is clearly documented
as human display only. `Display` cannot return errors and is a poor fit for
path quoting failures or invalid config state.

## CLI Rendering Semantics

The method should render an equivalent CLI call for the current filled config.
It is not responsible for remembering how the config was originally built.

Rules:

- Render command arguments in deterministic order.
- Include all fields that affect scientific computation or output location.
- Include reporting-related flags only when the config owns them. Current
  `RunOptions` are not stored in command configs, so they need a separate
  rendering path if provenance should include quiet or progress behavior.
- Prefer long flags over short flags.
- Use one spelling per enum value and keep it synchronized with Clap parsing.
- Render paths as `OsString` tokens in `to_cli_args()`.
- Shell-quote only in `to_cli_string()`.
- Validate unsupported config states and error clearly instead of omitting
  fields.

Open decision:

- Decide whether rendered commands should include every resolved default or only
  non-default values plus required arguments.

The safer first implementation is to render a fully explicit command for fields
that affect computation. It is longer, but it avoids ambiguity when defaults
change before a command is re-run.

## Implementation Strategy

Do not hand-code argument formatting independently from Clap parsing if it can
be avoided. That creates two command vocabularies that can drift.

Preferred approach:

- Define a small command-argument rendering trait for leaf config structs and
  shared config structs.
- Reuse the same field names, enum value parsers, and validation helpers used
  by config construction.
- Add command-level implementations that assemble the subcommand and delegate to
  shared config renderers.
- Keep rendering code close to config definitions, because config structs own
  the mapping from Rust fields to CLI flags.

If repetitive implementations become too large, consider a small internal macro
or derive helper. Do not start with a procedural macro unless the manual pattern
has become a real maintenance problem.

## Testing

Add public API tests under `tests/` because downstream crates are the target
caller.

Loader tests should cover:

- real small outputs produced by command runners
- malformed or wrong-schema files
- compressed TSV reading
- Zarr schema validation
- selector bounds and duplicate selector errors
- no eager load path for large positional outputs

CLI rendering tests should cover:

- each command config renders a deterministic argument vector
- rendered arguments can be parsed by the existing CLI config path
- parsed config is equivalent to the original filled config for fields that
  affect computation and output paths
- paths with spaces are preserved as single tokens in `to_cli_args()`
- `to_cli_string()` shell-quotes paths and values safely enough for display

Do not require tests to run full scientific pipelines just to validate rendering.
Use config parse round-trips and small generated outputs where possible.

## Acceptance Criteria

- A Rust caller can run a command, load its primary output through public Rust
  APIs, and inspect the loaded object without parsing TSV, Zarr, or compressed
  files directly.
- Every public command config exported under `run_like_cli` has a documented way
  to render an equivalent CLI call.
- The rendered CLI arguments round-trip through the same config validation path
  used by the CLI.
- Loader modules and CLI rendering implementations are feature-gated with their
  backing commands.
- Current implemented behavior is eventually distilled into
  `.AI/docs/specs/rust_public_api.md` after the design is accepted and built.
