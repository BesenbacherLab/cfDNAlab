# Rust Output Loader Design

Status: future design, not an accepted implementation plan.

CLI command rendering was split out and implemented separately. The remaining
future work in this note is exported Rust loaders for files produced by
cfDNAlab command runners.

## Design Target

Downstream Rust tools should be able to run cfDNAlab commands, take the returned
output paths, and load those outputs without reimplementing TSV, zstd, Zarr, or
schema parsing. The loaders should make the output easy to work with in Rust,
not copy the R or Python package APIs.

That means:

- no dataframe dependency such as Polars or Arrow
- one public loader function per command output family
- typed structs, enums, native vectors, slices, and iterators
- explicit methods when data may be large
- schema validation before data access
- standalone path loaders that also work for files produced by previous runs
- convenience methods for command results only after the standalone loaders
  exist

## Non-Goals

- Do not expose private command internals as the loader API.
- Do not make the Rust crate a general dataframe layer.
- Do not eagerly load large positional tracks by default.
- Do not hide storage differences that matter for memory, such as dense versus
  sparse end-motif counts.
- Do not add backwards-compatibility adapters unless a release explicitly needs
  them.
- Do not treat plot/image outputs as loadable scientific data. Return paths and
  metadata for those.
- Do not add loaders for converter utility outputs such as frag files or BAM
  files in the first pass. Those can stay path-oriented unless a concrete Rust
  caller needs more.

## Public Module Shape

Add a curated public module:

```rust
pub mod output_loaders;
```

The public surface should be a facade, not a public mirror of the internal file
layout. The command-specific loader functions and their returned types are
public; parser helpers and storage-specific machinery stay `pub(crate)`.

Suggested internal layout:

```text
output_loaders/
  mod.rs
  common.rs
  compressed_text.rs
  zarr_store.rs
  lengths.rs
  midpoints.rs
  ends.rs
  fcoverage.rs
  scaling_weights.rs
  gc_bias.rs
```

`mod.rs` should re-export the public command loader functions and public loaded
types. The helper modules can stay private or `pub(crate)` unless a downstream
Rust caller has a direct reason to use them.

Command-dependent loader modules should be gated by the same cargo features as
the commands that produce the files:

```rust
#[cfg(feature = "cmd_lengths")]
mod lengths;

#[cfg(feature = "cmd_midpoints")]
mod midpoints;

#[cfg(feature = "cmd_lengths")]
pub use lengths::{load_lengths_output, LengthsOutput};

#[cfg(feature = "cmd_midpoints")]
pub use midpoints::{load_midpoints_output, MidpointsOutput};
```

Keep `run_like_cli` for running commands and `output_loaders` for reading
command outputs. Loader modules may reuse internal low-level parsers through
deliberate public wrappers, but the command modules themselves should stay
private.

## Public Loader Shape

Prefer one public loader function per command output family. If the command has
multiple output modes, return a command-specific enum whose variants contain
mode-specific structs:

```rust
pub fn load_lengths_output(
    path: impl AsRef<std::path::Path>,
) -> anyhow::Result<LengthsOutput>;

pub enum LengthsOutput {
    Global(GlobalLengthCounts),
    Windows(WindowLengthCounts),
    Groups(GroupLengthCounts),
}
```

The enum should provide common methods when they are truly shared, while the
variant structs own mode-specific accessors:

```rust
impl LengthsOutput {
    pub fn mode(&self) -> LengthOutputMode;
    pub fn length_bins(&self) -> &[LengthBin];
    pub fn row_count(&self) -> usize;
    pub fn count(&self, row_index: usize, length_bin_index: usize) -> Option<f64>;
}
```

Variant structs should expose mode-specific metadata without forcing every mode
through one generic row model:

```rust
impl WindowLengthCounts {
    pub fn windows(&self) -> &[WindowRow];
}

impl GroupLengthCounts {
    pub fn groups(&self) -> &[GroupRow];
}
```

This follows the useful part of the R and Python loaders: callers get an object
appropriate to the output mode. It should not copy their dataframe-oriented
return types.

The practical rule is:

- use methods on the output enum for common questions such as mode, length bins,
  row count, and direct count lookup
- use `match` or `try_as_*` accessors when the caller needs mode-specific
  metadata
- keep the full data on the variant structs so their APIs stay small and honest
  about the output mode

## Core API Conventions

Rust public loaders should use:

- zero-based indexing
- half-open genomic intervals
- explicit `Result<T>` errors with path and schema context
- existing crate vocabulary: fragment length, pos, end, and reference_end
- native `Vec<T>` storage and `&[T]` accessors for ordinary loaded data
- streaming iterators for large row-oriented outputs
- method names that make full materialization visible, such as
  `read_all_counts` or `read_dense_counts`

The public API should accept normal output paths. Callers should not need to
know whether an output is zstd-compressed or which Zarr metadata file contains a
schema field.

Dense multidimensional data should use small Rust-owned containers instead of
returning `ndarray` as the primary API:

```rust
pub struct DenseMatrix<T> {
    values: Vec<T>,
    rows: usize,
    columns: usize,
}

impl<T> DenseMatrix<T> {
    pub fn values_row_major(&self) -> &[T];
    pub fn shape(&self) -> (usize, usize);
    pub fn get(&self, row_index: usize, column_index: usize) -> Option<&T>;
    pub fn into_values_row_major(self) -> Vec<T>;
}
```

Use the same pattern for tensors, such as `DenseTensor3<T>`. This keeps the
default Rust API easy to use with standard collections while preserving shape
and indexing semantics. Optional `ndarray` view helpers can be added later if
there is a concrete caller need, but callers should not have to use `ndarray` to
load cfDNAlab outputs.

## Common Types

Start with small shared types that are stable and useful across output families:

```rust
pub struct WindowRow {
    pub index: usize,
    pub chrom: String,
    pub interval: cfdnalab::interval::Interval<u64>,
    pub blacklisted_fraction: Option<f64>,
}

pub type LengthBin = cfdnalab::interval::IndexedInterval<u32, usize>;

pub struct GroupRow {
    pub index: usize,
    pub name: String,
    pub eligible_items: Option<u64>,
    pub blacklisted_fraction: Option<f64>,
}
```

Do not add a `GenomicInterval` wrapper unless a later implementation needs one.
The crate already exposes checked half-open intervals as
`cfdnalab::interval::Interval`, so genomic rows should carry chromosome name and
interval separately. That keeps interval validation centralized and avoids
duplicating interval helpers.

Length bins can use the already exported `IndexedInterval` because they are
checked half-open intervals plus a stable output-column index.

Use output-specific accessor names for `eligible_items`, such as
`eligible_intervals()` for midpoint profiles and `eligible_windows()` for
grouped BED outputs. The shared struct should not force one biological label
onto all commands.

Selectors should exist only when they encode meaning beyond "all rows" versus
"these row indices". For simple row filtering, prefer ordinary optional index
slices:

```rust
pub fn rows(&self) -> &[WindowRow];
pub fn select_rows(&self, row_indices: Option<&[usize]>) -> Result<WindowLengthCountsView>;
```

Use `None` as the general "all rows" pattern. Do not add a `RowSelector::All`
enum only to wrap `None`.

Selectors are still useful where there are several distinct ways to select the
same axis:

```rust
pub enum LengthBinSelector {
    Indices(Vec<usize>),
    Lengths(Vec<u32>),
    Range(std::ops::Range<u32>),
}

pub enum GroupSelector {
    Indices(Vec<usize>),
    Names(Vec<String>),
}
```

Callers can pass `Option<LengthBinSelector>` or `Option<GroupSelector>` when
"all" is needed. Loaders should validate selector bounds and duplicate selector
entries where duplicates would make the returned data ambiguous.

## Length Counts

Provide one path loader for the command output:

```rust
pub fn load_lengths_output(
    path: impl AsRef<std::path::Path>,
) -> anyhow::Result<LengthsOutput>;
```

`LengthsOutput` should be mode-specific:

```rust
pub enum LengthsOutput {
    Global(GlobalLengthCounts),
    Windows(WindowLengthCounts),
    Groups(GroupLengthCounts),
}

pub struct GlobalLengthCounts {
    length_bins: Vec<LengthBin>,
    counts: DenseMatrix<f64>,
}
```

Expected methods:

- `length_bins(&self) -> &[LengthBin]`
- `counts(&self) -> &DenseMatrix<f64>`
- `counts_row_major(&self) -> &[f64]`
- `count(&self, row_index: usize, length_bin_index: usize) -> Option<f64>`
- `length_bin_for_length(&self, fragment_length_bp: u32) -> Option<usize>`
- `select_length_bins(&self, selector: Option<LengthBinSelector>) -> Result<LengthsOutputView>`

This output is small enough that eager loading is acceptable.

## Midpoint Profiles

Provide one opener for the command output, not an eager full-count loader:

```rust
pub fn load_midpoints_output(
    path: impl AsRef<std::path::Path>,
) -> anyhow::Result<MidpointsOutput>;
```

`load_midpoints_output` should validate the Zarr root schema immediately and
load lightweight metadata, but it should not read the full count tensor:

```rust
pub struct MidpointsOutput {
    path: std::path::PathBuf,
    groups: Vec<GroupRow>,
    length_bins: Vec<LengthBin>,
    position_bins: Vec<PositionBin>,
}
```

Expected methods:

- `groups(&self) -> &[GroupRow]`
- `length_bins(&self) -> &[LengthBin]`
- `position_bins(&self) -> &[PositionBin]`
- `read_all_counts(&self) -> Result<DenseTensor3<f32>>`
- `read_group_counts(&self, group_index: usize) -> Result<DenseMatrix<f32>>`
- `read_profile(&self, group_index: usize, length_bin_index: usize) -> Result<Vec<f32>>`

The name `read_all_counts` is intentionally explicit. A caller choosing it is
choosing to materialize the full `[group, length_bin, position]` tensor.
Document this clearly on `load_midpoints_output`: "load" means load the output
handle and metadata; `read_*` methods materialize count arrays.

## End-Motif Counts

Provide one opener for the command output:

```rust
pub fn load_ends_output(
    path: impl AsRef<std::path::Path>,
) -> anyhow::Result<EndsOutput>;
```

The loader should preserve the storage mode instead of always densifying. Users
should not need to know the storage mode before loading:

```rust
pub enum EndsOutput {
    Dense(DenseEndMotifCounts),
    SparseCoo(SparseEndMotifCounts),
}

pub enum EndMotifStorage {
    Dense,
    SparseCoo,
}

pub struct SparseEndMotifEntry {
    pub row_index: usize,
    pub motif_index: usize,
    pub count: f64,
}
```

Expected methods:

- `storage(&self) -> EndMotifStorage`
- `row_count(&self) -> usize`
- `motifs(&self) -> &[String]`
- `read_dense_counts(&self) -> Result<DenseMatrix<f64>>`
- `iter_sparse_counts(&self) -> Result<impl Iterator<Item = Result<SparseEndMotifEntry>>>`
- `read_sparse_counts(&self) -> Result<Vec<SparseEndMotifEntry>>`

For sparse stores, `read_dense_counts` should either return a clear error or be
named to make densification explicit, such as `read_dense_counts_allowing_densify`.
The first implementation should prefer the clear error unless a caller has a
concrete need for densification.

Callers can discover the storage mode either by matching `EndsOutput` or by
calling `storage()`:

```rust
match cfdnalab::output_loaders::load_ends_output(path)? {
    EndsOutput::Dense(output) => {
        let counts = output.counts();
    }
    EndsOutput::SparseCoo(output) => {
        for entry in output.iter_counts()? {
            // consume non-zero counts
        }
    }
}
```

## FCoverage Outputs

Provide one loader for the command output path:

```rust
pub fn load_fcoverage_output(
    path: impl AsRef<std::path::Path>,
) -> anyhow::Result<FCoverageOutput>;

pub enum FCoverageOutput {
    Average(FCoverageAverageTable),
    Total(FCoverageTotalTable),
    SummaryStats(FCoverageSummaryStatsTable),
}
```

Aggregate variants can load eagerly into typed rows. Each table should expose
rows as typed structs and provide `iter()` plus `rows()` accessors. Avoid a
generic "string column map" API unless the format is actually unstable.

Positional `fcoverage` outputs are out of scope for the first public Rust
loader API. `load_fcoverage_output` should reject `.bedgraph.zst` and
`per_position_per_window.tsv.zst` paths with a clear error instead of opening a
streaming reader or silently materializing genome-wide positional data. The first
pass should avoid bigWig or other indexing dependencies; callers that need
indexed positional queries can convert positional outputs outside the loader API.

## Scaling Weights

Wrap the existing scaling-factor parser rather than adding a second parser.

Expose two layers:

- an inspection reader that returns rows and metadata from the TSV alone
- a validated loader for applying factors that accepts the contig/chromosome
  context needed by the existing scaling-factor logic

Suggested public command loaders:

```rust
pub fn load_coverage_weights_output(
    path: impl AsRef<std::path::Path>,
) -> anyhow::Result<ScalingWeightRows>;

pub fn load_fragment_count_weights_output(
    path: impl AsRef<std::path::Path>,
) -> anyhow::Result<ScalingWeightRows>;
```

If callers need factors in the validated form used by counting code, expose a
separate explicit method or conversion that takes the required contig context:

```rust
pub fn into_loaded_scaling_factors(
    rows: ScalingWeightRows,
    chromosomes: &[String],
    contigs: &ContigCatalog,
) -> anyhow::Result<LoadedScalingFactors>;
```

The exact contig type should follow the current public reference/contig API
rather than inventing another representation.

## GC Bias Packages

GC packages already have production readers. The public loader work should add
stable wrappers around those package structs instead of duplicating Zarr logic.

Suggested names:

```rust
pub fn load_ref_gc_bias_output(
    path: impl AsRef<std::path::Path>,
) -> anyhow::Result<ReferenceGcPackage>;

pub fn load_gc_bias_output(
    path: impl AsRef<std::path::Path>,
) -> anyhow::Result<GcCorrectionPackage>;
```

These loaders should validate package role, schema version, coordinate
footprint, bin edges, and support masks through the same logic used by correction
code.

## Relationship To Command Results

Standalone path loaders are the primary API:

```rust
let result = cfdnalab::run_like_cli::lengths::run_lengths(&config, options)?;
let loaded = cfdnalab::output_loaders::load_lengths_output(
    result.primary_output().expect("lengths has a primary output"),
)?;
```

After path loaders exist, add convenience methods or extension traits for
command results:

```rust
pub trait LoadPrimaryOutput {
    type Loaded;

    fn load_primary_output(&self) -> anyhow::Result<Self::Loaded>;
}
```

This should be convenience only. It must not become the only supported route,
because workflow engines and previous command runs often hand callers a path
without a live `CommandRunResult`.

## Documentation

Document this library API in Rust docstrings on the exported
`cfdnalab::output_loaders` functions and returned public types. Those docstrings
should explain:

- what is loaded immediately
- which methods materialize large data
- how to match mode-specific output enums
- coordinate and indexing conventions
- dense versus sparse behavior

Do not add this to CLI user docs unless the CLI docs already discuss Rust
library usage. These loaders are a Rust API surface, not a CLI feature.

## Error Handling

Loaders should fail early and with context:

- missing path
- wrong file or directory shape
- unsupported schema name or version
- unsupported output family, such as a positional `fcoverage` bedGraph path
- missing required columns, arrays, attributes, or dimensions
- malformed UTF-8 where text labels are expected
- coordinate or selector values outside valid bounds
- dense-only API called on sparse output, or the reverse

Errors should include the path and the output family. For Zarr, include the
array or attribute name when possible.

## Implementation Order

1. Add `output_loaders::common`, `output_loaders::compressed_text`, and
   `output_loaders::zarr_store` as private or `pub(crate)` wrappers around
   existing internal helpers.
2. Implement `lengths` and scaling-weight loaders first. They validate the TSV
   reader shape and give downstream users immediate value.
3. Implement `midpoints` and `ends` with Zarr schema validation and typed
   metadata accessors before adding any convenience command-result methods.
4. Implement fcoverage aggregate table loaders. Positional fcoverage outputs
   stay out of scope for the first public Rust loader API.
5. Add GC package wrappers as needed by Rust callers.
6. Add command-result convenience methods once each output family has a stable
   standalone path loader.

## Testing

Public loader tests belong under `tests/` because downstream Rust crates are the
target caller.

Cover:

- tiny valid outputs for each supported output family
- zstd-compressed TSV handling
- Zarr schema name and version validation
- missing required TSV columns or Zarr arrays
- row metadata for global, window, and group modes
- length-bin and row selector bounds
- sparse end-motif iteration without silent densification
- command-result convenience methods, after those methods exist

Use small generated fixtures or existing public testing helpers. Do not require
full scientific pipeline tests just to prove loader parsing behavior.

## Acceptance Criteria

- A Rust caller can load the primary scientific output of each supported command
  through `cfdnalab::output_loaders`.
- Callers can inspect metadata and numeric values without parsing TSV, zstd, or
  Zarr directly.
- Large positional outputs and sparse stores are not silently materialized.
- Loader modules are feature-gated with their backing commands.
- The implemented public behavior is eventually distilled into
  `.AI/docs/specs/rust_public_api.md` and any schema-specific spec that applies.
