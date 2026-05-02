# `cfdna lengths` length-bin plan

Date: 2026-05-01

Scope: Replace the `lengths` command's `--min-fragment-length` and
`--max-fragment-length` output-axis configuration with `--length-bins`, while
preserving the current per-bp distribution as the default.

This is a planning document only. It describes the intended implementation and
the files likely to change.

## Goal

Make the length axis the single source of truth for `cfdna lengths`.

The default should be:

```text
--length-bins 30:1001:1
```

This preserves the current behavior because it creates one half-open bin per
integer fragment length:

```text
[30,31), [31,32), ..., [1000,1001)
```

Custom bins should reduce memory use linearly with the number of output columns:

```text
--length-bins 30:221:10
--length-bins 10 151 221
```

The effective fragment length filter should be derived from the resolved bin
edges:

```text
min_fragment_length = first_edge
max_fragment_length = last_edge - 1
```

Do not keep `--min-fragment-length`, `--max-fragment-length`, and
`--length-bins` as three independent knobs for `lengths`. That would create
unclear behavior around partial bin coverage, out-of-bin fragments, and which
argument owns the filtering range.

## Semantics

Length bins are half-open intervals:

```text
[edge_i, edge_i + 1)
```

For example:

```text
--length-bins 10 151 221
```

means:

```text
[10,151), [151,221)
```

The output matrix shape becomes:

```text
(# count vectors, # length bins)
```

The first dimension contains count vectors:

- One global vector in global mode.

- One vector per fixed-size or BED window in ordinary windowed modes.

- One vector per group in grouped BED mode.

Columns are length bins. In the default configuration, each length bin has width
1 bp, so the old per-position length distribution is still available without
special-case code.

Binned output remains fragment-count mass by length bin. It is not exact
base-weighted coverage. The old documentation that suggests multiplying counts
by length only remains exact for width-1 bins. For wider bins, exact coverage
requires either retaining the full per-bp spectrum or adding a separate
length-weighted output metric later.

## Shared parser refactor

Current state:

- `LengthBin`, `LengthBins`, and `parse_length_bins` live in
  `src/commands/cli_common.rs`.

- `parse_length_bins` handles compact specs such as `30:1001:1`.

- `MidpointsConfig::resolve_length_bins()` separately handles multiple explicit
  edge values such as `30 80 150 221`.

Refactor this into one shared resolver so `midpoints` and `lengths` cannot drift:

```rust
pub fn resolve_length_bin_edges(
    raw_values: &[String],
    min_allowed_length: u32,
    max_supported_fragment_length: u32,
) -> Result<Vec<u32>>
```

Rules:

- A single compact spec is accepted. Examples: `30:1001:1`, `30-80,80-150`
  only if the existing parser intentionally supports that form.

- The current code intentionally rejects start-end list specs like
  `30-80,80-150`. Preserve that rejection unless the format is deliberately
  reintroduced with tests and docs.

- Multiple integer values are interpreted as explicit edges.

- At least two edges are required.

- Edges must be strictly increasing.

- Every edge must be `>= min_allowed_length`.

- The final exclusive edge may be `MAX_SUPPORTED_FRAGMENT_LENGTH + 1`.

- No edge may exceed `MAX_SUPPORTED_FRAGMENT_LENGTH + 1`.

Then update `MidpointsConfig::resolve_length_bins()` to call the shared helper.

Do not create a shared clap argument struct for `midpoints` and `lengths` in
this pass. Their defaults are intentionally different:

- `lengths` should default to per-bp bins: `30:1001:1`.

- `midpoints` currently defaults to one broad bin: explicit edges `30 1001`.

Sharing the parser is an improvement. Sharing the CLI default would be a
behavior change and should not happen accidentally.

Also do not refactor `ProfileGroupsCounts` onto the new `LengthAxis` unless
there is a concrete benefit beyond reuse. Its existing public/artifact-facing
shape is already edge-based and tested. Parser unification is enough to remove
drift without disturbing midpoint counting internals.

## `LengthAxis`

Add a resolved length-axis type for `lengths`, probably in
`src/commands/lengths/counting.rs`:

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LengthAxis {
    edges: Vec<u32>,
    length_to_bin: Vec<usize>,
    single_bp_bins: bool,
}
```

Use the same LUT idea already used by `ProfileGroupsCounts` in
`src/commands/midpoints/counting_by_group.rs`:

- Allocate `length_to_bin` up to the final exclusive edge.

- Fill in-bounds lengths with the corresponding output column.

- Use `usize::MAX` as the invalid sentinel.

- Compute whether every bin is width 1 bp during construction and store it,
  rather than scanning the edges on every call.

Methods:

```rust
pub fn new(edges: Vec<u32>) -> Result<Self>
pub fn edges(&self) -> &[u32]
pub fn num_bins(&self) -> usize
pub fn min_fragment_length(&self) -> u32
pub fn max_fragment_length(&self) -> u32
pub fn contains(&self, length: u32) -> bool
pub fn bin_index(&self, length: usize) -> Option<usize>
pub fn is_single_bp_bins(&self) -> bool
```

Validation belongs in `LengthAxis::new()` as a second line of defense even if
the CLI resolver already validates edges.

The LUT must not be cloned per output row. `lengths` can have many thousands of
window rows, and a full supported-range LUT can have 100,001 entries. Store the
resolved axis behind `Arc<LengthAxis>` in per-row counters so cloning a
`LengthCounts` row only clones a pointer, not the lookup table.

## `LengthCounts`

Change `LengthCounts` from an integer-range counter to a length-axis counter.

Current shape:

```rust
pub struct LengthCounts {
    pub counts: Vec<f64>,
    pub length_min: usize,
    pub length_max: usize,
}
```

Planned shape:

```rust
use std::sync::Arc;

pub struct LengthCounts {
    pub counts: Vec<f64>,
    pub axis: Arc<LengthAxis>,
}
```

Keep the call-site shape stable, but do not silently drop an out-of-axis length
inside the command. The fragment filter should prevent this, so a counted
fragment with no output bin indicates an internal mismatch.

```rust
pub fn incr_weighted(&mut self, length: usize, weight: f64) -> Result<()>
```

Internally:

```rust
if let Some(bin_index) = self.axis.bin_index(length) {
    self.counts[bin_index] += weight;
} else {
    bail!("fragment length {length} did not map to any configured length bin");
}
```

This keeps the tile processing, window assignment, GC correction, scaling, and
reducers on one path. There should not be separate binned and non-binned
implementations of the `lengths` pipeline.

Compatibility checks should compare the full axis, not only vector width. This
prevents accidental merging of count vectors with the same number of columns
but different bin definitions.

If a non-failing lookup helper is still useful for tests or optional reads, name
it explicitly, for example `bin_index()` returning `Option<usize>`. The counting
path should use the checked increment.

`get(length)` becomes ambiguous with wider bins because `get(35)` would return
the same bin as `get(39)` for `[30,40)`. Either document this as "get by
absolute length's bin" or rename it to `get_by_length_bin_for_length()` before
the behavior becomes confusing.

## `LengthsConfig`

Replace the flattened `FragmentLengthArgs` field in
`src/commands/lengths/config.rs` with a command-specific `length_bins` argument:

```rust
/// Edges of fragment length bins to count in `[string(s)]`
#[cfg_attr(
    feature = "cli",
    clap(
        long,
        value_parser,
        num_args = 1..,
        default_values_t = [String::from("30:1001:1")],
        help_heading = "Core"
    )
)]
pub length_bins: Vec<String>,
```

Add:

```rust
pub fn set_length_bins(&mut self, edges: Vec<u32>)
pub fn set_length_bins_spec<S: Into<String>>(&mut self, spec: S)
pub fn set_per_bp_length_bins(&mut self, min_length: u32, max_length: u32)
pub fn resolve_length_bins(&self) -> Result<Vec<u32>>
```

Use the shared resolver for `resolve_length_bins()`.

`set_per_bp_length_bins()` is a programmatic helper only. It should not add
CLI min/max options back. It exists because the current `lengths` tests and
callers often mean "exact per-bp columns from min through max", and the helper
avoids duplicating `format!("{min}:{}:1", max + 1)` everywhere.

The command help should explain:

- Bins are half-open.

- The default preserves the old per-bp output.

- Memory use scales with the number of bins.

- The effective fragment length filter comes from the first and last edges.

- Wider bins do not preserve exact base-weighted coverage recoverability.

Also update the `--gc-length-range requested` help text. It currently describes
the requested range as `--min-fragment-length..--max-fragment-length`. After this
change it should describe the span covered by `--length-bins`.

## `lengths::run` and `process_tile`

Resolve the axis once near command startup:

```rust
let length_edges = opt.resolve_length_bins()?;
let length_axis = LengthAxis::new(length_edges)?;
```

Then replace the old `FragmentLengthArgs` uses with the resolved axis:

- GC correction loading uses `min_fragment_length()` and
  `max_fragment_length()`.

- Tile halo uses `max_fragment_length()`.

- Tile-window span right reach uses `max_fragment_length()`.

- `LengthCounts::new(...)` becomes `LengthCounts::new(length_axis.clone())`.

- BED fetch narrowing uses `max_fragment_length()`.

- Fragment filtering uses `length_axis.contains(adjusted_len)`.

- Overlap threshold epsilons use `max_fragment_length()`.

- Blacklist max-length reach uses `max_fragment_length()`.

- Scaling overlap max-length reach uses `max_fragment_length()`.

- Plotting uses the axis edges.

Avoid reading min/max back from `opt`. After resolution, the axis should own
the filtering and output range.

This does not fix L-001 from `major_reviews/04_lengths.md`. Length-bin max still
describes the maximum adjusted output length, while tiling and fetch halos still
need an aligned-coordinate reach. Do not present this refactor as solving
indel-adjusted fragments whose aligned span exceeds the output max. That remains
a separate aligned-span cap or fail-fast decision.

Indel and clip adjustment modes must not be treated as BAM-fetch geometry. Fetch
coordinates are aligned reference coordinates. Mode-specific behavior belongs in
fragment filtering, length binning, and window-assignment reach.

## Metadata and sidecars

The current `fragment_length_settings.json` is too thin and only writes min/max.
Replace it with a structured settings object. At minimum:

```json
{
  "length_axis": {
    "column_intervals": "half_open",
    "min_fragment_length": 30,
    "max_fragment_length": 1000,
    "n_bins": 971,
    "single_bp_bins": true,
    "bin_definition": {
      "kind": "stepped_range",
      "start": 30,
      "end": 1001,
      "step": 1
    }
  },
  "aggregation_level": "windows",
  "indel_mode": "ignore",
  "clip_mode": "aligned",
  "max_soft_clips": 256,
  "assign_by": "count-overlap",
  "gc_correction_used": false,
  "scaling_factors_used": false
}
```

Use a real serde struct rather than hand-writing JSON.

Do not rely on debug formatting for enum-like fields. Use stable lower-case
strings, like the `ends` sidecar helpers do for motif settings.

Include the fields already called out in L-003:

- `gc_length_weighting`

- `gc_length_range`

- Whether GC correction was used.

- Whether genomic scaling factors were used.

- Window mode and aggregation level.

- Assignment mode.

- Indel and clip settings.

Consider also writing a `length_bins.tsv` for ergonomic inspection:

```text
column_idx	start	end	label
0	30	31	30
1	31	32	31
```

For irregular bins, `bin_definition` should use `"kind": "explicit_edges"` and
store the full edge vector there. Axes that can be described by `start:end:step`
should use `"kind": "stepped_range"` metadata instead of dumping hundreds of
edges for the default per-bp axis.

This TSV is optional if the JSON contains the full bin definition, but it would
make the binned output easier to inspect without loading JSON.

## Plotting

The current overall line plot assumes every column is one integer length.

For single-bp bins, keep the existing x-axis behavior.

For wider bins, be explicit about what is plotted. Raw normalized bin counts are
mass per bin, not density per bp. If the existing line plot is retained, either:

- Divide by bin width and label the y-axis as density per bp.

- Or keep mass per bin and label the y-axis as relative counted mass per bin.

Do not keep the current `Density` label for unequal-width bins.

For wider bins, x-values can use bin midpoints:

```text
(start + end - 1) / 2
```

Change the x-axis label to `Fragment length bin midpoint (bp)` when bins are not
single-bp. Do not silently pretend that wide bins are exact integer lengths.

If a suitable plot cannot be made without misleading labels, skip the plot for
non-single-bp bins and log that decision.

## Additional pitfalls

The second code pass found these implementation traps:

- Do not clone a by-value LUT into every `LengthCounts` row. Use
  `Arc<LengthAxis>`.

- `parse_length_bins` currently only accepts colon ranges. The shared resolver
  must preserve explicit edge lists from `MidpointsConfig::resolve_length_bins()`
  and preserve the existing rejection of `30-80,80-150`.

- Do not let the `lengths` default leak into `midpoints`. Midpoints default must
  remain `[30,1001]`, which means one length bin.

- Avoid a broad "shared length-bin args" abstraction unless it can represent
  command-specific defaults without hiding them. A parsing helper is enough.

- Many integration tests construct `LengthsConfig` through
  `fragment_lengths_mut()`. Add a small test helper or config helper for exact
  single-bp ranges before migrating those tests.

- The CLI smoke test for `lengths` currently passes `--min-fragment-length` and
  `--max-fragment-length`. It must switch to `--length-bins 10:201:1` and assert
  the new settings structure.

- Cross-command artifact tests also call `LengthsConfig::fragment_lengths_mut()`.
  They should use `set_per_bp_length_bins()` or explicit edges.

- Temporary NPZ files only persist count matrices, not axis metadata. That is
  acceptable because the temp directory is per run and the reducer receives the
  template axis, but the reducer should continue to check matrix width.

- The final BED reordering still has the L-004 `zip()` truncation risk. Since the
  final write path will be touched for metadata anyway, fold in an explicit
  `bin_info.len() == all_bins.len()` check before zipping.

- `stack_length_counts()` should be hardened while touching counting output:
  accept a slice, return `Result<Array2<f64>>`, reject empty input, and reject
  rows with incompatible axes or widths.

## Tests

Do not run tests. Add focused tests and rely on `cargo check` for local
verification.

Shared parser tests:

- `resolve_length_bin_edges(["30:1001:1"])` returns dense edges.

- `resolve_length_bin_edges(["10", "151", "221"])` returns explicit edges.

- Reject non-increasing edges.

- Reject an edge below 10.

- Accept final edge `MAX_SUPPORTED_FRAGMENT_LENGTH + 1`.

- Reject final edge `MAX_SUPPORTED_FRAGMENT_LENGTH + 2`.

- Preserve rejection of `30-80,80-150` unless that syntax is intentionally
  reintroduced.

`LengthAxis` tests:

- `[30,40,50]` maps `30..39` to column 0 and `40..49` to column 1.

- Values below the first edge and at or above the final edge return `None`.

- `min_fragment_length()` and `max_fragment_length()` return first edge and last
  edge minus one.

- `is_single_bp_bins()` is true for `30:1001:1` and false for wider bins.

`LengthCounts` tests:

- `incr_weighted()` increments the expected column for dense bins.

- `incr_weighted()` increments the expected column for wider bins.

- `merge_from()` rejects count vectors with different axes.

- `incr_weighted()` returns an error rather than silently ignoring a length
  outside the configured axis.

- `get()` or its replacement has a test documenting wider-bin behavior.

Command tests:

- Default `lengths` output shape stays equivalent to the old default.

- `--length-bins 30:101:10` produces the expected number of columns and counts.

- Explicit short/long bins count fragments into expected columns.

- Settings JSON records compact `stepped_range` metadata for axes that can be
  described by `start:end:step`, and explicit edges for irregular axes.

- Grouped BED output records aggregation level as groups.

- CLI smoke uses `--length-bins 10:201:1` and no longer relies on removed
  min/max flags.

- Exact single-length ranges migrate cleanly, for example old min=max=10 becomes
  `[10,11)`.

- A non-single-bp bin that contains multiple observed lengths collapses both
  counts into one column. For example fragments of length 35 and 39 both count
  into `[30,40)`.

- A boundary fragment at the final exclusive edge is filtered out. For example
  `[10,20)` should include length 19 and exclude length 20.

- Plotting behavior for wider bins is covered enough to prevent the old density
  label from being reused incorrectly.

Midpoints regression tests:

- Existing `midpoints` length-bin behavior remains unchanged after moving the
  parser.

- Default `MidpointsConfig::resolve_length_bins()` still returns `[30,1001]`,
  not dense per-bp edges.

- `set_length_bins_spec("30:1001:1")` still works for users who explicitly ask
  midpoints for per-bp bins.

- The existing rejection of `30-80,80-150` remains covered unless that syntax is
  intentionally reintroduced.

## Expected files

Likely production files:

- `src/commands/cli_common.rs`

- `src/commands/midpoints/config.rs`

- `src/commands/lengths/config.rs`

- `src/commands/lengths/counting.rs`

- `src/commands/lengths/lengths.rs`

- `src/commands/lengths/tiling.rs`

- `src/commands/lengths/writer.rs`

Likely test files:

- `src/commands/cli_common_tests.rs`

- `src/commands/midpoints/config_tests.rs`

- `tests/test_lengths_command.rs`

- `tests/test_cli_smoke.rs`

- `tests/test_cross_command_artifact_matrix.rs`

## Verification

After implementation:

```text
cargo check --features cli,plotters
```

If test code is changed, the project instruction says to run:

```text
cargo check --tests --features cli,plotters
```

Do not run tests.
