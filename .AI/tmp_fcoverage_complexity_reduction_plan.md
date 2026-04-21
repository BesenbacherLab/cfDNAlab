# fcoverage Complexity Reduction Plan

This document updates the earlier plan after comparing it against
`.AI/fcoverage_refactor_analysis.md`. The earlier version had the right
high-level direction, but it missed some of the biggest low-risk complexity
reducers and proposed one unnecessary intermediate data shape.

The target stays the same:

- reduce reducer/writer complexity without changing behavior
- keep temp files narrow in basic mode
- preserve row identity and cross-tile correctness
- keep summary-stat derivation outside the reducers

What changes in this revision is the refactor boundary and the recommended
order of work.

## Hard Constraints

These are non-negotiable. Any refactor that violates one of them is wrong.

- Keep the temp-file width optimization
  Basic partial rows must stay narrow
  Summary partial rows may carry extra additive fields
  Do not switch back to always writing `nonzero_positions` and
  `coverage_sum_of_squares`

- Keep persisted row identity unchanged
  BED reduction must still key by original BED `orig_idx`
  Fixed-size reduction must still key by full bin `start`
  Partial rows must still preserve full window/bin identity, not tile-local
  overlap identity

- Keep cross-index semantics unchanged
  Missing from all cross-index files means exactly one expected contribution
  Presence count across cross-index sidecars means exact expected contribution
  count

- Keep grouped behavior unchanged
  Group folding must still add already-reduced segment rows into group
  accumulators
  Group span must still come from the grouped layout, not from reduced segment
  lengths

- Keep summary-stat math unchanged
  Negative variance handling, CV near-zero handling, and `>1e6` display
  behavior must stay as-is

- Keep fast path behavior unchanged
  Aligned `--by-size` outputs must still use concatenation instead of cross-tile
  reduction

- Do not silently change public output schema unless explicitly intended
  Any header rename or formatting change should be treated as a separate
  decision

- Do not silently change reducer parse errors unless intended
  A shared parse helper will likely standardize the error text
  That is probably acceptable because these are internal temp-file parse errors
  But it should be treated as an explicit consequence, not an unnoticed one

## Invariant Summary

This is the invariant summary that should survive any implementation.

- A partial row represents one tile's contribution to one logical BED window or
  one logical fixed-size bin
- Cross-tile correctness comes from:
  - grouping by stable row identity
  - accumulating additive raw fields
  - waiting until `seen_contributions == expected_contributions`
- The merge heap only chooses which open stream to read from next
  It does not prove that a row is complete
- Basic mode and summary-stats mode differ in temp-file schema, not in row
  identity rules
- BED and fixed-size reduction differ in key type and interval recovery
  They do not differ in the high-level reduction algorithm
- BED partial rows do not carry the final interval on disk
  The BED merge engine must still recover the interval from `coords_by_idx`
- Size partial rows do carry interval coordinates
  The size merge engine must still clip the last bin to chromosome end after
  reduction

## Current Complexity Problem

The original plan was right that the reducer has too many copies of the same
merge algorithm. But that was not the whole problem.

Today there are effectively four copies of the same merge pattern:

- BED basic reducer
- BED summary reducer
- size basic reducer
- size summary reducer

And there are several additional duplication layers around them:

- four stream structs with nearly identical line-reading, blank-line-skipping,
  and open logic
- four parsed row structs or row-like parsing paths
- four accumulators that differ by only two fields
- four copies of cross-index loading
- large amounts of repeated per-column parsing boilerplate
- grouped writer fold logic that branches on basic vs summary even though the
  accumulation pattern is the same

The biggest single source of line count is not only the merge loop. It is the
mechanical parser boilerplate around the merge loop.

## What Was Right In The Earlier Plan

These points should be preserved.

- Keep separate on-disk schemas for basic vs summary partials
- Keep raw additive row reduction as a separate stage from final value
  derivation
- Keep one BED reduction path and one fixed-size reduction path as the right
  high-level simplification boundary
- Keep BED vs size identity semantics explicit
- Keep grouped outputs based on already-reduced raw rows
- Keep the aligned `--by-size` fast path out of the generic reducer

## What Needed Correction

### 1. The parser boilerplate was underestimated

The earlier plan focused on merge-loop duplication and said to keep separate
parsers. That misses one of the highest-value low-risk cleanups.

Across the four partial schemas there are 24 parsed columns. Each column parse
currently repeats the same missing-field and invalid-field scaffolding with only
the field name changed.

That should be reduced first with a small helper such as:

```rust
fn parse_col<T: FromStr>(
    cols: &mut std::str::Split<'_, char>,
    field_name: &str,
    chr: &str,
    tile_index: u32,
    line_number: u64,
) -> Result<T>
where
    T::Err: std::fmt::Display,
```

This is low risk and immediately cuts a large amount of noise.

### 2. `ParsedContribution` was the wrong extra abstraction

The earlier plan proposed a `ParsedContribution` struct that was almost the same
as `ReducedAggregateRow`.

That is not a useful layer. It adds one more type without solving a real
problem.

The cleaner split is:

- `ParsedPartialRow` for what the parser reads from disk
- `ReducedAggregateRow` for the fully reduced row emitted by the reducer

For BED, the parsed row can carry `interval: None` because the interval is
recovered later from `coords_by_idx`.

For size, the parsed row can carry `interval: Some(...)`.

### 3. The four stream structs should be collapsed

The earlier plan noted duplicated stream structs but did not turn that into an
actual design decision.

There is no strong reason to keep four stream structs when their IO behavior is
identical. A better target is:

```rust
enum PartialsSchema {
    BedBasic,
    BedSummary,
    SizeBasic,
    SizeSummary,
}

struct ParsedPartialRow {
    key: u64,
    interval: Option<Interval<u64>>,
    coverage_sum: f64,
    eligible_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
}

struct PartialsStream {
    reader: BufReader<Box<dyn Read + Send>>,
    line_buf: String,
    chr: String,
    line_number: u64,
    tile_index: u32,
    schema: PartialsSchema,
}
```

With this shape:

- `open()` is shared
- blank-line skipping is shared
- line-number tracking is shared
- column parsing is still explicit through a `match self.schema`
- basic-mode rows fill summary-only fields with zeros

This does not widen temp files. It only unifies the read side.

### 4. The writer grouped fold opportunity was missing

The earlier plan discussed shared accumulators in reducers, but it did not
explicitly connect that to `write_grouped_bed_aggregate_output`.

That function currently has separate chromosome-iteration loops for basic and
summary grouped folding. Once the accumulation shape is unified, those loops can
also be unified and only the final output derivation needs to branch.

This is not the primary payoff, but it is a real and safe secondary cleanup.

### 5. Some public reducers are already trivial wrappers

Two of the public reducer entry points are thin finalize-and-write wrappers
around the `_rows` variants.

That means the plan should not treat all public reducer functions as equally
important API layers. After reducer unification, the non-`_rows` wrappers should
either stay intentionally tiny or be inlined into their single writer call
sites.

This is optional, but the plan should acknowledge it.

## Recommended Target Shape

The best target is now:

- one shared `parse_col` helper
- one `ParsedPartialRow`
- one `PartialsStream` with `PartialsSchema`
- one shared reducer accumulator shape
- one BED merge engine
- one size merge engine
- grouped writer folding simplified to reuse the same additive field set

What should still stay separate:

- BED vs size merge engines
- basic vs summary on-disk schemas
- raw-row reduction vs final-value derivation
- aligned size fast path vs reduction path

## Internal Data Shapes

### Parsed row

Use one parsed row shape that reflects disk inputs without pretending BED rows
already know their interval:

```rust
struct ParsedPartialRow {
    key: u64,
    interval: Option<Interval<u64>>,
    coverage_sum: f64,
    eligible_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
}
```

Semantics:

- BED rows
  `key = orig_idx`
  `interval = None`
- size rows
  `key = bin_start`
  `interval = Some(full_bin_interval)`
- basic rows
  `nonzero_positions = 0`
  `coverage_sum_of_squares = 0.0`

### Reduced row

Keep the reduced output shape explicit:

```rust
struct ReducedAggregateRow {
    idx: u64,
    interval: Interval<u64>,
    coverage_sum: f64,
    eligible_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
}
```

### Shared reducer accumulator

Use one additive accumulator for the reducer merge engines:

```rust
#[derive(Default)]
struct AggregateAccum {
    coverage_sum: f64,
    eligible_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
    seen_contributions: u32,
}
```

This replaces the four reducer accumulators:

- `BasicWindowAccum`
- `SummaryWindowAccum`
- `BasicBinAccum`
- `SummaryBinAccum`

It does not replace the grouped writer accumulator. Grouped folding has a
different responsibility:

- it folds already-reduced rows, not open K-way streams
- it does not need `seen_contributions`
- it does need `span_positions` from the grouped layout

So `GroupedAggregateAccum` should stay separate. The writer-side cleanup still
works because the additive field set matches once basic-mode reduced rows carry
zero-valued summary-only fields.

Also note the naming normalization:

- current reducer accumulators use `allowed_positions`
- reduced rows and grouped folding use `eligible_positions`

This plan intentionally standardizes on `eligible_positions` for the unified
reducer accumulator so the field names match the reduced-row and writer-side
terminology.

## Recommended Refactor Order

The earlier plan started too late in the stack. The best order is to remove the
lowest-risk duplication first, then unify the merge engines.

### Step 1: Freeze behavior with focused regression tests

Before changing reducer structure, pin the behavior that must survive:

- BED basic reducer output
- BED summary raw-row reduction output
- size basic reducer output
- size summary raw-row reduction output
- grouped fold behavior
- aligned size fast path
- final-bin clipping for size reduction
- cross-tile merge completeness based on expected contribution count

If coverage is missing, add direct reducer tests first instead of relying only
on end-to-end writer behavior.

Important missing tests from the earlier version:

- round-trip reducer tests that write synthetic partials files and verify exact
  reduced rows
- explicit cross-tile merge tests for the same `orig_idx`
- basic vs summary equivalence tests for shared additive fields

### Step 2: Add `parse_col`

This is the smallest high-value change.

Goals:

- remove repetitive missing/invalid-field scaffolding
- keep parse diagnostics explicit
- make each parser readable at a glance

Note:

- shared formatting will likely slightly change parse error strings
- if that matters, pin it in tests or keep the helper output as close as
  possible to current wording

### Step 3: Replace the four stream structs with one `PartialsStream`

Unify:

- open logic
- blank-line skipping
- line buffering
- line-number tracking
- schema selection

Do not over-genericize it. A `match self.schema` inside `next_row()` is clearer
than four near-identical structs.

This step should also replace the per-schema parsed row types with
`ParsedPartialRow`.

### Step 4: Extract shared expected-contribution loading

Create one helper for cross-index loading. It should own:

- zstd/plain opening
- line iteration
- blank-line skipping
- parse failures
- default contribution count of `1` when a key is absent from all cross-index
  sidecars

This should be one function:

```rust
fn load_expected_contributions(
    files_by_tile: &FxHashMap<u32, TileFiles>,
) -> Result<FxHashMap<u64, u32>>
```

There is no BED vs size behavioral difference here. Both cross-index formats are
just `u64` keys counted per line.

### Step 5: Introduce one shared accumulator

Replace the four reducer accumulators with one additive accumulator.

Reason:

- it removes repeated accumulation code
- it aligns the reducer field set with grouped fold accumulation
- it keeps summary-only fields zero in basic mode without changing disk schema

Keep comments explicit that zero summary-only fields in basic mode are
intentional.

### Step 6: Collapse BED basic and BED summary into one merge engine

Create one internal BED reducer that is parameterized by summary mode, not by a
full schema enum.

Pseudo-shape:

```rust
fn reduce_bed_rows_internal(
    ...,
    summary: bool,
    on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()>
```

The function should map `summary` to the correct stream schema internally:

- `false` -> `PartialsSchema::BedBasic`
- `true` -> `PartialsSchema::BedSummary`

Responsibilities:

- open the right partial streams
- K-way merge by `orig_idx`
- accumulate additive fields
- wait until `seen_contributions == expected_contributions`
- recover the interval from `coords_by_idx`
- emit one `ReducedAggregateRow`

This should replace:

- `reduce_bed_basic_with_cross_index_for_chr_rows`
- `reduce_bed_with_cross_index_for_chr_rows`

Keep thin public wrappers only if they still improve call-site clarity.

### Step 7: Collapse size basic and size summary into one merge engine

Use the same pattern for fixed-size reduction.

Pseudo-shape:

```rust
fn reduce_size_rows_internal(
    ...,
    summary: bool,
    on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()>
```

The function should map `summary` to:

- `false` -> `PartialsSchema::SizeBasic`
- `true` -> `PartialsSchema::SizeSummary`

Responsibilities:

- open the right size partial streams
- K-way merge by full bin `start`
- accumulate additive fields
- wait until `seen_contributions == expected_contributions`
- preserve full-bin identity
- clip the final bin to chromosome end after reduction

This should replace:

- `reduce_aggregates_by_size_basic_with_cross_index_for_chr_rows`
- `reduce_aggregates_by_size_with_cross_index_for_chr_rows`

### Step 8: Simplify grouped writer folding

After the reducer accumulator is unified, simplify grouped folding in
`writers.rs`.

The goal is not to genericize the whole writer. The goal is narrower:

- one grouped accumulation path
- summary vs basic branching only at final-value derivation/output

This does not require sharing the reducer accumulator struct with grouped
folding. `GroupedAggregateAccum` should stay separate and keep
`span_positions`. The simplification comes from using the same additive field
set in both basic and summary grouped loops.

This is a second-pass cleanup, not the primary refactor driver.

### Step 9: Remove trivial finalize wrappers if they no longer help

If the non-`_rows` reducers are only single-call-site wrappers, inline them into
the writers.

These are:

- `reduce_bed_with_cross_index_for_chr`
- `reduce_aggregates_by_size_with_cross_index_for_chr`

If they still improve readability by naming the mode-specific intent, keep them.

This is optional.

## Writer-Side Scope

The main payoff is still in the reducer. Do not turn `writers.rs` into a
generic callback maze.

Reasonable writer-side cleanup:

- simplify grouped folding after reducer accumulator unification
- share obvious header/output helpers where they are already nearly identical
- keep summary derivation explicit at output time

Not reasonable:

- hiding BED/grouped/size modes behind deeply generic callbacks
- moving summary-stat math back into reducers
- making control flow harder to audit than the current explicit branches

## What Should Not Be Unified

These remain bad ideas even after the revision.

### Do not widen basic temp-file schemas

All unification in this plan is read-side or in-memory.

Wrong:

- always writing `nonzero_positions`
- always writing `coverage_sum_of_squares`

Right:

- basic temp files stay narrow
- summary-only fields appear only in memory as zero-filled fields when needed

### Do not hide BED vs size key semantics

BED and fixed-size reduction are not the same identity model.

- BED uses stable `orig_idx`
- size uses full bin `start`

That difference should stay visible in the internal engine signatures and
comments.

### Do not pretend BED parsed rows know their interval

BED partial files do not carry the final interval. Any refactor that acts as if
the interval is available directly from disk is conceptually wrong.

### Do not move summary-stat derivation back into reducers

Keep:

- reducers = exact additive raw rows
- writers/finalization = derived values and presentation

### Do not replace explicit code with trait-heavy parser abstractions

One `PartialsStream` plus a small schema enum is fine.

Deep generic abstractions are not. The point is to make the code easier to
trust, not harder to follow.

## Risk Areas To Watch

These are the places where the refactor can silently go wrong.

### Basic rows accidentally gain wider temp-file schemas

Guardrail:

- unification happens after reading the already-written temp rows

### Wrong identity key for reduction

Examples of broken refactors:

- BED keyed by interval instead of `orig_idx`
- size keyed by clipped overlap start instead of full bin start

Guardrail:

- keep key extraction explicit per schema

### Wrong BED interval recovery

Guardrail:

- only the BED merge engine should recover intervals from `coords_by_idx`

### Summary fields accidentally influence basic outputs

Guardrail:

- basic output paths must continue to ignore summary-only fields
- keep exact-output tests for basic modes

### Final-bin clipping regressions

Guardrail:

- keep targeted tests for final bin clipping in both size modes

### Parse helper subtly changes diagnostics

Guardrail:

- treat helper-driven parse error text as an explicit reviewed change

## Suggested Test Matrix

### Reducer identity tests

- BED basic: cross-tile same `orig_idx` merges exactly once
- BED summary: cross-tile same `orig_idx` merges exact raw additive fields
- size basic: full bin start is used as identity, not clipped overlap identity
- size summary: full bin start is used as identity, not clipped overlap identity

### Direct reducer round-trip tests

- synthetic BED partials -> reducer -> exact reduced rows
- synthetic size partials -> reducer -> exact reduced rows

These should verify the reducer itself rather than only writer behavior.

### Schema preservation tests

- basic BED partials still write only basic columns
- basic size partials still write only basic columns
- summary variants still write summary columns

### Basic vs summary additive equivalence tests

- the shared additive fields from basic reduction match the same fields from
  summary reduction for equivalent underlying data

### Wrapper equivalence tests

- raw-row reducer + `finalize_value` matches current basic final output
- raw-row reducer + summary derivation matches current summary final output

### Grouped fold tests

- grouped total from raw segment rows
- grouped summary stats from raw segment rows
- grouped + blacklist
- grouped + cross-tile segment reduction

### Invariance tests

- tile-size invariance for summary stats
- tile-size invariance for grouped summary stats
- cross-tile boundary invariance for size summary stats
- final-bin clipping invariance for size reduction

## Recommended Sequence For Actual Implementation

If the goal is maximum complexity reduction per unit risk, the implementation
sequence should be:

1. Add reducer tests that directly pin merge behavior
2. Introduce `parse_col`
3. Collapse the four stream structs into `PartialsStream`
4. Extract expected-contribution loading
5. Introduce one shared accumulator
6. Unify BED merge engines
7. Unify size merge engines
8. Simplify grouped writer folding
9. Inline trivial wrappers only if they are clearly no longer useful

The earlier sequence started by unifying accumulators and merge engines first.
That leaves too much parser duplication untouched and makes the higher-risk
changes land before the low-risk cleanup.

## Bottom Line

The right simplification is now:

- keep separate on-disk schemas
- keep reducer output as exact additive raw rows
- first remove parser and stream duplication
- then collapse the four merge paths into one BED engine and one size engine
- then simplify grouped folding as a downstream benefit of the shared
  accumulator

The wrong simplification is still:

- widening temp files
- hiding BED vs size identity semantics
- pretending BED parsed rows already know their interval
- moving summary-stat derivation back into reducers
- replacing explicit code with generic abstractions that are harder to audit

This revised plan is stricter than the earlier one about where the real bloat
is. The merge loops matter, but the parser/stream duplication is too large to
leave out of the plan.
