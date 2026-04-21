# fcoverage Complexity Reduction Plan

This document proposes how to reduce the current `fcoverage` reducer and writer complexity
without changing behavior, weakening invariants, or writing more data than necessary to
temporary files.

The main idea is not to undo the recent architectural direction. The useful part is:

- reducers can produce exact raw aggregate rows
- summary-stat derivation can happen after reduction
- grouped outputs can reuse the same raw reduced rows
- basic outputs can still avoid summary-only temp-file columns

The part that should be reduced is the duplicated merge machinery.

## Hard Constraints

These are non-negotiable. Any refactor that violates one of them is wrong.

- Keep the temp-file width optimization
  Basic partial rows must stay narrow
  Summary partial rows may carry extra additive fields
  Do not switch back to always writing `nonzero_positions` and `coverage_sum_of_squares`

- Keep persisted row identity unchanged
  BED reduction must still key by original BED `orig_idx`
  Fixed-size reduction must still key by full bin `start`
  Partial rows must still preserve full window/bin identity, not tile-local overlap identity

- Keep cross-index semantics unchanged
  Missing from all cross-index files means exactly one expected contribution
  Presence count across cross-index sidecars means exact expected contribution count

- Keep grouped behavior unchanged
  Group folding must still add already-reduced segment rows into group accumulators
  Group span must still come from the grouped layout, not from reduced segment lengths

- Keep summary-stat math unchanged
  Negative variance handling, CV near-zero handling, and `>1e6` display behavior must stay as-is

- Keep fast path behavior unchanged
  Aligned `--by-size` outputs must still use concatenation instead of cross-tile reduction

- Do not silently change public output schema unless explicitly intended
  Any header rename or formatting change should be treated as a separate decision

## Invariant Summary

This is the invariant summary that should be carried into any actual code edits.

- A partial row represents one tile's contribution to one logical BED window or one logical
  fixed-size bin
- Cross-tile correctness comes from:
  - grouping by stable row identity
  - accumulating additive raw fields
  - waiting until `seen_contributions == expected_contributions`
- The merge heap only chooses which open stream to read from next
  It does not prove that a row is complete
- Basic mode and summary-stats mode differ in temp-file schema, not in row identity rules
- BED and fixed-size reduction differ in key type and interval recovery
  They do not differ in the high-level reduction algorithm

## Current Complexity Problem

The current code pays for the useful architecture with too much duplicated control flow.

Today there are effectively four copies of the same merge pattern:

- BED basic reducer
- BED summary reducer
- size basic reducer
- size summary reducer

And each copy brings along matching duplicate pieces:

- stream reader structs
- parsed row structs
- accumulators
- cross-index counting loops
- heap setup
- merge loops
- finalize-on-complete logic

Some duplication is justified. Some is not.

### Duplication That Is Justified

- Separate temp-file schemas for basic vs summary
- Separate BED vs fixed-size parsing
- Separate public wrapper functions for:
  - raw-row reduction
  - final numeric output writing

These separations reflect real behavioral differences.

### Duplication That Is Not Justified

- Four distinct merge loops
- Separate basic/summary accumulators when one shared accumulator can hold zero-valued
  summary-only fields
- Repeated cross-index counting logic
- Repeated heap setup and stream advancement logic

## Target Shape

The optimal target is:

- Keep separate parsers for each on-disk schema
- Keep raw-row public reducer entry points
- Collapse the actual reduction algorithm to:
  - one internal BED merge engine
  - one internal fixed-size merge engine

This preserves the temp-file optimization while removing most of the algorithm duplication.

### Shared Internal Data Shape

Use one internal additive row shape for reduced rows and one internal accumulator shape.

Suggested shape:

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

struct AggregateAccum {
    coverage_sum: f64,
    eligible_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
    seen_contributions: u32,
}
```

Basic parsers populate the summary-only fields with zeros.

That does not increase temp-file size. It only standardizes the in-memory reduced row shape.

### Separate Parsing, Shared Reduction

Keep these separate parsers:

- BED basic partial row parser
- BED summary partial row parser
- size basic partial row parser
- size summary partial row parser

Each parser should map on-disk columns into a common internal contribution shape:

```rust
struct ParsedContribution {
    key: u64,
    interval: Interval<u64>,
    coverage_sum: f64,
    eligible_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
}
```

For BED rows:

- `key = orig_idx`
- `interval` is recovered later from `coords_by_idx`, or supplied by a callback

For size rows:

- `key = interval.start()`
- `interval` is already present in the partial row

## Recommended Refactor Plan

### Step 1: Freeze Behavior Before Simplifying

Before changing the structure, explicitly pin the current behavior that must survive:

- basic BED reducer output
- summary BED raw-row output
- basic size reducer output
- summary size raw-row output
- grouped folding behavior
- aligned size fast path

If any of these are not already directly covered, add focused tests first.

This step matters because once the merge logic is unified, a regression can hit all modes at once.

### Step 2: Introduce One Shared Accumulator

Replace:

- `BasicWindowAccum`
- `SummaryWindowAccum`
- `BasicBinAccum`
- `SummaryBinAccum`

with one shared accumulator that always contains the full additive field set.

Reason:

- the in-memory cost is tiny compared to tile coverage arrays and file IO
- it removes repeated accumulation code
- it does not widen the persisted temp-file schema

Important:

- keep comments explicit that zero-valued summary-only fields in basic mode are intentional

### Step 3: Factor Cross-Index Counting Into Helpers

Create helpers such as:

```rust
fn load_expected_contributions_for_orig_idx(...)
fn load_expected_contributions_for_bin_start(...)
```

or, if it stays readable:

```rust
fn load_expected_contributions<T: FromStr>(...)
```

I recommend two explicit helpers rather than a generic one if the generic version starts
obscuring the diagnostics.

Goal:

- one source of truth for:
  - zstd/plain handling
  - blank-line skipping
  - parse errors
  - default `1` contribution rule

### Step 4: Introduce One BED Merge Engine

Create one internal function that performs the full BED merge algorithm independent of
basic vs summary mode.

It should take:

- a stream factory or iterator of BED partial streams
- the expected-contribution map
- the `coords_by_idx` lookup
- a callback for each completed `ReducedAggregateRow`

Pseudo-shape:

```rust
fn reduce_bed_rows_internal(
    ...,
    stream_kind: BedPartialsKind,
    on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()>
```

Where `BedPartialsKind` selects the parser, not the merge algorithm.

Behavior:

- basic mode parser fills zero summary-only fields
- summary mode parser fills real summary fields
- shared merge loop accumulates the same `ReducedAggregateRow` shape

This should replace both:

- `reduce_bed_basic_with_cross_index_for_chr_rows`
- `reduce_bed_with_cross_index_for_chr_rows`

while preserving both public entry points as thin wrappers if keeping them is useful.

### Step 5: Introduce One Fixed-Size Merge Engine

Do the same for fixed-size reduction.

One internal function should own:

- expected-contribution counting by full bin start
- K-way merge by full bin start
- accumulation by start
- clipping the final bin to chromosome end

Basic vs summary should again be decided only at parse time.

This should replace both:

- `reduce_aggregates_by_size_basic_with_cross_index_for_chr_rows`
- `reduce_aggregates_by_size_with_cross_index_for_chr_rows`

while preserving any public wrappers needed by callers.

### Step 6: Keep Thin Public Wrappers

After the shared engines exist, the public reducers should become intentionally small.

Recommended public surface:

- raw BED row reducer
- raw size row reducer
- final BED aggregate writer reducer for average/total
- final size aggregate writer reducer for average/total

The final numeric wrappers should only:

- validate mode
- call the raw reducer
- run `finalize_value(...)`
- round
- write the final row

They should not duplicate merge mechanics.

### Step 7: Reduce Writer-Orchestration Duplication

After reducer simplification, do a second pass over `writers.rs`.

There is currently repeated branching across:

- BED aggregate writing
- grouped aggregate writing
- size aggregate writing

Some duplication is still warranted because grouped mode has a real second-stage fold.
But some repeated shape can be reduced:

- shared header selection
- shared summary-stat derivation call sites
- shared chromosome iteration helpers

Do not over-abstract this step. The main payoff is in reducer simplification, not in turning
`writers.rs` into a generic callback maze.

## What Should Not Be Unified

To avoid overengineering, do not try to unify these too aggressively.

### Do Not Force One Universal Partials Parser

The file formats are genuinely different.

Trying to parse them all with:

- one enum-heavy parser
- one dynamic column switch
- one generic trait soup

is likely to make the code harder to read than the current duplication.

Separate parser entry points are fine.

### Do Not Hide BED vs Size Key Semantics

BED and fixed-size rows differ in the identity key used for reduction:

- BED uses stable `orig_idx`
- fixed-size uses full bin `start`

That difference is fundamental and should stay explicit.

### Do Not Move Summary Math Back Into Reducers

That would collapse responsibilities again and make grouped folding less clean.

Keep:

- reducers = exact additive raw rows
- writers/finalization = derived presentation-level values

## Risk Areas To Watch During Refactor

These are the places where simplification can silently go wrong.

### Basic Rows Accidentally Gain Wider Temp-File Schema

This is the main risk against the goal you stated.

Guardrail:

- parser unification must happen after read-time, not by widening the written temp rows

### Wrong Identity Key For Grouping

Examples of incorrect refactors:

- grouping size rows by clipped tile-local start instead of full bin start
- grouping BED rows by start/end instead of `orig_idx`

Guardrail:

- keep the grouping key explicit in the internal engine signature

### Wrong Interval Recovery For BED

BED partial rows do not carry their interval in the file. The reducer reconstructs it from
`coords_by_idx`.

Guardrail:

- do not blur the distinction between:
  - row identity key
  - interval lookup source

### Summary Fields Accidentally Influence Basic Outputs

If shared accumulators are introduced, basic rows will carry zero-valued summary fields in memory.
That is fine. But basic-mode output must still ignore them.

Guardrail:

- keep basic output wrappers explicit
- keep tests that compare exact output rows for basic modes

### Final Bin Clipping Regressions

The size reducer still needs to clip the last bin to chromosome end after reduction.

Guardrail:

- keep targeted tests for final-bin clipping across both basic and summary size reducers

## Suggested Test Matrix For This Refactor

Claude should review whether each of these is already covered well enough.

### Reducer Identity Tests

- BED basic: cross-tile same `orig_idx` merges exactly once
- BED summary: cross-tile same `orig_idx` merges exact raw moments
- size basic: full bin start is used as identity, not clipped overlap
- size summary: full bin start is used as identity, not clipped overlap

### Schema Preservation Tests

- basic BED partials still write only basic columns
- basic size partials still write only basic columns
- summary variants still write summary columns

### Wrapper Equivalence Tests

- raw-row reducer + `finalize_value` wrapper matches current basic final output
- raw-row reducer + summary derivation matches current summary final output

### Grouped Fold Tests

- grouped total from raw segment rows
- grouped summary-stats from raw segment rows
- grouped + blacklist
- grouped + cross-tile segment reduction

### Invariance Tests

- tile-size invariance for summary-stats
- tile-size invariance for grouped summary-stats
- cross-tile boundary invariance for size summary-stats

## Recommended Review Questions For Claude

When asking Claude to review this plan, the useful questions are:

1. Does the proposed target shape keep temp-file schemas narrow in basic mode?
2. Is the suggested shared accumulator safe, or does it create any hidden semantic risk?
3. Is one BED engine + one size engine the right simplification boundary?
4. Are there any invariants in the current reducer code that this plan forgot to preserve?
5. Are there tests missing from the suggested matrix that would make this refactor safer?

## Bottom Line

The right simplification is:

- keep separate on-disk schemas
- keep raw-row reduction architecture
- remove duplicated merge engines

The wrong simplification is:

- widening temp files to avoid parser differences
- hiding BED vs size identity semantics
- reintroducing final-value logic into reducers
- replacing explicit code with generic abstractions that are harder to trust

If this plan is accepted, the first implementation step should be the smallest high-value one:

- unify accumulators
- factor expected-contribution loading
- then unify BED merge logic

Only after that should the size reducers be unified in the same style.
