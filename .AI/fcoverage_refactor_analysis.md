# fcoverage reducer/writer refactor: independent analysis

This is an independent review of the code and the Codex-generated plan
(`tmp_fcoverage_complexity_reduction_plan.md`). It is based on reading every
line of `reducer.rs` (1619 lines), `writers.rs` (1144 lines), the relevant
parts of `fcoverage.rs` (1518 lines), `tiling.rs`, and `window_results.rs`.

## Where the plan is right

The plan correctly identifies:

- Four nearly-identical merge loops (BED basic, BED summary, size basic, size summary)
- Four accumulators that differ only by 2 fields
- Cross-index loading copy-pasted 4 times
- The hard constraints around temp-file width, row identity keys, and summary-stats math
- The right simplification boundary: one BED merge engine + one size merge engine
- The right things NOT to unify (BED vs size key semantics, on-disk schema, summary math placement)

## Where the plan is wrong or incomplete

### 1. The plan underestimates the parsing boilerplate problem

The plan says "keep separate parsers" and focuses on the merge loop. But the
biggest single source of bloat is the per-column parsing boilerplate. Each
column in each parser takes ~15 lines:

```rust
let field: T = cols
    .next()
    .ok_or_else(|| anyhow::anyhow!(
        "Missing FIELD in partials for chromosome '{}' tile {} line {}",
        self.chr, self.tile_index, self.line_number
    ))?
    .parse()
    .with_context(|| format!(
        "Invalid FIELD in chromosome '{}' tile {} line {}",
        self.chr, self.tile_index, self.line_number
    ))?;
```

There are **24 column parses** across 4 stream types (4+6+5+7 columns).
That is ~360 lines of mechanically identical error-wrapping code that differs
only in the field name string. A single `parse_col` helper would turn each
parse into one line:

```rust
let orig_idx: u64 = parse_col(&mut cols, "orig_idx", &self.chr, self.tile_index, self.line_number)?;
```

This alone would cut reducer.rs by ~300 lines without changing any behavior.

### 2. The plan's `ParsedContribution` duplicates `ReducedAggregateRow`

The plan proposes a new `ParsedContribution` struct with fields:
`key, interval, coverage_sum, eligible_positions, blacklisted_positions,
nonzero_positions, coverage_sum_of_squares`.

Compare with the existing `ReducedAggregateRow`:
`idx, interval, coverage_sum, eligible_positions, blacklisted_positions,
nonzero_positions, coverage_sum_of_squares`.

These are the same struct with one field renamed. For BED rows the interval
is not known at parse time (recovered from `coords_by_idx` later), but
that is already handled by the merge engine, not the parsed row type.

Recommendation: don't create `ParsedContribution`. Instead, use a
lightweight `ParsedPartialRow` that carries only what the parser reads
from disk, and let the merge engine build `ReducedAggregateRow` when
the window is complete:

```rust
struct ParsedPartialRow {
    key: u64,                       // orig_idx for BED, bin_start for size
    interval: Option<Interval<u64>>, // Some for size, None for BED
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,         // 0 when basic
    coverage_sum_of_squares: f64,   // 0.0 when basic
}
```

### 3. The plan doesn't address the 4 stream structs

The plan mentions separate stream reader structs as part of the duplication
problem but then says "keep separate parsers" without addressing the fact
that the 4 stream structs (`BasicPartialsStream`, `SummaryPartialsStream`,
`BasicSizePartialsStream`, `SummarySizePartialsStream`) share 100% of their
line-reading, blank-line-skipping, and file-opening logic.

Recommendation: one `PartialsStream` struct that handles IO, plus a schema
enum that drives column parsing:

```rust
enum PartialsSchema { BedBasic, BedSummary, SizeBasic, SizeSummary }

struct PartialsStream {
    reader: BufReader<Box<dyn Read + Send>>,
    line_buf: String,
    chr: String,
    line_number: u64,
    tile_index: u32,
    schema: PartialsSchema,
}
```

The `next_row(&mut self) -> Result<Option<ParsedPartialRow>>` method uses
`self.schema` to decide how many columns to parse and how to fill the
result struct. This replaces 4 struct definitions and 4 nearly-identical
`open()` + `next_row()` implementations with 1 struct, 1 `open()`, and
1 `next_row()` with a match on schema.

### 4. The plan misses an opportunity in the writer fold logic

`write_grouped_bed_aggregate_output` (writers.rs:570-635) has two
separate chromosome-iteration loops: one for summary, one for non-summary.
They differ only in which fields are accumulated (3 vs 5). With a unified
accumulator, both loops become identical and the branching moves entirely
to the output phase. The plan mentions "shared accumulator" for reducers
but doesn't explicitly note that it also collapses the grouped fold
branching in writers.rs.

### 5. The plan doesn't notice that 2 of the 6 public reducers are trivial wrappers

The call graph from writers.rs:

```
write_bed_aggregate_output
  summary   -> reduce_bed_with_cross_index_for_chr_rows
  basic     -> reduce_bed_with_cross_index_for_chr
                 -> reduce_bed_basic_with_cross_index_for_chr_rows + finalize + write

write_grouped_bed_aggregate_output
  summary   -> reduce_bed_with_cross_index_for_chr_rows
  basic     -> reduce_bed_basic_with_cross_index_for_chr_rows

write_size_aggregate_output
  aligned   -> concat_aligned_size_tile_finals (fast path, no reduction)
  summary   -> reduce_aggregates_by_size_with_cross_index_for_chr_rows
  basic     -> reduce_aggregates_by_size_with_cross_index_for_chr
                 -> reduce_aggregates_by_size_basic_with_cross_index_for_chr_rows + finalize + write
```

`reduce_bed_with_cross_index_for_chr` and
`reduce_aggregates_by_size_with_cross_index_for_chr` are thin wrappers
that call the basic `_rows` variant and apply `finalize_value` + `round_to`
+ `write_final_row` inside the callback. After unification these wrappers
become even thinner and could be inlined into the writer call sites. 
Whether to keep them is a judgment call, but the plan should acknowledge
they exist and are already trivial.

## Recommended concrete refactoring order

### Step 1: `parse_col` helper (saves ~300 lines, zero risk)

Add one small helper function at the top of reducer.rs:

```rust
fn parse_col<T: std::str::FromStr>(
    cols: &mut std::str::Split<'_, char>,
    field_name: &str,
    chr: &str,
    tile_index: u32,
    line_number: u64,
) -> Result<T>
where
    T::Err: std::fmt::Display,
{
    cols.next()
        .ok_or_else(|| anyhow::anyhow!(
            "Missing {} in partials for chromosome '{}' tile {} line {}",
            field_name, chr, tile_index, line_number
        ))?
        .parse::<T>()
        .map_err(|e| anyhow::anyhow!(
            "Invalid {} in chromosome '{}' tile {} line {}: {}",
            field_name, chr, tile_index, line_number, e
        ))
}
```

Then each column parse becomes one line. This immediately makes every
parser readable at a glance.

### Step 2: One `PartialsStream` + one `ParsedPartialRow` (saves ~400 lines)

Replace the 4 stream structs + 4 row structs with:
- `ParsedPartialRow` (as above)
- `PartialsStream` with `PartialsSchema`
- One `open()`, one `next_row()` that matches on schema

This also eliminates the 4 accumulator structs because the parsed row already
carries all fields (zeros for basic-mode summary fields).

### Step 3: `load_expected_contributions` helper (saves ~80 lines, zero risk)

Extract the cross-index loading loop into:

```rust
fn load_expected_contributions(
    files_by_tile: &FxHashMap<u32, TileFiles>,
) -> Result<FxHashMap<u64, u32>>
```

This is called identically in all 4 merge functions with zero behavioral
difference.

### Step 4: One unified accumulator (saves ~40 lines)

Replace `BasicWindowAccum`, `SummaryWindowAccum`, `BasicBinAccum`,
`SummaryBinAccum` with one `Accumulator`:

```rust
#[derive(Default)]
struct Accumulator {
    sum: f64,
    allowed_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
    seen_contributions: u32,
}
```

Basic-mode rows feed zeros into the last two fields. The accumulation
loop becomes one function. This does not widen temp files (they are
already written; this is read-side only).

### Step 5: One BED merge engine (saves ~150 lines)

Collapse `reduce_bed_basic_with_cross_index_for_chr_rows` and
`reduce_bed_with_cross_index_for_chr_rows` into:

```rust
pub(crate) fn reduce_bed_with_cross_index_for_chr_rows(
    chr: &str,
    temp_dir: &Path,
    partials_prefix: &str,
    windows_chr: &[IndexedInterval<u64>],
    summary: bool,
    on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()>
```

The `summary` flag selects `PartialsSchema::BedBasic` vs
`PartialsSchema::BedSummary` for the stream. Everything else is identical.

Keep `reduce_bed_with_cross_index_for_chr` as a thin public wrapper for
the non-summary finalize-and-write path, or inline it into the writer.

### Step 6: One size merge engine (same pattern, saves ~150 lines)

Collapse the two size reducers into one with a `summary: bool` parameter.

### Step 7: Simplify grouped fold in writers.rs

With the unified accumulator, the two chromosome-iteration loops in
`write_grouped_bed_aggregate_output` (lines 570-635) become one loop that
always accumulates all 5 fields. The summary vs non-summary branching
moves entirely to the output phase (which already has that branching).

### Step 8: Remove the non-`_rows` finalize wrappers (optional)

`reduce_bed_with_cross_index_for_chr` and
`reduce_aggregates_by_size_with_cross_index_for_chr` can be inlined
into their single call sites in `write_bed_aggregate_output` and
`write_size_aggregate_output`. This removes one level of indirection.

## Expected results

| Area | Before | After |
|------|--------|-------|
| reducer.rs | ~1619 lines | ~600-700 lines |
| writers.rs | ~1144 lines | ~1050-1100 lines (modest savings here) |
| Stream structs | 4 | 1 |
| Row structs | 4 | 1 |
| Accumulator structs | 4 | 1 |
| Merge loop copies | 4 | 2 (BED + size) |
| Cross-index loading copies | 4 | 1 |

Total savings: roughly 800-900 lines from reducer.rs with minimal
changes to writers.rs. The main complexity reduction is in the reducer.

## What NOT to do (agreeing with the plan)

- Don't widen on-disk temp file schemas
- Don't hide BED vs size key semantics behind a generic
- Don't move summary-stats derivation back into reducers
- Don't make the writer into a callback maze
- Don't create trait-heavy abstractions for the parsers

## Risks the plan correctly identifies

- Basic rows must not accidentally gain summary-stat columns in temp files
  (addressed: this refactor is read-side only)
- Wrong identity key for grouping
  (addressed: key extraction stays explicit per schema)
- Wrong interval recovery for BED
  (addressed: `coords_by_idx` lookup stays in the BED merge engine)
- Final bin clipping regression for size reducer
  (addressed: stays in the size merge engine)

## One additional risk the plan doesn't mention

The `parse_col` helper changes error messages from field-specific strings to
a templated format. If any downstream process pattern-matches on reducer error
messages, this would be a silent break. In practice this is unlikely (these are
internal temp-file parse errors), but worth noting.

## Test coverage assessment

The existing tests in `writers_tests.rs` cover summary-stats derivation math
thoroughly (variance repair, CV thresholds, zero-eligible rows, tiny positive
variance). This is the area most likely to regress from a refactor. Good.

What's missing and should be added before or during the refactor:

1. **Round-trip reducer test**: write synthetic partials files, run the reducer,
   verify exact output rows. This would catch regressions in the merge engine
   itself. Currently there's no unit test for the merge loop -- it's tested
   only through integration (full pipeline tests, if any).

2. **Cross-tile merge test**: synthetic scenario with 2 tile partials for the
   same `orig_idx`, verify they merge into exactly one output row with
   correct sums.

3. **Basic vs summary output equivalence**: verify that running the basic
   reducer on basic partials produces the same `coverage_sum`,
   `eligible_positions`, `blacklisted_positions` as running the summary
   reducer on summary partials for the same underlying data.

These tests should exist BEFORE the merge engine is unified so they can serve
as regression guards.
