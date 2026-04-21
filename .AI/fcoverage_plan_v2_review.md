# Review of updated fcoverage complexity reduction plan (v2)

The updated plan correctly incorporated every major point from the earlier
analysis. The refactor order, hard constraints, invariant summary, and
"what not to unify" sections are solid. Four issues remain.

## Issue 1: The shared accumulator claim is wrong for grouped folding

The plan says:

> Use one additive accumulator for both reducers and grouped folding

And shows one `AggregateAccum` struct. But the reducer accumulators and the
grouped accumulator serve different purposes and have different fields:

**Reducer accumulators** (4 copies today, 1 after unification):

```rust
struct AggregateAccum {
    coverage_sum: f64,
    eligible_positions: u64,
    blacklisted_positions: u64,
    nonzero_positions: u64,
    coverage_sum_of_squares: f64,
    seen_contributions: u32,       // <-- K-way merge completion tracking
}
```

**Grouped fold accumulator** (in writers.rs):

```rust
struct GroupedAggregateAccum {
    span_positions: u64,           // <-- comes from grouped layout, not reduction
    blacklisted_positions: u64,
    eligible_positions: u64,
    nonzero_positions: u64,
    coverage_sum: f64,
    coverage_sum_of_squares: f64,
    // no seen_contributions -- this folds already-reduced rows, not K-way streams
}
```

These are not the same struct. The reducer accumulator tracks K-way merge
completion (`seen_contributions`) and has no group-level span. The grouped
accumulator carries `span_positions` from the layout and doesn't do K-way
merging at all.

**What the plan should say instead:**

The 4 reducer accumulators unify into one `AggregateAccum`. That is the primary
win. `GroupedAggregateAccum` stays as a separate type because it serves a
different purpose.

The grouped fold simplification (Step 8) still works because the *additive
field set* is the same: once basic-mode reduced rows carry zero-valued
summary-only fields, the two grouped fold loops collapse into one loop that
always accumulates all 5 additive fields. The struct doesn't need to be shared
for this to work -- just the field names.

Fix: change the `AggregateAccum` description to say "replaces the 4 reducer
accumulators" instead of "for both reducers and grouped folding." Keep
`GroupedAggregateAccum` separate, or at most note that after the reducer
accumulator is settled, the grouped accumulator should use matching field names
for consistency.

## Issue 2: Merge engine should take `summary: bool`, not `PartialsSchema`

The plan proposes:

```rust
fn reduce_bed_rows_internal(
    ...,
    schema: PartialsSchema,
    on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()>
```

This accepts any `PartialsSchema`, including `SizeBasic` and `SizeSummary`,
which would be nonsensical for a BED engine. The type doesn't prevent misuse.

A simpler signature is more precise:

```rust
fn reduce_bed_rows_internal(
    ...,
    summary: bool,
    on_row: impl FnMut(ReducedAggregateRow) -> Result<()>,
) -> Result<()>
```

The function maps `summary` to the right schema internally
(`BedBasic` / `BedSummary`). Same for the size engine.

This is a small API design point, but it matters because passing a valid-typed
but semantically wrong enum variant silently produces garbage. A `bool` can't be
wrong in the same way.

Fix: change the merge engine signatures to `summary: bool` and derive the
schema inside the function.

## Issue 3: Cross-index loading doesn't need separate BED/size helpers

Step 4 says:

> Use explicit helpers if that stays clearer, for example one for BED keys and
> one for size keys

But the cross-index files contain one `u64` per line regardless of whether it's
an `orig_idx` (BED) or a `bin_start` (size). The loading logic is byte-for-byte
identical:

1. Open each tile's cross-index file (zstd or plain)
2. Read lines, skip blanks
3. Parse each line as `u64`
4. Count occurrences per key

One function suffices:

```rust
fn load_expected_contributions(
    files_by_tile: &FxHashMap<u32, TileFiles>,
) -> Result<FxHashMap<u64, u32>>
```

There is no BED/size behavioral difference in cross-index loading.

Fix: drop the hedge about separate helpers. Recommend one function.

## Issue 4: Minor clarity gaps

a) **Which wrappers does Step 9 refer to?** The plan says "inline trivial
   finalize wrappers" but doesn't name them. For someone implementing the plan,
   these are:
   - `reduce_bed_with_cross_index_for_chr` (calls the basic BED `_rows` reducer
     + `finalize_value` + `write_final_row`)
   - `reduce_aggregates_by_size_with_cross_index_for_chr` (same pattern for size)

b) **`allowed_positions` vs `eligible_positions`**: The current reducer
   accumulators use `allowed_positions` while `ReducedAggregateRow` and the
   grouped accumulator use `eligible_positions`. The plan's `AggregateAccum`
   uses `eligible_positions`, which normalizes the naming. This is a good
   change but should be noted explicitly so an implementer doesn't think it's
   an oversight.

## Summary

| Issue | Severity | Fix |
|-------|----------|-----|
| Shared accumulator claim wrong for grouped fold | Real design error | Separate the claim: reducer accum unifies 4 structs; grouped accum stays separate |
| Merge engine takes `PartialsSchema` | API misuse risk | Use `summary: bool`, derive schema internally |
| Cross-index helper hedge | Unnecessary complexity | One function, no hedge |
| Unnamed wrappers, implicit rename | Clarity | Name the functions, note the rename |

Everything else in the plan is correct and ready to implement.
