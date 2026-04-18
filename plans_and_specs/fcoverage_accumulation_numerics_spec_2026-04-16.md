# fcoverage Accumulation Numerics Spec

## Summary

`fcoverage` currently accumulates tile-local fragment support in a global `f32` delta array,
finalizes that into `f32` positional coverage, and then applies a residue cleanup floor before
genomic scaling. The current fixed cleanup floor is wrong because it can exceed the smallest real
positive support allowed by some argument combinations, for example GC correction combined with
`--normalize-by-length`.

This spec separates the work into two phases:

1. Immediate fix:
   Replace the fixed cleanup floor with a per-run threshold derived from the smallest real
   positive pre-scaling support allowed by the active arguments.

2. Proper numerical fix:
   Replace the current tile-local `f32` delta accumulation with a rolling local `f64`
   accumulation design so fake support is prevented at the source rather than cleaned up later.

The final positional coverage vector written by `fcoverage` must stay structurally identical to
today. Downstream reducers, indexing, scaling, masking, and output writers should continue to
consume the same finalized tile coverage layout.

## Current Code Shape

The current per-tile `fcoverage` path is:

1. Build a tile-local `Coverage` object with an `f32` delta array
2. Add clipped fragment or segment spans into that delta array
3. Finalize the delta array into an `f32` per-position coverage vector
4. Clamp tiny positive residue to zero with a fixed threshold
5. Apply genomic scaling in place
6. Reuse the finalized coverage vector for positional output or for prefix-sum indexing

The current cleanup threshold is a fixed constant:

```rust
const INTERNAL_RESIDUAL_COVERAGE_FLOOR: f32 = 1.0e-4;
```

That is invalid because the smallest real positive support depends on the arguments.

## Phase 1: Dynamic Theoretical Cleanup Floor

### Goal

Keep the cleanup step, but make it impossible for it to remove any real positive pre-scaling
support that the current argument combination could legitimately produce.

### Derivation

The cleanup runs before genomic scaling, so the relevant support bound is the smallest real
positive **pre-scaling** support at one counted position.

Current pre-scaling multiplicative components are:

- Base weight:
  - `1.0` by default
  - `1.0 / counted_length` when `--normalize-by-length`
- GC weight:
  - `1.0` when GC correction is off
  - minimum usable positive GC weight `1e-3` when GC correction is on

Safe lower bound for the smallest real positive pre-scaling support in a run:

```text
min_positive_support =
    min_positive_base_weight * min_positive_gc_weight
```

with:

```text
min_positive_base_weight =
    if normalize_by_length:
        1 / max_fragment_length
    else:
        1

min_positive_gc_weight =
    if gc_file or gc_tag:
        1e-3
    else:
        1
```

This is safe because counted span length can only reduce from the fragment span, and the fragment
span is already bounded by `max_fragment_length`.

### Floor Rule

Use a cleanup floor strictly below that theoretical minimum, for example:

```text
cleanup_floor = min_positive_support / 2
```

This keeps the cleanup below any real positive support while still removing impossible smaller
residue.

### Important Maintenance Rule

Any future multiplicative weighting or normalization that can lower pre-scaling per-position
support must update this helper. The helper and its tests should explicitly say so.

### Immediate Test Coverage

Add tests for:

- No GC, no length normalization
- GC enabled, no length normalization
- Length normalization enabled, no GC
- GC enabled and length normalization enabled

These tests should verify both:

- the theoretical minimum positive support
- the derived cleanup floor is strictly smaller

## Phase 2: Proper Numerical Fix With Rolling Local f64 Accumulation

### Goal

Prevent fake positional support from being created in the first place, instead of depending on
post-finalization cleanup.

### Preferred Direction

Use a rolling local `f64` accumulation scheme and finalize directly into the final tile coverage
vector.

The preferred end state after the recent design discussion is:

- keep the current external `fcoverage` behavior unchanged
- keep full tile-level parallelization unchanged
- keep the final positional coverage vector structurally identical to today
- change only how that final tile-local coverage vector is numerically produced

This means the rewrite should happen inside the per-tile accumulation path, not by changing
reducers, output formats, or downstream consumers.

### Why This Is Better

- Removes the main source of repeated mixed-sign `f32` delta accumulation
- Keeps real low support intact
- Still ends in the same final `Vec<f32>` tile coverage representation used today

### Ordering Prerequisite

Any rolling flush design depends on knowing when a prefix is settled, so the first technical
question is the ordering guarantee of fragments fed into `process_tile`.

The rewrite must start by making this explicit:

- either prove and document that the per-tile fragment stream is nondecreasing by fragment start
- or prove and document a bounded-disorder guarantee
- or add a small explicit reorder buffer keyed by fragment start before the rolling accumulator

This must not be left implicit. The safe flush boundary depends on it.

### Preferred Rewrite Shape: Direct Rolling `f64` Coverage

The current best direction is to accumulate directly into rolling `f64` positional coverage rather
than into a second delta representation.

State:

- `final_coverage: Vec<f32>` for the full tile core
- `local_coverage: Vec<f64>` for a rolling genomic window inside the tile core
- `buffer_start_abs`
- `buffer_end_abs`
- a small ordering helper if the fragment stream is not strictly start-sorted

Per fragment:

1. Clip the fragment or segments to the tile core exactly as today
2. Determine the counted tile-local spans exactly as today
3. Ensure the rolling buffer covers those spans
4. Add the fragment weight directly to `local_coverage[start..end]`
5. Flush any settled prefix of `local_coverage` into `final_coverage`

End of tile:

6. Flush the remaining tail into `final_coverage`
7. Convert the completed `final_coverage` into the same downstream `Coverage`-compatible
   representation the rest of the tile code expects

This removes the repeated mixed-sign `f32` add/subtract pattern that currently creates the fake
support.

### Why Direct Coverage Is Preferred

- It directly computes the quantity we actually want at the end of the tile
- It avoids carrying a second delta lifecycle inside the rolling buffer
- It makes the numerical story easier to reason about: only positive weighted additions happen in
  the local `f64` state
- Any eventual cast back to `f32` happens once per finalized position, not on every fragment edge

### Fallback Variant: Rolling `f64` Delta With Limited `f32` Writeback

If the direct-coverage flush logic proves too invasive, the fallback is:

- keep local state as `f64` delta over a rolling window
- finalize only the settled prefix locally
- write the finalized prefix to the final tile coverage output
- if any overlap region must be written back to a global `f32` delta, keep that overlap bounded and
  explicitly floor only impossible near-zero edge residue there

This is less clean than direct coverage, but still much better than letting the entire tile live in
global `f32` delta state.

### Flush Trigger

The flush is not "flush after every fragment". It means:

- after each fragment insertion, check whether a prefix of the rolling buffer can no longer be
  touched by any future fragment
- flush exactly that settled prefix
- then slide the buffer forward

The exact settled-prefix rule depends on the fragment-ordering guarantee above. That rule must be
written down explicitly in code comments before the rewrite is considered complete.

### Concrete Rewrite Steps

1. Audit and document the current per-tile fragment ordering guarantee.
2. Add targeted tests or assertions for that ordering guarantee, or add an explicit reorder buffer.
3. Extract the current tile-local fragment insertion logic behind a small helper boundary so the
   new accumulator can reuse the same clipped spans and weights.
4. Introduce a dedicated rolling-accumulator type for `fcoverage` tile counting.
5. Implement the direct rolling `f64` coverage variant inside that accumulator.
6. Keep the rest of `process_tile` unchanged: blacklist application, scaling, indexes, and output
   writing should continue to operate on the same final tile coverage shape.
7. Keep the dynamic theoretical cleanup floor in place initially as a safety net, but make it
   diagnostic rather than relied on. If it still fires materially after the rewrite, that is a
   bug to investigate.
8. Once equivalence is established, decide whether the cleanup floor should remain as a redundant
   guard or be reduced further.

### Invariants The Rewrite Must Preserve

- final coverage values should match today's exact counting semantics
- segment-aware counting must still count only the same reference positions as today
- `--ignore-gap` must still only affect which positions are counted, not the weight formula
- `--normalize-by-length` must still use the same counted-span denominator semantics as today
- GC correction and genomic scaling must still be applied in the same order as today
- blacklist masking must still happen after raw coverage construction, exactly as today
- downstream reducers must see the same tile-core positional coverage layout as today
- external file formats and file names must not change

### Why This Needs A Plan

This is a meaningful rewrite of the currently tested accumulation path. The output structure stays
the same, but the internals change in a way that needs careful equivalence testing.

### Required Equivalence Tests Before Rewrite

Before touching the accumulation internals, confirm we have strong coverage for:

- Plain fragment spans
- Segment-aware fragments
- `--ignore-gap`
- `--normalize-by-length`
- GC file correction
- GC tag correction
- Blacklist masking
- Positional output
- By-bed aggregates
- By-size aggregates
- Tile boundary edge cases

If any of these are weak, strengthen them first so the rewrite can be checked afterwards.

### What This Rewrite Is Not

- It is not a change to the public meaning of `fcoverage`
- It is not a change to tile shapes, output rows, or reducer semantics
- It is not a reason to keep more whole-tile state in memory than necessary
- It is not a justification for changing downstream code that already works

## Decision

Immediate next step:

- implement the dynamic theoretical cleanup floor now

Next planning step:

- audit missing test coverage around the current accumulation behavior

Then:

- implement the ordering prerequisite
- implement the rolling local `f64` direct-to-coverage rewrite inside the tile accumulator
- keep the bounded-overlap `f64` delta writeback variant only as a fallback if the direct coverage
  version turns out to be materially more invasive than expected
