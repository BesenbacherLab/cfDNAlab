# fcoverage Accumulation Numerics Spec

## Summary

This spec covers the correctness-focused numerical work for `fcoverage`.

The original problem was that `fcoverage` accumulated tile-local fragment support in a global
`f32` delta array, finalized that into `f32` positional coverage, and then applied a residue
cleanup floor before genomic scaling. The fixed cleanup floor was wrong because it could exceed
the smallest real positive support allowed by some argument combinations, for example GC
correction combined with `--normalize-by-length`.

The current correctness direction is:

1. Keep the theoretical cleanup floor, but derive it from the active arguments.
2. Keep the tile-level accumulator structure unchanged for now, but perform the internal weight
   and delta arithmetic in `f64`.

The final positional coverage vector written by `fcoverage` must stay structurally identical to
today. Downstream reducers, indexing, scaling, masking, and output writers should continue to
consume the same finalized tile coverage layout.

Larger performance or memory optimizations are future work and are listed in:

- `plans_and_specs/FUTURE_OPTIMIZATIONS_SPEC_2026-04-18.md`

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

## Phase 2: Current Numerical Fix

### Goal

Remove the large mixed-sign `f32` delta accumulation error source without changing tile shapes,
reducers, output formats, or downstream consumers.

### Current Direction

Keep the existing tile-local delta accumulation structure, but:

- store the internal delta array in `f64`
- accept fragment weights as `f64`
- keep theoretical minimum-support derivation in `f64`
- keep the final positional coverage vector as `Vec<f32>`

This keeps the current per-tile algorithm recognizable while removing the main early precision
losses.

### Why This Is The Current Fix

- Removes the main source of repeated mixed-sign `f32` accumulation residue
- Keeps the implementation local to `Coverage` and the per-tile counting path
- Preserves the same final `Vec<f32>` positional output shape used by downstream code
- Avoids introducing a second internal tile/window size just for the correctness fix

### Deferred Optimization Work

Future optimization ideas are tracked separately in:

- `plans_and_specs/FUTURE_OPTIMIZATIONS_SPEC_2026-04-18.md`

### Invariants The Rewrite Must Preserve

- final coverage values should match today's exact counting semantics
- segment-aware counting must still count only the same reference positions as today
- `--ignore-gap` must still only affect which positions are counted, not the weight formula
- `--normalize-by-length` must still use the same counted-span denominator semantics as today
- GC correction and genomic scaling must still be applied in the same order as today
- blacklist masking must still happen after raw coverage construction, exactly as today
- downstream reducers must see the same tile-core positional coverage layout as today
- external file formats and file names must not change

### Required Equivalence Tests

The accumulator internals changed, even though the output structure stayed the same. We therefore
need strong coverage for:

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

- implement the `f64` delta and `f64` weight path inside the existing tile-local accumulator

Next planning step:

- audit missing test coverage around the current accumulation behavior

Future optimization work:

- keep rolling local delta or coverage accumulation as a separate future optimization track
- evaluate those changes on reducer IO, tile-size pressure, and implementation complexity before
  taking them on
