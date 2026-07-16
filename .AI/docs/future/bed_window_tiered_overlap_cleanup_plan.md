# BED Window Tiered Overlap Cleanup Plan

## Purpose

The current broad/narrow BED overlap optimization works, but it is layered on top of the old BED
window flow. `lengths` and `ends` now pass both the original chromosome-local BED window slice and
the split narrow/broad slices, plus separate tile spans for each. That makes the command flow harder
to reason about and keeps forcing lookups back through the original slice to recover the identity
needed for counting.

This plan describes the cleaner implementation: build one BED window collection per chromosome,
store the counting identities directly with each tiered window entry, precompute one tile-span
object per tile, and let `lengths` and `ends` use the same tiered BED overlap path.

## Goals

- Apply the tiered BED overlap finder to all explicit BED modes used by `lengths` and `ends`,
  including grouped BED.
- Keep the exact same candidate window set as the current finder. Result order is not part of the
  contract and tests should compare order-insensitive signatures.
- Avoid per-fragment sorting or other per-fragment cleanup that would erase the runtime gain.
- Stop overloading `IndexedInterval::idx()` in the optimized BED path. It should be clear whether an
  index is the chromosome-local source-window position, the original BED row, or the grouped-BED
  group index.
- Keep scaling-bin selection separate from count-window identity. Scaling helpers may need to
  preserve count-window identity, but they must not treat scaling-bin indices as count-window
  indices.
- Make the implementation general enough to support more than two length tiers later, while starting
  with the current `100_000 bp` broad threshold.

## Non-Goals

- Do not add a user-facing threshold option.
- Do not add per-fragment sorting to match old overlap order.
- Do not replace the overlap finder with an interval tree.
- Do not change output semantics for BED, grouped BED, fixed-size, or global modes.
- Do not change how scaling bins are selected or averaged. The only scaling-related change should be
  preserving the selected count-window identity through the scaling helper output if needed.

## Current Flow To Replace

Both `lengths` and `ends` currently do this in BED-like modes:

1. Build `indexed_windows_map`, where `IndexedInterval::idx()` is the original BED row id for plain
   BED or the group index for grouped BED.
2. Build `source_indexed_bed_windows_map`, where the same interval type is reused but
   `IndexedInterval::idx()` is changed to mean chromosome-local source-window position.
3. Precompute one tile span over the original window slice.
4. Precompute separate tile spans over the narrow split and the broad split.
5. Pass all of these into `process_tile`.
6. Use the split lists for optimized overlap lookup, but still pass the original window slice for
   fetch narrowing, count allocation, grouped BED lookup, and output id lookup.

The cleanup should replace this with a single BED window input object per chromosome plus a single
per-tile span object.

## Proposed Data Model

Use explicit BED window entries instead of storing different meanings in `IndexedInterval::idx()`.

```rust
struct BedWindowEntry {
    interval: Interval<u64>,
    source_window_idx: usize,
    output_idx: u64,
}
```

Meanings:

- `source_window_idx`: chromosome-local position in the original start-sorted BED-like window list.
  This is the position needed for plain BED tile-local count allocation and crossing-window
  tracking.
- `output_idx`: downstream output id. For plain BED this is the original BED row id. For grouped BED
  this is the group index.

The chromosome-level window collection should own the source list and the tier lists:

```rust
struct ChromosomeBedWindows {
    source_windows: Vec<BedWindowEntry>,
    tiers: Vec<BedWindowTier>,
}

struct BedWindowTier {
    windows: Vec<BedWindowEntry>,
}
```

Initial thresholds:

```text
narrow: [0, 100_000)
broad:  [100_000, +inf)
```

This is intentionally represented as tiers rather than hard-coded `narrow` and `broad` fields so a
future `[0, 10_000)`, `[10_000, 100_000)`, `[100_000, +inf)` split can be added without changing the
command flow again.

Keep threshold bounds in the builder, not in every hot-path tier. The overlap context only needs the
tier order and the windows inside each tier. The initial builder can use a private tier spec such as
`[(0, Some(100_000)), (100_000, None)]` and produce the lean `BedWindowTier` values above.

`source_windows` is still useful even though the overlap scan uses tiers. It supports fetch-span
narrowing and plain BED count allocation in source-window order. The tier entries carry their own
identity, so overlap hits do not need to look back into `source_windows` just to find the output id.

The tier lists should duplicate `BedWindowEntry` values instead of storing source indices and reading
indirectly during overlap lookup. That keeps the hot path direct and makes the identity carried by a
hit explicit at the point where the hit is found.

## Per-Tile Inputs

Precompute one tile BED span object instead of parallel span arrays:

```rust
struct TileBedWindowSpans {
    source_span: Option<TileWindowSpan>,
    tier_spans: Vec<Option<TileWindowSpan>>,
}

struct TileBedWindowInputs<'a> {
    chromosome_windows: &'a ChromosomeBedWindows,
    spans: &'a TileBedWindowSpans,
}
```

`process_tile` should receive either `None` for non-BED modes or one `TileBedWindowInputs` for BED
modes. It should not receive `windows_chr`, `source_indexed_windows_for_chr`,
`narrow_tile_window_span`, and `broad_tile_window_span` as separate arguments.

The existing `windows_map` and `grouped_windows_map` can remain at the outer command level for
metadata writers and final output setup. The tile processor should not need them once the BED window
collection carries source and output identities explicitly.

If `precompute_tile_window_spans` needs to work with both `IndexedInterval<u64>` and
`BedWindowEntry`, prefer a small crate-private bounds trait or helper over copying the span
precomputation logic. Keep existing non-BED command call sites unchanged.

## Tiered Overlap Context

Replace the current narrow/broad-specific context with a tier-based context:

```rust
struct TileBedOverlapContext<'a> {
    chrom_len: u64,
    always_hit_windows: Vec<BedWindowEntry>,
    tiers: Vec<TileBedOverlapTier<'a>>,
}

struct TileBedOverlapTier<'a> {
    windows: &'a [BedWindowEntry],
    window_ptr: usize,
}
```

Construction:

1. Slice each tier by the tile's corresponding tier span.
2. For the broadest tier, move windows that contain the full tile assignment envelope into
   `always_hit_windows`.
3. Leave narrower tiers in the normal scan. If a test fixture has tile sizes below the broad
   threshold, a narrow window that contains the whole envelope will still be checked by the normal
   tier scan. That is correct and avoids extra classification work on the common path.

Lookup:

1. Add all `always_hit_windows` with overlap fraction `1.0`.
2. Run the same pointer-retirement and exact-overlap loop for each tier:
   `window.end <= query_start.saturating_sub(look_back)` retires windows, and
   `window.start < query_end` bounds the forward scan.
3. Return the union of all hits. The order is tier order plus scan order within each tier and is not
   guaranteed to match the old single-slice finder.

Each returned BED hit should preserve at least:

```rust
source_window_idx: usize
output_idx: u64
interval: Interval<u64>
overlap_fraction: f64
```

Use the existing `OverlappingWindow` shape with one optional output identity field, rather than
adding a separate BED-only overlap row type:

```rust
struct OverlappingWindow {
    idx: usize,
    output_idx: Option<u64>,
    interval: Interval<u64>,
    overlap_fraction: f64,
}
```

Meanings:

- `idx` remains the source-window or bin position used by existing plain BED, fixed-size, and
  scaling-bin code.
- `output_idx` is present only when the overlap row already knows the BED-like downstream row id.
  BED constructors set it from `BedWindowEntry::output_idx`. Generic/global/fixed-size/scaling-bin
  rows leave it as `None` and keep using the existing mapping logic.

This is the smallest change that removes grouped BED lookup through a separate chromosome slice,
while avoiding a false universal-output-id invariant for fixed-size windows and scaling bins. It
also keeps the existing scaling helper input type and most command code shape.

## Command Flow After Cleanup

### Shared Setup

For BED and grouped BED modes:

1. Load the existing `windows_map` or `grouped_windows_map`.
2. Build `BedWindowsByChr = FxHashMap<String, ChromosomeBedWindows>`.
3. Precompute `Vec<TileBedWindowSpans>` in one helper call. Internally this can call
   `precompute_tile_window_spans` once for `source_windows` and once per tier, but the command only
   handles one returned span object.
4. In the tile loop, look up the chromosome's BED windows and the tile span object and pass them as
   one BED input.

For non-BED modes, keep the existing global/fixed-size flow.

### `lengths`

Use `TileBedWindowInputs` for:

- BAM fetch narrowing through `source_span` and `source_windows`.
- Plain BED count allocation through `source_span` and `source_window_idx`.
- BED and grouped BED overlap lookup through `TileBedOverlapContext`.
- Grouped BED counting through `output_idx`, without looking back into `windows_chr`.

Scaling remains fragment-level or overlap-level as today. If the scaling helper returns
`WindowScaling`, that struct must preserve the selected count-window identity from the input overlap
row. It does not need to know anything about which scaling bin was used beyond the existing scaling
weight calculation.

### `ends`

Use the same `TileBedWindowInputs` shape for:

- BAM fetch narrowing.
- BED and grouped BED overlap lookup.
- Output row lookup in the sparse count records.

For BED-like modes, `count_fragment_in_window(...)` should receive `output_idx` directly from the
overlap or scaling row. `WindowContext` can still handle global and fixed-size modes, but BED-like
tile code should no longer need a separate `windows_chr` lookup just to map source index to output
id.

## Scaling Identity Rule

Scaling-bin indices and count-window indices are separate.

The scaling overlap finder still returns scaling-bin positions into `scaling_chr`. Those indices are
only for averaging the scaling factor over the aligned reference span.

Count-window overlap rows represent selected output windows. If a scaling helper transforms
count-window overlap rows into `WindowScaling` rows, it must copy the count-window identity from each
input row to the corresponding output row. That identity is still about the count window, not the
scaling bin.

Practical implication:

- `WindowScaling.window_idx` can remain the source-window or bin index used by existing fixed-size
  and plain BED count allocation.
- Add `WindowScaling.output_idx: Option<u64>` and copy it from the input count-window overlap row.
- `compute_per_window_scaling_over_overlap` should still validate that remapped scaling rows match
  count rows by `idx`, `output_idx`, and overlap fraction.
- BED-like output identity should never be recovered from another slice in `process_tile`. BED-like
  command paths should require `Some(output_idx)` and return a helpful error if it is missing.

## Implementation Steps

1. Add the BED window collection structs and builder in shared code near the existing
   overlap/windowing helpers.
2. Build `ChromosomeBedWindows` from both plain BED windows and grouped BED windows, preserving
   `source_window_idx` and `output_idx`.
3. Add a helper that precomputes `TileBedWindowSpans` for all tiles from `BedWindowsByChr`.
4. Refactor `TileBedOverlapContext` from hard-coded narrow/broad fields to tier fields, keeping
   always-hit classification broadest-tier-only.
5. Extend `OverlappingWindow` and `WindowScaling` with optional `output_idx`, and update generic
   constructors so existing fixed-size/global/scaling-bin paths leave it unset and keep their current
   behavior.
6. Update the BED overlap return path so hits carry explicit source and output identities.
7. Refactor `lengths` setup and `process_tile` to accept one optional BED input object instead of
   the current parallel window arguments.
8. Refactor `ends` setup and `process_tile` the same way.
9. Adjust scaled counting call sites to require `WindowScaling.output_idx` for BED-like output rows
   and to keep using `WindowScaling.window_idx` where they need source-window or bin positions.
10. Remove `SourceIndexedBedWindowSets`, `split_bed_windows_by_length_with_source_indices`, and the
   separate narrow/broad span variables once both commands use the new BED window collection.
11. Distill the lasting behavior into the relevant current spec after implementation. Keep this file
    as the migration plan.

## Tests To Add

Add focused tests for the shared tiered overlap code:

- Broad plus narrow nested windows returns the same candidate set as the generic finder.
- Grouped BED entries return the correct `output_idx` group index while preserving
  `source_window_idx`.
- Always-hit broad windows return overlap fraction `1.0`.
- Broad windows at the tile assignment envelope boundary are not incorrectly classified as
  always-hit.
- Pointer retirement in every tier uses `query_start.saturating_sub(look_back)`.
- Returned order is allowed to differ from the generic finder, and tests compare sorted signatures.

Add only targeted command-level tests where existing fixtures make this practical:

- One grouped BED path in `lengths` or `ends` proves that output rows use `Some(output_idx)`
  directly from overlaps.
- One scaled BED-like path proves that `WindowScaling.output_idx` survives the scaling helper and is
  required for output identity.
- Avoid a full plain/grouped/scaled matrix unless the shared overlap tests or the targeted command
  tests expose a gap. The cleanup should not force broad revalidation of counting behavior that did
  not change.

Verification commands after implementation:

```text
cargo check
cargo check --tests --features testing
cargo check --all-features
```

The project convention is that the user runs tests. The implementation should add the tests, but not
run them unless explicitly requested.

## Open Decisions

- Whether to add a middle tier now, for example `[0, 10_000)`, `[10_000, 100_000)`,
  `[100_000, +inf)`, or keep the initial two-tier split and leave middle tiers for measurement.
  - Answer: Keep two-tiered for now. Not clear what the thresholds would be in a mid-tier
- Whether the overlap row should be a generic `OverlappingWindow` with an explicit count-window
  identity or a BED-specific overlap row type converted only where shared scaling helpers require it.
  - Answer: Extend `OverlappingWindow` with optional `output_idx`. This is less churn because
    scaling helpers already consume `OverlappingWindows`, and `Option<u64>` avoids pretending that
    fixed-size/global/scaling-bin rows already know a BED-like output identity.
- Whether `source_windows` should duplicate interval coordinates already present in tiers, or whether
  tiers should store source positions and use indirect reads. The direct-entry tier is the simpler
  and faster hot path. The source list is still useful for fetch narrowing and count allocation.
  - Answer: Duplicate `BedWindowEntry` values in tiers.
- Whether always-hit classification should stay broadest-tier-only or scan every tier during context
  construction. Broadest-tier-only is the intended starting point because it avoids extra work for
  narrow windows and still fixes the real nested broad-window case.
  - ANSWER: broad only
