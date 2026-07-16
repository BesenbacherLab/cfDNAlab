# Lengths Overlap Finder Heuristics Plan

## Purpose

This note tracks a future optimization for `cfdna lengths` window assignment when a grouped BED
mixes ordinary local intervals with broad intervals, for example a chromosome-wide background
group plus short cell-type intervals.

The current start-sorted overlap scan is fast for the common case. It should stay the default
model. The problem is a special nested-window shape where a broad early interval prevents the
moving window pointer from advancing, causing repeated scans across local windows that cannot
overlap the current fragment.

The goal is not to replace the overlap finder with a general interval tree. The goal is a small
tile-local heuristic that keeps the common path simple and handles broad-window nesting cleanly.

## Current Constraint

`lengths` calls the shared BED overlap finder once per counted fragment. The finder keeps a moving
window pointer and advances it only past windows where:

```text
window.end <= query_start.saturating_sub(look_back)
```

That `look_back` bound is required for correctness. Fragment assignment intervals in a tile are
not guaranteed to arrive in strictly increasing query-start order, and a later query can still
need a window that ended before the current query start but within the configured reach.

Any optimized path must preserve the same retirement rule unless it proves a stronger ordering
invariant.

## Proposed Design

Build a tile-local overlap context for BED and grouped BED windows. The context classifies the
tile candidate windows into three categories:

- `always_hit`: windows that contain the full tile assignment envelope.
- `broad`: long windows that should not participate in the normal local-window pointer scan.
- `local`: all remaining windows, scanned with the existing pointer logic.

For `lengths`, the tile assignment envelope is:

```text
[tile.core_start - left_assignment_reach, tile.core_end + max_fragment_reach)
```

A window in `always_hit` is valid for every fragment owned by that tile. It can be added to every
fragment's overlap result with overlap fraction `1.0`.

Broad windows use a separate active sweep. They are not exact-checked as a full broad list for
every fragment. The broad active set should:

- add broad windows while `window.start < query_end`
- retire broad windows only when `window.end <= query_start.saturating_sub(look_back)`
- exact-check only the currently active broad windows against the current query interval

Local windows keep the existing moving-pointer scan, but broad and always-hit windows are removed
from that local stream so they cannot pin the pointer.

The per-fragment result is the union of:

```text
always_hit overlaps
active broad overlaps
local overlaps
```

Returned window indices must still refer to the original chromosome-local window slice, because
grouped BED counting maps overlap rows back through `IndexedInterval.idx()`.

## Broad Threshold

Use a fixed internal threshold. Relative thresholds based on fragment reach or window-length
distribution add unpredictability without encoding the actual cost model clearly.

Initial candidate:

```text
BROAD_WINDOW_MIN_BP = 500_000
```

Classification:

```text
if window contains tile_assignment_envelope:
    always_hit
else if window.len() >= BROAD_WINDOW_MIN_BP:
    broad
else:
    local
```

This should catch chromosome-scale background intervals and 1 Mb windows mixed with small feature
windows, while leaving ordinary 100 kb-scale windows on the current fast path.

If this turns out too aggressive, raise the threshold to 1 Mb. If it misses real broad-window
pathologies, lower it. Avoid adding a user-facing CLI knob unless measurements show that a single
internal threshold is not stable enough.

## Expected Behavior

Common non-nested BED files should behave like today. The only added work is a tile-local
classification pass over candidate windows and small per-fragment handling for active broad
windows when they exist.

In a grouped BED with one background row plus many local intervals:

- the background row becomes `always_hit`
- local intervals stay in the existing pointer scan
- the background row no longer pins the local pointer

In a BED with 1 Mb windows mixed with short intervals:

- 1 Mb windows become broad
- broad windows enter and leave a separate active set
- short intervals stay in the existing pointer scan

In a BED where most windows are broad, this heuristic may not help much. It should still be
correct, and the active broad count should make the performance cost visible if instrumentation is
added later.

## Implementation Notes

Keep the implementation local to the overlap/caller boundary. Downstream counting should still
consume an `OverlappingWindows`-style result and should not need to know whether a hit came from
`always_hit`, `broad`, or `local`.

Possible shape:

```text
TileWindowOverlapContext
  always_hit: Vec<IndexedInterval>
  broad: Vec<IndexedInterval>
  local: Vec<IndexedInterval>
  local_to_original_idx: Vec<usize>
  broad_active: Vec<usize>
  broad_start_ptr: usize
  local_wd_ptr: usize
```

The exact data structure can be simpler if local windows retain original indices directly. The
important invariant is that overlap rows keep the original chromosome-local index used by the
loaded window slice.

Before implementing, decide whether this should be:

- a `lengths`-specific helper first, because `lengths` has the fragment-reach and tile-envelope
  context needed for `always_hit`
- or a shared BED overlap helper with command-specific envelope inputs

Start with `lengths` unless another command has the same pathological broad/local grouped BED use
case.

## Testing Targets

Use focused internal tests for the overlap context, not end-to-end performance tests.

Cases to cover:

- a chromosome-wide background interval plus short local windows returns the same overlaps as the
  current finder
- `always_hit` windows produce overlap fraction `1.0` for all assignment modes that consume the
  returned overlap fraction
- broad windows at tile-envelope boundaries are not classified as `always_hit`
- active broad retirement uses `query_start.saturating_sub(look_back)`, not `query_start`
- returned indices map back to the original grouped BED window slice
- local pointer behavior matches the current finder after broad windows are removed

If performance instrumentation is added, keep it internal:

```text
always-hit windows per tile
broad windows per tile
max active broad windows per tile
local windows scanned per fragment
```

