# Ref-kmers Tiered Fixed-Width Overlap Plan

## Purpose

`ref-kmers` currently uses `FixedWidthOverlapCursor` for repeated k-mer window lookups. The cursor
keeps the useful fixed-width cache, but in BED mode it still owns one forward pointer over one
start-sorted BED slice. A chromosome-wide background row can keep that pointer pinned near the
chromosome start, so every cache rebuild can scan many retired nested windows to the left of the
current k-mer.

`lengths` and `ends` already solve this class of BED input with tiered BED windows. This plan
describes how to make tier-aware BED lookup the normal BED implementation inside
`FixedWidthOverlapCursor`, without reimplementing the existing tier model or the fixed-width cache
math.

## Invariants

- Keep `ref-kmers` count values unchanged for global, fixed-size, BED, and grouped BED windows.
- Keep BED row identity unchanged. Plain BED rows use original BED row ids. Grouped BED rows use
  group indices.
- Keep `OverlappingWindow.idx` as the chromosome-local source-window index for BED-like rows.
- Keep `OverlappingWindow.output_idx` as the optional downstream BED-like row id when the normal
  BED path already knows it.
- Keep `count-overlap` as `overlap_bp / kmer_size`.
- Keep `any`, `all`, `midpoint`, and `proportion=<threshold>` threshold behavior unchanged.
- Keep k-mers that cross `chrom_len` skipped before overlap lookup. Do not clip them into counts.
- Do not change global or fixed-size lookup paths.
- Do not keep a separate non-tiered BED implementation inside `FixedWidthOverlapCursor`.
- Do not add a `new_tiered_bed` constructor or any other opt-in path that treats tiered BED lookup
  as a special case.
- Do not add a user-facing tier threshold option.

## Existing Code To Reuse

Reuse these existing pieces from `src/shared/overlaps.rs`:

- `DEFAULT_BROAD_WINDOW_MIN_BP`
- `BedWindowEntry`
- `BedWindowTier`
- `ChromosomeBedWindows`
- `TileBedWindowSpans`
- `TileBedWindowInputs`
- `build_bed_windows_by_chr`
- `precompute_tile_bed_window_spans`
- `span_slice`
- `clamp_bed_window_to_chrom`
- `CachedFixedWidthWindow` range math
- `required_overlap_bases`
- `first_query_start_reaching_window`
- `later_query_start_or_never`

Do not add a second ref-kmers-only representation for broad and narrow windows.

## Design

### Shared BED Tile Split

Extract the shared part of `TileBedOverlapContext::new` into a helper that receives:

```rust
chrom_len: u64
chromosome_windows: &ChromosomeBedWindows
spans: &TileBedWindowSpans
tile_assignment_envelope: Interval<u64>
```

It should return:

```rust
struct TileBedWindowSplit {
    always_hit_windows: Vec<BedWindowEntry>,
    scanned_tiers: Vec<Vec<BedWindowEntry>>,
}
```

The helper should do the current work exactly once:

1. Slice each tier with `span_slice`.
2. Clamp each candidate with `clamp_bed_window_to_chrom`.
3. Move windows containing the full tile envelope into `always_hit_windows`.
4. Keep every other candidate in the tier it came from.

Then `TileBedOverlapContext::new` and `FixedWidthOverlapCursor` both use this helper.

### Fixed-Width Cache Object

Move the current BED-cache state out of `FixedWidthOverlapCursor` into a small internal object:

```rust
struct FixedWidthBedCache {
    wd_ptr: usize,
    cache_ready: bool,
    cached_windows: Vec<CachedFixedWidthWindow>,
    next_candidate_change_query_start: u64,
    refresh_count: usize,
}
```

BED mode keeps one `FixedWidthBedCache` per scanned tier. A broad chromosome-wide row can then live
in `always_hit_windows`, while narrow tiers advance their own pointers without being pinned behind
the broad row.

Each size tier must own independent cache state:

- its own `wd_ptr`
- its own `cached_windows`
- its own `next_candidate_change_query_start`
- its own refresh and inspection counters for test-only instrumentation

This is the point of carrying tiers into `FixedWidthOverlapCursor`. A cache rebuild in a narrow tier
must not scan broad rows, and a long broad row must not keep the narrow tier's pointer near the
chromosome start.

The owned shape should make that explicit:

```rust
struct FixedWidthBedTier {
    windows: Vec<BedWindowEntry>,
    cache: FixedWidthBedCache,
}

enum FixedWidthWindowSource<'a> {
    Global,
    FixedSize(u64),
    Bed {
        always_hit_windows: Vec<BedWindowEntry>,
        tiers: Vec<FixedWidthBedTier>,
        all_windows: &'a [BedWindowEntry],
    },
}
```

The BED source builder should turn each scanned size tier into a separate `FixedWidthBedTier`:

```rust
let tiers: Vec<FixedWidthBedTier> = split
    .scanned_tiers
    .into_iter()
    .map(|windows| FixedWidthBedTier {
        windows,
        cache: FixedWidthBedCache::new(),
    })
    .collect();
```

No cache field should be shared between entries in `tiers`.

### Cache Refresh Input

Generalize cache refresh over a small private trait or helper input type that exposes the fields
needed by the cache:

```rust
fn start(&self) -> u64;
fn end(&self) -> u64;
fn result_idx(&self) -> usize;
fn output_idx(&self) -> Option<u64>;
```

Implement the input behavior for `BedWindowEntry`. This avoids copying the current refresh logic
while making row identity explicit:

- `BedWindowEntry` returns `all_windows_idx` as `idx` and `Some(output_idx)`.

The cursor should not need an `IndexedInterval<u64>` BED refresh path after this refactor. Incoming
BED and grouped BED windows are converted into `ChromosomeBedWindows` before tile processing.

### Cached Window Identity

Change `CachedFixedWidthWindow::new` so it no longer requires `IndexedInterval<u64>`. It should
accept the pieces the cache actually needs:

```rust
idx: usize
output_idx: Option<u64>
window_start: u64
window_end: u64
chrom_len: u64
query_width: u64
required_overlap: u64
```

Keep the existing accepted-range and constant-overlap calculations. Add `output_idx` to
`CachedFixedWidthWindow`, then return it through `OverlappingWindow::new_with_output_idx`.

### Cursor Modes

Replace the current optional `windows` BED state with the `FixedWidthWindowSource` shape above.
The key point is that `Bed` owns `Vec<FixedWidthBedTier>`, not a single BED cache plus tier slices:

```rust
enum FixedWidthWindowSource<'a> {
    Global,
    FixedSize(u64),
    Bed {
        always_hit_windows: Vec<BedWindowEntry>,
        tiers: Vec<FixedWidthBedTier>,
        all_windows: &'a [BedWindowEntry],
    },
}

struct FixedWidthBedTier {
    windows: Vec<BedWindowEntry>,
    cache: FixedWidthBedCache,
}
```

Change the existing constructor instead of adding a special constructor. The shape should be closer
to:

```rust
FixedWidthOverlapCursor::new(
    chrom_len: u64,
    source: FixedWidthWindowSource<'a>,
    query_width: u64,
    min_overlap_fraction: f64,
) -> Result<Self>
```

Add a small source builder for BED-like modes if needed:

```rust
FixedWidthWindowSource::bed(
    chrom_len: u64,
    bed_window_inputs: TileBedWindowInputs<'a>,
    tile_query_envelope: Interval<u64>,
) -> Result<Self>
```

That builder should call the shared BED tile split helper and build one cache per scanned tier. It
is a source constructor, not a second cursor mode.

### Lookup Behavior

For BED lookups:

1. Build the checked query interval with the existing `query_interval(...)` helper.
2. For full-width queries, use the per-tier fixed-width caches.
3. For clipped queries, use the same tier inputs without the fixed-width cache and compute overlap
   fractions against `query_interval.len()`. Do not fall back to a separate non-tiered BED path.
   `ref-kmers` should not submit clipped k-mer spans because counting skips `kmer_end_abs >
   chrom_len`, but the cursor should still have one internally consistent BED implementation.
4. Add `always_hit_windows` to the result with overlap fraction `1.0`.
5. For each scanned tier, refresh only that tier's `FixedWidthBedCache` when its own
   `next_candidate_change_query_start` has passed.
6. Append cached hits from every tier.
7. Return `None` only if no always-hit or cached tier hit is present.

Result order is not a contract. Tests should compare order-insensitive signatures where tier order
can differ from the old single-slice finder.

## Ref-kmers Command Wiring

In `src/commands/ref_kmers/ref_kmers.rs`:

1. Build `bed_windows_by_chr` from `indexed_windows_map` using `build_bed_windows_by_chr` and
   `DEFAULT_BROAD_WINDOW_MIN_BP`.
2. Replace the current `precompute_tile_window_spans(...)` BED span with
   `precompute_tile_bed_window_spans(...)` for BED and grouped BED modes.
3. Keep non-BED modes on the existing fixed-size/global path.
4. In the tile loop, pass `TileBedWindowInputs` into `process_tile` for BED-like modes.
5. Do not pass `windows_chr` only to support old BED cursor behavior. If metadata or output setup
   still needs the original loaded maps, keep those at the outer command layer.

In `src/commands/ref_kmers/counting.rs`:

1. Add an optional `TileBedWindowInputs` argument to `count_kmers_by_window`.
2. Build one `FixedWidthWindowSource` for the current mode:
   - `Global` for global mode
   - `FixedSize(window_bp)` for fixed-size mode
   - `Bed { ... }` from `TileBedWindowInputs` for BED and grouped BED modes
3. Create the cursor with the single `FixedWidthOverlapCursor::new(...)` constructor.

## Tile Query Envelope

The BED source builder needs a tile envelope that contains every query submitted to the cursor.

For `count-overlap`, `any`, `all`, and `proportion`, query intervals are full k-mer spans. Use:

```text
[first_counted_kmer_start, last_counted_kmer_start + kmer_size)
```

For `midpoint`, query intervals are 1 bp midpoint spans. Use:

```text
[first_counted_midpoint, last_counted_midpoint + 1)
```

Derive these from the owned start range after applying `sequence_start`, `chrom_len`, `kmer_size`,
and assignment mode. If no countable start exists, return early before building the cursor.

## Test Sequence

Add characterization tests before implementation changes. These tests should compile and pass
against the current implementation first, then keep the same expected signatures and output values
through the refactor. After the cursor constructor changes, update only test setup code needed to
build the cursor source. Do not update the expected windows or counts to match the new code.

The tests must not use the current overlap implementation as the oracle. Expected windows per query
should be hand-derived from the fixture intervals. Command-level expected counts should be derived
from those same per-position assignments.

### Shared Overlap Tests

In `src/shared/overlaps_tests.rs`:

- Add a test that asserts exact per-query window signatures for the existing mixed broad/narrow
  fixture. The signatures should include `idx`, `output_idx`, window coordinates, and overlap
  fraction.
- Add a test with a chromosome-wide background window plus many narrow windows to the left of later
  query starts. It should assert exact per-query window signatures first, then use instrumentation
  after the refactor to prove tier caches do not keep rescanning retired narrow windows because a
  broad row is still active.
- Add threshold-specific cases for `any`, `all`, `proportion`, and `midpoint` where the same BED
  fixture can expose an off-by-one or wrong-denominator error.
- Keep expected overlap counts hand-derived. Do not generate expectations with another overlap
  implementation.

Use test-only instrumentation on `FixedWidthBedCache` if needed:

- refresh count per tier
- candidate windows inspected per refresh

Do not expose this outside tests.

### Ref-kmers Counting Tests

In `src/commands/ref_kmers/counting_tests.rs`:

- Add a grouped BED test that proves `all_windows_idx` still maps to the
  stored group index.
- Add a plain BED test with a background row plus nested rows. The expected row counts should be
  derived from explicit per-position assignments, not by comparing to the old cursor.

### Integration Tests

Add command-level regression coverage before the refactor in `tests/`:

- Use a tiny reference and BED fixture where every countable k-mer start can be enumerated in the
  test comments.
- Include a chromosome-wide background row and nested short rows.
- Include grouped BED if the public command can build the fixture cheaply.
- Assert final output rows and motif weights exactly enough to detect a changed window assignment
  per k-mer position.

The integration tests cannot observe the cursor directly, so they should document the per-position
assignment table in comments and assert the resulting row totals and motif weights.

## Verification

After code changes:

```text
cargo check
cargo check --tests --features testing,cmd_ref_kmers
```

Do not run tests unless explicitly requested.
