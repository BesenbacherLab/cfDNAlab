# Current State Review — Tile/Fetch/Window Refactor

**Review date**: 2026-04-02
**Commit reviewed**: 6f22e18c "refactors tile+fetch+refspan with more"
**Method**: Full manual read of every changed file; no grep shortcuts; derivation from first principles where relevant

---

## 1. What the commit does

The commit decomposes tile candidate-window selection from BAM fetch narrowing and introduces an explicit vocabulary for the two fundamentally different window-selection models used across commands.

### 1.1 New infrastructure in `src/shared/`

**`tiled_run.rs`** gains three additions:

- `span_bounds_without_cache(windows, left_bound, right_bound)` — private helper; performs the
  two-pointer scan underlying all span computations. Already existed; now also used by the two
  new public functions below.

- `candidate_window_span_for_tile_core_overlap(windows, tile)` — public, single-tile,
  model-named wrapper around `span_bounds_without_cache(core_start, core_end)`.
  Returns the indices of windows that could overlap the tile core.

- `candidate_window_span_for_tile_fragment_reach(windows, tile, left_reach_bp, right_reach_bp)` —
  public, single-tile, model-named wrapper around
  `span_bounds_without_cache(core_start - left_reach, core_end + right_reach)`.
  Returns indices of windows reachable from fragments owned by this tile under the
  "aligned fragment start in tile core" ownership rule.

**`window_fetch.rs`** gains:

- `BedFetchPolicy` enum with three variants:
  - `CoreOverlap` — fetch envelope derives from windows that truly overlap the core
  - `CandidateWindowExtent` — fetch envelope trusts the precomputed candidate span directly
  - `KeepTileFetch` — returns tile.fetch unchanged

- `fetch_span_for_tile` gains a `bed_fetch_policy: BedFetchPolicy` parameter; BED mode now
  dispatches through this enum.

- `aligned_window_extent_for_core_overlap_bed(windows_chr, tile, candidate_span)` — re-applies
  `overlapping_windows_for_tile` to filter the cached span to true core-overlap windows before
  computing the fetch envelope.

- `aligned_window_extent_for_bed_candidates(windows_chr, tile_window_span)` — trusts the
  precomputed candidate span directly; no re-filtering.

- `full_tile_fetch_span(tile, chrom_len)` — returns tile.fetch clamped to chromosome length.

### 1.2 Changed commands

**`src/commands/ends/ends.rs`**:
- `precompute_tile_window_spans` now uses conditional halos based on `ClipStrategy`:
  - Raw: `left_halo = max_soft_clips`, `right_halo = max_fragment_length + max_soft_clips`
  - Aligned/Drop: `left_halo = 0`, `right_halo = max_fragment_length`
- `bed_fetch_halo_bp` similarly conditional on raw mode
- `fetch_span_for_tile` called with `BedFetchPolicy::CandidateWindowExtent` for **all** clip
  strategies (raw, aligned, drop — no branching)

**`src/commands/lengths/lengths.rs`**:
- `precompute_tile_window_spans(0, max_fragment_length)` — no left reach
- `fetch_span_for_tile` with `BedFetchPolicy::CandidateWindowExtent`
- `process_tile` iterates the full candidate span; stores `None` for windows ending before
  `core_start` (these windows are in the span due to straddling the left boundary, but their
  count slot stays empty because a fragment starting in the core cannot reach them)

**`src/commands/gc_bias/windows.rs`**:
- `prepare_tile_windows` (BED mode) now iterates `span.first_idx..span.last_idx_exclusive`,
  using the full fragment-reach candidate span. Previously was restricted to core-overlap windows.
- `precompute_tile_window_spans(0, max_fragment_length)` — same as `lengths`

**`src/commands/ref_gc_bias/ref_gc_bias.rs`**:
- Explicit `precompute_tile_window_spans(0, 0)` (pure core-overlap candidate bounds)
- `process_tile` calls `overlapping_windows_for_tile` for the precise per-window filter

### 1.3 New test files

- `src/shared/tiled_run_tests.rs` (403 lines, 15 tests for `precompute_tile_window_spans`,
  `candidate_window_span_for_tile_core_overlap`, `candidate_window_span_for_tile_fragment_reach`)
- `src/shared/window_fetch_tests.rs` (309 lines, tests for all three new fetch helpers)
- `src/commands/ref_gc_bias/ref_gc_bias_tests.rs` (145 lines)
- `src/commands/fragment_kmers/tiling_tests.rs` (103 lines)

### 1.4 Key semantic contracts

The `precompute_tile_window_spans` streaming function computes per-tile candidate bounds using
the same boundary logic as `span_bounds_without_cache`:
- Left: drops windows where `end <= left_bound = core_start - left_halo`
- Right: includes windows where `start < right_bound = core_end + right_halo`

This is implemented **inline** in `precompute_tile_window_spans`; it does **not** call
`span_bounds_without_cache` internally. The two standalone single-tile helpers DO call
`span_bounds_without_cache`. They are parallel implementations of identical boundary logic.

---

## 2. Architecture issues

### 2.1 `BedFetchPolicy::CoreOverlap` and `BedFetchPolicy::KeepTileFetch` are dead production code

**What the code does**: Both variants exist in the enum and are dispatched in `fetch_span_for_tile`.
The `CoreOverlap` branch calls `aligned_window_extent_for_core_overlap_bed`.

**What the code never does**: No command passes either of these variants to `fetch_span_for_tile`.
- `src/commands/ends/ends.rs:536` uses `CandidateWindowExtent`
- `src/commands/lengths/lengths.rs:525` uses `CandidateWindowExtent`
- All other commands (`fragment_kmers`, `fcoverage`, `midpoints`, `wps`) do NOT call
  `fetch_span_for_tile` at all; they use their own code paths (`tile_window_min_max`,
  `adapt_fetch_to_extreme_windows`, `get_overlapping_sites_and_adapt_fetch_to_extremes`)

**Implication**: `aligned_window_extent_for_core_overlap_bed` is tested in `window_fetch_tests.rs`
but is unreachable from any production code path. The tests verify a real contract, but that
contract is not enforced in practice for any current command.

**Risk**: If a new CoreOverlap command is written and its author sees `BedFetchPolicy::CoreOverlap`
exists, they might correctly reach for it. But if an existing CoreOverlap command is refactored
to use `fetch_span_for_tile`, it would need to choose `CoreOverlap`, and the existing tests of
`aligned_window_extent_for_core_overlap_bed` do cover this path. The real risk is that the enum
implies these variants are currently used, when they are not.

### 2.2 `candidate_window_span_for_tile_core_overlap` and `candidate_window_span_for_tile_fragment_reach` are tested but never called from production code

**What the code does**: Two public single-tile functions; both call `span_bounds_without_cache`.
Tested by 17 tests in `tiled_run_tests.rs`.

**What the code never does**: Neither function is called from any command or from
`precompute_tile_window_spans`. The streaming precompute function implements identical boundary
logic inline.

**Implication**: The tests for these helpers verify `span_bounds_without_cache` behavior, which IS
used by `overlapping_windows_for_tile` (the CoreOverlap per-window filter). If `span_bounds_without_cache`
has a regression, those tests will catch it. But they will NOT catch a boundary-condition regression in
`precompute_tile_window_spans` itself, because that function does not call `span_bounds_without_cache`.

The production coverage for `precompute_tile_window_spans` comes from the six direct tests in
`tiled_run_tests.rs` and multiple tests in `test_tiling.rs`.

**Verdict**: The design intent appears to be "named, testable artifacts that document each model's
contract." This is valuable. The naming (`CoreOverlap`, `ReachableFromTileOwnedFragments`) makes
implicit halo-encoding explicit. The risk is that they look like building blocks that are wired in
production but are not.

### 2.3 Comment in `precompute_tile_window_spans` misleads about who does the precise check

Lines 191–192 say: _"The iterator helpers perform the precise overlap check (end > core_start &&
start < core_end) so the recorded range only needs to bound the candidate windows."_

This comment describes the CoreOverlap flow correctly: `overlapping_windows_for_tile` (a CoreOverlap
iterator) re-checks each window against the core.

But for fragment-reach commands (`lengths`, `ends`, `gc_bias`), the precomputed candidate span is
trusted directly — no second filter is applied. The comment misleads readers into thinking the span
is always a bounding box that gets re-filtered later. It is not.

---

## 3. Documentation inaccuracies

### 3.1 Workdoc contradicts both the rules spec and the implementation

`TILE_WINDOW_AND_FETCH_EXECUTION_WORKDOC_2026-04-01.md` marks the task "Wire `ends` to
aligned-vs-raw BED fetch logic" as completed with this justification:

> _"`ends` must branch between aligned fragment-reach narrowing and raw full-tile fetch policy"_

The rules spec says the production contract for raw mode is `CandidateWindowExtent` with a
conservative halo (`max_fragment_length + max_soft_clips`), NOT `KeepTileFetch`. The
implementation in `ends.rs` uses `CandidateWindowExtent` for all clip strategies — no branch.

The workdoc describes a design that was evaluated and rejected or superseded. The rules spec
documents the chosen design correctly. The workdoc is stale.

### 3.2 Most new tests carry "Human verification status: unverified"

Of 15 tests in `tiled_run_tests.rs`, only one (`precompute_tile_window_spans_preserves_or_extends_span_in_example_when_right_reach_grows`)
is marked "Verified". All tests in `window_fetch_tests.rs` appear to have no verification label
(neither verified nor explicitly unverified).

None of the derivations appear incorrect upon inspection, but none have been independently
confirmed either. For tests that pin critical boundary conditions (e.g., exclusive-bound
semantics), human verification is the spec's own stated requirement.

---

## 4. Test coverage gaps relative to the test spec

The test spec (`TILE_WINDOW_AND_FETCH_TEST_SPEC_2026-04-01.md`) requires specific test coverage
at five layers. The following gaps were identified.

### 4.1 Layer 2A — Missing leftward monotonicity test for `aligned_window_extent_for_core_overlap_bed`

The spec requires: _"Adding another core-overlap window to the left widens or preserves fetch."_

The test `aligned_window_extent_for_core_overlap_bed_widens_monotonically_when_more_core_windows_are_added`
tests only **rightward** widening. It adds window [18,19) to window [10,11) inside core [10,20) and
checks that the right fetch bound widens. There is no test for adding a window to the left of an
existing core window.

**Why it matters**: The same implementation path applies to both directions. The spec called out
both directions explicitly. A leftward-narrowing bug would not be caught.

### 4.2 Layer 4 — No dedicated `fcoverage` halo-only BED window guard test

The spec requires for `fcoverage`: _"Halo-only BED window is ignored for BED outputs."_

The new test `by_bed_total_mixed_core_and_downstream_windows_is_tile_size_invariant` has a
window [10,11) that produces zero counts, but for the wrong reason: the fragment is [19,29) which
simply does not cover position 10. It is not testing that the CoreOverlap model **blocked** a
halo-only window from being processed by a tile it falls outside of.

A true guard test for `fcoverage` would require a fragment that **does** physically cover a
halo-only window position, where the window stays zero because the CoreOverlap model prevents the
owning tile from assigning coverage to it (while the downstream tile correctly gives it a count).

Without this test, a regression that switches `fcoverage` to a fragment-reach model would cause
double-counting but would not be detected by the existing test.

**This is a double-counting protection gap.** For `fcoverage`, using the wrong model causes
incorrect output, not just inefficiency. The current test does not cover this.

### 4.3 Layer 4 — `midpoints` halo-only tests are correct but only at the command output level

The test `bed_sites_mixed_core_and_halo_rows_keep_only_the_core_midpoint_count_across_tile_sizes`
correctly shows the halo-only site [22,23) stays zero for a midpoint at 10. This is model-correct
and tile-size-invariant.

However, there is no test that verifies a tile with a halo-only site in its cached span does not
**widen its fetch** to include that halo-only site. The spec requires: _"Mixed cached span with
one core window and one halo-only window narrows fetch from the core window only."_ That is a
fetch-narrowing test, not a count test. The command-level test only catches the counting error,
not the fetch-over-fetching error.

The `fragment_kmers` tiling_tests.rs does test this at the `determine_fetch_span` helper level,
but `midpoints` uses a different code path (`get_overlapping_sites_and_adapt_fetch_to_extremes`).
There is no fetch-narrowing test for `midpoints` specifically.

### 4.4 Layer 5 — No BED vs fixed-size parity tests for `ends` (aligned/drop mode)

The spec requires: _"Equivalent BED and fixed-size windows produce the same motif counts when
their windows describe the same genomic intervals."_

The new `ends` tests cover:
- BED mode with right-halo and far-right raw windows (Layer 4)
- Tile-size invariance within BED mode (Layer 4 correctness cross-checks)
- Fixed-size mode individually (Layer 2D/4D)

But there is no test that directly compares BED mode output to fixed-size mode output for the
same genomic intervals under aligned or drop clip strategy.

### 4.5 Layer 5 — No BED vs fixed-size parity tests for `gc_bias`

Same gap as 4.4 but for `gc_bias`. The spec requires this comparison. The
`test_gc_bias_windows.rs` tests `prepare_tile_windows` in isolation but does not compare the
resulting output to a fixed-size run.

### 4.6 Layer 2D — Chromosome-edge cases for CoreOverlap are not explicitly tested

The spec requires: _"Chromosome-edge case: first and last tiles still obey the same half-open
rules."_ The `precompute_tile_window_spans` tests use tiles with halos that do not hit chromosome
boundaries. The `build_tiles_respects_chrom_end_and_halo` test in `test_tiling.rs` checks tile
construction near a chromosome edge but does not check window span selection near the edge.

---

## 5. Test correctness checks (manual derivation)

The following critical boundary tests were independently derived to confirm correctness.

### 5.1 `precompute_tile_window_spans_excludes_windows_too_far_left_for_raw_end_reach`

Tile core [10,20), `left_halo=2`, so `left_bound = 10 - 2 = 8`.
Window [7,8): `end = 8 <= left_bound = 8` → **dropped** (test expects drop) ✓
Window [8,9): `end = 9 > left_bound = 8` → **kept** (test expects keep) ✓

### 5.2 `precompute_tile_window_spans_excludes_windows_too_far_right_for_raw_end_reach`

Tile core [10,20), `right_halo=14`, so `right_bound = 20 + 14 = 34`.
Window [33,34): `start = 33 < right_bound = 34` → **kept** (test expects keep) ✓
Window [34,35): `start = 34 >= right_bound = 34` → **dropped** (test expects drop) ✓

### 5.3 `candidate_window_span_for_tile_core_overlap_drops_window_ending_at_core_start`

From reading the test setup: window ending exactly at core_start.
`span_bounds_without_cache` left filter: `end <= core_start` → **dropped** ✓
This is the exact half-open boundary the spec requires to be tested.

### 5.4 `raw_endpoint_assignment_does_not_count_a_window_ending_at_the_left_raw_boundary`

Fragment: 2S10M2S at pos 10 (aligned start=10, raw left endpoint = 10-2 = 8).
Window [7,8): `end = 8`. Raw left endpoint at 8. Endpoint assignment is: left terminal base is
`boundary_pos = 8`; window must contain `boundary_pos` for left endpoint to count.
Half-open window [7,8) contains position 8? No — [7,8) contains {7}. Position 8 is NOT in [7,8).
So zero counts. Test expects 0. ✓

### 5.5 `raw_endpoint_assignment_keeps_a_left_window_that_only_raw_reach_can_touch`

Fragment: 2S10M2S at pos 10, aligned start=10, raw left endpoint = 8.
Window [8,9): left terminal base is 8. Is 8 in [8,9)? Yes. Test expects count `_T`. ✓

### 5.6 Fetch halo for raw `ends` BED mode

With `bed_fetch_halo_bp = max_fragment_length + max_soft_clips`:

For a far-right raw window at position `W_start`, the farthest-left fragment whose raw right
endpoint can reach `W_start` has:
- raw right endpoint = `W_start` → `aligned_end + right_soft_clips = W_start`
- worst case: `right_soft_clips = max_soft_clips` → `aligned_end = W_start - max_soft_clips`
- `aligned_start = aligned_end - frag_len >= W_start - max_soft_clips - max_fragment_length`

The fetch left bound is `window_span.start - (max_fragment_length + max_soft_clips)`. Taking
`window_span.start <= W_start`:
`fetch_left <= W_start - max_fragment_length - max_soft_clips = aligned_start_minimum`

The halo is sufficient. ✓

### 5.7 `clamp_fetch_respects_halo_and_chrom` in test_tiling.rs

Tile: fetch [80,170). Halo=20 (inferred from fetch width: (100-20)=20 on left, (170-150)=20 right).
Window span [90,160). Left fetch: 90-20=70, clamped to fetch_start=80 → 80. Right fetch: 160+20=180, clamped to fetch_end=170 and chrom_len=155 → 155.
Expected: [80,155). Test expects [80,155). ✓

But wait — `clamp_fetch_to_window_span` receives explicit `halo_bp=0` in this test (fourth arg).
Looking at the test more carefully: `clamp_fetch_to_window_span(&tile, 155, window_span, 0)`.
With `halo_bp=0`, the left fetch is `window_span.start - 0 = 90`, clamped to `tile.fetch_start = 80` → 80.
Right fetch: `window_span.end + 0 = 160`, clamped to `min(tile.fetch_end=170, chrom_len=155) = 155`.
Result: [80,155). ✓ (The tile carries the halo via tile.fetch, not via explicit halo_bp.)

### 5.8 `midpoint_fetch_span_keeps_fragment_start_that_old_symmetric_halo_would_drop`

Tile: core [80,95), fetch [70,95). Last tile with chromosome end at 95 clips the right halo to 0.
Window [90,95). Required fragment start at 84. With correct per-side halo tracking, left fetch = 80.
Old symmetric-halo bug: `total_extra = 10 → per_side = 5 → fetch_start = 90-5 = 85 > 84`. Bug confirmed.
Current implementation gives [80,95). ✓

---

## 6. `ref_gc_bias` test isolation concern

`process_tile_skips_a_halo_only_bed_window_in_core_overlap_mode` manually constructs a
`TileWindowSpan` that includes window [22,32) for tile core [10,20). In production,
`precompute_tile_window_spans(0, 0)` would never produce this span: window [22,32) has
`start = 22 >= core_end = 20`, so it would be excluded by the `start < right_bound = 20` check.

The test is testing that `overlapping_windows_for_tile` correctly filters the artificially wide
span. This is a valid unit test for `overlapping_windows_for_tile`. But it does not test the
production path, where `precompute_tile_window_spans(0, 0)` would never include this window.

The concern is narrow: someone reading the test might assume `process_tile` can receive such a
span from the production precompute, which it cannot. The test would be clearer if it noted that
the wide span is artificial.

---

## 7. Summary table

| # | Category | File | Severity | Description |
|---|----------|------|----------|-------------|
| B.1 | Dead code | `window_fetch.rs` | Low | `BedFetchPolicy::CoreOverlap` and `KeepTileFetch` never passed by any command |
| B.2 | Dead code | `window_fetch.rs` | Low | `aligned_window_extent_for_core_overlap_bed` unreachable from production |
| B.3 | Dead code | `tiled_run.rs` | Low | `candidate_window_span_for_tile_*` tested but not in production path |
| D.1 | Documentation | Workdoc | Low | "ends must branch" description contradicts implementation and rules spec |
| D.2 | Documentation | `tiled_run.rs` comment | Low | "iterator helpers perform precise check" only true for CoreOverlap commands |
| D.3 | Documentation | `tiled_run_tests.rs` | Medium | 14/15 tests unverified despite spec requiring derivation |
| G.1 | Test gap | `window_fetch_tests.rs` | Low | Missing leftward monotonicity test for `aligned_window_extent_for_core_overlap_bed` |
| G.2 | Test gap | `test_fcoverage_command.rs` | **High** | No halo-only-window double-counting prevention test; wrong model would not be caught |
| G.3 | Test gap | `test_profile_groups_command.rs` | Low | No fetch-narrowing test for midpoints; only count correctness tested |
| G.4 | Test gap | `test_ends_command.rs` | Low | No BED vs fixed-size parity comparison for aligned/drop mode |
| G.5 | Test gap | `test_gc_bias_windows.rs` | Low | No BED vs fixed-size parity comparison for gc_bias |
| G.6 | Test gap | `tiled_run_tests.rs` | Low | No chromosome-edge window-selection test for CoreOverlap |
| T.1 | Test isolation | `ref_gc_bias_tests.rs` | Info | Test uses artificially wide span that production code never produces |

---

## 8. What the commit gets right

The following are correct and well-executed:

1. **`precompute_tile_window_spans` boundary logic is correct.** The inline streaming logic
   implements exactly `end > left_bound && start < right_bound`, which is the correct half-open
   overlap test for both the CoreOverlap model (halos=0) and the fragment-reach model
   (right_halo > 0).

2. **`ends` raw-mode halo derivation is correct.** `tile_span_left_halo = max_soft_clips` and
   `tile_span_right_halo = max_fragment_length + max_soft_clips` cover all reachable endpoints.
   `bed_fetch_halo_bp = max_fragment_length + max_soft_clips` is sufficient to reconstruct every
   eligible read (see §5.6 derivation).

3. **No branching on ClipStrategy for the fetch policy.** Using `CandidateWindowExtent` for all
   modes is simpler and conservative. The conservative halo for raw mode guarantees no reads are
   dropped. This matches the rules spec intent.

4. **`gc_bias` BED preparation now correctly uses the full fragment-reach span.** Before this
   commit, the production code for `gc_bias` BED mode was likely using the CoreOverlap span,
   which would have missed fragment contributions to right-halo-only windows. The commit fixes
   this by using `precompute_tile_window_spans(0, max_fragment_length)`.

5. **`lengths` correctly stores `None` for windows ending before `core_start`.** Even though these
   windows are in the candidate span (they straddle the left core boundary), a fragment with
   aligned start in the core cannot have its right endpoint inside such a window. Storing `None`
   for them prevents false counts.

6. **The `clamp_fetch_to_window_span` asymmetric-halo fix is correct.** The new code keeps the
   per-side halo from the tile rather than inferring a symmetric halo from the total extra fetched
   bases. The old symmetric-halo bug is documented and tested in `test_tiling.rs`.

7. **`fragment_kmers` tests verify CoreOverlap isolation.** The new tests in `tiling_tests.rs`
   show that `determine_fetch_span` returns `None` for halo-only windows and narrows fetch from
   core-only windows when a mixed span is present.

8. **`midpoints` CoreOverlap tests work correctly.** The `core_overlap_bed_site_is_kept_for_midpoints`
   and `bed_sites_mixed_core_and_halo_rows_keep_only_the_core_midpoint_count_across_tile_sizes`
   tests correctly verify that halo-only sites receive zero count and output is tile-size-invariant.

9. **`ref_gc_bias` correctly uses CoreOverlap.** `precompute_tile_window_spans(0, 0)` plus
   `overlapping_windows_for_tile` provides double filtering: the span bounding is CoreOverlap, and
   the per-window iterator further applies the precise `end > core_start && start < core_end` check.

---

## 9. Most critical actions

Listed by priority:

### P1 — Add `fcoverage` CoreOverlap double-counting guard test

The absence of a test where a fragment physically covers a halo-only window (which must stay zero
because CoreOverlap blocks the owning tile from counting it) leaves the most consequential failure
mode untested. A regression here produces silently incorrect output.

Required test: fragment spanning [9,23), windows [10,11) and [22,23), tile_size small enough that
the fragment's owning tile is NOT the same as the window [22,23)'s tile. With CoreOverlap: both
windows count once (from their respective owning tiles). With a fragment-reach model: [22,23)
would be double-counted.

### P2 — Mark dead code explicitly

All three dead code clusters (§2.1 and §2.2) should carry a note in their docstring stating they
are not in the current production path. This prevents future readers from assuming these tested
functions are exercised by commands.

### P3 — Verify all "unverified" test derivations

The test spec requires hand-derivation for every boundary case. The unverified status on 14 tests
in `tiled_run_tests.rs` means the most important boundary conditions (exact equality drops,
exclusive bounds, asymmetric halos) have not been confirmed by a human reader.

### P4 — Fix the workdoc task description (§3.1)

Remove or correct the "ends must branch between aligned narrowing and raw KeepTileFetch"
justification. Replace with the actual production contract: `CandidateWindowExtent` with
conservative halo covers both modes.

### P5 — Fix the `precompute_tile_window_spans` comment (§2.3)

Qualify the comment to say it applies to CoreOverlap commands. For fragment-reach commands, the
span IS trusted directly.

### P6 — Add BED vs fixed-size parity tests for `ends` (aligned/drop) and `gc_bias` (Layer 5)

These are required by the test spec and constitute the only end-to-end verification that the two
windowing modes agree.

---

## 10. Notes on test files not primarily related to the core refactor

**`src/commands/ends/motifs_tests.rs`**: Tests `count_fragment_in_window` (endpoint counting),
`encode_inside_code`, `encode_outside_code`, `motif_reference_span_for_tile`. These cover
edge-of-window boundary semantics (right terminal base at `boundary_pos - 1`, left at
`boundary_pos`) and reference-span expansion for raw/aligned modes. All look correct.

**`src/commands/ends/write_tests.rs`**: Tests dense matrix stacking and JSON settings sidecar.
Correct.

**`src/shared/fragment/ends_fragment.rs`**: Assignment boundary computation for raw mode
(`assignment_boundary_pos = aligned_start - left_soft_clip_bp` for left end). Correct and
consistent with how window assignment uses this field.

**`src/shared/fragment_iterators/ends.rs`**: Thin adapter; no logic change.

**`tests/test_tiling.rs` additions**: Layer 3 clamp tests are well-constructed. The
`clamp_fetch_uses_explicit_halo_when_tile_has_no_inferred_halo` test pinpoints the old symmetric-halo
bug precisely. The `adapt_fetch_keeps_left_halo_when_chrom_end_clips_last_tile_in_fcoverage` test is
the most realistic regression guard for the chromosome-end asymmetric-halo issue.
