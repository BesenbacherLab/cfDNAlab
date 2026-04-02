# Tile, window, and fetch test spec

Date: 2026-04-01

## Scope

This spec turns [TILE_WINDOW_AND_FETCH_RULES_SPEC_2026-04-01.md](/Users/au547627/Documents/Development/rust/cfDNAlab/plans_and_specs/TILE_WINDOW_AND_FETCH_RULES_SPEC_2026-04-01.md) into a concrete test plan.

The goal is not "some useful regressions". The goal is to define the minimum test set that would
fail clearly if the new decomposition is implemented with the wrong semantics.

This spec exists so implementation does not drift into:

- changing helpers without pinning down their contracts
- testing only the already-known bug cases
- missing off-by-one errors at half-open boundaries
- missing cross-command model mismatches

## Test philosophy

The tests in this area must do all of the following:

- derive expectations by hand
- use the strongest boundary that still isolates the rule being tested
- fail loudly if the wrong ownership model is used
- fail loudly if half-open semantics are wrong
- fail loudly if candidate-window reach and fetch narrowing drift apart

Fixtures should stay small enough to reason about exactly, but they must not be so small that a
wrong implementation can still pass accidentally.

Every major rule below needs:

- at least one positive case
- the nearest off-by-one case
- at least one mixed case that distinguishes the intended model from an easy wrong model

Ownership assumptions must also be tested explicitly.

If a helper or command assumes "fragments starting in tile", the tests should say that directly.
They should not rely on the reader inferring it from the fixture.

If a future helper or command uses a different ownership rule, such as "fragments ending in tile",
that must get its own test block rather than silently reusing the "starts in tile" cases.

The corresponding helper docstrings should use the same vocabulary as the tests. A review in this
area is incomplete if the tests make the ownership rule explicit but the helper docstring still
leaves it implicit.

## Test layers

The tests should be built in this order:

1. Helper-level window-selection tests
2. Helper-level fetch-envelope tests
3. Final clamp tests
4. Command classification tests
5. BED vs fixed-size consistency tests where the command contract requires it

Do not start by adding only end-to-end command tests. That makes failures much harder to interpret.

## Layer 1: candidate window selection

This layer tests which windows are relevant for a tile before fetch narrowing is considered.

These tests should live with the helper that precomputes candidate BED spans and with any helper
that iterates windows under an explicit model.

### A. `CoreOverlap` model

These tests prove that the model means exactly "window overlaps tile core".

Required cases:

- Keep a window fully inside the core
- Keep a window crossing the left core boundary
- Keep a window crossing the right core boundary
- Drop a window ending exactly at `core_start`
- Drop a window starting exactly at `core_end`
- Mixed case: left halo-only, core-overlap, right halo-only
  - Expect only the core-overlap windows to survive
- Cached-span guard case:
  - input candidate span includes halo-only windows
  - core-overlap iterator still yields only true core-overlap windows
- Multi-tile case:
  - a boundary-crossing BED window appears in both neighboring tiles
  - a fully internal BED window appears only in its owning tile
- Chromosome-edge case:
  - first and last tiles still obey the same half-open rules

These tests should fail if the implementation starts trusting a wider candidate span directly in a
core-overlap command.

### B. `ReachableFromTileOwnedFragments` model, aligned interval geometry

These tests cover commands like `lengths` and aligned/drop `ends`.

Assumed ownership rule:

- fragment is owned by the tile iff aligned fragment start lies inside the tile core

Required cases:

- No left reach
  - a BED window left of the core is dropped
- Keep one right halo-only window that an owned fragment can reach
- Drop a right window starting exactly at the exclusive reachable bound
- Mixed case: core-overlap BED window plus right-halo target BED window
  - expect both candidate windows to be kept
- Mixed case: left non-overlap window plus core window plus right-halo target
  - expect left to be dropped and the other two kept
- Monotonicity case:
  - increasing right reach grows or preserves the candidate span
  - it never shrinks the candidate span or moves the first candidate rightward
- Multi-tile ownership case:
  - downstream window is selected only for the tile that owns fragments reaching it

These tests should fail if the implementation falls back to plain core overlap.

They should also fail if the implementation silently switches to a different ownership rule.

### C. `ReachableFromTileOwnedFragments` model, raw endpoint geometry

These tests cover `ends` with `clip_strategy = raw`.

Required cases:

- Keep a left halo-only window reachable only because raw clipping moves the left endpoint left
- Drop a left window ending exactly at the exclusive raw-left bound
- Keep a far-right window reachable only because raw clipping moves the right endpoint right
- Drop a far-right window starting exactly at the exclusive raw-right bound
- Asymmetry case:
  - left reach and right reach differ
  - implementation must not collapse them into one symmetric rule
- Mixed case:
  - left raw-only window, core-overlap window, and far-right raw-only window are all present
  - candidate selection keeps all three

These tests should fail if the implementation treats raw mode as aligned mode with one symmetric halo.

They should also fail if the implementation stops using the documented "aligned start in tile"
ownership rule.

## Layer 2: aligned fetch-envelope derivation

This layer tests the helper that turns relevant windows into an aligned-coordinate support envelope
before the final clamp to `tile.fetch`.

This layer must be model-aware. It must not rely on ambiguous cached spans.

### A. `CoreOverlap` fetch-envelope derivation

Required cases:

- One core-overlap BED window narrows fetch around that window
- Mixed cached span with halo-only windows still derives the envelope from core-overlap windows only
- No core-overlap windows returns `None`
- Adding another core-overlap window to the left widens or preserves fetch
- Adding another core-overlap window to the right widens or preserves fetch
- Mixed case:
  - one core-overlap window plus one halo-only candidate
  - fetch envelope is identical to the one-core-window case

These tests should fail if the implementation trusts all cached candidates instead of true
core-overlap windows.

### B. Fragment-reach fetch-envelope derivation with aligned counting geometry

This covers commands where count geometry and aligned geometry are the same, such as `lengths`.

Required cases:

- One right halo-only candidate BED window produces a non-empty fetch envelope
- Mixed case:
  - one core-overlap BED window plus one right-halo target BED window
  - fetch envelope uses the full candidate extent, not only the core window
- Coupling rule case:
  - if a window is kept because reach is `R`, the left fetch bound must still allow reads starting
    `R` bases earlier
- Off-by-one case:
  - a candidate window that only touches fragment end stays candidate-relevant if the selection
    model keeps it, but later counting tests must still stay zero
- Monotonicity case:
  - adding more relevant candidate windows widens or preserves fetch, never narrows it

These tests should fail if fetch narrowing uses a smaller reach value than candidate-window
selection used.

### C. Raw `ends` fetch policy

This layer is intentionally different.

The spec says raw endpoint windows are count-relevant in raw-shifted geometry, but production fetch
narrowing should still derive an aligned BAM support envelope from the raw candidate set.

So the tests here should assert the chosen policy explicitly.

The chosen policy is:

- raw BED mode narrows fetch from the raw candidate-window extent using an aligned halo large
  enough not to lose any eligible reads
- `KeepTileFetch` remains a safe fallback and may still be useful in tests, but it is not the
  target production contract

Required cases are:

- left raw-only BED window produces a narrowed fetch interval that still preserves enough aligned
  context
- far-right raw-only BED window produces a narrowed fetch interval on the right
- mixed core and far-right raw windows widen or preserve fetch relative to the smaller candidate
  set
- if tile already carries a larger fetch halo on one side, the final clamp keeps that larger tile
  halo
- adding more raw-only candidate windows never narrows the aligned fetch

### D. Fixed-size mode is a required explicit test area

This is a required explicit test area in general, not an optional follow-up for one command.

Reason:

- the repo has already had real `--by-size` bugs
- so BED-mode coverage is not sufficient evidence that fixed-size mode is correct

Required rule:

- if a command supports both BED and fixed-size windows, and its contract says those two window
  modes represent the same logical partition, fixed-size mode needs its own explicit tests

Minimum required cases by command family:

- fragment-reach commands:
  - a count that lands outside the owning tile core but inside the correct fixed-size bin is kept
  - a mixed case with one non-target bin and one true target bin counts only the true target
  - changing `tile_size` does not change the final fixed-size output rows
- core-overlap commands:
  - a non-target neighboring bin stays zero
  - changing `tile_size` does not create or remove final fixed-size rows
- half-open boundary cases:
  - exact-boundary fixed-size bins must stay zero when the command's BED-mode analogue would stay
    zero

Important subcase:

- raw `ends` still needs its own dedicated fixed-size block because raw endpoint reach is
  asymmetric

Required raw `ends` fixed-size cases:

- a raw left-reachable endpoint that lands in the previous fixed-size bin is counted
- a raw far-right endpoint that lands in the next fixed-size bin is counted
- a mixed case with one non-target core bin and one true raw-endpoint bin counts only the true
  target
- exact half-open boundary cases stay zero
- changing `tile_size` does not change the final fixed-size output rows

These tests should be written at the command level, not on the current mixed fetch helpers.

## Layer 3: final clamp helper

This layer is the generic aligned clamp logic. It should stay independent from the window-selection
model.

Required cases:

- Preserves tile-carried left halo near chromosome start
- Preserves tile-carried right halo near chromosome end
- Uses explicit halo when tile has no inferred halo on one or both sides
- Wider support interval never narrows the clamped fetch interval
- Left clamp to `tile.fetch_start`
- Right clamp to `tile.fetch_end`
- Clamp to chromosome end
- Empty-after-clamp returns `None`
- Exact-boundary cases at `fetch_start` and `fetch_end`

These tests should fail if the clamp helper starts inferring one combined symmetric halo from
already-clipped tile fetch width.

## Layer 4: command classification tests

These tests prove each command is using the correct model.

They must also prove each command is using the documented ownership rule for that model.

They should be added only after the helper tests above exist.

### A. Representative `CoreOverlap` commands

#### `fragment_kmers`

Required cases:

- Halo-only BED window does not create fetch or output
- Mixed cached span with one core window and one halo-only window narrows fetch from the core
  window only

#### `midpoints`

Required cases:

- Halo-only site is ignored
- Core-overlap site is kept
- Mixed cached span keeps only the core-overlap site

#### `fcoverage`

Required cases:

- Halo-only BED window is ignored for BED outputs
- Core-overlap BED window is kept
- Mixed cached span still behaves like core overlap only

#### `ref_gc_bias`

Required cases:

- BED preparation clips windows to tile-owned core overlap only
- Halo-only BED window does not create a tile-local window

### B. Fragment-reach commands

#### `lengths`

Required cases:

- Right halo-only BED window receives a count
- Mixed core-overlap plus right-halo target case counts only the real target
- Window touching fragment end stays zero under half-open overlap

#### `ends`, aligned and drop modes

Required cases:

- Right halo-only aligned-reachable BED window can count
- Mixed core plus right-halo target case counts the real target

#### `ends`, raw mode

Required cases:

- Left raw-only BED window counts
- Far-right raw-only BED window counts
- Mixed core plus far-right raw case counts only the endpoint-containing target
- Left and right half-open boundary cases stay zero
- Fetch policy tests from Layer 2C still hold at command level

#### `gc_bias`

Required cases:

- BED preparation keeps the same fragment-reach windows that the fixed-size ownership model would
  logically count into
- Right halo-only BED window is prepared
- Mixed core plus right-halo BED windows are both prepared

These tests should fail if `gc_bias` BED mode falls back to strict core overlap while fixed-size
mode still counts by fragment reach.

## Layer 5: BED vs fixed-size consistency

These are cross-mode semantic tests.

They are required only for commands whose contract says BED and fixed-size windowing express the
same ownership model.

### `lengths`

Required cases:

- one BED window `[20,30)` and one fixed-size bin `[20,30)` on the same fragment set produce the
  same counts
- mixed multi-window case preserves row-wise agreement

### `ends`, aligned and drop modes

Required cases:

- equivalent BED and fixed-size windows produce the same motif counts when their windows describe
  the same genomic intervals

### `gc_bias`

Required cases:

- equivalent BED and fixed-size windows prepare and count the same owned fragment contributions

### Raw `ends`

Do not fake a BED vs fixed-size consistency test unless the command defines a true fixed-size
analogue for raw endpoint assignment.

If no such analogue exists, state that explicitly in the implementation notes and do not invent a
weak comparison.

## Boundary and off-by-one checklist

Every implementation in this area must be tested against all of the following boundaries.

- window ends exactly at fragment start
- window starts exactly at fragment end
- window ends exactly at `core_start`
- window starts exactly at `core_end`
- right reachable bound exact equality drops the window
- left reachable bound exact equality drops the window
- right endpoint counting uses `boundary_pos - 1`
- left endpoint counting uses `boundary_pos`
- first tile and last tile near chromosome boundaries keep the correct carried halo

If any new helper is added and these exact-equality cases are not covered, the test suite is incomplete.

## Mixed-case checklist

Every model needs at least one mixed case where:

- one window should be kept by the correct model
- one window should be dropped by a tempting wrong model
- and one additional window ensures the implementation cannot pass by accident with a trivial
  single-window shortcut

Examples:

- left halo-only, core-overlap, right halo-only
- core-overlap, reachable-right, unreachable-right
- core-overlap, endpoint-touching-only, true endpoint-containing

These cases are important because many wrong implementations still pass single-window fixtures.

## Expected failure signals

The test suite should make the common wrong implementations fail in obvious clusters.

### Wrong implementation: everything treated as `CoreOverlap`

Expected failures:

- `lengths` right halo-only tests
- `ends` raw left and far-right tests
- `gc_bias` BED fragment-reach preparation tests

### Wrong implementation: everything trusts the cached candidate span directly

Expected failures:

- `fragment_kmers` halo-only guard tests
- `midpoints` core-overlap guard tests
- `fcoverage` core-overlap guard tests
- core-overlap fetch-envelope mixed tests

### Wrong implementation: fetch narrowing uses a smaller reach than candidate selection

Expected failures:

- fragment-reach fetch-envelope coupling tests
- mixed core plus halo target tests in `lengths`
- corresponding aligned `ends` tests

### Wrong implementation: raw `ends` candidate windows directly define aligned fetch

Expected failures:

- raw `ends` fetch-policy tests
- far-right raw-only cases near the tile edge

### Wrong implementation: half-open boundaries are wrong

Expected failures:

- window touching fragment end
- left raw boundary [7,8) case
- right raw boundary [24,25) case
- exact `core_start` and `core_end` cases

## Suggested file placement

- Shared selection and clamp helper tests:
  - module-local `*_tests.rs`
- Command-specific helper tests:
  - module-local `*_tests.rs`
- Public command behavior regressions:
  - `tests/`

Do not widen helper visibility only to satisfy this plan.

## Acceptance criteria

The implementation work may begin only when the following are true:

- every layer above has an explicit test checklist entry
- every command that uses tile/window helpers is classified into the correct model
- raw `ends` fetch policy is stated explicitly before coding
- BED vs fixed-size consistency expectations are stated per command instead of assumed globally

The refactor is complete only when the implemented tests make the wrong-model substitutions fail
clearly and reproducibly.
