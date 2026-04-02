# Tile, window, and fetch rules spec

Date: 2026-04-01

## Scope

This spec defines the contract between:

- tile construction
- BAM fetch ranges
- per-tile window candidate selection
- optional fetch narrowing against windows

The immediate motivation is the current ambiguity between:

- aligned-coordinate fetch halos
- window selection halos used to find windows that can receive counts from tile-owned fragments

Those are different concepts. They must not share one implicit contract.

This spec is intentionally stricter than the current code. The goal is to remove ambiguity before more helpers are changed.

## Decision summary

### 1. BAM fetch geometry is always aligned-coordinate geometry

All BAM fetch intervals are expressed in aligned reference coordinates.

This includes:

- `Tile.fetch`
- any later narrowed fetch span
- any halo used when clamping a fetch span

Raw clipping, endpoint shifting, and other count-space adjustments do not change the coordinate system of BAM fetch.

### 2. Window selection is a separate problem from BAM fetch

Each command must define which windows are relevant for a tile.

That decision may depend on:

- tile core overlap
- fragment ownership by aligned start
- fragment assignment intervals
- raw clipping reach

But that window-selection rule does not by itself define a BAM fetch interval.

### 3. A helper that narrows fetch from windows may only consume aligned-space relevance

`clamp_fetch_to_window_span(...)` is an aligned fetch helper.

Its input `window_span` must mean:

- "an aligned-coordinate envelope that is safe to widen by aligned halo and then clamp to `tile.fetch`"

It must not be fed a span that only means:

- "some windows are count-relevant in another coordinate model"

### 4. The repo needs two explicit window-selection models

Window selection must be modeled explicitly as one of:

- `CoreOverlap`
- `ReachableFromTileOwnedFragments`

These models can share low-level interval code, but they must not share an ambiguous type or helper contract.

### 5. Adding more relevant windows must never narrow the final fetch interval

If two candidate window sets `A` and `B` satisfy `A subset B`, then the derived fetch interval for `B` must be:

- equal to the fetch interval for `A`
- or wider before the final clamp
- or the same after the final clamp to `tile.fetch`

It must never become narrower.

This is now an explicit invariant of the design, not just a desirable property.

## Definitions

### Tile core

The half-open interval `[tile.core_start(), tile.core_end())`.

This is the region whose ownership is unique per tile.

### Tile fetch band

The half-open interval `[tile.fetch_start(), tile.fetch_end())`.

This is the maximal aligned-coordinate BAM fetch region that the tile is allowed to read.

The tile fetch band must already be large enough for the command's aligned-context needs.

### Aligned fetch halo

The extra aligned bases carried by `tile.fetch` beyond the tile core.

Examples:

- fragment reconstruction halo for paired/unpaired fragment tools
- WPS context halo from fragment length plus window size
- any other aligned reference context needed before per-tile processing starts

### Candidate windows

The windows that the command considers relevant for this tile before counting starts.

This is a window-selection concept, not a BAM fetch concept.

### Fragment ownership rule

The rule that decides which tile owns a fragment for counting.

Examples:

- fragment aligned start lies in the tile core
- fragment aligned end lies in the tile core
- fragment midpoint lies in the tile core

This is a command-level semantic choice.

It must never be left implicit.

Any helper that depends on fragment ownership must state in its docstring which ownership rule it
assumes.

### Fetch narrowing

An optional step that takes the already valid `tile.fetch` band and intersects it with a smaller aligned-space envelope derived from the current tile's relevant windows.

Fetch narrowing is an optimization only. It must never change which counts are possible.

## Fundamental separation

The system has three separate layers:

1. `build_tiles(...)` decides the maximal aligned fetch band for each tile.
2. Per-command logic decides which windows are relevant for the tile.
3. Some commands may further narrow the aligned BAM fetch band from those relevant windows, but only when there is a proven aligned-space derivation.

The current bug comes from collapsing steps 2 and 3 into one ambiguous helper path.

## Window-selection models

### Model A: `CoreOverlap`

A window is relevant for a tile when the window interval overlaps the tile core.

Use this when:

- the command writes only core-overlapping windows or sites
- tile ownership is fundamentally about the tile core, not about fragments that start there and travel onward
- later outputs clip to the tile core or otherwise treat core overlap as the ownership boundary

Formal rule:

- a BED window `[w_start, w_end)` is relevant iff `w_end > core_start && w_start < core_end`

Properties:

- windows fully left of the core are irrelevant
- windows fully right of the core are irrelevant
- widening this candidate set with non-overlapping halo-only windows is semantically wrong for this model

### Model B: `ReachableFromTileOwnedFragments`

A window is relevant for a tile when at least one fragment owned by the tile can contribute to that window.

Use this when:

- fragment ownership is "aligned fragment start lies in the tile core"
- the tile may count into windows outside the core
- the command's output semantics are about fragment-to-window assignment, not only about core-overlap windows

Formal rule:

- define the command's fragment ownership rule
- define the command's counting or assignment interval
- a BED window is relevant iff there exists at least one owned fragment whose counting interval can contribute to that window under the command's assignment rule

Properties:

- windows fully right of the core may be relevant
- windows fully left of the core may be relevant if the counting interval can extend left of aligned start
- this model is allowed to use left and right reach bounds that differ

## Fetch narrowing rules

### Rule 1: fetch narrowing is optional

A command may always use the full `tile.fetch` band.

This is always semantically safe when `tile.fetch` already carries the command's required aligned halo.

### Rule 2: fetch narrowing may only consume aligned-space support envelopes

A narrowed fetch band may be derived from windows only when the command can prove that the windows define a safe aligned-coordinate envelope for all aligned reads that could produce the tile's counts.

This proof must hold for the command's actual fragment ownership and assignment semantics.

### Rule 3: `clamp_fetch_to_window_span(...)` stays an aligned clamp helper

`clamp_fetch_to_window_span(...)` should keep its current role:

- widen an already aligned-space support interval by aligned halo
- clamp that widened interval back to `tile.fetch`

It should not become a generic "turn any candidate windows into a fetch band" helper.

### Rule 4: count-space reach must not be passed to aligned fetch helpers without an explicit support-envelope conversion

Raw endpoint windows in `ends` are the clearest example.

Those windows are count-relevant, but they do not directly define where the aligned read must be fetched from.

So a raw endpoint window span may narrow aligned BAM fetch only after the code converts it into a
proven aligned support envelope and then widens by a halo that is at least as permissive as the
reach used to keep those windows.

### Rule 5: larger relevant window sets are monotone

When fetch narrowing is used, enlarging the relevant window set must not narrow the final fetch interval.

That monotonicity must hold before and after clamping to `tile.fetch`.

### Rule 6: if a reach value makes a window relevant, the same reach value must be respected when narrowing fetch

Some commands include extra windows because tile-owned fragments can reach beyond the tile core.

If fetch narrowing is derived from those windows, it must preserve the same reach.

Plainly:

- if a command says "a fragment can reach a window `R` bases away, so keep that window"
- then fetch narrowing must also assume "a read starting up to `R` bases earlier may still be needed"

It must not do this:

- include the window using reach `R`
- then narrow BAM fetch using a smaller value than `R`

Worked left-bound example:

- a command owns fragments by aligned start in the tile core
- a BED window starting at `min_ws` is kept because an owned fragment may extend `R` bases to the right and still reach it
- then the narrowed aligned BAM fetch must start no later than `min_ws - R` before the final clamp to `tile.fetch`

Why:

- a read that starts at `min_ws - R` is exactly the read that just barely makes that window relevant
- if fetch starts to the right of `min_ws - R`, that read is no longer fetched
- then the code keeps a window that its own narrowed fetch can no longer satisfy

This rule does not require one symmetric helper or one universal halo value. It only requires that
the fetch-narrowing logic uses a bound that is at least as permissive as the bound that was used to
decide the window was relevant.

## Per-command classification

### Commands that use `CoreOverlap`

- `midpoints`
- `fcoverage`
- `wps`
- `wps-peaks` through `wps`
- `fragment_kmers`
- `ref_gc_bias`

These commands fundamentally treat windows or sites as tile-core-owned outputs.

For them:

- BED candidate windows should be selected by core overlap
- fetch narrowing may use the core-overlap min/max window span in aligned coordinates

### Commands that use `ReachableFromTileOwnedFragments`

- `lengths`
- `ends`
- `gc_bias`

These commands count fragments owned by aligned start in the tile core and may legitimately contribute to windows outside the core.

For them:

- BED candidate windows must be selected by fragment reach, not by core overlap
- the precomputed span and the later counting logic must use the same reach model

## Exact reach rules for fragment-owned models

### `lengths`

Owned fragment rule:

- a fragment is owned by the tile iff its aligned start lies in the tile core

Count interval:

- aligned fragment interval

Safe BED reach envelope:

- left reach: `0`
- right reach: `max_fragment_length`

So a BED window is a candidate when it can overlap a fragment interval starting inside the tile core with aligned length up to `max_fragment_length`.

Fetch narrowing:

- allowed
- derive the min/max from the reachable BED candidate set
- clamp with aligned halo `max_fragment_length`

This is safe because the same aligned interval drives both counting and fetch.

### `ends` with `clip_strategy = aligned` or `drop`

Owned fragment rule:

- aligned fragment start lies in the tile core

Count interval:

- aligned assignment interval

Safe BED reach envelope:

- left reach: `0`
- right reach: `max_fragment_length`

Fetch narrowing:

- allowed
- same rule as `lengths`

This is safe because count-space and aligned-space interval geometry are the same.

### `ends` with `clip_strategy = raw`

Owned fragment rule:

- aligned fragment start lies in the tile core

Count interval:

- raw-shifted assignment interval

Safe BED candidate envelope:

- left reach: `max_soft_clips`
- right reach: `max_fragment_length + max_soft_clips`

This rule defines which raw BED windows must stay visible to the tile.

Fetch narrowing:

- allowed
- derive the aligned support interval from the full raw candidate-window extent
- widen that interval by an aligned halo large enough to preserve every aligned read that could
  make those raw windows relevant
- a conservative halo is allowed if it is easier to prove safe, for example using
  `max_fragment_length + max_soft_clips` symmetrically even though the raw left reach is smaller
- clamp the widened aligned interval back to `tile.fetch`
- `KeepTileFetch` remains a safe fallback and may still be useful in tests, but it is not the
  target production contract for raw `ends`

This is the key asymmetry:

- raw endpoint windows may be farther left or right than the aligned read span
- the windows are still count-relevant
- so fetch narrowing must use an aligned support-envelope conversion, not the raw window extrema
  as-is

### Fixed-size mode is an explicit regression area

This must be treated as an explicit regression area in general, not only for BED mode and not only
for raw `ends`.

Reason:

- this repo has already had real `--by-size` bugs
- so the implementation must not rely on "BED mode was fixed, therefore size mode is probably fine"

Required semantic rule:

- every command whose semantics depend on window ownership or fragment reach must test `--by-size`
  explicitly and independently of BED mode

Clarification:

- fixed-size counting does not use the same precomputed BED candidate-span path as BED mode
- fixed-size windows are implicit bins, and the overlap query is built directly against those bins
- so BED-mode correctness does not prove fixed-size correctness, even when the command's high-level
  ownership model is the same

But:

- different implementation path does not reduce the testing requirement
- tile decomposition must not change which fixed-size rows receive counts
- half-open boundary rules must hold the same way in fixed-size mode as in BED mode whenever the
  command contract says the two window modes represent the same partition

Important subcase:

- raw `ends` still needs its own fixed-size tests because raw endpoint reach is asymmetric and has
  already been a source of bugs

### `gc_bias`

Owned fragment rule:

- aligned fragment start lies in the tile core

Count interval:

- `fragment_assignment_interval(...)` for the active assignment mode

Required semantic rule:

- BED and fixed-size modes must use the same fragment-reach ownership model

This command already proves the fragment-reach idea in fixed-size streaming:

- a "next" window may lie wholly outside the current tile core and still receive counts from tile-owned fragments

So the BED path must not keep a stricter core-overlap-only interpretation if the fixed-size path uses fragment reach.

Fetch narrowing:

- currently unnecessary
- full `tile.fetch` is semantically safe and already used

## Required helper decomposition

The implementation should move toward four explicit helper roles.

## Required helper docstrings

Every new or newly-repurposed helper in this area must state its assumptions in the docstring.

This is not optional documentation polish. It is part of the contract.

At minimum, the docstring must say all of the following:

- what coordinate space the helper consumes and returns
  - aligned BAM fetch space
  - BED window candidate space
  - endpoint count space
- which window-selection model it implements
  - `CoreOverlap`
  - `ReachableFromTileOwnedFragments`
- which fragment ownership rule it assumes
  - for example: "fragment is owned iff aligned start lies in the tile core"
- which counting or assignment interval it assumes
  - aligned fragment interval
  - aligned assignment interval
  - raw-shifted assignment interval
- whether it is allowed to narrow aligned BAM fetch
- if it narrows fetch, which proven aligned-space reach bound makes that safe

This must be written plainly enough that a reader can answer, without reading the implementation:

- does this helper assume "starts in tile"?
- would it still be correct for an "ends in tile" command?
- does it work in aligned geometry, raw endpoint geometry, or both?
- is it selecting candidate windows, deriving an aligned support envelope, or clamping an already
  aligned envelope?

If the answer to "would it still be correct for an `ends in tile` command?" is "no", the docstring
should say so directly rather than implying generality it does not have.

Bad docstring pattern:

- "Returns relevant windows for this tile"

Why it is bad:

- it hides the ownership rule
- it hides the coordinate model
- it hides whether "relevant" means core overlap or fragment reach

Acceptable docstring pattern:

- "Returns BED candidate windows for fragments owned by aligned start in the tile core. This uses
  `ReachableFromTileOwnedFragments` in aligned fragment geometry and may include windows right of
  the core. It is not valid for commands owned by fragment end or midpoint."

### 1. Tile construction helper

Responsibility:

- build `tile.core`
- build `tile.fetch`

Inputs:

- aligned fetch halo only

It must not know anything about BED candidate windows.

### 2. Candidate window span helper

Responsibility:

- precompute per-tile BED candidate spans

Inputs:

- an explicit fragment ownership rule when the model is fragment-owned
- an explicit selection model
- model-specific left and right reach parameters

It must not imply anything about BAM fetch narrowing on its own.

### 3. Optional aligned fetch-envelope helper

Responsibility:

- convert a tile's relevant windows into an aligned-space support envelope when that conversion is valid for the command

This helper must be model-aware.

It is valid for:

- `CoreOverlap`
- fragment-reach models whose count interval is also aligned interval geometry

It is not valid for:

- `ends` raw without an explicit aligned-space derivation

### 4. Final clamp helper

Responsibility:

- widen an aligned-space support envelope by aligned halo
- clamp to `tile.fetch`
- return `None` for empty results

This is the current role of `clamp_fetch_to_window_span(...)`.

## Type and API rules

### Rule 1: `TileWindowSpan` must stop being semantically ambiguous

Today `TileWindowSpan` is only an index slice, but different call sites implicitly treat it as:

- a core-overlap window set
- a fragment-reach candidate set

That ambiguity is the root design problem.

The implementation must do one of:

- split it into separate types
- or add an explicit model tag wherever a span is passed across helper boundaries

### Rule 2: fetch helpers must know which model produced the span

A fetch helper must not accept only:

- `tile`
- `tile_window_span`
- `windows_chr`

without also knowing which model produced that span.

Otherwise it will eventually reinterpret a fragment-reach span as a core-overlap span again.

If the model is fragment-owned, the helper contract must also make the fragment ownership rule
explicit. "Starts in tile" is one valid rule, but it is not the only possible future rule.

If a future command uses a different ownership rule, such as "ends in tile", it must either:

- use a helper whose contract explicitly supports that ownership rule
- or introduce a new helper instead of silently reusing one built for "starts in tile"

### Rule 3: docs must stop calling index spans "min/max bounds"

An index slice is not a coordinate envelope.

The docs and names should distinguish:

- candidate index span
- aligned support envelope

### Rule 4: helper docstrings must state the semantic assumptions they rely on

For this area, docstrings must say the assumptions plainly.

At minimum, relevant helpers must state:

- whether they operate in aligned fetch space or candidate-window space
- whether they assume `CoreOverlap` or fragment-reach semantics
- if fragment-reach is used, what fragment ownership rule is assumed
- if fragment-reach is used, what counting or assignment interval is assumed

This is not optional documentation polish.

The same helper name must not quietly mean:

- "fragments starting in tile" for one command
- "fragments ending in tile" for another command

without that difference being made explicit in the type, arguments, or docstring

## Acceptance criteria

The eventual implementation should satisfy all of these.

- A core-overlap command never counts into a BED window that does not overlap the tile core.
- A fragment-reach command can count into a BED window outside the tile core when an owned fragment can reach it.
- Adding extra relevant windows never narrows the final fetch interval.
- `ends` raw can keep halo-only raw endpoint windows while still deriving a narrower aligned BAM
  fetch that preserves every eligible aligned read.
- BED and fixed-size modes of the same command agree on the underlying ownership model.
- Helper docs describe whether they operate in aligned fetch space or candidate-window space.

## Current inconsistencies this spec resolves

- `TileWindowSpan` docs still read like a core-overlap span even though some tests now keep halo-only windows.
- `window_fetch` treats a cached index span like a fetch-ready min/max envelope.
- `lengths` and `ends` precompute fragment-reach candidate spans but later reinterpret them with core-overlap helpers.
- `gc_bias` fixed-size logic already behaves like a fragment-reach model, while the BED preparation path still uses core-overlap helpers.

## Non-goals

This spec does not lock in final function names.

This spec does not require every command to narrow fetch.

This spec does not change the existing fragment ownership rule that uses aligned fragment start inside the tile core to prevent double counting.
