## `--normalize-by-length` restore-mean plan

Date: 2026-04-22

## Goal

Keep the current clean `--normalize-by-length` semantics by default, while
allowing an optional mode that restores the mean scale afterward using the
observed mean normalization length.

This solves a practical usability gap:

- users cannot rescale correctly downstream unless the command exposes the
  exact mean normalization length it used

## Current semantics

Today `--normalize-by-length` means:

- each fragment contributes total mass `1.0`
- per-base weight is `1 / counted_length`

Interpretation:

- `--per-window total` is count-like and approximates fragment counts in
  sufficiently large windows
- `--per-window average` is count-density-like and approximates fragment counts
  divided by window length

This is good and should remain the default.

## Desired design

Replace the boolean-style concept with a small mode family:

- default: no length normalization
- `--normalize-by-length`
  shorthand for the current unit-mass behavior
- `--normalize-by-length=restore-mean`
  normalize each fragment to total mass `1.0`, then multiply by the observed
  mean normalization length

This keeps one conceptual argument:

- how should length normalization behave?

instead of adding a second flag for rescaling.

## Why `restore-mean`

`restore-mean` is preferred over broader names like `restore-scale` because:

- it describes the actual operation more honestly
- it does not imply that GC correction or other scaling is being undone
- it points at the mean normalization denominator specifically, not generic
  scale recovery

## What "mean normalization length" means

This must use the same normalization-length semantics as the normalization
itself.

That means:

- use the full fragment span when no explicit counted segments are present
- use the summed counted reference-segment lengths when segment-aware counting
  is active
- calculate after the same fragment filters are applied
- do not ask users to guess this value from an external QC tool
- define it before later masking or window-exclusion steps

This statistic should be named explicitly in code, logs, and docs.

Suggested name:

- `mean_normalization_length`

Recommended docs wording:

- the normalization length is the denominator used by
  `--normalize-by-length`
- for ordinary fragments, it is the fragment span length
- for segment-aware fragments, it is the summed counted reference-segment
  length
- it is defined before later masking/window exclusion steps

## User-facing behavior

### `--normalize-by-length`

Meaning:

- current behavior
- one fragment -> total mass `1.0`

### `--normalize-by-length=restore-mean`

Meaning:

- first perform ordinary length normalization
- then multiply by the observed mean normalization length

Interpretation:

- restores the plain global mean level by multiplying with the observed mean
  normalization length
- removes the strongest length bias while keeping the output easier to compare
  with non-length-normalized coverage-like outputs

## Logging and metadata

The observed mean normalization length should always be reported when
length normalization is used, regardless of mode.

Minimum requirement:

- log `mean_normalization_length`

Optional but reasonable:

- save a tiny scalar sidecar file with the observed mean normalization length
  when persistent reuse is useful

If a sidecar is added, keep it trivial and machine-readable. Suggested content:

- header: `mean_normalization_length`
- one value row

Reason:

- users otherwise cannot reproduce the exact restore-mean rescaling
- even for plain unit-mass mode, it is valuable context for interpretation

## Code touchpoints

### 1. CLI/config shape

Change `FCoverageConfig.normalize_by_length: bool` into an enum-like mode.

Needs updates in:

- [config.rs](/Users/au547627/Documents/Development/rust/cfDNAlab/src/commands/fcoverage/config.rs)
- any builders/setters using `set_normalize_by_length(...)`
- tests that assume boolean behavior

Suggested conceptual modes:

- `off`
- `unit-mass`
- `restore-mean`

Bare `--normalize-by-length` should map to `unit-mass`.

Clap shape:

- `None` -> `off`
- bare `--normalize-by-length` -> `unit-mass`
- `--normalize-by-length=restore-mean` -> `restore-mean`

### 2. Base-weight calculation

Current base weight:

- `1.0`
- or `1.0 / counted_length`

New logic:

- `off` -> `1.0`
- `unit-mass` -> `1.0 / counted_length`
- `restore-mean` -> count exactly like `unit-mass` during fragment iteration,
  then apply the observed mean normalization length in a later finalization
  step

Reason:

- there is no first pass
- the mean normalization length is only known after the normal counting pass
- trying to inject `restore-mean` directly into per-fragment base weights would
  require a statistic that does not exist yet

### 3. Mean counted-length accumulation

Need one run-level statistic:

- sum of counted lengths across included fragments
- divided by counted fragments

Important:

- this should reflect the same fragment population that contributes to the
  output
- excluded fragments must not contribute
- the sum must only increase when a counted fragment is owned by the current
  tile for normalization statistics
- ownership should be assigned once per fragment, by fragment start landing in
  exactly one tile core under the current contiguous non-overlapping core
  layout
- this intentionally differs from `counted_fragments`, which can be larger
  because tile halos may make one fragment visible in multiple tiles for
  coverage statistics
- fragments that passed early filters but never contributed to the tile core
  must not affect the mean

Suggested names:

- `tile_owned_normalization_length_sum`
- `tile_owned_normalization_fragments`
- `mean_normalization_length`

Implementation note:

- `FCoverageCounters` already carries `counted_fragments`
- the normalization-length accumulator should keep its own tile-owned fragment
  count so the mean is defined from the same once-per-fragment sample it sums

### 3b. Zero-count edge case

If `counted_fragments == 0`, the mean is undefined.

Required behavior:

- do not attempt `restore-mean` multiplication
- keep the output as the already-counted zero-valued signal
- log a warning explaining that `mean_normalization_length` was undefined
  because no fragments contributed

## Finalization design

`restore-mean` should be implemented as a late scaling step on count-like raw
aggregate values, not as a rewrite of already-derived human-facing summary
statistics.

This matters most for `--by-size`, where the current aligned fast path writes
already-finalized tile outputs and later concatenates them. That path is fine
for ordinary `unit-mass`, but becomes brittle for `restore-mean`, especially
for summary statistics.

### Ordinary aggregate outputs

For `average` and `total`, `restore-mean` can multiply the raw count-like
values by `mean_normalization_length` before final output rounding.

### Summary-stats outputs

For `summary-stats`, the clean design is:

- keep exact raw additive fields until the final write stage
- apply `restore-mean` there
- derive `average_coverage`, `total_coverage`, `variance_coverage`,
  `sd_coverage`, and `coefficient_of_variation_coverage` afterward

This avoids schema-aware rewriting of already-finished summary-stat files, which
would be easy to get wrong.

### Aligned `--by-size` fast path

When tile and size-window boundaries align, there is currently a fast path that
writes per-tile final outputs and concatenates them without reopening the rows.

For `--normalize-by-length=restore-mean`, that fast path should be disabled in
its current form.

Instead:

- keep the existing alignment marker so the code still knows no cross-tile
  reduction is needed
- write aligned per-bin raw rows rather than already-finalized rows
- run a lightweight finalization pass over those raw rows
- apply `mean_normalization_length` there
- then derive final summary stats and round once

This is intentionally slower than pure concatenation. That is acceptable for an
explicit opt-in mode where correctness and trustworthiness matter more than the
last fast-path optimization.

### 4. Run statistics/logging

Add one line to the run summary when length normalization is enabled:

- `Mean normalization length: <value>`

Keep it concise.

### 5. Docs/help strings

Update the `fcoverage` help text so the modes are clear without turning into a
long tutorial.

Main points to preserve:

- unit-mass mode keeps one fragment at total mass `1.0`
- `total` is count-like
- `average` is count-density-like
- `restore-mean` multiplies the final signal by the observed mean counted
  normalization length
- this restores the plain global mean level in the ordinary length-normalized
  case, but not local length-driven variation

### 6. Output scope

`restore-mean` should be a true `fcoverage` mode, not an aggregate-only
special case.

That means the design must cover all output families:

- positional outputs
- BED-window aggregate outputs
- grouped BED aggregate outputs
- size-window aggregate outputs

The implementation details can differ by output family, but the user-facing
meaning must stay the same:

- first count in ordinary `unit-mass` space
- then apply the observed mean normalization length before final output
  rounding/derivation

This keeps the mode coherent across the command instead of making it an
aggregate-only quirk.

### 7. `--by-size` implementation touchpoints

Rework the current aligned-size behavior only for `restore-mean`.

Main points:

- do not rely on concatenating already-finalized aligned tile outputs
- preserve the existing "no reduction needed" alignment marker
- add a raw-row finalization path that can apply `restore-mean` before final
  row derivation
- keep the ordinary aligned concatenation fast path for modes that do not need
  the extra scaling step

### 8. Other output-family touchpoints

The plan should stay explicit that `restore-mean` is not just a `--by-size`
concern.

For the other output families:

- positional outputs need an explicit scaled-merge path:
  write the usual unit-mass tile files during counting, then during final merge
  reopen the tile rows, multiply the value column by
  `mean_normalization_length`, and write the scaled merged output instead of
  using plain file concatenation
- BED and grouped BED aggregates should apply the scalar before final row
  derivation/rounding in the existing reducer/finalizer flow

The exact plumbing can differ, but the mode semantics should not.

### 9. Downstream command docs

Re-read and adjust any docs that rely on current boolean wording, especially:

- `fragment-count-weights`
- any internal comments/specs describing `normalize_by_length`

`fragment-count-weights` should remain explicitly `unit-mass`.

Reason:

- it normalizes the resulting stride-bin values into scaling factors
- a global `restore-mean` multiplier would cancel during that normalization
- exposing `restore-mean` there would add confusion without changing the final
  factors

## Testing plan

Add targeted tests, not a large matrix.

### Semantics tests

- `off` gives base weight `1.0`
- `unit-mass` gives `1.0 / counted_length`
- `restore-mean` uses the same intrinsic base weight as `unit-mass` during
  counting and differs only at the finalization stage

### Normalization-length consistency tests

- segment-aware fragments use summed counted segment length
- plain fragments use full fragment span length

### Logging/stat reporting tests

- when normalization mode is active, mean counted length is reported
- when no fragments are counted, the warning path is used

### Regression-level interpretation tests

Use a tiny synthetic example where the expected mean normalization length is
obvious, then verify:

- `unit-mass` and `restore-mean` differ by exactly that mean for ordinary
  aggregate outputs

### Positional-output tests

- positional `restore-mean` uses the scaled-merge path rather than plain
  concatenation
- the merged positional values equal the unit-mass values multiplied by
  `mean_normalization_length`

### `fragment-count-weights` tests

- internal `fcoverage` still uses `unit-mass`
- no restore-mean-specific behavior is exposed or needed there

### Aligned `--by-size` tests

- when alignment holds and mode is `unit-mass`, the ordinary fast concat path is
  still used
- when alignment holds and mode is `restore-mean`, the code takes the slower
  raw-row finalization path instead
- aligned `summary-stats` with `restore-mean` matches the same result that would
  be obtained by rescaling raw count-like fields before deriving the summary
  columns

## Non-goals

This plan does not propose:

- changing the default semantics
- silently rescaling by mean counted length
- redefining `fragment-count-weights`
- exposing `restore-mean` through `fragment-count-weights`
- broad changes to GC correction or scaling-factor logic

## Recommendation

Implement this if and only if you want `fcoverage` to support both:

- mathematically clean unit-mass fragment support
- and a coverage-like restored-mean variant that users can compare more
  directly to ordinary coverage outputs

The design is coherent as long as:

- `unit-mass` stays the default for `--normalize-by-length`
- `restore-mean` is explicit
- `mean_normalization_length` is logged and reproducible
- `restore-mean` is applied before final row derivation, not by trying to patch
  already-finished summary-stat files

## Detailed restore-mean test matrix

The implementation should be considered incomplete until all of the following
behaviors are covered explicitly.

### Config and bookkeeping

- bare `--normalize-by-length` selects `unit-mass`
- `--normalize-by-length=restore-mean` selects `restore-mean`
- run statistics report `mean_normalization_length`
- the optional scalar sidecar, if implemented, writes the expected key/value
- `mean_normalization_length` uses the same denominator semantics as
  `--normalize-by-length`
- the accumulator advances only for tile-owned counted fragments, so each
  fragment contributes once even when tile halos make it visible more than once
- zero counted fragments leave output unchanged and emit the documented warning

### Denominator semantics

- ordinary paired fragments use fragment span length
- unpaired `reads_are_fragments` uses read span length
- gapped and ref-skip fragments use summed counted reference-segment length
- `--ignore-gap` changes the normalization denominator to the remaining counted
  spans
- GC correction and genomic scaling still multiply on top of the restored
  signal

### Positional outputs

- global positional output restores the expected per-base values for mixed
  normalization lengths
- global positional output is basewise tile-size invariant even when the merged
  bedGraph run boundaries differ
- `keep_zero_runs` still behaves correctly under `restore-mean`
- blacklist masking still removes masked positions after restore-mean scaling
- scaling factors and GC correction still act on genomic coordinates before the
  late restore-mean multiply
- whole-genome positional output across three chromosomes keeps chromosome order
  and expected restored values

### Windowed positional outputs

- `OnlyIncludeThesePositionsUnique` restores the selected positional values
- `OnlyIncludeThesePositionsIndexed` restores the selected positional values and
  preserves original window indices
- overlapping indexed BED windows still duplicate overlap intentionally under
  `restore-mean`
- tile-size changes do not alter the basewise restored signal for positional BED
  outputs

### By-size aggregate outputs

- `total` restores each fragment to `mean_normalization_length` total mass
- `average` restores the plain global mean level in full-chromosome windows
- aligned fast path and general reducer path agree exactly under
  `restore-mean`
- bin coordinates remain the full requested bin coordinates when bins cross tile
  boundaries
- three-chromosome by-size output keeps chromosome order and restored values

### By-size summary stats

- raw additive fields are restored before final derivation:
  `coverage_sum`, `coverage_sum_of_squares`
- derived fields are recomputed from restored raw fields:
  `average_coverage`, `total_coverage`, `variance_coverage`, `sd_coverage`,
  `coefficient_of_variation_coverage`, `covered_fraction`
- aligned and non-aligned tile paths agree exactly
- zero-support windows still keep the documented summary-stat conventions

### BED aggregate outputs

- `total` restores expected window sums
- `average` restores expected window means
- BED windows crossing one tile, two tiles, and more than two tiles still give
  identical final rows
- same-start windows remain coordinate-sorted after reduction under
  `restore-mean`
- three-chromosome BED-average output still respects global window indices
- chromosomes without BED windows are still skipped cleanly

### BED summary stats

- raw and derived summary-stat fields are restored correctly
- overlapping BED windows crossing tiles stay invariant under tile-size changes
- halo-only windows are not double-counted after restore-mean scaling

### Grouped BED outputs

- plain grouped `total` keeps site-weighted semantics under `restore-mean`
- `total-on-unique-bases` merges same-group overlaps once under
  `restore-mean`
- `average-on-unique-bases` restores expected group means
- grouped sidecar ordering remains stable
- grouped outputs stay invariant when grouped segments cross tile boundaries
- filtered-out chromosomes still remove their grouped rows cleanly

### Grouped summary stats

- plain grouped summary stats restore both raw and derived fields correctly
- unique-base grouped summary stats restore both raw and derived fields
  correctly
- grouped summary-stat tile-size invariance still holds
- three-chromosome grouped summary-stat workflows still support downstream
  Pearson-style derivations with restored values

### Related command expectations

- `fragment-count-weights` keeps using `unit-mass`
- no restore-mean-specific behavior leaks into normalized scaling-factor output
