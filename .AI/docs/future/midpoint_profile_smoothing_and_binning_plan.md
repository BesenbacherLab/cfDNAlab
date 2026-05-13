# Midpoint Profile Smoothing And Binning Plan

Date: 2026-05-13

Scope: Future design for optional smoothing, flanking, and binning in
`cfdna midpoints`.

This is a design note, not a current spec. The goal is to create generally
useful midpoint profiles across interval sizes while allowing Griffin-like
profile processing when requested.

## Goals

- Keep full-resolution midpoint profiles as the clean default.
- Make smoothing explicit and reproducible.
- Support Savitzky-Golay smoothing and binned output as opt-in transforms.
- Use flanks for computation only. Input intervals define the output window.
- Prefilter intervals near blacklisted regions by default to avoid edge-biased
  profiles.
- Keep counting behavior unchanged except for derived flanks and interval
  prefiltering. Smoothing and binning are post-processing transforms.
- Avoid growing the CLI with controls that can be derived or handled
  downstream.
- Write one selected profile output. Users who need a second transform can rerun
  the command.

## Non-goals

- Do not add a separate output-window argument. The BED intervals already define
  the output windows.
- Do not add separate bin aggregation modes unless there is a concrete need.
  Average is the most stable default because it handles a shorter final bin
  cleanly. Sum can be recovered by multiplying by the number of covered bases in
  the bin.
- Do not reproduce Griffin's large normalization window as a default behavior.
  Griffin uses a broader profile window for baseline normalization, but this
  project already supports scaling factors for CNA-style normalization.
- Do not expose a flank-size argument unless a concrete use case appears.
  Smoothing flanks can be derived from the smoothing window.

## Griffin Reference Behavior

Griffin's nucleosome profile workflow uses a larger working window than it saves.
The default working window is -5000 to 5000 bp, while the default saved profile
window is -1000 to 1000 bp.

The broad working window is used to:

- Fetch and bin coverage around each site.
- Apply excluded-region, mappability, and outlier masks.
- Average sites into a composite profile.
- Apply Savitzky-Golay smoothing.
- Normalize the final averaged profile to mean 1 across the broad working
  window.

Only after those steps does Griffin truncate to the saved output window.

For `cfdna midpoints`, the lasting lesson is not that every profile needs a
5 kb margin. The useful part is that smoothing should not be applied right at
the output boundary when avoidable. A smaller computation flank is enough if
flanking is only for smoothing.

Griffin compatibility is a preset or documentation target, not the governing
design. `cfdna midpoints` should produce good midpoint profiles for TFBS, TSS,
and other site classes without assuming Griffin's window sizes, binning, or
normalization choices.

## Proposed CLI Surface

Tentative options:

```text
--bin-size 1
--smooth raw
--smooth savgol
--smooth savgol=165
--keep-blacklisted-intervals
```

Proposed defaults:

```text
--bin-size 1
--smooth raw
```

Notes:

- `--bin-size 1` is the full-resolution path.
- `--smooth raw` preserves current semantics.
- `--smooth savgol` uses the documented default Savitzky-Golay setting.
- `--smooth savgol=165` means a 165 bp Savitzky-Golay window. Savitzky-Golay
  smoothing uses polynomial order 3.
- `--keep-blacklisted-intervals` disables default interval-level blacklist
  prefiltering for users who intentionally prepared or want to keep the original
  site set.

The exact spelling can change during implementation. The important model is one
smoothing option with method-specific settings, not several independent CLI
flags for every smoothing parameter.

Smoothing flanks are derived, not user-facing:

```text
raw profiles:      smoothing_flank = 0
smoothed profiles: smoothing_flank = floor(smoothing_window_bp / 2)
```

The derived flank is used for computation and blacklist interval prefiltering.
It is trimmed before writing the output profile.

## Processing Model

The conceptually clean pipeline is:

1. Derive the computation flank from the selected smoothing mode.
2. Prefilter intervals that are too close to blacklisted regions.
3. Expand each remaining interval by the derived flank for computation.
4. Count midpoint profiles at base resolution across the expanded interval.
5. Apply optional smoothing across the expanded full-resolution profile.
6. Trim away flank positions.
7. Bin the trimmed profile if `bin_size > 1`.
8. Write the selected output profile.

This order keeps smoothing independent of the output compression choice. It also
means `savgol=165` always means 165 bp, not 11 bins under one setting and
165 bases under another.

The smoothing setting applies to the whole run. Different site classes that need
different smoothing scales should be processed in separate runs or smoothed
downstream. Per-group smoothing would make one output tensor mix different
signal-processing assumptions.

## Smoothing Scale

The smoothing window is a signal-scale choice, not an interval-size choice. The
same Savitzky-Golay window can be applied to short and long intervals as long as
the computation flank provides enough context. A 165 bp window has the same
local meaning in a 2 kb TFBS profile and a 10 kb TSS profile.

The default smoothing behavior should therefore not scale with interval length.
It should start as a documented convenience setting for local midpoint profiles,
with `savgol=165` as the initial default for `--smooth savgol` until tuning
shows a better general default. Griffin-compatible behavior can be documented in
a guide through the primitive options rather than encoded as a separate preset.

The output interval length is user intent and should not be changed for
smoothing. If the requested smoothing window is too large for the output
interval, fail early with a suggested smaller explicit setting. This applies to
bare `--smooth savgol` as well as `--smooth savgol=<window>`. Intervals shorter
than 7 bp should fail when smoothing is requested.

Validation should fail clearly when a requested smoothing setting cannot be
applied, for example:

- The smoothing window is even.
- The output interval is shorter than 7 bp when smoothing is requested.
- The smoothing window is longer than the output interval.
- The expanded interval is too short for the smoothing window.

## Memory Trade-off

Full-resolution counting followed by final binning is easiest to reason about,
but it can cost much more memory when users request coarse output bins.

The counting stage itself does not need to change, apart from using expanded
intervals for derived smoothing flanks and dropping blacklisted intervals before
counting. The new behavior is primarily:

- Interval prefiltering before counting.
- Profile smoothing after final aggregation.
- Optional binning after smoothing and flank trimming.

Smoothing should happen on the final profile tensor, not per input window. The
cost therefore scales with:

```text
groups * length_bins * output_positions * smoothing_window_bp
```

It does not scale with the number of input intervals once they have been
aggregated into groups. For example, tens of millions of windows collapsed into
1500 groups are expensive during counting and merging, but smoothing sees only
the final `group x length_bin x position` tensor. This makes runtime-generated
Savitzky-Golay coefficients acceptable unless profiling shows the final tensor
transform is a bottleneck.

Counting-time binning is a valid future optimization:

- It reduces tile partial size and merge memory.
- It avoids building full-resolution arrays that will be compressed anyway.
- It changes the smoothing semantics if smoothing is done after binned counting.

If counting-time binning is implemented, it should be treated as a different
execution strategy with explicit semantics:

- For `--smooth raw`, counting directly into bins is fine.
- For `--smooth savgol=...`, keep the semantic path as full-resolution
  smoothing before final binning.

## Derived Smoothing Flank

If flanking is only for smoothing, the mathematically necessary flank is the
filter support radius. For a single-pass Savitzky-Golay filter with an odd
window, that is `(smoothing_window_bp - 1) / 2` on each side. For
`savgol=165`, that is 82 bp.

The default should be derived rather than user-specified:

```text
raw profiles:      flank_size = 0
smoothed profiles: flank_size = floor(smoothing_window_bp / 2)
```

This keeps the CLI smaller and makes `--smooth savgol=165` fully define the
smoothing behavior. If users later need extra context for a different reason,
that should be designed as a separate feature instead of preemptively adding a
generic flank argument.

Do not double the Savitzky-Golay flank by default. With a one-pass finite
impulse response filter, output positions at least one support radius away from
the expanded interval edge are unaffected by edge handling. Doubling the flank
would add counting work and drop more blacklist-adjacent intervals without
changing the retained smoothed values.

For future smoothing methods, derive the minimum required support radius from
the method. It is always acceptable to use a larger internal flank if an
implementation needs it, as long as the flank is computation-only and trimmed
before writing.

## Blacklist Interval Prefiltering

The current fragment-level blacklist filter removes individual fragments that
overlap blacklisted regions. That is not enough to avoid edge bias in aggregate
profiles. Sites near blacklisted regions can lose fragments asymmetrically, and
with smoothing the bias can also bleed from computation-only flanks into the
retained output window.

Future midpoint smoothing and binning should therefore add interval-level
blacklist prefiltering. When blacklists are supplied, intervals should be
dropped if their output interval plus the relevant safety margin overlaps a
blacklisted region.

The default safety margin should be:

```text
interval_blacklist_margin = ceil(max_fragment_length / 2) + smoothing_flank
```

Rationale:

- Half the maximum fragment length protects against fragments whose midpoint is
  inside the output interval but whose aligned fragment span reaches a
  blacklist.
- The smoothing flank protects against blacklist-induced artifacts outside the
  output interval bleeding into the retained profile.
- This is a site-set policy, not a fragment filter. Dropped intervals should be
  reported clearly in run statistics.

This should be the default policy when blacklists are supplied. Users can still
pre-clean their interval file manually, but the command should not rely on that
for unbiased profiles.

`--keep-blacklisted-intervals` should disable this interval-level prefiltering
while preserving the existing fragment-level blacklist filter. This is an
escape hatch for deliberate site-set choices, not the recommended path.

## Output Naming

Preserve the current main output name for the selected profile:

```text
<prefix>.midpoint_profiles.npy
```

Write a settings sidecar so binned and smoothed arrays are self-describing:

```text
<prefix>.midpoint_profile_settings.json
```

The sidecar should record the array axes, fragment length bins, original output
interval length, final position bin size, bin aggregation, smoothing method and
parameters, correction flags, and interval blacklist prefilter settings.

Do not add a built-in "both raw and smoothed" output mode. Users who need both
can rerun the command with a different `--smooth` setting and output prefix.

## Decisions

- Smoothing flanks are derived from the selected smoothing method. Future
  methods may use larger internal flanks if needed.
- `--smooth savgol` starts with the same setting as Griffin, currently
  `savgol=165` with fixed polynomial order 3, but Griffin compatibility is
  documented through primitive options rather than a separate preset.
- Smoothing happens at base resolution before final binning.
- The command writes one selected profile output. There is no built-in `both`
  mode.
- The command writes a settings sidecar for interpreting length bins, position
  bins, smoothing, correction flags, and interval blacklist prefiltering.
- Interval-level blacklist prefiltering is enabled by default when blacklists
  are supplied. `--keep-blacklisted-intervals` disables that prefilter while
  preserving fragment-level blacklist filtering.

## Implementation Checklist

### Configuration And CLI

- [x] Add `--bin-size`, default `1`.
- [x] Add `--smooth`, default `raw`.
- [x] Parse `--smooth raw`, `--smooth savgol`, and `--smooth savgol=<odd_bp>`.
- [x] Fix Savitzky-Golay polynomial order at 3. Do not expose order as a CLI
  argument unless a real use case appears.
- [x] Add `--keep-blacklisted-intervals`.
- [x] Validate `--bin-size >= 1`.
- [x] Validate explicit Savitzky-Golay windows:
  - [x] Window is odd.
  - [x] Window fits the output interval after computation flanks are derived.
- [x] For bare `--smooth savgol`, use the documented default window `165`.
- [x] If the smoothing window does not fit the output interval, fail early and
  suggest the largest odd `--smooth savgol=<window>` value that would fit.
- [x] Fail when smoothing is requested for output intervals shorter than 7 bp.
- [x] Include clear config docstrings that explain:
  - [x] Raw profiles remain the default.
  - [x] Savitzky-Golay order is fixed at 3.
  - [x] Smoothing window is in base pairs, not bins.
  - [x] Binning happens after smoothing and flank trimming.
  - [x] `--keep-blacklisted-intervals` disables only interval-level prefiltering,
    not fragment-level blacklist filtering.

### Interval Preparation

- [x] Preserve the original interval start/end as the output span.
- [x] Derive smoothing flanks from the selected smoothing method.
- [x] Expand intervals by the derived flank only for computation.
- [x] Require expanded intervals to fit chromosome bounds so smoothed output
  positions have full filter support.
- [x] Keep output positions indexed relative to the original interval, not the
  expanded interval.
- [x] When blacklists are supplied and `--keep-blacklisted-intervals` is false,
  drop intervals whose original output span plus safety margin overlaps a
  blacklist.
- [x] Compute the interval-level blacklist safety margin as:

```text
ceil(max_fragment_length / 2) + smoothing_flank
```

- [x] Report interval prefiltering counts in run statistics:
  - [x] Loaded intervals after chromosome filtering.
  - [ ] Intervals dropped by chromosome filtering.
  - [x] Intervals dropped by blacklist prefiltering.
  - [x] Intervals retained for counting.
- [x] Fail clearly if no intervals remain after prefiltering.

### Counting

- [x] Keep existing midpoint counting semantics unchanged.
- [x] Count over expanded intervals when smoothing flanks are active.
- [x] Preserve current fragment-level blacklist filtering.
- [x] Preserve GC and scaling weight semantics.
- [x] Preserve sparse tile partial and dense merge behavior where possible.
- [x] Confirm tile fetch halos still use maximum fragment length and are not
  accidentally replaced by smoothing flanks.

### Profile Post-processing

- [x] Add a post-merge profile transform stage before writing.
- [x] Apply smoothing to the final profile tensor, not per input interval.
- [x] Smooth along the position axis independently for every group and length
  bin.
- [x] Trim computation-only flanks after smoothing.
- [x] Apply final binning after smoothing and flank trimming.
- [x] Store binned values as averages.
- [x] Handle a shorter final bin by dividing by the actual covered positions in
  that bin.
- [x] Preserve the current selected-profile output name:

```text
<prefix>.midpoint_profiles.npy
```

- [x] Write `<prefix>.midpoint_profile_settings.json` with the resolved length
  axis, position axis, binning, smoothing, correction flags, and interval
  blacklist prefilter settings.
- [x] Do not add a built-in mode that writes both raw and smoothed profiles.

### Savitzky-Golay Implementation

- [x] Prefer a small internal implementation with runtime-generated
  coefficients unless a maintained dependency is clearly better.
- [x] If adapting MIT-licensed code, keep the original copyright and MIT notice
  in a NOTICE/license file or source header. Not applicable, no external code was
  adapted.
- [x] Credit any adapted implementation in concise module-level documentation.
  Not applicable, no external code was adapted.
- [x] The implementation documentation must explain the scientific basis:
  - [x] Coefficients are the center row of the least-squares polynomial fit over
    centered integer positions.
  - [x] The order-3 smoother preserves constants, linear trends, quadratics, and
    cubics at the window center.
  - [x] For coefficients `c_i` over centered offsets `x_i`, the expected
    moment constraints are:

```text
sum(c_i) = 1
sum(c_i * x_i) = 0
sum(c_i * x_i^2) = 0
sum(c_i * x_i^3) = 0
```

- [x] Include a hand-derived 7 bp order-3 coefficient test as the smallest
  useful example:

```text
[-2, 3, 6, 7, 6, 3, -2] / 21
```

- [x] Add coefficient invariant tests for larger windows, including 165 bp:
  - [x] Symmetry.
  - [x] Sum equals 1 within tolerance.
  - [x] Moments 1, 2, and 3 equal 0 within tolerance.
  - [x] Constant input remains constant.
  - [x] Linear, quadratic, and cubic input are preserved at valid centers.
- [x] Keep mentally derived tests as the scientific proof of correctness.
  Expectations must be derived from the polynomial-preservation properties or
  tiny hand-derived examples. A supplemental hard-coded SciPy coefficient
  regression is allowed as a drift guard because it was explicitly requested,
  but it must not replace the hand-derived and invariant tests.

### Unit Tests

- [x] Test smoothing parser accepts `raw`, `savgol`, and `savgol=<odd_bp>`.
- [x] Test smoothing parser rejects even windows, malformed values, and
  unsupported methods.
- [x] Test smoothing validation rejects output intervals shorter than 7 bp.
- [x] Test bare `savgol` fails early with a suggested explicit smaller window
  when the default window does not fit.
- [x] Test explicit `savgol=<window>` fails with a suggested smaller window when
  it does not fit.
- [x] Test derived smoothing flank is `floor(window / 2)` for odd Savitzky-Golay
  windows.
- [x] Test expanded intervals are trimmed back to the original output span.
- [x] Test output position indexing remains relative to the original interval.
- [x] Test interval-level blacklist prefiltering drops an interval when a
  blacklist overlaps output span plus safety margin.
- [x] Test `--keep-blacklisted-intervals` keeps the interval while fragment-level
  blacklist filtering still removes blacklisted fragments.
- [x] Test prefiltering fails clearly when all intervals are dropped.
- [x] Test smoothing is applied after group aggregation, not per input window.
- [x] Test smoothing is applied independently per group and length bin.
- [x] Test bin-size 1 is identity.
- [x] Test final bin averaging for exact and shorter final bins.
- [x] Test output shapes for raw, smoothed, binned, and smoothed-plus-binned
  profiles.
- [ ] Test run statistics include interval prefiltering counts.
- [x] Add an end-to-end `midpoints` command test that exercises smoothing and
  binning together through counting, merge, post-processing, and output writing.
- [x] Test the settings sidecar records the resolved smoothing and binning
  choices needed to interpret the written `.npy` file.
- [x] Add a margin-composition test where both `ceil(max_fragment_length / 2)`
  and `smoothing_flank` contribute to interval blacklist prefiltering.
- [x] Strengthen polynomial-preservation coverage with a profile-level test at
  non-zero retained positions, not only coefficient moments around offset 0.

All numeric expectations should be written before any test is run. Derive them
by hand from tiny fixtures, coefficient identities, or conservation properties.

### Documentation

- [x] Update `.AI/docs/specs/midpoints_spec.md` only after behavior is
  implemented and finalized.
- [x] Update `.AI/docs/command_pipelines/midpoints.md` to show interval
  prefiltering and post-processing.
- [x] Document the settings sidecar as part of the midpoint output contract.
- [x] Add concise, precise, pedagogically readable CLI docstrings for smoothing,
  binning, and interval blacklist prefiltering.
- [x] Mention in docs where users need to know:
  - [x] Raw full-resolution profiles remain the default.
  - [x] `--smooth savgol` uses order-3 Savitzky-Golay smoothing.
  - [x] Smoothing happens on final profiles, not per input interval.
  - [x] Binning happens after smoothing and trimming.
  - [x] Smoothing flanks must fit chromosome bounds.
  - [x] Blacklisted intervals are prefiltered by default, and
    `--keep-blacklisted-intervals` disables only that prefilter.
  - [x] Griffin-compatible behavior is available through primitive options and a
    guide, not a separate preset.
- [x] Do not rewrite unrelated docs while implementing this feature.

### Verification

- [x] Run `cargo check` after implementation changes.
- [x] Run `cargo check --tests` after adding or changing tests.
- [x] Do not run tests unless explicitly asked. Test expectations must be
  mentally derived before any test execution.

### Final Simplification Pass

- [x] Re-read every changed file before answering.
- [x] Remove unnecessary wrapper functions, types, and indirection.
- [x] Inline helpers that are only used once unless the helper isolates a real
  invariant or substantially improves readability.
- [x] Check that new names are explicit but not verbose.
- [x] Keep the public CLI surface limited to the planned options.
- [x] Keep comments and docstrings tied to behavior that users or maintainers
  need to understand.
- [x] Remove any speculative abstraction for future smoothing methods unless the
  current implementation actually needs it.
- [x] Confirm the final diff is scoped to interval prefiltering, derived flanks,
  profile smoothing, final binning, tests, and necessary documentation.

## Claude review

This review covers the staged diff for `midpoints` smoothing, blacklist
prefiltering, and binning. Every claim was checked against the source.

### Correctness

What checks out:

- **Order-3 Savitzky-Golay center coefficient formula**
  ([smoothing.rs:100-123](../../../src/commands/midpoints/smoothing.rs#L100-L123))
  matches the standard closed form. The denominator
  `(2m + 1) * (4m^2 + 4m - 3)` expands to `(2m-1)(2m+1)(2m+3)`, which is the
  textbook order-2 normalizer. Using the order-2 form for order 3 is justified
  by symmetry: for an odd, symmetric window the odd-order Gram polynomial
  contributions vanish at the center, so order 3 equals order 2 there. The
  hand-derived 7 bp case `[-2, 3, 6, 7, 6, 3, -2] / 21` verified independently.
- **Blacklist safety margin** in
  [midpoints.rs:156-158](../../../src/commands/midpoints/midpoints.rs#L156-L158).
  `(max_fragment_length + 1) / 2` is exactly `ceil(max_fragment_length / 2)`
  for all `u32` inputs (verified for even, odd, and 0).
- **Streaming blacklist pointer** in
  [windows.rs:150-175](../../../src/commands/midpoints/windows.rs#L150-L175)
  is sound: input windows are sorted by start, length is uniform, and the
  margin is constant — so the margin-expanded start and end are both monotonic
  per chromosome. Merged-sorted blacklists let `bl_ptr` advance without
  lookback.
- **Output position centering after smoothing** is correct.
  [postprocess.rs:163-175](../../../src/commands/midpoints/postprocess.rs#L163-L175)
  sets `out[..., 0] = sum_{k=0..W} input[k] * c_k`, whose effective center is
  `input[W/2]`. Counting was done into an interval expanded by exactly `W/2`
  on the left, so `out[..., 0]` lands on the original interval's first base.
  No off-by-one.
- **Final partial bin averaging** divides by actual covered positions
  ([postprocess.rs:184-211](../../../src/commands/midpoints/postprocess.rs#L184-L211)),
  matching the spec.
- **Sort-order preservation after expansion**: a uniform flank is added to
  every retained window, so sort-by-start of the original BED becomes
  sort-by-start of the expanded vector. Downstream tile helpers
  (`precompute_tile_window_spans`, `overlapping_windows_for_tile`) continue
  to work unchanged.
- **`is_identity` fast path**: when `smoothing.is_none() && bin_size == 1`,
  `flanked_length == output_len`, so the merged counting tensor is written
  directly without an extra copy.

What is questionable:

- **Error message lies in raw mode**.
  [windows.rs:112-135](../../../src/commands/midpoints/windows.rs#L112-L135)
  always frames a chromosome-bounds failure as a smoothing error
  ("Cannot smooth interval..." / "Use a smaller --smooth window..."), even
  when `flank == 0`. With raw profiles, a BED row with `end > chrom_len`
  would now produce a misleading smoothing-flavored error. Pre-feature
  behavior did not have this validation at this point in the pipeline.
- **All-or-nothing bound failure vs. silent blacklist drop** is an
  intentional design choice from the plan, but it is worth surfacing in the
  user-facing docs: one telomere-adjacent interval can kill an otherwise
  large run when `--smooth savgol=165` is set. The blacklist prefilter
  silently drops; the chromosome-bounds check bails. Users will likely hit
  this on real TSS/TFBS sets without a clear up-front mitigation.
- **`interval_blacklist_margin` is computed even when no prefilter runs**
  ([midpoints.rs:156-159](../../../src/commands/midpoints/midpoints.rs#L156-L159)).
  Harmless and cheap; just noting that the value is unused when
  `use_blacklist_prefilter` is false.

### Optimality

- **Smoothing is single-threaded and O(G · L · N · W)** in
  [postprocess.rs:147-178](../../../src/commands/midpoints/postprocess.rs#L147-L178).
  For the documented use case (≈1500 groups × ≈100 length bins × ≈2 kb
  output × 165 bp window ≈ 5 × 10^10 multiply-adds) this is not free. The
  inner loop is a textbook target for `rayon::par_iter` over the
  `(group, length_bin)` dimensions — every output cell is independent and
  the coefficient vector is shared. The plan acknowledges that "runtime
  coefficients are acceptable unless profiling shows the final tensor
  transform is a bottleneck", so leaving this serial is defensible; just
  flagging that this is the natural place to parallelize when it shows up
  in profiles.
- **Peak-memory advisory understates reality**.
  [midpoints.rs:189-222](../../../src/commands/midpoints/midpoints.rs#L189-L222)
  warns based on `final_profile_bytes` (the binned/trimmed output). Real
  peak is `counting_tensor + transformed_tensor`, both alive simultaneously
  during the write
  ([midpoints.rs:369-374](../../../src/commands/midpoints/midpoints.rs#L369-L374)).
  For large binning factors the gap is significant — e.g. a 50 GB counting
  tensor with `--bin-size 50` shows a 1 GB warning while peak is ~51 GB.
- **`bin_profile`'s `bin_size == 1` early return is dead code on the actual
  call paths.** `is_identity` catches raw+bin=1, and the smoothing+bin=1
  branch returns `smoothed` directly without calling `bin_profile`. Harmless
  defensive code; could be removed.
- **`ProfileLayout::output_positions` duplicates what `bin_profile` recomputes
  from `positions.div_ceil(bin_width)`**. Same result, two sources of truth.
  Minor.
- The f64 accumulator in `smooth_trimmed_profile` with an f32 store is a
  good numerical choice for 165 bp windows.

### Test coverage

The strong tests:

- Hand-derived 7 bp coefficients
  ([smoothing_tests.rs:101-119](../../../src/commands/midpoints/smoothing_tests.rs#L101-L119))
  with the derivation written out in the comment.
- Moments 0..3 conservation at 165 bp
  ([smoothing_tests.rs:122-156](../../../src/commands/midpoints/smoothing_tests.rs#L122-L156)).
- Polynomial preservation at the center for 21 bp.
- Constant-input preservation per-cell, plus shape, in
  [postprocess_tests.rs:119-148](../../../src/commands/midpoints/postprocess_tests.rs#L119-L148).
- Smoothing-before-binning end-to-end check with hand-derived expectations
  ([postprocess_tests.rs:98-116](../../../src/commands/midpoints/postprocess_tests.rs#L98-L116)).
- Blacklist margin half-open coordinate boundary
  ([windows_tests.rs:122-141](../../../src/commands/midpoints/windows_tests.rs#L122-L141)).

The problems:

- **`order3_coefficients_match_scipy_regression_values`
  ([smoothing_tests.rs:57-98](../../../src/commands/midpoints/smoothing_tests.rs#L57-L98))
  is a direct violation of an explicit plan rule.** The plan checklist line
  "Do not use SciPy, Python, or another implementation as an oracle in tests"
  is marked `[x]` but this test is exactly that: constants generated by
  `scipy.signal.savgol_coeffs(window, 3, use="conv")` and hard-coded. The
  in-test comment acknowledges this and reframes it as a drift guard, but
  drift-guarding via an external oracle is what the rule was forbidding.
  The 5 bp and 9 bp values used here are derivable from the closed-form
  formula in `order3_coefficients` itself; rewriting the test to compute
  expectations from the formula derivation (or from `(2m+1)(2m+3)(2m-1)`
  factorization tables) would satisfy the rule without losing coverage.
- **Polynomial preservation only checked at center offset 0**.
  [smoothing_tests.rs:172-200](../../../src/commands/midpoints/smoothing_tests.rs#L172-L200)
  evaluates each polynomial at offset 0 and asserts it equals 5. By
  construction every polynomial used evaluates to 5 at offset 0, so this
  test is mostly checking moment-0 = 1. Checking preservation at a
  non-zero center (i.e., shift the polynomial origin) would actually
  exercise that the linear/quadratic/cubic terms cancel out.
- **Three plan checkboxes are openly unchecked** in the implementation
  checklist:
  - "Intervals dropped by chromosome filtering" stat is not collected
    ([midpoints.rs:117-179](../../../src/commands/midpoints/midpoints.rs#L117-L179)
    only reports the post-filter count).
  - No test that `--keep-blacklisted-intervals` retains the interval
    *while* fragment-level filtering still removes the fragment. The
    closest is `keep_blacklisted_intervals_disables_interval_prefilter`,
    which only exercises the prefilter half. The integration test in
    `test_profile_groups_command.rs` was *defanged* by setting
    `keep_blacklisted_intervals=true` — that is the right fix to keep the
    old test meaningful, but it leaves zero coverage of "both filters
    active at once with a prefilter-surviving but fragment-overlapping
    blacklist".
  - No test that prefiltering fails clearly when every interval is dropped
    (the `ensure!` at
    [midpoints.rs:168-171](../../../src/commands/midpoints/midpoints.rs#L168-L171)
    is unexercised).
  - No test asserting interval prefilter counts appear in run statistics
    output.
- **No end-to-end integration test** exercises the new path (smoothing or
  binning) through `cfdna midpoints` against a small BAM/BED. Unit tests
  cover the pieces, but nothing wires expanded counting → merge →
  postprocess → write together with a known input.
- **The `prepare_count_windows_drops_blacklist_margin_overlaps` test**
  passes `smoothing_flank = 0` and `blacklist_margin = 15`, so it does not
  cover the case where the margin is constructed from
  `ceil(max_fragment_length/2) + smoothing_flank` with both terms nonzero.
  A targeted test for that composition would protect against future
  refactors of the margin formula in `run()`.

### Critique of the approach

- **`MidpointSmoothing::Raw` as an enum variant with `Default` is fine, but
  `--smooth raw` is also accepted as user input.** This is documented as
  intentional in the parser. It just means `--help` shows the flag without
  a visible default, and round-tripping configs from tests works. No issue,
  noted because it might surprise reviewers.
- **The `ProfileLayout::resolve` function does both layout math and
  user-facing validation.** It is the right shape for the job, but it
  emits CLI-suggestion strings (`--smooth savgol=99`) from a layer that
  has no other CLI awareness. This is a small layering smell — if the
  layout struct ever gets reused outside the `midpoints` CLI, the
  suggestion text will look strange.
- **Failing fast on smoothing-window-too-large is good UX**, and the
  suggested largest fitting odd window in the error is genuinely helpful.
- **The chromosome-edge expansion error and the blacklist prefilter have
  inconsistent semantics** (bail vs. silent drop) for what is conceptually
  the same situation: "this interval cannot produce a clean smoothed
  profile". The plan picked these explicitly, so this is design intent
  rather than a bug — just worth mentioning in user docs because a single
  bad interval can cost a long run.
- **The post-merge transform allocates a new tensor of `output_len`** even
  when only binning is requested. For pure binning that allocation is
  about `output_len / bin_size`-fold smaller than the counting tensor, so
  it is fine. For smoothing + binning, you transiently hold counting +
  smoothed + (during construction of) binned. Not a regression but a
  factor to remember when sizing memory.
- **Loop ordering in `smooth_trimmed_profile`** iterates
  `(group, length_bin, output_pos, offset)`. For `ndarray` row-major
  storage with shape `(G, L, P)`, the innermost stride is on `P`, so the
  current iteration order is already cache-friendly along the position
  axis. Good.

### Summary

Math and counting semantics are correct. The plan's explicit "no SciPy
oracle" rule is violated by one test, despite the checklist being marked
done; that one is worth fixing before this lands, both to respect the rule
and because it is the single most likely-to-be-flagged item in code
review. A handful of plan checkboxes are honestly left `[ ]` and should
either be addressed or removed from scope. The remaining items above are
optimization and UX polish rather than blocking issues.

## Follow-up After Review

The SciPy-regression criticism is intentionally not accepted. The project keeps
the SciPy constants as a supplemental regression check because this was an
explicit design request. The correctness proof must still come from the
hand-derived coefficient case and the polynomial moment/preservation tests.

Changes to make:

- [x] Fix chromosome-bound error messages in raw mode. When `smoothing_flank ==
  0`, an interval outside chromosome bounds should report an invalid BED/contig
  bound, not "Cannot smooth" or "Use a smaller --smooth window".
- [x] Document the chromosome-edge smoothing behavior where users will see it:
  smoothing flanks must fit chromosome bounds, so telomere-adjacent intervals can
  fail the run unless the smoothing window is smaller or those intervals are
  removed.
- [x] Add an integration test for `--keep-blacklisted-intervals` where the
  interval survives prefiltering but fragment-level blacklist filtering still
  removes the fragment.
- [x] Add an integration test for the "all intervals dropped by blacklist
  prefiltering" error.
- [x] Add an end-to-end smoothing-plus-binning command test with tiny fixtures
  and hand-derived expected output.
- [x] Add a test that protects the composed blacklist margin:
  `ceil(max_fragment_length / 2) + smoothing_flank`.
- [x] Strengthen polynomial-preservation testing at the profile level by
  smoothing a quadratic or cubic profile and checking retained non-zero
  positions.
- [x] Remove dead defensive code in `bin_profile` if it remains unreachable after
  the follow-up edits.

Deferred or no action for now:

- Do not remove the SciPy regression test.
- Do not change the dense output-size warning into a peak-memory warning. That
  warning exists to tell users when the final output file/tensor is going to be
  very large. Peak memory has more contributors and will fail naturally if the
  system cannot allocate it.
- Do not parallelize `smooth_trimmed_profile` until profiling shows the
  post-processing stage is a real bottleneck.
- Do not move CLI-flavored smoothing suggestions out of `ProfileLayout::resolve`
  unless the layout code becomes shared outside the command.
- Do not add "intervals dropped by chromosome filtering" unless the BED loader
  can report the pre-filter total without adding another large interval copy.
