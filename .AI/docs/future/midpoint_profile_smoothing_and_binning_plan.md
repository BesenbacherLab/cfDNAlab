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
--smooth savgol=165,3
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
- `--smooth savgol=165,3` means a 165 bp Savitzky-Golay window with polynomial
  order 3. This is Griffin-like, but it is not the only supported smoothing
  scale.
- `--keep-blacklisted-intervals` disables default interval-level blacklist
  prefiltering for users who intentionally prepared or want to keep the original
  site set.

The exact spelling can change during implementation. The important model is one
smoothing option with method-specific settings, not several independent CLI
flags for every smoothing parameter.

Smoothing flanks are derived, not user-facing:

```text
raw profiles:      smoothing_flank = 0
smoothed profiles: smoothing_flank = ceil(smoothing_window_bp / 2)
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
means `savgol=165,3` always means 165 bp, not 11 bins under one setting and
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
with `savgol=165,3` as the initial default for `--smooth savgol` until tuning
shows a better general default. Griffin-compatible behavior can be documented in
a guide through the primitive options rather than encoded as a separate preset.

Validation should fail clearly when a requested smoothing setting cannot be
applied, for example:

- The smoothing window is even.
- The smoothing window is shorter than `poly_order + 2`.
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

If flanking is only for smoothing, the mathematically necessary flank is about
half the smoothing window on each side. For `savgol=165,3`, that is roughly
82 bp.

The default should be derived rather than user-specified:

```text
raw profiles:      flank_size = 0
smoothed profiles: flank_size = ceil(smoothing_window_bp / 2)
```

This keeps the CLI smaller and makes `--smooth savgol=165,3` fully define the
smoothing behavior. If users later need extra context for a different reason,
that should be designed as a separate feature instead of preemptively adding a
generic flank argument.

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

Do not add a built-in "both raw and smoothed" output mode. Users who need both
can rerun the command with a different `--smooth` setting and output prefix.

## Decisions

- Smoothing flanks are derived from the selected smoothing method. Future
  methods may use larger internal flanks if needed.
- `--smooth savgol` starts with the same setting as Griffin, currently
  `savgol=165,3`, but Griffin compatibility is documented through primitive
  options rather than a separate preset.
- Smoothing happens at base resolution before final binning.
- The command writes one selected profile output. There is no built-in `both`
  mode.
- Interval-level blacklist prefiltering is enabled by default when blacklists
  are supplied. `--keep-blacklisted-intervals` disables that prefilter while
  preserving fragment-level blacklist filtering.

## Implementation Steps

1. Add config parsing for smoothing mode and bin size without changing current
   defaults.
2. Derive computation-only flanks from the smoothing mode while preserving the
   original output span.
3. Add default interval-level blacklist prefiltering with clear run statistics.
4. Add a profile transform stage after merge and before write.
5. Add final binning as an output transform.
6. Preserve selected-profile output naming.
7. Add focused tests for boundary trimming, blacklist interval prefiltering,
   Savitzky-Golay window validation, last-bin averaging, and output shape.
8. Distill implemented behavior into `.AI/docs/specs/midpoints_spec.md` only
   after the design is finalized.
