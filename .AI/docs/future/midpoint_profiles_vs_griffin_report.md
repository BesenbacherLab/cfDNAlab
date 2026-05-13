# Midpoint Profiles Compared With Griffin

Analysis by Codex 5.5 on May 13th 2026. Last edit: May 13th 2026

Status: future analysis, not a finalized user-facing specification.

This report compares the planned midpoint-profile command with Griffin's nucleosome
profiling workflow. It assumes the next midpoint-profile work adds strand-aware
orientation for site-centered inputs, including TSS sites. It also assumes the
smoothing documentation clearly states that `savgol=165` is a nucleosome-scale
preset, not a universal profile-smoothing default.

The goal is not to make a Griffin clone. The goal is to keep the midpoint command a
general profile generator while making it possible to explain which settings are
Griffin-like and which results are intentionally different.

## Source Points Checked

- Griffin upstream repository:
  `https://github.com/adoebley/Griffin`

- Griffin coverage counting:
  `scripts/griffin_coverage.py`

- Griffin site merging, binning, masking, smoothing, normalization, and feature
  calculation:
  `scripts/griffin_merge_sites.py`

- Griffin site orientation and chromosome-edge clamping:
  `scripts/griffin_functions.py`

- Griffin nucleosome profiling defaults:
  `snakemakes/griffin_nucleosome_profiling/config/config.yaml`

- Current midpoint implementation:
  `src/commands/midpoints/`

## Smoothing Window Interpretation

Griffin does not simply apply a base-resolution 165 bp Savitzky-Golay filter to
the final profile.

In `griffin_merge_sites.py`, Griffin first sums base-resolution values into
`step`-sized bins. With the nucleosome profiling defaults, `step = 15`. It then
rounds `smoothing_length` to a multiple of `step` and converts the smoothing
length to a Savitzky-Golay window measured in binned positions. With the default
`smoothing_length = 165`, this becomes an 11-bin Savitzky-Golay window over
15 bp bins.

So Griffin's default keeps the same nominal physical support, roughly 165 bp, but
it smooths a lower-resolution binned signal. That is not mathematically identical
to smoothing the full-resolution profile with a 165 bp kernel and only then
binning.

For our command, a 165 bp order-3 Savitzky-Golay preset remains the natural
Griffin-like nucleosome-scale setting. There is no reason to change that number
just because we smooth at base resolution. It should still not be treated as
universally optimal for TFBS, TSS, very short windows, or non-nucleosomal fragment
classes. The clean policy is:

- raw profiles remain the default output

- `--smooth savgol` is a documented nucleosome-scale preset

- `--smooth savgol=<odd_bp>` lets users tune the biological scale

- smoothing happens at full resolution before any final binning

This makes the default conservative. It avoids irreversible smoothing unless the
user asks for it, while keeping a Griffin-familiar smoothing scale available.

## Where Our Approach Is Better For A General Tool

### Exact Output Windows

Our intervals define the output window exactly. Griffin has separate normalization,
save, center, and FFT windows, and it rounds profile windows inward to multiples of
`step`.

For a general command, exact interval-defined output is better. Users can prepare
TFBS, enhancer, TSS, promoter, or arbitrary windows and know that the array axis
matches their input definition.

### Base Resolution Until The End

Our planned processing keeps profiles at base resolution through counting and
smoothing, then optionally bins at the end. Griffin bins first and smooths the
binned signal.

For a general tool, base resolution until the end is the better default. It avoids
making the smoothing result depend on an early binning choice, and it lets users
choose final binning as an output-size or downstream-analysis decision.

### Raw By Default

Griffin's final profiles are smoothed by default. Our command writes raw profiles
unless smoothing is requested.

Raw by default is better for a reusable profile generator. Smoothing is a modeling
choice, not a counting requirement. Users can reproduce smoothed outputs from raw
profiles, but they cannot reconstruct raw high-frequency signal from a smoothed
profile.

### Smoothing Flanks Instead Of Edge Modes

Our smoothing counts extra flank positions, smooths with real observed context,
and trims back to the requested output interval. Griffin uses SciPy's default
Savitzky-Golay edge behavior on the binned profile.

For biological aggregate profiles, computation-only flanks are preferable when
available. They make the reported edge positions depend on nearby observed signal
rather than an edge interpolation rule.

### Balanced Even-Fragment Midpoints

Griffin places even-length fragment midpoints at `floor((start + end) / 2)`, which
chooses the right center for half-open intervals. Our midpoint code assigns
even-length fragments reproducibly to left or right center using a coordinate-
derived tie-break.

Our behavior is better for unbiased midpoint profiles because it avoids a
systematic one-base shift for the large fraction of even-length fragments.

### Length-Resolved Output

Griffin filters to one size range for a profile run. Our midpoint array keeps a
length-bin axis. A single run can retain broad fragment classes, narrow
nucleosome-sized classes, and finer fragment-length structure.

That is more useful for exploratory fragmentomics and for methods where the
length distribution is part of the signal.

### Stricter Input Validation

Griffin clamps windows at chromosome edges and pads missing positions with `NaN`.
Our current policy is to fail if the requested interval or derived smoothing flank
does not fit the chromosome.

For a command that consumes user-defined intervals, failing early is better. A BED
row outside the chromosome often means the wrong assembly or a coordinate bug.
Silently padding can hide that error.

### Blacklist Policy Is Less Edge-Biased

Griffin masks excluded bins during site merging. Our blacklist inputs can include
low- or zero-mappability regions, ENCODE exclusions, gaps, centromeres, and other
regions that should not contribute to profiles. The planned default removes
intervals whose output span plus fragment and smoothing safety margin overlaps a
blacklisted region, unless users explicitly keep them.

For aggregate profile shape, the interval-level prefilter is better. It avoids
one-sided fragment loss near profile edges, which can otherwise produce artificial
dips or shoulders.

### Metadata Sidecar For Arrays

Our `.npy` profile values are accompanied by a settings JSON that records axis
order, length bins, bin size, smoothing, flanks, and blacklist prefiltering.

That is better for array reuse. Griffin TSVs are human-readable and include useful
metadata columns, but the final profile values are already transformed by Griffin's
workflow choices.

### No Mappability Correction In The Core Profile

Griffin disables mappability correction by default in the Snakemake workflow, and
its own config says it was not recommended because it did not improve signals.
Griffin still defaults to excluding zero-mappability bins during merge.

For our midpoint command, not doing mappability correction in the core profile is
the better default. Mappability correction can easily become a model that invents
coverage rather than a transparent count adjustment.

This does not mean ignoring mappability. Mappability exclusion belongs in the
blacklist layer. That keeps the policy explicit: low-mappability regions are
removed or cause nearby intervals to be removed, rather than receiving inflated
coverage weights.

## Where Griffin Is Better Or More Complete

### Mature Site Orientation

Griffin already handles forward, reverse, and undirected sites and reverses
reverse-strand profiles before merging.

This is the main real feature gap for TSS-style analyses. Once our command has
strand-aware orientation, this difference should mostly close.

### End-To-End Griffin Outputs

Griffin produces final normalized TSVs, plots, central coverage, FFT amplitude,
and metadata in one workflow.

Our command currently focuses on producing midpoint profile arrays. That is the
right scope for the core command, but Griffin is in some ways currently more complete 
as a ready-made nucleosome-profiling analysis pipeline.

### Built-In Partial-Site Masking

Griffin's `NaN` masking allows a site to contribute to some bins but not others.
This can be convenient when retaining as many sites as possible matters more than
keeping a uniform site set across the whole profile.

Our interval-level prefilter is better for avoiding edge bias in aggregate shape,
but it is stricter. That is an intentional trade-off, not a universal win.

### Direct Griffin Feature Compatibility

Griffin's central coverage and FFT amplitude are available immediately in the
final TSV.

For us, these should stay downstream for now. A future Griffin guide can document
the equivalent calculations from our profile arrays. Adding them to the command
should wait until there is clear user demand, because features are analysis
choices rather than profile-counting requirements.

## Differences That Are Mostly Policy Choices

### Normalization

Griffin normalizes final profiles to mean 1 over the normalization window, with an
optional CNA-style per-site normalization step. Our command can apply GC and
genome scaling factors but does not force Griffin-style profile normalization.

This is a deliberate separation. The midpoint command should count and optionally
apply explicit correction weights. Shape normalization can stay downstream unless
we add a clearly named profile-normalization option.

### Output Scale

Griffin averages sites by default. Our groups are summed across intervals.

The summed representation is more primitive and therefore more reusable. The
trade-off is that users need the group index and interval counts to interpret
profile scale.

For midpoint profiles, `group_index.tsv` exposes `eligible_intervals`, the
number of output intervals retained per group after chromosome filtering and
interval-level blacklist prefiltering. That count is independent of fragment
overlap. An eligible interval still counts even if no fragment midpoint overlaps
it.

For groups with at least one eligible interval, this lets users convert summed
group profiles to mean profiles without reconstructing the command's internal
interval filtering decisions.

### Outlier Handling

Griffin masks bins with extreme coverage outliers before smoothing and
normalization.

This may be useful in a closed workflow, but it is too opinionated for the core
midpoint profile command. Outlier handling should remain a downstream analysis or
a separately named option.

We have plans of adding an outlier scaling approach (similar to fragment-count-weights
but on positional level) that could also reduce outlier issues.

### Chromosome Edges

Griffin pads truncated edge windows with `NaN`. Our command should fail when
intervals or required flanks fall outside the chromosome.

This is stricter, but it makes mistakes obvious. Users who intentionally want
edge-padding can preprocess intervals or use a future explicit mode.

## Remaining Work For A Griffin-Like Guide

- Add strand-aware interval orientation that works naturally for TSS and other
  directional sites.

- Document that `savgol=165` is the Griffin-like nucleosome-scale preset. Any
  future tuning should be driven by biology and interval class, not by a mistaken
  idea that Griffin's 165 bp setting used a different physical support.

- Document a Griffin-like recipe:

  - use nucleosome-sized fragment bins such as `100 201`

  - use MAPQ 20 if reproducing Griffin defaults

  - use `--smooth savgol` for a Griffin-familiar smoothing scale

  - use `--bin-size 15` when comparing to Griffin-style output tables

  - normalize downstream to mean 1 over the intended comparison window

  - compute central coverage and FFT amplitude downstream if needed

- Be explicit that exact Griffin reproduction would still require Griffin's
  right-center even-fragment convention, bin-first smoothing, `NaN` masks,
  site-list averaging, and final TSV feature calculations.

## Bottom Line

After adding strand-aware orientation and making the smoothing preset clear, the
midpoint command is better suited than Griffin for general midpoint-profile
generation. It is less opinionated, keeps higher-resolution information longer,
has clearer interval semantics, and avoids several hidden transformations.

Griffin remains useful as a reference implementation for a specific
nucleosome-profiling workflow. It should inform our guide and defaults, but it
should not define the command's core behavior.
