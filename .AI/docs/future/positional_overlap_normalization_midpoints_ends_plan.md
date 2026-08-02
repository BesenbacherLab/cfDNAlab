# Positional Overlap-Length Normalization for Midpoints and End Motifs

Status: future experimental plan. This document does not describe current command behavior or an
established fragmentomics method.

## Purpose

Evaluate whether the sample-specific average overlapping fragment length model can reduce dataset
shift in midpoint profiles and selected end-motif features.

The model is already fitted from raw positional fragment coverage. Its current validated
application is `fcoverage`, where positional lookup weights are averaged across each fragment's
counted span before the fragment contributes coverage. Midpoints and ends require different
application rules because their observations occur at individual genomic positions.

The research question is not whether LIONHEART applied these transformations. It did not. The
question is whether a positional normalization learned without phenotype labels can improve
fragmentomics model performance on datasets, laboratories, or library protocols that were absent
from model training.

This work must remain opt-in and experimental until external-dataset testing demonstrates a useful
tradeoff between reduced dataset effects and preserved biological signal.

## Motivation

For a covered genomic position `x`, define:

```text
depth(x) = number of accepted fragments covering x

average_length(x) =
    sum of full directional fragment lengths covering x / depth(x)

overlap_weight(x) = fitted_lookup(average_length(x))
```

The fitted lookup homogenizes the sample-specific relationship between positional coverage and
average overlapping fragment length toward the model's symmetric 166 bp target. Fragmentation,
library preparation, and sequencing differences can make this relationship dataset-specific.
Those differences can also affect where midpoint and fragment-end events are observed.

For machine learning, the intended benefit is removal of an unstable shortcut. A classifier may
otherwise use dataset-specific fragment length and coverage structure that predicts the label in
its training cohorts but fails in an external cohort.

The cost is equally important. If phenotype-associated fragmentation changes the same local
average-length context, this normalization can remove useful biological signal. Making samples
more alike is not sufficient evidence of success. The transformation is useful only when it
improves held-out dataset performance or calibration while preserving the feature structure needed
by the model.

## Scientific Boundary

The overlap lookup was fitted with positional coverage as the observed signal. Applying it to
midpoint or end counts is therefore a transfer of a learned positional normalization, not a fitted
inverse observation probability for those event types.

That distinction has practical consequences:

- The transformation has a clear deterministic meaning at an event position.
- It is not guaranteed to flatten or otherwise calibrate midpoint or end density conditional on
  average overlapping fragment length.
- Evaluation must focus on domain generalization rather than agreement with an assumed unbiased
  midpoint or end estimator.
- Standard command behavior must remain unchanged when no overlap model is supplied.

Do not describe these outputs as LIONHEART midpoint or end features. Describe them as
overlap-length-normalized midpoint profiles or overlap-length-normalized end-motif counts.

## Shared Counting Invariants

The following invariants apply to both commands.

- Positional overlap context is calculated from raw accepted fragments before GC correction,
  genomic scaling, motif filtering, profile smoothing, or downstream normalization.
- Paired fragment length remains `forward.pos` to `reverse.reference_end`. Unpaired fragment length
  remains `read.pos` to `read.reference_end`.
- The positional context must follow the model package's inter-mate-gap behavior and must preserve
  mapped-reference segments around deletions and skipped regions.
- Lookup uses the package's existing clipped length-bin behavior. Do not add lower-resolution
  arithmetic or float conversion merely to resemble another implementation.
- The model package must be loaded and structurally validated before BAM processing begins.
- A mismatch that changes positional context cannot be ignored silently. Initial experiments
  should require matching fragment length limits, mapping quality, paired or unpaired mode, proper
  pair policy, and inter-mate-gap behavior.
- The fragments defining positional context and the events counted into the feature output are
  related but distinct populations. Command-specific event filters must not retroactively change
  `depth(x)` or `average_length(x)`.
- Blacklisted positions are not valid lookup positions. The initial implementation must choose and
  document an explicit skip or error policy for a surviving event whose lookup position is
  blacklisted. It must not silently substitute a neutral weight.
- All enabled multiplicative weights must be finite and non-negative. A material negative weight
  is an error.
- GC, genomic scaling, and overlap-length normalization remain independently optional. Their
  multiplication order does not change the result, but their coordinate definitions must remain
  distinct and documented.

The requirement to separate context fragments from feature events is particularly important for
`ends`. A motif can be removed by clipping, indel, base-quality, missing-reference, or motifs-file
selection logic. That removal must not make the same physical fragment disappear from the raw
overlap context used to weight neighboring events.

## Midpoint Application

### Definition

For a fragment assigned the deterministic counted midpoint `m`, define:

```text
midpoint_overlap_weight = overlap_weight(m)

final_midpoint_weight =
    midpoint_overlap_weight
    * fragment_gc_weight
    * existing_genomic_scaling_weight
```

The midpoint event receives the positional lookup at its selected midpoint. Do not reuse the
current fcoverage scalar that averages positional weights across the complete fragment span.

Even-length fragments must use the command's existing reproducible left-or-right midpoint choice
before the positional lookup. Window assignment and profile position must use that same selected
midpoint.

The positional weight is applied before grouped-window aggregation, smoothing, flank trimming, and
final position binning. Post-processing therefore acts on the weighted midpoint profile exactly as
it currently acts on GC- or scaling-weighted counts.

### Interpretation

The corrected profile estimates midpoint event intensity after reweighting genomic positions by
their local average overlapping fragment length context.

At a fixed genomic position, all midpoint events receive the same overlap weight. Their relative
fragment length composition at that exact position is unchanged. Across different positions and
different input intervals, however, weights vary. Aggregated length-bin profiles can therefore
change.

At low depth, the counted fragment contributes to its own `average_length(m)`. At depth one, the
average is the fragment's own length. This self-dependence is also present in the positional model
but becomes relevant when interpreting length-resolved midpoint output. Do not introduce a
leave-one-out variant in the initial experiment because it would use a different quantity from the
fitted model and is undefined at depth one.

### Expected Use

The most plausible use is a grouped site-centered profile where the downstream signal is expected
to reside in positional shape, occupancy, or a selected set of fragment length strata, while
sample-wide differences in average overlapping fragment length are treated as nuisance variation.

The correction is less attractive when the downstream model is intended to use global fragment
length shifts. In that case, forcing samples toward the 166 bp target can remove part of the
intended signal.

## End-Motif Application

### Definition

Use independent positional weights for the two aligned terminal bases of a fragment spanning
`[start, end)`:

```text
left_overlap_weight = overlap_weight(start)
right_overlap_weight = overlap_weight(end - 1)
```

For each end that survives motif extraction and window selection:

```text
final_end_weight =
    assignment_overlap_mass
    * end_overlap_weight
    * fragment_gc_weight
    * existing_genomic_scaling_weight
```

The two ends of the same fragment may receive different overlap weights. A shared fragment scalar
would answer a different question and must not be used as a shortcut.

The end-specific overlap weight does not depend on which output window receives the motif. Existing
window assignment mass and genomic scaling remain row-specific where they are row-specific today.
For example, `count-overlap` can produce different assignment and scaling factors for different
rows while reusing the same left or right positional overlap weight for that end.

### Shifted Boundaries

The model defines context on aligned reference coverage. The initial implementation should reject
the combination of positional overlap normalization with
`--clip-strategy include-at-shifted-boundary`.

The shifted motif boundary can lie outside the aligned fragment and may have no defined overlap
context. Looking up the nearest aligned base or the shifted boundary would introduce a new rule
without validation. The aligned `skip`, `aligned`, and `include-at-aligned-boundary` modes have
unambiguous terminal-base coordinates and are suitable for the initial experiment.

### Interpretation

The weighted count for a motif is the sum of its end events after positions have been reweighted by
local average overlapping fragment length context. For downstream motif proportions, the overlap
weight must be applied to counts before dividing by the total weighted motif mass:

```text
weighted_proportion(motif) =
    weighted_count(motif) / sum of weighted counts across the chosen motif set
```

This can reduce dataset differences caused by libraries sampling different genomic contexts. It
can also change motif composition because sequence context, chromatin state, copy number, and
fragment length are geographically structured. That change is the proposed normalization effect,
not an implementation artifact, and must be evaluated rather than assumed beneficial.

Reference k-mer correction remains a separate downstream operation. If both are used, first form
the overlap-weighted observed end-motif counts, then apply the chosen reference denominator. Do not
alter reference k-mer frequencies with a sample-derived overlap model.

## Positional Inference Design

The existing `OverlappingLengthWeightIterator` is specialized for fcoverage. It delays fragments
until their complete counted spans have been finalized, then stores the span-average lookup weight
on `FragmentWithSegments`.

Do not add midpoint and end options by reading that span-average scalar. Refactor the shared
rolling calculation so it can answer two explicit query types:

- Mean lookup weight over counted segments for fcoverage
- Lookup weight at a requested covered position for midpoint and end events

The shared component should own raw depth and fragment length sum events, safe-position
finalization, blacklist eligibility, bounded prefix storage, and package lookup. Command-specific
code should own midpoint selection, resolved ends, window assignment, motif extraction, GC, and
genomic scaling.

The design must allow a context fragment to be ingested even when no feature event is eventually
returned for it. Possible implementation approaches are:

1. Extend the paired-fragment construction path to return raw overlap-context segments together
   with an optional command-specific event representation.
2. Use a separate context pass and feature pass per tile.

Prefer a single BAM pass if it can keep the code understandable. A second pass is acceptable for an
initial research implementation if it substantially reduces correctness risk. Measure memory and
runtime before selecting the more complex design. Do not merge motif-specific filtering into the
context fragment filter to avoid the second pass.

Tile fetch spans must include every fragment that can affect a queried position. Targeted BED fetch
narrowing must retain this overlap-context halo even when only a small number of midpoint or end
events can enter output rows. Derive and test the command-specific halo rather than assuming the
existing counting halo is sufficient. A conservative halo equal to twice the maximum fragment
length is acceptable for the initial implementation if memory remains bounded.

## Configuration And Metadata

Use the existing option name:

```text
--overlap-length-file <model.zarr>
```

The option remains conditional on `cmd_overlap_length_model`. Do not add another Cargo feature for
the experimental applications unless the existing dependency structure causes a concrete build
problem.

CLI help must state:

- The transformation is experimental for midpoint and end-motif features.
- The lookup is evaluated at the midpoint or aligned terminal base, respectively.
- Raw positional overlap context is calculated before GC and genomic scaling.
- The initial supported fragment and filtering settings must match the package.
- Shifted-boundary end motifs are unsupported initially.

Settings outputs must record more than a generic correction flag. At minimum record:

```text
overlap_length_normalization_used
overlap_length_application = midpoint_position | aligned_end_position
overlap_length_package_schema_version
overlap_length_target_mean_fragment_length
overlap_length_package_minimum_fragment_length
overlap_length_package_maximum_fragment_length
```

Do not place the model path in the public Zarr schema merely because it was supplied on the command
line. The settings JSON may record it if that is consistent with the command's existing path
provenance policy. Persisted scientific meaning must not depend on the path remaining valid.

The main midpoint and ends Zarr array schemas do not need to change solely because count values are
weighted differently. Existing `f32` midpoint output and `f64` sparse end counts remain their
current storage types.

## Correctness Tests

Write hand-derived expectations before running tests.

### Shared Positional Inference

- Construct a small fragment set where adjacent positions have different raw depth and average
  overlapping fragment lengths. Assert exact lookup-bin selection at each queried position.
- Assert that a point query differs from a whole-fragment mean when the fragment crosses multiple
  lookup regions.
- Assert equivalent values across tile boundaries and prefix chunk sizes.
- Assert that deletions, skipped regions, and optional inter-mate gaps affect positional context in
  the same way as model fitting.
- Assert that a context fragment with no returned feature event still affects later point queries.
- Assert errors for queries before finalization, missing package state, invalid weights, and
  incompatible context settings.
- Assert the chosen explicit behavior for blacklisted query positions.

### Midpoints

- Assert that odd-length fragments use the exact center base for lookup.
- Assert that even-length fragments use the same reproducible selected midpoint for lookup and
  profile placement.
- Assert that all length bins at the same genomic midpoint receive the same positional overlap
  weight.
- Assert multiplication with GC and genomic scaling using nontrivial hand-derived factors.
- Assert that smoothing and final position binning operate on weighted counts without changing
  existing aggregation rules.
- Assert that output without an overlap package is byte-for-byte or value-for-value identical to
  the existing path where deterministic serialization allows it.
- Add CLI round-trip and settings JSON assertions.

### Ends

- Assert left lookup at `start` and right lookup at `end - 1`.
- Construct a fragment whose two ends have different positional lookup values and assert separate
  motif weights.
- Assert that dropping a motif at an end does not remove the fragment from positional context.
- Assert each window assignment mode preserves its existing assignment mass while adding the
  correct end-specific overlap factor.
- Assert multiplication with GC and genomic scaling using nontrivial hand-derived factors.
- Assert an explicit error for `include-at-shifted-boundary` with overlap normalization.
- Assert that motifs-file selection changes counted targets but not positional context.
- Assert unchanged output when no overlap package is supplied.
- Add CLI round-trip and settings JSON assertions.

Public command behavior and output metadata belong in integration tests under `tests/`. Private
rolling inference and event-weight logic belong in sibling `*_tests.rs` files without widening the
public API.

## Domain-Generalization Evaluation

Do not promote either application based on visual similarity, within-dataset cross-validation, or
reduced dataset separation alone.

For every feature family, generate paired outputs from the same samples:

```text
existing feature path
existing feature path plus positional overlap normalization
```

When GC correction or genomic scaling is part of the intended downstream workflow, evaluate a
small fixed factorial comparison rather than changing several corrections at once. At minimum,
compare the established pipeline with and without positional overlap normalization.

### Dataset Splits

- Hold out complete datasets, laboratories, or library protocols.
- Fit predictive model hyperparameters without using the held-out dataset.
- Fit each sample's overlap model from that sample only, without phenotype labels.
- Keep any feature selection and downstream normalization inside the training split.
- Repeat with several held-out datasets when available. A gain on a single cohort is not sufficient
  evidence of generalization.

### Outcomes

Measure:

- External discrimination metrics appropriate to the task
- Calibration on each held-out dataset
- Technical replicate and downsampling concordance
- Dataset-of-origin predictability from the final features
- Association between final features and overlap-model fit parameters
- Weight distributions and the fraction of feature mass dominated by extreme weights

For midpoints, also inspect profile shape, amplitude, and each fragment length bin separately. For
ends, inspect the selected feature values, complete motif proportions where available, and whether
the weighting concentrates counts in a small set of genomic regions.

### Graduation Criterion

An experimental application may become documented supported behavior when it reproducibly:

- Improves or preserves performance on multiple held-out datasets
- Improves calibration or technical reproducibility
- Reduces dataset-specific structure without collapsing biologically expected feature variation
- Avoids unstable dependence on a small number of extreme positional weights
- Has a clear supported configuration and complete provenance metadata

Failure to improve the available end-motif features is evidence against that application for those
features. It is not evidence that the positional implementation is incorrect. Likewise, success on
selected end motifs does not automatically justify enabling the correction for every motif output.

## Proposed Order

1. Extract and test a positional point-lookup capability without changing current fcoverage
   results.
2. Implement end-specific weighting first because selected end-motif features already have a
   downstream evaluation path.
3. Run leave-dataset-out comparisons on those fixed end-motif features.
4. Keep the midpoint design documented but defer its research conclusion until a concrete midpoint
   modeling setup exists.
5. Implement midpoint weighting when that evaluation can be run immediately.
6. Promote neither option to default behavior without the graduation evidence above.

## Non-Goals

- Do not change the established default output of `midpoints` or `ends`.
- Do not apply the model to fragment length distributions.
- Do not refit the overlap model from GC-corrected or genomically scaled coverage.
- Do not make GC or scaling commands depend on an overlap model.
- Do not change reference k-mer frequencies with a sample-derived model.
- Do not add a generic universal fragment weight and reuse it for all feature types.
- Do not add shifted-boundary behavior until its lookup coordinate has a defensible definition and
  validation data.

## Open Decisions Before Implementation

- Whether the first implementation should use a combined single-pass fragment representation or a
  separate positional-context pass.
- Whether corrected commands should require exact model-package filter equality or initially allow
  a narrowly defined subset relationship.
- Whether a surviving event at a blacklisted lookup position should be skipped or rejected. Neutral
  weighting is not acceptable without an explicit scientific justification.
- Whether package paths belong in command settings JSON or only package-derived scientific
  attributes should be stored.
- Which existing selected end-motif features and external datasets define the first fixed
  evaluation benchmark.
