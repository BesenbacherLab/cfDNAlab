# Fragment-level nuisance outlier handling plan

Date: 2026-08-10

Status: Future design. This is not an implementation specification yet.

This plan develops a fragment-level interpretation of sample-specific coverage
outliers. It supersedes the positional application direction in the earlier
[`cfdna outliers` research spec](outliers_command_spec.md) without replacing
that document or its broader review of possible detection methods.

## Modeling assumptions

### The nuisance process acts on fragments

The primary model is that an extreme coverage region may contain an excess
population of nuisance fragments. Examples include fragments that were mapped
to the wrong reference location or fragments that are overrepresented after
ordinary duplicate filtering.

Under this model, a nuisance fragment is not problematic only at the positions
where aggregate coverage crosses an outlier threshold. Its complete aligned
fragment span and all features derived from it are affected.

The correction should therefore attach a single keep weight to each fragment
during downstream processing. That weight should apply consistently to:

- positional fragment coverage
- fragment mass
- midpoints
- fragment ends
- fragment length counts
- GC-bias estimation and correction
- other fragment-derived features

Paired fragment spans keep the existing project semantics from `pos` to
`reference_end`. This plan does not introduce special handling of the gap
between paired alignments.

### Coverage identifies excess support, not its cause

For a position `x`, model the observed raw fragment coverage as:

```text
observed_coverage(x) = expected_coverage(x) + nuisance_coverage(x)
```

A single sample can show that coverage is extreme. It generally cannot prove
whether the cause is mismapping, residual duplication, reference structure,
or real narrow biology.

The command should therefore describe the result as a nuisance weight or
outlier keep weight. It should not claim to identify erroneous fragments or
specific artifact mechanisms.

### Detection and correction are different decisions

The detection threshold answers:

> Is this region extreme enough to correct?

The replacement target answers:

> What fraction of its fragment support should remain?

These values need not be equal. A strict upper-tail threshold may be useful
for calling a region, while the expected mean or another baseline may be the
appropriate correction target after the region has been called.

Supported target concepts should include:

- retain coverage equivalent to the expected local mean from the underlying
  ZIP fitted with the right-truncated second likelihood, as the default
- retain coverage equivalent to the detection threshold
- retain zero contribution from fragments associated with the region

The detection method and the target mode must be recorded separately in the
output metadata.

### Fragment weights are expected removal, not classification

If a region receives a keep weight of `0.4`, this should be interpreted as a
fractional correction of the fragment population associated with that region.
It does not mean the command has identified each fragment as 40 percent
artifactual.

The output should therefore store `keep_weight` directly. A field named
`prob_nuisance` would imply a level of probability calibration that the
initial method does not establish.

### Broad copy-number changes are not the primary target

The command targets sparse extreme support. Broad copy-number changes occur
on a much larger scale and are ordinarily far smaller than the extreme
pileups of interest.

Local calibration is intended to let the expected coverage shift with a broad
baseline change while still detecting sparse extremes within it. Full CNA
segmentation is not part of this command.

## Why positional clipping is not the primary application model

Multiplying positional coverage by a capped positional factor controls the
final value at those positions. It does not remove the influence of the
fragments that caused the excess.

In particular, positional clipping leaves those fragments represented in:

- nearby positions outside the called interval
- midpoint and end counts outside the interval
- fragment length distributions
- GC-bias estimation
- any other fragment-level statistic

It can also make the same fragment fully trusted at some positions and only
partly trusted at others, even when the assumed nuisance mechanism acts on the
fragment as a whole.

The primary downstream behavior should instead be:

1. Determine whether the current fragment overlaps any outlier region.
2. Resolve a single outlier keep weight for the fragment.
3. Multiply that weight with the fragment's other independent correction
   factors.
4. Use the resulting fragment weight for every contribution from that
   fragment.

This naturally reduces the surrounding coverage according to the actual
fragment starts, ends, and lengths. It avoids constructing an approximate
positional decay curve from an average fragment length.

## Proposed command boundary

`cfdna outliers` should:

- read a single BAM or CRAM
- use the ordinary fragment inclusion filters
- calculate raw positional fragment coverage before GC correction or genomic
  scaling
- build configurable local coverage histograms in a first BAM pass
- fit local zero-inflated Poisson models and upper-tail thresholds
- detect sparse extreme regions in a second BAM pass
- assign a regional fragment keep weight
- write a sparse interval TSV with implicit weight `1.0` outside the listed
  intervals
- write exact and flanked BED files for exclusion workflows
- write the observed histograms, fitted distributions, thresholds, and a
  combined diagnostic plot

The command should not rewrite the input alignment file by default. Downstream
commands already encounter each fragment and can resolve its weight from the
sparse interval file at that point.

## Initial detector

The initial detector should be a zero-inflated Poisson model fitted to raw
positional fragment coverage. This is the first model to implement and
evaluate, not a permanent claim that other upper-tail models can never be
useful.

The model must be evaluated from its output plots across samples with
different sequencing depths. If the observed positive-count distribution has
substantial overdispersion that ZIP does not capture, another count model can
be considered later.

Upper-tail decisions must use the survival probability:

```text
P(coverage >= observed_coverage)
```

They must not use the probability mass at only the observed coverage value.

The automatic probability threshold should be:

```text
1 / number_of_eligible_positions_in_the_selected_reference
```

Although the ZIP is fitted locally, the denominator remains global. The
default therefore means that fewer than a single eligible position this
extreme is expected in the analyzed sample under its local null model.

Users should be able to override the automatic probability directly or apply
a multiplier to it. Exact CLI names remain open. A multiplier greater than
`1.0` makes calling less stringent.

## Fixed-size core coverage histograms

The first BAM pass should calculate a raw positional coverage histogram for
each contiguous fixed-size model core and write those histograms to disk.

Initial settings should match the scale used by the scaling-weight commands:

- model core size of 500 kb
- configurable core size
- ordinary eligible-position and blacklist handling
- raw coverage before GC correction or genomic scaling

The fixed 500 kb size is a starting default, not part of the scientific
definition. `cfDNAlab` must work across species and cfDNA-like analyses with
different reference sizes and contig structures.

Each core histogram should record how many eligible positions have each
integer coverage. Its accompanying summary must contain enough additive
support information to decide whether a wider fitting context contains enough
fragments. The exact support measure remains open. Candidates include a
fragment count assigned once per fragment or eligible fragment coverage mass
converted to an effective fragment count.

Histograms from adjacent cores are additive. The genome-wide histogram is
therefore obtained by summing the core histograms after the first pass. No
separate global coverage calculation is required.

The global histogram should also receive the same diagnostic ZIP fits and
threshold calculation. The global result provides a sample-level fit plot and
reference point, while the local `T2` assigned to each core controls calling.

The histogram files are not a temporary positional coverage track. Their size
depends on the observed coverage values per core rather than on a row for
every position or constant-coverage run.

## Local fitting contexts

Each model core should receive its own local ZIP model and final discrete
coverage threshold.

Build its fitting context as follows:

1. Start with the core's coverage histogram.
2. Expand by a complete model core on the left and a complete model core on the right.
3. Continue paired expansion until the context covers at least the configured
   minimum span and contains the configured minimum fragment support.
4. Sum the selected core histograms and fit the core's ZIP model from the
   combined histogram.

With a 500 kb core and paired expansion, a minimum span of 5 Mb resolves to a
5.5 Mb context because the possible centered spans are 0.5, 1.5, 2.5, 3.5,
4.5, and 5.5 Mb.

At a contig boundary, expansion should continue on the available side after
the other side is exhausted. If the complete contig still lacks enough
eligible support, the core should use a clearly recorded fallback model built
from the global histogram. A fixed physical span must never be assumed to
guarantee a valid model because references can contain short contigs and
heavily blacklisted regions.

Adjacent cores will usually have strongly overlapping fitting contexts. The
histograms are combined after the first BAM pass, so this running-window
calculation does not recalculate positional coverage or reread the BAM.

Whether distance weighting within the fitting context improves the local ZIP
fit remains open. The first implementation should not silently copy the
triangular smoothing semantics of scaling weights unless weighted ZIP fitting
is explicitly defined and validated.

## Two-stage ZIP fitting

Each local context should use two ZIP fits. The second fit protects the normal
coverage model from values already identified as implausible by the first
fit.

The order is:

1. Fit an initial ZIP to the complete, unmodified local histogram.
2. Calculate the initial threshold `T1` from the selected upper-tail
   probability.
3. Exclude histogram observations with coverage greater than or equal to `T1`
   from the second fit.
4. Fit the parameters of an underlying ZIP to the remaining observations with
   a right-truncated likelihood conditioned on coverage being below `T1`.
5. Calculate the final calling threshold `T2` from the untruncated survival
   distribution of that fitted underlying ZIP.
6. Call positions from the original coverage using `T2`.

No mass is moved into `T1`, and the procedure does not iterate. The raw
histogram remains unchanged for diagnostics, global aggregation, and final
calling.

The initial design should use `T1` itself as the fit-exclusion boundary. A
future option could place the exclusion boundary at a still more extreme tail
probability, but that adds another parameter and should only be introduced if
the diagnostics show a concrete need.

Both thresholds have lasting meanings:

- `T1` is the fit-exclusion threshold
- `T2` is the final outlier threshold

## Two BAM passes

The command should read the BAM twice.

The first pass should:

- calculate the per-core coverage histograms and additive support summaries
- write the histogram data to disk
- sum the core histograms into the global histogram
- build adaptive local contexts
- perform the two-stage ZIP fits
- assign a final threshold to every model core

The second pass should:

- recalculate raw positional coverage
- apply the finalized threshold for the current model core
- combine adjacent qualifying positions into outlier regions
- calculate observed and target regional coverage mass
- calculate regional keep weights
- write exact and flanked BED intervals

Two BAM passes avoid a genome-wide temporary positional coverage track,
provisional candidate thresholds, species-specific pilot assumptions, and
candidate recall failures in regions with unusually low local coverage.

## Building outlier regions

The basic region caller should combine adjacent positions that meet the high
detection threshold into a contiguous outlier region.

For each region, calculate:

```text
observed_mass = sum of observed positional coverage in the region
target_mass = sum of the selected target coverage in the region
keep_weight = min(1, target_mass / observed_mass)
```

If the target mode is exclusion, `keep_weight` is `0.0`.

Using a regional mass ratio avoids making the single highest position dictate
the weight of every fragment touching the region. It also gives a clear
regional invariant before fragments overlap multiple outlier regions:

> Applying the same regional keep weight to every fragment contributing inside
> an isolated region reduces the region's total coverage mass to its target
> mass.

It does not guarantee an exact cap at every position. That is intentional.
The correction models a fraction of the contributing fragments as nuisance
rather than enforcing a positional output constraint.

## Optional hysteresis region growth

A high threshold may identify only the center of an extreme pileup even when
the same nuisance process produces high but sub-threshold coverage around it.
Full-fragment weighting already propagates correction beyond the called
region for fragments that touch the high-threshold core. It cannot affect
nuisance fragments that occur only in the surrounding elevated region and do
not touch the core.

Hysteresis region growth is an optional way to include those fragments:

1. Use a strict high threshold to seed a region.
2. Grow left and right while a lower extension criterion remains satisfied.
3. Calculate the regional keep weight over the resulting grown interval.

The extension criterion could be:

- a lower absolute coverage threshold
- a lower tail-score threshold from the same detection model
- coverage above a multiple of the estimated expected baseline

The correct criterion, gap tolerance, and safeguards against growing through
broad elevated regions are unresolved. Hysteresis should remain disabled or
experimental until those choices are evaluated on representative data.

This option is preferable to applying a generic smoothing curve because it
changes which fragments are associated with the outlier. It does not create a
second positional correction layer that downstream fragment averaging could
smooth again.

## Downstream fragment application

For a fragment `f`, let `overlapping_regions(f)` contain the outlier regions
intersecting its complete fragment span.

Resolve the fragment's outlier keep weight as:

```text
1.0                                      when there are no overlaps
minimum regional keep weight             otherwise
```

Equivalently, if nuisance scores are discussed internally, this selects the
maximum nuisance score.

### Why minimum is the default

Outlier regions overlapped by the same fragment are not independent evidence.
They were detected from the same coverage signal and may be different parts
of the same pileup.

- Multiplication would assume independent evidence and can overcorrect.
- Averaging would dilute the strongest evidence.
- The minimum keep weight lets the strongest associated region determine the
  fragment's reliability without counting correlated evidence repeatedly.

This is a deterministic combination rule, not a probabilistic independence
model.

The final fragment contribution should be:

```text
final_fragment_weight =
    outlier_keep_weight
    * gc_weight
    * genomic_scaling_weight
    * any other independent fragment correction
```

Within the outlier channel, combine overlaps by minimum. Between independent
correction channels, combine factors by multiplication.

### Required consistency across commands

All participating commands should use the same fragment interval lookup and
minimum-over-overlaps rule. A fragment should receive the same outlier keep
weight regardless of whether the command later counts its span, midpoint,
ends, fragment length, or another feature.

For positional coverage, add the fragment's resolved scalar weight across its
complete span. Do not first accumulate unweighted coverage and then multiply
only the called outlier positions.

GC-bias estimation should accept the outlier intervals so nuisance-weighted
fragments have reduced influence while the correction model is fitted. This
places outlier handling before GC correction and downstream genomic scaling in
the processing hierarchy.

## Output format

The primary correction output should be a sparse, bedGraph-like TSV:

```text
chrom  start  end  keep_weight
```

Required semantics:

- coordinates are 0-based and end-exclusive
- rows are sorted and non-overlapping within each chromosome
- omitted positions have implicit keep weight `1.0`
- weights are finite and in the closed interval from `0.0` to `1.0`
- adjacent intervals with identical weights may be combined
- the file is interpreted using fragment-overlap semantics, not ordinary
  positional scaling semantics

The command should also write:

- an exact BED containing the called outlier regions
- a BED containing the same regions with a configurable flank
- the per-core coverage histograms
- a model summary containing the fitting context and thresholds for each core
- global observed and fitted histogram data
- a combined diagnostic plot of the observed histograms, fitted ZIP
  distributions, and thresholds

The exact BED is the sample-specific blacklist without expansion. The flanked
BED supports workflows where removing the surrounding fragment-associated
region is preferable to fractional correction. Both should be retained so
the chosen flank does not replace the original calls.

Diagnostic tables should contain only information needed to reproduce or
assess the result. The core histogram data need core coordinates, coverage,
and observed eligible-position count. The model summary needs core and
context coordinates, eligible support, model source or fallback, fitted ZIP
parameters, `T1`, and `T2`. Zero inflation belongs in this model summary, not
in the outlier interval TSV.

Metadata should record at least:

- command and package version
- fragment inclusion filters
- core size, minimum context span, and minimum support
- detection method, upper-tail probability, and automatic-threshold
  denominator
- target mode and target parameters
- whether hysteresis growth was enabled
- blacklist flank
- interval combination rule
- coordinate convention

The outlier loader should be separate from the current genomic scaling loader.
The current scaling semantics average positional factors over spans and assume
full genomic coverage. This outlier format is sparse and assigns a regional
weight to every fragment that overlaps the interval.

Combining exact or flanked BED files across samples can be explored as a
separate cohort workflow. The single-sample command should not silently turn
sample-specific calls into a shared blacklist.

## Important edge cases

### A fragment barely overlaps an outlier

The initial rule gives the fragment the full regional weight even if the
overlap is short. This follows the fragment-level nuisance assumption, but it
creates sensitivity to region boundaries.

Do not silently introduce overlap-fraction dilution. That would weaken the
regional correction and return toward positional scaling semantics. Evaluate
boundary behavior explicitly before adding an alternative.

### A fragment overlaps multiple regions

Use the minimum regional keep weight. If this regularly causes unexpected
correction, inspect whether nearby regions should have been called as a single
connected outlier rather than changing the combination rule first.

### Very long called regions

Very long regions may represent a different process, including broad baseline
change, and can affect a large fragment population. The initial design should
report their lengths and make them easy to audit. A future implementation may
need a warning or an optional maximum region length, but should not silently
discard them without evidence.

### Regions crossing model-core boundaries

Thresholds may differ between adjacent cores. Qualifying positions should
still be joined across a core boundary when they are contiguous. Their
observed and target mass should be accumulated using the threshold or target
assigned to each position's core. A core boundary must not create an
artificial region boundary.

### Short or poorly supported contigs

A contig may be shorter than the minimum context span or contain too few
eligible positions after masking. The command should use the recorded
fallback model rather than fitting an unstable local ZIP or silently omitting
the contig.

### Zero keep weights

A zero regional keep weight excludes every fragment touching that region from
all participating downstream statistics. The CLI and metadata should state
this directly because it is stronger than masking only the interval itself.

## Validation plan

Validation should distinguish detection quality from application behavior.

### Detection validation

- Inspect observed histograms, fitted ZIP distributions, and both thresholds
  across samples with different depths.
- Compare the observed positive-count tail with the fitted tail. Zero
  inflation alone does not establish that the positive counts are Poisson.
- Verify the first fit, the right-truncated second fit, and the final call
  against the original histogram independently.
- Check how `T1` and `T2` change with sequencing depth and downsampling.
- Validate adaptive context expansion, minimum support, contig-edge handling,
  and global fallback on short and heavily masked contigs.
- Verify that the global histogram is exactly the sum of eligible core
  histograms.
- Verify that broad amplifications and deep deletions shift their local models
  instead of being interpreted through a genome-wide mean.
- Evaluate the 500 kb core and minimum 5 Mb context defaults across reference
  sizes and cfDNA-like assays.
- Evaluate whether hysteresis growth joins only clearly associated elevated
  shoulders.

### Application validation

- Construct fragments that generate a narrow extreme core and verify that
  full-fragment weighting reduces their complete spans.
- Verify that an isolated region reaches its requested total target mass.
- Show explicitly that exact per-position clipping is not promised.
- Verify implicit weight `1.0` for fragments with no outlier overlap.
- Verify minimum keep weight for fragments overlapping multiple regions.
- Verify zero-weight exclusion.
- Verify that the same fragment receives the same outlier weight in coverage,
  midpoint, end, fragment length, and GC-related paths.
- Compare the data-derived surrounding reduction with a positional cap and a
  fragment-length smoothing heuristic.

Tests should use valid fragment lengths and preserve the project's directional
`pos` to `reference_end` fragment semantics.

## Implementation stages

### Stage 1. Histogram and model diagnostics

- Calculate and persist configurable per-core coverage histograms and additive
  support summaries in a BAM pass.
- Build adaptive local contexts and the global histogram from core
  histograms.
- Implement the initial and right-truncated ZIP fits.
- Calculate `T1` and `T2` from the configured survival probability.
- Write model summaries, fitted distribution data, and the combined diagnostic
  plot before treating ZIP as an accepted default.

### Stage 2. Outlier calling and standalone outputs

- Recalculate positional coverage in the second BAM pass.
- Call and join regions using each core's finalized threshold.
- Support threshold, expected-mean, and zero target modes.
- Write sparse regional keep weights plus exact and flanked BED files.
- Verify output and metadata invariants, including fallback reporting.

### Stage 3. Shared downstream application

- Add a dedicated sparse outlier interval loader.
- Add shared fragment-overlap weight resolution.
- Use the minimum keep weight across overlaps.
- Integrate the resolved fragment weight into positional and fragment-level
  commands.

### Stage 4. Earliest correction layers

- Allow GC-bias estimation and correction workflows to consume outlier
  intervals.
- Confirm that outlier weights combine multiplicatively with independent GC
  and genomic scaling factors.

### Stage 5. Optional improvements

- Evaluate hysteresis region growth.
- Evaluate distance weighting within local fitting contexts.
- Evaluate overdispersed or empirical alternatives if ZIP diagnostics show a
  concrete lack of fit.
- Evaluate principled cohort-level combination of exact or flanked blacklist
  intervals.
- Consider more complex fragment-weight fitting only if regional weights show
  a concrete failure that simpler rules cannot address.

## Decisions retained by this plan

- The primary correction is a full-fragment keep weight.
- Detection uses raw positional coverage before GC correction and genomic
  scaling.
- The initial detector is a local two-stage ZIP using survival probabilities.
- Coverage histograms default to configurable 500 kb model cores.
- Every core receives an adaptive local model using at least the configured
  context span and fragment support, with a recorded fallback when necessary.
- Automatic tail probability is based on the global number of eligible
  positions.
- Coverage is calculated in two BAM passes, without a temporary genome-wide
  positional coverage track.
- The output is sparse with implicit keep weight `1.0`.
- The default replacement target is the underlying local ZIP mean estimated by the
  right-truncated second fit.
- Exact and flanked BED outputs support exclusion workflows.
- Histograms, fitted distributions, thresholds, and a combined plot are
  required diagnostics.
- Regional severity is based on observed and target coverage mass.
- Multiple outlier overlaps combine by minimum keep weight.
- Independent correction channels combine by multiplication.
- Fixed 10 bp bins are not part of the generalized design.
- Positional smoothing based on average fragment length is not the default.
- Hysteresis is an optional region-building improvement with unresolved
  defaults.
- The command remains single-sample and does not attempt causal artifact
  classification or CNA segmentation.

## Open questions

- What additive support statistic and minimum value should decide whether a
  local context contains enough fragments?
- Should the configured minimum local span default to 5 Mb, acknowledging
  that paired expansion from a 500 kb core produces a 5.5 Mb context?
- Should local context histograms be unweighted, or can distance weighting be
  defined without making ZIP fitting and threshold interpretation unclear?
- Should the first-fit threshold `T1` remain the second-fit exclusion boundary,
  or do diagnostics justify a distinct, more extreme boundary?
- What concise CLI options should expose the automatic probability multiplier
  and a manual probability?
- Should regional mass use every called position equally or account for masked
  and otherwise ineligible bases?
- How should nearby regions be merged before fragment application?
- Is a minimum seed-region length needed, or does fragment coverage already
  provide enough spatial persistence?
- Which lower criterion and stopping rule make hysteresis useful without
  absorbing broad elevated regions?
- Should hard exclusion require an explicit flag distinct from a target value
  of zero?
- What should the default blacklist flank be?
- Which commands must support the outlier file before the feature is considered
  usable?
