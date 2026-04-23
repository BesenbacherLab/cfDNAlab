## `cfdna outliers` system spec

Date: 2026-04-23

## Purpose

`cfdna outliers` detects extremely high local cfDNA support that is more
consistent with mapping or technical artifacts than with interpretable
fragmentomic biology, then writes a multiplicative scaling track that caps
their influence in downstream commands.

The command is upper-tail only. It reduces unusually high signal. It does not
try to correct low-support regions, broad copy-number structure, or general
local weirdness.

## Core Contract

The command should:

- read one sample
- count raw support in fixed genomic scoring bins
- identify bins with implausibly high upper-tail support
- convert those bins into scaling factors in `[0, 1]`
- write a full-coverage scaling TSV compatible with existing cfDNAlab scaling
  loaders
- write diagnostics explaining what was capped, what was not capped, and why

The command should not:

- replace blacklists or mappability resources
- call CNAs
- suppress broad regions just because they look unusual
- fit local residual models whose main effect is to remove biology
- build cohort-derived resources
- combine sample evidence with cohort priors in the core detector
- infer anything from already normalized or smoothed output unless a separate
  raw-count reconstruction path is explicitly designed

## Artifact Target

The target event is a narrow high-support pileup.

Typical causes:

- collapsed or repetitive sequence in the reference
- poor or ambiguous mappability
- missing decoy or alternative sequence
- alignment artifacts
- library-specific technical pileups
- residual nuisance loci missed by ordinary blacklist or mappability filters

Non-target events:

- chromosome-arm shifts
- broad copy-number changes
- broad tissue-of-origin or nucleosome-organization signal
- moderate plateaus spread across many bins
- low-support regions

Broad regions are treated as biology or upstream normalization problems unless
there is explicit evidence that they are technical artifacts. The command should
not silently flatten them.

## Input Signal

The detector operates on raw count-like support produced directly from the
sample, before GC correction, length normalization, genomic smoothing, or
outlier scaling.

The primary signal is coverage-style support:

- divide each selected chromosome into fixed scoring bins
- count raw fragment-base support per bin
- use integer support values for the tail model

Coverage-style support is the cleanest target because a per-bin cap maps
directly to positional coverage scaling. Fragment-averaged downstream commands
will dilute narrow low-weight intervals across the full fragment span, so the
same track does not have identical numerical meaning in every command family.

## Ordinary Filtering

Ordinary read and fragment filters are applied before scoring:

- duplicate, secondary, supplementary, and QC-failed reads are excluded
- paired-end consistency filters follow the same project conventions as other
  fragment commands
- mapping-quality filters are applied before counting
- fragment-length filters are applied before counting
- user-provided blacklist or mappability exclusions can be provided

Excluded genomic regions do not enter model fitting. By default, external
exclusions should be emitted with neutral weight `1.0` rather than converted to
outlier weights. This keeps `cfdna outliers` from becoming a hidden blacklist
application step.

Scoring bins that overlap excluded regions need explicit exposure handling.
The simplest correct behavior is:

- bins fully covered by excluded regions are unscored and emitted with neutral
  weight
- bins partially overlapping excluded regions are unscored unless an
  exposure-adjusted scoring mode is implemented
- unscored bins are reported in diagnostics

If users want hard exclusion, that should be explicit and separately labeled in
the output diagnostics.

## Scoring Bins

The scoring bin size is a first-class parameter.

The bin size must be:

- fine enough to localize narrow technical pileups
- coarse enough that the selected count model has useful support
- aligned with the intended downstream use case

The default scoring bin size is an open question and should be calibrated
against the artifact target.

Terminal bins shorter than the scoring bin size should not distort the fit.
Reasonable handling options:

- exclude terminal short bins from model fitting and emit neutral weight unless
  they are long enough to score reliably
- score terminal bins with an exposure-adjusted threshold

The implementation should choose one behavior and report it in metadata.

The detector cannot localize signal within a scoring bin. A bin-level weight
caps the total support assigned to that bin and applies uniformly across its
span in the scaling track.

## Tail Model

The core detector is chromosome-wise upper-tail clipping on raw counts.

For each chromosome:

1. collect kept scoring-bin counts
2. fit an upper-tail count model
3. find the smallest clip count `k_clip` whose upper-tail probability is below
   the configured threshold
4. mark bins with `observed_count > k_clip` as candidate caps

The default model should be a simple zero-aware count model:

- `mu = mean(counts)` over kept bins
- `p_nonzero = number_of_nonzero_bins / number_of_kept_bins`
- `P(X > k) = p_nonzero * (1 - PoissonCDF(k; mu))`

This is not a full maximum-likelihood zero-inflated Poisson. It is a pragmatic
upper-tail model that keeps the detector interpretable and stable.

## Fit Robustness

The detector should not let the artifacts it is trying to cap define the
background.

The model fit should therefore include one of these safeguards:

- fit once, identify extreme preliminary candidates, then refit after excluding
  those candidates
- fit the mean on a high-tail-trimmed count set
- fit on all bins but fail or warn when capped bins carry enough mass to
  materially change the chromosome mean

The exact safeguard is an implementation decision, but the invariant is not:
extreme candidate bins must not be allowed to silently inflate the cap threshold
enough to hide themselves.

## Tail Threshold

The threshold should be expressed in probability space, not as a fixed count.

The default threshold should control how many bins are expected to cross the
tail by chance under the fitted model. A natural rule is:

`tail_threshold = expected_false_caps / total_scored_bins`

where `total_scored_bins` is computed across all chromosomes being processed.

The default value for `expected_false_caps` is an open question. It should be
small because the command is meant to cap only extreme technical pileups.

The same probability threshold can be used for every chromosome while the
model parameters are fit per chromosome.

## Candidate Runs

Candidate bins are merged into contiguous candidate runs after thresholding.

Each candidate run should be classified as one of:

- `capped`: narrow enough to match the artifact target
- `broad_high_support`: too broad to cap automatically
- `excluded`: ignored because it overlaps user-provided excluded regions
- `unscored`: ignored because the chromosome or bin did not meet scoring
  requirements

Broad candidate runs should not be capped by default. They should be reported
as diagnostics because they may represent biology, copy-number structure,
sample contamination, or a failed model assumption.

The maximum run width eligible for automatic capping is an open question. It
should be specified in bases, not only in number of bins.

## Weight Mapping

For capped bins:

`scaling_factor = min(1, k_clip / observed_count)`

For bins that are not capped:

`scaling_factor = 1`

For broad high-support candidate runs:

`scaling_factor = 1` by default, with the run reported in diagnostics.

For external blacklist or mappability exclusions:

`scaling_factor = 1` by default, unless the user explicitly asks for hard
exclusion.

If `k_clip = 0`, the formula yields zero for positive observed counts. The
implementation should decide whether this is acceptable for the selected mode
or whether a minimum soft floor is needed. That decision must be reported in
metadata.

## Output Scaling Track

The main output is a full-coverage TSV with the columns required by the current
scaling-factor loader:

- `chromosome`
- `start`
- `end`
- `scaling_factor`

Coordinates are 0-based, half-open `[start, end)`.

The output must satisfy the existing scaling-factor contract:

- one set of rows per processed chromosome
- sorted by chromosome and start coordinate
- no gaps
- no overlaps
- every processed chromosome starts at `0`
- every processed chromosome ends at the chromosome length
- scaling factors are finite and non-negative

Adjacent bins with identical scaling factors should be merged to keep the file
compact. Most rows should have `scaling_factor = 1`.

## Metadata

The scaling TSV should include leading comment metadata lines.

Required metadata:

- command name and version
- source sample identifier when available
- counting model
- scoring bin size
- tail model
- tail threshold rule
- total scored bins
- expected false caps setting
- mapq filter
- fragment-length filter
- chromosome list
- external exclusions used
- hard-exclusion behavior
- broad-run behavior

Metadata is required because the scaling factor alone is not enough to know
what scientific assumption produced the track.

## Diagnostics

The command should write diagnostics in addition to the scaling TSV.

Per-chromosome model diagnostics:

- chromosome
- number of scoring bins
- number of kept bins
- number of excluded bins
- number of zero bins
- `mu`
- `p_nonzero`
- `k_clip`
- tail threshold
- number of candidate bins
- number of capped bins
- number of broad high-support bins
- capped base pairs
- capped support mass before capping
- estimated support mass removed by capping

Run-level diagnostics:

- chromosome
- start
- end
- class
- number of bins
- max observed count
- mean observed count
- `k_clip`
- minimum scaling factor
- mean scaling factor
- minimum tail probability
- overlap with external exclusions, if any

Summary diagnostics:

- total capped base pairs
- total capped runs
- fraction of scored genome capped
- fraction of support mass removed
- chromosomes with no usable fit
- warnings and fail-fast reasons

## Failure And Warning Conditions

The command should fail or warn loudly rather than silently producing a
misleading track.

Fail:

- input BAM or reference metadata cannot be read
- requested chromosomes are missing
- scoring bin size is invalid
- no scoring bins remain after filtering
- output scaling track cannot satisfy the full-coverage contract
- raw support cannot be represented as non-negative integer counts for the
  selected detector

Warn:

- a chromosome has too few kept bins to fit reliably
- a chromosome is all zero or nearly all zero
- a large fraction of a chromosome is classified as candidate high support
- capped bins remove a large fraction of total support mass
- broad high-support runs are found but not capped
- sex chromosomes are included without an explicit chromosome choice
- external exclusions remove a large fraction of scoring bins

Warnings should be visible in logs and summarized in the diagnostics output.

## Command Shape

Likely CLI shape:

```text
cfdna outliers \
  --bam sample.bam \
  --output-dir outliers \
  --output-prefix sample \
  --score-bin-size <bp> \
  --mapq <int> \
  --min-fragment-length <bp> \
  --max-fragment-length <bp> \
  --tail-expected-false-caps <float> \
  --max-capped-run-bp <bp> \
  --chromosomes <list> \
  --exclude-bed <path>
```

The exact names can follow project CLI conventions during implementation.

The important public concepts are:

- raw support is counted from the sample
- scoring happens in fixed bins
- the detector is upper-tail only
- broad candidate runs are reported, not silently capped
- output is a scaling-factor track

## Interaction With Existing Scaling Factors

Outlier weights answer a different question than broad smoothing weights.

- broad smoothing weights ask: what is the expected large-scale baseline here
- outlier weights ask: is this local observation too extreme to trust

These corrections are orthogonal enough to compose by multiplication.

The command should output an outlier track only. Users can combine it with
other scaling factors by multiplication. A helper command for multiplying
scaling tracks may be useful, but it is separate from outlier detection.

## Hard Masking

Soft capping is the default behavior.

Hard masking can be useful for inspection or for workflows that truly want
exclusion, but it changes the product from "cap extreme influence" to "remove
regions". If hard masking exists, it should be explicit and reflected in:

- output metadata
- diagnostics
- run-level class labels

Hard masking should not be triggered implicitly by external blacklist input.

## Open Questions

Default scoring bin size:

- Needs calibration against real artifact examples and typical cfDNA depth.
- Candidate values should be judged by whether they catch narrow technical
  pileups without turning broad biology into caps.

Tail false-cap budget:

- The default `expected_false_caps` should be small.
- Needs calibration against clean samples and samples with known mapping
  artifacts.

Fit robustness method:

- Choose between iterative refit, high-tail trimming, or fit-all plus strict
  diagnostics.
- The chosen method must keep true extreme bins from inflating their own cap.

Maximum capped run width:

- Needs calibration in base pairs.
- The default should prevent broad biological regions from being capped.

Terminal-bin handling:

- Choose neutral terminal bins or exposure-adjusted scoring.

Minimum soft floor:

- Decide whether `scaling_factor = 0` is allowed when `k_clip = 0`, or whether
  a minimum positive floor is required for soft mode.

Fragment-count mode:

- Decide whether a separate fragment-count support mode is needed.
- If added, it must be documented as a different support definition, not as the
  same coverage clip applied differently.

Sex chromosomes:

- Decide whether autosomes should be the default or whether selected
  chromosomes should exactly follow user input.

Hard mask option:

- Decide whether it belongs in the command or should remain a downstream use of
  diagnostics.

Scaling-track multiplication helper:

- Decide whether users need a built-in helper to multiply smoothing and
  outlier tracks.

## Acceptance Criteria

The system is behaving correctly when:

- known narrow high-coverage technical pileups receive weights below `1`
- ordinary regions receive weight `1`
- broad high-support regions are reported but not automatically flattened
- external blacklist input does not silently become an outlier track
- the output scaling TSV loads through the existing scaling-factor loader
- diagnostics make every capped region auditable
- changing the scoring bin size changes resolution, not the scientific target
- the command fails or warns when the model is being used outside its valid
  regime
