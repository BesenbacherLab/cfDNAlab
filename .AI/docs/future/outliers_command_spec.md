## `cfdna outliers` spec

Date: 2026-04-21

## Scope

This spec defines a research-backed design direction for a possible
`cfdna outliers` command.

The command is intended to produce a genome-wide multiplicative weighting
track that leaves most of the genome at `1.0` and selectively downweights a
small number of clearly aberrant loci or short regions.

This is a concept and product spec only. It does not propose implementation
work yet.

## Status and assumptions

This revision is based on:

- the current `cfDNAlab` scaling-factor and `fcoverage` architecture
- the LIONHEART Nature Communications paper
- literature on blacklists, mappability, CNV normalization, and robust
  outlier handling

Important limitation:

- I could not inspect the local `/Users/au547627/Downloads/cfdna_outlier_design.pdf`
  from the shell because the sandbox could not access `Downloads`
- the public LIONHEART paper was accessible online

So this spec reflects:

- your written notes
- the public LIONHEART methods
- the current codebase
- external literature

It does not yet reflect the contents of the local design PDF unless they
overlap with your notes.

## Why this command exists

The current scaling commands solve a different problem:

- `cfdna coverage-weights`
- `cfdna fragment-count-weights`

These commands normalize broad baseline shifts by generating multiplicative
weights across the whole genome. They are designed for slow, large-scale
variation such as CNA-like background changes.

That is not the same as sparse local extremes.

For many downstream fragmentomics analyses, a few highly aberrant loci can
dominate statistics even when the rest of the genome behaves well. In those
cases, a targeted soft-capping track is conceptually cleaner than pretending
the problem is just another broad smoothing problem.

The realistic scope for this package is narrower than "make a better
blacklist from one BAM".

The more defensible scope is:

- assume standard blacklist-style resources may already be available
- then detect clear residual nuisance signal those resources missed
- and softly downweight those residual high-signal loci

## Core distinction

The word "outlier" hides multiple targets that should not be collapsed into
one method:

1. Assembly or mappability artifact
2. Recurrently high-signal locus across many samples
3. Sample-specific extreme local spike
4. Broad copy-number shift
5. True narrow biology that is inconvenient for one downstream analysis

These are not the same thing.

If one command tries to optimize for all of them simultaneously, it will
either become incoherent or quietly bias the science.

The strongest conclusion from the literature review is therefore:

- a sample-specific soft-capping command and
- a cohort-derived reusable blacklist resource

should be treated as related but distinct products.

For this package, the first product should also be framed as residual
correction rather than blacklist replacement.

## Fit with the current codebase

The current scaling-factor loader already expects:

- per chromosome
- sorted bins
- no gaps
- no overlaps
- non-negative multiplicative factors
- full chromosome coverage

That means an outlier track already fits the existing application path if it
is written as a full-coverage scaling TSV.

Crucially, the bins do not need to be fixed-width.

So an efficient output is possible:

- long runs with `scaling_factor = 1.0`
- short intervals with `0.0 < scaling_factor < 1.0`

This means the file can be bedGraph-like in spirit while still satisfying the
current scaling-factor contract.

## Relationship to existing scaling factors

Conceptually:

- downstream smoothing weights handle broad baseline variation
- outlier weights handle sparse local extremes

If both are present, the correct numerical relationship is simple:

- `combined_weight = smoothing_weight * outlier_weight`

There are two architectural ways to do that:

- pre-multiply the tracks into one scaling file
- let downstream commands apply two scaling tracks during counting

Recommendation:

- v1 should not add two-track runtime logic
- prefer pre-multiplication into one final scaling file

Reason:

- the mathematics is trivial
- the runtime code stays simpler
- the semantics stay easier to audit

If track combination becomes common, a separate helper command could be added
later.

Important distinction:

- `cfdna outliers` may internally estimate a broad baseline or smooth trend in
  order to define relative extremeness
- but that internal baseline is part of detection methodology, not part of the
  downstream command contract

So:

- downstream genomic smoothing and outlier correction should be independent
- the emitted outlier track should work whether or not the user also applies a
  separate smoothing track downstream

## The main design axes

Before choosing a method, it helps to separate the design axes that are easy
to conflate.

### Axis 1. Where the problem is handled

Possible layers:

- reference or aligner level
- static exclusion resource level
- cohort-derived recurrent exclusion level
- sample-specific post-alignment weighting level

Examples:

- sponge or decoy sequences and T2T-like assemblies attack the problem before
  read pileups are turned into features
- ENCODE blacklist and Umap attack it with reusable exclusion resources
- LIONHEART-style recurrent empirical bin removal attacks it at the cohort
  level
- a future `cfdna outliers` command would attack it at the per-sample
  weighting level

Implication:

- `cfdna outliers` should be framed as one layer in a stack, not as the whole
  artifact strategy

### Axis 2. Which baseline is removed first

Possible spaces for scoring:

- raw support
- GC-corrected support
- support after broad baseline normalization
- residuals after local trend removal

This matters because the same local spike can look:

- mild on the raw scale
- severe after broad normalization
- or invisible after an overaggressive local detrending step

Implication:

- the scoring space is a first-class design decision
- broad baseline normalization should be explicit, not implicit
- and that baseline can be internal to the outlier detector rather than a
  required external preprocessing step

### Axis 3. How anomaly evidence is computed

Main families:

- annotation-derived artifact priors
- recurrent control or cohort evidence
- global count-model tail probabilities
- local robust residual scores
- segmentation or change-point scores
- multiscale denoising or decomposition scores
- higher-dimensional anomaly scores

Implication:

- different scores answer different questions
- combining them blindly is not principled unless each score has a clear role
- in this package, annotation-aware information is more plausible as an
  optional secondary prior or annotation layer than as the main command
  identity

### Axis 4. How scores become weights

Possible policies:

- hard exclusion
- fixed soft cap
- monotone score-to-weight mapping
- rank-based or quantile-based capping
- winsorization of support followed by reconstruction of weights

Implication:

- two methods with the same detector can behave very differently downstream if
  the capping policy differs

### Axis 5. What unit is actually being penalized

This is easy to miss.

Possible semantics:

- penalize genomic positions
- penalize fragments by averaging over the fragment span
- penalize counted overlap only
- reject a fragment entirely if it touches a bad region

These are not minor implementation details.

They define what "outlier correction" means to downstream commands.

Implication:

- a genomic scaling track is only one way to express outlier handling
- it is attractive because the repo already supports it
- but it is not semantically equivalent to fragment rejection or per-window
  gating

## Research summary

## What existing methods actually do

### 1. Static annotation-driven exclusion

This is the classical blacklist and mappability route.

Key references:

- ENCODE blacklist
- Umap / Bismap
- QDNAseq-style exclusion of problematic bins

What they do well:

- identify known problematic reference regions
- remove many repeat-driven or low-mappability artifacts before any modeling
- provide stable, reusable resources across studies

What they do badly:

- they are not sample specific
- they do not detect new sample-specific extremes
- they can become reference-build and read-length specific
- they risk silently becoming outdated or population-biased

Evidence:

- ENCODE blacklist was built from many input datasets using read depth and
  mappability in overlapping 1 kb windows, deliberately using many samples to
  avoid blacklisting cell-type-specific CNVs
- Umap quantifies uniquely mappable regions for a chosen read length, directly
  addressing short-read ambiguity
- QDNAseq developed an additional empirical blacklist from 1000 Genomes data,
  with most flagged regions beyond the earlier ENCODE set

Interpretation:

- this is best viewed as baseline exclusion, not as the core logic of a new
  sample-specific `cfdna outliers` command
- the command should assume these resources may already have removed many of
  the obvious mappability and repeat-driven artifacts
- its main value is more likely to be catching residual nuisance spikes rather
  than rediscovering blacklist logic from one shallow sample

### 1b. Alignment-guided alternatives to blacklists

This family does not try to score outliers after alignment. It tries to reduce
 the artifact load before it appears.

Examples:

- sponge or decoy sequences
- more complete assemblies such as T2T-style references

What they do well:

- reduce alignment-driven pileups without needing a large post hoc exclusion
  set
- may be scientifically cleaner than blacklisting when the true issue is
  missing or collapsed sequence in the reference

What they do badly for this project:

- require changing upstream mapping assumptions
- are not a drop-in extension of the current scaling-factor infrastructure
- do not solve sample-specific nuisance spikes by themselves

Interpretation:

- these are important alternatives when thinking about future hg38 resources
- they are upstream of `cfdna outliers`, not part of the command contract

### 2. Cohort-derived empirical recurrent nuisance-bin resources

This is an important separate family.

Instead of using only annotation, these methods flag bins that repeatedly
behave abnormally across many samples or datasets.

Examples:

- ENCODE blacklist itself is effectively this kind of method
- LIONHEART built a recurrent empirical bin resource by thresholding 10 bp bins
  within each sample, then keeping bins that crossed the threshold in a chosen
  fraction of samples within datasets and unioning those bins across datasets
- WisecondorX and related CNA tools rely heavily on blacklists, and showed
  that different tools can mask materially different fractions of the genome

What they do well:

- find recurring problem loci missed by static annotation
- support building reusable hg38 resources
- reduce dependence on any single study's hand-picked blacklist

What they do badly:

- require many samples
- depend on assay, aligner, read length, reference, and preprocessing details
- can accidentally encode cohort-specific biology if the cohort mix is wrong

Interpretation:

- this is valuable, but it is a resource-generation problem
- it should not be the only or default behavior of a sample-specific command
- it is distinct from annotation-driven mappability logic
- it is also distinct from single-sample clipping, even when the recurrent
  resource is built by aggregating the outputs of a single-sample detector

### 2b. Dynamic greylists

This is a useful conceptual middle ground between fixed blacklists and
per-sample soft weights.

The best-known analogy is GreyListChIP:

- tile the genome
- count support in an input or control sample
- fit a negative-binomial background
- derive a threshold
- merge nearby flagged intervals into a sample-specific or experiment-specific
  exclusion set

What they do well:

- adapt to a specific sample or batch
- remain simple to explain
- produce explicit intervals rather than opaque scores
- pair naturally with existing static blacklists

What they do badly for this project:

- they are intrinsically hard-thresholded
- they fit better to "exclude these regions" than to "softly reduce their
  influence"
- they still need a choice of tiling, distributional model, and merge logic

Interpretation:

- a greylist is not the same product as an outlier scaling track
- but it is a strong conceptual reference for a future hard-mask mode or
  cohort-resource workflow

### 3. Sample-specific global count-model tails

This is the family closest to the LIONHEART outlier removal step.

Typical models:

- Poisson
- zero-inflated Poisson
- negative binomial
- overdispersed Poisson
- zero-inflated or hurdle extensions

What they do well:

- define an explicit tail probability
- work naturally on sparse count-like bins
- make it easy to explain why a bin was considered extreme

What they do badly:

- global thresholds do not distinguish local spikes from broad CNAs
- misspecification is easy when overdispersion is strong
- chromosome-wise fits can still be distorted by large aberrant segments

Interpretation:

- useful as one score
- insufficient as the only score for a general outlier command

### 3b. Empirical-rank and quantile methods

This family deliberately avoids a strong parametric distributional claim.

Examples:

- top-x-percent capping
- chromosome-wise extreme-rank thresholds
- empirical tail cutoffs after baseline normalization

What they do well:

- are robust to model misspecification
- are simple to audit
- can be surprisingly stable when the only real question is "how extreme is
  this compared with the rest of the chromosome?"

What they do badly:

- the thresholds are partly arbitrary
- they do not naturally extrapolate out-of-sample
- they can behave badly when chromosome-wide distributions differ strongly
  between runs

Interpretation:

- these are reasonable fallback or benchmarking methods
- they are weak as the sole scientific story for the main command

### 4. Sample-specific local robust residual methods

This is the strongest conceptual match for a new soft-capping command.

Representative ideas:

- rolling median baseline
- local quantile baseline
- local MAD or Hampel-style residual scoring
- Winsorization of residuals instead of binary masking

Evidence from related genomics literature:

- the `copynumber` package explicitly recommends outlier handling before
  segmentation and describes median-filter plus MAD-based Winsorization of the
  residuals
- LMADCNV builds local features and MAD-based anomaly scores to exploit the
  positional dependence of read-depth data

What they do well:

- ask the right question for local spikes:
  "is this bin extreme relative to its neighborhood?"
- are robust to a modest number of extremes
- naturally support soft downweighting rather than hard exclusion

What they do badly:

- can suppress real narrow biology if the command is too aggressive
- are sensitive to neighborhood width and chromosome-edge handling
- still need some protection against broad CNA baselines

Interpretation:

- this should be the center of a sample-specific `cfdna outliers` design

### 4b. Panel-of-normals and negative-control normalization

This family uses controls to estimate background bias rather than trying to
infer everything from the case sample alone.

Examples:

- panel-of-normals workflows in production CNV callers
- CODEX2 negative-control samples or negative-control regions

What they do well:

- subtract workflow-specific systematic bias
- help separate technical background from true biology
- are much better than pure single-sample methods when a well-matched control
  resource exists

What they do badly for a generic cfDNA outlier command:

- require carefully matched library prep, aligner, panel or genome binning,
  and often sex balance
- controls can themselves contain CNVs or workflow artifacts
- the result is less portable than it first appears

Interpretation:

- this is compelling for future resource building
- it is not a strong default dependency for a single-sample cfDNA CLI unless
  the project commits to maintaining a matched control resource

### 5. Segmentation and state-space methods

Examples:

- CBS
- HMMs
- piecewise-constant segmentation
- negative-binomial HMMs

What they do well:

- detect broad copy-number states
- model persistent shifts well
- produce interpretable segments

What they do badly for this use case:

- they are designed for state changes, not sparse local capping
- they are relatively heavy for a command whose desired output is
  mostly `1.0`
- they solve the wrong primary problem if the target is local anomalies

Interpretation:

- useful as context
- useful for broad baseline normalization
- not the core engine for `cfdna outliers`

### 5b. Multiscale denoising and decomposition

Examples:

- wavelet denoising
- total variation denoising
- cumulative-sum style persistent-shift detection

What they do well:

- separate broad structure from local noise
- preserve breakpoints better than naive moving averages in some settings
- suggest useful ways to split the signal into broad and local components

What they do badly for this use case:

- they are usually designed to recover cleaner copy-number profiles, not to
  produce interpretable soft-capping weights
- the resulting "denoised" signal still needs another rule to decide which
  loci are downweighted
- they can become hard to explain in a command that should be transparent

Interpretation:

- this family is valuable mainly as inspiration for decomposition
- broad component versus local residual is the useful lesson here, more than
  the exact denoiser

### 6. Feature-aware or multivariate fragmentomic scores

Until now, most of the spec has implicitly discussed one-dimensional support.
That is too narrow if the long-term goal is a genuinely fragmentomics-aware
 outlier command.

Examples:

- CRAG integrated fragmentation score, which combines local fragment count and
  average fragment size
- joint use of coverage, fragment size, end patterns, GC residual, or
  mappability context

What they do well:

- may better isolate loci that are strange in fragmentomic terms rather than
  only high in raw support
- can separate "boring high coverage" from "high coverage plus abnormal size
  profile"
- are closer to the eventual downstream scientific features in some workflows

What they do badly:

- are harder to apply consistently across commands that do not all use the same
  fragmentomic dimensions
- reduce portability of the output scaling track
- can create a mismatch between the signal used to build the weights and the
  signal modified downstream

Interpretation:

- this family is scientifically important
- but if used, it should probably be explicit that the resulting weights are
  feature-aware and not universally meaningful

### 7. Density-based or general anomaly-learning methods

Examples:

- Local Outlier Factor
- covariance-based anomaly scoring
- broader unsupervised anomaly methods

What they do well:

- capture complex local structure
- can combine multiple features

What they do badly:

- are harder to interpret and tune
- add substantial methodological complexity
- are weak defaults for a scientific CLI that should be auditable

Interpretation:

- not suitable as the default v1 method
- maybe later for offline research, not for the main command contract

## What the cfDNA literature implies specifically

### LIONHEART

The LIONHEART paper is directly relevant because it already uses both kinds of
normalization that matter here:

- recurrent extreme-bin handling
- broad CNA-like baseline normalization

The reported method conceptually mixes separate layers:

- excludes ENCODE blacklist v2 and `umap k100` non-mappable bins
- counts coverage in 10 bp bins
- fits chromosome-wise zero-inflated Poisson distributions per sample
- derives per-sample upper-tail thresholds and truncates sample-specific
  extreme coverages
- separately removes recurrent empirical extreme bins across datasets
- then divides coverage by a broad 5 Mb / 500 kb overlapping-average baseline

That is not a pure local method.
It is already a hybrid of:

- annotation/resource-based exclusion
- empirical recurrent nuisance-bin removal
- single-sample upper-tail clipping
- broad baseline normalization

Important implication:

- broad baseline normalization and local extreme handling were separated
- that separation should be preserved in a new command design

### Other cfDNA and shallow-WGS workflows

Related cfDNA and sWGS methods repeatedly rely on:

- blacklist filtering
- mappability correction or exclusion
- GC correction
- broad-bin normalization for CNAs

This is consistent with the idea that a good outlier command should not try to
replace those layers. It should sit after or on top of them.

### Reference bias warning

Recent work on reference-based cfDNA fragmentomics showed that fragmentomic
biases can vary with the reference genome and ancestry background.

Implication:

- a future reusable hg38 "outlier resource" should be treated as conditional
  on assay, reference, and preprocessing
- the command should not oversell any cohort-derived blacklist as universally
  correct

### What the broader CNV literature adds

The non-cfDNA CNV literature reinforces several points:

- a large amount of apparent local extremeness is really background bias or
  broad state change viewed at the wrong scale
- panel-of-normals and negative-control regions can improve technical
  background estimation when matched controls exist
- segmentation-oriented denoisers such as total variation are mainly useful
  because they preserve persistent structure and edges, not because they
  directly define good soft-capping weights

The main transferable lesson is:

- decompose first
- score second
- cap third

not:

- pick a complicated detector and hope it solves everything at once

### What dynamic greylists add conceptually

Greylist-style thinking adds one useful design lesson:

- there may be value in supporting both
  - a continuous soft-weight output
  - and a discrete interval output

These are not interchangeable, but they solve related problems.

The soft-weight output is better for:

- gradual influence reduction
- reuse through current scaling infrastructure
- avoiding brittle yes/no cutoff behavior

The interval output is better for:

- explicit exclusion workflows
- inspection in genome browsers
- future resource generation

This does not mean v1 needs both.
It means the spec should not assume the scaling-track output is the only
possible representation.

### What multivariate fragmentomic scores add conceptually

CRAG is useful here not because it solves the same problem, but because it
shows a different scoring philosophy:

- coverage alone is not the only meaningful local signal
- joint signals such as count plus fragment size can produce more biologically
  targeted loci
- local background plus global background can be combined in one detector

The catch is important:

- a weight track derived from a combined fragmentation statistic is less
  neutral than a coverage-only track

That may be good or bad depending on downstream intent.

The same caution applies to annotation-aware priors:

- they can be useful for ranking confidence or annotating outputs
- but if they dominate the score, the command starts to become a
  blacklist-generator surrogate rather than a residual outlier detector

## Downstream command impact

This section matters because the same scaling track has different practical
effects depending on how each command consumes it.

### Family A. Per-position signal multiplication

Commands in this family multiply a per-base array directly:

- `fcoverage`
- `wps`
- positional parts of `fragment-kmers`

Implications:

- narrow low-weight bins have immediate positional effect
- hard zeroing creates literal holes in the positional signal
- even very short downweighted regions can change downstream segment structure
  or correlations

This family is the most sensitive to aggressive narrow masking.

### Family B. Fragment-averaged scaling over the full fragment

Commands in this family average scaling over the whole fragment span before
counting:

- `midpoints`
- many paths in `ends`
- many paths in `lengths`
- `bam_to_bam`
- `bam_to_frag`

Implications:

- a very narrow downweighted region gets diluted across the full fragment
- single-bin spikes may have little effect unless the bins are fine enough or
  the low-weight interval is long relative to fragment length
- broad low-weight regions are much more influential than point-like ones

This means:

- a track optimized for per-position analyses may be too narrow to matter much
  for fragment-level feature extraction

### Family C. Window-overlap-weighted scaling

Some `ends` and `lengths` modes use overlap-based scaling rather than
whole-fragment averaging.

Implications:

- the same outlier interval can affect only the overlapping part of the count
  geometry
- window assignment policy now matters
- the difference between "scale over fragment" and "scale over counted
  overlap" can become meaningful near narrow outlier regions

### Family D. Separate coverage and fragment-count weight channels

`bam_to_bam` and `bam_to_frag` already distinguish:

- coverage-like scaling
- count-like scaling

Implications:

- the repo already has a semantic precedent for two weight meanings
- an outlier command that ignores this distinction may be scientifically muddy
  for downstream tagged-fragment workflows

### Family E. Commands that could have been handled by fragment rejection

Some downstream tasks could conceptually respond better to:

- "drop fragments touching problematic loci"

than to:

- "reduce their average weight"

Implications:

- the current scaling-track infrastructure biases the design toward continuous
  weighting
- that is good for reuse and composability
- but it also means the command should acknowledge that it is choosing a soft
  weighting semantics rather than a fragment-filtering semantics

### Consequence for command design

The above behavior argues against a single naive objective such as:

- "find tiny local spikes and downweight them"

because that objective will:

- strongly affect per-position commands
- weakly affect fragment-averaged commands

A more honest framing is:

- choose the counting model according to the downstream use case
- choose the scoring resolution with fragment-averaging dilution in mind
- keep hard masking limited because its effects are uneven across command
  families

It also argues for being explicit about whether the command is:

- coverage-centric
- fragment-centric
- or feature-aware

It also argues for separating:

- internal baseline estimation used to define relative outliers
- downstream genomic smoothing used as an independent normalization layer

## Product conclusion

The most coherent product split is:

### Product A: sample-specific soft outlier capping

This is the main candidate for `cfdna outliers`.

Core behavior:

- take one sample
- compute one or more local anomaly scores
- turn extreme loci into multiplicative downweights
- write a full-coverage scaling TSV

This directly serves downstream feature extraction.

### Product B: cohort-derived recurrent exclusion resource

This should be a separate workflow, command, or documented pipeline.

Core behavior:

- analyze many controls or many datasets
- identify recurrently problematic loci
- export a reusable blacklist or default exclusion resource

This could be built quite simply from the outputs of a good single-sample
downweighting command:

- run the command per sample
- collect bins or intervals that are repeatedly downweighted above some chosen
  severity threshold
- aggregate recurrence within and across datasets

This is scientifically valuable, but it is not the same user contract.

Recommendation:

- do not make `cfdna outliers` try to be both products in v1
- do not promise that it can recover most of what mature blacklist resources
  already catch from one shallow cfDNA BAM

## Scope for this package

The command should be scoped as:

- residual nuisance detection after ordinary filtering

not as:

- blacklist replacement
- one-sample mappability discovery
- a universal problematic-region finder

In practical terms that means:

- existing blacklist and mappability resources remain the first line of
  defense
- `cfdna outliers` operates on what is left
- success should be judged by whether it catches obvious residual spikes, not
  by whether it reconstructs ENCODE blacklist or Umap-like resources
- if its per-sample downweighting turns out to be useful across many samples,
  those outputs can later be aggregated into a recurrent nuisance-bin resource
- the command should emit one standalone weight track rather than requiring
  users to also run a specific smoothing workflow for coherence

## Annotation-aware priors

Annotation-aware information is still useful, but only in a demoted role.

Reasonable roles:

- annotate detected loci as being in known suspicious contexts
- slightly lower the evidence threshold in already suspicious regions
- help prioritize output review

Unreasonable role for v1:

- make annotation the main signal used to "discover" problematic regions from
  one shallow BAM

That broader promise is not realistic under the intended use case.

## Ranked ideas for this package

This ranking is not about what is most sophisticated in general.
It is about what is most useful for this package under the actual constraints:

- one BAM file
- short-read cfDNA
- usually low coverage
- existing scaling-factor infrastructure
- downstream commands that already consume genome weights in different ways

### Tier 1. Most useful for this package

#### 1. Broad-baseline normalization plus local robust residual scoring

Summary:

- remove broad baseline first
- then score local residual extremeness with a robust method such as rolling
  median plus MAD
- convert only upper-tail extremes into soft downweights

Why it ranks first:

- works with a single BAM
- does not require a maintained control cohort
- is robust to low-coverage noise if the local bins are not too fine
- maps naturally to the current scaling-track infrastructure
- remains interpretable enough for CLI users

Main strengths:

- best balance of practicality and scientific honesty
- handles local spikes without pretending to solve all artifact sources
- naturally supports soft capping instead of brittle exclusion
- aligns well with the narrower "residual nuisance after standard filtering"
  scope
- can use an internal broad baseline for relative scoring while staying
  independent of downstream smoothing choices

Main weaknesses:

- still depends strongly on bin size and local window width
- can suppress real narrow biology if used too aggressively
- needs broad baseline normalization to avoid confusing CNAs with local spikes

Package fit:

- excellent

#### 2. Broad-baseline normalization plus one simple global tail score

Summary:

- after broad normalization, compute a chromosome-wise extreme-tail score
- use it either alone or as a companion to the local score

Why it ranks second:

- still single-sample compatible
- still easy to explain
- gives a complementary notion of absolute extremeness

Main strengths:

- explicit statistical interpretation
- useful cross-check against local-only scoring
- especially relevant in shallow data where random local noise can still
  produce odd neighborhoods

Main weaknesses:

- model choice is fragile under low coverage and zero inflation
- global tails alone are not enough
- more sensitive than empirical methods to misspecification

Package fit:

- strong as part of a hybrid
- weaker as a standalone default

#### 3. Coverage-vs-fragment-count modes for the same command

Summary:

- one outlier command
- at least two counting semantics:
  - coverage-centric
  - fragment-count-centric

Why it ranks third:

- the repo already has this semantic distinction elsewhere
- downstream commands do not all respond to weights in the same way

Main strengths:

- keeps the design honest about downstream meaning
- lets the same command serve both positional and fragment-level workflows

Main weaknesses:

- increases conceptual burden in docs
- users may choose the wrong mode without clear guidance

Package fit:

- high

This is not a detector by itself, but it is one of the most important design
choices for package usefulness.

### Tier 2. Useful, but secondary or optional

#### 4. Empirical-rank or quantile tail methods

Summary:

- use rank or quantile thresholds rather than a strong parametric model

Why this ranks here:

- very compatible with single-sample shallow data
- but scientifically less satisfying as the main story

Main strengths:

- robust benchmark
- good fallback when parametric fitting is unstable
- easy to audit

Main weaknesses:

- threshold choice is arbitrary
- poor out-of-sample story
- less elegant than local robust residuals

Package fit:

- good as a benchmark mode or sanity-check mode

#### 4b. Annotation-aware prior as a secondary layer

Summary:

- use mappability, repeat burden, or related context as an annotation prior or
  annotation
- do not let it define the main command contract

Why this ranks here:

- useful in moderation
- unrealistic as the primary engine for one shallow cfDNA BAM

Main strengths:

- improves interpretability
- may help separate "residual nuisance in suspicious context" from unexplained
  spikes

Main weaknesses:

- easy to overstate what one sample can infer
- drifts toward blacklist reconstruction if made too central

Package fit:

- moderate as a secondary layer
- poor as the main method

#### 5. Greylist-style hard-threshold interval output

Summary:

- identify extreme intervals and emit a discrete exclusion set

Why this ranks here:

- conceptually useful
- may be handy for inspection and later resource-building
- but less aligned with the existing soft scaling pathway

Main strengths:

- easy to inspect
- easy to compare across samples
- good complement to static blacklists

Main weaknesses:

- hard-thresholded
- less reusable across downstream commands than soft weights
- more brittle in shallow data

Package fit:

- moderate as an optional separate output or future mode
- not the main v1 output

#### 6. Feature-aware fragmentomic scoring

Summary:

- score loci using a combined statistic such as count plus fragment size

Why this ranks here:

- scientifically interesting
- may eventually be better for some fragmentomics tasks
- but much less universal across downstream commands

Main strengths:

- can capture loci that are weird in fragmentomic terms rather than merely
  enriched in support
- may align better with some downstream biological questions

Main weaknesses:

- harder to explain
- lower portability of the resulting weight track
- easy to mismatch the score definition with the downstream usage

Package fit:

- plausible future research branch
- not the best default package contract

### Tier 3. Interesting, but low priority for this package

#### 7. Dynamic greylists as the main product

Why it ranks low:

- the package already has a natural soft-weight application path
- a hard interval product is less aligned with that path

Still useful for:

- separate interval outputs
- future recurrent-resource workflows

#### 8. Panel-of-normals or negative-control normalization as a dependency

Why it ranks low:

- not compatible with the default one-BAM use case
- maintaining a matched control resource is a different project

Still useful for:

- later cohort-resource building
- offline research

#### 9. Segmentation-first or CNV-caller-like approaches

Why it ranks low:

- they solve the wrong primary problem
- they are broader and heavier than needed for sparse local capping

Still useful for:

- broad baseline estimation
- contextual comparison

#### 10. LOF, general anomaly learning, or deep models

Why it ranks low:

- too hard to justify and tune under shallow single-sample cfDNA
- poor fit for a transparent CLI

Still useful for:

- offline method exploration
- publications, if they truly outperform simpler methods

### Ranked recommendation

If the package needs one strongest direction, it is:

1. broad baseline normalization
2. local robust residual score
3. optional global tail companion score
4. soft downweighting
5. explicit coverage-vs-fragment-count mode
6. optional annotation-aware labeling or mild prior, but not annotation-led
   discovery

That is the best package fit I can currently see.

## Ranked ideas beyond this package

These are interesting, but they should not dominate the package design unless
the project scope changes.

### Most interesting beyond-package directions

#### 1. Cohort-derived recurrent outlier resources

This is the strongest next research branch if you eventually want:

- better hg38 nuisance resources
- robust defaults for many datasets

Why interesting:

- it could produce genuinely reusable knowledge
- it could improve more than this one command

Why not the main package focus now:

- it depends on many samples and workflow consistency
- but the barrier is lower if the single-sample command already emits
  interpretable downweight tracks that can be compared across samples

#### 2. Upstream alignment-guided strategies

Examples:

- sponge sequences
- decoy-aware references
- more complete assemblies

Why interesting:

- if the problem is fundamentally reference collapse, this attacks the root
  cause rather than the symptom

Why not the main package focus now:

- it lives upstream of the current package abstractions

#### 3. Feature-aware fragmentomic outlier maps

Why interesting:

- this may eventually be more biologically meaningful than coverage-only
  capping

Why not the main package focus now:

- it is much harder to make one weight track that is honest across all current
  downstream commands

#### 4. Joint sample-plus-resource models

Examples:

- single sample plus panel-of-normals
- single sample plus negative-control regions
- single sample plus recurrent outlier prior

Why interesting:

- this may be the long-term best scientific solution

Why not the main package focus now:

- it requires new infrastructure and curated resources

## Recommended v1 philosophy

The v1 command should be:

- sample specific
- soft by default
- auditable
- composable with existing smoothing factors
- explicit about which score produced each downweight
- explicit that it targets residual nuisance after ordinary blacklist-style
  filtering
- independent of whether downstream genomic smoothing is also used

The command should not:

- pretend to build a universal hg38 blacklist from one sample
- pretend to replace existing blacklist and mappability resources
- use a single global threshold as the whole method
- silently hard-mask bins unless explicitly requested

It should also avoid pretending that one default is equally optimal for:

- positional coverage correlations
- midpoint profile extraction
- end-motif counting
- fragment length distributions

because the scaling semantics differ materially across those command families

And it should be explicit about asymmetry:

- high-support outliers are the primary target
- low-support anomalies are scientifically real, but they do not map naturally
  to a "downweight to cap them" contract

So a general anomaly detector and a capping command are not identical things.

## Recommended v1 method

The strongest design is a hybrid, but a narrow hybrid.

### Step 1. Start from already filtered data

The command should expect the usual baseline filtering layers to remain
available:

- mapping-quality filtering
- duplicate filtering
- existing blacklist support
- optional GC correction

For a future cohort workflow, fixed mappability exclusions should also be used.

This is not an implementation footnote.
It is part of the command identity.

The intended order is:

1. ordinary exclusion and correction layers
2. choose the signal space for scoring
3. estimate any internal baseline needed for relative outlier detection
4. residual outlier scoring
5. soft capping

This should not be read as:

- users must run downstream genomic smoothing before or after `cfdna outliers`

Instead:

- the outlier command should be coherent on its own
- downstream smoothing remains optional and independent

### Step 2. Count support in smaller bins than smoothing weights

The command needs a local scoring resolution.

This should be much smaller than the broad smoothing windows used in
`coverage-weights`.

The exact default is open, but the important point is conceptual:

- scoring bins should represent local structure
- output bins can later be merged aggressively wherever the final scaling
  factor is unchanged

Additional consequence:

- if the scoring bins are too fine, fragment-averaged downstream commands may
  barely feel the weights
- if the scoring bins are too wide, per-position commands may lose the local
  specificity that motivated the command in the first place

So the scoring-bin choice is not just a detector tuning parameter. It is a
cross-command semantics choice.

If feature-aware scoring is ever added, this step becomes:

- choose the local summarization unit
- choose the feature vector per unit

not just:

- choose the bin width

### Step 3. Normalize broad baseline before local scoring

This is critical.

Without this step, the command will confuse:

- broad CNAs
- broad library-depth shifts
- arm-level biases

with the local spikes it is supposed to capture.

Best current direction:

- divide local support by a much broader internal baseline before computing
  local anomaly scores

Possible sources of that baseline:

- an internal broad smoother inside `cfdna outliers`
- a robust broad trend estimate computed directly in the command
- in special cases, an external precomputed track if the user explicitly wants
  that

My recommendation:

- conceptually reuse the idea of broad-trend estimation
- do not make the local score operate on raw coverage by default
- do not make the command depend on an external smoothing workflow for default
  coherence

On your note about larger scales than `50 kb` to avoid CNA switching:

- that reasoning is directionally correct
- the broad baseline should be clearly larger than the local anomaly scale
- if the local scoring resolution becomes small, the baseline should remain
  much broader

Inference, not a direct literature rule:

- a broad baseline somewhere in the `250 kb` to multi-megabase range with a
  smaller stride is plausible, but the right default depends on the scoring
  bin size and on whether the goal is cfDNA open-chromatin work or a more
  generic artifact detector

What matters most is not whether this baseline matches the user's downstream
genomic smoothing exactly.
What matters is that it is broad and robust enough to define relative
extremeness without forcing actual segmentation into scope.

### Step 3b. Alternative internal baseline candidate: coarse-bin normalization
with breakpoint-bin handling

There is a more specific candidate that may be better than a single broad
smoother in samples with visible CNA-like shifts.

Concept:

1. partition the genome into coarse bins such as `50 kb`
2. calculate one mean per coarse bin
3. normalize the finer local signal within each coarse bin by that bin's own
   mean
4. perform outlier detection in the normalized space using neighboring coarse
   bins as the context

The main appeal is that broad local level shifts are mostly absorbed by the
coarse-bin normalization itself. The neighboring bins then become more
comparable without requiring full segmentation.

Potential advantage:

- a `100 kb` shift can largely be absorbed into two neighboring `50 kb` means
  rather than contaminating an entire wider pooled local window

Main weakness:

- the `50 kb` bin containing the actual shift boundary is mixed and therefore
  mis-normalized
- this creates a breakpoint-bin blind spot where outlier calling may be
  unreliable

That suggests an extension:

- compare adjacent coarse-bin means
- when the left-right difference exceeds a threshold consistent with a likely
  real shift rather than ordinary shallow-noise fluctuation, mark the boundary
  bin as breakpoint-like
- then either skip normal outlier calling there or use a more conservative
  rule

This is essentially segmentation-lite:

- more CNA-aware than a naive pooled local baseline
- much simpler than full segmentation

It should only activate when the mean shift is large enough to matter.
In many samples there will not be a strong enough signal for breakpoint
handling to be worth invoking.

The main conceptual tradeoff versus a global or chromosome-wise model is:

- coarse-bin normalization can reveal local-relative outliers that are not
  globally extreme
- that may be desirable for residual nuisance detection
- but it also changes the meaning of "outlier" away from pure absolute
  extremeness

So if this route is chosen, the docs should state explicitly that the command
is detecting local-relative excess after coarse-bin level-setting, not just
global upper-tail bins.

### Step 4. Compute two scores, not ten

The method should stay scientifically legible.

The best v1 score pair is:

- a global tail score on the broad-normalized support
- a local robust residual score on the same broad-normalized support

Why this pair:

- the global tail score catches bins that are extreme in absolute terms
- the local score catches bins that are extreme relative to nearby bins
- both are interpretable

Good v1 candidates:

- global: Poisson, ZIP, or negative-binomial-like tail score
- local: rolling median residual divided by local MAD

### Step 4b. When Poisson or ZIP is scientifically defensible

This needs to be stated explicitly because not all score spaces are equally
compatible with discrete count models.

#### True count-like space

Poisson or ZIP is most defensible when the scored quantity is still close to a
non-negative integer count, for example:

- raw count-like support in fixed bins
- raw coverage counts in sufficiently fine bins
- raw fragment-count-like support before strong smoothing or ratio
  transformations

In this space, Poisson- or ZIP-like tail scoring is conceptually cleanest.

#### Pseudo-count space after rescaling and discretization

Some transformed signals may still be usable for a Poisson-like tail fit if
they can be put back onto a count-like scale before fitting.

Examples:

- GC-corrected support, if the correction is roughly scale-preserving
- `--normalize-by-length` support after rescaling by a representative fragment
  length, ideally sample-specific when feasible

In these cases, discretization can be a practical heuristic:

- rescale to a count-like mean scale
- round or otherwise discretize
- then fit the global tail model

This is no longer a literal Poisson generative model.
It is a pseudo-count heuristic for thresholding.

That can still be useful, but the docs should say so plainly.

#### Smoothed or ratio-like residual space

This is the weakest case for Poisson or ZIP.

Examples:

- broad-smoothed support
- broad-baseline-divided support
- local-relative residuals near a baseline of `1`
- strongly averaged or smoothed signals

Even if such signals are rescaled back to the genome-wide or tile-wide mean and
then discretized, smoothing has already changed the object:

- neighboring bins become correlated
- local variance is reduced
- the tail becomes less count-like
- the signal becomes more ratio-like or residual-like than event-count-like

So in this space, Poisson or ZIP should be viewed as much more heuristic.
Robust continuous-space scoring is usually the more honest default.

#### Practical implication for this command

This suggests a three-way distinction:

- true count-like branch
- pseudo-count branch after rescaling and discretization
- continuous residual branch

Current best interpretation:

- raw count-like support -> ZIP or Poisson is fully plausible
- GC-corrected support -> ZIP or Poisson may still be plausible after
  discretization
- `normalize-by-length` support -> plausible after rescaling by representative
  fragment length and discretization
- smoothed or broad-baseline-divided support -> much weaker justification for
  Poisson or ZIP, robust continuous residual scoring is preferable

Rounding or discretization policy should be fixed and documented because it can
slightly affect extreme-tail thresholds.

Recommendation:

- if bins are very sparse and many zeros are expected, ZIP is defensible in
  the count-like or pseudo-count branches
- if overdispersion dominates, a negative-binomial-like model is more honest
- if the scored space is smoothed or baseline-divided, robust continuous
  residual scoring is the safer centerpiece

Weakness to keep in mind:

- a global ZIP or Poisson-like tail model may miss bins that are only extreme
  relative to their local neighborhood
- a local-relative method may flag bins that are not globally extreme at all

Neither notion is universally correct.
That is why the meaning of the score space and baseline must be explicit in
the output metadata and docs.

Benchmark alternatives that should stay in view even if they are not defaults:

- an empirical-rank tail score
- a greylist-style hard-threshold score
- a feature-aware combined score such as count plus fragment-size deviation

Optional non-core support:

- annotation-aware ranking or confidence annotation for already detected spikes

### Step 5. Convert scores to soft downweights

This step must stay separate from detection.

The command should compute raw scores first and then apply a capping policy.

Examples of sensible policies:

- hard mask above a strict threshold
- soft cap with a monotone penalty
- piecewise penalty with a mild default and a harsher optional mode

Recommendation:

- default to soft downweighting
- reserve hard zeroing for an explicit flag

Additional recommendation:

- allow the capping curve to be more conservative for very narrow intervals
  than for broader intervals

Reason:

- narrow intervals are diluted in fragment-averaged commands but harsh in
  per-position commands
- interval width is therefore part of the downstream effect, not just the
  detector output

Another important asymmetry:

- the command contract currently suggests weights mostly in `[0, 1]`
- that makes sense for capping enriched or nuisance-high loci
- it does not naturally express "this locus is abnormally depleted"

Recommendation:

- v1 should state clearly that it is primarily an upper-tail suppression
  command, not a symmetric anomaly normalizer

For relative detectors, the implied semantics are:

- the weight does not cap to one global absolute depth
- it caps to a threshold relative to the command's internal baseline

That is acceptable and likely preferable, as long as the docs state it
clearly.

### Step 6. Merge neighboring bins

After final scaling factors are produced:

- merge adjacent bins with the same or numerically equivalent final weight

This keeps the output small while preserving full genome coverage.

## Recommended output contract

The output should be a full-coverage TSV with at least:

- `chromosome`
- `start`
- `end`
- `raw_support`
- `broad_baseline`
- `broad_normalized_support`
- `global_outlier_score`
- `local_outlier_score`
- `baseline_mode`
- `outlier_class`
- `scaling_factor`

Possible future output files:

- sparse interval BED for hard-thresholded greylist-like regions
- summary table of the most extreme loci before and after capping

Possible extra metadata:

- `# counting_model=coverage|fragment_count`
- `# gc_mode=...`
- `# broad_baseline=...`
- `# score_space=...`
- `# fit_space=true_count|pseudo_count|continuous_residual`
- `# baseline_mode=broad_smoother|coarse_bin_norm|...`
- `# method=...`

The metadata should make it clear that:

- `broad_baseline` describes the internal detection baseline
- it does not imply that downstream smoothing must use the same baseline or be
  used at all
- `baseline_mode` affects the meaning of "outlier", especially when a
  local-relative or coarse-bin-normalized method is used
- `fit_space` describes whether any global discrete tail score was fit on a
  true count-like signal, a rescaled/discretized pseudo-count signal, or not

Recommendation:

- keep the final application column named `scaling_factor` so the file can fit
  the current infrastructure

## Counting model question

This command probably does need two scientific interpretations:

- coverage-based support
- fragment-count-based support

That mirrors the current smoothing split and matters for downstream meaning.

Recommendation:

- do not create two public commands immediately
- prefer one `cfdna outliers` command with an explicit counting mode

Reason:

- the user goal is one conceptual product
- the counting choice is a scientific mode, not a separate workflow

This mode should likely be tied to downstream intent in the docs:

- use `coverage` when the goal is to modify positional or coverage-like
  signals
- use `fragment-count` when the goal is to modify fragment-level support more
  directly

A future third family could exist:

- feature-aware

but only if the project is comfortable with weights that are no longer
universally meaningful across all downstream commands

## Candidate CLI direction

This is intentionally provisional.

Possible shape:

- `cfdna outliers`
- `--counting-model coverage|fragment-count`
- `--score-bin-size <bp>`
- `--broad-bin-size <bp>`
- `--broad-stride <bp>`
- `--score-space raw|gc-corrected|...`
- `--fit-space auto|true-count|pseudo-count|continuous`
- `--baseline-mode broad-smoother|coarse-bin-norm`
- `--method local-mad|global-tail|hybrid`
- `--cap-policy soft|hard`
- `--max-weight-reduction <float>`
- `--min-soft-interval-width <bp>`
- `--min-run-length <bp>`
- `--merge-distance <bp>`

Recommendation:

- keep the public surface small in v1
- the main default should be `hybrid`
- do not expose too many research knobs initially

If the CLI needs one extra knob beyond the basics, the most defensible one is
probably interval-width sensitivity rather than a long menu of detectors.

If a second extra knob is needed, the next most defensible one is probably the
scoring space, not another detector variant.

## What should be deferred

These ideas are real, but should probably not be in v1:

- density-based anomaly models such as LOF
- large score ensembles
- automatic universal hg38 blacklist generation from arbitrary sample sets
- runtime application of multiple scaling files simultaneously
- deep-learning-style anomaly detectors

Those add complexity faster than they add trustworthy science.

## Use cases

The command makes sense for:

- positional coverage analyses where a few loci dominate correlations
- non-positional analyses that already support genome scaling factors
- future reproducible downweighting of recurrent cfDNA nuisance loci
- exploratory work toward improved hg38 fragmentomics exclusion resources

The command is especially well matched to:

- feature extraction where extreme local support can distort correlation-like
  summaries
- workflows where soft downweighting is preferable to fully removing
  fragments or positions

The command is less naturally matched to:

- workflows that fundamentally want hard fragment rejection
- workflows where the nuisance signal is only visible in a multivariate
  fragmentomic space rather than in coverage-like support alone
- workflows whose main goal is rebuilding blacklist-like resources from one
  sample

The command does not make sense as:

- a replacement for broad CNA normalization
- a replacement for blacklist handling in general
- a CNV caller
- a claim that the flagged loci are always technical artifacts

## Risks and anti-goals

The biggest scientific risk is suppressing real biology.

That risk increases when:

- the local window is too small
- the capping is too aggressive
- broad baseline normalization is skipped
- coarse-bin breakpoint bins are treated as ordinary bins
- the method is framed as "artifact removal" rather than "downweighting
  potentially problematic extremes"

Another real risk is semantic drift:

- a local-relative method can create "outliers" that are not globally unusual
- a global tail method can ignore residual nuisance that is obvious in local
  context

The command should not hide that tradeoff.

The strongest anti-goal is therefore:

- do not present the output as truth
- present it as a controllable weighting layer

## Best current recommendation

If this command is built, the most coherent v1 is:

- one sample-specific command
- one explicit counting mode
- broad baseline normalization first
- one global tail score plus one local MAD-style score
- soft downweighting by default
- full-coverage scaling TSV output
- easy composition with existing smoothing factors by multiplication

And separately:

- a future workflow for recurrent cohort-derived blacklists or hg38 outlier
  resources

That split matches both the literature and the current architecture better
than trying to force one command to be both.

## Remaining open questions

- What should the default scoring bin size be for general use?
- Should the broad baseline be computed internally or supplied externally?
- Should the default internal baseline be a broad smoother or coarse-bin
  normalization with breakpoint-bin handling?
- Should the default global model be ZIP or something more overdispersion-aware?
- At what adjacent-bin difference should breakpoint handling activate in noisy
  low-coverage samples?
- Should interval width explicitly influence the score-to-weight mapping?
- Should positional and fragment-averaged downstream use cases get different
  defaults even if they share one command?
- Should v1 expose only a continuous scaling output, or also support a
  greylist-like interval file?
- Is coverage-only scoring the right scope for v1, or should feature-aware
  scoring already be reserved in the design?
- How much of the output should be hard masked versus softly capped?
- Should the first version support only one score mode internally even if the
  spec recommends a hybrid endpoint?
- Do we want a separate helper for pre-multiplying scaling tracks?
- Do we want a future `make-blacklist` or `recurrent-outliers` command rather
  than overloading `cfdna outliers`?

## References

- LIONHEART, Nature Communications 2025:
  https://www.nature.com/articles/s41467-025-66503-3
- ENCODE blacklist, Scientific Reports 2019:
  https://www.nature.com/articles/s41598-019-45839-z
- Umap and Bismap, Nucleic Acids Research 2018:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6237805/
- QDNAseq blacklist and correction paper, Genome Research 2014:
  https://genome.cshlp.org/content/early/2014/11/03/gr.175141.114
- WisecondorX, Nucleic Acids Research 2019:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6393301/
- Beyond blacklists and sponge-sequence alternative, Bioinformatics 2026:
  https://academic.oup.com/bioinformatics/article/42/3/btag110/8519621
- CODEX2 negative-control normalization, Genome Biology 2018:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6260772/
- GreyListChIP package manual, dynamic greylist thresholding:
  https://bioconductor.org/packages/release/bioc/manuals/GreyListChIP/man/GreyListChIP.pdf
- `copynumber` outlier Winsorization and segmentation, BMC Genomics 2012:
  https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-591
- Total variation denoising for CNV signals, BMC Bioinformatics 2018:
  https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2332-x
- LOF-based CNV detection, IEEE/ACM TCBB 2021:
  https://pubmed.ncbi.nlm.nih.gov/31880558/
- LMADCNV local-feature MAD approach, IEEE/ACM TCBB 2025:
  https://pubmed.ncbi.nlm.nih.gov/41082433/
- CRAG integrated fragmentation score and local/global NB backgrounds,
  Genome Medicine 2022:
  https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01141-8
- Reference-bias warning for cfDNA fragmentomics, Cell Reports Methods 2024:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11228372/
