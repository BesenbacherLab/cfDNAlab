# Outliers: Claude review findings and open decisions

Current implementation update (2026-10-15): the first BAM-inspection choices and remaining work
are recorded in [section 8](#8-what-still-needs-to-be-settled). Earlier discussion and simulation
settings below remain historical. The production defaults are now ZINB and 5%, with selectable ZIP.

Date: 2026-10-12. Updated: 2026-10-15.

Status: Temporary discussion document. This records the review and subsequent discussion.
The original review discussion remains historical. The first BAM-inspection version is now
implemented, with current choices and remaining decisions recorded in section 8.

Related documents: [rectification plan](outliers_review_rectification_plan.md),
[original review](../major_reviews/outliers.md), and
[later real-data investigations](outliers_real_data_validation.md).

Sources are Claude's first response supplied as an attachment, its follow-up pasted into the
conversation, and the user's replies to both. The findings below preserve the substance of
those responses and distinguish it from Codex's assessment. Relevant code was re-read on
branch `outliers2`, HEAD `f5a75211`. Claude's simulation code, seeds, and complete results were
not supplied. The user subsequently authorized a fragment simulation study. Its methods and
results are recorded in [the simulation report](outliers_fragment_simulations/results_report.md)
and [reproduction notes](outliers_fragment_simulations/README.md). The first study completed
399 independent fragment replicates plus reanalysis of three primary replicates. The subsequent
[fitted ZIP/ZINB comparison](outliers_fragment_simulations/zinb_comparison_report.md) completed
684 paired fragment replicates. Those studies left production unchanged. The subsequent first
BAM-inspection implementation is recorded separately above. No project tests were run by the agent.

## Constraints carried forward

- The command supplements existing blacklists by handling residual extreme coverage.
- The required operating range remains 0.1x through 200x. Current data are about 0.5x.
  This range is a requirement to establish, not a result already demonstrated.
- Use the active-component mean as the agreed local-mean replacement target. This is `lambda`
  under ZIP and the negative-binomial component mean under ZINB. Keep reporting the component
  mean and the unconditional mean rather than changing their interpretation with depth.
- ZINB is now the default, with ZIP selectable. Both follow model-consistent initial fitting,
  T1 exclusion, right-truncated refitting, and T2 calculation. Numerical failure does not silently
  substitute ZIP.
- The current protection proposal starts a 100 kb genomic window at every exact call start.
  Every overlapping call receives evidence from all evaluated windows, including earlier
  calls' windows. The initial implemented cutoff is 5%. The earlier 10% setting remains a
  historical simulation setting, not the production default.
- Windows must retain physical locality. Do not collect 100,000 eligible positions by
  walking across large blacklisted blocks. Minimum support and the denominator floor are
  initially 50 kb, with the supported backward fallback described in section 8.
- Calculate extent from exact called bases before flanking. Flank all candidates and merge
  overlapping and touching intervals before filtering. A preservation decision applies to
  the entire merged region and all its member calls.
- Keep the original exact intervals and individual weights for selected calls. Do not pool
  weights, subtract preserved bases from flanks, or add fragment-level protection overrides.
- Outlier weights must work with `fcoverage --gc-file`, GC tags, and no GC correction.
- Zero outlier weight means skip before GC correction and counting. Keep normalization
  formulas and existing normalization ownership unchanged. Weighted restore-mean statistics
  are a user-owned later consideration only.
- Do not compare blacklists or identify/compare BAM files. Do not add BAM paths, hashes, or
  sample-name comparisons. The agreed reference check uses contig names and lengths, with a
  separate list of analysed chromosomes. It does not identify a BAM or sample.
- Keep input checks focused. Accept custom tables with a clear warning when the relevant
  reference or analysis information is missing. No broad checking framework is requested.
- CLI help is the primary documentation. Existing feature-branch outputs do not justify
  extra behavior solely for compatibility.

## 1. GC correction: estimate and apply from the same weighted signal

**Claude's finding:** It questioned whether multiplication restores the GC-corrected
background when the detector's target comes from raw coverage.

**Code evidence:** [fcoverage's three counting branches](../../../src/commands/fcoverage/fcoverage.rs#L1062)
already apply `base_weight * gc_weight * outlier_keep_weight`, or the corresponding expression
without GC.

**Working decision:** The useful workflow recommendation is to apply outlier scaling during
`gc-bias` estimation and use the same outlier scaling during fcoverage counting. Both then
operate on the same weighted fragment signal. Producer support remains part of the later
rollout, rather than an extra input requirement for the current fcoverage consumer.

Multiplication is already supported and does not require raw positional coverage inside
fcoverage. Keep the combined application in the simulations. Do not turn this review point
into a new consumer-side reconstruction, matching-file check, or normalization change.

## 2. Window support near masks and chromosome ends

**Claude's finding:** Dividing by a very small eligible span can preserve a narrow technical
call simply because its denominator is small. Its examples were a call just before a large
blacklisted block and the boundary fixture's final 30 bp call, with only 40 bp remaining.

**Codex assessment:** The small-denominator problem is real under that proposed formula.
It is not an argument that all rightward windows miss calls to their left. Earlier evaluated
windows already contribute to every overlapping call's maximum.

**Rejected proposal:** Claude first suggested gathering the first 100,000 eligible positions,
skipping masks and extending left at chromosome ends. The user rejected this because it can
join unrelated events across megabases of blacklisting. Keep genomic spans fixed.

**Claude's follow-up proposal:** Ignore windows below a minimum eligible count. Only for a
call with no usable overlapping window, evaluate a backward fallback ending at the call's
end. Use 50% of the configured window as the minimum. Neither the fallback nor 50% is agreed.
This fallback is unnecessary for calls that already have usable earlier windows. It also
introduces a second window-placement rule that needs justification.

**User's latest proposal:** An undersized window should not supply its inflated percentage.
However, it can still support preservation if it contains enough called bases to exceed the
cutoff at the smallest acceptable window size. An isolated tiny call should not be preserved
automatically merely because little chromosome remains to its right.

### A concrete interpretation to discuss

Define `window_bp`, `minimum_eligible_bases`, and `cutoff_fraction`. For each fixed genomic
window, calculate the union of exact called bases and the eligible unblacklisted base count.

1. When `eligible_bases >= minimum_eligible_bases`, use the ordinary fraction
   `called_bases / eligible_bases`.
2. When `eligible_bases < minimum_eligible_bases`, do not use that observed percentage.
   Instead, it can provide sufficient called-base evidence for preservation when
   `called_bases > cutoff_fraction * minimum_eligible_bases`.
3. Each call receives the supported percentages and any sufficient called-base evidence
   from every evaluated window it overlaps. Do not limit it to its own window.
4. Preserve a merged flanked region if any member meets a preservation rule. A low-support
   window alone would no longer automatically protect every member under this candidate.

For decision purposes, the two numeric comparisons can be expressed as
`called_bases / max(eligible_bases, minimum_eligible_bases) > cutoff_fraction`.
That is a denominator floor, not the actual observed percentage in an undersized window.
The audit should keep the actual counts and support status clear.

This captures the user's 15 kb of calls in 20 kb example without allowing the 30 bp/40 bp
boundary example to become a broad-event flag solely through its denominator. The sufficient
count is judged at the minimum acceptable scale, not necessarily at the full 100 kb scale.
Using `cutoff_fraction * window_bp` instead would be a stricter alternative that guarantees
the observed called count alone exceeds the cutoff even with a full-window denominator.
These are different policies and must not be conflated.

**Working direction:** Keep candidates correctable when this rule supplies insufficient
evidence for preservation. A small eligible denominator does not itself imply preservation.
This avoids separate rules for hundreds versus thousands of called bases. The minimum eligible
count remains unresolved. Compare candidate values in simulations without presenting a value
as an agreed CLI default. The backward fallback is not part of this working direction.

For a homogeneous amplification, the copy-number increase is the same throughout the event.
With the same ordinary rate and threshold, the chance of a T2 call is also the same throughout.
Sampling fluctuations and differences in locally fitted thresholds affect observed calls.
Absence of upstream calls is useful negative evidence when an event continuation should
generate many calls. Do not turn the lack of a logical proof into automatic preservation
of every small boundary call. Define the minimum extent the safeguard intends to protect
instead of relying on an assumed rarity of amplifications near blacklisted regions.

## 3. Called-base density does not bound lost fragment signal

**Claude's finding:** Its simulation used random fragment starts, 167 bp fragments, a 1 Mb
gain, a fixed background `lambda`, automatic T2, the lambda target, and the minimum-weight
rule. It reported that every scenario stayed below 10% called bases and was corrected.
A 2.35-fold gain at background coverage 17 reportedly fell to approximately 1.7-fold.

**Codex assessment:** The mechanism is valid. Short calls can affect many fragments and
remove coverage outside the exact called bases because a positive overlap applies the full
regional weight. Called-base fraction measures spatial extent, not affected fragment mass.

The initial review identified these gaps in treating that number as a prediction of the full
proposed behavior. Our subsequent reproduction is summarized in section 9 and the linked report:

- A mean called fraction across a 1 Mb event does not establish that every evaluated 100 kb
  call-start window falls below 10%. The maximum across all overlapping windows matters.
- Preservation propagates through merged flanked regions. A protected member can withhold
  correction from other members. It is unclear whether the simulation included that step.
- A fixed background fit is a stress case, not the actual two-stage fit. For a uniform
  1 Mb elevation of 2.35 times 17 within 5.5 Mb, the complete-context mean is about 21.2.
  Observations below T1 remain in the second fit. The fitted baseline need not return to 17.
- Context placement, chromosome-edge contexts, support fallback, and the run-wide eligible
  denominator all affect the actual thresholds and decisions.

**Claude's proposed remedies:** Add an intensity requirement, measure affected fragment mass,
or lower the density cutoff. It also proposed recording mean and maximum coverage per call.

**Working direction:** Positional coverage supplies evidence for preservation. Autocorrelation
affects uncertainty, not the validity of measuring a called-base fraction. Do not add an
affected-fragment-mass screen. Measure fragment signal loss in simulations to assess whether
the coverage-based safeguard actually preserves biological elevations. Mean and maximum
coverage remain useful audit candidates, not an agreed new classifier.

For a homogeneous amplified window and a fixed threshold, the expected called fraction is
`P(coverage >= T2 under the amplified rate)`. This depends on depth, amplification strength,
dispersion, and the threshold. Autocorrelation does not change that marginal expectation.
For an event covering only part of a window, combine the amplified and ordinary probabilities
according to their eligible spans. Calculate these relationships when selecting the density
cutoff. A gain generating no calls needs no protection from correction. The concerning cases
are sparse calls that cause substantial signal loss through whole-fragment weighting.

Claude also observed that a contiguous call exceeding 10% of a full window is preserved by
construction. The user accepts this tradeoff: protecting broad elevations necessarily leaves
some broad technical events uncorrected. Do not add an independent length rule that defeats
that protection. If mitochondrial mismapping and an amplification produce the same coverage
pattern, coverage alone cannot distinguish their causes. Document which spatial patterns
the tool corrects and preserves rather than claiming biological classification.

## 4. Poisson at high coverage: the concern is a varying baseline

**Claude's follow-up explanation:** If position-specific rates vary, sampling variance grows
with the mean while variation among rates adds a term growing with the squared mean. It
assumed relative rate variation of 20-30%, without measuring it. It suggested negative
binomial or zero-inflated negative binomial fitting, or a warning based on existing diagnostics.

**Codex assessment:** Poisson probabilities do not become mathematically invalid at 200x.
Under a Poisson fragment-start process, coverage at a position conditional on its expected
rate is Poisson. Neighboring positions are dependent because fragments span multiple bases.
Pooling positions with different rates generally gives a broader marginal distribution.

Let `position_rate` vary among ordinary active positions with mean `mean_rate` and relative
standard deviation `rate_cv`. Conditional on that rate, let coverage be Poisson. The law of
total variance gives:

```text
E[coverage] = mean_rate
Var(coverage) = E[position_rate] + Var(position_rate)
              = mean_rate + mean_rate^2 * rate_cv^2
variance / Poisson variance = 1 + mean_rate * rate_cv^2
```

For Claude's hypothetical `rate_cv = 0.2`, the total variance is 1.02 times Poisson variance
at mean 0.5 and 9 times at mean 200. This follows from the assumptions. It does not measure
the user's data or establish a universal failure depth. `mean_rate` here is the command's
positional fragment coverage, which must not be silently equated with nominal sequencing depth.

Claude's statement that the tail is approximately right at 0.5x also needs qualification.
A small change in bulk variance does not guarantee a well-calibrated extreme tail. Local
fitting can account for baseline differences between contexts, but a shared rate still cannot
describe variation among positions within a context.

Claude's first response reported ordinary-position call fractions of 0.4-2.9% at 100x,
2.2-7.6% at 200x, and at most 0.15% at 30x or below. Its follow-up quoted `T2 = 1.47 * lambda`
at 200x and described 1.5-fold GC elevations as routine. Our subsequent population calculations
reproduce threshold 294 at mean 200 and gamma-mixture ordinary tail fractions of 2.23-7.63%
with the explicit genome-scale probability used in our study. This supports those numbers
under stated assumptions. It does not establish a universal 30x boundary or verify the claim
about routine GC elevations. Claude's original simulation code remains unavailable.

### Why multiplying the probability threshold by squared depth is not the fix

The squared term describes extra variance, not a probability rescaling rule. In the current
algorithm, a position is called when its inclusive fitted upper-tail probability is sufficiently
small. Increasing the permitted tail probability lowers the coverage threshold and can make
ordinary elevations easier to call. Decreasing it raises the threshold.

Under a normal approximation, widening the variance would widen the excess over the mean
by a factor resembling `sqrt(1 + mean_rate * rate_cv^2)`. This is an explanation of direction,
not an exact tail calculation or a proposed implementation. Gamma and lognormal rate mixtures
can have the same mean and variance but different extreme tails. A universal multiplier on
the probability cutoff cannot be derived from the variance formula alone.

### Model and threshold candidates

| Candidate | What it addresses | What remains to establish |
| --- | --- | --- |
| Current ZIP | Sampling variation plus structural zeros | Whether positive coverage tails are adequately calibrated across the required range |
| Calibrated threshold adjustment with ZIP fitting | Raises the calling threshold using an explicit dispersion or empirical calibration rule | Tail calibration, sensitivity, and the consequences of still fitting/excluding with a misspecified distribution |
| Negative binomial | A broader positive-count distribution, with a fitted dispersion parameter | Robust fitting without letting pileups or amplifications define ordinary variation |
| Zero-inflated negative binomial, preferred working model | Broader counts plus a separate structural-zero component | Robust dispersion estimation, weaker-artifact sensitivity, tail assumptions, and stable fitting at low depth |

A negative binomial with mean `mean_rate` can have variance
`mean_rate + dispersion * mean_rate^2`. A gamma mixture of Poisson rates gives this form.
It provides a candidate probability distribution rather than merely changing a cutoff.
[Stan's distribution reference](https://mc-stan.org/docs/functions-reference/unbounded_discrete_distributions.html#negative-binomial-distribution-alternative-parameterization)
documents the equivalent mean/precision parameterization. This mathematical model does not
establish that a gamma rate distribution describes the residual coverage in this application.

For ZINB, T1, the right-truncated refit, T2, structural zeros, and target reporting must be
treated coherently. The research comparison now does this. Production integration must use
the appropriate likelihood and survival function rather than the truncated-Poisson equations.
At low depth, dispersion and structural-zero variation can be difficult to distinguish.
At high depth, a flexible fit can absorb genuine artifacts and reduce detection sensitivity.

### Existing diagnostics are evidence, not an automatic solution

[diagnose_zip_fit](../../../src/commands/outliers/model.rs#L168) compares the variance among
original positive observations with the fitted zero-truncated Poisson variance. At sufficiently
high component coverage, with negligible Poisson zeros and a correctly identified ordinary
mean, its ratio can resemble `1 + mean_rate * rate_cv^2`. That approximation is unsuitable
as an exact low-depth dispersion estimate.

The observed variance includes candidate pileups and biological elevations. Local contexts
overlap, so an event can affect multiple model rows. A median across rows may be useful, but
does not isolate GC variation or prove a particular cause. Observed-versus-fitted tails should
also be inspected. A warning reports a problem but does not make the default correction valid.

**Working direction:** Establish a general operating range under explicit assumptions.
The goal is not calibration to a particular user's dataset. Combine theoretical derivations
with simulations across 0.1x-200x, dispersion, event strengths, and event widths. A universal
failure depth cannot be inferred from an unmeasured CV, but that does not remove the need
for general theoretical assurances. ZINB is now the preferred working model, with the fitting
and preservation decisions below still open.

Two assurances should be assessed separately:

- Under correctly calibrated ordinary-position tail probabilities, the expected number of
  falsely called positions is bounded by their sum. Independence is not needed for this
  expectation. Estimated parameters, truncation, and misspecification need explicit treatment.
- For biological elevations, derive the preservation probability and signal loss as functions
  of depth, event strength and extent, ordinary rate variation, and fragment spans.

These establish an operating range rather than a promise to distinguish observationally
identical technical and biological patterns. A model change requires rethinking the complete
fitting and correction process, not just exchanging the T2 calculation.

## 5. Proposed simulation study, without a 200x BAM

Simulation can expose theoretical failures without a deep real sample. It establishes behavior
under stated generating assumptions, not empirical validity for every 200x library. This is a
separate investigation authorized by the user on 2026-10-13. Keep source, settings, and
results in `outliers_fragment_simulations/`. The study does not change production behavior.

Use two complementary levels:

1. **Histogram-level distributions:** Poisson/ZIP as correctly specified controls, then gamma
   and lognormal Poisson mixtures with optional structural zeros. Check fit recovery and tail
   thresholds against known generating distributions. Independent draws at this level do not
   establish spatial call behavior or full-fragment correction.
2. **Fragment-level simulations:** Generate fragment starts and accepted lengths, derive raw
   coverage from their `pos` to `reference_end` spans, and run actual detection and weighting.
   Use spatial rate variation so neighboring coverage and fragmented calls arise from shared
   fragments. Include GC-associated rate variation and corresponding fragment GC weights.

Suggested scenarios, with values treated as sensitivity cases rather than measured biology:

- Positional component means spanning 0.1, 0.5, 1, 5, 17, 30, 100, and 200.
- Rate CV including zero and values below and around Claude's assumed 0.2-0.3. Vary spatial
  scale from hundreds of bases through kilobases and wider blocks. Compare more than a gamma
  mixture so a candidate is not assessed only against its own assumed generating model.
- Ordinary coverage without artifacts, narrow technical additions, moderate and strong copy
  gains, broad gains, and a sharp spike within an elevated region. Revisit the reported
  2.35-fold/1 Mb/mean-17 example with both a fixed oracle baseline and the actual two-stage fit.
- Event positions relative to the 5.5 Mb fitting contexts, model cores, tile boundaries,
  chromosome ends, large masks, and short eligible islands. Include isolated and clustered calls.
- The reported 167 bp fragments, then varying accepted lengths and overlapping calls. Respect
  actual command length limits. Do not use fragment lengths below 10 bp, or below an enabled
  inclusion minimum such as the usual 30 bp.
- No GC weights, constant GC weights as a control, and varying GC weights associated with
  simulated rate bias. Compare no-outlier, GC-only, outlier-only, and combined outputs.
  Define the desired GC-corrected baseline from the generating process rather than assuming
  the combined output should numerically equal raw `lambda`.

For the full workflow, evaluate every call-start window, propagate each window's evidence
to all overlapping calls, flank and merge all candidates, then apply the shared decision.
Do not substitute the mean percentage of an event for those maxima. Use the actual fitting
and lookup implementation where possible rather than an unverified rewrite of its formulas.

Report fitted parameters, T1/T2, exact call fractions, window maxima and support, preserved
merged regions, affected fragments, signal removed inside and outside calls, biological gain
retention, technical excess removed, and background distortion. Retain spatial plots so the
mechanism can be explained, not just a final call percentage.

Use reproducible seeds and replicates, with variability across replicates. Keep simulation
code and settings with the results. Select tolerable background distortion and amplification
loss before comparing models. Do not infer acceptable performance from a low called fraction.

The automatic cutoff depends on the complete eligible-position count. A small simulated
chromosome gives a different cutoff from a whole-genome run. Report that count and either
use the intended run's explicit probability or explain why the smaller automatic cutoff is
the scenario being studied. Ordinary Monte Carlo over a few million positions cannot validate
an extreme probability around `1e-10`. Use analytic tail calculations for known distributions
alongside simulations of the practically important failure modes.

Upsampling or duplicating a 0.5x BAM would reuse molecules and can create artificial pileups.
It should not be presented as an independent 200x validation dataset.

## 6. Fragment ownership and input checks: keep the scope small

**Claude's proposal:** Count excluded and downweighted fragments by fragment-start ownership
to avoid duplicate counts across processing tiles. It also proposed additional metadata and
sample-name checks.

**Codex clarification:** Ownership here was a reporting convention, not detection of fragment
starts, a correction rule, or a new biological criterion. Fetch halos can expose a fragment
to multiple tiles. Existing `counted_fragments` increments when that fragment contributes to
a tile core. The existing start-based helper is specifically for normalization statistics.

Reusing it indiscriminately would be wrong in BED-restricted modes, where a fragment's start
tile can be skipped while a neighboring requested tile counts its overlapping contribution.
This follows from [window-based tile skipping](../../../src/commands/fcoverage/tiling.rs#L324)
and [the normalization-specific helper](../../../src/commands/fcoverage/fcoverage.rs#L1728).
It is a reason to avoid promising new unique-fragment summaries, not a request to redesign
tiling or normalization.

**Direction for the revised plan:** Keep existing accounting conventions and label any new
summaries accurately. Do not add a fragment-start ownership scheme or a unique-counting
refactor. Implement early zero-weight skipping in each existing GC/no-GC branch. Fragments
skipped before GC handling will naturally not increment GC-failure counters afterward.

For input checks, retain existing required columns, coordinate/order validation, and finite
weights in `[0, 1]`. Add the agreed reference names/lengths and analysed-chromosome checks,
with a concise warning for missing information in custom tables and useful filter-setting
warnings. Do not compare blacklists, BAM identity/content, `@RG SM`, or private sample names.
Do not turn hypothetical combinations of metadata into a large new checking system.

Recording only relevant analysed-contig names and lengths in the weight file is a reasonable
way to keep this information small. The representation already exists as
[ContigFootprintEntry](../../../src/shared/reference.rs#L13). It can be populated from header
information without adding a reference-sequence input. Exact metadata layout remains to be
settled when the plan is revised.

## 7. Other findings and corrections to Claude's response

| Finding or proposal | Codex assessment and disposition |
| --- | --- |
| O-001 conditional-expectation argument | Claude acknowledged the error. `E[X given X > 0]` is the zero-truncated mean, not lambda. Keep the agreed active-component baseline argument for lambda. |
| f32 boundary guard | Claude's correction is valid. Although `2^24` is representable, `2^24 + 1` can round to it. The stored boundary must be rejected because the original integer is ambiguous. Clarify the explanation without changing the guard. |
| Narrowing O-012 through O-014 | Claude accepted the earlier scope corrections. Retain the small terminology, lookup-table, and missing-behavior coverage work unless another decision changes it. |
| Assembly gaps counted as zeros | Eligible positions are determined by masking, not by inspecting reference bases. Recommend adding an assembly-gap BED to the existing blacklist when appropriate. Do not require a new reference input each run. The actual supplied hg38 blacklist has not been inspected for gap coverage. |
| One-time gap BED preparation | Confirmed: UCSC `twoBitInfo -nBed` outputs N intervals from a 2bit file. This can be preparation per assembly, not a new per-run dependency. [UCSC source](https://github.com/ucscGenomeBrowser/kent/blob/master/src/utils/twoBitInfo/twoBitInfo.c) |
| Female chrY and other biologically absent regions | A lambda target does not classify the cause of zero coverage or establish that an artifact is at an active locus. Explain this interpretation limit when users include such chromosomes. Do not infer biological presence from fitted pi alone. |
| Fragment can overlap preserved and corrected calls with small flanks | The span argument is valid. A fragment can have a positive overlap with both calls across a gap no larger than `maximum_fragment_length - 2`. Default flanking covers this when detection and counting length bounds agree. |
| Group across `max(2 * flank, maximum_fragment_length - 2)` | This is an additional grouping proposal, not agreed behavior. It changes what a small or zero flank means. Do not add fragment-level overrides. If grouping is reconsidered, distinguish sharing a decision from filling extra genomic gaps in a blacklist. |
| No warning-capture support exists | Incorrect. [test_run_options.rs](../../../tests/test_run_options.rs#L88) captures tracing output with a subscriber that also admits warnings. A new dependency or checking framework is not required by that claim. |
| Every integration fixture will hit the safeguard | All fixtures need review, but automatic preservation does not follow from their sizes alone. The simple 30 bp call has 900 bp to the contig end, giving about 3.33%. The 150 bp roundtrip call has 5,600 bp, giving about 2.68%. A future minimum-support rule can still affect these fixtures. Neither percentage exceeds 10%. |
| Final boundary call is preserved | Under the naive clipped-denominator rule, yes: 30 called bases out of 40 yields 75%. Its disposition changes with the support policy discussed above. |
| Roundtrip expectation changes | The local-mean weight and resulting coverage change from the unconditional target near 2.695 to lambda near 2.82. The diagnostic `underlying_mean` retains its definition and is not replaced by lambda. |
| Audit layout is underspecified | Valid. The revised plan should fix a filename, row identity/order, decision reasons, proposed weight, actual window support, and merged-region membership. Existing masses and candidate mean/max coverage can support later analysis. A separate preserved-region BED is optional, not decided. |
| Per-call lambda/pi/T2 audit columns | A connected exact call can cross model cores and contain multiple parameter triples. Reference the relevant model rows or define a summary. Do not imply every call has a single fit. [Call merging](../../../src/commands/outliers/outliers.rs#L675) retains connected intervals and summed masses. |
| Error action publishes no diagnostics | Failure before final publication would publish no new diagnostic files. That does not invalidate an explicit strict pipeline action. Decide its reporting contract rather than automatically dropping it or broadening output publication. |
| Target-definition metadata | Useful: make the meaning of the local-mean target explicit, alongside preservation settings and withheld-call summaries. This explains the output mathematically and is not a request for compatibility behavior. |
| Argument vocabulary and CLI help | Valid housekeeping. Explain what the percentage counts, window support, contiguous-call preservation, pipeline defaults, and full-fragment weighting. Settle names after the decisions. |
| Generated documentation and result paths | When implemented, regenerate CLI documentation, remove the existing unwanted metadata term from help, correct the claim that multiplicative factors must be independent, and include audit artifacts in output/result paths. Do not edit them ahead of settled behavior. |

## 8. What still needs to be settled

The first BAM-inspection version now makes the initial choices: ZINB by default with a ZIP option,
active-component targeting, call-start 100 kb windows, 50 kb minimum support/denominator floor,
a strict 5% preservation cutoff, and shared preservation after flanking/merging. Unsupported
boundary calls try a backward window ending at the exact call's end. A sufficiently supported
fallback replaces the forward score. Otherwise retain the forward maximum and floor.
These choices allow real inspection rather than claiming the scientific questions are resolved.

1. **BAM findings and priority:** Inspect actionable/preserved calls and matched raw/weighted
   coverage. Measure analysed genomic extent separately from removed coverage. Decide whether
   the observed benefit warrants further work before the application-note deadline.
2. **Biological retention and technical sensitivity:** Establish acceptable loss and useful event
   scope across depths, including narrow gains, moderate artifacts, and spikes within gains.
   The earlier 10% simulation results remain historical. They are not measurements of the new 5% rule.
3. **Preservation evidence:** Assess the initial density heuristic and boundary fallback. Broader
   regional coverage evidence remains a later algorithm decision, rather than part of this version.
4. **Fit and probability calibration:** Inspect fitted dispersion/boundary use, low-depth parameter
   uncertainty, convergence, and observed-versus-fitted tails. Compiled code and a favorable BAM
   example do not settle the extreme-tail assumption across all libraries and depths.
5. **CLI and documentation review:** The initial arguments, model/audit fields, summaries, and
   local inspection workflow exist. The user will finalize names and review documentation.
6. **Merge readiness:** Evaluate the first BAM run before broader validation or rollout. Midpoints,
   calibration producers, and normalization changes remain outside this milestone.

The existing whole-fragment minimum-weight semantics, focused reference/chromosome checks,
no BAM/blacklist comparisons, and user-owned normalization decisions remain unchanged.

## 9. Simulation findings now available

The [full report](outliers_fragment_simulations/results_report.md) contains the methods,
assumptions, tables, and figures. Rates below are positional fragment coverage, not nominal
sequenced-base depth.

- The frozen-baseline mean-17, 2.35-fold gain reproduces approximately 1.70-fold output without
  preservation. With complete window maxima and merging it averages 2.09-fold. Actual context
  fitting gives about 2.22-fold, with 5.48% gain signal removed. The original claimed 1.7 result
  must not be presented as the complete proposed workflow's prediction.
- Sparse calls still remove meaningful biological signal. In that fitted case, about 83% of
  the removed coverage lies outside exact called bases. A 10 kb gain at mean 17/CV 0.2 loses
  about 31%. A strong shallow-depth stress case also loses signal. The tested 10% cutoff is
  not yet an established protection rule across the requested operating range.
- With ordinary lognormal rate variation at mean 200/CV 0.2, fitted ZIP removes about 2.82%
  of ordinary raw background signal despite the screen. This supports the general model
  mismatch concern under explicit assumptions, rather than relying on a particular dataset.
- The denominator floor corrects the short-window preservation bias in the tested narrow
  additions. The 30/50/100 kb alternatives behave differently for a few thousand called bases.
  Selecting the minimum still requires the intended spatial scope.
- Broader oracle thresholds improve some outcomes. They do not establish fitted NB/ZINB
  behavior or calibrated tails for all rate distributions. Population calculations show
  gamma and lognormal mixtures with the same mean and variance have different extreme tails.
- Multiplicative GC weighting behaves as intended in this research application. Consistent
  estimation and downstream weighting remain the useful workflow recommendation. There is
  no new raw-coverage requirement or normalization change.

The fitted follow-up below supplies the coherent broader-model comparison. The protection
scope, signal-loss criteria, and remaining fitting policies are listed in section 8.
No new fragment-mass screen is selected.

## 10. Fitted ZINB and smaller-CV follow-up

The [comparison report](outliers_fragment_simulations/zinb_comparison_report.md) describes
684 paired fragment simulations. ZINB estimates its own component mean, dispersion, and
zero inflation, including a correctly normalized right-truncated refit. The comparison
adds CV 0.05 and 0.1 and selected 200 bp and 10 kb rate blocks.

At ordinary mean 200 with 1 kb lognormal rates, ZIP removes about 0.0004%, 0.09%, and 2.8%
of ordinary signal at CV 0.05, 0.1, and 0.2. ZINB removes none in those finite replicates.
The CV scenarios remain unmeasured assumptions, so the larger loss is not a typical-library
estimate. Both methods preserve the same full-fragment and multiplicative GC semantics.

ZINB preserves the tested 1 Mb gains but still loses about 18% of a 10 kb, 2.35-fold gain at
mean 17/CV 0.1. The larger regressions are at mean 200: ZIP preserves a 10 kb gain while
ZINB removes about 58% at CV 0.1 and 26% at CV 0.2. ZINB also removes about 5.4% from a
100 kb gain at mean 200/CV 0.2 that ZIP preserves. In the first 10 kb/CV 0.1 replicate,
window scores are 10.063% under ZIP and 9.983% under ZINB. The provisional safeguard turns
that small difference into preservation versus correction. These are model/safeguard
interactions, not evidence that ZIP is a better model for ordinary rate variation.

At mean 200/CV 0.1 inside a broad gain, ZINB removes about 95% of a 5,000-fragment technical
spike while losing about 0.032% of biological gain signal. ZIP preserves that spike with
the gain. Both leave a smaller 500-fragment spike untouched. Independent population-tail
calculations still expose the NB tail assumption for lognormal rates.

**Overall preference:** Use ZINB as the working model for the next design iteration. It
includes the Poisson limit, handles simulated ordinary variation across the required range,
and retains broadly comparable correction of the tested isolated spikes. The genuine costs
are more complicated fitting and possible loss of sensitivity when technical excess defines
the estimated dispersion. Low-depth separation of dispersion and zero inflation also needs
a fitting policy. The simulations do not establish universal tail calibration or sensitivity
over all artifact strengths. Production code is unchanged, and the final preservation rule
and fitting details remain open.
