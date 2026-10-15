# Outliers review verification and rectification plan

Date: 2026-10-12. Updated: 2026-10-15.

Status: The first BAM-inspection milestone is implemented. Project tests are authored and
compiled, not executed by the agent. The historical review verification below remains the
record of the original branch state, not a claim that all 14 items have been completed.

## First BAM-inspection milestone

Current choices supersede earlier open defaults below: ZINB is the default, `--model zip` retains
the analytic ZIP fitter, local-mean targets the active component, and both means are reported.
Preservation uses call-start 100 kb physical windows, a 50 kb support/denominator floor, a strict
5% cutoff, and preserve action. An unsupported boundary call tries the backward window ending
at its exact end. Known exact intervals supply called bases. A sufficiently supported backward
window replaces undersupported forward evidence. Otherwise retain the forward maximum and floor.
Flank all candidates, merge, and apply a shared decision to every member.

The command now writes complete audit and genomic-summary TSVs, and fcoverage skips zero weights
before correction/counting while retaining normalization formulas. Recorded reference and analysed
chromosome checks are focused. BAM identities and blacklists are not compared. The guide includes
matched raw/weighted commands and `scripts/plot_outliers.R` for local PNG/PDF inspection.

This milestone stops at BAM inspection. Remaining scientific validation, smarter regional evidence,
other consumers, calibration producers, normalization changes, and merge readiness remain later work.
Numerical compilation does not establish biological validity on the user's BAM.

Verified branch: `outliers2`, HEAD `f5a75211`.

Source review: [Outliers Command Review](../major_reviews/outliers.md).

This plan checks O-001 through O-014 against the current implementation and separates
implementation defects from scientific decisions, future integration work, and wording issues.
The original review is preserved. The proposed scientific changes below are not current specs.

The [working notes](outliers_claude_review_working_notes.md#8-what-still-needs-to-be-settled)
record the remaining validation and prioritization decisions. The
[fitted comparison](outliers_fragment_simulations/zinb_comparison_report.md) supplies the earlier
simulation evidence, including failures under its historical 10% cutoff. The current first
inspection version uses the initial choices stated above. It does not establish universal
calibration or complete all remaining review work.

The intended workflow uses existing blacklists first, then detects the remaining sample-specific
extremes. This command supplements blacklisting. Current data are about 0.5x, and the required
operating range is 0.1x through 200x. This is a validation requirement, not an established result.
Nominal sequencing depth and the command's counted positional fragment coverage are different
quantities and must both be recorded in comparisons.

## Working implementation scope

These decisions reflect the discussion with the user. They supersede the earlier suggestions
in the review. The first milestone above is implemented. Other review items below remain a
record of proposed work and must not be inferred complete from this milestone.

- Use the fitted active-component mean for `--target local-mean`. This is `lambda` under ZIP
  and the negative-binomial component mean under ZINB. Continue reporting the component and
  unconditional means.
- Use ZINB by default and keep ZIP selectable. Keep the two-stage structure with its own initial fit,
  T1 exclusion, correctly normalized right-truncated refit, and T2 calculation. Preservation
  is decided after final calling and does not put observations back into the baseline fit.
- Start a 100 kb window at each exact call's start. Assign each call the maximum called-base
  score across evaluated forward windows it overlaps, with the supported backward fallback.
  The initial cutoff is 5%, with a 50 kb support/denominator floor.
- Calculate percentages before flanking. Flank all candidate calls, then merge overlapping
  and touching intervals. If any member requires preservation, preserve the entire merged
  region and omit all its member calls from correction and actionable blacklist outputs.
- Resolve a fragment's outlier keep weight before counting. A zero weight skips the fragment
  before GC correction, coverage insertion, and counted-fragment accounting.
- Keep normalization formulas and ownership rules unchanged. Weighting restore-mean statistics
  is excluded from this work and may be considered later only when the user explicitly asks.
- Reuse the existing encoded reference contig names and lengths, alongside analysed chromosome
  names. Known reference mismatches and unanalysed requested chromosomes are errors. Differences
  in detection/counting settings are warnings. Accept custom tables without recorded settings
  and warn that those checks could not be performed. No override option is required.
- CLI help is the primary documentation. Explain the targeted residual high-coverage pileups,
  possible biological causes, automatic preservation, complete-fragment weighting, and controls.
- Inspect the initial BAM outputs before broader real-data validation. The initial model and
  preservation choices are made. Their usefulness and biological retention remain to assess.
- Keep changes focused on outliers and its current fcoverage integration. Producer rollout and
  wider consumers remain future work. Broader refactors require separate user consent.

## What still needs to be settled

The [full remaining list](outliers_claude_review_working_notes.md#8-what-still-needs-to-be-settled)
now starts with BAM inspection and prioritization. Initial defaults and the backward fallback are
implemented. Remaining work concerns biological retention, moderate-artifact sensitivity,
calibration and uncertainty, possible smarter preservation evidence, CLI naming/documentation
review, and eventual merge readiness. These are not blockers for the initial inspection run.

## Assessment of all 14 findings

### O-001: The default target includes the structural-zero fraction

**Confirmed behavior, qualified review argument. The lambda target is now agreed.**

[Calling](../../../src/commands/outliers/outliers.rs#L630) uses
`underlying_fit.mean()`, which [the model](../../../src/commands/outliers/model.rs#L65)
defines as `(1 - pi) * lambda`. It is the correct unconditional ZIP mean, not an algebra bug.
For fixed `lambda`, adding structural zeros to the context lowers this target by exactly
`1 - pi`. A called mappable locus can therefore be corrected below its active-locus baseline.

The review's conditional-expectation argument is incorrect as stated. Under an uncontaminated
ZIP, observing `X > 0` identifies the Poisson component, but also conditions its count on being
positive. Thus:

```text
E[X] = (1 - pi) * lambda
E[X | Poisson component] = lambda
E[X | X > 0] = lambda / (1 - exp(-lambda))
```

At a threshold-crossing position, `E[X | X >= T2]` is another truncated mean. None of these
conditioning facts alone specifies the desired replacement for a contaminated observation.
An artifact can also add positive coverage at a location whose true baseline is structurally
zero. Fitted `pi` does not identify assembly gaps, unmappable positions, or deletions individually.

The distinction matters at shallow depth. At `lambda = 0.5`, Poisson zeros have probability
about 0.607 and the positive conditional mean is about 1.271. Those ordinary sampling zeros
are already part of the Poisson component and must not be described as structural zeros.
Using `lambda` to restore the baseline at an active locus does not mean its observed positive
counts have mean `lambda`. Validate both fitting and correction from 0.1x through 200x, with
0.5x as the immediate priority. Do not rely on a high-depth approximation to justify a target.

The two means keep these definitions at every depth. Do not choose lambda at high coverage
and the unconditional ZIP mean at low coverage solely because sampling zeros become common.
When `pi = 0`, both targets already equal lambda, including at shallow depth. For `lambda = 0.5`
and `pi = 0.2`, the active-component baseline remains 0.5 while the whole-context mean is 0.4.
The reduction to 0.4 represents the extra structural-zero component, not ordinary Poisson zeros.
At low depth, pi and lambda can be less clearly separated by the observed distribution because
multi-count observations are rarer. Assess that uncertainty and model fit instead of changing
the biological target definition according to depth.

The [roundtrip fixture](../../../tests/test_outliers_command.rs#L429) supports the numerical
claim. After exclusion it has 1,000 zeros and 8,850 values of 3. Its positive mean identifies
`lambda` near 2.82, `pi` is near 0.045, and the unconditional target is near 2.695. This verifies
the implemented target, not which target is biologically appropriate.

**Agreed change:** Use `underlying_fit.lambda` for the local-mean correction target, with the
same definition from 0.1x through 200x. Restore the Poisson-component baseline rather than
diluting it by the fitted structural-zero fraction. This is a theoretical decision. The
real-data comparisons are later validation, not a prerequisite for changing this target.

Do not redefine `ZipParameters::mean()`. It must remain the unconditional ZIP mean.
Both metrics are already written as `underlying_lambda` and `underlying_mean` in the model
TSV and global-fit metadata. Keep both, with their existing mathematical meanings.

Update target selection, target metadata/help, affected regional-weight expectations, and
detector documentation together. CLI help must explain that local-mean now uses lambda.
Update current specs only when the behavior is implemented. Keep fits, both thresholds,
the regional mass ratio, and the fragment overlap rule unchanged.

The mixture definition is supported by [Lambert (1992)](https://www.stat.cmu.edu/technometrics/90-00/vol-34-01/v3401001.pdf).
The conditional means above follow directly from that definition.

### O-002: Calibration commands do not consume the outlier weights

**Confirmed future rollout dependency, not a defect in the intentionally limited rollout.**

[Scaling configuration](../../../src/commands/coverage_weights/scaling_weights_config.rs)
has no outlier input, and [its internal fcoverage configuration](../../../src/commands/coverage_weights/coverage_weights.rs#L425)
never sets one. The commands use average coverage or length-normalized total fragment mass.
[Scaling normalization](../../../src/commands/coverage_weights/striding.rs#L226) uses
`global_mean / smoothed_value`. Both local smoothing and the global mean can be contaminated.
[GC counting](../../../src/commands/gc_bias/gc_bias.rs#L1170) also lacks this input.

Only fcoverage integration has been implemented so far. Missing inputs in the other commands
are expected at this stage and should not be counted as a current implementation bug. Existing
blacklists are intended to remove known problematic positions before residual-outlier detection.
The automatic denominator already excludes input-blacklisted bases, as do the fit histograms.
Whether the remaining outliers materially affect calibration needs measurement after that masking.

The illustrative 500 kb coverage calculation is sound. Replacing background 17 with coverage
5,000 over 2 kb gives `17 + (5000 - 17) * 2000 / 500000 = 36.932` before smoothing. Such an
excess can lower the scaling factor and depress the corrected bin and its neighbors after the
outlier fragments are downweighted. The direction and magnitude of GC distortion depend on
the excess fragments' length/GC distribution, window normalization, and matrix processing.
It is not necessarily a second uniform reduction at the pileup. Multiplicative application
does not require statistical independence, and matching outlier inputs alone does not prove it.

**Plan:** Keep support for `gc-bias`, `coverage-weights`, and `fragment-count-weights` in the
future rollout. First compare factors estimated with and without the residual outlier correction,
using the same existing blacklist. Measure the excess mass relative to each smoothing context
and length/GC stratum, including less sparse alternative use cases. Do not make this a blocker
for discussing or validating the current detector, or automatically reject existing combinations
of corrections without evidence about their effect.

When weighted calibration is added, the proposed order is existing blacklist -> raw residual
detection -> outlier-weighted GC estimation -> outlier- and GC-weighted scaling estimation
-> feature counting with matching corrections. Record weight-file identity, including a content
digest, in calibration outputs. Decide which recorded correction differences warrant errors or warnings
from the supported workflows and measurements. GC-tag inputs need explicit treatment because
an AUX tag does not carry the calibration settings themselves.

For `gc-bias`, apply the keep weight to observed fragment counts in both fixed-size streaming
and BED/global branches, in addition to the existing window-assignment fraction. Recalculate
derived length frequencies from those counts. Do not weight the reference distribution as
though it were observed fragments. Review zero-support windows and reference masking separately.

### O-003: Focal biological amplifications can be suppressed

**Confirmed scientific risk and missing run summaries. The proposed length limit is incomplete.**

[Context expansion](../../../src/commands/outliers/striding.rs#L204) normally produces 5.5 Mb
centered contexts with the defaults, while an edge context can be 5 Mb. The two-stage fit has
no causal or copy-number classification. The idealized plateau example is consistent with it:
for a 5.5 Mb context containing a 1 Mb plateau at 93 and background 17, the initial mean is
about 30.82 and `T1` is about 72. Removing the plateau gives a baseline near 17 and `T2` near
49 at `alpha = 3.5e-10`. The plateau can receive a weight near `17 / 93 = 0.183`.

With 3 Mb at 93 in the same context, the initial mean is about 58.45 and the passing threshold
is near 112, rather than exactly 110. The plateau remains in the second fit and can avoid
calling. These are hand-derived idealized comparisons, not BAM experiments. Actual calls
depend on context placement, variation in coverage, and the eligible-position denominator.
The automatic probability is not a fixed human-genome constant.

Coverage 93 versus 17 is a roughly 5.47-fold elevation. In a simplified diploid mixture with
tumor contribution `f` and tumor locus copy number `C`, the expected relative coverage is
`R = 1 + f * (C / 2 - 1)`. A ratio of 93/17 corresponds to about 11 copies at `f = 1`,
20 copies at `f = 0.5`, or 91 copies at `f = 0.1`. It is a strong-amplification example,
not a typical event at every tumor fraction. High copy numbers are biologically possible,
including the EGFR examples in [Turner et al. (2017)](https://www.nature.com/articles/nature21356).
That evidence does not establish the frequency of this exact plasma scenario.

The earlier 20-copy example means about 20 copies in the tumor component, not 20-fold plasma
coverage. Direct plasma evidence also exists. The selected HER2-amplified
[HERALD trial](https://pubmed.ncbi.nlm.nih.gov/39088783/) reported baseline plasma HER2 copy
numbers with median 8.55 and range 2.4 through 73.9. This supports the possibility of strong
plasma signals, not their prevalence in unselected samples or a flat 1 Mb coverage plateau.

The first fit has already moved lambda from 17 to about 31, but 93 still exceeds `T1` near
72. Removing the second stage would therefore not protect this idealized plateau. It would
retain a higher threshold and target, reducing severity rather than eliminating the call.
The second-stage exclusion assumes that the tail should not define the background. Accounting
for truncation makes the fit mathematically coherent under that assumption, but does not prove
that a removed component is technical noise rather than a genuine alternative local baseline.
For stochastic amplicon coverage, the first exclusion removes only observations at or above T1.
Lower amplicon observations can remain, so recovery of exactly lambda 17 is only an idealization.

At shallow depth, the same relative amplification must be examined with discrete observed
counts and depth-specific thresholds. A 5.47-fold increase over a mean of 0.5 does not reproduce
the deterministic 17/93 example. Compare one-stage fitting, the present two-stage fitting, and
spatial evidence for a distinct biological baseline before deciding whether exclusion is right.

[Run statistics](../../../src/commands/outliers/outliers.rs#L337) report the number of calls,
but not total called bases or longest region. Exact intervals are already available in the TSV
and BED outputs, so lengths are inspectable. They are not summarized or highlighted.

Real focal amplifications are observed in plasma sequencing, including in
[Ulz et al. (2016)](https://www.nature.com/articles/ncomms12008.pdf) and
[Chabon et al. (2016)](https://www.nature.com/articles/ncomms11815).
Those studies establish biological relevance. They do not establish how often this detector
would suppress such events across assays or depths.

**Working change:** Add called-base totals, longest exact call, corrected/preserved call counts,
and merged-region summaries. Retain all candidates in an audit output. Settle preservation
evidence and its cutoff with ZINB rather than assuming the tested 10% setting is final.
Do not add an independent maximum-region-length policy that defeats broad-event protection.

**Two-stage structure retained, model details pending:** Favor ZINB with its own likelihood and
survival probabilities, and decide preservation after final calling. A biological amplification
can be preserved without defining ordinary background coverage. Exclude values at or above T1,
not entire amplified regions. Observations below T1 still contribute. The right-truncated
normalization must match the chosen model. Do not selectively add preserved high values back
or introduce an iterative baseline refit.

The implementation order is:

1. Fit the complete eligible histograms and calculate global and local T1 thresholds.
2. Exclude values at or above each fit's T1, refit with the chosen model's right-truncated likelihood,
   and calculate T2.
3. Call against T2 on original coverage and join exact calls across cores and tiles.
4. Evaluate all call-start windows on the complete exact call set. Assign each exact call its
   maximum window score and actual support. The final preservation evidence remains to be settled.
5. Flank every candidate call with the configured flank, clip to chromosome bounds, and merge
   overlapping and touching flanked intervals. Retain which exact calls belong to each region.
6. Apply the ambiguity action to each merged region. With the default preserve action, any
   member requiring preservation protects the entire merged region and all of its member calls.
7. Write all candidates and decisions to the audit table. Write only selected exact calls to
   the correction TSV and exact BED, and only selected merged regions to the flanked BED.

**Call-start window calculation:** Start a window at every exact T2 call's `start`. The agreed
default is 100 kb. Calculate:

```text
called_score(window) = 100 * T2-called unblacklisted bases / max(unblacklisted bases, minimum_eligible_bases)
outlier_score(call) = maximum called_score of windows overlapping the call
merged_region_score = maximum outlier_score among the member calls
```

Earlier calls' windows also update later overlapping calls. A call does not use only the
window it starts. The maximum is over these evaluated call-start windows, not every possible
genomic window. Evaluate all windows before filtering or flanking. Count the union of exact
called bases. Flanks and the gaps joined by merging never inflate the numerator.

Include unblacklisted zero-coverage positions in the denominator. Clip windows at chromosome
ends and record their actual span and eligible-position count. The working short-window rule
uses `called_bases / max(eligible_bases, minimum_eligible_bases)`, with the numeric minimum
still to be settled. Low support alone does not mark a call for preservation. This density
calculation uses exact calls and masking, but any additional coverage-based evidence needs
its own definition before implementation.

The tested percentage cutoff is 10%, exposed as a configurable candidate. It is not a settled
CLI default. Under the candidate comparison, a percentage strictly greater than the cutoff
marks a call as ambiguous, while equality does not. The final evidence and threshold must
be assessed with ZINB, including near-cutoff regressions. High percentages indicate broad or
clustered coverage elevations, not proven amplifications. Technical excess and model mismatch
can produce the same pattern.

Keep the current fitting-context definition. With the current defaults, centered fits use
5.5 Mb. The 100 kb window assesses extent and does not supply a correction baseline. A separate
baseline for correcting spikes within an amplification is additional fitting work and remains
outside this implementation.

**Preservation after flanking and merging:** Build merged regions from all candidate calls,
including calls individually marked for preservation. Do not discard those calls before
merging, since they may connect neighboring calls into the same region. A merged region is
ambiguous if any member meets the final preservation criterion. Low support alone is not
that criterion. This propagates
preservation through the complete merged region, including calls with a lower individual score.

The default action is to omit the whole ambiguous merged region from the flanked BED and all
of its exact member calls from the keep-weight TSV and exact BED. Keep them in the audit table
as preserved detection candidates. Do not subtract protected bases from neighboring flanks,
split a merged region around protected intervals, or introduce fragment-level protection
overrides. The existing whole-fragment minimum-weight rule remains unchanged.

For merged regions selected for correction, keep the original exact intervals and their
individual proposed weights in the correction TSV. Merging determines a shared action. It
does not pool regional masses or create a new weight over flanks and intervening gaps.

**Automatic handling and outputs:** Manual editing must not be required for pipeline use.
Write a full audit table with each exact interval, `proposed_keep_weight`, its own maximum
window score and location, actual eligible support, merged-region bounds and maximum score,
selected action, and decision reason. Explain when a call is preserved because another member
of its merged region is ambiguous. Using `proposed_keep_weight` rather than `keep_weight`
prevents the audit table from accidentally satisfying the correction loader's required columns.

Keep the correction TSV as the input consumed by `fcoverage`. Exact and flanked actionable
BED files must use the same merged-region decisions. Report raw candidates, selected calls,
preserved calls, and merged-region outcomes separately. Preserved candidates are not operational
outliers to remove. Their biological or technical cause remains unresolved.

**Arguments and CLI documentation:** Proposed names are implementation choices. The existing
flank argument also controls merged-region grouping and therefore preservation propagation.

| Argument | Purpose | Default |
| --- | --- | --- |
| `--protection-window-size` | Window span in bp, anchored at each exact T2 call's start | 100000 bp |
| `--max-outlier-percentage` | Candidate maximum T2-called-base percentage before ambiguity | Provisional 10%, final value to settle |
| `--min-protection-positions` | Minimum denominator for the candidate short-window score | Still to agree |
| `--ambiguous-region-action` | Handle an ambiguous merged region with preserve, correct, or error | preserve |
| `--blacklist-flank` | Existing flank used before merging candidate calls | Maximum accepted fragment length |

CLI help is the primary documentation. The command description and argument help must explain
what users need to know where they run the command:

- The intended workflow uses existing blacklists first, then detects residual sample-specific
  regions with extreme high raw positional fragment coverage. It targets isolated pileups.
- Extreme coverage can also be biological. The command does not identify causes from coverage
  alone. Broad or densely called regions are preserved by default.
- Zeros and low coverage influence the fit but are not themselves called as outliers.
- The local-mean correction target is the active-component mean. Both it and the unconditional
  mean are reported, together with the chosen model and its dispersion when applicable.
- Listed exact intervals apply the full regional weight to complete overlapping fragments.
- Under the density candidate, a call's score comes from all overlapping evaluated windows,
  including earlier calls' windows. It uses exact called bases before flanking and the
  configured denominator floor. Explain any additional preservation evidence when selected.
- Flanking and merging propagate preservation to all member calls in an ambiguous merged region.
- The fitting context, the protection window, and the flank have different purposes.
- The keep-weight TSV and actionable BED files are ready for pipeline use. The audit table
  records every candidate and explains withheld corrections.

Explain the effect of changing each argument. A lower percentage cutoff preserves more calls.
Changing the flank can join more calls into a region with a shared action. The preserve action
withholds the merged region's corrections, correct uses the existing proposed weights despite
ambiguity, and error stops rather than publishing actionable correction files. These actions
leave the baseline fit unchanged. The guide can illustrate the CLI explanations.

Use percentages on the 0 through 100 scale. A maximum of 100 disables percentage-based
preservation under that candidate rule. Validate finite percentages
in range, positive window sizes/support, and minimum positions no larger than the window.
Reject invalid settings before scanning with errors naming the argument and relevant values.

Expose matching programmatic configuration and CLI serialization. Record actual settings and
actions in metadata. Include hand-derived behavioral cases for earlier-window maxima, cutoff
equality, blacklisted denominators, chromosome-end clipping, overlapping/touching flanks,
transitive merges, and a low-percentage call sharing a region with a preserved call.

### O-004: The loader ignores recorded reference and detection settings

**Confirmed. The input/checking policy is now agreed.**

[The loader](../../../src/shared/outlier_weights.rs#L45) skips comments and validates only
selected chromosomes' interval rows and weights. A selected chromosome absent from the rows
uses identity weight. That is valid if it was analysed and had no calls. The missing metadata
check cannot distinguish that case from a chromosome never analysed.

**Agreed change:** Reuse the existing reference representation:
[`ContigFootprintEntry`](../../../src/shared/reference.rs#L15) records a contig's name and length,
and GC packages encode the sorted entries as `reference_contig_footprint_json`. Populate the
same representation from the BAM header for outliers. No new reference input is required.
Record the analysed chromosome names separately, using the existing `selected_chromosomes`
metadata. Contig names and lengths are the existing check, not a sequence-content checksum.

When the relevant information is supplied, fail on a reference mismatch or a requested
chromosome absent from the analysed list. An analysed chromosome with no weight rows retains
identity. Continue to reject malformed rows, invalid coordinates/weights, and malformed or
contradictory recognized metadata.

Differences in fragment length bounds, mapping quality, paired/unpaired mode, proper-pair
policy, or gap handling should produce clear warnings. Do not compare detection blacklists,
BAM files, or sample identifiers. Downstream users
may intentionally choose different filters. Do not require a population-expansion override,
a matching BAM path, or an additional unchecked-input option.

Custom tables produced by another outlier method are supported. Accept tables without
reference/chromosome/settings metadata while validating the required interval and weight
columns as usual. Check any supplied recognized information, and warn once about checks that
could not be performed. Use plain language, for example: "This table does not record the
reference, analysed chromosomes, or detection settings. These checks could not be performed."
Do not use "provenance" in CLI help or warnings.

### O-005: Zero-weight fragments are counted and weighting is invisible in statistics

**Confirmed zero-weight accounting defect. Normalization changes are explicitly out of scope.**

[Segment insertion](../../../src/commands/fcoverage/fcoverage.rs#L1863) reports success when a
fragment overlaps a tile even if its scalar weight is zero. All three GC/no-GC branches then
increment `counted_fragments` and the existing normalization sums. The code multiplies by zero
without skipping the fragment. This is a missing early exclusion, not a design limitation.

**Agreed change:** Resolve the outlier keep weight once the filtered fragment's original
`pos` to `reference_end` span is available. For a zero keep weight, record its exclusion and
skip before GC correction, length-normalization work, optional trimming, coverage insertion,
and counted-fragment accounting. Apply this to all three existing counting branches.

Add focused summaries for downweighted and excluded fragments. Keep the current tiling and
ownership design. Label the counting convention accurately and do not introduce a general
unique-fragment accounting refactor. Zero-weight fragments do not become counted contributions
or enter the existing normalization sums. Cover the all-zero case and cross-tile reporting.

Keep normalization formulas, ownership rules, and the treatment of partially downweighted
fragments in restore-mean statistics unchanged. Do not weight restore-mean sums, change GC
normalization, or broaden this patch into other correction channels. Skipping zero-weight
fragments naturally removes excluded fragments from existing accounting, which is the intended
consequence of early exclusion.

The user's working expectation is that these residual outliers will not materially change the
relevant global fragment length distribution. A possible investigation of weighted restore-mean
statistics is recorded under [user-owned later consideration](outliers_real_data_validation.md#user-owned-later-consideration).
Only the user may choose to revisit it. It is not part of this rectification plan and must not
be initiated automatically.

### O-006: Local-fit errors omit the genomic context

**Confirmed.** [Per-core fitting and diagnostics](../../../src/commands/outliers/striding.rs#L171)
propagate their numerical errors without chromosome/core/context information.

**Plan:** Add chromosome, core `start`/`end`, context `start`/`end`, eligible positions, and
fragment support to errors at the local fit and diagnostic boundaries. Preserve the underlying
cause and the deliberate rule that complete local contexts cannot silently use global fallback.
Add a focused internal error-context regression.

### O-007: An accepted tail probability can leave no second-fit support

**Confirmed, with a broader cause than shallow sequencing alone.**

[Configuration](../../../src/commands/outliers/config.rs#L44) accepts `(0, 0.5]`.
[The second fit](../../../src/commands/outliers/model.rs#L320) rejects `T1 = 1`.
For a Poisson mean 0.17, `P(X >= 1) = 1 - exp(-0.17)`, about 0.156, so alpha 0.2 gives this
failure. High zero inflation can also cause it, including with a larger Poisson-component mean.
The global fit fails after the first scan and already has a global-fit error prefix.

**Plan:** Report the effective alpha, fitted positive-position probability, and `T1`, explaining
that the probability leaves no retained positive observations and advising a smaller alpha.
Also contextualize failures where `T1 > 1` removes every observed positive value. Keep them as
explicit errors. Do not silently force `T1 = 2` or bypass the second stage. Config-only validation
cannot determine the sample-dependent condition before the coverage scan.

### O-008: Fragment-support lookup has a silent zero fallback

**Confirmed invariant-handling defect, not a demonstrated current index miss.**

[Histogram construction](../../../src/commands/outliers/outliers.rs#L392) uses
`.get(support_offset).copied().unwrap_or(0)`. The first pass always requests support and its
allocated range covers the core indices being processed, so the lookup should succeed.

**Plan:** Replace the fallback with an error containing chromosome, tile interval/index, core
index, requested offset, and available support length. Preserve a genuine zero support value.
Exercise nonaligned tiles at the integration boundary. Do not add an otherwise unnecessary
wrapper solely to make the impossible-index branch testable.

### O-009: Documentation incorrectly describes downstream support as absent

**Confirmed in the original verification, partially addressed in the workflow clarification.**
Stale statements were found in outliers config/run docs, the public guide, and the current
outliers spec. `fcoverage` already consumes the weights. The public guide and current outliers
spec have now been updated while documenting the intended use after existing blacklisting.
The stale statements in Rust config/run documentation remain to be fixed.

**Plan:** Update the remaining statements to name `fcoverage --outlier-weights` and describe
other consumers as planned. The guide now includes a usable weighting example and the
full-fragment, any-positive-overlap, minimum-over-regions rule. Preserve correct surrounding
explanations.

### O-010: Guide/spec descriptions omit or misdescribe relevant behavior

**Mostly confirmed. One wording criticism is stronger than necessary.**

[The guide](../../../website/docs/guides/handle_outliers_guide.md#L130) calls
`eligible_positions` data used to estimate the model. For fallback rows it instead describes
the chromosome context compared with an independently fitted global model. Its
`fragment_support` wording omits the unblacklisted-midpoint requirement.
[The fcoverage spec](../specs/fcoverage_spec.md#L16) omits outlier weighting.
The fcoverage help's phrase "retains a fraction" can mathematically include zero, so it is not
strictly false. It should explicitly mention complete exclusion and avoid implying that
weighting guarantees protection of every biological signal.

**Plan:** Correct fallback/context language, the midpoint definition, the weighting spec,
and zero-weight help. Document early zero-weight exclusion and the agreed reference/settings
checks when implemented. Keep normalization definitions unchanged. Update current specs only
to implemented behavior.

### O-011: The overlap resolver is implemented inside fcoverage

**Confirmed reuse concern, not a current lookup defect.**

[The closure](../../../src/commands/fcoverage/fcoverage.rs#L1043) owns the pointer, lookback,
and positive-overlap threshold. The loader/minimum rule are shared, but the module is gated
only by `cmd_fcoverage`. The current fraction threshold is valid because the query is the
original filtered `pos` to `reference_end` span, before trimming.

**Plan:** Before adding consumers, move this state and lookup into a small crate-private
resolver in `shared::outlier_weights`. Preserve original span lookup, minimum weight,
half-open boundaries, and maximum-fragment-length lookback. The existing BED overlap helper
already checks positive interval intersection before accepting an overlap, so a zero minimum
fraction in this specific BED-only lookup can express any positive overlap without a new
general overlap API. Alternatively use `Interval::intersects` in a small dedicated sweep.
Choose the simplest implementation and test the ordering invariant. Widen gates only for
commands actually wired in each patch. Public API visibility must not grow for test access.

### O-012: Two outlier meanings and different gap metadata keys

**Confirmed terminology inconsistency, not a CLI option collision.**

[GC matrix outlier options](../../../src/commands/gc_bias/config.rs#L293) have distinct names
from the proposed `--outlier-weights` input. They can coexist mechanically. Coverage-fragment
weights and correction-matrix clipping nevertheless need clearly separate help descriptions.
Outliers writes `ignore_inter_mate_gap`, whereas scaling writes `ignore_gap`.

**Plan:** Use "Coverage Outlier Weighting" for the new GC input's help and explain that the
existing outlier-method options act on the correction matrix. Use the common `ignore_gap`
metadata key when implementing the agreed reference/settings checks.

### O-013: Repeated log-factorial work

**Confirmed inefficiency, qualified complexity and impact.**

[Conditional means and CDFs](../../../src/commands/outliers/model.rs#L501) repeatedly call
`log_factorial` in two passes. Each call below 256 resums logarithms. This gives quadratic
work while the coverage range remains below 256. At or above 256 the helper uses a fixed-cost
Stirling approximation, so growth above that boundary is linear plus the bounded lower-range
cost. The bisection always performs 96 iterations after bracketing. No timing evidence was
provided, so its share of total runtime is unmeasured.

**Plan:** Make a small table for log factorials 0 through 255, initialized once. This is about
2 KiB of f64 values and directly removes the repeated summation. Keep the existing Stirling
branch for larger counts and preserve stable log-space calculations and threshold boundaries.
Do not gate this small improvement on a profiling project. Verify low counts, 255/256, shallow
means, and extreme tails. Bracket-stagnation stopping and broader numerical work are optional
later optimizations. No new distribution dependency is needed.

### O-014: Remaining behavioral coverage gaps

**Mostly confirmed. The missing-chromosome claim needs narrowing.**

Current outlier-weight integration tests exercise no-GC and GC-package products, 1 bp overlaps,
overlap in an ignored gap, lookup before trimming, and zero output. They do not exercise the
GC-tag product, out-of-start-order pairing lookback, outlier-weight application across tiles,
zero-weight counters, or declared-but-unanalysed chromosome handling in `fcoverage`.

[A loader test](../../../src/shared/outlier_weights_tests.rs#L26) already selects chr2 without
any chr2 rows and asserts that no weight entry exists. It verifies sparse identity behavior,
though not its use in `fcoverage` or the distinction from missing analysis information.
Outliers config tests cover incompatible gap/unpaired flags, but there are no detection
end-to-end tests for `--ignore-gap` or `--reads-are-fragments`.

**Plan:** Add the missing integration cases while retaining and extending the existing loader
test. The review's A/B lookback fixture is valid: A spans `[100, 300)`, with mates starting at
100 and 270. B spans `[240, 280)`, with mates starting at 240 and 260. Coordinate-sorted reads
complete B before A. With max fragment length 200, `[110, 120)` must still affect A. These
fragment lengths meet the default minimum inclusion filter. Add policy/error tests alongside
O-004 through O-007 changes. Include custom tables, missing/partial recorded settings, filter
warnings, reference/chromosome errors, merged-region preservation, and early zero-weight
exclusion in all GC branches. Derive expected weights, coverage, counters, and masses before
implementation and keep internal helpers private. Do not add weighted-normalization cases.

## Checks of claims listed as non-findings

- The minimum rule uses slice indices returned by BED overlap lookup, not the TSV line-number
  identifiers. The current indexing is correct.
- Pairing completes on the later encountered mate. Coordinate-sorted BAM input is enforced.
  A maximum-fragment-length lookback bounds the possible backward jump in original fragment
  starts. The proposed A/B fixture is useful regression coverage, not evidence of a current bug.
- The current `1 / (max_fragment_length + 1)` fraction admits 1 bp overlap for every accepted
  original fragment span. Trimming does not change the query used by the resolver.
- BED parsing skips the generated metadata comments. The listed receiving commands default to
  a blacklist minimum size of 1 bp. Masking semantics still depend on the receiving command.
- Raw coverage accumulation is in f64 and retained values below the guard are exact integers.
  The review's boundary explanation is off by one: 16,777,216 is exactly representable in f32,
  while 16,777,217 is not and can round to 16,777,216. Rejecting stored values at or above
  16,777,216 is therefore deliberately conservative and should remain. Correct explanatory
  wording, including the test comment calling that boundary value itself inexact.
- Histogram support is assigned by unblacklisted midpoint within the tile core. Partial-core
  histograms and support are added across tiles. Ordered tile collection and call merging
  preserve connected calls and their masses. A blacklisted base interrupts a call.
- The untruncated and right-truncated ZIP equations match their likelihoods. With retained
  positive observations and zero fraction `r < 1`, the stated truncated pi denominator is
  positive. The conditional means are increasing in identifiable cases. Degenerate retained
  supports still need explicit errors rather than a blanket guarantee that every input fits.
- Inclusive survival and predecessor checks in threshold search match `P(X >= threshold)`.
  Direct tail summation begins at or above lambda. No threshold off-by-one defect was found.
- The simple fixture has 880 zeros, 90 ones, and 30 values of 20. Its first positive mean is
  5.75, giving `T1 = 7`. Retention leaves 970 values and the boundary Poisson fit gives `T2 = 2`.
  Threshold-target weight is `2 / 20 = 0.1`. The boundary fixture has corresponding threshold
  weights of 0.1. The coverage-3 roundtrip has initial mean 3.5 and final underlying mean near
  2.695, consistent with thresholds 13 and 12. These were checked by reading and hand reasoning,
  not by executing tests.

## Implementation order

### 1. Small mechanical and CLI documentation fixes

Address O-006, O-007, O-008, O-009, and the unambiguous parts of O-010. Correct the f32
explanation. Add error and supported-mode regression cases from O-014. O-013's small
log-factorial table can be a separate patch, with explicit numerical-boundary verification.
These mechanical changes preserve fitting and correction results for valid inputs.

CLI help is the primary documentation throughout. Add the agreed descriptions when each
behavior is implemented, including useful explanations for the command and its arguments.

### 2. Reference/settings checks and shared weight lookup

Address O-004 with the agreed encoded reference, analysed chromosomes, warnings for setting
differences, and support for custom tables. Move the small O-011 resolver without changing
lookup behavior or public API visibility. Apply O-012's plain terminology and common gap key.
Cover lookback, positive overlap, touching boundaries, multiple intervals, original-span lookup,
and tiles. Do not undertake a broad overlaps or command-plumbing refactor.

### 3. Early zero-weight exclusion and focused reporting

Address O-005 by skipping zero-weight fragments before counting and correction work. Add
downweighted/excluded summaries using the existing accounting conventions. Preserve normalization
formulas and ownership rules. Do not implement weighted restore-mean statistics. Broader changes
to counters, normalization, or other correction channels require separate user consent.

### 4. Working ZINB model, active-component target, and preservation policy

Implemented for the first BAM-inspection milestone, with the defaults and supported backward
fallback stated at the top. The current spec and CLI define implemented behavior. Historical
10% formulations in the original review discussion below are superseded. Biological validation
and any later changes to the preservation evidence remain separate work.

### 5. Later validation and future rollout

The [real-data validation list](outliers_real_data_validation.md) records later investigations
across 0.1x through 200x, including 0.5x, known amplifications, and tiny fragment overlaps.
It does not authorize changes to normalization.

Keep O-002's GC/scaling producer integration and wider consumer rollout outside this
rectification scope. Measure calibration effects after ordinary blacklisting before deciding
that work. Use the shared resolver when rollout is requested. Further performance work may
use profiling. None of this authorizes an unrequested broader refactor.

## Implementation invariants

- Paired fragment spans remain directional from `forward.pos` to `reverse.reference_end`.
  Coverage segments, gap policy, and weight lookup span retain their separate meanings.
- The local-mean target remains the active-component mean. The current ZIP unconditional mean
  retains its definition. ZINB integration must use its own fitted parameters and likelihood,
  with coherent initial/refitted means, dispersion, zero inflation, and survival probabilities.
- Original raw histograms and exact call identities remain available for auditing.
- Scores use unflanked exact calls and eligible unblacklisted bases, including zeros. The
  candidate short-window rule applies a denominator floor and records actual support.
  Earlier evaluated windows contribute to later overlapping calls' maxima.
- Flank and merge all candidates before filtering. Ambiguity in any member applies to the whole
  merged region. Every actionable output follows the same region decision.
- No subtraction of protected bases, pooled merged-region weights, or fragment-level protection
  overrides are introduced. Selected exact weights retain the minimum-over-overlaps rule.
- Zero outlier keep weight means early exclusion from counting. Preserve normalization formulas
  and ownership rules, including the current treatment of partial weights.
- Weight lookup remains consistent across tiles and GC branches. Missing overlap is identity.
  Missing analysis information produces a warning, while known reference/chromosome mismatches
  produce errors.
- Reporting uses accurately named existing accounting conventions. Do not describe tile
  insertions as unique fragments or add a general accounting refactor.
- Weighted restore-mean statistics belong only to a user-owned later consideration. They must
  not be implemented or investigated unless the user explicitly requests that work.
- Current specs describe implemented behavior. This future plan does not make its proposed
  behavior current.

## Verification completed and required

Read-only verification completed:

- `cargo fmt --all -- --check`: passed.
- `cargo check --offline --tests --features "testing cmd_outliers"`: passed.
- `cargo check --offline --no-default-features --features "cli cmd_outliers"`: passed with
  six dead-code warnings in shared I/O and heatmap helpers, matching the source review.
- No tests were run. No numerical replicas or scripts were used to derive test expectations.

For implementation, run `cargo check` after code changes and
`cargo check --tests --features "testing cmd_outliers"` after test changes. Include non-default
consumer features when wiring those commands. Run `cargo check --all-features` after a major
refactor. Test execution remains with the user. Larger refactors and unresolved scientific
policy choices must be discussed before their implementation.
