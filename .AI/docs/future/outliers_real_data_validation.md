# Outliers real-data validation list

Status: Pending. These are investigations, not claims that the current correction is validated.

Related work: [remaining decisions and code work](outliers_remaining_work.md),
[current spec](../specs/outliers.md), and the
[short biological discussion report](outliers_biological_release_criteria_report.md).

The intended workflow starts with an existing blacklist and detects residual sample-specific
extremes. Current data are about 0.5x. Validate the required range of 0.1x through 200x, recording
both nominal sequencing depth and counted positional fragment coverage. Support over this range
has not yet been established by these investigations. These comparisons are later work, not
prerequisites for the active-component target or initial ZINB default. The first BAM-inspection
version uses a 5% cutoff, 100 kb physical windows, and a 50 kb support/denominator floor, with the
supported backward fallback. These initial choices still need biological assessment.
This list does not authorize normalization changes.

## Decisions for the owner and supervisor

- [ ] The owner and supervisor choose the biological event sizes, acceptable losses, downstream
  measurements, sample types and depth range that define a useful result. The discussion report
  proposes questions without selecting tolerances. The investigation list below remains open.

## Depth range and baseline estimation

- [ ] Evaluate 0.1x, 0.2x, 0.5x, 1x, 5x, 30x, 100x, and 200x, using compatible real samples or
  downsampling where possible. Do not extrapolate a shallow sample upward as evidence for 200x.
- [ ] Check Poisson sampling zeros separately from fitted structural-zero inflation. Record
  numerical fit failures, support-driven context expansion, thresholds, and correction severity.
- [ ] Verify that the automatic probability denominator and all fitted histograms exclude
  blacklisted bases while retaining ordinary eligible positions with zero coverage.
- [ ] Compare the current two-stage fit with a single fit on the same contexts. Show first-stage
  and second-stage calls, excluded position counts and mass, fitted baselines, and changed calls.
- [ ] Examine whether excluded tails are isolated pileups or parts of coherent biological
  elevations. Decide whether exclusion needs additional spatial evidence or another baseline model.
- [ ] Estimate ordinary expected-rate variation within actual fitting contexts after masking
  and filtering. Separate sampling noise, structural zeros, technical excess, and copy-number
  changes. The simulated CV values are not measured library estimates. Assess spatial scale
  and tail shape rather than substituting observed count CV for underlying rate CV.
- [ ] Compare fitting minima of 5, 10 and 20 million eligible bases at the same 500 kb stride,
  masks and tail probability. Reuse saved core histograms for fitting comparisons. Include the control's chr1
  region near 16.6 Mb and chr21's first usable cores at 8.5-10 Mb, whose supplied models have
  T2=464 and contexts containing only 719,053 eligible bases.
  Context expansion now requires the configured eligible-base minimum and midpoint-assigned
  fragment support. Reassess these previously undersupported regions. Track eligible support,
  actual physical context spans, fitted means, dispersion, optimizer
  agreement and T1/T2. Check whether a localized excess raises cutoffs
  across neighboring cores, and whether larger contexts dilute that effect or spread it farther.
  Include broad copy-number changes before choosing a default. Smoother cutoffs alone do not
  establish better detection or biological retention.
- [ ] Compare ZIP and the working ZINB model on identical contexts, first examining detection
  and then the final preservation/correction policy. Record fitted dispersion, pi, both means,
  T1/T2, and observed-versus-fitted tails. Assess weakly determined dispersion at low depth.
- [ ] Assess weaker and moderate technical elevations across heights, widths, and depths.
  Check whether excess inflates the initial dispersion enough to escape T1 exclusion and
  final calling. Strong-spike results alone do not establish sensitivity over this range.

## Full-fragment correction after any positive overlap

- [ ] Inspect the current rule on real loci: even 1 bp overlap of a fragment's complete `pos`
  to `reference_end` span applies the full regional weight. Multiple overlaps use the minimum.
- [ ] Show fragments fully inside a call, crossing either boundary, barely touching it, and
  spanning two regions with different weights. Include overlaps confined to the inter-mate
  gap or other uncounted reference segments.
- [ ] Build aligned locus panels with raw coverage, local fitted baseline and `T2`, exact called
  intervals and weights, fragment spans colored by resolved weight, corrected coverage, and
  coverage removed. Include the fragment endpoints so the effect outside calls is explainable.
  The current coverage plots do not show individual fragments, so those panels remain later work.
- [ ] Plot overlap bases and overlap fraction against assigned weight. Summarize how much removed
  coverage lies inside calls and in surrounding positions, grouped by distance from boundaries
  and fragment length.
  Compare the expected ordinary-background contribution with total observed corrected coverage.
  The linear recovery example assumes stationary starts, a fixed contiguous counted fragment span,
  an isolated constant-weight call, and no other corrections. Variable lengths, alignment gaps,
  artifact fragments, and nearby calls can change the observed shape.
- [ ] Compare full correction, overlap-fraction-diluted whole-fragment correction, and positional
  correction on the same detected regions. Define each candidate mathematically before comparing
  it. Assess changes to regional mass and neighboring signal, rather than assuming all choices
  can satisfy the same target with the current regional weights.
- [ ] Check whether short calls downweight ordinary neighboring fragments disproportionately,
  and whether affected fragments have independent evidence of mismapping or residual duplication.
  Record when that evidence is unavailable. Extreme coverage alone cannot identify bad fragments.
- [ ] Produce annotated, exportable plots suitable for explaining the chosen rule in the guide.
  Show examples where the rule is useful and examples where it removes unwanted neighboring
  signal. Do not select only examples that support the current design.

## Replacement target and reference gaps

- [ ] Compare unconditional ZIP mean `(1 - pi) * lambda`, Poisson-component mean `lambda`,
  and positive conditional mean `lambda / (1 - exp(-lambda))`. Label their meanings explicitly.
  For ZINB, compare the corresponding active mean `mu`, unconditional mean `(1 - pi) * mu`,
  and positive conditional mean `mu / (1 - P(NB = 0))`. The active mean remains the working target.
- [ ] Inspect loci with differing structural-zero fractions at comparable active coverage,
  including reference gaps, centromeres, and deletion candidates. Compare masked and unmasked
  context definitions and annotate reference/mappability evidence where available.
- [ ] Include shallow data, where lambda and the positive conditional mean differ substantially.
- [ ] Keep the target's biological definition consistent across depths. Assess uncertainty in
  lambda and pi at shallow depth rather than choosing a different mean solely because Poisson
  sampling zeros are common. Check that pi near zero makes both candidate baseline targets agree.
- [ ] Measure corrected coverage against surrounding comparable loci, and quantify the dependence
  of correction severity on fitted pi. Assess the agreed lambda target against the previous
  unconditional target. The default target choice is already made theoretically.

## Biological amplifications and spatial calling

- [ ] Include independently annotated focal amplifications and broad copy-number changes,
  alongside likely technical pileups and healthy/control samples.
- [ ] Plot the event span, local fitting contexts, `T1`, `T2`, calls, weights, and correction.
  Look for fragmented calls within a broader elevated event as well as long connected calls.
  Prioritize sparse calls within moderate genuine gains and technical spikes within gains. Record
  biological signal lost separately from technical signal removed, using independent evidence
  or known simulated fragment origins where available. Compare events near 2.5 kb and 5 kb
  preservation boundaries with wider events. The current called-base percentage is not the event span.
- [ ] Compare tumor fractions, event widths, coverage depths, and event positions relative to
  model cores. Include strong and moderate copy gains. State the mixture assumptions connecting
  copy number to observed relative coverage. Evaluate whether contamination of the fitting context
  hides larger events and whether exclusion makes genuine biological events more severely corrected.
- [ ] Assess the final preservation safeguard on all actionable outputs, including exact/flanked
  blacklists. Verify that every member of a preserved merged region is omitted consistently.
- [ ] Evaluate the current 100 kb call-start windows and strict 5% cutoff on final T2 calls.
  Use the maximum only among adequately supported overlapping forward windows when any exist.
  Otherwise try the supported backward fallback before retaining the floored forward maximum.
  Include eligible sampling zeros and record actual support. Calculate these metrics before
  flanking. Compare fragmented events with isolated spikes and repeat-associated technical elevations.
- [ ] Plot the fraction tracks with calls, raw coverage, baselines, and known copy-number spans.
  Assess percentage thresholds and minimum eligible support across depths, including 0.5x.
  Compare alternatives against the initial 5% cutoff and 50 kb minimum support.
- [ ] Evaluate withholding correction from ambiguous final calls while keeping the existing
  two-stage baseline fit. Verify that preservation leaves histograms, fit parameters, T1, and
  T2 unchanged, while retaining all final detection evidence.
  Flank all candidates and merge overlapping/touching intervals before filtering. Verify that
  ambiguity in any member preserves the entire merged region, including lower-percentage calls.
  Apply that decision consistently to weights and actionable blacklist outputs. Do not subtract
  protected bases from flanks or introduce fragment-level protection overrides. Keep the cause
  of a flag unresolved until independent biological or technical evidence is available.
- [ ] Include calls appearing only at T2 when calculating window percentages. Check that the
  preservation decisions use the complete final call set before any output filtering.
- [ ] Assess the chosen preservation cutoff and minimum eligible-position setting across
  depths and window sizes. Measure the actual fragment-weighted signal removed from known
  biological events, including sparse calls at shallow depth and fragments crossing boundaries.
- [ ] Include cases near the current preservation boundary, where a small change in called-base
  percentage can change whether a biological gain is corrected. Measure signal loss using the
  working model and current settings. Also compare sustained coverage evidence if selected.
- [ ] Verify that user changes to window size, percentage cutoff, minimum support, and ambiguity
  action have the effects explained in CLI help. Include threshold equality and chromosome-end
  windows in the hand-derived behavioral cases.
- [ ] Evaluate automatic preservation of ambiguous calls with a full audit table and a correction
  file usable directly in a pipeline. Report how much correction is withheld and how often true
  technical excess is left in the output. Manual editing must not be required for this default.
- [ ] Consider larger screening windows or full-coverage event separation only if the 100 kb
  comparisons expose a failure they can address. A separate correction baseline inside an
  amplification is additional fitting work, not part of the current window-percentage proposal.
- [ ] Investigate narrow technical spikes inside an amplified region and fragments crossing
  protected boundaries. Distinguish preservation of the regional baseline from preservation of
  every fragment contribution at neighboring corrected calls. Keep the existing fragment lookup
  rule and evaluate the effects of the configured flank on merged-region membership.
- [ ] Measure uncertainty and spatial resolution at each depth. A fixed window and percentage
  definition do not imply that evidence is equally informative at 0.1x and 200x.
- [ ] Decide whether independent protected-region annotations or another explicit policy are
  needed. Do not describe corrected candidate loci as proven technical artifacts.

## Correction order and recorded settings

Integration across applicable commands is required before release, as recorded in the
[integration plan](outliers_remaining_work.md#required-integration-before-release). GC/scaling
estimation must use outlier weights even when the final measurement uses blacklisting.
The comparisons below assess effects and limitations. They do not gate this design decision.

- [ ] Compare matched raw, outlier-only, GC-only, and combined coverage at selected loci with the
  same fragment filters and masks. Check GC content, mappability, promoter/CpG annotations where
  useful, and independent copy-gain evidence against matched background. Multiplying weights alone
  does not demonstrate harmful double correction. Distinguish consumer effects from producer
  estimation and the separately deferred restore-mean question.
- [ ] Compare scaling and GC corrections estimated before and after applying residual outlier
  weights, with the same existing blacklist. Show surrounding bins as well as the pileup itself.
  Quantify residual excess mass relative to each smoothing context before judging its importance.
- [ ] Inspect length/GC strata to identify when nuisance fragments distort GC estimation and
  when their contribution has little effect. Check supported and unsupported package bins.
- [ ] Compare the full raw -> weighted GC -> weighted scaling -> feature-counting workflow with
  workflows that mix corrections built with different settings. Record factor identities,
  filters, reference, and sample selection.

Implementation tasks belong in the integration plan. Keep zero-support windows explicit when
assessing weighted GC estimation and scaling.

## Model fit, reporting, and performance

- [ ] Audit local/global fallback use, context support, positive-count variance ratios, and
  observed-versus-fitted tails across assays and depths. Separate localized excess from broad
  positive-count overdispersion when assessing the working ZINB model. Check the assumed extreme
  tail and parameter uncertainty separately from favorable finite-run signal retention.
- [ ] Show unchanged, downweighted, excluded, and effective kept fragment support. Verify early
  zero-weight exclusion while preserving the existing normalization formulas and ownership rules.
- [ ] Profile BAM scans, context aggregation, numerical fitting, and output generation at shallow
  and deeper coverage for broader performance work.

## User-owned later consideration

- [ ] Consider whether outlier keep weights should affect restore-mean normalization statistics.
  This is deferred to the user, not scheduled agent work. Do not investigate or implement it unless
  the user explicitly requests it. The working expectation is that residual outliers do not
  materially change the relevant global fragment length distribution. Normalization changes
  require a separate explicit request.

Record the sample/locus selection, settings, annotations, plots, observations, and unresolved
counterexamples for each investigation. Keep this list pending until those artifacts exist.
