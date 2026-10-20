# Archived annotated outliers review, 2026-10-19

This preserves the annotated review and proposals at revision `9e9bdac5`, including Claude's
inline responses. It is historical evidence, not an active checklist. The owner's later decisions
and completed implementation are recorded in [the current checklist](outliers_review_rectification_plan.md).
Later direction supersedes the first-offending-region proposal: report all ambiguity errors after
the second BAM pass. GC-counter work and naming remain deferred.

# Outliers review verification record

Date: 2026-10-19.
Status: Dated verification record, including Claude's annotations. Implementation and scientific validation remain pending.
Baseline: `outliers2` at `9e9bdac5dedb16103817c1b88fdc986004fcae0c`.
Input: [Claude's review](../../../major_reviews/outliers_new.md).

## Owner response and current work tracking

The owner's 2026-10-19 responses supersede the original implementation proposals below.
The [rectification checklist](outliers_review_rectification_plan.md#owner-direction-2026-10-19)
now holds the consolidated eight decisions and remaining implementation work.
The [validation list](outliers_real_data_validation.md#short-report-for-the-owner-and-supervisor)
holds the requested short, accessible report on biological risks and possible release criteria.
This file preserves the derivations, evidence, and review exchange rather than serving as a third
active checklist. Condense it after the remaining ridge-rule and GC-counter decisions are settled.

In particular, the older proposals to publish a complete failure audit and to attempt an unused
global fit for diagnostics are superseded. Failures must publish no new files under any normal
success-output filename. Error diagnostics should be human-readable and focused on the actual
failure. The preferred global-fit direction is to skip it when no core needs fallback, provided
the changes remain contained. Feature/default decisions and renaming remain with the owner.

## Original verification and annotated proposals

This is a future work plan, not a change to the implemented specification. No production code or tests
were changed during this verification. The supplied review was left untouched. Read the source again
when implementing, since the locations below describe this revision.

The consequential findings are loss of failure diagnostics, an unnecessary global-fit dependency,
unidentified truncated fits, missing positive-dispersion command coverage, and biological signal loss.
Several recommendations need correction before implementation. In particular, counting occupied
positive coverage levels is not an identifiability criterion. ZIP also has an unidentified case.
Neither an explicit alpha nor preservation of a called region guarantees preservation of neighboring
fragment contributions.

> **Claude note (2026-10-19):** Reviewed after Codex wrote this plan. The derivations below were
> re-checked by hand and hold. The corrections to the original review are accepted. Inline notes mark
> the few places where this plan needs a concrete fix or a decision. Answers to the questions, the
> consolidated owner decisions, and the remaining issues are in
> [Claude review of this plan](#claude-review-of-this-plan) at the end.

Related documents:

- [Current implementation invariants](../specs/outliers.md).
- [Existing remaining-work checklist](outliers_review_rectification_plan.md).
- [Scientific validation list](outliers_real_data_validation.md).

Main source locations: [runner and calling](../../../src/commands/outliers/outliers.rs),
[ZIP fitting](../../../src/commands/outliers/model.rs),
[ZINB fitting](../../../src/commands/outliers/zinb.rs),
[two-stage fitting and tails](../../../src/commands/outliers/distribution.rs),
[preservation](../../../src/commands/outliers/preservation.rs),
[detector command tests](../../../tests/test_outliers_command.rs), and
[fcoverage application](../../../src/commands/fcoverage/fcoverage.rs).

This plan adds detail to those pending items. It does not authorize additional consumers, GC/scaling
producer integration, changes to restore-mean statistics, or a broad refactor.

## Evidence and verification limits

Re-read the detector, fitting, context selection, preservation, inspection outputs, plotting, CLI help,
generated README source, website guide, command tests, relevant private tests, weight loader, and
the downstream application paths. Traced their relevant shared interval, tile, read-filter, fragment
assembly, coverage, overlap, and output-publication helpers. Also read the locally installed
argmin 0.11 Nelder-Mead and statrs 0.18 beta implementation.

The branch HEAD matches the revision reviewed by Claude. The initial working-tree change was the
staged review file. Evidence below distinguishes code behavior and algebra from data-dependent risks.
No BAM analysis, simulation, benchmark, or project test was run. Test compilation does not establish
test success.

| Command | Observed result |
| --- | --- |
| `cargo check --offline --locked --tests --features "testing cmd_outliers"` | Passed with no warnings reported. |
| `cargo check --offline --locked --no-default-features --features cmd_outliers` | Passed with the six reported dead-code warnings. |
| `cargo check --offline --locked --no-default-features --features cmd_fcoverage` | Passed with no warnings reported. |
| `cargo clippy --offline --locked --lib --features cmd_outliers` | Failed with a denied comparison lint in `src/shared/plotters/heatmap.rs:476`. Reported 97 warnings, including the three cited model lints. That shared file has no changes in `main...HEAD`. |
| `git check-ignore -v .AI/docs/future/outliers_fragment_simulations/zinb_comparison_report.md` | The report is ignored by `.gitignore:56` and is not tracked. It nevertheless exists and is readable in this workspace. |

The clippy invocation was repeated only to capture its truncated diagnostics. The error is the
`RIGHT_PANEL_GAP > 0` comparison with a `u32` constant equal to zero. Do not describe the clippy run as
passing after fixing three outlier warnings, and do not bundle unrelated lint cleanup into this work.

## Findings register

Each numbered finding in review sections 1-6 appears below. "Confirmed" means the behavior or missing
coverage follows from the source. It does not mean every proposed remedy or severity is accepted.

### Correctness and behavior

| Review ID | Verdict and evidence | Handling |
| --- | --- | --- |
| 1.1 | Confirmed usability gap. `run_outliers` writes the histogram at `outliers.rs:297`, but publishes at `:376`. Any intervening error drops the staging guard. `evaluate_regions` returns at the first flagged merged region in error mode, before a complete audit exists. This obeys the current no-results-on-error contract. | Change the failure-diagnostics contract deliberately. Preserve the completed histogram and contextual failure report. For ambiguity errors, evaluate all decisions and publish diagnostic evidence while withholding weights/BEDs. See the failure-output design below. |
| 1.2 | Confirmed dependency. `outliers.rs:306` always fits the pooled histogram. `fit_core_models` uses it only for incomplete chromosome contexts, and also obtains the model family from it. Whether extreme pileups actually drive dispersion to the upper bound is unproven. Cross-chromosome rate heterogeneity is another possible source of poor pooled fit. | Determine fallback need before fitting. Keep failures fatal when a fallback model is required. Recommend an explicit warning and unavailable global-fit status when only the diagnostic fit fails and all local fits succeed. Preserve observed global counts and local diagnostics. |
| 1.3 | Confirmed statistical issue, with material corrections. ZINB has exact likelihood ridges for retained supports below 2 or 3 in the stated interior cases. ZIP is also unidentified when T1 is 2. Conversely, an occupied-level count would reject valid boundary fits. `distribution.rs:143` and `model.rs:226` only guard T1 greater than 1. | Prioritize a model-specific identifiability rule and regression cases. Do not implement the proposed three-level/two-level rule. See the derivations below. |
| 1.4 | Confirmed optimizer stopping rule. argmin stops on vertex-cost standard deviation, without a parameter-spread condition. `zinb.rs:206` selects improvements above a fixed average-cost tolerance. The frequency and practical threshold spread remain unmeasured. | Record converged start outcomes and distinguish near-optimal alternatives from inferior local optima. Assess mean, dispersion, pi, and T2 sensitivity. Agreement across three starts is a diagnostic, not proof of identification. |
| 1.5 | Confirmed latent error path. argmin initialization calls `cost(...).unwrap()`. In-domain errors from `ProfiledHistogram::evaluate` can therefore panic. Later trial errors propagate through `Executor::run`. No realistic triggering histogram was demonstrated. | Preflight initial vertices with contextual errors. Decide explicitly which numerical trial failures are recoverable. Do not indiscriminately convert every error to infinity. Keep strict, finite re-evaluation of any accepted parameters. |
| 1.6 | Reporting omission confirmed. Both thresholds already exist in the table, but no relation flag is provided. T2 greater than T1 is mathematically reachable, even for an identified ZIP fit. See the example below. | Explain that refitting need not lower T2. Add a count or flag only if it makes inspection clearer. Never clamp T2 to T1. This is not intrinsically an incorrect result. |
| 1.7 | Confirmed selection dependence at `outliers.rs:277`. It is intentional selected-base accounting. "T2 drops in every core" is too strong because thresholds are discrete and the retained second-fit data also change. | Document split-run consequences. Prefer detecting once over the intended chromosome set and consuming subsets. A shared manual alpha fixes one source of differences, but cannot reproduce a different pooled fallback model. |
| 1.8 | Confirmed counter discrepancy at `fcoverage.rs:1103` versus `:1159` and `:1216`. File-based GC rejects full fragment spans outside the tile fetch before lookup. Other branches count those zero-weight encounters. | Use the same encounter eligibility in all branches, with a boundary fixture. Preserve zero-weight exclusion before GC/normalization. This concerns reporting, not evidence of changed counted coverage. |
| 1.9 | Confirmed plotting limitations. `write_preserved_plots` and its chunk loop are sequential. Only 2,000 context bins are retained and reused for the close-up. Only group-member calls are shaded. | First improve close-up resolution and show all relevant candidates with their actual decisions. Report group count and requested BAM span/progress. Profile before adding concurrency. Parallel BAM readers can increase I/O contention. Any plot cap must report omissions explicitly. |
| 1.10 | Confirmed misleading labels at `plotting.rs:480` and unlabeled chromosome separators. The unconditional mean is diagnostic. Even the active mean is not the correction target when `--target threshold` or `zero` is selected. | Use fitted-mean labels and state the selected target separately. Add chromosome labels or a genomic coordinate mapping. |
| 1.11 | Metadata contents confirmed at `outputs.rs:112-215`. Their presence is not itself a correctness failure. `detector` also encodes the two-stage procedure, so it is not wholly redundant with distribution family. | Document the CLI mapping, make automatic-only fields explicitly inactive in manual mode, and remove obsolete hysteresis metadata if unused. If simplifying `detector`, retain the algorithm identity, for example `two_stage`. Coordinate any metadata-key changes with consumers. |
| 1.12 | Confirmed cosmetic parser/validator difference. `parse_positive_f64` accepts 0.5 for the plot multiplier, while `validate` rejects it. Percentage validation is also in `validate`. `run_outliers` validates before opening the BAM. | Low priority. Align CLI parsing with validation when touching these arguments, while retaining validation for programmatic configs. |

### Biological and methodological findings

| Review ID | Verdict and evidence | Handling |
| --- | --- | --- |
| B1 | The loss of neighboring background contributions follows directly from whole-fragment minimum weighting. The linear formula requires stationary starts, a fixed fragment length, contiguous counted spans, an isolated constant-weight call, and no other corrections. A literal V-shaped trough beside every real call is not established. | Treat surrounding-signal loss as a release-level scientific concern. Document the conditional calculation and inspect real loci with fragment endpoints. Do not promise a universal 5% observed coverage level or a universal shape. |
| B2 | Reachable failure mechanism. Neither calling nor density preservation recognizes a subthreshold broad gain whose occasional excursions form sparse calls. The particular 30x/3x/0.06 example and the claim that this is the most likely harmful mode are hypotheses. The context itself can also absorb a broad gain. | Prioritize independently supported gains, sparse excursions, and spikes within gains. Measure biological mass loss and downstream feature distortion. Do not select a new preservation algorithm solely from the approximation. |
| B3 | Raw detection followed by multiplicative GC weighting is confirmed in `build_tile_coverage` and `fcoverage::process_tile`. Whether this causes harmful excess correction at GC-rich elements is unmeasured. Multiplication alone is not proof of double correction. | Compare matched raw, outlier-only, GC-only, and combined tracks. Record GC/mappability/CpG annotations and filters. Keep producer integration and restore-mean investigations outside this plan. |
| B4 | Confirmed default upper bound for isolated exact calls under `preserve`. Above 5,000 called bases, the selected window exceeds 5% when its denominator cannot exceed 100 kb. Equality is allowed. Ends/masks can lower the effective limit, and nearby calls can preserve shorter calls. | Explain this conditional bound, not a general 5 kb detector limit. `correct`, a 100% cutoff, or different window settings change it. Include numerical boundary cases. |
| B5 | Confirmed at `src/shared/read.rs:5-29`. Both paired and unpaired paths reject duplicate flags, secondary/supplementary alignments, QC failures, unmapped reads, and low MAPQ. Paired reads also have mate/orientation/reference requirements. | Document the accepted fragment population. Say duplicate-flagged alignments are excluded. Avoid claiming all PCR artifacts are absent from every marked BAM. Residual unflagged pileups remain detectable. |

### Documentation findings

The help, `readme.rs`, and website guide were checked separately. Some requested information already
exists in one of them. Update each document for its audience without copying the entire spec into help.

| Review ID | Verdict | Handling |
| --- | --- | --- |
| D1 | Confirmed lack of troubleshooting in help/guide. ZIP truncated errors still omit probability advice. | Give stage-specific guidance, described below. A smaller alpha cannot repair an initial-fit failure, and a larger local context cannot repair a failed pooled fit. |
| D2 | Confirmed omission of the two detection passes in help/guide. The optional local reread is already documented. | State two tiled BAM passes plus optional context rereads. Do not promise exactly three genome-wide passes or invent runtime estimates. |
| D3 | Confirmed missing read-filter explanation. | Add the exclusions in B5, including MAPQ and the paired-read requirements where appropriate. |
| D4 | Partly confirmed, partly incorrect. Short contigs use the pooled model when span/support is insufficient. With `--chromosomes all`, the pool includes all selected contigs, not just autosomes. A chromosome being almost entirely called is data-dependent. If 16,569 eligible bases were all called, the default floor would give 33.138%. | Explain mixed-reference populations and fallback source precisely. Recommend a deliberate selection of comparable contigs. Blacklisting a contig removes it from eligible histograms but does not turn the remaining pool into a guaranteed autosomal model. |
| D5 | Confirmed eligibility rule. The detector never reads reference sequence. Unblacklisted N-gap coordinates remain eligible and contribute zeros when uncovered. They are not explicitly forced to zero. | Recommend assembly-matched gap masks when those coordinates should be excluded. Explain effects on histograms, alpha, and preservation scores. Do not remove ordinary sampling zeros. |
| D6 | Confirmed absent numeric bound. | Document B4 beside preservation settings with strict equality and boundary/mask qualifications. |
| D7 | Confirmed missing split-run warning. | Document alpha and fallback differences together. Demonstrate the preferred detect-once/consume-subset workflow. |
| D8 | Confirmed avoidable vagueness for the explicit `fcoverage` example. Positional outputs mask blacklist bases, rather than excluding whole overlapping fragments. BED `#` comments are accepted. | State that exact masking leaves coverage outside the call. Matching maximum-fragment-length flanks can mask contributions from fragments spanning the call, at the cost of masking neighboring signal too. Do not claim this removes their contributions from normalization accounting. |
| D9 | Confirmed missing quantitative shoulder explanation. | Add the qualified B1 derivation and an exportable real example showing raw, corrected, removed coverage, and endpoints. |
| D10 | Confirmed stale optional PNGs on prefix reuse. `FinalOutputFiles` moves registered new outputs and performs no old-artifact cleanup. | Document current-table authority and fresh prefixes. If adding cleanup, restrict it to known artifacts of the same successful output set. Do not delete files by a broad wildcard or remove prior valid results on a failed rerun. |
| D11 | Confirmed missing metadata-to-CLI mapping in the output README. | Add a compact mapping after CLI names are settled. Preserve the consumer's recorded-filter checks. |
| D12 | Confirmed missing core-index explanation in the README. The complete-histogram panel has no legend, unlike the detailed panel. | Add a legend or a clear shared key, explain the core-index axis, and label chromosome spans. Coordinate with optional/unavailable global-fit reporting. |

### Naming and CLI surface

These are design choices, not numerical defects. Proposed names below are recommendations for a
single pre-release cleanup, not accepted API decisions.

| Review ID | Verdict | Handling |
| --- | --- | --- |
| N1 | Confirmed `protection` versus `preservation` terminology. | Recommend `--preservation-window-size` and `--min-preservation-positions` if renaming. Keep consistent help, config serialization, metadata, and docs. |
| N2 | Confirmed `--bin-size` is a minimum fitting-context span, unlike output bins elsewhere. | Recommend `--min-context-span`. Keep `--stride` unless there is a concrete reason to rename it too. |
| N3 | Confirmed model-heading placement of correction target and output/preservation flank controls. | Move target under correction and flank under preservation/output. Explain that the flank also propagates preservation decisions. |
| N4 | Confirmed separate outlier concepts in `src/commands/gc_bias/config.rs`. The future `gc-bias --outlier-weights` combination is outside the current implementation. | Record the naming distinction for producer integration. Do not rename the existing GC estimator controls as incidental work. Explain histogram/estimator outliers versus sample-specific fragment keep weights. |
| N5 | Confirmed mix of ambiguous, flagged, density, and percentage vocabulary. | Prefer user-facing "preservation" and "called-base percentage". Consider `--preservation-action` and `--max-called-percentage` together. Keep audit reasons explicit. Finalize names with the owner before changing them. |

No compatibility aliases are needed by default. Do not add them for speculative backwards
compatibility. Ask only if compatibility becomes an actual requirement.

### Tests and assertion quality

| Review ID | Verdict and existing coverage | Handling |
| --- | --- | --- |
| T1 | Confirmed wrong column in `local_mean_and_zero_targets_write_their_regional_mass_weights` at lines 1014/1023. This fixture has pi zero, so merely switching columns does not strengthen discrimination. The roundtrip test already distinguishes active mean near 2.82 from unconditional mean near 2.70. | Correct the column and wording. Use header-based lookup for semantic assertions. Keep a fixture with nonzero fitted pi that explicitly rejects the unconditional target. Avoid a duplicate expensive roundtrip solely for this fix. |
| T2 | Confirmed missing complete two-stage ZINB test and command run retaining positive dispersion. The existing `truncated_geometric_histogram_recovers_its_underlying_component` does test a positive-dispersion truncated fit. Guard failures lack targeted coverage. | Add a hand-derived full-path fixture and private numerical failure/identifiability cases. Preserve the existing valid Poisson-boundary command case. |
| T3 | Confirmed missing partially split core in the command fixtures. Existing calls cross core and tile boundaries, but the cores themselves align with tiles. `build_tiles` rounds to stride multiples when at least ten fit. | Add a 600 kb core over 1 Mb tiles and assert partial histograms, midpoint ownership, and output invariance. Concrete counts are below. |
| T4 | Confirmed no command error/correct ambiguity coverage. The private flank-propagation test covers both policy branches. | Add publication/audit assertions at the command boundary after defining the new diagnostics contract. Test an entire flagged merged group, not only the policy enum. |
| T5 | Confirmed missing command-level alpha/multiplier assertion. A private metadata test writes a supplied probability and cannot verify the automatic calculation. | Assert eligible denominator, recorded alpha, multiplier, and equivalent manual settings. Include masked bases and eligible zeros. |
| T6 | Confirmed no successful detector run using `--reads-are-fragments`. Config tests do reject incompatible paired-only settings. | Add a small unpaired command fixture. Test wiring and persisted results, without duplicating all shared gap/segment tests. |
| T7 | Confirmed plot chunk accumulation is only simulated through bins. The command plotting fixture is a 1 kb chromosome under a 1 Mb tile. | Cover real chunk construction/fetch halos across a plot boundary and compare numerical bins. Do not assert only a PNG signature. |
| T8 | Confirmed no command-level preserved audit action/reason assertion. Existing command coverage verifies the preserved table and omission from weights/exact BED. | Extend that fixture to assert `action=preserve`, `density_above_cutoff`, and consistent flanked BED omission. Add `merged_with_preserved_call` at an appropriate boundary. |

### Build and repository readiness

| Review ID | Verdict | Handling |
| --- | --- | --- |
| M1 | Confirmed default feature, release-doc command list, docs script, and sidebar exposure. Help recommends inspection but does not explicitly label the command experimental. | Owner decision before release: gate it, or clearly label the command and consumer option experimental. Default availability alone does not demonstrate validation. Align the README's claim about validation of released commands with the decision. |
| M2 | Confirmed by both focused compile commands. `writes_text_outputs` unnecessarily includes outliers, and enabling plotters exposes a heatmap entry point used by midpoints. | Remove outliers from the writer feature group and gate the unused heatmap entry point to actual consumers, retaining its private tests. Do not hide warnings with broad allowances. Check affected feature combinations. |
| M3 | The three lints reproduce. The review's abbreviated list does not describe the full clippy outcome. | Apply the small equivalent expressions when touching the fitter. Track the unchanged shared heatmap failure separately. Do not claim a clean clippy result. |
| M4 | Confirmed portability problem. The linked report is ignored/untracked. The claim that Codex cannot read it is false in this workspace. The report itself describes historical settings and a different optimizer workflow, not validation of the current Rust fitter. | Publish a concise, tracked historical methods/results note with settings and provenance, or deliberately track the minimal reproducible evidence. Do not force-add the entire ignored results directory. Repair both existing links. |
| M5 | Confirmed missing Unreleased notes and main command-list entry. | Add them in the release/documentation pass after the experimental status decision. Describe the actual supported consumer only. |
| M6 | Confirmed public result exposes `FCoverageCounters` with unused GC/outlier/normalization fields. Only its base counters are updated. `add_fragment_clipped_to_core` and Cargo features also couple commands. | Low priority cleanup. Reuse `BaseCounters` directly if sufficient instead of inventing another wrapper. Consider moving the existing insertion helper with unchanged behavior when decoupling is justified. A broad shared-module refactor is not a prerequisite for correcting results. |

## Statistical decisions that need more care

### Identification depends on retained model support and constraints

An empty empirical bin still contributes information through normalization of the model's retained
support. Therefore "three occupied positive levels for ZINB, two for ZIP" is not a valid general rule.
It would reject the existing ZIP fixture with only retained coverage 1 even though its T1 is 7 and
its constrained pi-zero fit is finite and identified.

For T1 equal to 2, the retained distribution has just two categories, 0 and 1. For ZIP,

```text
q0 = [pi + (1 - pi) exp(-mu)] /
     [pi + (1 - pi) exp(-mu) (1 + mu)]
```

One observed zero fraction cannot identify both mu and pi. For example, both these underlying ZIPs
give retained probabilities `q0=0.9, q1=0.1`:

- `mu=1/9, pi=0`.
- `mu=1, pi=8/(exp(1)+8)`.

They have different untruncated tails. The current truncated ZIP solver chooses the pi-zero boundary
in this case. That is a selection from a tied family, not evidence that ZIP is identified.
The problem also applies to ZINB. Without retained zeros, the corresponding maximum can instead be
an unattained parameter limit. T1 equal to 2 must not be accepted as a generally identified free-pi fit.

For T1 equal to 3 and retained counts `n0,n1,n2`, let `ratio=n2/n1` when both positive counts exist.
The NB2 positive ratio fixes

```text
P(2) / P(1) = (1 + dispersion) * mu / [2 * (1 + dispersion * mu)]
mu = 2 * ratio / [1 + dispersion * (1 - 2 * ratio)]
```

This is a curve in mean/dispersion. Along it, the component's retained zero probability is

```text
q0_component = (1 + dispersion) /
               (1 + dispersion + 2 * ratio + 2 * ratio^2)
```

If the observed retained zero fraction is strictly above the Poisson value
`1 / (1 + 2 * ratio + 2 * ratio^2)`, a range of small positive dispersions can match the same
positive ratio and zero fraction by changing pi. This is an exact likelihood ridge.
When the observed zero fraction is below that value, the constrained boundary problem is different.
It is not correct to declare every T1-equal-to-3 fit unidentified.

An exact counterexample already resembles `truncated_profile_accounts_for_retention_when_solving_zero_inflation`.
Retained counts `12,2,1` at `0,1,2` have conditional probabilities `12/15,2/15,1/15`.
Both of these full models match them exactly:

- NB component `mu=1, dispersion=1` with `pi=1/2`.
- Poisson component `mu=1` with `pi=5/(exp(1)+5)`.

The first has `P(X>=3)=1/16=0.0625`. The second has
`P(X>=3)=(exp(1)-2.5)/(exp(1)+5)`, about 0.0283.
At alpha 0.04, their full-distribution thresholds are respectively 4 and 3 despite equal retained
likelihood. This is enough to show why accepting a tied extrapolation needs a stated policy.
This algebraic example fixes the truncation level. A complete two-stage regression must also
establish its first-fit T1 independently.

Recommended implementation direction:

1. Reject T1 equal to 2 for both free-pi families with a contextual insufficient-identification error.
2. For ZINB at T1 equal to 3, detect the interior ridge rather than rejecting every sparse histogram.
   The retained-count criterion above, equivalently a positive-pi truncated ZIP candidate in the
   regular two-positive-bin case, gives a concrete starting point. Handle missing bins and
   numerical boundary equality explicitly.

   > **Claude note:** This check costs almost nothing. `fit_zinb` already computes the ZIP boundary
   > candidate with `fit_right_truncated_zip(histogram, threshold)` (`zinb.rs:163-167`). Its lambda
   > matches the retained positive ratio, and its pi is positive exactly when the observed zero
   > fraction exceeds `1 / (1 + 2 * ratio + 2 * ratio^2)`. So the ridge test is
   > `truncation == Some(3) && n1 > 0 && n2 > 0 && boundary.zero_inflation > tolerance`. With `n2 = 0`,
   > the interior candidate's lambda hits `MIN_POSITIVE_LAMBDA` and its pi is negative. The returned
   > candidate is then the truncated Poisson with pi = 0, which matches the identified boundary case. Use a small explicit tolerance (for example `1e-12`), because at exact equality the
   > ridge has zero length and the fit is identified. The existing
   > `default_zinb_bam_run_reports_a_valid_empty_call_set` fixture has T1 = 3 with retained counts
   > `150, 150, 0`, so it stays valid under this rule. A blanket "T1 >= 4" rule would reject it
   > (changing its alpha to 0.01 would give T1 = 4).
3. Preserve identifiable pi-zero/dispersion-zero boundary cases, including the existing default
   ZINB no-call fixture, unless a deliberately stricter supported-use policy is chosen.
4. For larger retained support, retain explicit numerical-bound and convergence errors. Record
   occupied levels and near-optimal-start spread as diagnostics rather than treating them as proof.
5. Include model family, stage, T1, alpha, positive support, core/context, and the suggested remedy
   in the failure evidence. Do not raise T1 or switch models silently.

A blanket minimum T1 of 4 for ZINB would be a conservative policy, but would reject some valid
boundary fits. Discuss that tradeoff before adopting it. It is not required by the algebra above.

### Numerical stability diagnostics must concern comparable candidates

The objective is an average negative log likelihood. Neither the `1e-12` selection tolerance nor the
`1e-13` Nelder-Mead cost-spread tolerance is a confidence level. Positional observations also share
fragments, so nominal position counts do not provide independent sample sizes.

For each fitting stage, retain the boundary candidate and each start's convergence status, finite
cost, parameters, and rejected-evaluation count if recovery is added. Compute T2 spread for valid
near-optimal candidates, using an explicit recorded cost criterion. A poor converged local optimum
should not generate the same warning as a tied alternative tail. Identical starts' answers cannot
rule out an unsampled ridge.

Start with initial-vertex preflight and strict accepted-result evaluation. Only recover numerical
trial-domain failures that are understood and reported. Invalid histograms, overflow, and invariant
violations should still return errors. Avoid an opaque "all errors become infinity" adapter.
Tests for optimizer-outcome selection belong with private numerical logic. Do not expose internals
or add broad dependency injection just to manufacture an optimizer failure.

### T2 above T1 can be correct

A hand-derived identified ZIP example uses `n0=810, n1=180, n2=10` and alpha `0.001175`.
The complete fit is the Poisson boundary with mean 0.2. Its inclusive tail at 3 is about
0.00114848, giving T1 equal to 3.

No observed values are removed, but the second likelihood still conditions on X below 3.
Its pi-zero mean solves

```text
mu * (1 + mu) / (1 + mu + mu^2/2) = 0.2
0.9 * mu^2 + 0.8 * mu - 0.2 = 0
mu = (sqrt(1.36) - 0.8) / 1.8, approximately 0.203439
```

Its full Poisson tail at 3 exceeds 0.001175, while the tail at 4 is below it. Thus T2 is 4.
The positive counts at 1 and 2 and the zero constraint identify this boundary solution.
This example also separates proper truncation from "just fit the retained data as complete".

## Failure diagnostics and optional global fitting

These changes must be designed together. Current output staging protects individual files from
being half written, but `FinalOutputFiles::move_into_place` is explicitly best effort across files.
Do not claim that publishing the whole output set is atomic.

Recommended failure contract:

- A fitting or ambiguity-policy failure remains a nonzero exit. It never publishes newly computed
  keep weights or actionable BEDs.
- Preserve completed diagnostics in a clearly named, unique failed-run location. Minimum contents
  are the complete histogram, actual settings, failed stage/context, and full error chain.
  Identify this location in the error message. Avoid mixing new failure files with an old successful prefix.
- Preserve completed model results when practical, but do not require a large parallel-fit result
  redesign merely to save the histogram. Partial model tables must identify missing/failed contexts.
- For ambiguity error mode, collect decisions over every selected chromosome before returning the
  policy failure. Write the full audit and summary with an explicit failed/withheld-publication status.
  Settle how the audit represents blocked actions. It must not imply that correction files exist.
- Keep actionable and diagnostic output registration separate. Flushing, closing, and finishing
  compressed writers must succeed before treating a diagnostic as complete.
- If diagnostic publication also fails, report that failure in addition to the original cause.
  Do not hide the fitting error or silently leave users without the promised evidence.
- Leave previous successful outputs intact on failure, and make their different run identity clear.
  Test reused prefixes as well as empty output directories. Do not solve stale-output ambiguity
  by deleting unrelated or previously valid files.

> **Claude note:** Two details need a decision before implementation.
>
> - A "unique" failed-run location (for example with a random or time suffix) piles up one directory
>   per failed attempt, and users usually retry several times. A fixed per-prefix name, such as
>   `<prefix>.outliers.failed_run/`, replaced as a whole on each new failure and never deleted by a
>   successful run, keeps one current failure record without deleting valid results. The error
>   message should print its path. Either way, document the rule.
> - For the ambiguity error, a concrete audit representation: write every call with
>   `action=withheld`, keep `decision_reason` as the policy would have produced it, and add a
>   `# publication=withheld_by_ambiguous_region_action_error` comment. Then no row implies that weights
>   or BEDs exist.

For global fitting, determine whether any selected chromosome cannot meet the minimum physical span
or total midpoint fragment support. Since a local context can expand to the whole chromosome,
these totals establish whether any fallback is required without an additional BAM pass.
Pass the selected model family explicitly instead of obtaining it from an obligatory global model.

Attempt the pooled model for diagnostics even when unused. If only that fit fails, recommend a
prominent warning plus a structured unavailable status and recorded reason. Keep observed counts
and local tracks. Do not construct a `CoverageParameters` full of NaNs, draw nonexistent fitted curves,
or disguise failure with zero expected counts. If any core needs fallback, failure remains fatal.
This is an explicit change to the current required-global-diagnostics contract and should be
confirmed before implementation.

For manual probability split jobs, successful local models should agree when their chromosome
inputs and contexts agree. Global fallback remains dependent on the pooled population. Detecting
once across the desired chromosomes and letting `fcoverage` consume a subset avoids both issues
with the existing loader. Do not add a new denominator flag or global-model exchange format merely
to document the current behavior.

## Biological validation priorities

### Background contributions outside calls

For stationary fragment starts, a fixed contiguous counted length `fragment_length`, and an isolated
call with weight `keep_weight`, a base at distance `distance` outside the edge has
`fragment_length - distance` of its `fragment_length` possible covering starts touching the call,
for integer distances from 1 through the fragment length. Its expected background contribution is

```text
background_mean *
  [1 - (1 - keep_weight) * max(fragment_length - distance, 0) / fragment_length]
```

The edge approaches `background_mean * keep_weight` and recovery takes a fragment length. This
describes the background component, not necessarily the total observed corrected coverage.
Artifact fragments extending outside the call, variable lengths, nearby calls, alignment gaps,
nonuniform starts, trimming, GC weights, and normalization can change the observed shape.

For an ordinary molecule length distribution with constant start density, the affected background
fraction is `E[max(fragment_length - distance, 0)] / E[fragment_length]`.
The coverage weights are therefore length biased. A simple average over molecules of
`max(fragment_length - distance,0)/fragment_length` would use the wrong coverage weighting.
This is a mechanistic concern for local fragment length distributions in future consumers.
The `lengths` integration itself remains deferred.

Required comparison artifacts:

- Raw, outlier-weighted, and removed positional coverage with exact calls, weights, active means,
  T1/T2, masks, and original fragment endpoints.
- Loss of known ordinary fragment contributions inside calls and by distance outside calls,
  split by fragment length. Separate this from removed technical-fragment contributions when truth exists.
- Features used by the intended analysis, including whether a new dip appears at a promoter or
  other profiled element. Coverage plots alone cannot validate a nucleosome or footprint interpretation.
- Calls with a single weight, several nearby weights, gap-only span overlaps, and preserved calls
  next to actionable calls. Omission of a preserved call does not exempt its fragments from a
  different actionable overlap.

Do not change whole-fragment weighting in this rectification pass. If validation rejects it,
compare mathematically defined alternatives separately. Diluting a whole-fragment weight by overlap
fraction changes the achieved target and cannot inherit the current mass-ratio derivation unchanged.

### Sparse calls in genuine gains

Use independently annotated gains and controls. Cover narrower gains, broad moderate gains with
sparse excursions, gains straddling model cores, and technical spikes within gains. Report the
elevation's physical span separately from the called-base fraction and the flanked merged span.
Determine whether failure occurs in initial fitting, tail exclusion, final calling, preservation,
or whole-fragment application.

Record nominal sequencing depth and counted positional fragment coverage separately. The review's
30x example is not a calibrated prediction unless mean fragment coverage, dispersion, pi, context
contamination, and alpha are specified. A gamma approximation at a probability near 1e-10 is not
sufficient evidence for a quantitative operating boundary.

Use biologically coherent mixtures. For a diploid background and tumor population with fraction
`tumor_fraction` and local copy number `copy_number`, under equal per-copy sampling and no ploidy
rescaling, the expected relative coverage is `1 + tumor_fraction * (copy_number - 2) / 2`.
Arbitrary gains cannot be described as high tumor fraction or a particular copy number without those
assumptions. Record independent evidence and uncertainty when real locus truth is unavailable.

Start with a small set of decisive loci and hand-defined fragment examples. Then expand through the
existing 0.1x to 200x validation requirement. Include lengths around 2.5 kb and 5 kb preservation
boundaries as well as 10 kb, 100 kb, and 1 Mb events. Vary elevation strength and fragment lengths.
Historical 10% preservation simulations are useful counterexamples, but their outcomes do not
establish what the current 5% rule or current Rust optimizer will do.

Measure biological mass lost, distortion of local features, technical signal removed, preserved
technical excess, fit failure rate, and threshold stability. Decide acceptable tradeoffs with the owner
before using these comparisons as a release gate. Do not invent a universal acceptable percentage
or infer broad validation from a favorable BAM example.

### GC-associated risk and reference populations

At matched loci compare raw, outlier-only, GC-only, and combined tracks using identical accepted
fragments and masks. Stratify by library bias, fragment length, GC content, reference gaps,
mappability, and independent copy-gain evidence. Matched backgrounds are necessary to interpret
CpG-island overlap. Selection of a called region by its high raw coverage already biases naive
enrichment comparisons.

Include deliberately selected short/alternate/mitochondrial contigs only when testing that use case.
The pooled fallback can be biologically inappropriate even when fitted successfully. Alpha adjustment
does not fix a mismatched reference population. Do not treat preservation of an almost entirely called
mitochondrial contig as evidence that the pooled autosomal-like baseline was appropriate.

## Concrete test work

All new expectations must be derived before execution. The user runs tests. Keep public behavior
under `tests/` and private numerical/aggregation checks in sibling `*_tests.rs` files. No visibility
widening, duplicate implementation in another language, or fragments below 10 bp. Prefer 30 bp or
longer so fixtures satisfy command inclusion filters.

1. **Target semantics and preserved audit.** Fix test finding T1's column/comment. Extend the
   existing nonzero-pi roundtrip with an explicit inequality against the unconditional target if
   needed. Extend the preserved command fixture for action/reason and both BED omissions.
2. **Exact identification cases.** Test the tied ZIP and ZINB examples above through private fitting
   logic. Test a valid sparse boundary case and near-equality of the retained zero constraint.
   Check that errors identify model, stage, T1, and alpha. Add the identified T2-above-T1 ZIP example.
3. **Positive-dispersion two-stage path.** Use a deterministic stepped coverage fixture with exact
   geometric/NB ratios and realistic fragment spans. Existing geometric private fixtures provide
   analytic PMF/truncation expectations, but do not by themselves establish a complete initial fit.
   Independently derive or bound the initial T1, then verify the truncated refit, T2 and predecessor
   tails, positive dispersion, calls, masses, and output columns. Do not choose expected thresholds
   by observing the Rust fitter. If the complete initial MLE cannot be established by hand, state
   that limit and choose invariant assertions or discuss an independently justified reference.

   > **Claude note:** The complete initial MLE can be set by hand with a private histogram that is
   > exactly proportional to a ZINB. Use component NB with `mu = 1, dispersion = 1` (geometric,
   > `P(Y = k) = 2^-(k+1)`) and `pi = 1/3`. Then `P(X = 0) = 2/3` and `P(X = k) = 2^-k / 3` for k >= 1.
   > Counts `c_0 = 2^41` and `c_k = 2^(40-k)` for `k = 1..40` match `N * P(X = k)` exactly for
   > `N = 3 * 2^40`. The only deviation is the missing tail above 40 (relative mass about 2^-41) plus
   > the one-count difference in N. So the complete fit must return `(1, 1, 1/3)` within numerical
   > tolerance, with interior pi (zero fraction 2/3 above the component zero 1/2).
   >
   > - `S(k) = (2/3) * 2^-k`. With alpha `1e-6`, `S(19) = 1.27e-6` fails and `S(20) = 6.36e-7` passes,
   >   so T1 = 20. The margins of about 1.27x and 1.57x tolerate parameter error far above fitter noise.
   > - The retained part `k = 0..19` is exactly the truncated model, so the refit returns the same
   >   parameters and T2 = 20.
   > - `second_fit_retained_positions = 3 * 2^40 - 2^21 = 3,298,532,786,176`.
   > - The observed initial tail is `2^21 - 1 = 2,097,151` positions, expected about `2,097,152`.
   > - The active mean is 1 and the unconditional mean is 2/3, so the same fixture also separates the
   >   `local-mean` target from the unconditional mean at the private level.
   >
   > Call `fit_two_stage(&histogram, 1e-6, OutlierModel::Zinb)` from `distribution_tests.rs`. No BAM
   > is needed. Keep the command-level positive-dispersion run, but with invariant assertions only.
4. **Partially split model core.** Use ZIP, 1 Mb processing tiles, a 600 kb stride, a 2,000,030 bp
   chromosome, a whole-chromosome fitting context, minimum fragment support 1, manual alpha 1e-6,
   threshold targeting, and no flank. Accept exactly 30 bp fragments. Include baseline fragments
   at 100 and 700,000 and twenty fragments at 999,990, all with 30 bp spans. The model core
   `[600000,1200000)` must have 599,940 zeros, 30 ones, 30 twenties, and 21 midpoint-owned fragments.
   Its left tile piece has 399,960 zeros, 30 ones, 10 twenties, and one owned fragment. Its right
   piece has 199,980 zeros, 20 twenties, and twenty owned fragments. The pileup midpoint is 1,000,005.
   Assert these counts and the joined call `[999990,1000020)`, then compare histogram/model/call content with a
   different legal tile size. Both stride and bin-size constraints must remain valid.

   > **Claude note:** The counts check out. Specify the bin size, because "whole-chromosome fitting
   > context" depends on it. Cores are `[0,600000)`, `[600000,1200000)`, `[1200000,1800000)` and
   > `[1800000,2000030)`. Any `--bin-size` in `(1800000, 2000030]` forces every core's context to the
   > whole chromosome. A value above 2,000,030 makes every core fall back to the global model. Use
   > `--bin-size 2000030`. Hand-derived models for this fixture: the initial ZIP has a positive mean of
   > 660/90 = 7.33, so T1 is far above 2. The refit keeps only the 60 positions at coverage 1, so it
   > lands on the Poisson boundary with `lambda = 60 / (retained positions)`, about 3e-5. Then T2 = 2,
   > because `P(X >= 1)` is about 3e-5 > 1e-6 and `P(X >= 2)` is about 4.5e-10. So the coverage-1
   > positions are not called. For the comparison tile size, 2 Mb is legal and also unaligned with
   > the 600 kb stride (fewer than ten strides fit), so the last core straddles a different boundary.
5. **Automatic alpha.** Clear the base config's manual probability and mask `[900,1000)` in the
   existing 1 kb paired fixture, without masking its fragments. There are 900 eligible positions,
   including uncovered bases. Expect alpha `1/900`, or `2/900` with multiplier 2. Compare the latter
   with manual `2/900`, allowing probability-mode and inactive-setting metadata differences.
   The private writer test does not replace this command-level check.
6. **Ambiguity actions and failures.** Use known call percentages for preserve/correct/error and
   shared flank membership. Assert every audit reason and withheld/actionable row. On fitting
   failures and ambiguity errors assert non-success, complete promised diagnostics, no new weights
   or BEDs, cleanup of incomplete staging, and integrity of earlier successful-prefix outputs.
   Use deterministic internal failure selection checks for difficult optimizer outcomes.
7. **Unused versus required global model.** At the private orchestration boundary, supply an explicit
   failed global-fit outcome. All-local contexts must succeed with recorded warning/status. A context
   requiring fallback must fail. A failing supported local fit must still fail rather than use global
   parameters. Add the strongest feasible public fixture without searching for an unstable optimizer failure.
8. **Unpaired and plotting boundaries.** Add an unpaired successful run with hand-counted coverage.
   For plotting, cross a real chunk boundary with fragments and a mask, and compare mean, peak,
   peak position, and eligible counts across chunkings. Test fine close-up bins from actual positions,
   rather than pretending a coarse-bin mean is fine-resolution evidence.
9. **Shoulder calculation.** A custom keep-weight interval `[100,101)` of weight 0.2 and stationary
   30 bp fragments starting once per base gives background coverage 30. At distances 1, 15, and 30
   outside either edge, expect 6.8, 18, and 30. Keep away from chromosome edges and use no other
   corrections. A 1 bp call interval is valid even though 1 bp fragment fixtures are forbidden.
   This tests application, while existing detector roundtrips test file production.
10. **GC-mode counter consistency.** Place a zero-weight fragment beyond the tile core whose full
    span extends outside the fetch, plus a zero-weight fragment contributing to the core. Compare
    plain, GC-tag, and neutral file-GC modes after adopting a common encounter rule. Counted coverage
    must agree and zero weights must still bypass invalid-GC and normalization accounting.

Keep runtime appropriate to each boundary. Tiny algebraic tests should not require a BAM.
Do not turn every documentation improvement into a snapshot test.

## Implementation sequence and acceptance

### First pass: expose failures and protect statistical behavior

- Correct the weak target assertion and extend the preserved audit test.
- Define the failure-diagnostic contract, including reused prefixes and ambiguity audit status.
- Add the exact identifiability regressions, then implement the agreed guard without rejecting
  valid boundary fits by occupied-level count.
- Decouple global fallback need from diagnostic fitting, retaining strict local and required-fallback failures.
- Save failure evidence before pursuing optimizer tuning.
- Add positive-dispersion two-stage/command coverage and the split-core regression.

Acceptance: the source makes identification failures explicit, successful local inference does not
depend on an unused pooled fit, and a failed expensive run leaves usable diagnostics without new
actionable files. Author and compile the behavioral tests. Test execution remains with the user.

### Second pass: inspect scientific effects with adequate reporting

- Complete split-alpha and unpaired coverage.
- Fix stage-specific error guidance, plot labels, chromosome labels, and the missing legend.
- Improve close-up data resolution and draw all relevant candidate decisions.
- Measure B1-B3 on the decisive loci before proposing a new weighting or preservation algorithm.
- Record near-optimal-start sensitivity and numerical limitations in inspection evidence.

Acceptance: plots and tables explain the decisions without implying biological truth, and the owner
can assess measured signal loss and technical benefit. Statistical fit quality alone is insufficient.

### Third pass: prepare the release surface

- Finalize CLI terminology with the owner, then update serialization, metadata mapping, help,
  README generation, guide, and generated CLI docs consistently.
- Add experimental status or feature gating as decided. Update release notes and command listing.
- Fix focused build warnings and the three equivalent-expression lints.
- Repair the historical evidence links without treating old results as current validation.
- Document optional-plot cost and stale-file behavior. Do performance/concurrency work only if measured.
- Consider the counter/helper coupling separately if it interferes with the chosen release boundary.

Acceptance: documentation says what is actually supported and scientifically established. A clean
compile is not a declaration of merge readiness. Preserve the existing deferred normalization scope.

After production changes, run `cargo check --offline --locked`. After test changes, run
`cargo check --offline --locked --tests --features "testing cmd_outliers"`. Keep the focused outliers
and fcoverage feature checks. Use `cargo check --all-features` after any major refactor, as required
by repository instructions. Do not run tests. Re-read every changed file before handing work back.
Only distill implemented, lasting invariants into `.AI/docs/specs/`.

## Troubleshooting content to write

- **Initial complete local fit fails:** inspect its saved histogram and context. Larger
  `--bin-size` or `--min-context-fragments` changes local pooling but is not guaranteed to help.
  Check technical extremes, mixed biological rates, masks, and the actual bound/convergence reason.
- **Initial pooled fit fails:** changing local context settings does not change that histogram.
  Explain whether it was required for fallback or only unavailable for diagnostics.
- **T1 leaves no support or an unidentified second fit:** a smaller manual tail probability makes
  T1 larger for the fixed initial fit. Report current alpha and T1 so this advice is actionable.
  It does not guarantee sufficient evidence or repair every numerical failure.
- **Refit reaches a bound or has unresolved convergence:** inspect retained support and recorded
  candidate outcomes. Do not automatically increase numerical bounds or hide the failure with global reuse.
- **Using `--model zip`:** describe this as a different statistical assumption and comparison,
  with less allowance for ordinary overdispersion. It is not a generally validated recovery mode.
- **Ambiguity policy error:** point to the full diagnostic audit and explain the withheld correction.
  Switching to `correct` accepts a biological signal-loss tradeoff. Do not recommend it as a routine fix.

## Questions and objections for Claude

These are prepared for discussion. They have not been sent through another application.

1. **ZIP at T1 equal to 2:** do you agree that the Bernoulli retained distribution also leaves mu/pi
   unidentified, as the two exact ZIP examples above show? Your statement that ZIP is fine needs
   qualification beyond the particular T1-equal-to-7 integration fixture.
2. **Occupied levels:** what argument justifies requiring two occupied positive levels for ZIP or
   three for ZINB? Empty levels in the allowed retained support enter through normalization.
   Can we instead use the exact ridge criterion while preserving valid constrained boundary fits?
3. **Global-fit risk:** can you supply a concrete histogram/BAM where the pooled fit hits the bound
   while every required local fit succeeds? The unnecessary dependency is proven, but the claimed
   pileup-driven failure frequency is not. Which failure modes should be diagnostic-only?
4. **Split-run equivalence:** please qualify both "drops in every core" and the manual-alpha remedy.
   Discrete thresholds can stay equal, refitting changes with T1, and pooled fallbacks differ by selection.
5. **Trough shape:** can you distinguish expected ordinary-background loss from total observed
   corrected coverage, and state the assumptions behind the V-shape and 5% example? Do you have
   data showing the proposed promoter/footprint distortion under the current filters?
6. **Gain approximation:** what component mean, structural-zero fraction, actual fragment coverage,
   tumor mixture, context contamination, and tail calculation support the 3x/30x example?
   Treat it as a scenario to examine until those are specified.
7. **Short contigs:** please replace "sample-wide autosomal model" with the all-selected-contig pool
   and make the mitochondrial call/preservation example conditional on observed eligible coverage.
8. **Trial errors:** which actual in-domain evaluations fail? Can we prevent initialization panics
   and classify numerical trial failures without suppressing invalid-data or implementation errors?
9. **Build evidence:** which exact clippy invocation and exit status did you observe? The invocation
   recorded here fails at unchanged `heatmap.rs:476`, in addition to the cited outlier lints.
10. **Historical report:** it is readable locally but untracked. What minimal methods/settings/results
    should be tracked for review? Its historical 10% preservation and alternate optimizer outcomes
    cannot be presented as evidence that the current Rust default is validated.

## Recheck of the review's "found consistent" section

These are reasons to preserve the existing invariants, not blanket numerical or biological approval.

| Area | Independent check and remaining limit |
| --- | --- |
| Complete ZIP MLE | Positive conditional mean and the pi-zero Poisson boundary follow from `fit_zip`. The complete finite fit is distinct from the unidentified T1-equal-to-2 truncated case. |
| Right-truncated ZIP | Re-derived the retained normalization and `pi=(q0*F-P0)/[(1-P0)-q0*(1-F)]`. The algebra is correct. Identifiability of the chosen support is a separate requirement. |
| Profiled ZINB | Re-derived the pi constraint and factorization into the retained zero fraction and conditional positive likelihood in `ProfiledHistogram::evaluate`. Boundary and interior costs use the same per-retained-position scale. |
| NB2 log mass | Re-derived the rising-factorial/Stirling difference and checked the coefficients and small-argument series in `component_log_probability`. This is not an empirical error bound for the entire parameter domain. |
| NB survival | Confirmed the inclusive beta mapping with shape `1/dispersion` and argument `dispersion*mu/(1+dispersion*mu)`. statrs uses direct beta evaluation in the relevant upper-tail branch and silently returns its current approximation after 140 iterations. The claimed relative-error magnitude and unreachable iteration cap were not independently established here. |
| Positive variance | From `E[X^2 given X>0]=(mu+(1+dispersion)*mu^2)/(1-P0)`, obtained `positive_mean*(1+dispersion*mu+mu-positive_mean)`, matching the code. |
| Threshold search | Confirmed inclusive survival and minimum-integer search. Numerical monotonicity/accuracy remain assumptions of the evaluated tails. Even a small relative error can change an integer cutoff near equality, so "harmless" is not a universal guarantee. |
| Preservation heap | Confirmed lazy expiry is valid for sorted exact calls and monotonically ending forward windows. Supported evidence outranks unsupported scores, with the documented earliest-window tie break and supported backward fallback. |
| Preserved groups | Confirmed seeds retain merged membership, sort by window span, and merge intersecting spans rather than touching spans. Plot context does not alter grouping or decisions. |
| Tile order and joined calls | `build_tiles` produces chromosome/position order. Indexed parallel collection preserves that order when successful. Joining touching pieces sums masses and takes maximum peak coverage. The partial-core regression remains necessary. |
| Fragment support | Midpoints in the tile core and outside blacklists are assigned once. Partial histogram merges sum counts/support. The once-per-fragment claim assumes complete accepted-pair retrieval through the configured halo, which the split-core fixture should exercise. |
| Coverage representation | Coverage uses f64 deltas and casts to f32. Eligible histogram values at or above 2^24 are rejected. Blacklisted values do not enter fitting. |
| Contigs and metadata | `bam_contigs_info` returns selected contigs only, and output metadata sorts their names/lengths. The histogram pool is consequently selection-dependent. |
| fcoverage application | Verified all output modes share `process_tile`, lookup uses original `pos` to `reference_end` bounds before trimming, one-base overlap is admitted by `1/(maximum_fragment_length+1)`, returned scan indices match the weight vector, and minimum weights combine. |
| Pairing and gaps | Both commands use segment trigger 1 and the same inter-mate-gap policy. Pairs complete when the second mate is read. With accepted spans bounded by the maximum fragment length, the overlap look-back preserves the required start-order slack. Zero weights precede normalization/GC counting, subject to the reporting discrepancy in 1.8. |

Do not replace statrs or the optimizer solely because an unmeasured numerical concern exists.
Extend analytic geometric/integer-shape and near-boundary tail checks, then investigate independent
high-precision evidence if those checks or real diagnostics expose a failure.

## Claude review of this plan

Date: 2026-10-19. Reviewed against `9e9bdac5`. The working tree contains only this plan and the
original review. No code was changed.

### Verdict

Use this plan as the working follow-up. Its corrections to the original review are right, and its
statistical sections improve on the review. The derivations below were re-checked by hand:

- **ZIP at T1 = 2:** both examples give `q0 = 0.9`. For `mu = 1, pi = 8/(e+8)`, numerator and
  denominator are `9/(e+8)` and `10/(e+8)`.
- **ZINB at T1 = 3:** the ratio curve `mu = 2 * ratio / [1 + dispersion * (1 - 2 * ratio)]` and the
  component zero probability `(1 + dispersion) / (1 + dispersion + 2 * ratio + 2 * ratio^2)` are
  correct.
- **The `12, 2, 1` counterexample:** both models reproduce `12/15, 2/15, 1/15`. Their tails at 3 are
  0.0625 and 0.0283, and their thresholds at alpha 0.04 are 4 and 3.
- **T2 above T1:** the complete fit is the Poisson boundary 0.2 with T1 = 3. The refit gives
  `mu = 0.203439`, with tails of about 0.001207 at 3 and 6.2e-5 at 4, so T2 = 4.
- **Fixtures:** the split-core counts (test 4), the automatic alpha `1/900` (test 5), and the
  shoulder values 6.8, 18 and 30 (test 9) all hold. Test 5 still calls the pileup at alpha `1/900`
  and `2/900`, since T2 = 3 in both cases.
- **Codex's other claims checked:** `write_heatmap` is used by `midpoints` (M2), and
  `move_into_place` is documented as best effort across files.

### Corrections to the original review, accepted

- **1.3:** "ZIP is fine here" was derived for the fixture with T1 = 7 and only coverage 1 observed. At
  T1 = 2 the ZIP positive conditional is degenerate and `(mu, pi)` is not identified. The
  occupied-level rule was a heuristic with no identifiability argument and is withdrawn.
- **1.7:** "T2 drops in every core" is wrong as stated. Integer thresholds can stay equal, and T1
  changes the retained data. A manual alpha removes only the alpha difference, not the
  pooled-fallback difference.
- **B1:** "V-shaped troughs next to every call" was too strong. The formula gives the expected
  ordinary-background contribution under its stated assumptions. With a fragment length
  distribution, the affected fraction `E[max(l - d, 0)] / E[l]` is convex in `d`, so the recovery is
  not a literal V.
- **D4:** "sample-wide autosomal model" should be "the pool of all selected contigs". chrM being almost
  entirely called depends on its observed coverage.
- **M3:** the original clippy run (`cargo clippy --offline --locked --features "testing cmd_outliers"
  --tests`) was piped through `grep` and its exit status was not checked. Rerunning it gives exit
  101: a denied lint at `src/shared/plotters/heatmap.rs:476` (lib) and 7 lib-test errors, including
  `src/commands/midpoints/counting_by_group_tests.rs:364`. None are in outliers code. "3 trivial
  lints" misdescribed the outcome.
- **M4:** the report is readable in this workspace. The link is broken only for anyone working from
  the tracked repository.

### Remaining issues with this plan

1. **Three overlapping documents.** This plan, `outliers_review_rectification_plan.md` and
   `outliers_real_data_validation.md` now repeat each other. This is worst for B1-B3 versus the
   validation list, and for implementation items versus the checklist. Once the owner decisions
   below are made, move the accepted implementation items into the rectification checklist and the
   biological items into the validation list. Then shrink this file to a dated review record, or
   delete it.
2. **Owner decisions are scattered** across "owner decision", "confirm before implementation" and
   "discuss". They are collected below.
3. **Effort calibration.** 1.8 and test item 10 fix an informational counter and need a GC-package
   fixture. Documenting "tile encounters, not comparable across GC modes" is enough for the first two
   passes. Move the code fix to the third pass, or drop it.
4. **Test items 3 and 4** were underspecified. See the inline notes for an exactly hand-derivable ZINB
   two-stage histogram and for the missing `--bin-size`.
5. **T2 above T1 in practice.** The example uses `alpha * positions` of about 1.2. A 5.5 Mb context at
   the default automatic alpha (about 3.6e-10) has `alpha * positions` of about 0.002. The first-fit
   tail mass above T1 is then negligible, so the truncation barely changes the likelihood. Expect
   T2 > T1 to be rare in real contexts. Not clamping it is still right.

### Owner decisions needed

1. **Failure diagnostics.** Change the "nothing published on error" contract. Decide the location
   rule (unique per attempt or fixed per prefix, see the inline note) and how the audit represents
   withheld actions.
2. **Global fit.** When no core needs a fallback, make a pooled-fit failure a warning with an
   explicit unavailable status.
3. **Identifiability policy.** Either the ridge rule (recommended, cheap with the existing boundary
   candidate) or a blanket ZIP T1 >= 3 and ZINB T1 >= 4. The blanket rule needs the default ZINB
   no-call fixture's alpha changed.
4. **Ambiguity error output.** Should `--ambiguous-region-action error` publish the full audit and
   summary, with actions shown as withheld?
5. **Release surface.** Feature gate the command, or label it experimental, before any release.
6. **Naming.** CLI renames (N1-N5) and any metadata key changes, settled together in one pass.
7. **Release-gate criteria** for the B1-B3 measurements (acceptable signal loss versus technical
   benefit).
8. **1.8.** Fix the counter or only document it.

### Answers to the questions for Claude

1. **ZIP at T1 = 2:** agreed. See the corrections above.
2. **Occupied levels:** there is no valid argument. Empty levels inside the retained support enter
   through normalization. Use the support-based rule, which keeps identified boundary fits such as
   the T1 = 7 single-level fixture and the default ZINB no-call fixture.
3. **Global-fit risk:** there is no concrete histogram. Pileup-driven failure was an inference. The
   defect is the unnecessary dependency, whatever the failure frequency. Make every failure of the
   pooled fit diagnostic-only when no core needs a fallback: initial fit, truncated refit, bounds and
   convergence alike. Keep all of them fatal when any core does need a fallback.
4. **Split runs:** qualified as in the corrections. The recommended workflow is to detect once over
   the intended chromosomes and let `fcoverage` consume a subset. The existing loader already supports
   that.
5. **Trough shape:** the formula covers the ordinary-background component only. It assumes stationary
   starts, a fixed fragment length, contiguous counted spans (no `--ignore-gap`, no D/N gaps), an
   isolated single-weight call, and no GC, normalization or trimming. Total observed coverage also
   contains the artifact fragments' own downweighted contributions outside the call and the effect of
   other calls. The 5% figure is the hypothetical `w = 0.05`. There are no data on promoter or
   footprint distortion. It is a mechanistic concern to measure.
6. **Gain example:** the assumptions were component mean 30 (counted positional fragment coverage,
   not nominal depth), pi = 0, dispersion 0.06, an uncontaminated context, and alpha about 3.4e-10.
   The calculation used a Wilson-Hilferty gamma quantile (z of about 6.17) and ignored Poisson
   sampling, giving T2 of about 3.35 times the mean. Under the plan's mixture formula, relative
   coverage 3 needs `tumor_fraction * (copy_number - 2) = 4`. That is copy number 10 at tumor fraction
   0.5, or copy number 22 at 0.2. So the scenario concerns high-level amplification at high tumor
   fraction. Lower fitted dispersion lowers T2 relative to the mean and widens the risk band to more
   moderate gains. Treat it as a scenario, not a calibration.
7. **Short contigs:** agreed. See the corrections above.
8. **Trial errors:** I know of no in-domain input that fails. The `Err` sites in `evaluate` are the
   finiteness check on `log_retained`/`log_positive` and the `zero_inflation` range check
   (`zinb.rs:94-97`, `:118-121`). Agree with preflighting the three initial vertices with contextual
   errors before `Executor::run`, and with keeping invalid-data and invariant errors fatal.
9. **Build evidence:** see the M3 correction. My invocation also exits 101 on unchanged shared and
   test files.
10. **Historical report:** track a short note with:
    - the generator model and parameters (means, CV values, event heights and widths)
    - replicates and seeds
    - the fitting implementation used (the Python harness and optimizer, not the Rust fitter)
    - the preservation settings (10% cutoff)
    - an outcome table
    - date and commit
    - an explicit statement that these results do not validate the current Rust defaults or the 5%
      rule

    Describe the harness rather than force-adding the ignored directory.
