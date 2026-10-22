# Outliers remaining work

Current behavior belongs in the [spec](../specs/outliers.md). This file keeps unresolved decisions
and deferred work only. Scientific investigations are in the
[validation plan](outliers_real_data_validation.md), with a
[short biological discussion report](outliers_biological_release_criteria_report.md) for the owner
and supervisor.

## Owner review

Finalize CLI terminology during manual code verification. These are proposals, not agreed names:

| Current argument | Proposed argument |
| --- | --- |
| `--protection-window-size` | `--preservation-window-size` |
| `--min-protection-positions` | `--min-preservation-positions` |
| `--bin-size` | `--min-context-bases` |
| `--ambiguous-region-action` | `--preservation-action` |
| `--max-outlier-percentage` | `--max-called-percentage` |

Review the help headings for `--target` and `--blacklist-flank` at the same time. Keep `--stride`
unless there is a concrete reason to rename it. Apply approved names together across config
serialization, metadata, producer/consumer checks, help, the generated README, guides and generated
references. Preserve the distinction between GC-estimator outlier controls and fragment keep weights
during producer integration.
The context option and its metadata must describe eligible bases, including zeros, rather than
physical span. Review `minimum_context_span` in the same naming pass.
Clarify whole-core rounding in that review. A 5 Mb eligible-base minimum gives 5.5 Mb centered
contexts on an unmasked 500 kb grid because expansion retains a central core and equal numbers
of neighbors. The owner flagged this distinction as unintuitive. Keep the current expansion
rule until the intended option semantics are settled.

Feature/default/release-status changes remain the owner's decision. The owner and supervisor also
need to choose biological retention criteria, supported event sizes, sample types and depths.
The biological report proposes discussion questions rather than selecting those criteria.

## Required integration before release

Support optional outlier keep weights across fragment-counting commands, including `gc-bias`,
`coverage-weights` and `fragment-count-weights`. Complete integration across all applicable
commands before public release. There is no phased public rollout.
Blacklisting remains an alternative chosen for the analysis. Further empirical comparisons are
not prerequisites for this design decision.

In an outlier-handling workflow, `gc-bias`, `coverage-weights` and `fragment-count-weights` must
estimate corrections using the outlier keep weights, even when the final measurement blacklists
the called regions. Do not substitute sample-specific outlier blacklisting for weights in these
estimators. Ordinary blacklists still apply. The final measurement can use weighting or exclusion.

- **GC bias:** Excluding regions only from observed counts changes their genomic representation
  relative to the precomputed reference distribution and can introduce an artificial GC bias.
  Weighting reduces excess observed mass while retaining the intended regional contribution.
  Expected reference counts remain unweighted.
- **Scaling weights:** Unweighted peaks can distort the broad smoothed background and therefore
  the corrections applied to neighboring regions. Replacing weighting with exclusion changes
  the represented genomic area and can lower stride totals, particularly for fragment-count
  weights. Estimate the broad background from outlier-weighted signal regardless of whether the
  final measurement uses weighting or exclusion.

- [ ] In `gc-bias`, weight observed fragment counts in both streaming and BED/global branches,
  multiplying existing window-assignment fractions. Derive fragment length frequencies from
  weighted counts. Keep expected reference counts unweighted and retain existing window
  normalization. Keep correction-matrix outlier controls distinct from fragment keep weights.
- [ ] In both scaling-weight commands, forward keep weights to internal `fcoverage` before
  stride aggregation, smoothing and global normalization. Preserve existing count and support
  definitions. Document the effects of partly and fully masked strides separately from weighting.
- [ ] Extend other fragment-counting consumers with the same original `pos` to `reference_end`
  overlap rule, minimum overlapping weight and early zero-weight exclusion. Share the existing
  resolver when another consumer needs it. Keep encounter counters distinct from weighted mass
  and leave restore-mean normalization statistics unchanged.
- [ ] Complete user-facing examples for the integrated release. Explain when exclusions are
  suitable for coverage profiles and regional estimates, and when fixed-track counts or distinct
  regional motif/length composition motivate retaining weighted contributions. Weighting assumes
  the retained composition remains useful. It cannot recover true composition from technical
  contamination. Explain matching reference masks for GC exclusion and usable-area effects on
  count contrasts. A small genome-wide excluded fraction need not be small within a track.

## Deferred work

- **GC-mode exclusion counters, low priority.** Zero-weight exclusion counts are tile encounters.
  Padding-only fragments can be counted differently with file-GC, GC-tag and plain counting,
  without establishing different counted coverage. If addressed later, agree on encounter
  eligibility and derive a boundary fixture with neutral GC weights. Preserve early zero-weight
  exclusion before GC and normalization work. The owner rarely sees GC exclusions and does not
  regard this as a current focus.
- **Numerical calibration.** Assess optimizer sensitivity, parameter uncertainty, larger retained
  supports, near-boundary fits and extreme-tail accuracy. Ridge guards and agreement between starts
  do not establish identification or calibrated probabilities. Use analytic cases and investigate
  independent numerical evidence when diagnostics expose a problem.
- **Plotting and performance.** Consider plot concurrency, caps or changes to BAM rereads only
  after a measured need. Any cap must make omitted plots explicit.
- **Shared test feature gates.** Minimal command test builds reference helpers from disabled
  commands in shared CLI, loader, overlap, window, scaling and BAM tests. Coverage-weights tests
  also reference the fragment-count variant when its command is disabled. Discuss a package-wide
  cleanup of those test cfgs rather than adding unrelated product dependencies to outliers.
- **Shared lint failure.** The recorded clippy check failed at the existing zero-valued
  `RIGHT_PANEL_GAP > 0` comparison in `src/shared/plotters/heatmap.rs`. Review this separately if
  requiring a clean package clippy check. It is not an outliers implementation task.
- **Output publication rollback.** Retain the existing shared publisher. A possible replacement
  requires the package-wide discussion in the
  [nice-to-have plan](output_publication_rollback_plan.md).

Restore-mean statistics and normalization changes require an explicit request from the owner.
They are not scheduled work in this plan.
