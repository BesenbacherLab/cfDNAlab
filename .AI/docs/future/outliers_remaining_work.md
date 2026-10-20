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
| `--bin-size` | `--min-context-span` |
| `--ambiguous-region-action` | `--preservation-action` |
| `--max-outlier-percentage` | `--max-called-percentage` |

Review the help headings for `--target` and `--blacklist-flank` at the same time. Keep `--stride`
unless there is a concrete reason to rename it. Apply approved names together across config
serialization, metadata, producer/consumer checks, help, the generated README, guides and generated
references. Preserve the distinction between GC-estimator outlier controls and fragment keep weights
when producer integration is considered.

Feature/default/release-status changes remain the owner's decision. The owner and supervisor also
need to choose biological retention criteria, supported event sizes, sample types and depths.
The biological report proposes discussion questions rather than selecting those criteria.

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
- **Additional consumers and calibration producers.** Consider GC/scaling producer integration and
  consumers beyond fcoverage separately. Estimate and apply corrections from consistently weighted
  fragments. The validation plan records the proposed comparisons. Move overlap-resolver state
  out of fcoverage only when another consumer needs it. Shared counter restructuring needs a
  separate reason.
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
