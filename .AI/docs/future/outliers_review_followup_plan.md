# Outliers review record, 2026-10-19

Reviewed Claude's findings against `9e9bdac5dedb16103817c1b88fdc986004fcae0c`, then re-read
its annotations and applied the owner's decisions. This file is a dated record, not another
implementation plan.

- [Implementation checklist and deferred work](outliers_review_rectification_plan.md).
- [Scientific validation list](outliers_real_data_validation.md).
- [Short supervisor discussion report](outliers_biological_release_criteria_report.md).
- [Historical simulation methods and selected outcomes](outliers_historical_simulation_evidence.md).
- [Original annotated verification and mathematical derivations](outliers_review_annotated_archive_2026-10-19.md).

## Decisions applied

1. Failed runs may publish readable diagnostics only under distinct error filenames. Normal
   success filenames are withheld, including histograms, models, PNGs and README. Earlier
   successful files stay intact after analysis failures.
2. Skip the pooled model fit unless some chromosome needs fallback. Keep the pooled histogram.
3. Use the ridge-specific rule, retaining valid constrained boundary fits. T1=2 fails in ZIP
   and ZINB. At T1=3, the ZINB interior dispersion/zero-inflation ridge fails explicitly.
4. After the second BAM pass, collect every ambiguity-policy failure across selected chromosomes.
   Save focused failed-region evidence and contributing calls with `withheld_error` actions.
   No extra BAM pass, success summary or plots are needed to explain those failures.
5. Feature/default status stays with the owner. No feature or experimental-label change.
6. CLI and metadata renaming remains a single later task during the owner's manual verification.
7. The requested biological report is written. It proposes measurements and discussion questions,
   and does not claim new real-BAM results or settled release thresholds.
8. GC-mode exclusion counters are documented and deferred as low priority.

Update 2026-10-20: the owner deferred publication rollback to a
[package-wide nice-to-have plan](output_publication_rollback_plan.md). Outliers uses the existing
shared publication helper. A final rename error can leave some completed outputs in place.

## Corrections retained from the exchange

Occupied observed levels are not a valid general identification rule. Empty levels still enter
through normalization. A smaller alpha changes retained support but does not repair initial fitting.
T2 can exceed T1 and must not be clamped. Shared manual alpha does not reproduce a different pooled
fallback model. The pool contains all selected chromosomes, not only autosomes. Whole-fragment
weights affect ordinary neighboring coverage, but do not predict a literal V for every observed
track. Sparse calls inside a genuine gain can still be corrected. GC calibration contamination
remains a hypothesis requiring measurement.

## Claude review of this plan

Claude accepted these corrections. Its [original annotations and responses](outliers_review_annotated_archive_2026-10-19.md#claude-review-of-this-plan),
fixture derivations, open scientific questions and build evidence are preserved in the annotated
archive rather than duplicated here. This heading keeps the original review's incoming link usable.

## Verification boundary

Project tests are authored and compiled only. The owner runs them. Compile checks cannot establish
that numerical or command assertions pass. No new BAM analysis, simulation, benchmark or biological
validation was performed in this implementation task. The current checklist records final compile
commands and known unrelated clippy failures.
