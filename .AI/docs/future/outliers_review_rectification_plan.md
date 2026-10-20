# Outliers review checklist and remaining work

Updated 2026-10-20. This is the implementation checklist. The [current spec](../specs/outliers.md)
contains lasting invariants. The [validation list](outliers_real_data_validation.md) contains
scientific investigations, and the [dated review record](outliers_review_followup_plan.md) links
the original annotations and mathematical derivations.

## Owner decisions

| Item | Direction |
| --- | --- |
| Failure output | Distinct readable error filenames only. No new normal files after analysis or report-writing errors. Publication uses the existing package helper. |
| Global fitting | Skip entirely unless a chromosome needs fallback. Preserve pooled observations and alpha accounting. |
| Small-T1 fits | Apply the ridge-specific rule, preserving constrained boundary fits. Never silently change family or T1. |
| Ambiguity errors | After the second BAM pass, discover all failures across selected chromosomes. Report focused failed-region evidence without additional BAM passes or plots. |
| Feature/default status | Reserved to the owner. No changes in this update. |
| Naming | Defer combined CLI and metadata renames until the owner's manual verification. |
| Biological criteria | Report written for discussion. The owner and supervisor choose tolerances and scope. |
| GC-mode counter | Documented. Implementation deferred as low priority. |

## Implementation

- [x] Separate completed analysis from normal-file publication. Withhold weights, BEDs, histograms,
  models, PNGs, audit, summary and README on analysis or report-writing failure. Use the existing
  shared publication helper and document that a rename error can leave partial results.
- [x] Save `<prefix>.outliers.error.txt`, with stage, original errors, settings, attempt identifier
  and current companion paths. Save plain `error_histogram.tsv` for failing contexts. Include
  first-fit parameters and T1 in second-fit failures. Global rows have `NA` genomic coordinates.
- [x] Collect independent local fitting failures across chromosomes from the first-pass histograms.
  A required pooled fit failure is fatal. Complete local contexts never silently borrow that fit.
- [x] Evaluate all second-pass candidates before enforcing ambiguity error policy. Save every failed
  merged region's member evidence, selected windows, scores, cutoff and contributing exact calls
  in `error_audit.tsv`, with `withheld_error` actions. Omit unrelated successful decisions.
- [x] Check fallback need from full chromosome span and summed accepted midpoint support. Make the
  pooled fit optional and family explicit. Keep observed pooled counts. Use `global_fit_status`
  and `NA` expected counts when no global fit was needed. Adapt plotting and documentation.
- [x] Reject T1=2 for both free-zero mixtures. At ZINB T1=3 with retained counts at 1 and 2,
  reject positive-pi truncated ZIP candidates above `1e-12`. Preserve equality, numerical
  near-equality and constrained boundary cases. Add stage-specific alpha/support guidance.
- [x] Preflight optimizer vertices and strictly re-evaluate accepted candidates. Keep likelihood
  and invariant errors fatal. Report converged/unfinished starts and mean, dispersion, pi and
  threshold ranges across valid near-optimal candidates, with the cost criterion recorded.
  Do not present those ranges as confidence intervals.
- [x] Keep T2>T1 valid and add `final_threshold_exceeds_initial`. Never clamp thresholds.
- [x] Use independent close-up bins accumulated from the same optional plot reread. Show all
  visible candidates with actual decisions. Log each plot's group, count and requested span.
  Add chromosome labels and complete-histogram legend. Label fitted means by their meaning.
- [x] Align plot-multiplier and percentage parser bounds with programmatic validation. Remove
  unused hysteresis metadata, mark automatic-only fields `NA` with manual alpha, and document
  existing metadata-to-CLI mappings. Preserve existing CLI/key names for the later naming pass.
- [x] Document filters, two detection passes, optional rereads, reference/gap eligibility,
  selected-chromosome fallback, split-run alpha limitations, positional masks versus whole-fragment
  weights, qualified neighboring-background effects, default 5 kb preservation boundary and equality,
  stage-specific troubleshooting, and stale optional/error files on reused prefixes.
- [x] Remove outliers from the unused text-writer feature grouping, gate the heatmap entry point
  to its consumers, and replace the three equivalent-expression model lints.
- [x] Update changelog, main command listing, output README source and website guide. Track concise
  [historical methods/outcomes](outliers_historical_simulation_evidence.md), explicitly retaining
  the old 10% cutoff and Python fitter. Preserve annotations in a historical archive and shrink
  the old follow-up plan to a dated record.
- [x] Write the [short biological discussion report](outliers_biological_release_criteria_report.md).
  It contains proposals, not new real-BAM measurements or selected release criteria.

## Regression coverage authored

Expectations are hand-derived. These tests are compiled only. The owner runs them.

- [x] Positive-dispersion two-stage ZINB with geometric-mixture counts, T1=T2=20 and exact retained
  count `3298532786176`. Also a command-level geometric-block ZINB example with positive dispersion.
- [x] ZIP T1=2, ZINB T1=3 ridge, equality/near-equality, valid sparse boundary and T2>T1 fixtures.
- [x] Correct local-mean column and explicitly discriminate active versus unconditional targets
  in the existing nonzero-pi roundtrip example.
- [x] A core crossing 1 Mb tiles, with stride 600 kb and context 2,000,030 bp. Assert counts,
  midpoint ownership, joined call and equality with 2 Mb processing tiles.
- [x] Automatic alpha and multiplier against 900 eligible masked positions and equivalent manual
  settings. Successful reads-as-fragments detection, including per-read coverage and support.
- [x] Command-level consumer checks distinguish analysed chromosomes with no calls from unanalysed
  chromosomes, reject reference-length mismatches, allow a subset, and warn on filter differences.
- [x] Preservation, correction and error actions/reasons. Four ambiguity errors across two
  chromosomes. Empty error-prefix destinations and preserved bytes for previous normal outputs.
- [x] Multiple local zero-coverage fit failures with sparse readable histograms and no normal
  files. Blocked publication destinations are named in the terminal error and readable report.
  A separate regression keeps the original analysis error visible when saving diagnostics fails.
- [x] Actual optional plotting chunk loop with a peak split at 1 Mb, a two-base mask, and
  independently resolved close-up bins. The existing accumulator fixture checks means/peaks.
- [x] Whole-fragment shoulders for stationary 30 bp fragments and a 1 bp call at weight 0.2.
  Expected ordinary coverage is 6 at the call, then 6.8, 18 and 30 at distances 1, 15 and 30.

## Deferred work

- [Package-wide publication rollback](output_publication_rollback_plan.md) is a nice to have.
  The owner deferred it on 2026-10-20. Discuss general behavior before implementation.
- GC-mode zero-weight exclusion counters. They count tile encounters and can differ because of
  padding-only fragments, without demonstrating different counted coverage. The owner does not
  regard this reporting issue as a current focus. No GC-package fixture or counter fix is included.
- Combined CLI/metadata renames N1-N5 during manual verification. Include config serialization,
  producer/consumer filter keys, help, README, website and generated references in that pass.
- Biological release criteria and the empirical work in the validation list, including ordinary
  neighboring signal, genuine gains, moderate artifacts, GC interaction and the 0.1x-200x range.
- Feature/default/release-status decisions remain owner-controlled.
- Broader optimizer calibration and parameter uncertainty. The cheap ridge guard and start ranges
  do not prove every fit is well determined or that extreme tails are calibrated.
- GC/scaling producer integration O-002 and consumers beyond fcoverage. Move the resolver only
  when another consumer needs it. Shared counter restructuring M6 is not needed for these fixes.
- Plot concurrency, caps and broader performance work require measured need. No new limits or
  parallel BAM readers are introduced.
- Restore-mean statistics and normalization changes remain outside scope unless requested.

## Verification

Final checks on 2026-10-19:

| Check | Result |
| --- | --- |
| `cargo check --offline --locked` | Passed without warnings. |
| `cargo check --offline --locked --tests --features "testing cmd_outliers"` | Passed without warnings. Tests were compiled, not run. |
| `cargo check --offline --locked --no-default-features --features cmd_outliers` | Passed without warnings. The six focused warnings are resolved. |
| `cargo check --offline --locked --no-default-features --features cmd_fcoverage` | Passed without warnings. |
| `cargo check --offline --locked --all-features` | Passed with 22 warnings in unchanged fragment-kmer and position-visualization code. |
| `cargo fmt --all -- --check` | Passed. |
| `git diff --check` | Passed. |
| CLI reference generation with `gen_cli_docs`, features `cli docs_gen cmd_outliers`, scope `release` | Passed, 13 command pages generated under the existing ignored destination. |
| Local links in active outlier documentation | Resolved, including the original incoming Claude-review anchor. |
| `cargo clippy --offline --locked --lib --no-default-features --features cmd_outliers` | Exit 101 at the existing `RIGHT_PANEL_GAP > 0` comparison, now `heatmap.rs:477`. Reported 51 warnings. The three numerical-expression warnings and newly introduced needless borrows/type-complexity warning are resolved. Writer argument-count and module-naming style warnings remain. |

No project tests, new simulations, real-BAM analyses or benchmarks were run. The owner runs tests.
Compilation cannot establish that the numerical assertions pass. Historical clippy runs with
lib-test targets also found unrelated existing errors. This update does not claim a clean clippy run.
