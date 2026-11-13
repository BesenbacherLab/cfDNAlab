# Ideas outside the readability refactor

These came up during the refactor in [readability_plan.md](readability_plan.md) but were deliberately not done, because they change behavior, results or files outside the plan. Each is a suggestion to decide on, not a pending task.

## Possibly worth doing

- **Copy the scripts into the results.** `run.tsv` records the benchmark commit and a dirty flag, but a dirty checkout leaves the exact code unknown. Copying the job script, `commands.sh`, `measure_command.sh` and the helpers into the results directory would make each result self-describing. It adds a few small files per job.
- **Replace the SQLite residual counter.** All approved per-base commands write integer depths, so their residuals have few distinct values and a plain `Counter` would bound memory too. `ResidualCounts` in `compare_coverage.py` would shrink to a dictionary. This is the open D7 decision in `auditability_plan.md`, so it needs your call on fractional per-base tracks first.
- **Say why uncorrected `cfdnalab` is run.** `experiment_commands.md`'s "cfDNAlab without GC correction" section reads like a speed comparison. One sentence could say it supplies the uncorrected concordance against `griffin-gc`'s raw table, and that its timing has no Griffin counterpart.
- **Move the FinaleToolkit header check out of the ends job.** The 14-line Python heredoc that checks BAM/2bit contigs and the 1 Mb window guard could become a small script like `midpoints/prepare_sites.py`, and then be tested.
- **One autosome rule.** The regex `^(chr)?([1-9]|1[0-9]|2[0-2])$` appears in all four job scripts and in coverage's `save_profile`. A shared function in `common.sh` would keep them identical.

## Cosmetic or organizational

- **One worker-count column name.** Coverage and lengths use `tool_threads`, while end motifs and midpoints use `tool_workers` for the same idea. Renaming changes `timings.tsv` columns, which matters for any existing results and summaries.
- **Remove the old FrEIA profile branch.** `compare_end_motifs.py` skips `T__sample.pq` profiles saved before the left-end fix. New jobs never write them, so the branch only serves old results directories.
- **Split the comparator tests by feature.** `test_benchmarks.py` and `test_review_regressions.py` mix coverage, length, motif and midpoint comparator tests with installer and submission tests. Files such as `test_compare_coverage.py` would be easier to find.
- **Fewer copies of the run rules.** Timing, failure handling and saved results are described in `README.md`, `benchmarks_spec.md` and the feature READMEs, in slightly different words. A single description with links would drift less.
- **Tidy the plan and review files.** `codex_review.md`, `new_claude_review.md`, `benchmarks_review_plan.md`, `new_review_plan.md` and `auditability_plan.md` are historical. A `notes/` folder would leave the benchmark folder showing only what runs.
- **The midpoint status table.** `compare_midpoints.py` decides each group's status with an if/elif chain whose order matters. A short table of cases with one comment each would be easier to check.
