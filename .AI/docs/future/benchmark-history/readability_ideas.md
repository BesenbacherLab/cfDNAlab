# Ideas outside the readability refactor

These came up during the refactor in [readability_plan.md](readability_plan.md) but were deliberately not done. Each is a suggestion to decide on, not a pending task.

- **Replace the SQLite residual counter.** All approved per-base commands write integer depths, so their residuals have few distinct values and a plain `Counter` would bound memory too. `ResidualCounts` in `compare_coverage.py` would shrink to a dictionary. This is the open D7 decision in `auditability_plan.md`, so it needs your call on fractional per-base tracks first.
- **Move the FinaleToolkit header check out of the ends job.** The 14-line Python heredoc that checks BAM/2bit contigs and the 1 Mb window guard could become a small script like `midpoints/prepare_sites.py`, and then be tested.
- **Tidy the plan and review files.** `codex_review.md`, `new_claude_review.md`, `benchmarks_review_plan.md`, `new_review_plan.md` and `auditability_plan.md` are historical. A `notes/` folder would leave the benchmark folder showing only what runs.
