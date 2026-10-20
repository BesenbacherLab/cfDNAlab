# Review: `cfdna outliers` on branch `outliers2`

> **Errata (2026-10-19, after Codex verification).** These claims below are wrong or overstated.
> The corrected versions and derivations are in
> [the follow-up plan](../.AI/docs/future/outliers_review_followup_plan.md#claude-review-of-this-plan).
>
> - **1.3:** ZIP is also unidentified at T1 = 2. The proposed occupied-level rule is withdrawn. Use a
>   support-based ridge rule instead.
> - **1.7:** "T2 drops in every core" is too strong, and a manual alpha does not reproduce pooled
>   fallbacks.
> - **B1 / Summary item 5:** the formula is the expected ordinary-background loss under stated
>   assumptions, not a V-shaped trough "next to every call".
> - **D4:** the fallback pool is all selected contigs, not an autosomal model. chrM behavior depends
>   on its coverage.
> - **M3:** the clippy run exits 101 on pre-existing denied lints in unchanged shared and test files.
>   It is not just "3 trivial lints".
> - **M4:** the linked report is readable locally. It is broken only for the tracked repository.
> - **T1:** switching to column 12 fixes the label but adds no discrimination, because pi = 0 in that
>   fixture. A nonzero-pi fixture is needed.

Date: 2026-10-19. Scope: `git diff main...HEAD` (commits `5354a196` to `9e9bdac5`), about 10k added lines.
Integration into `gc-bias`, `lengths` and other commands is out of scope, as requested. Forward-looking
notes are marked as such.

## Method and evidence labels

Read in full: every file under `src/commands/outliers/`, `src/shared/outlier_weights.rs`, the
`fcoverage` diff, `cli_common.rs`/`counters.rs`/`build.rs`/`Cargo.toml` diffs, `tests/test_outliers_command.rs`,
the module-local tests, the spec, both `.AI/docs/future/outliers_*` files, the website guide and
`scripts/plot_outliers.R`. Also read the shared helpers the command relies on (`build_tiles`, `Tile`,
`Coverage`, `bam_contigs_info`, the BED blacklist loader, `find_overlapping_windows`, segment fragment
assembly, read filters) and the two external numerical dependencies (argmin 0.11 Nelder-Mead,
statrs 0.18 `checked_beta_reg`).

Commands run: `cargo check --offline --locked --tests --features "testing cmd_outliers"` (clean),
`cargo check --offline --locked --no-default-features --features cmd_outliers` (6 warnings, see M2),
`cargo check ... --features cmd_fcoverage` for comparison (0 warnings), and `cargo clippy` (see M3).
No tests were run.

- **[Verified]**: read in the code and traced, derived by hand, or observed from a command.
- **[Inference]**: reasoned from the code but not demonstrated on data. Treat as a hypothesis to check.

Severity: **High** = wrong results or lost work in realistic use. **Medium** = wrong or misleading
results in reachable cases, or a significant usability gap. **Low** = edge cases, cosmetics, hygiene.

## Summary

No crash-level bug was found in the main path. The ZIP/ZINB likelihood algebra, the truncated-fit
normalization, the NB tail mapping, the preservation heap, the tile merge and the `fcoverage` lookup
all check out (section 7 lists what was verified). The main problems are:

1. Fit failures and `--ambiguous-region-action error` delete every staged output, including the
   histograms that would explain the failure (1.1, Medium).
2. The sample-wide ZINB fit must succeed even when no core uses it, on the histogram most exposed to
   extreme pileups (1.2, Medium, partly inference).
3. ZINB accepts second fits whose mean and dispersion are not identified (T1 <= 3) and then silently
   picks the lightest-tailed member, giving the most permissive T2 (1.3, Medium for an edge case).
4. Automatic alpha changes with chromosome selection, so per-chromosome jobs give different calls
   than a genome-wide run. Undocumented (1.7, Medium).
5. Whole-fragment weighting carves V-shaped troughs next to every call, at nucleosome scale (B1).
   Documented qualitatively but not its size or shape.
6. User documentation lacks troubleshooting, runtime expectations, read-filter details, the
   effective maximum correctable call length, and short-contig behavior (section 3).
7. One integration test asserts the local-mean weight against the wrong column and cannot detect a
   regression to the unconditional mean (T1).
8. The unvalidated command is in default features and in the release CLI-docs scope (M1).

---

## 1. Correctness and behavior

### 1.1 Medium. Failed runs discard the diagnostics needed to fix them [Verified]

`run_outliers` writes `outliers.histograms.tsv.zst` into the final-output staging directory
(`src/commands/outliers/outliers.rs:297-303`) before fitting (`:305-315`). Any later error returns
before `final_outputs.move_into_place()` (`:376`), and `TempDirGuard::drop` removes the staging
directory (`src/shared/tiled_run.rs:789-798`).

Consequences:

- After a full BAM pass, a failed local fit (for example "fitted ZINB dispersion ... reached its
  numerical upper bound") leaves only the error text. The per-core histogram that would show why the
  context failed is deleted.
- With `--ambiguous-region-action error`, the check fires inside `evaluate_regions`
  (`preservation.rs:255-261`) for the first flagged merged region found. The audit and summary that
  would list all offending regions are deleted. The message gives one merged interval (with the
  chromosome in its context) and the cutoff.

Suggestion: keep the "no actionable outputs" contract but publish diagnostics. For example, on a fit
error move the histogram (and models for cores that did fit) into place under a clearly marked
name, and for `error` publish audit and summary while withholding weights and BEDs.

### 1.2 Medium. The sample-wide model must fit even when it is unused [Verified path, Inference on risk]

`fit_two_stage(&global_histogram.counts, ...)` always runs and its failure aborts the command
(`outliers.rs:306-307`). With the default `chr1..chr22` and `--bin-size 5000000`, every chromosome can
meet the context requirements, so no core uses the global model. It then only feeds
`outliers.global_fit.tsv/png`.

This pooled histogram is also the one most likely to fail. Its first (complete) fit sees every
unblacklisted extreme pileup in the genome with no T1 exclusion, which pushes the NB dispersion up
toward the `MAXIMUM_DISPERSION * 0.999` guard (`zinb.rs:228-232`) [Inference].

Suggestion: fit the global model only when a fallback is needed, or make a diagnostic-only failure
a warning with `NaN` diagnostics. If strictness is intended, say so in the help.

### 1.3 Medium (edge case). ZINB accepts unidentified second fits [Verified algebra]

`fit_two_stage` only requires `initial_threshold > 1` (`distribution.rs:143-146`).

- T1 = 2: the retained positive support is {1}, so the positive term in
  `ProfiledHistogram::evaluate` is `log(P(1) / P(1)) = 0` for every (mu, dispersion)
  (`zinb.rs:122-133`). Every pair with `P(0) / P(X < T1) <= zero_fraction` has the same cost.
- T1 = 3: one observed ratio n2/n1 constrains two parameters, leaving a curve of equal-likelihood
  solutions.

No start beats the ZIP boundary by more than `LIKELIHOOD_TOLERANCE = 1e-12` (`zinb.rs:206`), so the
Poisson member is kept. Within the tied set it has the lightest tail and therefore the lowest T2,
the most permissive calling. Nothing is reported.

ZIP is fine here. With a single positive level its constrained maximum lies well defined on the
`pi = 0` boundary. The `paired_outlier_fixture` integration test exercises exactly that path. The
retained positives are all coverage 1, and the second ZIP fit moves from (lambda 5.73, pi 0.88) to
(lambda about 0.093, pi 0) with T2 = 2 (hand-derived, matches the test's expected `initial_threshold=7`,
`second_fit_retained_positions=970`, `final_threshold=2`).

Reachability: needs `P(X >= 2)` or `P(X >= 3) <= alpha` under the first fit, so very shallow contexts
or a permissive manual `--tail-probability` (accepted up to 0.5). With automatic alpha at 0.1x,
T1 is about 7 [Inference].

The spec says "Fail explicitly if T1 leaves insufficient positive support". Suggest requiring at least
three non-empty positive retained levels for ZINB (two for ZIP), with the existing "choose a smaller
tail probability" advice.

### 1.4 Low. Nelder-Mead "convergence" is a cost-spread test only [Verified source, Inference on impact]

argmin's `NelderMead::terminate` stops when the standard deviation of the three vertex costs is below
`sd_tolerance` (`argmin-0.11.0/src/solver/neldermead/mod.rs:418-431`). Parameter spread is not
checked. On weakly identified surfaces (shallow depth, few positive levels) different starts can stop
at different points of a ridge with costs within 1e-12. The kept start is the first that beats the
running best by more than 1e-12 (`zinb.rs:206`).

T2 extrapolates beyond T1, so the ridge position matters more for T2 than for the fit quality.
Suggestion: record the spread of (mu, dispersion, T2) across converged starts in `models.tsv`. It is a
cheap stability diagnostic and directly supports validation item 4 of the review checklist.

### 1.5 Low. Optimizer error handling can panic or abort on a trial vertex [Verified]

- argmin evaluates the initial simplex with `problem.cost(p).unwrap()`
  (`neldermead/mod.rs:336`). An `Err` from `evaluate` at a start vertex panics inside a rayon worker
  instead of returning the contextual fitting error.
- During iterations, an `Err` at any trial vertex propagates through `run()` and aborts the whole
  command (`zinb.rs:187-191`), rather than rejecting that vertex.

I found no in-domain parameters that make `evaluate` fail on realistic histograms [Inference]. Option:
return `+inf` for trial-vertex evaluation errors (and count them), then keep the strict `evaluate` of
the selected parameters. That is not a silent fallback because the accepted result is still checked.

### 1.6 Low. T2 above T1 is not reported [Inference]

Positions in `[T1, T2)` are excluded from the second fit and also not called. It can happen when the
refit has a slightly heavier tail than the first fit, for example in a clean context where both fits
estimate nearly the same parameters. Suggest a `models.tsv` flag or a summary count.

### 1.7 Medium. Automatic alpha depends on the chromosome selection [Verified]

`tail_probability = multiplier / eligible_positions` over the selected chromosomes
(`outliers.rs:277-279`). A chr21-only job has roughly 40 M eligible positions against roughly 2.8 G
for chr1-22, so alpha is about 70 times larger and T2 drops in every core. Splitting by chromosome is
a common HPC pattern. The help says "across the selected chromosomes" but does not warn that split
runs are not comparable. Suggest documenting this and recommending an explicit `--tail-probability`
(copied from a genome-wide run's `# tail_probability`) for split runs.

### 1.8 Low. `fcoverage` exclusion counter depends on GC mode [Verified]

The GC-file branch skips fragments outside the fetch interval before the outlier lookup
(`src/commands/fcoverage/fcoverage.rs:1103-1107`, counter at `:1114`). The GC-tag and plain branches
look up every fetched fragment (`:1161`, `:1218`). The same input therefore gives different
`outlier_excluded_fragments` values with and without `--gc-file`. The spec calls it "tile encounters",
which stays true, but the number is not comparable across modes.

### 1.9 Low. Preserved-region plots: unbounded, sequential, coarse close-up [Verified]

- Groups and chunks are processed sequentially (`outliers.rs:443-486`, `:499-523`) with no cap on
  the number of groups and no progress bar. Plot width is 10 times the group span (minimum 100 kb), so a
  20 Mb preserved group rereads 200 Mb or the whole chromosome. Overlapping contexts are reread.
- The close-up reuses the context bins (`preserved_plotting.rs:136-143`, `:156-161`). With the default
  100 kb minimum span, bins are 50 bp. For a 30 bp call the close-up spans 230 bp, about five bins.
- Only member calls are shaded (`:214-228`). Actionable calls inside the same window, which count in
  `maximum_window_percentage`, are invisible in the figure.

Suggestions: parallelize per group, warn or cap above a number of groups, compute separate fine bins
for the close-up, and draw all candidates colored by action.

### 1.10 Low. Global-fit plot labels [Verified]

- `plotting.rs:480` labels `(1 - pi) * mu` as "Used unconditional mean". Only the active-component mean
  is used (as the `local-mean` target). The unconditional mean is a diagnostic.
- The bottom panel's x-axis is the core index across all chromosomes (`:441`). Chromosome boundaries
  are faint unlabeled lines, so a user cannot tell which chromosome a spike belongs to.

### 1.11 Low. Recorded metadata [Verified]

In `src/commands/outliers/outputs.rs`:

- `# automatic_probability_denominator` is written in manual mode too (`:183-186`).
- `# hysteresis=false` (`:150`) refers to an option that does not exist.
- `# detector=two_stage_zinb` (`:112`) duplicates `# model`.
- Keys differ from CLI flags: `model_core_size` is `--stride`, `minimum_context_span` is `--bin-size`,
  `minimum_context_fragment_support` is `--min-context-fragments`, `minimum_mapping_quality` is
  `--min-mapq` (`:192-196`, `:215`). The README does not give this mapping, so a user cannot rebuild
  the command line from a table header.

### 1.12 Low. Validation timing [Verified]

`--preserved-plot-context-multiplier` parses with `parse_positive_f64` (accepts 0.5) while
`validate()` rejects values below 1 (`config.rs:400` vs `:616-620`). `--max-outlier-percentage` has no
parser and relies on `validate()`. Both fail before BAM IO, so this is cosmetic.

---

## 2. Biological and methodological consistency

### B1. Whole-fragment weighting creates troughs next to calls [Verified math]

Background fragments that overlap a call by at least 1 bp receive the call's weight `w`. For a
background of fixed fragment length `l` and mean `mu`, a position at distance `d` (1 to `l`) outside
the call edge has a fraction `(l - d) / l` of its fragments reaching into the call. Its corrected
coverage is

```text
mu * (d / l + (1 - d / l) * w)
```

So coverage drops to about `mu * w` right at the edge and recovers linearly over one fragment length
on each side. For a strong pileup with `w = 0.05` and mono-nucleosomal fragments, that is a dip to
about 5% of background spanning roughly 167 bp on each side. This is the scale of nucleosome-depleted
regions and footprints, so it can create NDR-like features in TSS, WPS-like or footprint profiles
computed from `fcoverage`. The README and help say "Even a 1 bp overlap applies the full weight", but
not what that does to neighboring coverage. The validation list covers inspection of this. The user
guide should show the expected shape.

Forward-looking: overlap probability scales with `call length + fragment length`, so near short calls
longer fragments are downweighted more often than short ones. When `lengths` consumes these weights,
this skews local fragment length distributions toward shorter fragments around calls.

### B2. Sparse calls inside moderate broad gains can be corrected [Inference]

When a gain's mean coverage sits below T2, only its upper tail crosses T2. If those positions cover
less than 5% of the 100 kb window, the calls are corrected toward the context's background mean,
which lies below the elevated level, and all overlapping fragments are downweighted. That carves holes
into a genuine gain. Rough scale: at 30x with dispersion near 0.06, T2 is about 3.4 times the mean
(Wilson-Hilferty approximation of the gamma-Poisson tail at alpha 3.4e-10), so focal gains at about
3x relative coverage are in the risky band. This is the most likely biologically harmful failure mode
at high depth and high tumor fraction. It is listed under O-003 in the validation plan, but deserves
priority there.

### B3. Detection on raw coverage before GC correction [Inference]

Short GC-extreme elements can be locally elevated in GC-biased libraries. If called, they are
downweighted here and then GC-corrected again by `fcoverage --gc-file`. CpG islands overlap TSSs, so
B1 and B3 can combine exactly where TSS profiles are measured. The fitted dispersion absorbs part of
the GC variation and raises T2, so this may be self-limiting. Suggested check: GC content and CpG-island
overlap of calls versus matched background, at 0.5x and at 30x or higher.

### B4. Effective maximum correctable call length [Verified]

With defaults and no nearby blacklist, an isolated call's own forward window holds at least
`min(L, 100000)` called bases out of at most 100,000 eligible positions. Its score is therefore above 5%
whenever `L > 5000`. Technical pileups longer than about 5 kb are never corrected by default. This is a
valid trade-off but users should be told the number.

### B5. Read filters define which pileups exist [Verified]

`default_include_read_paired_end` drops duplicate-flagged, secondary, supplementary, QC-failed and
low-MAPQ reads (`src/shared/read.rs:5-20`). Whether PCR-duplicate pileups can be detected therefore
depends on duplicate marking. The help's coverage definition does not say this (see D3).

---

## 3. User-facing documentation (help, `outliers.README.txt`, guide)

From the perspective of a user who does not read the code.

- **D1. No troubleshooting.** The help and guide do not say what to change when a fit fails
  ("dispersion ... reached its numerical upper bound", "an unconverged ZINB start improved materially",
  "initial ZINB threshold ... leaves no positive second-fit support", "fitted ZINB mean reached its
  numerical upper bound"). Options to document: larger `--bin-size` or `--min-context-fragments`,
  blacklisting the region, `--model zip`, a smaller `--tail-probability`. The ZIP truncation errors
  still lack the tail-probability advice (known O-007).
- **D2. Runtime.** The command reads the BAM twice, plus once more around preserved groups with
  `--plot-preserved-regions`. Not stated in the help or guide.
- **D3. Read filters.** Add one sentence that duplicate-flagged, secondary, supplementary and QC-failed
  reads are excluded, so PCR-duplicate pileups are only visible in unmarked BAMs.
- **D4. Short contigs and `--chromosomes all`.** Chromosomes shorter than `--bin-size` use the
  sample-wide autosomal model. chrM (16,569 bp) is then almost entirely called and preserved through the
  50 kb denominator floor (16,569 / 50,000 = 33%). Unplaced, alt or decoy contigs can be corrected
  against an autosomal background. Advise against `all` or recommend blacklisting such contigs.
- **D5. Reference gaps.** Unblacklisted N-gaps count as eligible zero-coverage positions. They enter the
  automatic alpha denominator and the preservation window denominators (lowering nearby scores).
  Recommend including gap regions in the blacklist.
- **D6. Maximum correctable length.** State the about 5 kb isolated-call limit from B4 next to
  `--max-outlier-percentage`.
- **D7. Split runs.** Document 1.7.
- **D8. "Use a blacklist instead" in the guide.** The example passes `outliers.flanked.bed` to
  `fcoverage --blacklist`, then says "The receiving command's blacklist rules determine what is
  excluded". For `fcoverage` the rule is known: positions are masked from outputs. Say it, and say that
  `exact.bed` therefore leaves the pileup fragments' shoulders in place while `flanked.bed` masks them.
  The BED loader skips `#` lines (`src/shared/bed.rs:81-83`), so the metadata headers are fine.
- **D9. Trough shape.** Add the B1 formula or a figure to the guide's "Check the effect on coverage".
- **D10. Stale plots.** PNGs from an earlier run with the same prefix remain in the directory when a
  later run has fewer groups or no plots. The table lists the current files only.
- **D11. Metadata mapping.** See 1.11.
- **D12. `global_fit.png`.** The README does not explain that the bottom panel's x-axis is a core
  index. The second panel ("Complete original histogram") draws two fitted curves and eight vertical
  threshold lines with no legend (`plotting.rs:319-414` never calls `configure_series_labels`), so its
  lines can only be decoded by comparing colors with the first panel.

Claims checked against code and found accurate: window and backward-fallback rules, strict cutoff
equality, flank merging, target definitions, histogram format including the fully blacklisted row,
summary columns and `NaN` rule, preserved-group joining (overlap, not touching), and the `fcoverage`
overlap rule.

---

## 4. Naming and CLI surface

Relevant to checklist item 5 ("Finalize CLI names").

- **N1.** "protection" (`--protection-window-size`, `--min-protection-positions`) versus "preserve" and
  "preservation" everywhere else, including the help heading.
- **N2.** `--bin-size` means the minimum context span. "Bin" means output bins elsewhere in cfDNAlab
  (`fcoverage --by-size`) and plot bins in this command. The metadata already says
  `minimum_context_span`.
- **N3.** `--blacklist-flank` sits under the "Outlier Model" heading but controls BED output and
  preservation merging. `--target` is a correction setting under the same heading.
- **N4.** "Outlier" already names a different concept in `gc-bias`: `--outlier-method`,
  `--outlier-scope`, `--outlier-quantiles`, `--outlier-k` (`src/commands/gc_bias/config.rs:295-338`).
  When `gc-bias` gains `--outlier-weights`, one command will expose two unrelated "outlier" concepts.
- **N5.** `--max-outlier-percentage` is the preservation cutoff. `--ambiguous-region-action` uses
  "ambiguous" while the audit uses "flagged" and `density_above_cutoff`.

---

## 5. Tests: gaps and weak assertions

- **T1. Weak assertion [Verified].** `local_mean_and_zero_targets_write_their_regional_mass_weights`
  compares the local-mean weight with column 14, `underlying_mean`
  (`tests/test_outliers_command.rs:1014`, `:1023`). The target uses `underlying_component_mean`
  (column 12, `outliers.rs:808`). In this fixture the second fit has `pi = 0` (1.3), so both columns are
  equal and the test cannot catch a regression to the unconditional mean. The comment ("the same
  underlying fitted mean") repeats the mix-up. The fcoverage roundtrip test does distinguish 2.82 from
  2.70, so coverage exists elsewhere, but this test should use column 12.
- **T2.** No end-to-end ZINB run with positive dispersion. The only ZINB command test uses an
  underdispersed fixture that selects the Poisson boundary. No unit test covers `fit_two_stage` for ZINB
  (first fit, T1, truncated refit, T2), and none covers the failure guards in `fit_zinb` (dispersion
  bound, mean bound, unconverged better start) or the ZINB `initial_threshold > 1` error.
- **T3.** No test with a model core straddling two processing tiles. Every fixture either uses a
  stride that divides the 1 Mb tile (250, 1,000, 5,000, 250,000) or has a single tile (stride 300 on a
  300 bp chromosome). The partial-core histogram merge and split fragment support
  (`outliers.rs:542-565`, `:720-751`) are untested. Because `build_tiles` rounds the tile down to a
  multiple of the stride whenever at least 10 strides fit, a straddling core needs a stride above
  tile/10 that does not divide the tile (for example a 1.5 Mb stride with the 10 Mb default tile) or a
  stride larger than the tile.
- **T4.** No command-level test for `--ambiguous-region-action error` (nothing published) or `correct`
  (flagged calls written to the weights with `flagged_region_correction_requested`).
- **T5.** No command-level assertion of the automatic alpha. The roundtrip test relies on
  `1 / 10,000` but never checks `# tail_probability`. The multiplier is untested end to end.
- **T6.** No `--reads-are-fragments` run of `outliers`.
- **T7.** `collect_preserved_plot_coverage`'s chunk loop is untested across a `--tile-size` boundary.
  The bins unit test simulates chunks directly.
- **T8.** No command-level check that a preserved call appears in the audit with `action=preserve` and
  the right `decision_reason`.

---

## 6. Build, repository and merge readiness

- **M1.** `cmd_outliers` is in the default features (`Cargo.toml:42`), in the release command list of
  `examples/gen_cli_docs.rs:161-162`, in `website/scripts/generate_cli_docs.sh`, and the guide is in the
  sidebar. The command is explicitly unvalidated. Decide whether to gate it or label it experimental in
  the help before any release.
- **M2. [Verified]** `--no-default-features --features cmd_outliers` builds with 6 dead-code warnings
  (`TextWriter` family in `src/shared/io.rs` under `writes_text_outputs`, `write_heatmap` under the
  `plotters` feature). The `cmd_fcoverage`-only build has none. Cause: `cmd_outliers` was added to the
  `writes_text_outputs` group (`build.rs:76-85`, which does not list `cmd_fcoverage`) and enables
  `plotters`, but uses neither item.
- **M3. [Verified]** clippy: `manual_range_contains` at `model.rs:200` and `:268`,
  `manual_saturating_arithmetic` at `model.rs:501-503`.
- **M4. [Verified]** The checklist and validation docs link
  `.AI/docs/future/outliers_fragment_simulations/zinb_comparison_report.md`, which is gitignored
  (`.gitignore:56`). The link is broken for anyone else, including Codex.
- **M5.** No `CHANGELOG.md` "Unreleased" entry for `cfdna outliers` or `fcoverage --outlier-weights`,
  and the README command list is unchanged.
- **M6.** `OutliersRunResult.counters` is the public `FCoverageCounters`, so the outliers result exposes
  `outlier_excluded_fragments`, `gc_failed_fragments` and `tile_owned_normalization_*`, which are always
  zero for this command. Together with `add_fragment_clipped_to_core`, this makes `cmd_outliers` depend
  on `cmd_fcoverage` internals.

---

## 7. Checked and found consistent

Listed so the verification pass can target the same derivations.

- **ZIP MLE** (`model.rs:169-214`): lambda from the zero-truncated positive mean, pi from the zero
  fraction, Poisson boundary `lambda = S/N` when pi would be negative.
- **Right-truncated ZIP** (`model.rs:222-291`): lambda from `E[X | 1 <= X < T1]`. Solving
  `q0 = (pi + (1 - pi) e0) / (pi + (1 - pi) F)` gives `pi = (q0 F - e0) / ((1 - e0) - q0 (1 - F))`, as
  coded. The boundary uses `E[X | X < T1]`.
- **ZINB profile** (`zinb.rs:76-134`): the likelihood factorizes into `q0` and the positive conditional.
  The interior `pi = (zf R - P0) / (1 - P0 - zf (1 - R))` matches the code's numerator and
  denominator. The boundary condition is `zf <= P0 / R`. Both branches compute the same per-position
  objective, so comparisons with the ZIP boundary are valid.
- **NB2 log mass** (`zinb.rs:245-267`): `log P(0) = -log(1 + dispersion * mu) / dispersion`. The Stirling
  branch equals `r (log1p(x) - x) + (k - 1/2) log1p(x) + c(r + k) - c(r)` with `x = k / r`. The Stirling
  coefficients (1/12, -1/360, 1/1260, -1/1680, 1/1188) and the `log1p(x) - x` series are correct.
- **NB tail** (`distribution.rs:75-100`): `P(X >= k) = I_x(k, 1/dispersion)` with
  `x = dispersion * mu / (1 + dispersion * mu)`, matching scipy's parameterization. statrs returns
  `bt * h / a` without subtracting from 1 when `k` is in the upper tail. For near-Poisson dispersion
  (shape up to 4.9e8) the log-gamma cancellation costs about 1e-5 relative accuracy, which is harmless for
  thresholds. The 140-iteration continued-fraction cap (silent) is not reached for tail evaluations at
  cfDNA depths [Inference].
- **Positive variance** (`distribution.rs:231-235`): `m+ (1 + dispersion * mu + mu - m+)`.
- **Threshold searches** (`distribution.rs:171-201`, `model.rs:449-526`) return the smallest `k >= 1`
  with `S(k) <= alpha` under monotone survival.
- **Preservation heap** (`preservation.rs:166-215`): lazy deletion is valid because a window that ends
  at or before one call start ends before every later call start. Key order is supported first, then
  score, then the earliest window.
- **Preserved groups** (`preserved_regions.rs`): groups cannot interleave. A group's span runs from its
  first member's window start (before that member's end) to its last member's window end (after that
  member's start), so it contains any call in between, whose window then overlaps and joins.
- **Tile order**: `par_iter().collect()` preserves the order of `build_tiles` (chromosome, then position),
  which `merge_calling_pass_results` relies on.
- **Fragment support** is counted once per fragment (midpoint inside the tile core), excludes
  blacklisted midpoints, and does not depend on tile size.
- **Coverage exactness**: `Coverage` accumulates `f64` deltas and casts to `f32`, so integers are exact
  below 2^24, and the histogram pass rejects larger eligible values.
- **Contigs and metadata**: `bam_contigs_info` returns only the selected chromosomes, so
  `analysed_chromosomes` is correct.
- **`fcoverage` integration**: one `process_tile` path serves all tile modes. Zero weights skip before
  normalization and GC accounting. `1 / (max_fragment_length + 1)` admits any 1 bp overlap for accepted
  fragments. The overlap index is the scan index, which matches `keep_weights` order. The look-back of
  `max_fragment_length` covers the pairing stash's start-order slack. The fragment iterator settings
  (gap trigger 1, `!ignore_gap`) match between detection and counting.

---

## 8. Suggested order of work

1. Fix T1 and add T2-T4 tests. They protect the default model and the paths most likely to break
   during the coming integrations.
2. Publish diagnostics on failure (1.1) and decide on the global fit (1.2).
3. Add the ZINB identifiability guard (1.3) and the start-spread diagnostic (1.4).
4. Write the documentation items D1-D9, especially troubleshooting, split runs, the about 5 kb limit and
   the trough shape.
5. Run the real-data checks for B1-B3 before deciding merge readiness. B2 is the main biological risk.
6. Settle names (N1-N5), then M1, M2, M4 and M5.
