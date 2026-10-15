# Outliers Command Review

Review date: 2026-09-29

Branch: `outliers2` compared with `main` (commits `aeb4af30`, `5354a196`, `f5a75211`).

This review covers the new `cfdna outliers` command and the first downstream consumer of its keep
weights, `cfdna fcoverage --outlier-weights`. The intended workflow is to learn which regions to
downweight with `outliers` and then apply those fragment weights during counting in every feature
command. Only `fcoverage` consumes the weights so far.

I reviewed:

- `src/commands/outliers/` (config, model, outliers, outputs, plotting, striding, and their sibling
  test files)
- `src/shared/outlier_weights.rs`, `src/shared/outlier_weights_tests.rs`, and
  `src/testing/outlier_weights.rs`
- The `fcoverage` changes in `src/commands/fcoverage/config.rs`,
  `src/commands/fcoverage/fcoverage.rs`, and `src/commands/fcoverage/fcoverage_tests.rs`
- `OutlierWeightsArgs` in `src/commands/cli_common.rs`, plus `src/cli_command/helpers.rs`,
  `src/cli_app.rs`, and `src/lib.rs`
- Shared helpers used directly by the new code: `src/shared/overlaps.rs`,
  `src/shared/coverage.rs`, `src/shared/tiled_run.rs`, `src/shared/blacklist/load.rs`,
  `src/shared/bed.rs`, `src/shared/fragment_iterators/core.rs`,
  `src/shared/fragment_iterators/segments.rs`, and `src/shared/scale_genome.rs`
- `src/commands/coverage_weights/coverage_weights.rs` and
  `src/commands/coverage_weights/striding.rs`, for the interaction between scaling factors and
  outlier weights
- `tests/test_outliers_command.rs`, the new outlier tests in `tests/test_fcoverage_command.rs`, and
  `src/cli_app_roundtrip_tests.rs`
- `.AI/docs/specs/outliers.md`, `.AI/docs/specs/fcoverage_spec.md`,
  `.AI/docs/future/fragment_nuisance_outliers_plan.md`,
  `.AI/docs/future/fragment_nuisance_outliers_review.md`, and
  `website/docs/guides/handle_outliers_guide.md`

I did not run tests, per project instruction. I ran:

- `cargo fmt --all -- --check`: clean
- `cargo check --tests --features "testing cmd_outliers"`: clean
- `cargo check --no-default-features --features "cli cmd_outliers"`: compiles, with 6 dead-code
  warnings in `src/shared/io.rs` and `src/shared/plotters/heatmap.rs`, which this branch does not
  touch

Notation: `pi` is the fitted zero inflation, `lambda` is the Poisson-component mean, and `T1` and
`T2` are the initial and final thresholds.

## Release Triage

No mechanical defect was found in the `fcoverage` weight wiring (overlap lookup, lookback, minimum
rule, weight product, cleanup floor) or in the ZIP algebra. The active findings are scientific
defaults, pipeline consistency, silent behavior, stale documentation, and rollout design.

Pre-release correctness and safety:

- O-001: `--target local-mean` dilutes the target by the structural-zero fraction.
- O-002: Scaling factors and GC corrections are derived without outlier weights, so combining them
  with `--outlier-weights` double-corrects.
- O-003: High-level focal amplifications are called and flattened, and long regions are not
  reported.
- O-004: The keep-weight loader ignores provenance metadata. Chromosomes that `outliers` never
  analyzed silently keep weight `1.0`.
- O-005: Zero-weight and downweighted fragments are invisible in `fcoverage` statistics and are
  still reported as counted.
- O-006: Local-fit errors do not identify the failing chromosome or core.
- O-007: Legal manual `--tail-probability` values abort the run on shallow data.
- O-008: A silent `unwrap_or(0)` hides an impossible fragment-support index miss.
- O-014: Remaining high-value coverage gaps.

Pre-release docs and API polish:

- O-009: Several docs still say that no command consumes the keep weights.
- O-010: Inaccurate or missing guide and spec statements.
- O-011: Weight resolution is local to `fcoverage`, which invites drift once other commands are
  wired.
- O-012: Naming collisions in CLI options and metadata keys.

Post-release performance optimization:

- O-013: Quadratic log-factorial work inside every ZIP fit.

## Findings

### O-001 - High - `local-mean` dilutes the target by the structural-zero fraction

Evidence:

- `src/commands/outliers/outliers.rs:635` sets the per-position target to
  `model.zip.underlying_fit.mean()`.
- `src/commands/outliers/model.rs:69-70` defines `mean()` as `(1 - zero_inflation) * lambda`.
- A called position has coverage `>= T2 >= 1`. Under the fitted ZIP it therefore belongs to the
  Poisson component, and its expected coverage is `lambda`, not `(1 - pi) * lambda`.
- `tests/test_outliers_command.rs:519-523` and `:575` show the effect on a small scale. The fitted
  `lambda` is about 2.82, `pi` is about 0.045 (from the uncovered chromosome ends), and the asserted
  target is `(1 - 0.045) * 2.82`, about 2.695.

Impact:

At ordinary WGS depth, Poisson zeros are negligible (`lambda = 17` gives `P(X = 0)` of about 4e-8).
`pi` is then essentially the unmappable, N, or deleted fraction of the fitting context, which has
nothing to do with the called position. A corrected pileup is pushed below its mappable neighbors by
exactly that fraction. A context that is half assembly gap gives a target of about `lambda / 2`.
Contexts next to centromeres and assembly gaps are where satellite pileups concentrate, so the
strongest over-correction lands where most calls are expected.

Recommendation:

- Use `lambda` as the `local-mean` target.
- Update the definitions in `.AI/docs/future/fragment_nuisance_outliers_plan.md`,
  `.AI/docs/specs/outliers.md`, the `OutlierTarget::LocalMean` and `--target` help in
  `src/commands/outliers/config.rs`, and the guide.
- Keep reporting `underlying_mean` in the model TSV as a diagnostic, but do not use it as the
  correction target.

### O-002 - High - Corrections derived from coverage do not see outlier weights

Evidence:

- `src/commands/coverage_weights/coverage_weights.rs:425-465` builds the internal `FCoverageConfig`
  for `coverage-weights` without outlier weights. `ScalingWeightsArgs` has no outlier-weight field.
- `src/commands/coverage_weights/coverage_weights.rs:99-105` uses per-bin `Average` coverage for
  `coverage-weights` and `Total` mass for `fragment-count-weights`. Neither is robust to a pileup.
- `src/commands/coverage_weights/striding.rs:226-310` sets each scaling factor to
  `global_mean / smoothed_bin_value`.
- `gc-bias` has no outlier-weight input either.
- `tests/test_fcoverage_command.rs:1090` tests `--outlier-weights` together with `--scaling-factors`
  and a GC package as a plain product.

Impact:

Scaling factors built from raw coverage are already shrunk by any pileup in their bin.
`fcoverage --outlier-weights W --scaling-factors S` then removes the pileup twice: once through `W`,
and again through the too-small factor in `S`, which applies to the whole bin and, after smoothing,
to its neighbors. For example, a 5,000x pileup over 2 kb in a 500 kb bin at mean coverage 17 raises
that bin's average from 17 to about 37 before smoothing. The same issue applies to GC correction
fitted from counts that include the pileup. `src/commands/fcoverage/config.rs:94` describes the
outlier weights and the other corrections as independent. That only holds if the other corrections
were estimated from outlier-weighted data.

Recommendation:

- Wire `--outlier-weights` into `gc-bias`, `coverage-weights`, and `fragment-count-weights` before
  or together with the feature commands, as in stage 4 of the plan.
- Record the weights file (at minimum its path and a digest) in their output metadata.
- In `fcoverage` and later consumers, warn or fail when `--outlier-weights` is combined with scaling
  factors or GC corrections whose metadata shows they were built without outlier weights.

### O-003 - Medium - High-level focal amplifications are called and flattened

Evidence:

Worked example with the defaults (5.5 Mb contexts, automatic tail probability about 3.5e-10 for
chr1-chr22, ordinary coverage `lambda = 17`):

- A 1 Mb amplicon at about 93x (5.5 x 17) inside a 5.5 Mb context gives a first-fit positive mean
  of about `(4.5 * 17 + 1 * 93) / 5.5`, about 31, and `T1` about 72. That excludes the amplicon but
  keeps all ordinary positions.
- The second fit recovers `lambda` of about 17, so `T2` is about 49. The whole amplicon is called
  with a keep weight of about 0.18.
- A 3 Mb amplicon at the same level raises the first-fit mean to about 58 and `T1` to about 110. It
  stays in the second fit and is not called. Whether an amplicon is flattened is therefore not
  monotone in its size.
- Amplicons whose coverage sits near `T2` should receive many scattered partial calls.
- `.AI/docs/future/fragment_nuisance_outliers_plan.md` requires very long called regions to be easy
  to audit. The run statistics only report the number of regions
  (`src/commands/outliers/outliers.rs:340`).

Impact:

MYCN-, ERBB2-, or EGFR-type amplicons in high tumor fraction plasma are a common cfDNA signal.
`fcoverage --outlier-weights`, or blacklisting with `outliers.flanked.bed`, would erase them without
any warning.

Recommendation:

- Report the total called bases, the longest region, and the number of regions above a length
  threshold in the run statistics and in the keep-weight metadata.
- Warn when regions exceed a length threshold.
- Consider a `--max-region-length` option that keeps longer calls in the BED outputs but excludes
  them from the keep weights.

### O-004 - Medium - The keep-weight loader ignores provenance metadata

Evidence:

- `src/shared/outlier_weights.rs:68` skips every `#` line without parsing it.
- `src/shared/outlier_weights.rs:40` documents missing selected chromosomes as implicit weight
  `1.0`.
- `outliers` defaults to chr1-chr22 (`ChromosomeArgs` in `src/commands/cli_common.rs`), and every
  output records `# selected_chromosomes`.
- The loader never compares `ignore_inter_mate_gap`, `reads_are_fragments`,
  `fragment_length_range`, `minimum_mapping_quality`, or `require_proper_pair` with the consuming
  run.
- The scaling loader already parses its metadata and fails or warns on mismatches:
  `src/shared/scale_genome.rs:674`, `:744`, and `:790`.

Impact:

A sparse file cannot tell "no outliers here" from "not analyzed". `fcoverage` on chrX with a default
`outliers` run silently leaves chrX uncorrected. The regional mass relationship behind each keep
weight, where weighting every overlapping fragment brings the region to its target mass, only holds
if both runs use the same fragment filters.

Recommendation:

- Parse `selected_chromosomes` and fail, or at least warn, when a consuming command processes a
  chromosome that `outliers` did not analyze.
- Warn on fragment filter and gap handling mismatches, following the scaling-factor pattern.
- Keep accepting files without metadata, with a warning that compatibility could not be checked.

### O-005 - Medium - Zero-weight and downweighted fragments are invisible in statistics

Evidence:

- `src/commands/fcoverage/fcoverage.rs:1879-1887` marks a fragment as counted after
  `add_fragment_weighted`, including when the weight is `0.0`.
- `src/commands/fcoverage/fcoverage.rs:1102-1119`, `:1152-1169`, and `:1179-1196` then increment
  `counted_fragments`. With length normalization they also increment
  `tile_owned_normalization_fragments` and `tile_owned_normalization_length_sum`, which set the
  restore-mean multiplier.
- `src/commands/fcoverage/fcoverage.rs:165-192` reports a GC summary but has no outlier summary.
- `tests/test_fcoverage_command.rs:1028` checks only that the output is empty for a zero-weight
  fragment. It does not check counters.

Impact:

Users cannot see how many fragments were downweighted or excluded. Excluded fragments are reported
as counted and influence the restore-mean fragment length.

Recommendation:

- Add counters for fragments with a keep weight in `(0, 1)` and fragments with keep weight `0`, and
  print them in the run statistics.
- Decide explicitly whether weight-0 fragments count as counted. Skipping them before insertion is
  the simplest consistent choice.
- Decide whether the restore-mean length should be weighted by the keep weight.

### O-006 - Low - Local-fit errors do not identify the failing context

Evidence:

- `src/commands/outliers/striding.rs:180-186` propagates errors from `fit_two_stage_zip` and
  `diagnose_zip_fit` without context.
- `src/commands/outliers/outliers.rs:281-288` calls `fit_core_models` without context.
- A single failing complete local context aborts the whole run by design.

Impact:

A user who hits a fit failure sees messages such as "cannot fit ZIP because every eligible position
has zero coverage", but no chromosome, core, or context interval to investigate.

Recommendation:

Add `.with_context(...)` with the chromosome, core interval, and context interval at the per-core
fit and diagnostic calls.

### O-007 - Low - Legal manual `--tail-probability` values abort shallow runs

Evidence:

- `--tail-probability` accepts values up to 0.5 (`parse_probability` and `validate` in
  `src/commands/outliers/config.rs`).
- At about 0.1x (`lambda` about 0.17), `--tail-probability 0.2` gives a fitted `P(X >= 1)` of about
  0.15, which is `<= 0.2`, so `T1 = 1`.
- `src/commands/outliers/model.rs:326-330` then fails the global fit with "initial ZIP threshold 1
  leaves no positive coverage values for the second fit".

Impact:

The run fails after the whole first BAM pass, with a message that does not tell the user that the
probability is too large for the sequencing depth.

Recommendation:

Detect `T1 = 1` in `fit_two_stage_zip` and report it as a user-facing error that names the tail
probability and explains that it is too permissive for the coverage level.

### O-008 - Low - Silent zero fallback for fragment support

Evidence:

`src/commands/outliers/outliers.rs:398-402` reads fragment support with
`.get(support_offset).copied().unwrap_or(0)`. By construction the index is always in range during
the histogram pass.

Impact:

If a future change breaks that invariant, cores would silently get zero fragment support. That
changes context expansion and fallback decisions without any error.

Recommendation:

Replace `unwrap_or(0)` with an error that names the model-core index and the tile.

### O-009 - Low - Documentation still says no command consumes the keep weights

Evidence:

`fcoverage --outlier-weights` now consumes the file, but these statements still say that no command
does:

- `src/commands/outliers/config.rs:98`, `:109`, and `:138`
- `src/commands/outliers/outliers.rs:58` and `:143-145`
- `website/docs/guides/handle_outliers_guide.md:54` and `:85`
- `.AI/docs/specs/outliers.md:4`, `:147`, `:164`, `:181`, and `:302`

Recommendation:

Name `fcoverage --outlier-weights` as the current consumer and keep the other commands as planned.

### O-010 - Low - Inaccurate or missing guide and spec statements

Evidence:

- `website/docs/guides/handle_outliers_guide.md:136` says that `eligible_positions` is the number
  of unmasked positions used to estimate the model. For `global_fallback` rows it is the count for
  the chromosome context, while the model came from the global histogram. The diagnostics in those
  rows compare the chromosome's data with the global model.
- `website/docs/guides/handle_outliers_guide.md:138` omits that `fragment_support` only counts
  fragments with an unblacklisted midpoint.
- `.AI/docs/specs/fcoverage_spec.md:16-21` (Weighting) does not mention outlier weights, how they
  are resolved, or how they combine with the other weights.
- `src/commands/fcoverage/config.rs:85` says that outlier weighting retains a fraction of the
  fragment support. That is not true for files written with `--target zero`.

Recommendation:

Correct the guide wording for fallback rows and fragment support, add an outlier-weight bullet to
the fcoverage spec, and qualify the `fcoverage` help text for zero weights.

### O-011 - Low - Weight resolution is local to `fcoverage`

Evidence:

- `src/commands/fcoverage/fcoverage.rs:1043-1062` implements the overlap pointer, the
  `max_fragment_length` lookback, and the `1 / (max_fragment_length + 1)` any-overlap fraction as a
  closure inside `process_tile`.
- `src/shared/mod.rs:43-44` compiles `shared::outlier_weights` only with `cmd_fcoverage`.
- The plan requires all participating commands to use the same fragment interval lookup and
  minimum-over-overlaps rule.

Impact:

Wiring the remaining commands (ends, midpoints, lengths, wps, wps-peaks, fragment-kmers,
transitions, bam-to-frag, bam-to-bam, gc-bias, coverage-weights, and fragment-count-weights) would
copy this closure about ten times, and the copies can drift. The any-overlap fraction is correct
today but depends on the fragment length filter matching the lookup interval.

Recommendation:

- Move the lookup into `src/shared/outlier_weights.rs` as a small cursor type that owns the pointer
  and lookback and returns the keep weight for a fragment span.
- Replace the fraction threshold with an explicit any-overlap query, or document it next to the
  cursor.
- Widen the feature gate to every consuming command.

### O-012 - Low - Naming collisions in CLI options and metadata keys

Evidence:

- `gc-bias` already has `--outlier-method`, `--outlier-scope`, `--outlier-quantiles`, and
  `--outlier-k` for correction-matrix outliers (`src/commands/gc_bias/config.rs:307` and the
  surrounding fields).
- `src/commands/outliers/outputs.rs:138` writes `# ignore_inter_mate_gap=`, while scaling files use
  the `ignore_gap` metadata key.

Impact:

Once `gc-bias` accepts `--outlier-weights`, two unrelated meanings of "outlier" share a help screen.
A future metadata check (O-004) would need to translate between the two gap keys.

Recommendation:

Decide the `gc-bias` naming before wiring it, for example by using "coverage outlier weights" in help
text and headings. Align the gap metadata key with `ignore_gap`.

### O-013 - Low - Quadratic log-factorial work in every ZIP fit

Evidence:

- `src/commands/outliers/model.rs:673` computes `ln(value!)` by summing logs for every value below
  256.
- `src/commands/outliers/model.rs:501-541` calls it for every coverage value in the range, twice, in
  `truncated_poisson_mean` and in `poisson_log_cdf`. Each evaluation is therefore quadratic in `T1`.
- `fit_lambda_from_conditional_mean` evaluates the truncated mean about 100 times per fit
  (`LAMBDA_BISECTION_ITERATIONS = 96` at `src/commands/outliers/model.rs:49`, plus the bracket
  search), and every model core runs its own fit.

Impact:

This does not change results. It adds avoidable CPU time that grows with sequencing depth, because
deeper samples have larger `T1`.

Recommendation:

Use a precomputed `ln(k!)` table, or a running sum inside the loops. Optionally stop the bisection
once the bracket stops shrinking.

### O-014 - Low - Remaining high-value coverage gaps

- The GC-tag branch of `fcoverage` with outlier weights is not exercised. Only the no-GC and
  GC-file branches are tested.
- No test covers the lookback when fragments arrive out of start order. With
  `--max-fragment-length 200`, take fragment A = `[100, 300)` with reads `[100, 130)` and
  `[270, 300)`, fragment B = `[240, 280)` with reads `[240, 260)` and `[260, 280)`, and an outlier
  interval `[110, 120)`. The iterator returns B first because its later mate starts at 260. With a
  lookback of 0, A would miss the interval.
- No `fcoverage` test places an outlier interval or its overlapping fragments across a tile
  boundary.
- The zero-weight test (`tests/test_fcoverage_command.rs:1028`) does not assert counters (O-005).
- No test covers a selected chromosome that is absent from the weights file (O-004).
- `outliers` still lacks end-to-end tests for `--ignore-gap` and `--reads-are-fragments`, as already
  noted in `.AI/docs/future/fragment_nuisance_outliers_review.md`.

## Checked And Not Currently Findings

- The keep-weight lookup is correct. `minimum_overlapping_keep_weight` indexes `keep_weights` by
  `OverlappingWindow.idx`. In BED mode that is the position in the interval slice
  (`src/shared/overlaps.rs:850-852`), not the line-number index stored in each `IndexedInterval`.
- A lookback of `max_fragment_length` is sufficient. The pairing iterator returns a fragment when it
  reads the later mate (`src/shared/fragment_iterators/core.rs:262`). If fragment B is returned
  before fragment A, then
  `B.start <= B.later_mate_pos <= A.later_mate_pos <= A.start + A.length - 1 < A.start + max_fragment_length`.
  Intervals skipped while processing B therefore cannot overlap A.
- The `1 / (max_fragment_length + 1)` threshold (`src/commands/fcoverage/fcoverage.rs:1046-1047`)
  admits a 1 bp overlap for every fragment. The fragment filter at
  `src/commands/fcoverage/fcoverage.rs:1012-1015` guarantees that no fragment is longer than
  `max_fragment_length`.
- The exact and flanked BED files can be used as cfDNAlab blacklists. The BED parser skips `#` lines
  (`src/shared/bed.rs:82`), and `--blacklist-min-size` defaults to 1 in ends, midpoints, lengths,
  fragment-kmers, bam-to-frag, and bam-to-bam, so short calls are kept.
- Raw coverage is exact. `Coverage` accumulates in `f64` and stores `f32`
  (`src/shared/coverage.rs:70-71` and `:257-268`). The guard at
  `src/commands/outliers/striding.rs:12` and `:64` rejects exactly the first value that `f32` cannot
  store as an exact integer.
- Fragment support is counted exactly once. The midpoint must lie in the tile core
  (`src/commands/outliers/outliers.rs:477`), and the halo equals `max_fragment_length`. Tiles that
  do not align with the stride are handled by merging partial core histograms
  (`src/commands/outliers/outliers.rs:560-591`).
- Called regions are merged across tile boundaries in chromosome and tile order, with masses summed
  (`src/commands/outliers/outliers.rs:664-703`). A blacklisted position splits a region and is
  never bridged.
- ZIP algebra:
  - The first fit matches the zero-truncated Poisson mean and then solves the zero fraction for
    `pi`. It falls back to `pi = 0` when the solution would be negative
    (`src/commands/outliers/model.rs:269-314`).
  - The truncated solution `pi = (r * F - p0) / ((1 - p0) - r * (1 - F))`
    (`src/commands/outliers/model.rs:361-367`) follows from the right-truncated likelihood. Its
    denominator is positive whenever positive values are retained, and `pi < 1` holds whenever
    `r < 1`.
  - The bisection bracket is valid because every conditional mean used is increasing in `lambda`.
- The threshold search (`src/commands/outliers/model.rs:549-626`) uses the inclusive survival
  `P(X >= k)` and only calls the direct tail sum at `k >= lambda`. Both the walk-left and the search
  branches return the smallest passing `k`.
- Hand-checked test expectations in `tests/test_outliers_command.rs`:
  - `paired_bam_produces_the_hand_derived_models_calls_and_all_outputs`: `T1 = 7`, 970 retained
    positions, `T2 = 2`, and keep weight 0.1
  - `boundaries_blacklist_and_short_contig_preserve_calls_and_record_fallback`: `T2 = 2` and keep
    weight 0.1 for every call
  - `detected_outlier_weights_flatten_an_extreme_pileup_in_fcoverage`: `T1 = 13`, `T2 = 12`, and an
    underlying mean of about 2.695
