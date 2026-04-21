# Post-refactoring review of fcoverage reducer and writers

Reviewed: commit `23246a91` (naming normalisation) plus staged structural refactoring.
Baseline: commit `f3c3ef32` (added `--by-grouped-bed` and `--per-window summary-stats`).

## Verdict

**No functional regressions found.** All reducer semantics, writer outputs, accumulator logic,
and cross-index bookkeeping are preserved. The code is substantially cleaner and easier to audit.

---

## 1. Breaking output changes (committed in 23246a91)

These changes are deliberate naming normalisation, not bugs, but they break any existing
downstream parser or script that references the old column or filename.

| Location | Old | New |
|---|---|---|
| fcoverage aggregate header (average mode) | `avg_coverage` (pre-f3c3ef32) / `mean_coverage` (f3c3ef32) | `average_coverage` |
| fcoverage summary-stats header column | `mean_coverage` | `average_coverage` |
| fcoverage grouped summary-stats header | `mean_coverage` | `average_coverage` |
| WPS aggregate header (average mode) | `avg_coverage` | `average_coverage` |
| WPS output filename | `wps.avg.tsv.zst` | `wps.average.tsv.zst` |
| coverage_weights expected input header | `avg_coverage` | `average_coverage` |
| coverage_weights output header | `avg_pos_cov`, `avg_overlapping_pos_cov` | `average_pos_coverage`, `average_overlapping_pos_coverage` |
| coverage_weights struct fields | `avg_coverage`, `avg_overlap_coverage` | `average_coverage`, `average_overlap_coverage` |
| coverage_weights normalisation fn | `normalize_avg_overlap_by_global_mean` | `normalize_average_overlap_by_global_mean` |

Note: `f3c3ef32` introduced `mean_coverage` in fcoverage but did not update coverage_weights,
which still expected `avg_coverage`. So coverage_weights was already broken on this branch.
`23246a91` fixes that by aligning everything to `average_coverage`.

**Action needed:** if any external tooling or documentation references `mean_coverage`,
`avg_coverage`, or `wps.avg.tsv.zst`, update it before merging.

---

## 2. Structural refactoring (staged changes)

### What was unified

| Before | After | Lines saved |
|---|---|---|
| 4 stream structs (`BedBasicStream`, etc.) | 1 `PartialsStream` with `PartialsSchema` | ~200 |
| 4 row structs | 1 `ParsedPartialRow` | ~80 |
| 4 accumulators | 1 `AggregateAccum` (reducer only) | ~60 |
| 4 merge loops | 2 internal engines + shared helpers | ~200 |
| 4 cross-index loading copies | 1 `load_expected_contributions` | ~120 |
| Inline column parsing everywhere | 1 `parse_col<T>` helper | many |
| Total: reducer.rs ~1619 lines | 959 lines | ~660 |

### What was correctly kept separate

- `GroupedAggregateAccum` stays as its own type in writers.rs (has `span_positions`, no
  `seen_contributions`).
- BED vs size merge engines stay as two functions. They differ in key type (`orig_idx` vs
  `bin_start`), interval source (lookup table vs parsed row), and post-reduction clipping.
- Summary-stats derivation stays in writers, not in reducers.
- The internal engines take `summary: bool` and derive the schema themselves, preventing callers
  from passing e.g. `SizeBasic` to a BED engine.

### Grouped fold extracted cleanly

`fold_reduced_segment_into_group` is now an explicit helper in writers.rs. The fold logic is
identical for basic and summary sources because basic reducers zero the summary-only fields in
the `ParsedPartialRow`. This means both grouped code paths produce the same accumulation, which
is correct.

---

## 3. Line-by-line observations

### 3a. `clip_upper` vs manual `Interval::new(start, end.min(chrom_len))`

reducer.rs:809 now uses `full_interval.clip_upper(chrom_len)`. The original created a new
`Interval::new(start, end.min(chrom_len))`. Both clamp the end to `chrom_len`. The difference is:
- `Interval::new` fails if `start >= end` after clipping
- `clip_upper` returns `None` in that case

The reducer converts `None` to a descriptive `anyhow::anyhow!` error. Both paths produce an
error for the same degenerate input. In practice, `start >= chrom_len` cannot happen because
partial rows only exist for bins that overlap the chromosome.

**Verdict:** equivalent. No risk.

### 3b. `\r` trimming added to line parsing

reducer.rs:175 trims `\r` from line endings. The old code did not. Since all files are written
by the same Rust binary on the same OS, lines will never actually contain `\r`. But the addition
is harmless and improves robustness if files are ever transferred across platforms.

**Verdict:** benign.

### 3c. `allowed_positions` renamed to `eligible_positions`

The old accumulator fields used `allowed_positions` while `ReducedAggregateRow` and grouped
accumulators already used `eligible_positions`. The refactoring normalises the naming. This is
correct and removes a source of confusion.

**Verdict:** improvement.

### 3d. Error diagnostics improved

`parse_col<T>` produces error messages that include the field name, chromosome, tile index, and
line number. The old inline parsing often had less context. This makes debugging partial-file
corruption much easier.

**Verdict:** improvement.

### 3e. `open_maybe_zstd_reader` shared

The helper replaces identical zstd-opening logic that was previously duplicated across each
stream constructor and cross-index loader. Correct and safe.

### 3f. `write_reduced_value_row` in reducer

This helper owns the `finalize_value -> round_to -> write_final_row` chain that was previously
inlined in each reducer public function. The `unmasked_span_bp = interval.len()` derivation
matches the original code for both BED (where interval comes from `coords_by_idx`) and size
(where interval is the clipped full bin).

**Verdict:** equivalent.

---

## 4. Potential future concerns

### 4a. Two merge loops still share ~80% of logic

`reduce_bed_rows_internal` and `reduce_size_rows_internal` both do: open streams -> build heap ->
pop min -> accumulate -> check contribution count -> emit. They differ in:
- key source: `orig_idx` (BED) vs `bin_start` (size)
- interval recovery: lookup table (BED) vs parsed row field (size)
- post-reduction clipping: none (BED) vs `clip_upper(chrom_len)` (size)

A further unification into one generic merge engine is possible (extract key/interval via
closures), but the current two-function split is defensible because the differences are
load-bearing and the code is already clear.

If a merge-loop bug is ever found, remember to fix **both** engines.

### 4b. `ReducedAggregateRow` is the new stable contract

All communication between reducer and writers now flows through `ReducedAggregateRow`. Any future
field (e.g., min/max coverage) must be added to this struct. This is cleaner than the old ad-hoc
approach where each reducer inlined its own output format, but it means the struct must evolve
carefully.

### 4c. Missing drain check is already present

Both engines `ensure!(accum_by_*.is_empty(), ...)` after the merge loop, catching any window or
bin that received fewer contributions than expected. This is a good invariant check.

### 4d. Non-summary grouped fold reads basic temp files

`write_grouped_bed_aggregate_output` calls `reduce_bed_basic_with_cross_index_for_chr_rows` for
non-summary grouped actions. This reads the narrow basic temp files. The `ReducedAggregateRow`
then has `nonzero_positions = 0` and `coverage_sum_of_squares = 0.0`. The grouped fold adds
these zeroed fields to the group accumulator, which is harmless. The final `finalize_grouped_value`
only uses `coverage_sum` and `eligible_positions`, so the zeroed fields never affect the output.

**Verdict:** correct.

---

## 5. Test coverage

The three regression tests added before the refactoring cover:

1. **Cross-tile BED average invariance** (`by_bed_average_is_invariant_when_overlapping_windows_cross_tiles`):
   overlapping BED windows with `tile_size=33` vs `tile_size=1000`. Catches merge-loop bugs that
   produce different results depending on tile decomposition.

2. **Cross-tile BED summary-stats invariance** (`by_bed_summary_stats_is_invariant_when_windows_cross_tiles`):
   same tile-size comparison but for the full summary-stats row including variance, SD, CV,
   and covered_fraction.

3. **Cross-tile size summary-stats invariance** (`by_size_summary_stats_is_invariant_when_bin_crosses_tiles`):
   fixed-size bins that straddle tile boundaries. Verifies that `coverage_sum_of_squares` sums
   correctly across tiles and that the derived variance matches.

These tests are well-designed: they compare two runs with different tile sizes against each other
and against hand-calculated expected values. They would catch any regression in the merge loop,
accumulator arithmetic, or `clip_upper` clipping.

---

## 6. Summary

| Category | Count | Details |
|---|---|---|
| Functional regressions | 0 | All reducer semantics preserved |
| Breaking output changes | 7+ | All in the naming commit, not in the structural refactoring |
| Code improvements | 5 | Unified stream, shared parser, better diagnostics, consistent naming, extracted fold |
| Potential future risks | 2 | Dual merge loops, `ReducedAggregateRow` as contract |
| Tests guarding the refactoring | 3 | Cross-tile invariance for BED average, BED summary, size summary |
