# `fcoverage` grouped BED implementation plan

Date: 2026-04-20

## Goal

Implement `--by-grouped-bed` for `cfdna fcoverage` in a way that:

- supports both existing grouped semantics and grouped unique-base aggregate variants
- supports downstream binary-mask correlation and variance summaries through `summary-stats`
- reuses the current tile and reducer design where possible
- keeps docs and CLI help compact

## Phase 1: CLI and mode plumbing

### 1. Switch `fcoverage` to grouped-capable window args

- replace `WindowsArgs` with `DistributionWindowsArgs` in `FCoverageConfig`
- update `run_inner` to resolve `DistributionWindowSpec`
- keep current no-window, `--by-size`, and `--by-bed` behavior unchanged

### 2. Extend `CoverageWindowAction`

- add `SummaryStats`
- update clap parsing and help text
- validate the allowed support matrix early in `run_inner`

Validation rules to add:

- `summary-stats` is valid for BED, size, and grouped aggregate modes
- grouped BED only allows `average`, `total`, and `summary-stats`
- grouped BED additionally allows:
  - `average-on-unique-bases`
  - `total-on-unique-bases`
  - `summary-stats-on-unique-bases`
- `--by-size` still only allows `average`, `total`, and `summary-stats`

## Phase 2: Grouped BED preprocessing

### 4. Reuse the existing grouped merge path where possible

Before adding new preprocessing code, explicitly reuse or extract the existing merge logic:

- `prepare_windows::mergers::merge_within_groups` is the closest semantic match for grouped
  `*-on-unique-bases`
- `shared::interval::{push_merged_interval, merge_sorted_intervals, TouchingMergePolicy}` already
  contains the lower-level interval-collapse logic
- `shared::bed::Windows::into_flattened_reindexed` is the ordinary-window precedent for
  overlap/touch collapse plus reindexing

Preferred implementation order:

1. reuse an existing shared helper if one already fits cleanly
2. otherwise extract the interval-collapse core from `prepare_windows` into a shared helper
3. only duplicate logic inside `fcoverage` if reuse/extraction turns out to be clearly worse

Do not introduce a third independent definition of "merge within group" if it can be avoided.

### 5. Add a grouped unique-base preprocessing helper

Add a helper in `shared/bed.rs` or a closely related shared module that takes:

- grouped BED windows
- chromosome order

and returns:

- per chromosome, non-overlapping grouped unique-base segments with stable internal `segment_idx`
- `segment_idx -> group_idx`
- per-group unique-base `span_positions`

Required semantics:

- merge overlapping intervals within the same group
- merge touching intervals within the same group
- do not merge across groups
- preserve deterministic group index assignment from the grouped BED loader

This helper is only used for grouped `*-on-unique-bases` actions. Plain grouped actions should
keep the existing grouped interval view as loaded.

The helper should be implemented on top of the reused/extracted merge semantics above, not with a
fresh ad hoc overlap loop.

### 6. Decide internal metadata shape

Use explicit small structs rather than tuples.

Suggested shapes:

```rust
pub struct GroupedCoverageSegment {
    pub interval: IndexedInterval<u64>, // idx is internal segment_idx
    pub group_idx: u64,
}

pub struct GroupedCoverageLayout {
    pub segments_by_chr: FxHashMap<String, Vec<IndexedInterval<u64>>>,
    pub segment_idx_to_group_idx: FxHashMap<u64, u64>,
    pub group_span_positions: FxHashMap<u64, u64>,
    pub group_idx_to_name: FxHashMap<u64, String>,
}
```

The exact struct names can change, but keep the concepts separate.

## Phase 3: Tile counting integration

### 7. Add shared summary-stats support for all aggregate window modes

Implement `summary-stats` for:

- `--by-size`
- `--by-bed`
- `--by-grouped-bed`

Raw per-row stats required:

- `span_positions`
- `blacklisted_positions`
- `eligible_positions`
- `nonzero_positions`
- `coverage_sum`
- `coverage_sum_of_squares`

Derived per-row stats required:

- `mean_coverage`
- `total_coverage`
- `variance_coverage`
- `sd_coverage`
- `coefficient_of_variation_coverage`
- `covered_fraction`

The final reducer should derive these columns from the raw finalized stats, not by rescanning
coverage.

Numeric rules to implement:

- when a derived metric cannot be computed, write `NaN`
- keep integer count columns exact and unclamped
- clamp floating-point values sufficiently close to zero to exact `0.0`
- if roundoff makes variance slightly negative, clamp it to `0.0` before `sqrt`

### 8. Feed unique-base grouped segments into tile-span precomputation

- reuse `precompute_tile_window_spans`
- treat grouped segments like ordinary BED windows for fetch narrowing and overlap lookup
- keep the expensive fragment iteration path unchanged

### 9. Reuse `AggregatesByBed` for grouped aggregates

- run grouped aggregate counting over the unique-base segment view when needed
- write tile partials keyed by `segment_idx`
- reuse cross-index sidecars for boundary-crossing segments

This is the critical simplification. Grouping should happen after segment-level correctness is
settled, not inside the tile hot path.

For plain grouped actions, reuse the grouped intervals as loaded instead of the unique-base segment
view.

### 10. Add tile-level support for `summary-stats`

Add summary-stats tile support that also tracks:

- `nonzero_positions`
- `coverage_sum_of_squares`

Two reasonable implementation choices:

- extend the existing aggregate partial row format with optional `sum_of_squares`
- add dedicated summary-stats partial formats and reducers

Preferred choice:

- use dedicated summary-stats partial formats

Reason:

- it avoids making the existing `average` and `total` reducer path more conditional than it needs
- it keeps the summary-stats semantics explicit

## Phase 4: Final reducers and writers

### 11. Add grouped final output writers

Use the final filenames and grouped headers exactly as specified in the spec.

Important:

- plain grouped actions and `*-on-unique-bases` actions must not share the same output filename
- use the `--per-window` action name in the filename, replacing hyphens with underscores

### 12. Add summary-stats reducers for all aggregate window modes

Implement reducers for:

- `--by-size` summary-stats
- `--by-bed` summary-stats
- `--by-grouped-bed` summary-stats

### 13. Add grouped aggregate reducers

Implement grouped reducers in two steps:

1. reduce per-tile partials into exact per-segment totals
2. fold per-segment totals into per-group totals

For `average` and `total`, per-group accumulators need:

- `span_positions`
- `blacklisted_positions`
- `eligible_positions`
- `coverage_sum`

For `summary-stats`, per-group accumulators additionally need:

- `nonzero_positions`
- `coverage_sum_of_squares`

### 14. Derive final statistics from raw summary-stats output

After raw row stats are finalized, compute:

- `mean_coverage`
- `total_coverage`
- `variance_coverage`
- `sd_coverage`
- `coefficient_of_variation_coverage`
- `covered_fraction`

These must be derived from the raw finalized row values, not from any other output mode.

Implementation detail:

- apply close-to-zero cleanup only to floating-point derived values
- never clamp integer count columns
- emit `NaN` whenever the metric is undefined rather than silently substituting `0.0`

### 15. Keep group row order deterministic

- iterate grouped output rows by ascending `group_idx`
- emit zero rows for groups that received no contributions
- always write `group_index.tsv`

### 16. Use a plain group mapping sidecar

Do not use `write_group_index_with_blacklist_tsv`.

Instead:

- write `group_idx\tgroup_name`
- keep any grouped blacklist-derived metadata out of the sidecar

That avoids publishing interval-weighted metadata that contradicts the unique-base counting
semantics.

## Phase 5: Tests

### 17. Add preprocessing tests

Test the grouped union helper for:

- same-group overlap merge
- same-group touching merge
- different-group overlap remains separate
- stable `group_idx` mapping
- deterministic `segment_idx` assignment
- plain grouped actions bypass the unique-base helper
- parity with the reused/extracted `prepare_windows` within-group merge semantics

### 18. Add reducer tests

Test grouped final reduction for:

- one group with multiple segments on one chromosome
- one group split across chromosomes
- segments spanning multiple tiles
- zero-contribution groups
- zero-eligible groups
- by-size summary-stats reduction
- by-bed summary-stats reduction
- exact derived stats from raw row values
- NaN behavior for undefined derived metrics
- close-to-zero cleanup behavior for floating-point derived metrics
- filename separation between plain grouped and `*-on-unique-bases` outputs

### 19. Add end-to-end `fcoverage` tests

Cover:

- by-size `summary-stats`
- by-bed `summary-stats`
- grouped `average`
- grouped `total`
- grouped `summary-stats`
- grouped `average-on-unique-bases`
- grouped `total-on-unique-bases`
- grouped `summary-stats-on-unique-bases`
- grouped mode with blacklist
- grouped mode with `global` as just another group
- invalid CLI combinations
- coexistence of plain grouped and `*-on-unique-bases` outputs in one output directory
- at least one direct implementation test beyond the standalone mathematical proof, including
  downstream Pearson derivation from `summary-stats-on-unique-bases`

Keep fixtures small and use exact expected outputs.

## Phase 6: Docs

### 20. Keep CLI help short

In `fcoverage/config.rs`:

- document `--by-grouped-bed` as grouped aggregate mode
- document the grouped-only `*-on-unique-bases` modes
- document that same-group intervals are only collapsed in the `*-on-unique-bases` modes
- document that grouped mode does not support positional outputs
- document `summary-stats` in one short block, including which windowing modes support it

### 21. Keep website docs minimal

Update only `website/docs/guides/fragment_coverage_guide.md`.

Add one compact subsection:

- what `summary-stats` means in `fcoverage`
- how plain grouped actions and grouped `*-on-unique-bases` actions differ
- how to include a full-universe `global` group
- what `summary-stats` contains and which columns are derived from raw outputs
- how undefined metrics appear as `NaN`

Do not add a separate guide unless a later command consumes the output directly.

NOTE: Only humans update guides.

## Suggested implementation order

1. CLI enum and validation
2. BED and size `summary-stats`
3. grouped unique-base preprocessing helper
4. grouped average and total output path
5. grouped `summary-stats`
6. tests
7. docs

## Risks and watchpoints

### Semantic risk

The biggest risk is accidentally mixing two grouped meanings:

- interval-weighted grouped distributions
- grouped unique-base aggregates

For `fcoverage`, both may be useful, but they must never share one hidden semantic path.

### Numeric risk

The next biggest risk is inconsistent edge-case handling:

- one reducer writing `0.0` while another writes `NaN`
- clamping integer counts by mistake
- deriving one metric from rounded output instead of raw finalized values

### Reducer risk

Do not aggregate to group level before tile-boundary reconciliation is complete.

If that happens, boundary-crossing segments can be double-counted or partially counted in ways
that are hard to detect.

### Metadata risk

Do not reuse grouped sidecar helpers whose blacklist metrics are based on raw loaded intervals.

That metadata would silently disagree with the grouped rows.

### Docs risk

Do not let grouped support expand into a long matrix of every possible mode. The docs should guide
users toward grouped aggregate coverage and stop there.

## Corrective follow-up

The first implementation pass introduced a few shortcuts that should be removed rather than
accepted as permanent structure.

### 1. Unify summary-stats with the existing per-window computation path

`SummaryStats` for ordinary BED and fixed-size windows should not live on an unnecessary parallel
path when the same per-window coverage math already exists.

Follow-up work:

- extend `WindowValue` with a summary-stats payload
- make `compute_window_outputs()` support `SummaryStats`
- keep `AverageOnUniqueBases` and `TotalOnUniqueBases` as aliases of `Average` and `Total` at that
  layer once the input windows have already been merged upstream
- reserve the grouped-only extra logic for grouped row construction and grouped writing, not for a
  separate copy of the window math

### 2. Add the aligned fast path for `--by-size summary-stats`

When tile boundaries and fixed-size window boundaries align, `summary-stats` should use the same
“write final rows directly per tile” pattern as `average` and `total`, just with the wider summary
payload.

Follow-up work:

- add a size-aligned direct writer for raw and derived summary-stats rows
- keep the general reducer path only for non-aligned `--by-size` runs
- add a dedicated regression test that aligned and non-aligned `summary-stats` produce identical
  output under blacklist, GC correction, and genomic scaling

### 3. Keep `fcoverage` correlation-free

Pearson should stay outside `fcoverage`. The only contract inside `fcoverage` is that
`summary-stats-on-unique-bases` writes the exact raw inputs needed for downstream correlation.

Follow-up work:

- keep the `global` row as ordinary grouped output with no special correlation column
- document the downstream formula separately from the `fcoverage` output contract
- add one direct regression test that derives Pearson from written summary rows and matches the
  explicit positional formula

### 4. Reduce duplicated aggregate output plumbing where possible

The final aggregate writers should differ because the output shapes differ, but the raw summary
statistic computation and derivation rules should be shared rather than copied.

Follow-up work:

- keep one shared summary-stats derivation path
- keep one shared definition of close-to-zero cleanup and undefined-metric handling
- only branch where row semantics genuinely differ:
  - ordinary per-window rows
  - grouped site-weighted rows
  - grouped unique-base rows

### 5. Remove any clutter

Follow-up work:

- remove stale names from the plan and spec so they match the code:
  - `coverage_sum_of_squares`
  - `coefficient_of_variation_coverage`
- remove duplicated section numbering and other stale plan text
- remove dead code and helpers left behind by refactors instead of carrying unused parallel paths
- keep `fcoverage.rs` as orchestration, `writers.rs` as final-output plumbing, and avoid letting
  helper placement drift again

## Done criteria

This work is done when:

- `fcoverage` accepts `--by-grouped-bed`
- `summary-stats` works for BED, size, and grouped aggregate runs
- grouped `average` and `total` work in both plain and unique-base grouped variants
- grouped `summary-stats` writes exact first and second moments
- grouped `summary-stats` also writes `nonzero_positions` and the requested derived metrics
- undefined derived metrics appear as `NaN`
- close-to-zero cleanup applies only to floating-point derived values
- `group_index.tsv` is written for grouped runs
- invalid grouped positional combinations fail clearly
- tests cover grouped preprocessing, reduction, and end-to-end output
- the fragment coverage guide documents grouped coverage in one compact section
