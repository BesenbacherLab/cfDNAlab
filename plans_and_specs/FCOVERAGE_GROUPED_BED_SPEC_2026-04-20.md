# `fcoverage` grouped BED spec

Date: 2026-04-20

## Scope

This spec defines grouped BED support for `cfdna fcoverage`.

The main goals are:

- support `--by-grouped-bed` in `fcoverage`
- define lightweight and summary aggregate modes for grouped BED under both site-weighted and
  unique-base semantics
- define a `summary-stats` aggregate mode that also works for `--by-size` and `--by-bed`
- keep grouped behavior exact for binary-mask correlation workflows
- keep grouped implementation aligned with the existing `fcoverage` tile and reducer model
- avoid cluttering `fcoverage` docs with modes that are not actually useful here

Initial scope includes:

- grouped aggregate outputs
- grouped sidecar row metadata
- grouped lightweight aggregate modes under unique-base semantics
- an exact grouped summary-statistics mode for downstream analysis, including Pearson correlation to
  binary masks
- matching summary-statistics support for `--by-size` and `--by-bed`

Explicitly out of scope for this spec:

- grouped positional outputs
- a direct `cfdna correlation` command
- any "outside all groups" synthetic mode
- changing existing grouped-BED semantics in other commands

## Summary

### Core meaning

`--by-grouped-bed` means:

- input intervals come from a BED-like file with columns
  `chromosome, start, end, group_name`
- output rows correspond to `group_name`, not to individual BED intervals
- the grouped row semantics depend on an explicit grouped mode selector
- overlaps across different groups are allowed and counted independently

This spec should not silently redefine `--by-grouped-bed` for `fcoverage` while `ends` and
`lengths` already use site-weighted grouped semantics. If `fcoverage` needs unique-base grouped
counting for binary-mask work, that must be a deliberate grouped action, not an accidental
semantic fork hidden behind the same flag.

### Grouped semantics live in the aggregate action

When `--by-grouped-bed` is used, the grouped semantic split should be expressed through the
`--per-window` action, not through an extra grouped-BED flag.

Reason:

- `fcoverage` already has many arguments
- the semantic difference only matters for grouped aggregate outputs
- the special unique-base meaning should be attached directly to the action that needs it

Definitions:

- default grouped actions
  - intervals are counted exactly as loaded
  - same-group overlaps count multiple times
  - grouped rows represent collections of sites, consistent with `ends` and `lengths`
- `*-on-unique-bases` grouped actions
  - same-group overlaps and touches are merged before counting
  - grouped rows represent unique bases within each group
  - these are the only grouped actions that support exact binary-mask interpretation and Pearson
    correlation to a universe row

### Global is just another group

There is no special global code path at the API level.

If the grouped BED contains full-chromosome windows assigned to a group such as `global` or
`__global__`, that group is treated exactly like every other group:

- same blacklist handling
- same GC correction
- same genomic smoothing
- same fragment filters
- same chromosome subset
- same normalization rules

This keeps the semantics auditable and removes any risk of drift between "global" and
"non-global" counting.

In grouped summary outputs, Pearson correlation is only meaningful when:

- the selected action is `summary-stats-on-unique-bases`
- one grouped row is designated as the full analysis universe

### Why grouped positional output is out of scope

Grouped positional output would either:

- duplicate positions when groups overlap
- require group-indexed positional rows that are large and awkward to consume
- expand CLI help and docs for a mode that is not needed for the motivating workflow

That is unnecessary complexity. Grouped support in `fcoverage` should stay aggregate-only in v1.

## CLI shape

### Window selection

`FCoverageConfig` should switch from `WindowsArgs` to `DistributionWindowsArgs`, so
`fcoverage` gets:

- `--by-size <bp>`
- `--by-bed <path>`
- `--by-grouped-bed <path>`

These remain mutually exclusive.

Working enum shape:

```rust
pub enum DistributionWindowSpec {
    Global,
    Size(u64),
    Bed(PathBuf),
    GroupedBed(PathBuf),
}
```

### `--per-window` values

Keep the existing values:

- `average`
- `total`
- `unique-positions`
- `indexed-positions`

Add one new value:

- `summary-stats`
- `average-on-unique-bases`
- `total-on-unique-bases`
- `summary-stats-on-unique-bases`

### Support matrix

- no windows
  - positional genome-wide output as today
- `--by-size`
  - `average`
  - `total`
  - `summary-stats`
- `--by-bed`
  - `average`
  - `total`
  - `unique-positions`
  - `indexed-positions`
  - `summary-stats`
- `--by-grouped-bed`
  - `average`
  - `total`
  - `summary-stats`
  - `average-on-unique-bases`
  - `total-on-unique-bases`
  - `summary-stats-on-unique-bases`

Invalid combinations must fail with direct error messages:

- `--by-grouped-bed --per-window unique-positions`
- `--by-grouped-bed --per-window indexed-positions`
- `--per-window average-on-unique-bases` without `--by-grouped-bed`
- `--per-window total-on-unique-bases` without `--by-grouped-bed`
- `--per-window summary-stats-on-unique-bases` without `--by-grouped-bed`

### Universe row selection for Pearson

If grouped `summary-stats` should include Pearson correlation to a binary group mask, the run must
explicitly identify which grouped row is the full analysis universe.

Working CLI shape:

```text
--summary-stats-universe-group <group_name>
```

Rules:

- only valid with `--by-grouped-bed --per-window summary-stats-on-unique-bases`
- Pearson values are derived from the raw stats written for that designated group row
- if the option is omitted, the run still writes full summary stats but does not compute valid
  Pearson values

## Summary-stats as a shared aggregate mode

### Purpose

`summary-stats` should be the machine-friendly aggregate output mode for `fcoverage`.

It should exist for:

- `--by-size`
- `--by-bed`
- `--by-grouped-bed`

It should write both:

- exact raw statistics that are composable and auditable
- derived statistics computed only from those raw output values

This keeps `average` and `total` as simple convenience modes while making `summary-stats` the
fuller downstream-friendly output.

### Raw statistics required

Every `summary-stats` row should include enough exact raw statistics to derive the requested final
metrics without revisiting the coverage track:

- `span_positions`
- `blacklisted_positions`
- `eligible_positions`
- `nonzero_positions`
- `coverage_sum`
- `coverage_sumsq`

Definitions:

- `nonzero_positions` counts eligible positions with final coverage strictly greater than zero
- `coverage_sum = Σx`
- `coverage_sumsq = Σx²`

### Derived statistics required

Every `summary-stats` row should also include these derived columns, computed from the raw stats
in the final output:

- `mean_coverage`
- `total_coverage`
- `variance_coverage`
- `sd_coverage`
- `cv_coverage`
- `covered_fraction`

Definitions:

- `mean_coverage = coverage_sum / eligible_positions`
- `total_coverage = coverage_sum`
- `variance_coverage = coverage_sumsq / eligible_positions - mean_coverage^2`
- `sd_coverage = sqrt(variance_coverage)`
- `cv_coverage = sd_coverage / mean_coverage`
- `covered_fraction = nonzero_positions / eligible_positions`

Numerics:

- derived columns must be computed after the exact raw stats are finalized
- do not compute derived columns from rounded `average` or `total` outputs
- after derivation, clamp finite floating-point values that are sufficiently close to zero to
  exact `0.0` (we already have implementations of these, reuse those)
- if floating-point roundoff yields a tiny negative variance, clamp it to `0.0` before `sqrt`
- when a derived metric cannot be computed, it should be `NaN`
- when `eligible_positions == 0`, derived coverage metrics should be `NaN` except
  `total_coverage`, which stays `0.0`
- when `mean_coverage == 0`, `cv_coverage` should be `NaN`

Clamping rules:

- only apply the close-to-zero cleanup to floating-point coverage-derived values
- do not clamp integer count columns such as `span_positions`, `eligible_positions`, or
  `nonzero_positions`
- do not clamp nonzero finite values to one another, only to exact zero

### Pearson correlation in grouped summary-stats

Grouped `summary-stats` may additionally include:

- `pearson_r_to_universe_binary_mask`

This column is valid only when:

- the selected action is `summary-stats-on-unique-bases`
- `--summary-stats-universe-group` identifies a row that represents the full analysis universe

The column must be computed from the raw stats already present in the final output:

- from the designated universe row:
  - `n = eligible_positions`
  - `S = coverage_sum`
  - `Q = coverage_sumsq`
- from each target grouped row:
  - `n1 = eligible_positions`
  - `S1 = coverage_sum`

Then:

```text
r = (n*S1 - S*n1) / sqrt((n*Q - S^2) * (n*n1 - n1^2))
```

Behavior:

- if the selected action is plain `summary-stats`, `pearson_r_to_universe_binary_mask` must be
  `NaN`
- if no universe row is designated, `pearson_r_to_universe_binary_mask` must be `NaN`
- the universe row itself should get `NaN`, because the binary mask is all ones and the variance is
  zero

## Grouped BED semantics

### Site-weighted grouped actions

Plain grouped aggregate actions match the existing grouped-BED semantics already specified for
`ends` and `lengths`.

That means:

- intervals are counted exactly as loaded
- overlapping intervals in the same group contribute multiple times
- touching intervals stay separate
- grouped rows represent site collections, not merged territory

This is the default grouped mode for consistency across commands.

### Grouped actions on unique bases

The `*-on-unique-bases` grouped actions define grouped rows over the unique bases of each group.

That means:

- overlapping intervals in the same group must be merged before counting
- touching intervals in the same group should also be merged
- merged segments inside one group are counted once
- same-group overlaps must not double-count coverage or eligible positions

These are the only grouped actions where the grouped row means unique bases and where Pearson to a
binary mask is valid.

Implementation note:

- reuse the existing "merge within group" semantics where possible rather than reimplementing them
- `prepare_windows::mergers::merge_within_groups` already expresses the correct grouped behavior:
  merge within one group, do not merge across groups, and treat touching intervals as merged when
  the gap is zero
- if that code is too coupled to `prepare_windows::Window`, extract the shared interval-collapse
  core into a shared helper and let both commands call it
- `shared::interval::{push_merged_interval, merge_sorted_intervals, TouchingMergePolicy}` and
  `shared::bed::Windows::into_flattened_reindexed` are the existing lower-level precedents to
  build on

### Cross-group overlap semantics

Intervals from different groups remain independent.

If two groups overlap the same base:

- that base contributes to both groups
- this is expected
- no cross-group deduplication is performed

The grouped BED represents many binary masks on the same genome, not one partition of the genome.

### Blacklist and eligible bases

For each grouped row:

- in plain grouped actions:
  - `span_positions` is the total loaded interval width for the group
  - `blacklisted_positions` is the total blacklisted bp across the loaded intervals
  - same-group overlaps contribute multiple times
- in `*-on-unique-bases` grouped actions:
  - `span_positions` is the number of unique bases in the group before blacklist masking
  - `blacklisted_positions` is the number of blacklisted bases inside those unique bases
  - same-group overlaps contribute once

In both cases:

- `eligible_positions = span_positions - blacklisted_positions`

All grouped aggregates must be based on `eligible_positions`, using the same blacklist semantics as
the rest of `fcoverage`.

### Zero-eligible groups

All groups from the grouped BED file must appear in the output, even if:

- they receive no fragment coverage
- all their bases are blacklisted
- they are on selected chromosomes but have zero eligible positions

This keeps row identity fixed across samples.

For grouped `average`, a zero-eligible group should follow the current masked `fcoverage`
convention and output `0.0`.

## Output files

Summary-stats should follow the same top-level filename across aggregate windowing modes, with row
schema depending on the selected windowing mode. This keeps the mode discoverable without creating
one filename family per window source.

### Main output names

- `<prefix>.fcoverage.avg.tsv.zst`
- `<prefix>.fcoverage.total.tsv.zst`
- `<prefix>.fcoverage.summary_stats.tsv.zst`
- `<prefix>.fcoverage.avg_on_unique_bases.tsv.zst`
- `<prefix>.fcoverage.total_on_unique_bases.tsv.zst`
- `<prefix>.fcoverage.summary_stats_on_unique_bases.tsv.zst`

Filename rule:

- the `--per-window` action name should be reflected in the output filename
- use underscores in filenames where the CLI action uses hyphens
- plain grouped actions reuse the standard aggregate filenames
- `*-on-unique-bases` actions must write distinct filenames so plain grouped outputs and
  unique-base grouped outputs can coexist safely in the same output directory

For grouped runs, the output header and `group_index.tsv` sidecar define the row identity.

### Sidecar name

- `<prefix>.group_index.tsv`

This file maps stable row ids to group names.

## Output schemas

### Row meaning by action

| Windowing | Action family | One row means |
| --- | --- | --- |
| `--by-size` | plain actions | one fixed genomic window |
| `--by-bed` | plain actions | one BED interval |
| `--by-grouped-bed` | plain actions | one site-weighted grouped collection |
| `--by-grouped-bed` | `*-on-unique-bases` actions | one grouped unique-base aggregate |

### Grouped `average`

Header:

```text
group_idx	span_positions	blacklisted_positions	eligible_positions	avg_coverage
```

Definitions:

- `group_idx`: stable grouped row id
- `span_positions`, `blacklisted_positions`, and `eligible_positions` follow the selected grouped
  mode semantics described above
- `avg_coverage`: `coverage_sum / eligible_positions`, or `0.0` when `eligible_positions == 0`

### Grouped `total`

Header:

```text
group_idx	span_positions	blacklisted_positions	eligible_positions	total_coverage
```

`total_coverage` is the sum of per-base coverage across eligible bases in the selected grouped
mode semantics.

### BED or size `summary-stats`

Header:

```text
chromosome	start	end	span_positions	blacklisted_positions	eligible_positions	nonzero_positions	coverage_sum	coverage_sumsq	mean_coverage	total_coverage	variance_coverage	sd_coverage	cv_coverage	covered_fraction
```

### Grouped `summary-stats`

Header:

```text
group_idx	span_positions	blacklisted_positions	eligible_positions	nonzero_positions	coverage_sum	coverage_sumsq	mean_coverage	total_coverage	variance_coverage	sd_coverage	cv_coverage	covered_fraction	pearson_r_to_universe_binary_mask
```

For plain grouped `summary-stats`, `pearson_r_to_universe_binary_mask` must be `NaN`.

For grouped `summary-stats-on-unique-bases`, the global row does not get special columns. It is
just the row for whichever grouped BED label represents the full analysis universe.

### Group index sidecar

Header:

```text
group_idx	group_name
```

Do not use the existing grouped blacklist-fraction sidecar for this mode. That helper in
`shared/windowing.rs` explicitly counts intervals exactly as loaded, including same-group overlap.
That metadata would be inconsistent with the required unique-base semantics here.

## Internal design

### Reuse the existing bed aggregate counting path

The expensive part of `fcoverage` is already solved:

- tile-local counting over a set of indexed BED windows
- partial rows written per tile
- cross-index sidecars for windows that cross tile boundaries
- a final reducer that merges per-window tile contributions

Grouped support should reuse that path through one of two internal views:

- plain grouped actions use grouped intervals as loaded
- `*-on-unique-bases` grouped actions convert grouped input into a unique-base segment view

### Required grouped preprocessing

The grouped BED loader currently provides:

- intervals keyed by chromosome
- `group_idx -> group_name`

For `*-on-unique-bases` grouped actions, preprocessing must additionally build:

- per chromosome, a list of non-overlapping unique-base segments
- each segment carrying a unique `segment_idx`
- a mapping `segment_idx -> group_idx`
- per group, the unique-base `span_positions`

This preprocessing should be deterministic and preserve stable `group_idx` assignment from the BED
loader.

### Counting model

Tile counting should treat the selected grouped representation exactly like ordinary BED windows:

- in plain grouped actions:
  - overlap lookup uses grouped intervals as loaded
  - tile partials are keyed by interval identity
- in `*-on-unique-bases` grouped actions:
  - overlap lookup uses unique-base grouped segments
  - tile partials are keyed by `segment_idx`
  - segments that cross tile cores use the existing cross-index mechanism

No grouped aggregation should happen inside the hot counting loop.

### Final grouped reduction

Grouped finalization needs one extra fold after per-segment reduction:

1. merge tile contributions for each interval or unique-base segment identity
2. finalize each segment's aggregate values exactly once
3. accumulate finalized segment totals into a per-group accumulator
4. write one row per `group_idx`

This keeps cross-tile correctness local to the existing reducer model.

### `summary-stats` accumulation

Per-segment counting must retain enough information to compute:

- `nonzero_positions`
- `coverage_sum`
- `coverage_sumsq`
- `blacklisted_positions`
- `eligible_positions`

The simplest exact model is to accumulate the necessary raw moments and counts at tile level when
`--per-window summary-stats` is selected.

That should not be bolted onto the current average/total finalizer by inference. `Σx²` and
`nonzero_positions` are real counting outputs and must be tracked explicitly.

The same applies to `summary-stats-on-unique-bases`, except the grouped geometry first collapses to
unique-base segments before those raw statistics are accumulated.

## Validation rules

The implementation should error if:

- the grouped BED has no valid rows on selected chromosomes
- `--summary-stats-universe-group` is requested outside grouped `summary-stats-on-unique-bases`
  mode
- grouped positional output is requested
- the reducer sees duplicate internal segment ids within one chromosome mapping
- a grouped row expected from the BED metadata is missing from the final output universe

It should not error merely because:

- groups overlap each other
- a group has zero coverage
- a group has zero eligible positions

## Docs strategy

The docs should stay compact.

CLI help should:

- mention that `summary-stats` is available for global, BED, size, and grouped aggregate runs
- mention that grouped BED rows are aggregated by group name
- mention that grouped BED supports plain grouped actions and `*-on-unique-bases` grouped actions
- mention that grouped mode is aggregate-only

Website docs should only add one compact subsection to the fragment coverage guide:

- grouped and non-grouped `summary-stats`
- one short example with a synthetic `global` group
- one short note explaining that Pearson is only valid for `summary-stats-on-unique-bases` with a
  designated universe row

Do not create a separate long-form guide unless a later command consumes this output directly.

## Test requirements

Tests must cover:

- union of overlapping intervals within one group
- merging of touching intervals within one group
- independence of overlapping intervals across different groups
- distinction between plain grouped actions and `*-on-unique-bases` grouped actions
- grouped average denominator with blacklist masking
- grouped total over unique bases
- grouped reduction across tile boundaries
- exact `coverage_sum`, `coverage_sumsq`, and `nonzero_positions` for `summary-stats`
- exact derived statistics from raw summary-stats output
- Pearson derivation from the designated universe row in grouped `summary-stats-on-unique-bases`
- NaN behavior for undefined derived metrics
- close-to-zero clamping for floating-point derived metrics only
- distinct filenames for plain grouped outputs versus `*-on-unique-bases` outputs
  mode
- zero-eligible group output
- stable `group_idx` mapping
- invalid mode combinations
