# Grouped BED distributions spec

Date: 2026-04-18

## Scope

This spec defines a `--by-grouped-bed` windowing mode for commands whose outputs are already
"one distribution per row".

Initial target commands:

- `cfdna ends`
- `cfdna lengths`

Explicitly out of scope for this spec:

- any "outside all groups" mode
- union or deduplicated within-group semantics
- positional outputs such as `midpoints`
- `fragment-kmers` implementation

The intent is to match the grouped-BED mental model from `midpoints`, but without a positional
axis in the output. Each output row is one group label, and the row stores one collapsed
distribution for that group.

## Summary

### Core meaning

`--by-grouped-bed` means:

- input intervals come from a BED-like file with columns
  `chromosome, start, end, group_name`
- output rows correspond to `group_name`, not to individual BED intervals
- a fragment contributes to every grouped-BED interval that it qualifies for under the command's
  existing interval-assignment semantics
- the group's output distribution is defined by all contributions made to intervals with that
  group label

This keeps grouped mode conceptually aligned with `midpoints`: many genomic sites are treated as
members of one biological group, and the command reports one aggregate distribution per group.

### Multi-hit same-group semantics

Version 1 uses site-weighted semantics.

That means:

- if a fragment qualifies for two intervals from the same group, both contributions count
- there is no automatic deduplication within a group
- there is no implicit unioning of same-group intervals

This is deliberate. A grouped BED in this mode represents a collection of sites, not a merged
territory.

Consequences:

- grouped counts can exceed `1.0` for one fragment within one group
- this is expected when same-group intervals overlap each other
- the exact amount depends on the command's existing assignment mode

This must be documented clearly. Users should not infer union semantics from the grouped output.

### Assignment semantics stay command-native

`--by-grouped-bed` does not introduce a new assignment rule. It only changes row identity.

Each command keeps using its current assignment options exactly as today for deciding whether a
fragment contributes to an interval and with what weight.

Grouped mode only changes which output row receives those contributions:

- plain `--by-bed`: row identity is the BED interval
- `--by-grouped-bed`: row identity is the interval's `group_idx`

This is the simplest model that stays predictable across commands.

## Why this is the right first version

The alternative would be to define grouped mode as a union of all intervals in each group.

That would be a different biological question:

- site-weighted grouped BED asks "what is the aggregate distribution across this set of sites?"
- unioned grouped BED asks "what is the aggregate distribution across the merged genomic territory
  of this group?"

Those are not interchangeable.

The site-weighted version is the closer analogue to current `midpoints` behavior and to the
expected "many binding sites for one transcription factor" use case.

Union semantics can be added later as a separate mode if needed. It should not be silently folded
into the first grouped-BED implementation.

## CLI shape

Add a new mutually exclusive window flag:

- `--by-grouped-bed <path>`

The BED-like file format is:

- no header
- columns: `chromosome, start, end, group_name`

Shared CLI constraints:

- `--by-grouped-bed` is mutually exclusive with `--by-bed` and `--by-size`
- grouped BED uses the existing grouped loader in `shared/bed.rs`
- group labels are preserved in a sidecar mapping file

This should be represented as a new window mode, not as a hidden interpretation of `--by-bed`.

Working enum shape:

```rust
pub enum WindowSpec {
    Global,
    Size(u64),
    Bed(PathBuf),
    GroupedBed(PathBuf),
}
```

This keeps the semantics explicit and avoids command code having to guess whether a BED file is
plain or grouped.

## Shared behavior

### Row identity

Grouped mode changes the output row key from "window index" to "group index".

That means:

- the grouped BED loader provides `group_idx -> group_name`
- interval overlap lookup still operates on the full list of intervals
- once an interval qualifies, the command counts into that interval's `group_idx`

### Row ordering

Grouped output row order is the grouped BED loader's `group_idx` order.

That order is:

- deterministic
- defined by first occurrence of each `group_name` during the top-to-bottom BED scan
- not alphabetical unless the BED file itself happens to make it alphabetical
- not genomic-sort order of the group's full interval set unless first occurrence happens to imply
  that

So output row order should be documented as:

- row `i` corresponds to `group_idx == i`
- `group_index.tsv` is the authoritative mapping from row index to group label

### Zero-contribution groups

All groups loaded from the grouped BED file must be represented in the output, even when no
fragment qualifies for them.

That means:

- output shape is always `# unique groups`, not `# groups with observed counts`
- groups with no qualifying contributions get an all-zero row
- `group_index.tsv` must still contain those groups

This matches the current "fixed row universe" behavior of the windowed commands and keeps outputs
stable across samples.

### Output metadata

`bins.tsv` is not the right row metadata file in grouped mode.

In plain BED mode, one output row corresponds to one interval, so `bins.tsv` is sufficient.
In grouped mode, one output row corresponds to many intervals, so a file named `bins.tsv` would
misrepresent row identity.

Grouped mode should instead write:

- `group_index.tsv`
  - columns: `group_idx`, `group_name`, `blacklisted_fraction`

`group_index.tsv` is the row metadata file in grouped mode.

`blacklisted_fraction` is aggregated across all loaded intervals in the group as:

- `sum(interval_blacklisted_bp) / sum(interval_bp)`

That is, it is width-weighted across the group's loaded intervals. Intervals are counted exactly as
loaded, so overlapping intervals in the same group contribute separately to both the numerator and
denominator. This matches the grouped counting semantics better than writing interval-level
metadata that no longer lines up with the collapsed output rows.

### Mixed interval widths

Grouped mode does not require equal interval widths for `ends` or `lengths`.

Because these commands output distributions without a positional axis:

- interval width affects which fragments qualify and with what weight
- interval width does not require output reshaping

So mixed widths are allowed in grouped mode for these two commands.

This is different from `midpoints`, where the positional axis requires a uniform window width.

## `ends` fit

### Why `ends` fits well

`ends` already produces one motif distribution per row.

Today:

- global mode produces one row
- `--by-size` and `--by-bed` produce one row per window

Grouped BED is a natural extension:

- `--by-grouped-bed` produces one motif distribution per group

No positional axis or window-width normalization is involved.

### `ends` grouped semantics

`ends` should keep its existing `--assign-by` modes unchanged.

Grouped mode only changes the destination row:

- `endpoint`
  - count the end motif in every grouped interval whose window contains that endpoint
- `count-overlap`
  - use the existing overlap-fraction semantics for each qualifying interval
- `any`
  - count full weight in every qualifying interval
- `all`
  - count full weight in every qualifying interval that fully contains the fragment assignment
    interval
- `midpoint`
  - count full weight in every grouped interval hit by the fragment midpoint
- `proportion=<threshold>`
  - count full weight in every qualifying interval

Because version 1 is site-weighted:

- one endpoint may count multiple times in one group if multiple same-group intervals contain it
- `count-overlap` contributions from same-group intervals may sum to more than `1.0`

This is acceptable and should be documented as expected behavior.

### `ends` outputs in grouped mode

Grouped `ends` should write:

- `end_motifs.npy` or `end_motifs.sparse.npz`
  - shape becomes `(# groups, # motifs)`
- `end_motifs.txt`
- `end_motif_settings.json`
- `group_index.tsv`

Grouped `ends` should not write `bins.tsv`.

### `ends` docs that need explicit wording

The grouped-BED help/docs must state:

- output rows are groups, not intervals
- assignment rules still apply per interval
- same-group overlaps are counted separately

Suggested wording:

"Count motifs in grouped BED intervals and report one motif distribution per group label.
Fragments are evaluated against each grouped interval using the selected `--assign-by` rule.
Intervals with the same group label contribute to the same output row. Overlapping intervals from
the same group are counted separately."

That describes user-visible semantics without explaining reducer internals.

## `lengths` fit

### Why `lengths` fits well

`lengths` is structurally an even cleaner fit than `ends`.

It already writes one fragment-length distribution per row, with shape:

- global: `(# 1, # lengths)`
- windowed: `(# windows, # lengths)`

Grouped BED would simply make that:

- grouped: `(# groups, # lengths)`

No positional normalization is involved, and the current reducer model already thinks in terms of
"one counter per output row".

### `lengths` grouped semantics

`lengths` should keep the current `--assign-by` behavior exactly as today:

- `count-overlap`
  - interval contributions use the current overlap-fraction semantics
- `any`
  - full contribution for each qualifying interval
- `all`
  - full contribution for each qualifying interval
- `midpoint`
  - full contribution for each qualifying interval hit by the midpoint
- `proportion=<threshold>`
  - full contribution for each qualifying interval

As with `ends`, grouped mode is site-weighted:

- same-group intervals are not merged
- same-group overlaps count separately
- grouped mass can exceed `1.0` for one fragment in one group

This is acceptable because the output is an aggregate distribution over sites, not a partition of
genome territory.

### `lengths` outputs in grouped mode

Grouped `lengths` should write:

- `length_counts.npy`
  - shape becomes `(# groups, # lengths)`
- `fragment_length_settings.json`
- `group_index.tsv`

Grouped `lengths` should not write `bins.tsv`.

The optional overall QC plot can stay unchanged. It already plots the sum across all rows, which
still makes sense for grouped mode.

## Shared implementation strategy

This should be implemented as a focused extension, not as a large shared-window refactor.

Recommended strategy:

1. Extend the shared window enum with `GroupedBed`
2. Reuse the existing grouped BED loader
3. Let each command map qualifying intervals to `group_idx` in its own hot path
4. Add grouped-mode output metadata helpers
5. Keep plain-BED and grouped-BED code paths explicit where row identity differs

Do not try to force `ends`, `lengths`, `midpoints`, and `fragment-kmers` into one generic grouped
window abstraction immediately. The commands do not share the same output shape.

### Implementation hazards to call out explicitly

Adding `WindowSpec::GroupedBed` will require touching every exhaustive `match` on `WindowSpec`,
not just the primary loader path.

That includes:

- command-local loader branches
- tile and fetch setup branches
- output-row count and row-id mapping branches
- metadata-writing branches
- shared helpers such as `WindowContext`, `compute_window_offsets`, `build_bin_info`, and
  `fetch_span_for_tile`

This matters because grouped mode is not just "BED with another parser". Some sites should reuse
plain-BED behavior, while others must switch row identity from interval index to `group_idx`.

Related note:

- `GCWindowsArgs::resolve_windows()` returns the same shared `WindowSpec` type, even though it
  should never produce `GroupedBed`

So commands using `GCWindowsArgs` may still need an explicit `GroupedBed` arm after the enum is
extended. In those commands, that arm should usually be an explicit unreachable case or an error,
not silent fallthrough.

### Metadata-writing rule

Grouped metadata writing should be driven by an explicit `WindowSpec::GroupedBed(_)` branch.

Do not implement it as:

- "write `bins.tsv` for every non-global mode"

because that would incorrectly include grouped mode.

The intended rule is:

- `Global`
  - write no per-window or per-group metadata table
- `Size(_)` or `Bed(_)`
  - write `bins.tsv`
- `GroupedBed(_)`
  - write `group_index.tsv`

## Validation requirements

At minimum, grouped mode should be tested for:

- grouped BED parsing failure on missing fourth column
- grouped BED row count equals number of unique groups, not number of intervals
- same-group non-overlapping intervals aggregate into one output row
- same-group overlapping intervals count separately
- plain BED and grouped BED produce identical per-interval qualification behavior before grouping
- grouped metadata files are written and stable
- mixed interval widths are accepted for `ends` and `lengths`

Command-specific tests:

### `ends`

- `endpoint` grouped counting with two same-group intervals
- `count-overlap` grouped counting where same-group overlap yields > 1 grouped weight
- grouped dense output shape under `--all-motifs`

### `lengths`

- grouped `count-overlap` aggregation over multiple same-group intervals
- grouped output ordering follows `group_idx`
- grouped mode still writes valid overall QC plot input

## `fragment-kmers` later, not now

Grouped BED should be considered for `fragment-kmers` later, but it should not be part of this
first implementation.

Reason:

- `fragment-kmers` has a row-per-window model today
- it also has a positional selection model layered on top of that row identity
- grouped aggregation raises extra design questions about how positional outputs should collapse
  across windows in the same group

This is not the same problem as `ends` or `lengths`, where each row is already just a
distribution.

So the package-level rule should be:

- grouped BED is first added only to commands whose output rows are already non-positional
  distributions
- `fragment-kmers` gets its own later design once the desired grouped positional semantics are
  clear

## Non-goals

This spec does not define:

- grouped union semantics
- deduplicated within-group counting
- group-vs-outside comparison modes
- automatic normalization by number of windows per group
- `fragment-kmers` grouped outputs

Those may all be reasonable later, but they should not be bundled into the first grouped-BED
distribution feature.
