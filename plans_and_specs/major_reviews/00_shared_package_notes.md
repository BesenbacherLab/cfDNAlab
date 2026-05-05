# Shared package review notes

Date: 2026-04-24

Scope: package-level findings discovered while reviewing released commands. Do not duplicate these findings in command-specific review files unless there is a command-specific consequence that needs separate handling.

This file now tracks only active shared findings. Findings that were already implemented have been removed from the active queue. Points that were reviewed and deliberately rejected as findings are listed separately so they are not re-added in later command reviews.

Implemented findings removed from active tracking:

- G-001: shared fragment length ranges can be inverted without a direct error.
- G-003: tiled commands clean temporary directories only on successful completion.
- G-004: `--gc-file` lacks shared fail-fast validation for the required `--ref-2bit`.
- G-005: shared even-fragment midpoint tie-breaking is not reproducible.
- G-007: plain BED modes need a shared no-surviving-windows guard.
- G-010: GC correction packages cannot identify the sample or inputs they were built from.
- G-011: scaling-factor TSV compatibility metadata is minimal but sufficient for first release.
- G-012: short final stride bins are length-weighted only in the numerator.
- G-013: smoothing-weight docs claim the inverted scaling factors have mean 1.0.
- G-014: smoothing-weight TSV writes do not explicitly flush the final writer.
- G-015: all-zero smoothing runs fail with a misleading normalization error.
- G-016: smoothing-weight commands cannot match `fcoverage --ignore-gap` segmentation.
- G-022: output prefixes can escape the output directory.

Reviewed and deliberately not tracked:

- G-020: release builds with `plotters` can create QC plots or other auxiliary files by default. This is not a review finding by itself. Commands may legitimately have multiple outputs, and automated pipelines should depend on the specific required input/output artifacts they consume rather than failing because extra files are present. Do not re-add this as an issue unless a command overwrites a primary output, hides a required output contract, or makes a requested primary artifact unavailable.

## Pre-release docs/API polish

### G-002 - Low - README OPTIONS blocks need clearer alternative-choice labeling

The README examples are intentionally written as command skeletons with an `OPTIONS` section, not as many separate runnable snippets. That structure is acceptable and does not need a radical rewrite. The remaining issue is narrower: some OPTIONS groups list alternatives that a new user might read as cumulative unless the alternative-choice contract is stated explicitly.

The `fcoverage` OPTIONS block lists mutually exclusive window selectors together: `--by-size`, `--by-bed`, and `--by-grouped-bed` ([README.md](../../README.md#L275-L294)).

The `ends` and `lengths` OPTIONS blocks also list mutually exclusive window selectors together ([README.md](../../README.md#L316-L345), [README.md](../../README.md#L361-L385)). That is fine as long as the text makes clear that these are alternatives, not one command to run unchanged.

The `midpoints` block shows two alternative `--length-bins` forms inside the OPTIONS section ([README.md](../../README.md#L397-L417)). The intent is clear to a careful reader, but a short note would make the contract explicit.

This is documentation polish, not a code-correctness issue and not a reason to replace the OPTIONS format.

Recommended fix:

- Keep the current command-plus-OPTIONS structure.
- Add one explicit sentence before the OPTIONS blocks: choose one option from each alternative group; do not run every OPTIONS line together.
- Where alternatives are shown, label them as alternatives in prose or comments rather than implying they are cumulative arguments.

## Post-release performance optimizations

### G-006 - Sparse-window reference sequence reads happen before no-window pruning

Several tiled commands build GC prefixes from the full tile fetch span before asking the windowing code whether a sparse BED/grouped run should skip or narrow the tile. The pattern appears in `fcoverage` ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L780-L807)), `midpoints` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L368-L415)), `ends` ([ends.rs](../../src/commands/ends/ends.rs#L612-L645)), and `lengths` ([lengths.rs](../../src/commands/lengths/lengths.rs#L633-L668)).

`ref-gc-bias` has the same ordering in its BED mode: `process_tile()` reads the tile sequence and builds GC prefixes before deriving `tile_windows`, and only then returns early for empty tile-window sets ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L416-L475)).

Impact: targeted runs can spend substantial I/O and CPU reading and prefix-building reference sequence for tiles that cannot contribute any output. This directly hurts the use case where sparse windowing should be fastest. This is a performance optimization, not a first-release correctness blocker.

Recommended future fix:

- Move fetch-span adaptation before GC-prefix construction.
- Build GC prefixes over the narrowed fetch span and shift fragment coordinates relative to that narrowed span.
- Add a shared helper-level regression, or one command-level regression per fetch helper family, proving no reference sequence is requested for a no-window tile.
