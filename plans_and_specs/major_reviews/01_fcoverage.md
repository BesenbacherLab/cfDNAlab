# `cfdna fcoverage` review

Date: 2026-04-24

Scope: `src/commands/fcoverage/*`, the CLI dispatch for `fcoverage`, and directly used shared helpers. I also read the existing `tests/test_fcoverage_command.rs` function list and targeted sections to understand coverage shape. I did not run tests.

Shared findings that affect this command:

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Post-release performance:

- G-006: sparse-window GC reference pruning.

## Existing coverage notes

The command already has unusually broad end-to-end coverage in `tests/test_fcoverage_command.rs`, including normalization, restore-mean, gap handling, blacklist masking, tile-boundary invariance, by-size/by-BED/grouped-BED modes, GC file/tag paths, and invalid mode combinations. The deferred sparse-window GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.

## Re-review additions (2026-05-04)

Shared findings that affect this command:

- G-019 in `00_shared_package_notes.md`: tiled temporary files use raw chromosome names as path components.

### Release triage additions

Pre-release correctness/safety:

- FCOV-001: grouped `average` and `average-on-unique-bases` report `0` instead of `NaN` when no positions are eligible.
- FCOV-002: grouped BED mode can succeed with header-only outputs when no grouped windows survive chromosome selection.
- G-019: raw chromosome names in per-tile temporary filenames.

Post-release performance:

- G-006 remains the only fcoverage-specific performance item from this pass.

### Findings

#### FCOV-001 - Medium - Grouped average rows with no eligible bases report `0` instead of `NaN`

The command-level docs state that window averages are `NaN` when no positions are eligible, for example when a window is fully blacklisted ([config.rs](../../src/commands/fcoverage/config.rs#L43-L48)). The ordinary BED/by-size reducer path follows that rule: `finalize_value()` returns `f64::NAN` for average modes with no allowed positions ([tiling.rs](../../src/commands/fcoverage/tiling.rs#L575-L610)), and the reduced row writer uses that helper for non-summary aggregate output ([reducer.rs](../../src/commands/fcoverage/reducer.rs#L526-L553)).

Grouped non-summary output has a separate finalization helper. For `Average` and `AverageOnUniqueBases`, `finalize_grouped_value()` returns `0.0` when `accum.eligible_positions == 0` ([writers.rs](../../src/commands/fcoverage/writers.rs#L483-L501)), and `write_grouped_bed_aggregate_output()` uses that helper for final grouped rows ([writers.rs](../../src/commands/fcoverage/writers.rs#L760-L779)).

Impact: a grouped row that is fully blacklisted is reported as true zero coverage instead of an undefined average. That makes it indistinguishable from a valid, unmasked group with zero supporting fragments. Existing tests cover the ordinary no-eligible average case ([test_fcoverage_command.rs](../../tests/test_fcoverage_command.rs#L3105-L3144)) and grouped average cases with nonzero eligible positions ([test_fcoverage_command.rs](../../tests/test_fcoverage_command.rs#L5978-L6035), [test_fcoverage_command.rs](../../tests/test_fcoverage_command.rs#L6182-L6235)), but I did not find a grouped average test where `eligible_positions == 0`.

Recommended fix:

- Make grouped average finalization return `NaN` when `eligible_positions == 0`.
- Add grouped BED tests for `average` and `average-on-unique-bases` with a fully blacklisted group.
- Keep total modes returning `0` for empty eligible denominators, because totals do not divide by the denominator.

#### FCOV-002 - Medium - Grouped BED mode can silently produce header-only output after chromosome filtering

Plain BED mode explicitly fails when no windows survive chromosome selection by calling `ensure_plain_bed_windows_not_empty()` after loading windows ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L209-L238), [windowing.rs](../../src/shared/windowing.rs#L116-L127)). Grouped BED mode does not make the equivalent check: it loads grouped windows, builds a grouped coverage layout, and continues even when the selected chromosome set leaves no grouped rows ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L240-L257)).

The grouped BED loader skips disallowed chromosomes before registering their group names ([bed.rs](../../src/shared/bed.rs#L563-L587)), so a grouped BED whose rows are all on filtered-out chromosomes yields an empty grouped layout and empty group-name map. The final grouped writer iterates `grouped_layout.group_idx_to_name.keys()` for output row order ([writers.rs](../../src/commands/fcoverage/writers.rs#L737-L779)), which means the command can successfully write only headers plus an empty group index sidecar.

Impact: a chromosome mismatch, typo, or over-restrictive chromosome selection can look like a successful zero-row analysis. That differs from plain BED behavior and is easy to miss in automated pipelines. Existing grouped tests confirm groups on filtered-out chromosomes are ignored when at least one selected group remains ([test_fcoverage_command.rs](../../tests/test_fcoverage_command.rs#L2373-L2414), [test_fcoverage_command.rs](../../tests/test_fcoverage_command.rs#L5856-L5917)), but I did not find a test for the all-filtered grouped case.

Recommended fix:

- Add a grouped equivalent of `ensure_plain_bed_windows_not_empty()` and call it after grouped BED loading/layout construction.
- Add an end-to-end test where `--chromosomes` excludes every grouped BED row and assert a direct error.
- Preserve the current behavior where groups on unselected chromosomes are omitted when other selected groups remain.
