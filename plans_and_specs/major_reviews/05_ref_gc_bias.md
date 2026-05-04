# `cfdna ref-gc-bias` review

Date: 2026-04-24

Scope: `src/commands/ref_gc_bias/*`, the reference GC package writer/loader boundary, sampling helpers, support-mask helpers, README usage for the GC correction pipeline, and existing `ref-gc-bias` tests in `tests/test_ref_gc_bias.rs` plus module tests in `src/commands/ref_gc_bias/ref_gc_bias_tests.rs`. I did not run tests.

Shared findings that affect this command:

- G-002 in `00_shared_package_notes.md`: README OPTIONS blocks need clearer alternative-choice labeling.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release semantic/docs:

- G-002: README OPTIONS blocks should keep their current structure but clarify alternative choices.

Post-release performance:

- G-006: sparse-window reference pruning.

## Findings

No command-specific findings remain open in this review.

## Existing coverage notes

The command already has good coverage for written package shapes and scalar metadata, exact distributions, blacklist masking, end offsets, smoothing, interpolation, BED flattening, full-chromosome BED vs global mode, multiple blacklist files, rejection when sampling density exceeds 1.0, and fixed-seed determinism for thread count and identical tile size.

The deferred sparse-window reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.

## Released-command re-review additions (2026-05-04)

Shared findings that affect this command:

- G-022 in `00_shared_package_notes.md`: `--output-prefix` can escape the output directory for the `.ref_gc_package.npz` output.
- G-023 in `00_shared_package_notes.md`: reference GC packages do not record or check the reference identity later used by `gc-bias`.

Shared findings reviewed and not applied:

- G-019 does not affect `ref-gc-bias`: this command tiles in memory and does not write per-chromosome temporary tile files.
- G-021 does not affect `ref-gc-bias`: this command does not read BAM GC tags.

### Release triage additions

Pre-release correctness/safety:

- G-023: reference GC package/reference genome identity mismatch can pass silently into `gc-bias`.
- G-022: unchecked output prefix for the package output path.

Post-release performance:

- G-006: sparse-window reference sequence pruning.

### Command-specific findings

No new command-specific findings beyond the shared artifact and prefix findings above.
