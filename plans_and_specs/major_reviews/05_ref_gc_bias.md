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
