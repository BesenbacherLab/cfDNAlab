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
