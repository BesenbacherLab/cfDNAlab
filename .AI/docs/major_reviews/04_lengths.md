# `cfdna lengths` review

Date: 2026-04-24

Scope: `src/commands/lengths/*`, the `lengths` CLI configuration, and directly used shared helpers for indel/clip-aware fragment construction, BED/grouped windows, tile spans, GC correction, scaling factors, blacklist checks, midpoint assignment, and length-count reduction. I also read the existing `lengths` tests in `tests/test_lengths_command.rs`. I did not run tests.

Shared findings that affect this command:

- None active. The README OPTIONS labeling issue originally noted here as G-002 has since been implemented.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release docs/API polish:

- None active.

Post-release performance:

- G-006: sparse-window GC reference pruning.

## Existing coverage notes

The command already has broad integration coverage: global, fixed-size, ordinary BED, grouped BED, multi-chromosome runs, tile-boundary invariance, unpaired mode, default MAPQ, GC correction, GC weighting modes, scaling factors, blacklist behavior, indel and clip modes, soft-clip filtering, all window assignment modes, grouped metadata, and reducer helper behavior are represented.

The deferred sparse-window GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.

## Released-command re-review additions (2026-05-04)

The shared raw-chromosome temporary filename issue (G-019) and unchecked output-prefix issue (G-022) originally noted here have since been implemented. The overlong `--gc-tag` issue (G-021) did not affect `lengths` because this command uses file-only GC correction arguments.

### Release triage additions

Pre-release correctness/safety:

- None active from this re-review.

Post-release performance:

- G-006: sparse-window GC reference pruning.

### Command-specific findings

All solved.
