# `cfdna gc-bias` review

Date: 2026-04-24

Scope: `src/commands/gc_bias/*`, `gc-bias` CLI configuration, correction-package writing/loading boundary, cross-tile window reducers, directly used tiled/window/fragment helpers, README GC-correction usage, and existing `gc-bias` tests in `tests/test_gc_bias.rs`. I did not run tests.

Shared findings that affect this command:

- Original review note: no active shared correctness findings were tracked in `00_shared_package_notes.md` at the time. Current re-review additions below list shared findings that now affect this command.

## Release triage

Post-release performance:

- GB-006: empty BED tiles still open a BAM reader before skipping.

## Findings

### GB-006 - Post-release performance - Empty BED tiles still open a BAM reader before skipping

The tile worker opens a fresh chromosome reader before preparing the tile-window state ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L784-L825)). `prepare_tile_windows()` can identify BED tiles with no candidate windows and return `skip_tile` before reference sequence reads or BAM fetches happen ([windows.rs](../../src/commands/gc_bias/windows.rs#L372-L448)), but the BAM reader has already been opened. The heavier reference sequence read and BAM fetch happen later ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L831-L847), [gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L880-L882)).

Impact: sparse BED runs avoid most per-tile work, but still pay one BAM-reader open per no-window tile. On very sparse targeted designs over whole-genome tiling, that can be avoidable overhead.

Recommended fix:

- For BED mode, prepare/skip windows before opening the BAM reader, using contig lengths already resolved in `run()`.
- Keep the current early return before reference sequence reads and `fetch()`.
- Add a small helper-level test or instrumentation-friendly regression proving no reader is opened for no-window BED tiles.


## Existing coverage notes

The command already has broad coverage for default MAPQ, global/fixed/BED windows, no-count failure, reference end-offset propagation, overlapping/touching BED behavior, aligned vs misaligned single-chromosome tiling, aligned multi-chromosome accumulation, empty middle tiles, fixed-size vs BED cross-tile equivalence on one chromosome, real `ref-gc-bias` integration, saved intermediate sequence/content, minimum window ACGT filtering, outlier methods, hard-clamp behavior, and greedy binning.

The missing coverage from this review is avoiding BAM-reader opens for no-window BED tiles.

## Released-command re-review additions (2026-05-04)

### Shared findings that affect this command

The shared unchecked output-prefix issue (G-022) and reference package identity issue (G-023) originally noted here have since been implemented.

Reviewed shared findings that do not currently apply:

- G-019: `gc-bias` crossing-window temp files are named by tile index, not raw chromosome name.
- G-021 did not affect `gc-bias` because this command does not expose `--gc-tag`.

### Release triage additions

Pre-release correctness/safety:

- None active from this re-review.

Post-release performance:

- GB-006 remains a post-release performance item; no new correctness issue was found in that path.
