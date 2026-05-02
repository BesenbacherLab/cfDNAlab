# `cfdna gc-bias` review

Date: 2026-04-24

Scope: `src/commands/gc_bias/*`, `gc-bias` CLI configuration, correction-package writing/loading boundary, cross-tile window reducers, directly used tiled/window/fragment helpers, README GC-correction usage, and existing `gc-bias` tests in `tests/test_gc_bias.rs`. I did not run tests.

Shared findings that affect this command:

- No active shared correctness findings from `00_shared_package_notes.md`; remaining items below are command-specific.

## Release triage

Pre-release correctness/safety:

- GB-002: `--save-intermediates` should respect `--output-prefix`.

Post-release performance:

- GB-006: empty BED tiles still open a BAM reader before skipping.

## Findings

### GB-002 - Medium - `--save-intermediates` ignores `--output-prefix`

The output-prefix help says the prefix enables writing multiple calls to the same output directory and documents `<prefix>.gc_bias_correction.npz` ([config.rs](../../src/commands/gc_bias/config.rs#L112-L124)). The main correction package and plots use the trimmed prefix ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L243-L246), [gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L715-L738)), but `IntermediateFileSaver` is constructed with only `save_intermediates` and `out_dir` ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L243-L246)). It writes fixed names like `gc_bias.avg_cfdna_counts.0.npy` without the prefix ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L1380-L1417)).

Impact: two prefixed runs in the same output directory can overwrite or mix their `.npy` diagnostic intermediates even though the primary package names remain distinct. The existing intermediate-file test pins the unprefixed names, so it would need updating when this is fixed ([test_gc_bias.rs](../../tests/test_gc_bias.rs#L1907-L1967)).

Recommended fix:

- Add `prefix` to `IntermediateFileSaver` and write names through `dot_join(&[prefix, "gc_bias", file_tag, ...])`, or put intermediates under a per-prefix subdirectory.
- Add a regression with two different prefixes and `save_intermediates = true` in one output directory.

### GB-006 - Post-release performance - Empty BED tiles still open a BAM reader before skipping

The tile worker opens a fresh chromosome reader before preparing the tile-window state ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L784-L825)). `prepare_tile_windows()` can identify BED tiles with no candidate windows and return `skip_tile` before reference sequence reads or BAM fetches happen ([windows.rs](../../src/commands/gc_bias/windows.rs#L372-L448)), but the BAM reader has already been opened. The heavier reference sequence read and BAM fetch happen later ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L831-L847), [gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L880-L882)).

Impact: sparse BED runs avoid most per-tile work, but still pay one BAM-reader open per no-window tile. On very sparse targeted designs over whole-genome tiling, that can be avoidable overhead.

Recommended fix:

- For BED mode, prepare/skip windows before opening the BAM reader, using contig lengths already resolved in `run()`.
- Keep the current early return before reference sequence reads and `fetch()`.
- Add a small helper-level test or instrumentation-friendly regression proving no reader is opened for no-window BED tiles.


## Existing coverage notes

The command already has broad coverage for default MAPQ, global/fixed/BED windows, no-count failure, reference end-offset propagation, overlapping/touching BED behavior, aligned vs misaligned single-chromosome tiling, aligned multi-chromosome accumulation, empty middle tiles, fixed-size vs BED cross-tile equivalence on one chromosome, real `ref-gc-bias` integration, saved intermediate sequence/content, minimum window ACGT filtering, outlier methods, hard-clamp behavior, and greedy binning.

The missing coverage from this review is prefix-safe intermediates and avoiding BAM-reader opens for no-window BED tiles.
