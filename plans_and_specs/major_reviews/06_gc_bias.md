# `cfdna gc-bias` review

Date: 2026-04-24

Scope: `src/commands/gc_bias/*`, `gc-bias` CLI configuration, correction-package writing/loading boundary, cross-tile window reducers, directly used tiled/window/fragment helpers, README GC-correction usage, and existing `gc-bias` tests in `tests/test_gc_bias.rs`. I did not run tests.

Shared findings that affect this command:

- No active shared correctness findings from `00_shared_package_notes.md`; remaining items below are command-specific.

## Release triage

Pre-release correctness/safety:

- GB-001: cross-tile spill files can collide or merge incorrectly across chromosomes.
- GB-002: `--save-intermediates` should respect `--output-prefix`.

Pre-release docs/API polish:

Post-release performance:

- GB-006: empty BED tiles still open a BAM reader before skipping.

## Findings

### GB-001 - High - Cross-tile spill files are keyed by chromosome-local tile/window ids

`Tile.index` resets to zero for each chromosome ([tiled_run.rs](../../src/shared/tiled_run.rs#L488-L516)). `gc-bias` passes that chromosome-local index to `write_crossing_parts()` ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L1160-L1169)), and the writer creates filenames as `cross.<tile_idx>.npz` with no chromosome or global tile id ([cross_tile_parts.rs](../../src/commands/gc_bias/cross_tile_parts.rs#L20-L30)). The crossing file also stores only `idx` plus counts/support arrays ([cross_tile_parts.rs](../../src/commands/gc_bias/cross_tile_parts.rs#L14-L18), [cross_tile_parts.rs](../../src/commands/gc_bias/cross_tile_parts.rs#L54-L57)), and the streaming reducer sorts by that filename index and merges active parts by `idx` only ([cross_tile_parts.rs](../../src/commands/gc_bias/cross_tile_parts.rs#L96-L106), [cross_tile_parts.rs](../../src/commands/gc_bias/cross_tile_parts.rs#L145-L155)).

Impact: any multi-chromosome run that enters the crossing-file path can corrupt reduction. For BED windows this can happen when windows cross tile boundaries; for fixed-size windows it can happen when tiles are not guaranteed aligned to the window grid. Tiles with the same chromosome-local index on different chromosomes write the same temp path, so parallel workers can overwrite each other or return duplicate paths. For fixed-size windows, even after filename collisions are fixed, per-chromosome window indices can still merge across chromosomes when the reducer keys only by `idx`.

Recommended fix:

- Use the global `par_iter().enumerate()` tile index, not `tile.index`, in crossing filenames and reducer sorting.
- Include chromosome/tid or a globally unique window id in every crossing row key. BED windows already carry global BED row ids, but fixed-size windows need a chromosome offset or compound `(tid, window_idx)` key.
- Add a regression with two chromosomes, crossing fixed-size windows or crossing BED windows, and a tile size/window layout that forces the crossing reducer on both chromosomes.

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

The missing coverage from this review is multi-chromosome crossing-file reduction, prefix-safe intermediates, avoiding BAM-reader opens for no-window BED tiles, and copy-pasteable CLI help.
