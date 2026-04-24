# `cfdna gc-bias` review

Date: 2026-04-24

Scope: `src/commands/gc_bias/*`, `gc-bias` CLI configuration, correction-package writing/loading boundary, cross-tile window reducers, directly used tiled/window/fragment helpers, README GC-correction usage, and existing `gc-bias` tests in `tests/test_gc_bias.rs`. I did not run tests.

Shared findings that affect this command:

- G-003 in `00_shared_package_notes.md`: tiled temp directories need cleanup guards.
- G-005 in `00_shared_package_notes.md`: `--assign-by midpoint` uses shared non-reproducible even-fragment midpoint tie-breaking.
- G-008 in `00_shared_package_notes.md`: feature-gated QC plots are default command side effects.
- G-010 in `00_shared_package_notes.md`: GC correction packages cannot identify the sample or inputs they were built from.

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

### GB-003 - Medium - The documented correction clamp does not bound final package weights

The outlier help says correction values are clipped at `[0.1, 10.0]` after outlier detection ([config.rs](../../src/commands/gc_bias/config.rs#L288-L303)). In the implementation, the hard clamp is applied to the normalized bias ratio before the command re-normalizes each row and then inverts values into multiplicative weights ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L672-L693)). The existing clamp regression explicitly documents and expects a final package value of `10.49`, outside the nominal clamp range ([test_gc_bias.rs](../../tests/test_gc_bias.rs#L3499-L3513), [test_gc_bias.rs](../../tests/test_gc_bias.rs#L3559-L3565)). Also, the summary prints the hard-clamp count only when `outlier_method != none`, even though the hard clamp runs regardless ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L630-L685), [gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L746-L768)).

Impact: users who expect the written multiplicative weights to be bounded by `[0.1, 10.0]` can get larger or smaller final weights, and `--outlier-method none` can still clamp without reporting the clamp count.

Recommended fix:

- Decide the contract: clamp the pre-inversion bias ratio, or clamp the final multiplicative weights.
- If the current math is intended, update CLI help and the printed summary to say the safety clamp is pre-recentering/pre-inversion and final weights can move outside the nominal range.
- If the final weights should be bounded, apply a final finite positive clamp after inversion and re-check row centering expectations.

### GB-004 - Medium - The command does not validate that the written correction matrix is finite

The correction matrix starts with elementwise division of normalized cfDNA counts by normalized reference counts ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L640-L647)). Non-finite values can survive the rest of the pipeline: outlier handling ignores non-finite values, and the hard clamp uses ordinary comparisons, so `NaN` passes through unchanged ([outliers.rs](../../src/commands/gc_bias/outliers.rs#L154-L161), [outliers.rs](../../src/commands/gc_bias/outliers.rs#L195-L198), [gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L672-L684)). The package is then written without a final finite/positive check ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L687-L719)). Downstream correction sanitizes `NaN` and infinity as unusable weights, usually causing fragments to be skipped unless the user neutralizes invalid GC weights ([gc_tag.rs](../../src/shared/gc_tag.rs#L86-L113)).

Impact: a bad reference package, all-zero supported row, or another degenerate normalization case can produce a `.gc_bias_correction.npz` that loads successfully but later drops GC-corrected fragments. This should fail at package creation with a direct explanation.

Recommended fix:

- After the final inversion, `ensure!` every correction entry is finite and non-negative/positive according to the intended weight contract.
- Include row/column/bin-edge context in the error, so users can identify the problematic length/GC bin.
- Add a malformed or zero-supported reference-package regression that currently produces non-finite output and should instead fail before writing.

### GB-005 - Medium - Reference GC package structural validation is incomplete and can fail late

The loader validates scalar array lengths and that `counts`, support masks, and GC-percent widths share a shape ([load_reference_bias.rs](../../src/commands/gc_bias/load_reference_bias.rs#L69-L159)). It does not validate that `length_range = [min, max]` implies exactly `max - min + 1` rows, and it casts the stored `u32` `end_offset` to `u8` without checking range ([load_reference_bias.rs](../../src/commands/gc_bias/load_reference_bias.rs#L72-L80), [load_reference_bias.rs](../../src/commands/gc_bias/load_reference_bias.rs#L132-L140)). The run then uses this metadata to size its cfDNA `GCCounts` template before tiled counting ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L319-L324)).

Impact: malformed or stale `.ref_gc_package.npz` files can pass initial loading and fail only after tile counting, for example when GC-width correction discovers a shape mismatch. A large `end_offset` in a malformed package can also be silently truncated.

Recommended fix:

- Validate `min_fragment_length <= max_fragment_length`, row count equals the inclusive length range, and `end_offset <= u8::MAX` before constructing `ReferenceGCMetadata`.
- Validate the minimum effective length contract that `ref-gc-bias` enforces when it writes packages.
- Add loader-level regressions for row-count mismatch, inverted length range, and out-of-range end offset.

### GB-006 - Low - Empty BED tiles still open a BAM reader before skipping

The tile worker opens a fresh chromosome reader before preparing the tile-window state ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L784-L825)). `prepare_tile_windows()` can identify BED tiles with no candidate windows and return `skip_tile` before reference sequence reads or BAM fetches happen ([windows.rs](../../src/commands/gc_bias/windows.rs#L372-L448)), but the BAM reader has already been opened. The heavier reference sequence read and BAM fetch happen later ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L831-L847), [gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L880-L882)).

Impact: sparse BED runs avoid most per-tile work, but still pay one BAM-reader open per no-window tile. On very sparse targeted designs over whole-genome tiling, that can be avoidable overhead.

Recommended fix:

- For BED mode, prepare/skip windows before opening the BAM reader, using contig lengths already resolved in `run()`.
- Keep the current early return before reference sequence reads and `fetch()`.
- Add a small helper-level test or instrumentation-friendly regression proving no reader is opened for no-window BED tiles.

### GB-007 - Low - The CLI help example is not a clean runnable shell snippet

The command-level example in `GCConfig` includes blank lines inside the continued command and comments after line-continuation backslashes ([config.rs](../../src/commands/gc_bias/config.rs#L76-L88)). The README `gc-bias` example is cleaner, but the help text is likely what users see first from `cfdna gc-bias --help`.

Impact: users copying from CLI help can get shell-dependent behavior or a broken command, especially on the `--ref-2bit ... \ # Or some other assembly` line.

Recommended fix:

- Make the help example a single runnable snippet with no comments after `\`.
- Put optional notes before or after the code block.
- Keep the README and CLI help examples consistent.

## Existing coverage notes

The command already has broad coverage for default MAPQ, global/fixed/BED windows, no-count failure, reference end-offset propagation, overlapping/touching BED behavior, aligned vs misaligned single-chromosome tiling, aligned multi-chromosome accumulation, empty middle tiles, fixed-size vs BED cross-tile equivalence on one chromosome, real `ref-gc-bias` integration, saved intermediate sequence/content, minimum window ACGT filtering, outlier methods, hard-clamp behavior, and greedy binning.

The missing coverage from this review is multi-chromosome crossing-file reduction, prefix-safe intermediates, the final clamp/summary contract, finite correction-package validation, malformed reference-package metadata validation, avoiding BAM-reader opens for no-window BED tiles, and copy-pasteable CLI help.
