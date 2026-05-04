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

- G-022: `gc-bias` uses the shared unchecked output prefix for the correction package, optional intermediates, and plot files.
- G-023: `gc-bias` consumes reference GC packages but currently checks only chromosome-name set equality against the current run.

Reviewed shared findings that do not currently apply:

- G-019: `gc-bias` crossing-window temp files are named by tile index, not raw chromosome name.
- G-021: `gc-bias` does not expose `--gc-tag`.

### Release triage additions

Pre-release correctness/safety:

- G-023: wrong reference GC package can be accepted and repackaged under the current 2bit signature.
- GB-007: fixed-size window counting can miss fragment contributions when valid window sizes are smaller than fragment spans.
- G-022: unchecked output prefixes can write files outside the requested output directory.

Post-release performance:

- GB-006 remains a post-release performance item; no new correctness issue was found in that path.

### GB-007 - Medium - Fixed-size windowing can miss windows when fragments span more than two bins

The fixed-size `gc-bias` path uses a streaming implementation with only two live window buffers: the current fixed window and the next fixed window ([windows.rs](../../src/commands/gc_bias/windows.rs#L196-L263), [windows.rs](../../src/commands/gc_bias/windows.rs#L425-L432)). During fragment counting, it computes overlap against only those two buffers and can increment only `current` and `next` ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L965-L1084)).

That assumption is not enforced by the CLI. `--by-size` accepts any `u64` and defaults to 100000 only when absent ([cli_common.rs](../../src/commands/cli_common.rs#L343-L405)). The default assignment is `count-overlap` ([cli_common.rs](../../src/commands/cli_common.rs#L456-L489)), so a valid command such as a small-bin `gc-bias --by-size 50` run can have fragments that overlap three or more fixed windows. In that case, fixed-size counting silently drops the third and later window contributions. The BED/global path does not share this limitation because it calls `find_overlapping_windows()` and iterates every returned overlap ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L1146-L1185)).

Impact: default 100 kb windows are not affected for the current fragment-length defaults, but user-requested small fixed windows can produce undercounted sample GC tables and therefore wrong correction factors.

Recommended fix:

- Keep the two-buffer fixed-size streaming branch. That design is the right performance invariant for this command; small `gc-bias` windows are not an important use case.
- Add a validation guard before counting that rejects fixed window sizes too small to preserve the two-buffer invariant. A conservative rule would require `--by-size >= max_fragment_length`; an exact coordinate proof could allow the tight half-open boundary case, but it should still guarantee that any assignment interval can overlap only the current and next window.
- Add a focused validation test proving an invalid small `--by-size` run fails fast with a direct message, instead of expanding the streaming implementation for an edge case.
