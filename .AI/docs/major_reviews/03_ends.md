# `cfdna ends` review

Date: 2026-04-24

Scope: `src/commands/ends/*`, the `ends` CLI configuration, and directly used shared helpers for end-aware fragment construction, BED/grouped windows, tile spans, sparse motif output, GC correction, scaling factors, blacklist checks, and midpoint assignment. I also read the existing `ends` tests in `tests/test_ends_command.rs` plus the module-level tests under `src/commands/ends/*_tests.rs`. I did not run tests.

Shared findings that affect this command:

- None active.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Post-release performance:

- G-006: sparse-window GC reference pruning.
- E-003: motif reference preload is full-tile even when BED windowing narrows the BAM fetch.

## Findings

### E-003 - Post-release performance - Motif reference preload is full-tile even when BED windowing narrows the BAM fetch

The GC-prefix part of this issue is now covered by G-006 in `00_shared_package_notes.md`. The ends-specific remaining issue is motif reference preparation: `motif_reference_span_for_tile()` always expands the tile fetch span, not the narrowed fetch span ([motifs.rs](../../src/commands/ends/motifs.rs#L153-L168)), and `build_tile_motif_context()` then reads that whole span and precomputes reference k-mer codes when reference bases are needed ([motifs.rs](../../src/commands/ends/motifs.rs#L239-L261)). The current tests pin this full-tile behavior ([motifs_tests.rs](../../src/commands/ends/motifs_tests.rs#L434-L485)).

Impact: targeted window runs that need reference motifs can spend a large fraction of runtime reading and encoding reference sequence that cannot contribute to output.

Recommended fix:

- Consider narrowing motif reference preload to the reachable assignment/motif span for the tile's relevant windows, while keeping enough padding for outside bases and raw-shifted soft clips.
- Add regression coverage with a sparse BED where no-window tiles do not request motif reference sequence.

## Existing coverage notes

The command already has broad coverage: read-backed and reference-backed motifs, inside-only and outside motifs, dense and sparse outputs, prefix handling, base-quality filters, unpaired mode, hard/soft clipping strategies, indel policy interactions, fragment and motif blacklisting, grouped BED aggregation, all window assignment modes, raw endpoint tile-boundary behavior, GC-file and GC-tag weighting, scaling factors, grouped metadata, settings formatting, and required-reference checks for motif extraction.

Motif reference narrowing for skipped sparse-window tiles and shared GC reference pruning are deferred performance optimizations tracked here and in G-006 in `00_shared_package_notes.md`.

## Released-command re-review additions (2026-05-04)

The shared raw-chromosome temporary filename issue (G-019), overlong `--gc-tag` issue (G-021), and unchecked output-prefix issue (G-022) originally noted here have since been implemented through shared temporary chromosome tokens, shared AUX-tag validation, and shared output-prefix validation.

No new `ends`-only pre-release correctness issue was found in this re-review. The previously tracked `ends` performance issue E-003 still stands, and the shared sparse-window GC preload issue remains tracked as G-006.

Additional coverage checked during this pass: the command now has explicit guards for zero-length motif definitions, base-quality filter incompatibilities, missing `--ref-2bit` for reference-backed motifs/outside bases/blacklists, raw-aligned clipping with reference-backed inside bases, and grouped BED inputs with no surviving selected-chromosome windows.

# Claude

Date: 2026-05-04

Scope: deep dive into `src/commands/ends/{ends,counting,motifs,output,tiling,write,config,config_structs,mod}.rs` plus the directly-used shared helpers (`shared/scale_genome.rs`, `shared/window_fetch.rs`, `shared/tiled_run.rs`, `shared/overlaps.rs`, `shared/bam.rs`, `shared/fragment/ends_fragment.rs`, `commands/cli_common.rs`). For comparison I cross-checked the analogous fragment-streaming loops in `lengths.rs` and `midpoints.rs`. I also scanned `tests/test_ends_command.rs` and the per-module `*_tests.rs` files for behaviour the current code is pinning. I did not run tests.

## Findings

## Items I checked and ruled out

These were investigated but I did not find evidence of an actual problem:

- **Window check order vs GC weight**: `ends` runs `find_overlapping_windows` *before* computing `gc_weight` ([ends.rs:830-888](../../src/commands/ends/ends.rs#L830-L888)) and explicitly comments why ([ends.rs:807-810](../../src/commands/ends/ends.rs#L807-L810)). This is the order I flagged as missing in `midpoints` (M-C-002); `ends` already gets it right.
- **Tile-core ownership vs blacklist counter**: `ends` checks `fragment.start() < tile.core_start()` *before* the blacklist filter ([ends.rs:790-805](../../src/commands/ends/ends.rs#L790-L805)), so `blacklisted_fragments` is incremented exactly once per fragment. This is the order I flagged as inverted in `midpoints` (M-C-001); `ends` is consistent with `lengths` here.
- **`scaling_with_bin_idx` placeholder index of 0**: same pattern as midpoints, and the inline comment at [ends.rs:690-697](../../src/commands/ends/ends.rs#L690-L697) explicitly documents that `find_overlapping_windows` returns the scan position rather than the carried `idx`. Verified against [overlaps.rs:269-300](../../src/shared/overlaps.rs#L269-L300).
- **`raw_shifted_gc_length_warning_issued` race**: the `load(Relaxed)` early-out followed by `swap(true, Relaxed)` correctly emits the warning at most once per run ([ends.rs:755-768](../../src/commands/ends/ends.rs#L755-L768)). Multiple threads can pass the load, but only one can observe `swap(...) == false`.
- **`right_endpoint_pos = boundary_pos.checked_sub(1).expect(...)`**: `boundary_pos` is the half-open right-end position, so it must be `> 0` for any non-empty fragment, including shifted-right boundaries (which only move outward). Safe.
- **`outside_kmer_clip_strategy_warning` firing for `Aligned`**: even though the aligner has placed the boundary, `Aligned` mode still exposes the user to the alignment's choice when soft clipping was present, which the strategy doc acknowledges ([config_structs.rs:300-306](../../src/commands/ends/config_structs.rs#L300-L306)). The warning is consistent with the docs.
- **`tile_window_span` right halo missing `max_soft_clips`**: this is purely an optimization knob — `find_overlapping_windows` advances `wd_ptr` through the full chromosome window slice and is not bounded by the cached span, so windows past the span end are still found. The actual correctness gap on the right shift comes from the BAM fetch halo (E-C-001), not from `tile_window_span`.

## Codex comments on Claude findings
