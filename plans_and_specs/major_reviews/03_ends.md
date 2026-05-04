# `cfdna ends` review

Date: 2026-04-24

Scope: `src/commands/ends/*`, the `ends` CLI configuration, and directly used shared helpers for end-aware fragment construction, BED/grouped windows, tile spans, sparse motif output, GC correction, scaling factors, blacklist checks, and midpoint assignment. I also read the existing `ends` tests in `tests/test_ends_command.rs` plus the module-level tests under `src/commands/ends/*_tests.rs`. I did not run tests.

Shared findings that affect this command:

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

Shared findings that affect this command:

- G-019 in `00_shared_package_notes.md`: `ends` writes per-tile temporary count files using the raw chromosome name in the filename ([ends.rs](../../src/commands/ends/ends.rs#L607-L612)).
- G-021 in `00_shared_package_notes.md`: `ends` passes `--gc-tag` directly as bytes to the end-fragment iterator ([ends.rs](../../src/commands/ends/ends.rs#L726-L735)), so overlong tag names follow the shared first-two-byte BAM AUX lookup behavior.
- G-022 in `00_shared_package_notes.md`: `ends` uses the raw output prefix for final outputs and sidecars ([write.rs](../../src/commands/ends/write.rs#L47-L64), [ends.rs](../../src/commands/ends/ends.rs#L460-L477)) and for the per-run temporary directory path ([ends.rs](../../src/commands/ends/ends.rs#L291-L294)).

No new `ends`-only pre-release correctness issue was found in this re-review. The previously tracked `ends` performance issue E-003 still stands, and the shared sparse-window GC preload issue remains tracked as G-006.

Additional coverage checked during this pass: the command now has explicit guards for zero-length motif definitions, base-quality filter incompatibilities, missing `--ref-2bit` for reference-backed motifs/outside bases/blacklists, raw-aligned clipping with reference-backed inside bases, and grouped BED inputs with no surviving selected-chromosome windows.

# Claude

Date: 2026-05-04

Scope: deep dive into `src/commands/ends/{ends,counting,motifs,output,tiling,write,config,config_structs,mod}.rs` plus the directly-used shared helpers (`shared/scale_genome.rs`, `shared/window_fetch.rs`, `shared/tiled_run.rs`, `shared/overlaps.rs`, `shared/bam.rs`, `shared/fragment/ends_fragment.rs`, `commands/cli_common.rs`). For comparison I cross-checked the analogous fragment-streaming loops in `lengths.rs` and `midpoints.rs`. I also scanned `tests/test_ends_command.rs` and the per-module `*_tests.rs` files for behaviour the current code is pinning. I did not run tests.

## Findings

### E-C-001 - Low - BAM fetch halo and tile halo do not account for `max_soft_clips` under `--clip-strategy raw-shifted-boundary`

In `RawShiftedBoundary` mode the assignment interval can extend up to `max_soft_clips` bp beyond the aligned reference span on either side ([config_structs.rs:316-327](../../src/commands/ends/config_structs.rs#L316-L327), [ends_fragment.rs:39-48](../../src/shared/fragment/ends_fragment.rs#L39-L48)). The tile build and the BED-narrowed BAM fetch in `ends` only use `max_fragment_length` for halos, with no `+max_soft_clips` term:

- [ends.rs:237](../../src/commands/ends/ends.rs#L237): `let halo_bp = opt.fragment_lengths.max_fragment_length;` — passed to `build_tiles`.
- [ends.rs:254](../../src/commands/ends/ends.rs#L254): `let tile_span_right_halo = opt.fragment_lengths.max_fragment_length as u64;`.
- [ends.rs:634](../../src/commands/ends/ends.rs#L634): `let bed_fetch_halo_bp = opt.fragment_lengths.max_fragment_length as u64;` — fed into `fetch_span_for_tile(... CandidateWindowExtent ...)`.

`lengths` solves the same problem with `configured_max_fragment_reach_bp = max_fragment_length + max_soft_clips` (or `+max_deletion_bases`) when clip-adjustment is on, and threads that value through both `build_tiles` and `fetch_span_for_tile` ([lengths.rs:87-102](../../src/commands/lengths/lengths.rs#L87-L102), [lengths.rs:277](../../src/commands/lengths/lengths.rs#L277), [lengths.rs:293-298](../../src/commands/lengths/lengths.rs#L293-L298), [lengths.rs:702](../../src/commands/lengths/lengths.rs#L702)).

Impact: when the leftmost relevant window in a tile sits more than `max_fragment_length` to the right of the tile core start (for example, the BED windows for this tile are all at the far right of the core), the BED-narrowed BAM fetch starts at `min_window - max_fragment_length`, which can be strictly greater than `tile.core_start()`. Owned fragments with aligned start in `[tile.core_start(), min_window - max_fragment_length)` are therefore not fetched at all. Without the shift those fragments cannot reach any window, so dropping them is correct. Under `RawShiftedBoundary`, however, their right boundary can shift up to `max_soft_clips` bp further right, so they can reach `min_window` via the shifted assignment interval — and they are silently missed.

Practically this is a narrow correctness gap (`max_soft_clips ≤ MAX_MAX_SOFT_CLIPS = 256` per [shared/constants.rs:11](../../src/shared/constants.rs#L11), so at most a 256 bp band of fragments at one edge of the affected tile core), but it is a real divergence from the equivalent `lengths` geometry and is unguarded by tests. The pinned `raw_shifted_boundary_endpoint_assignment_keeps_a_left_window_that_only_raw_reach_can_touch` and `..._tile_size_invariant` cases ([test_ends_command.rs:3834-3939](../../tests/test_ends_command.rs#L3834-L3939)) only exercise *left*-shifted reach into a window that lies before the tile core, where the fetch already has full left coverage.

Recommended fix: introduce an `ends`-side equivalent of `configured_max_fragment_reach_bp` that adds `max_soft_clips` when `clip_strategy.uses_shifted_boundary()`, and feed it into both `build_tiles` and `bed_fetch_halo_bp`. Add a regression test where the leftmost BED window for a tile is more than `max_fragment_length` to the right of `core_start()` and the only reachers are owned fragments whose right shift carries them into the window.

### E-C-002 - Low - Scaling-aware path silently drops candidate windows that overlap only the shifted assignment

For `--clip-strategy raw-shifted-boundary`, the candidate set is built from `fragment.assignment_interval` ([ends.rs:824](../../src/commands/ends/ends.rs#L824)). The no-scaling branch then iterates *every* candidate window and counts it ([ends.rs:962-982](../../src/commands/ends/ends.rs#L962-L982)). The scaling-enabled, non-`CountOverlap` branch instead routes through `compute_window_scaling_over_fragment` and passes `fragment.interval` (the aligned span, *not* the assignment span) as the "fragment" used for filtering ([ends.rs:928-934](../../src/commands/ends/ends.rs#L928-L934)). That helper then drops any candidate window whose start/end does not overlap the aligned span:

```rust
for window in &count_overlaps.windows {
    if window.end() > fragment_start_bp && window.start() < fragment_end_bp {
        per_window_scaling.push((window.idx, avg_over_fragment, 1.0));
    }
}
```

at [scale_genome.rs:209-216](../../src/shared/scale_genome.rs#L209-L216).

Impact: with `RawShiftedBoundary` enabled and `--scaling-factors` provided (and any assigner other than `CountOverlap`), candidate windows that overlap only the *shifted* portion of the assignment — exactly the windows that pinned tests like `raw_shifted_boundary_endpoint_assignment_keeps_a_left_window_that_only_raw_reach_can_touch` rely on without scaling — are silently filtered out. The same input therefore produces fewer counted windows when scaling is added than when it is not. The `CountOverlap` branch escapes this because `compute_window_scaling_over_overlap` re-uses `count_overlaps.query_start()/query_end()` ([scale_genome.rs:149-150](../../src/shared/scale_genome.rs#L149-L150)), which carries the assignment-interval start/end through `OverlappingWindows.interval`.

There is no test covering scaling combined with `RawShiftedBoundary` (`grep` for `set_scaling_factors` near `RawShifted*` in `tests/test_ends_command.rs` returns no matches), so no test pins either contract.

Recommended fix: pick one explicit semantics and document it. Either pass `fragment.assignment_interval` to `compute_window_scaling_over_fragment` so scaling matches the candidate-window geometry, or apply the same aligned-span filter in the no-scaling branch (and document that shifted-only overlaps are not counted at all). The current behaviour is "scaling silently changes the set of counted windows," which is the worst of the three options.

### E-C-003 - Low - `gc_failed_fragments` is incremented in both arms of the file-based GC `if/else`

The file-based GC branch increments `gc_failed_fragments` separately inside the `neutralize_invalid_gc` and `else` arms ([ends.rs:875-883](../../src/commands/ends/ends.rs#L875-L883)):

```rust
(None, true) => {
    if opt.gc.neutralize_invalid_gc {
        counter.gc_failed_fragments += 1;
        1.0
    } else {
        counter.gc_failed_fragments += 1;
        continue;
    }
}
```

The gc-tag branch right above it follows the cleaner pattern of incrementing once *before* the `if`, which is also what `midpoints` does ([ends.rs:850-857](../../src/commands/ends/ends.rs#L850-L857)). The current shape is behaviourally equivalent today but makes it easy for one of the two `+= 1` updates to drift out of step in a future edit.

Recommended fix: hoist `counter.gc_failed_fragments += 1;` above the `if opt.gc.neutralize_invalid_gc` so both branches share a single increment, mirroring the gc-tag path.

### E-C-004 - Low - `overlapping_window_intervals` is rebuilt as an FxHashMap per fragment in the scaling path

For every fragment with a non-empty scaling map and at least one candidate window, the scaling branch builds a fresh `FxHashMap<usize, Interval<u64>>` keyed by chromosome-local window index just to look up `window_interval` inside the inner loop ([ends.rs:938-961](../../src/commands/ends/ends.rs#L938-L961)). Both `compute_window_scaling_over_fragment` ([scale_genome.rs:210-216](../../src/shared/scale_genome.rs#L210-L216)) and `compute_window_scaling_over_overlap` ([scale_genome.rs:155-174](../../src/shared/scale_genome.rs#L155-L174)) emit their tuples in the same order as `count_overlaps.windows`, just possibly with some entries skipped, so the lookup table is mostly redundant.

Impact: a hash-map allocation and population on every counted, scaled fragment in the hot path. Not a correctness issue, but it scales linearly with fragment count.

Recommended fix: have the scaling helpers carry the window's `Interval<u64>` (or the original `OverlappingWindow` reference) alongside `(idx, weight, fraction)`, or change `ends` to walk `count_overlaps.windows` and `overlap_weights` together in lockstep (skipping the same entries that the helper skipped). Either avoids the per-fragment hash map.

### E-C-005 - Informational - `max_soft_clips.try_into::<u32>` cannot fail

[ends.rs:711-715](../../src/commands/ends/ends.rs#L711-L715):

```rust
let max_soft_clips: u32 = opt
    .clip
    .max_soft_clips
    .try_into()
    .context("max_soft_clips does not fit in u32 for fragment iteration")?;
```

`max_soft_clips` is `u16` ([config_structs.rs:347-356](../../src/commands/ends/config_structs.rs#L347-L356)), and `u16::try_into::<u32>` is statically infallible, so the `.context(...)` is unreachable. A plain `u32::from(opt.clip.max_soft_clips)` documents the conversion and removes the misleading error path.

### E-C-006 - Informational - `_tid_check == tile.tid as u32` wraps on negative tids

[ends.rs:600-601](../../src/commands/ends/ends.rs#L600-L601) compares the BAM-resolved tid (`u32` per [bam.rs:10](../../src/shared/bam.rs#L10)) to `tile.tid as u32`, where `tile.tid: i32` ([tiled_run.rs:21](../../src/shared/tiled_run.rs#L21)). For a negative `tile.tid` the cast wraps and the assert can pass spuriously. Same defensive concern as flagged in `02_midpoints.md` (M-C-008); not a real bug today because tids are never negative for valid contigs, but worth normalizing across commands.

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

### E-C-001 comment

I do not think E-C-001 is a confirmed bug as written. The argument misses that `ends` filters fragments by `assignment_len()` after clip-strategy boundary resolution ([ends.rs](../../src/commands/ends/ends.rs#L703-L707), [ends_fragment.rs](../../src/shared/fragment/ends_fragment.rs#L351-L359)). In `RawShiftedBoundary`, `assignment_len = aligned_span + left_soft_clip + right_soft_clip`; therefore a kept fragment's shifted assignment end is still at most `aligned_start + max_fragment_length`.

Given a BED-narrowed fetch start of `min_window - max_fragment_length`, an owned fragment with `aligned_start < min_window - max_fragment_length` cannot have `assignment_end > min_window` without violating the `assignment_len <= max_fragment_length` filter. So the proposed counterexample does not survive the current geometry.

The existing tile-window span setup also does account for shifted reach in the two directions that matter for owned fragments: left reach is `max_soft_clips`, and right reach is `max_fragment_length` ([ends.rs](../../src/commands/ends/ends.rs#L249-L267)). That is different from `lengths`, where adjusted output length can hide additional aligned reference reach from deletions or clip adjustment, so `lengths` needs a broader `configured_max_fragment_reach_bp` ([lengths.rs](../../src/commands/lengths/lengths.rs#L81-L102)).

Tile ownership is also central here: `ends` owns a fragment by its aligned `fragment.start()` being inside the tile core before assignment-window counting begins ([ends.rs](../../src/commands/ends/ends.rs#L786-L792)). The right-reach question is therefore not "how far can a shifted right boundary move from the aligned end?", but "how far can a kept assignment interval extend from an owned aligned start?" That distance is bounded by `max_fragment_length`.

I would not add `max_soft_clips` to the `ends` BAM fetch halo without a failing example that respects both tile ownership and `assignment_len()`. At most, add a regression or code comment that pins this invariant: `max_fragment_length` bounds the right shifted assignment reach from an owned aligned start.

### E-C-002 comment

I agree this one looks real. Candidate windows are found from `fragment.assignment_interval` ([ends.rs](../../src/commands/ends/ends.rs#L810-L838)), but the scaling-enabled non-`CountOverlap` path asks `compute_window_scaling_over_fragment()` to emit only windows overlapping `fragment.interval`, the aligned span ([ends.rs](../../src/commands/ends/ends.rs#L922-L934), [scale_genome.rs](../../src/shared/scale_genome.rs#L189-L216)). That means a shifted-only endpoint/window can be counted without scaling but skipped when scaling factors are present.

The fix should preserve the documented scaling span: file-based GC and scaling weighting still use the aligned reference span ([config_structs.rs](../../src/commands/ends/config_structs.rs#L321-L327)). So I would not pass the shifted assignment interval as the span to average over. Better: compute the average scaling over `fragment.interval`, but emit that same value for every already-selected candidate window in non-`CountOverlap` modes. Add a regression combining `RawShiftedBoundary`, a scaling TSV, and a shifted-only endpoint/window.

### E-C-003 comment

True, but not a correctness issue. Both arms increment exactly once today ([ends.rs](../../src/commands/ends/ends.rs#L873-L883)). Hoisting the increment would reduce drift risk and match the GC-tag branch, but this should be treated as cleanup.

### E-C-004 comment

This is a plausible hot-path allocation cleanup, not a correctness finding. If optimized, prefer changing the scaling helpers to carry each selected window interval with the emitted tuple. A lockstep walk is easy to get subtly wrong because the helpers may skip windows.

### E-C-005 comment

Confirmed cleanup only. `max_soft_clips` is `u16`, so `u32::from(opt.clip.max_soft_clips)` is the clearer conversion ([config_structs.rs](../../src/commands/ends/config_structs.rs#L340-L356), [ends.rs](../../src/commands/ends/ends.rs#L711-L715)).

### E-C-006 comment

Confirmed as defensive consistency cleanup, not a runtime bug for valid tiles. If touched, use `u32::try_from(tile.tid)` before comparing to the BAM-resolved `u32` tid ([ends.rs](../../src/commands/ends/ends.rs#L599-L601), [tiled_run.rs](../../src/shared/tiled_run.rs#L17-L25)).
