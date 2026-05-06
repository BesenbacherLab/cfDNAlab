# `cfdna midpoints` review

Date: 2026-04-24

Scope: `src/commands/midpoints/*`, the CLI dispatch for `midpoints`, and directly used shared helpers for BED loading, tiling, overlap lookup, GC correction, scaling, blacklist checks, midpoint placement, and grouped count merging. I also read the existing midpoint-focused tests in `tests/test_profile_groups_command.rs`, `tests/test_profile_groups_counts.rs`, `tests/test_tiling.rs`, `tests/test_cli_smoke.rs`, and cross-command roundtrip/artifact tests. I did not run tests.

Shared findings that affect this command:

- None active.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release semantic/docs:

- M-005: midpoint blacklist semantics should be either aligned with counted midpoint placement or documented explicitly.

Post-release performance/scalability:

- G-006: sparse-window GC reference pruning.
- M-001B: optional sparse final output would help very large group, width, and length-bin shapes.

## Findings

### M-001B - Medium - Very large final profile shapes need an optional sparse output format

The final output is intentionally dense today: `<prefix>.midpoint_profiles.npy` stores counts shaped `(group, length_bin, position)` after collapsing all intervals in each group. Dense output is straightforward for downstream numpy workflows and should remain the default for ordinary profile sizes.

Impact: dense final output becomes impractical when users request high positional resolution, many groups, and near-base-pair fragment length bins. For example, `2001 * 3000 * 971` cells is about 5.83 billion `f32` values, or roughly 21.7 GiB for the final counts alone. Sparse internal tile partials solve per-tile memory and temp disk pressure, but they do not solve the final artifact size when the requested output tensor itself is too large.

Recommended tasks:

- Add an explicit sparse final output option rather than changing the default dense artifact.
- Choose and document a stable sparse layout. A practical first version is a SciPy-compatible sparse `.npz` with rows representing `group_idx * num_length_bins + length_bin_idx` and columns representing `position`, plus sidecar metadata for group names and length-bin edges.
- Preserve zero-count groups in the sparse matrix shape so group indices remain stable even when a group has no observed counts.
- Keep dense plotting tied to dense output, or require densifying only selected groups for plotting. Do not accidentally densify the full sparse output just to make plots.
- Make the CLI estimate compare dense final size with expected sparse payload size when possible. If expected density is unknown before counting, report the dense final size and explain that sparse output avoids storing zeros in the final artifact.
- Add roundtrip tests that load the sparse final output, densify it in test code, and compare it against the dense `.npy` output for the same small fixture.
- Document downstream loading in Python, including how to reconstruct `(group, length_bin, position)` from the sparse two-dimensional representation.

### M-005 - Medium - `--blacklist-strategy midpoint` uses a different even-fragment midpoint than counting

Blacklist filtering happens before the command samples the counted midpoint ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L495-L510)). The shared blacklist midpoint strategy uses `start + (end - start) / 2`, which is the right-center base for even-length half-open fragments ([overlaps.rs](../../src/shared/blacklist/overlaps.rs#L74-L96)). The count placement then randomly selects left or right center for the same even-length fragment ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L508-L510)).

The existing test suite pins this current behavior: an even fragment can be deterministically blacklisted by the right-center base even though the profile placement would otherwise be randomized between two central bases ([test_profile_groups_command.rs](../../tests/test_profile_groups_command.rs#L586-L662)).

Impact: users may read `midpoint` blacklist strategy as "filter by the midpoint that would be counted", but that is not true for even-length fragments at blacklist edges.

Recommended fix:

- Decide whether `midpoints` should reuse the sampled counted midpoint for midpoint-blacklist filtering.
- If the deterministic blacklist midpoint is preferred, state that explicitly in the `--blacklist-strategy midpoint` help for commands that also randomize midpoint placement.

## Sparse final output design note

The urgency of M-001B depends heavily on the intended axis sizes. In one recent benchmark, `--length-bins 30:1001:50` produced a dense midpoint profile file of about 191 MB. That range resolves to 20 length bins. A full 1 bp length axis for 80-600 bp inclusive would use `--length-bins 80:601:1`, which resolves to 521 bins. Holding groups and profile width fixed, that scales the dense output by about `521 / 20 = 26.05x`, or roughly 5.0 GB for the same sample. That is large but still manageable on cluster storage and memory for many single-sample workflows.

Sparse final output is therefore less urgent for one sample at this shape, but it still matters for cohort workflows with thousands of samples side by side, higher group counts, wider windows, or full 30-1000 bp resolution. Dense should remain the default because it is easy to load and slice in NumPy. A sparse option should optimize downstream ergonomics, not only file size.

Recommended sparse-output shape:

- Add an explicit option such as `--save-as-sparse`.
- Write a SciPy-compatible CSR `.npz`, e.g. `<prefix>.midpoint_profiles.sparse.npz`.
- Store rows as `group_idx * num_length_bins + length_bin_idx` and columns as `position`.
- Preserve empty groups and empty length bins in the matrix shape so group and bin indices remain stable.
- Keep sidecar metadata for group names and length-bin edges, likely reusing `<prefix>.group_index.tsv` and adding a small length-bin sidecar if needed.
- Document Python loading with `scipy.sparse.load_npz`, plus reconstruction of `(group, length_bin, position)` from sparse row and column indices.
- Keep plotting dense-only by default, or densify only the selected rows needed for plotting. Do not densify the full sparse matrix just to make plots.

# Claude

Date: 2026-05-03

Scope: deep dive into `src/commands/midpoints/{midpoints,counting_by_group,config,windows,plotting,mod}.rs` plus the directly-used shared helpers (`shared/midpoint.rs`, `shared/overlaps.rs`, `shared/scale_genome.rs`, `shared/window_fetch.rs`, `shared/tiled_run.rs`, `shared/length_axis.rs`, `shared/bed.rs`, `commands/cli_common.rs`). For comparison I cross-checked the analogous fragment-streaming loops in `lengths.rs` and `ends.rs`. I did not run tests.

## Findings

### M-C-001 - Low - Blacklist counter is incremented before tile-core ownership check, unlike `lengths`/`ends`

In `midpoints.rs` the blacklist filter runs **before** the midpoint-in-core ownership check ([midpoints.rs:548-569](../../src/commands/midpoints/midpoints.rs#L548-L569)):

```
is_blacklisted(...) -> bl_count += 1; continue   // line 548-559
midpoint = ...                                   // line 563-564
if midpoint < core_start || >= core_end { continue }  // line 567-569
```

The sister commands invert this order: tile-core check first, then blacklist ([lengths.rs:914-931](../../src/commands/lengths/lengths.rs#L914-L931), [ends.rs:790-805](../../src/commands/ends/ends.rs#L790-L805)).

Impact: because tile fetch bands overlap by the pairing halo, the same fragment can be visited by two adjacent tile workers. In `lengths`/`ends` the second tile drops the fragment at the core check, so `blacklisted_fragments` is incremented exactly once per blacklisted fragment. In `midpoints` both tiles increment the counter on boundary fragments before any ownership check, inflating the reported "Blacklist-excluded fragments" by a small amount near tile borders. Actual midpoint counts are unaffected (they sit behind the in-core filter at line 567).

This is partially absorbed by `TILE_DOUBLE_COUNT_NOTE`, but the note is a generic disclaimer and the divergence from `lengths`/`ends` is avoidable.

Recommended fix: move the midpoint computation and core-ownership check above the blacklist call so the loop body order matches `lengths` and `ends`. The existing comment "Must come after filtering for midpoints lying within the core!" at [midpoints.rs:572](../../src/commands/midpoints/midpoints.rs#L572) (about GC correction) already implies this ordering principle for stat counters, so applying it to `blacklisted_fragments` is the natural extension.

### M-C-002 - Low - `gc_failed_fragments` is incremented for fragments whose midpoint hits no count window

The current order inside the fragment loop is:

1. Compute `gc_weight` (and possibly `gc_failed_fragments += 1` / `gc_missing_tags += 1` / `gc_out_of_range_tags += 1`) at [midpoints.rs:573-613](../../src/commands/midpoints/midpoints.rs#L573-L613).
2. Call `find_overlapping_windows` for the 1bp midpoint at [midpoints.rs:616-624](../../src/commands/midpoints/midpoints.rs#L616-L624).
3. If no window is hit, `continue` at line 627-629.

This means a fragment whose midpoint sits in the tile core but does not overlap any grouped BED window still pays for GC correction *and* is permanently recorded as a GC failure if its tag/file lookup fails. Those fragments would never have contributed counts.

Impact: small inflation of `gc_failed_fragments`, `gc_missing_tags`, and `gc_out_of_range_tags` for sparse BED inputs (the use case where window coverage of the tile core is partial), plus avoidable per-fragment GC work when the midpoint falls outside every BED window. `ends.rs` already documents and applies the opposite ordering to "all later work only happens for windows that can actually receive counts" ([ends.rs:807-810](../../src/commands/ends/ends.rs#L807-L810)).

Recommended fix: move the `find_overlapping_windows` call above the GC weight block so the GC computation only runs when the midpoint actually has a candidate count window. The comment at [midpoints.rs:572](../../src/commands/midpoints/midpoints.rs#L572) should be tightened to say "after filtering for midpoints lying within the core *and* overlapping at least one count window".

### M-C-003 - Low - `get_overlapping_sites_and_adapt_fetch_to_extremes` walks tile windows twice

The helper calls `overlapping_windows_for_tile` once to materialize `overlapping_sites` and then again, indirectly, inside `window_derived_fetch_extent_for_core_overlap` which itself calls `overlapping_windows_for_tile` to find the min/max window edges ([midpoints.rs:773-789](../../src/commands/midpoints/midpoints.rs#L773-L789), [window_fetch.rs:171-189](../../src/shared/window_fetch.rs#L171-L189)).

Impact: per-tile overhead for re-iterating the candidate window slice. Small per-tile, but multiplied by the tile count it is pure waste because the second pass produces information that the first pass could have computed inline.

Recommended fix: fold the min-start / max-end accumulation into the single pass that builds `overlapping_sites`, and drop the second `window_derived_fetch_extent_for_core_overlap` call (or specialize a variant that accepts a pre-collected slice).

### M-C-004 - Low - `--plot-groups` validation runs after the dense `.npy` and TSV are already on disk

`plot_midpoint_profiles` is the only place that checks whether each requested plot index exists in `group_idx_to_name` ([plotting.rs:47-57](../../src/commands/midpoints/plotting.rs#L47-L57)). It is invoked at [midpoints.rs:301-318](../../src/commands/midpoints/midpoints.rs#L301-L318), after `write_npy` and `write_group_idx_to_name_tsv` have already written `<prefix>.midpoint_profiles.npy` and `<prefix>.group_index.tsv`.

Impact: a typo in `--plot-groups` (e.g. `--plot-groups 0,17` when only 5 groups exist) currently aborts the run only after the entire counting pipeline finishes and the public outputs are already written. The user gets an error and no `print_fragment_run_statistics` summary, but the artifacts are present. This is recoverable but jarring.

Recommended fix: validate `plot_groups` against `group_idx_to_name` (and against `num_groups`) right after the BED has been loaded — the same place that already logs `Num. groups: ...`. Behaviorally equivalent for valid inputs, fail-fast for invalid ones.

### M-C-005 - Low - `0.99` midpoint overlap threshold and its rationale comment

[midpoints.rs:622](../../src/commands/midpoints/midpoints.rs#L622) passes `0.99` as `min_overlap_fraction` for a 1bp midpoint query and labels it `// "Full" 1bp overlap but avoid rounding error`. For a 1bp query, `fraction_overlap_of_a` always evaluates to exactly `0.0` or exactly `1.0` in IEEE-754 (`overlap_bp` is 0 or 1, divided by an `interval_a_len` of 1) — there is no rounding to defend against ([overlaps.rs:336-340](../../src/shared/overlaps.rs#L336-L340)).

The numerical effect of `0.99` here is correct, but the comment is misleading and the constant is inconsistent with the analogous `Midpoint` branch in `ends.rs` which uses `1.0 - 1.0/(max_fragment_length + 1.0)` ([ends.rs:684-687](../../src/commands/ends/ends.rs#L684-L687)).

Impact: documentation/maintenance only; behavior is unchanged.

Recommended fix: either drop the misleading rationale and use `1.0` directly for a 1bp query, or align with the formula used in `ends.rs` so the two commands explain the threshold the same way.

### M-C-006 - Low - Stale "profile-groups" reference in `run` doc comment

The `run` doc says the parameter is "Fully resolved configuration for the `profile-groups` command" ([midpoints.rs:76-77](../../src/commands/midpoints/midpoints.rs#L76-L77)), but the command was renamed to `midpoints` (the file, the module path `commands::midpoints::config::MidpointsConfig`, and `COMMAND_TARGET = "midpoints"` all reflect the current name).

Impact: documentation only. Worth fixing while touching this file.

### M-C-007 - Low - Stray `// Handle deletions?` TODO at top of `midpoints.rs`

[midpoints.rs:56](../../src/commands/midpoints/midpoints.rs#L56) carries an isolated `// Handle deletions?` comment with no surrounding code. It is unclear whether the author intended this as a question about CIGAR D operations affecting `Fragment.interval` (which currently uses the aligned reference span and includes deletions) or about something else.

Impact: ambiguous follow-up that future readers may misinterpret. If it is a real open question about midpoint semantics for fragments containing reference deletions, it deserves a tracked review entry; otherwise it should be removed.

Recommended fix: resolve to either an explicit review item with concrete semantics, or delete the comment.

### M-C-008 - Informational - `_tid_check == tile.tid as u32` silently wraps on negative tids

[midpoints.rs:416](../../src/commands/midpoints/midpoints.rs#L416) compares the BAM-resolved tid (`u32`) to `tile.tid as u32`, where `tile.tid: i32` per [tiled_run.rs:21](../../src/shared/tiled_run.rs#L21). For a negative `tile.tid`, the cast wraps and the assert can pass spuriously. In practice, tids for valid chromosomes are non-negative, so this is not currently exploitable, but the assertion is weaker than it appears.

Impact: defensive sanity check only; flagged for completeness.

Recommended fix: use `u32::try_from(tile.tid)` and assert equality on the unwrapped value, or change the field type if no command ever uses negative tids.

## Items I checked and ruled out

These were investigated but I did not find evidence of an actual problem:

- **`scaling_with_bin_idx` placeholder index of 0 plus `OverlappingWindow.idx` reuse**: works correctly because `find_overlapping_windows` returns the scan position (not the carried `idx`) in BED mode ([overlaps.rs:269-300](../../src/shared/overlaps.rs#L269-L300)), and the scaling vector is built in the same order as `scaling_chr`. The contract is documented in the inline comment at [midpoints.rs:438-444](../../src/commands/midpoints/midpoints.rs#L438-L444) and matches the same pattern used in `ends.rs` and `lengths.rs`.
- **`wd_ptr` and midpoints going backward across consecutive fragments**: `find_overlapping_windows` uses `query.start.saturating_sub(look_back)` with `look_back = max_fragment_length`, which strictly upper-bounds the maximum decrease in midpoint between two start-sorted fragments (`(max_len - min_len)/2 ≤ max_len`). No off-by-one in the streaming pointer.
- **`compute_per_window_scaling_over_fragment` filter pruning midpoint-overlap windows**: cannot occur. Any window returned by the midpoint query satisfies `window.end > fragment_start` and `window.start < fragment_end` because the midpoint lies inside both intervals.
- **Sparse merge correctness**: `read_sparse_profile_partial_file` validates length agreement, ascending sort, `usize` fit, in-bounds indices, and shape match before any merge work begins ([counting_by_group.rs:553-620](../../src/commands/midpoints/counting_by_group.rs#L553-L620)). The chunked-Mutex merge in `merge_sparse_profile_partial_file` correctly partitions by `flat_idx / chunk_size` and only holds one chunk lock at a time, with the wrap-around split via `partition_point` preserving full coverage of every entry.
- **`gc_corrector.clone()` per tile**: clones an `Array2<f64>`, which is non-trivial, but the same pattern is used in `ends`, `lengths`, and `fcoverage` (see explicit comment "Quite small memory footprint" at [fcoverage.rs:562](../../src/commands/fcoverage/fcoverage.rs#L562)). Not a midpoints-specific concern.
- **`incr_weighted` casting `f64 -> f32`**: realistic GC × scaling products stay far below `f32::MAX`; no overflow in practice and consistent with the dense `.npy` output type.
- **Even-fragment blacklist midpoint vs counted midpoint divergence**: already tracked as M-005 above; not duplicated here.

## Codex comments on Claude findings

Overall: the Claude findings are mostly low-severity cleanup and statistics-consistency issues rather than output-count correctness issues. The most useful near-term items are M-C-001, M-C-002, M-C-004, and M-C-006 because they make stats and failure behavior easier to reason about.

- **M-C-001**: Agree. Moving midpoint-core ownership before blacklist would make `blacklisted_fragments` less tile-size-sensitive and more consistent with `lengths` and `ends`. This should not change final midpoint counts, only reported statistics and possibly a little wasted blacklist work. It needs a focused stats test near a tile boundary.
- **M-C-002**: Agree in principle. GC failure counters should probably reflect fragments that could have contributed to a count. Moving midpoint-window overlap before GC also avoids unnecessary GC work for sparse BED inputs. The main implementation detail is to keep the existing ordering invariant that file-based GC correction only happens after tile-core ownership is known.
- **M-C-003**: Agree, but this is a small optimization. If changed, keep the helper readable because the fetch-narrowing semantics are more important than saving one short per-tile pass.
- **M-C-004**: Agree. This is a straightforward fail-fast improvement. Plot group validation should happen immediately after grouped BED loading, while `num_groups` and `group_idx_to_name` are already available.
- **M-C-005**: Mostly agree. For a 1 bp midpoint query, `1.0` is clearer than `0.99`, but this should not be framed as a general "no roundoff" argument. Float-heavy counting code always has numerical noise. The practical point is narrower: a 1 bp query produces only complete-overlap or no-overlap cases, so `1.0` expresses the intended predicate better than `0.99`. If we keep a non-1.0 threshold for consistency with other commands, the comment should describe command consistency rather than floating-point safety.
- **M-C-006**: Agree. This is just stale wording and should be fixed when touching `midpoints.rs`.
- **M-C-007**: Agree that the comment is too vague. If deletion/CIGAR semantics matter, they need a real review item. Otherwise the comment should go away in a cleanup pass.
- **M-C-008**: Agree as defensive cleanup, but it is not operationally important if all real `Tile` values come from contig tids. A `try_from` assert would make the invariant explicit.
- **Items ruled out**: Mostly agree. The note about `wd_ptr` only rules out backward movement across fragment order inside one window slice. It does not cover the separate tile-local versus chromosome-wide pointer bug found during the tile-size comparison, where `core_overlapping_windows` is a compact vector and must use a tile-local pointer starting at zero.

## Existing coverage notes

The command already has broad integration coverage: length-bin parsing, default MAPQ, paired/unpaired parity, group index ordering, even-length midpoint edge placement, blacklist midpoint behavior, real GC packages, GC tags, GC/scaling multiplication, scaling TSV validation, tile-boundary behavior, chromosome-end fetch narrowing, CLI smoke output, and cross-command BAM/fragment roundtrips are all represented.

The most important midpoint-specific missing tests from this review are sparse final-output roundtrips if M-001B is implemented and the chosen even-midpoint blacklist semantics if that contract changes. The deferred skipped-tile GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.

## Re-review additions (2026-05-04)

The shared raw-chromosome temporary filename issue (G-019) and overlong `--gc-tag` issue (G-021) originally noted here have since been implemented. `midpoints` now uses shared temporary chromosome tokens and shared AUX-tag validation.

### Release triage additions

Pre-release correctness/safety:

- None active from this re-review.

Post-release performance:

- G-006 remains the only midpoint-specific performance item from this pass.

### Command-specific findings

No new midpoints-only counting correctness finding was added in this pass. The shared input/output contract issues found during the re-review have since been implemented.
