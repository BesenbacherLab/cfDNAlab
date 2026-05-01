# `cfdna lengths` review

Date: 2026-04-24

Scope: `src/commands/lengths/*`, the `lengths` CLI configuration, and directly used shared helpers for indel/clip-aware fragment construction, BED/grouped windows, tile spans, GC correction, scaling factors, blacklist checks, midpoint assignment, and length-count reduction. I also read the existing `lengths` tests in `tests/test_lengths_command.rs`. I did not run tests.

Shared findings that affect this command:

- G-002 in `00_shared_package_notes.md`: README OPTIONS blocks need clearer alternative-choice labeling.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release correctness/safety:

- L-001: indel-adjusted output length and aligned fetch span share one max-length setting.
- L-003: settings sidecar is too thin to interpret the matrix.
- L-004: BED row reordering should fail loudly if metadata and counts drift apart.

Pre-release docs/API polish:

- G-002: README OPTIONS blocks should keep their current structure but clarify alternative choices.

Post-release performance:

- G-006: sparse-window GC reference pruning.

## Findings

### L-001 - High - Indel-adjusted output length and aligned fetch span share one max-length setting

The command help says that when indel or clip adjustment is enabled, fragment length filtering is based on the adjusted length ([config.rs](../../src/commands/lengths/config.rs#L28-L30)). The implementation does that in the fragment filter by calling `fragment.adjusted_len(indel_mode, clip_mode)` ([lengths.rs](../../src/commands/lengths/lengths.rs#L793-L807)). For `IndelMode::Adjust`, adjusted length starts from the aligned reference span and subtracts deletion-like events ([indel_counting_fragment.rs](../../src/shared/fragment/indel_counting_fragment.rs#L124-L140)).

The adjusted length does not directly set fetch coordinates; fetching is still aligned-coordinate based. The issue is that `lengths` uses the same configured `max_fragment_length` for two different meanings: the maximum output length after adjustment, and the aligned-coordinate halo used for tile construction and BED fetch narrowing. `aligned_fetch_halo_bp = opt.fragment_lengths.max_fragment_length` is passed into `build_tiles()` ([lengths.rs](../../src/commands/lengths/lengths.rs#L215-L233)), and BED fetch narrowing also uses `max_fragment_length` as the aligned halo ([lengths.rs](../../src/commands/lengths/lengths.rs#L650-L668)). The paired-end iterator only emits a fragment after both mates have been seen in the fetched records ([core.rs](../../src/shared/fragment_iterators/core.rs#L300-L320)).

Impact: with `--indel-mode adjust`, a fragment containing a long deletion or skipped region can have an adjusted length inside `[min, max]` while its aligned reference span is larger than `max_fragment_length`. If a mate lands outside the aligned fetch halo, the pair is not emitted in the owning tile. If the pair is emitted but its aligned interval extends beyond preloaded reference sequence, GC correction can fail or skip it. This is most plausible with `N`/large deletion CIGAR operations, which the help explicitly treats as deletions ([config.rs](../../src/commands/lengths/config.rs#L104-L107)).

Recommended fix:

- Separate "maximum adjusted output length" from "maximum aligned fetch span" when `--indel-mode adjust` is allowed.
- Either add an explicit aligned-span cap used for tiling/fetch/pairing, or fail fast for indel adjustment unless such a cap is configured.
- Add a regression with an adjusted-in-range fragment whose aligned span exceeds `max_fragment_length`; it should either count correctly or produce the new explicit configuration error.

### L-003 - Medium - The settings sidecar is too thin to interpret the matrix

The command can change the meaning of both length columns and row assignment through `--indel-mode`, `--clip-mode`, `--max-soft-clips`, `--assign-by`, GC weighting, scaling, and grouped BED mode ([config.rs](../../src/commands/lengths/config.rs#L104-L178), [config.rs](../../src/commands/lengths/config.rs#L202-L323)). The sidecar currently writes only `min_fragment_length` and `max_fragment_length` by hand ([lengths.rs](../../src/commands/lengths/lengths.rs#L490-L505)).

There is also a small documentation mismatch: the top-level command doc says the `.npy` shape is `(# windows, # lengths)` ([config.rs](../../src/commands/lengths/config.rs#L18-L20)), but grouped BED mode writes one row per group, not one row per window ([lengths.rs](../../src/commands/lengths/lengths.rs#L408-L425), [lengths.rs](../../src/commands/lengths/lengths.rs#L557-L579)).

Impact: downstream consumers can recover the numeric length range, but not whether those lengths are aligned, indel-adjusted, clip-adjusted, or assigned by overlap/midpoint/proportion. For grouped outputs, the sidecar also does not state that rows correspond to groups.

Recommended fix:

- Replace the manual JSON string with a small serde settings struct.
- Include at least `min_fragment_length`, `max_fragment_length`, `indel_mode`, `clip_mode`, `max_soft_clips`, `assign_by`, window mode, grouped-vs-window row semantics, GC length weighting, and whether GC/scaling inputs were used.
- Clarify the command help so grouped BED output is described as `(# groups, # lengths)`.

### L-004 - Low - BED row reordering silently truncates if metadata and counts drift apart

After BED reduction, the command pairs `bin_info` with `all_bins` using `zip()` and sorts the pairs by original BED index ([lengths.rs](../../src/commands/lengths/lengths.rs#L430-L464)). `zip()` stops at the shorter iterator, so any upstream mismatch between metadata rows and count rows would silently drop the extras before writing `length_counts.npy` and `bins.tsv`.

Impact: current invariants probably keep these lengths equal in normal non-empty runs, but this is a fragile final-output boundary. If an earlier window-filtering or reducer bug appears, the final files can become self-consistent-looking while missing rows.

Recommended fix:

- Add `ensure!(bin_info.len() == all_bins.len(), ...)` immediately before zipping.
- Make `stack_length_counts()` accept a slice and return `Result<Array2<f64>>`, including an explicit empty-input error.
- Add a small unit test for the defensive mismatch check rather than relying only on end-to-end happy paths.

## Existing coverage notes

The command already has broad integration coverage: global, fixed-size, ordinary BED, grouped BED, multi-chromosome runs, tile-boundary invariance, unpaired mode, default MAPQ, GC correction, GC weighting modes, scaling factors, blacklist behavior, indel and clip modes, soft-clip filtering, all window assignment modes, grouped metadata, and reducer helper behavior are represented.

The most important lengths-specific missing tests from this review are deletion-adjusted fragments whose aligned span exceeds the configured max adjusted output length, sidecar completeness for indel/clip/grouped outputs, and defensive metadata/count mismatch handling. The deferred sparse-window GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.
