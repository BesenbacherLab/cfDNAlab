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

## Existing coverage notes

The command already has broad integration coverage: global, fixed-size, ordinary BED, grouped BED, multi-chromosome runs, tile-boundary invariance, unpaired mode, default MAPQ, GC correction, GC weighting modes, scaling factors, blacklist behavior, indel and clip modes, soft-clip filtering, all window assignment modes, grouped metadata, and reducer helper behavior are represented.

The most important lengths-specific missing test from this review covers deletion-adjusted fragments whose aligned span exceeds the configured max adjusted output length. The deferred sparse-window GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.
