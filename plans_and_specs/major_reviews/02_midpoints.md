# `cfdna midpoints` review

Date: 2026-04-24

Scope: `src/commands/midpoints/*`, the CLI dispatch for `midpoints`, and directly used shared helpers for BED loading, tiling, overlap lookup, GC correction, scaling, blacklist checks, midpoint placement, and grouped count merging. I also read the existing midpoint-focused tests in `tests/test_profile_groups_command.rs`, `tests/test_profile_groups_counts.rs`, `tests/test_tiling.rs`, `tests/test_cli_smoke.rs`, and cross-command roundtrip/artifact tests. I did not run tests.

Shared findings that affect this command:

- G-002 in `00_shared_package_notes.md`: README OPTIONS blocks need clearer alternative-choice labeling.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release semantic/docs:

- G-002: README OPTIONS blocks should keep their current structure but clarify alternative choices.
- M-005: midpoint blacklist semantics should be either aligned with counted midpoint placement or documented explicitly.

Post-release performance/scalability:

- G-006: sparse-window GC reference pruning.
- M-001: dense per-tile profile storage can be too large for sparse targeted runs.

## Findings

### M-001 - High - Per-tile profile storage is dense across all groups, positions, and length bins

Every active tile allocates a full `ProfileGroupsCounts` for the complete output shape ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L432-L434)). The flattened allocation is `num_groups * window_size * num_length_bins` `f32`s ([counting_by_group.rs](../../src/commands/midpoints/counting_by_group.rs#L63-L69)), and each tile writes that full dense vector to a temp `.npy` even when only a few groups or positions were touched ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L660-L663)). The final merge then allocates the same full shape again and reads every dense tile file back ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L263-L266)).

Impact: targeted midpoint profiles can become memory- and disk-bound long before BAM iteration is the bottleneck. With many groups, a 2001 bp profile, and range-style length bins, a single tile can require hundreds of MB or more even if that tile contains sparse sites. The help text warns only that memory increases with the number of bins ([config.rs](../../src/commands/midpoints/config.rs#L84-L98)), but the real scaling also includes groups, window size, thread count, and number of active tile temp files.

Recommended fix:

- Consider a sparse per-tile accumulator keyed by `(group_idx, length_bin_idx, position)` and reduce sparse triplets into the final dense NPY.
- If dense output remains the final format, add an upfront memory/disk estimate and fail fast above a practical threshold.
- At minimum, document the full shape cost in `--length-bins` or command help.

### M-005 - Medium - `--blacklist-strategy midpoint` uses a different even-fragment midpoint than counting

Blacklist filtering happens before the command samples the counted midpoint ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L495-L510)). The shared blacklist midpoint strategy uses `start + (end - start) / 2`, which is the right-center base for even-length half-open fragments ([overlaps.rs](../../src/shared/blacklist/overlaps.rs#L74-L96)). The count placement then randomly selects left or right center for the same even-length fragment ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L508-L510)).

The existing test suite pins this current behavior: an even fragment can be deterministically blacklisted by the right-center base even though the profile placement would otherwise be randomized between two central bases ([test_profile_groups_command.rs](../../tests/test_profile_groups_command.rs#L586-L662)).

Impact: users may read `midpoint` blacklist strategy as "filter by the midpoint that would be counted", but that is not true for even-length fragments at blacklist edges.

Recommended fix:

- Decide whether `midpoints` should reuse the sampled counted midpoint for midpoint-blacklist filtering.
- If the deterministic blacklist midpoint is preferred, state that explicitly in the `--blacklist-strategy midpoint` help for commands that also randomize midpoint placement.

## Existing coverage notes

The command already has broad integration coverage: length-bin parsing, default MAPQ, paired/unpaired parity, group index ordering, even-length midpoint edge placement, blacklist midpoint behavior, real GC packages, GC tags, GC/scaling multiplication, scaling TSV validation, tile-boundary behavior, chromosome-end fetch narrowing, CLI smoke output, and cross-command BAM/fragment roundtrips are all represented.

The most important midpoint-specific missing tests from this review are memory-shape guardrails for large group/bin configurations and the chosen even-midpoint blacklist semantics if that contract changes. The deferred skipped-tile GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.
