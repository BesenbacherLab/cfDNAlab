# `cfdna midpoints` review

Date: 2026-04-24

Scope: `src/commands/midpoints/*`, the CLI dispatch for `midpoints`, and directly used shared helpers for BED loading, tiling, overlap lookup, GC correction, scaling, blacklist checks, midpoint placement, and grouped count merging. I also read the existing midpoint-focused tests in `tests/test_profile_groups_command.rs`, `tests/test_profile_groups_counts.rs`, `tests/test_tiling.rs`, `tests/test_cli_smoke.rs`, and cross-command roundtrip/artifact tests. I did not run tests.

Shared findings that affect this command:

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

## Existing coverage notes

The command already has broad integration coverage: length-bin parsing, default MAPQ, paired/unpaired parity, group index ordering, even-length midpoint edge placement, blacklist midpoint behavior, real GC packages, GC tags, GC/scaling multiplication, scaling TSV validation, tile-boundary behavior, chromosome-end fetch narrowing, CLI smoke output, and cross-command BAM/fragment roundtrips are all represented.

The most important midpoint-specific missing tests from this review are sparse final-output roundtrips if M-001B is implemented and the chosen even-midpoint blacklist semantics if that contract changes. The deferred skipped-tile GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.
