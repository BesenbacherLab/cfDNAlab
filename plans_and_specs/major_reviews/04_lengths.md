# `cfdna lengths` review

Date: 2026-04-24

Scope: `src/commands/lengths/*`, the `lengths` CLI configuration, and directly used shared helpers for indel/clip-aware fragment construction, BED/grouped windows, tile spans, GC correction, scaling factors, blacklist checks, midpoint assignment, and length-count reduction. I also read the existing `lengths` tests in `tests/test_lengths_command.rs`. I did not run tests.

Shared findings that affect this command:

- G-002 in `00_shared_package_notes.md`: README OPTIONS blocks need clearer alternative-choice labeling.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release docs/API polish:

- G-002: README OPTIONS blocks should keep their current structure but clarify alternative choices.

Post-release performance:

- G-006: sparse-window GC reference pruning.

## Existing coverage notes

The command already has broad integration coverage: global, fixed-size, ordinary BED, grouped BED, multi-chromosome runs, tile-boundary invariance, unpaired mode, default MAPQ, GC correction, GC weighting modes, scaling factors, blacklist behavior, indel and clip modes, soft-clip filtering, all window assignment modes, grouped metadata, and reducer helper behavior are represented.

The deferred sparse-window GC reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.

## Released-command re-review additions (2026-05-04)

Shared findings that affect this command:

- G-019 in `00_shared_package_notes.md`: tiled temporary files use raw chromosome names as filename components.
- G-022 in `00_shared_package_notes.md`: `--output-prefix` can escape the output directory for final outputs and temporary directories.

Shared findings reviewed and not applied:

- G-021 does not affect `lengths`: this command uses file-only GC correction arguments and does not expose `--gc-tag`.

### Release triage additions

Pre-release correctness/safety:

- L-001: the length reducer can associate temporary tile files with the wrong chromosome when contig names overlap by dotted substring.
- G-019: raw chromosome names in tiled temporary filenames.
- G-022: unchecked output prefixes for final outputs and temp directories.

Post-release performance:

- G-006: sparse-window GC reference pruning.

### Command-specific findings

#### L-001 - Medium - Length reducer matches chromosome temp files by ambiguous substring

`lengths` writes per-tile partial files as `{prefix}.{chr}.{tile_idx}.npz` and crossing-window files as `{prefix}.{chr}.{tile_idx}.npy` ([tiling.rs](../../src/commands/lengths/tiling.rs#L52-L85)). During reduction, it scans the temp directory and keeps files whose names both start with the partial/cross prefix and contain `.{chr}.` ([tiling.rs](../../src/commands/lengths/tiling.rs#L119-L124), [tiling.rs](../../src/commands/lengths/tiling.rs#L152-L157)).

That substring match is ambiguous for legal custom contig names. For example, reducing `chr1` would also match files for `chr1.extra`, because `partials.chr1.extra.0.npz` contains `.chr1.`. Depending on the window counts and crossing metadata, this can either add counts to the wrong chromosome's output rows or fail with a confusing contribution or bounds error. The existing reducer regression proves that `chr1` ignores clearly distinct `chr2` files ([test_lengths_command.rs](../../tests/test_lengths_command.rs#L4830-L4844)), but it does not cover overlapping names such as `chr1` and `chr1.extra`.

Recommended fix:

- Stop discovering chromosome-specific tile files by raw chromosome-name substring. Have tile workers return the exact paths and chromosome/tile identity they wrote, or write safe ordinal-based filenames with a manifest.
- Reuse the same safe per-contig temp-name helper proposed in G-019 so the reducer and writer share one identity model.
- Add a regression with overlapping contig names, for example `chr1` and `chr1.extra`, proving the reducer only reads files for the requested chromosome.
