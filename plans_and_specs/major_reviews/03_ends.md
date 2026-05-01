# `cfdna ends` review

Date: 2026-04-24

Scope: `src/commands/ends/*`, the `ends` CLI configuration, and directly used shared helpers for end-aware fragment construction, BED/grouped windows, tile spans, sparse motif output, GC correction, scaling factors, blacklist checks, and midpoint assignment. I also read the existing `ends` tests in `tests/test_ends_command.rs` plus the module-level tests under `src/commands/ends/*_tests.rs`. I did not run tests.

Shared findings that affect this command:

- G-002 in `00_shared_package_notes.md`: README OPTIONS blocks need clearer alternative-choice labeling.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release correctness/safety:

- E-005: raw-clipping geometry needs a visible contract for blacklist/scaling behavior.

Pre-release docs/API polish:

- G-002: README OPTIONS blocks should keep their current structure but clarify alternative choices.

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

### E-005 - Medium - Raw-clipping geometry is split across counting features without a visible contract

The raw-shifted help says the shifted boundary is used for outside-base lookup, window assignment, and motif-level blacklist validation ([config_structs.rs](../../src/commands/ends/config_structs.rs#L315-L321)). In the implementation, window queries use `assignment_interval` ([ends.rs](../../src/commands/ends/ends.rs#L791-L805)), but full-fragment blacklist filtering uses the aligned interval ([ends.rs](../../src/commands/ends/ends.rs#L775-L780)), scaling-bin lookup uses the aligned interval ([ends.rs](../../src/commands/ends/ends.rs#L871-L879)), and non-overlap scaling also computes fragment scaling over the aligned interval ([ends.rs](../../src/commands/ends/ends.rs#L899-L914)).

Impact: users can reasonably expect a raw-shifted boundary choice to apply consistently to fragment-level filters and weights, especially near blacklist or scaling-bin boundaries. The code may be intentionally using aligned genomic span for those fragment-level operations, but that contract is not explicit in the CLI help or settings sidecar.

Recommended fix:

- Document the geometry matrix for `ends`: motif bases, endpoint/window assignment, length filtering, full-fragment blacklist filtering, GC correction, and scaling.
- If aligned geometry is intentional for blacklist/scaling, add tests that pin that behavior at raw-shifted soft-clip boundaries.
- If the adjusted assignment interval is intended instead, update blacklist/scaling calls to use it and add boundary regressions.

## Existing coverage notes

The command already has broad coverage: read-backed and reference-backed motifs, inside-only and outside motifs, dense and sparse outputs, prefix handling, base-quality filters, unpaired mode, hard/soft clipping strategies, indel policy interactions, fragment and motif blacklisting, grouped BED aggregation, all window assignment modes, raw endpoint tile-boundary behavior, GC-file and GC-tag weighting, scaling factors, grouped metadata, settings formatting, and required-reference checks for motif extraction.

The most important ends-specific missing tests from this review are explicit raw-clipping geometry contracts for fragment-level blacklist/scaling behavior. Motif reference narrowing for skipped sparse-window tiles and shared GC reference pruning are deferred performance optimizations tracked here and in G-006 in `00_shared_package_notes.md`.
