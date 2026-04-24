# `cfdna ends` review

Date: 2026-04-24

Scope: `src/commands/ends/*`, the `ends` CLI configuration, and directly used shared helpers for end-aware fragment construction, BED/grouped windows, tile spans, sparse motif output, GC correction, scaling factors, blacklist checks, and midpoint assignment. I also read the existing `ends` tests in `tests/test_ends_command.rs` plus the module-level tests under `src/commands/ends/*_tests.rs`. I did not run tests.

Shared findings that affect this command:

- G-001 in `00_shared_package_notes.md`: shared fragment length ranges can be inverted without a direct error.
- G-002 in `00_shared_package_notes.md`: README examples are not clean runnable snippets.
- G-003 in `00_shared_package_notes.md`: tiled temp directories need cleanup guards.
- G-004 in `00_shared_package_notes.md`: `--gc-file` needs shared fail-fast `--ref-2bit` validation.
- G-005 in `00_shared_package_notes.md`: `--assign-by midpoint` uses shared non-reproducible even-fragment midpoint tie-breaking.
- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.
- G-007 in `00_shared_package_notes.md`: ordinary BED modes need a shared no-surviving-windows guard.

## Findings

### E-001 - High - Covered by shared plain-BED empty-window finding

The ends-specific evidence for ordinary `--by-bed` producing no primary sparse output is now tracked in G-007 in `00_shared_package_notes.md`, together with the equivalent `lengths` failure mode.

### E-002 - High - Raw-shifted GC correction filters by adjusted length but corrects by aligned length

`FragmentWithEnds` keeps two geometries: `interval` is the aligned fragment span and `assignment_interval` is the boundary-adjusted span ([ends_fragment.rs](../../src/shared/fragment/ends_fragment.rs#L37-L44)). `ends` filters fragment lengths using `assignment_len()` ([ends.rs](../../src/commands/ends/ends.rs#L701-L705)), but file-based GC correction passes the aligned `fragment.interval` into `GCCorrector::correct_fragment()` ([ends.rs](../../src/commands/ends/ends.rs#L739-L750)).

The corrector then derives the correction length from the interval it receives ([correct.rs](../../src/commands/gc_bias/correct.rs#L64-L70)) and indexes the GC package by subtracting `length_min` ([correct.rs](../../src/commands/gc_bias/correct.rs#L99-L116)). Package compatibility is validated against the configured min/max fragment lengths ([correct.rs](../../src/commands/gc_bias/correct.rs#L274-L282), [correct.rs](../../src/commands/gc_bias/correct.rs#L332-L340)), which are the same settings that `ends` applies to the adjusted assignment length.

Impact: in `--clip-strategy raw-shifted-boundary`, a soft-clipped fragment can pass length filtering because its adjusted length is inside the configured range, but then GC correction can request a weight for the shorter aligned length. If the package covers only the configured adjusted range, this can error during counting; if the package covers both, the GC weight is still based on a different length/span than the endpoint assignment.

Recommended fix:

- Decide the intended GC geometry for raw-shifted `ends`: adjusted assignment interval, aligned interval, or unsupported with `--gc-file`.
- If aligned GC is intentional, validate that the GC package covers possible aligned lengths after soft-clip adjustment and document the split geometry.
- Add a focused regression with a soft-clipped raw-shifted fragment, an adjusted length different from the aligned length, and a GC package whose length range exposes the mismatch.

### E-003 - High - Motif reference preload is full-tile even when BED windowing narrows the BAM fetch

The GC-prefix part of this issue is now covered by G-006 in `00_shared_package_notes.md`. The ends-specific remaining issue is motif reference preparation: `motif_reference_span_for_tile()` always expands the tile fetch span, not the narrowed fetch span ([motifs.rs](../../src/commands/ends/motifs.rs#L153-L168)), and `build_tile_motif_context()` then reads that whole span and precomputes reference k-mer codes when reference bases are needed ([motifs.rs](../../src/commands/ends/motifs.rs#L239-L261)). The current tests pin this full-tile behavior ([motifs_tests.rs](../../src/commands/ends/motifs_tests.rs#L434-L485)).

Impact: targeted window runs that need reference motifs can spend a large fraction of runtime reading and encoding reference sequence that cannot contribute to output.

Recommended fix:

- Consider narrowing motif reference preload to the reachable assignment/motif span for the tile's relevant windows, while keeping enough padding for outside bases and raw-shifted soft clips.
- Add regression coverage with a sparse BED where no-window tiles do not request motif reference sequence.

### E-004 - Medium - The settings sidecar omits core interpretation fields

`write_end_settings_json()` says it records motif-definition settings and fragment-length filter basis ([write.rs](../../src/commands/ends/write.rs#L70-L73)), but it currently writes only `source_inside`, `clip_strategy`, `window_assignment`, optional base-quality filters, and optional complement collapse ([write.rs](../../src/commands/ends/write.rs#L88-L109)). It does not write the required motif sizes `k_inside` and `k_outside` ([config.rs](../../src/commands/ends/config.rs#L121-L127)), the fragment length range ([config.rs](../../src/commands/ends/config.rs#L225-L228), [cli_common.rs](../../src/commands/cli_common.rs#L118-L129)), `max_soft_clips`, `indel_filter`, or whether dense all-motif output was requested.

Impact: most non-empty motif label files let users infer motif dimensions, but zero-observed sparse outputs have an empty motif file. In those cases the sidecar is the natural place to recover the run's motif and filtering definition, and it is currently incomplete.

Recommended fix:

- Serialize the sidecar through a small serde struct instead of manual string assembly.
- Include at least `k_inside`, `k_outside`, `min_fragment_length`, `max_fragment_length`, `clip_strategy`, `max_soft_clips`, `source_inside`, `indel_filter`, `window_assignment`, `all_motifs`, and configured base-quality filters.
- Add a regression for a sparse run with zero observed motifs that still writes enough settings to interpret the empty result.

### E-005 - Medium - Raw-clipping geometry is split across counting features without a visible contract

The raw-shifted help says the shifted boundary is used for outside-base lookup, window assignment, and motif-level blacklist validation ([config_structs.rs](../../src/commands/ends/config_structs.rs#L315-L321)). In the implementation, window queries use `assignment_interval` ([ends.rs](../../src/commands/ends/ends.rs#L791-L805)), but full-fragment blacklist filtering uses the aligned interval ([ends.rs](../../src/commands/ends/ends.rs#L775-L780)), scaling-bin lookup uses the aligned interval ([ends.rs](../../src/commands/ends/ends.rs#L871-L879)), and non-overlap scaling also computes fragment scaling over the aligned interval ([ends.rs](../../src/commands/ends/ends.rs#L899-L914)).

Impact: users can reasonably expect a raw-shifted boundary choice to apply consistently to fragment-level filters and weights, especially near blacklist or scaling-bin boundaries. The code may be intentionally using aligned genomic span for those fragment-level operations, but that contract is not explicit in the CLI help or settings sidecar.

Recommended fix:

- Document the geometry matrix for `ends`: motif bases, endpoint/window assignment, length filtering, full-fragment blacklist filtering, GC correction, and scaling.
- If aligned geometry is intentional for blacklist/scaling, add tests that pin that behavior at raw-shifted soft-clip boundaries.
- If the adjusted assignment interval is intended instead, update blacklist/scaling calls to use it and add boundary regressions.

## Existing coverage notes

The command already has broad coverage: read-backed and reference-backed motifs, inside-only and outside motifs, dense and sparse outputs, prefix handling, base-quality filters, unpaired mode, hard/soft clipping strategies, indel policy interactions, fragment and motif blacklisting, grouped BED aggregation, all window assignment modes, raw endpoint tile-boundary behavior, GC-file and GC-tag weighting, scaling factors, grouped metadata, settings formatting, and required-reference checks for motif extraction.

The most important ends-specific missing tests from this review are raw-shifted `--gc-file` with adjusted length different from aligned length, motif reference narrowing for skipped sparse-window tiles, sidecar completeness for zero-observed sparse outputs, and explicit raw-clipping geometry contracts for fragment-level blacklist/scaling behavior. Shared missing coverage for ordinary `--by-bed` with zero surviving windows and GC reference pruning is tracked in `00_shared_package_notes.md`.
