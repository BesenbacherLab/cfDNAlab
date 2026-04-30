# `cfdna ref-gc-bias` review

Date: 2026-04-24

Scope: `src/commands/ref_gc_bias/*`, the reference GC package writer/loader boundary, sampling helpers, support-mask helpers, README usage for the GC correction pipeline, and existing `ref-gc-bias` tests in `tests/test_ref_gc_bias.rs` plus module tests in `src/commands/ref_gc_bias/ref_gc_bias_tests.rs`. I did not run tests.

Shared findings that affect this command:

- G-002 in `00_shared_package_notes.md`: README OPTIONS blocks need clearer alternative-choice labeling.
- G-009 in `00_shared_package_notes.md`: `--chromosomes all` is BAM-header only, so reference-only commands cannot use it.

Post-release performance optimizations that affect this command:

- G-006 in `00_shared_package_notes.md`: sparse-window reference sequence reads happen before no-window pruning.

## Release triage

Pre-release correctness/safety:

- RGC-001: reject reference packages with no usable sampled counts.
- RGC-003: add enough reference-package metadata to prevent unsafe reuse.

Pre-release semantic/docs:

- G-002: README OPTIONS blocks should keep their current structure but clarify alternative choices.
- G-009: decide whether reference-only `--chromosomes all` should work for the first release.
- RGC-002: document approximate/tile-size-dependent `--n-positions` behavior or make quotas exact.
- RGC-004: decide whether `--skip-smoothing` should validate unused smoothing parameters.

Post-release performance:

- G-006: sparse-window reference pruning.

## Findings

### RGC-001 - High - The command can write a reference package with no sampled counts

The sampling helpers treat "zero requested positions" and "no valid starts exist" as a zero sampling density rather than an error ([sampling.rs](../../src/shared/sampling.rs#L11-L31)). Per tile, `sample_starts_in_core()` also returns no starts when the selected chromosome is shorter than `max_fragment_length`, when the tile has no eligible start range, or when rounded per-tile sampling chooses zero positions ([sampling.rs](../../src/shared/sampling.rs#L39-L63)).

`ref-gc-bias` only checks that the global sampling density is not above 1.0 ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L143-L152)). It then counts tiles, computes `used_start_positions`, but only logs that value after the package has been written ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L200-L259), [ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L323-L350)). If no selected ACGT positions remain, the support threshold becomes zero and every zero-valued bin is considered supported (`value >= threshold`) ([support_masking.rs](../../src/commands/gc_bias/support_masking.rs#L112-L128)).

Impact: configurations such as `--n-positions 0`, all selected contigs shorter than `--max-fragment-length`, sparse sampling rounded to zero per tile, or a fully masked/empty selected reference can produce a plausible `.ref_gc_package.npz` with no usable empirical information. This finding is about the post-sampling/post-counting guard the command still needs.

Recommended fix:

- Reject `--n-positions 0` at startup.
- Fail if the selected reference has zero valid starts for `max_fragment_length`.
- After reduction, `ensure!(used_start_positions > 0.0)` and `ensure!(total_covered_acgt_positions > 0)` before smoothing, interpolation, support masking, or writing the package.
- Add regressions for zero requested positions, `max_fragment_length` longer than all selected contigs, and all selected sequence masked to non-ACGT.

### RGC-002 - Medium - `--n-positions` is approximate and tile-size-dependent

The CLI describes `--n-positions` as the number of genomic starting positions to sample, uniformly across chromosomes ([config.rs](../../src/commands/ref_gc_bias/config.rs#L56-L71)). The implementation computes one global density, then each tile samples `round(density * possible_in_tile)` positions independently ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L143-L147), [ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L217-L224), [sampling.rs](../../src/shared/sampling.rs#L58-L63)).

Impact: the actual number of sampled starts can differ from `--n-positions`, and changing `--tile-size` can change both the count and the sampled set even with the same seed. For small requested sample sizes across many tiles, per-tile rounding can sample zero starts in every tile. Existing tests cover fixed-seed determinism for the same tile layout and thread-count invariance, but not invariance or documented variance across tile sizes.

Recommended fix:

- Either document `--n-positions` as an approximate density target and print requested vs actual sampled starts before writing, or allocate exact per-chromosome/per-tile quotas.
- Consider reusing the existing Hamilton-apportioned `sample_starts_per_chrom()` helper as the basis for an exact quota plan.
- Add a regression showing the intended behavior when a fixed seed is run with different tile sizes.

### RGC-003 - Medium - Reference package metadata is too thin for safe reuse

The writer stores arrays for counts, support masks, GC-percent widths, schema version, length range, end offset, and smoothing/interpolation settings ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L355-L391)). The downstream loader reconstructs only those fields into `ReferenceGCMetadata` ([load_reference_bias.rs](../../src/commands/gc_bias/load_reference_bias.rs#L16-L25), [load_reference_bias.rs](../../src/commands/gc_bias/load_reference_bias.rs#L69-L140)).

Impact: a package advertised as reusable per assembly cannot identify the reference assembly/path/digest, selected chromosomes, BED inclusion set, blacklist inputs, requested/actual `n_positions`, seed, or tile size. Downstream `gc-bias` can validate shape and length compatibility, but it cannot detect accidental reuse of a package built from the wrong reference, chromosome subset, or inclusion mask.

Recommended fix:

- Add a machine-readable metadata entry, for example a JSON string array, with reference identity, chromosome list, BED/blacklist provenance, requested and actual sampled starts, seed, tile size, and cfDNAlab version.
- Have the loader expose this metadata and let downstream commands log it.
- Keep numeric arrays for fast compatibility checks, but do not rely on filename prefixes for provenance.

### RGC-004 - Low - `--skip-smoothing` still validates unused smoothing parameters

Startup always calls `check_smoothing_settings()` ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L70-L75)). That helper rejects `smoothing_sigma <= 0.0` and `smoothing_sigma > 10.0` even when `skip_smoothing` is true ([config.rs](../../src/commands/ref_gc_bias/config.rs#L166-L177)).

Impact: an ignored smoothing parameter can fail a run. This is low severity because the defaults are valid and the CLI parser handles the common path, but the behavior is still surprising for a flag named `--skip-smoothing`.

Recommended fix:

- Validate smoothing parameters only when smoothing is enabled, or rename/document the check as whole-config validation that applies even to skipped stages.
- Add a small config-level regression for `skip_smoothing = true` with an otherwise invalid sigma.

## Existing coverage notes

The command already has good coverage for written package shapes and scalar metadata, exact distributions, blacklist masking, end offsets, smoothing, interpolation, BED flattening, full-chromosome BED vs global mode, multiple blacklist files, rejection when sampling density exceeds 1.0, and fixed-seed determinism for thread count and identical tile size.

The important missing coverage from this review is zero usable sampled starts, zero covered ACGT positions, tile-size-dependent sampling behavior, `--chromosomes all` on a reference-only command, richer package provenance, and the `skip_smoothing` validation contract. The deferred sparse-window reference pruning optimization is tracked in G-006 in `00_shared_package_notes.md`.
