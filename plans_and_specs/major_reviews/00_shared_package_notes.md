# Shared package review notes

Date: 2026-04-24

Scope: package-level findings discovered while reviewing released commands. Do not duplicate these findings in command-specific review files unless there is a command-specific consequence that needs separate handling.

This file now tracks only active shared findings. Findings that were already implemented have been removed from the active queue.

Implemented findings removed from active tracking:

- G-001: shared fragment length ranges can be inverted without a direct error.
- G-003: tiled commands clean temporary directories only on successful completion.
- G-004: `--gc-file` lacks shared fail-fast validation for the required `--ref-2bit`.
- G-005: shared even-fragment midpoint tie-breaking is not reproducible.
- G-007: plain BED modes need a shared no-surviving-windows guard.
- G-008: feature-gated QC plots are default command side effects.
- G-010: GC correction packages cannot identify the sample or inputs they were built from.
- G-013: smoothing-weight docs claim the inverted scaling factors have mean 1.0.
- G-014: smoothing-weight TSV writes do not explicitly flush the final writer.
- G-015: all-zero smoothing runs fail with a misleading normalization error.
- G-016: smoothing-weight commands cannot match `fcoverage --ignore-gap` segmentation.

## Pre-release correctness and release-safety

### G-011 - Medium - Scaling-factor TSV metadata is too thin for safe reuse

The smoothing-weight writer emits only one metadata comment, `# gc_mode=...`, before the scaling TSV header ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L149-L163)). The downstream loader models only that GC mode in `ScalingFactorsMetadata` ([scale_genome.rs](../../src/shared/scale_genome.rs#L47-L55)), ignores every other metadata key ([scale_genome.rs](../../src/shared/scale_genome.rs#L502-L530)), and downstream command setup checks only GC-mode compatibility before using the map ([cli_common.rs](../../src/commands/cli_common.rs#L987-L1004), [scale_genome.rs](../../src/shared/scale_genome.rs#L546-L587)).

Impact: downstream commands cannot detect a scaling TSV made from the wrong BAM, wrong smoothing source (`coverage-weights` vs `fragment-count-weights`), chromosome set, bin size, stride, fragment-length range, MAPQ threshold, pairing mode, blacklist set, reference, GC source/package, or cfDNAlab version. The README tells users to build the factors from the same BAM, but the artifact itself cannot prove or even summarize that provenance.

Recommended fix:

- Add machine-readable metadata comments, preferably one JSON line, with at least source command, BAM identity/fingerprint, contigs/chromosome selection, `bin_size`, `stride`, fragment-length range, MAPQ/proper-pair/unpaired settings, blacklist provenance, GC mode and GC source identity, reference identity when relevant, and cfDNAlab version.
- Have `load_scaling_factors_tsv()` parse and expose the metadata rather than ignoring all non-`gc_mode` keys.
- Decide which mismatches should be hard errors versus warnings. Source type, BAM identity when available, contig identity, and GC/source incompatibilities are the highest-value checks.

## Pre-release docs/API polish

### G-002 - Low - README OPTIONS blocks need clearer alternative-choice labeling

The README examples are intentionally written as command skeletons with an `OPTIONS` section, not as many separate runnable snippets. That structure is acceptable and does not need a radical rewrite. The remaining issue is narrower: some OPTIONS groups list alternatives that a new user might read as cumulative unless the alternative-choice contract is stated explicitly.

The `fcoverage` OPTIONS block lists mutually exclusive window selectors together: `--by-size`, `--by-bed`, and `--by-grouped-bed` ([README.md](../../README.md#L275-L294)).

The `ends` and `lengths` OPTIONS blocks also list mutually exclusive window selectors together ([README.md](../../README.md#L316-L345), [README.md](../../README.md#L361-L385)). That is fine as long as the text makes clear that these are alternatives, not one command to run unchanged.

The `midpoints` block shows two alternative `--length-bins` forms inside the OPTIONS section ([README.md](../../README.md#L397-L417)). The intent is clear to a careful reader, but a short note would make the contract explicit.

This is documentation polish, not a code-correctness issue and not a reason to replace the OPTIONS format.

Recommended fix:

- Keep the current command-plus-OPTIONS structure.
- Add one explicit sentence before the OPTIONS blocks: choose one option from each alternative group; do not run every OPTIONS line together.
- Where alternatives are shown, label them as alternatives in prose or comments rather than implying they are cumulative arguments.

### G-009 - Low - `--chromosomes all` is BAM-header only, so reference-only commands cannot use it

The shared chromosome resolver treats `--chromosomes all` as "read contigs from the BAM header" and errors when no BAM path is available ([cli_common.rs](../../src/commands/cli_common.rs#L504-L525)). `ref-gc-bias` is reference-only and calls that resolver with `None` before opening the `.2bit` file ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L70-L75)), even though the reference helper can read `.2bit` contig names and sizes ([reference.rs](../../src/shared/reference.rs#L38-L51)).

Impact: the default remains `chr1` through `chr22`, but users cannot ask `ref-gc-bias` to use all contigs in a reference with the same `--chromosomes all` spelling used by BAM-backed commands. They must know and provide an explicit chromosome list or file, which is easy to miss for non-human references or assemblies with non-autosomal contigs.

Recommended fix:

- Split chromosome resolution by source type, or add a resolver that accepts a contig-provider callback.
- For reference-only commands, make `all` resolve from the reference file rather than requiring a BAM.
- Add one regression for `ref-gc-bias --chromosomes all` on a tiny multi-contig `.2bit`.

## Post-release performance optimizations

### G-006 - Sparse-window reference sequence reads happen before no-window pruning

Several tiled commands build GC prefixes from the full tile fetch span before asking the windowing code whether a sparse BED/grouped run should skip or narrow the tile. The pattern appears in `fcoverage` ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L780-L807)), `midpoints` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L368-L415)), `ends` ([ends.rs](../../src/commands/ends/ends.rs#L612-L645)), and `lengths` ([lengths.rs](../../src/commands/lengths/lengths.rs#L633-L668)).

`ref-gc-bias` has the same ordering in its BED mode: `process_tile()` reads the tile sequence and builds GC prefixes before deriving `tile_windows`, and only then returns early for empty tile-window sets ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L416-L475)).

Impact: targeted runs can spend substantial I/O and CPU reading and prefix-building reference sequence for tiles that cannot contribute any output. This directly hurts the use case where sparse windowing should be fastest. This is a performance optimization, not a first-release correctness blocker.

Recommended future fix:

- Move fetch-span adaptation before GC-prefix construction.
- Build GC prefixes over the narrowed fetch span and shift fragment coordinates relative to that narrowed span.
- Add a shared helper-level regression, or one command-level regression per fetch helper family, proving no reference sequence is requested for a no-window tile.
