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
- G-011: scaling-factor TSV compatibility metadata is minimal but sufficient for first release.
- G-012: short final stride bins are length-weighted only in the numerator.
- G-013: smoothing-weight docs claim the inverted scaling factors have mean 1.0.
- G-014: smoothing-weight TSV writes do not explicitly flush the final writer.
- G-015: all-zero smoothing runs fail with a misleading normalization error.
- G-016: smoothing-weight commands cannot match `fcoverage --ignore-gap` segmentation.

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

## Post-release performance optimizations

### G-006 - Sparse-window reference sequence reads happen before no-window pruning

Several tiled commands build GC prefixes from the full tile fetch span before asking the windowing code whether a sparse BED/grouped run should skip or narrow the tile. The pattern appears in `fcoverage` ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L780-L807)), `midpoints` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L368-L415)), `ends` ([ends.rs](../../src/commands/ends/ends.rs#L612-L645)), and `lengths` ([lengths.rs](../../src/commands/lengths/lengths.rs#L633-L668)).

`ref-gc-bias` has the same ordering in its BED mode: `process_tile()` reads the tile sequence and builds GC prefixes before deriving `tile_windows`, and only then returns early for empty tile-window sets ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L416-L475)).

Impact: targeted runs can spend substantial I/O and CPU reading and prefix-building reference sequence for tiles that cannot contribute any output. This directly hurts the use case where sparse windowing should be fastest. This is a performance optimization, not a first-release correctness blocker.

Recommended future fix:

- Move fetch-span adaptation before GC-prefix construction.
- Build GC prefixes over the narrowed fetch span and shift fragment coordinates relative to that narrowed span.
- Add a shared helper-level regression, or one command-level regression per fetch helper family, proving no reference sequence is requested for a no-window tile.

## Converter review additions (2026-05-04)

### Pre-release correctness/safety

#### G-017 - High - Converter AUX tag names are longer than the BAM tag field

`bam-to-bam` documents output tags named `COV`, `CNT`, and `FLEN` ([config.rs](../../src/commands/bam_to_bam/config.rs#L23-L35)), and its writer constants use those same byte strings ([sorted_writer.rs](../../src/commands/bam_to_bam/sorted_writer.rs#L9-L12)) before passing them to `Record::push_aux()` ([sorted_writer.rs](../../src/commands/bam_to_bam/sorted_writer.rs#L201-L218)). `frag-to-bam` uses the same `COV`, `CNT`, and `FLEN` byte strings when converting frag rows back to BAM records ([frag_to_bam.rs](../../src/commands/frag_to_bam/frag_to_bam.rs#L501-L523)).

The package currently depends on `rust-htslib = "0.50.0"` ([Cargo.toml](../../Cargo.toml#L71)). The `rust-htslib` `Record::aux()` contract says only the first two bytes of a tag are used for lookup ([rust-htslib 0.50.0 docs](https://docs.rs/rust-htslib/0.50.0/rust_htslib/bam/record/struct.Record.html#method.aux)). That matches the BAM AUX field shape: the tag key is two bytes. By inspection, `COV`, `CNT`, and `FLEN` therefore do not create three- or four-character BAM tags. They are effectively written and read as `CO`, `CN`, and `FL`.

Impact: downstream tools that inspect normal BAM/SAM optional fields will not see the documented `COV`, `CNT`, or `FLEN` names. The current tests can miss this because they also query with overlong byte strings such as `b"COV"` and `b"CNT"` ([test_bam_to_bam_command.rs](../../tests/test_bam_to_bam_command.rs#L643-L647), [test_bam_to_bam_command.rs](../../tests/test_bam_to_bam_command.rs#L682-L685), [test_bam_to_bam_command.rs](../../tests/test_bam_to_bam_command.rs#L727-L729)), which exercises the same first-two-byte lookup behavior instead of proving the serialized tag names.

Recommended fix:

- Pick explicit two-byte tag names for GC, coverage scaling, count scaling, and fragment length, then update CLI docs, README/docs, converter writers, and downstream consumers to use those names.
- Add tests that inspect the actual serialized AUX tag keys, for example through `aux_iter()` or a SAM text roundtrip, rather than querying overlong names.
- Decide and document what should happen when an input record already has one of the chosen tags. Silent replacement would be risky; a clear fail-fast error or an explicit overwrite flag would be easier to reason about.

#### G-018 - Medium - Converter temporary files use raw chromosome names as path components

`bam-to-frag` resolves chromosome names from the BAM header or user selection ([cli_common.rs](../../src/commands/cli_common.rs#L570-L606), [cli_common.rs](../../src/commands/cli_common.rs#L994-L1002)) and then creates per-chromosome temporary files with `temp_dir.join(format!("{chr}.frag.tsv.zst"))` ([bam_to_frag.rs](../../src/commands/bam_to_frag/bam_to_frag.rs#L256-L267)). `frag-to-bam` accepts chromosome names from `--chrom-sizes` without filename sanitization ([reference.rs](../../src/shared/reference.rs#L149-L187)) and creates temporary files with `temp_dir.join(format!("{}.frag.tmp", frag.chrom))` ([frag_to_bam.rs](../../src/commands/frag_to_bam/frag_to_bam.rs#L289-L295)).

Because those names are interpolated into filesystem paths rather than encoded as filenames, chromosome names containing `/`, `..`, or an absolute-looking path are interpreted as path components. Normal human genome names will not trigger this, but custom references, malformed inputs, or hostile files can cause temp-file creation to fail, create unexpected subdirectories, or write outside the per-run temp directory.

Impact: a converter can be made to write temporary files outside the intended temp directory before the final output is produced. In the common non-hostile case this is more likely to show up as a confusing failure for unusual contig names with slashes, but the safer invariant is that reference names should never control filesystem paths directly.

Recommended fix:

- Name per-chromosome temp files by a trusted ordinal or generated token, for example `chrom.000001.frag.tmp`, and keep the real chromosome name only in an in-memory map.
- Add regression coverage for a chromosome name containing `/` and one containing `..`, proving the command either handles the name safely or rejects it with a direct validation error before writing temp files.
- Prefer one shared helper for per-contig temp filenames so future converter and tiled-command code does not reintroduce raw path components.
