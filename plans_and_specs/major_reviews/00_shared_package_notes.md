# Shared package review notes

Date: 2026-04-24

Scope: package-level findings discovered while reviewing released commands. Do not duplicate these findings in command-specific review files unless there is a command-specific consequence that needs separate handling.

This file now tracks only active shared findings. Findings that were already implemented have been removed from the active queue. Points that were reviewed and deliberately rejected as findings are listed separately so they are not re-added in later command reviews.

Implemented findings removed from active tracking:

- G-001: shared fragment length ranges can be inverted without a direct error.
- G-003: tiled commands clean temporary directories only on successful completion.
- G-004: `--gc-file` lacks shared fail-fast validation for the required `--ref-2bit`.
- G-005: shared even-fragment midpoint tie-breaking is not reproducible.
- G-007: plain BED modes need a shared no-surviving-windows guard.
- G-010: GC correction packages cannot identify the sample or inputs they were built from.
- G-011: scaling-factor TSV compatibility metadata is minimal but sufficient for first release.
- G-012: short final stride bins are length-weighted only in the numerator.
- G-013: smoothing-weight docs claim the inverted scaling factors have mean 1.0.
- G-014: smoothing-weight TSV writes do not explicitly flush the final writer.
- G-015: all-zero smoothing runs fail with a misleading normalization error.
- G-016: smoothing-weight commands cannot match `fcoverage --ignore-gap` segmentation.

Reviewed and deliberately not tracked:

- G-020: release builds with `plotters` can create QC plots or other auxiliary files by default. This is not a review finding by itself. Commands may legitimately have multiple outputs, and automated pipelines should depend on the specific required input/output artifacts they consume rather than failing because extra files are present. Do not re-add this as an issue unless a command overwrites a primary output, hides a required output contract, or makes a requested primary artifact unavailable.

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

## Released-command re-review additions (2026-05-04)

### Pre-release correctness/safety

#### G-021 - Medium - `--gc-tag` accepts overlong BAM AUX tag names and silently reads the first two bytes

`ApplyGCArgs::validate()` rejects `--gc-file`/`--gc-tag` conflicts and missing `--ref-2bit` for GC files, but it does not validate the supplied `gc_tag` string itself ([cli_common.rs](../../src/commands/cli_common.rs#L705-L713)). Commands then pass the raw string bytes into fragment iteration. In `midpoints`, for example, `gc_tag.as_bytes()` is passed through to `fragments_from_bam()` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L495-L515)), and the shared read mappers pass those bytes directly to `read_gc_tag_from_record()` ([minimal_fragment.rs](../../src/shared/fragment/minimal_fragment.rs#L49-L58), [segment_fragment.rs](../../src/shared/fragment/segment_fragment.rs#L167-L178), [ends_fragment.rs](../../src/shared/fragment/ends_fragment.rs#L114-L124)).

`read_gc_tag_from_record()` calls `Record::aux(tag)` without checking the tag length ([gc_tag.rs](../../src/shared/gc_tag.rs#L171-L193)). The same `rust-htslib` behavior noted in G-017 applies here: BAM AUX tags are two-byte keys, and overlong lookup strings are effectively interpreted by their first two bytes. A user typo such as `--gc-tag GCP` can therefore read `GC` instead of failing fast.

Impact: GC-tag correction can silently use a different tag than the user requested. In the common case users pass `GC`, so this does not affect normal GCParagon/GCfix usage. It matters for typo detection and for pipelines that expose tag names as configuration.

Recommended fix:

- Add shared validation for `--gc-tag`: exactly two ASCII bytes, with a direct error for empty, one-byte, overlong, or non-ASCII values.
- Reuse the same validator anywhere tag names are written or read.
- Add `ApplyGCArgs` tests for invalid tag lengths and at least one command-level regression proving an overlong tag is rejected before reading the BAM.

#### G-022 - Medium - Output prefixes can escape the output directory

Output prefixes are described as sample/file prefixes, but the shared `dot_join()` helper only trims and dot-joins text; it does not reject path separators, `..`, or absolute-looking path components ([io.rs](../../src/shared/io.rs#L14-L27)). Commands then join those names directly under the requested output directory. Confirmed examples include final `ends` output paths ([write.rs](../../src/commands/ends/write.rs#L47-L64), [ends.rs](../../src/commands/ends/ends.rs#L460-L477)), `midpoints` output paths ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L189-L200)), `fcoverage` output paths ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L349-L359)), `lengths` output paths ([lengths.rs](../../src/commands/lengths/lengths.rs#L527-L563), [writer.rs](../../src/commands/lengths/writer.rs#L80-L83)), the `ref-gc-bias` package output path ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L355-L357)), `gc-bias` correction, intermediate, and plot paths ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L742-L746), [gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L1460-L1467), [plotting.rs](../../src/commands/gc_bias/plotting.rs#L140-L180)), and scaling-weight final output/internal source directories used by `fragment-count-weights` and `coverage-weights` ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L135-L178)). The same prefix also feeds temporary directory names in tiled runs through `TempDirGuard` ([ends.rs](../../src/commands/ends/ends.rs#L291-L294), [lengths.rs](../../src/commands/lengths/lengths.rs#L323-L327), [tiled_run.rs](../../src/shared/tiled_run.rs#L642-L656)).

Because the prefix string is interpolated into filenames before `Path::join()`, values containing `/` or `..` are interpreted as path components. A normal sample name is fine, but a mistaken or hostile prefix can write final outputs or temporary directories outside the intended output directory.

Impact: commands do not fully enforce the user-facing invariant that outputs land inside `--output-dir`. This is lower risk than raw chromosome names from input files because the prefix is user-supplied CLI/config text, but automated pipelines commonly pass sample names from metadata, so it should still fail fast.

Recommended fix:

- Add one shared validator for output prefixes and other filename-stem inputs that rejects empty path components, `/`, platform path separators, `..`, and absolute paths.
- Use the validator before any `dot_join()` result is passed to `Path::join()`.
- Add command-level or shared helper regressions showing that prefixes such as `../sample` and `nested/sample` are rejected before output or temporary directories are created.
