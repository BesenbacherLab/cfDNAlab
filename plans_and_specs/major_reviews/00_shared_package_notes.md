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
- G-008: feature-gated QC plots were previously listed as resolved; the active remaining issue is now tracked as G-020.
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

## Released-command re-review additions (2026-05-04)

### Pre-release correctness/safety

#### G-019 - Medium - Tiled command temporary files use raw chromosome names as path components

Several released tiled commands interpolate BAM/reference chromosome names directly into per-tile temporary filenames. Confirmed sites include `fcoverage` ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L457-L482)), `midpoints` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L234-L240)), `ends` ([ends.rs](../../src/commands/ends/ends.rs#L607-L612)), and `lengths` ([tiling.rs](../../src/commands/lengths/tiling.rs#L52-L85)). `prepare-windows` already shows the safer direction by sanitizing `/` before creating chromosome temp files ([writers.rs](../../src/commands/prepare_windows/writers.rs#L105-L108)), though a more complete shared helper should handle `..` and other path-component edge cases as well.

Because chromosome names are not filesystem-safe by definition, names containing `/`, `..`, or other path separators can turn one intended filename into nested path components. Normal references will not trigger this, but custom references or hostile inputs can cause confusing write failures or writes outside the intended per-run temp directory.

Impact: tiled commands can fail or write intermediate files in unexpected locations before final output is produced. This is the same class of issue as G-018, but it affects the tiled command path rather than the converter-specific temp files.

Recommended fix:

- Introduce one shared per-contig temp-name helper that maps chromosome names to safe ordinals or escaped tokens.
- Use that helper in `fcoverage`, `midpoints`, `ends`, `lengths`, and converter temp-file creation so the same invariant is enforced everywhere.
- Add regression coverage with chromosome names containing `/` and `..`, proving the commands either encode names safely or reject them before creating temp files.

#### G-021 - Medium - `--gc-tag` accepts overlong BAM AUX tag names and silently reads the first two bytes

`ApplyGCArgs::validate()` rejects `--gc-file`/`--gc-tag` conflicts and missing `--ref-2bit` for GC files, but it does not validate the supplied `gc_tag` string itself ([cli_common.rs](../../src/commands/cli_common.rs#L705-L713)). Commands then pass the raw string bytes into fragment iteration. In `midpoints`, for example, `gc_tag.as_bytes()` is passed through to `fragments_from_bam()` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L495-L515)), and the shared read mappers pass those bytes directly to `read_gc_tag_from_record()` ([minimal_fragment.rs](../../src/shared/fragment/minimal_fragment.rs#L49-L58), [segment_fragment.rs](../../src/shared/fragment/segment_fragment.rs#L167-L178), [ends_fragment.rs](../../src/shared/fragment/ends_fragment.rs#L114-L124)).

`read_gc_tag_from_record()` calls `Record::aux(tag)` without checking the tag length ([gc_tag.rs](../../src/shared/gc_tag.rs#L171-L193)). The same `rust-htslib` behavior noted in G-017 applies here: BAM AUX tags are two-byte keys, and overlong lookup strings are effectively interpreted by their first two bytes. A user typo such as `--gc-tag GCP` can therefore read `GC` instead of failing fast.

Impact: GC-tag correction can silently use a different tag than the user requested. In the common case users pass `GC`, so this does not affect normal GCParagon/GCfix usage. It matters for typo detection and for pipelines that expose tag names as configuration.

Recommended fix:

- Add shared validation for `--gc-tag`: exactly two ASCII bytes, with a direct error for empty, one-byte, overlong, or non-ASCII values.
- Reuse the same validator anywhere tag names are written or read.
- Add `ApplyGCArgs` tests for invalid tag lengths and at least one command-level regression proving an overlong tag is rejected before reading the BAM.

### Pre-release docs/API polish

#### G-020 - Low - Release builds with `plotters` still create QC plots by default

The README install/build commands explicitly enable `plotters` ([README.md](../../README.md#L40-L48)). In that build, several released commands still write plot files as a default side effect. `midpoints` defaults `plot_groups` to `[0]` in both CLI and programmatic config ([config.rs](../../src/commands/midpoints/config.rs#L201-L216), [config.rs](../../src/commands/midpoints/config.rs#L240-L247)) and invokes plotting after writing the main outputs ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L301-L318)). `lengths` always writes an overall fragment-length PNG when `plotters` is enabled and bins are available ([lengths.rs](../../src/commands/lengths/lengths.rs#L568-L627)). `gc-bias` likewise plots correction summaries under the same feature gate ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L748-L765)).

Impact: users following the documented installation path get extra PNG outputs even when they only requested machine-readable artifacts. This is not a scientific correctness issue, but it is still a default side effect and can surprise automated pipelines.

Recommended fix:

- Make plot output opt-in, or add a shared `--no-plots`/`--plots` policy and document it consistently.
- If plots remain enabled by default in `plotters` builds, list the PNG outputs explicitly in command docs and README examples.
- Add smoke coverage for the chosen default so this does not drift silently again.
