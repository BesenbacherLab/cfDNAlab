# Shared package review notes

Date: 2026-04-24

Scope: package-level findings discovered while reviewing released commands. Do not duplicate these findings in command-specific review files unless there is a command-specific consequence that needs separate handling.

## Findings

### G-001 - High - Shared fragment length ranges can be inverted without a direct error [IMPLEMENTED]

`FragmentLengthArgs` validates each bound independently as `>= 10`, but there is no shared validation that `min_fragment_length <= max_fragment_length` ([cli_common.rs](../../src/commands/cli_common.rs#L118-L129)). The shared `contains()` helper then returns false for every fragment when the range is inverted ([cli_common.rs](../../src/commands/cli_common.rs#L139-L142)).

This affects multiple released commands that flatten `FragmentLengthArgs`, including `fcoverage`, `ends`, `lengths`, conversion commands, and the smoothing-weight commands. `lengths` has a sharper failure mode: `LengthCounts::new()` computes `length_max - length_min + 1` without a guard ([counting.rs](../../src/commands/lengths/counting.rs#L25-L27)), and the command builds that template from the CLI range before counting ([lengths.rs](../../src/commands/lengths/lengths.rs#L253-L258)).

`ref-gc-bias` has a different sharp failure mode: it calls `gc_percent_widths()` from the raw CLI range before constructing `GCCounts` ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L137-L142)), and that helper uses an `assert!` rather than returning a user-facing error for inverted ranges ([counting.rs](../../src/commands/gc_bias/counting.rs#L765-L771)).

Impact: a user typo like `--min-fragment-length 500 --max-fragment-length 100` can silently produce empty or misleading outputs instead of a direct configuration error. In `lengths`, it can also panic in debug builds or wrap into an impractical allocation in optimized builds.

Recommended fix:

- Add a `FragmentLengthArgs::validate()` method that errors when `min_fragment_length > max_fragment_length`.
- Call it at the start of every command that accepts `FragmentLengthArgs`.
- Add one shared unit test for the helper and at least one command-level regression test for a released command.

### G-002 - Medium - README examples mix runnable commands with mutually exclusive option variants

The `fcoverage` README block is written as one shell command, but it includes all three mutually exclusive window selectors in the same continued command: `--by-size`, `--by-bed`, and `--by-grouped-bed` ([README.md](../../README.md#L275-L294)). It also places comments after line-continuation backslashes, which makes copy-paste behavior shell-sensitive ([README.md](../../README.md#L276-L281)).

The `midpoints` block has the same copy-paste problem: comments follow continued lines, alternative `--length-bins` forms appear inside one apparent command, and one continued command line is followed by a blank/comment block before more arguments ([README.md](../../README.md#L397-L417)).

The `ends` block repeats the pattern: comments follow continued lines, optional arguments appear after a blank line that terminates the continued shell command, and all three mutually exclusive window selectors are shown in one apparent command ([README.md](../../README.md#L316-L345)).

The `lengths` block repeats it too: comments follow continued lines, a blank/comment block splits the apparent command, and all three mutually exclusive window selectors are shown together ([README.md](../../README.md#L361-L385)).

The `ref-gc-bias` block has the same copy-paste risk: the `--output-prefix` line puts a comment after a continuation backslash, so the example is not a clean shell command as shown ([README.md](../../README.md#L136-L142)).

This is not limited to one implementation module. It is a documentation pattern risk for a CLI with many mutually exclusive or alternative options.

Impact: users may copy examples that cannot run as written, or they may miss which options are alternatives.

Recommended fix:

- Split examples into separate runnable snippets.
- Use prose before each snippet instead of inline comments after `\`.
- Keep one short "option variants" list outside code blocks if needed.

### G-003 - High - Tiled commands clean temporary directories only on successful completion

Several released tiled commands create a temp directory inside the output directory and rely on an explicit `remove_dir_all()` on the normal path instead of a cleanup guard. `midpoints` creates its temp directory before tiled counting ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L155-L157)) and removes it only after merge, final writes, and optional plotting have all succeeded ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L302-L311)). The same unguarded cleanup pattern exists in `fcoverage` ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L295-L296), [fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L745-L752)), `ends` ([ends.rs](../../src/commands/ends/ends.rs#L435-L443)), `lengths` ([lengths.rs](../../src/commands/lengths/lengths.rs#L467-L475)), and `gc-bias` ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L327-L329), [gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L399-L408)). `lengths` removes its temp directory before final matrix/metadata writes, so failures before that cleanup leak temp files while failures after cleanup leave partial final outputs without the temp files available for inspection.

Impact: any error after temp directory creation and before explicit cleanup can leave `tmp.<prefix>.<random>` directories behind. This includes validation failures discovered inside tile workers, interrupted writes, merge errors, and plotting errors in commands that plot before cleanup.

Recommended fix:

- Introduce a small temp-dir guard that removes the directory on drop unless explicitly preserved.
- Keep the current warning behavior for cleanup failures, but make cleanup run on error paths too.
- Add one helper-level test for the guard and one command-level regression that forces a post-temp error.

### G-004 - High - `--gc-file` lacks shared fail-fast validation for the required `--ref-2bit` [IMPLEMENTED]

`ApplyGCArgs::validate()` checks only that `--gc-file` and `--gc-tag` are not combined ([cli_common.rs](../../src/commands/cli_common.rs#L629-L637)). Commands with file-based GC correction then perform their own `ref_2bit` checks later. In `midpoints`, startup calls `opt.gc.validate()` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L75-L81)), but the missing-reference error is raised inside `process_tile()` after output setup and temp directory creation ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L368-L372)). `ends` has the same delayed check: startup only requires `--ref-2bit` for motif extraction cases ([ends.rs](../../src/commands/ends/ends.rs#L132-L140)), while the GC-file missing-reference error is raised inside tile processing after the temp directory is created ([ends.rs](../../src/commands/ends/ends.rs#L281-L284), [ends.rs](../../src/commands/ends/ends.rs#L612-L617)). `lengths` loads the GC package at startup ([lengths.rs](../../src/commands/lengths/lengths.rs#L204-L213)), but it does not reject a missing `--ref-2bit` until tile processing after temp directory creation ([lengths.rs](../../src/commands/lengths/lengths.rs#L260-L305), [lengths.rs](../../src/commands/lengths/lengths.rs#L633-L638)).

`fcoverage` also validates only the mutually exclusive GC source pair at startup ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L134-L140)); the missing-reference error is raised later inside tile processing ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L780-L784)). The smoothing-weight commands inherit that delayed failure through their internal `fcoverage` call: they validate `ApplyGCArgs`, create an internal output directory, build an `fcoverage` config, and only then call `fcoverage_run_inner()` ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L102-L123), [coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L247-L252)). Those commands wrap the internal output directory in a cleanup guard, but the configuration error is still discovered after output-side effects.

Impact: a simple configuration error does unnecessary setup/work and can combine with G-003 to leak temp directories. The help text already says `--ref-2bit` is required when `--gc-file` is used, so this should fail before output-side effects.

Recommended fix:

- Add a shared helper for optional-reference GC validation, for example `validate_gc_reference(&ApplyGCArgs, Option<&Path>)`.
- Call it at the start of every command that accepts `ApplyGCArgs` plus an optional `ref_2bit`.
- Keep deeper checks only as defensive assertions/context, not the first user-visible validation.

### G-005 - Medium - Shared even-fragment midpoint tie-breaking is not reproducible

The shared midpoint helper uses thread-local randomness for even-length fragments ([midpoint.rs](../../src/shared/midpoint.rs#L9-L35)). Released commands that call it inherit run-to-run variation for exact boundary placement, including `midpoints` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L508-L510)), `ends` ([ends.rs](../../src/commands/ends/ends.rs#L791-L799)), `lengths` ([lengths.rs](../../src/commands/lengths/lengths.rs#L890-L897)), and `gc-bias` ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L104-L112)).

Impact: repeated runs over identical inputs can produce different exact counts at the two central bases of even-length fragments. The stochastic effect may average out at scale, but it complicates bit-for-bit reproducibility and can matter in small targeted fixtures or edge windows.

Recommended fix:

- Replace thread-local random tie-breaking with a deterministic hash of stable fragment identity and coordinates, or expose an explicit seed.
- Document the chosen reproducibility contract in the shared midpoint helper and command help.
- Keep tests that assert both centers are reachable, but add deterministic-output tests for command-level runs.

### G-006 - High - Sparse-window reference sequence reads happen before no-window pruning

Several tiled commands build GC prefixes from the full tile fetch span before asking the windowing code whether a sparse BED/grouped run should skip or narrow the tile. The pattern appears in `fcoverage` ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L780-L807)), `midpoints` ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L368-L415)), `ends` ([ends.rs](../../src/commands/ends/ends.rs#L612-L645)), and `lengths` ([lengths.rs](../../src/commands/lengths/lengths.rs#L633-L668)).

`ref-gc-bias` has the same ordering in its BED mode: `process_tile()` reads the tile sequence and builds GC prefixes before deriving `tile_windows`, and only then returns early for empty tile-window sets ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L416-L475)).

Impact: targeted runs can spend substantial I/O and CPU reading and prefix-building reference sequence for tiles that cannot contribute any output. This directly hurts the use case where sparse windowing should be fastest.

Recommended fix:

- Move fetch-span adaptation before GC-prefix construction.
- Build GC prefixes over the narrowed fetch span and shift fragment coordinates relative to that narrowed span.
- Add a shared helper-level regression, or one command-level regression per fetch helper family, proving no reference sequence is requested for a no-window tile.

### G-007 - High - Plain BED modes need a shared no-surviving-windows guard

The shared BED offset helper returns `total = 0` when no BED windows survive the selected chromosome filter ([windowing.rs](../../src/shared/windowing.rs#L99-L111)). Grouped BED paths in reviewed commands often guard this case, but ordinary BED paths do not have a shared validation point. `ends` loads ordinary BED windows without checking for at least one surviving window ([ends.rs](../../src/commands/ends/ends.rs#L161-L170)); `lengths` has the same ordinary BED load path ([lengths.rs](../../src/commands/lengths/lengths.rs#L150-L159)). `ref-gc-bias` also loads, filters, and flattens ordinary BED windows without checking that any selected positions remain ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L97-L118)).

Impact: a wrong chromosome selection, empty BED, or wrong-assembly BED can fail late and unclearly. In `ends`, the default sparse writer can return without writing the primary sparse output for an empty bin set. In `lengths`, `all_bins` can remain empty and then `stack_length_counts()` indexes `all_counts[0]` ([lengths.rs](../../src/commands/lengths/lengths.rs#L378-L405), [counting.rs](../../src/commands/lengths/counting.rs#L235-L238)). In `ref-gc-bias`, this can proceed into a zero-count reference package instead of telling the user that no reference positions were selected.

Recommended fix:

- Add a shared helper that validates ordinary BED windows after chromosome filtering and returns a command-neutral error like "BED file did not contain any valid windows on the selected chromosomes".
- Call it immediately after ordinary BED loading, mirroring the grouped BED checks.
- Keep output writers and stackers defensive by returning clear errors on empty primary output collections.

### G-008 - Medium - Feature-gated QC plots are default command side effects

When the `plotters` feature is enabled, some commands do plotting by default as part of counting. `midpoints` defaults `plot_groups` to `[0]` ([config.rs](../../src/commands/midpoints/config.rs#L197-L212)) and always calls plotting under the feature ([midpoints.rs](../../src/commands/midpoints/midpoints.rs#L283-L300)). `lengths` has no plot option and always writes an overall PNG when bins exist ([lengths.rs](../../src/commands/lengths/lengths.rs#L507-L553)); that plot happens after `length_counts.npy` and settings are written but before `bins.tsv` or `group_index.tsv` metadata ([lengths.rs](../../src/commands/lengths/lengths.rs#L481-L585)). `gc-bias` writes the correction package first, then always writes four QC PNGs under the feature, with no plot opt-out ([gc_bias.rs](../../src/commands/gc_bias/gc_bias.rs#L706-L738)).

Impact: production counting can spend time on PNG generation and can fail after primary outputs have already been written. In `lengths`, that can leave the primary matrix without the metadata file that tells users what rows mean. In `gc-bias`, a plot failure can make the command return an error after the reusable correction package already exists.

Recommended fix:

- Make plotting opt-in, or add an explicit `--no-plots`/empty CLI value that is easy to use.
- Write required machine-readable metadata before optional QC plots.
- Consider treating plot failures as warnings when primary scientific outputs have already succeeded.

### G-009 - Medium - `--chromosomes all` is BAM-header only, so reference-only commands cannot use it

The shared chromosome resolver treats `--chromosomes all` as "read contigs from the BAM header" and errors when no BAM path is available ([cli_common.rs](../../src/commands/cli_common.rs#L504-L525)). `ref-gc-bias` is reference-only and calls that resolver with `None` before opening the `.2bit` file ([ref_gc_bias.rs](../../src/commands/ref_gc_bias/ref_gc_bias.rs#L70-L75)), even though the reference helper can read `.2bit` contig names and sizes ([reference.rs](../../src/shared/reference.rs#L38-L51)).

Impact: the default remains `chr1` through `chr22`, but users cannot ask `ref-gc-bias` to use all contigs in a reference with the same `--chromosomes all` spelling used by BAM-backed commands. They must know and provide an explicit chromosome list or file, which is easy to miss for non-human references or assemblies with non-autosomal contigs.

Recommended fix:

- Split chromosome resolution by source type, or add a resolver that accepts a contig-provider callback.
- For reference-only commands, make `all` resolve from the reference file rather than requiring a BAM.
- Add one regression for `ref-gc-bias --chromosomes all` on a tiny multi-contig `.2bit`.

### G-010 - High - GC correction packages cannot identify the sample or inputs they were built from

The correction package stores only schema version, end offset, length/GC edges, the correction matrix, and length-bin frequencies ([package.rs](../../src/commands/gc_bias/package.rs#L11-L19), [package.rs](../../src/commands/gc_bias/package.rs#L46-L57)). The shared CLI help says `--gc-file` should come from the same BAM file ([cli_common.rs](../../src/commands/cli_common.rs#L581-L586)), but downstream loading validates only package shape, schema version, length-bin coverage, and end-offset compatibility ([package.rs](../../src/commands/gc_bias/package.rs#L75-L118), [correct.rs](../../src/commands/gc_bias/correct.rs#L274-L342)).

Impact: released commands that accept `--gc-file` cannot detect a correction package made from the wrong BAM, reference 2bit, reference GC package, blacklist set, window settings, pairing mode, MAPQ threshold, or cfDNAlab version. A swapped or stale package can produce plausible weighted outputs with no warning beyond the user remembering the file provenance.

Recommended fix:

- Add a machine-readable metadata entry to the correction package, for example JSON containing BAM identity/fingerprint, reference identity, reference GC package metadata, blacklist provenance, chromosome/window settings, filtering settings, binning/outlier settings, and cfDNAlab version.
- Have downstream loaders expose and log this metadata.
- Decide which mismatches should be hard errors versus warnings; at minimum, detect known-incompatible settings such as reference identity, BAM identity when available, and end-offset/length-range conflicts.

### G-011 - Medium - Scaling-factor TSV metadata is too thin for safe reuse

The smoothing-weight writer emits only one metadata comment, `# gc_mode=...`, before the scaling TSV header ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L149-L163)). The downstream loader models only that GC mode in `ScalingFactorsMetadata` ([scale_genome.rs](../../src/shared/scale_genome.rs#L47-L55)), ignores every other metadata key ([scale_genome.rs](../../src/shared/scale_genome.rs#L502-L530)), and downstream command setup checks only GC-mode compatibility before using the map ([cli_common.rs](../../src/commands/cli_common.rs#L987-L1004), [scale_genome.rs](../../src/shared/scale_genome.rs#L546-L587)).

Impact: downstream commands cannot detect a scaling TSV made from the wrong BAM, wrong smoothing source (`coverage-weights` vs `fragment-count-weights`), chromosome set, bin size, stride, fragment-length range, MAPQ threshold, pairing mode, blacklist set, reference, GC source/package, or cfDNAlab version. The README tells users to build the factors from the same BAM, but the artifact itself cannot prove or even summarize that provenance.

Recommended fix:

- Add machine-readable metadata comments, preferably one JSON line, with at least source command, BAM identity/fingerprint, contigs/chromosome selection, `bin_size`, `stride`, fragment-length range, MAPQ/proper-pair/unpaired settings, blacklist provenance, GC mode and GC source identity, reference identity when relevant, and cfDNAlab version.
- Have `load_scaling_factors_tsv()` parse and expose the metadata rather than ignoring all non-`gc_mode` keys.
- Decide which mismatches should be hard errors versus warnings. Source type, BAM identity when available, contig identity, and GC/source incompatibilities are the highest-value checks.

### G-012 - Medium - Short final stride bins are length-weighted only in the numerator

`fill_triangular_overlap()` accounts for a short final stride bin by multiplying its contribution by `len_ratio = bin_slice[j].size() / stride` in the weighted sum, but it still adds the full integer kernel weight to the denominator ([striding.rs](../../src/commands/coverage_weights/striding.rs#L156-L168)). The command docs describe chromosome-edge kernels as truncated and normalized by the remaining weights, not as padded with missing sequence that behaves like zero support ([config.rs](../../src/commands/coverage_weights/config.rs#L79-L82)).

Impact: when a contig length is not an exact multiple of `--stride`, smoothed support near the right chromosome edge is depressed by the missing part of the final short stride bin. That inflates the reciprocal scaling factors near contig ends in both `coverage-weights` and `fragment-count-weights`. A constant-support genome should remain constant after edge smoothing; with the current denominator, the tail can drift below the interior solely because the final bin is shorter.

Recommended fix:

- Decide the intended contract: length-normalized truncation at contig edges, or zero-padding outside the contig.
- If truncation is intended, accumulate the denominator using the same length factor as the numerator, for example `sum_w += weight * len_ratio`, or compute the numerator and denominator in base-pair units.
- Add a helper-level regression with a short final stride bin where every bin has identical average support; the smoothed values should stay identical if edge truncation is the contract.

### G-013 - Medium - Smoothing-weight docs claim the inverted scaling factors have mean 1.0

Both smoothing-weight command docs say the factors are inverted and then claim non-zero factors have `mean == 1.0` ([fragment_count_weights/config.rs](../../src/commands/fragment_count_weights/config.rs#L19-L20), [coverage_weights/config.rs](../../src/commands/coverage_weights/config.rs#L19-L20)). The implementation first computes a global mean of `average_overlap_coverage`, then writes each non-zero scaling factor as `1.0 / (value / mean)` when `invert` is true ([striding.rs](../../src/commands/coverage_weights/striding.rs#L196-L252)).

Impact: the pre-inversion normalized support has mean 1.0, but the reciprocal scaling factors generally do not. Users doing QC on the `scaling_factor` column can be told by the help text to expect a property the output does not mathematically preserve.

Recommended fix:

- Update both command docs to say that non-zero smoothed support is normalized to mean 1.0 before inversion, and `scaling_factor` is the reciprocal multiplier.
- If a factor-column mean of 1.0 is actually desired, the normalization algorithm needs a second factor-space normalization step and downstream expectations should be reviewed.

### G-014 - Medium - Smoothing-weight TSV writes do not explicitly flush the final writer

The shared smoothing-weight writer creates a `BufWriter<File>` for the final scaling TSV, writes metadata, header, and rows, then logs the output path and returns without calling `flush()` ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L144-L185)). Any write error that surfaces only during the final buffered flush would be lost because the writer is dropped after the function has already decided to return `Ok(())`.

Impact: `coverage-weights` and `fragment-count-weights` can report success after a late filesystem error, especially on full disks, network filesystems, or interrupted writes where the final buffered bytes are still pending.

Recommended fix:

- Call `tsv_writer.flush().context("flushing scaling-factors TSV")?` before logging success.
- Consider `File::sync_all()` only if you want a stronger durability guarantee than ordinary command-line tools usually provide.
- Add a small writer-helper regression if the writer logic is extracted; otherwise, keep this as a targeted code-review fix.

### G-015 - Low - All-zero smoothing runs fail with a misleading normalization error

`normalize_average_overlap_by_global_mean()` skips non-finite and effectively zero overlap values when computing the global mean ([striding.rs](../../src/commands/coverage_weights/striding.rs#L204-L223)). If every stride bin is zero after filtering, blacklisting, or GC correction, it returns `no bins to normalize or all had length 0` ([striding.rs](../../src/commands/coverage_weights/striding.rs#L225-L227)). The caller invokes this after internal `fcoverage` has already completed and the stride bins have been smoothed ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L126-L137)).

Impact: a common user problem such as the wrong chromosome selection, an overly strict MAPQ/length filter, a fully blacklisted selection, or all GC-corrected fragments receiving zero weight produces an error that points at missing bins or zero-length bins rather than "no non-zero fragment support after filtering".

Recommended fix:

- Track total bins separately from bins with finite positive support and return a specific error when all bins are zero.
- Include the selected chromosomes and key filters in the error context if practical.
- Add one command-level regression with a tiny BAM and a filter combination that removes all support.

### G-016 - Medium - Smoothing-weight commands cannot match `fcoverage --ignore-gap` segmentation [IMPLEMENTED]

`fcoverage` exposes `--ignore-gap` to exclude the inter-mate gap from paired-end fragment coverage ([config.rs](../../src/commands/fcoverage/config.rs#L242-L249)), and the tile worker passes `!opt.ignore_gap` into fragment segmentation ([fcoverage.rs](../../src/commands/fcoverage/fcoverage.rs#L840-L845)). The smoothing-weight shared args do not expose an equivalent setting ([scaling_weights_config.rs](../../src/commands/coverage_weights/scaling_weights_config.rs#L12-L109)), and the internal `fcoverage` config builder sets filtering, blacklist, GC, and stride windows, but never sets `ignore_gap`, leaving the `fcoverage` default of full-span counting in place ([coverage_weights.rs](../../src/commands/coverage_weights/coverage_weights.rs#L233-L252), [config.rs](../../src/commands/fcoverage/config.rs#L312-L320)).

Impact: users can run downstream `fcoverage --ignore-gap` while the smoothing factors were still built from full-span coverage. That makes the genomic smoothing profile include inter-mate gap mass that the downstream coverage command intentionally omits. The same option-parity problem applies to fragment-count smoothing when users want count support over read-covered segments only.

Recommended fix:

- Add `--ignore-gap` to the shared smoothing-weight args and pass it through to the internal `fcoverage` config.
- Reuse the existing `fcoverage` validation that rejects `--ignore-gap` with `--reads-are-fragments`, but run it before creating output-side effects.
- Add a regression with a paired fragment containing an inter-mate gap that proves smoothing weights differ when `--ignore-gap` is enabled.
