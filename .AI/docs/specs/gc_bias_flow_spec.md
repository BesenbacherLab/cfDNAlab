# GC Bias Flow Spec

This spec describes the current `ref-gc-bias` -> `gc-bias` -> downstream `--gc-file` flow. It is based on the implementation in `src/commands/ref_gc_bias`, `src/commands/gc_bias`, shared fragment construction, and shared GC-correction loading.

## Source Anchors

- `src/commands/ref_gc_bias/config.rs`
- `src/commands/ref_gc_bias/ref_gc_bias.rs`
- `src/commands/gc_bias/config.rs`
- `src/commands/gc_bias/gc_bias.rs`
- `src/commands/gc_bias/counting.rs`
- `src/commands/gc_bias/package.rs`
- `src/commands/gc_bias/load_reference_bias.rs`
- `src/commands/gc_bias/correct.rs`
- `src/commands/gc_bias/binning.rs`
- `src/commands/gc_bias/support_masking.rs`
- `src/commands/gc_bias/interpolation.rs`
- `src/commands/gc_bias/windows.rs`
- `src/shared/fragment/minimal_fragment.rs`
- `src/shared/read.rs`
- `src/shared/constants.rs`

## Overall Contract

- `cfdna ref-gc-bias` builds an assembly-level expected GC table. Its output is reusable across samples aligned to the same selected reference contigs, chromosome set, fragment length range, end offset, and blacklist/windowing decisions.
- `cfdna gc-bias` builds a sample-specific correction package by comparing a BAM-derived observed GC table to the reference table.
- Downstream commands load the sample-specific correction package with `--gc-file` and multiply fragment contributions by the looked-up correction weight.
- The shared schema version for both reference and sample GC packages is `GC_CORRECTION_SCHEMA_VERSION = 2`.
- The shared minimum ACGT support for GC-fraction calculation is 10 bases.

## Fragment Geometry

- Paired-end fragments span `[forward.pos, reverse.reference_end)`.
- The forward and reverse reads must be on the same `tid`, have opposite strand flags, and be inward-oriented with `forward.pos <= reverse.pos`.
- Unpaired read-as-fragment mode spans `[read.pos, read.reference_end)`.
- `ref-gc-bias` does not read fragments from BAM. It simulates fragments by sampled start position and candidate fragment length.
- `gc-bias` does not expose its own fragment length range. It inherits the inclusive fragment length range from the loaded reference GC package.
- Fragment GC is counted after trimming `end_offset` bases from each fragment end. Effective GC length is `fragment_length - 2 * end_offset`.
- `ref-gc-bias` rejects configurations where the minimum effective GC length is below 10 bp.

## Reference GC Command

`cfdna ref-gc-bias` samples reference start positions and counts, for every configured fragment length, how many sampled reference fragments fall into each length by GC-percent cell.

### Inputs And Validation

- `--ref-2bit` is required.
- `--output-dir` is required.
- Fragment length defaults come from shared `FragmentLengthArgs`: 30..1000 bp inclusive.
- `--end-offset` defaults to 10 and is validated so the minimum effective GC length is at least 10 bp.
- `--n-positions` is a positive approximate sampling target.
- `--seed` makes tile-level sampling deterministic for a fixed thread-independent tile set.
- `--smoothing-sigma` must be finite, positive, and `<= 10.0` when smoothing is enabled.
- `--smoothing-radius` is parsed in the CLI as 1..9.
- `--tile-size` is parsed in the CLI as at least 1,000,000 bp.
- Selected chromosomes are resolved from the 2bit reference and optional chromosome arguments before counting.

### Sampling

- Sampling density is `n_positions / total_possible_starts`.
- `total_possible_starts` sums `chrom_len - max_fragment_length + 1` for each selected chromosome long enough to fit the maximum fragment length.
- The command fails if no selected chromosome is long enough for `max_fragment_length`.
- The command fails if the sampling density exceeds 1.0.
- Each tile samples from starts in its core, capped so the maximum fragment length stays within the chromosome.
- Per-tile sampling draws `round(density * possible_in_tile)` unique starts without replacement, sorts them, and returns all valid starts only when the rounded estimate reaches the number of possible starts.

### Windows

- With no BED, global mode is used.
- With `--by-bed`, BED intervals are loaded, selected by chromosome, and flattened by merging overlapping or touching intervals into unique positions.
- BED windowing limits reference counting to sampled starts inside the merged intervals, and each simulated fragment must fit inside the tile-local window used for counting.
- There is no fixed-size reference-window mode in `ref-gc-bias`.

### Blacklists

- Optional blacklist BED files are loaded into a chromosome map.
- For each tile, blacklisted positions are masked to `N` in the loaded reference sequence before GC prefixes are built.
- A simulated reference fragment is counted only when its end-trimmed GC interval has full ACGT support. Because the reference counter is called with `min_acgt_fraction = 1.0`, any `N` or blacklisted base in the trimmed interval excludes that fragment-length/start combination.
- Reference ACGT support used for support-mask scaling is counted only in the tile core, clipped to the active reference windows, to avoid double-counting tile halos.

### Counting Shape

- Raw `GCCounts` stores ragged rows by absolute fragment length.
- Within a row, raw columns are absolute GC counts from `0..=effective_length`.
- After optional smoothing, rows are collapsed to integer GC percent bins `0..=100`.
- GC percent rounding uses integer half-up rounding:

```text
gc_percent = min(100, (100 * gc_count + acgt_count / 2) / acgt_count)
```

- `gc_percent_widths` records, for each length and GC percent, how many raw GC counts round to that percent.
- After conversion to GC percent bins, both `ref-gc-bias` and `gc-bias` divide each populated percent bin by its width and rescale the row to preserve the row sum. This corrects uneven integer percent-bin widths.

### Reference Processing

- Raw GC-count rows are smoothed before conversion to GC percent bins unless `--skip-smoothing` is set.
- Smoothing is one-dimensional per fragment length row over raw GC-count bins, using a Gaussian kernel with the configured sigma and radius and reflected boundary indexing.
- After width correction, `support_mask_outliers` is built before interpolation.
- The reference support threshold is:

```text
threshold_per_mb = 1 + n_positions / 100000000
threshold = total_covered_acgt_positions / 1,000,000 * threshold_per_mb
```

- A cell in `support_mask_outliers` is true when its reference count is at least that threshold.
- If interpolation is enabled, unsupported cells in each length row are filled with a degree-2 local polynomial using up to 3 neighboring anchors per side and at least 3 neighbors. The support mask is not updated by this reference interpolation.
- `support_mask_unobservables` is theoretical. It is true for length by GC-percent cells that can be produced by at least one integer GC count for that effective length and false for theoretically impossible cells.

### Reference Package Output

The output path is `<prefix>.ref_gc_package.zarr`.

The Zarr store contains:

- Root attributes include `cfdnalab_schema = "reference_gc_package"`,
  `cfdnalab_schema_version`, `end_offset`, `skip_interpolation`, `smoothing_radius`,
  `smoothing_sigma`, `skip_smoothing`, and the GC-percent rounding/minimum-ACGT settings.
- `counts`: `float64[length, gc_percent]` matrix after smoothing, GC-percent width correction, and optional interpolation.
- `support_mask_unobservables`: `bool[length, gc_percent]` theoretical observable-cell mask.
- `support_mask_outliers`: `bool[length, gc_percent]` low-reference-support mask from before interpolation.
- `gc_percent_widths`: `uint16[length, gc_percent]` width table for GC percent bins.
- `length`: `int32[length]` concrete fragment lengths.
- `gc_percent`: `int32[gc_percent]` concrete integer GC-percent values.
- `chromosome`: `int32[chromosome]` selected chromosome indices with JSON chromosome labels in array attributes.
- `reference_contig_footprint_json`: `uint8[json_byte]` 2bit contig footprint encoded as JSON bytes.

## Sample GC Command

`cfdna gc-bias` loads the reference package, counts observed cfDNA fragment GC by length and GC percent from a BAM, normalizes the observed table against the expected reference table, and writes multiplicative correction weights.

### Inputs And Validation

- `--bam`, `--output-dir`, `--ref-2bit`, and `--ref-gc-file` are required through the command input structs.
- The reference GC package is loaded before sample counting.
- The selected chromosome set must equal the reference package chromosome set. The check is set-based, not order-based.
- The current `--ref-2bit` contig footprint must exactly match the footprint stored in the reference package.
- Fixed-size `gc-bias --by-size` windows must be at least the reference package maximum fragment length.
- `--require-proper-pair` cannot be combined with `--reads-are-fragments`.
- The default sample GC windowing is fixed 100,000 bp windows unless `--global`, `--by-size`, or `--by-bed` is set.

### Read And Fragment Filters

- Paired-end read inclusion removes unmapped reads, unmapped mates, cross-tid mates, same-orientation pairs, secondary reads, supplementary reads, duplicates, QC-failed reads, reads below `--min-mapq`, and non-proper pairs when `--require-proper-pair` is set.
- Unpaired read-as-fragment inclusion removes unmapped, secondary, supplementary, duplicate, QC-failed, and low-MAPQ reads.
- Fragment construction then requires inward-oriented paired fragments.
- Fragment length is filtered to the reference package length range.
- Tile ownership in `gc-bias` is by fragment start in the tile core. Fragments whose start is outside the tile core are not counted in that tile.

### Sample GC Counting

- For each tile, the command loads the reference sequence for the tile fetch span, masks blacklist intervals to `N`, and builds GC/ACGT prefixes.
- Sample `gc-bias` counts a fragment only when the contracted GC interval is inside the fetched sequence.
- The contracted GC interval must contain at least 10 ACGT bases.
- The contracted GC interval must be 100% ACGT for `gc-bias` estimation because `get_fragment_gc` is called with `min_acgt_fraction = 1.0`.
- The counted GC value is the absolute GC count in the end-trimmed interval. It is converted to GC percent later using the same half-up rounding as the reference command.

### Window Assignment

- `gc-bias` supports global, fixed-size, and BED windows.
- Fixed-size windows are contiguous from chromosome coordinate 0, with the last window clipped to chromosome length.
- BED windows are loaded as plain windows and are not flattened by `gc-bias`.
- `--assign-by` controls how fragments map to windows:
  - `count-overlap`: contributes the fraction of the assignment interval that overlaps each selected window.
  - `any`: contributes 1.0 to windows with any overlap.
  - `all`: contributes 1.0 only when nearly the full assignment interval overlaps.
  - `midpoint`: assignment interval is the deterministic midpoint base.
  - `proportion=<threshold>`: contributes 1.0 when the overlap fraction meets the threshold.
- For `any` and `count-overlap`, the effective minimum overlap fraction is `1 / (max_fragment_length + 1)`.
- For `all` and `midpoint`, the effective threshold is `1 - 1 / (max_fragment_length + 1)`.
- For all non-midpoint assignment modes, the assignment interval is the full fragment span.

### Window Scaling

- Global mode does not apply per-window scaling. Tile-level counts merge directly.
- BED and fixed-size modes scale each accepted window before averaging windows together.
- A window is ignored if its ACGT percentage is below `--min-window-acgt-pct`, if it has zero fragment counts, or if its counted ACGT support is zero.
- For accepted windows, counts are scaled by:

```text
(1 / mean_count_across_gc_count_cells) * (usable_acgt_bases_in_window / average_window_span)
```

- `average_window_span` is the mean BED span, the whole selected genome span in global mode, or the mean fixed-size chromosome-bin span including clipped chromosome-end bins.
- For windows crossing tile boundaries, partial counts are spilled to temporary `.npz` sidecars and merged by stable window index before window scaling.

### Sample Processing

- Scaled window counts are averaged by dividing the merged scaled sum by the number of accepted windows. In global mode, the effective scaling weight is 1 when at least one tile contributed.
- Sample raw GC-count rows are smoothed with the reference package smoothing settings unless the reference package has `skip_smoothing = true`.
- Sample counts are converted to length by GC-percent bins and corrected by `gc_percent_widths` from the reference package.
- The command verifies that cells marked false by `support_mask_unobservables` have zero sample counts.
- The full-resolution sample GC-percent matrix can be saved as an intermediate named `gc_bias.avg_cfdna_counts.0.npy` when `--save-intermediates` is set.
- The sample matrix is globally mean-scaled using only cells true in `support_mask_outliers` to estimate the mean.
- If reference interpolation was enabled, sample cells false in `support_mask_outliers` are interpolated per length row with the same degree-2 local polynomial settings used by the reference flow. In this sample interpolation, the mutable mask is updated for interpolated cells.
- A disabled code path exists for extra 2D Gaussian smoothing of normalized sample counts. In current behavior `do_smoothing` is hard-coded to `false`.

### Binning And Correction Matrix

- Length bins are built greedily from the normalized sample matrix along the length axis using `--min-length-bin-mass` and `--min-length-bin-width`.
- GC bins are built greedily from the normalized sample matrix along the GC-percent axis using `--min-gc-bin-mass` and width 1.
- Greedy binning accumulates contiguous indices until both minimum mass and minimum width are met. A trailing partial bin is merged into the previous bin when possible.
- Reference and sample matrices are first collapsed across GC bins with sum aggregation.
- They are then collapsed across length bins with unweighted mean aggregation.
- Extreme correction support is false for the configured number of GC bins at both GC tails and for the configured number of shortest length bins.
- Binned sample and reference counts are mean-scaled per length row while excluding those extreme-support false cells from the mean.
- Extreme-support false cells are set to 1.0 in both normalized binned matrices before division.
- Raw bias is computed as:

```text
normalized_binned_sample_counts / normalized_binned_reference_counts
```

- The raw bias matrix is mean-scaled per length row while excluding extreme-support false cells from the mean.
- Extreme-support false correction cells are interpolated first across rows and then across columns using the shared unsupported-bin polynomial interpolation.
- Outlier handling is then applied unless `--outlier-method none` is selected.
- The default outlier method is global Tukey IQR with `k = 3.0`.
- `quantile`, `iqr`, `stddev`, and `mad` modes estimate bounds from supported cells. Unsupported cells are also winsorized using those supported-cell bounds.
- After outlier handling, all correction values are hard-clamped to `[0.1, 10.0]`.
- The matrix is mean-scaled per length row again without a mask.
- The matrix is inverted elementwise, keeping exact zeros as zero. The final stored values are multiplicative correction weights.
- Length-bin frequencies are computed from binned sample counts as each length-bin total divided by the total binned count.

### Sample Package Output

The output path is `<prefix>.gc_bias_correction.zarr`.

The Zarr store contains:

- Root attributes include `cfdnalab_schema = "gc_correction_package"`,
  `cfdnalab_schema_version`, `end_offset`, and the GC-percent rounding/minimum-ACGT settings.
- `correction_matrix`: `float64[length_bin, gc_bin]` multiplicative weights.
- `length_edges`: `uint32[length_edge]` inclusive/exclusive bin edges with the final edge treated as inclusive on readback.
- `gc_edges`: `uint32[gc_edge]` inclusive/exclusive GC-percent bin edges with the final edge treated as inclusive on readback.
- `length_bin_frequencies`: `float64[length_bin]` normalized sample length-bin frequencies.
- `reference_contig_footprint_json`: `uint8[json_byte]` 2bit contig footprint inherited from the reference package, encoded as JSON bytes.

With the `plotters` feature, `gc-bias` also writes:

- `avg_gc_bias_across_lengths_unweighted.png`
- `avg_gc_bias_across_lengths_weighted.png`
- `gc_bias_by_selected_lengths_80_220bp.png` when at least one selected length from 80, 100, ..., 220 bp is covered by the correction package.
- `gc_bias_heatmap.png`
- `gc_bias_heatmap.bins.png`

With `--save-intermediates`, staged `.npy` files use names of the form `<prefix>.gc_bias.<tag>.<index>.npy`.

## Downstream `--gc-file` Application

- Downstream commands load `gc_bias_correction.zarr` through `GCCorrectionPackage::from_file`.
- The public package path must exist and have a `.zarr` extension.
- The package schema version must match `GC_CORRECTION_SCHEMA_VERSION`.
- `--gc-file` requires `--ref-2bit` at CLI/config validation time.
- When a downstream command supplies `--ref-2bit`, the current reference contig footprint must exactly match the footprint stored in the correction package.
- The downstream requested fragment length range must lie inside the correction package length range.
- The downstream minimum fragment length must exceed twice the package end offset.
- Downstream file-based correction reads reference sequence for the fetch span needed by that command and builds GC prefixes without applying command blacklist masks to those prefixes.
- `GCCorrector::correct_fragment` contracts the fragment interval by the package end offset, computes integer GC percent with at least 10 ACGT bases and `min_acgt_fraction = 0.0`, maps fragment length and GC percent into package bins, and returns the stored multiplicative weight.
- A fragment outside the package length range returns no correction weight.
- Missing, invalid, out-of-range, non-finite, or otherwise unusable GC weights skip fragments by default in downstream counting commands.
- `--neutralize-invalid-gc` keeps those fragments with weight 1.0 instead.
- `--gc-file` and `--gc-tag` are mutually exclusive in commands that expose both.
- `--gc-tag` uses a two-character SAM/BAM AUX tag. For paired fragments, the combined fragment tag value comes from the mate tags and downstream commands use the classified tag weight instead of file-based GC lookup.

## Practical Invariants

- Use the same reference assembly for `ref-gc-bias`, `gc-bias`, and downstream `--gc-file` consumers. The code validates this by contig footprint.
- Use the same selected chromosome set for `ref-gc-bias` and `gc-bias`. The code validates set equality.
- Keep blacklist choices aligned between `ref-gc-bias` and `gc-bias`. The code does not prove the blacklist inputs are the same, but both commands mask blacklisted reference bases before estimating bias.
- The correction package is sample-specific. CLI help and shared args state it should be made from the same BAM file that downstream commands correct.
- Correction weights are multipliers. A fragment's existing contribution is multiplied by the GC weight.
