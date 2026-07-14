# cfDNAlab | R Loaders <img src="https://raw.githubusercontent.com/BesenbacherLab/cfDNAlab/refs/heads/main/cfdnalab_logo_little_guy_172x200_144dpi.png" align="right" height="155" />

R helpers for loading [cfDNAlab](https://github.com/BesenbacherLab/cfDNAlab) analysis outputs.

This package does not install or run the cfDNAlab command-line tool. The CLI is distributed separately as the Rust `cfdna` binary. Use this R package after running cfDNAlab to load, inspect, and reshape output files in R.

This package supports midpoint, end-motif, and reference k-mer Zarr outputs plus length-count TSV outputs: `<prefix>.midpoint_profiles.zarr`, `<prefix>.end_motifs.zarr`, `<prefix>.ref_kmers.zarr`, and `<prefix>.length_counts.tsv.zst`.

The helpers return base `data.frame` objects, R arrays, and `Matrix` sparse matrices. Convert data frames with `tibble::as_tibble()` or `data.table::as.data.table()` when you want those workflows.

Numeric indices returned by this R package are one-based, matching ordinary R indexing.

NOTE: While the main CLI tool is highly tested and validated, this R package is currently being built and may have bugs or use too AI'ish language in the documentation. The core functions should work and we are actively improving it over the coming weeks. We decided to share it early to help you use the outputs of the main tool.

<br>

## Install

Install the current development version from GitHub:

```r
install.packages("pak")
pak::pak("cfdnalab=github::BesenbacherLab/cfDNAlab/r-cfdnalab")
```

<br>

## Midpoint Profiles

Midpoint profile stores contain a 3D count array with axes:

```text
group x length_bin x position
```

Use the metadata helpers to inspect the axes, then extract one profile for a group and length bin.

```r
library(cfdnalab)

midpoints <- read_midpoints("sample.midpoint_profiles.zarr")

group_metadata(midpoints)
length_bins(midpoints)
positions(midpoints)

profile <- midpoint_data_frame(
  midpoints,
  groups = "LYL1",
  with_lengths = 167
)

head(profile)
```

Use `with_length_range = c(start, end)` to select all whole length bins that overlap a half-open bp range.

Use `profile_array()` when you only need the count vector for one profile. Use `midpoint_array()` only when you want the full 3D count array in memory.

<br>

## End Motifs

End-motif stores can be dense or sparse. Check the storage mode before choosing which count helper to use.

```r
ends <- read_end_motifs("sample.end_motifs.zarr")

storage_mode(ends)
row_mode(ends)
motifs(ends)
has_motif(ends, "_AA")
```

`storage_mode(ends)` tells you how counts are stored on disk. Either as sparse or dense arrays.

`row_mode(ends)` tells you what each row in the count table represents:

- `"global"`: one row for the whole input file. Use the global count helpers.
- `"size"`: rows are fixed-size genomic windows from `--window-size`. Use `window_metadata(ends)` and the window count helpers.
- `"bed"`: rows are BED intervals. Use `window_metadata(ends)` and the window count helpers.
- `"grouped_bed"`: rows are BED groups. Use `group_metadata(ends)`, `group_idx()`, and the group count helpers.

`window_metadata(ends)` returns `window_idx`, `chrom`, `start`, `end`, and `blacklisted_fraction`.

Sparse output keeps only non-zero counts in memory:

```r
counts <- sparse_counts_matrix(ends)
motif_counts <- end_motif_data_frame(ends, motifs = "_AA")
selected_motifs <- end_motif_data_frame(ends, motifs = c("_AA", "_CC"))
selected_motif_idxs <- end_motif_data_frame(ends, motif_idxs = c(1L, 4L))
```

Dense output can be read as a matrix or data frame:

```r
counts <- dense_counts_matrix(ends)
motif_counts <- end_motif_data_frame(ends, motifs = "_AA")
```

Dense helpers do not silently convert sparse stores. If you want a dense matrix from sparse output, pass `allow_densify = TRUE`.

```r
counts <- dense_counts_matrix(ends, allow_densify = TRUE)
```

Windowed output supports `window_idxs` and blacklist filtering:

```r
windows <- window_metadata(ends)
window_counts <- end_motif_data_frame(
  ends,
  window_idxs = c(1L, 3L, 4L),
  motifs = c("_AA", "_CC"),
  densify = TRUE,
  max_blacklisted_fraction = 0.1
)
```

Grouped output supports group names and group indices:

```r
groups <- group_metadata(ends)
group_idx(ends, "t-cells")
group_counts <- sparse_counts_matrix(
  ends,
  groups = c("t-cells", "b-cells"),
  motifs = c("_AA", "_CC")
)
group_idx_counts <- end_motif_data_frame(ends, group_idxs = c(1L, 3L), motif_idxs = 2L)
```

<br>

## Reference K-Mers

Reference k-mer outputs contain row-wise frequencies. Each row is the whole reference, a genomic window, or a BED group depending on how `cfdna ref-kmers` was run. Counts are reconstructed by multiplying each row by its `row_scaling_factor`.

```r
ref_kmers <- read_ref_kmers("sample.ref_kmers.zarr")

storage_mode(ref_kmers)
row_mode(ref_kmers)
motifs(ref_kmers)
row_scaling_factors(ref_kmers)

frequencies <- sparse_frequencies_matrix(ref_kmers)
counts <- sparse_counts_matrix(ref_kmers)
rows <- ref_kmer_data_frame(ref_kmers)
```

Sparse output stores only non-zero values. Dense helpers do not silently create a zero-filled matrix from sparse output. If you want a dense matrix from sparse output, pass `allow_densify = TRUE`. For observed-only output, zero filling covers the motif axis returned by `motifs(ref_kmers)`: the combined set of motifs or motifs-file targets observed anywhere in the output. It does not add every possible k-mer unless `all_motifs(ref_kmers)` is `TRUE`.

```r
dense_frequencies_matrix(ref_kmers, allow_densify = TRUE)
dense_counts_matrix(ref_kmers, allow_densify = TRUE)
```

Use `motifs` for k-mer or k-mer-group labels and `motif_idxs` for one-based
motif indices:

```r
selected_kmers <- ref_kmer_data_frame(ref_kmers, motifs = c("ACGT", "TGCA"))
selected_kmer_counts <- dense_counts_matrix(
  ref_kmers,
  motifs = c("ACGT", "TGCA"),
  allow_densify = TRUE
)
```

Windowed output supports `window_metadata()` and `window_idxs`:

```r
windows <- window_metadata(ref_kmers)
window_rows <- ref_kmer_data_frame(
  ref_kmers,
  window_idxs = c(1L, 6L, 7L),
  motifs = c("ACGT", "TGCA"),
  densify = TRUE,
  max_blacklisted_fraction = 0.1
)
```

Grouped output supports `group_metadata()`, `group_idx()`, `groups`, and
`group_idxs`:

```r
groups <- group_metadata(ref_kmers)
group_idx(ref_kmers, "promoters")
group_counts <- sparse_counts_matrix(
  ref_kmers,
  groups = c("promoters", "enhancers"),
  motifs = c("ACGT", "TGCA")
)
group_rows <- ref_kmer_data_frame(ref_kmers, groups = "promoters", densify = TRUE)
```

`ref_kmer_data_frame()` returns both `frequency` and reconstructed `count`.

<br>

## Correct End-Motif Counts

Pass a matching reference k-mer output to correct end-motif counts for reference sequence composition. The correction is normalized so that a uniform reference composition leaves counts unchanged. Motifs that are common in the reference are scaled down, while rare motifs are scaled up. Corrected data frames include both `corrected_count` and `corrected_frequency`.

When both `--k-inside` and `--k-outside` was used, creating motif labels such as `"AC_GT"`, specify how the two sides are used in the correction via `two_sided_correction`:

- `"joint"` keeps the full `"AC_GT"` label and corrects its count using the frequency of the exact reference k-mer `"ACGT"`. Use this when the full pairing of outside and inside bases is the quantity of interest.

- `"split"` also keeps the full `"AC_GT"` label, but calculates separate correction factors for outside label `"AC"` and inside label `"GT"`, then multiplies them. Use this when full two-sided sample motifs should remain separate, but outside and inside reference composition should be modeled independently or exact full reference k-mers are too sparse.

- `"outside"` combines sample counts that share the same outside bases before correction. For example, `"AC_AA"` and `"AC_GT"` both contribute to `"AC_"`. The result contains outside labels such as `"AC_"` and uses the summed reference frequency of full k-mers beginning with `"AC"`.

- `"inside"` combines sample counts that share the same inside bases before correction. For example, `"AA_GT"` and `"AC_GT"` both contribute to `"_GT"`. The result contains inside labels such as `"_GT"` and uses the summed reference frequency of full k-mers ending with `"GT"`.

For `"split"`, `"outside"`, and `"inside"`, these side frequencies are calculated from the loaded full-length reference k-mers. If the reference output was restricted by a motifs file, only k-mers in that file contribute to the correction.

```r
ends <- read_end_motifs("sample.end_motifs.zarr")
ref_kmers <- read_ref_kmers("reference.ref_kmers.zarr")

corrected_rows <- end_motif_data_frame(
  ends,
  ref_kmers = ref_kmers,
  two_sided_correction = "joint"
)
```

The choice also determines the motif axis of corrected matrices. `"joint"` and `"split"` retain the selected full-motif axis. `"outside"` and `"inside"` create a new axis after combining counts by side. The matrix column names report the resulting motif labels in column order.

```r
corrected_matrix <- dense_corrected_counts_matrix(
  ends,
  ref_kmers,
  two_sided_correction = "outside"
)
corrected_motifs <- colnames(corrected_matrix)
```

<br>

## Length Counts

Length-count outputs are wide TSV files from `cfdna lengths`.

```r
lengths <- read_lengths("sample.length_counts.tsv.zst")

length_bins(lengths)
length_counts_matrix(lengths)
length_data_frame(lengths, value = "fraction")
length_data_frame(lengths, with_length_range = c(100L, 221L))
length_data_frame(
  lengths,
  with_length_range = c(100L, 221L),
  value = "fraction",
  denominator = "selected_bins"
)
```

- Global outputs also support `length_counts_vector(lengths)`.
- Windowed outputs support `window_metadata(lengths)` and optional `window_idxs` selection in `length_data_frame()`.
- Grouped outputs support `group_metadata(lengths)`, `group_idx()`, and optional `groups` or `group_idxs` selection.
- Length-bin selection supports `with_lengths`, `with_length_range`, and `length_bin_idxs`.
- For `fraction` and `density`, `denominator = "all_bins"` uses all length bins and `denominator = "selected_bins"` uses only the returned bins.

For windowed or grouped outputs, `max_blacklisted_fraction` filters rows by `blacklisted_fraction`. Outputs without that column only accept the default keep-all cutoff.
