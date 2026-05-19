# cfDNAlab | R Loaders <img src="https://raw.githubusercontent.com/BesenbacherLab/cfDNAlab/refs/heads/main/cfdnalab_logo_little_guy_172x200_144dpi.png" align="right" height="155" />

R helpers for loading [cfDNAlab](https://github.com/BesenbacherLab/cfDNAlab) analysis outputs.

This package does not install or run the cfDNAlab command-line tool. The CLI is distributed separately as the Rust `cfdna` binary. Use this R package after running cfDNAlab to load, inspect, and reshape output files in R.

The first supported output types are midpoint and end-motif Zarr outputs plus length-count TSV outputs: `<prefix>.midpoint_profiles.zarr`, `<prefix>.end_motifs.zarr`, and `<prefix>.length_counts.tsv.zst`.

The helpers return base `data.frame` objects, R arrays, and `Matrix` sparse matrices. Convert data frames with `tibble::as_tibble()` or `data.table::as.data.table()` when you want those workflows.

Numeric indices returned by this R package are one-based, matching ordinary R indexing.

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

length_bin <- length_bin_idx(midpoints, 167)

profile <- midpoint_data_frame(
  midpoints,
  groups = "LYL1",
  length_bin_idxs = length_bin
)

head(profile)
```

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

<br>

## Length Counts

Length-count outputs are wide TSV files from `cfdna lengths`.

```r
lengths <- read_lengths("sample.length_counts.tsv.zst")

length_bins(lengths)
length_counts_matrix(lengths)
length_data_frame(lengths, value = "fraction")
```

- Global outputs also support `length_counts_vector(lengths)`.
- Windowed outputs support `window_metadata(lengths)` and optional `window_idxs` selection in `length_data_frame()`.
- Grouped outputs support `group_metadata(lengths)`, `group_idx()`, and optional `groups` or `group_idxs` selection.

For windowed or grouped outputs, `max_blacklisted_fraction` filters rows by `blacklisted_fraction`. Outputs without that column only accept the default keep-all cutoff.
