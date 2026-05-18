# cfDNAlab | R Loaders <img src="https://raw.githubusercontent.com/BesenbacherLab/cfDNAlab/refs/heads/main/cfdnalab_logo_little_guy_172x200_144dpi.png" align="right" height="155" />

R helpers for loading [cfDNAlab](https://github.com/BesenbacherLab/cfDNAlab) analysis outputs.

This package does not install or run the cfDNAlab command-line tool. The CLI is distributed separately as the Rust `cfdna` binary. Use this R package after running cfDNAlab to load, inspect, and reshape output files in R.

The first supported output types are midpoint and end-motif Zarr outputs: `<prefix>.midpoint_profiles.zarr` and `<prefix>.end_motifs.zarr`.

The helpers return base `data.frame` objects, R arrays, and `Matrix` sparse matrices. Convert data frames with `tibble::as_tibble()` or `data.table::as.data.table()` when you want those workflows.

Numeric indices returned by this R package are one-based, matching ordinary R indexing.

<br>

## Install

Install the current development version from GitHub:

```r
install.packages("pak")
pak::pak("cfdnalab=github::BesenbacherLab/cfDNAlab/r-cfdnalab")
```

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

profile <- profile_data_frame(
  midpoints,
  group = "LYL1",
  length_bin_idx = length_bin
)

head(profile)
```

Use `profile_array()` when you only need the count vector for one profile. Use `midpoint_array()` only when you want the full 3D count array in memory.

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
- `"size"`: rows are fixed-size genomic windows from `--window-size`. Use `windows(ends)` and the window count helpers.
- `"bed"`: rows are BED intervals. Use `windows(ends)` and the window count helpers.
- `"grouped_bed"`: rows are BED groups. Use `group_metadata(ends)`, `group_idx()`, and the group count helpers.

Sparse output keeps only non-zero counts in memory:

```r
counts <- sparse_counts_matrix(ends)
motif_counts <- sparse_data_frame_for_motif(ends, "_AA")
```

Dense output can be read as a matrix or data frame:

```r
counts <- dense_counts_matrix(ends)
motif_counts <- dense_data_frame_for_motif(ends, "_AA")
```

Dense helpers do not silently convert sparse stores. If you want a dense matrix from sparse output, pass `allow_densify = TRUE`.

```r
counts <- dense_counts_matrix(ends, allow_densify = TRUE)
```
