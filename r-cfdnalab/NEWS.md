# Changelog

**NOTE**: This is the changelog for the **R** package that provide output loaders for the main output files from the main `cfDNAlab` command line tool. The changelog for the CLI tool is found [here](https://github.com/BesenbacherLab/cfDNAlab/blob/main/CHANGELOG.md).

<br />

## r-cfDNAlab 0.4.0

- Adds zarr loader for the outputs of `cfdna ref-kmers`.
- Adds reference correction in the loader for outputs of `cfdna ends` using the `ref-kmers` frequencies.

<br />

## r-cfDNAlab 0.3.0

 - Adds support for schema v2 end-motif Zarr stores written by `cfdna ends --motifs-file`, including outputs where the motif axis contains user-defined motif groups instead of concrete motif labels.
 - Adds the common `motifs`/`motif_idxs`, `groups`/`group_idxs`, and `window_idxs` selectors to the `sparse_counts_matrix()` and `dense_counts_matrix` functions.
 - `ends`: Improves error message when sparse data has no counts.

<br />

## r-cfDNAlab 0.2.0

This version cleans up the API and introduces `read_lengths` for reading the output of `cfdna lengths`.

**BREAKING CHANGES**:

 - Renames `groups()` to `group_metadata()` to remove the masking with tidyverse.
 - Adds `window_metadata()` as the R window metadata helper and removes the
   `windows()` helper.
 - Standardizes public genomic window metadata columns to `window_idx`,
   `chrom`, `start`, `end`, and `blacklisted_fraction`.
 - Adds `read_lengths()` for `cfdna lengths` TSV outputs, including length-bin
   metadata, matrix/vector getters, and long or wide data frame reshaping.
 - Adds row selectors to `length_counts_matrix()` for windowed and grouped
   length-count outputs.
 - Adds length-bin selectors `with_lengths`, `with_length_range`, and
   `length_bin_idxs` to length-count outputs, and `with_length_range` to
   midpoint outputs.
 - Adds `denominator` to length-count data frames so `fraction` and `density`
   can be normalized over all bins or only selected bins.
 - Removes the internal `count_column` field from `length_bins()` output. Wide
   length data frames still preserve the source count column labels.
 - Adds `max_blacklisted_fraction` filtering to row-based end-motif data frame
   helpers.
 - Replaces the R midpoint and end-motif data frame families with
   `midpoint_data_frame()` and `end_motif_data_frame()`, using plural selectors
   such as `groups`, `group_idxs`, `window_idxs`, `motifs`, and `motif_idxs`.

<br />

## r-cfDNAlab 0.1.0

 - Adds zarr loaders for the outputs of `cfdna midpoints` and `cfdna ends`.
