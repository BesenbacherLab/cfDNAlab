# Changelog

**NOTE**: This is the changelog for the **R** package that provide output loaders for the main output files from the main `cfDNAlab` command line tool. The changelog for the CLI tool is found [here](https://github.com/BesenbacherLab/cfDNAlab/blob/main/CHANGELOG).

"x.x.x.9000" is the current development version.

## r-cfDNAlab 0.1.0.9000

 - Renames `groups()` to `group_metadata()` to remove the masking with tidyverse. 
 - Adds `window_metadata()` as the R window metadata helper and removes the
   `windows()` helper.
 - Standardizes public genomic window metadata columns to `window_idx`,
   `chrom`, `start`, `end`, and `blacklisted_fraction`.
 - Adds `read_lengths()` for `cfdna lengths` TSV outputs, including length-bin
   metadata, matrix/vector getters, and long or wide data-frame reshaping.
 - Adds `max_blacklisted_fraction` filtering to row-based end-motif data-frame
   helpers.

## r-cfDNAlab 0.1.0

 - Adds zarr loaders for the outputs of `cfdna midpoints` and `cfdna ends`.
