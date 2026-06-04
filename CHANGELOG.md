# Changelog

This is the changelog for the main CLI tool. You can find the changelog for the downstream loader packages here:

 - [Python](https://github.com/BesenbacherLab/cfDNAlab/blob/main/py-cfdnalab/CHANGELOG.md)
 - [R](https://github.com/BesenbacherLab/cfDNAlab/blob/main/r-cfdnalab/NEWS.md)

<br />

## Unreleased

 - Adds `--motifs-file` to `cfdna ends` for pre-specifying the motifs to count. This allows counting larger motifs without exploding the memory.
 - In `cfdna lengths`, the `max_soft_clips` cap is only applied when `--clip-mode adjust`.
 - In `cfdna ends`, fixes theoretical bug that wrongly filtered fragments when kmer-size was larger than read sequence but reference was used as source.

<br />

## cfDNAlab 0.2.0

**BREAKING CHANGES**:

The output formats are changed to simplify downstream work in R and python. While the previous NumPy `.npy` format was easy to load in python, the user had to infer the column names, etc. We are moving towards more self-contained formats for a better and safer experience.

 - `cfdna lengths` now outputs a wide-format Zstandard-compressed TSV (`.tsv.zst`) file with one row per window/group. This is readable by `pandas` in Python and `fread(cmd = "zstd -dc file.tsv.zst")` in R.
 - `cfdna lengths` writes command settings to `<prefix>.length_settings.json`.
 - `cfdna lengths` adds `--decimals` to control written count precision.
 - `cfdna midpoints` now outputs profiles as a self-contained Zarr store at `<prefix>.midpoint_profiles.zarr`.
 - `cfdna midpoints` writes command settings to `<prefix>.midpoint_settings.json`.
 - `cfdna ends` now outputs dense and sparse motif counts as a self-contained Zarr store at `<prefix>.end_motifs.zarr`.
 - `cfdna ends` writes command settings to `<prefix>.end_settings.json`.
 - `cfdna ref-gc-bias` now writes the public reference GC package as `<prefix>.ref_gc_package.zarr`.
 - `cfdna gc-bias` now writes the public sample correction package as `<prefix>.gc_bias_correction.zarr`.

Any GC-correction packages previously computed are no longer compatible. Please recompute.

**Python and R packages for downstream loading**:

While the new `.zarr` format has many benefits, such as being self-contained and working in both R and Python, it's not the easiest to work with in its raw form. We added both an R package and a Python package with friendly loaders and data frame + array getters.

 - Adds the **Python** `cfDNAlab` package for loading the `.zarr` outputs of the CLI tool and extracting data frames, meta data, etc. See the `cfDNAlab/py-cfdnalab/` directory.
 - Adds the **R** `cfDNAlab` package for loading the `.zarr` outputs of the CLI tool and extracting data frames, meta data, etc. See the `cfDNAlab/r-cfdnalab/` directory.
 
**Other changes**:

- Adds `zarrs=0.23.11` dependency.
- Adds downstream Zarr compatibility checks for Python and R reader packages.
- We no longer expose the former `gen_cli_docs` binary that is for internal use only.

<br />

## cfDNAlab 0.1.0

Adds the following public commands:

 - `cfdna fcoverage`
 - `cfdna lengths`
 - `cfdna ends`
 - `cfdna midpoints`
 - `cfdna gc-bias`
 - `cfdna coverage-weights`
 - `cfdna fragment-count-weights`
 - `cfdna bam-to-bam`
 - `cfdna bam-to-frag` 
 - `cfdna frag-to-bam`
 - `cfdna ref-gc-bias`

Additional commands are in development and will be added in future releases.
