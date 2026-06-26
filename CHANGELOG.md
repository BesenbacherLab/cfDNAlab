# Changelog

This is the changelog for the main CLI tool. You can find the changelog for the downstream loader packages here:

 - [Python](https://github.com/BesenbacherLab/cfDNAlab/blob/main/py-cfdnalab/CHANGELOG.md)
 - [R](https://github.com/BesenbacherLab/cfDNAlab/blob/main/r-cfdnalab/NEWS.md)

<br />

## Unreleased

**BREAKING CHANGES**:
- Grouped motifs-file outputs now order motif groups alphabetically by group name.

**Other changes**:

 - Adds `cfdna ref-kmers` command for counting k-mers in the reference assembly for downstream normalization of k-mer counts.
   - Library: Adds Rust loaders for `ref-kmers` output.
 - Adds minimal validation of BED files to catch obviously non-BED formats.
 - Library: The Rust `fcoverage` run result now includes grouped `group_index.tsv` filepaths in `output_files()` when they are written.
 - Library: Adds missing blacklist setters on the Rust `LengthsConfig` API.
 
 
<br />

## cfDNAlab 0.5.0

**BREAKING CHANGES**:

 - `cfdna fcoverage --normalize-by-length` is now actually called `--normalize-by-length` instead of `--normalize-by-length-mode`, matching the guides and documentation.
 - The Rust crate feature `cmd_ref_gc_bias` has been removed. Enable `cmd_gc_bias` to compile both the `gc-bias` and `ref-gc-bias` command APIs.

**Other changes**:

 - Pins the build-time Clang/LLVM conda packages to version 21.x in installation instructions to avoid broken `rust-htslib` bindings with libclang version 22.
 - `cfdna fragment-count-weights` and `cfdna coverage-weights` shows progress bar by default in the internal call to `cfdna fcoverage`.
 - CLI runs now log the *equivalent* full `cfdna <COMMAND> <OPTION>` command call for reproducibility and transparency.
 - The Rust library can render *equivalent* command calls from exported configs without enabling the `cli` feature. This adds the `RunOptions.log_equivalent_cli` field.
 - Adds initial output loaders for `lengths`, `fcoverage`, `ends`, and `midpoints` to the Rust library.
 - General feature gating improvements for an improved library experience.

<br />

## cfDNAlab 0.4.0

**BREAKING CHANGES**:

 - `cfdna fcoverage --per-window summary-stats*` output TSV columns are renamed and deduplicated.
 - For `cfdna fcoverage --normalize-by-length`, aggregate TSV headers now use `fragment_mass` instead of `coverage` for the length-normalized signal columns.

**Other changes**:

 - Allows reading BAM files from URLs via the `curl` feature in `rust-htslib`. This adds `url` as dependency.
 - Updates `rust-htslib` dependency to `1.0.0` and reduces its feature set for fewer dependencies.
 - Enables `libdeflate-sys` for faster BAM reading.
 - `bam-to-bam` and `frag-to-bam` writes `.bam.bai` index files alongside the BAM outputs.
 - Requires `bindgen >=0.69.5, <0.70` to remove issue with `0.69.4` when users don't install with lockfile.
 - Adds `--locked` to install commands in installation instructions.

<br />

## cfDNAlab 0.3.0

**BREAKING CHANGES**:

 - Major refactor of the rust library. Makes a clearer boundary of what is public/private. While cfDNAlab is primarily a CLI tool, the library side needed a clean up. This does not affect CLI usage.

**Other changes**:

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
