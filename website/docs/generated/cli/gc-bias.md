<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->

# cfdna gc-bias

Calculate a multiplicative GC correction matrix based on the GC fraction and length of fragments in a BAM-file.

The observed distribution of cfDNA fragments is corrected to a precomputed reference bias.

Requirements: Please precompute the reference GC bias with `cfdna ref-gc-bias`. This file can be reused for all samples aligned to the same assembly.

**NOTE**: This command is highly flexible, enabling experimentation. The default values have been tuned and should be useful in most use cases. Start with the example below.

## Interpolations

The most extreme GC and shortest-length bins get interpolated corrections based on neighbours to avoid extreme corrections due to sparsity.

The combinations of GC fractions and fragment lengths that are either theoretically unobservable or *very* rarely observed in the **reference genome** are interpolated from surrounding counts. Other combinations with post-smoothing zero counts in the *cfDNA* remains zero in the correction matrix. The final correction matrix thus works for all possible GC x Length combinations.

## Fragment length definition

**Paired-end**: `end(reverse) - start(forward)`.

**Unpaired** where each read is a fragment: `end(read) - start(read)`.

## Windowing

Technical GC bias is assumed to be a "global" bias. To control how each region of the genome (which may have amplified/reduced coverage) contributes to the calculation of this global bias, we can calculate the bias in genomic windows and combine them via weighted averaging: The counts of each window are divided by their window-mean and scaled by the number of valid ACGT positions in the window. The windows are then averaged.

## Example

```bash

cfdna gc-bias --bam {BAM_FILE} --output-dir {PATH}/gc_bias \

--ref-2bit {PATH}/hg38.2bit \ # Or some other assembly

--ref-gc-dir {REFERENCE_GC_DIRECTORY} \

--min-fragment-length 30 --max-fragment-length 1000 \

--blacklist {PATH}/encode_blacklist.bed # Or some other blacklist(s)

```

Besides these arguments, the default values should work in most cases.

## Always-on exclusion criteria

The following criteria always exclude a read:

The read is secondary, supplementary or duplicate. The read failed quality check.

**Paired-end input only**: The read or mate read is unmapped. The read is mapped to a different `tid` than the mate. The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).


<hr class="cli-usage-separator" />

## Usage

`cfdna gc-bias [OPTIONS] --bam <BAM> --output-dir <OUTPUT_DIR> --ref-2bit <REF_2BIT> --ref-gc-dir <REF_GC_DIR>`

## Options

- `-h, --help`

  Print help (see a summary with '-h')
  

## Core

- `-i, --bam <BAM>`

  Indexed, coordinate-sorted BAM input file `[path]`
  
  Can be either **paired-end** or **unpaired** (set `--reads-are-fragments`). Unpaired assumes the reads span their fragments exactly (so read size is fragment size).
  

- `-o, --output-dir <OUTPUT_DIR>`

  Output directory for results `[path]`
  

- `-t, --n-threads <N_THREADS>`

  Number of threads to use (increases RAM usage) `[integer]`
  
  Defaults to the number of available CPU cores (-1).
  
  [default: 7]
  

- `--reads-are-fragments`

  The input has one read per fragment and the **read spans exactly the full fragment** (e.g. Nanopore) `[flag]`
  
  Each aligned read is treated as a fragment spanning its aligned reference interval `[pos, reference_end)`. This uses the mapped span only (soft clips excluded).
  
  Cannot be combined with `--require-proper-pair` (when available).
  

- `-r, --ref-2bit <REF_2BIT>`

  2bit reference genome file [path]
  
  E.g., "hg38.2bit" from UCSC ( https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit ).
  

- `--ref-gc-dir <REF_GC_DIR>`

  Path to directory with reference GC bias files to correct against `[path]`
  
  Precompute with `cfdna ref-gc-bias`. The directory must include all files created by `cfdna ref-gc-bias`.
  

- `--tile-size <TILE_SIZE>`

  Size of tiles to parallelize over `[integer]`
  
  Chromosomes are processed in tiles of this size to reduce memory usage.
  
  [default: 10000000]
  

- `--save-intermediates`

  Whether to save key intermediate files for inspecting the correction process `[flag]`
  

## Windows (select max. one arg.)

- `--by-size <BY_SIZE>`

  Window definition: a fixed window size `[integer]`
  
  Default window definition is `--by-size 100000` window.
  

- `--by-bed <BY_BED>`

  Window definition: a BED file of windows `[path]`
  

- `--global`

  Window definition: one global window `[flag]`
  

## Window Assignment

- `--assign-by <ASSIGN_BY>`

  The **fragment positions** that should overlap a window for it to be counted in that window, OR the option to count the fraction of overlapping bases `[string]`
  
  Possible values: `"count-overlap"`, `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
  
  `'count-overlap'`: Count up the fraction of overlapping fragment bases.
  
  Example of proportion: `--assign-by proportion=0.2` (no space around `=`)
  
  Midpoints for even-sized fragments are randomly selected as either the left or right base to avoid bias.
  
  **NOTE**: In the rare case where windows are smaller than fragments, it's still the proportion of the fragment positions that overlap that is considered. If the window size is 30% of the fragment size, that fragment cannot overlap more than 30%.
  
  **NOTE**: Ignored when no windows are specified.
  
  [default: count-overlap]
  

## Chromosome Selection (select max. one arg.)

- `--chromosomes <CHROMOSOMES>...`

  Names of chromosomes to process (comma-separated or repeated). E.g. `'chr1,chr2,chr3'`.
  
  When no chromosomes are specified, it defaults to `chr1..chr22`.
  
  Specify `"all"` *as the only string* to use all present chromosomes. For BAM-backed commands this uses the BAM header order. For commands that read chromosome order from their input, this may use the input order or some other order.
  

- `--chromosomes-file <CHROMOSOMES_FILE>`

  File with chromosome names to process (one per line)
  

## Binning

- `--min-length-bin-mass <MIN_LENGTH_BIN_MASS>`

  Minimum percentage of counts to have in each length bin `[float]`
  
  Greater than 0, lower than 100. Default is 0.5% (i.e., a max. of 200 bins).
  
  [default: 0.5]
  

- `--min-length-bin-width <MIN_LENGTH_BIN_WIDTH>`

  Minimum number of fragment lengths per fragment length bin `[float]`
  
  Reduces sparsity-related issues in ultra low-coverage samples.
  
  [default: 3]
  

- `--min-gc-bin-mass <MIN_GC_BIN_MASS>`

  Minimum percentage of counts to have in each GC contents bin `[float]`
  
  Greater than 0, lower than 100. Default is 1% (i.e., a max. of 100 bins).
  
  [default: 1]
  

- `--num-extreme-gc-bins <NUM_EXTREME_GC_BINS>`

  Number of extreme GC bins (`--min_gc_bin_mass`) from each side to interpolate from neighbouring corrections `[float]`
  
  The most extreme GC fractions are very sparsely observed. This can lead to extreme corrections. Set the number of bins from each side where we interpolate a correction based on the neighbouring corrections. The default of 1 should be fine but this can be tuned via visualization of the created correction matrix and intermediate files (`--save-intermediates`).
  
  [default: 1]
  

- `--num-short-length-bins <NUM_SHORT_LENGTH_BINS>`

  Number of the **shortest** fragment length bins (`--min_length_bin_mass`) to interpolate from neighbouring corrections `[float]`
  
  The shortest fragment lengths can be very sparsely observed. This can lead to extreme corrections. Set the number of short-length bins where we interpolate a correction based on the neighbouring corrections. With the default minimum fragment length setting in `cfdna ref-gc-bias` (30bp), the default of 1 should be fine. This can be tuned via visualization of the created correction matrix and intermediate files (`--save-intermediates`).
  
  [default: 1]
  

## Filtering

- `-b, --blacklist <BLACKLIST>...`

  Optional BED file(s) with blacklisted regions `[path]`
  
  Masking: Blacklisted positions are set to 'N' in the reference sequence the GC fraction is calculated from. See the `Minimum ACGT` options for when to ignore a fragment with too few ACGT (non-'N' and non-blacklisted) bases.
  
  NOTE: Ensure the same positions were blacklisted when calculating the reference bias (`cfdna ref-gc-bias`).
  

- `--min-mapq <MIN_MAPQ>`

  Minimum mapping quality to include `[integer]`
  
  [default: 30]
  

- `--require-proper-pair`

  Only count properly paired reads `[flag]`
  
  This is **NOT** recommended by default as it trims the tails of the length distribution.
  

## Minimum ACGT

- `--min-window-acgt-pct <MIN_WINDOW_ACGT_PCT>`

  Minimum percentage of ACGT positions in a **window** to consider it in the bias estimation `[integer]`
  
  If you believe windows that are mostly blacklisted may be too noisy in their remaining positions, use this to threshold to remove them from the analysis.
  
  [default: 10]
  

## Outliers

- `--outlier-method <OUTLIER_METHOD>`

  Handle extreme correction factors to avoid unstable weights `[string]`
  
  Options:
  
    - `none`: Disable outlier handling.
  
    - `quantile`: Clamp using `--outlier-quantiles` (one symmetric value or two explicit values).
  
    - `iqr`, `stddev`, `mad`: Use the corresponding rule with multiplier `--outlier-k`.
  
  **NOTE**: After outlier detection, correction values are further clipped at `[0.1, 10.0]`.
  
  [default: iqr]
  [possible values: none, quantile, iqr, stddev, mad]
  

- `--outlier-scope <OUTLIER_SCOPE>`

  Whether to detect outliers per fragment length or across the full matrix `[string]`
  
    - `per-length`: Detect separately per fragment length.
  
    - `global`: Detect from the full correction matrix.
  
  [default: global]
  [possible values: per-length, global]
  

- `--outlier-quantiles <OUTLIER_QUANTILES>...`

  Quantiles for `quantile` outlier detection `[float or float,float]`
  
  Used when `--outlier-method quantile`. Provide one value to apply symmetrically (`q` -> lower=`q`, upper=`1-q`) or two values for explicit `lower,upper`.
  
  [default: 0.03 0.97]
  

- `--outlier-k <OUTLIER_K>`

  Multiplier `k` for `iqr`, `stddev`, or `mad` outlier detection `[float]`
  
  Used when `--outlier-method` is one of `iqr`, `stddev`, or `mad`.
  
  [default: 3]

