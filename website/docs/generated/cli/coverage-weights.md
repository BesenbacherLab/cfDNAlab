<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->

# cfdna coverage-weights

````text
Extract fragment coverage in large genomic bins ("megabins") with a rolling window and calculate normalizing scaling factors for smoothing the genome.

Outputs scaling factors per stride to allow other methods to apply the normalization (by weighting fragment counts).

The scaling factors are *inverted*, so normalization becomes multiplication. Zero-valued coverages lead to zero-valued scaling factors. Non-zero factors have `mean == 1.0`.

## Coverage

The full fragment span is counted without consideration of deletions and gaps. This is fine for genome-scale normalization that reduces relative changes in coverage across the genome.

## Fragment span definition

**Paired-end**: `[forward.pos, reverse.end)`.

**Unpaired** where each read is a fragment: `[read.pos, read.end)`.

## Smoothing

Smoothing is performed as a triangular moving average, calculating a weighted average of coverages from all bins overlapping a stride.

### Example

Assuming a bin-size of 6 and stride size of 2 (normally defaults to 5Mb and 0.5Mb respectively).

**Stride bins** (fixed along genome, each with an average coverage):

`[A] [B] [C] [D] [E] [F] [G] ...`

**Overlapping megabins** (`MB*`) (each covers 3 stride-bins). **`W_D`**, the number of overlapping megabins, is the (unnormalized) weight of each stride-bin in the weighted-average coverage for stride-bin `D`:

```text

MB1: [A][B][C]

MB2:    [B][C][D]

MB3:       [C][D][E]

MB4:          [D][E][F]

MB5:             [E][F][G]

W_D: [0][1][2][3][2][1][0]

```

At chromosome edges, the weights are truncated (e.g., `W_D: [2][3][2][1][0]`).

The weights are normalized by their sum (after potential truncation at edges).

## Always-on exclusion criteria

The following criteria always exclude a read:

The read is secondary, supplementary or duplicate. The read failed quality check.

**Paired-end input only**: The read or mate read is unmapped. The read is mapped to a different `tid` than the mate. The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).

Usage: coverage-weights [OPTIONS] --bam <BAM> --output-dir <OUTPUT_DIR>

Options:
  -h, --help
          Print help (see a summary with '-h')

Core:
  -i, --bam <BAM>
          Indexed, coordinate-sorted BAM input file `[path]`
          
          Can be either **paired-end** or **unpaired** (set `--reads-are-fragments`). Unpaired assumes the reads span their fragments exactly (so read size is fragment size).

  -o, --output-dir <OUTPUT_DIR>
          Output directory for results `[path]`

  -t, --n-threads <N_THREADS>
          Number of threads to use (increases RAM usage) `[integer]`
          
          Defaults to the number of available CPU cores (-1).
          
          [default: 7]

      --reads-are-fragments
          The input has one read per fragment and the **read spans exactly the full fragment** (e.g. Nanopore) `[flag]`
          
          Each aligned read is treated as a fragment spanning its aligned reference interval `[pos, reference_end)`. This uses the mapped span only (soft clips excluded).
          
          Cannot be combined with `--require-proper-pair` (when available).

  -x, --output-prefix <OUTPUT_PREFIX>
          Prefix for output files (e.g., a sample name) `[string]`
          
          E.g., specify to enable writing to the same output directory from multiple calls to this software.
          
          Examples produce files like: `<prefix>.scaling_factors.tsv`
          
          [default: normalize_genome]

Filtering:
      --bin-size <BIN_SIZE>
          Size (bp) of large genomic bins to calculate coverage in [integer]
          
          Larger values lead to a more smooth coverage across the genome.
          
          **NOTE**: The normalizing scaling factors are calculated per stride-sized overlap of these bins. Technically, we only count the coverage per stride-sized bin and then calculate the overlap with a triangular weighting scheme.
          
          [default: 5000000]

      --stride <STRIDE>
          Size (bp) of stride [integer]
          
          **NOTE**: `--bin-size` must be divisible by `stride`. I.e., `--bin-size % stride` == 0`.
          
          A normalizing scaling factor is calculated per stride as the (inverse) weighted average coverage of the overlapping large-scale bins.
          
          Smaller values lead to a higher precision in the downstream normalization but also require saving a larger BED file in the end (one line per stride-bin) and take longer to compute.
          
          [default: 500000]

      --min-fragment-length <MIN_FRAGMENT_LENGTH>
          Minimum fragment length to include `[integer]`
          
          [default: 30]

      --max-fragment-length <MAX_FRAGMENT_LENGTH>
          Maximum fragment length to include `[integer]`
          
          [default: 1000]

      --min-mapq <MIN_MAPQ>
          Minimum mapping quality to include [integer]
          
          [default: 30]

      --require-proper-pair
          Only count properly paired reads [flag]

  -b, --blacklist <BLACKLIST>...
          Optional BED file(s) with blacklisted regions [path]

      --blacklist-min-size <BLACKLIST_MIN_SIZE>
          Minimum size of blacklist intervals to load (bp) [integer]
          
          [default: 1]

      --blacklist-strategy <BLACKLIST_STRATEGY>
          The fragment positions that should overlap blacklisted regions for it to be excluded [string]
          
          Possible values: `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"` [string]
          
          Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
          
          [default: any]

Chromosome Selection (select max. one arg.):
      --chromosomes <CHROMOSOMES>...
          Names of chromosomes to process (comma-separated or repeated). E.g. `'chr1,chr2,chr3'`.
          
          When no chromosomes are specified, it defaults to `chr1..chr22`.
          
          Specify `"all"` *as the only string* to use all present chromosomes. For BAM-backed commands this uses the BAM header order. For commands that read chromosome order from their input, this may use the input order or some other order.

      --chromosomes-file <CHROMOSOMES_FILE>
          File with chromosome names to process (one per line)
````
