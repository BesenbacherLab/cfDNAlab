<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->

# cfdna fcoverage

```text
Count positional **fragment** coverage across the genome.

In paired-end mode, only fragments with both reads present are considered. By default, the entire fragment span is counted, except for deletions and skipped regions that are not covered by the other read.

## Fragment span definition

**Paired-end**: `[forward.pos, reverse.end)`.

**Unpaired** where each read is a fragment: `[read.pos, read.end)`.

## Windowing

When specifying windows (`--by-bed` or `--by-size`), one of the following outputs is possible:

- Get the average (default) or total coverage per window.

- Get the positional coverage for the included windows only (`--by-bed` *only*). Excludes all positions that do not overlap a window from the output. Choose between: 1) Indexed: Adds the original window index as an output column and keeps duplicate positions. 2) Unique: Overlapping windows are merged to avoid duplicate positions.

Without windowing, positional coverage are outputted for the selected chromosomes.

## Blacklisting

Positions in blacklisted regions are set to `f32::NaN` (and thus not included in sums or averages).

## GC correction

Reduce the global GC bias (common technically-induced bias) in the coverage by weighting the contribution of fragments. Two options:

`--gc-file`: Weight the contribution of each fragment by its length and GC content using a precomputed correction matrix from `cfdna gc-bias`. The GC correction matrix should be calculated from the same BAM file, as the bias is sample-specific.

`--gc-tag`: Weight the contribution of each fragment by a weight saved as an aux tag in the BAM reads. Allows using external GC packages like `GCParagon` and `GCfix` (both use the tag "GC").

## Temporary files

We write temporary files to a `<output-dir>/tmp.<output-prefix>.<random>` directory to reduce memory. This directory is deleted at the end of the run. If the software is disrupted, the directory may be left behind.

## Always-on exclusion criteria

The following criteria always exclude a read:

The read is secondary, supplementary or duplicate. The read failed quality check.

**Paired-end input only**: The read or mate read is unmapped. The read is mapped to a different `tid` than the mate. The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).

Usage: fcoverage [OPTIONS] --bam <BAM> --output-dir <OUTPUT_DIR>

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
          
          Examples produce files like: `<prefix>.per_position.bedgraph.zst`, `<prefix>.per_position_per_window.tsv.zst`, `<prefix>.avg.tsv.zst`, or `<prefix>.total.tsv.zst`.
          
          [default: coverage]

      --decimals <DECIMALS>
          Decimals to round coverage to when writing `[integer]`
          
          **NOTE**: When floating point precision is not needed, all coverages are integers, we remove all decimal points!
          
          [default: 2]

      --keep-zero-runs
          Output zero-coverage runs in positional coverage outputs `[flag]`
          
          By default, only covered positions are written to the output.

      --tile-size <TILE_SIZE>
          Size of tiles to parallelize over `[integer]`
          
          Chromosomes are processed in tiles of this size to reduce memory usage.
          
          [default: 20000000]

      --per-window <PER_WINDOW>
          What to return per window `[string]`
          
          Possible values:
          
          - `"average"`: Get the average coverage per window (default).
          
          - `"total"`: Get the total coverage per window.
          
          - `"unique-positions"`: Get the positional coverage for the included windows only (`--by-bed` *only*). Overlapping windows are merged to avoid duplicate positions. Excludes all positions that do not overlap a window from the output.
          
          - `"indexed-positions"`: Get the positional coverage for the included windows only (`--by-bed` *only*). Adds the original window index as an output column and keeps duplicate positions. Excludes all positions that do not overlap a window from the output.
          
          **NOTE**: Ignored when no windows are specified.
          
          [default: average]

      --ignore-gap
          Ignore inter-mate gap `[flag]`
          
          Disable counting of the gap between reads (i.e., `[forward.end, reverse.start)`) when the two reads do not overlap.

Windows (select max. one arg.):
      --by-size <BY_SIZE>
          Window definition: a fixed window size `[integer]`
          
          Default is one global window.

      --by-bed <BY_BED>
          Window definition: a BED file of windows `[path]`

Chromosome Selection (select max. one arg.):
      --chromosomes <CHROMOSOMES>...
          Names of chromosomes to process (comma-separated or repeated). E.g. `'chr1,chr2,chr3'`.
          
          When no chromosomes are specified, it defaults to `chr1..chr22`.
          
          Specify `"all"` *as the only string* to use all present chromosomes. For BAM-backed commands this uses the BAM header order. For commands that read chromosome order from their input, this may use the input order or some other order.

      --chromosomes-file <CHROMOSOMES_FILE>
          File with chromosome names to process (one per line)

Normalization:
      --scaling-factors <SCALING_FACTORS>
          Optional path to *non-negative* scaling factors for normalizing/smoothing the genome `[path]`
          
          `.tsv` file as produced by `cfdna coverage-weights` containing a scaling factor to *multiply* by per **scaling-bin**.
          
          The scaling-bin-overlapping parts of the fragments are counted as the scaling factor of the bin (`w=sf`).
          
          ## File Requirements
          
          The TSV file **must** have a header. Column names are matched **case-insensitively**.
          
          Required columns: `chromosome`, `start`, `end`, `scaling_factor`.
          
          Coordinates are 0-based, half-open `[start, end)`.
          
          `scaling_factor` must be finite and strictly >= 0.
          
          Bins are filtered to the provided `chromosomes`.
          
          For every chromosome in `chromosomes`, bins must:
          
          - start at 0
          
          - be perfectly contiguous (no gaps, no overlaps)
          
          - end exactly at that chromosome’s length

Filtering:
      --min-fragment-length <MIN_FRAGMENT_LENGTH>
          Minimum fragment length to include `[integer]`
          
          [default: 30]

      --max-fragment-length <MAX_FRAGMENT_LENGTH>
          Maximum fragment length to include `[integer]`
          
          [default: 1000]

      --min-mapq <MIN_MAPQ>
          Minimum mapping quality to include `[integer]`
          
          [default: 30]

      --require-proper-pair
          Only count properly paired reads `[flag]`
          
          Not recommended, as we already select only inward-directed read pairs within fragment length bounds.

  -b, --blacklist <BLACKLIST>...
          Optional BED file(s) with blacklisted regions `[path]`

GC Correction (select max. one source):
      --gc-file <GC_FILE>
          Optional path to GC correction file *made from the same BAM file* with `cfdna gc-bias` `[path]`
          
          The file is usually called `gc_bias_correction.npz`.
          
          **NOTE**: Requires specifying the reference genome 2bit file as well.

      --gc-tag <GC_TAG>
          Optional aux tag to get GC weight from when using external GC correction packages `[string]`
          
          Packages like `GCParagon` and `GCfix` allow saving GC weights directly to the reads in a BAM file. They often assign a "GC" aux tag.
          
          The average per-read weight is used to count the fragment. When any of the reads have a zero-weight, the fragment gets a zero-weight.

      --drop-invalid-gc
          Whether to drop fragments where the GC correction could not be calculated `[flag]`
          
          If a GC correction weight could not be computed/retrieved for a fragment, the default is to weight it as `1.0` (no correction). If you prefer to exclude it instead, set this flag.

GC Correction:
  -r, --ref-2bit <REF_2BIT>
          Optional 2bit reference genome file [path]
          
          NOTE: Required for GC correction, otherwise ignored.
          
          E.g., "hg38.2bit" from UCSC ( https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit ).
```
