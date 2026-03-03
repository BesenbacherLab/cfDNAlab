<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->

# cfdna midpoints

```text
Count positional fragment **midpoint** coverage in groups of genomic windows.

**Midpoints**: The center of the fragment span, with ties (in even-sized fragments) randomly assigned to either the left or right mid-position to reduce rounding bias.

**Groups**: The coverage profiles are "collapsed" (summed per position) for all windows in a group. E.g., groups can be transcription factors with windows being binding sites. We then get the overall midpoint profile per transcription factor.

## Fragment span definition

**Paired-end**: `[forward.pos, reverse.end)`.

**Unpaired** where each read is a fragment: `[read.pos, read.end)`.

## Always-on exclusion criteria

The following criteria always exclude a read:

The read is secondary, supplementary or duplicate. The read failed quality check.

**Paired-end input only**: The read or mate read is unmapped. The read is mapped to a different `tid` than the mate. The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).

Usage: midpoints [OPTIONS] --bam <BAM> --output-dir <OUTPUT_DIR> --intervals <INTERVALS>

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
          
          E.g., specify to enable writing to the same output directory from multiple calls of the command.
          
          Examples produce files like: `<prefix>.midpoint_profiles.npy`.
          
          [default: sites]

  -w, --intervals <INTERVALS>
          The grouped fixed-size intervals to count within `[path]`
          
          A BED-like file of genomic intervals and their respective group names.
          
          Must be sorted by the `chromosome` and `start` coordinates, and all intervals must have the same size.
          
          Sites with the same group name are collapsed to a single profile.
          
          Columns: `chromosome, start, end, group_name`. No header.

      --length-bins <LENGTH_BINS>...
          Edges of fragment length bins to count in `[string(s)]`
          
          Accepted forms:
          
          - A single value with `start:end:step`: Creates contiguous bins from `start` to `end` (end-exclusive) in `step` increments. Example: `30:1000:10` -> bins `[30,40), [40,50), ..., [990,1000)`.
          
          - Multiple integer values interpreted as bin edges: Example: `--length-bins 30 80 150 220 500 1001` -> bins `[30,80), [80,150), ..., [500,1001)`.
          
          **NOTE**: Memory consumption increases linearly with the number of bins.
          
          [default: 30 1001]

      --tile-size <TILE_SIZE>
          Size of tiles to parallelize over `[integer]`
          
          Chromosomes are processed in tiles of this size to reduce memory usage.
          
          [default: 63000000]

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
      --min-mapq <MIN_MAPQ>
          Minimum mapping quality to include `[integer]`
          
          [default: 30]

      --require-proper-pair
          Only count properly paired reads `[flag]`
          
          This is **NOT** recommended by default as it trims the tails of the length distribution.

  -b, --blacklist <BLACKLIST>...
          Optional BED file(s) with blacklisted regions `[path]`
          
          **NOTE**: It may be an advantage to instead remove intervals that lie within half the maximum fragment length of blacklisted regions from the `--intervals` file.

      --blacklist-min-size <BLACKLIST_MIN_SIZE>
          Minimum size of blacklist intervals to load (bp) `[integer]`
          
          [default: 1]

      --blacklist-strategy <BLACKLIST_STRATEGY>
          The fragment positions that should overlap blacklisted regions for it to be excluded `[string]`
          
          Possible values: `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
          
          Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
          
          [default: any]

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

Plotting:
      --plot-groups <PLOT_GROUPS>
          Group indices to plot as midpoint profiles `[integers]`
          
          Comma separated list of zero-based group indices to plot after counting.
          
          This plotting step is intended for quick QC of the outputs. It's not optimized for publication etc. (although feel free!)
          
          [default: 0]
```
