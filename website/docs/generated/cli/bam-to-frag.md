<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->

# cfdna bam-to-frag

```text
Write the fragments from a BAM file to a finaleDB-style frag file.

Information in the `.frag.tsv` file:

- **Chromosome**

- **Start**: forward.pos

- **End**: reverse.end

- **MapQ**: Minimum mapping quality for the two reads

- **Strand**: The strand alignment of read1

AND, when either/both `--gc-file` and `--scaling-factors` are specified:

- **GC Weight**: The multiplicative weight needed to correct for GC bias.

- **Scaling Weight**: The multiplicative weight needed to perform genomic smoothing.

Note: When GC correction is not specified but genomic scaling is, the sixth column is the scaling weight.

The accompanying `*.frag.header.tsv` file has the matching column names.

Fragments are sorted by `(chromosome, start, end)`, using the chromosome order in `--chromosomes`.

## Always-on exclusion criteria

The following criteria always exclude a read:

The read is secondary, supplementary or duplicate. The read failed quality check.

**Paired-end input only**: The read or mate read is unmapped. The read is mapped to a different `tid` than the mate. The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).

Usage: bam-to-frag [OPTIONS] --bam <BAM> --output-dir <OUTPUT_DIR>

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
          Prefix for output file (e.g., a sample name) `[string]`
          
          E.g., specify to enable writing to the same output directory from multiple calls to this software.
          
          Examples produce files like: `<prefix>.frag.tsv.gz`,
          
          [default: fragments]

Windows:
      --by-bed <BY_BED>
          Intervals to keep overlapping fragments from `[path]`

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
          
          Defaults to 0 to allow making filtering decisions downstream.
          
          [default: 0]

      --require-proper-pair
          Only count properly paired reads `[flag]`
          
          This is **NOT** recommended by default, as it trims the tails of the length distribution. It may be useful to match the files in FinaleDB.
          
          Note, that we only keep inward-directed fragments within a specified length range, so there's no real need for proper-pair filtering.

  -b, --blacklist <BLACKLIST>...
          Optional BED file(s) with blacklisted regions `[path]`

      --blacklist-min-size <BLACKLIST_MIN_SIZE>
          Minimum size of blacklist intervals to load (bp) `[integer]`
          
          [default: 1]

      --blacklist-strategy <BLACKLIST_STRATEGY>
          The fragment positions that should overlap blacklisted regions for it to be excluded `[string]`
          
          Possible values: `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
          
          Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
          
          [default: any]

GC Correction:
      --gc-file <GC_FILE>
          Optional path to GC correction file *made from the same BAM file* with `gc-bias` `[path]`
          
          The file is usually called `gc_bias_correction.npz`.
          
          **NOTE**: Requires specifying the reference genome 2bit file as well.

      --drop-invalid-gc
          Whether to drop fragments where the GC correction could not be calculated `[flag]`
          
          If a GC correction weight could not be computed for a fragment, the default is to weight it as `1.0` (no correction). If you prefer to exclude it instead, set this flag.

  -r, --ref-2bit <REF_2BIT>
          Optional 2bit reference genome file [path]
          
          NOTE: Required for GC correction, otherwise ignored.
          
          E.g., "hg38.2bit" from UCSC ( https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit ).
```
