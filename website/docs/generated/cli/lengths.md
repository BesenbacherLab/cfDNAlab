<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->

# cfdna lengths

```text
Count fragment lengths in a BAM-file.

Writes an `.npy` file with shape (# windows, # lengths).

## Fragment length definition

**Paired-end**: `end(reverse) - start(forward)`.

**Unpaired** where each read is a fragment: `end(read) - start(read)`.

See also `--indel-mode` for adjusting the length to present indels.

## GC correction

Weight the contribution of each fragment based on their GC contents.

Note: The GC percentage is calculated from the full genomic coordinates (does not consider `--indel-mode`).

The length-dimension of the original correction matrix is averaged out with a specifiable weighting scheme (`--gc-length-weighting`).

## Genomic smoothing (--scaling-factors)

Weight how genomic regions contribute to the length distribution, e.g. to reduce the influence of copy number alterations. This weights the contribution of each fragment by region-wise precomputed scaling factors.

Can be precomputed with `cfdna coverage-weights`.

## Window assignment

By default, fragments are counted by their window-overlap fraction. That is, most fragments are counted as `1.0` (before correction/scaling), while fragments overlapping the edge of a window are counted as the fraction it overlaps the window (`< 1.0`).

For consecutive non-overlapping windows, this conserves the total mass, as an edge-overlapping fragment will count `f` in one window and `1-f` in the other window.

To get base-weighted counts (i.e. coverage in the window), you can multiply the output counts by their lengths (`C'[L] = L * C[L]`; Remember to account for the minimum fragment length offset).

Other options include counting the full fragment if the *fragment midpoint* or a given *proportion* of positions overlaps the window.

## Blacklisting

Ignores fragments that overlap blacklisted regions with a given proportion.

## Always-on exclusion criteria

The following criteria always exclude a read:

The read is secondary, supplementary or duplicate. The read failed quality check.

**Paired-end input only**: The read or mate read is unmapped. The read is mapped to a different `tid` than the mate. The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).

Usage: lengths [OPTIONS] --bam <BAM> --output-dir <OUTPUT_DIR>

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
          
          Examples produce files like: `<prefix>.length_counts.npy`
          
          [default: lengths]

      --indel-mode <INDEL_MODE>
          How to handle insertions and deletions in fragments `[string]`
          
          Deletions: Both 'D' and 'N' in the cigar string are considered deletions.
          
          Possible values:
          
          - `"ignore"`: Ignore whether indels are present or not.
          
          Lengths are calculated from the reference coordinates `end(reverse) - start(forward)`.
          
          - `"adjust"`: Adjust the reference length by the observed insertions and deletions (we cannot adjust in the mate-gap).
          
          For bases only covered by a single read, all insertions and deletions are adjusted for.
          
          In the mate-overlap, only adjust when both reads show the indel at the same reference position.
          
          Deletions: subtract the reference bases deleted in both reads.
          
          Insertions: add the shortest insertion length per position.
          
          **NOTE**: Blacklist exclusion and calculation of scaling weights (--scaling-factors) use the full reference span.
          
          - `"skip"`: Skip fragments with any insertion or deletion present.
          
          [default: ignore]

      --tile-size <TILE_SIZE>
          Size of tiles to parallelize over `[integer]`
          
          Chromosomes are processed in tiles of this size to reduce memory usage.
          
          [default: 20000000]

Windows (select max. one arg.):
      --by-size <BY_SIZE>
          Window definition: a fixed window size `[integer]`
          
          Default is one global window.

      --by-bed <BY_BED>
          Window definition: a BED file of windows `[path]`

Window Assignment:
      --assign-by <ASSIGN_BY>
          The **fragment positions** that should overlap a window for it to be counted in that window, OR the option to count the fraction of overlapping bases `[string]`
          
          Possible values: `"count-overlap"`, `"any"`, `"all"`, `"midpoint"`, or `"proportion=<threshold>"`
          
          `'count-overlap'`: Count up the fraction of overlapping fragment bases.
          
          Example of proportion: `--assign-by proportion=0.2` (no space around `=`)
          
          Midpoints for even-sized fragments are randomly selected as either the left or right base to avoid bias.
          
          **NOTE**: In the rare case where windows are smaller than fragments, it's still the proportion of the fragment positions that overlap that is considered. If the window size is 30% of the fragment size, that fragment cannot overlap more than 30%.
          
          **NOTE**: Ignored when no windows are specified.
          
          [default: count-overlap]

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
          
          This is **NOT** recommended by default as it trims the tails of the length distribution.

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

      --gc-length-weighting <GC_LENGTH_WEIGHTING>
          How to weight the fragment length bins when estimating the global GC bias correction `[string]`
          
          To GC correct a fragment length distribution, the correction weights should be **length-agnostic**.
          
          The default `fragment-length x GC` matrix has one correction curve per length bin, so using it would preserve the original length distribution (assuming we're correcting the same fragments seen by `cfdna gc-bias`).
          
          We therefore average out the fragment length dimension to get a single, length-agnostic GC bias curve.
          
          We have three weighting options when averaging the fragment-length-wise correction curves:
          
          - `"equal"` weighting (default): Weight every length bin the same.
          
          Keeps the correction independent of the distribution we are trying to estimate.
          
          Downside: Rare fragment length bins contribute the same as the most present fragment lengths.
          
          For low-coverage BAM files, this *could* make the correction more volatile to outliers.
          
          - `"coverage"`-based weighting: Weight lengths by how often they were observed in `cfdna gc-bias`.
          
          This should work better for the majority of the observed fragments **BUT**:
          
          Downside: **Biases** the correction based on the length distribution we are trying to estimate.
          
          - `"max-coverage"` weighting: Use the GC curve for the most-observed fragment length bin.
          
          [default: equal]

  -r, --ref-2bit <REF_2BIT>
          2bit reference genome file [path]
          
          NOTE: Required when specifying `--gc-file`.
          
          E.g., "hg38.2bit" from UCSC (https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit).
```
