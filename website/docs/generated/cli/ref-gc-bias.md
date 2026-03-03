<!-- AUTO-GENERATED FILE - DO NOT EDIT -->
<!-- Source: cfdna Clap config and command tree -->

# cfdna ref-gc-bias

```text
Build a reference GC bias table for cfDNA correction.

Samples `n_positions` across all chromosomes and counts GC for every fragment length in range (optionally trimmed in ends). Creates one genome-wide GC-by-length table that downstream GC bias correction uses as the expected bias. If you provide a BED file via `--by-bed`, overlapping intervals are merged and counting is limited to those bases. Problematic regions can be excluded via a blacklist. Otherwise, the full genome is used.

This command never produces per-window outputs. Use `ref-gc-counts` if you need window-level counts. After counting, the table is smoothed length-wise and converted to GC percentages. A support mask flags bins with too few counts per megabase (including theoretically unobservable GC-by-length combinations), and the sparse bins are interpolated using neighbours.

Usage: ref-gc-bias [OPTIONS] --ref-2bit <REF_2BIT> --output-dir <OUTPUT_DIR>

Options:
  -h, --help
          Print help (see a summary with '-h')

Core:
  -r, --ref-2bit <REF_2BIT>
          2bit reference genome file [path]
          
          E.g., "hg38.2bit" from UCSC ( https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit ).

  -o, --output-dir <OUTPUT_DIR>
          Output directory for results [path]

  -t, --n-threads <N_THREADS>
          Number of threads to use (increases RAM usage) [integer]
          
          Defaults to the minimum of 22 (one thread per chromosome) and the number of available CPU cores (-1).
          
          [default: 7]

      --n-positions <N_POSITIONS>
          Number of genomic starting positions to sample [integer]
          
          The positions are uniformly sampled across the chromosomes with the GC of each fragment length being counted from those same starting positions.
          
          **NOTE**: Sampling is independent of windowing and blacklisting! The per-length-sum of the output counts may thus be significantly lower than the specified `n_positions` and different between lengths. **TIP**: Add 20% extra starting positions than you think you need, since blacklisting likely removes a big chunk of them.
          
          [default: 500000000]

      --seed <SEED>
          Seed for sampling of start positions `[integer]`
          
          Use this to reproduce identical reference GC outputs across runs.

      --end-offset <END_OFFSET>
          Number of bases to exclude from each fragment end `[integer]`
          
          The nucleotides in the cfDNA fragment ends can reflect biological biases (e.g., DNase activity). This argument allows isolating the GC correction from this signal.
          
          The default of `10 bp` is based on the GCfix paper by Rahman et al. 2025.
          
          [default: 10]

      --skip-interpolation
          Whether to skip the interpolation of zero-counts `[flag]`
          
          By default, `0`s are interpolated **independently per fragment length**. The assumption is that 0s are caused due to the GC content not being possible to observe with a given fragment length (e.g., a fragment length of 47 can never achieve a 99% GC). To avoid errors from this in downstream use, we use polynomial interpolation based on the neighbourhood of non-zero counts.

      --tile-size <TILE_SIZE>
          Size of tiles to process the reference in `[integer]`
          
          Chromosomes are processed in tiles of this size to reduce memory usage.
          
          [default: 10000000]

Windows:
      --by-bed <BY_BED>
          BED file with regions to include `[path]`
          
          We count at the **unique positions** included in the specified intervals.

Chromosome Selection (select max. one arg.):
      --chromosomes <CHROMOSOMES>...
          Names of chromosomes to process (comma-separated or repeated). E.g. `'chr1,chr2,chr3'`.
          
          When no chromosomes are specified, it defaults to `chr1..chr22`.
          
          Specify `"all"` *as the only string* to use all present chromosomes. For BAM-backed commands this uses the BAM header order. For commands that read chromosome order from their input, this may use the input order or some other order.

      --chromosomes-file <CHROMOSOMES_FILE>
          File with chromosome names to process (one per line)

Filtering:
  -b, --blacklist <BLACKLIST>...
          Optional BED file(s) with blacklisted regions [path]
          
          We count no fragment intervals that overlap a blacklisted base. This results in a lower count for long fragment lengths, which is not a problem due to length-wise normalization in the downstream `cfdna gc-bias` command.

      --min-fragment-length <MIN_FRAGMENT_LENGTH>
          Minimum fragment length to include `[integer]`
          
          [default: 30]

      --max-fragment-length <MAX_FRAGMENT_LENGTH>
          Maximum fragment length to include `[integer]`
          
          [default: 1000]

Smoothing:
      --smoothing-sigma <SMOOTHING_SIGMA>
          Standard deviation for Gaussian kernel that smoothes raw GC counts for each fragment length `[float]`
          
          Before converting to discrete GC percentages, we apply smoothing to the raw GC counts separately for each fragment length. For a fragment length of 150, we thus have counts of fragments with GCs ranging from 0..=150, and smoothing happens on this scale so the distance between elements are the same for all fragment lengths.
          
          Note: The same smoothing parameters (sigma and radius) are used for downstream `cfdna gc-bias` calls.
          
          [default: 0.55]

      --smoothing-radius <SMOOTHING_RADIUS>
          Radius of Gaussian kernel that smoothes raw GC counts for each fragment length `[integer]`
          
          Kernel size is `2 * radius + 1`.
          
          [default: 2]

      --skip-smoothing
          Whether to skip the smoothing of raw GC counts `[flag]`
```
