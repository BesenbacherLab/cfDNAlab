# Extract Fragment Coverage

Fragment coverage measures how many fragments overlap each genomic position. In contrast to many non-cfDNA-tools, we (optionally) count the gap between paired reads along with the aligned bases of the reads. We avoid double counting when reads overlap.

When no GC correction or genomic smoothing is applied, each fragment is counted as `1.0` in the overlapping (aligned / gap) positions. Using GC correction and/or genomic smoothing changes this to a weight (floating point).

For a fragment *counts*-like signal, see ["Normalize by countable bases"](#normalize-by-countable-bases).

## Examples

The following examples show different aspects of the `cfdna fcoverage` command. They can of course be combined in a multitude of ways, but for simplification we just show one aspect at a time.

### Positional fragment coverage

Extract the per-position fragment coverage into BedGraph files: 

```bash
cfdna fcoverage --help

cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
```

BedGraph files are commonly converted to BigWig files to allow indexed lookups in downstream analyses. E.g., for the UCSC hg38 chromosome:

```bash
# Install UCSC tools
conda install -c bioconda ucsc-bedgraphtobigwig ucsc-twobitinfo

# Get the chromosome sizes from the 2bit reference genome (used in --ref-2bit)
twoBitInfo hg38.2bit hg38.chrom.sizes

# Decompress the bedgraph
zstd -d <prefix>.fcoverage.per_position.bedgraph.zst

# Convert bedGraph to bigWig
bedGraphToBigWig \
  <prefix>.fcoverage.per_position.bedgraph \
  hg38.chrom.sizes \
  <prefix>.fcoverage.per_position.bw
```

### Average coverage in windows

Extract the average positional coverage in 1Mb bins:

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Average coverage in 1Mb windows
  --by-size 1000000 \
  --per-window 'average'
```

### GC-bias correction example

To correct the sample-specific GC-bias, you need to precompute the correction matrix (see the [GC-bias guide](./correct_gc_bias_guide.md)). Then you provide that correction file as:

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Precomputed GC correction file and reference genome
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit
```

### Genomic smoothing example

When you're interested in local, relative coverage changes instead of large-scale changes from CNVs etc., you can use genomic smoothing to weight contributions more similarly across genomic regions. In some analyses, this can be thought of as a copy-number normalization. See more in the [genomic smoothing guide](./genomic_smoothing_guide.md).

When calling `cfdna fcoverage` **WITHOUT** `--normalize-by-length` (this is the default), we suggest using the `cfdna coverage-weights` command for calculating scaling factors, as it has the same fragment-length-weighting.

When calling **WITH** `--normalize-by-length`, we suggest using `cfdna fragment-count-weights` instead.

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Precomputed scaling factors
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.coverage.scaling_factors.tsv
```

### GC-bias correction + genomic smoothing

**NOTE**: *Requires* using GC-bias correction when calculating the scaling factors to avoid double-correction.

See ["Genomic smoothing example"](#genomic-smoothing-example) for when to use `cfdna coverage-weights` or `cfdna fragment-count-weights`.

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Precomputed GC correction file and reference genome
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit \
  # Precomputed scaling factors
  --scaling-factors <sample_directory>/gc_corrected_coverage_weights/<sample_id>.coverage.scaling_factors.tsv
```

### Normalize by countable bases

In the default coverage mode, longer fragments weight more than shorter fragments simply because they cover more positions (each counted with `1.0` before correction/scaling). When you instead want all fragments to contribute the same mass, you can enable `--normalize-by-length` and count them as `1.0/num_countable_bases` (the fragment length in most cases). 

On a positional-level, this might be a bit funky to interpret, but once you're working with larger windows, this represents **fragment counts**, where fragments not fully contained by a window are weighted by their overlap fraction.

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Get aggregated "counts" in 1Mb windows
  --by-size 1000000 \
  # Normalize by countable bases
  --normalize-by-length \
  # "total" gives us the fragment count scale in larger windows
  --per-window 'total' \
  # Additional precision in the output
  --decimals 3
```

If you want positional values on the same scale as the regular coverage tracks, you can use `--normalize-by-length=restore-mean`. This multiplies each coverage value by the mean `num_countable_bases` from all the counted fragments. This makes interpretation of positional values *closer* to "how many fragments overlap this position". For windows, this pairs better with `--per-window average` (**NOTE**: rescaling can lead to longer runtime, so for "fragment counts per window", the `--normalize-by-length --per-window 'total'` approach is better).

### Windowed summary statistics

So, you're not satisfied with a mere average coverage metric per window, ay? Me neither, at least not in all cases. Hence, I implemented the option to extract the following summary statistics per window (or per group when supplying grouped windows with `--by-grouped-bed`):

Position stats:

- `span_positions`: Number of reference positions spanned by the window.

- `blacklisted_positions`: Number of positions removed by blacklist masking.

- `eligible_positions`: Number of positions that can contribute to coverage after masking.

- `nonzero_positions`: Number of eligible positions with coverage above zero.

- `covered_fraction`: Fraction of eligible positions with coverage above zero.

Together, these columns describe the denominator and coverage breadth used by the coverage statistics: `span_positions` is the full reference span, `blacklisted_positions` is the masked part of that span, `eligible_positions` is the unmasked part used for coverage calculations, and `nonzero_positions` plus `covered_fraction` describe the subset and fraction of eligible positions that actually have coverage.

Coverage stats:

- `total_coverage`: Sum of coverage over eligible positions.

- `total_squared_coverage`: Sum of squared coverage over eligible positions.

- `average_coverage`: Mean coverage across eligible positions.

- `variance_coverage`: Coverage variance across eligible positions.

- `sd_coverage`: Coverage standard deviation across eligible positions.

- `coefficient_of_variation_coverage`: `sd_coverage / average_coverage`.

When `--normalize-by-length` is active, the same aggregate columns are named with `fragment_mass` instead of `coverage`.

Calculate these summary statistics for groups of windows:

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Summary statistics per group
  --by-grouped-bed <path>/grouped_windows.bed \
  --per-window 'summary-stats'
```