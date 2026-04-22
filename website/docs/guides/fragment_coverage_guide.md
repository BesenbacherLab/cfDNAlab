# Extract Fragment Coverage

Fragment coverage measures how many fragments overlap each genomic position. In contrast to many non-cfDNA-tools, we (optionally) count the gap between paired reads along with the aligned bases of the reads. We avoid double counting when reads overlap.

When no GC correction or genomic smoothing is applied, each fragment is counted as `1.0` in the overlapping (aligned / gap) positions. Using GC correction and/or genomic smoothing changes this to a weight (floating point).

## Base command

```bash
cfdna fcoverage --help

cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  # Get average coverage in 1Mb windows (remove for positional output)
  --by-size 1000000 \
  --per-window 'average'
```

## Normalize by countable bases

In the default coverage mode, longer fragments weight more than shorter fragments simply because they cover more positions (each counted with `1.0` before correction/scaling). When you instead want all fragments to contribute the same mass, you can enable `--normalize-by-length` and count them as `1.0/num_countable_bases` (the fragment length in most cases). 

On a positional-level, this might be a bit funky to interpret, but once you're working with larger windows, this represents **fragment counts**, where fragments not fully contained by a window are weighted by their overlaps.

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --by-size 1000000 \
  # Normalize by countable bases
  --normalize-by-length \
  # "total" gives us the fragment count scale in larger windows
  --per-window 'total' \
  # Additional precision in the output
  --decimals 3
```

If you want positional values on the same scale as the regular coverage tracks, you can use `--normalize-by-length=restore-mean`. This multiplies each coverage value by the mean `num_countable_bases` from all the counted fragments. This makes interpretation of positional values *closer* to "how many fragments overlap this position". For windows, this pairs better with `--per-window average` (**NOTE**: rescaling can lead to longer runtime, so for "fragment counts per window", the `--normalize-by-length --per-window 'total'` path is better).

## GC-bias correction example

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --by-size 1000000 \
  --per-window 'average' \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit
```

## Genomic smoothing example

When you're interested in local, relative coverage changes instead of large-scale changes from CNVs etc., you can use genomic smoothing to weight contributions more similarly across genomic regions.

When calling `cfdna fcoverage` **WITHOUT** `--normalize-by-length` (this is the default), we suggest using the `cfdna coverage-weights` command for calculating scaling factors, as it has the same fragment-length-weighting.

When calling **WITH** `--normalize-by-length`, we suggest using `cfdna fragment-count-weights` instead.

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --by-size 1000000 \
  --per-window 'average' \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```

## GC-bias correction + genomic smoothing

**NOTE**: *Requires* using GC-bias correction when calculating the scaling factors to avoid double-correction.

See "Genomic smoothing example" for when to use `cfdna coverage-weights` or `cfdna fragment-count-weights`.

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --by-size 1000000 \
  --per-window 'average' \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors <sample_directory>/gc_corrected_coverage_weights/<sample_id>.scaling_factors.tsv
```
