# Extract Fragment Coverage

Fragment coverage measures how many fragments overlap each genomic position. In contrast to many non-cfDNA-tools, we (optionally) count the gap between paired reads along with the aligned bases of the reads. We avoid double counting when reads overlap.

When no GC correction or genomic smoothing is applied, each fragment is counted as `1` in the overlapping (aligned / gap) positions. Using GC correction and/or genomic smoothing changes this to a weight (floating point).

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
  --by-size 1000000 \
  --per-window 'average'
```

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
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```
