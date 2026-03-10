# Fragment Midpoint Recipes

Multiple studies have profiled midpoint coverage around e.g. transcription factor binding sites (summed per transcription factor, per position) [REFS]. This can inform about the binding activity of different transcription factors related to cancer.

## Base command

```bash
cfdna midpoints --help

cfdna midpoints \
  --bam <sample>.bam \
  --output-dir <sample_directory>/midpoints \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --intervals <fixed_size_intervals>.tsv \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --length-bins {30..1000..10}
```

## GC-bias correction example

```bash
cfdna midpoints \
  --bam <sample>.bam \
  --output-dir <sample_directory>/midpoints \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --intervals <fixed_size_intervals>.tsv \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --length-bins {30..1000..10} \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit
```

## Genomic smoothing example

```bash
cfdna midpoints \
  --bam <sample>.bam \
  --output-dir <sample_directory>/midpoints \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --intervals <fixed_size_intervals>.tsv \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --length-bins {30..1000..10} \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```

## GC-bias correction + genomic smoothing

```bash
cfdna midpoints \
  --bam <sample>.bam \
  --output-dir <sample_directory>/midpoints \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --intervals <fixed_size_intervals>.tsv \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --length-bins {30..1000..10} \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```

The intervals must have the same fixed size. The expected columns are: `chromosome, start, end, group_name` (where `group_name` is the group to collapse profiles by, e.g., the transcription factor ID). The intervals should be sorted by chromosome and start coordinates.
