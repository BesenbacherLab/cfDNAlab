# Genomic Smoothing

For some commands, like `cfdna midpoints`, you may want all genomic regions to have approximately the same contribution to the features, for example to reduce the effect of copy number alterations.

The workflow uses coverage-based scaling factors from `cfdna coverage-weights`. In downstream feature extraction, each fragment contribution is multiplied by a scaling factor for its genomic bin.

## Step 1. Build per-sample scaling factors

```bash
cfdna coverage-weights --help

cfdna coverage-weights \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage_weights \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed
```

Tip: use the same blacklist set you use in your other per-sample commands.

## Step 2. Apply smoothing in feature extraction commands

```bash
cfdna midpoints \
  --bam <sample>.bam \
  ... \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```

The same `--scaling-factors` input pattern works for `fcoverage` and `lengths`.

## Optional: combine with GC-bias correction

Genomic smoothing and GC-bias correction can be used together.

```bash
cfdna midpoints \
  --bam <sample>.bam \
  ... \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```
