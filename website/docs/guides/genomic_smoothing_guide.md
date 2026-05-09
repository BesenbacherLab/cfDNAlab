# Genomic Smoothing

For some commands, like `cfdna midpoints`, you may want all genomic regions to contribute approximately the same "mass" to the features, for example to reduce the effect of copy number alterations.

`cfdna` comes with two scaling weight calculation commands: 

 - `cfdna fragment-count-weights`
 - `cfdna coverage-weights`

**Conceptually**, the main difference is that coverage weights longer fragments higher than shorter fragments just due to them covering more positions. When that is not specifically beneficial, we recommend using `cfdna fragment-count-weights` where all fragments contribute the same mass independently of length. 

**Technically**, the only difference in their calculations is that where `cfdna coverage-weights` counts each fragment as `1.0` in all covered positions (before the optional GC-bias correction). `cfdna fragment-count-weights` instead counts `1.0 / num_countable_bases` (usually just the fragment length) in the covered positions. When counting in large windows, this approximates fragment counts very closely (the few fragments not fully contained in their windows are weighted by their overlap).

When applied in downstream feature extractions, each fragment contribution is multiplied by the scaling factor of its genomic scaling bin.

## Step 1. Build per-sample scaling factors

```bash
cfdna fragment-count-weights --help
cfdna coverage-weights --help

cfdna fragment-count-weights \
  --bam <sample>.bam \
  --output-dir <sample_directory>/count_weights \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed

cfdna coverage-weights \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage_weights \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed
```

**Tip**: use the same blacklist set you use in your other per-sample commands.

## Step 2. Apply smoothing in feature extraction commands

For midpoints, where the signal is specifically fragment counts, not a length-weighted coverage, we prefer the fragment count-based weights:

```bash
cfdna midpoints \
  --bam <sample>.bam \
  ... \
  --scaling-factors <sample_directory>/count_weights/<sample_id>.scaling_factors.tsv
```

The same `--scaling-factors` input pattern works for `fcoverage` and `lengths`.

## Combine with GC-bias correction

Genomic smoothing and GC-bias correction can be used together. 

The only *requirement* is that they are used together consistently, meaning that the GC-bias correction must be used **when calculating the scaling factors**:

```bash

cfdna fragment-count-weights \
  --bam <sample>.bam \
  --output-dir <sample_directory>/gc_corrected_count_weights \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit

cfdna midpoints \
  --bam <sample>.bam \
  ... \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors <sample_directory>/gc_corrected_count_weights/<sample_id>.scaling_factors.tsv
```
