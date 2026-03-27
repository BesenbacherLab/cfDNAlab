# Correct GC-Bias

Fragmentomics features are vulnerable to biases from various sample-handling and sequencing processes, such as PCR amplification. `cfDNAlab` commands thus allow the correction of the commonly observed GC-bias.

This requires only a few steps.

## Step 1. Build reference GC bias once per assembly

Calculate the expected GC bias in the reference genome assembly (for example hg38). This output can be reused for all samples aligned to that assembly.

```bash
cfdna ref-gc-bias --help

# Run once per assembly
cfdna ref-gc-bias \
  --ref-2bit <path>/hg38.2bit \
  --output-dir <ref_gc_directory> \
  --output-prefix hg38 \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed
```

## Step 2. Build sample-specific GC correction

```bash
cfdna gc-bias --help

cfdna gc-bias \
  --bam <sample>.bam \
  --output-dir <sample_directory>/gc_bias \
  --n-threads 12 \
  --ref-2bit <path>/hg38.2bit \
  --ref-gc-file <ref_gc_directory>/hg38.ref_gc_package.npz \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed
```

Use the same blacklist inputs as in step 1.

## Step 3. Apply correction in feature extraction commands

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  ... \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit
```

The same pattern works for `lengths` and `midpoints`.

## Alternative: read GC weights from BAM aux tags

If you prefer a different or custom GC-bias tool, the feature extraction commands also accept reading a GC weight from a BAM aux tag.

```bash
cfdna fcoverage \
  --bam <sample>.bam \
  ... \
  --gc-tag 'GC'
```
