# Extract Fragment Lengths

Multiple studies have used fragment lengths (count distributions) to detect cancer [REFS].

## Base command

```bash
cfdna lengths --help

cfdna lengths \
  --bam <sample>.bam \
  --output-dir <sample_directory>/lengths \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --by-size 1000000
```

### Use length bins

...

The default length bins are specified as `--length-bins 30:1001:1`, which gives us the full single-bp resolution fragment length distribution between 30bp and 1000bp.

```bash
cfdna lengths \
  ... \
  # Bin every 10bp from 30-500bp
  --length-bins 30:500:10
  # OR Specify edges directly (last end is exclusive)
  --length-bins 100 151 221

```

## GC-bias correction example

```bash
cfdna lengths \
  --bam <sample>.bam \
  --output-dir <sample_directory>/lengths \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --by-size 1000000 \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit
```

## Genomic smoothing example

```bash
cfdna lengths \
  --bam <sample>.bam \
  --output-dir <sample_directory>/lengths \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --by-size 1000000 \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```

## GC-bias correction + genomic smoothing

```bash
cfdna lengths \
  --bam <sample>.bam \
  --output-dir <sample_directory>/lengths \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --by-size 1000000 \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```

## Adjusting for indels

The default `lengths` behavior is to use the fragment span on the reference genome and ignore whether insertions or deletions (InDels) are present in the reads.

By specifying `--indel-mode adjust`, the fragment lengths are adjusted for indels, which should be closer to the size of the original DNA molecule.

```bash
cfdna lengths \
  ... \
  --indel-mode adjust
```

This is an analysis choice, not a requirement. If you are unsure, start without it.
