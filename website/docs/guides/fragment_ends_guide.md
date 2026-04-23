# Extract Fragment End Motifs

Multiple studies have used fragment end- and breakpoint-motifs to study cfDNA fragmentation biology [REFS]. These motif frequencies can capture sequence preferences around where fragments start and end.

## Base command

```bash
cfdna ends --help

cfdna ends \
  --bam <sample>.bam \
  --output-dir <sample_directory>/ends \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --by-size 1000000 \
  --k-inside 2 \
  --k-outside 2
```

## GC-bias correction example

```bash
cfdna ends \
  --bam <sample>.bam \
  --output-dir <sample_directory>/ends \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --by-size 1000000 \
  --k-inside 2 \
  --k-outside 2 \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit
```

## Genomic smoothing example

```bash
cfdna ends \
  --bam <sample>.bam \
  --output-dir <sample_directory>/ends \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --by-size 1000000 \
  --k-inside 2 \
  --k-outside 2 \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```

## GC-bias correction + genomic smoothing

```bash
cfdna ends \
  --bam <sample>.bam \
  --output-dir <sample_directory>/ends \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \
  --by-size 1000000 \
  --k-inside 2 \
  --k-outside 2 \
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```

## Handling clipped ends

The default `ends` behavior is conservative around soft clipping. With `--clip-strategy skip`, motifs are discarded when the relevant fragment end is soft-clipped.

If you want to keep using the aligned fragment boundaries, you can switch to `--clip-strategy aligned`.

```bash
cfdna ends \
  ... \
  --clip-strategy aligned
```

The `raw-aligned-boundary` and `raw-shifted-boundary` modes are stronger analysis choices. Use them only when you specifically want raw read bases, including soft-clipped sequence, to contribute to the motif.
