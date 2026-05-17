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
  --gc-file <sample_directory>/gc_bias/<sample_id>.gc_bias_correction.zarr \
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
  --scaling-factors <sample_directory>/scaling_factors/<sample_id>.fragment_counts.scaling_factors.tsv
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
  --gc-file <sample_directory>/gc_bias/<sample_id>.gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors <sample_directory>/scaling_factors/<sample_id>.fragment_counts.scaling_factors.tsv
```

## Downstream usage

`cfdna ends` writes `<sample_id>.end_motifs.zarr`. The helper packages give a smaller user-facing API than working with the Zarr arrays directly.

In Python:

```python
import cfdnalab as cfl

ends = cfl.read_end_motifs("<sample_id>.end_motifs.zarr")

print(ends.storage_mode())
print(ends.row_mode())
print(ends.motif_metadata().head())

motif = ends.motifs()[0]
if ends.storage_mode() == "sparse_coo":
    nonzero_counts = ends.sparse_coo_data_frame()
    motif_counts = nonzero_counts[nonzero_counts["motif"] == motif]
else:
    motif_counts = ends.dense_data_frame_for_motif(motif)
```

In R:

```r
library(cfdnalab)

ends <- read_end_motifs("<sample_id>.end_motifs.zarr")

storage_mode(ends)
row_mode(ends)
head(motifs(ends))

motif <- motifs(ends)$motif[[1]]
if (storage_mode(ends) == "sparse_coo") {
  motif_counts <- sparse_data_frame_for_motif(ends, motif)
} else {
  motif_counts <- dense_data_frame_for_motif(ends, motif)
}
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
