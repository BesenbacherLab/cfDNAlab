# Extract Fragment End Motifs

Multiple studies have used fragment end- and breakpoint-motifs to study cfDNA fragmentation biology. These motif frequencies can capture sequence preferences around where fragments start and end.

With `cfdna ends`, we can count the bases just inside (`--k-inside`) and just outside (`--k-outside`) the fragment ends. The outside bases are taken from the reference genome, while the inside bases can be taken from either the read sequence (default) or the reference genome.

The following shows the counting for aligned fragment ends in paired-end sequencing data:

For `--k-inside 2 --k-outside 2`:

```text
Reference 5' >>>>>>>>>>>>>>>>>>>> 3'
             ATCGTTTTTTTTTTTTCATC
Fragment     --|--------------|--
Forward     5' |>>>>>>>| 3'
  Outside    AT
  Inside       CG
Reverse            3' |<<<<<<<| 5'
  Inside                     CA
  Outside                      TC
```

For Nanopore-like single-end data, where each read represents a fragment, the aligned 5' and 3' ends are used.

`cfdna ends` writes motifs as `"<outside>_<inside>"`. To achieve this, the right-side motifs are reverse complemented together. For the above example, we would get the count: `AT_CG: 1, GA_TG: 1`.

## Handling clipped ends

When the ends of fragments have been soft-clipped, it indicates that something couldn't be aligned properly. `cfdna ends` defaults to *not* counting motifs in soft-clipped ends. This reduces the risk of counting motifs that do not correspond to the actual ends of the DNA molecules. 

If you want to keep using the aligned fragment boundaries, you can switch to `--clip-strategy aligned`.

When you are only counting the inside bases, you can also choose to *include* the soft-clipped bases with either `--clip-strategy include-at-aligned-boundary` or `--clip-strategy include-at-shifted-boundary`. The "boundary" part of the option names refers to where the software should consider the fragment ends to lie on the reference genome. Since the clipped part could technically stem from a different region of the reference genome, this choice is not trivial. If you want to include the clipped bases, we recommend you read the `cfdna ends` [documentation](https://cfdnalab.tools/docs/generated/cli/ends).

## Sparse or dense outputs

**Dense**: When you work with small combined `--k-inside --k-outside` settings, you can set `--all-motifs` to ensure all possible motifs are in the file even when they are not observed in the data. This makes the output file dense, meaning that we also store the zero counts.

**Sparse**: For larger `k`-settings, we recommend keeping the output sparse (default). This allows us to only store the nonzero counts. If some motifs are not observed at all though, they will not be in the list of motif labels. This should be handled downstream when performing cross-sample analyses where the sets of motifs might differ.

## Examples

The following examples show different aspects of the `cfdna ends` command. They can of course be combined in a multitude of ways, but for simplification we just show one aspect at a time.

### Base command

The following example counts 2 bases outside plus 2 bases inside the fragment:

```bash
cfdna ends --help

cfdna ends \
  --bam <sample>.bam \
  --output-dir <sample_directory>/ends \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --ref-2bit <path>/hg38.2bit \
  --k-inside 2 \
  --k-outside 2
```

### GC-bias correction example

If you haven't already, please start by [computing the GC correction matrix](https://cfdnalab.tools/docs/guides/correct_gc_bias_guide). Then pass in the sample-specific GC-bias correction package and the reference genome:

```bash
cfdna ends \
  --bam <sample>.bam \
  --output-dir <sample_directory>/ends \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --ref-2bit <path>/hg38.2bit \
  --k-inside 2 \
  --k-outside 2 \
  --gc-file <sample_directory>/gc_bias/<sample_id>.gc_bias_correction.zarr
```

### Genomic smoothing example

If you want each region of the genome to contribute approximately the same to the motif counts, you can pass in the [precomputed scaling factors](https://cfdnalab.tools/docs/guides/genomic_smoothing_guide). Since end motifs are fragment-count features, we recommend computing the scaling factors with `cfdna fragment-count-weights`.

```bash
cfdna ends \
  --bam <sample>.bam \
  --output-dir <sample_directory>/ends \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --ref-2bit <path>/hg38.2bit \
  --k-inside 2 \
  --k-outside 2 \
  --scaling-factors <sample_directory>/scaling_factors/<sample_id>.fragment_counts.scaling_factors.tsv
```

Note: You can easily use GC correction and genomic smoothing together. Just ensure you enable GC correction when calculating the scaling factors.

### Filter on base qualities

Depending on the question you are asking, it can make sense to filter your counts based on the base qualities in the **inside** read bases of the motifs. `--bq-filter` lets you filter either the whole fragment or individual ends with a syntax like:

 - `--bq-filter "min in end >= 30"` (for "all bases have decent quality")

 - `--bq-filter "mean in fragment >= 30"` (for "average bases have decent quality")

 - `--bq-filter "max in fragment < 20"` (for "no bases have decent quality")

Repeat `--bq-filter` to count only ends that pass all **end filters** and belong to fragments that pass all **fragment filters**.

For the full set of options, see the `cfdna ends` [documentation](https://cfdnalab.tools/docs/generated/cli/ends).

```bash
cfdna ends \
  --bam <sample>.bam \
  --output-dir <sample_directory>/ends \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --ref-2bit <path>/hg38.2bit \
  --k-inside 2 \
  --k-outside 2 \
  --bq-filter "min in end >= 30"
```

## Load motif counts in Python

Our [`cfdnalab` **Python** package](https://cfdnalab.tools/docs/generated/loaders/python-api) can help you load the `<sample_id>.end_motifs.zarr` output and extract data frames, arrays and metadata.

The methods differ when the data was built with `--by-size` / `--by-bed`, `by-grouped-bed` or the default "global" window. The below examples are a small subset of the options for subsetting the counts. See more in the Python package [documentation](https://cfdnalab.tools/docs/generated/loaders/python-api).

```python
import cfdnalab as cfl

ends = cfl.read_end_motifs("<sample_id>.end_motifs.zarr")

print(ends.storage_mode())
print(ends.row_mode())
motifs = ends.motifs_metadata()
print(motifs.head())

# Get the first motif
motif = motifs.loc[0, "motif"]

# Extract counts as a long-format data frame
all_counts = ends.data_frame()
motif_counts = ends.data_frame(motifs=motif)

# Extract the counts in a sparse matrix or dense array
# depending on how the data was stored
if ends.storage_mode() == "sparse_coo":
  motif_count_matrix = ends.sparse_counts_matrix(motifs=motif)
else:
  motif_count_array = ends.dense_counts_array(motifs=motif)
```

## Load motif counts in R

Our [`cfdnalab` **R** package](https://cfdnalab.tools/docs/generated/loaders/r-api) can help you load the `<sample_id>.end_motifs.zarr` output and extract data frames, arrays and metadata.

The methods differ when the data was built with `by-size/by-bed`, `by-grouped-bed` or the default "global" window. The below examples are a small subset of the options for subsetting the counts. See more in the R package [documentation](https://cfdnalab.tools/docs/generated/loaders/r-api).

```r
library(cfdnalab)

ends <- read_end_motifs("<sample_id>.end_motifs.zarr")

storage_mode(ends)
row_mode(ends)
head(motifs(ends))

# Get the first motif
motif <- motifs(ends)$motif[[1]]

# Extract counts as a long-format data frame
all_counts <- end_motif_data_frame(ends)
motif_counts <- end_motif_data_frame(ends, motifs = motif)

# Extract the counts in a sparse/dense matrix 
# depending on how the data was stored
if (storage_mode(ends) == "sparse_coo") {
  counts <- sparse_counts_matrix(ends, motifs=motif)
} else {
  counts <- dense_counts_matrix(ends, motifs=motif)
}
```

## Correct reference k-mers bias

The end-motif frequencies reflect the reference sequence composition in the regions you count. When a region is high in GC, you're more likely to observe end motifs with GC in them. To make regions more comparable in enrichment analyses, you can remove (to some degree) this reference bias from the end-motif frequencies.

Reference correction divides each observed end-motif count by a reference-based correction factor for the matched row. Motifs that are common in the reference row are scaled down. Motifs that are rare in the reference row are scaled up.

First precompute matching reference k-mer frequencies with `cfdna ref-kmers`, then load the sample end motifs and reference k-mers in Python or R and pass the reference k-mers to the usual end-motif output loaders.

You only need to run `cfdna ref-kmers` once per reference setup. That means once per reference genome, blacklist, k-mer size, etc.

### Two-sided correction

When both `--k-outside` and `--k-inside` are specified, choose a `two_sided_correction` option:

- `"joint"` keeps full labels such as `AC_TG` and corrects each count with the matching reference k-mer, `ACTG`.
- `"split"` keeps full labels such as `AC_TG`, but calculates the correction factor from the two sides separately. For `AC_TG`, separate correction factors are calculated for outside label `AC` and inside label `TG`. Those two correction factors are multiplied and applied to the observed `AC_TG` count. Use this when the full two-sided motif labels should remain in the result, but the exact full reference k-mers are too sparse or the correction should treat outside and inside sequence composition separately.
- `"outside"` returns outside labels such as `AC_`. For each outside label, all full motif counts with that outside label are summed first. For example, `AC_AA` and `AC_TG` both contribute to the `AC_` count. That summed count is corrected using the outside label `AC`.
- `"inside"` returns inside labels such as `_TG`. For each inside label, all full motif counts with that inside label are summed first. For example, `AA_TG` and `AC_TG` both contribute to the `_TG` count. That summed count is corrected using the inside label `TG`.

For `"split"`, `"outside"`, and `"inside"`, side-specific reference frequencies are calculated from the loaded full-length reference k-mers. For example, the outside frequency for `AC` is the sum of frequencies for loaded k-mers with prefix `AC`, such as `ACTG` and `ACAA`. The inside frequency for `TG` is the corresponding sum over loaded k-mers with suffix `TG`. Separate shorter `ref-kmers` runs are not supported.

A motifs file used for `ref-kmers` restricts these sums to the k-mers in that file. Without a motifs file, all k-mers in the reference output can contribute, including k-mers absent from the sample end-motif output.

One-sided outputs do not use the `two_sided_correction` argument.

### Extracting reference k-mers

Use `--kmer-size` equal to `--k-outside + --k-inside`. For `--k-outside 2 --k-inside 2`, use `--kmer-size 4`:

```bash
cfdna ref-kmers \
  --ref-2bit <path>/hg38.2bit \
  --output-dir <reference_directory>/ref_kmers \
  --output-prefix hg38.k4 \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --kmer-size 4 \
  --all-motifs
```

If the `cfdna ends` call used `--by-size`, `--by-bed`, or `--by-grouped-bed`, pass the same windowing or grouping option to `cfdna ref-kmers`. Use the same motifs file for both commands when you want correction to use exactly that motif subset. Omit the motifs file from `ref-kmers` when the complete reference motif set should contribute to correction. `ref-kmers` accepts end-motif labels such as `AT_CG` and treats them as the reference k-mer `ATCG`.

Both outputs use **forward-oriented** motif labels only. `cfdna ends` reverse-complements right-end motifs before writing them, and `cfdna ref-kmers` counts reference k-mers left-to-right. Ungrouped motif labels are matched by removing the underscore, for example `AT_CG -> ATCG`, `_AA -> AA`, and `AA_ -> AA`. Grouped motifs-file output are matched directly on the group names.

### Apply correction in Python

In Python:

```python
import cfdnalab as cfl

ends = cfl.read_end_motifs("<sample_id>.end_motifs.zarr")
ref_kmers = cfl.read_ref_kmers("<reference>.ref_kmers.zarr")

corrected = ends.data_frame(ref_kmers=ref_kmers, two_sided_correction="joint")

# Drop or NA motifs that are only present in the sample (can happen with large k-mers)
corrected_without_unsupported = ends.data_frame(
  ref_kmers=ref_kmers,
  two_sided_correction="split",
  unsupported_motifs="drop",
)

# Matrix columns can differ from the stored motif axis in outside/inside mode.
corrected_matrix = ends.corrected_counts_array(
  ref_kmers,
  two_sided_correction="outside",
)
corrected_motifs = ends.corrected_motifs_metadata(
  ref_kmers,
  two_sided_correction="outside",
)
# corrected_motifs row 0 describes corrected_matrix column 0, and so on.

```

### Apply correction in R

In R:

```r
library(cfdnalab)

ends <- read_end_motifs("<sample_id>.end_motifs.zarr")
ref_kmers <- read_ref_kmers("<reference>.ref_kmers.zarr")

corrected <- end_motif_data_frame(
  ends,
  ref_kmers = ref_kmers,
  two_sided_correction = "joint"
)

# Drop or NA motifs that are only present in the sample (can happen with large k-mers)
corrected_without_unsupported <- end_motif_data_frame(
  ends,
  ref_kmers = ref_kmers,
  two_sided_correction = "split",
  unsupported_motifs = "drop"
)

# Corrected R matrices carry their corrected motif axis as column names.
corrected_matrix <- dense_corrected_counts_matrix(
  ends,
  ref_kmers,
  two_sided_correction = "outside"
)
corrected_motifs <- colnames(corrected_matrix)

```

### Corrected output

Corrected data frames include `corrected_count` and `corrected_frequency`. For each output row, `corrected_frequency` is the corrected count divided by the sum of corrected counts over the full motif axis produced by the correction mode. Motif selection is applied afterward and does not renormalize the frequencies, so a selected subset may sum to less than 1. If the corrected total is zero, finite frequencies are zero.

For Python corrected arrays and sparse matrices, call `corrected_motifs_metadata()` with the same reference, correction mode, and motif selection. It returns `matrix_column`, `motif_index`, and `motif`, in matrix-column order. This is especially important for `"outside"` and `"inside"`, which replace the stored full-motif axis with a shorter side-motif axis. R corrected matrices expose the same corrected axis through `colnames()`.

The loaders handle dense and sparse outputs without manually joining the data frames. By default they stop if an observed sample motif has no positive correction factor under the selected mode. Use `unsupported_motifs = "drop"` to omit unsupported rows before frequencies are calculated, or `unsupported_motifs = "keep_na"` to keep them with undefined corrected counts. With `keep_na`, one undefined positive corrected count makes every corrected frequency in that output row undefined because the full corrected total is unknown.
