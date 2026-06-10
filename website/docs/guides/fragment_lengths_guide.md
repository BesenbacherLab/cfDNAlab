# Extract Fragment Lengths

Multiple studies have used fragment lengths (count distributions) to detect cancer (Renaud et al. 2022, Cristiano et al. 2019).

## Examples

The following examples show different aspects of the `cfdna lengths` command. They can of course be combined in a multitude of ways, but for simplification we just show one aspect at a time.

### Global length distribution

Count the global fragment length distribution:

```bash
cfdna lengths --help

cfdna lengths \
  --bam <sample>.bam \
  --output-dir <sample_directory>/lengths \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed
```

### Specify length bins

By default, the command counts each fragment length from 30-1000bp. To change those limits or to count in bins of fragment lengths, you can specify the `--length-bins` as `start:stop:step` or as space-separated edges like `--length-bins 100 151 221`. The default is `--length-bins 30:1001:1`.

```bash
cfdna lengths \
  ... \
  # Bin every 10bp from 30-500bp
  --length-bins 30:500:10
  # OR Specify edges directly (last end is exclusive)
  --length-bins 100 151 221

```

### GC-bias correction example

To correct the sample-specific GC-bias, you need to precompute the correction matrix (see the [GC-bias guide](./correct_gc_bias_guide.md)). Then you provide that correction file as:

```bash
cfdna lengths \
  --bam <sample>.bam \
  --output-dir <sample_directory>/lengths \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Precomputed GC correction file and reference genome
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit
```

### Genomic smoothing example

When you're interested in local, relative coverage changes instead of large-scale changes from CNVs etc., you can use genomic smoothing to weight contributions more similarly across genomic regions. In some analyses, this can be thought of as a copy-number normalization.

Since we're counting fragments per lengths, we suggest using the `cfdna fragment-count-weights` command for calculating scaling factors (see the [genomic smoothing guide](./genomic_smoothing_guide.md)).

```bash
cfdna lengths \
  --bam <sample>.bam \
  --output-dir <sample_directory>/lengths \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Precomputed scaling factors
  --scaling-factors <sample_directory>/count_weights/<sample_id>.fragment_counts.scaling_factors.tsv
```

### GC-bias correction + genomic smoothing

**NOTE**: *Requires* using GC-bias correction when calculating the scaling factors to avoid double-correction.

```bash
cfdna lengths \
  --bam <sample>.bam \
  --output-dir <sample_directory>/lengths \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Precomputed GC correction file and reference genome
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit \
  # Precomputed scaling factors
  --scaling-factors <sample_directory>/gc_corrected_count_weights/<sample_id>.fragment_counts.scaling_factors.tsv
```

### Adjusting for indels

The default `lengths` behavior is to use the fragment span on the reference genome and ignore whether insertions or deletions (InDels) are present in the reads.

By specifying `--indel-mode adjust`, the fragment lengths are adjusted for indels, which should be closer to the size of the original DNA molecule.

```bash
cfdna lengths \
  ... \
  --indel-mode adjust
```

This is an analysis choice, not a requirement. If you are unsure, start without it.

## Load length counts in Python

Our [`cfdnalab` **Python** package](https://cfdnalab.tools/docs/generated/loaders/python-api) can help you load the `<sample_id>.length_counts.tsv.zst` output and extract data frames, arrays and metadata.

The methods differ when the data was built with `--by-size` / `--by-bed`, `by-grouped-bed` or the default "global" window. The below examples are a small subset of the options for subsetting the counts. See more in the Python package [documentation](https://cfdnalab.tools/docs/generated/loaders/python-api).

```python
import cfdnalab as cfl

lengths = cfl.read_lengths("sample.length_counts.tsv.zst")

# Extract metadata (depends on applied windowing and grouping)
length_bins = lengths.length_bins()
groups = lengths.group_metadata()
windows = lengths.window_metadata()

# Extract long data frame with all data
length_data = lengths.data_frame()

# Extract length bins containing a given length range
length_range_data = lengths.data_frame(with_length_range=(100, 220))

# Normalize to fractions across all length bins
fraction_data = lengths.data_frame(value="fraction")

# Normalize to fractions of only the selected bins
selected_fraction_data = lengths.data_frame(
    with_length_range=(100, 220),
    value="fraction",
    denominator="selected_bins",
)

# For windowed counts, filter highly blacklisted windows/groups
filtered = lengths.data_frame(max_blacklisted_fraction=0.3)

# Select counts for specific group(s) or window(s)
group_data = lengths.data_frame(groups="t-cells", value="fraction")
window_data = lengths.data_frame(window_idxs=0, value="fraction")

# Extract counts as an array
counts = lengths.counts_array()
```

## Load length counts in R

Our [`cfdnalab` **R** package](https://cfdnalab.tools/docs/generated/loaders/r-api) can help you load the wide `<sample_id>.length_counts.tsv.zst` output and extract data frames, arrays, and metadata.

The methods differ when the data was built with `by-size/by-bed`, `by-grouped-bed` or the default "global" window. The below examples are a small subset of the options for subsetting the counts. See more in the R package [documentation](https://cfdnalab.tools/docs/generated/loaders/r-api).

```r
library(cfdnalab)

lengths <- read_lengths("sample.length_counts.tsv.zst")

# Extract metadata (depends on applied windowing and grouping)
length_bins(lengths)
window_metadata(lengths)
group_metadata(lengths)

# Extract long data frame with all data
length_data_frame(lengths)

# Extract length bins containing a given length range
length_data_frame(lengths, with_length_range = c(100L, 220L))

# Normalize to fractions across all length bins
length_data_frame(lengths, value = "fraction")

# Extract length bins containing a given length range
# And normalize to fractions of only those selected bins
length_data_frame(
  lengths,
  with_length_range = c(100L, 220L),
  value = "fraction",
  denominator = "selected_bins"
)

# For windowed counts, filter highly blacklisted windows/groups
length_data_frame(lengths, max_blacklisted_fraction = 0.3)

# Select counts for specific group(s) or window(s)
length_data_frame(lengths, groups = "t-cells", value = "fraction")
length_data_frame(lengths, window_idxs = 1L, value = "fraction")

# Extract counts as a matrix
length_counts_matrix(lengths)
```

## References

- Cristiano et al. Genome-wide cell-free DNA fragmentation in patients with cancer. Nature. 2019 Jun;570(7761):385-389. [https://doi.org/10.1038/s41586-019-1272-6](https://doi.org/10.1038/s41586-019-1272-6)

- Renaud et al. Elife. 2022 Jul 27;11:e71569. [https://doi.org/10.7554/eLife.71569](https://doi.org/10.7554/eLife.71569)
