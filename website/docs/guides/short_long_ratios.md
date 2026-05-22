# Short/Long Ratios

Multiple studies have suggested that count ratios of short/long fragments across the genome reflect tumor presence (Cristiano et al. 2019).

With cfDNAlab, we can calculate these ratios in two steps:

1)  Count fragments of sizes 100-150bp and 151-220bp in 5Mb bins with `cfdna lengths`.

2)  Load the output in either Python or R and calculate the ratios.

## Count fragments

First set the paths to your BAM file, blacklist, and output directory. This is just an example of a setup, you can of course just specify these directly in the `cfdna lengths` call.

```bash
# Working directory
PROJECT_DIR="$HOME/short_long_ratios"

# Sample level
SAMPLE_ID="sample_01"
BAM="$PROJECT_DIR/inputs/$SAMPLE_ID.bam"

# File with regions to exclude
# You can specify multiple blacklists by repeating `--blacklist <path>`
BLACKLIST="$PROJECT_DIR/hg38-blacklist.bed"

# Output directory
LEN_DIR="$PROJECT_DIR/output/$SAMPLE_ID/lengths"

# Available CPU cores on your system 
# Higher is faster but uses more memory
N_CORES=12
```

Count the short (100-150bp) and long (151-220bp) fragments in 5Mb bins across the genome:

```bash
cfdna lengths \
  --bam "$BAM" \
  --output-dir "$LEN_DIR" \
  --output-prefix "$SAMPLE_ID" \
  --n-threads $N_CORES \
  --by-size 5000000 \
  --length-bins 100 151 221 \
  --blacklist "$BLACKLIST"
```

This creates a file like `"sample_01.length_counts.tsv.zst"`.

Next, we will see how to load the output in either Python or R with our loader packages and then calculate the ratios.

## Calculate ratios in R

First, we load in the TSV file and convert it to a wide-format data frame. For demonstrational purposes, we skip positional bins where more than 60% (arbitrary) of bases were blacklisted.

```r
library(cfdnalab)

# Specify path to the output file
path <- ".../sample_01.length_counts.tsv.zst"

# Load the TSV file
lengths <- read_lengths(path)

# If more than 60% of the 5Mb bin was blacklisted, exclude it
# (Arbitrary choice to show the filtering option)
max_blacklisted_fraction <- 0.6

# Extract as a data.frame
length_counts <- length_data_frame(
  lengths, 
  max_blacklisted_fraction = max_blacklisted_fraction, 
  keep_wide = TRUE
)

head(length_counts)
```

Now, we can calculate the ratio per bin as such:

```r
length_counts <- length_counts |>
  dplyr::mutate(short_long_ratio = count_100_151 / count_151_221)
```

## Calculate ratios in Python

First, we load in the TSV file and convert it to a wide-format data frame. For demonstrational purposes, we skip positional bins where more than 60% (arbitrary) of bases were blacklisted.

```python
import cfdnalab as cfl

# Specify path to the output file
path = ".../sample_01.length_counts.tsv.zst"

# Load the TSV file
lengths = cfl.read_lengths(path)

# If more than 60% of the 5Mb bin was blacklisted, exclude it
# (arbitrary choice to show the filtering option)
max_blacklisted_fraction = 0.6

# Extract as a Pandas DataFrame
length_counts = lengths.data_frame(
  max_blacklisted_fraction = max_blacklisted_fraction,
  keep_wide = True
)

length_counts.head()
```

Now, we can calculate the ratio per bin as such:

```python
length_counts["short_long_ratio"] = (
  length_counts["count_100_151"] / length_counts["count_151_221"]
)
```

## References

- Cristiano et al. Genome-wide cell-free DNA fragmentation in patients with cancer. Nature. 2019 Jun;570(7761):385-389. [https://doi.org/10.1038/s41586-019-1272-6](https://doi.org/10.1038/s41586-019-1272-6)
