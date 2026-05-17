# Extract Fragment Midpoint Profiles

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
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
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
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv
```

The intervals must have the same fixed size. The expected columns are: `chromosome, start, end, group_name` (where `group_name` is the group to collapse profiles by, e.g., the transcription factor ID). The intervals should be sorted by chromosome and start coordinates.

## Load midpoint profiles in Python

The main output is a Zarr store:

```text
<sample_id>.midpoint_profiles.zarr/
```

```python
import numpy as np
import pandas as pd
import zarr

profiles = zarr.open_group(
    "<sample_id>.midpoint_profiles.zarr",
    mode="r",
    zarr_format=3,
)

group_index = 0
counts = profiles["counts"][group_index, :, :]
length_index, position_index = np.indices(counts.shape)

df = pd.DataFrame({
    "group_name": profiles["group_name"][group_index],
    "eligible_intervals": int(profiles["eligible_intervals"][group_index]),
    "length_start_bp": profiles["length_start_bp"][:][length_index.ravel()],
    "length_end_bp": profiles["length_end_bp"][:][length_index.ravel()],
    "position_bin_start_bp": profiles["position_bin_start_bp"][:][position_index.ravel()],
    "position_bin_end_bp": profiles["position_bin_end_bp"][:][position_index.ravel()],
    "count": counts.ravel(),
})
```

Subset before converting to a dataframe. Expanding all groups, length bins, and
positions can create a much larger object than the Zarr store.

## Load midpoint profiles in R

Use `Rarr` to read the arrays directly:

```r
library(Rarr)

store <- "<sample_id>.midpoint_profiles.zarr"

counts <- read_zarr_array(
  file.path(store, "counts"),
  index = list(NULL, NULL, NULL)
)
group_name <- read_zarr_array(file.path(store, "group_name"), index = list(NULL))
eligible_intervals <- read_zarr_array(file.path(store, "eligible_intervals"), index = list(NULL))
length_start_bp <- read_zarr_array(file.path(store, "length_start_bp"), index = list(NULL))
length_end_bp <- read_zarr_array(file.path(store, "length_end_bp"), index = list(NULL))
position_bin_start_bp <- read_zarr_array(file.path(store, "position_bin_start_bp"), index = list(NULL))
position_bin_end_bp <- read_zarr_array(file.path(store, "position_bin_end_bp"), index = list(NULL))

group_index <- 1L
profile <- counts[group_index, , ]

df <- do.call(rbind, lapply(seq_along(length_start_bp), function(length_index) {
  do.call(rbind, lapply(seq_along(position_bin_start_bp), function(position_index) {
    data.frame(
      group_name = group_name[group_index],
      eligible_intervals = eligible_intervals[group_index],
      length_start_bp = length_start_bp[length_index],
      length_end_bp = length_end_bp[length_index],
      position_bin_start_bp = position_bin_start_bp[position_index],
      position_bin_end_bp = position_bin_end_bp[position_index],
      count = profile[length_index, position_index]
    )
  }))
}))
```
