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
  --scaling-factors <sample_directory>/count_weights/<sample_id>.fragment_counts.scaling_factors.tsv
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
  --scaling-factors <sample_directory>/gc_corrected_count_weights/<sample_id>.fragment_counts.scaling_factors.tsv
```

The intervals must have the same fixed size. The expected columns are: `chromosome, start, end, group_name` (where `group_name` is the group to collapse profiles by, e.g., the transcription factor ID). The intervals should be sorted by chromosome and start coordinates.

## Load midpoint profiles in Python

The main output is a Zarr store:

```text
<sample_id>.midpoint_profiles.zarr/
```

The Python helper package gives the shortest route to a plotting table:

```python
import cfdnalab as cfl

profiles = cfl.read_midpoints("<sample_id>.midpoint_profiles.zarr")

df = profiles.data_frame(groups="LYL1", with_lengths=167)
```

Select groups or length bins before converting to a data frame. Expanding all groups, length bins, and
positions can create a much larger object than the Zarr store.

## Load midpoint profiles in R

The R helper package follows the same idea but uses ordinary R functions:

```r
library(cfdnalab)

profiles <- read_midpoints("<sample_id>.midpoint_profiles.zarr")

df <- midpoint_data_frame(
  profiles,
  groups = "LYL1",
  with_lengths = 167
)
```
