# Extract Fragment Midpoint Profiles

Multiple studies have used midpoint coverage profiles around e.g. transcription factor binding sites (summed per transcription factor, per position)(Doebley et al., 2022). This can inform about the binding activity of different transcription factors related to cancer.

## Examples

The following examples show different aspects of the `cfdna midpoints` command. They can of course be combined in a multitude of ways, but for simplification we just show one aspect at a time.

### Base command

Extract midpoint profiles given a set of fixed-size intervals:

```bash
cfdna midpoints --help

cfdna midpoints \
  --bam <sample>.bam \
  --output-dir <sample_directory>/midpoints \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --intervals <fixed_size_intervals>.tsv \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed
```

### Count in length bins

By default, a single profile is created per interval group for all fragments with a length of 30-1000bp. To change those limits or to count in bins of fragment lengths, you can specify the `--length-bins` as `start:stop:step` or as space-separated edges like `--length-bins 100 151 221`. The default is `--length-bins 30 1001`.

```bash
cfdna lengths \
  ... \
  # Bin every 10bp from 80-499bp
  --length-bins 80:500:10
  # OR Specify edges directly (last end is exclusive)
  --length-bins 100 151 221

```

### GC-bias correction example

To correct the sample-specific GC-bias, you need to precompute the correction matrix (see the [GC-bias guide](./correct_gc_bias_guide.md)). Then you provide that correction file as:

```bash
cfdna midpoints \
  --bam <sample>.bam \
  --output-dir <sample_directory>/midpoints \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --intervals <fixed_size_intervals>.tsv \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Precomputed GC correction file and reference genome
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit
```

### Genomic smoothing example

When you're interested in local, relative coverage changes instead of large-scale changes from CNVs etc., you can use genomic smoothing to weight contributions more similarly across genomic regions. In some analyses, this can be thought of as a copy-number normalization.

Since we're counting fragments, we suggest using the `cfdna fragment-count-weights` command for calculating scaling factors (see the [genomic smoothing guide](./genomic_smoothing_guide.md)).

```bash
cfdna midpoints \
  --bam <sample>.bam \
  --output-dir <sample_directory>/midpoints \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --intervals <fixed_size_intervals>.tsv \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Precomputed scaling factors
  --scaling-factors <sample_directory>/count_weights/<sample_id>.fragment_counts.scaling_factors.tsv
```

### GC-bias correction + genomic smoothing

**NOTE**: *Requires* using GC-bias correction when calculating the scaling factors to avoid double-correction.

```bash
cfdna midpoints \
  --bam <sample>.bam \
  --output-dir <sample_directory>/midpoints \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --intervals <fixed_size_intervals>.tsv \
  --blacklist <path>/hg38-blacklist.v2.bed \
  # Precomputed GC correction file and reference genome
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.zarr \
  --ref-2bit <path>/hg38.2bit \
  # Precomputed scaling factors
  --scaling-factors <sample_directory>/gc_corrected_count_weights/<sample_id>.fragment_counts.scaling_factors.tsv
```

The intervals must have the same fixed size. The expected columns are: `chromosome, start, end, group_name` (where `group_name` is the group to collapse profiles by, e.g., the transcription factor ID). The intervals should be sorted by chromosome and start coordinates.

## Load midpoint profiles in Python

The main output is a Zarr store:

```text
<sample_id>.midpoint_profiles.zarr/
```

Load it with the Python helper package and extract a profile as a long data frame:

```python
import cfdnalab as cfl

profiles = cfl.read_midpoints("<sample_id>.midpoint_profiles.zarr")

df = profiles.data_frame(groups="LYL1", with_lengths=167)
```

We suggest selecting subsets of groups or length bins before converting to a data frame. Expanding all groups, length bins, and positions can create a very large object.

## Load midpoint profiles in R

The main output is a Zarr store:

```text
<sample_id>.midpoint_profiles.zarr/
```

Load it with the R helper package and extract a profile as a long data frame:

```r
library(cfdnalab)

profiles <- read_midpoints("<sample_id>.midpoint_profiles.zarr")

df <- midpoint_data_frame(
  profiles,
  groups = "LYL1",
  with_lengths = 167
)
```

We suggest selecting subsets of groups or length bins before converting to a data frame. Expanding all groups, length bins, and positions can create a very large object.

## References

- Doebley et al. A framework for clinical cancer subtyping from nucleosome profiling of cell-free DNA. Nat Commun 13, 7475 (2022). [https://doi.org/10.1038/s41467-022-35076-w](https://doi.org/10.1038/s41467-022-35076-w)
