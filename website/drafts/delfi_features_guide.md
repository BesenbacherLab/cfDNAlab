# Calculating DELFI-like features

Multiple papers by authors related to DELFI Diagnostics have suggested that features such as count ratios of short/long fragments and chromosome arm-level coverage reflect tumor presence [REFS].

While reproducing their results is beyond the scope of this guide, let's create those two feature sets with cfDNAlab and a bit of Python. This will not be a complete recreation of DELFI's features, since we use fragment-level GC-bias correction, but you can adjust the Python code if that is your goal. The cfDNAlab parts should be the same (expect for the addition of fragment-level GC-bias correction).

[TODO: Make it clear that we don't do exactly what they do (not their GC correction or weird normalization)]

## File paths

We will start by assigning variables with paths to the needed input files. The guide will use hg38 but feel free to use another assembly. The below paths should be modified to point to their location on your disk.

```bash
# Project / assembly level
PROJECT_DIR="$HOME/delfi_features"
REF_2BIT="$PROJECT_DIR/refs/hg38.2bit"
ARM_WINDOWS="$PROJECT_DIR/regions/chromosome_arms.bed"
BLACKLIST="$PROJECT_DIR/refs/blacklist/hg38-blacklist.bed"
REF_GC_DIR="$PROJECT_DIR/refs/ref_gc" # May not exist yet

# Sample level
SAMPLE_ID="sample_01"
BAM="$PROJECT_DIR/inputs/$SAMPLE_ID.bam"

# Available CPU cores on your system 
# Higher is faster but uses more memory
N_CORES=12

# Output directories (sample-level)
GC_DIR="$PROJECT_DIR/output/$SAMPLE_ID/gc_bias"
LEN_DIR="$PROJECT_DIR/output/$SAMPLE_ID/lengths"
COV_DIR="$PROJECT_DIR/output/$SAMPLE_ID/coverage"
```

## 1. GC correction

The old Cristiano et al. () paper does not use fragment-level GC correction, but since this is one of the benefits of cfDNAlab, we will include it.

The first computational step is to calculate the GC-bias correction matrix for our sample. This requires knowing the expected GC-bias in the reference genome (`cfdna ref-gc-bias`). Since we only need to calculate this once for a given assembly, you may already have computed it. If not, here is the command to do so:

```bash
# Run once per assembly
cfdna ref-gc-bias \
  --ref-2bit "$REF_2BIT" \
  --output-dir "$REF_GC_DIR" \
  --output-prefix hg38 \
  --n-threads $N_CORES \
  --blacklist "$BLACKLIST"
```

With the reference bias ready, we can calculate the gc-bias correction matrix for our BAM file:

```bash
GC_FILE="$GC_DIR/gc_bias_correction.npz"

cfdna gc-bias \
  --bam "$BAM" \
  --output-dir "$GC_DIR" \
  --n-threads $N_CORES \
  --ref-2bit "$REF_2BIT" \
  --ref-gc-file "$REF_GC_DIR/hg38.ref_gc_package.npz" \
  --blacklist "$BLACKLIST"
```

## 2. Fragment lengths

With the GC-bias correction ready, it is time to count the fragment lengths in 5Mb windows. DELFI defines short and long fragment lengths as 100-150bp and 151-220bp, respectively, so we can count those two bins directly:

[TODO: Figure out exactly what the DELFI settings are]
[TODO: How do DELFI actually use GC correction on lengths? Do they?]

```bash
cfdna lengths \
  --bam "$BAM" \
  --output-dir "$LEN_DIR" \
  --output-prefix "$SAMPLE_ID" \
  --n-threads $N_CORES \
  --by-size 5000000 \
  --length-bins 100 151 221 \
  --ref-2bit "$REF_2BIT" \
  --blacklist "$BLACKLIST" \
  --gc-file "$GC_FILE"
```

**Note**: The calculation of the short/long ratios are shown later in the guide.

## 3. Fragment coverage

The arm-level Z-scores features start by counting the coverage in 50kb windows and then later aggregate per chromosome arm. In this example, we won't limit the fragment lengths to 100-220bp, but you could choose to do that.

[TODO: Do they use all fragment lengths?]

```bash
cfdna fcoverage \
  --bam "$BAM" \
  --output-dir "$COV_DIR" \
  --output-prefix "$SAMPLE_ID" \
  --n-threads $N_CORES \
  --by-size 50000 \
  --per-window 'average' \
  --ref-2bit "$REF_2BIT" \
  --blacklist "$BLACKLIST" \
  --gc-file "$GC_FILE"
```

## 4. Short/long feature calculation with Python

To calculate the short/long fragment length ratios, we take the following steps:

1) Load the length counts from the saved TSV table.
2) Read short (100-150bp) and long (151-220bp) counts.
3) Calculate the ratios of those counts per 5Mb window.


```bash
LEN_COUNTS="$LEN_DIR/$SAMPLE_ID.length_counts.tsv.zst"
DELFI_LENGTHS="$LEN_DIR/$SAMPLE_ID.delfi_short_long.npy"
```

```python
from pathlib import Path
import numpy as np
import pandas as pd

# Path to files 
# Change these to fit with your paths or pass through the command line
sample_id = "sample_01"
len_dir = Path.home() / "delfi_features" / "output" / sample_id / "lengths"

# Load lengths and extract short/long columns.
# count_100_151 is [100,151), count_151_221 is [151,221).
length_counts = pd.read_csv(len_dir / f"{sample_id}.length_counts.tsv.zst", sep="\t")
short_counts = length_counts["count_100_151"].to_numpy()
long_counts = length_counts["count_151_221"].to_numpy()

# Divide the two count vectors across all windows
# Some windows might have little coverage due to the blacklisting
# so we ensure 0-divisions become NaN instead of raising an error
short_long_ratio = np.divide(
    short_counts,
    long_counts,
    out=np.full(short_counts.shape, np.nan),
    where=long_counts != 0,
)

# Save the features to disk
np.save(len_dir / f"{sample_id}.delfi_short_long.npy", short_long_ratio)
```

## 5. Arm-level depth calculation with Python



To get the arm-level coverage, DELFI first counts the coverage in 50kb windows, normalize them and 

[TODO: The actual DELFI approach needs to be found in their code]




To calculate the arm-level Z-scores, we first count the fragment coverage in 50kb windows, then normalize them by the sample's read-depth, log2 transforms them, and calculate the mean value per chromosome arm. 

**NOTE 1**: DELFI uses a different sample-level read-depth normalization as part of their bin-level GC-bias correction.

**NOTE 2**: DELFI standardizes the final scores by the mean and standard deviation of 50 healthy samples, but since we don't have these, we will stop before doing this. 

In this step, we count the fragment coverage in 50kb windows:



Since we already have the fragment-level GC correction, we will 
