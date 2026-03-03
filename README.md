# cfDNAlab <a href='https://github.com/ludvigolsen/cfdnalab'><img src='https://raw.githubusercontent.com/ludvigolsen/cfdnalab/main/cfdnalab_logo_257x285_250dpi.png' align="right" height="160" /></a>

Ultra-fast command-line tools for analysis of cell-free DNA. Extract **fragment coverage**, **midpoint coverage**, and **fragment lengths** across the whole genome (or in windows) in mere seconds or minutes. Apply sample-specific GC correction and large-scale genomic smoothing.

Works on cfDNA **fragments** from either *paired-end* sequencing data or unpaired data where each read represents a full fragment. Written in rust for *speed*.

The commands are **highly flexible** with many options and good default settings. See [recipes](#recipes) in the end of this README for usage examples.

The package is in alpha-stage (being developed). Multiple additional commands are currently being built.
Suggest a tool or feature [here](https://github.com/LudvigOlsen/cfDNAlab/issues/new/choose)!

---

## Installation

### Compile from source

You may need a few dependencies that can be installed as a conda environment with:

```bash
conda create -n cfdnalab rust=1.87.0 zstandard perl fontconfig conda-forge::llvmdev conda-forge::clangdev
conda activate cfdnalab
```

Compile and install:

```bash
cargo install --git https://github.com/ludvigolsen/cfDNAlab --release --features cli,plotters
cfdna --help
# or clone + build
git clone https://github.com/ludvigolsen/cfDNAlab
cd cfdnalab && cargo build --release --features cli,plotters
target/release/cfdna --help
```

---

## Commands

The following commands are currently available:

| Command                              | Description                                                                                                                                                                                                            |
| ------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Feature extraction**               | Extract fragmentomics features                                                                                                                                                                                         |
| `cfdna fcoverage`                    | Count *fragment* coverage per position or aggregated in windows                                                                                                                                                        |
| `cfdna midpoints`                    | Count fragment *midpoint* coverage in fixed-size intervals, collapsed by groups across the genome<br />E.g. transcription factor binding sites, aggregated per transcription factor<br />Fast alternative to *Griffin* |
| `cfdna lengths`                      | Count fragment lengths<br />Defined as: `end(reverse) - start(forward)` for inwardly directed pairs only                                                                                                               |
| **Normalization**                    | Precompute normalization/correction factors to enable their use in the feature extraction commands                                                                                                                     |
| `cfdna gc-bias`, `cfdna ref-gc-bias` | Calculate GC-bias for correcting a sample in the main commands                                                                                                                                                         |
| `cfdna coverage-weights`             | Calculate scaling factors for normalizing/smoothing coverage across the genome                                                                                                                                         |
| **Conversion**                       | Convert BAM > frag > BAM or BAM > BAM                                                                                                                                                                                  |
| `cfdna bam-to-bam`                   | Apply our filters and/or write GC correction and coverage weight tags to a BAM file                                                                                                                                    |
| `cfdna bam-to-frag`                  | Write fragment coordinates to a "frag" file (bed-like tsv file)                                                                                                                                                        |
| `cfdna frag-to-bam`                  | Convert fragment coordinates to a single-read unpaired BAM file                                                                                                                                                        |

Planned: `cfdna ends` (end-motifs, breakpoint motifs), `cfdna fragment-kmers` (count kmers within fragments), `cfdna wps-peaks` (call windowed protection score peaks). Let us know what other fragmentomics features you would like to extract with `cfDNAlab`.

### Common options

- **GC bias correction**: Perform GC bias correction by weighting the contribution of each fragment by their GC content.

- **Blacklist filtering**: Supply BED files with regions to exclude. The implementation is specific to each tool (filtering of full fragments or just the overlapping positions).

- **Windowing**: Perform the command in genomic windows. Either a single global window (default), windows specified in a BED file, or via a fixed window size. Assign fragments to windows by how they overlap.

- **Genomic smoothing**: Scale the contribution of fragments by their coverage in megabase-scale overlapping bins. This reduces the effect of amplifications and deletions.

---

## FAQ

- How is *fragment* coverage different from the outputs of similar tools like `mosdepth` and `samtools`?
  - `mosdepth` counts the coverage of aligned bases per *read* [TODO: Not that simple]. `fcoverage` instead first collects the paired reads into a fragment and then counts the coverage of the aligned bases and (optionally) the gap between mate reads. (TODO on samtools!).

- How do you define a "fragment" in paired-end sequencing data?
  - We define the *fragment* as the bases from the start of the forward read till the end of the reverse read (`[start(forward), end(reverse))`) for *inwardly directed* pairs only (i.e., where `start(forward) <= start(reverse)`), as suggested by Wang, H. et al. 2025. Some methods exclude deletions and skipped-regions.

  Fragment visualization:

  ```text
  Reference 5' >>>>>>>>>>>>>>> 3'
  Fragment     |-------------|
  Forward   5' |>>>>>>>| 3'     
  Reverse        3' |<<<<<<<<| 5' 
  ``` 

- Should I order the BAM files differently to allow pairing of reads into fragments?
  - No, we expect BAM files to be *coordinate-sorted* and indexed.

- How do I run the command for unpaired data?
  - Most commands accept `--reads-are-fragments`. Each read is then assumed to represent a full fragment.

- How did you use LLMs (AI) in this project?
  - OpenAI's codex models were used for pair programming to speed up development and testing. All code for the released commands have been designed and validated by us.

---

## Recipes

We aim for high flexibility to make the commands useful for both established and novel use cases. This leads to commands having many options. The following recipes (examples) will get you quickly up and running with common cfDNA analyses.

The final example is a full pipeline for running everything (but without the explanations from the separate examples).

**Assembly**: The below examples use file names specific to the hg38 assembly, but any assembly (hg19 etc.) should work. Just be consistent, of course. Note that most commands use only the autosomes (`chr1-chr22`) by default (see `--chromosomes` in help files).

**Fragment length range**: The min/max fragment length range defaults to `30-1000bp`. This can be specified via `--min-fragment-length` and `--max-fragment-length`. We suggest keeping this range for the `ref-gc-bias`, `gc-bias`, and `coverage-weights` commands (unless you want to support longer fragments than 1000bp). In the downstream feature extraction commands you can then narrow the range, if you want.

### GC correction pipeline

Fragmentomics features are vulnerable to biases from various sample-handling and sequencing processes, such as PCR amplification. `cfDNAlab` commands thus allow the correction of the commonly observed **GC-bias**.

This requires only a few steps:

1) Calculate the "expected" GC bias in the reference genome assembly (e.g., hg38). This can be **REUSED for all samples** aligned to that assembly:

```bash

cfdna ref-gc-bias --help

# Run once per assembly
cfdna ref-gc-bias \
  --ref-2bit <path>/hg38.2bit \
  --output-dir <ref_gc_directory> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed  # As many as you want

```

2) Calculate the GC-bias correction factors per sample:

```bash

cfdna gc-bias --help

cfdna gc-bias \
  --bam <sample>.bam \
  --output-dir <sample_directory>/gc_bias \
  --n-threads 12 \
  --ref-2bit <path>/hg38.2bit \
  --ref-gc-dir <ref_gc_directory> \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed  # Should match those specified in ref-gc-bias!

```

3) Provide the correction factors when running the feature extraction commands **on the same BAM file**:

```bash

cfdna fcoverage \
  --bam <sample>.bam \
  ... \  # See fcoverage example
  --gc-file <sample_directory>/gc_bias/gc_bias_correction.npz \
  --ref-2bit <path>/hg38.2bit

```

If you prefer a different/custom GC-bias tool, the feature extraction commands also accept reading a GC weight (how much a fragment should contribute) from an aux tag in the BAM file:

```bash

cfdna fcoverage \
  --bam <sample>.bam \
  ... \  # See fcoverage example
  --gc-tag 'GC'

```

### Genomic smoothing pipeline

For some commands, like `cfdna midpoints`, you may want all genomic regions to have approximately the same contribution to the features. E.g., to reduce the effect of copy number alterations.

**Simplified**, this can be achieved by calculating the fragment coverage in a kilo/megabase resolution and dividing the contribution of each fragment (`1.0` or the gc-weight) by the coverage value.

**More detailed**, for a more smooth scaling, `cfdna coverage-weights` builds a smoothed normalization map using a sliding window:

**A**) It splits the genome into "stride-bins" (default: 500kb) and counts the average positional fragment coverage in each bin.

**B**) It smoothes each bin with a triangular weighting kernel, that weights the coverage of the neighbouring stride-bins by how many overlapping megabins (default: 5Mb) they are part of. E.g.:

Using a megabin-size of `6` and stride size of `2` for demonstrational purposes:

**Stride bins** (fixed along genome, each with an average positional coverage):

`[A] [B] [C] [D] [E] [F] [G] ...`

**Overlapping megabins** (`MB*`) each cover 3 stride-bins.
**`W_D`** weights each stride-bin by how many `D`-overlapping megabins it is part of.
Stride-bin `B` is only part of one megabin that overlaps `D`, so its (unnormalized) weight is 1.
In contrast, stride-bin `D` is naturally part of all three megabins, so its weight is 3:

<pre>

<i>MB1</i>: [A][B][C]

MB2:    [B][C][<b>D</b>]

MB3:       [C][<b>D</b>][E]

MB4:          [<b>D</b>][E][F]

<i>MB5</i>:             [E][F][G]

W_D: [0][1][2][3][2][1][0]

</pre>

$$smoothCoverage_{D} = (0A + 1B + 2C + 3D + 2E + 1F + 0G) / (1+2+3+2+1)$$

**C**) Finally, the values are **inverted** with $1/smoothCoverage$ to become multiplicative scaling factors (one per stride-bin). A fragment's contribution (`1.0` or the gc-weight) can then be scaled by multiplying by the scaling factor of the stride-bin it's located in.

You can think of this approach as a very fast alternative to e.g. Gaussian smoothing.

The genomic smoothing can be achieved in two steps:

1) Calculate the coverage-based scaling factors:

```bash

cfdna coverage-weights --help

cfdna coverage-weights \
  --bam <sample>.bam \
  --output-dir <sample_directory>/coverage_weights \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed  # Tip: Use the same blacklists everywhere / per sample!

```

2) Provide the scaling factors when running the feature extraction commands **on the same BAM file**:

```bash

cfdna midpoints \
  --bam <sample>.bam \
  ...  # See midpoints example
  --scaling-factors <sample_directory>/coverage_weights/<sample_id>.scaling_factors.tsv

```

### Fragment coverage

Fragment coverage measures how many fragments overlap each genomic position. In contrast to many non-cfDNA-tools, we (optionally) count the gap between paired reads along with the aligned bases of the reads. We avoid double counting when reads overlap.

When no GC correction or genomic smoothing is applied, each fragment is counted as `1` in the overlapping (aligned / gap) positions. Using GC correction and/or genomic smoothing changes this to a weight (floating point).

```bash

cfdna fcoverage --help

cfdna fcoverage \
  --bam <sample>.bam \                          # Coordinate-sorted bam file with cfDNA
  --output-dir <sample_directory>/coverage \    # Where to write files
  --output-prefix <sample_id> \                 # A file prefix to identify the sample (optional)
  --n-threads 12 \                              # Use 12 CPU cores (speed vs. RAM tradeoff)
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \  # Tip: Use the same blacklists everywhere / per sample!

  # OPTIONS:

  # Average per 1Mb positions
  --by-size 1000000 \
  --per-window 'average' \
  # OR sums per interval in a bed file
  --by-bed <path>/<some_intervals>.bed \
  --per-window 'total' \

  # Add GC correction and / or genomic smoothing (see above)
  --gc-file ... \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors ...

```

### Fragment lengths

Multiple studies have used fragment lengths (count distributions) to detect cancer [REFS].

```bash

cfdna lengths --help

cfdna lengths \
  --bam <sample>.bam \                          # Coordinate-sorted bam file with cfDNA
  --output-dir <sample_directory>/lengths \     # Where to write files
  --output-prefix <sample_id> \                 # A file prefix to identify the sample (optional)
  --n-threads 12 \                              # Use 12 CPU cores (speed vs. RAM tradeoff)
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \  # Tip: Use the same blacklists everywhere / per sample!

  # OPTIONS:

  # Adjust lengths to indels
  --indel-mode 'adjust' \

  # Separate counts per 1Mb positions
  --by-size 1000000 \
  # OR supply a bed file
  --by-bed <path>/<some_intervals>.bed \

  # Add GC correction and / or genomic smoothing (see above)
  --gc-file ... \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors ...

```

### Fragment midpoint profiles

Multiple studies have used profiled the midpoint coverage around e.g. transcription factor binding sites (summed per transcription factor, per position) [REFS]. This can inform about the binding activity of different transcription factors related to cancer.

```bash

cfdna midpoints --help

cfdna midpoints \
  --bam <sample>.bam \                          # Coordinate-sorted bam file with cfDNA
  --output-dir <sample_directory>/midpoints \   # Where to write files
  --output-prefix <sample_id> \                 # A file prefix to identify the sample (optional)
  --n-threads 12 \                              # Use 12 CPU cores (speed vs. RAM tradeoff)
  --intervals <fixed_size_intervals>.tsv \      # The grouped fixed-size intervals (see --help)
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <path>/<another_blacklist>.bed \  # Tip: Use the same blacklists everywhere / per sample!

  # OPTIONS:

  # Separate counts per 10bp lengths (last edge is exclusive, 1000bp is excluded)
  --length-bins {30..1000..10} \

  # Add GC correction and / or genomic smoothing (see above)
  --gc-file ... \
  --ref-2bit <path>/hg38.2bit \
  --scaling-factors ...

```

[TODO: Note on how to get griffin-like profiles]

The **intervals** must have the same fixed size. A common binding site window size is `2001bp`, centered around the binding site center. The expected columns are: `chromosome, start, end, group_name` (where `group_name` is the group to collapse profiles by, e.g., the transcription factor ID). The intervals should be sorted by chromosome and start-coordinates.

Consider removing intervals that lie closer than half the maximum fragment length from any blacklisted region, to reduce mappability biases.

### Everything combined

[TODO: Add output-prefix for remaining commands]

```bash

BAM="..."                  # The cfDNA sample
OUT="..."                  # Sample specific directory
SAMPLE_NAME="..."          # Sample ID for filename prefix
BLACKLIST="..."            # One or more (repeat argument per file)
ASSEMBLY="..."             # E.g. hg38.2bit
REF_GC="..."               # Precompute with `cfdna ref-gc-bias`
MIDPOINT_INTERVALS="..."   # Fixed-size intervals BED-like tsv-file
THREADS=12
MINLENGTH=30
MAXLENGTH=1000

# GC bias correction matrix
cfdna gc-bias --bam $BAM --output-dir $OUT/gc_bias --ref-2bit $ASSEMBLY --ref-gc-dir $REF_GC --blacklist $BLACKLIST --min-fragment-length $MINLENGTH --max-fragment-length $MAXLENGTH --n-threads $THREADS 

# Coverage weights for genomic smoothing
cfdna coverage-weights --bam $BAM --output-dir $OUT/coverage_weights --output-prefix $SAMPLE_NAME --blacklist $BLACKLIST --min-fragment-length $MINLENGTH --max-fragment-length $MAXLENGTH --n-threads $THREADS 

# Fragment coverage per position
cfdna fcoverage --bam $BAM --output-dir $OUT/coverage --output-prefix $SAMPLE_NAME --blacklist $BLACKLIST --min-fragment-length $MINLENGTH --max-fragment-length $MAXLENGTH --gc-file $OUT/gc_bias/gc_bias_correction.npz --ref-2bit $ASSEMBLY --n-threads $THREADS 

# Fragment coverage in 5Mb bins (averaged)
cfdna fcoverage --bam $BAM --output-dir $OUT/coverage_per_5Mb --output-prefix $SAMPLE_NAME --blacklist $BLACKLIST --min-fragment-length $MINLENGTH --max-fragment-length $MAXLENGTH --gc-file $OUT/gc_bias/gc_bias_correction.npz --ref-2bit $ASSEMBLY --n-threads $THREADS --by-size 5000000 --per-window 'average'

# Fragment lengths (global)
cfdna lengths --bam $BAM --output-dir $OUT/lengths_$MINLENGTH_$MAXLENGTH --blacklist $BLACKLIST  --min-fragment-length $MINLENGTH --max-fragment-length $MAXLENGTH --gc-file $OUT/gc_bias/gc_bias_correction.npz --ref-2bit $ASSEMBLY --scaling-factors $OUT/coverage_weights/$SAMPLE_NAME.scaling_factors.tsv --n-threads $THREADS 

# Fragment lengths in 5Mb bins 
# E.g., to calculate short/long ratios from (100-150bp, 151-220bp)
cfdna lengths --bam $BAM --output-dir $OUT/lengths_per_5mb_100_220 --by-size 5000000 --blacklist $BLACKLIST --min-fragment-length 100 --max-fragment-length 220 --gc-file $OUT/gc_bias/gc_bias_correction.npz --ref-2bit $ASSEMBLY --scaling-factors $OUT/coverage_weights/$SAMPLE_NAME.scaling_factors.tsv --n-threads $THREADS

# Midpoint profiles (very fast alternative to Griffin)
cfdna midpoints --bam $BAM --output-dir $OUT/midpoints --intervals $MIDPOINT_INTERVALS --blacklist $BLACKLIST --length-bins $MINLENGTH $(($MAXLENGTH+1)) --gc-file $OUT/gc_bias/gc_bias_correction.npz --ref-2bit $ASSEMBLY --scaling-factors $OUT/coverage_weights/$SAMPLE_NAME.scaling_factors.tsv --n-threads $THREADS

```

---

## Unpaired (--reads-are-fragments)

If you have Nanopore-sequenced cell-free DNA (or similar) where each read represents the full fragment, you can supply the `--reads-are-fragments` flag. This will consider each read a full fragment.

Simplest example:

```bash

cfdna fcoverage --help

cfdna fcoverage \
  --bam <sample>.bam \                          # Coordinate-sorted bam file with cfDNA
  --output-dir <sample_directory>/coverage \    # Where to write files
  --output-prefix <sample_id> \                 # A file prefix to identify the sample (optional)
  --reads-are-fragments                         # Consider each read a fragment

```

---

## Output formats

The various commands produce either Numpy arrays (`.npy`, `.npz`), `.tsv.zst`, `.bedgraph.zst`, or text/json files.

Numpy files can be read with python:

```python

import numpy as np

x = np.load("path_to/file.npy")
y = np.load("path_to/file.npz")

```

zstd-compressed files can be decompressed with:

```bash

zstd -d path_to/file.tsv.zst

# Read with cat
zstdcat path_to/file.tsv.zst

```

Tip: Bedgraph files can be converted to bigwig files for indexed lookup.

File columns: [TODO: Check correctness]

- Bedgraph: `chromosome  start  end  value`
- Scaling factors TSV: `chromosome  start  end  scaling_factor`

---

## References

- Wang, H., Mennea, P.D., Chan, Y.K.E. et al. A standardized framework for robust fragmentomic feature extraction from cell-free DNA sequencing data. Genome Biol 26, 141 (2025). <https://doi.org/10.1186/s13059-025-03607-5>
