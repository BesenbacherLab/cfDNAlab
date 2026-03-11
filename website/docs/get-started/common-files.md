# Common Files

The `cfDNAlab` commands use a common set of reference files, specific to your chosen assembly (e.g. `hg38`). This section describes the most common files and shows where to find them for `hg38`. The commands should work with any assembly.

The following external files are used by the main commands (some are optional):

| File              | Format                                 | Argument      | Where to get it                                                                                                                                    |
| ----------------- | -------------------------------------- | ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| Sample alignment  | Coordinate-sorted BAM + index (`.bai`) | `--bam`       | From your alignment pipeline or preprocessing workflow                                                                                             |
| Reference genome  | `.2bit`                                | `--ref-2bit`  | Download the exact assembly that matches your sample alignment. E.g. from: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit       |
| Blacklist regions | BED                                    | `--blacklist` | Download one or more assembly-matched blacklist BEDs. E.g. from: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz |
| Windows           | BED-like file                          | `--by-bed`    | From your own analysis design or external source.                                                                                                  |
| Intervals         | BED-like file                          | `--intervals` | From your own analysis design or external source.                                                                                                  |

In general, we recommend reusing the same blacklists across all steps in the analyses. Exceptions can be sample-specific blacklists that are only relevant for the last feature extraction steps.

The following files are made with cfDNAlab and passed to the main feature extraction commands (except the `Reference GC directory` which is necessary for calculating the GC-bias):

| File                      | Format                                    | Argument            | How to make it                                                                      |
| ------------------------- | ----------------------------------------- | ------------------- | ----------------------------------------------------------------------------------- |
| Sample scaling factors    | `.scaling_factors.tsv`                    | `--scaling-factors` | Create with `cfdna coverage-weights` per sample                                     |
| Sample GC correction file | `gc_bias_correction.npz`                  | `--gc-file`         | Create with `cfdna gc-bias` per sample                                              |
| Reference GC directory    | Directory produced by `cfdna ref-gc-bias` | `--ref-gc-dir`      | Create once per assembly, then reuse it across samples when calling `cfdna gc-bias` |

## Quick explanation

The **BAM** file (`--bam`) contains the actual sample-specific sequencing data, that you wish to extract features from. It has been aligned to a **reference genome** (`--ref-2bit`).
**Blacklists** show which regions of that genome are hard to map to. Those regions add noise to our features, so we usually exclude them from the analysis (`--blacklist`).

**Windows** are genomic intervals that we are specifically interested in (`--by-bed`). Depending on the command and settings, we either get separate features per window _or_ include only the specified window positions.

**Intervals** are fixed-size genomic intervals (`--intervals`) used by `cfdna midpoints` to create coverage profiles. 

**GC-bias**: Fragmentation patterns in cfDNA are vulnerable to GC-bias (see the GC-bias guide). `cfdna gc-bias` calculates this bias for a given sample BAM file, which can then be passed to the feature extraction commands to correct the bias in the features. Before we can calculate the sample bias though, we need to know the GC-bias in _reference genome_. This should be done _once_ per assembly using `cfdna ref-gc-bias`.

**Genomic smoothing**: When we specifically care about local changes in e.g. coverage, we might want all genomic regions to have roughly the same contributions to the features. This can be done by counting the coverage in large Mb bins and dividing the contributions of fragments by their overlapping bin coverage. `cfdna coverage-weights` calculates such scaling factors (although using a running-window, triangular weighting scheme for a smoother effect; similar to Gaussian smoothing but very efficient). Again, you calculate these scaling factors once per sample BAM file and pass them into feature extraction commands.

## Store paths in bash variables

To avoid writing the filepaths again and again, you can assign them to shell variables. Of course, some of the files need to be created first. Adjust these to your paths.

Tip: Use quotes around variable values and around variable expansions in commands. That protects you from path names with spaces.

```bash
# Project / assembly level
PROJECT_DIR="$HOME/cfdna_project"
REF_2BIT="$PROJECT_DIR/refs/hg38.2bit"
REF_GC_DIR="$PROJECT_DIR/refs/ref_gc"
WINDOWS="$PROJECT_DIR/regions/windows.bed"
INTERVALS="$PROJECT_DIR/regions/intervals.bed"

BLACKLIST_PRIMARY="$PROJECT_DIR/refs/blacklist/hg38-blacklist.bed"
BLACKLIST_EXTRA="$PROJECT_DIR/refs/blacklist/custom-mask.bed"
BLACKLIST_ARGS=(
  --blacklist "$BLACKLIST_PRIMARY"
  --blacklist "$BLACKLIST_EXTRA"
)

# Sample level
SAMPLE_ID="sample_01"
BAM="$PROJECT_DIR/inputs/$SAMPLE_ID.bam"
GC_FILE="$PROJECT_DIR/outputs/$SAMPLE_ID/gc_bias/gc_bias_correction.npz"
SCALING_FACTORS="$PROJECT_DIR/outputs/$SAMPLE_ID/coverage_weights/$SAMPLE_ID.scaling_factors.tsv"
```

If you only use one blacklist file, you can skip the array and keep a single variable:

```bash
BLACKLIST="$PROJECT_DIR/refs/blacklist/hg38-blacklist.bed"
```

## Use the variables in commands

**Note**: *The below examples only show the arguments for the shell variables, they are not full examples.*

Here is a sample-specific `gc-bias` call that reuses the variables above:

```bash
cfdna gc-bias \
  --bam "$BAM" \
  --output-dir "$PROJECT_DIR/outputs/$SAMPLE_ID/gc_bias" \
  --ref-2bit "$REF_2BIT" \
  --ref-gc-dir "$REF_GC_DIR" \
  "${BLACKLIST_ARGS[@]}"
```

Here is a downstream feature extraction call using the same shared variables plus the derived sample files:

```bash
cfdna midpoints \
  --bam "$BAM" \
  --output-dir "$PROJECT_DIR/outputs/$SAMPLE_ID/midpoints" \
  --intervals "$INTERVALS" \
  --ref-2bit "$REF_2BIT" \
  --gc-file "$GC_FILE" \
  --scaling-factors "$SCALING_FACTORS" \
  "${BLACKLIST_ARGS[@]}"
```

## Conceptual folder layout

While you can use any folder structure you want, the below layout conceptualizes the various types of files.

Keeping shared reference files separate from sample-specific outputs makes the pipeline easier to understand:

```text
project/
├── refs/
│   ├── hg38.2bit
│   ├── blacklist/
│   │   └── hg38-blacklist.bed
│   └── ref_gc/
├── inputs/
│   └── sample_01.bam
├── regions/
│   └── intervals.bed
└── outputs/
    └── sample_01/
        ├── gc_bias/
        │   └── gc_bias_correction.npz
        ├── coverage_weights/
        │    └── sample_01.scaling_factors.tsv
        └── lengths/
             └── sample_01.length_counts.npy
```

## Next step

Once you have these files in place, continue with [Quickstart](./quickstart.md) to find the right command and guide for your workflow.
