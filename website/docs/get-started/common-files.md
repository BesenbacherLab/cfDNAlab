# Common Files

The `cfDNAlab` commands use a common set of reference files specific to your chosen assembly (e.g. `hg38`). This section describes the most common files and shows where to find them for `hg38`. The commands should work with any assembly.

The following external files are used by the main commands (some are optional):

| File              | Format                                 | Arguments                      | Where to get it                                                                                                                                    |
| ----------------- | -------------------------------------- | ------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| Sample alignment  | Coordinate-sorted BAM + index (`.bai`) | `--bam`                        | From your alignment pipeline or preprocessing workflow                                                                                             |
| Reference genome  | `.2bit`                                | `--ref-2bit`                   | Download the exact assembly that matches your sample alignment. E.g. from: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit       |
| Blacklist regions | BED                                    | `--blacklist`                  | Download one or more assembly-matched blacklist BEDs. E.g. from: https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz |
| Windows           | BED-like file                          | `--by-bed`, `--by-grouped-bed` | From your own analysis design or external source.                                                                                                  |
| Intervals         | BED-like file                          | `--intervals`                  | From your own analysis design or external source.                                                                                                  |

In general, we recommend reusing the same blacklists across all steps in the analyses. Exceptions can be sample-specific blacklists that are only relevant for the last feature extraction steps.

The following files are made with cfDNAlab and passed to the main feature extraction commands (except the `Reference GC file` which is necessary for calculating the GC-bias):

| File                      | Format                                  | Argument            | How to make it                                                                        |
| ------------------------- | --------------------------------------- | ------------------- | ------------------------------------------------------------------------------------- |
| Sample scaling factors    | `.scaling_factors.tsv`                  | `--scaling-factors` | Create with `cfdna fragment-count-weights` and/or `cfdna coverage-weights` per sample |
| Sample GC correction file | `gc_bias_correction.npz`                | `--gc-file`         | Create with `cfdna gc-bias` per sample                                                |
| Reference GC file         | Package produced by `cfdna ref-gc-bias` | `--ref-gc-file`     | Create once per assembly, then reuse it across samples when calling `cfdna gc-bias`   |

## Quick explanation

The **BAM** file (`--bam`) contains the actual sample-specific sequencing data, that you wish to extract features from. It has been aligned to a **reference genome** (`--ref-2bit`).
**Blacklists** show which regions of that genome are hard to map to. Those regions add noise to our features, so we usually exclude them from the analysis (`--blacklist`).

**Windows** are genomic intervals that we are specifically interested in (`--by-bed`, `--by-grouped-bed`). Depending on the command and settings, we either get separate features per window/group _or_ include only the specified window positions.

**Intervals** are fixed-size genomic intervals (`--intervals`) used by `cfdna midpoints` to create midpoint coverage profiles. 

**GC-bias**: Fragmentation patterns in cfDNA are vulnerable to GC-bias (see the GC-bias guide). `cfdna gc-bias` calculates this bias for a given sample BAM file, which can then be passed to the feature extraction commands to correct the bias in the features. Before we can calculate the sample bias though, we need to know the GC-bias in _reference genome_. This should be done _once_ per assembly using `cfdna ref-gc-bias`.

**Genomic smoothing**: Use genomic smoothing when you care about **local** changes in fragment counts or coverage. This makes all non-blacklisted genomic regions contribute roughly the same total weight to the features. There are two related modes: **coverage**, where longer fragments count more because they cover more positions, and **fragment counts**, where each fragment has the same total weight regardless of length. We calculate local fragment counts or coverage in large genomic windows, then divide each fragment’s contribution to features by the count or coverage of the windows it overlaps. `cfdna fragment-count-weights` and `cfdna coverage-weights` calculate these scaling factors once per sample BAM file. They use a running-window, triangular weighting scheme, which gives a smooth effect similar to Gaussian smoothing. Pass the resulting scaling factors into feature extraction commands.

## Store paths in bash variables

To avoid writing the filepaths again and again, you can assign them to shell variables. Of course, some of the files need to be created first. Adjust these to your paths.

Tip: Use quotes around variable values and around variable expansions in commands. That protects you from path names with spaces.

```bash
# Project / assembly level
PROJECT_DIR="$HOME/cfdna_project"
REF_2BIT="$PROJECT_DIR/refs/hg38.2bit"
REF_GC_FILE="$PROJECT_DIR/refs/ref_gc/hg38.ref_gc_package.npz"
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
SCALING_FACTORS="$PROJECT_DIR/outputs/$SAMPLE_ID/scaling_factors/$SAMPLE_ID.scaling_factors.tsv"
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
  --ref-gc-file "$REF_GC_FILE" \
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
        │   └── sample_01.gc_bias_correction.npz
        ├── scaling_factors/
        │    ├── sample_01.coverage.scaling_factors.tsv
        │    └── sample_01.fragment_counts.scaling_factors.tsv
        ├── lengths/
        │    ├── sample_01.length_counts.tsv.gz
        │    └── sample_01.length_settings.json
        └── midpoints/
             ├── sample_01.midpoint_profiles.zarr/
             ├── sample_01.group_index.tsv
             └── sample_01.midpoint_settings.json
```

## Next step

Once you have these files in place, continue with the **Guides** section to find the right workflow for your analysis.
