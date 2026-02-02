# cfDNAlab

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
$ conda create -n cfdnalab rust=1.87.0 zstandard perl fontconfig conda-forge::llvmdev conda-forge::clangdev
$ conda activate cfdnalab
```

Compile and install:

```bash
$ cargo install --git https://github.com/ludvigolsen/cfDNAlab
$ cfdna --help
# or clone + build
$ git clone https://github.com/ludvigolsen/cfDNAlab
$ cd cfdnalab && cargo build --release --features cli,plotters
$ target/release/cfdna --help
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
| `cfdna bam-to-bam`                   | Apply our read filters and write GC correction and coverage weight tags to a BAM file                                                                                                                                  |
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
   - `mosdepth` counts the coverage of aligned bases per *read* [TODO: Not that simple]. `fcoverage` instead first collects the paired reads into a fragment and then counts the coverage of the aligned bases and (optionally) the gap between mate reads.  (TODO on samtools!).
 
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
 
 - How did you use LLMs (AI) in this project?
   - OpenAI's codex models were used for pair programming to speed up development and testing. All released code have been designed and validated by us.

---

## Recipes

We aimed for high flexibility so the commands are useful for both established and more novel usecases. This led to the commands having many options. The following recipes (examples) will get you quickly up and running with common cfDNA analyses. Then you can dive into the options when needed. The final example is a full pipeline for running everything (but without the explanations in the separate examples).

### GC correction pipeline

Fragmentomics features are vulnerable to biases from various sample-handling and sequencing processes, such as PCR amplification. `cfDNAlab` commands thus allow the correction of the commonly observed **GC-bias**.

This requires only a few steps: 

1) Calculate the "expected" GC bias in the reference genome assembly (e.g., hg38). This can be **reused** for all samples aligned to that assembly:

```bash

cfdna ref-gc-bias \
  ...

```

2) Calculate the GC-bias correction factors per sample:

```bash

cfdna gc-bias \
  ...
  --ref-gc-dir <path>/

```

3) Provide the correction factors when running the feature extraction commands:


```bash

cfdna fcoverage \
  ...
  --gc-file <path>/gc_bias_correction.npz

```

If you prefer a different/custom GC-bias tool, the feature extraction commands also accepts reading a GC weight (how much a fragment should contribute) from an aux tag in the BAM file:

```bash

cfdna fcoverage \
  ...
  --gc-tag 'GC'

```

### Genomic smoothing pipeline

For some commands, like `cfdna midpoints`, you may want all genomic regions to have a similar contribution to the features. E.g. to reduce the effect of copy number alterations. 

**Simplified**, this can be achieved by calculating the fragment coverage in a kilo/megabase resolution and dividing the contribution of each fragment (`1.0` or the gc-weight) with the coverage value.

**More detailed**, for a more smooth scaling, `cfdna coverage-weights` builds a smoothed normalization map using a sliding window:

A) It splits the genome into "stride-bins" (default: 500kb) and counts the average positional fragment coverage in each bin.

B) It smoothes each bin with a triangular weighting kernel, that weights the coverage of the neighbouring stride-bins by how many overlapping megabins (default: 5Mb) they are part of. 

E.g.: 

Using a megabin-size of `6` and stride size of `2` for demonstrational purposes:

**Stride bins** (fixed along genome, each with an average positional coverage):

`[A] [B] [C] [D] [E] [F] [G] ...`

**Overlapping megabins** (`MB*`) (each covers 3 stride-bins). 
**`W_D`** weights each stride-bin by how many `D`-overlapping megabins they are part of. 
Stride-bin `B` is only part of one of the megabins overlapping `D`, so it has an (unnormalized) weight of 1:

<pre>

<i>MB1</i>: [A][B][C]

MB2:    [B][C][<b>D</b>]

MB3:       [C][<b>D</b>][E]

MB4:          [<b>D</b>][E][F]

<i>MB5</i>:             [E][F][G]

W_D: [0][1][2][3][2][1][0]

</pre>

The further away a stride-bin is from the center stride-bin, the less it contributes to the smoothed average coverage. The MB1 and MB5 bins do not contribute since they don't overlap the center bin.

C) Finally, the values are *inverted* to `1/coverage` to become multiplicative scaling factors (one per stride-bin). A fragment can be scaled by multiplying its contribution (`1.0` or the gc-weight) with the scaling factor of the stride-bin it's located in. 

You can think of this approach as a very fast alternative to e.g. Gaussian smoothing.

This can be achieved with two steps:

1) Calculate the coverage-based scaling factors that leads to such genomic smoothing:

```bash

cfdna coverage-weights \
  ...

```

2) Provide the scaling factors when running the feature extraction commands:

```bash

cfdna midpoints \
  ...
  --scaling-factors <path>/<prefix>.scaling_factors.tsv

```

### Fragment coverage

Fragment coverage measures how many fragments overlap each genomic position. In contrast to many non-cfDNA-tools, we (optionally)count the gap between paired reads along with the aligned bases of the reads. We avoid double counting when reads overlap. When no GC correction or genomic smoothing is applied, each fragment counts `1` in the overlapping (aligned / gap) positions. GC correction and/or genomic smoothing changes this to a weight (floating point).

```bash

cfdna fcoverage \
  --bam sample.bam \                      # coordinate-sorted bam file with paired-end cfDNA
  --output-dir results \                  # where to write files
  --n-threads 12 \                        # use 12 CPU cores (max. one per chromosome)
  --blacklist encode_blacklist.bed        # exclude ENCODE blacklist intervals
  
# Add GC correction and / or genomic smoothing
  --gc ... \
  --scale-genome ...

```

### Fragment lengths

Multiple studies have used fragment lengths (count distributions) to detect cancer [REFS].

NOTE: For fragment lengths, we use the same GC correction for all lengths (based only on GC contents). 

```bash

cfdna lengths \
  --bam sample.bam \                      # coordinate-sorted bam file with paired-end cfDNA
  --output-dir results \                  # where to write files
  --n-threads 12 \                        # use 12 CPU cores (max. one per chromosome)
  --blacklist encode_blacklist.bed        # exclude ENCODE blacklist intervals

```

### Fragment midpoint profiles

```bash

# Optional preparation of intervals to count midpoints in
cfdna prepare_windows ...

# Count midpoints, summed per group and position
cfdna midpoints ...


```

### Everything combined

The below does not show the midpoint profiles. See the separate examples above.

```bash

BAM="..." # ?
OUT="..." # sample specific
BLACKLIST="..." # ?
THREADS=12
MINLENGTH=30
MAXLENGTH=600

# Coverage weights for genomic smoothing
cfdna coverage-weights --bam $BAM --output-dir $OUT/coverage_weights --min-fragment-length $MINLENGTH --max-fragment-length $$MAXLENGTH --n-threads $THREADS 

# GC bias correction matrix
cfdna gc-bias --bam $BAM --output-dir $OUT/gc_bias --min-fragment-length $MINLENGTH --max-fragment-length $$MAXLENGTH --scale-genome $OUT/coverage_weights --n-threads $THREADS 

# Fragment coverage
cfdna fcoverage --bam $BAM --output-dir $OUT/coverage --min-fragment-length $MINLENGTH --max-fragment-length $$MAXLENGTH --gc $OUT/gc_bias --scale-genome $OUT/coverage_weights --blacklist $BLACKLIST --n-threads $THREADS 

# Fragment lengths (global)
cfdna lengths --bam $BAM --output-dir $OUT/lengths_$(MINLENGTH)_$(MAXLENGTH) --min-fragment-length $MINLENGTH --max-fragment-length $$MAXLENGTH --gc $OUT/gc_bias --scale-genome $OUT/coverage_weights --blacklist $BLACKLIST --n-threads $THREADS 

# Fragment lengths in 5Mb bins 
# E.g., to calculate short/long ratios from (100-150bp, 151-220bp)
cfdna lengths --bam $BAM --output-dir $OUT/lengths_per_5mb_100_220 --by-size 5000000 --min-fragment-length 100 --max-fragment-length 220 --gc $OUT/gc_bias --scale-genome $OUT/coverage_weights --blacklist $BLACKLIST --n-threads $THREADS

# Transition probabilities in the first 10bp from each end
cfdna transitions --bam $BAM --output-dir $OUT --orders 1 --frame nearest --positions '..10' --min-fragment-length $MINLENGTH --max-fragment-length $$MAXLENGTH --gc $OUT/gc_bias --scale-genome $OUT/coverage_weights --blacklist $BLACKLIST --n-threads $THREADS

# End motifs
cfdna ends ... # NOT IMPLEMENTED YET (breakpoint motif example also?)

# Nucleosome peaks from windowed protection scores
cfdna wps-peaks --bam $BAM --output-dir $OUT/wps_peaks --min-fragment-length 120 --max-fragment-length 180 --window-size 120 --gc $OUT/gc_bias --scale-genome $OUT/coverage_weights --blacklist $BLACKLIST --n-threads $THREADS

# Statistics on nucleosome peaks per 5Mb
cfdna wps-peaks --bam $BAM --output-dir $OUT/wps_peaks_statistics_per_5mb --by-size 5000000 --per-window stats --min-fragment-length 120 --max-fragment-length 180 --window-size 120 --gc $OUT/gc_bias --scale-genome $OUT/coverage_weights --blacklist $BLACKLIST --n-threads $THREADS

```

---

## TODO

    - Figure out --output-prefix (default to remove prefix?) and use consistently across commmands!
    - Bin chromosomes for higher parallelization where meaningful.
    - Allow input BED files to be compressed.
    - Check / optimize RAM usage in `cfdna coverage-weights`.
    - Find way to handle deletion of temp dirs when command fails (and not in debug mode)
    - Fix double-counting of reads and fragments in stats counters around tile edges (fetch halo-related)

---

## References

 - Wang, H., Mennea, P.D., Chan, Y.K.E. et al. A standardized framework for robust fragmentomic feature extraction from cell-free DNA sequencing data. Genome Biol 26, 141 (2025). https://doi.org/10.1186/s13059-025-03607-5
