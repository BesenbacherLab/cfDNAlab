# cfDNAlab

Ultra-fast command-line tools for analysis of cell-free DNA. Extract the *fragment coverage*, *midpoint coverage*, *fragment lengths*, *fragment k-mers*, or *nucleosome peaks* across the whole genome (or in windows) in mere seconds or minutes. Apply sample-specific GC correction and large-scale genomic smoothing.

Written in rust for *speed*. Built for *paired-end* sequencing data.

To enable a wide set of usecases, the commands are highly flexible with many options and good default settings. See [recipes](#recipes) for examples.

Suggest a tool or feature [here](https://github.com/LudvigOlsen/cfDNAlab/issues/new/choose)!

---

## Installation

### Compile from source

You may need a few dependencies that can be installed as a conda environment with:

```bash
$ conda create -n cfdnalab rust=1.87.0 zstandard perl conda-forge::llvmdev conda-forge::clangdev
$ conda activate cfdnalab
```

Compile and install:

```bash
$ cargo install --git https://github.com/ludvigolsen/cfDNAlab
$ cfdna --help
# or clone + build
$ git clone https://github.com/ludvigolsen/cfDNAlab
$ cd cfdnalab && cargo build --release
$ target/release/cfdna --help
```

---

## Commands
The following commands are currently available:

| Command                  | Description                                                                                                                                                                                                            |
| ------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cfdna fcoverage`        | Count *fragment* coverage per position or aggregated in windows                                                                                                                                                        |
| `cfdna profile-groups`   | Count fragment *midpoint* coverage in fixed-size intervals, collapsed by groups across the genome<br />E.g. transcription factor binding sites, aggregated per transcription factor<br />Fast alternative to *Griffin* |
| `cfdna lengths`          | Count fragment lengths<br />Defined as: `end(reverse) - start(forward)` for inwardly directed pairs only                                                                                                               |
| `cfdna fragment-kmers`   | Count fragment k-mers with highly flexible selection of positions                                                                                                                                                      |
| `cfdna transitions`      | Extract nth-order transition probabilities in specifiable parts of the fragments                                                                                                                                       |
| `cfdna wps-peaks`        | Estimate nucleosome peaks from windowed protection scores                                                                                                                                                              |
| `cfdna wps`              | Calculate windowed protection scores per position (independent of `wps-peaks`)                                                                                                                                         |
| **Normalization**        | Precompute normalization/correction factors to enable their use in the main commands                                                                                                                                   |
| `cfdna coverage-weights` | Calculate scaling factors for normalizing/smoothing coverage across the genome                                                                                                                                         |


### Common options

 - **Windowing**: Perform the command in genomic windows. Either a single global window (default), windows specified in a BED file, or via a fixed window size. Assign fragments to windows by how they overlap.
 
 - **Blacklist filtering**: Supply BED files with regions to exclude. The implementation is specific to each tool (filtering of full fragments or just the overlapping positions).
 
 - **GC bias correction**: Perform GC bias correction by weighting the contribution of each fragment by their GC content.
 
 - **Genomic smoothing**: Scale the contribution of fragments by their coverage in large-scale overlapping bins. This reduces the effect of amplifications and deletions.

---

## FAQ

 - How is *fragment* coverage different from the outputs of similar tools like `mosdepth` and `samtools`?
   - `mosdepth` counts the coverage of aligned bases per *read*, independently. `fcoverage` instead first collects the paired reads into a fragment and then counts the coverage of the aligned bases and (optionally) the gap between mate reads.  (TODO on samtools!).
 
 - How do you define a "fragment"?
   - We define the *fragment* as the bases from the start of the forward read till the end of the reverse read (`[start(forward), end(reverse))`) for *inwardly directed* pairs only (i.e., where `start(forward) <= start(reverse)`), as suggested by Wang, H. et al. 2025. Some methods exclude deletions and skipped-regions.

  Fragment visualization:

  ```text
  Reference 5' >>>>>>>>>>>>>>> 3'
  Fragment     |-------------|
  Forward   5' |>>>>>>>| 3'     
  Reverse        3' |<<<<<<<<| 5' 
  ``` 
 
 - Should I order the BAM files differently to allow pairing of reads into fragments?
   - We expect BAM files to be *coordinate-sorted* and indexed.
 
 - How did you use LLMs (AI) in this project?
   - OpenAI's GPT5 thinking models were used for pair programming to speed up development and testing. All released code have been validated by us.

---

## Recipes

To allow high flexibility, the commands have many options. The following recipes (examples) will get you quickly up and running with common cfDNA analyses. Then you can dive into the options when needed. The final example is a full pipeline for running everything (but without the explanations in the separate examples).

### Correction pipeline

To allow correcting downstream features for GC bias and performing genomic smoothing, it makes sense to first extract the coverage weights and GC correction matrix required per sample. [TODO: Explain genomic smoothing first] NOTE: To run both, you should calculate the GC bias on the already smoothed genome (see below).

Calculate coverage weights:

```bash

cfdna coverage-weights \
  ...

```

Calculate the GC bias correction matrix (for using both, see next example instead):

```bash

cfdna gc-bias \
  ...

```

To combine both transformations, supply the coverage weights to the GC bias estimator:

```bash

cfdna gc-bias \
  ... \
  --scale-genome <>

```

NOTE: When GC bias is calculated on "smoothed" fragments, be conscious about using both transformations together or not in the feature extraction. Usually, we would be consistent with using both together or not throughout the pipeline, but technically you *could* use smoothing only during GC bias estimation to not let large-scale amplifications and deletions affect the GC bias estimation. Whereas the order is always 1) genomic smoothing, 2) GC bias correction (when both are specified of course), the combination of the two is up to the user.

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

$BAM = ... # ?
$OUT = ... # sample specific
$BLACKLIST = ... # ?
$THREADS = 12 #?
$MINLENGTH = 20
$MAXLENGTH = 600

# Coverage weights for genomic smoothing
cfdna coverage-weights --bam $BAM --output-dir $OUT/coverage_weights --min-fragment-length $MINLENGTH --max-fragment-length $$MAXLENGTH --n-threads $THREADS 

# GC bias correction matrix
cfdna gc-bias --bam $BAM --output-dir $OUT/gc_bias --min-fragment-length $MINLENGTH --max-fragment-length $$MAXLENGTH --scale-genome $OUT/coverage_weights --n-threads $THREADS 

# Fragment coverage
cfdna fcoverage --bam $BAM --output-dir $OUT/coverage --min-fragment-length $MINLENGTH --max-fragment-length $$MAXLENGTH --gc $OUT/gc_bias --scale-genome $OUT/coverage_weights --blacklist $BLACKLIST --n-threads $THREADS 

# Fragment lengths (global)
cfdna lengths --bam $BAM --output-dir $OUT/lengths_$MINLENGTH_$MAXLENGTH --min-fragment-length $MINLENGTH --max-fragment-length $$MAXLENGTH --gc $OUT/gc_bias --scale-genome $OUT/coverage_weights --blacklist $BLACKLIST --n-threads $THREADS 

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

    - Bin chromosomes for higher parallelization where meaningful.
    - Add GC correction tools and implementations.
    - Allow input BED files to be compressed.
    - Check / optimize RAM usage in `cfdna coverage-weights`.

---

## References

 - Wang, H., Mennea, P.D., Chan, Y.K.E. et al. A standardized framework for robust fragmentomic feature extraction from cell-free DNA sequencing data. Genome Biol 26, 141 (2025). https://doi.org/10.1186/s13059-025-03607-5
