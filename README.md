# cfDNAlab

Incredibly fast command-line tools for analysis of cell-free DNA. Count the fragment coverage, midpoint coverage, or fragment lengths across the whole genome (or selected windows) in mere seconds or minutes. Apply sample-specific GC correction and large-scale genomic smoothing.

Written in rust for *speed*. Built for *paired-end* sequencing data.

To enable a wide set of usecases, the commands are highly flexible with many options and good default settings.

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

| Command                  | Description                                                                                                                                                                                                               |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cfdna fcoverage`        | Count *fragment* coverage per position or aggregated in windows                                                                                                                                                           |
| `cfdna profile-groups`   | Count fragment *midpoint* coverage in fixed-size intervals, collapsed by groups across the genome.<br />E.g. transcription factor binding sites, aggregated per transcription factor.<br />Fast alternative to *Griffin*. |
| `cfdna lengths`<br />    | Count fragment lengths<br />Defined as: `end(reverse) - start(forward)` for inwardly directed pairs only                                                                                                                  |
| **Normalization**        | Precompute these normalization/correction factors to enable their use in the main commands                                                                                                                                |
| `cfdna normalize-genome` | Calculate scaling factors for normalizing/smoothing coverage across the genome                                                                                                                                            |
 

### Common options

 - **Windowing**: Perform the command in genomic windows. Either a single global window (default), windows specified in a BED file, or via a fixed window size. Assign fragments to windows by how they overlap.
 - **Blacklist filtering**: Supply BED files with regions to exclude. The implementation is specific to each tool (filtering of full fragments or just the overlapping positions).

---

## Quick‑start example

Runtime depends on the size of the bam file. The below example ran in < 3min with 12 cores for a ~28x WGS file:

```bash
cfdna lengths \
  --bam sample.bam \                      # coordinate-sorted bam file with paired-end cfDNA
  --output-dir results \                  # where to write files
  --n-threads 12 \                        # use 12 CPU cores (max. one per chromosome)
  --blacklist encode_blacklist.bed        # exclude ENCODE blacklist intervals
```

---

## FAQ

 - How is *fragment* coverage different from the outputs of similar tools like `mosdepth` and `samtools`?
   - `mosdepth` counts the coverage of aligned bases per *read* independently. `fcoverage` instead first collects the paired reads into a fragment and then counts the coverage of the aligned bases and (optionally) the gap between mate reads.  (TODO on samtools!).
 - How do you define a "fragment"?
   - We define the *fragment* as the bases from the start of the forward read till the end of the reverse read: `[start(forward), end(reverse))` (inwardly directed pairs only), as suggested by Wang, H. et al. 2025. Some methods exclude deletions and skipped-regions.
 - Should I order the BAM files differently to allow pairing of reads into fragments?
   - We expect BAM files to be *coordinate-sorted* and indexed.

## TODO

    - Bin chromosomes for higher parallelization where meaningful.
    - Add GC correction tools and implementations.
    - Allow input BED files to be compressed.

## References

 - Wang, H., Mennea, P.D., Chan, Y.K.E. et al. A standardized framework for robust fragmentomic feature extraction from cell-free DNA sequencing data. Genome Biol 26, 141 (2025). https://doi.org/10.1186/s13059-025-03607-5