# cfDNAlab

Incredibly fast command-line tools for analysis of cell-free DNA.

Written in rust for *speed*.

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

| Command                  | Description                                                                    | Output                                                                   |
| ------------------------ | ------------------------------------------------------------------------------ | ------------------------------------------------------------------------ |
| `cfdna fcoverage`        | Count *fragment* coverage per position or aggregated in windows                | TODO                                                                     |
| `cfdna lengths`          | Count fragment lengths<br />Defined as: `end(reverse) - start(forward)`        | `all_length_counts.npy`: Count array<br />`bins.bed`: Window coordinates |
| `cfdna normalize-genome` | Calculate scaling factors for normalizing/smoothing coverage across the genome | TODO                                                                     |
 

### Common options

 - **Windowing**: Perform the command in genomic windows. Either a single global window (default), windows specified in a BED file, or via a fixed window size. Assign fragments to windows by how they overlap.
 - **Blacklist filtering**: Supply BED files with regions to exclude. The implementation is specific to each tool (filtering of full fragments or just the overlapping positions).

---

## Quick‑start example

Runtime depends on the size of the bam file. The below example ran in ~3min with 12 cores for a ~25x WGS file:

```bash
cfdna lengths \
  --bam sample.bam \                      # bam file with paired-end cfDNA
  --output-dir results \                  # where to write files
  --n-threads 12 \                        # use 12 CPU cores (max. one per chromosome)
  --blacklist encode_blacklist.bed        # exclude ENCODE blacklist intervals
```

---

## TODO

    - Bin chromosomes for higher parallelization where meaningful.
    - Add GC correction tools and implementations.