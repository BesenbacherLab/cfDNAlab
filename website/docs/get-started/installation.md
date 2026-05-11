# Installation

## Environment 

Create the following conda environment to allow building the software.

```bash
conda create -n cfdnalab rust=1.94.0 zstandard perl fontconfig conda-forge::llvmdev conda-forge::clangdev
conda activate cfdnalab
```

## Build from GitHub

```bash
cargo install --git https://github.com/BesenbacherLab/cfDNAlab
```

Verify installation:

```bash
cfdna --help
```

## Build from source

```bash
# Once downloaded, enter the directory
cd cfDNAlab
# Then build it as so:
cargo build --release
```

The binary is available at `target/release/cfdna`.

Verify installation:

```bash
./target/release/cfdna --help
```
