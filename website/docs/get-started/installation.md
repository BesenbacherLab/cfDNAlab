# Installation

## Environment 

Create the following conda environment to allow building the software.

```bash
conda create -n cfdnalab rust=1.94.0 zstandard perl fontconfig conda-forge::llvmdev conda-forge::clangdev
conda activate cfdnalab
```

## Build from GitHub

```bash
cargo install --git https://github.com/ludvigolsen/cfDNAlab --release --features cli,plotters
```

Verify installation:

```bash
cfdna --help
```

## Build from source

```bash
cargo build --release --features cli,plotters
```

The binary is available at `target/release/cfdna`.

Verify installation:

```bash
./target/release/cfdna --help
```
