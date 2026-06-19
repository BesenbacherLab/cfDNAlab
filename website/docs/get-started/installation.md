# Installation

## Environment 

Create the following conda environment to allow building the software.

```bash
conda create -n cfdnalab \
  rust=1.94.0 zstandard perl fontconfig \
  'conda-forge::clang=21.*' 'conda-forge::clangdev=21.*' \
  'conda-forge::libclang=21.*' 'conda-forge::llvmdev=21.*'
conda activate cfdnalab
export LIBCLANG_PATH="$CONDA_PREFIX/lib"
```

## Latest release

```bash
cargo install cfdnalab --locked
cfdna --help
```

## Build from GitHub

```bash
cargo install --git https://github.com/BesenbacherLab/cfDNAlab --locked
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
cargo build --release --locked
```

The binary is available at `target/release/cfdna`.

Verify installation:

```bash
./target/release/cfdna --help
```
