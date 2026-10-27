#!/usr/bin/env bash

# Install the five coverage benchmark tools on Linux x86_64
# Also install an isolated Python environment for Snyder BAM preparation
# Usage from the repository root:
#   bash benchmarks/coverage/install_tools.sh /faststorage/project/benchmark-tools
#
# Use a new installation directory. Each tool gets its own directory or Conda
# environment. The cfDNAlab build uses this repository's current checkout
# Conda, wget, tar, GNU make and a C compiler must already be available
# For samtools, load the development libraries listed in its INSTALL guide:
# zlib, bzip2, liblzma, libcurl, OpenSSL, ncurses and libdeflate
# Include libdeflate so samtools can use its faster BAM decompression

set -euo pipefail

install_root=${1:?Usage: bash benchmarks/coverage/install_tools.sh /path/to/benchmark-tools}
repository_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
mkdir -p "$install_root"
install_root=$(cd "$install_root" && pwd)

# Enable conda activate inside this Bash script
source "$(conda info --base)/etc/profile.d/conda.sh"

# cfDNAlab
# Guide: https://github.com/BesenbacherLab/cfDNAlab#installation
# Use the documented dependencies, then install the checkout with Cargo
# --root keeps cfdna inside this environment instead of the shared ~/.cargo/bin
# Parentheses keep activation and LIBCLANG_PATH local to this section
(
    conda create --yes --prefix "$install_root/cfdnalab" --channel conda-forge \
        rust=1.94.0 zstandard perl fontconfig \
        'clang=21.*' 'clangdev=21.*' 'libclang=21.*' 'llvmdev=21.*'
    conda activate "$install_root/cfdnalab"
    export LIBCLANG_PATH="$CONDA_PREFIX/lib"
    cargo install --path "$repository_root" --locked --root "$CONDA_PREFIX"
    "$CONDA_PREFIX/bin/cfdna" --version
)

# mosdepth
# Guide: https://github.com/brentp/mosdepth#installation
# The guide recommends the static release binary, which needs no Conda environment
mkdir -p "$install_root/mosdepth-0.3.14"
wget -O "$install_root/mosdepth-0.3.14/mosdepth" \
    https://github.com/brentp/mosdepth/releases/download/v0.3.14/mosdepth
chmod +x "$install_root/mosdepth-0.3.14/mosdepth"
"$install_root/mosdepth-0.3.14/mosdepth" --version

# samtools
# Guide: https://www.htslib.org/download/
# Dependencies: https://github.com/samtools/samtools/blob/1.24/INSTALL
# Build the release tarball with its bundled HTSlib, using the documented
# configure -> make -> make install sequence and a private installation prefix
wget -O "$install_root/samtools-1.24.tar.bz2" \
    https://github.com/samtools/samtools/releases/download/1.24/samtools-1.24.tar.bz2
tar -xjf "$install_root/samtools-1.24.tar.bz2" -C "$install_root"
(
    cd "$install_root/samtools-1.24"
    ./configure --prefix="$install_root/samtools-1.24/install"
    make
    make install
    "$install_root/samtools-1.24/install/bin/samtools" --version
)

# PanDepth
# Guide: https://github.com/HuiyangYu/PanDepth#1-install
# Use the Linux x86_64 release archive shown in the guide
wget -O "$install_root/PanDepth-2.26-Linux-x86_64.tar.gz" \
    https://github.com/HuiyangYu/PanDepth/releases/download/v2.26/PanDepth-2.26-Linux-x86_64.tar.gz
tar -xzf "$install_root/PanDepth-2.26-Linux-x86_64.tar.gz" -C "$install_root"
"$install_root/PanDepth-2.26-Linux-x86_64/pandepth" -h

# FinaleToolkit
# Guide: https://epifluidlab.github.io/FinaleToolkit/documentation/user_guide/installation.html
# Follow the clean Python 3.11 Conda environment + pip installation recipe
# conda-pack is included to make the archive required by benchmark_coverage.sh
# Archive guide: https://conda.github.io/conda-pack/
(
    conda create --yes --prefix "$install_root/finaletoolkit" --channel conda-forge \
        python=3.11 pip conda-pack=0.8.1
    conda activate "$install_root/finaletoolkit"
    # Ignore Python packages and import paths from other environments
    export PYTHONNOUSERSITE=1
    unset PYTHONPATH PYTHONHOME
    python -m pip install finaletoolkit==1.1.0
    finaletoolkit --version
    conda-pack -p "$CONDA_PREFIX" -o "$install_root/finaletoolkit.tar.gz"
)

# Snyder BAM preparation
# Guide: https://pysam.readthedocs.io/en/latest/installation.html
# Use pip in its own environment, following pysam's documented wheel installation
(
    conda create --yes --prefix "$install_root/snyder-preparation" --channel conda-forge \
        python=3.11 pip
    conda activate "$install_root/snyder-preparation"
    python -m pip install pysam==0.24.1
)

# Print the installed paths for copying into the benchmark submission shell
cat <<EOF
export BENCH_CFDNA="$install_root/cfdnalab/bin/cfdna"
export BENCH_MOSDEPTH="$install_root/mosdepth-0.3.14/mosdepth"
export BENCH_SAMTOOLS="$install_root/samtools-1.24/install/bin/samtools"
export BENCH_PANDEPTH="$install_root/PanDepth-2.26-Linux-x86_64/pandepth"
export BENCH_FINALE_ARCHIVE="$install_root/finaletoolkit.tar.gz"
export BENCH_PREPARE_PYTHON="$install_root/snyder-preparation/bin/python"
EOF
