#!/usr/bin/env bash

# Install all benchmark tools on Linux x86_64, with a section per tool
# Also install isolated environments for BAM preparation and output comparisons
# Usage from the repository root:
#   bash benchmarks/install_tools.sh /faststorage/project/benchmark-tools
#
# Use a new installation directory. Conda environments go in envs/<tool>
# Standalone tools, source checkouts and packed archives stay at the root
# The cfDNAlab build uses this repository's current checkout
# Tool versions and installation routes were checked on 2026-10-29
# Pin those versions here so the benchmark can be reproduced later
# Conda, Apptainer, Git, wget, tar, GNU make and C/C++ compilers must be available
# Only cfDNAlab itself needs a local native build
# See benchmarks/README.md for the installation audit and deliberate adaptations

set -euo pipefail

install_root=${1:?Usage: bash benchmarks/install_tools.sh /path/to/benchmark-tools}
repository_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
mkdir -p "$install_root"
install_root=$(cd "$install_root" && pwd)
mkdir -p "$install_root/envs"

# Python normally also searches ~/.local for packages installed with pip --user
# Ignore those packages so each tool imports its own environment's versions
# PYTHONPATH adds outside import directories. PYTHONHOME can select another
# Python installation. Clear both for the same reason
# These changes affect this script and its children, not the calling shell
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME

# Enable conda activate inside this Bash script
source "$(conda info --base)/etc/profile.d/conda.sh"

# Both Bioconda and FrEIA recommend strict channel priority to avoid incompatible
# packages from lower-priority channels. Apply it here without editing ~/.condarc
export CONDA_CHANNEL_PRIORITY=strict

# Archive preparation, used outside all benchmark timers
# Guide: https://conda.github.io/conda-pack/
# Install conda-pack once, separately from the tools whose environments it packs
# Complete archives let the runners relocate tools and libraries to local scratch
# The runners use an already-staged Python to execute each conda-unpack script
conda create --yes --prefix "$install_root/envs/packaging" --channel conda-forge \
    python=3.11 conda-pack=0.8.1

# cfDNAlab
# Guide: https://github.com/BesenbacherLab/cfDNAlab#installation
# Use the documented dependencies, then install the checkout with Cargo
# --root keeps cfdna inside this environment instead of the shared ~/.cargo/bin
# Parentheses keep activation and LIBCLANG_PATH local to this section
(
    conda create --yes --prefix "$install_root/envs/cfdnalab" --channel conda-forge \
        rust=1.94.0 zstandard perl fontconfig \
        'clang=21.*' 'clangdev=21.*' 'libclang=21.*' 'llvmdev=21.*'
    conda activate "$install_root/envs/cfdnalab"
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
# Guide: https://bioconda.github.io/recipes/samtools/README.html
# Install the pinned Bioconda package in its own environment
# Pack the complete environment so its HTSlib and other libraries run locally
(
    conda create --yes --prefix "$install_root/envs/samtools" \
        --channel conda-forge --channel bioconda samtools=1.25
    conda activate "$install_root/envs/samtools"
    samtools --version-only
    "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/samtools.tar.gz"
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
# Conda alternative: https://bioconda.github.io/recipes/finaletoolkit/README.html
# Use the same release through Bioconda so compiled dependencies arrive prebuilt
# Coverage, end-motif and length benchmarks use this same environment archive
(
    conda create --yes --prefix "$install_root/envs/finaletoolkit" \
        --channel conda-forge --channel bioconda python=3.11 finaletoolkit=1.1.0
    conda activate "$install_root/envs/finaletoolkit"
    finaletoolkit --version
    "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/finaletoolkit.tar.gz"
)

# cfDNAPro
# Guide: https://github.com/hw538/cfDNAPro#installation
# Use the authors' recommended Apptainer image, including its R dependencies
# and built-in reference packages. The benchmark copies this image locally
# The guide recommends an immutable digest. This identifies their 1.7.4 image
# published from commit ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c
apptainer pull "$install_root/cfdnapro-1.7.4.sif" \
    docker://ghcr.io/hw538/cfdnapro@sha256:26a28821a272856e77fdf3d13d902ce8ca5c1e60475b1d4f2f076ef2ac3e2b45
apptainer exec --cleanenv "$install_root/cfdnapro-1.7.4.sif" \
    Rscript --vanilla -e 'cat(as.character(packageVersion("cfDNAPro")), "\n")'

# Picard
# Guide: https://broadinstitute.github.io/picard/
# Conda: https://bioconda.github.io/recipes/picard-slim/README.html
# picard-slim contains the official JAR and Java, without the R dependency
# Our CollectInsertSizeMetrics call omits the optional PDF, so it does not need R
(
    conda create --yes --prefix "$install_root/envs/picard" \
        --channel conda-forge --channel bioconda picard-slim=3.5.0 openjdk=17
    conda activate "$install_root/envs/picard"
    # The packaged launcher honors JAVA_HOME. Keep it on this environment's Java
    export JAVA_HOME="$CONDA_PREFIX"
    picard CollectInsertSizeMetrics --version
    "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/picard.tar.gz"
)

# FrEIA
# Guide: https://github.com/mouliere-lab/FrEIA#installation-guide
# There are no published releases. Pin the latest source commit instead
# Use the authors' complete environment file without changing dependency pins
# Only the end-extraction and abundance scripts are needed for an existing BAM
git clone https://github.com/mouliere-lab/FrEIA.git "$install_root/FrEIA-source"
git -C "$install_root/FrEIA-source" checkout 2c32fb9d9fab774f4ae57931df7fd7d75654d343
(
    conda env create --prefix "$install_root/envs/freia" \
        --file "$install_root/FrEIA-source/workflow/envs/FrEIA_env.yaml"
    conda activate "$install_root/envs/freia"
    # Keep the unmodified scripts inside the archive so they also run locally
    mkdir -p "$CONDA_PREFIX/share/freia"
    cp "$install_root/FrEIA-source/workflow/scripts/1_FrEIA/1_extract_fragment_ends.py" \
        "$install_root/FrEIA-source/workflow/scripts/1_FrEIA/2_data_transformation.py" \
        "$install_root/FrEIA-source/workflow/scripts/1_FrEIA/FrEIA_tools.py" \
        "$CONDA_PREFIX/share/freia/"
    git -C "$install_root/FrEIA-source" rev-parse HEAD > "$CONDA_PREFIX/share/freia/version.txt"
    "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/freia.tar.gz"
)

# Griffin
# Guide: https://github.com/adoebley/Griffin/wiki
# v0.2.0 is the newest published version, marked as a prerelease upstream
# Follow the wiki's Python and numerical-library versions in an isolated env
# Install pybedtools and pysam through Bioconda at the wiki's versions to avoid
# compiling them locally. PyYAML is also imported by the native scripts
# Omit Snakemake, plotting and samtools because these direct script calls use none
git clone --branch v0.2.0 --depth 1 https://github.com/adoebley/Griffin.git "$install_root/Griffin-source"
(
    conda create --yes --prefix "$install_root/envs/griffin" \
        --channel conda-forge --channel bioconda \
        python=3.7.4 nomkl numpy=1.21.2 pandas=1.3.2 scipy=1.7.1 \
        pybigwig=0.3.17 pyyaml bedtools=2.29.2 pybedtools=0.8.0 pysam=0.19.0
    conda activate "$install_root/envs/griffin"
    # Put unmodified source scripts in the archive for node-local execution
    mkdir -p "$CONDA_PREFIX/share/griffin"
    cp "$install_root/Griffin-source/scripts/griffin_coverage.py" \
        "$install_root/Griffin-source/scripts/griffin_merge_sites.py" \
        "$install_root/Griffin-source/scripts/griffin_functions.py" \
        "$CONDA_PREFIX/share/griffin/"
    printf 'v0.2.0 %s\n' "$(git -C "$install_root/Griffin-source" rev-parse HEAD)" \
        > "$CONDA_PREFIX/share/griffin/version.txt"
    "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/griffin.tar.gz"
)

# cfDNAlab Python reader, used only after timing to compare saved outputs
# Guide: ../py-cfdnalab/README.md for installing the Python loader package
# Keep these dependencies separate from the timed tools
(
    conda create --yes --prefix "$install_root/envs/comparison" --channel conda-forge \
        python=3.11 pip fastparquet
    conda activate "$install_root/envs/comparison"
    # fastparquet reads FrEIA's native output. Install our loader from this checkout
    python -m pip install "$repository_root/py-cfdnalab"
)

# Snyder BAM preparation
# Guide: https://pysam.readthedocs.io/en/latest/installation.html
# The pysam guide recommends Conda/Bioconda, which supplies compiled dependencies
(
    conda create --yes --prefix "$install_root/envs/snyder-preparation" \
        --channel conda-forge --channel bioconda python=3.11 pysam=0.24.1
    conda activate "$install_root/envs/snyder-preparation"
    python -c 'import pysam; print(pysam.__version__)'
)

# Print the installed paths for copying into the benchmark submission shell
cat <<EOF
export BENCH_CFDNA="$install_root/envs/cfdnalab/bin/cfdna"
export BENCH_MOSDEPTH="$install_root/mosdepth-0.3.14/mosdepth"
export BENCH_SAMTOOLS_ARCHIVE="$install_root/samtools.tar.gz"
export BENCH_PANDEPTH="$install_root/PanDepth-2.26-Linux-x86_64/pandepth"
export BENCH_FINALE_ARCHIVE="$install_root/finaletoolkit.tar.gz"
export BENCH_CFDNAPRO_IMAGE="$install_root/cfdnapro-1.7.4.sif"
export BENCH_PICARD_ARCHIVE="$install_root/picard.tar.gz"
export BENCH_FREIA_ARCHIVE="$install_root/freia.tar.gz"
export BENCH_GRIFFIN_ARCHIVE="$install_root/griffin.tar.gz"
export BENCH_COMPARE_PYTHON="$install_root/envs/comparison/bin/python"
export BENCH_PREPARE_PYTHON="$install_root/envs/snyder-preparation/bin/python"
EOF
