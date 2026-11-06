#!/usr/bin/env bash

# This script installs benchmark tools on Linux x86_64, with a section per tool
# It also installs separate environments for BAM preparation and output comparisons
# Usage from the repository root:
#   bash benchmarks/install_tools.sh /faststorage/project/benchmark-tools
#
# Rerunning resumes completed versions. Use --force TOOL to reinstall
# Conda environments go in envs/<tool>
# Standalone tools, source checkouts and packed archives stay at the root
# The cfDNAlab build uses this repository's current checkout
# Tool versions and installation routes were checked on 2026-10-29
# Those versions are pinned here so the benchmark can be reproduced later
# Conda, Apptainer, Git, wget, tar, GNU make and C/C++ compilers must be available
# Only cfDNAlab itself needs a local native build
# See benchmarks/README.md for the installation audit and deliberate adaptations
#
# The installer reads the requested tool list and checks manifest.tsv for completed versions
# It installs or replaces selected tools as needed, then prints paths for the job scripts
# The functions below define the installation commands and dependencies for each tool

set -euo pipefail
if ((BASH_VERSINFO[0] < 4)); then echo 'This installer needs Bash 4 or later.' >&2; exit 1; fi

install_root=${1:?Usage: bash benchmarks/install_tools.sh ROOT [--only TOOLS] [--force TOOLS] [--mark-outdated TOOLS]}
shift
repository_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$repository_root/benchmarks/common.sh"
only_tools=all
force_tools=
outdated_tools=
while (($#)); do
    case "$1" in
        --only|--force|--mark-outdated)
            option=$1
            selection=${2:?Supply a comma-separated tool list}
            case "$option" in
                --only) only_tools=$selection ;;
                --force) force_tools=$selection ;;
                --mark-outdated) outdated_tools=$selection ;;
            esac
            shift 2
            ;;
        *) echo "Unknown installer option: $1" >&2; exit 1 ;;
    esac
done
mkdir -p "$install_root"
install_root=$(cd "$install_root" && pwd -P)
if [[ "$install_root" == / ]]; then echo 'Choose a dedicated installation directory.' >&2; exit 1; fi
mkdir -p "$install_root/envs" "$install_root/logs" "$install_root/records"
manifest="$install_root/manifest.tsv"
if [[ ! -f "$manifest" ]]; then printf 'tool\tversion\tstatus\tartifacts\n' > "$manifest"; fi

# The tables below specify each tool's version, expected files and BENCH_* paths
# All file paths are relative to ROOT. A matching version with complete files is reused
# User control: source/dependency changes that keep the same version require an explicit
# reinstallation request with --force or --mark-outdated
# all_tools lists installation order, starting with the Conda packaging helper
declare -A versions artifacts tool_exports
versions[packaging]=0.9.2
versions[cfdnalab]=$(sed -n 's/^version = "\([^"]*\)"/\1/p' "$repository_root/Cargo.toml" | head -n 1)
versions[mosdepth]=0.3.14
versions[samtools]=1.25
versions[pandepth]=2.26
versions[finaletoolkit]=1.1.0
versions[cfdnapro]=1.7.4
versions[picard]=3.5.0
versions[freia]=2c32fb9d9fab774f4ae57931df7fd7d75654d343
versions[griffin]=v0.2.0
versions[comparison]=checkout
versions[snyder-preparation]=0.24.1
artifacts[packaging]='envs/packaging/bin/conda-pack'
artifacts[cfdnalab]='envs/cfdnalab/bin/cfdna cfdnalab-build.tsv'
artifacts[mosdepth]='mosdepth-0.3.14/mosdepth'
artifacts[samtools]='envs/samtools/bin/samtools samtools.tar.gz'
artifacts[pandepth]='PanDepth-2.26-Linux-x86_64/pandepth'
artifacts[finaletoolkit]='envs/finaletoolkit/bin/finaletoolkit finaletoolkit.tar.gz'
artifacts[cfdnapro]='cfdnapro-1.7.4.sif'
artifacts[picard]='envs/picard/bin/picard picard.tar.gz'
artifacts[freia]='envs/freia/share/freia/version.txt freia.tar.gz'
artifacts[griffin]='envs/griffin/share/griffin/version.txt envs/griffin/share/griffin/CTCF.hg38.1000.txt griffin.tar.gz'
artifacts[comparison]='envs/comparison/bin/python'
artifacts[snyder-preparation]='envs/snyder-preparation/bin/python'
tool_exports[cfdnalab]='BENCH_CFDNA=envs/cfdnalab/bin/cfdna'
tool_exports[mosdepth]='BENCH_MOSDEPTH=mosdepth-0.3.14/mosdepth'
tool_exports[samtools]='BENCH_SAMTOOLS_ARCHIVE=samtools.tar.gz'
tool_exports[pandepth]='BENCH_PANDEPTH=PanDepth-2.26-Linux-x86_64/pandepth'
tool_exports[finaletoolkit]='BENCH_FINALE_ARCHIVE=finaletoolkit.tar.gz'
tool_exports[cfdnapro]='BENCH_CFDNAPRO_IMAGE=cfdnapro-1.7.4.sif'
tool_exports[picard]='BENCH_PICARD_ARCHIVE=picard.tar.gz'
tool_exports[freia]='BENCH_FREIA_ARCHIVE=freia.tar.gz'
tool_exports[griffin]='BENCH_GRIFFIN_ARCHIVE=griffin.tar.gz'
tool_exports[comparison]='BENCH_COMPARE_PYTHON=envs/comparison/bin/python'
tool_exports[snyder-preparation]='BENCH_PREPARE_PYTHON=envs/snyder-preparation/bin/python'
all_tools=(packaging cfdnalab mosdepth samtools pandepth finaletoolkit cfdnapro picard freia griffin comparison snyder-preparation)
for selection in "$only_tools" "$force_tools" "$outdated_tools"; do
    for tool in ${selection//,/ }; do
        if [[ "$tool" != all && -z ${versions[$tool]:-} ]]; then echo "Unknown tool: $tool" >&2; exit 1; fi
    done
done

record_status() {
    # This function updates one tool's record, preserving the records for other tools
    # It writes a complete replacement file, then renames it over manifest.tsv
    # If writing is interrupted, the previous installation records remain readable
    # Tools are installed sequentially, so only this script writes this file
    local tool=$1 status=$2
    awk -F '\t' -v tool="$tool" '$1 != tool' "$manifest" > "$manifest.tmp"
    printf '%s\t%s\t%s\t%s\n' "$tool" "${versions[$tool]}" "$status" "${artifacts[$tool]}" >> "$manifest.tmp"
    mv "$manifest.tmp" "$manifest"
}

completed_tool() {
    # This check permits reuse only if the version matches, installation finished
    # and all required files still exist. Otherwise the main loop installs it again
    # This does not inspect changed source or dependencies. Those require an explicit
    # reinstallation request when the recorded tool version stays the same
    local tool=$1 artifact
    [[ $(awk -F '\t' -v tool="$tool" -v version="${versions[$tool]}" '$1 == tool && $2 == version {print $3}' "$manifest") == complete ]] || return 1
    for artifact in ${artifacts[$tool]}; do [[ -e "$install_root/$artifact" ]] || return 1; done
}

for tool in "${all_tools[@]}"; do
    if [[ ",$outdated_tools," == *,all,* || ",$outdated_tools," == *,$tool,* ]]; then record_status "$tool" outdated; fi
done
if [[ -n "$outdated_tools" ]]; then
    # --mark-outdated only updates the installation record, then exits
    # A later installer run replaces the tool. No Conda commands or downloads occur here
    echo "Marked outdated: $outdated_tools"
    exit 0
fi

# Python imports and Conda channels stay inside the requested environments
export PYTHONNOUSERSITE=1 CONDA_CHANNEL_PRIORITY=strict
unset PYTHONPATH PYTHONHOME
source "$(conda info --base)/etc/profile.d/conda.sh"

install_packaging() {
    # Archive preparation, used outside all benchmark timers
    # Guide: https://conda.github.io/conda-pack/
    # This function installs conda-pack once, apart from the tools it packages
    # Complete archives let the runners relocate tools and libraries to local scratch
    # The runners use BENCH_SETUP_PYTHON for each untimed conda-unpack script
    conda create --yes --override-channels --prefix "$install_root/envs/packaging" --channel conda-forge \
        python=3.11 conda-pack=0.9.2
    conda list --prefix "$install_root/envs/packaging" --explicit > "$install_root/records/packaging.explicit.txt"
}

install_cfdnalab() {
    # cfDNAlab
    # Guide: https://github.com/BesenbacherLab/cfDNAlab#installation
    # The documented dependencies are installed before Cargo builds this checkout
    # --root keeps cfdna inside this environment instead of the shared ~/.cargo/bin
    # Parentheses keep activation and LIBCLANG_PATH local to this section
    (
        conda create --yes --override-channels --prefix "$install_root/envs/cfdnalab" --channel conda-forge \
            rust=1.94.0 zstandard perl fontconfig \
            'clang=21.*' 'clangdev=21.*' 'libclang=21.*' 'llvmdev=21.*'
        conda activate "$install_root/envs/cfdnalab"
        export LIBCLANG_PATH="$CONDA_PREFIX/lib"
        export CARGO_TARGET_DIR="$install_root/build/cfdnalab"
        cargo install --path "$repository_root" --locked --root "$CONDA_PREFIX"
        "$CONDA_PREFIX/bin/cfdna" --version
    )
    {
        printf "commit\t%s\n" "$(git -C "$repository_root" rev-parse HEAD)"
        if [[ -n $(git -C "$repository_root" status --porcelain) ]]; then dirty=true; else dirty=false; fi
        printf "dirty\t%s\n" "$dirty"
        cargo_version=$("$install_root/envs/cfdnalab/bin/rustc" --version)
        printf "rustc\t%s\n" "$cargo_version"
    } > "$install_root/cfdnalab-build.tsv"
    conda list --prefix "$install_root/envs/cfdnalab" --explicit > "$install_root/records/cfdnalab.explicit.txt"
}

install_mosdepth() {
    # mosdepth
    # Guide: https://github.com/brentp/mosdepth#installation
    # The guide recommends the static release binary, which needs no Conda environment
    mkdir -p "$install_root/mosdepth-0.3.14"
    wget -O "$install_root/mosdepth-0.3.14/mosdepth" \
        https://github.com/brentp/mosdepth/releases/download/v0.3.14/mosdepth
    chmod +x "$install_root/mosdepth-0.3.14/mosdepth"
    "$install_root/mosdepth-0.3.14/mosdepth" --version
}

install_samtools() {
    # samtools
    # Guide: https://bioconda.github.io/recipes/samtools/README.html
    # The pinned Bioconda package is installed in its own environment
    # The archive includes HTSlib and other libraries needed on the compute node
    (
        conda create --yes --override-channels --prefix "$install_root/envs/samtools" \
            --channel conda-forge --channel bioconda samtools=1.25
        conda activate "$install_root/envs/samtools"
        samtools --version-only
        "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/samtools.tar.gz"
    )
    conda list --prefix "$install_root/envs/samtools" --explicit > "$install_root/records/samtools.explicit.txt"
}

install_pandepth() {
    # PanDepth
    # Guide: https://github.com/HuiyangYu/PanDepth#1-install
    # This function downloads the Linux x86_64 release archive shown in the guide
    wget -O "$install_root/PanDepth-2.26-Linux-x86_64.tar.gz" \
        https://github.com/HuiyangYu/PanDepth/releases/download/v2.26/PanDepth-2.26-Linux-x86_64.tar.gz
    tar -xzf "$install_root/PanDepth-2.26-Linux-x86_64.tar.gz" -C "$install_root"
    "$install_root/PanDepth-2.26-Linux-x86_64/pandepth" -h
}

install_finaletoolkit() {
    # FinaleToolkit
    # Guide: https://epifluidlab.github.io/FinaleToolkit/documentation/user_guide/installation.html
    # Conda alternative: https://bioconda.github.io/recipes/finaletoolkit/README.html
    # Bioconda supplies the selected release and its prebuilt compiled dependencies
    # The required loess version comes from PyPI, after Conda has finished installing
    # Coverage, end-motif and length benchmarks use this same environment archive
    (
        conda create --yes --override-channels --prefix "$install_root/envs/finaletoolkit" \
            --channel conda-forge --channel bioconda python=3.11 pip finaletoolkit=1.1.0 \
            'numpy>=2' 'pandas>=2.3.0' \
            pysam=0.24.0 pybigwig=0.3.25 py2bit=0.3.3 numba=0.63.1 \
            scipy=1.17.0 statsmodels=0.14.6 tqdm=4.67.3 matplotlib=3.10.8
        conda activate "$install_root/envs/finaletoolkit"
        # Pip supplies loess and any missing dependencies inside this environment
        # The check rejects inconsistent requirements before the environment is packed
        python -m pip install loess==2.1.2
        python -m pip check
        python -m pip freeze > "$install_root/records/finaletoolkit.pip.txt"
        finaletoolkit --version
        "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/finaletoolkit.tar.gz"
    )
    conda list --prefix "$install_root/envs/finaletoolkit" --explicit > "$install_root/records/finaletoolkit.explicit.txt"
}

install_cfdnapro() {
    # cfDNAPro
    # Guide: https://github.com/hw538/cfDNAPro#installation
    # This function downloads the authors' recommended Apptainer image, including
    # R dependencies and built-in references. The benchmark copies this image locally
    # The guide recommends an immutable digest. This identifies their 1.7.4 image
    # published from commit ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c
    apptainer pull "$install_root/cfdnapro-1.7.4.sif" \
        docker://ghcr.io/hw538/cfdnapro@sha256:26a28821a272856e77fdf3d13d902ce8ca5c1e60475b1d4f2f076ef2ac3e2b45
    apptainer exec --cleanenv "$install_root/cfdnapro-1.7.4.sif" \
        Rscript --vanilla -e 'cat(as.character(packageVersion("cfDNAPro")), "\n")'
    apptainer inspect --json "$install_root/cfdnapro-1.7.4.sif" > "$install_root/records/cfdnapro-image.json"
    apptainer exec --cleanenv "$install_root/cfdnapro-1.7.4.sif" Rscript --vanilla -e 'for (package in c("cfDNAPro", "GenomicAlignments", "Rsamtools")) cat(package, as.character(packageVersion(package)), "\n")' > "$install_root/records/cfdnapro-packages.txt"
}

install_picard() {
    # Picard
    # Guide: https://broadinstitute.github.io/picard/
    # Conda: https://bioconda.github.io/recipes/picard-slim/README.html
    # picard-slim contains the official JAR and Java, without the R dependency
    # Our CollectInsertSizeMetrics call omits the optional PDF, so it does not need R
    (
        conda create --yes --override-channels --prefix "$install_root/envs/picard" \
            --channel conda-forge --channel bioconda picard-slim=3.5.0 openjdk=17
        conda activate "$install_root/envs/picard"
        # JAVA_HOME selects the Java in this environment for Picard's version check
        export JAVA_HOME="$CONDA_PREFIX"
        benchmark_picard_version "$CONDA_PREFIX/bin/picard"
        printf "\n"
        "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/picard.tar.gz"
    )
    conda list --prefix "$install_root/envs/picard" --explicit > "$install_root/records/picard.explicit.txt"
}

install_freia() {
    # FrEIA
    # Guide: https://github.com/mouliere-lab/FrEIA#installation-guide
    # This installation uses a pinned source commit rather than a release archive
    # The authors' dependency pins are retained, with outside Conda channels excluded
    # Only the end-extraction and abundance scripts are needed for an existing BAM
    git clone https://github.com/mouliere-lab/FrEIA.git "$install_root/FrEIA-source"
    git -C "$install_root/FrEIA-source" checkout 2c32fb9d9fab774f4ae57931df7fd7d75654d343
    (
        awk '/^dependencies:/ {print "  - nodefaults"} {print}' "$install_root/FrEIA-source/workflow/envs/FrEIA_env.yaml" \
            > "$install_root/freia-environment.yaml"
        conda env create --yes --prefix "$install_root/envs/freia" \
            --file "$install_root/freia-environment.yaml"
        conda activate "$install_root/envs/freia"
        # The archive includes unmodified scripts for execution on the compute node
        mkdir -p "$CONDA_PREFIX/share/freia"
        cp "$install_root/FrEIA-source/workflow/scripts/1_FrEIA/1_extract_fragment_ends.py" \
            "$install_root/FrEIA-source/workflow/scripts/1_FrEIA/2_data_transformation.py" \
            "$install_root/FrEIA-source/workflow/scripts/1_FrEIA/FrEIA_tools.py" \
            "$CONDA_PREFIX/share/freia/"
        git -C "$install_root/FrEIA-source" rev-parse HEAD > "$CONDA_PREFIX/share/freia/version.txt"
        "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/freia.tar.gz"
    )
    conda list --prefix "$install_root/envs/freia" --explicit > "$install_root/records/freia.explicit.txt"
    "$install_root/envs/freia/bin/python" -m pip freeze > "$install_root/records/freia.pip.txt"
}

install_griffin() {
    # Griffin
    # Guide: https://github.com/adoebley/Griffin/wiki
    # v0.2.0 is the newest published version, marked as a prerelease upstream
    # This environment uses the wiki's libraries with Python 3.7.12 instead of 3.7.4
    # Bioconda supplies the wiki's pybedtools and pysam versions without a local build
    # PyYAML is included because the scripts import it. Snakemake, plotting and
    # samtools are omitted because these direct profiling calls do not use them
    git clone --branch v0.2.0 --depth 1 https://github.com/adoebley/Griffin.git "$install_root/Griffin-source"
    (
        conda create --yes --override-channels --prefix "$install_root/envs/griffin" \
            --channel conda-forge --channel bioconda \
            python=3.7.12 nomkl numpy=1.21.2 pandas=1.3.2 scipy=1.7.1 \
            pybigwig=0.3.17 pyyaml=5.3.1 bedtools=2.29.2 pybedtools=0.8.0 pysam=0.19.0
        conda activate "$install_root/envs/griffin"
        # The archive includes unmodified scripts for execution on the compute node
        mkdir -p "$CONDA_PREFIX/share/griffin"
        cp "$install_root/Griffin-source/scripts/griffin_coverage.py" \
            "$install_root/Griffin-source/scripts/griffin_merge_sites.py" \
            "$install_root/Griffin-source/scripts/griffin_functions.py" \
            "$CONDA_PREFIX/share/griffin/"
        # The release's 1000-site demo is included for the short, untimed warmup
        cp "$install_root/Griffin-source/demo/griffin_nucleosome_profiling_demo_files/sites/CTCF.hg38.1000.txt" \
            "$CONDA_PREFIX/share/griffin/"
        printf 'v0.2.0 %s\n' "$(git -C "$install_root/Griffin-source" rev-parse HEAD)" \
            > "$CONDA_PREFIX/share/griffin/version.txt"
        "$install_root/envs/packaging/bin/conda-pack" -p "$CONDA_PREFIX" -o "$install_root/griffin.tar.gz"
    )
    conda list --prefix "$install_root/envs/griffin" --explicit > "$install_root/records/griffin.explicit.txt"
}

install_comparison() {
    # cfDNAlab Python reader and zstd, used only after timing to compare outputs
    # Guide: ../py-cfdnalab/README.md for installing the Python loader package
    # These dependencies have their own environment, separate from the timed tools
    (
        conda create --yes --override-channels --prefix "$install_root/envs/comparison" --channel conda-forge \
            python=3.11 pip zstd pyarrow
        conda activate "$install_root/envs/comparison"
        # pyarrow reads FrEIA's intact left-only native Parquet abundance table
        # The Python output reader is installed from this checkout
        python -m pip install "$repository_root/py-cfdnalab"
    )
    conda list --prefix "$install_root/envs/comparison" --explicit > "$install_root/records/comparison.explicit.txt"
    "$install_root/envs/comparison/bin/python" -m pip freeze > "$install_root/records/comparison.pip.txt"
}

install_snyder_preparation() {
    # Snyder BAM preparation
    # Guide: https://pysam.readthedocs.io/en/latest/installation.html
    # The pysam guide recommends Conda/Bioconda, which supplies compiled dependencies
    (
        conda create --yes --override-channels --prefix "$install_root/envs/snyder-preparation" \
            --channel conda-forge --channel bioconda python=3.11 pysam=0.24.1
        conda activate "$install_root/envs/snyder-preparation"
        python -c 'import pysam; print(pysam.__version__)'
    )
    conda list --prefix "$install_root/envs/snyder-preparation" --explicit > "$install_root/records/snyder-preparation.explicit.txt"
}

# Conda tools need an archive containing their programs and libraries so the job
# can copy them to a compute node. conda-pack is installed first when a selected
# tool needs such an archive. It is the only tool added automatically to --only
need_packaging=false
for tool in samtools finaletoolkit picard freia griffin; do
    if [[ "$only_tools" == all || ",$only_tools," == *,$tool,* ]]; then need_packaging=true; fi
done
# The error message uses this name to locate the tool's installation log
current_tool=""
trap 'status=$?; if ((status != 0)); then echo "Installation failed for $current_tool. See $install_root/logs/$current_tool.log" >&2; tail -n 30 "$install_root/logs/$current_tool.log" >&2; fi' EXIT
for tool in "${all_tools[@]}"; do
    # --only chooses the tools to consider. --force reinstalls among those tools
    if [[ "$only_tools" != all && ",$only_tools," != *,$tool,* && !( "$tool" == packaging && "$need_packaging" == true ) ]]; then continue; fi
    if [[ ",$force_tools," != *,all,* && ",$force_tools," != *,$tool,* ]] && completed_tool "$tool"; then
        echo "Keeping $tool ${versions[$tool]}"
        continue
    fi
    current_tool=$tool
    echo "Installing $tool ${versions[$tool]} (log: $install_root/logs/$tool.log)"
    # This tool is marked unfinished before its files are replaced. If installation
    # stops, the next run retries it while reusing tools that completed successfully
    record_status "$tool" outdated
    # Maintenance requirement: this block must stay outside an if/|| error check
    # That lets any installation failure stop the script. Such a check would make
    # Bash ignore automatic stopping and could mark a failed installation complete
    (
        # Conda activation here affects this installation, not the next tool's setup
        # Only this tool's known files are removed. cfDNAlab's compiled dependency
        # files stay in its separate Cargo build folder to avoid unnecessary rebuilds
        cleanup_paths=("envs/$tool")
        case "$tool" in
            samtools|finaletoolkit|picard|freia|griffin) cleanup_paths+=("$tool.tar.gz") ;;
            mosdepth) cleanup_paths+=(mosdepth-0.3.14) ;;
            pandepth) cleanup_paths+=(PanDepth-2.26-Linux-x86_64 PanDepth-2.26-Linux-x86_64.tar.gz) ;;
            cfdnapro) cleanup_paths+=(cfdnapro-1.7.4.sif) ;;
        esac
        if [[ "$tool" == freia ]]; then cleanup_paths+=(FrEIA-source); fi
        if [[ "$tool" == griffin ]]; then cleanup_paths+=(Griffin-source); fi
        for artifact in "${cleanup_paths[@]}"; do
            if [[ -L "$install_root/$artifact" || -L "$install_root/envs" ]]; then
                echo "Refusing cleanup through a symlink: $install_root/$artifact" >&2
                exit 1
            fi
            rm -rf "$install_root/$artifact"
        done
        "install_${tool//-/_}"
    ) > "$install_root/logs/$tool.log" 2>&1
    # Success is recorded only after the expected programs/archives are found
    # A command that exits successfully without creating its output is not enough
    for artifact in ${artifacts[$tool]}; do
        if [[ ! -e "$install_root/$artifact" ]]; then echo "Missing artifact: $artifact" >&2; exit 1; fi
    done
    record_status "$tool" complete
done
# Installation is finished, so the tool name is cleared before paths are printed
current_tool=""

# This loop prints the paths needed by the job scripts, only for completed installations
# %q adds shell quoting so these export commands also work for folders with spaces
for tool in "${all_tools[@]}"; do
    if completed_tool "$tool" && [[ -n ${tool_exports[$tool]:-} ]]; then
        assignment=${tool_exports[$tool]}
        printf 'export %s=%q\n' "${assignment%%=*}" "$install_root/${assignment#*=}"
    fi
done
printf 'export BENCH_INSTALL_MANIFEST=%q\n' "$manifest"
if completed_tool cfdnalab; then printf 'export BENCH_CFDNA_BUILD=%q\n' "$install_root/cfdnalab-build.tsv"; fi
