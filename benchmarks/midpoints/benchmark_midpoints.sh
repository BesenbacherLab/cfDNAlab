#!/usr/bin/env bash
#SBATCH --job-name=cfdna-midpoints
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --no-requeue

# This job produces native midpoint profiles with 165 bp smoothing at 15 bp resolution
#
# Submit from the repository root: bash benchmarks/submit_benchmark.sh midpoints
# The launcher requests an exclusive node for full runs and allows shared pilots
# See README.md for required paths, concordance and expected differences
#
# Every midpoint job runs a single method (cfdnalab, cfdnalab-gc or griffin-gc). It:
#   1. Checks the settings
#   2. Copies the BAM, references, sites and selected tools to the compute node, and
#      converts the sites to each tool's format with prepare_sites.py (untimed)
#   3. For Griffin, runs the release's short demo as an untimed warmup in its own step
#   4. Runs the warmups and measured repetitions. Every run is its own Slurm step,
#      bound to BENCH_CORES physical cores. The commands are in commands.sh
# Full-input runs get timing rows. Griffin's short demo warmup is untimed
# Results: timings.tsv, runs/<run_order>/ (each run's command and logs) and
# profiles/<method>/ (first measured output). Concordance needs profiles from
# separate jobs, so it is run afterwards (see README.md)
# All BAM/reference copies, environments and temporary files use local scratch

set -euo pipefail
export LC_ALL=C
# Slurm runs a copy of this script from elsewhere, so the repository is found
# through the submission directory. Submit from the repository root
benchmarks="${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks"
source "$benchmarks/common.sh"
source "$benchmarks/midpoints/commands.sh"

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

: "${BENCH_INTERVALS:?Set the grouped fixed-size BED path}"
BENCH_WARMUPS=${BENCH_WARMUPS:-1}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-3}
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-20}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-100}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-200}
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM:-}.bai}
BENCH_FASTA_INDEX=${BENCH_FASTA_INDEX:-${BENCH_REF_FASTA:-}.fai}
check_job_settings midpoints

if ((${#methods[@]} != 1)); then
    if [[ " ${methods[*]} " == *" griffin-gc "* ]]; then
        echo "Run Griffin in its own job: set BENCH_METHODS=griffin-gc." >&2
    else
        echo "Submit each midpoint workload in a separate job and results directory." >&2
    fi
    exit 1
fi
method=${methods[0]}
case "$method" in
    cfdnalab|cfdnalab-gc|griffin-gc) ;;
    *) echo "Unknown midpoint tool: $method" >&2; exit 1 ;;
esac
if [[ "$method" == griffin-gc ]] && ((BENCH_MIN_LENGTH < 15 || BENCH_MAX_LENGTH > 500)); then
    echo "Griffin GC estimation uses 15-500 bp. Keep the profiling range within those bounds." >&2
    exit 1
fi

# Check the paths the method needs
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the packed samtools environment for BAM header preparation}"
: "${BENCH_GRIFFIN_REF_DIR:?Set the Griffin Ref directory containing the shared exclusion BEDs}"
if [[ "$method" == cfdnalab* ]]; then : "${BENCH_CFDNA:?Set the cfDNAlab executable}"; fi
if [[ "$method" == cfdnalab-gc ]]; then
    : "${BENCH_REF_2BIT:?Set the matching reference 2bit file}"
    : "${BENCH_CFDNA_REF_GC:?Set the precomputed cfDNAlab reference GC Zarr package}"
fi
if [[ "$method" == griffin-gc ]]; then
    : "${BENCH_GRIFFIN_ARCHIVE:?Set the packed Griffin environment}"
    : "${BENCH_REF_FASTA:?Set the matching reference FASTA}"
    if [[ "$BENCH_REF_FASTA" == *.gz ]]; then
        echo "Use an uncompressed FASTA and its existing .fai index." >&2
        exit 1
    fi
    : "${BENCH_GRIFFIN_GC_REGIONS:?Set the autosomal BED matching the Griffin reference GC tables}"
    : "${BENCH_GRIFFIN_REF_GC:?Set the precomputed Griffin reference GC-frequency directory}"
    BENCH_GRIFFIN_MAPPABILITY_BW=${BENCH_GRIFFIN_MAPPABILITY_BW:-$BENCH_GRIFFIN_REF_DIR/k100.Umap.MultiTrackMappability.bw}
    # The release demo replaces full-input warmups, including in pilot jobs
    BENCH_WARMUPS=0
fi

# ---------------------------------------------------------------------------
# 2. Copy inputs and tools to the compute node, and prepare the sites (untimed)
# ---------------------------------------------------------------------------

BENCH_WORK=$(mktemp -d "$TMPDIR/midpoints-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)
export BENCH_WORK BENCH_NODE

# Results are created before staging, so a setup failure leaves an inspectable record
# The columns of timings.tsv are listed in commands.sh
create_results_directory
log_progress "Initializing midpoint benchmark on $BENCH_NODE"
record_job_setup midpoints
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib"

# The same public exclusion BEDs are passed to each tool's native masking code
# cfDNAlab can drop intervals near them. Griffin masks individual positions
# Those native differences are kept rather than forcing identical site counts
mkdir -p "$BENCH_WORK/exclusions"
exclusions=()
for exclusion_file in encode_unified_GRCh38_exclusion_list.bed hg38_centromeres.bed hg38_gaps.bed hg38_fix_patches.bed hg38_alternative_haplotypes.bed; do
    cp -L "$BENCH_GRIFFIN_REF_DIR/$exclusion_file" "$BENCH_WORK/exclusions/"
    exclusions+=("$BENCH_WORK/exclusions/$exclusion_file")
done
printf 'exclusion_reference_directory\t%s\n' "$BENCH_GRIFFIN_REF_DIR" >> "$BENCH_RESULTS/run.tsv"
if [[ "$method" == griffin-gc ]]; then
    log_progress "Copying Griffin mappability mask to node-local scratch"
    cp -L "$BENCH_GRIFFIN_MAPPABILITY_BW" "$BENCH_WORK/mappability.bw"
    printf 'griffin_mappability_bw\t%s\n' "$BENCH_GRIFFIN_MAPPABILITY_BW" >> "$BENCH_RESULTS/run.tsv"
fi

# Python packages, font caches and temporary files stay with this job's inputs
# Removing inherited Python paths keeps imports inside the selected environments
export TMPDIR="$BENCH_WORK/tmp" XDG_CACHE_HOME="$BENCH_WORK/cache"
export MPLCONFIGDIR="$BENCH_WORK/matplotlib" MPLBACKEND=Agg PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
if [[ "$method" == cfdnalab* ]]; then copy_executable "$BENCH_CFDNA" cfdna; fi
unpack_environment "$BENCH_SAMTOOLS_ARCHIVE" "$BENCH_WORK/samtools-env"
if [[ "$method" == griffin-gc ]]; then
    unpack_environment "$BENCH_GRIFFIN_ARCHIVE" "$BENCH_WORK/griffin-env"
    for script in griffin_GC_counts.py griffin_GC_bias.py; do
        if [[ ! -f "$BENCH_WORK/griffin-env/share/griffin/$script" ]]; then
            echo "The Griffin archive lacks $script. Rerun install_tools.sh ROOT --only griffin and use the new archive." >&2
            exit 1
        fi
    done
    # Griffin's GC bias step draws plots, so check that matplotlib loads
    "$BENCH_WORK/griffin-env/bin/python" -c 'import matplotlib.pyplot'
    # pybedtools must find the packed bedtools, rather than a shared installation
    export PATH="$BENCH_WORK/griffin-env/bin:$PATH"
fi

log_progress "Checking tool versions"
cfdna="$BENCH_WORK/bin/cfdna"
cfdna_version=NA
griffin_version=NA
if [[ "$method" == cfdnalab* ]]; then cfdna_version=$("$cfdna" --version 2>&1 | tr '\t\n' '  '); fi
if [[ "$method" == griffin-gc ]]; then griffin_version=$(cat "$BENCH_WORK/griffin-env/share/griffin/version.txt"); fi

# Reference GC inputs are prepared once, independently of these per-sample runs
# Their paths are recorded so corrected runs cannot be mistaken for the baseline
if [[ "$method" == cfdnalab ]]; then
    printf 'gc_workload\tuncorrected-profiles\n' >> "$BENCH_RESULTS/run.tsv"
else
    printf 'gc_workload\tsample-bias-and-profiles\n' >> "$BENCH_RESULTS/run.tsv"
    cp "$benchmarks/midpoints/gc_corrected_profile.sh" "$BENCH_WORK/gc_corrected_profile.sh"
fi
if [[ "$method" == cfdnalab-gc ]]; then
    # The default reference package is reused. Native GC loading checks reference
    # compatibility, and midpoint correction checks that its range is covered
    log_progress "Copying cfDNAlab reference GC inputs to node-local scratch"
    cp -L "$BENCH_REF_2BIT" "$BENCH_WORK/reference.2bit"
    cp -RL "$BENCH_CFDNA_REF_GC" "$BENCH_WORK/reference-gc.zarr"
    printf 'reference_2bit\t%s\ncfdnalab_reference_gc\t%s\n' "$BENCH_REF_2BIT" "$BENCH_CFDNA_REF_GC" >> "$BENCH_RESULTS/run.tsv"
fi
if [[ "$method" == griffin-gc ]]; then
    cp -L "$BENCH_GRIFFIN_GC_REGIONS" "$BENCH_WORK/gc-regions.bed"
    sha256sum "$BENCH_GRIFFIN_GC_REGIONS" > "$BENCH_RESULTS/gc-regions.sha256"
    cp "$BENCH_WORK/gc-regions.bed" "$BENCH_RESULTS/gc-regions.bed"
    printf 'griffin_gc_regions\t%s\n' "$BENCH_GRIFFIN_GC_REGIONS" >> "$BENCH_RESULTS/run.tsv"
    # Griffin names its reference GC tables after the estimation BED, without its extension
    export BENCH_GC_REFERENCE_NAME
    BENCH_GC_REFERENCE_NAME=${BENCH_GRIFFIN_GC_REGIONS##*/}
    BENCH_GC_REFERENCE_NAME=${BENCH_GC_REFERENCE_NAME%.*}
    # Griffin estimates bias over its release config's 15-500 bp range
    # These tables include the neighbours used to smooth the 100-200 bp profiles
    printf 'griffin_gc_estimation_min_length\t15\ngriffin_gc_estimation_max_length\t500\n' >> "$BENCH_RESULTS/run.tsv"
    log_progress "Copying Griffin reference GC frequencies to node-local scratch"
    mkdir -p "$BENCH_WORK/reference-gc"
    for ((fragment_length = 15; fragment_length <= 500; fragment_length++)); do
        reference_file="$BENCH_GRIFFIN_REF_GC/$BENCH_GC_REFERENCE_NAME.${fragment_length}bp.GC_frequency.txt"
        if [[ ! -s "$reference_file" ]]; then echo "Missing Griffin reference GC frequency: $reference_file" >&2; exit 1; fi
        cp -L "$reference_file" "$BENCH_WORK/reference-gc/"
    done
    printf 'griffin_reference_gc\t%s\n' "$BENCH_GRIFFIN_REF_GC" >> "$BENCH_RESULTS/run.tsv"
fi

log_progress "Copying BAM and index to node-local scratch"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.${BENCH_INDEX##*.}"
cp -L "$BENCH_INTERVALS" "$BENCH_WORK/input.bed"
if [[ "$method" == griffin-gc ]]; then
    cp -L "$BENCH_REF_FASTA" "$BENCH_WORK/reference.fa"
    cp -L "$BENCH_FASTA_INDEX" "$BENCH_WORK/reference.fa.fai"
    cp "$benchmarks/midpoints/griffin_profile.sh" "$BENCH_WORK/griffin_profile.sh"
fi
# The BAM's contig names and lengths come from its index. The unmapped '*' row is skipped
"$BENCH_WORK/samtools-env/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/bam.chrom.sizes"

# Convert the same windows to cfDNAlab's BED and Griffin's site tables. See
# prepare_sites.py for the files it writes. Its summary line gives the window
# width, half width, number of sites and groups, and the chromosomes with sites
site_arguments=("$BENCH_WORK" --min-length "$BENCH_MIN_LENGTH" --max-length "$BENCH_MAX_LENGTH")
if [[ "$method" == griffin-gc ]]; then site_arguments+=(--griffin); fi
site_summary=$("$BENCH_SETUP_PYTHON" "$benchmarks/midpoints/prepare_sites.py" "${site_arguments[@]}")
export BENCH_WINDOW_WIDTH BENCH_HALF_WIDTH BENCH_SITES BENCH_GROUPS BENCH_CHROMOSOMES
IFS=$'\t' read -r BENCH_WINDOW_WIDTH BENCH_HALF_WIDTH BENCH_SITES BENCH_GROUPS BENCH_CHROMOSOMES <<< "$site_summary"
cp "$BENCH_WORK/sites.bed" "$BENCH_RESULTS/sites.bed"
sha256sum "$BENCH_INTERVALS" > "$BENCH_RESULTS/sites.sha256"

# The remaining names used by commands.sh
cores=$BENCH_CORES
bam="$BENCH_WORK/sample.bam"
output="$BENCH_WORK/output"
min_mapq=$BENCH_MIN_MAPQ
min_length=$BENCH_MIN_LENGTH
max_length=$BENCH_MAX_LENGTH
length_end=$((max_length + 1))

# Staged files finish writing before the first run. This does not evict the page cache
sync -f "$BENCH_WORK"
log_progress "Setup complete. Starting warmups and measurements"

# ---------------------------------------------------------------------------
# 3. Griffin's untimed demo warmup
# ---------------------------------------------------------------------------

# A full Griffin run can take more than a day, so Griffin gets no full-input warmup
# Instead it profiles the release's 1000-site CTCF demo on the same BAM, reference
# and node, with unit GC weights. This runs in its own Slurm step, so its memory is
# not counted in a measured step, and it has no timing row
if [[ "$method" == griffin-gc ]]; then
    printf 'griffin_warmup\trelease-demo-1000-sites\n' >> "$BENCH_RESULTS/run.tsv"
    sha256sum "$BENCH_WORK/griffin-env/share/griffin/CTCF.hg38.1000.txt" > "$BENCH_RESULTS/griffin-warmup.sha256"
    warmup_chromosomes=$(cut -f 1 "$BENCH_WORK/warmup.chrom.sizes" | paste -sd, -)
    if ! srun --jobid="$SLURM_JOB_ID" --exact --nodes=1 --ntasks=1 \
        --cpus-per-task="$cores" --threads-per-core=1 --cpu-bind=threads \
        --export=ALL bash "$BENCH_WORK/griffin_profile.sh" "$cores" \
        "$BENCH_WORK/warmup-sites.yaml" "$BENCH_WORK/warmup.chrom.sizes" \
        "$warmup_chromosomes" "$BENCH_WORK/warmup-output" \
        > "$BENCH_RESULTS/griffin-warmup.log" 2> "$BENCH_RESULTS/griffin-warmup-error.log"; then
        cat "$BENCH_RESULTS/griffin-warmup-error.log" "$BENCH_RESULTS/griffin-warmup.log" >&2
        echo "Griffin's demo warmup failed. Measurements have not started." >&2
        exit 1
    fi
    rm -rf "$BENCH_WORK/warmup-output"
    sync -f "$BENCH_WORK"
fi

# ---------------------------------------------------------------------------
# 4. Warmups and measured repetitions
# ---------------------------------------------------------------------------

# A failed run is recorded and the remaining repetitions are skipped. The job then
# exits with status 1
run_order=0
job_exit_code=0
for ((repetition = 1; repetition <= BENCH_WARMUPS + BENCH_REPETITIONS; repetition++)); do
    prepare_command "$method"
    run_in_slurm_step "$method" "$repetition"
    save_run_record
    # The first measured repetition's output is kept for concordance
    if ((step_exit_code == 0 && repetition == BENCH_WARMUPS + 1)); then
        save_profile "$method"
    fi
    clear_run_output
    if ((step_exit_code != 0)); then
        echo "Skipping the remaining repetitions of $method." >&2
        job_exit_code=1
        break
    fi
done

# Concordance compares profiles from separate jobs, so it is not run here
# See README.md for combining the saved profiles and running compare_midpoints.py

# The job reports failure if any run failed, while keeping completed results
exit "$job_exit_code"
