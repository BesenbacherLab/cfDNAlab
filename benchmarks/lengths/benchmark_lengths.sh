#!/usr/bin/env bash
#SBATCH --job-name=cfdna-lengths
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --no-requeue

# This job benchmarks global fragment length histograms
#
# Submit from the repository root: bash benchmarks/submit_benchmark.sh lengths
# The launcher requests an exclusive node for full runs and allows shared pilots
# See README.md for installation, a concordance pilot and the full timing run
#
# The job:
#   1. Checks the settings
#   2. Copies the BAM, index, R helper and selected tools to the compute node (untimed)
#   3. Runs each method's warmups, then its measured repetitions. Every run is its own
#      Slurm step, bound to BENCH_CORES physical cores. The commands are in commands.sh
#   4. Makes Picard's filtered histogram for concordance (untimed)
#   5. Compares the saved histograms with compare_lengths.py (untimed)
# Results: timings.tsv (every run, including warmups), runs/<run_order>/ (each run's
# command and logs), profiles/<method>/ (first measured output) and concordance.tsv
# Host caches are not evicted

set -euo pipefail
export LC_ALL=C
# Slurm runs a copy of this script from elsewhere, so the repository is found
# through the submission directory. Submit from the repository root
benchmarks="${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks"
source "$benchmarks/common.sh"
source "$benchmarks/lengths/commands.sh"

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

BENCH_WARMUPS=${BENCH_WARMUPS:-1}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-3}
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-30}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-1000}
BENCH_GENOME=${BENCH_GENOME:-hg38}
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM:-}.bai}
check_job_settings lengths
if [[ "$BENCH_INDEX" != *.bai ]]; then
    echo "Supply a BAI index. FinaleToolkit requires one." >&2
    exit 1
fi

# Check the method names and the tool paths they need
needs_cfdna=false
needs_cfdnapro=false
needs_finale=false
needs_picard=false
for method in "${methods[@]}"; do
    case "$method" in
        cfdnalab|cfdnalab-autosomes|cfdnalab-proper) needs_cfdna=true ;;
        cfdnapro|cfdnapro-outer-span) needs_cfdnapro=true ;;
        finaletoolkit) needs_finale=true ;;
        picard) needs_picard=true ;;
        # samtools is always copied, because setup also uses it to read the BAM header
        samtools) ;;
        *) echo "Unknown length benchmark method: $method" >&2; exit 1 ;;
    esac
done
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the packed samtools environment for BAM header preparation}"
if [[ "$needs_cfdna" == true ]]; then : "${BENCH_CFDNA:?Set the cfDNAlab executable}"; fi
if [[ "$needs_cfdnapro" == true ]]; then : "${BENCH_CFDNAPRO_IMAGE:?Set the cfDNAPro image}"; fi
if [[ "$needs_finale" == true ]]; then : "${BENCH_FINALE_ARCHIVE:?Set the packed FinaleToolkit environment}"; fi
if [[ "$needs_picard" == true ]]; then : "${BENCH_PICARD_ARCHIVE:?Set the packed Picard environment}"; fi

# ---------------------------------------------------------------------------
# 2. Copy inputs and tools to the compute node (untimed)
# ---------------------------------------------------------------------------

BENCH_WORK=$(mktemp -d "$TMPDIR/lengths-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)
export BENCH_WORK BENCH_NODE

# Results are created before staging, so a setup failure leaves an inspectable record
# The columns of timings.tsv are listed in commands.sh
create_results_directory
log_progress "Initializing lengths benchmark on $BENCH_NODE"
record_job_setup lengths

mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib"
if [[ "$needs_cfdna" == true ]]; then copy_executable "$BENCH_CFDNA" cfdna; fi
if [[ "$needs_cfdnapro" == true ]]; then
    apptainer=$(command -v apptainer) || { echo "cfDNAPro requires Apptainer." >&2; exit 1; }
    cp -L "$BENCH_CFDNAPRO_IMAGE" "$BENCH_WORK/cfdnapro.sif"
    # The container can only see files under $BENCH_WORK, mounted as /work
    cp "$benchmarks/lengths/cfdnapro_lengths.R" "$BENCH_WORK/cfdnapro_lengths.R"
fi
if [[ "$needs_picard" == true ]]; then
    cp "$benchmarks/lengths/picard_filtered_histogram.sh" "$BENCH_WORK/picard_filtered_histogram.sh"
fi
export TMPDIR="$BENCH_WORK/tmp" APPTAINER_TMPDIR="$BENCH_WORK/tmp"
export XDG_CACHE_HOME="$BENCH_WORK/cache" MPLCONFIGDIR="$BENCH_WORK/matplotlib"
# Python user packages and outside import paths are disabled so imports use
# the copied environment. Outside Java options are cleared for the same reason
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
# JAVA_HOME selects the Java copied with Picard, rather than another Java on the node
export JAVA_HOME="$BENCH_WORK/picard-env"

unpack_environment "$BENCH_SAMTOOLS_ARCHIVE" "$BENCH_WORK/samtools-env"
if [[ "$needs_finale" == true ]]; then unpack_environment "$BENCH_FINALE_ARCHIVE" "$BENCH_WORK/finale-env"; fi
if [[ "$needs_picard" == true ]]; then unpack_environment "$BENCH_PICARD_ARCHIVE" "$BENCH_WORK/picard-env"; fi

# Tool paths, as used by commands.sh and experiment_commands.md
cfdna="$BENCH_WORK/bin/cfdna"
samtools="$BENCH_WORK/samtools-env/bin/samtools"
# A packed Python launcher can use env python3.11, which searches the node's PATH
# FinaleToolkit therefore runs with its own Python, here and in the timed command
finale_python="$BENCH_WORK/finale-env/bin/python"
finale="$BENCH_WORK/finale-env/bin/finaletoolkit"
picard="$BENCH_WORK/picard-env/bin/picard"
cfdnapro_image="$BENCH_WORK/cfdnapro.sif"

# Versions are checked before the large BAM is copied, so a broken tool fails early
# Tabs and newlines become spaces, so each version fits in a timings.tsv cell
log_progress "Checking tool versions"
cfdna_version=NA
finale_version=NA
picard_version=NA
cfdnapro_version=NA
if [[ "$needs_cfdna" == true ]]; then cfdna_version=$("$cfdna" --version | tr '\t\n' '  '); fi
if [[ "$needs_finale" == true ]]; then
    finale_version=$("$finale_python" "$finale" --version | tr '\t\n' '  ')
    finale_version+=" commit:$(cat "$BENCH_WORK/finale-env/share/finaletoolkit/revision.txt")"
fi
samtools_version=$("$samtools" --version-only | tr '\t\n' '  ')
if [[ "$needs_picard" == true ]]; then picard_version=$(benchmark_picard_version "$picard"); fi
if [[ "$needs_cfdnapro" == true ]]; then
    cfdnapro_version=$("$apptainer" exec --cleanenv --containall "$cfdnapro_image" Rscript --vanilla -e 'cat(as.character(packageVersion("cfDNAPro")))')
    cfdnapro_version+=" image-sha256:$(sha256sum "$cfdnapro_image" | cut -d ' ' -f 1)"
fi

log_progress "Copying BAM and index to node-local scratch"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.bai"
# The BAM's contig names and lengths come from its index. The unmapped '*' row is skipped
"$samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/chrom.sizes"
# Autosomes use the BAM's own names, chr1..chr22 or 1..22
autosomes=$(autosome_sizes "$BENCH_WORK/chrom.sizes" | cut -f 1 | paste -sd, -)
if [[ " ${methods[*]} " == *" cfdnalab-autosomes "* && -z "$autosomes" ]]; then
    echo "The BAM has no autosomes named chr1..chr22 or 1..22." >&2
    exit 1
fi

# The remaining names used by commands.sh
cores=$BENCH_CORES
bam="$BENCH_WORK/sample.bam"
output="$BENCH_WORK/output"
min_mapq=$BENCH_MIN_MAPQ
min_length=$BENCH_MIN_LENGTH
max_length=$BENCH_MAX_LENGTH
length_end=$((max_length + 1))
genome=$BENCH_GENOME

# Staged files finish writing before the first run. This does not evict the page cache
sync -f "$BENCH_WORK"
log_progress "Setup complete. Starting warmups and measurements"

# ---------------------------------------------------------------------------
# 3. Warmups and measured repetitions
# ---------------------------------------------------------------------------

# Each method runs as a block: its warmups first, then its measured repetitions
# A failed run is recorded, the method's remaining repetitions are skipped and the
# next method starts. The job then exits with status 1 at the end
run_order=0
job_exit_code=0
for method in "${methods[@]}"; do
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
            echo "Skipping the remaining repetitions of $method and continuing with the other methods." >&2
            job_exit_code=1
            break
        fi
    done
done

# ---------------------------------------------------------------------------
# 4. Filtered Picard histogram for concordance (untimed)
# ---------------------------------------------------------------------------

# Picard's timed run keeps all MAPQs. Concordance needs the shared MAPQ and flag
# filters, so a filtered run is made once here, in its own Slurm step with the same
# CPU binding. It adds job time but no timing row. A failure keeps the completed
# timings and the timed Picard histogram
if [[ -s "$BENCH_RESULTS/profiles/picard/metrics.txt" ]]; then
    log_progress "Generating the untimed filtered Picard concordance histogram"
    if ! srun --jobid="$SLURM_JOB_ID" --exact --nodes=1 --ntasks=1 \
        --cpus-per-task="$cores" --threads-per-core=1 --cpu-bind=threads \
        --export=ALL bash "$BENCH_WORK/picard_filtered_histogram.sh" "$cores" \
        "$BENCH_RESULTS/profiles/picard/concordance.metrics.txt"; then
        echo "Untimed Picard concordance failed. Completed timing results are retained." >&2
        rm -f "$BENCH_RESULTS/profiles/picard/concordance.metrics.txt"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# 5. Concordance (untimed)
# ---------------------------------------------------------------------------

compare_saved_profiles "$benchmarks/lengths/compare_lengths.py" \
    --min-length "$min_length" --max-length "$max_length"

# The job reports failure if any method failed, while keeping completed results
exit "$job_exit_code"
