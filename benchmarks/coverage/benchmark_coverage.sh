#!/usr/bin/env bash
#SBATCH --job-name=cfdna-coverage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --no-requeue

# This job benchmarks per-base coverage and 500 kb regional coverage
#
# Submit from the repository root: bash benchmarks/submit_benchmark.sh coverage
# The launcher requests BENCH_CORES CPUs for this job's single core setting
# Full runs also request --exclusive and --mem=0 to reserve the node and its memory
# Pilots allow shared nodes and request 50 GB unless --mem is supplied
# User setup: the chosen partition must honor exclusive allocations for full runs
# Its policy takes precedence over the job's sharing request
# https://slurm.schedmd.com/sbatch.html#OPT_exclusive
#
# The job:
#   1. Checks the settings
#   2. Copies the BAM, index and selected tools to the compute node, and writes the
#      shared 500 kb windows (untimed)
#   3. Runs each method's warmups, then its measured repetitions. Every run is its own
#      Slurm step, bound to BENCH_CORES physical cores. The commands are in commands.sh
#   4. Compares saved profiles that measure the same quantity with
#      compare_coverage.py (untimed)
# Methods are grouped by the quantity they measure (see coverage_quantity in
# commands.sh): per-base read depth, per-base fragment coverage, regional read depth
# and regional fragment counts. Each quantity is compared separately
# Results: timings.tsv (every run, including warmups), runs/<run_order>/ (each run's
# command and logs), profiles/<method>/ (first measured output), concordance.tsv,
# residuals.tsv and coverage_jaccard.tsv
#
# Example submission from the repository root:
#   bash benchmarks/install_tools.sh /path/to/benchmark-tools
#   export BENCH_BAM=/faststorage/project/sample.bam
#   export BENCH_CFDNA=/path/to/cfdna
#   export BENCH_MOSDEPTH=/path/to/mosdepth
#   export BENCH_SAMTOOLS_ARCHIVE=/path/to/benchmark-tools/samtools.tar.gz
#   export BENCH_PANDEPTH=/path/to/pandepth
#   export BENCH_FINALE_ARCHIVE=/path/to/benchmark-tools/finaletoolkit.tar.gz
#   The batch job uses TMPDIR supplied on the compute node
#   export BENCH_RESULTS=/faststorage/project/coverage-results
#   export BENCH_CORES=4
#   export BENCH_METHODS="cfdnalab mosdepth"
#   bash benchmarks/submit_benchmark.sh coverage --constraint=gen6 --time=YOUR_PILOT_BASED_BUDGET
# User setup: adjust the generation selector, partition, time and paths for the cluster

# Shell errors stop the job, and LC_ALL=C gives decimal numbers a dot in the TSV
set -euo pipefail
export LC_ALL=C
# Slurm runs a copy of this script from elsewhere, so the repository is found
# through the submission directory. Submit from the repository root
benchmarks="${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks"
source "$benchmarks/common.sh"
source "$benchmarks/coverage/commands.sh"

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

# Input and tool paths are supplied in the environment when submitting the job
# BENCH_INDEX can override the usual sample.bam.bai path. FinaleToolkit requires
# a BAI index, so CSI indices are not supported by this combined benchmark
# Copied programs must have the libraries they need available on the compute node
# BENCH_FINALE_ARCHIVE contains a Linux Conda environment with FinaleToolkit
# installed, packed with conda-pack so Python and its packages can run locally
# BENCH_SAMTOOLS_ARCHIVE similarly includes samtools and its shared libraries
BENCH_WARMUPS=${BENCH_WARMUPS:-1}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-3}
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-30}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-1000}
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM:-}.bai}
check_job_settings coverage
if [[ "$BENCH_INDEX" != *.bai ]]; then
    echo "Supply a BAI index. FinaleToolkit requires one." >&2
    exit 1
fi

# Check the method names and the tool paths they need
needs_cfdna=false
needs_mosdepth=false
needs_pandepth=false
needs_finale=false
for method in "${methods[@]}"; do
    case "$method" in
        cfdnalab|cfdnalab-autosomes|cfdnalab-fragment|cfdnalab-fragment-proper|cfdnalab-regions|cfdnalab-region-counts|cfdnalab-region-counts-proper)
            needs_cfdna=true ;;
        mosdepth|mosdepth-fragment|mosdepth-regions) needs_mosdepth=true ;;
        pandepth|pandepth-regions) needs_pandepth=true ;;
        finaletoolkit-coverage) needs_finale=true ;;
        # samtools is always copied, because setup also uses it to read the BAM header
        samtools-depth|samtools-regions) ;;
        *) echo "Unknown coverage method: $method" >&2; exit 1 ;;
    esac
done
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the path of the packed samtools environment}"
if [[ "$needs_cfdna" == true ]]; then : "${BENCH_CFDNA:?Set the cfDNAlab executable}"; fi
if [[ "$needs_mosdepth" == true ]]; then : "${BENCH_MOSDEPTH:?Set the mosdepth executable}"; fi
if [[ "$needs_pandepth" == true ]]; then : "${BENCH_PANDEPTH:?Set the PanDepth executable}"; fi
if [[ "$needs_finale" == true ]]; then : "${BENCH_FINALE_ARCHIVE:?Set the packed FinaleToolkit environment}"; fi

# ---------------------------------------------------------------------------
# 2. Copy inputs and tools to the compute node (untimed)
# ---------------------------------------------------------------------------

BENCH_WORK=$(mktemp -d "$TMPDIR/coverage-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)
export BENCH_WORK BENCH_NODE

# Results are created before staging, so a setup failure leaves an inspectable record
# The columns of timings.tsv are listed in commands.sh
create_results_directory
log_progress "Initializing coverage benchmark on $BENCH_NODE"
record_job_setup coverage
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib"
export TMPDIR="$BENCH_WORK/tmp"
# Python user packages and outside import paths are disabled so imports use
# the copied environment
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
export XDG_CACHE_HOME="$BENCH_WORK/cache" MPLCONFIGDIR="$BENCH_WORK/matplotlib"

# Only selected tools are copied
unpack_environment "$BENCH_SAMTOOLS_ARCHIVE" "$BENCH_WORK/samtools-env"
if [[ "$needs_cfdna" == true ]]; then copy_executable "$BENCH_CFDNA" cfdna; fi
if [[ "$needs_mosdepth" == true ]]; then copy_executable "$BENCH_MOSDEPTH" mosdepth; fi
if [[ "$needs_pandepth" == true ]]; then copy_executable "$BENCH_PANDEPTH" pandepth; fi
if [[ "$needs_finale" == true ]]; then unpack_environment "$BENCH_FINALE_ARCHIVE" "$BENCH_WORK/finale-env"; fi

# Tool paths, as used by commands.sh and experiment_commands.md
cfdna="$BENCH_WORK/bin/cfdna"
mosdepth="$BENCH_WORK/bin/mosdepth"
pandepth="$BENCH_WORK/bin/pandepth"
samtools="$BENCH_WORK/samtools-env/bin/samtools"
# A packed Python launcher can use env python3.11, which searches the node's PATH
# FinaleToolkit therefore runs with its own Python, here and in the timed command
finale_python="$BENCH_WORK/finale-env/bin/python"
finale="$BENCH_WORK/finale-env/bin/finaletoolkit"

# Versions are checked before the large BAM is copied, so a broken tool fails early
# Tabs and newlines become spaces, so each version fits in a timings.tsv cell
log_progress "Checking tool versions"
cfdna_version=NA
mosdepth_version=NA
pandepth_version=NA
finale_version=NA
if [[ "$needs_cfdna" == true ]]; then cfdna_version=$("$cfdna" --version | tr '\t\n' '  '); fi
if [[ "$needs_mosdepth" == true ]]; then mosdepth_version=$("$mosdepth" --version | tr '\t\n' '  '); fi
samtools_version=$("$samtools" --version-only | tr '\t\n' '  ')
if [[ "$needs_finale" == true ]]; then
    finale_version=$("$finale_python" "$finale" --version | tr '\t\n' '  ')
    finale_version+=" commit:$(cat "$BENCH_WORK/finale-env/share/finaletoolkit/revision.txt")"
fi
if [[ "$needs_pandepth" == true ]]; then
    # PanDepth has no version option. Its help text contains "show this help [vX.Y]"
    pandepth_version=$("$pandepth" -h 2>&1 | sed -n 's/.*show this help \[\(v[^]]*\)\].*/\1/p')
    if [[ -z "$pandepth_version" ]]; then echo "Could not read PanDepth version." >&2; exit 1; fi
fi

log_progress "Copying BAM and index to node-local scratch"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.bai"
# The BAM's contig names and lengths come from its index. The unmapped '*' row is skipped
"$samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/chrom.sizes"

# The shared BED4 of disjoint 500 kb windows over every BAM contig, with a shorter
# final window per contig. FinaleToolkit and samtools bedcov read this BED
# cfDNAlab, mosdepth and PanDepth make the same windows with their own options
awk 'BEGIN {OFS="\t"} {
    for (window_start = 0; window_start < $2; window_start += 500000) {
        window_end = window_start + 500000
        if (window_end > $2) window_end = $2
        print $1, window_start, window_end, "."
    }
}' "$BENCH_WORK/chrom.sizes" > "$BENCH_WORK/regions.bed"
if [[ ! -s "$BENCH_WORK/regions.bed" ]]; then
    echo "No coverage windows could be generated from the BAM index." >&2
    exit 1
fi
# The comparator reads the windows from the results directory
cp "$BENCH_WORK/regions.bed" "$BENCH_RESULTS/regions.bed"

# Autosomes use the BAM's own names. Snyder's reference can use 1..22, while
# cfDNAlab defaults to chr1..chr22. The file's chromosome names remain unchanged
autosomes=$(awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ {print $1}' "$BENCH_WORK/chrom.sizes" | paste -sd, -)
if [[ " ${methods[*]} " == *" cfdnalab-autosomes "* && -z "$autosomes" ]]; then
    echo "The BAM has no autosomes named chr1..chr22 or 1..22." >&2
    exit 1
fi

# The remaining names used by commands.sh
cores=$BENCH_CORES
bam="$BENCH_WORK/sample.bam"
output="$BENCH_WORK/output"
regions="$BENCH_WORK/regions.bed"
min_mapq=$BENCH_MIN_MAPQ
min_length=$BENCH_MIN_LENGTH
max_length=$BENCH_MAX_LENGTH

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
        # The first measured repetition's output is kept for concordance. This
        # includes every processed chromosome, also for samtools and PanDepth
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
# 4. Concordance (untimed)
# ---------------------------------------------------------------------------

# Only profiles of the same quantity are compared, so concordance needs two saved
# profiles with the same quantity
if saved_profiles_share_a_quantity; then
    # All-contig profiles keep the full scope. The optional autosome configuration
    # is compared on the chromosomes it shares with the all-contig profiles
    compare_saved_profiles "$benchmarks/coverage/compare_coverage.py" --common-scope
else
    echo 'Concordance skipped: no pair of saved profiles measures the same quantity.'
fi

# The job reports failure if any method failed, while keeping completed results
exit "$job_exit_code"
