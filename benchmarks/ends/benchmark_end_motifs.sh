#!/usr/bin/env bash
#SBATCH --job-name=cfdna-end-motifs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --no-requeue

# This job runs the tools' end-motif commands on the same node
#
# Submit from the repository root: bash benchmarks/submit_benchmark.sh ends
# The launcher requests an exclusive node for full runs and allows shared pilots
# See README.md for the settings, known differences and pilot invocation
#
# The job:
#   1. Checks the settings
#   2. Copies the BAM, index, reference, helpers and selected tools to the compute
#      node (untimed). Source files are never modified
#   3. Runs each method's warmups, then its measured repetitions. Every run is its own
#      Slurm step, bound to BENCH_CORES physical cores. The commands are in commands.sh
#   4. Compares the saved motif profiles with compare_end_motifs.py (untimed)
# Results: timings.tsv (every run, including warmups), runs/<run_order>/ (each run's
# command and logs), profiles/<method>/ (first measured output) and concordance.tsv
# No host caches are explicitly evicted. Only this job's scratch directory is deleted

set -euo pipefail
export LC_ALL=C
# Slurm runs a copy of this script from elsewhere, so the repository is found
# through the submission directory. Submit from the repository root
benchmarks="${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks"
source "$benchmarks/common.sh"
source "$benchmarks/ends/commands.sh"

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

BENCH_WARMUPS=${BENCH_WARMUPS:-1}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-3}
export BENCH_K=${BENCH_K:-4}
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-30}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-1000}
BENCH_GENOME=${BENCH_GENOME:-hg38}
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM:-}.bai}
check_job_settings ends
if [[ "$BENCH_INDEX" != *.bai ]]; then
    echo "FinaleToolkit requires a BAI index." >&2
    exit 1
fi
if [[ "$BENCH_K" != 3 && "$BENCH_K" != 4 ]]; then
    echo "Use a motif length of 3 or 4 (BENCH_K)." >&2
    exit 1
fi

# Check the method names and the tool paths they need
needs_cfdna=false
needs_finale=false
needs_cfdnapro=false
needs_reference=false
for method in "${methods[@]}"; do
    case "$method" in
        cfdnalab-reference|cfdnalab-reference-proper|cfdnalab-reference-autosomes)
            needs_cfdna=true
            needs_reference=true
            ;;
        cfdnalab-read|cfdnalab-read-aligned|cfdnalab-read-include) needs_cfdna=true ;;
        finaletoolkit)
            needs_finale=true
            needs_reference=true
            ;;
        cfdnapro) needs_cfdnapro=true ;;
        *) echo "Unknown end-motif method: $method" >&2; exit 1 ;;
    esac
done
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the packed samtools environment for BAM header preparation}"
if [[ "$needs_cfdna" == true ]]; then : "${BENCH_CFDNA:?Set the cfDNAlab executable}"; fi
if [[ "$needs_finale" == true ]]; then : "${BENCH_FINALE_ARCHIVE:?Set the packed FinaleToolkit environment}"; fi
if [[ "$needs_cfdnapro" == true ]]; then : "${BENCH_CFDNAPRO_IMAGE:?Set the cfDNAPro image}"; fi
if [[ "$needs_reference" == true ]]; then : "${BENCH_REF_2BIT:?Set the matching 2bit reference}"; fi

# ---------------------------------------------------------------------------
# 2. Copy inputs and tools to the compute node (untimed)
# ---------------------------------------------------------------------------

BENCH_WORK=$(mktemp -d "$TMPDIR/end-motifs-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)
export BENCH_WORK BENCH_NODE

# Results are created before staging, so a setup failure leaves an inspectable record
# The columns of timings.tsv are listed in commands.sh
create_results_directory
log_progress "Initializing end-motif benchmark on $BENCH_NODE"
record_job_setup ends

mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib"
log_progress "Copying BAM and index to node-local scratch"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.bai"
if [[ "$needs_reference" == true ]]; then cp -L "$BENCH_REF_2BIT" "$BENCH_WORK/reference.2bit"; fi
if [[ "$needs_cfdna" == true ]]; then copy_executable "$BENCH_CFDNA" cfdna; fi
if [[ "$needs_cfdnapro" == true ]]; then
    apptainer=$(command -v apptainer) || { echo "cfDNAPro requires Apptainer." >&2; exit 1; }
    cp -L "$BENCH_CFDNAPRO_IMAGE" "$BENCH_WORK/cfdnapro.sif"
    # The container can only see files under $BENCH_WORK, mounted as /work
    cp "$benchmarks/ends/cfdnapro_end_motifs.R" "$BENCH_WORK/cfdnapro_end_motifs.R"
fi
export TMPDIR="$BENCH_WORK/tmp" APPTAINER_TMPDIR="$BENCH_WORK/tmp"
# Python user packages and outside import paths are disabled so imports use
# the copied environments. The installation notes in ../README.md explain these settings
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
export XDG_CACHE_HOME="$BENCH_WORK/cache"
export MPLCONFIGDIR="$BENCH_WORK/matplotlib"

# samtools reads the BAM's contig names and lengths from its index. The unmapped
# '*' row is skipped. Autosomes use the BAM's own names, chr1..chr22 or 1..22
unpack_environment "$BENCH_SAMTOOLS_ARCHIVE" "$BENCH_WORK/samtools-env"
"$BENCH_WORK/samtools-env/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/chrom.sizes"
autosomes=$(autosome_sizes "$BENCH_WORK/chrom.sizes" | cut -f 1 | paste -sd, -)
if [[ " ${methods[*]} " == *" cfdnalab-reference-autosomes "* && -z "$autosomes" ]]; then
    echo "The BAM has no autosomes named chr1..chr22 or 1..22." >&2
    exit 1
fi
if [[ "$needs_finale" == true ]]; then unpack_environment "$BENCH_FINALE_ARCHIVE" "$BENCH_WORK/finale-env"; fi

# FinaleToolkit visits the reference's contigs, while cfDNAlab --chromosomes all uses
# the BAM's contigs. This check requires matching names and lengths. It reads
# headers only, with no filtering or changes to the BAM
if [[ "$needs_finale" == true ]]; then
    "$BENCH_WORK/finale-env/bin/python" - <<'PYHEADER'
import os
import py2bit
import pysam
work = os.environ['BENCH_WORK']
with pysam.AlignmentFile(work + '/sample.bam', 'rb') as bam:
    bam_contigs = dict(zip(bam.references, bam.lengths))
with py2bit.open(work + '/reference.2bit') as reference:
    if bam_contigs != reference.chroms():
        raise SystemExit('BAM and 2bit contig names/lengths differ.')
for chromosome, length in bam_contigs.items():
    if length % 1000000 == 0:
        raise SystemExit('The pinned FinaleToolkit window builder omits its last window for this contig length: ' + chromosome)
PYHEADER
fi

# Tool paths, as used by commands.sh and experiment_commands.md
cfdna="$BENCH_WORK/bin/cfdna"
# A packed Python launcher can use env python3.11, which searches the node's PATH
# FinaleToolkit therefore runs with its own Python, here and in the timed command
finale_python="$BENCH_WORK/finale-env/bin/python"
finale="$BENCH_WORK/finale-env/bin/finaletoolkit"
cfdnapro_image="$BENCH_WORK/cfdnapro.sif"

# Tabs and newlines become spaces, so each version fits in a timings.tsv cell
log_progress "Checking tool versions"
cfdna_version=NA
finale_version=NA
cfdnapro_version=NA
if [[ "$needs_cfdna" == true ]]; then cfdna_version=$("$cfdna" --version | tr '\t\n' '  '); fi
if [[ "$needs_finale" == true ]]; then
    finale_version=$("$finale_python" "$finale" --version | tr '\t\n' '  ')
    finale_version+=" commit:$(cat "$BENCH_WORK/finale-env/share/finaletoolkit/revision.txt")"
fi
if [[ "$needs_cfdnapro" == true ]]; then
    cfdnapro_version=$("$apptainer" exec --cleanenv --containall "$cfdnapro_image" Rscript --vanilla -e 'cat(as.character(packageVersion("cfDNAPro")))')
    cfdnapro_version+=" image-sha256:$(sha256sum "$cfdnapro_image" | cut -d ' ' -f 1) genome:$BENCH_GENOME"
fi

# The remaining names used by commands.sh
cores=$BENCH_CORES
bam="$BENCH_WORK/sample.bam"
output="$BENCH_WORK/output"
k=$BENCH_K
min_mapq=$BENCH_MIN_MAPQ
min_length=$BENCH_MIN_LENGTH
max_length=$BENCH_MAX_LENGTH
genome=$BENCH_GENOME
reference_2bit="$BENCH_WORK/reference.2bit"

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
# 4. Concordance (untimed)
# ---------------------------------------------------------------------------

compare_saved_profiles "$benchmarks/ends/compare_end_motifs.py" --k "$k"

# The job reports failure if any method failed, while keeping completed results
exit "$job_exit_code"
