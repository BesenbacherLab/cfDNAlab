#!/usr/bin/env bash
#SBATCH --job-name=cfdna-reference-gc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --no-requeue

# Calculate a GC reference once, with the midpoint benchmark's five exclusion BEDs
# Submit the separate preparation job using benchmarks/README.md
# BENCH_RESULTS names a new directory, independent of any existing GC package
# Input staging and copying the completed package are outside the measured step

set -euo pipefail
export LC_ALL=C
: "${SLURM_JOB_ID:?Submit reference preparation through sbatch}"
benchmarks="${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks"
source "$benchmarks/common.sh"
: "${BENCH_CFDNA:?Load the installed tool paths}"
: "${BENCH_REF_2BIT:?Set the alignment reference 2bit path}"
: "${BENCH_GRIFFIN_REF_DIR:?Load the shared exclusion BED directory}"
: "${BENCH_RESULTS:?Choose a new reference-results directory}"
: "${BENCH_CORES:?Set the core count for reference preparation}"
if [[ ! "$BENCH_CORES" =~ ^[1-9][0-9]*$ ]] || ((BENCH_CORES > ${SLURM_CPUS_PER_TASK:-0})); then
    echo 'BENCH_CORES must be positive and fit within SLURM_CPUS_PER_TASK.' >&2
    exit 1
fi
if [[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 ]]; then
    echo 'Reference preparation requires a single allocated node.' >&2
    exit 1
fi
if [[ -e "$BENCH_RESULTS" || -L "$BENCH_RESULTS" ]]; then
    echo "Reference results directory already exists: $BENCH_RESULTS. Choose a new directory." >&2
    exit 1
fi

# Resolve scratch on the compute node, using the same setting as benchmark jobs
LOCAL_SCRATCH_VAR=${LOCAL_SCRATCH_VAR-TMPDIR}
if [[ ! "$LOCAL_SCRATCH_VAR" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]]; then
    echo 'LOCAL_SCRATCH_VAR must name the variable supplying node-local scratch.' >&2
    exit 1
fi
scratch_root=${!LOCAL_SCRATCH_VAR:-}
if [[ -z "$scratch_root" || ! -d "$scratch_root" || ! -w "$scratch_root" ]]; then
    echo "The compute node must provide a writable scratch directory through $LOCAL_SCRATCH_VAR." >&2
    exit 1
fi
BENCH_NODE=$(hostname)
if [[ -n ${BENCH_EXPECTED_NODE:-} && "$BENCH_NODE" != "$BENCH_EXPECTED_NODE" ]]; then
    echo "Expected node $BENCH_EXPECTED_NODE, received $BENCH_NODE." >&2
    exit 1
fi
BENCH_WORK=$(mktemp -d "$scratch_root/reference-gc-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
export BENCH_WORK BENCH_NODE
export TMPDIR="$BENCH_WORK/tmp" XDG_CACHE_HOME="$BENCH_WORK/cache"
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/exclusions"

# One preparation run, without a warmup. Its records stay outside the feature results
BENCH_WARMUPS=0
export BENCH_PILOT=0 BENCH_RUN_TIMEOUT_SECONDS=0
timing_columns=(
    run_order repetition warmup method tool_version cores tool_threads
    node cpu_model cpu_ids physical_core_ids job_id step_id started_utc
    wall_seconds user_seconds system_seconds max_rss_kib tool_exit_code
    flush_seconds process_plus_flush_seconds output_directory_bytes pilot step_exit_code
)
# Create the chosen parent path, then claim a new results directory without reuse
mkdir -p "$(dirname "$BENCH_RESULTS")"
mkdir "$BENCH_RESULTS"
BENCH_RESULTS=$(cd "$BENCH_RESULTS" && pwd)
mkdir "$BENCH_RESULTS/runs"
(IFS=$'\t'; printf '%s\n' "${timing_columns[*]}") > "$BENCH_RESULTS/timings.tsv"
benchmark_commit=$(git -C "$SLURM_SUBMIT_DIR" rev-parse HEAD)
benchmark_changes=$(git -C "$SLURM_SUBMIT_DIR" status --porcelain)
{
    printf 'job_id\t%s\ncores\t%s\nnode\t%s\n' "$SLURM_JOB_ID" "$BENCH_CORES" "$BENCH_NODE"
    printf 'input_reference_2bit\t%s\n' "$BENCH_REF_2BIT"
    printf 'exclusion_reference_directory\t%s\n' "$BENCH_GRIFFIN_REF_DIR"
    printf 'local_scratch_variable\t%s\nlocal_scratch_directory\t%s\n' "$LOCAL_SCRATCH_VAR" "$scratch_root"
    printf 'output_prefix\tgriffin_blacklisted\n'
    printf 'benchmark_commit\t%s\n' "$benchmark_commit"
    if [[ -n "$benchmark_changes" ]]; then printf 'benchmark_dirty\ttrue\n'; else printf 'benchmark_dirty\tfalse\n'; fi
} > "$BENCH_RESULTS/run.tsv"
if ! scontrol show job "$SLURM_JOB_ID" > "$BENCH_RESULTS/slurm-job.txt" 2> "$BENCH_RESULTS/slurm-job-error.log"; then
    echo 'Slurm job details were unavailable. See slurm-job-error.log.' >&2
fi
if [[ -n ${BENCH_INSTALL_MANIFEST:-} ]]; then cp "$BENCH_INSTALL_MANIFEST" "$BENCH_RESULTS/installation.tsv"; fi
if [[ -n ${BENCH_CFDNA_BUILD:-} ]]; then cp "$BENCH_CFDNA_BUILD" "$BENCH_RESULTS/cfdnalab-build.tsv"; fi
cp "$benchmarks/measure_command.sh" "$BENCH_WORK/measure_command.sh"

log_progress 'Staging the executable, reference and exclusion BEDs on node-local scratch'
copy_executable "$BENCH_CFDNA" cfdna
cfdna_version=$("$BENCH_WORK/bin/cfdna" --version | tr '\t\n' '  ')
cp -L "$BENCH_REF_2BIT" "$BENCH_WORK/reference.2bit"
blacklist_arguments=()
for filename in encode_unified_GRCh38_exclusion_list.bed hg38_centromeres.bed hg38_gaps.bed hg38_fix_patches.bed hg38_alternative_haplotypes.bed; do
    cp -L "$BENCH_GRIFFIN_REF_DIR/$filename" "$BENCH_WORK/exclusions/"
    blacklist_arguments+=(--blacklist "$BENCH_WORK/exclusions/$filename")
done
sha256sum "$BENCH_WORK/reference.2bit" > "$BENCH_RESULTS/reference.sha256"
sha256sum "$BENCH_WORK/exclusions/"*.bed > "$BENCH_RESULTS/blacklists.sha256"
sync -f "$BENCH_WORK"

# Use the same measurement and CPU-placement checks as the sample benchmarks
output="$BENCH_WORK/output"
expected_output="$output/griffin_blacklisted.ref_gc_package.zarr"
stdout_file="$BENCH_WORK/tool-stdout.log"
empty_output_allowed=false
command=("$BENCH_WORK/bin/cfdna" ref-gc-bias
    --ref-2bit "$BENCH_WORK/reference.2bit"
    --output-dir "$output" --output-prefix griffin_blacklisted
    --n-threads "$BENCH_CORES" --temp-dir "$BENCH_WORK/tmp"
    "${blacklist_arguments[@]}")
row=([tool_version]="$cfdna_version" [tool_threads]="$BENCH_CORES")
run_order=0
run_in_slurm_step ref-gc-bias 1
save_run_record
if ((step_exit_code != 0)); then exit "$step_exit_code"; fi

# Publish only a successful package, after the measured Slurm step has finished
log_progress 'Saving the completed GC reference package'
cp -R "$expected_output" "$BENCH_RESULTS/"
sync -f "$BENCH_RESULTS"
package_path="$BENCH_RESULTS/griffin_blacklisted.ref_gc_package.zarr"
printf 'reference_gc_package\t%s\n' "$package_path" >> "$BENCH_RESULTS/run.tsv"
printf 'export BENCH_CFDNA_REF_GC=%q\n' "$package_path"
