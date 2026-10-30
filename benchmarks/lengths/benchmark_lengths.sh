#!/usr/bin/env bash
#SBATCH --job-name=cfdna-lengths
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --no-requeue

# Benchmark global fragment length histograms on an exclusively allocated node
# See README.md for installation, a concordance pilot and the full timing run
# Save every run in timings.tsv, including warmups, plus a native profile per method
# Stage inputs and dependencies locally before warmups. Do not evict host caches

set -euo pipefail
export LC_ALL=C
source "${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks/common.sh"
BENCH_FEATURE=lengths

: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${BENCH_BAM:?Set the prepared input BAM path}"
: "${BENCH_LOCAL_ROOT:?Set an existing node-local scratch directory}"
: "${BENCH_RESULTS:?Set a new results directory}"
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM}.bai}
: "${BENCH_CORES:?Set a single core budget for this job: 1, 4, 12 or 32}"
BENCH_WARMUPS=${BENCH_WARMUPS:-2}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-5}
: "${BENCH_METHODS:?Set the methods for this job explicitly in BENCH_METHODS}"
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-30}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-1000}
export BENCH_GENOME=${BENCH_GENOME:-hg38}

benchmark_controls
read -r -a methods_to_run <<< "$BENCH_METHODS"
if [[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 || -e "$BENCH_RESULTS" ]]; then
    echo "Use a single-node allocation and a new results directory." >&2
    exit 1
fi
if [[ "$BENCH_INDEX" != *.bai || ${#methods_to_run[@]} == 0 ]]; then
    echo "Supply a BAI index, at least a method and at least a core budget." >&2
    exit 1
fi
for value in "$BENCH_WARMUPS" "$BENCH_REPETITIONS" "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH"; do
    if [[ ! "$value" =~ ^(0|[1-9][0-9]*)$ ]]; then
        echo "Core counts, repetitions and filters must be integers: $value" >&2
        exit 1
    fi
done
if ((BENCH_REPETITIONS < 1 || BENCH_MIN_MAPQ > 255 || BENCH_MIN_LENGTH < 10 || BENCH_MAX_LENGTH < BENCH_MIN_LENGTH || BENCH_MAX_LENGTH > 50000)); then
    echo "Use MAPQ 0-255, lengths within 10-50000 bp, and at least a measured repetition." >&2
    exit 1
fi
if ((BENCH_CORES > ${SLURM_CPUS_PER_TASK:-0})); then
    echo "Core budget $BENCH_CORES exceeds allocated CPUs: ${SLURM_CPUS_PER_TASK:-0}" >&2
    exit 1
fi
for method in "${methods_to_run[@]}"; do
    case "$method" in
        cfdnalab|cfdnalab-autosomes|cfdnalab-proper|cfdnapro|cfdnapro-outer-span|finaletoolkit|samtools|picard) ;;
        *) echo "Unknown length benchmark method: $method" >&2; exit 1 ;;
    esac
done

run_cfdna=false
run_finale=false
run_cfdnapro=false
run_picard=false
for method in "${methods_to_run[@]}"; do
    case "$method" in
        cfdnalab*) run_cfdna=true ;;
        finaletoolkit) run_finale=true ;;
        cfdnapro*) run_cfdnapro=true ;;
        picard) run_picard=true ;;
    esac
done
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the packed samtools environment for BAM header preparation}"
if [[ "$run_cfdna" == true ]]; then : "${BENCH_CFDNA:?Set the cfDNAlab executable}"; fi
if [[ "$run_finale" == true ]]; then : "${BENCH_FINALE_ARCHIVE:?Set the packed FinaleToolkit environment}"; fi
if [[ "$run_cfdnapro" == true ]]; then : "${BENCH_CFDNAPRO_IMAGE:?Set the cfDNAPro image}"; fi
if [[ "$run_picard" == true ]]; then : "${BENCH_PICARD_ARCHIVE:?Set the packed Picard environment}"; fi

# sbatch runs a spool copy of this script. Submit from the repository root so
# the R helper can be staged from SLURM_SUBMIT_DIR
export BENCH_WORK BENCH_NODE BENCH_RESULTS BENCH_AUTOSOMES BENCH_APPTAINER
export BENCH_CFDNA_VERSION BENCH_CFDNAPRO_VERSION BENCH_FINALE_VERSION BENCH_SAMTOOLS_VERSION BENCH_PICARD_VERSION
BENCH_WORK=$(mktemp -d "$BENCH_LOCAL_ROOT/lengths-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)
if [[ "$run_cfdnapro" == true ]]; then
    BENCH_APPTAINER=$(command -v apptainer) || { echo "cfDNAPro requires Apptainer." >&2; exit 1; }
fi
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.bai"
if [[ "$run_cfdna" == true ]]; then benchmark_stage_binary "$BENCH_CFDNA" cfdna; fi
if [[ "$run_cfdnapro" == true ]]; then cp -L "$BENCH_CFDNAPRO_IMAGE" "$BENCH_WORK/cfdnapro.sif"; fi
cp "$SLURM_SUBMIT_DIR/benchmarks/lengths/cfdnapro_lengths.R" "$BENCH_WORK/cfdnapro_lengths.R"
export TMPDIR="$BENCH_WORK/tmp" APPTAINER_TMPDIR="$BENCH_WORK/tmp"
export XDG_CACHE_HOME="$BENCH_WORK/cache" MPLCONFIGDIR="$BENCH_WORK/matplotlib"
# Ignore ~/.local Python packages and inherited import paths so imports stay
# inside the staged environment. Clear injected Java options for the same reason
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
# Picard's packaged launcher honors JAVA_HOME. Select its node-local Java
export JAVA_HOME="$BENCH_WORK/picard-env"

# Stage selected tools and reuse samtools for untimed BAM header preparation
benchmark_unpack "$BENCH_SAMTOOLS_ARCHIVE" "$BENCH_WORK/samtools-env"
if [[ "$run_finale" == true ]]; then benchmark_unpack "$BENCH_FINALE_ARCHIVE" "$BENCH_WORK/finale-env"; fi
if [[ "$run_picard" == true ]]; then benchmark_unpack "$BENCH_PICARD_ARCHIVE" "$BENCH_WORK/picard-env"; fi
"$BENCH_WORK/samtools-env/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/chrom.sizes"
BENCH_AUTOSOMES=$(awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ {print $1}' "$BENCH_WORK/chrom.sizes" | paste -sd, -)
if [[ " $BENCH_METHODS " == *" cfdnalab-autosomes "* && -z "$BENCH_AUTOSOMES" ]]; then echo "No autosomes in BAM." >&2; exit 1; fi
if [[ "$run_cfdna" == true ]]; then BENCH_CFDNA_VERSION=$("$BENCH_WORK/bin/cfdna" --version | tr '\t\n' '  '); fi
if [[ "$run_finale" == true ]]; then BENCH_FINALE_VERSION=$("$BENCH_WORK/finale-env/bin/finaletoolkit" --version | tr '\t\n' '  '); fi
BENCH_SAMTOOLS_VERSION=$("$BENCH_WORK/samtools-env/bin/samtools" --version-only | tr '\t\n' '  ')
if [[ "$run_picard" == true ]]; then BENCH_PICARD_VERSION=$(benchmark_picard_version "$BENCH_WORK/picard-env/bin/picard"); fi
if [[ "$run_cfdnapro" == true ]]; then
    BENCH_CFDNAPRO_VERSION=$("$BENCH_APPTAINER" exec --cleanenv --containall "$BENCH_WORK/cfdnapro.sif" Rscript --vanilla -e 'cat(as.character(packageVersion("cfDNAPro")))')
    image_hash=$(sha256sum "$BENCH_WORK/cfdnapro.sif" | cut -d ' ' -f 1)
    BENCH_CFDNAPRO_VERSION+=" image-sha256:$image_hash"
fi

mkdir -p "$BENCH_RESULTS/profiles"
header=$'run_order\trepetition\twarmup\tmethod\ttool_version\tcores\ttool_threads'
header+=$'\tmin_mapq\tmin_length\tmax_length\tscope\tlength_definition\tproper_pair'
header+=$'\tnode\tcpu_model\tcpu_ids\tphysical_core_ids\tjob_id\tstep_id\tstarted_utc'
header+=$'\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kib\ttool_exit_code'
header+=$'\tflush_seconds\tprocess_plus_flush_seconds\toutput_bytes\tpilot\tstep_exit_code'
printf '%s\n' "$header" > "$BENCH_RESULTS/timings.tsv"

# Run each method's warmups and measurements together on the recorded node
benchmark_record_setup
benchmark_schedule > "$BENCH_WORK/order.tsv"
if [[ ! -s "$BENCH_WORK/order.tsv" ]]; then echo "The requested run matrix is empty." >&2; exit 1; fi
run_order=0
job_exit_code=0
failed_method=
while IFS=$'\t' read -r method cores repetition is_warmup; do
    if [[ "$method" == "$failed_method" ]]; then continue; fi
    run_order=$((run_order + 1))
    step_exit_code=0
    rm -f "$BENCH_WORK/time.tsv" "$BENCH_WORK/io.tsv" "$BENCH_WORK/profile.path" "$BENCH_WORK/step.id" "$BENCH_WORK/stdout.path"

    # Each step attaches to this existing single-node allocation
    # Check its actual node and CPU mask below. Hard confinement also depends
    # on the cluster's Slurm cgroup configuration
    srun --jobid="$SLURM_JOB_ID" --exact --nodes=1 --ntasks=1 \
        --cpus-per-task="$cores" --threads-per-core=1 --cpu-bind=threads \
        --export=ALL bash -s -- "$run_order" "$repetition" "$is_warmup" "$method" "$cores" \
        <<'STEP' || step_exit_code=$?
set -euo pipefail
source "$BENCH_WORK/common.sh"
printf '%s\n' "$SLURM_STEP_ID" > "$BENCH_WORK/step.id"
run_order=$1
repetition=$2
is_warmup=$3
method=$4
cores=$5

# Allow numerical libraries the full allocation even when a tool has no thread
# flag. CPU binding still limits the complete process tree to this step's cores
export OMP_NUM_THREADS="$cores" OPENBLAS_NUM_THREADS="$cores"
export MKL_NUM_THREADS="$cores" NUMEXPR_NUM_THREADS="$cores"

node=$(hostname)
cpu_model=$(awk -F ': ' '/^model name/ {print $2; exit}' /proc/cpuinfo)
cpu_ids=$(awk '/Cpus_allowed_list/ {print $2}' /proc/self/status)
started_utc=$(date -u +%FT%TZ)
physical_cores=()
for cpu_range in ${cpu_ids//,/ }; do
    first_cpu=${cpu_range%-*}
    last_cpu=${cpu_range#*-}
    for ((cpu_id=first_cpu; cpu_id<=last_cpu; cpu_id++)); do
        topology="/sys/devices/system/cpu/cpu${cpu_id}/topology"
        socket_id=$(cat "$topology/physical_package_id")
        core_id=$(cat "$topology/core_id")
        physical_cores+=("$socket_id:$core_id")
    done
done
physical_core_ids=$(printf '%s\n' "${physical_cores[@]}" | sort -u | paste -sd, -)
physical_core_count=$(printf '%s\n' "${physical_cores[@]}" | sort -u | wc -l)

mkdir -p "$BENCH_WORK/output"
profile="$BENCH_WORK/output/lengths.tsv"
scope=all
length_definition=abs-TLEN
proper_pair=not-required
tool_threads=NA
case "$method" in
    cfdnalab|cfdnalab-autosomes|cfdnalab-proper)
        tool_version=$BENCH_CFDNA_VERSION
        tool_threads=$cores
        length_definition=directional-aligned-span
        command=("$BENCH_WORK/bin/cfdna" lengths
            --bam "$BENCH_WORK/sample.bam" --output-dir "$BENCH_WORK/output"
            --length-bins "$BENCH_MIN_LENGTH:$((BENCH_MAX_LENGTH + 1)):1"
            --indel-mode ignore --clip-mode aligned --min-mapq "$BENCH_MIN_MAPQ"
            --n-threads "$cores")
        if [[ "$method" == cfdnalab-autosomes ]]; then
            scope=autosomes
            command+=(--chromosomes "$BENCH_AUTOSOMES")
        else
            command+=(--chromosomes all)
        fi
        if [[ "$method" == cfdnalab-proper ]]; then
            proper_pair=required
            command+=(--require-proper-pair)
        fi
        profile="$BENCH_WORK/output/length_counts.tsv.zst"
        ;;
    cfdnapro|cfdnapro-outer-span)
        tool_version=$BENCH_CFDNAPRO_VERSION
        length_definition=directional-aligned-span
        curate=TRUE
        if [[ "$method" == cfdnapro-outer-span ]]; then
            length_definition=outer-aligned-span
            curate=FALSE
        fi
        # The image contains R, cfDNAPro and the authors' reference packages
        # Containment and --vanilla avoid reading home-directory R settings
        # Pass the full library thread budget because --cleanenv drops host settings
        command=("$BENCH_APPTAINER" exec --cleanenv --containall
            --bind "$BENCH_WORK:/work" --pwd /work
            --env "OMP_NUM_THREADS=$cores,OPENBLAS_NUM_THREADS=$cores,MKL_NUM_THREADS=$cores"
            "$BENCH_WORK/cfdnapro.sif" Rscript --vanilla /work/cfdnapro_lengths.R
            /work/sample.bam /work/output/lengths.tsv "$BENCH_MIN_LENGTH"
            "$BENCH_MAX_LENGTH" "$BENCH_MIN_MAPQ" "$BENCH_GENOME" "$curate")
        ;;
    finaletoolkit)
        tool_version=$BENCH_FINALE_VERSION
        proper_pair=required
        command=("$BENCH_WORK/finale-env/bin/finaletoolkit" frag-length-bins
            "$BENCH_WORK/sample.bam" --bin-size 1 --min-length "$BENCH_MIN_LENGTH"
            --max-length "$BENCH_MAX_LENGTH" --min-mapq "$BENCH_MIN_MAPQ" --output "$profile")
        ;;
    samtools)
        tool_version=$BENCH_SAMTOOLS_VERSION
        # -@ counts additional decoding threads, leaving a core for calculation
        # The HTSlib input filter applies MAPQ without an extra BAM/filtering pass
        # -m 1 prevents tail truncation. The final bin collects overflow, so place
        # it above our range. The comparator selects inward pairs after timing
        tool_threads=$((cores - 1))
        profile="$BENCH_WORK/output/stats.txt"
        command=("$BENCH_WORK/samtools-env/bin/samtools" stats -@ "$tool_threads"
            --input-fmt-option "filter=mapq >= $BENCH_MIN_MAPQ"
            -f 1 -F 3852 -m 1 -i "$((BENCH_MAX_LENGTH + 1))" "$BENCH_WORK/sample.bam")
        ;;
    picard)
        tool_version=$BENCH_PICARD_VERSION
        profile="$BENCH_WORK/output/metrics.txt"
        # Picard has no MAPQ option for this command. Include its input filter
        # in timing, streaming uncompressed BAM directly into the native command
        # Both processes share the step's binding. Omit PDF plotting
        # Explicit width prevents Picard from trimming the histogram tail
        command=(bash -c '
            set -euo pipefail
            cores=$1
            "$BENCH_WORK/samtools-env/bin/samtools" view -@ "$((cores - 1))" -u -q "$BENCH_MIN_MAPQ" -f 1 -F 3852 "$BENCH_WORK/sample.bam" |
                "$BENCH_WORK/picard-env/bin/picard" -Xmx4g -XX:ActiveProcessorCount="$cores" \
                -Djava.io.tmpdir="$BENCH_WORK/tmp" \
                CollectInsertSizeMetrics \
                INPUT=/dev/stdin OUTPUT="$BENCH_WORK/output/metrics.txt" \
                HISTOGRAM_WIDTH="$BENCH_MAX_LENGTH" MINIMUM_PCT=0 INCLUDE_DUPLICATES=false
        ' -- "$cores")
        ;;
esac

# Twenty identifying fields precede GNU time's five measurement fields
printf -v prefix '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "$run_order" "$repetition" "$is_warmup" "$method" "$tool_version" "$cores" "$tool_threads" \
    "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" "$scope" "$length_definition" \
    "$proper_pair" "$node" "$cpu_model" "$cpu_ids" "$physical_core_ids" "$SLURM_JOB_ID" "$SLURM_STEP_ID" "$started_utc"
logical_cpu_count=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
if [[ "$node" != "$BENCH_NODE" ]] || ((logical_cpu_count != cores || physical_core_count != cores)); then
    printf '%s\tNA\tNA\tNA\tNA\t1\n' "$prefix" > "$BENCH_WORK/time.tsv"
    echo "Unexpected placement: node=$node, CPUs=$cpu_ids, physical cores=$physical_core_ids" >&2
    exit 1
fi

# samtools prints its full native report to stdout. The other commands write
# their profiles directly. Keep success logs temporary and print failures
tool_stdout="$BENCH_WORK/tool.log"
if [[ "$method" == samtools ]]; then tool_stdout=$profile; fi
tool_exit_code=0
benchmark_time_command "$prefix" "$BENCH_WORK/output" \
    "$tool_stdout" "$BENCH_WORK/tool-error.log" "${command[@]}" || tool_exit_code=$?
if ((tool_exit_code != 0)); then
    cat "$BENCH_WORK/tool-error.log" "$tool_stdout" >&2
    exit "$tool_exit_code"
fi
if [[ ! -s "$profile" ]]; then
    cat "$BENCH_WORK/tool-error.log" "$tool_stdout" >&2
    echo "The expected length profile is missing or empty: $profile" >&2
    exit 1
fi
printf '%s\n' "$profile" > "$BENCH_WORK/profile.path"
STEP

    # Retain known identifiers even when Slurm fails before GNU time can run
    if [[ ! -s "$BENCH_WORK/time.tsv" ]]; then
        fallback=()
        for ((column=1; column<=25; column++)); do fallback+=(NA); done
        fallback[0]=$run_order
        fallback[1]=$repetition
        fallback[2]=$is_warmup
        fallback[3]=$method
        fallback[5]=$cores
        fallback[17]=$SLURM_JOB_ID
        if [[ -s "$BENCH_WORK/step.id" ]]; then fallback[18]=$(cat "$BENCH_WORK/step.id"); fi
        (IFS=$'\t'; printf '%s\n' "${fallback[*]}") > "$BENCH_WORK/time.tsv"
    fi
    benchmark_finish_run "$method" "$cores" "$repetition" "$is_warmup" "$step_exit_code"
    if ((step_exit_code != 0)); then
        failed_method=$method
        job_exit_code=1
        echo "Skipping remaining repetitions of $method after failure; continuing other selected methods." >&2
    fi
done < "$BENCH_WORK/order.tsv"

# Keep successful methods usable while marking an incomplete job as failed
exit "$job_exit_code"
