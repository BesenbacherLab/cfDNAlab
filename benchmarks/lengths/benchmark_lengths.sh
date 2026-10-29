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

: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${BENCH_BAM:?Set the prepared input BAM path}"
: "${BENCH_CFDNA:?Set the cfDNAlab executable path}"
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the packed samtools environment}"
: "${BENCH_FINALE_ARCHIVE:?Set the packed FinaleToolkit environment}"
: "${BENCH_CFDNAPRO_IMAGE:?Set the cfDNAPro Apptainer image from its authors}"
: "${BENCH_PICARD_ARCHIVE:?Set the packed Picard/Java environment}"
: "${BENCH_LOCAL_ROOT:?Set an existing node-local scratch directory}"
: "${BENCH_RESULTS:?Set a new results directory}"
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM}.bai}
BENCH_CORES=${BENCH_CORES:-1 4 8 12 32}
BENCH_WARMUPS=${BENCH_WARMUPS:-3}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-10}
BENCH_METHODS=${BENCH_METHODS:-cfdnalab cfdnalab-autosomes cfdnapro finaletoolkit samtools picard}
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-30}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-1000}
export BENCH_GENOME=${BENCH_GENOME:-hg38}

read -r -a cores_to_run <<< "$BENCH_CORES"
read -r -a methods_to_run <<< "$BENCH_METHODS"
if [[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 || -e "$BENCH_RESULTS" ]]; then
    echo "Use a single-node allocation and a new results directory." >&2
    exit 1
fi
if [[ "$BENCH_INDEX" != *.bai || ${#cores_to_run[@]} == 0 || ${#methods_to_run[@]} == 0 ]]; then
    echo "Supply a BAI index, at least a method and at least a core budget." >&2
    exit 1
fi
for value in "$BENCH_WARMUPS" "$BENCH_REPETITIONS" "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" "${cores_to_run[@]}"; do
    if [[ ! "$value" =~ ^[0-9]+$ ]]; then
        echo "Core counts, repetitions and filters must be integers: $value" >&2
        exit 1
    fi
done
if ((BENCH_REPETITIONS < 1 || BENCH_MIN_MAPQ > 255 || BENCH_MIN_LENGTH < 10 || BENCH_MAX_LENGTH < BENCH_MIN_LENGTH || BENCH_MAX_LENGTH > 50000)); then
    echo "Use MAPQ 0-255, lengths within 10-50000 bp, and at least a measured repetition." >&2
    exit 1
fi
for cores in "${cores_to_run[@]}"; do
    if ((cores < 1 || cores > ${SLURM_CPUS_PER_TASK:-0})); then
        echo "Core budget $cores is outside this job's CPUs per task." >&2
        exit 1
    fi
done
for method in "${methods_to_run[@]}"; do
    case "$method" in
        cfdnalab|cfdnalab-autosomes|cfdnalab-proper|cfdnapro|cfdnapro-outer-span|finaletoolkit|samtools|picard) ;;
        *) echo "Unknown length benchmark method: $method" >&2; exit 1 ;;
    esac
done

# sbatch runs a spool copy of this script. Submit from the repository root so
# the R helper and memory sampler can be staged from SLURM_SUBMIT_DIR
export BENCH_WORK BENCH_NODE BENCH_RESULTS BENCH_AUTOSOMES BENCH_CAPTURE_CORES BENCH_APPTAINER
export BENCH_CFDNA_VERSION BENCH_CFDNAPRO_VERSION BENCH_FINALE_VERSION BENCH_SAMTOOLS_VERSION BENCH_PICARD_VERSION
BENCH_WORK=$(mktemp -d "$BENCH_LOCAL_ROOT/lengths-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)
BENCH_CAPTURE_CORES=${cores_to_run[0]}
BENCH_APPTAINER=$(command -v apptainer)
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.bai"
cp -L "$BENCH_CFDNA" "$BENCH_WORK/bin/cfdna"
cp -L "$BENCH_CFDNAPRO_IMAGE" "$BENCH_WORK/cfdnapro.sif"
cp "$SLURM_SUBMIT_DIR/benchmarks/lengths/cfdnapro_lengths.R" "$BENCH_WORK/cfdnapro_lengths.R"
cp "$SLURM_SUBMIT_DIR/benchmarks/measure_memory.py" "$BENCH_WORK/measure_memory.py"
export TMPDIR="$BENCH_WORK/tmp" APPTAINER_TMPDIR="$BENCH_WORK/tmp"
export XDG_CACHE_HOME="$BENCH_WORK/cache" MPLCONFIGDIR="$BENCH_WORK/matplotlib"
# Ignore ~/.local Python packages and inherited import paths so imports stay
# inside the staged environment. Clear injected Java options for the same reason
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME JAVA_TOOL_OPTIONS JDK_JAVA_OPTIONS _JAVA_OPTIONS
# Picard's packaged launcher honors JAVA_HOME. Select its node-local Java
export JAVA_HOME="$BENCH_WORK/picard-env"

# Relocate complete environments before timing, rather than copying their launchers
mkdir -p "$BENCH_WORK/finale-env" "$BENCH_WORK/picard-env" "$BENCH_WORK/samtools-env"
tar -xzf "$BENCH_FINALE_ARCHIVE" -C "$BENCH_WORK/finale-env"
"$BENCH_WORK/finale-env/bin/python" "$BENCH_WORK/finale-env/bin/conda-unpack"
tar -xzf "$BENCH_PICARD_ARCHIVE" -C "$BENCH_WORK/picard-env"
# Relocation uses only Python's standard library, so neither Java nor samtools
# needs its own Python installation. Reuse FinaleToolkit's staged Python
"$BENCH_WORK/finale-env/bin/python" "$BENCH_WORK/picard-env/bin/conda-unpack"
tar -xzf "$BENCH_SAMTOOLS_ARCHIVE" -C "$BENCH_WORK/samtools-env"
"$BENCH_WORK/finale-env/bin/python" "$BENCH_WORK/samtools-env/bin/conda-unpack"
BENCH_AUTOSOMES=$("$BENCH_WORK/finale-env/bin/python" - <<'PY'
import os
import re
import pysam

with pysam.AlignmentFile(os.environ['BENCH_WORK'] + '/sample.bam', 'rb') as bam:
    autosomes = [name for name in bam.references if re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2])', name)]
if not autosomes:
    raise SystemExit('No autosomes named chr1..chr22 or 1..22 were found.')
print(','.join(autosomes))
PY
)

BENCH_CFDNA_VERSION=$("$BENCH_WORK/bin/cfdna" --version 2>&1 | tr '\t\n' '  ')
BENCH_FINALE_VERSION=$("$BENCH_WORK/finale-env/bin/finaletoolkit" --version 2>&1 | tr '\t\n' '  ')
BENCH_SAMTOOLS_VERSION=$("$BENCH_WORK/samtools-env/bin/samtools" --version | sed -n '1p')
BENCH_PICARD_VERSION=$("$BENCH_WORK/picard-env/bin/picard" CollectInsertSizeMetrics --version 2>&1 | tr '\t\n' '  ')
BENCH_CFDNAPRO_VERSION=$("$BENCH_APPTAINER" exec --cleanenv --containall "$BENCH_WORK/cfdnapro.sif" \
    Rscript --vanilla -e 'cat(as.character(packageVersion("cfDNAPro")))')
# The image hash identifies the actual container, even if its upstream tag changes
image_hash=$(sha256sum "$BENCH_WORK/cfdnapro.sif" | cut -d ' ' -f 1)
BENCH_CFDNAPRO_VERSION+=" image-sha256:$image_hash"

mkdir -p "$BENCH_RESULTS/profiles"
header=$'run_order\trepetition\twarmup\tmethod\ttool_version\tcores\ttool_threads'
header+=$'\tmin_mapq\tmin_length\tmax_length\tscope\tlength_definition\tproper_pair'
header+=$'\tnode\tcpu_model\tcpu_ids\tphysical_core_ids\tjob_id\tstep_id\tstarted_utc'
header+=$'\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kib\ttool_exit_code'
header+=$'\tsampled_peak_rss_kib\tsampled_mean_rss_kib\tmemory_samples\tmemory_interval_seconds\tstep_exit_code'
printf '%s\n' "$header" > "$BENCH_RESULTS/timings.tsv"

# Give every tool the same core budgets, including commands without a thread flag
# Their CPU times will show whether they actually use the allocated cores
for method in "${methods_to_run[@]}"; do
    for cores in "${cores_to_run[@]}"; do printf '%s\t%s\n' "$method" "$cores"; done
done > "$BENCH_WORK/combinations.tsv"

run_order=0
total_repetitions=$((BENCH_WARMUPS + BENCH_REPETITIONS))
for ((repetition=1; repetition<=total_repetitions; repetition++)); do
    is_warmup=false
    if ((repetition <= BENCH_WARMUPS)); then is_warmup=true; fi
    shuf "$BENCH_WORK/combinations.tsv" > "$BENCH_WORK/order.tsv"
    while IFS=$'\t' read -r method cores; do
        run_order=$((run_order + 1))
        step_exit_code=0
        rm -f "$BENCH_WORK/time.tsv" "$BENCH_WORK/memory.tsv"

        # Each step attaches to this existing single-node allocation
        # Check its actual node and CPU mask below. Hard confinement also depends
        # on the cluster's Slurm cgroup configuration
        srun --jobid="$SLURM_JOB_ID" --exact --nodes=1 --ntasks=1 \
            --cpus-per-task="$cores" --threads-per-core=1 --cpu-bind=threads \
            --export=ALL bash -s -- "$run_order" "$repetition" "$is_warmup" "$method" "$cores" \
            <<'STEP' || step_exit_code=$?
set -euo pipefail
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
"$BENCH_WORK/finale-env/bin/python" "$BENCH_WORK/measure_memory.py" "$BENCH_WORK/memory.tsv" \
    /usr/bin/time --quiet --format "$prefix\t%e\t%U\t%S\t%M\t%x" \
    --output "$BENCH_WORK/time.tsv" "${command[@]}" \
    > "$tool_stdout" 2> "$BENCH_WORK/tool-error.log" || tool_exit_code=$?
if ((tool_exit_code != 0)); then
    cat "$BENCH_WORK/tool-error.log" "$tool_stdout" >&2
    exit "$tool_exit_code"
fi
if [[ ! -s "$profile" ]]; then
    cat "$BENCH_WORK/tool-error.log" "$tool_stdout" >&2
    echo "The expected length profile is missing or empty: $profile" >&2
    exit 1
fi
if ((repetition == 1 && cores == BENCH_CAPTURE_CORES)); then
    mkdir -p "$BENCH_RESULTS/profiles/$method"
    cp "$profile" "$BENCH_RESULTS/profiles/$method/"
fi
STEP

        # Append immediately so completed runs survive a later failure
        # Scheduler failures without timing measurements also receive a row
        if [[ ! -s "$BENCH_WORK/time.tsv" ]]; then
            printf '%s\t%s\t%s\t%s\tNA\t%s' \
                "$run_order" "$repetition" "$is_warmup" "$method" "$cores" > "$BENCH_WORK/time.tsv"
            for ((column=7; column<=25; column++)); do printf '\tNA' >> "$BENCH_WORK/time.tsv"; done
            printf '\n' >> "$BENCH_WORK/time.tsv"
        fi
        memory_values=$'NA\tNA\tNA\tNA'
        if [[ -s "$BENCH_WORK/memory.tsv" ]]; then memory_values=$(cat "$BENCH_WORK/memory.tsv"); fi
        printf '%s\t%s\t%s\n' "$(cat "$BENCH_WORK/time.tsv")" "$memory_values" "$step_exit_code" >> "$BENCH_RESULTS/timings.tsv"
        if ((step_exit_code != 0)); then
            echo "Run $run_order failed. Completed rows are in $BENCH_RESULTS/timings.tsv" >&2
            exit "$step_exit_code"
        fi
        sync
        rm -rf "$BENCH_WORK/output" "$BENCH_WORK/tool.log" "$BENCH_WORK/tool-error.log"
    done < "$BENCH_WORK/order.tsv"
done
