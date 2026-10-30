#!/usr/bin/env bash
#SBATCH --job-name=cfdna-coverage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --no-requeue

# Reserve the whole node for this job's method blocks: --exclusive requests all
# CPUs, and --mem=0 requests all node memory. The 32 CPUs per task ensure the node
# supports the largest run. Each srun binds this job's single selected core count
# Use a partition that honors exclusive allocations, as partition policy takes
# precedence over the job's sharing request
# https://slurm.schedmd.com/sbatch.html#OPT_exclusive
#
# Results contain timings.tsv, retained profiles and setup provenance
# Keep the first measured profile for concordance. Other outputs are temporary
# Diagnostics from a failed command go to the normal Slurm job log
# cfdnalab, mosdepth and samtools-depth process all BAM contigs for per-base depth
# cfdnalab-autosomes selects autosomes named chr1..chr22 or 1..22
# This preserves the default biological scope for either BAM naming convention
# pandepth also processes all BAM contigs, counting overlapping mates twice
# Its workload is labeled per-base-read-depth to retain that distinction
# cfdnalab-region-counts and finaletoolkit-coverage count fragments overlapping
# identical 500 kb windows over all BAM contigs. The workload column separates
# these regional counts from per-base depth. They measure different quantities
#
# Example submission from the repository root:
#   bash benchmarks/install_tools.sh /path/to/benchmark-tools
#   export BENCH_BAM=/faststorage/project/sample.bam
#   export BENCH_CFDNA=/path/to/cfdna
#   export BENCH_MOSDEPTH=/path/to/mosdepth
#   export BENCH_SAMTOOLS_ARCHIVE=/path/to/benchmark-tools/samtools.tar.gz
#   export BENCH_PANDEPTH=/path/to/pandepth
#   export BENCH_FINALE_ARCHIVE=/path/to/benchmark-tools/finaletoolkit.tar.gz
#   export BENCH_LOCAL_ROOT=/path/to/node-local-scratch
#   export BENCH_RESULTS=/faststorage/project/coverage-results
#   export BENCH_CORES=4
#   export BENCH_METHODS="cfdnalab mosdepth"
#   sbatch --constraint=gen6 --exclusive --mem=0 --time=<pilot-based-budget> \
#     benchmarks/coverage/benchmark_coverage.sh
# Adjust the generation selector, partition, time and paths for the cluster

# Stop on errors and unset variables. Use a dot for decimal numbers in the TSV
set -euo pipefail
export LC_ALL=C
source "${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks/common.sh"
BENCH_FEATURE=coverage

# These paths are supplied in the environment when submitting the job
# BENCH_INDEX can override the usual sample.bam.bai path. FinaleToolkit requires
# a BAI index, so CSI indices are not supported by this combined benchmark
# The native executables must work with node-local or system shared libraries
# BENCH_FINALE_ARCHIVE contains a Linux Conda environment with FinaleToolkit
# installed, packed with conda-pack so Python and its packages can run locally
# BENCH_SAMTOOLS_ARCHIVE similarly includes samtools and its shared libraries
: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${BENCH_BAM:?Set the input BAM path}"
: "${BENCH_LOCAL_ROOT:?Set an existing node-local scratch directory}"
: "${BENCH_RESULTS:?Set a new results directory}"
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the path of the packed samtools environment}"
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM}.bai}
: "${BENCH_CORES:?Set a single core budget for this job: 1, 4, 12 or 32}"
BENCH_WARMUPS=${BENCH_WARMUPS:-2}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-5}
: "${BENCH_METHODS:?Set the methods for this job explicitly in BENCH_METHODS}"
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-30}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-1000}
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
    if [[ ! "$value" =~ ^(0|[1-9][0-9]*)$ ]]; then echo "Invalid integer setting: $value" >&2; exit 1; fi
done
if ((BENCH_REPETITIONS < 1 || BENCH_MIN_MAPQ > 255 || BENCH_MIN_LENGTH < 10 || BENCH_MAX_LENGTH < BENCH_MIN_LENGTH || BENCH_MAX_LENGTH > 50000)); then
    echo "Use MAPQ 0-255, lengths within 10-50000 bp and at least a measured repetition." >&2
    exit 1
fi
if ((BENCH_CORES > ${SLURM_CPUS_PER_TASK:-0})); then
    echo "Core budget $BENCH_CORES exceeds allocated CPUs: ${SLURM_CPUS_PER_TASK:-0}" >&2
    exit 1
fi
for method in "${methods_to_run[@]}"; do
    case "$method" in
        cfdnalab|cfdnalab-autosomes|cfdnalab-fragment|cfdnalab-fragment-proper|mosdepth|mosdepth-fragment|samtools-depth|pandepth|cfdnalab-region-counts|cfdnalab-region-counts-proper|finaletoolkit-coverage) ;;
        *) echo "Unknown coverage method: $method" >&2; exit 1 ;;
    esac
done

# Stage the BAM, index and executables once before warmup
# Every srun below uses this job's node and these same local files
export BENCH_WORK BENCH_NODE BENCH_CFDNA_VERSION BENCH_MOSDEPTH_VERSION BENCH_AUTOSOMES
export BENCH_SAMTOOLS_VERSION BENCH_FINALE_VERSION
export BENCH_PANDEPTH_VERSION
BENCH_WORK=$(mktemp -d "$BENCH_LOCAL_ROOT/coverage-${SLURM_JOB_ID}.XXXXXX")
BENCH_NODE=$(hostname)
trap 'rm -rf "$BENCH_WORK"' EXIT
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp"
export TMPDIR="$BENCH_WORK/tmp"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.${BENCH_INDEX##*.}"
# sbatch uses a spool copy of this script. Submit from the repository root

# Stage selected competitors and the samtools header-reading helper
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
export XDG_CACHE_HOME="$BENCH_WORK/cache" MPLCONFIGDIR="$BENCH_WORK/matplotlib"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
benchmark_unpack "$BENCH_SAMTOOLS_ARCHIVE" "$BENCH_WORK/samtools-env"
for method in "${methods_to_run[@]}"; do
    case "$method" in
        cfdnalab*)
            : "${BENCH_CFDNA:?Set the cfDNAlab executable}"
            if [[ ! -f "$BENCH_WORK/bin/cfdna" ]]; then benchmark_stage_binary "$BENCH_CFDNA" cfdna; fi ;;
        mosdepth*)
            : "${BENCH_MOSDEPTH:?Set the mosdepth executable}"
            if [[ ! -f "$BENCH_WORK/bin/mosdepth" ]]; then benchmark_stage_binary "$BENCH_MOSDEPTH" mosdepth; fi ;;
        pandepth)
            : "${BENCH_PANDEPTH:?Set the PanDepth executable}"
            benchmark_stage_binary "$BENCH_PANDEPTH" pandepth ;;
        finaletoolkit-coverage)
            : "${BENCH_FINALE_ARCHIVE:?Set the packed FinaleToolkit environment}"
            benchmark_unpack "$BENCH_FINALE_ARCHIVE" "$BENCH_WORK/finale-env" ;;
    esac
done
"$BENCH_WORK/samtools-env/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/chrom.sizes"

# Generate a shared BED4 of 500 kb windows, including each shorter final window
# idxstats obtains contig lengths from the BAM index. Skip its unmapped '*' row
# This setup is outside every timer. Both regional methods read the same BED
"$BENCH_WORK/samtools-env/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk 'BEGIN {OFS="\t"} $1 != "*" {
        for (window_start=0; window_start<$2; window_start+=500000) {
            window_end=window_start+500000
            if (window_end>$2) window_end=$2
            print $1, window_start, window_end, "."
        }
    }' > "$BENCH_WORK/regions.bed"
if [[ ! -s "$BENCH_WORK/regions.bed" ]]; then
    echo "No coverage windows could be generated from the BAM index." >&2
    exit 1
fi

# Use the BAM's actual autosome names. Snyder's reference can use 1..22, while
# cfDNAlab's default names are chr1..chr22. Keep the file's original names intact
BENCH_AUTOSOMES=$(cut -f 1 "$BENCH_WORK/regions.bed" | uniq |
    awk '/^(chr)?([1-9]|1[0-9]|2[0-2])$/' | paste -sd, -)
if [[ " $BENCH_METHODS " == *" cfdnalab-autosomes "* && -z "$BENCH_AUTOSOMES" ]]; then
    echo "The BAM has no autosomes named chr1..chr22 or 1..22." >&2
    exit 1
fi

# Record only versions of the selected native tools
if [[ -f "$BENCH_WORK/bin/cfdna" ]]; then BENCH_CFDNA_VERSION=$("$BENCH_WORK/bin/cfdna" --version | tr '\t\n' '  '); fi
if [[ -f "$BENCH_WORK/bin/mosdepth" ]]; then BENCH_MOSDEPTH_VERSION=$("$BENCH_WORK/bin/mosdepth" --version | tr '\t\n' '  '); fi
BENCH_SAMTOOLS_VERSION=$("$BENCH_WORK/samtools-env/bin/samtools" --version-only | tr '\t\n' '  ')
if [[ -d "$BENCH_WORK/finale-env" ]]; then BENCH_FINALE_VERSION=$("$BENCH_WORK/finale-env/bin/finaletoolkit" --version | tr '\t\n' '  '); fi
if [[ -f "$BENCH_WORK/bin/pandepth" ]]; then
    BENCH_PANDEPTH_VERSION=$("$BENCH_WORK/bin/pandepth" -h 2>&1 | sed -n 's/.*show this help \[\(v[^]]*\)\].*/\1/p')
    if [[ -z "$BENCH_PANDEPTH_VERSION" ]]; then echo "Could not read PanDepth version." >&2; exit 1; fi
fi
# Native text outputs are large. Preserve an explicitly recorded small scope
if [[ -n ${BENCH_CAPTURE_CHROMOSOMES:-} ]]; then
    awk -v names="$BENCH_CAPTURE_CHROMOSOMES" '
        BEGIN {split(names, entries, ","); for (entry_index in entries) selected[entries[entry_index]]=1}
        $1 in selected {found[$1]=1; print}
        END {for (chromosome in selected) if (!(chromosome in found)) {
            print "Capture chromosome absent from BAM: " chromosome > "/dev/stderr"
            exit 1
        }}' "$BENCH_WORK/chrom.sizes" > "$BENCH_WORK/capture.chrom.sizes"
else
    awk '$1 ~ /^(chr)?(21|22)$/' "$BENCH_WORK/chrom.sizes" > "$BENCH_WORK/capture.chrom.sizes"
    if [[ ! -s "$BENCH_WORK/capture.chrom.sizes" ]]; then head -n 1 "$BENCH_WORK/chrom.sizes" > "$BENCH_WORK/capture.chrom.sizes"; fi
fi
if [[ ! -s "$BENCH_WORK/capture.chrom.sizes" ]]; then echo "No capture chromosomes found." >&2; exit 1; fi
mkdir -p "$BENCH_RESULTS/profiles"
cp "$BENCH_WORK/regions.bed" "$BENCH_RESULTS/regions.bed"
export BENCH_RESULTS BENCH_MIN_MAPQ BENCH_MIN_LENGTH BENCH_MAX_LENGTH

# run_order is global execution order. repetition includes warmups per method/core
# cpu_ids are Linux logical CPU IDs. physical_core_ids are socket:core pairs
# tool_threads means cfDNAlab threads, FinaleToolkit worker processes, or EXTRA
# decompression threads for mosdepth/samtools. It is not the allocated core count
# For PanDepth it records the requested -t setting, within the same CPU binding
# GNU time reports elapsed/user/system seconds, maximum RSS in KiB and exit code
# step_exit_code also catches placement failures or a missing output track
header=$'run_order\trepetition\twarmup\tmethod\tworkload\ttool_version\tcores\ttool_threads'
header+=$'\tnode\tcpu_model\tcpu_ids\tphysical_core_ids\tjob_id\tstep_id\tstarted_utc'
header+=$'\tmin_mapq\tmin_length\tmax_length\tscope\tproper_pair\tlength_filter'
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
    workload=per-base-read-depth
    if [[ "$method" == *fragment* ]]; then workload=per-base-fragment-coverage; fi
    if [[ "$method" == cfdnalab-region-counts* || "$method" == finaletoolkit-coverage ]]; then
        workload=region-counts-500kb
    elif [[ "$method" == pandepth ]]; then
        workload=per-base-read-depth
    fi

    # --jobid attaches to this allocation instead of requesting another node
    # --exact and --cpus-per-task restrict the resources requested by this step
    # Binding uses one hardware thread per core. Hard confinement additionally
    # depends on the cluster's Slurm cgroup configuration
    # The quoted STEP delimiter sends the following Bash block unchanged to
    # the compute node. Its six positional arguments are supplied here
    srun --jobid="$SLURM_JOB_ID" --exact --nodes=1 --ntasks=1 \
        --cpus-per-task="$cores" --threads-per-core=1 --cpu-bind=threads \
        --export=ALL bash -s -- "$run_order" "$repetition" "$is_warmup" "$method" "$cores" "$workload" \
        <<'STEP' || step_exit_code=$?
set -euo pipefail
source "$BENCH_WORK/common.sh"
printf '%s\n' "$SLURM_STEP_ID" > "$BENCH_WORK/step.id"
run_order=$1
repetition=$2
is_warmup=$3
method=$4
cores=$5
workload=$6

# Inspect placement inside the step, where the tool will actually execute
node=$(hostname)
cpu_model=$(awk -F ': ' '/^model name/ {print $2; exit}' /proc/cpuinfo)
cpu_ids=$(awk '/Cpus_allowed_list/ {print $2}' /proc/self/status)
started_utc=$(date -u +%FT%TZ)

# Expand a logical CPU list such as 0-3,8 and look up its physical socket/core IDs
# A single-number range works too, because both parameter expansions return it
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

mkdir -p "$BENCH_WORK/coverage"
scope=all
proper_pair=not-required
length_filter=applied
if [[ "$method" == cfdnalab || "$method" == cfdnalab-autosomes || "$method" == cfdnalab-fragment* ]]; then
    tool_threads=$cores
    tool_version=$BENCH_CFDNA_VERSION
    # Both cfDNAlab variants use the same coverage, filtering and pairing settings
    command=("$BENCH_WORK/bin/cfdna" fcoverage
        --bam "$BENCH_WORK/sample.bam" --output-dir "$BENCH_WORK/coverage"
        --min-mapq "$BENCH_MIN_MAPQ"
        --min-fragment-length "$BENCH_MIN_LENGTH" --max-fragment-length "$BENCH_MAX_LENGTH" --n-threads "$tool_threads")
    # Explicit all matches mosdepth's whole-BAM scope. The autosome variant uses
    # the equivalent autosome selection with this BAM's chromosome names
    if [[ "$method" != cfdnalab-autosomes ]]; then
        command+=(--chromosomes all)
    else
        scope=autosomes
        command+=(--chromosomes "$BENCH_AUTOSOMES")
    fi
    if [[ "$method" != *fragment* ]]; then command+=(--ignore-gap); fi
    if [[ "$method" == cfdnalab-fragment-proper ]]; then command+=(--require-proper-pair); proper_pair=required; fi
    track="$BENCH_WORK/coverage/fcoverage.per_position.bedgraph.zst"
elif [[ "$method" == mosdepth || "$method" == mosdepth-fragment ]]; then
    # mosdepth has a main calculation thread plus extra decompression threads
    # Give decompression the same available budget as the other tools
    tool_threads=$((cores - 1))
    tool_version=$BENCH_MOSDEPTH_VERSION
    # Ordinary mode excludes inter-mate gaps and retains normal CIGAR processing
    # The flag mask removes unmapped, secondary, supplementary, failed and duplicate
    # alignments, plus unmapped mates. Proper-pair selection is not required
    command=("$BENCH_WORK/bin/mosdepth" --threads "$tool_threads"
        --include-flag 1 --flag 3852 --mapq "$BENCH_MIN_MAPQ"
        --min-frag-len "$BENCH_MIN_LENGTH" --max-frag-len "$BENCH_MAX_LENGTH")
    if [[ "$method" == mosdepth-fragment ]]; then
        command+=(--fragment-mode)
        proper_pair=required
    fi
    command+=("$BENCH_WORK/coverage/coverage" "$BENCH_WORK/sample.bam")
    track="$BENCH_WORK/coverage/coverage.per-base.bed.gz"
elif [[ "$method" == samtools-depth ]]; then
    # Leave a core for the main calculation, requesting cores - 1 decoding threads
    # -s removes mate overlap, -aa includes zero-depth bases on all contigs
    # -q is base quality and -Q is MAPQ. No -J means deletions are excluded
    # Match basic flags without imposing proper pairs. samtools depth has no
    # fragment length filter and does not require both mates to pass filtering
    # Its native output is uncompressed text with a 1-based position per row
    # https://www.htslib.org/doc/samtools-depth.html
    tool_threads=$((cores - 1))
    tool_version=$BENCH_SAMTOOLS_VERSION
    length_filter=unsupported
    track="$BENCH_WORK/coverage/samtools.depth.tsv"
    command=("$BENCH_WORK/samtools-env/bin/samtools" depth -@ "$tool_threads"
        -aa -s -q 0 -Q "$BENCH_MIN_MAPQ" -G 3852 --require-flags 1
        -o "$track" "$BENCH_WORK/sample.bam")
elif [[ "$method" == pandepth ]]; then
    # -a writes every base, including zeros, to a compressed SiteDepth table
    # Match MAPQ and exclusion flags. PanDepth has no fragment length filter,
    # paired-read inclusion flag, complete-pair check or mate-overlap correction
    # Overlapping mates therefore count twice, unlike the paired coverage runs
    # Its native position column is 0-based. It also writes chromosome statistics
    # https://github.com/HuiyangYu/PanDepth/blob/main/src/PanDepth.cpp
    tool_threads=$cores
    tool_version=$BENCH_PANDEPTH_VERSION
    length_filter=unsupported
    command=("$BENCH_WORK/bin/pandepth"
        -i "$BENCH_WORK/sample.bam" -o "$BENCH_WORK/coverage/pandepth"
        -a -q "$BENCH_MIN_MAPQ" -x 3852 -t "$tool_threads")
    track="$BENCH_WORK/coverage/pandepth.SiteDepth.gz"
elif [[ "$method" == cfdnalab-region-counts* ]]; then
    # A single length bin gives a fragment count per region, with any overlap
    # Each included fragment contributes 1 to each window it overlaps
    # This matches FinaleToolkit coverage's quantity. fcoverage window means
    # measure base depth instead. Keep cfDNAlab's normal pairing policy
    tool_threads=$cores
    tool_version=$BENCH_CFDNA_VERSION
    command=("$BENCH_WORK/bin/cfdna" lengths
        --bam "$BENCH_WORK/sample.bam" --output-dir "$BENCH_WORK/coverage"
        --chromosomes all --by-bed "$BENCH_WORK/regions.bed" --assign-by any
        --length-bins "$BENCH_MIN_LENGTH" "$((BENCH_MAX_LENGTH + 1))" --indel-mode ignore --clip-mode aligned
        --min-mapq "$BENCH_MIN_MAPQ" --n-threads "$tool_threads")
    if [[ "$method" == cfdnalab-region-counts-proper ]]; then command+=(--require-proper-pair); proper_pair=required; fi
    track="$BENCH_WORK/coverage/length_counts.tsv.zst"
else
    # Use raw fragment counts, with no normalization and scale factor 1
    # --threads controls worker processes in the current FinaleToolkit CLI
    # Its BAM reader requires proper pairs and uses read1/TLEN. Its indexed query
    # also needs read1 to overlap the window, which can miss boundary fragments
    # These differences remain, so assess regional agreement with r/MAE/RMSE
    # https://epifluidlab.github.io/FinaleToolkit/documentation/cli_reference/index.html#coverage
    tool_threads=$cores
    tool_version=$BENCH_FINALE_VERSION
    proper_pair=required
    track="$BENCH_WORK/coverage/finaletoolkit.coverage.bed.gz"
    command=("$BENCH_WORK/finale-env/bin/finaletoolkit" coverage
        "$BENCH_WORK/sample.bam" "$BENCH_WORK/regions.bed"
        --output "$track" --intersect-policy any --scale-factor 1
        --min-length "$BENCH_MIN_LENGTH" --max-length "$BENCH_MAX_LENGTH" --min-mapq "$BENCH_MIN_MAPQ" --threads "$tool_threads")
fi

# This prefix supplies the identifying columns. GNU time supplies the measurements
printf -v prefix '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "$run_order" "$repetition" "$is_warmup" "$method" "$workload" "$tool_version" "$cores" "$tool_threads" \
    "$node" "$cpu_model" "$cpu_ids" "$physical_core_ids" "$SLURM_JOB_ID" "$SLURM_STEP_ID" "$started_utc" \
    "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" "$scope" "$proper_pair" "$length_filter"

# Check both logical and physical counts. Ignore OpenMP settings when using nproc
# A failed placement check still produces a row, with unavailable timings as NA
logical_cpu_count=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
if [[ "$node" != "$BENCH_NODE" ]] || ((logical_cpu_count != cores || physical_core_count != cores)); then
    printf '%s\tNA\tNA\tNA\tNA\t1\n' "$prefix" > "$BENCH_WORK/time.tsv"
    echo "Unexpected placement: node=$node, CPUs=$cpu_ids, physical cores=$physical_core_ids" >&2
    exit 1
fi

# Time only the native command, including creation and compression of its track
# Scheduler startup, placement checks and cleanup are outside this timer
tool_exit_code=0
benchmark_time_command "$prefix" "$BENCH_WORK/coverage" \
    "$BENCH_WORK/tool.log" "$BENCH_WORK/tool-error.log" "${command[@]}" || tool_exit_code=$?
if ((tool_exit_code != 0)); then
    cat "$BENCH_WORK/tool-error.log" "$BENCH_WORK/tool.log" >&2
    exit "$tool_exit_code"
fi
if [[ ! -s "$track" ]]; then
    echo "The command succeeded but its expected track is missing: $track" >&2
    exit 1
fi
printf '%s\n' "$track" > "$BENCH_WORK/profile.path"
STEP

    # Retain known identifiers even when Slurm fails before GNU time can run
    if [[ ! -s "$BENCH_WORK/time.tsv" ]]; then
        fallback=()
        for ((column=1; column<=26; column++)); do fallback+=(NA); done
        fallback[0]=$run_order
        fallback[1]=$repetition
        fallback[2]=$is_warmup
        fallback[3]=$method
        fallback[6]=$cores
        fallback[12]=$SLURM_JOB_ID
        if [[ -s "$BENCH_WORK/step.id" ]]; then fallback[13]=$(cat "$BENCH_WORK/step.id"); fi
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
