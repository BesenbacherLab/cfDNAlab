#!/usr/bin/env bash
#SBATCH --job-name=cfdna-coverage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --no-requeue

# Reserve the whole node for the full benchmark series: --exclusive requests all
# CPUs, and --mem=0 requests all node memory. The 12 CPUs per task ensure the node
# supports the largest run. Each srun below still binds its own 1/4/8/12 CPUs
# Use a partition that honors exclusive allocations, as partition policy takes
# precedence over the job's sharing request
# https://slurm.schedmd.com/sbatch.html#OPT_exclusive
#
# This script writes a single TSV, with a row for every attempted run
# Coverage files and tool logs are temporary and are removed after each run
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
#   conda pack -n finaletoolkit -o /faststorage/project/finaletoolkit.tar.gz
#   export BENCH_BAM=/faststorage/project/sample.bam
#   export BENCH_CFDNA=/path/to/cfdna
#   export BENCH_MOSDEPTH=/path/to/mosdepth
#   export BENCH_SAMTOOLS=/path/to/samtools
#   export BENCH_PANDEPTH=/path/to/pandepth
#   export BENCH_FINALE_ARCHIVE=/faststorage/project/finaletoolkit.tar.gz
#   export BENCH_LOCAL_ROOT=/path/to/node-local-scratch
#   export BENCH_OUTPUT=/faststorage/project/coverage.tsv
#   sbatch --constraint=gen6 --exclusive --mem=0 --time=12:00:00 \
#     benchmarks/coverage/benchmark_coverage.sh
# Adjust the generation selector, partition, time and paths for the cluster

# Stop on errors and unset variables. Use a dot for decimal numbers in the TSV
set -euo pipefail
export LC_ALL=C

# These paths are supplied in the environment when submitting the job
# BENCH_INDEX can override the usual sample.bam.bai path. FinaleToolkit requires
# a BAI index, so CSI indices are not supported by this combined benchmark
# The native executables must work with node-local or system shared libraries
# BENCH_FINALE_ARCHIVE contains a Linux Conda environment with FinaleToolkit
# installed, packed with conda-pack so Python and its packages can run locally
: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${BENCH_BAM:?Set the input BAM path}"
: "${BENCH_LOCAL_ROOT:?Set an existing node-local scratch directory}"
: "${BENCH_OUTPUT:?Set the path of a new output TSV}"
: "${BENCH_FINALE_ARCHIVE:?Set the path of the packed FinaleToolkit environment}"
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM}.bai}
BENCH_CFDNA=${BENCH_CFDNA:-$(command -v cfdna)}
BENCH_MOSDEPTH=${BENCH_MOSDEPTH:-$(command -v mosdepth)}
BENCH_SAMTOOLS=${BENCH_SAMTOOLS:-$(command -v samtools)}
BENCH_PANDEPTH=${BENCH_PANDEPTH:-$(command -v pandepth)}
if [[ "$BENCH_INDEX" != *.bai ]]; then
    echo "FinaleToolkit requires a BAI index. Set BENCH_INDEX to its .bai path." >&2
    exit 1
fi

# Repetition 1-3 is warmup for every method/core combination
# Repetition 4-13 is measured. All rows remain in the TSV for downstream filtering
cores_to_run=(1 4 8 12)
warmups=3
measured_repetitions=10
total_repetitions=$((warmups + measured_repetitions))

if [[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 ]]; then
    echo "The benchmark must run inside a single-node allocation." >&2
    exit 1
fi
if [[ -e "$BENCH_OUTPUT" ]]; then
    echo "Refusing to overwrite an existing TSV: $BENCH_OUTPUT" >&2
    exit 1
fi

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
cp -L "$BENCH_CFDNA" "$BENCH_WORK/bin/cfdna"
cp -L "$BENCH_MOSDEPTH" "$BENCH_WORK/bin/mosdepth"
cp -L "$BENCH_SAMTOOLS" "$BENCH_WORK/bin/samtools"
cp -L "$BENCH_PANDEPTH" "$BENCH_WORK/bin/pandepth"

# Unpack the complete Python environment before timing, fixing its old paths
# Copying only the finaletoolkit launcher would still use the original Python
# See https://conda.github.io/conda-pack/ for the archive preparation procedure
mkdir -p "$BENCH_WORK/finale-env" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib"
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
export XDG_CACHE_HOME="$BENCH_WORK/cache"
export MPLCONFIGDIR="$BENCH_WORK/matplotlib"
# Limit numerical libraries to a thread each. FinaleToolkit's --threads controls
# its worker processes, which share the CPU allocation with their parent
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
tar -xzf "$BENCH_FINALE_ARCHIVE" -C "$BENCH_WORK/finale-env"
"$BENCH_WORK/finale-env/bin/python" "$BENCH_WORK/finale-env/bin/conda-unpack"

# Generate a shared BED4 of 500 kb windows, including each shorter final window
# idxstats obtains contig lengths from the BAM index. Skip its unmapped '*' row
# This setup is outside every timer. Both regional methods read the same BED
"$BENCH_WORK/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
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
if [[ -z "$BENCH_AUTOSOMES" ]]; then
    echo "The BAM has no autosomes named chr1..chr22 or 1..22." >&2
    exit 1
fi

# Record versions in each TSV row, replacing any tabs/newlines with spaces
BENCH_CFDNA_VERSION=$("$BENCH_WORK/bin/cfdna" --version 2>&1 | tr '\t\n' '  ')
BENCH_MOSDEPTH_VERSION=$("$BENCH_WORK/bin/mosdepth" --version 2>&1 | tr '\t\n' '  ')
BENCH_SAMTOOLS_VERSION=$("$BENCH_WORK/bin/samtools" --version-only 2>&1 | tr '\t\n' '  ')
BENCH_FINALE_VERSION=$("$BENCH_WORK/finale-env/bin/finaletoolkit" --version 2>&1 | tr '\t\n' '  ')
# PanDepth reports its version in the -h help text, without a --version option
BENCH_PANDEPTH_VERSION=$("$BENCH_WORK/bin/pandepth" -h 2>&1 |
    sed -n 's/.*show this help \[\(v[^]]*\)\].*/\1/p' | tr '\t\n' '  ')
if [[ -z "$BENCH_PANDEPTH_VERSION" ]]; then
    echo "Could not find the PanDepth version in its -h output." >&2
    exit 1
fi

# run_order is the global execution order. repetition counts 1-13 per combination
# cpu_ids are Linux logical CPU IDs. physical_core_ids are socket:core pairs
# tool_threads means cfDNAlab threads, FinaleToolkit worker processes, or EXTRA
# decompression threads for mosdepth/samtools. It is not the allocated core count
# For PanDepth it records the requested -t setting, within the same CPU binding
# GNU time reports elapsed/user/system seconds, maximum RSS in KiB and exit code
# step_exit_code also catches placement failures or a missing output track
header=$'run_order\trepetition\twarmup\tmethod\tworkload\ttool_version\tcores\ttool_threads'
header+=$'\tnode\tcpu_model\tcpu_ids\tphysical_core_ids\tjob_id\tstep_id\tstarted_utc'
header+=$'\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kib\ttool_exit_code\tstep_exit_code'
printf '%s\n' "$header" > "$BENCH_OUTPUT"

# Shuffle the method/core combinations within each repetition block
# The actual order is recorded in the TSV, including all warmups
for method in cfdnalab cfdnalab-autosomes mosdepth samtools-depth pandepth cfdnalab-region-counts finaletoolkit-coverage; do
    for cores in "${cores_to_run[@]}"; do
        printf '%s\t%s\n' "$method" "$cores"
    done
done > "$BENCH_WORK/combinations.tsv"

run_order=0
for ((repetition=1; repetition<=total_repetitions; repetition++)); do
    is_warmup=false
    if ((repetition <= warmups)); then
        is_warmup=true
    fi
    shuf "$BENCH_WORK/combinations.tsv" > "$BENCH_WORK/order.tsv"

    while IFS=$'\t' read -r method cores; do
        run_order=$((run_order + 1))
        step_exit_code=0
        rm -f "$BENCH_WORK/time.tsv"
        workload=per-base
        if [[ "$method" == cfdnalab-region-counts || "$method" == finaletoolkit-coverage ]]; then
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
if [[ "$method" == cfdnalab || "$method" == cfdnalab-autosomes ]]; then
    tool_threads=$cores
    tool_version=$BENCH_CFDNA_VERSION
    # Both cfDNAlab variants use the same coverage, filtering and pairing settings
    command=("$BENCH_WORK/bin/cfdna" fcoverage
        --bam "$BENCH_WORK/sample.bam" --output-dir "$BENCH_WORK/coverage"
        --ignore-gap --min-mapq 30
        --min-fragment-length 30 --max-fragment-length 1000 --n-threads "$tool_threads")
    # Explicit all matches mosdepth's whole-BAM scope. The autosome variant uses
    # the equivalent autosome selection with this BAM's chromosome names
    if [[ "$method" == cfdnalab ]]; then
        command+=(--chromosomes all)
    else
        command+=(--chromosomes "$BENCH_AUTOSOMES")
    fi
    track="$BENCH_WORK/coverage/fcoverage.per_position.bedgraph.zst"
elif [[ "$method" == mosdepth ]]; then
    # mosdepth has a main calculation thread plus extra decompression threads
    # Use 0,3,4,4 extra threads for 1,4,8,12 cores, capping decompression at four
    tool_threads=$((cores - 1))
    if ((tool_threads > 4)); then
        tool_threads=4
    fi
    tool_version=$BENCH_MOSDEPTH_VERSION
    # Ordinary mode excludes inter-mate gaps and retains normal CIGAR processing
    # The flag mask removes unmapped, secondary, supplementary, failed and duplicate
    # alignments, plus unmapped mates. Proper-pair selection is not required
    command=("$BENCH_WORK/bin/mosdepth" --threads "$tool_threads"
        --include-flag 1 --flag 3852 --mapq 30
        --min-frag-len 30 --max-frag-len 1000
        "$BENCH_WORK/coverage/coverage" "$BENCH_WORK/sample.bam")
    track="$BENCH_WORK/coverage/coverage.per-base.bed.gz"
elif [[ "$method" == samtools-depth ]]; then
    # Leave a core for the main calculation, requesting 0,3,7,11 decoding threads
    # -s removes mate overlap, -aa includes zero-depth bases on all contigs
    # -q is base quality and -Q is MAPQ. No -J means deletions are excluded
    # Match basic flags without imposing proper pairs. samtools depth has no
    # fragment length filter and does not require both mates to pass filtering
    # Its native output is uncompressed text with a 1-based position per row
    # https://www.htslib.org/doc/samtools-depth.html
    tool_threads=$((cores - 1))
    tool_version=$BENCH_SAMTOOLS_VERSION
    track="$BENCH_WORK/coverage/samtools.depth.tsv"
    command=("$BENCH_WORK/bin/samtools" depth -@ "$tool_threads"
        -aa -s -q 0 -Q 30 -G 3852 --require-flags 1
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
    command=("$BENCH_WORK/bin/pandepth"
        -i "$BENCH_WORK/sample.bam" -o "$BENCH_WORK/coverage/pandepth"
        -a -q 30 -x 3852 -t "$tool_threads")
    track="$BENCH_WORK/coverage/pandepth.SiteDepth.gz"
elif [[ "$method" == cfdnalab-region-counts ]]; then
    # A single length bin gives a fragment count per region, with any overlap
    # Each fragment of length 30-1000 bp contributes 1 to each window it overlaps
    # This matches FinaleToolkit coverage's quantity. fcoverage window means
    # measure base depth instead. Keep cfDNAlab's normal pairing policy
    tool_threads=$cores
    tool_version=$BENCH_CFDNA_VERSION
    command=("$BENCH_WORK/bin/cfdna" lengths
        --bam "$BENCH_WORK/sample.bam" --output-dir "$BENCH_WORK/coverage"
        --chromosomes all --by-bed "$BENCH_WORK/regions.bed" --assign-by any
        --length-bins 30 1001 --indel-mode ignore --clip-mode aligned
        --min-mapq 30 --n-threads "$tool_threads")
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
    track="$BENCH_WORK/coverage/finaletoolkit.coverage.bed.gz"
    command=("$BENCH_WORK/finale-env/bin/finaletoolkit" coverage
        "$BENCH_WORK/sample.bam" "$BENCH_WORK/regions.bed"
        --output "$track" --intersect-policy any --scale-factor 1
        --min-length 30 --max-length 1000 --min-mapq 30 --threads "$tool_threads")
fi

# This prefix supplies the identifying columns. GNU time supplies the measurements
prefix="$run_order\t$repetition\t$is_warmup\t$method\t$workload\t$tool_version\t$cores\t$tool_threads"
prefix+="\t$node\t$cpu_model\t$cpu_ids\t$physical_core_ids\t$SLURM_JOB_ID\t$SLURM_STEP_ID\t$started_utc"

# Check both logical and physical counts. Ignore OpenMP settings when using nproc
# A failed placement check still produces a row, with unavailable timings as NA
logical_cpu_count=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
if [[ "$node" != "$BENCH_NODE" ]] || ((logical_cpu_count != cores || physical_core_count != cores)); then
    printf '%b\tNA\tNA\tNA\tNA\t1\n' "$prefix" > "$BENCH_WORK/time.tsv"
    echo "Unexpected placement: node=$node, CPUs=$cpu_ids, physical cores=$physical_core_ids" >&2
    exit 1
fi

# Time only the native command, including creation and compression of its track
# Scheduler startup, placement checks and cleanup are outside this timer
tool_exit_code=0
/usr/bin/time --quiet --format "$prefix\t%e\t%U\t%S\t%M\t%x" \
    --output "$BENCH_WORK/time.tsv" "${command[@]}" \
    > "$BENCH_WORK/tool.log" 2>&1 || tool_exit_code=$?
if ((tool_exit_code != 0)); then
    cat "$BENCH_WORK/tool.log" >&2
    exit "$tool_exit_code"
fi
if [[ ! -s "$track" ]]; then
    echo "The command succeeded but its expected track is missing: $track" >&2
    exit 1
fi
STEP

        # Append immediately so completed rows survive a later failure
        # If the step failed before producing a timing row, retain the attempted
        # run's known identifiers and mark the unavailable fields NA
        if [[ ! -s "$BENCH_WORK/time.tsv" ]]; then
            printf '%s\t%s\t%s\t%s\t%s\tNA\t%s\tNA\tNA\tNA\tNA\tNA\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n' \
                "$run_order" "$repetition" "$is_warmup" "$method" "$workload" "$cores" "$SLURM_JOB_ID" \
                > "$BENCH_WORK/time.tsv"
        fi
        printf '%s\t%s\n' "$(cat "$BENCH_WORK/time.tsv")" "$step_exit_code" >> "$BENCH_OUTPUT"
        if ((step_exit_code != 0)); then
            echo "Run $run_order failed. Completed rows are in $BENCH_OUTPUT" >&2
            exit "$step_exit_code"
        fi

        # Flush normal output writes, then remove the temporary tracks and logs
        # The next invocation gets fresh outputs, with no explicit input-cache reset
        sync
        rm -rf "$BENCH_WORK/coverage" "$BENCH_WORK/tool.log"
    done < "$BENCH_WORK/order.tsv"
done
