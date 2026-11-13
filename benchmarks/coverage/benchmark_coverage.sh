#!/usr/bin/env bash
#SBATCH --job-name=cfdna-coverage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --no-requeue

# The submission launcher requests BENCH_CORES CPUs for this job's single core setting
# Full runs also request --exclusive and --mem=0 to reserve the node and its memory
# Pilots allow shared nodes and request 50 GB unless --mem is supplied
# Each srun binds this job's single selected core count
# User setup: the chosen partition must honor exclusive allocations for full runs
# Its policy takes precedence over the job's sharing request
# https://slurm.schedmd.com/sbatch.html#OPT_exclusive
#
# Results contain timings.tsv, saved profiles and a record of inputs/settings
# The first measured profile is kept for numerical comparison. Other outputs are temporary
# Diagnostics from a failed command go to the normal Slurm job log
# cfdnalab, mosdepth and samtools-depth process all BAM contigs for per-base depth
# cfdnalab-autosomes selects autosomes named chr1..chr22 or 1..22
# This preserves the default biological scope for either BAM naming convention
# pandepth also processes all BAM contigs, counting overlapping mates twice
# Its workload is labeled per-base-read-depth to retain that distinction
# cfdnalab-region-counts sums length-normalized fragment coverage in 500 kb windows
# finaletoolkit-coverage counts fragments overlapping the same window coordinates
# cfDNAlab splits a boundary fragment's contribution between windows. FinaleToolkit
# assigns it to the window containing its midpoint. Regional counts are separate from base depth
# cfdnalab-regions, mosdepth-regions, samtools-regions and pandepth-regions instead
# summarize read depth in the same windows, without writing per-base tracks
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
source "${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks/common.sh"
BENCH_FEATURE=coverage

# Input and tool paths are supplied in the environment when submitting the job
# BENCH_INDEX can override the usual sample.bam.bai path. FinaleToolkit requires
# a BAI index, so CSI indices are not supported by this combined benchmark
# Copied programs must have the libraries they need available on the compute node
# BENCH_FINALE_ARCHIVE contains a Linux Conda environment with FinaleToolkit
# installed, packed with conda-pack so Python and its packages can run locally
# BENCH_SAMTOOLS_ARCHIVE similarly includes samtools and its shared libraries
: "${SLURM_JOB_ID:?Submit with bash benchmarks/submit_benchmark.sh coverage}"
: "${BENCH_BAM:?Set the input BAM path}"
# TMPDIR is supplied on the compute node when the batch job starts
: "${TMPDIR:?The compute node must provide TMPDIR for node-local scratch}"
: "${BENCH_RESULTS:?Set a new results directory}"
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the path of the packed samtools environment}"
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM}.bai}
: "${BENCH_CORES:?Set a single core budget for this job: 1, 4, 12 or 32}"
BENCH_WARMUPS=${BENCH_WARMUPS:-1}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-3}
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
        cfdnalab|cfdnalab-autosomes|cfdnalab-fragment|cfdnalab-fragment-proper|mosdepth|mosdepth-fragment|samtools-depth|pandepth|cfdnalab-region-counts|cfdnalab-region-counts-proper|finaletoolkit-coverage|cfdnalab-regions|mosdepth-regions|samtools-regions|pandepth-regions) ;;
        *) echo "Unknown coverage method: $method" >&2; exit 1 ;;
    esac
done

# The BAM, index and programs are copied to the compute node before warmup
# Every srun below uses this job's node and these same local files
export BENCH_WORK BENCH_NODE BENCH_CFDNA_VERSION BENCH_MOSDEPTH_VERSION BENCH_AUTOSOMES
export BENCH_SAMTOOLS_VERSION BENCH_FINALE_VERSION
export BENCH_PANDEPTH_VERSION
BENCH_WORK=$(mktemp -d "$TMPDIR/coverage-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)

# Create result files before staging so setup failures leave an inspectable record
mkdir -p "$BENCH_RESULTS/profiles"
# run_order counts every attempted run in this job, while repetition counts runs
# of the current method, including warmups. CPU records identify jobs run on
# different hardware when their timings are interpreted
# tool_threads records the native option, set to the allocated core count
# It means cfDNAlab threads, FinaleToolkit workers, PanDepth's -t setting, or
# additional BAM decoding threads for mosdepth/samtools. Each tool's threads
# share the CPUs allocated by Slurm, including its main calculation thread
# Rows contain elapsed time, CPU time, RAM reported by GNU time and the tool's exit code
# step_exit_code also records failure to get the expected CPUs or output file
header=$'run_order\trepetition\twarmup\tmethod\tworkload\ttool_version\tcores\ttool_threads'
header+=$'\tnode\tcpu_model\tcpu_ids\tphysical_core_ids\tjob_id\tstep_id\tstarted_utc'
header+=$'\tmin_mapq\tmin_length\tmax_length\tscope\tproper_pair\tlength_filter'
header+=$'\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kib\ttool_exit_code'
header+=$'\tflush_seconds\tprocess_plus_flush_seconds\toutput_bytes\tpilot\tstep_exit_code'
printf '%s\n' "$header" > "$BENCH_RESULTS/timings.tsv"
benchmark_progress "Initializing $BENCH_FEATURE benchmark on $BENCH_NODE"
benchmark_record_setup
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp"
export TMPDIR="$BENCH_WORK/tmp"
# User setup: submit from the repository root so SLURM_SUBMIT_DIR locates helpers
# Slurm keeps its own copy of this job script elsewhere, away from those helpers

# Only selected tools are copied. Samtools also reads BAM chromosome names/lengths
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
export XDG_CACHE_HOME="$BENCH_WORK/cache" MPLCONFIGDIR="$BENCH_WORK/matplotlib"
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
        pandepth|pandepth-regions)
            : "${BENCH_PANDEPTH:?Set the PanDepth executable}"
            benchmark_stage_binary "$BENCH_PANDEPTH" pandepth ;;
        finaletoolkit-coverage)
            : "${BENCH_FINALE_ARCHIVE:?Set the packed FinaleToolkit environment}"
            benchmark_unpack "$BENCH_FINALE_ARCHIVE" "$BENCH_WORK/finale-env" ;;
    esac
done
benchmark_progress "Checking tool versions"
# Version records are collected only for the selected tools
# FinaleToolkit uses its packed Python so its launcher does not search the node's PATH
if [[ -f "$BENCH_WORK/bin/cfdna" ]]; then BENCH_CFDNA_VERSION=$("$BENCH_WORK/bin/cfdna" --version | tr '\t\n' '  '); fi
if [[ -f "$BENCH_WORK/bin/mosdepth" ]]; then BENCH_MOSDEPTH_VERSION=$("$BENCH_WORK/bin/mosdepth" --version | tr '\t\n' '  '); fi
BENCH_SAMTOOLS_VERSION=$("$BENCH_WORK/samtools-env/bin/samtools" --version-only | tr '\t\n' '  ')
if [[ -d "$BENCH_WORK/finale-env" ]]; then
    BENCH_FINALE_VERSION=$("$BENCH_WORK/finale-env/bin/python" "$BENCH_WORK/finale-env/bin/finaletoolkit" --version | tr '\t\n' '  ')
    BENCH_FINALE_VERSION+=" commit:$(cat "$BENCH_WORK/finale-env/share/finaletoolkit/revision.txt")"
fi
if [[ -f "$BENCH_WORK/bin/pandepth" ]]; then
    BENCH_PANDEPTH_VERSION=$("$BENCH_WORK/bin/pandepth" -h 2>&1 | sed -n 's/.*show this help \[\(v[^]]*\)\].*/\1/p')
    if [[ -z "$BENCH_PANDEPTH_VERSION" ]]; then echo "Could not read PanDepth version." >&2; exit 1; fi
fi
benchmark_progress "Copying BAM and index to node-local scratch"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.${BENCH_INDEX##*.}"
"$BENCH_WORK/samtools-env/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/chrom.sizes"

# The script generates a shared BED4 of 500 kb windows, including shorter final windows
# idxstats supplies contig lengths from the BAM index. Its unmapped '*' row is excluded
# This setup is outside every timer. FinaleToolkit and samtools bedcov read this BED
# Native fixed-size windows in cfDNAlab/mosdepth use the same coordinates
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
# PanDepth generates its own 500 kb windows with -w, like cfDNAlab and mosdepth
# Its native window table is converted to BED coordinates only during comparison

# Autosome selection uses the BAM's actual names. Snyder's reference can use 1..22,
# while cfDNAlab defaults to chr1..chr22. The file's chromosome names remain unchanged
BENCH_AUTOSOMES=$(cut -f 1 "$BENCH_WORK/regions.bed" | uniq |
    awk '/^(chr)?([1-9]|1[0-9]|2[0-2])$/' | paste -sd, -)
if [[ " $BENCH_METHODS " == *" cfdnalab-autosomes "* && -z "$BENCH_AUTOSOMES" ]]; then
    echo "The BAM has no autosomes named chr1..chr22 or 1..22." >&2
    exit 1
fi

# Every processed chromosome is retained for concordance, including positional
# samtools/PanDepth outputs. Saving those profiles happens after the tool timer
cp "$BENCH_WORK/regions.bed" "$BENCH_RESULTS/regions.bed"
export BENCH_RESULTS BENCH_MIN_MAPQ BENCH_MIN_LENGTH BENCH_MAX_LENGTH

# The BAM and programs are ready on the compute node. Preparation time above
# is excluded from tool timings. The script now lists each method's warmup and measurements
# srun starts a separate program for each run. After it finishes, this main script
# saves that run's results. common.sh explains the files used to pass them back
# Staged files finish writing before warmups. This does not evict the page cache
sync -f "$BENCH_WORK"
benchmark_progress "Setup complete. Starting warmups and measurements"
benchmark_schedule > "$BENCH_WORK/order.tsv"
if [[ ! -s "$BENCH_WORK/order.tsv" ]]; then echo "The requested run matrix is empty." >&2; exit 1; fi
run_order=0
job_exit_code=0
failed_method=
while IFS=$'\t' read -r method cores repetition is_warmup; do
    if [[ "$method" == "$failed_method" ]]; then continue; fi
    run_order=$((run_order + 1))
    benchmark_progress "Starting run $run_order: $method, repetition $repetition (warmup=$is_warmup)"
    step_exit_code=0
    # Old timings and output paths are removed before each run. If this run fails
    # before writing its own files, it cannot accidentally reuse the preceding results
    rm -f "$BENCH_WORK/time.tsv" "$BENCH_WORK/io.tsv" "$BENCH_WORK/profile.path" "$BENCH_WORK/step.id" "$BENCH_WORK/stdout.path"
    workload=per-base-read-depth
    if [[ "$method" == *fragment* ]]; then workload=per-base-fragment-coverage; fi
    if [[ "$method" == cfdnalab-region-counts* || "$method" == finaletoolkit-coverage ]]; then
        workload=region-counts-500kb
    elif [[ "$method" == *-regions ]]; then
        workload=regional-read-depth-500kb
    elif [[ "$method" == pandepth ]]; then
        workload=per-base-read-depth
    fi

    # srun starts the tool on the reserved node with exactly the requested cores
    # --cpu-bind keeps the program on those CPUs. Enforcement also depends on the
    # cluster's settings, so the program below checks the CPUs it actually received
    # The program between <<'STEP' and the closing STEP line performs this run
    # Quoting STEP delays expansion of its variables until that program starts
    # --export=ALL passes the input/tool paths, and the six arguments identify the run
    # It writes results under BENCH_WORK for this main script to read afterward
    # Its exit code is saved so failure can be recorded before another method is tried
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

# Node and CPU assignment are read from the program that will run the tool
node=$(hostname)
cpu_model=$(awk -F ': ' '/^model name/ {print $2; exit}' /proc/cpuinfo)
cpu_ids=$(awk '/Cpus_allowed_list/ {print $2}' /proc/self/status)
started_utc=$(date -u +%FT%TZ)

# Hyperthreading can give a physical core two CPU IDs. The script looks up
# the actual core behind each allowed ID to verify the requested core count
# Core numbers can repeat across CPU chips, so chip (socket) and core number
# together identify a core. Lists such as 0-3,8 are expanded to inspect each ID
# A single ID such as 8 is handled as the range 8-8
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
if [[ "$method" == cfdnalab || "$method" == cfdnalab-autosomes || "$method" == cfdnalab-fragment* || "$method" == cfdnalab-regions ]]; then
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
    if [[ "$method" == cfdnalab-regions ]]; then
        # Regional read depth is the sum of covered bases, without length normalization
        command+=(--by-size 500000 --per-window total)
        track="$BENCH_WORK/coverage/fcoverage.total.tsv.zst"
    fi
elif [[ "$method" == mosdepth || "$method" == mosdepth-fragment || "$method" == mosdepth-regions ]]; then
    # --threads controls BAM decoding, while coverage calculation runs separately
    # Pass the allocated core count directly, as for every native thread option
    tool_threads=$cores
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
    if [[ "$method" == mosdepth-regions ]]; then
        # The documented precision setting limits rounding when means are later
        # multiplied by window widths. It changes output precision, not counting
        export MOSDEPTH_PRECISION=5
        command+=(--by 500000 --no-per-base)
    fi
    command+=("$BENCH_WORK/coverage/coverage" "$BENCH_WORK/sample.bam")
    track="$BENCH_WORK/coverage/coverage.per-base.bed.gz"
    if [[ "$method" == mosdepth-regions ]]; then track="$BENCH_WORK/coverage/coverage.regions.bed.gz"; fi
elif [[ "$method" == samtools-depth ]]; then
    # -@ receives the allocated core count and controls additional decoding threads
    # -s removes mate overlap. Native sparse output omits zero-depth positions
    # -q is base quality and -Q is MAPQ. No -J means deletions are excluded
    # Basic flags match without requiring proper pairs. The HTSlib input expression
    # filters both TLEN signs with inclusive bounds during the same BAM read pass
    # This does not require both mates to pass filtering or reconstruct their span
    # Its native output is uncompressed text with a 1-based position per row
    # https://www.htslib.org/doc/samtools-depth.html
    tool_threads=$cores
    tool_version=$BENCH_SAMTOOLS_VERSION
    track="$BENCH_WORK/coverage/samtools.depth.tsv"
    command=("$BENCH_WORK/samtools-env/bin/samtools" depth -@ "$tool_threads"
        --input-fmt-option "filter=(tlen >= $BENCH_MIN_LENGTH && tlen <= $BENCH_MAX_LENGTH) || (tlen <= -$BENCH_MIN_LENGTH && tlen >= -$BENCH_MAX_LENGTH)"
        -s -q 0 -Q "$BENCH_MIN_MAPQ" -G 3852 --require-flags 1
        -o "$track" "$BENCH_WORK/sample.bam")
elif [[ "$method" == samtools-regions ]]; then
    # bedcov sums read depth directly over the BED. -j excludes D/N positions
    # Unlike depth -s, bedcov counts overlapping mates twice and has no thread
    # option. samtools' global input filter supplies paired selection and TLEN bounds
    tool_threads=NA
    tool_version=$BENCH_SAMTOOLS_VERSION
    track="$BENCH_WORK/coverage/samtools.bedcov.tsv"
    command=("$BENCH_WORK/samtools-env/bin/samtools" bedcov
        --input-fmt-option "filter=flag.paired && ((tlen >= $BENCH_MIN_LENGTH && tlen <= $BENCH_MAX_LENGTH) || (tlen <= -$BENCH_MIN_LENGTH && tlen >= -$BENCH_MAX_LENGTH))"
        -j -Q "$BENCH_MIN_MAPQ" -G 3852
        "$BENCH_WORK/regions.bed" "$BENCH_WORK/sample.bam")
elif [[ "$method" == pandepth || "$method" == pandepth-regions ]]; then
    # -a writes every base, including zeros, to a compressed SiteDepth table
    # MAPQ and exclusion flags match. PanDepth has no fragment length filter,
    # paired-read inclusion flag, complete-pair check or mate-overlap correction
    # Overlapping mates therefore count twice, unlike the paired coverage runs
    # Its native position column is 0-based. It also writes chromosome statistics
    # https://github.com/HuiyangYu/PanDepth/blob/main/src/PanDepth.cpp
    tool_threads=$cores
    tool_version=$BENCH_PANDEPTH_VERSION
    length_filter=unsupported
    command=("$BENCH_WORK/bin/pandepth"
        -i "$BENCH_WORK/sample.bam" -o "$BENCH_WORK/coverage/pandepth"
        -q "$BENCH_MIN_MAPQ" -x 3852 -t "$tool_threads")
    track="$BENCH_WORK/coverage/pandepth.SiteDepth.gz"
    if [[ "$method" == pandepth-regions ]]; then
        # Native regional statistics include TotalDepth, avoiding a per-base table
        command+=(-w 500000)
        track="$BENCH_WORK/coverage/pandepth.win.stat.gz"
    else
        command+=(-a)
    fi
elif [[ "$method" == cfdnalab-region-counts* ]]; then
    # Length normalization gives each fragment a total weight of 1. Summing this
    # coverage estimates fragment counts in each fixed-size window. A fragment
    # crossing a boundary splits its contribution, unlike FinaleToolkit's full
    # count in its midpoint's window. Native pairing and coverage rules remain
    tool_threads=$cores
    tool_version=$BENCH_CFDNA_VERSION
    command=("$BENCH_WORK/bin/cfdna" fcoverage
        --bam "$BENCH_WORK/sample.bam" --output-dir "$BENCH_WORK/coverage"
        --chromosomes all --by-size 500000 --per-window total --normalize-by-length
        --min-fragment-length "$BENCH_MIN_LENGTH" --max-fragment-length "$BENCH_MAX_LENGTH"
        --min-mapq "$BENCH_MIN_MAPQ" --n-threads "$tool_threads")
    if [[ "$method" == cfdnalab-region-counts-proper ]]; then command+=(--require-proper-pair); proper_pair=required; fi
    track="$BENCH_WORK/coverage/length_normalized.fcoverage.total.tsv.zst"
else
    # This command calculates raw fragment counts with the native midpoint rule
    # A boundary fragment belongs to the window containing its midpoint
    # No normalization is requested, and the native scale factor is 1
    # --threads controls worker processes in the current FinaleToolkit CLI
    # Its BAM reader requires proper pairs and uses read1/TLEN. Its indexed query
    # also needs read1 to overlap the window, which can miss boundary fragments
    # The comparator reports r/MAE/RMSE to assess agreement despite these differences
    # https://epifluidlab.github.io/FinaleToolkit/documentation/cli_reference/index.html#coverage
    tool_threads=$cores
    tool_version=$BENCH_FINALE_VERSION
    proper_pair=required
    track="$BENCH_WORK/coverage/finaletoolkit.coverage.bed.gz"
    command=("$BENCH_WORK/finale-env/bin/python" "$BENCH_WORK/finale-env/bin/finaletoolkit" coverage
        "$BENCH_WORK/sample.bam" "$BENCH_WORK/regions.bed"
        --output "$track"
        --min-length "$BENCH_MIN_LENGTH" --max-length "$BENCH_MAX_LENGTH" --min-mapq "$BENCH_MIN_MAPQ" --threads "$tool_threads")
fi

# The timing row starts with method, settings and run IDs, followed by measurements
printf -v prefix '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "$run_order" "$repetition" "$is_warmup" "$method" "$workload" "$tool_version" "$cores" "$tool_threads" \
    "$node" "$cpu_model" "$cpu_ids" "$physical_core_ids" "$SLURM_JOB_ID" "$SLURM_STEP_ID" "$started_utc" \
    "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" "$scope" "$proper_pair" "$length_filter"

# This check requires the intended node and CPU counts. OpenMP settings are ignored
# when calling nproc so it counts CPUs actually available to the program
# A wrong assignment produces a failed row with NA measurements
logical_cpu_count=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
if [[ "$node" != "$BENCH_NODE" ]] || ((logical_cpu_count != cores || physical_core_count != cores)); then
    printf '%s\tNA\tNA\tNA\tNA\t1\n' "$prefix" > "$BENCH_WORK/time.tsv"
    echo "Unexpected placement: node=$node, CPUs=$cpu_ids, physical cores=$physical_core_ids" >&2
    exit 1
fi

# Timing covers the tool command, including writing and compressing its output
# Waiting for Slurm, checking CPUs and later saving/cleaning output are excluded
tool_exit_code=0
# bedcov prints its native table to stdout rather than accepting an output path
tool_stdout="$BENCH_WORK/tool.log"
if [[ "$method" == samtools-regions ]]; then tool_stdout=$track; fi
benchmark_time_command "$prefix" "$BENCH_WORK/coverage" \
    "$tool_stdout" "$BENCH_WORK/tool-error.log" "${command[@]}" || tool_exit_code=$?
if ((tool_exit_code != 0)); then
    cat "$BENCH_WORK/tool-error.log" "$tool_stdout" >&2
    exit "$tool_exit_code"
fi
# With sparse samtools depth output, a real all-zero result is an empty file
# Other native formats retain headers/compression frames even with no coverage
if [[ ! -f "$track" || ( ! -s "$track" && "$method" != samtools-depth ) ]]; then
    echo "The command succeeded but its expected track is missing: $track" >&2
    exit 1
fi
printf '%s\n' "$track" > "$BENCH_WORK/profile.path"
STEP

    # A run that failed before timing still gets a row naming the attempted
    # method and repetition. Unknown measurements are recorded as NA
    # Column numbers below start at 0 and must match the timing-table header above
    # benchmark_finish_run adds disk-write measurements and the srun exit code
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
    # benchmark_finish_run saves the timing row and useful output before cleanup
    # It is called on its own line so Bash stops the job if saving fails
    # A tool failure is recorded first and then handled by the following if block
    benchmark_finish_run "$method" "$repetition" "$is_warmup" "$step_exit_code"
    if ((step_exit_code != 0)); then
        # All repetitions of this method are listed together. After failure, its
        # remaining runs are skipped until the next method name is reached
        failed_method=$method
        job_exit_code=1
        echo "Skipping remaining repetitions of $method after failure; continuing other selected methods." >&2
    fi
done < "$BENCH_WORK/order.tsv"

# Concordance runs after all measurements, in both pilots and full benchmark jobs
# All-contig profiles retain the full scope. Optional autosome diagnostics are
# compared on the chromosomes shared with their whole-reference counterparts
benchmark_compare_results "$SLURM_SUBMIT_DIR/benchmarks/coverage/compare_coverage.py" --common-scope

# The job reports failure if any method failed, while keeping completed results
exit "$job_exit_code"
