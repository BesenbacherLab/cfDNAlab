#!/usr/bin/env bash
#SBATCH --job-name=cfdna-end-motifs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --no-requeue

# Run native end-motif commands on the same exclusively allocated node
# See README.md for the settings, known differences and pilot invocation
# Results: timings.tsv, plus one native motif profile per method in profiles/
# All warmups remain in the timing TSV. No host caches are explicitly evicted
# Input BAM, BAI, reference and packed environments are copied to local scratch
# Source files are never modified. Only this script's scratch directory is deleted

set -euo pipefail
export LC_ALL=C
source "${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks/common.sh"
BENCH_FEATURE=ends

: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${BENCH_BAM:?Set the prepared input BAM path}"
: "${BENCH_LOCAL_ROOT:?Set an existing node-local scratch directory}"
: "${BENCH_RESULTS:?Set a new results directory}"
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM}.bai}
: "${BENCH_CORES:?Set a single core budget for this job: 1, 4, 12 or 32}"
BENCH_WARMUPS=${BENCH_WARMUPS:-2}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-5}
BENCH_METHODS=${BENCH_METHODS:-cfdnalab-reference cfdnalab-reference-proper finaletoolkit}
BENCH_FREIA_ARCHIVE=${BENCH_FREIA_ARCHIVE:-}
export BENCH_K=${BENCH_K:-4}
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-30}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-1000}

benchmark_controls
read -r -a cores_to_run <<< "$BENCH_CORES"
read -r -a methods_to_run <<< "$BENCH_METHODS"
if ((${#cores_to_run[@]} == 0 || ${#methods_to_run[@]} == 0)); then
    echo "Supply at least a run label in BENCH_METHODS and a core budget in BENCH_CORES." >&2
    exit 1
fi
if [[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 ]]; then
    echo "A single-node allocation is required." >&2
    exit 1
fi
if [[ -e "$BENCH_RESULTS" ]]; then
    echo "The results directory already exists: $BENCH_RESULTS" >&2
    exit 1
fi
if [[ "$BENCH_INDEX" != *.bai ]]; then
    echo "FinaleToolkit requires a BAI index." >&2
    exit 1
fi
for value in "$BENCH_WARMUPS" "$BENCH_REPETITIONS" "$BENCH_K" "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" "${cores_to_run[@]}"; do
    if [[ ! "$value" =~ ^(0|[1-9][0-9]*)$ ]]; then
        echo "Core counts, repetitions and motif/filter settings must be integers: $value" >&2
        exit 1
    fi
done
if ((BENCH_REPETITIONS < 1 || (BENCH_K != 3 && BENCH_K != 4) || BENCH_MIN_MAPQ > 255 || BENCH_MIN_LENGTH < 10 || BENCH_MAX_LENGTH < BENCH_MIN_LENGTH || BENCH_MAX_LENGTH > 50000)); then
    echo "Use k=3 or 4, MAPQ 0-255, lengths within 10-50000 bp, and at least one measured repetition." >&2
    exit 1
fi
for cores in "${cores_to_run[@]}"; do
    if ((cores < 1 || cores > ${SLURM_CPUS_PER_TASK:-0})); then
        echo "Core budget $cores is outside this job's CPUs per task." >&2
        exit 1
    fi
done
run_freia=false
for method in "${methods_to_run[@]}"; do
    case "$method" in
        cfdnalab-read|cfdnalab-read-aligned|cfdnalab-read-include|cfdnalab-reference|cfdnalab-reference-proper|cfdnalab-reference-proper-wide|cfdnalab-reference-autosomes|finaletoolkit) ;;
        freia)
            : "${BENCH_FREIA_ARCHIVE:?Set the packed FrEIA environment for the freia comparison}"
            if ((BENCH_K != 3)); then
                echo "FrEIA's native abundance command produces 3-mers. Set BENCH_K=3 for this series." >&2
                exit 1
            fi
            run_freia=true
            echo "FrEIA timing only: pinned both-end aggregation is incorrect. Concordance is disabled." >&2
            ;;
        *) echo "Unknown end-motif method: $method" >&2; exit 1 ;;
    esac
done

run_cfdna=false
run_finale=false
need_reference=false
for method in "${methods_to_run[@]}"; do
    case "$method" in
        cfdnalab*) run_cfdna=true ;;
        finaletoolkit) run_finale=true ;;
    esac
    if [[ "$method" == cfdnalab-reference* || "$method" == finaletoolkit ]]; then need_reference=true; fi
done
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the packed samtools environment for BAM header preparation}"
if [[ "$run_cfdna" == true ]]; then : "${BENCH_CFDNA:?Set the cfDNAlab executable}"; fi
if [[ "$run_finale" == true ]]; then : "${BENCH_FINALE_ARCHIVE:?Set the packed FinaleToolkit environment}"; fi
if [[ "$need_reference" == true ]]; then : "${BENCH_REF_2BIT:?Set the matching 2bit reference}"; fi

# Stage once before every tool's warmups, following the coverage benchmark policy
export BENCH_WORK BENCH_NODE BENCH_AUTOSOMES BENCH_CONTIGS BENCH_CAPTURE_CORES BENCH_RESULTS
export BENCH_CFDNA_VERSION BENCH_FINALE_VERSION BENCH_FREIA_VERSION
BENCH_WORK=$(mktemp -d "$BENCH_LOCAL_ROOT/end-motifs-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)
BENCH_CAPTURE_CORES=${cores_to_run[0]}
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.bai"
if [[ "$need_reference" == true ]]; then cp -L "$BENCH_REF_2BIT" "$BENCH_WORK/reference.2bit"; fi
if [[ "$run_cfdna" == true ]]; then benchmark_stage_binary "$BENCH_CFDNA" cfdna; fi
# sbatch uses a spool copy of this script. Submit from the repository root
export TMPDIR="$BENCH_WORK/tmp"
# Ignore ~/.local Python packages and inherited import paths so imports stay
# inside the staged environments. See the installation notes in ../README.md
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
export XDG_CACHE_HOME="$BENCH_WORK/cache"
export MPLCONFIGDIR="$BENCH_WORK/matplotlib"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# Relocate only selected competitors. Samtools supplies the shared BAM header
benchmark_unpack "$BENCH_SAMTOOLS_ARCHIVE" "$BENCH_WORK/samtools-env"
"$BENCH_WORK/samtools-env/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/chrom.sizes"
BENCH_AUTOSOMES=$(awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/ {print $1}' "$BENCH_WORK/chrom.sizes" | paste -sd, -)
BENCH_CONTIGS=$(cut -f 1 "$BENCH_WORK/chrom.sizes" | paste -sd, -)
if [[ "$run_finale" == true ]]; then benchmark_unpack "$BENCH_FINALE_ARCHIVE" "$BENCH_WORK/finale-env"; fi
if [[ "$run_freia" == true ]]; then
    benchmark_unpack "$BENCH_FREIA_ARCHIVE" "$BENCH_WORK/freia"
    BENCH_FREIA_VERSION=$(cat "$BENCH_WORK/freia/share/freia/version.txt")
fi

# FinaleToolkit visits the reference's contigs. cfDNAlab --chromosomes all uses
# the BAM's contigs. Check their names and lengths to give both the same scope
# The check reads headers only, with no filtering or changes to the BAM
if [[ "$run_finale" == true ]]; then
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
        raise SystemExit('FinaleToolkit 1.1.0 omits its last window for this contig length: ' + chromosome)
PYHEADER
fi
if [[ " $BENCH_METHODS " == *" cfdnalab-reference-autosomes "* && -z "$BENCH_AUTOSOMES" ]]; then echo "No autosomes in BAM." >&2; exit 1; fi
if [[ "$run_cfdna" == true ]]; then BENCH_CFDNA_VERSION=$("$BENCH_WORK/bin/cfdna" --version | tr '\t\n' '  '); fi
if [[ "$run_finale" == true ]]; then BENCH_FINALE_VERSION=$("$BENCH_WORK/finale-env/bin/finaletoolkit" --version | tr '\t\n' '  '); fi

mkdir -p "$BENCH_RESULTS/profiles"

# The requested length bounds and whether a method applies them are separate
# columns. This makes FinaleToolkit 1.1.0's limitation visible in every run
header=$'run_order\trepetition\twarmup\tmethod\ttool_version\tcores\ttool_workers\tk\tmin_mapq'
header+=$'\trequested_min_length\trequested_max_length\tlength_filter\tscope\tsource\tclip\tproper_pair'
header+=$'\tnode\tcpu_model\tcpu_ids\tphysical_core_ids\tjob_id\tstep_id\tstarted_utc'
header+=$'\teffective_min_length\teffective_max_length'
header+=$'\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kib\ttool_exit_code'
header+=$'\tflush_seconds\tprocess_plus_flush_seconds\toutput_bytes\tpilot\tstep_exit_code'
printf '%s\n' "$header" > "$BENCH_RESULTS/timings.tsv"

# Run each method's warmups and measurements together on the recorded node
benchmark_record_setup
benchmark_schedule > "$BENCH_WORK/order.tsv"
if [[ ! -s "$BENCH_WORK/order.tsv" ]]; then echo "The requested run matrix is empty." >&2; exit 1; fi
run_order=0
while IFS=$'\t' read -r method cores repetition is_warmup; do
    run_order=$((run_order + 1))
    step_exit_code=0
    rm -f "$BENCH_WORK/time.tsv" "$BENCH_WORK/io.tsv" "$BENCH_WORK/profile.path" "$BENCH_WORK/step.id"

    # Attach each sequential step to this existing job's node
    # Request binding and check the actual CPU mask below. Hard confinement
    # additionally depends on the cluster's Slurm cgroup configuration
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
scope=all
source=reference
clip=aligned
proper_pair=not-required
length_filter=applied
effective_min_length=$BENCH_MIN_LENGTH
effective_max_length=$BENCH_MAX_LENGTH
if [[ "$method" == cfdnalab-reference-proper-wide ]]; then effective_min_length=10; effective_max_length=50000; fi
tool_workers=$cores
if [[ "$method" == cfdnalab-* ]]; then
    tool_version=$BENCH_CFDNA_VERSION
    command=("$BENCH_WORK/bin/cfdna" ends
        --bam "$BENCH_WORK/sample.bam" --output-dir "$BENCH_WORK/output"
        --k-inside "$BENCH_K" --k-outside 0 --all-motifs
        --min-mapq "$BENCH_MIN_MAPQ" --min-fragment-length "$effective_min_length"
        --max-fragment-length "$effective_max_length" --n-threads "$cores")
    if [[ "$method" == cfdnalab-read* ]]; then
        source=read
        clip=skip
        if [[ "$method" == cfdnalab-read-aligned ]]; then clip=aligned; fi
        if [[ "$method" == cfdnalab-read-include ]]; then clip=include-at-aligned-boundary; fi
        command+=(--source-inside read --clip-strategy "$clip")
    else
        command+=(--source-inside reference --clip-strategy aligned
            --ref-2bit "$BENCH_WORK/reference.2bit")
    fi
    if [[ "$method" == cfdnalab-reference-autosomes ]]; then
        scope=autosomes
        command+=(--chromosomes "$BENCH_AUTOSOMES")
    else
        command+=(--chromosomes all)
    fi
    if [[ "$method" == cfdnalab-reference-proper* ]]; then
        proper_pair=required
        command+=(--require-proper-pair)
    fi
    profile="$BENCH_WORK/output/end_motifs.zarr"
elif [[ "$method" == freia ]]; then
    tool_version=$BENCH_FREIA_VERSION
    source=read
    clip=aligned
    profile="$BENCH_WORK/output/4_FrEIA/3_Abundances/benchmark/Genome_Lvl/T__sample.pq"
    # Follow the authors' BAM filtering, extraction and abundance commands
    # Count both ends and avoid bootstrap sampling or intentional downsampling
    # The upper length bound is exclusive in FrEIA, hence max_length + 1
    # FrEIA adds coordinate deduplication, and removes a whole pair if either
    # extracted 3-base sequence contains N. These native choices are retained
    # The samtools exclusion flags are unmapped, secondary, duplicate and
    # supplementary, respectively, as in the authors' preprocessing rule
    printf 'sample_name group\nsample benchmark\n' > "$BENCH_WORK/output/samples.tsv"
    command=(bash -c '
        set -euo pipefail
        cores=$1
        environment="$BENCH_WORK/freia"
        scripts="$environment/share/freia"
        "$environment/bin/samtools" view -@ "$cores" -q "$BENCH_MIN_MAPQ" \
            -F 4 -F 256 -F 1024 -F 2048 \
            -b -o "$BENCH_WORK/output/filtered.bam" "$BENCH_WORK/sample.bam"
        "$environment/bin/samtools" index -@ "$cores" "$BENCH_WORK/output/filtered.bam"
        "$environment/bin/python" "$scripts/1_extract_fragment_ends.py" \
            -b "$BENCH_WORK/output/filtered.bam" -o "$BENCH_WORK/output/sample.pq" \
            -c "$BENCH_CONTIGS" -t "$cores" -n 3
        "$environment/bin/python" "$scripts/2_data_transformation.py" \
            -i "$BENCH_WORK/output/sample.pq" -o "$BENCH_WORK/output/" \
            -st "$BENCH_WORK/output/samples.tsv" --both \
            -fsmin "$BENCH_MIN_LENGTH" -fsmax "$((BENCH_MAX_LENGTH + 1))" \
            -subs 1 --bootstrap_sample 0 --subchroff
    ' -- "$cores")
    # Extraction forces one worker above 8000 MiB of filtered BAM input
    # Its actual choice is only known after filtering. Record it below
else
    # FinaleToolkit reconstructs coordinates from read1/TLEN and requires 0x2
    # Its current end-motif reader ignores the requested length bounds
    tool_version=$BENCH_FINALE_VERSION
    length_filter=ignored-by-release
    effective_min_length=NA
    effective_max_length=NA
    proper_pair=required
    clip=unrestricted
    profile="$BENCH_WORK/output/motifs.tsv"
    command=("$BENCH_WORK/finale-env/bin/finaletoolkit" end-motifs
        "$BENCH_WORK/sample.bam" "$BENCH_WORK/reference.2bit" -k "$BENCH_K"
        --min-length "$BENCH_MIN_LENGTH" --max-length "$BENCH_MAX_LENGTH"
        --min-mapq "$BENCH_MIN_MAPQ" --threads "$cores" --strand both --output "$profile")
fi

# There are 25 identifying fields before GNU time's five measurement fields
printf -v prefix '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "$run_order" "$repetition" "$is_warmup" "$method" "$tool_version" "$cores" "$tool_workers" \
    "$BENCH_K" "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" "$length_filter" \
    "$scope" "$source" "$clip" "$proper_pair" "$node" "$cpu_model" "$cpu_ids" \
    "$physical_core_ids" "$SLURM_JOB_ID" "$SLURM_STEP_ID" "$started_utc"
printf -v prefix '%s\t%s\t%s' "$prefix" "$effective_min_length" "$effective_max_length"
logical_cpu_count=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
if [[ "$node" != "$BENCH_NODE" ]] || ((logical_cpu_count != cores || physical_core_count != cores)); then
    printf '%s\tNA\tNA\tNA\tNA\t1\n' "$prefix" > "$BENCH_WORK/time.tsv"
    echo "Unexpected placement: node=$node, CPUs=$cpu_ids, physical cores=$physical_core_ids" >&2
    exit 1
fi

# Time only the native CLI. Profile copying and comparisons are outside timing
tool_exit_code=0
benchmark_time_command "$prefix" "$BENCH_WORK/output" \
    "$BENCH_WORK/tool.log" "$BENCH_WORK/tool-error.log" "${command[@]}" || tool_exit_code=$?
if ((tool_exit_code != 0)); then
    cat "$BENCH_WORK/tool-error.log" "$BENCH_WORK/tool.log" >&2
    exit "$tool_exit_code"
fi
if [[ "$method" == freia ]]; then
    # Correct the worker field using the exact native size threshold
    # This records extraction workers, while samtools still uses the core budget
    if (( $(stat -c '%s' "$BENCH_WORK/output/filtered.bam") > 8000 * 1024 * 1024 )); then
        awk 'BEGIN {FS=OFS="\t"} {$7=1; print}' "$BENCH_WORK/time.tsv" > "$BENCH_WORK/updated-time.tsv"
        mv "$BENCH_WORK/updated-time.tsv" "$BENCH_WORK/time.tsv"
    fi
fi
if [[ ! -e "$profile" ]]; then
    echo "The expected motif profile is missing: $profile" >&2
    exit 1
fi

# Keep one native profile per method. All core settings still recompute outputs
printf '%s\n' "$profile" > "$BENCH_WORK/profile.path"
STEP

    # Retain known identifiers even when Slurm fails before GNU time can run
    if [[ ! -s "$BENCH_WORK/time.tsv" ]]; then
        fallback=()
        for ((column=1; column<=30; column++)); do fallback+=(NA); done
        fallback[0]=$run_order
        fallback[1]=$repetition
        fallback[2]=$is_warmup
        fallback[3]=$method
        fallback[5]=$cores
        fallback[20]=$SLURM_JOB_ID
        if [[ -s "$BENCH_WORK/step.id" ]]; then fallback[21]=$(cat "$BENCH_WORK/step.id"); fi
        (IFS=$'\t'; printf '%s\n' "${fallback[*]}") > "$BENCH_WORK/time.tsv"
    fi
    benchmark_finish_run "$method" "$cores" "$repetition" "$is_warmup" "$step_exit_code"
done < "$BENCH_WORK/order.tsv"
