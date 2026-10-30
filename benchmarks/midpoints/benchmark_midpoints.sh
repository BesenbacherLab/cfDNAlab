#!/usr/bin/env bash
#SBATCH --job-name=cfdna-midpoints
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --no-requeue

# Compare native, unsmoothed 1 bp midpoint profiles on one exclusive node
# Every attempted run, including warmups, gets a row in timings.tsv
# Keep a native profile from the first measured repetition for each tool
# All BAM/reference copies, environments and temporary files use local scratch
# See README.md for required paths, concordance and expected differences

set -euo pipefail
export LC_ALL=C
source "${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks/common.sh"
BENCH_FEATURE=midpoints

: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${BENCH_BAM:?Set the prepared coordinate-sorted BAM path}"
: "${BENCH_INTERVALS:?Set the grouped fixed-size BED path}"
: "${BENCH_LOCAL_ROOT:?Set an existing node-local scratch directory}"
: "${BENCH_RESULTS:?Set a new results directory}"
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM}.bai}
BENCH_FASTA_INDEX=${BENCH_FASTA_INDEX:-${BENCH_REF_FASTA:-}.fai}
: "${BENCH_CORES:?Set a single core budget for this job: 1, 4, 12 or 32}"
BENCH_WARMUPS=${BENCH_WARMUPS:-2}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-5}
BENCH_METHODS=${BENCH_METHODS:-cfdnalab griffin}
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-100}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-200}

benchmark_controls
read -r -a cores_to_run <<< "$BENCH_CORES"
read -r -a methods_to_run <<< "$BENCH_METHODS"
if [[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 || -e "$BENCH_RESULTS" ]]; then
    echo "Use a single-node allocation and a new BENCH_RESULTS directory." >&2
    exit 1
fi
if ((${#cores_to_run[@]} == 0 || ${#methods_to_run[@]} == 0)); then
    echo "Supply at least a tool in BENCH_METHODS and a core budget in BENCH_CORES." >&2
    exit 1
fi
for value in "$BENCH_WARMUPS" "$BENCH_REPETITIONS" "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" "${cores_to_run[@]}"; do
    if [[ ! "$value" =~ ^(0|[1-9][0-9]*)$ ]]; then
        echo "Core counts, repetitions and filter settings must be integers: $value" >&2
        exit 1
    fi
done
if ((BENCH_REPETITIONS < 1 || BENCH_MIN_MAPQ > 255 || BENCH_MIN_LENGTH < 10 || BENCH_MAX_LENGTH < BENCH_MIN_LENGTH || BENCH_MAX_LENGTH > 50000)); then
    echo "Use MAPQ 0-255, lengths within 10-50000 bp, and at least one measured repetition." >&2
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
        cfdnalab|griffin) ;;
        *) echo "Unknown midpoint tool: $method" >&2; exit 1 ;;
    esac
done
if [[ " $BENCH_METHODS " == *" griffin "* && ${BENCH_REF_FASTA:-} == *.gz ]]; then
    echo "Use an uncompressed FASTA and its existing .fai index." >&2
    exit 1
fi

run_cfdna=false
run_griffin=false
for method in "${methods_to_run[@]}"; do
    case "$method" in
        cfdnalab*) run_cfdna=true ;;
        griffin) run_griffin=true ;;
    esac
done
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the packed samtools environment for BAM header preparation}"
if [[ "$run_cfdna" == true ]]; then : "${BENCH_CFDNA:?Set the cfDNAlab executable}"; fi
if [[ "$run_griffin" == true ]]; then : "${BENCH_GRIFFIN_ARCHIVE:?Set the packed Griffin environment}"; : "${BENCH_REF_FASTA:?Set the matching reference FASTA}"; fi

# Stage once before warmups. Never rewrite or index the original inputs
export BENCH_WORK BENCH_NODE BENCH_HALF_WIDTH BENCH_SITES BENCH_GROUPS BENCH_CHROMOSOMES
export BENCH_CAPTURE_CORES BENCH_RESULTS BENCH_CFDNA_VERSION BENCH_GRIFFIN_VERSION
BENCH_WORK=$(mktemp -d "$BENCH_LOCAL_ROOT/midpoints-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)
BENCH_CAPTURE_CORES=${cores_to_run[0]}
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib" "$BENCH_WORK/griffin-env"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.${BENCH_INDEX##*.}"
cp -L "$BENCH_INTERVALS" "$BENCH_WORK/input.bed"
if [[ "$run_griffin" == true ]]; then
    cp -L "$BENCH_REF_FASTA" "$BENCH_WORK/reference.fa"
    cp -L "$BENCH_FASTA_INDEX" "$BENCH_WORK/reference.fa.fai"
fi
if [[ "$run_cfdna" == true ]]; then benchmark_stage_binary "$BENCH_CFDNA" cfdna; fi
# sbatch uses a spool copy of this script. Submit from the repository root
export TMPDIR="$BENCH_WORK/tmp" XDG_CACHE_HOME="$BENCH_WORK/cache"
# Ignore ~/.local Python packages and inherited import paths so imports stay
# inside Griffin's staged environment. See the installation notes in ../README.md
export MPLCONFIGDIR="$BENCH_WORK/matplotlib" PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
benchmark_unpack "$BENCH_SAMTOOLS_ARCHIVE" "$BENCH_WORK/samtools-env"
"$BENCH_WORK/samtools-env/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/chrom.sizes"
if [[ "$run_griffin" == true ]]; then
    benchmark_unpack "$BENCH_GRIFFIN_ARCHIVE" "$BENCH_WORK/griffin-env"
    # pybedtools must find the packed bedtools, rather than a shared installation
    export PATH="$BENCH_WORK/griffin-env/bin:$PATH"
fi
export BENCH_RUN_GRIFFIN=$run_griffin

# Prepare two representations of exactly the same oriented windows
# Keep all sites, including repeated and overlapping rows. Fail on invalid input
# This setup reads headers and BED rows, without scanning or filtering BAM reads
site_summary=$("$BENCH_SETUP_PYTHON" - <<'PY'
import csv
import json
import os
from pathlib import Path

work = Path(os.environ['BENCH_WORK'])
with (work / 'chrom.sizes').open() as source:
    bam_contigs = {fields[0]: int(fields[1]) for fields in csv.reader(source, delimiter='\t')}
reference_contigs = bam_contigs
if os.environ['BENCH_RUN_GRIFFIN'] == 'true':
    with (work / 'reference.fa.fai').open() as source:
        reference_contigs = {fields[0]: int(fields[1]) for fields in csv.reader(source, delimiter='\t')}
rows = []
window_width = None
with (work / 'input.bed').open() as source:
    for line_number, fields in enumerate(csv.reader(source, delimiter='\t'), start=1):
        if not fields or fields[0].startswith('#'):
            continue
        if len(fields) < 4:
            raise ValueError('BED row {} needs chromosome, start, end and group'.format(line_number))
        chromosome, start, end, group = fields[:4]
        start, end = int(start), int(end)
        if not group or any(character.isspace() for character in group) or start < 0 or end <= start:
            raise ValueError('Invalid interval or group name on BED row {}'.format(line_number))
        if chromosome not in bam_contigs or reference_contigs.get(chromosome) != bam_contigs[chromosome]:
            raise ValueError('BAM/FASTA contig missing or length mismatch: ' + chromosome)
        if end > bam_contigs[chromosome]:
            raise ValueError('BED row {} extends beyond its chromosome'.format(line_number))
        if window_width is None:
            window_width = end - start
        if end - start != window_width or window_width % 2:
            raise ValueError('Use equal, even-width BED intervals for symmetric Griffin windows')
        # Follow cfDNAlab's BED5/BED6 strand convention, including ambiguity checks
        strand = fields[5] if len(fields) >= 6 else fields[4] if len(fields) == 5 else '.'
        if len(fields) >= 6 and fields[4] in ('+', '-', '.') and strand not in ('+', '-', '.'):
            raise ValueError('Ambiguous strand columns on BED row {}'.format(line_number))
        if strand not in ('+', '-', '.'):
            raise ValueError('Invalid strand on BED row {}'.format(line_number))
        rows.append((chromosome, start, end, group, 0, strand))
if not rows:
    raise ValueError('No midpoint sites in the BED')
rows.sort(key=lambda row: (row[0], row[1], row[2]))
with (work / 'sites.bed').open('w') as destination:
    csv.writer(destination, delimiter='\t', lineterminator='\n').writerows(rows)
# Griffin takes a table of anchor positions per group, referenced by YAML
# With symmetric even-width windows these anchors recover the same BED bounds
groups = {}
for chromosome, start, end, group, score, strand in rows:
    groups.setdefault(group, []).append((chromosome, start + window_width // 2, strand))
(work / 'sites').mkdir()
site_paths = {}
for group_index, (group, sites) in enumerate(groups.items()):
    if os.environ['BENCH_RUN_GRIFFIN'] == 'true' and all(site[0].isdigit() for site in sites):
        raise ValueError('Griffin infers numeric chromosome types for group ' + group + '. Use the selected hg38/UCSC input.')
    path = work / 'sites' / '{}.tsv'.format(group_index)
    with path.open('w') as destination:
        writer = csv.writer(destination, delimiter='\t', lineterminator='\n')
        writer.writerow(('Chrom', 'position', 'Strand'))
        writer.writerows(sites)
    site_paths[group] = str(path)
with (work / 'sites.yaml').open('w') as destination:
    # JSON is valid YAML and avoids requiring PyYAML in the setup interpreter
    json.dump({'site_lists': site_paths}, destination)
chromosomes = sorted(set(row[0] for row in rows))
# Its counting script sorts intervals lexicographically. Use that same contig
# order in the bigWig header so its exported entries remain correctly ordered
with (work / 'chrom.sizes').open('w') as destination:
    for chromosome in chromosomes:
        destination.write('{}\t{}\n'.format(chromosome, bam_contigs[chromosome]))
# Griffin always calculates GC-weighted output alongside raw counts. Unit bias
# supplies the mandatory table without performing a sample-specific correction
with (work / 'unit_gc_bias.tsv').open('w') as destination:
    writer = csv.writer(destination, delimiter='\t', lineterminator='\n')
    writer.writerow(('length', 'num_GC', 'smoothed_GC_bias'))
    for fragment_length in range(int(os.environ['BENCH_MIN_LENGTH']), int(os.environ['BENCH_MAX_LENGTH']) + 1):
        for gc_count in range(fragment_length + 1):
            writer.writerow((fragment_length, gc_count, 1.0))
print(window_width // 2, len(rows), len(groups), ','.join(chromosomes), sep='\t')
PY
)
IFS=$'\t' read -r BENCH_HALF_WIDTH BENCH_SITES BENCH_GROUPS BENCH_CHROMOSOMES <<< "$site_summary"
if [[ "$run_cfdna" == true ]]; then BENCH_CFDNA_VERSION=$("$BENCH_WORK/bin/cfdna" --version 2>&1 | tr '\t\n' '  '); fi
if [[ "$run_griffin" == true ]]; then BENCH_GRIFFIN_VERSION=$(cat "$BENCH_WORK/griffin-env/share/griffin/version.txt"); fi
mkdir -p "$BENCH_RESULTS/profiles"
cp "$BENCH_WORK/sites.bed" "$BENCH_RESULTS/sites.bed"
sha256sum "$BENCH_INTERVALS" > "$BENCH_RESULTS/sites.sha256"

# One row per tool run, with the actual node and allowed CPU IDs
header=$'run_order\trepetition\twarmup\tmethod\ttool_version\tcores\ttool_workers'
header+=$'\tmin_mapq\tmin_length\tmax_length\tsites\tgroups\twindow_bp\tresolution_bp\tmidpoint_policy'
header+=$'\tnode\tcpu_model\tcpu_ids\tphysical_core_ids\tjob_id\tstep_id\tstarted_utc'
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
    # Attach to this job's existing node, then check the binding in the step
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
mkdir -p "$BENCH_WORK/output/tmp" "$BENCH_WORK/output/results"
if [[ "$method" == cfdnalab ]]; then
    tool_version=$BENCH_CFDNA_VERSION
    midpoint_policy=coordinate-random-even
    command=("$BENCH_WORK/bin/cfdna" midpoints
        --bam "$BENCH_WORK/sample.bam" --output-dir "$BENCH_WORK/output"
        --intervals "$BENCH_WORK/sites.bed" --chromosomes "$BENCH_CHROMOSOMES"
        --length-bins "$BENCH_MIN_LENGTH" "$((BENCH_MAX_LENGTH + 1))"
        --min-mapq "$BENCH_MIN_MAPQ" --smoothing none --bin-size 1 --n-threads "$cores")
    profile="$BENCH_WORK/output/midpoint_profiles.zarr"
else
    tool_version=$BENCH_GRIFFIN_VERSION
    midpoint_policy=right-central-even
    # Time both native scripts, including bigWig writing and group aggregation
    # No masking, outlier removal, site-wise CNA normalization or smoothing
    # Griffin's final group-wise mean-one normalization cannot be disabled
    command=(bash -c '
        set -euo pipefail
        cores=$1
        environment="$BENCH_WORK/griffin-env"
        scripts="$environment/share/griffin"
        IFS=, read -r -a chromosomes <<< "$BENCH_CHROMOSOMES"
        common=(--sample_name sample --tmp_dir "$BENCH_WORK/output/tmp"
            --sites_yaml "$BENCH_WORK/sites.yaml" --griffin_scripts_dir "$scripts"
            --chrom_sizes_path "$BENCH_WORK/chrom.sizes" --chroms "${chromosomes[@]}"
            --chrom_column Chrom --position_column position --strand_column Strand
            --norm_window "-$BENCH_HALF_WIDTH" "$BENCH_HALF_WIDTH"
            --mappability_correction False --mappability_bw none --CPU "$cores")
        "$environment/bin/python" "$scripts/griffin_coverage.py" "${common[@]}" \
            --bam "$BENCH_WORK/sample.bam" --GC_bias "$BENCH_WORK/unit_gc_bias.tsv" \
            --reference_genome "$BENCH_WORK/reference.fa" --mappability_bias none \
            --size_range "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" --map_quality "$BENCH_MIN_MAPQ"
        "$environment/bin/python" "$scripts/griffin_merge_sites.py" "${common[@]}" \
            --uncorrected_bw_path "$BENCH_WORK/output/tmp/sample/tmp_bigWig/sample.uncorrected.bw" \
            --GC_corrected_bw_path "$BENCH_WORK/output/tmp/sample/tmp_bigWig/sample.GC_corrected.bw" \
            --GC_map_corrected_bw_path none --results_dir "$BENCH_WORK/output/results" \
            --save_window "-$BENCH_HALF_WIDTH" "$BENCH_HALF_WIDTH" \
            --center_window "-$BENCH_HALF_WIDTH" "$BENCH_HALF_WIDTH" \
            --fft_window "-$BENCH_HALF_WIDTH" "$BENCH_HALF_WIDTH" --fft_index 1 \
            --step 1 --exclude_paths none --exclude_zero_mappability False \
            --exclude_outliers False --CNA_normalization False --individual False --smoothing False
    ' -- "$cores")
    profile="$BENCH_WORK/output/results/sample/sample.uncorrected.coverage.tsv"
fi

# The first 22 fields identify the run. GNU time supplies five measurements
printf -v prefix '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "$run_order" "$repetition" "$is_warmup" "$method" "$tool_version" "$cores" "$cores" \
    "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" "$BENCH_SITES" "$BENCH_GROUPS" \
    "$((2 * BENCH_HALF_WIDTH))" 1 "$midpoint_policy" "$node" "$cpu_model" "$cpu_ids" \
    "$physical_core_ids" "$SLURM_JOB_ID" "$SLURM_STEP_ID" "$started_utc"
logical_cpu_count=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
if [[ "$node" != "$BENCH_NODE" ]] || ((logical_cpu_count != cores || physical_core_count != cores)); then
    printf '%s\tNA\tNA\tNA\tNA\t1\n' "$prefix" > "$BENCH_WORK/time.tsv"
    echo "Unexpected placement: node=$node, CPUs=$cpu_ids, physical cores=$physical_core_ids" >&2
    exit 1
fi
tool_exit_code=0
benchmark_time_command "$prefix" "$BENCH_WORK/output" \
    "$BENCH_WORK/tool.log" "$BENCH_WORK/tool-error.log" "${command[@]}" || tool_exit_code=$?
if ((tool_exit_code != 0)); then
    cat "$BENCH_WORK/tool-error.log" "$BENCH_WORK/tool.log" >&2
    exit "$tool_exit_code"
fi
if [[ ! -e "$profile" ]]; then
    echo "Expected midpoint profile is missing: $profile" >&2
    exit 1
fi
# Keep a native output outside the timer. Every run still recomputes its output
printf '%s\n' "$profile" > "$BENCH_WORK/profile.path"
STEP

    # Retain known identifiers even when Slurm fails before GNU time can run
    if [[ ! -s "$BENCH_WORK/time.tsv" ]]; then
        fallback=()
        for ((column=1; column<=27; column++)); do fallback+=(NA); done
        fallback[0]=$run_order
        fallback[1]=$repetition
        fallback[2]=$is_warmup
        fallback[3]=$method
        fallback[5]=$cores
        fallback[19]=$SLURM_JOB_ID
        if [[ -s "$BENCH_WORK/step.id" ]]; then fallback[20]=$(cat "$BENCH_WORK/step.id"); fi
        (IFS=$'\t'; printf '%s\n' "${fallback[*]}") > "$BENCH_WORK/time.tsv"
    fi
    benchmark_finish_run "$method" "$cores" "$repetition" "$is_warmup" "$step_exit_code"
done < "$BENCH_WORK/order.tsv"
