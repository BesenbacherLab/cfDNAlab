#!/usr/bin/env bash
#SBATCH --job-name=cfdna-midpoints
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --no-requeue

# This job produces native, unsmoothed 1 bp midpoint profiles
# The submission launcher requests an exclusive node for full runs and allows shared pilots
# Full-input runs get timing rows. Griffin's short demo warmup is untimed
# A profile is kept from the first measured repetition of each tool
# All BAM/reference copies, environments and temporary files use local scratch
# See README.md for required paths, concordance and expected differences

set -euo pipefail
export LC_ALL=C
source "${SLURM_SUBMIT_DIR:?Submit from the repository root}/benchmarks/common.sh"
BENCH_FEATURE=midpoints

: "${SLURM_JOB_ID:?Submit with bash benchmarks/submit_benchmark.sh midpoints}"
: "${BENCH_BAM:?Set the prepared coordinate-sorted BAM path}"
: "${BENCH_INTERVALS:?Set the grouped fixed-size BED path}"
# TMPDIR is supplied on the compute node when the batch job starts
: "${TMPDIR:?The compute node must provide TMPDIR for node-local scratch}"
: "${BENCH_RESULTS:?Set a new results directory}"
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM}.bai}
BENCH_FASTA_INDEX=${BENCH_FASTA_INDEX:-${BENCH_REF_FASTA:-}.fai}
: "${BENCH_CORES:?Set a single core budget for this job: 1, 4, 12 or 32}"
BENCH_WARMUPS=${BENCH_WARMUPS:-1}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-3}
: "${BENCH_METHODS:?Set the methods for this job explicitly in BENCH_METHODS}"
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-100}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-200}

benchmark_controls
read -r -a methods_to_run <<< "$BENCH_METHODS"
if [[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 || -e "$BENCH_RESULTS" ]]; then
    echo "Use a single-node allocation and a new BENCH_RESULTS directory." >&2
    exit 1
fi
if ((${#methods_to_run[@]} == 0)); then
    echo "Supply at least a tool in BENCH_METHODS and a core budget in BENCH_CORES." >&2
    exit 1
fi
for value in "$BENCH_WARMUPS" "$BENCH_REPETITIONS" "$BENCH_MIN_MAPQ" "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH"; do
    if [[ ! "$value" =~ ^(0|[1-9][0-9]*)$ ]]; then
        echo "Core counts, repetitions and filter settings must be integers: $value" >&2
        exit 1
    fi
done
if ((BENCH_REPETITIONS < 1 || BENCH_MIN_MAPQ > 255 || BENCH_MIN_LENGTH < 10 || BENCH_MAX_LENGTH < BENCH_MIN_LENGTH || BENCH_MAX_LENGTH > 50000)); then
    echo "Use MAPQ 0-255, lengths within 10-50000 bp, and at least one measured repetition." >&2
    exit 1
fi
if ((BENCH_CORES > ${SLURM_CPUS_PER_TASK:-0})); then
    echo "Core budget $BENCH_CORES exceeds allocated CPUs: ${SLURM_CPUS_PER_TASK:-0}" >&2
    exit 1
fi
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
if [[ "$run_griffin" == true ]]; then
    if ((${#methods_to_run[@]} != 1)); then
        echo "Run Griffin in its own job: set BENCH_METHODS=griffin." >&2
        exit 1
    fi
    # The release demo replaces full-input warmups, including in pilot jobs
    BENCH_WARMUPS=0
fi
: "${BENCH_SAMTOOLS_ARCHIVE:?Set the packed samtools environment for BAM header preparation}"
if [[ "$run_cfdna" == true ]]; then : "${BENCH_CFDNA:?Set the cfDNAlab executable}"; fi
if [[ "$run_griffin" == true ]]; then : "${BENCH_GRIFFIN_ARCHIVE:?Set the packed Griffin environment}"; : "${BENCH_REF_FASTA:?Set the matching reference FASTA}"; fi

# Inputs are copied once before warmups, including the supplied BAM/reference indexes
export BENCH_WORK BENCH_NODE BENCH_HALF_WIDTH BENCH_SITES BENCH_GROUPS BENCH_CHROMOSOMES
export BENCH_RESULTS BENCH_CFDNA_VERSION BENCH_GRIFFIN_VERSION
BENCH_WORK=$(mktemp -d "$TMPDIR/midpoints-${SLURM_JOB_ID}.XXXXXX")
trap 'rm -rf "$BENCH_WORK"' EXIT
BENCH_NODE=$(hostname)
mkdir -p "$BENCH_WORK/bin" "$BENCH_WORK/tmp" "$BENCH_WORK/cache" "$BENCH_WORK/matplotlib" "$BENCH_WORK/griffin-env"
cp -L "$BENCH_BAM" "$BENCH_WORK/sample.bam"
cp -L "$BENCH_INDEX" "$BENCH_WORK/sample.bam.${BENCH_INDEX##*.}"
cp -L "$BENCH_INTERVALS" "$BENCH_WORK/input.bed"
if [[ "$run_griffin" == true ]]; then
    cp -L "$BENCH_REF_FASTA" "$BENCH_WORK/reference.fa"
    cp -L "$BENCH_FASTA_INDEX" "$BENCH_WORK/reference.fa.fai"
    cp "$SLURM_SUBMIT_DIR/benchmarks/midpoints/griffin_profile.sh" "$BENCH_WORK/griffin_profile.sh"
fi
if [[ "$run_cfdna" == true ]]; then benchmark_stage_binary "$BENCH_CFDNA" cfdna; fi
# User setup: submit from the repository root so SLURM_SUBMIT_DIR locates helpers
# Slurm keeps its own copy of this job script elsewhere, away from those helpers
export TMPDIR="$BENCH_WORK/tmp" XDG_CACHE_HOME="$BENCH_WORK/cache"
# Python user packages and outside import paths are disabled so imports use
# Griffin's copied environment. The installation notes in ../README.md explain why
export MPLCONFIGDIR="$BENCH_WORK/matplotlib" PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
benchmark_unpack "$BENCH_SAMTOOLS_ARCHIVE" "$BENCH_WORK/samtools-env"
"$BENCH_WORK/samtools-env/bin/samtools" idxstats "$BENCH_WORK/sample.bam" |
    awk '$1 != "*" {print $1 "\t" $2}' > "$BENCH_WORK/chrom.sizes"
if [[ "$run_griffin" == true ]]; then
    benchmark_unpack "$BENCH_GRIFFIN_ARCHIVE" "$BENCH_WORK/griffin-env"
    # pybedtools must find the packed bedtools, rather than a shared installation
    export PATH="$BENCH_WORK/griffin-env/bin:$PATH"
fi
export BENCH_RUN_GRIFFIN=$run_griffin

# Setup converts the same oriented windows to cfDNAlab and Griffin input formats
# All sites are retained, including repeated and overlapping rows. Invalid input fails
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
        # Strand parsing follows cfDNAlab's BED5/BED6 convention and ambiguity checks
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
# The fixed release demo gives Griffin a short warmup without changing the
# supplied production windows. It uses the same copied BAM and reference
if os.environ['BENCH_RUN_GRIFFIN'] == 'true':
    demo_path = work / 'griffin-env/share/griffin/CTCF.hg38.1000.txt'
    demo_rows = []
    with demo_path.open() as source:
        for row in csv.DictReader(source, delimiter='\t'):
            chromosome, position = row['Chrom'], int(row['position'])
            if chromosome not in bam_contigs or reference_contigs.get(chromosome) != bam_contigs[chromosome]:
                raise ValueError('Griffin demo contig is absent or has mismatched lengths: ' + chromosome)
            if position - window_width // 2 < 0 or position + window_width // 2 > bam_contigs[chromosome]:
                raise ValueError('Griffin demo window extends beyond its chromosome: ' + chromosome)
            demo_rows.append((chromosome, position, '.'))
    if len(demo_rows) != 1000:
        raise ValueError('Expected exactly 1000 sites in the pinned Griffin demo')
    with (work / 'warmup-sites.tsv').open('w') as destination:
        writer = csv.writer(destination, delimiter='\t', lineterminator='\n')
        writer.writerow(('Chrom', 'position', 'Strand'))
        writer.writerows(demo_rows)
    with (work / 'warmup-sites.yaml').open('w') as destination:
        json.dump({'site_lists': {'CTCF-demo': str(work / 'warmup-sites.tsv')}}, destination)
    with (work / 'warmup.chrom.sizes').open('w') as destination:
        for chromosome in sorted(set(row[0] for row in demo_rows)):
            destination.write('{}\t{}\n'.format(chromosome, bam_contigs[chromosome]))
chromosomes = sorted(set(row[0] for row in rows))
# Griffin sorts intervals lexicographically. The bigWig header uses that same
# contig order so exported entries remain correctly ordered
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

# The BAM, sites and programs are ready on the compute node. Preparation time
# above and Griffin's short demo below are excluded from tool timings
# Each full-input run starts through srun. Once it finishes, this main script saves
# its results. common.sh explains the temporary files used to pass them back
benchmark_record_setup
if [[ "$run_griffin" == true ]]; then
    printf 'griffin_warmup\trelease-demo-1000-sites\n' >> "$BENCH_RESULTS/run.tsv"
    sha256sum "$BENCH_WORK/griffin-env/share/griffin/CTCF.hg38.1000.txt" > "$BENCH_RESULTS/griffin-warmup.sha256"
    warmup_chromosomes=$(cut -f 1 "$BENCH_WORK/warmup.chrom.sizes" | paste -sd, -)
    # A separate untimed step keeps demo work out of measured Slurm peak memory
    # This gives Griffin a practical warmup without another full production run
    srun --jobid="$SLURM_JOB_ID" --exact --nodes=1 --ntasks=1 \
        --cpus-per-task="$BENCH_CORES" --threads-per-core=1 --cpu-bind=threads \
        --export=ALL bash "$BENCH_WORK/griffin_profile.sh" "$BENCH_CORES" \
        "$BENCH_WORK/warmup-sites.yaml" "$BENCH_WORK/warmup.chrom.sizes" \
        "$warmup_chromosomes" "$BENCH_WORK/warmup-output" \
        > "$BENCH_RESULTS/griffin-warmup.log" 2> "$BENCH_RESULTS/griffin-warmup-error.log" || {
            cat "$BENCH_RESULTS/griffin-warmup-error.log" "$BENCH_RESULTS/griffin-warmup.log" >&2
            echo "Griffin's demo warmup failed. Measurements have not started." >&2
            exit 1
        }
    rm -rf "$BENCH_WORK/warmup-output"
    sync -f "$BENCH_WORK"
fi
benchmark_schedule > "$BENCH_WORK/order.tsv"
if [[ ! -s "$BENCH_WORK/order.tsv" ]]; then echo "The requested run matrix is empty." >&2; exit 1; fi
run_order=0
job_exit_code=0
failed_method=
while IFS=$'\t' read -r method cores repetition is_warmup; do
    if [[ "$method" == "$failed_method" ]]; then continue; fi
    run_order=$((run_order + 1))
    step_exit_code=0
    # Old timings and output paths are removed before each run. If this run fails
    # before writing its own files, it cannot accidentally reuse the preceding results
    rm -f "$BENCH_WORK/time.tsv" "$BENCH_WORK/io.tsv" "$BENCH_WORK/profile.path" "$BENCH_WORK/step.id" "$BENCH_WORK/stdout.path"
    # srun starts the next run on the reserved node with the requested cores. The code between
    # <<'STEP' and the closing STEP line does this run's work, then this main script
    # resumes and saves its results. Quoting STEP delays variable expansion until
    # that program starts. --export=ALL passes paths/settings and the arguments
    # identify the run. Results come back in files under BENCH_WORK
    # The code below checks the assigned CPUs because enforcement depends on cluster settings
    # srun's exit code is saved so a failed run still appears in the final timing table
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
# Hyperthreading can give a physical core two CPU IDs. The script looks up
# the actual core behind each allowed ID to verify the requested core count
# Core numbers can repeat across CPU chips, so chip (socket) and core number
# together identify a core
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
    # Timing includes both Griffin scripts, bigWig writing and group aggregation
    command=(bash "$BENCH_WORK/griffin_profile.sh" "$cores" "$BENCH_WORK/sites.yaml"
        "$BENCH_WORK/chrom.sizes" "$BENCH_CHROMOSOMES" "$BENCH_WORK/output")
    profile="$BENCH_WORK/output/results/sample/sample.uncorrected.coverage.tsv"
fi

# The first 22 fields identify the run, followed by four GNU time values and status
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
# The output path is passed back for preservation after timing. Every run recomputes it
printf '%s\n' "$profile" > "$BENCH_WORK/profile.path"
STEP

    # A run that failed before timing still gets a row naming the attempted
    # method and repetition. Unknown measurements are recorded as NA
    # Column numbers below start at 0 and must match the timing-table header above
    # benchmark_finish_run adds disk-write measurements and the srun exit code
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

# The job reports failure if any method failed, while keeping completed results
exit "$job_exit_code"
