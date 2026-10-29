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
# Keep a native profile from the first repetition for each tool
# All BAM/reference copies, environments and temporary files use local scratch
# See README.md for required paths, concordance and expected differences

set -euo pipefail
export LC_ALL=C

: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${BENCH_BAM:?Set the prepared coordinate-sorted BAM path}"
: "${BENCH_INTERVALS:?Set the grouped fixed-size BED path}"
: "${BENCH_REF_FASTA:?Set the matching uncompressed reference FASTA path}"
: "${BENCH_CFDNA:?Set the cfDNAlab executable path}"
: "${BENCH_GRIFFIN_ARCHIVE:?Set the packed Griffin environment}"
: "${BENCH_LOCAL_ROOT:?Set an existing node-local scratch directory}"
: "${BENCH_RESULTS:?Set a new results directory}"
BENCH_INDEX=${BENCH_INDEX:-${BENCH_BAM}.bai}
BENCH_FASTA_INDEX=${BENCH_FASTA_INDEX:-${BENCH_REF_FASTA}.fai}
BENCH_CORES=${BENCH_CORES:-1 4 8 12 32}
BENCH_WARMUPS=${BENCH_WARMUPS:-3}
BENCH_REPETITIONS=${BENCH_REPETITIONS:-10}
BENCH_METHODS=${BENCH_METHODS:-cfdnalab griffin}
export BENCH_MIN_MAPQ=${BENCH_MIN_MAPQ:-30}
export BENCH_MIN_LENGTH=${BENCH_MIN_LENGTH:-100}
export BENCH_MAX_LENGTH=${BENCH_MAX_LENGTH:-200}

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
    if [[ ! "$value" =~ ^[0-9]+$ ]]; then
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
if [[ "$BENCH_REF_FASTA" == *.gz ]]; then
    echo "Use an uncompressed FASTA and its existing .fai index." >&2
    exit 1
fi

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
cp -L "$BENCH_REF_FASTA" "$BENCH_WORK/reference.fa"
cp -L "$BENCH_FASTA_INDEX" "$BENCH_WORK/reference.fa.fai"
cp -L "$BENCH_CFDNA" "$BENCH_WORK/bin/cfdna"
# sbatch uses a spool copy of this script. Submit from the repository root
cp "$SLURM_SUBMIT_DIR/benchmarks/measure_memory.py" "$BENCH_WORK/measure_memory.py"
export TMPDIR="$BENCH_WORK/tmp" XDG_CACHE_HOME="$BENCH_WORK/cache"
# Ignore ~/.local Python packages and inherited import paths so imports stay
# inside Griffin's staged environment. See the installation notes in ../README.md
export MPLCONFIGDIR="$BENCH_WORK/matplotlib" PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
tar -xzf "$BENCH_GRIFFIN_ARCHIVE" -C "$BENCH_WORK/griffin-env"
"$BENCH_WORK/griffin-env/bin/python" "$BENCH_WORK/griffin-env/bin/conda-unpack"
# pybedtools must find the packed bedtools, rather than a shared installation
export PATH="$BENCH_WORK/griffin-env/bin:$PATH"

# Prepare two representations of exactly the same oriented windows
# Keep all sites, including repeated and overlapping rows. Fail on invalid input
# This setup reads headers and BED rows, without scanning or filtering BAM reads
site_summary=$("$BENCH_WORK/griffin-env/bin/python" - <<'PY'
import csv
import os
from pathlib import Path

import pysam
import yaml

work = Path(os.environ['BENCH_WORK'])
with pysam.AlignmentFile(str(work / 'sample.bam'), 'rb') as bam:
    bam_contigs = dict(zip(bam.references, bam.lengths))
with pysam.FastaFile(str(work / 'reference.fa')) as reference:
    reference_contigs = dict(zip(reference.references, reference.lengths))
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
        strand = strand if strand in ('+', '-', '.') else '.'
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
    path = work / 'sites' / '{}.tsv'.format(group_index)
    with path.open('w') as destination:
        writer = csv.writer(destination, delimiter='\t', lineterminator='\n')
        writer.writerow(('Chrom', 'position', 'Strand'))
        writer.writerows(sites)
    site_paths[group] = str(path)
with (work / 'sites.yaml').open('w') as destination:
    yaml.safe_dump({'site_lists': site_paths}, destination)
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
BENCH_CFDNA_VERSION=$("$BENCH_WORK/bin/cfdna" --version 2>&1 | tr '\t\n' '  ')
BENCH_GRIFFIN_VERSION=$(cat "$BENCH_WORK/griffin-env/share/griffin/version.txt")
mkdir -p "$BENCH_RESULTS/profiles"

# One row per tool run, with the actual node and allowed CPU IDs
header=$'run_order\trepetition\twarmup\tmethod\ttool_version\tcores\ttool_workers'
header+=$'\tmin_mapq\tmin_length\tmax_length\tsites\tgroups\twindow_bp\tresolution_bp\tmidpoint_policy'
header+=$'\tnode\tcpu_model\tcpu_ids\tphysical_core_ids\tjob_id\tstep_id\tstarted_utc'
header+=$'\twall_seconds\tuser_seconds\tsystem_seconds\tmax_rss_kib\ttool_exit_code'
header+=$'\tsampled_peak_rss_kib\tsampled_mean_rss_kib\tmemory_samples\tmemory_interval_seconds\tstep_exit_code'
printf '%s\n' "$header" > "$BENCH_RESULTS/timings.tsv"
for method in "${methods_to_run[@]}"; do
    for cores in "${cores_to_run[@]}"; do
        printf '%s\t%s\n' "$method" "$cores"
    done
done > "$BENCH_WORK/combinations.tsv"

# Shuffle within each repetition block. All warmup rows remain in the TSV
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
        # Attach to this job's existing node, then check the binding in the step
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
"$BENCH_WORK/griffin-env/bin/python" "$BENCH_WORK/measure_memory.py" "$BENCH_WORK/memory.tsv" \
    /usr/bin/time --quiet --format "$prefix\t%e\t%U\t%S\t%M\t%x" \
    --output "$BENCH_WORK/time.tsv" "${command[@]}" \
    > "$BENCH_WORK/tool.log" 2>&1 || tool_exit_code=$?
if ((tool_exit_code != 0)); then
    cat "$BENCH_WORK/tool.log" >&2
    exit "$tool_exit_code"
fi
if [[ ! -e "$profile" ]]; then
    echo "Expected midpoint profile is missing: $profile" >&2
    exit 1
fi
# Keep a native output outside the timer. Every run still recomputes its output
if ((repetition == 1 && cores == BENCH_CAPTURE_CORES)); then
    mkdir -p "$BENCH_RESULTS/profiles/$method"
    cp -R "$profile" "$BENCH_RESULTS/profiles/$method/"
fi
STEP
        # Preserve completed runs immediately, and retain a row for failed steps
        if [[ ! -s "$BENCH_WORK/time.tsv" ]]; then
            printf '%s\t%s\t%s\t%s\tNA\t%s' \
                "$run_order" "$repetition" "$is_warmup" "$method" "$cores" > "$BENCH_WORK/time.tsv"
            for ((column=7; column<=27; column++)); do printf '\tNA' >> "$BENCH_WORK/time.tsv"; done
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
        rm -rf "$BENCH_WORK/output" "$BENCH_WORK/tool.log"
    done < "$BENCH_WORK/order.tsv"
done
