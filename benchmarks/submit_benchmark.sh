#!/usr/bin/env bash

# Submit a benchmark from the repository root after sourcing the input/tool settings
# Usage: bash benchmarks/submit_benchmark.sh lengths [sbatch options]
# BENCH_PILOT=1 allows a shared node and requests 50 GB
# Full runs reserve the whole node and its memory
# Extra arguments go to sbatch, for example --partition, --time or a pilot's --mem

set -euo pipefail

feature=${1:?Choose coverage, lengths, ends or midpoints, followed by any sbatch options}
shift
case "$feature" in
    coverage|lengths|midpoints) job_script="benchmarks/$feature/benchmark_${feature}.sh" ;;
    ends) job_script=benchmarks/ends/benchmark_end_motifs.sh ;;
    *) echo "Unknown benchmark: $feature. Choose coverage, lengths, ends or midpoints." >&2; exit 1 ;;
esac
if [[ ! -f "$job_script" || ! -f benchmarks/common.sh ]]; then
    echo 'Run this launcher from the repository root.' >&2
    exit 1
fi
: "${BENCH_CORES:?Set a single core budget for this job: 1, 4, 12 or 32}"
if [[ ! "$BENCH_CORES" =~ ^[1-9][0-9]*$ ]]; then
    echo 'Set BENCH_CORES to a single positive core count for this job.' >&2
    exit 1
fi
export BENCH_PILOT=${BENCH_PILOT:-0}
if [[ "$BENCH_PILOT" != 0 && "$BENCH_PILOT" != 1 ]]; then
    echo 'BENCH_PILOT must be 0 or 1.' >&2
    exit 1
fi

# Catch an existing midpoint results directory before queueing. Failed attempts
# can leave output behind, and the runner preserves it rather than overwriting it
if [[ "$feature" == midpoints ]]; then
    # Check required reference files before queueing a long midpoint job
    # These are existence checks only. The native tools validate their contents
    : "${BENCH_GRIFFIN_REF_DIR:?Set the Griffin Ref directory containing the shared exclusion BEDs}"
    for exclusion_file in encode_unified_GRCh38_exclusion_list.bed hg38_centromeres.bed hg38_gaps.bed hg38_fix_patches.bed hg38_alternative_haplotypes.bed; do
        if [[ ! -f "$BENCH_GRIFFIN_REF_DIR/$exclusion_file" ]]; then
            echo "Missing midpoint exclusion BED: $BENCH_GRIFFIN_REF_DIR/$exclusion_file" >&2
            exit 1
        fi
    done
    if [[ " ${BENCH_METHODS:-} " == *" griffin"* ]]; then
        griffin_mappability=${BENCH_GRIFFIN_MAPPABILITY_BW:-$BENCH_GRIFFIN_REF_DIR/k100.Umap.MultiTrackMappability.bw}
        if [[ ! -f "$griffin_mappability" ]]; then
            echo "Missing Griffin mappability track: $griffin_mappability. Set BENCH_GRIFFIN_MAPPABILITY_BW to the hg38 k100 Umap bigWig." >&2
            exit 1
        fi
    fi
    : "${BENCH_RESULTS:?Set a new results directory}"
    if [[ -e "$BENCH_RESULTS" ]]; then
        echo "BENCH_RESULTS already exists: $BENCH_RESULTS" >&2
        echo "Choose a new path without creating it first. The benchmark creates the directory." >&2
        exit 1
    fi
    # Check only the first five lines for the expected grouped BED format
    # The tools remain responsible for input validation
    : "${BENCH_INTERVALS:?Set the grouped fixed-size BED path}"
    "${BENCH_SETUP_PYTHON:-python3}" - "$BENCH_INTERVALS" <<'PY'
import csv
from itertools import islice
import sys

with open(sys.argv[1]) as source:
    for line_number, fields in enumerate(csv.reader(islice(source, 5), delimiter='\t'), 1):
        if not fields or fields[0].startswith('#'):
            continue
        if len(fields) < 4:
            sys.exit('BED row {} needs tab-separated chromosome, start, end and group'.format(line_number))
        try:
            int(fields[1])
            int(fields[2])
        except ValueError:
            sys.exit('BED row {} needs integer start/end coordinates and no header'.format(line_number))
PY
fi

# Slurm reads resource requests before the batch script starts, so the pilot choice
# must be made here. Conditional #SBATCH lines inside the job would not work
# Pilots request 50 GB unless --mem is supplied. User options follow the defaults
# so an explicit memory request takes precedence
sbatch_options=(--export=ALL --cpus-per-task="$BENCH_CORES")
if [[ "$BENCH_PILOT" == 0 ]]; then
    sbatch_options+=(--exclusive --mem=0)
else
    sbatch_options+=(--mem=50G)
fi
exec sbatch "${sbatch_options[@]}" "$@" "$job_script"
