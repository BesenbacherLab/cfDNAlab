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
