#!/usr/bin/env bash

# FrEIA's end-motif workflow for a single sample, as timed by the benchmark
#
# This runs the authors' Snakemake rules filter_and_index_FrEIA, FragmEndSeq and
# data_transformation directly, without Snakemake, and stops before their group
# comparison and scoring. The benchmark times this whole script as FrEIA's command
# https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/rules/2_FrEIA/1_preprocessing.smk
#
# Usage: bash freia_motifs.sh CORES
# Reads BENCH_WORK, BENCH_MIN_MAPQ, BENCH_MIN_LENGTH, BENCH_MAX_LENGTH and BENCH_CONTIGS
# The two-line samples table ($BENCH_WORK/freia-samples.tsv) is written during setup

set -euo pipefail
cores=${1:?Supply the allocated core count}
freia_environment="$BENCH_WORK/freia"
freia_scripts="$freia_environment/share/freia"
bam="$BENCH_WORK/sample.bam"
output="$BENCH_WORK/output"

# Record each native call below with its variables filled in. Bash writes the trace
# to $BENCH_WORK/command-trace.log, which the job saves in runs/<run_order>/
exec {trace_descriptor}>> "$BENCH_WORK/command-trace.log"
BASH_XTRACEFD=$trace_descriptor
set -x

# 1. The authors' BAM filter. The exclusion flags are unmapped (4), secondary (256),
#    duplicate (1024) and supplementary (2048). The filtered BAM streams through tee
#    into the indexer, as in the upstream rule, so indexing does not read the
#    finished file again. Both samtools thread options receive the allocated cores
"$freia_environment/bin/samtools" view -@ "$cores" -h -q "$BENCH_MIN_MAPQ" \
    -F 4 -F 256 -F 1024 -F 2048 \
    -b "$bam" |
    tee "$output/filtered.bam" |
    "$freia_environment/bin/samtools" index -@ "$cores" - "$output/filtered.bam.bai"

# 2. Extract each fragment's end sequences (3 bases) from every BAM contig
#    The extractor itself switches to one worker when the filtered BAM is larger
#    than 8000 MiB. The ends job records the worker count it actually used
"$freia_environment/bin/python" "$freia_scripts/1_extract_fragment_ends.py" \
    -b "$output/filtered.bam" -o "$output/sample.pq" \
    -c "$BENCH_CONTIGS" -t "$cores" -n 3

# 3. Count motifs per fragment length. FrEIA's upper length bound is exclusive,
#    hence BENCH_MAX_LENGTH + 1. No subsampling or bootstrap resampling is done
#    FrEIA adds coordinate deduplication, and removes a whole pair if either
#    extracted 3-base sequence contains N. These native choices are kept
"$freia_environment/bin/python" "$freia_scripts/2_data_transformation.py" \
    -i "$output/sample.pq" -o "$output/" \
    -st "$BENCH_WORK/freia-samples.tsv" \
    -fsmin "$BENCH_MIN_LENGTH" -fsmax "$((BENCH_MAX_LENGTH + 1))" \
    -subs 1 --bootstrap_sample 0
