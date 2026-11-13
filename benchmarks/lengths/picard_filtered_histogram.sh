#!/usr/bin/env bash

# Makes Picard's fragment length histogram with the shared MAPQ and flag filters
#
# The timed Picard run reads the BAM directly, as Picard documents it, and Picard
# has no MAPQ option. Concordance needs the same filters as the other tools, so
# benchmark_lengths.sh runs this once after all measurements, in its own untimed
# Slurm step. samtools keeps paired reads (-f 1) with MAPQ >= BENCH_MIN_MAPQ and
# without the excluded flags (-F 3852), and streams them to Picard as uncompressed BAM
#
# Usage: bash picard_filtered_histogram.sh CORES OUTPUT_FILE
# Reads BENCH_WORK, BENCH_MIN_MAPQ and BENCH_MAX_LENGTH

set -euo pipefail
cores=${1:?Supply the allocated core count}
output_file=${2:?Supply the output metrics path}

"$BENCH_WORK/samtools-env/bin/samtools" view -@ "$cores" -u -q "$BENCH_MIN_MAPQ" -f 1 -F 3852 "$BENCH_WORK/sample.bam" |
    "$BENCH_WORK/picard-env/bin/picard" -Xmx2g CollectInsertSizeMetrics \
        I=/dev/stdin "O=$output_file" \
        "HISTOGRAM_WIDTH=$BENCH_MAX_LENGTH" "TMP_DIR=$BENCH_WORK/tmp"
