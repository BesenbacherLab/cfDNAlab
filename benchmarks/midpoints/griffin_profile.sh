#!/usr/bin/env bash

# This script runs Griffin's two profiling scripts for demo or production sites
# Demo and production calls share filters and window settings. Sites, chromosomes,
# output paths and the supplied GC bias can vary
# Usage: bash griffin_profile.sh CORES SITES_YAML CHROM_SIZES CHROMOSOMES OUTPUT [GC_BIAS]

set -euo pipefail
cores=${1:?Supply the allocated core budget}
sites_yaml=${2:?Supply the Griffin site-list YAML}
chromosome_sizes=${3:?Supply chromosome sizes}
chromosome_names=${4:?Supply comma-separated chromosome names}
output_directory=${5:?Supply the output directory}
# The untimed demo uses unit weights. Production passes freshly estimated bias
gc_bias=${6:-$BENCH_WORK/unit_gc_bias.tsv}
environment="$BENCH_WORK/griffin-env"
scripts="$environment/share/griffin"
# End coordinates are exclusive. The supplied wide window defines
# normalization. Saved output follows Griffin's official central window, rounded
# inward from [-1000, 1000) to [-990, 990) at its 15 bp step
window_end=$((BENCH_WINDOW_WIDTH - BENCH_HALF_WIDTH))
IFS=, read -r -a chromosomes <<< "$chromosome_names"
mkdir -p "$output_directory/tmp" "$output_directory/results"

# Keep the profiling workflow's native masks, outlier removal and 165 bp smoothing
# Site-wise CNA correction stays off, as in the release config. FFT features are
# compulsory in the merge script and remain timed with the official feature settings
exclusion_paths=(
    "$BENCH_WORK/exclusions/encode_unified_GRCh38_exclusion_list.bed"
    "$BENCH_WORK/exclusions/hg38_centromeres.bed"
    "$BENCH_WORK/exclusions/hg38_gaps.bed"
    "$BENCH_WORK/exclusions/hg38_fix_patches.bed"
    "$BENCH_WORK/exclusions/hg38_alternative_haplotypes.bed"
)
common=(--sample_name sample --tmp_dir "$output_directory/tmp"
    --sites_yaml "$sites_yaml" --griffin_scripts_dir "$scripts"
    --chrom_sizes_path "$chromosome_sizes" --chroms "${chromosomes[@]}"
    --chrom_column Chrom --position_column position --strand_column Strand
    --norm_window "-$BENCH_HALF_WIDTH" "$window_end"
    --mappability_correction False --mappability_bw "$BENCH_WORK/mappability.bw" --CPU "$cores")
"$environment/bin/python" "$scripts/griffin_coverage.py" "${common[@]}" \
    --bam "$BENCH_WORK/sample.bam" --GC_bias "$gc_bias" \
    --reference_genome "$BENCH_WORK/reference.fa" --mappability_bias none \
    --size_range "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" --map_quality "$BENCH_MIN_MAPQ"
"$environment/bin/python" "$scripts/griffin_merge_sites.py" "${common[@]}" \
    --uncorrected_bw_path "$output_directory/tmp/sample/tmp_bigWig/sample.uncorrected.bw" \
    --GC_corrected_bw_path "$output_directory/tmp/sample/tmp_bigWig/sample.GC_corrected.bw" \
    --GC_map_corrected_bw_path none --results_dir "$output_directory/results" \
    --save_window -1000 1000 \
    --center_window -30 30 --fft_window -960 960 --fft_index 10 \
    --step 15 --exclude_paths "${exclusion_paths[@]}" --exclude_zero_mappability True \
    --exclude_outliers True --CNA_normalization False --individual False \
    --smoothing True --smoothing_length 165
