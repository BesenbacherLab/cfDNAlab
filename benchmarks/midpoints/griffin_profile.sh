#!/usr/bin/env bash

# Run Griffin's two native profiling scripts for either demo or production sites
# Both calls share settings. Only sites, chromosomes and output paths vary
# Usage: bash griffin_profile.sh CORES SITES_YAML CHROM_SIZES CHROMOSOMES OUTPUT

set -euo pipefail
cores=${1:?Supply the allocated core budget}
sites_yaml=${2:?Supply the Griffin site-list YAML}
chromosome_sizes=${3:?Supply chromosome sizes}
chromosome_names=${4:?Supply comma-separated chromosome names}
output_directory=${5:?Supply the output directory}
environment="$BENCH_WORK/griffin-env"
scripts="$environment/share/griffin"
IFS=, read -r -a chromosomes <<< "$chromosome_names"
mkdir -p "$output_directory/tmp" "$output_directory/results"

# No masking, outlier removal, site-wise CNA normalization or smoothing
# Griffin's final group-wise mean-one normalization cannot be disabled
common=(--sample_name sample --tmp_dir "$output_directory/tmp"
    --sites_yaml "$sites_yaml" --griffin_scripts_dir "$scripts"
    --chrom_sizes_path "$chromosome_sizes" --chroms "${chromosomes[@]}"
    --chrom_column Chrom --position_column position --strand_column Strand
    --norm_window "-$BENCH_HALF_WIDTH" "$BENCH_HALF_WIDTH"
    --mappability_correction False --mappability_bw none --CPU "$cores")
"$environment/bin/python" "$scripts/griffin_coverage.py" "${common[@]}" \
    --bam "$BENCH_WORK/sample.bam" --GC_bias "$BENCH_WORK/unit_gc_bias.tsv" \
    --reference_genome "$BENCH_WORK/reference.fa" --mappability_bias none \
    --size_range "$BENCH_MIN_LENGTH" "$BENCH_MAX_LENGTH" --map_quality "$BENCH_MIN_MAPQ"
"$environment/bin/python" "$scripts/griffin_merge_sites.py" "${common[@]}" \
    --uncorrected_bw_path "$output_directory/tmp/sample/tmp_bigWig/sample.uncorrected.bw" \
    --GC_corrected_bw_path "$output_directory/tmp/sample/tmp_bigWig/sample.GC_corrected.bw" \
    --GC_map_corrected_bw_path none --results_dir "$output_directory/results" \
    --save_window "-$BENCH_HALF_WIDTH" "$BENCH_HALF_WIDTH" \
    --center_window "-$BENCH_HALF_WIDTH" "$BENCH_HALF_WIDTH" \
    --fft_window "-$BENCH_HALF_WIDTH" "$BENCH_HALF_WIDTH" --fft_index 1 \
    --step 1 --exclude_paths none --exclude_zero_mappability False \
    --exclude_outliers False --CNA_normalization False --individual False --smoothing False
