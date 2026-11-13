#!/usr/bin/env bash

# Each call estimates sample GC bias, then applies it to the midpoint profiles
# The benchmark times this entire script. Only reference tables are precomputed
# Usage: bash gc_corrected_profile.sh METHOD CORES
set -euo pipefail
method=${1:?Supply cfdnalab-gc or griffin-gc}
cores=${2:?Supply the allocated core budget}
bias_directory="$BENCH_WORK/output/gc-bias"
mkdir -p "$bias_directory"

# Record each native call below with its variables filled in. Bash writes the trace
# to $BENCH_WORK/command-trace.log, which the job saves in runs/<run_order>/
# griffin_profile.sh adds its own calls to the same file
exec {trace_descriptor}>> "$BENCH_WORK/command-trace.log"
BASH_XTRACEFD=$trace_descriptor
set -x

case "$method" in
    cfdnalab-gc)
        # The reference package is calculated once beforehand with ref-gc-bias
        # Sample estimation keeps native windows and chromosome selection, and
        # inherits the reference's full length range. Profiling uses the agreed
        # midpoint length range below
        "$BENCH_WORK/bin/cfdna" gc-bias \
            --bam "$BENCH_WORK/sample.bam" --output-dir "$bias_directory" \
            --ref-2bit "$BENCH_WORK/reference.2bit" \
            --ref-gc-file "$BENCH_WORK/reference-gc.zarr" \
            --min-mapq "$BENCH_MIN_MAPQ" --n-threads "$cores"
        "$BENCH_WORK/bin/cfdna" midpoints \
            --bam "$BENCH_WORK/sample.bam" --output-dir "$BENCH_WORK/output" \
            --intervals "$BENCH_WORK/sites.bed" --chromosomes "$BENCH_CHROMOSOMES" \
            --length-bins "$BENCH_MIN_LENGTH" "$((BENCH_MAX_LENGTH + 1))" \
            --min-mapq "$BENCH_MIN_MAPQ" --bin-size 15 --n-threads "$cores" \
            --blacklist "$BENCH_WORK/exclusions/encode_unified_GRCh38_exclusion_list.bed" \
            "$BENCH_WORK/exclusions/hg38_centromeres.bed" "$BENCH_WORK/exclusions/hg38_gaps.bed" \
            "$BENCH_WORK/exclusions/hg38_fix_patches.bed" "$BENCH_WORK/exclusions/hg38_alternative_haplotypes.bed" \
            --ref-2bit "$BENCH_WORK/reference.2bit" \
            --gc-file "$bias_directory/gc_bias_correction.zarr"
        ;;
    griffin-gc)
        environment="$BENCH_WORK/griffin-env"
        scripts="$environment/share/griffin"
        # These are the two native sample-bias steps from Griffin's workflow
        # GC counting has a CPU argument. Bias fitting itself is serial
        # The release config estimates GC bias over 15-500 bp, independently of
        # the profiling range. Neighbouring lengths are needed for native smoothing
        "$environment/bin/python" "$scripts/griffin_GC_counts.py" \
            --bam_file "$BENCH_WORK/sample.bam" --bam_file_name sample \
            --mappable_regions_path "$BENCH_WORK/gc-regions.bed" \
            --ref_seq "$BENCH_WORK/reference.fa" --chrom_sizes "$BENCH_WORK/bam.chrom.sizes" \
            --out_dir "$bias_directory" --map_q "$BENCH_MIN_MAPQ" \
            --size_range 15 500 --CPU "$cores"
        "$environment/bin/python" "$scripts/griffin_GC_bias.py" \
            --bam_file_name sample --mappable_name "$BENCH_GC_REFERENCE_NAME" \
            --genome_GC_frequency "$BENCH_WORK/reference-gc" \
            --out_dir "$bias_directory" --size_range 15 500
        bash "$BENCH_WORK/griffin_profile.sh" "$cores" "$BENCH_WORK/sites.yaml" \
            "$BENCH_WORK/chrom.sizes" "$BENCH_CHROMOSOMES" "$BENCH_WORK/output" \
            "$bias_directory/GC_bias/sample.GC_bias.txt"
        ;;
    *) echo "Unknown GC-corrected midpoint method: $method" >&2; exit 1 ;;
esac
