#!/usr/bin/env bash

# The native midpoint profile commands, one per method
#
# benchmark_midpoints.sh loads this file. Before every run it calls prepare_command,
# and after the first measured repetition it calls save_profile. Each method makes
# a profile of fragment midpoints around the supplied sites, in 15 bp bins with
# native 165 bp smoothing. Every midpoint job runs a single method
#
# The commands use these variables, which the job script sets after copying the
# inputs and tools to the compute node. experiment_commands.md uses the same names
#   cores              CPU cores allocated to every run (BENCH_CORES)
#   bam                The copied input BAM
#   output             The run's empty output directory, $BENCH_WORK/output
#   min_mapq, min_length, max_length
#                      Profiling filters. Length bounds are inclusive
#   length_end         max_length + 1, for options with an exclusive upper bound
#   exclusions         The five shared exclusion BEDs (an array)
#   cfdna              cfDNAlab's path
#   cfdna_version, griffin_version
#                      Tool versions recorded at setup
#   BENCH_CHROMOSOMES  Comma-separated chromosomes with sites
#   BENCH_SITES, BENCH_GROUPS, BENCH_WINDOW_WIDTH
#                      Site count, group count and window width, from prepare_sites.py
# The GC-corrected methods run gc_corrected_profile.sh, and Griffin's profiling
# runs griffin_profile.sh. Both read the BENCH_* settings from the environment

# The columns of timings.tsv, in order. common.sh fills the run's identity and
# measurements, and prepare_command below fills the method's settings
#   tool_workers       Value passed to the tool's thread/CPU option
#   sites, groups      Rows and groups in the supplied BED
#   window_bp          Width of every site window
#   resolution_bp      Bin size of the profile
#   midpoint_policy    Which central base gets an even-length fragment's midpoint
#   gc_correction      none, or corrected for the GC-corrected workloads
#   cpu_ids, physical_core_ids  The CPUs the Slurm step was bound to
#   wall_seconds ... max_rss_kib  GNU time's measurements of the native command
#   flush_seconds      Time spent waiting for the output to reach the disk
#   tool_exit_code     The native command's exit status
#   step_exit_code     The Slurm step's exit status, which is also nonzero for a
#                      wrong CPU placement or a missing output
timing_columns=(
    run_order repetition warmup method tool_version cores tool_workers
    min_mapq min_length max_length sites groups window_bp resolution_bp midpoint_policy
    node cpu_model cpu_ids physical_core_ids job_id step_id started_utc
    gc_correction
    wall_seconds user_seconds system_seconds max_rss_kib tool_exit_code
    flush_seconds process_plus_flush_seconds output_bytes pilot step_exit_code
)

# Sets the native command for a method and its settings for timings.tsv
#
# Parameters:
#   $1  Method name
# Sets (see common.sh for how they are used):
#   command, expected_output, empty_output_allowed, stdout_file, row[...]
prepare_command() {
    local method=$1
    # Defaults shared by all methods. The branches below change only what differs
    #   tool_workers     Value passed to the tool's thread/CPU option
    #   sites, groups    Rows and groups in the supplied BED
    #   window_bp        Width of every site window
    #   midpoint_policy  Which central base gets an even-length fragment's midpoint
    row=(
        [tool_workers]=$cores
        [min_mapq]=$min_mapq
        [min_length]=$min_length
        [max_length]=$max_length
        [sites]=$BENCH_SITES
        [groups]=$BENCH_GROUPS
        [window_bp]=$BENCH_WINDOW_WIDTH
        [resolution_bp]=15
        [gc_correction]=none
    )
    stdout_file="$BENCH_WORK/tool-stdout.log"
    empty_output_allowed=false

    case "$method" in
        cfdnalab)
            # Uncorrected profiles. They are compared with the raw table that the
            # griffin-gc job also writes. --blacklist drops whole intervals near the
            # exclusion BEDs, so cfDNAlab can report fewer eligible sites than Griffin
            command=("$cfdna" midpoints --bam "$bam" --output-dir "$output"
                --intervals "$BENCH_WORK/sites.bed" --chromosomes "$BENCH_CHROMOSOMES"
                --length-bins "$min_length" "$length_end"
                --min-mapq "$min_mapq" --bin-size 15 --n-threads "$cores"
                --blacklist "${exclusions[@]}")
            expected_output="$output/midpoint_profiles.zarr"
            row[tool_version]=$cfdna_version
            row[midpoint_policy]=coordinate-random-even
            ;;
        cfdnalab-gc)
            # Estimates the sample's GC bias, then makes GC-corrected profiles
            # Both cfDNAlab calls are in gc_corrected_profile.sh and share one timer
            command=(bash "$BENCH_WORK/gc_corrected_profile.sh" cfdnalab-gc "$cores")
            expected_output="$output/midpoint_profiles.zarr"
            row[tool_version]=$cfdna_version
            row[midpoint_policy]=coordinate-random-even
            row[gc_correction]=corrected
            ;;
        griffin-gc)
            # Griffin's presented workflow: sample GC estimation, then profiling
            # All four Griffin calls are in gc_corrected_profile.sh and
            # griffin_profile.sh and share one timer. Profiling writes both the
            # uncorrected and the GC-corrected table
            command=(bash "$BENCH_WORK/gc_corrected_profile.sh" griffin-gc "$cores")
            expected_output="$output/results/sample/sample.GC_corrected.coverage.tsv"
            row[tool_version]=$griffin_version
            row[midpoint_policy]=right-central-even
            row[gc_correction]=corrected
            ;;
        *)
            echo "Unknown midpoint tool: $method" >&2
            return 1
            ;;
    esac
}

# Copies a method's output into profiles/METHOD, for compare_midpoints.py
#
# Only the first measured repetition is saved. Every repetition still computes its
# full output, including a fresh sample GC estimate, so saving one copy reduces
# stored data without changing the timed work
#
# Parameters:
#   $1  Method name
save_profile() {
    local method=$1
    local destination="$BENCH_RESULTS/profiles/$method"
    mkdir -p "$destination"
    case "$method" in
        cfdnalab)
            cp -R "$output/midpoint_profiles.zarr" "$destination/"
            ;;
        cfdnalab-gc)
            # The corrected profiles are kept with the sample GC estimate behind them
            cp -R "$output/midpoint_profiles.zarr" "$destination/"
            cp -R "$output/gc-bias" "$destination/sample-gc-bias"
            ;;
        griffin-gc)
            # Both of Griffin's tables are kept. The uncorrected table supplies the
            # uncorrected concordance without another Griffin job
            cp "$output/results/sample/sample.GC_corrected.coverage.tsv" "$destination/"
            cp "$output/results/sample/sample.uncorrected.coverage.tsv" "$destination/"
            cp -R "$output/gc-bias" "$destination/sample-gc-bias"
            ;;
        *)
            echo "No saved profile defined for midpoint method: $method" >&2
            return 1
            ;;
    esac
}
