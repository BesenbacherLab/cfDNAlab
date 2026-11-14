#!/usr/bin/env bash

# The native fragment length commands, one per method
#
# benchmark_lengths.sh loads this file. Before every run it calls prepare_command,
# and after a method's first measured repetition it calls save_profile. Each method
# counts fragments per integer length from min_length through max_length, inclusive
#
# The commands use these variables, which the job script sets after copying the
# inputs and tools to the compute node. experiment_commands.md uses the same names
#   cores              CPU cores allocated to every run (BENCH_CORES)
#   bam                The copied input BAM
#   output             The run's empty output directory, $BENCH_WORK/output
#   min_mapq, min_length, max_length
#                      Shared filters. Length bounds are inclusive
#   length_end         max_length + 1, for options with an exclusive upper bound
#   genome             cfDNAPro's genome label (BENCH_GENOME)
#   cfdna, samtools, picard, finale_python, finale, apptainer, cfdnapro_image
#                      Tool paths
#   cfdna_version, samtools_version, picard_version, finale_version, cfdnapro_version
#                      Tool versions recorded at setup

# The columns of timings.tsv, in order. common.sh fills the run's identity and
# measurements, and prepare_command below fills the method's settings
#   tool_threads       Value passed to the tool's thread option, NA if it has none
#   min_mapq           MAPQ filter the tool applies (NA for timed Picard)
#   length_definition  How the tool measures a fragment's length
#   cpu_ids, physical_core_ids  The CPUs the Slurm step was bound to
#   wall_seconds ... max_rss_kib  GNU time's measurements of the native command
#   flush_seconds      Time spent waiting for the output to reach the disk
#   tool_exit_code     The native command's exit status
#   step_exit_code     The Slurm step's exit status, which is also nonzero for a
#                      wrong CPU placement or a missing output
timing_columns=(
    run_order repetition warmup method tool_version cores tool_threads
    min_mapq min_length max_length scope length_definition proper_pair
    node cpu_model cpu_ids physical_core_ids job_id step_id started_utc
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
    # Defaults shared by most methods. The branches below change only what differs
    #   tool_threads       Value passed to the tool's thread option, NA if it has none
    #   min_mapq           MAPQ filter the tool applies, NA if it has none
    #   length_definition  How the tool measures a fragment's length
    row=(
        [tool_threads]=NA
        [min_mapq]=$min_mapq
        [min_length]=$min_length
        [max_length]=$max_length
        [scope]=all
        [length_definition]=abs-TLEN
        [proper_pair]=not-required
    )
    stdout_file="$BENCH_WORK/tool-stdout.log"
    empty_output_allowed=false

    case "$method" in
        cfdnalab)
            command=("$cfdna" lengths --bam "$bam" --output-dir "$output"
                --length-bins "$min_length:$length_end:1"
                --indel-mode ignore --clip-mode aligned --min-mapq "$min_mapq"
                --n-threads "$cores" --chromosomes all)
            expected_output="$output/length_counts.tsv.zst"
            row[tool_version]=$cfdna_version
            row[tool_threads]=$cores
            row[length_definition]=directional-aligned-span
            ;;
        cfdnalab-proper)
            # Additional configuration: the same call, requiring properly paired reads
            command=("$cfdna" lengths --bam "$bam" --output-dir "$output"
                --length-bins "$min_length:$length_end:1"
                --indel-mode ignore --clip-mode aligned --min-mapq "$min_mapq"
                --n-threads "$cores" --chromosomes all --require-proper-pair)
            expected_output="$output/length_counts.tsv.zst"
            row[tool_version]=$cfdna_version
            row[tool_threads]=$cores
            row[length_definition]=directional-aligned-span
            row[proper_pair]=required
            ;;
        cfdnapro)
            # cfdnapro_lengths.R reads the BAM and calls cfDNAPro's callLength(). The
            # image contains R, cfDNAPro and the authors' reference packages
            # --cleanenv, --containall and --vanilla keep home-directory R settings out
            # The container sees $BENCH_WORK as /work, so /work/output is $output
            # Native coordinate curation measures directional fragment spans
            # cfDNAPro has no thread option here
            command=("$apptainer" exec --cleanenv --containall
                --bind "$BENCH_WORK:/work" --pwd /work
                "$cfdnapro_image" Rscript --vanilla /work/cfdnapro_lengths.R
                /work/sample.bam /work/output/lengths.tsv
                "$min_length" "$max_length" "$min_mapq" "$genome")
            expected_output="$output/lengths.tsv"
            row[tool_version]=$cfdnapro_version
            row[length_definition]=directional-aligned-span
            ;;
        finaletoolkit)
            # The packed environment's Python runs the FinaleToolkit launcher, so the
            # launcher does not look for its Python on the compute node's PATH
            # FinaleToolkit only counts properly paired reads. This command has no
            # thread option
            command=("$finale_python" "$finale" frag-length-bins "$bam"
                --bin-size 1 --min-length "$min_length" --max-length "$max_length"
                --min-mapq "$min_mapq" --output "$output/lengths.tsv")
            expected_output="$output/lengths.tsv"
            row[tool_version]=$finale_version
            row[proper_pair]=required
            ;;
        samtools)
            # samtools stats prints its full QC report, which is its native output
            # The input filter applies MAPQ while reading, without an extra BAM pass
            # -f 1 keeps paired reads and -F 3852 removes the excluded flags (see
            # experiment_commands.md). -m 1 keeps every length instead of trimming the
            # tail, and -i places the overflow bin above max_length. The comparator
            # later reads the inward-pair column of the IS section
            command=("$samtools" stats -@ "$cores"
                --input-fmt-option "filter=mapq >= $min_mapq"
                -f 1 -F 3852 -m 1 -i "$length_end" "$bam")
            stdout_file="$output/stats.txt"
            expected_output="$output/stats.txt"
            row[tool_version]=$samtools_version
            row[tool_threads]=$cores
            ;;
        picard)
            # Picard's documented standalone call, with its recommended 2 GB heap
            # Picard has no MAPQ option, so this timed run keeps low-MAPQ and QC-fail
            # pairs. An explicit histogram width keeps the tail. No PDF is requested
            # The filtered histogram used for concordance is made after all
            # measurements, by picard_filtered_histogram.sh
            command=("$picard" -Xmx2g CollectInsertSizeMetrics
                "I=$bam" "O=$output/metrics.txt" "HISTOGRAM_WIDTH=$max_length"
                "TMP_DIR=$BENCH_WORK/tmp")
            expected_output="$output/metrics.txt"
            row[tool_version]=$picard_version
            row[min_mapq]=NA
            ;;
        *)
            echo "Unknown length benchmark method: $method" >&2
            return 1
            ;;
    esac
}

# Copies a method's output into profiles/METHOD, for compare_lengths.py
#
# Only the first measured repetition is saved. Every repetition still computes its
# full output, so saving one copy reduces stored data without changing the timed work
#
# Parameters:
#   $1  Method name
save_profile() {
    local method=$1
    local destination="$BENCH_RESULTS/profiles/$method"
    mkdir -p "$destination"
    case "$method" in
        cfdnalab|cfdnalab-proper)
            cp "$output/length_counts.tsv.zst" "$destination/"
            ;;
        cfdnapro|finaletoolkit)
            cp "$output/lengths.tsv" "$destination/"
            ;;
        samtools)
            cp "$output/stats.txt" "$destination/"
            ;;
        picard)
            # The timed, unfiltered histogram. The filtered histogram for
            # concordance is added after all measurements
            cp "$output/metrics.txt" "$destination/"
            ;;
        *)
            echo "No saved profile defined for length method: $method" >&2
            return 1
            ;;
    esac
}
