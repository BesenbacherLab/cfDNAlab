#!/usr/bin/env bash

# The native end-motif commands, one per method
#
# benchmark_end_motifs.sh loads this file. Before every run it calls prepare_command,
# and after a method's first measured repetition it calls save_profile. Each method
# counts k-base motifs at fragment ends, inside the fragment
#
# The commands use these variables, which the job script sets after copying the
# inputs and tools to the compute node. experiment_commands.md uses the same names
#   cores              CPU cores allocated to every run (BENCH_CORES)
#   bam                The copied input BAM
#   output             The run's empty output directory, $BENCH_WORK/output
#   k                  Motif length (BENCH_K)
#   min_mapq, min_length, max_length
#                      Shared filters. Length bounds are inclusive
#   autosomes          Comma-separated autosome names found in the BAM
#   genome             cfDNAPro's genome label (BENCH_GENOME)
#   reference_2bit     The copied 2bit reference
#   cfdna, finale_python, finale, apptainer, cfdnapro_image
#                      Tool paths
#   cfdna_version, finale_version, cfdnapro_version, freia_version
#                      Tool versions recorded at setup

# The columns of timings.tsv, in order. common.sh fills the run's identity and
# measurements, and prepare_command below fills the method's settings. The
# requested length bounds and how a method applies them are separate columns, which
# tells native filters apart from cfDNAPro's subsetting in R
#   tool_threads       Value passed to the tool's thread or worker option, NA if it has none
#   source             Whether motif bases come from the reference or the read
#   clip               How soft-clipped ends are handled
#   effective_*_length The length bounds actually passed to the tool
#   cpu_ids, physical_core_ids  The CPUs the Slurm step was bound to
#   wall_seconds ... max_rss_kib  GNU time's measurements of the native command
#   flush_seconds      Time spent waiting for the output to reach the disk
#   tool_exit_code     The native command's exit status
#   step_exit_code     The Slurm step's exit status, which is also nonzero for a
#                      wrong CPU placement or a missing output
timing_columns=(
    run_order repetition warmup method tool_version cores tool_threads k min_mapq
    requested_min_length requested_max_length length_filter scope source clip proper_pair
    node cpu_model cpu_ids physical_core_ids job_id step_id started_utc
    effective_min_length effective_max_length
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
    #   tool_threads          Value passed to the tool's thread or worker option, NA if none
    #   requested_*_length    The job's length bounds
    #   length_filter         How the bounds are applied: natively, or in R for cfDNAPro
    #   source                Whether motif bases come from the reference or the read
    #   clip                  How soft-clipped ends are handled
    #   effective_*_length    The bounds actually passed to the tool
    row=(
        [tool_threads]=$cores
        [k]=$k
        [min_mapq]=$min_mapq
        [requested_min_length]=$min_length
        [requested_max_length]=$max_length
        [length_filter]=applied
        [scope]=all
        [source]=reference
        [clip]=aligned
        [proper_pair]=not-required
        [effective_min_length]=$min_length
        [effective_max_length]=$max_length
    )
    stdout_file="$BENCH_WORK/tool-stdout.log"
    empty_output_allowed=false

    case "$method" in
        cfdnalab-reference)
            # Reference bases at the aligned fragment ends, on all BAM contigs
            command=("$cfdna" ends --bam "$bam" --output-dir "$output"
                --k-inside "$k" --k-outside 0 --all-motifs
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --source-inside reference --clip-strategy aligned
                --ref-2bit "$reference_2bit" --chromosomes all)
            expected_output="$output/end_motifs.zarr"
            row[tool_version]=$cfdna_version
            ;;
        cfdnalab-reference-proper)
            # The reference call, requiring properly paired reads
            command=("$cfdna" ends --bam "$bam" --output-dir "$output"
                --k-inside "$k" --k-outside 0 --all-motifs
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --source-inside reference --clip-strategy aligned
                --ref-2bit "$reference_2bit" --chromosomes all --require-proper-pair)
            expected_output="$output/end_motifs.zarr"
            row[tool_version]=$cfdna_version
            row[proper_pair]=required
            ;;
        cfdnalab-reference-proper-wide)
            # Additional configuration: the proper-pair reference call with the widest
            # length bounds cfDNAlab accepts, 10-50000 bp
            command=("$cfdna" ends --bam "$bam" --output-dir "$output"
                --k-inside "$k" --k-outside 0 --all-motifs
                --min-mapq "$min_mapq"
                --min-fragment-length 10 --max-fragment-length 50000
                --n-threads "$cores" --source-inside reference --clip-strategy aligned
                --ref-2bit "$reference_2bit" --chromosomes all --require-proper-pair)
            expected_output="$output/end_motifs.zarr"
            row[tool_version]=$cfdna_version
            row[proper_pair]=required
            row[effective_min_length]=10
            row[effective_max_length]=50000
            ;;
        cfdnalab-reference-autosomes)
            # Additional configuration: the reference call on autosomes only
            command=("$cfdna" ends --bam "$bam" --output-dir "$output"
                --k-inside "$k" --k-outside 0 --all-motifs
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --source-inside reference --clip-strategy aligned
                --ref-2bit "$reference_2bit" --chromosomes "$autosomes")
            expected_output="$output/end_motifs.zarr"
            row[tool_version]=$cfdna_version
            row[scope]=autosomes
            ;;
        cfdnalab-read-aligned)
            # Read bases at the aligned fragment ends. Compared with FrEIA, which
            # also takes motif bases from the reads
            command=("$cfdna" ends --bam "$bam" --output-dir "$output"
                --k-inside "$k" --k-outside 0 --all-motifs
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --source-inside read --clip-strategy aligned
                --chromosomes all)
            expected_output="$output/end_motifs.zarr"
            row[tool_version]=$cfdna_version
            row[source]=read
            ;;
        cfdnalab-read)
            # Additional configuration: read bases, skipping soft-clipped ends
            command=("$cfdna" ends --bam "$bam" --output-dir "$output"
                --k-inside "$k" --k-outside 0 --all-motifs
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --source-inside read --clip-strategy skip
                --chromosomes all)
            expected_output="$output/end_motifs.zarr"
            row[tool_version]=$cfdna_version
            row[source]=read
            row[clip]=skip
            ;;
        cfdnalab-read-include)
            # Additional configuration: read bases, keeping clipped bases at the
            # aligned boundary
            command=("$cfdna" ends --bam "$bam" --output-dir "$output"
                --k-inside "$k" --k-outside 0 --all-motifs
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --source-inside read --clip-strategy include-at-aligned-boundary
                --chromosomes all)
            expected_output="$output/end_motifs.zarr"
            row[tool_version]=$cfdna_version
            row[source]=read
            row[clip]=include-at-aligned-boundary
            ;;
        finaletoolkit)
            # FinaleToolkit reconstructs fragment coordinates from read1 and TLEN and
            # requires properly paired reads. The pinned upstream commit applies
            # inclusive length bounds and assigns each fragment to the window
            # containing its midpoint. --strand both counts both fragment ends
            # The packed environment's Python runs the launcher (see the job script)
            command=("$finale_python" "$finale" end-motifs "$bam" "$reference_2bit" -k "$k"
                --min-length "$min_length" --max-length "$max_length"
                --min-mapq "$min_mapq" --threads "$cores" --strand both
                --output "$output/motifs.tsv")
            expected_output="$output/motifs.tsv"
            row[tool_version]=$finale_version
            row[proper_pair]=required
            row[clip]=unrestricted
            ;;
        cfdnapro)
            # cfdnapro_end_motifs.R calls cfDNAPro's callMotif(type="s"), which counts
            # reference-derived left ends after native coordinate curation. cfDNAPro
            # has no motif length filter, so the script subsets fragment lengths in R
            # It keeps the BAM contigs present in cfDNAPro's reference (all but chrEBV
            # for this BAM). cfDNAPro has no thread option for this calculation
            # The container sees $BENCH_WORK as /work, so /work/output is $output
            command=("$apptainer" exec --cleanenv --containall
                --bind "$BENCH_WORK:/work" --pwd /work
                "$cfdnapro_image" Rscript --vanilla /work/cfdnapro_end_motifs.R
                /work/sample.bam /work/output
                "$min_length" "$max_length" "$min_mapq" "$genome" "$k")
            expected_output="$output/cfdnapro-motifs.tsv"
            row[tool_version]=$cfdnapro_version
            row[tool_threads]=NA
            row[scope]=bam-contigs-in-cfdnapro-reference
            row[length_filter]=applied-in-r
            row[clip]=unrestricted
            ;;
        freia)
            # FrEIA's filtering, end extraction and motif counting, in freia_motifs.sh
            # Its intact output is the left-end table. The broken native both-end
            # aggregate is not requested
            command=(bash "$BENCH_WORK/freia_motifs.sh" "$cores")
            expected_output="$output/4_FrEIA/3_Abundances/benchmark/Genome_Lvl/T__sample.pq"
            row[tool_version]=$freia_version
            row[source]=read
            ;;
        *)
            echo "Unknown end-motif method: $method" >&2
            return 1
            ;;
    esac
}

# Copies a method's output into profiles/METHOD, for compare_end_motifs.py
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
        cfdnalab-*)
            cp -R "$output/end_motifs.zarr" "$destination/"
            ;;
        finaletoolkit)
            cp "$output/motifs.tsv" "$destination/"
            ;;
        cfdnapro)
            # The motif counts are saved with the reference cfDNAPro used and the
            # list of BAM contigs it kept
            cp "$output/cfdnapro-motifs.tsv" "$output/reference.tsv" "$output/contig_scope.tsv" "$destination/"
            ;;
        freia)
            # The new name says these are FrEIA's intact left-end counts, not its
            # incorrect both-end abundance table
            cp "$output/4_FrEIA/3_Abundances/benchmark/Genome_Lvl/T__sample.pq" "$destination/left-motifs.pq"
            ;;
        *)
            echo "No saved profile defined for end-motif method: $method" >&2
            return 1
            ;;
    esac
}
