#!/usr/bin/env bash

# The native coverage commands, one per method
#
# benchmark_coverage.sh loads this file. Before every run it calls prepare_command,
# and after a method's first measured repetition it calls save_profile
#
# The commands use these variables, which the job script sets after copying the
# inputs and tools to the compute node. experiment_commands.md uses the same names
#   cores              CPU cores allocated to every run (BENCH_CORES)
#   bam                The copied input BAM
#   output             The run's empty output directory, $BENCH_WORK/output
#   regions            Shared BED4 of 500 kb windows over every BAM contig
#   min_mapq, min_length, max_length
#                      Shared filters. Length bounds are inclusive
#   autosomes          Comma-separated autosome names found in the BAM
#   cfdna, mosdepth, pandepth, samtools, finale_python, finale
#                      Tool paths
#   cfdna_version, mosdepth_version, pandepth_version, samtools_version, finale_version
#                      Tool versions recorded at setup
#
# The exclusion mask 3852 combines unmapped (4), mate unmapped (8), secondary (256),
# QC-fail (512), duplicate (1024) and supplementary (2048) flags. It does not
# require the proper-pair flag

# The columns of timings.tsv, in order. common.sh fills the run's identity and
# measurements, and prepare_command below fills the method's settings
#   run_order          Counts every attempted run in this job
#   repetition         Counts runs of the current method, including warmups
#   workload           The quantity the method measures (see coverage_quantity)
#   tool_threads       Value passed to the tool's thread option, NA if it has none
#   cpu_model, cpu_ids, physical_core_ids
#                      Identify the hardware, and the CPUs the Slurm step was bound to
#   wall_seconds ... max_rss_kib  GNU time's measurements of the native command
#   flush_seconds      Time spent waiting for the output to reach the disk
#   tool_exit_code     The native command's exit status
#   step_exit_code     The Slurm step's exit status, which is also nonzero for a
#                      wrong CPU placement or a missing output
timing_columns=(
    run_order repetition warmup method workload tool_version cores tool_threads
    node cpu_model cpu_ids physical_core_ids job_id step_id started_utc
    min_mapq min_length max_length scope proper_pair length_filter
    wall_seconds user_seconds system_seconds max_rss_kib tool_exit_code
    flush_seconds process_plus_flush_seconds output_bytes pilot step_exit_code
)

# Prints the quantity a coverage method measures
#
# Only methods that measure the same quantity are compared with each other. The
# value is also the workload column of timings.tsv
#
# Parameters:
#   $1  Method name
coverage_quantity() {
    local method=$1
    case "$method" in
        cfdnalab|cfdnalab-autosomes|mosdepth|samtools-depth|pandepth)
            echo per-base-read-depth ;;
        cfdnalab-fragment|cfdnalab-fragment-proper|mosdepth-fragment)
            echo per-base-fragment-coverage ;;
        cfdnalab-regions|mosdepth-regions|samtools-regions|pandepth-regions)
            echo regional-read-depth-500kb ;;
        cfdnalab-region-counts|cfdnalab-region-counts-proper|finaletoolkit-coverage)
            echo region-counts-500kb ;;
        *)
            echo "Unknown coverage method: $method" >&2
            return 1
            ;;
    esac
}

# Succeeds when two saved profiles measure the same quantity
#
# The job only runs compare_coverage.py when there is such a pair, and otherwise
# explains that concordance was skipped
saved_profiles_share_a_quantity() {
    local directory quantity saved_quantities=' '
    for directory in "$BENCH_RESULTS"/profiles/*/; do
        if [[ ! -d "$directory" ]]; then continue; fi
        quantity=$(coverage_quantity "$(basename "$directory")")
        if [[ "$saved_quantities" == *" $quantity "* ]]; then return 0; fi
        saved_quantities+="$quantity "
    done
    return 1
}

# Sets the native command for a method and its settings for timings.tsv
#
# Parameters:
#   $1  Method name
# Sets (see common.sh for how they are used):
#   command, expected_output, empty_output_allowed, stdout_file, row[...]
prepare_command() {
    local method=$1 workload
    workload=$(coverage_quantity "$method")
    # Defaults shared by most methods. The branches below change only what differs
    #   workload       The quantity the method measures (see coverage_quantity)
    #   tool_threads   Value passed to the tool's thread option, NA if it has none.
    #                  It means cfDNAlab threads, FinaleToolkit workers, PanDepth's -t,
    #                  or additional BAM decoding threads for mosdepth and samtools.
    #                  All of a tool's threads share the CPUs Slurm allocated
    #   length_filter  applied, or unsupported when the tool has no length filter
    row=(
        [workload]=$workload
        [tool_threads]=$cores
        [min_mapq]=$min_mapq
        [min_length]=$min_length
        [max_length]=$max_length
        [scope]=all
        [proper_pair]=not-required
        [length_filter]=applied
    )
    stdout_file="$BENCH_WORK/tool-stdout.log"
    empty_output_allowed=false

    case "$method" in
        # ---- Read depth at individual bases ----
        cfdnalab)
            # --ignore-gap counts the bases covered by the reads, not the gap between
            # mates. Explicit "all" matches mosdepth's whole-BAM scope
            command=("$cfdna" fcoverage --bam "$bam" --output-dir "$output"
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --chromosomes all --ignore-gap)
            expected_output="$output/fcoverage.per_position.bedgraph.zst"
            row[tool_version]=$cfdna_version
            ;;
        cfdnalab-autosomes)
            # Additional configuration: the same call on autosomes only
            command=("$cfdna" fcoverage --bam "$bam" --output-dir "$output"
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --chromosomes "$autosomes" --ignore-gap)
            expected_output="$output/fcoverage.per_position.bedgraph.zst"
            row[tool_version]=$cfdna_version
            row[scope]=autosomes
            ;;
        mosdepth)
            # --threads sets BAM decompression threads, while the coverage calculation
            # runs separately. The allocated core count is passed directly, as for every
            # native thread option. Ordinary mode excludes the gap between mates and
            # keeps normal CIGAR processing. Proper pairs are not required
            command=("$mosdepth" --threads "$cores"
                --include-flag 1 --flag 3852 --mapq "$min_mapq"
                --min-frag-len "$min_length" --max-frag-len "$max_length"
                "$output/coverage" "$bam")
            expected_output="$output/coverage.per-base.bed.gz"
            row[tool_version]=$mosdepth_version
            ;;
        samtools-depth)
            # -s removes mate overlap. -q 0 is the base-quality threshold and -Q is MAPQ
            # Without -J, deletions are not counted. The input filter applies inclusive
            # absolute TLEN bounds to both signs while reading the BAM. It does not
            # require both mates to pass or reconstruct their span
            # The native output is sparse: zero-depth positions are omitted, so an empty
            # file is valid and means no coverage. Positions are 1-based
            # https://www.htslib.org/doc/samtools-depth.html
            command=("$samtools" depth -@ "$cores"
                --input-fmt-option "filter=(tlen >= $min_length && tlen <= $max_length) || (tlen <= -$min_length && tlen >= -$max_length)"
                -s -q 0 -Q "$min_mapq" -G 3852 --require-flags 1
                -o "$output/samtools.depth.tsv" "$bam")
            expected_output="$output/samtools.depth.tsv"
            empty_output_allowed=true
            row[tool_version]=$samtools_version
            ;;
        pandepth)
            # -a writes every base, including zeros, to a compressed SiteDepth table
            # with 0-based positions. MAPQ and exclusion flags match. PanDepth has no
            # fragment length filter, paired-read inclusion flag or mate-overlap
            # correction, so overlapping mates count twice
            # https://github.com/HuiyangYu/PanDepth/blob/main/src/PanDepth.cpp
            command=("$pandepth" -i "$bam" -o "$output/pandepth"
                -q "$min_mapq" -x 3852 -t "$cores" -a)
            expected_output="$output/pandepth.SiteDepth.gz"
            row[tool_version]=$pandepth_version
            row[length_filter]=unsupported
            ;;

        # ---- Fragment coverage at individual bases ----
        cfdnalab-fragment)
            # Without --ignore-gap, each fragment covers its full span, including the
            # gap between mates
            command=("$cfdna" fcoverage --bam "$bam" --output-dir "$output"
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --chromosomes all)
            expected_output="$output/fcoverage.per_position.bedgraph.zst"
            row[tool_version]=$cfdna_version
            ;;
        cfdnalab-fragment-proper)
            # Additional configuration: fragment coverage from properly paired reads only
            command=("$cfdna" fcoverage --bam "$bam" --output-dir "$output"
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --chromosomes all --require-proper-pair)
            expected_output="$output/fcoverage.per_position.bedgraph.zst"
            row[tool_version]=$cfdna_version
            row[proper_pair]=required
            ;;
        mosdepth-fragment)
            # Fragment mode requires proper pairs and fills the complete TLEN span
            command=("$mosdepth" --threads "$cores"
                --include-flag 1 --flag 3852 --mapq "$min_mapq"
                --min-frag-len "$min_length" --max-frag-len "$max_length"
                --fragment-mode "$output/coverage" "$bam")
            expected_output="$output/coverage.per-base.bed.gz"
            row[tool_version]=$mosdepth_version
            row[proper_pair]=required
            ;;

        # ---- Regional read depth: covered bases summed per 500 kb window ----
        cfdnalab-regions)
            # The sum of covered bases per window, without length normalization
            command=("$cfdna" fcoverage --bam "$bam" --output-dir "$output"
                --min-mapq "$min_mapq"
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --n-threads "$cores" --chromosomes all --ignore-gap
                --by-size 500000 --per-window total)
            expected_output="$output/fcoverage.total.tsv.zst"
            row[tool_version]=$cfdna_version
            ;;
        mosdepth-regions)
            # mosdepth reports mean depth per window, which the comparator multiplies by
            # the window width. The documented MOSDEPTH_PRECISION=5 limits rounding in
            # that product. It changes output precision, not counting
            command=(env MOSDEPTH_PRECISION=5 "$mosdepth" --threads "$cores"
                --include-flag 1 --flag 3852 --mapq "$min_mapq"
                --min-frag-len "$min_length" --max-frag-len "$max_length"
                --by 500000 --no-per-base "$output/coverage" "$bam")
            expected_output="$output/coverage.regions.bed.gz"
            row[tool_version]=$mosdepth_version
            ;;
        samtools-regions)
            # bedcov sums read depth over the shared BED and prints the table. -j
            # excludes deletions and skipped bases. Unlike depth -s, bedcov counts
            # overlapping mates twice, and it has no thread option. samtools' input
            # filter selects paired reads within the TLEN bounds
            command=("$samtools" bedcov
                --input-fmt-option "filter=flag.paired && ((tlen >= $min_length && tlen <= $max_length) || (tlen <= -$min_length && tlen >= -$max_length))"
                -j -Q "$min_mapq" -G 3852
                "$regions" "$bam")
            stdout_file="$output/samtools.bedcov.tsv"
            expected_output="$output/samtools.bedcov.tsv"
            row[tool_version]=$samtools_version
            row[tool_threads]=NA
            ;;
        pandepth-regions)
            # -w makes PanDepth's own 500 kb windows. Its window table includes
            # TotalDepth, so no per-base table is written. Coordinates are 1-based and
            # inclusive, and are converted only when the comparator reads them
            command=("$pandepth" -i "$bam" -o "$output/pandepth"
                -q "$min_mapq" -x 3852 -t "$cores" -w 500000)
            expected_output="$output/pandepth.win.stat.gz"
            row[tool_version]=$pandepth_version
            row[length_filter]=unsupported
            ;;

        # ---- Regional fragment counts per 500 kb window ----
        cfdnalab-region-counts)
            # Length normalization gives each fragment a total weight of 1, so the
            # window sum estimates its fragment count. A fragment crossing a window
            # boundary splits its weight between the windows, unlike FinaleToolkit
            command=("$cfdna" fcoverage --bam "$bam" --output-dir "$output"
                --chromosomes all --by-size 500000 --per-window total --normalize-by-length
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --min-mapq "$min_mapq" --n-threads "$cores")
            expected_output="$output/length_normalized.fcoverage.total.tsv.zst"
            row[tool_version]=$cfdna_version
            ;;
        cfdnalab-region-counts-proper)
            # Additional configuration: the same counts from properly paired reads only
            command=("$cfdna" fcoverage --bam "$bam" --output-dir "$output"
                --chromosomes all --by-size 500000 --per-window total --normalize-by-length
                --min-fragment-length "$min_length" --max-fragment-length "$max_length"
                --min-mapq "$min_mapq" --n-threads "$cores" --require-proper-pair)
            expected_output="$output/length_normalized.fcoverage.total.tsv.zst"
            row[tool_version]=$cfdna_version
            row[proper_pair]=required
            ;;
        finaletoolkit-coverage)
            # Raw fragment counts in the shared BED windows. A boundary fragment
            # belongs to the window containing its midpoint. No normalization is
            # requested, and the native scale factor is 1. --threads sets worker
            # processes. FinaleToolkit requires proper pairs and uses read1 and TLEN
            # Its indexed query also needs read1 to overlap the window, which can
            # miss boundary fragments
            # https://epifluidlab.github.io/FinaleToolkit/documentation/cli_reference/index.html#coverage
            command=("$finale_python" "$finale" coverage "$bam" "$regions"
                --output "$output/finaletoolkit.coverage.bed.gz"
                --min-length "$min_length" --max-length "$max_length" --min-mapq "$min_mapq"
                --threads "$cores")
            expected_output="$output/finaletoolkit.coverage.bed.gz"
            row[tool_version]=$finale_version
            row[proper_pair]=required
            ;;
        *)
            echo "Unknown coverage method: $method" >&2
            return 1
            ;;
    esac
}

# Copies a method's output into profiles/METHOD, for compare_coverage.py
#
# Only the first measured repetition is saved. Every repetition still computes its
# full output, so saving one copy reduces stored data without changing the timed work
#
# Per-base profiles are saved with chrom.sizes, the contigs the method processed. A
# contig missing from a sparse per-base table then means zero coverage, not an
# unprocessed contig. Regional profiles are compared on the shared regions.bed
#
# Parameters:
#   $1  Method name
save_profile() {
    local method=$1
    local destination="$BENCH_RESULTS/profiles/$method"
    mkdir -p "$destination"
    case "$method" in
        cfdnalab|cfdnalab-fragment|cfdnalab-fragment-proper)
            cp "$output/fcoverage.per_position.bedgraph.zst" "$destination/"
            cp "$BENCH_WORK/chrom.sizes" "$destination/chrom.sizes"
            ;;
        cfdnalab-autosomes)
            # This run only processed the autosomes, so only they are listed
            cp "$output/fcoverage.per_position.bedgraph.zst" "$destination/"
            autosome_sizes "$BENCH_WORK/chrom.sizes" > "$destination/chrom.sizes"
            ;;
        mosdepth|mosdepth-fragment)
            cp "$output/coverage.per-base.bed.gz" "$destination/"
            cp "$BENCH_WORK/chrom.sizes" "$destination/chrom.sizes"
            ;;
        samtools-depth)
            # The uncompressed native table is compressed here, after timing
            gzip -c "$output/samtools.depth.tsv" > "$destination/positions.tsv.gz"
            cp "$BENCH_WORK/chrom.sizes" "$destination/chrom.sizes"
            ;;
        pandepth)
            # PanDepth's compressed native table is copied unchanged
            cp "$output/pandepth.SiteDepth.gz" "$destination/positions.tsv.gz"
            cp "$BENCH_WORK/chrom.sizes" "$destination/chrom.sizes"
            ;;
        cfdnalab-regions)
            cp "$output/fcoverage.total.tsv.zst" "$destination/"
            ;;
        mosdepth-regions)
            cp "$output/coverage.regions.bed.gz" "$destination/"
            ;;
        samtools-regions)
            cp "$output/samtools.bedcov.tsv" "$destination/"
            ;;
        pandepth-regions)
            cp "$output/pandepth.win.stat.gz" "$destination/"
            ;;
        cfdnalab-region-counts|cfdnalab-region-counts-proper)
            cp "$output/length_normalized.fcoverage.total.tsv.zst" "$destination/"
            ;;
        finaletoolkit-coverage)
            cp "$output/finaletoolkit.coverage.bed.gz" "$destination/"
            ;;
        *)
            echo "No saved profile defined for coverage method: $method" >&2
            return 1
            ;;
    esac
}
