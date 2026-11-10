#!/usr/bin/env bash

# Functions shared by the coverage, fragment length, motif and midpoint benchmarks
# Each feature script chooses the tool commands and biological filters. This file
# handles repeated runs, timing, saving results and cleaning up temporary output
# The scripts load these functions with "source", then call them as needed
# They use set -euo pipefail to stop on errors, missing variables or failed pipelines
# The installer also uses the Picard version check defined here
#
# What happens to the sample and its results:
# 1. The feature script copies the BAM, index, required references and programs to
#    the compute node's local disk. These copies finish before any tool is timed
# 2. benchmark_schedule lists a tool's warmups followed by its measured repetitions
#    A "method" can also mean a labeled configuration of a tool
# 3. The feature script uses srun to start each run with the requested CPU cores
#    Each run has its own Slurm step ID for looking up memory usage after the job
# 4. benchmark_finish_run adds the timing row to timings.tsv and keeps the first
#    measured output for numerical comparison. Other output is deleted before
#    the next repetition
# 5. After tool attempts finish, benchmark_compare_results writes concordance.tsv
#    from comparable saved profiles. This work is outside every tool's measurements
#
# srun starts a separate Bash program. Setting a variable there cannot update a
# variable in the main job script. Both can read BENCH_WORK, the job's temporary
# directory on the compute node, so that program leaves its results in these files:
#   time.tsv     Which run this was, its run time, CPU time, memory and exit code
#   io.tsv       Time spent finishing disk writes, total time and output file sizes
#   profile.path Path to the coverage, length, motif or midpoint output to keep
#   stdout.path  Path to the tool's printed output, useful when investigating failure
#   step.id      Slurm step ID, available even if the tool fails before timing starts
# The main script reads these files after srun finishes. Each attempt starts with
# fresh files so a failed run cannot accidentally reuse the preceding run's results
#
# If a tool fails, the script keeps its available measurements and error messages,
# skips its remaining repetitions and tries the next selected method. The job is
# reported as failed at the end. Setup or result-saving errors stop it immediately,
# because subsequent measurements could otherwise be lost or treated as successful

benchmark_picard_version() {
    # This function checks that Picard is the version selected for this benchmark
    # Picard prints a valid --version result but exits with code 1, normally an
    # error code. This is accepted only when the printed version is also correct
    # Tabs/newlines become spaces so the version fits in a timing-table cell
    local executable=$1 version_output version_status=0
    version_output=$("$executable" CollectInsertSizeMetrics --version 2>&1) || version_status=$?
    if ((version_status > 1)) || [[ ! "$version_output" =~ (^|[[:space:]])Version:[[:space:]]*3\.5\.0($|[[:space:]]) ]]; then
        printf 'Picard version check failed (status %s): %s\n' "$version_status" "$version_output" >&2
        return 1
    fi
    printf '%s\n' "$version_output" | tr '\t\n' '  '
}

benchmark_controls() {
    # This function checks cores, methods, pilot mode and the optional run time limit
    # It runs before the BAM is copied so an invalid setting fails without a long wait
    # Each feature script separately checks its tool names and biological filters
    export BENCH_PILOT=${BENCH_PILOT:-0}
    if [[ "$BENCH_PILOT" != 0 && "$BENCH_PILOT" != 1 ]]; then
        echo 'BENCH_PILOT must be 0 or 1.' >&2
        return 1
    fi
    if [[ ! ${BENCH_CORES:-} =~ ^[1-9][0-9]*$ ]]; then
        echo 'Set BENCH_CORES to a single positive core count for this job.' >&2
        return 1
    fi
    # Duplicate method names are rejected because they would double the repetitions
    # Spaces around each name distinguish, for example, cfdnalab from cfdnalab-proper
    local selected_method seen_methods=' ' selected_methods=()
    read -r -a selected_methods <<< "${BENCH_METHODS:-}"
    for selected_method in "${selected_methods[@]}"; do
        if [[ "$seen_methods" == *" $selected_method "* ]]; then
            echo "Duplicate method in BENCH_METHODS: $selected_method" >&2
            return 1
        fi
        seen_methods+="$selected_method "
    done
    export BENCH_RUN_TIMEOUT_SECONDS=${BENCH_RUN_TIMEOUT_SECONDS:-0}
    if [[ ! "$BENCH_RUN_TIMEOUT_SECONDS" =~ ^(0|[1-9][0-9]*)$ ]]; then
        echo 'BENCH_RUN_TIMEOUT_SECONDS must be zero (disabled) or a positive integer.' >&2
        return 1
    fi
    if [[ "$BENCH_RUN_TIMEOUT_SECONDS" != 0 ]]; then
        command -v timeout >/dev/null || {
            echo 'A run time limit requires GNU timeout on the compute node.' >&2
            return 1
        }
    fi
    if [[ -n ${BENCH_EXPECTED_NODE:-} && $(hostname) != "$BENCH_EXPECTED_NODE" ]]; then
        echo "Expected node $BENCH_EXPECTED_NODE, received $(hostname)." >&2
        return 1
    fi
    if [[ "$BENCH_PILOT" == 1 ]]; then
        # A pilot runs each selected method once to estimate time and memory needs
        # These counts take effect before the main script lists runs or records settings
        BENCH_WARMUPS=0
        BENCH_REPETITIONS=1
    fi
    export BENCH_SETUP_PYTHON=${BENCH_SETUP_PYTHON:-python3}
    command -v "$BENCH_SETUP_PYTHON" >/dev/null || {
        echo 'Set BENCH_SETUP_PYTHON to an existing Python 3 executable for untimed setup.' >&2
        return 1
    }
    # Multiple methods may produce a comparable pair. Check the comparison Python
    # before staging inputs so a missing setting does not waste a full benchmark job
    if ((${#selected_methods[@]} >= 2)); then
        : "${BENCH_COMPARE_PYTHON:?Set the comparison Python printed by install_tools.sh}"
        command -v "$BENCH_COMPARE_PYTHON" >/dev/null || {
            echo 'BENCH_COMPARE_PYTHON must name an existing Python executable.' >&2
            return 1
        }
    fi
}

benchmark_compare_results() {
    # This function compares saved outputs after every native tool run has finished
    # Arguments: the feature's comparison script, followed by its biological settings
    # It runs in the batch shell, outside the timed srun steps and their memory records
    # Call it on its own line so an error saving the comparison stops the job
    local comparator=$1
    shift
    local directory method quantity seen_quantities=' ' comparable=false
    for directory in "$BENCH_RESULTS"/profiles/*; do
        [[ -d "$directory" ]] || continue
        method=${directory##*/}
        quantity=$BENCH_FEATURE
        if [[ "$BENCH_FEATURE" == coverage ]]; then
            # Read depth, fragment coverage and regional counts need separate pairs
            case "$method" in
                *region-counts*|finaletoolkit-coverage) quantity=regional-counts ;;
                *fragment*) quantity=fragment-coverage ;;
                *) quantity=read-depth ;;
            esac
        fi
        if [[ "$seen_quantities" == *" $quantity "* ]]; then
            comparable=true
            break
        fi
        seen_quantities+="$quantity "
    done
    if [[ "$comparable" == false ]]; then
        echo 'Concordance skipped: no pair of saved profiles measures the same quantity.'
        return 0
    fi
    : "${BENCH_COMPARE_PYTHON:?Set the comparison Python printed by install_tools.sh}"
    echo "Calculating concordance from $BENCH_RESULTS/profiles"
    # The comparison environment also supplies zstd for reading compressed tracks
    # Write to a temporary name so a failed comparison cannot replace valid metrics
    if ! PATH="$(dirname -- "$BENCH_COMPARE_PYTHON"):$PATH" \
        "$BENCH_COMPARE_PYTHON" "$comparator" "$BENCH_RESULTS/profiles" "$@" \
        > "$BENCH_RESULTS/.concordance.tsv" 2> "$BENCH_RESULTS/concordance-error.log"; then
        cat "$BENCH_RESULTS/concordance-error.log" >&2
        echo "Concordance failed. See $BENCH_RESULTS/concordance-error.log" >&2
        rm -f "$BENCH_RESULTS/.concordance.tsv"
        return 1
    fi
    mv "$BENCH_RESULTS/.concordance.tsv" "$BENCH_RESULTS/concordance.tsv"
    echo "Concordance saved to $BENCH_RESULTS/concordance.tsv"
}

benchmark_unpack() {
    # This function unpacks a tool's Conda environment onto the compute node's disk
    # Programs in the archive can still refer to the original installation folder
    # conda-unpack updates those paths to the folder just created here
    # It runs with the setup Python, which needs no extra packages. The Python in
    # the unpacked environment may depend on paths that still need updating
    local archive=$1 destination=$2
    benchmark_progress "Unpacking environment archive: $archive"
    mkdir -p "$destination"
    tar -xzf "$archive" -C "$destination"
    "$BENCH_SETUP_PYTHON" "$destination/bin/conda-unpack"
}

benchmark_stage_binary() {
    # This function copies a standalone program and checks for the libraries it needs
    # Arguments: path to the program and the name to give its local copy
    # A copied executable may still need libraries installed on the compute node
    # ldd checks them before the script attempts repeated runs that cannot start
    # ldd reports an error for a program with its libraries built in. That case is
    # accepted. Other errors and output naming a missing library stop setup
    local executable=$1 name=$2 linkage
    benchmark_progress "Copying executable: $name"
    cp -L "$executable" "$BENCH_WORK/bin/$name"
    [[ -x "$BENCH_WORK/bin/$name" ]] || {
        echo "Not an executable: $executable" >&2
        return 1
    }
    linkage=$(ldd "$BENCH_WORK/bin/$name" 2>&1) || {
        if [[ "$linkage" != *'not a dynamic executable'* && "$linkage" != *'statically linked'* ]]; then
            printf 'Could not inspect %s linkage: %s\n' "$name" "$linkage" >&2
            return 1
        fi
    }
    if [[ "$linkage" == *'not found'* ]]; then
        printf 'Missing shared libraries for %s:\n%s\n' "$name" "$linkage" >&2
        return 1
    fi
}

benchmark_progress() {
    # Save a short, dated message in the results directory and in Slurm's output
    # Setup can spend a long time copying a BAM or unpacking a tool. These messages
    # distinguish that work from measured runs and identify the last completed run
    local message=$1
    printf '%s\t%s\n' "$(date -u +%FT%TZ)" "$message" | tee -a "$BENCH_RESULTS/progress.log"
}

benchmark_record_setup() {
    # This function saves input paths, run settings, node name and code version
    # beside the results before copying inputs or unpacking tool environments
    # These records help explain differences between jobs run at different times
    # It also copies this file so each srun program loads the same helper functions
    # SLURM_SUBMIT_DIR locates it because Slurm keeps its copy of the submitted
    # job script elsewhere, away from the other repository files
    cp "$SLURM_SUBMIT_DIR/benchmarks/common.sh" "$BENCH_WORK/common.sh"
    {
        printf 'node\t%s\npilot\t%s\n' "$BENCH_NODE" "$BENCH_PILOT"
        printf 'methods\t%s\ncores\t%s\n' "$BENCH_METHODS" "$BENCH_CORES"
        printf 'warmups\t%s\nmeasured_repetitions\t%s\n' "$BENCH_WARMUPS" "$BENCH_REPETITIONS"
        printf 'run_timeout_seconds\t%s\n' "$BENCH_RUN_TIMEOUT_SECONDS"
        printf 'benchmark_commit\t%s\n' "$(git -C "$SLURM_SUBMIT_DIR" rev-parse HEAD)"
        if [[ -n $(git -C "$SLURM_SUBMIT_DIR" status --porcelain) ]]; then
            printf 'benchmark_dirty\ttrue\n'
        else
            printf 'benchmark_dirty\tfalse\n'
        fi
        printf 'feature\t%s\n' "$BENCH_FEATURE"
        printf 'input_bam\t%s\ninput_index\t%s\n' "$BENCH_BAM" "$BENCH_INDEX"
    } > "$BENCH_RESULTS/run.tsv"
    if [[ -n ${BENCH_INSTALL_MANIFEST:-} ]]; then
        cp "$BENCH_INSTALL_MANIFEST" "$BENCH_RESULTS/installation.tsv"
    fi
    if [[ -n ${BENCH_CFDNA_BUILD:-} ]]; then
        cp "$BENCH_CFDNA_BUILD" "$BENCH_RESULTS/cfdnalab-build.tsv"
    fi
    # Slurm's settings help interpret its reported peak memory after the job
    # If those settings are unavailable, the script saves an error log and warns
    # Tool timings and outputs can still be collected, so this failure is not fatal
    scontrol show config > "$BENCH_RESULTS/slurm-config.txt" 2> "$BENCH_RESULTS/slurm-config-error.log" || {
        echo 'Slurm configuration was unavailable. See slurm-config-error.log.' >&2
    }
    # The feature script waits for staged files to finish writing after setup
    # That wait happens before warmups and does not evict the input from memory
}

benchmark_schedule() {
    # This function lists each method's warmups immediately before its measurements
    # It reads methods_to_run, BENCH_CORES, BENCH_WARMUPS and BENCH_REPETITIONS from
    # the main script. Each line gives method, cores, repetition and warmup flag
    # With 1 warmup + 3 measurements, repetition 1 is marked true and 2-4 are false
    # The summary excludes warmup rows. Griffin instead gets its separate short demo,
    # so its runner sets full-input warmups to zero before asking for this list
    local method repetition is_warmup
    for method in "${methods_to_run[@]}"; do
        for ((repetition=1; repetition<=BENCH_WARMUPS+BENCH_REPETITIONS; repetition++)); do
            is_warmup=false
            if ((repetition <= BENCH_WARMUPS)); then is_warmup=true; fi
            printf '%s\t%s\t%s\t%s\n' "$method" "$BENCH_CORES" "$repetition" "$is_warmup"
        done
    done
}

benchmark_time_command() {
    # This function times the tool, then measures the wait for its disk writes
    # The first argument is the start of its timing row: method, settings and run ID
    # Next come its output folder, printed-output file and error file, followed by
    # the complete tool command. The command is stored as an array so a file path
    # with spaces stays a single argument
    # printf adds the identifying columns after timing. GNU time's format string
    # could misread '%' or backslashes in a name as formatting instructions
    #
    # The main script catches errors from this function so it can record failed runs
    # In that situation Bash does not automatically stop this function on an error,
    # even with set -e. Each required file operation therefore has "|| return 1" to
    # prevent a failed read/write from being treated as a successful measurement
    local prefix=$1 output_directory=$2 stdout_path=$3 stderr_path=$4
    shift 4
    local command_status=0 measurements flush_seconds output_bytes wall_seconds combined_seconds
    local timed_command=("$@")
    if [[ ${BENCH_RUN_TIMEOUT_SECONDS:-0} != 0 ]]; then
        # GNU time records the command even if the time limit stops it
        # timeout asks it to stop at the limit and forces it after another 60 seconds
        # Available measurements are retained, with the run counted as failed
        timed_command=(timeout --verbose --kill-after=60 "$BENCH_RUN_TIMEOUT_SECONDS" "${timed_command[@]}")
    fi
    # Old measurements are removed in case the timing program itself fails to start
    rm -f "$BENCH_WORK/measurements.tsv" "$BENCH_WORK/flush.tsv" || return 1
    # samtools stats prints its results instead of writing a named output file
    # This file records where that text goes so it can be kept if the run fails
    printf '%s\n' "$stdout_path" > "$BENCH_WORK/stdout.path" || return 1
    # GNU time records elapsed seconds, CPU seconds spent by the tool and operating
    # system, and maximum RAM use (RSS, in KiB). The exit code is recorded separately,
    # including termination by a signal or time limit
    /usr/bin/time --quiet --format $'%e\t%U\t%S\t%M' \
        --output "$BENCH_WORK/measurements.tsv" "${timed_command[@]}" \
        > "$stdout_path" 2> "$stderr_path" || command_status=$?
    measurements=$(cat "$BENCH_WORK/measurements.tsv") || return 1
    # Measurements are saved before a failure is returned, so unsuccessful attempts
    # remain visible in timings.tsv
    printf '%s\t%s\t%s\n' "$prefix" "$measurements" "$command_status" > "$BENCH_WORK/time.tsv" || return 1
    if ((command_status != 0)); then return "$command_status"; fi
    # A tool may finish while Linux is still holding output in memory to write later
    # This command waits for those writes and records the wait separately from run time
    # sync -f also waits for other pending writes on the same filesystem, so this
    # number can include more than the time spent writing this tool's output alone
    /usr/bin/time --quiet --format '%e' --output "$BENCH_WORK/flush.tsv" \
        sync -f "$output_directory" || return 1
    flush_seconds=$(cat "$BENCH_WORK/flush.tsv") || return 1
    # Output size includes intermediate BAMs/tables required by the tool's workflow
    output_bytes=$(du -sb "$output_directory" | cut -f 1) || return 1
    wall_seconds=$(cut -f 1 "$BENCH_WORK/measurements.tsv") || return 1
    combined_seconds=$(awk -v wall="$wall_seconds" -v flush="$flush_seconds" 'BEGIN {printf "%.3f", wall+flush}') || return 1
    # Disk-write time, combined time and output size are saved together once all are known
    # If the tool or disk-write check failed, the main script uses NA for these fields
    printf '%s\t%s\t%s\n' "$flush_seconds" "$combined_seconds" "$output_bytes" > "$BENCH_WORK/io.tsv" || return 1
}

benchmark_finish_run() {
    # This function adds an attempt to timings.tsv, keeps useful output and clears space
    # It runs after srun finishes so copying a large coverage track is excluded from
    # the tool's measured time and its Slurm step's memory usage
    # Arguments: method name, repetition number, warmup flag and srun exit code
    # The job's BENCH_* settings and run_order locate and label its results
    #
    # A failed tool is recorded here and the main script can continue with other tools
    # Call-site requirement: this function must be called on its own line, because
    # saving errors must stop the job. An if/|| condition would prevent Bash from
    # automatically stopping on errors inside the function
    local method=$1 repetition=$2 is_warmup=$3 step_exit_code=$4
    local io_values=$'NA\tNA\tNA' profile stdout_path timing_values
    if [[ -s "$BENCH_WORK/io.tsv" ]]; then io_values=$(cat "$BENCH_WORK/io.tsv"); fi
    # A failed tool may have no disk-write measurements, which is represented by NA
    # Its timing row must exist. Reading it separately makes a missing file stop the
    # job instead of printf succeeding with an empty row
    timing_values=$(cat "$BENCH_WORK/time.tsv")
    printf '%s\t%s\t%s\t%s\n' "$timing_values" "$io_values" "$BENCH_PILOT" "$step_exit_code" >> "$BENCH_RESULTS/timings.tsv"
    if ((step_exit_code != 0)); then
        # The tool's messages are saved before temporary files are deleted, along
        # with printed result tables such as samtools stats that may explain failure
        mkdir -p "$BENCH_RESULTS/failed-run-$run_order"
        for profile in "$BENCH_WORK/tool.log" "$BENCH_WORK/tool-error.log"; do
            if [[ -f "$profile" ]]; then cp "$profile" "$BENCH_RESULTS/failed-run-$run_order/"; fi
        done
        if [[ -s "$BENCH_WORK/stdout.path" ]]; then
            stdout_path=$(cat "$BENCH_WORK/stdout.path")
            if [[ -f "$stdout_path" ]]; then cp "$stdout_path" "$BENCH_RESULTS/failed-run-$run_order/tool-stdout.log"; fi
        fi
        echo "Run $run_order failed. Completed rows and diagnostics are in $BENCH_RESULTS" >&2
    fi
    # Output from the first measured repetition is kept for numerical comparisons
    # Every repetition still recomputes its full output, so saving a single copy
    # reduces stored data without reducing the work included in the timing
    if [[ "$is_warmup" == false ]] && ((step_exit_code == 0 && repetition == BENCH_WARMUPS+1)); then
        profile=$(cat "$BENCH_WORK/profile.path")
        mkdir -p "$BENCH_RESULTS/profiles/$method"
        case "$method" in
            freia)
                # The saved name distinguishes intact left-end counts from FrEIA's
                # old, incorrect both-end abundance tables
                cp "$profile" "$BENCH_RESULTS/profiles/$method/left-motifs.pq"
                ;;
            samtools-depth|pandepth)
                # Whole-genome per-base tables are very large, so only the selected
                # chromosomes are kept for comparison. The timed run still produced
                # the full table. The saved chrom.sizes file records which chromosomes
                # and lengths are present in this subset
                if [[ "$profile" == *.gz ]]; then
                    gzip -dc "$profile" |
                        awk 'NR==FNR {selected[$1]=1; next} $1 in selected' "$BENCH_WORK/capture.chrom.sizes" - |
                        gzip > "$BENCH_RESULTS/profiles/$method/positions.tsv.gz"
                else
                    awk 'NR==FNR {selected[$1]=1; next} $1 in selected' "$BENCH_WORK/capture.chrom.sizes" "$profile" |
                        gzip > "$BENCH_RESULTS/profiles/$method/positions.tsv.gz"
                fi
                cp "$BENCH_WORK/capture.chrom.sizes" "$BENCH_RESULTS/profiles/$method/chrom.sizes"
                ;;
            *)
                cp -R "$profile" "$BENCH_RESULTS/profiles/$method/"
                if [[ "$method" == cfdnapro && "$BENCH_FEATURE" == ends ]]; then
                    # The motif counts are saved with a record of the reference used
                    cp "$BENCH_WORK/output/reference.tsv" "$BENCH_RESULTS/profiles/$method/"
                fi
                if [[ -f "$BENCH_WORK/chrom.sizes" ]]; then
                    # An absent chromosome in a sparse coverage file could mean zero
                    # coverage or an unprocessed chromosome. The saved chromosome list
                    # and lengths distinguish those cases. An autosome-only run gets
                    # a restricted list matching its actual selection
                    if [[ "$method" == *autosomes* ]]; then
                        awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/' "$BENCH_WORK/chrom.sizes" > "$BENCH_RESULTS/profiles/$method/chrom.sizes"
                    else
                        cp "$BENCH_WORK/chrom.sizes" "$BENCH_RESULTS/profiles/$method/chrom.sizes"
                    fi
                fi
                ;;
        esac
    fi
    # Pending writes on the temporary disk finish before the next run starts
    # This also covers saved results if they use that disk. The timer has stopped
    sync -f "$BENCH_WORK"
    # This run's output and temporary messages are removed to make space for the next
    # run. The copied inputs/programs remain for reuse. The main script removes the whole
    # temporary directory when the job ends, including when it stops after an error
    rm -rf "$BENCH_WORK/output" "$BENCH_WORK/coverage" "$BENCH_WORK/tool.log" "$BENCH_WORK/tool-error.log" "$BENCH_WORK/stdout.path"
    benchmark_progress "Saved run $run_order: $method, repetition $repetition (step exit $step_exit_code)"
}
