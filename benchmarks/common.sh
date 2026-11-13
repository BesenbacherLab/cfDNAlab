#!/usr/bin/env bash

# Steps shared by the coverage, fragment length, end motif and midpoint benchmarks
#
# Each feature has a Slurm job script (for example lengths/benchmark_lengths.sh) that
# loads this file with "source". The job script reads from top to bottom:
#   1. Check the settings, before any slow work
#   2. Copy the BAM, references and tools to the compute node's local disk (untimed)
#   3. Run every method's warmups and measured repetitions. Each repetition is its own
#      srun step on the job's exclusive node, bound to BENCH_CORES physical cores.
#      The step runs measure_command.sh, which checks the CPUs and times the command
#   4. Compare the saved outputs (untimed)
# The native commands, and the files kept from them, are in <feature>/commands.sh
#
# This file only holds the mechanical parts shared by the four jobs. It does not
# choose commands, filters or which outputs are kept
#
# Variables shared between the job script, commands.sh and this file:
#   Set by the job script:
#     BENCH_*          Job settings from the submission environment
#     BENCH_WORK       The job's temporary directory on the compute node
#     BENCH_NODE       The compute node's name
#     output           The run's output directory, $BENCH_WORK/output
#     run_order        Counts the runs in this job. run_in_slurm_step increases it
#   Set when commands.sh is loaded:
#     timing_columns   The timings.tsv columns, in order
#   Set by prepare_command in commands.sh, before each run:
#     command               The native command, as an array of arguments
#     expected_output       The file or directory the command must create
#     empty_output_allowed  true when an empty output file is a valid result
#     stdout_file           Where the command's standard output goes
#     row[...]              The method's settings for its timings.tsv row
#   Set by run_in_slurm_step:
#     row[...]              The run's identity and measurements
#     step_exit_code        srun's exit status. 0 means the run succeeded
#
# The job scripts use set -euo pipefail. A failed tool is recorded and the job moves
# on to the next method. Failures in setup or in saving results stop the job
# immediately, because later results could otherwise be lost or mislabeled

# Holds the current run's timings.tsv row, as column name -> value. -g keeps it
# global even when this file is loaded from inside a function
declare -gA row=()

# ---------------------------------------------------------------------------
# Settings checks
# ---------------------------------------------------------------------------

# Checks the settings that all four jobs share, before any input is copied
#
# The job script sets its defaults first and checks its own methods and
# feature-specific settings afterwards
#
# Parameters:
#   $1  Feature name, used in error messages (coverage, lengths, ends or midpoints)
# Sets:
#   methods          The methods from BENCH_METHODS, as an array
#   BENCH_WARMUPS, BENCH_REPETITIONS
#                    Set to 0 and 1 for a pilot (BENCH_PILOT=1)
check_job_settings() {
    local feature=$1
    : "${SLURM_JOB_ID:?Submit with bash benchmarks/submit_benchmark.sh $feature}"
    : "${BENCH_BAM:?Set the input BAM path}"
    # Slurm sets TMPDIR to node-local scratch when the job starts
    : "${TMPDIR:?The compute node must provide TMPDIR for node-local scratch}"
    : "${BENCH_RESULTS:?Set a new results directory}"
    : "${BENCH_CORES:?Set a single core budget for this job: 1, 4, 12 or 32}"
    : "${BENCH_METHODS:?Set the methods for this job explicitly in BENCH_METHODS}"

    if [[ ! "$BENCH_CORES" =~ ^[1-9][0-9]*$ ]]; then
        echo 'Set BENCH_CORES to a single positive core count for this job.' >&2
        return 1
    fi
    if ((BENCH_CORES > ${SLURM_CPUS_PER_TASK:-0})); then
        echo "Core budget $BENCH_CORES exceeds allocated CPUs: ${SLURM_CPUS_PER_TASK:-0}" >&2
        return 1
    fi
    if [[ ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-0}} != 1 ]]; then
        echo "Expected one allocated node, but Slurm reports ${SLURM_JOB_NUM_NODES:-${SLURM_NNODES:-none}}." >&2
        return 1
    fi
    if [[ -e "$BENCH_RESULTS" ]]; then
        echo "BENCH_RESULTS already exists: $BENCH_RESULTS" >&2
        echo 'Choose a new path without creating it first. The benchmark creates the directory.' >&2
        return 1
    fi
    if [[ -n ${BENCH_EXPECTED_NODE:-} && $(hostname) != "$BENCH_EXPECTED_NODE" ]]; then
        echo "Expected node $BENCH_EXPECTED_NODE, received $(hostname)." >&2
        return 1
    fi

    export BENCH_PILOT=${BENCH_PILOT:-0}
    if [[ "$BENCH_PILOT" != 0 && "$BENCH_PILOT" != 1 ]]; then
        echo 'BENCH_PILOT must be 0 or 1.' >&2
        return 1
    fi
    if [[ "$BENCH_PILOT" == 1 ]]; then
        # A pilot runs each method once to estimate the time and memory it needs
        BENCH_WARMUPS=0
        BENCH_REPETITIONS=1
    fi

    local name
    for name in BENCH_WARMUPS BENCH_REPETITIONS BENCH_MIN_MAPQ BENCH_MIN_LENGTH BENCH_MAX_LENGTH; do
        if [[ ! "${!name}" =~ ^(0|[1-9][0-9]*)$ ]]; then
            echo "$name must be a nonnegative integer, not: ${!name}" >&2
            return 1
        fi
    done
    if ((BENCH_REPETITIONS < 1)); then
        echo 'Measure at least one repetition (BENCH_REPETITIONS).' >&2
        return 1
    fi
    if ((BENCH_MIN_MAPQ > 255 || BENCH_MIN_LENGTH < 10 || BENCH_MAX_LENGTH < BENCH_MIN_LENGTH || BENCH_MAX_LENGTH > 50000)); then
        echo 'Use a MAPQ within 0-255 and fragment lengths within 10-50000 bp.' >&2
        return 1
    fi

    read -r -a methods <<< "$BENCH_METHODS"
    if ((${#methods[@]} == 0)); then
        echo 'Set at least one method in BENCH_METHODS.' >&2
        return 1
    fi
    # A repeated name would double that method's repetitions. The spaces around each
    # name tell, for example, cfdnalab apart from cfdnalab-proper
    local method seen_methods=' '
    for method in "${methods[@]}"; do
        if [[ "$seen_methods" == *" $method "* ]]; then
            echo "Duplicate method in BENCH_METHODS: $method" >&2
            return 1
        fi
        seen_methods+="$method "
    done

    export BENCH_RUN_TIMEOUT_SECONDS=${BENCH_RUN_TIMEOUT_SECONDS:-0}
    if [[ ! "$BENCH_RUN_TIMEOUT_SECONDS" =~ ^(0|[1-9][0-9]*)$ ]]; then
        echo 'BENCH_RUN_TIMEOUT_SECONDS must be zero (disabled) or a positive integer.' >&2
        return 1
    fi
    if [[ "$BENCH_RUN_TIMEOUT_SECONDS" != 0 ]] && ! command -v timeout > /dev/null; then
        echo 'A run time limit requires GNU timeout on the compute node.' >&2
        return 1
    fi

    # The setup Python only runs untimed helpers that need its standard library
    export BENCH_SETUP_PYTHON=${BENCH_SETUP_PYTHON:-python3}
    if ! command -v "$BENCH_SETUP_PYTHON" > /dev/null; then
        echo 'Set BENCH_SETUP_PYTHON to an existing Python 3 executable for untimed setup.' >&2
        return 1
    fi
    # Two or more methods can be compared at the end of the job. Check the comparison
    # Python now, so a missing setting does not waste a full benchmark job
    if ((${#methods[@]} >= 2)); then
        : "${BENCH_COMPARE_PYTHON:?Set the comparison Python printed by install_tools.sh}"
        if ! command -v "$BENCH_COMPARE_PYTHON" > /dev/null; then
            echo 'BENCH_COMPARE_PYTHON must name an existing Python executable.' >&2
            return 1
        fi
    fi
}

# Checks that Picard is the version selected for this benchmark and prints it
#
# Picard prints a valid --version result but exits with status 1, normally an error
# status. That status is accepted only when the printed version is also correct
# Tabs and newlines become spaces, so the version fits in a timings.tsv cell
# install_tools.sh also uses this check
#
# Parameters:
#   $1  The Picard executable
benchmark_picard_version() {
    local executable=$1 version_output version_status=0
    version_output=$("$executable" CollectInsertSizeMetrics --version 2>&1) || version_status=$?
    if ((version_status > 1)) || [[ ! "$version_output" =~ (^|[[:space:]])Version:[[:space:]]*3\.5\.0($|[[:space:]]) ]]; then
        printf 'Picard version check failed (status %s): %s\n' "$version_status" "$version_output" >&2
        return 1
    fi
    printf '%s\n' "$version_output" | tr '\t\n' '  '
}

# ---------------------------------------------------------------------------
# Setup records and staging
# ---------------------------------------------------------------------------

# Saves a short, dated message in progress.log and prints it to Slurm's output
#
# Copying a BAM or unpacking a tool can take a long time. These messages tell that
# work apart from measured runs, and identify the last run that was started or saved
log_progress() {
    local message=$1
    printf '%s\t%s\n' "$(date -u +%FT%TZ)" "$message" | tee -a "$BENCH_RESULTS/progress.log"
}

# Creates the results directory with a header-only timings.tsv
#
# This happens before any input is copied, so a setup failure still leaves an
# inspectable record. A header-only timings.tsv means no run has finished yet
create_results_directory() {
    mkdir -p "$BENCH_RESULTS/profiles" "$BENCH_RESULTS/runs"
    (IFS=$'\t'; printf '%s\n' "${timing_columns[*]}") > "$BENCH_RESULTS/timings.tsv"
}

# Saves the input paths, run settings, node name and benchmark code version
#
# These records explain differences between jobs run at different times. This also
# copies measure_command.sh to the compute node. Every Slurm step runs that copy, so
# an edit to the repository during a long job cannot change its later runs
#
# Parameters:
#   $1  Feature name (coverage, lengths, ends or midpoints)
record_job_setup() {
    local feature=$1 commit uncommitted_changes
    cp "$SLURM_SUBMIT_DIR/benchmarks/measure_command.sh" "$BENCH_WORK/measure_command.sh"
    # Separate assignments, so a failing git stops the job instead of recording an
    # empty commit or a clean checkout
    commit=$(git -C "$SLURM_SUBMIT_DIR" rev-parse HEAD)
    uncommitted_changes=$(git -C "$SLURM_SUBMIT_DIR" status --porcelain)
    {
        printf 'node\t%s\npilot\t%s\n' "$BENCH_NODE" "$BENCH_PILOT"
        printf 'methods\t%s\ncores\t%s\n' "$BENCH_METHODS" "$BENCH_CORES"
        printf 'warmups\t%s\nmeasured_repetitions\t%s\n' "$BENCH_WARMUPS" "$BENCH_REPETITIONS"
        printf 'run_timeout_seconds\t%s\n' "$BENCH_RUN_TIMEOUT_SECONDS"
        printf 'benchmark_commit\t%s\n' "$commit"
        if [[ -n "$uncommitted_changes" ]]; then
            printf 'benchmark_dirty\ttrue\n'
        else
            printf 'benchmark_dirty\tfalse\n'
        fi
        printf 'feature\t%s\n' "$feature"
        printf 'input_bam\t%s\ninput_index\t%s\n' "$BENCH_BAM" "$BENCH_INDEX"
    } > "$BENCH_RESULTS/run.tsv"
    if [[ -n ${BENCH_INSTALL_MANIFEST:-} ]]; then
        cp "$BENCH_INSTALL_MANIFEST" "$BENCH_RESULTS/installation.tsv"
    fi
    if [[ -n ${BENCH_CFDNA_BUILD:-} ]]; then
        cp "$BENCH_CFDNA_BUILD" "$BENCH_RESULTS/cfdnalab-build.tsv"
    fi
    # Slurm's settings help interpret the peak memory it reports after the job
    # Timings and outputs can still be collected without them, so a failure here
    # only saves an error log and prints a warning
    if ! scontrol show config > "$BENCH_RESULTS/slurm-config.txt" 2> "$BENCH_RESULTS/slurm-config-error.log"; then
        echo 'Slurm configuration was unavailable. See slurm-config-error.log.' >&2
    fi
}

# Prints the autosome rows (name and length) of a chrom.sizes file
#
# Autosomes are the contigs named chr1..chr22 or 1..22, so both naming conventions
# work. The rows keep the file's names and order. A job that needs a comma-separated
# list uses: autosome_sizes FILE | cut -f 1 | paste -sd, -
#
# Parameters:
#   $1  A chrom.sizes file with contig name and length columns
autosome_sizes() {
    local chromosome_sizes=$1
    awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/' "$chromosome_sizes"
}

# Unpacks a tool's packed Conda environment onto the compute node's disk
#
# Programs in the archive can still refer to the original installation folder
# conda-unpack updates those paths to the folder created here. It runs with the
# setup Python, which needs no extra packages. The unpacked environment's own
# Python may depend on paths that are not yet updated
#
# Parameters:
#   $1  The packed archive (.tar.gz) made by install_tools.sh
#   $2  The folder to unpack into
unpack_environment() {
    local archive=$1 destination=$2
    log_progress "Unpacking environment archive: $archive"
    mkdir -p "$destination"
    tar -xzf "$archive" -C "$destination"
    "$BENCH_SETUP_PYTHON" "$destination/bin/conda-unpack"
}

# Copies a standalone program to $BENCH_WORK/bin and checks its shared libraries
#
# A copied program may still need libraries installed on the compute node. ldd checks
# them before any run, instead of every repetition failing to start. ldd reports an
# error for a program with its libraries built in, which is accepted
#
# Parameters:
#   $1  The program to copy
#   $2  The name of the local copy
copy_executable() {
    local executable=$1 name=$2 linkage
    log_progress "Copying executable: $name"
    cp -L "$executable" "$BENCH_WORK/bin/$name"
    if [[ ! -x "$BENCH_WORK/bin/$name" ]]; then
        echo "Not an executable: $executable" >&2
        return 1
    fi
    if ! linkage=$(ldd "$BENCH_WORK/bin/$name" 2>&1); then
        if [[ "$linkage" != *'not a dynamic executable'* && "$linkage" != *'statically linked'* ]]; then
            printf 'Could not inspect %s linkage: %s\n' "$name" "$linkage" >&2
            return 1
        fi
    fi
    if [[ "$linkage" == *'not found'* ]]; then
        printf 'Missing shared libraries for %s:\n%s\n' "$name" "$linkage" >&2
        return 1
    fi
}

# ---------------------------------------------------------------------------
# Warmups and measured repetitions
# ---------------------------------------------------------------------------

# Runs a repetition of the prepared command in its own Slurm step
#
# Call prepare_command from commands.sh first. srun starts measure_command.sh on the
# job's node, bound to BENCH_CORES physical cores. It receives the native command as
# arguments, so the timed command is exactly the one saved in runs/<run_order>/command.sh
# The step gets its own Slurm step ID, used to look up its peak memory with sacct
#
# Parameters:
#   $1  Method name
#   $2  Repetition number. Repetitions up to BENCH_WARMUPS are warmups
# Sets:
#   run_order, step_exit_code, and the run's identity and measurements in `row`
run_in_slurm_step() {
    local method=$1
    local repetition=$2
    local warmup=false
    if ((repetition <= BENCH_WARMUPS)); then warmup=true; fi
    run_order=$((run_order + 1))
    log_progress "Starting run $run_order: $method, repetition $repetition (warmup=$warmup)"

    # The command is saved before the run starts, so it is known even if the
    # allocation ends during the run. %q quotes each argument for a shell
    local run_directory="$BENCH_RESULTS/runs/$run_order"
    mkdir "$run_directory"
    printf '%q ' "${command[@]}" > "$run_directory/command.sh"
    printf '\n' >> "$run_directory/command.sh"

    # Remove the previous run's record and logs, so a failed step can never be
    # reported with values or messages from an earlier run. The output directory must
    # not exist yet, because clear_run_output deleted the previous run's
    rm -f "$BENCH_WORK/measurement.tsv" "$BENCH_WORK/tool-stdout.log" "$BENCH_WORK/tool-stderr.log" \
        "$BENCH_WORK/command-trace.log"
    mkdir "$output"

    # --exact and --cpus-per-task give the step exactly BENCH_CORES CPUs of the
    # allocation. --threads-per-core=1 makes those CPUs separate physical cores, and
    # --cpu-bind=threads keeps every process of the step on them
    step_exit_code=0
    srun --jobid="$SLURM_JOB_ID" --exact --nodes=1 --ntasks=1 \
        --cpus-per-task="$BENCH_CORES" --threads-per-core=1 --cpu-bind=threads \
        --export=ALL bash "$BENCH_WORK/measure_command.sh" \
        "$stdout_file" "$expected_output" "$empty_output_allowed" -- "${command[@]}" \
        || step_exit_code=$?

    row[run_order]=$run_order
    row[repetition]=$repetition
    row[warmup]=$warmup
    row[method]=$method
    row[cores]=$BENCH_CORES
    row[job_id]=$SLURM_JOB_ID
    row[pilot]=$BENCH_PILOT
    row[step_exit_code]=$step_exit_code

    # Copy the step's measurements into the row. A value the step did not reach,
    # for example the timing of a run with the wrong CPUs, is recorded as NA
    local -A measured=()
    local name value
    if [[ -f "$BENCH_WORK/measurement.tsv" ]]; then
        while IFS=$'\t' read -r name value; do
            measured[$name]=$value
        done < "$BENCH_WORK/measurement.tsv"
    fi
    for name in step_id node cpu_model cpu_ids physical_core_ids started_utc \
        wall_seconds user_seconds system_seconds max_rss_kib tool_exit_code \
        flush_seconds process_plus_flush_seconds output_bytes; do
        row[$name]=${measured[$name]:-NA}
    done
}

# Saves a finished run: its timings.tsv row, and its logs in runs/<run_order>/
#
# This runs after the Slurm step, outside the run's measured time and memory
# Failed runs are saved in the same way, so their row and messages stay available
save_run_record() {
    local run_directory="$BENCH_RESULTS/runs/$run_order"
    # A step that failed before running the command has no logs
    if [[ -f "$stdout_file" ]]; then cp "$stdout_file" "$run_directory/stdout.log"; fi
    if [[ -f "$BENCH_WORK/tool-stderr.log" ]]; then cp "$BENCH_WORK/tool-stderr.log" "$run_directory/stderr.log"; fi
    # Methods that run a helper script (freia, cfdnalab-gc, griffin-gc) also trace
    # each native call the helper made, with its variables filled in
    if [[ -f "$BENCH_WORK/command-trace.log" ]]; then cp "$BENCH_WORK/command-trace.log" "$run_directory/command-trace.log"; fi
    write_timing_row
    if ((step_exit_code != 0)); then
        if [[ -f "$run_directory/stderr.log" ]]; then cat "$run_directory/stderr.log" >&2; fi
        echo "Run $run_order (${row[method]}, repetition ${row[repetition]}) failed with step exit code $step_exit_code." >&2
        echo "Its command and logs are in $run_directory" >&2
    fi
}

# Appends `row` to timings.tsv, in the column order of timing_columns
#
# Every column needs a value and every value needs a column. A misspelled or
# forgotten setting therefore stops the job instead of leaving a wrong cell
write_timing_row() {
    local column values=()
    for column in "${timing_columns[@]}"; do
        # ${row[x]+present} is "present" only when row has the key x
        if [[ -z "${row[$column]+present}" ]]; then
            echo "No value for timing column '$column'. This is a bug in the benchmark scripts." >&2
            return 1
        fi
        values+=("${row[$column]}")
    done
    for column in "${!row[@]}"; do
        if [[ " ${timing_columns[*]} " != *" $column "* ]]; then
            echo "'$column' is not a timing column. This is a bug in the benchmark scripts." >&2
            return 1
        fi
    done
    (IFS=$'\t'; printf '%s\n' "${values[*]}") >> "$BENCH_RESULTS/timings.tsv"
}

# Deletes the run's output, after its row and any saved profile are in the results
#
# sync first lets pending writes on the scratch disk finish, so they cannot overlap
# the next measurement. The copied inputs and tools stay for the next run. The job's
# exit trap removes the whole scratch directory at the end
clear_run_output() {
    sync -f "$BENCH_WORK"
    rm -rf "$output"
    log_progress "Saved run $run_order: ${row[method]}, repetition ${row[repetition]} (step exit $step_exit_code)"
}

# ---------------------------------------------------------------------------
# Concordance
# ---------------------------------------------------------------------------

# Compares the saved profiles with the feature's comparison script
#
# This runs in the batch shell after all runs, outside every measured Slurm step
# It writes concordance.tsv in the results directory. A failed comparison keeps its
# error in concordance-error.log and stops the job. The timings and profiles are
# already saved by then
#
# Parameters:
#   $1   The comparison script
#   $2-  Its options, such as the job's length bounds
compare_saved_profiles() {
    local comparator=$1
    shift
    local saved_profile_count
    saved_profile_count=$(find "$BENCH_RESULTS/profiles" -mindepth 1 -maxdepth 1 -type d | wc -l)
    if ((saved_profile_count < 2)); then
        echo 'Concordance skipped: fewer than two methods saved a profile.'
        return 0
    fi
    : "${BENCH_COMPARE_PYTHON:?Set the comparison Python printed by install_tools.sh}"
    echo "Calculating concordance from $BENCH_RESULTS/profiles"
    # The comparison environment also supplies zstd, used to read cfDNAlab's tables
    # The table is first written under a temporary name, so a failed comparison
    # cannot replace an earlier valid concordance.tsv
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
