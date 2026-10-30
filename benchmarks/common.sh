#!/usr/bin/env bash

# Shared setup and measurement details for the four benchmark scripts
# Source this file. Native tool commands and their settings stay in each runner

benchmark_picard_version() {
    # Picard deliberately returns 1 for --version. Validate its output as well
    local executable=$1 version_output version_status=0
    version_output=$("$executable" CollectInsertSizeMetrics --version 2>&1) || version_status=$?
    if ((version_status > 1)) || [[ ! "$version_output" =~ (^|[[:space:]])Version:[[:space:]]*3\.5\.0($|[[:space:]]) ]]; then
        printf 'Picard version check failed (status %s): %s\n' "$version_status" "$version_output" >&2
        return 1
    fi
    printf '%s\n' "$version_output" | tr '\t\n' '  '
}

benchmark_controls() {
    # A pilot keeps this job's selected methods and core count, changing repetitions only
    export BENCH_PILOT=${BENCH_PILOT:-0}
    if [[ "$BENCH_PILOT" != 0 && "$BENCH_PILOT" != 1 ]]; then
        echo 'BENCH_PILOT must be 0 or 1.' >&2
        return 1
    fi
    if [[ ! ${BENCH_CORES:-} =~ ^[1-9][0-9]*$ ]]; then
        echo 'Set BENCH_CORES to a single positive core count for this job.' >&2
        return 1
    fi
    # A method is a single repetition block, not an accidentally repeated block
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
        BENCH_WARMUPS=0
        BENCH_REPETITIONS=1
    fi
    export BENCH_SETUP_PYTHON=${BENCH_SETUP_PYTHON:-python3}
    command -v "$BENCH_SETUP_PYTHON" >/dev/null || {
        echo 'Set BENCH_SETUP_PYTHON to an existing Python 3 executable for untimed setup.' >&2
        return 1
    }
}

benchmark_unpack() {
    # Relocation uses only the standard library, including for non-Python tools
    local archive=$1 destination=$2
    mkdir -p "$destination"
    tar -xzf "$archive" -C "$destination"
    "$BENCH_SETUP_PYTHON" "$destination/bin/conda-unpack"
}

benchmark_stage_binary() {
    # Copy a selected native executable and detect unavailable shared libraries
    local executable=$1 name=$2 linkage
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

benchmark_record_setup() {
    # Keep verbose provenance once per job, outside every native command timer
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
    scontrol show config > "$BENCH_RESULTS/slurm-config.txt" 2> "$BENCH_RESULTS/slurm-config-error.log" || {
        echo 'Slurm configuration was unavailable. See slurm-config-error.log.' >&2
    }
    # Clear staging write-back before measuring the first native command
    sync -f "$BENCH_WORK"
}

benchmark_schedule() {
    # Run every selected method at this job's core count, with warmups first
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
    # Time the native process, then filesystem flush. Identifiers are plain TSV
    # Keep them out of GNU time's format, which interprets '%' and backslashes
    local prefix=$1 output_directory=$2 stdout_path=$3 stderr_path=$4
    shift 4
    local command_status=0 measurements flush_seconds output_bytes wall_seconds combined_seconds
    local timed_command=("$@")
    if [[ ${BENCH_RUN_TIMEOUT_SECONDS:-0} != 0 ]]; then
        timed_command=(timeout --verbose --kill-after=60 "$BENCH_RUN_TIMEOUT_SECONDS" "${timed_command[@]}")
    fi
    rm -f "$BENCH_WORK/measurements.tsv" "$BENCH_WORK/flush.tsv" || return 1
    printf '%s\n' "$stdout_path" > "$BENCH_WORK/stdout.path" || return 1
    /usr/bin/time --quiet --format $'%e\t%U\t%S\t%M\t%x' \
        --output "$BENCH_WORK/measurements.tsv" "${timed_command[@]}" \
        > "$stdout_path" 2> "$stderr_path" || command_status=$?
    measurements=$(cat "$BENCH_WORK/measurements.tsv") || return 1
    # Keep the shell-observed status, including signal exits and timeout statuses
    measurements=${measurements%$'\t'*}$'\t'"$command_status"
    printf '%s\t%s\n' "$prefix" "$measurements" > "$BENCH_WORK/time.tsv" || return 1
    if ((command_status != 0)); then return "$command_status"; fi
    /usr/bin/time --quiet --format '%e' --output "$BENCH_WORK/flush.tsv" \
        sync -f "$output_directory" || return 1
    flush_seconds=$(cat "$BENCH_WORK/flush.tsv") || return 1
    output_bytes=$(du -sb "$output_directory" | cut -f 1) || return 1
    wall_seconds=$(cut -f 1 "$BENCH_WORK/measurements.tsv") || return 1
    combined_seconds=$(awk -v wall="$wall_seconds" -v flush="$flush_seconds" 'BEGIN {printf "%.3f", wall+flush}') || return 1
    printf '%s\t%s\t%s\n' "$flush_seconds" "$combined_seconds" "$output_bytes" > "$BENCH_WORK/io.tsv" || return 1
}

benchmark_finish_run() {
    # Append outside srun and preserve failed diagnostics before cleanup
    # A recorded tool failure returns normally so the parent can try other methods
    # Results-bookkeeping failures still abort the parent under set -e
    local method=$1 cores=$2 repetition=$3 is_warmup=$4 step_exit_code=$5
    local io_values=$'NA\tNA\tNA' profile stdout_path timing_values
    if [[ -s "$BENCH_WORK/io.tsv" ]]; then io_values=$(cat "$BENCH_WORK/io.tsv"); fi
    timing_values=$(cat "$BENCH_WORK/time.tsv")
    printf '%s\t%s\t%s\t%s\n' "$timing_values" "$io_values" "$BENCH_PILOT" "$step_exit_code" >> "$BENCH_RESULTS/timings.tsv"
    if ((step_exit_code != 0)); then
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
    if [[ "$method" != freia && "$is_warmup" == false ]] && ((step_exit_code == 0 && repetition == BENCH_WARMUPS+1)); then
        profile=$(cat "$BENCH_WORK/profile.path")
        mkdir -p "$BENCH_RESULTS/profiles/$method"
        case "$method" in
            samtools-depth|pandepth)
                # Preserve only the declared small concordance scope for huge tables
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
                if [[ -f "$BENCH_WORK/chrom.sizes" ]]; then
                    if [[ "$method" == *autosomes* ]]; then
                        awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2])$/' "$BENCH_WORK/chrom.sizes" > "$BENCH_RESULTS/profiles/$method/chrom.sizes"
                    else
                        cp "$BENCH_WORK/chrom.sizes" "$BENCH_RESULTS/profiles/$method/chrom.sizes"
                    fi
                fi
                ;;
        esac
    fi
    sync -f "$BENCH_WORK"
    rm -rf "$BENCH_WORK/output" "$BENCH_WORK/coverage" "$BENCH_WORK/tool.log" "$BENCH_WORK/tool-error.log" "$BENCH_WORK/stdout.path"
}
