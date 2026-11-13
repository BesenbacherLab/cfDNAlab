#!/usr/bin/env bash

# Times a native command inside a Slurm step and records what was measured
#
# Every warmup and measured repetition runs this script in its own srun step.
# Slurm binds that step to BENCH_CORES physical cores on the job's exclusive node,
# so the native command cannot use more CPUs than its budget. The script:
#   1. Checks that the step really got the job's node and BENCH_CORES physical cores
#   2. Runs the native command under GNU time (elapsed time, CPU time, peak memory)
#   3. Times the wait for the command's output to be written to disk
#   4. Measures the size of the output directory
#   5. Checks that the expected output exists
#
# Each value is appended to $BENCH_WORK/measurement.tsv as a "name<TAB>value" line
# as soon as it is known. After the step ends, the job script copies these values
# into the run's timings.tsv row. A value is missing if this script stopped before
# measuring it, and the job script then records NA
#
# Usage:
#   measure_command.sh STDOUT_FILE EXPECTED_OUTPUT EMPTY_OUTPUT_ALLOWED -- COMMAND [ARGUMENTS...]
#     STDOUT_FILE           Receives the command's standard output. This is a log for most
#                           tools, and the native result for tools that print their results
#     EXPECTED_OUTPUT       File or directory that the command must create
#     EMPTY_OUTPUT_ALLOWED  true when an empty output file is a valid result, otherwise false
#     COMMAND               The native command and its arguments, exactly as they are timed
# The command's standard error goes to $BENCH_WORK/tool-stderr.log. Its output
# directory is $BENCH_WORK/output, which the job script creates empty for every run
#
# Job settings read from the environment: BENCH_WORK, BENCH_CORES, BENCH_NODE and
# BENCH_RUN_TIMEOUT_SECONDS. srun sets SLURM_STEP_ID
#
# Exit status: the native command's exit status when it fails, 1 for a wrong CPU
# placement or a missing output, and 0 when everything succeeded

set -euo pipefail

if (($# < 5)) || [[ "$4" != -- ]]; then
    echo 'Usage: measure_command.sh STDOUT_FILE EXPECTED_OUTPUT EMPTY_OUTPUT_ALLOWED -- COMMAND [ARGUMENTS...]' >&2
    exit 2
fi
stdout_file=$1
expected_output=$2
empty_output_allowed=$3
shift 4
native_command=("$@")
if [[ "$empty_output_allowed" != true && "$empty_output_allowed" != false ]]; then
    echo "EMPTY_OUTPUT_ALLOWED must be true or false, not: $empty_output_allowed" >&2
    exit 2
fi

record_file="$BENCH_WORK/measurement.tsv"
output_directory="$BENCH_WORK/output"
stderr_file="$BENCH_WORK/tool-stderr.log"

# Appends a value to the measurement record. Tabs and newlines become spaces, so
# every value stays a single cell in timings.tsv
record() {
    local name=$1
    local value=$2
    local tab_and_newline=$'\t\n'
    printf '%s\t%s\n' "$name" "${value//[$tab_and_newline]/ }" >> "$record_file"
}

# ---------------------------------------------------------------------------
# 1. CPU placement
# ---------------------------------------------------------------------------

# The step ID comes first, so even a step that fails immediately can be found in
# Slurm's accounting
record step_id "$SLURM_STEP_ID"
node=$(hostname)
record node "$node"
record cpu_model "$(awk -F ': ' '/^model name/ {print $2; exit}' /proc/cpuinfo)"
cpu_ids=$(awk '/Cpus_allowed_list/ {print $2}' /proc/self/status)
record cpu_ids "$cpu_ids"

# Hyperthreading can give a physical core two CPU IDs. To count physical cores,
# look up the core behind each CPU ID this step may use. Core numbers repeat across
# CPU sockets, so a core is identified by its socket and core number together
# The allowed list looks like "0-3,8". A single ID such as 8 is read as the range 8-8
physical_cores=()
for cpu_range in ${cpu_ids//,/ }; do
    first_cpu=${cpu_range%-*}
    last_cpu=${cpu_range#*-}
    for ((cpu_id = first_cpu; cpu_id <= last_cpu; cpu_id++)); do
        topology="/sys/devices/system/cpu/cpu${cpu_id}/topology"
        socket_id=$(cat "$topology/physical_package_id")
        core_id=$(cat "$topology/core_id")
        physical_cores+=("$socket_id:$core_id")
    done
done
physical_core_ids=$(printf '%s\n' "${physical_cores[@]}" | sort -u | paste -sd, -)
physical_core_count=$(printf '%s\n' "${physical_cores[@]}" | sort -u | wc -l)
record physical_core_ids "$physical_core_ids"

# Slurm's binding depends on the cluster configuration, so check it here. The step
# must run on the job's node and see exactly BENCH_CORES CPUs, which must be
# BENCH_CORES different physical cores. nproc counts the CPUs this process may use
# OpenMP variables are removed for that call, because nproc would report them instead
logical_cpu_count=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)
if [[ "$node" != "$BENCH_NODE" ]] || ((logical_cpu_count != BENCH_CORES || physical_core_count != BENCH_CORES)); then
    echo "Unexpected CPU placement: node $node, CPUs $cpu_ids, physical cores $physical_core_ids." >&2
    echo "Expected $BENCH_CORES physical cores on $BENCH_NODE. The command was not run." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 2. Timed native command
# ---------------------------------------------------------------------------

timed_command=("${native_command[@]}")
if [[ "$BENCH_RUN_TIMEOUT_SECONDS" != 0 ]]; then
    # timeout asks the command to stop at the limit and kills it 60 seconds later
    # It runs inside GNU time, so a stopped run still gets its measurements
    timed_command=(timeout --verbose --kill-after=60 "$BENCH_RUN_TIMEOUT_SECONDS" "${timed_command[@]}")
fi

# GNU time records elapsed seconds, user and system CPU seconds, and peak memory
# (maximum resident set size, in KiB). It writes them to its own file, so the
# command's standard output and error stay unchanged. --quiet keeps GNU time from
# adding a line to that file when the command fails
gnu_time_file="$BENCH_WORK/gnu-time.tsv"
rm -f "$gnu_time_file"
record started_utc "$(date -u +%FT%TZ)"
tool_exit_code=0
/usr/bin/time --quiet --format $'%e\t%U\t%S\t%M' --output "$gnu_time_file" \
    "${timed_command[@]}" > "$stdout_file" 2> "$stderr_file" || tool_exit_code=$?
IFS=$'\t' read -r wall_seconds user_seconds system_seconds max_rss_kib < "$gnu_time_file"
record wall_seconds "$wall_seconds"
record user_seconds "$user_seconds"
record system_seconds "$system_seconds"
record max_rss_kib "$max_rss_kib"
# The exit status includes termination by a signal (128 + signal number) and an
# expired time limit (124 from timeout)
record tool_exit_code "$tool_exit_code"
if ((tool_exit_code != 0)); then
    echo "The native command failed with exit status $tool_exit_code." >&2
    exit "$tool_exit_code"
fi

# ---------------------------------------------------------------------------
# 3. Wait for the output to reach the disk
# ---------------------------------------------------------------------------

# A tool can finish while Linux still holds some of its output in memory, to write
# later. sync -f waits for those writes on the output's filesystem. Its duration is
# recorded separately from the command's elapsed time. It can include other pending
# writes on the same filesystem, not only this tool's output
flush_time_file="$BENCH_WORK/flush-time.txt"
/usr/bin/time --quiet --format '%e' --output "$flush_time_file" sync -f "$output_directory"
flush_seconds=$(cat "$flush_time_file")
record flush_seconds "$flush_seconds"
record process_plus_flush_seconds "$(awk -v wall="$wall_seconds" -v flush="$flush_seconds" 'BEGIN {printf "%.3f", wall + flush}')"

# ---------------------------------------------------------------------------
# 4. Output size
# ---------------------------------------------------------------------------

# The size includes intermediate files that a tool's workflow writes, such as
# FrEIA's filtered BAM
record output_bytes "$(du -sb "$output_directory" | cut -f 1)"

# ---------------------------------------------------------------------------
# 5. Expected output
# ---------------------------------------------------------------------------

# A command that exits with status 0 but leaves no output, or an empty output file,
# has not produced a usable result. Some outputs can validly be empty, for example
# a sparse depth table when nothing is covered. The caller says so explicitly
if [[ ! -e "$expected_output" ]]; then
    echo "The command succeeded but its expected output is missing: $expected_output" >&2
    exit 1
fi
if [[ -f "$expected_output" && ! -s "$expected_output" && "$empty_output_allowed" == false ]]; then
    echo "The command succeeded but its expected output is empty: $expected_output" >&2
    exit 1
fi
