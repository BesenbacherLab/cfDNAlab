#!/usr/bin/env python3
"""Run GNU time's command while sampling the tool's complete process tree.

Usage: python measure_memory.py MEMORY_TSV /usr/bin/time [options] COMMAND ...
Write sampled peak RSS, time-weighted mean RSS, sample count and nominal sampling
interval to a temporary TSV row. The Bash runner appends it to its timing row.
Use Linux ps and the standard library, with no additional Python dependencies.
"""

from pathlib import Path
import subprocess
import sys
import time


def process_tree_rss(time_process_id):
    """Sum resident KiB for GNU time's descendants in a single ps snapshot.

    Include worker processes and piped commands. Exclude the sampler and GNU
    time itself. RSS counts shared resident pages in each process that maps them,
    so this is a sum of process RSS rather than unique physical memory usage.
    """
    snapshot = subprocess.check_output(
        ["ps", "-eo", "pid=,ppid=,rss="], universal_newlines=True
    )
    processes = [tuple(map(int, line.split())) for line in snapshot.splitlines() if line.strip()]
    descendants = {time_process_id}
    while True:
        children = {process_id for process_id, parent_id, rss in processes if parent_id in descendants}
        if children.issubset(descendants):
            break
        descendants.update(children)
    return sum(rss for process_id, parent_id, rss in processes if process_id in descendants and process_id != time_process_id)


def main():
    """Sample once per second and preserve the native command's exit status.

    Weight each sample by the elapsed time until the next snapshot or command
    exit. Waiting on the child with a timeout stops promptly when it finishes,
    rather than adding a full sampling interval after every benchmark call.
    """
    if len(sys.argv) < 3:
        raise SystemExit("Usage: measure_memory.py MEMORY_TSV /usr/bin/time [options] COMMAND ...")
    memory_path = Path(sys.argv[1])
    interval_seconds = 1.0
    process = subprocess.Popen(sys.argv[2:])
    started = time.monotonic()
    previous_time = started
    previous_rss = 0
    peak_rss = 0
    rss_integral = 0.0
    sample_count = 0
    while True:
        rss = process_tree_rss(process.pid)
        sampled_at = time.monotonic()
        rss_integral += previous_rss * (sampled_at - previous_time)
        previous_time, previous_rss = sampled_at, rss
        peak_rss = max(peak_rss, rss)
        sample_count += 1
        try:
            exit_code = process.wait(timeout=interval_seconds)
            break
        except subprocess.TimeoutExpired:
            continue
    finished = time.monotonic()
    rss_integral += previous_rss * (finished - previous_time)
    mean_rss = rss_integral / (finished - started)
    memory_path.write_text("{}\t{:.3f}\t{}\t{}\n".format(peak_rss, mean_rss, sample_count, interval_seconds))
    raise SystemExit(exit_code if exit_code >= 0 else 128 - exit_code)


if __name__ == "__main__":
    main()
