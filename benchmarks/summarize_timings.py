#!/usr/bin/env python3
"""Summarize measured runs and optionally join completed Slurm step accounting.

Pass timing TSVs from the same experiment. Warmups, pilots and failed runs are
excluded from production summaries. --estimate uses successful pilot rows to
budget every observed method/core combination for the requested repetitions.
"""

import argparse
import csv
import math
from pathlib import Path
import re
import statistics
import sys


def read_accounting(path):
    """Index a single sacct export by JobIDRaw, preserving missing memory as NA."""
    records = {}
    if path is None:
        return records
    with path.open() as source:
        for row in csv.DictReader(source, delimiter="|"):
            identifier = row["JobIDRaw"]
            if identifier in records and records[identifier] != row:
                raise ValueError(f"Conflicting accounting rows for {identifier}")
            records[identifier] = row
    return records


def memory_kib(value):
    """Interpret sacct --units=K output, including explicit binary unit suffixes."""
    if not value or value == "NA":
        return "NA"
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([KMGT]?)", value.strip())
    if match is None:
        raise ValueError(f"Unrecognized sacct MaxRSS: {value}")
    return float(match[1]) * {"": 1, "K": 1, "M": 1024, "G": 1024**2, "T": 1024**3}[match[2]]


def read_runs(paths, accounting):
    """Read complete TSV rows and attach accounting to their exact job/step IDs."""
    runs = []
    identifiers = set()
    for path in paths:
        with path.open() as source:
            for row in csv.DictReader(source, delimiter="\t"):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"Malformed timing row in {path}")
                identifier = f"{row['job_id']}.{row['step_id']}"
                # Failed placement can have no step ID, so use its run order
                unique_identifier = (row["job_id"], row["run_order"])
                if unique_identifier in identifiers:
                    raise ValueError(f"Repeated run {unique_identifier}. Do not include a TSV twice")
                identifiers.add(unique_identifier)
                record = accounting.get(identifier, {})
                if record.get("NTasks", "") not in ("", "1"):
                    raise ValueError(f"Expected a single-task Slurm step: {identifier}")
                row["slurm_max_rss_kib"] = memory_kib(record.get("MaxRSS", ""))
                row["slurm_state"] = record.get("State") or "NA"
                row["slurm_exit_code"] = record.get("ExitCode") or "NA"
                row["source_tsv"] = str(path)
                runs.append(row)
    if not runs:
        raise ValueError("No timing rows found")
    return runs


def successful(row):
    """Require both native success and step success, including known accounting."""
    return (
        row["tool_exit_code"] == "0" and row["step_exit_code"] == "0"
        and row["slurm_state"] in ("NA", "COMPLETED")
        and row["slurm_exit_code"] in ("NA", "0:0")
    )


def distribution(values):
    """Return sample size, mean/SD, median and inclusive quartiles."""
    if not values:
        return [0, "NA", "NA", "NA", "NA", "NA"]
    lower, _, upper = statistics.quantiles(values, n=4, method="inclusive") if len(values) > 1 else [values[0]] * 3
    return [len(values), statistics.mean(values), statistics.stdev(values) if len(values) > 1 else "NA",
            statistics.median(values), lower, upper]


def summarize(runs, estimate=False):
    """Keep versions, nodes and command settings separate when pooling repetitions."""
    measurements = ("process_plus_flush_seconds", "wall_seconds", "flush_seconds",
                    "slurm_max_rss_kib", "max_rss_kib", "output_bytes", "effective_cores")
    # These fields are common identifiers or measurements, not command settings
    excluded = set(measurements) | {
        "run_order", "repetition", "warmup", "pilot", "tool_exit_code", "step_exit_code",
        "user_seconds", "system_seconds", "cpu_ids", "physical_core_ids", "job_id", "step_id",
        "started_utc", "slurm_state", "slurm_exit_code", "source_tsv",
    }
    setting_columns = sorted(set().union(*(row.keys() for row in runs)) - excluded)
    groups = {}
    failed = sum(not successful(row) for row in runs)
    if failed:
        print(f"Excluded {failed} failed runs. Their rows remain in the joined output.", file=sys.stderr)
    for row in runs:
        if not successful(row) or row["warmup"] != "false" or (row["pilot"] == "1") != estimate:
            continue
        key = tuple(row.get(column, "NA") for column in setting_columns)
        wall = float(row["wall_seconds"])
        row["effective_cores"] = (float(row["user_seconds"]) + float(row["system_seconds"])) / wall if wall > 0 else "NA"
        groups.setdefault(key, []).append(row)
    if not groups:
        raise ValueError("No successful pilot rows" if estimate else "No successful measured production rows")
    columns = setting_columns + [f"{metric}_{statistic}" for metric in measurements
                                 for statistic in ("n", "mean", "sd", "median", "q1", "q3")]
    output_rows = []
    for key, rows in sorted(groups.items()):
        result = dict(zip(setting_columns, key))
        for metric in measurements:
            values = [float(row[metric]) for row in rows if row[metric] != "NA"]
            if any(not math.isfinite(value) or value < 0 for value in values):
                raise ValueError(f"Invalid {metric} measurement")
            result.update(zip((f"{metric}_{statistic}" for statistic in ("n", "mean", "sd", "median", "q1", "q3")), distribution(values)))
        output_rows.append(result)
    return columns, output_rows


def main():
    """Write summary TSV to stdout, with optional joined rows and pilot budgeting."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("timings", nargs="+", type=Path)
    parser.add_argument("--accounting", type=Path, help="sacct pipe-delimited export in K units")
    parser.add_argument("--runs-output", type=Path, help="Write every original row with joined memory/status")
    parser.add_argument("--estimate", action="store_true", help="Summarize pilots and print the total production budget")
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--margin", type=float, default=0.25, help="Fraction added to pilot-based budget")
    arguments = parser.parse_args()
    if arguments.warmups < 0 or arguments.repetitions < 1 or not math.isfinite(arguments.margin) or arguments.margin < 0:
        parser.error("Use nonnegative warmups/margin and at least one measured repetition")
    runs = read_runs(arguments.timings, read_accounting(arguments.accounting))
    columns, rows = summarize(runs, arguments.estimate)
    writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    if arguments.runs_output:
        with arguments.runs_output.open("x") as destination:
            columns = list(dict.fromkeys(column for row in runs for column in row))
            writer = csv.DictWriter(destination, fieldnames=columns, delimiter="\t", lineterminator="\n", restval="NA")
            writer.writeheader()
            writer.writerows(runs)
    if arguments.estimate:
        attempted = {(row["method"], row["cores"]) for row in runs if row["pilot"] == "1"}
        completed = {(row["method"], row["cores"]) for row in runs
                     if row["pilot"] == "1" and row["warmup"] == "false" and successful(row)}
        if attempted != completed:
            raise ValueError("The pilot has failed configurations. Its total budget is incomplete")
        seconds = sum(row["process_plus_flush_seconds_mean"] for row in rows) * (arguments.warmups + arguments.repetitions)
        print(f"Observed matrix: {len(rows)} configurations, {seconds / 3600:.2f} h raw, "
              f"{seconds * (1 + arguments.margin) / 3600:.2f} h with margin. Add staging/cleanup. "
              "This excludes configurations absent from the supplied pilots.", file=sys.stderr)
    if any(row["slurm_max_rss_kib"] == "NA" for row in runs):
        print("Some Slurm memory is NA. Collect completed sacct records and rerun the join.", file=sys.stderr)


if __name__ == "__main__":
    main()
