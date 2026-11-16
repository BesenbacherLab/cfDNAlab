#!/usr/bin/env python3
"""Summarize measured runs and optionally join completed Slurm step accounting.

Pass timing TSVs from the same experiment. By default, summaries include only
successful measured runs, excluding warmups and shared-node debugging runs
(pilot=1). --estimate instead uses successful debugging runs to estimate the time
needed for each observed method/core combination and the requested repetitions.
"""

import argparse
import csv
from itertools import groupby
import math
from pathlib import Path
import re
import statistics
import sys


# Summarize each measurement independently. A missing memory value must not
# discard the same run's valid time measurements or become a zero-byte result
MEASUREMENT_COLUMNS = (
    "process_plus_flush_seconds",
    "wall_seconds",
    "flush_seconds",
    "slurm_max_rss_kib",
    "max_rss_kib",
    "output_directory_bytes",
    "effective_cores",
)
STATISTIC_NAMES = ("n", "mean", "sd", "median", "q1", "q3")

# These columns describe individual runs or their measurements. Everything left
# over is a configuration setting, including each feature's filters, tool version
# and node. Keeping those settings in the grouping prevents unlike runs being pooled
NON_SETTING_COLUMNS = set(MEASUREMENT_COLUMNS) | {
    "run_order",
    "repetition",
    "warmup",
    "pilot",
    "tool_exit_code",
    "step_exit_code",
    "user_seconds",
    "system_seconds",
    "cpu_ids",
    "physical_core_ids",
    "job_id",
    "step_id",
    "started_utc",
    "slurm_state",
    "slurm_exit_code",
    "source_tsv",
}


def read_accounting(path):
    """Read Slurm's records so each run can look up its exact job and step.

    Accounting is optional. Repeated identical records are accepted, but conflicting
    records for the same step stop the summary rather than choosing one silently.
    """
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
    """Convert sacct's MaxRSS value to KiB, keeping missing memory as NA.

    With sacct --units=K, a number without a suffix is already in KiB. Explicit
    K/M/G/T suffixes use powers of 1024 and are converted to the same unit.
    """
    if not value or value == "NA":
        return "NA"
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)([KMGT]?)", value.strip())
    if match is None:
        raise ValueError(f"Unrecognized sacct MaxRSS: {value}")
    number, unit = match.groups()
    unit_in_kib = {"": 1, "K": 1, "M": 1024, "G": 1024**2, "T": 1024**3}
    return float(number) * unit_in_kib[unit]


def read_runs(paths, accounting):
    """Read all timing rows and attach the matching Slurm measurements.

    Keep warmups, debugging runs and failures here so the optional detailed table
    retains the full record. Selection for averages happens later. A missing Slurm
    record stays NA, because accounting can arrive after the timing files.
    """
    runs = []
    identifiers = set()
    for path in paths:
        with path.open() as source:
            for row in csv.DictReader(source, delimiter="\t"):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"Malformed timing row in {path}")
                # Match a tool's step, not the containing job's aggregate record
                step_identifier = f"{row['job_id']}.{row['step_id']}"
                # Failed placement can have no step ID, so use its run order to
                # reject duplicate input rows without merging separate failures
                unique_identifier = (row["job_id"], row["run_order"])
                if unique_identifier in identifiers:
                    raise ValueError(
                        f"Repeated run {unique_identifier}. Do not include a TSV twice"
                    )
                identifiers.add(unique_identifier)
                record = accounting.get(step_identifier, {})
                if record.get("NTasks", "") not in ("", "1"):
                    raise ValueError(f"Expected a single-task Slurm step: {step_identifier}")
                row["slurm_max_rss_kib"] = memory_kib(record.get("MaxRSS", ""))
                row["slurm_state"] = record.get("State") or "NA"
                row["slurm_exit_code"] = record.get("ExitCode") or "NA"
                row["source_tsv"] = str(path)
                runs.append(row)
    if not runs:
        raise ValueError("No timing rows found")
    return runs


def successful(row):
    """Require success from the command, its Slurm step and available accounting.

    Missing accounting does not invalidate an otherwise successful run. A known
    scheduler failure, such as an out-of-memory kill, always excludes it.
    """
    return (
        row["tool_exit_code"] == "0"
        and row["step_exit_code"] == "0"
        and row["slurm_state"] in ("NA", "COMPLETED")
        and row["slurm_exit_code"] in ("NA", "0:0")
    )


def distribution(values):
    """Describe one measurement using named statistics.

    A singleton has no sample standard deviation. For an empty measurement, keep
    n=0 and all statistics as NA. Inclusive quartiles interpolate within the
    observed values, matching the documented summary definition.
    """
    if not values:
        return {
            "n": 0,
            "mean": "NA",
            "sd": "NA",
            "median": "NA",
            "q1": "NA",
            "q3": "NA",
        }

    if len(values) == 1:
        lower_quartile = upper_quartile = values[0]
        standard_deviation = "NA"
    else:
        lower_quartile, _, upper_quartile = statistics.quantiles(values, n=4, method="inclusive")
        standard_deviation = statistics.stdev(values)
    return {
        "n": len(values),
        "mean": statistics.mean(values),
        "sd": standard_deviation,
        "median": statistics.median(values),
        "q1": lower_quartile,
        "q3": upper_quartile,
    }


def group_runs_by_settings(runs, setting_columns):
    """Group runs with identical settings, keeping the groups in a stable order.

    A setting absent from a feature's TSV is treated as NA. Each group contains
    the setting values and an iterator over the matching runs.
    """
    def setting_values(run):
        return tuple(run.get(column, "NA") for column in setting_columns)

    # itertools.groupby only combines adjacent matches, so sort by the same key
    ordered_runs = sorted(runs, key=setting_values)
    return groupby(ordered_runs, key=setting_values)


def select_measured_runs(runs, estimate):
    """Select successful measurements and calculate each run's average busy cores.

    Normal summaries exclude debugging runs. Estimate mode selects them instead.
    Add the derived core measurement to each selected row so it is also available
    in the detailed output, while retaining every original row in runs.
    """
    selected_runs = []
    for run in runs:
        if not successful(run) or run["warmup"] != "false":
            continue
        is_debugging_run = run["pilot"] == "1"
        if is_debugging_run != estimate:
            continue

        # Calculate this ratio per run, before taking averages across repetitions
        wall_seconds = float(run["wall_seconds"])
        if wall_seconds > 0:
            cpu_seconds = float(run["user_seconds"]) + float(run["system_seconds"])
            run["effective_cores"] = cpu_seconds / wall_seconds
        else:
            run["effective_cores"] = "NA"
        selected_runs.append(run)
    return selected_runs


def summarize(runs, estimate=False):
    """Select measurements, group matching configurations and describe each group.

    Return the setting columns, complete output columns and summary rows. Tool
    versions, nodes and every feature-specific setting remain separate groups.
    Missing values are skipped separately for each measurement, giving it its own n.
    """
    all_columns = {column for run in runs for column in run}
    setting_columns = sorted(all_columns - NON_SETTING_COLUMNS)
    failed_runs = sum(not successful(run) for run in runs)
    if failed_runs:
        print(
            f"Excluded {failed_runs} failed runs. Their rows remain in the joined output.",
            file=sys.stderr,
        )

    selected_runs = select_measured_runs(runs, estimate)
    if not selected_runs and not estimate:
        raise ValueError(
            "No successful measured runs remain after excluding "
            "warmups and debugging runs (pilot=1)"
        )

    summary_columns = setting_columns + [
        f"{measurement}_{statistic}"
        for measurement in MEASUREMENT_COLUMNS
        for statistic in STATISTIC_NAMES
    ]
    summaries = []
    for settings, matching_runs in group_runs_by_settings(selected_runs, setting_columns):
        # Reuse the same repetitions for each measurement without exhausting the iterator
        repetitions = list(matching_runs)
        summary = dict(zip(setting_columns, settings))
        for measurement in MEASUREMENT_COLUMNS:
            values = [
                float(run[measurement])
                for run in repetitions
                if run[measurement] != "NA"
            ]
            if any(not math.isfinite(value) or value < 0 for value in values):
                raise ValueError(f"Invalid {measurement} measurement")
            for statistic, value in distribution(values).items():
                summary[f"{measurement}_{statistic}"] = value
        summaries.append(summary)
    return setting_columns, summary_columns, summaries


def pilot_budgets(runs, setting_columns, summaries, warmups, repetitions, margin):
    """Estimate time for each observed debugging configuration, including failures.

    Successful measurements provide a mean runtime. All attempts determine whether
    an estimate is complete, so an entirely failed configuration still appears
    with an unknown budget. No estimate is invented for an unobserved configuration.
    """
    # Look up each configuration's successful summary while also retaining those
    # with only failed attempts, which have no entry in summaries
    summary_by_settings = {
        tuple(row[column] for column in setting_columns): row for row in summaries
    }
    debugging_runs = [run for run in runs if run["pilot"] == "1"]
    if not debugging_runs:
        raise ValueError("No pilot configurations found")

    budgets = []
    for settings, matching_runs in group_runs_by_settings(debugging_runs, setting_columns):
        attempts = list(matching_runs)
        summary = summary_by_settings.get(settings)
        failures = sum(not successful(run) for run in attempts)
        budget = dict(zip(setting_columns, settings))
        if summary is not None:
            budget.update(summary)

        # Griffin gets a short untimed demo warmup instead of a full-input run
        planned_runs = repetitions
        if budget["method"] != "griffin-gc":
            planned_runs += warmups
        budget.update(
            {
                "planned_runs": planned_runs,
                "pilot_failures": failures,
                "budget_minutes": "NA",
                "budget_status": "incomplete",
            }
        )
        if summary is not None:
            budget["budget_minutes"] = (
                summary["process_plus_flush_seconds_mean"]
                * planned_runs
                * (1 + margin)
                / 60
            )
            if not failures:
                budget["budget_status"] = "complete"
        budgets.append(budget)
    return budgets


def parse_arguments():
    """Read input/output paths and reject invalid runtime-estimate settings."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("timings", nargs="+", type=Path)
    parser.add_argument(
        "--accounting", type=Path, help="sacct pipe-delimited export in K units"
    )
    parser.add_argument(
        "--runs-output",
        type=Path,
        help="Write every original row with joined memory/status",
    )
    parser.add_argument(
        "--estimate",
        action="store_true",
        help="Report per-configuration pilot budgets and flag incomplete totals",
    )
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument(
        "--margin",
        type=float,
        default=0.25,
        help="Fraction added to pilot-based budget",
    )
    arguments = parser.parse_args()
    if (
        arguments.warmups < 0
        or arguments.repetitions < 1
        or not math.isfinite(arguments.margin)
        or arguments.margin < 0
    ):
        parser.error(
            "Use nonnegative warmups/margin and at least one measured repetition"
        )
    return arguments


def write_tsv(destination, columns, rows):
    """Write a table in the supplied column order, filling absent fields with NA.

    Both the summary and detailed tables use the same TSV formatting. The caller
    selects stdout or opens an output file and decides whether overwriting is allowed.
    """
    writer = csv.DictWriter(
        destination,
        fieldnames=columns,
        delimiter="\t",
        lineterminator="\n",
        restval="NA",
    )
    writer.writeheader()
    writer.writerows(rows)


def report_budget_status(budgets, margin):
    """Print available runtime totals and return whether any estimate is incomplete.

    Unknown runtimes are excluded from the subtotal and named in the warnings.
    This diagnostic runs after the tables are saved, so incomplete estimates still
    leave usable records even though the command returns a failure status.
    """
    incomplete_budget = False
    for budget in budgets:
        if budget["budget_status"] != "complete":
            incomplete_budget = True
            print(
                f"Incomplete pilot budget: method={budget['method']}, cores={budget['cores']}, "
                f"node={budget.get('node', 'NA')}, failed attempts={budget['pilot_failures']}.",
                file=sys.stderr,
            )

    total_seconds = sum(
        budget["process_plus_flush_seconds_mean"] * budget["planned_runs"]
        for budget in budgets
        if budget["budget_minutes"] != "NA"
    )
    print(
        f"Observed matrix: {len(budgets)} configurations, {total_seconds / 3600:.2f} h raw, "
        f"{total_seconds * (1 + margin) / 3600:.2f} h with margin. "
        "Add staging/cleanup and any Griffin demo warmup. "
        "This excludes configurations absent from the supplied pilots.",
        file=sys.stderr,
    )
    if incomplete_budget:
        print(
            "The total is an incomplete successful-pilot subtotal, not a full job budget.",
            file=sys.stderr,
        )
    return incomplete_budget


def main():
    """Read and join runs, summarize selected measurements, then write the tables."""
    arguments = parse_arguments()
    accounting = read_accounting(arguments.accounting)
    runs = read_runs(arguments.timings, accounting)
    setting_columns, summary_columns, summaries = summarize(runs, arguments.estimate)

    # The default output is the measurement summary. Estimate mode adds runtime
    # budgets, including configurations whose attempts all failed
    output_rows = summaries
    if arguments.estimate:
        output_rows = pilot_budgets(
            runs,
            setting_columns,
            summaries,
            warmups=arguments.warmups,
            repetitions=arguments.repetitions,
            margin=arguments.margin,
        )
        summary_columns += ["planned_runs", "pilot_failures", "budget_minutes", "budget_status"]
    write_tsv(sys.stdout, summary_columns, output_rows)

    # Keep every original row, including warmups and failures, in the optional
    # detailed table. Preserve column order and refuse to overwrite an existing file
    if arguments.runs_output:
        detailed_columns = list(dict.fromkeys(column for run in runs for column in run))
        with arguments.runs_output.open("x") as destination:
            write_tsv(destination, detailed_columns, runs)

    incomplete_budget = False
    if arguments.estimate:
        incomplete_budget = report_budget_status(output_rows, arguments.margin)
    if any(row["slurm_max_rss_kib"] == "NA" for row in runs):
        print(
            "Some Slurm memory is NA. Collect completed sacct records and rerun the join.",
            file=sys.stderr,
        )
    return 1 if incomplete_budget else 0


if __name__ == "__main__":
    raise SystemExit(main())
