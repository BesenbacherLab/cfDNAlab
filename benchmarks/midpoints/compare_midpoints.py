#!/usr/bin/env python3
"""Compare the native midpoint profiles saved by benchmark_midpoints.sh.

Write a TSV row per group to stdout. Griffin always writes mean-one profiles, so
normalize cfDNAlab's native profiles over the same complete window before comparison.
Keep every saved bin, including zeros and native smoothing undershoots.
Do not shift, smooth or rebin profiles during comparison.
Use --gc-corrected to compare the separate cfdnalab-gc and griffin-gc workloads.
"""

import argparse
import csv
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from cfdnalab import read_midpoints

# The shared checks live in benchmarks/, above this feature's script directory
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_metrics import checked_metric, pearson_correlation


def read_cfdnalab(path):
    """Read native profiles, eligible site counts, window width and bin size.

    Require a single fragment length bin and a contiguous interval-relative
    position axis. Group metadata supplies the actual number of eligible sites.
    """
    result = read_midpoints(path)
    if len(result.length_bins()) != 1:
        raise ValueError("Expected a single fragment length bin in cfDNAlab output")
    positions = result.positions()
    if len(positions) == 0:
        raise ValueError("Expected a nonzero profile width")
    width = int(positions["position_bin_end_bp"].iloc[-1])
    bin_size = int(positions["position_bin_end_bp"].iloc[0])
    if width <= 0 or bin_size <= 0:
        raise ValueError("Expected positive midpoint window width and bin size")
    expected = np.arange(0, width, bin_size)
    if not np.array_equal(
        positions["position_bin_start_bp"], expected
    ) or not np.array_equal(
        positions["position_bin_end_bp"], np.minimum(expected + bin_size, width)
    ):
        raise ValueError("cfDNAlab positions must be contiguous, regularly spaced bins")
    counts = result.counts_array()[:, 0, :].astype(float)
    # Savitzky-Golay smoothing can give small negative profile values even though
    # underlying midpoint counts are nonnegative. Keep the native output unchanged
    if not np.isfinite(counts).all():
        raise ValueError("Expected finite midpoint profile values")
    profiles = {}
    for row, profile in zip(result.group_metadata().itertuples(index=False), counts):
        group = str(row.group_name)
        if group in profiles or row.eligible_intervals < 0:
            raise ValueError(f"Duplicate group or invalid eligible site count: {group}")
        profiles[group] = (profile, int(row.eligible_intervals))
    return profiles, width, bin_size


def read_griffin(path, width, bin_size):
    """Read native mean-one profiles by their signed position labels.

    Griffin's other columns contain features and metadata. Select only the exact
    shared position range and check that the file contains no extra positions.
    """
    first_position = -(width // 2)
    # Griffin rounds both signed bounds inward to multiples of its step. Its
    # bins match cfDNAlab's BED-relative bins only when both bounds already align
    # Check this before reading values rather than comparing shifted profiles
    if first_position % bin_size or (first_position + width) % bin_size:
        raise ValueError(
            "The BED window does not align with Griffin's signed bin boundaries. "
            "For 15 bp bins, use matching bounds such as +/-4995 bp (9990 bp total). "
            "The comparator does not trim or shift profiles."
        )
    table = pd.read_csv(path, sep="\t", dtype={"site_name": str}, keep_default_na=False)
    position_columns = [
        str(position) for position in range(first_position, first_position + width, bin_size)
    ]
    native_positions = {
        column for column in table.columns if column.lstrip("-").isdigit()
    }
    if native_positions != set(position_columns):
        raise ValueError("Griffin position labels do not match the cfDNAlab window")
    profiles = {}
    for row_index, row in table.iterrows():
        group = row["site_name"]
        site_count = float(row["number_of_sites"])
        if (
            group in profiles
            or not np.isfinite(site_count)
            or site_count <= 0
            or not site_count.is_integer()
        ):
            raise ValueError(f"Duplicate group or invalid site count: {group}")
        if str(row["CNA_normalization"]).lower() != "false":
            raise ValueError("Expected Griffin's default without site-wise CNA normalization")
        # Griffin writes an empty field for a NaN, including an empty group's
        # undefined mean-one profile. Preserve that distinction from a real zero
        profile = pd.to_numeric(row[position_columns].replace("", np.nan)).to_numpy(
            dtype=float
        )
        mean_per_site = float(row["mean_reads_per_bp_in_normalization_window"] or "nan")
        if np.isinf(mean_per_site) or mean_per_site < 0:
            raise ValueError(f"Invalid Griffin mean-count metadata for {group}")
        if np.isinf(profile).any():
            raise ValueError(f"Infinite Griffin profile values for {group}")
        if np.isfinite(profile).all():
            # Its native TSV rounds values to five decimal places
            # A nonempty sparse profile can have mean-count metadata rounded
            # to zero. Determine emptiness from the profile, not that metadata
            if not np.isclose(profile.mean(), 1.0, rtol=0, atol=1e-5):
                raise ValueError(f"Griffin profile does not have mean one: {group}")
        # Native masking can leave undefined positions. Keep them as NaN and
        # report unavailable metrics for that group rather than choosing a subset
        if not np.isfinite(mean_per_site) and not np.isnan(profile).all():
            raise ValueError(f"Invalid Griffin mean-count metadata for {group}")
        profiles[group] = (profile, int(site_count), mean_per_site)
    return profiles


def main():
    """Check shared groups/windows and compare native profiles over all positions."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "profiles", type=Path, help="The benchmark results/profiles directory"
    )
    parser.add_argument(
        "--gc-corrected", action="store_true",
        help="Compare cfdnalab-gc with griffin-gc, using Griffin's GC-corrected output",
    )
    arguments = parser.parse_args()
    # Separate method directories prevent pairing a corrected profile with an
    # uncorrected profile. Their native count readers and normalization are shared
    cfdnalab_method = "cfdnalab-gc" if arguments.gc_corrected else "cfdnalab"
    griffin_method = "griffin-gc" if arguments.gc_corrected else "griffin"
    griffin_filename = (
        "sample.GC_corrected.coverage.tsv"
        if arguments.gc_corrected
        else "sample.uncorrected.coverage.tsv"
    )
    cfdnalab, width, bin_size = read_cfdnalab(
        arguments.profiles / cfdnalab_method / "midpoint_profiles.zarr"
    )
    griffin = read_griffin(
        arguments.profiles / griffin_method / griffin_filename, width, bin_size
    )
    if set(cfdnalab) != set(griffin):
        raise ValueError("The tools did not report the same group set")
    columns = [
        "group",
        "cfdnalab_sites",
        "griffin_sites",
        "positions",
        "window_bp",
        "resolution_bp",
        "gc_correction",
        "cfdnalab_mean_midpoints_per_site_bp",
        "griffin_mean_midpoints_per_site_bp",
        "pearson_r",
        "pearson_r_squared",
        "mae_mean_one",
        "rmse_mean_one",
        "max_abs_mean_one",
        "status",
    ]
    writer = csv.DictWriter(
        sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for group in sorted(cfdnalab):
        counts, site_count = cfdnalab[group]
        griffin_profile, griffin_sites, griffin_mean = griffin[group]
        # Native blacklist handling differs. Report both eligible site counts
        # rather than silently changing either tool's site selection
        count_mean = checked_metric(float(counts.mean()), f"midpoint mean for {group}", minimum=None)
        row = {
            "group": group,
            "cfdnalab_sites": site_count,
            "griffin_sites": griffin_sites,
            "positions": len(counts),
            "window_bp": width,
            "resolution_bp": bin_size,
            "gc_correction": "corrected" if arguments.gc_corrected else "none",
            "cfdnalab_mean_midpoints_per_site_bp": count_mean / site_count if site_count else "NA",
            "griffin_mean_midpoints_per_site_bp": griffin_mean if np.isfinite(griffin_mean) else "NA",
            "pearson_r": "NA",
            "pearson_r_squared": "NA",
            "mae_mean_one": "NA",
            "rmse_mean_one": "NA",
            "max_abs_mean_one": "NA",
            "status": "both-zero-coverage",
        }
        if np.isnan(griffin_profile).all() and not np.isfinite(griffin_mean):
            row["status"] = "griffin-undefined-profile"
        elif count_mean == 0 and np.isfinite(griffin_profile).all():
            row["status"] = "cfdnalab-zero-coverage"
        elif count_mean < 0:
            row["status"] = "cfdnalab-nonpositive-profile-mean"
        elif count_mean > 0 and np.isnan(griffin_profile).all():
            row["status"] = "griffin-zero-coverage"
        elif np.isnan(griffin_profile).any() and not np.isnan(griffin_profile).all():
            row["status"] = "griffin-undefined-positions"
        if count_mean > 0 and np.isfinite(griffin_profile).all():
            normalized = counts / count_mean
            difference = normalized - griffin_profile
            row["status"] = "ok"
            row["pearson_r"] = pearson_correlation(normalized, griffin_profile)
            if row["pearson_r"] != "NA":
                row["pearson_r_squared"] = checked_metric(
                    row["pearson_r"] ** 2, "midpoint Pearson r squared", maximum=1.0
                )
            row["mae_mean_one"] = checked_metric(
                float(np.mean(np.abs(difference))), "midpoint MAE"
            )
            row["rmse_mean_one"] = checked_metric(
                float(np.sqrt(np.mean(difference**2))), "midpoint RMSE"
            )
            row["max_abs_mean_one"] = checked_metric(
                float(np.max(np.abs(difference))), "maximum midpoint difference"
            )
        writer.writerow(row)


if __name__ == "__main__":
    main()
