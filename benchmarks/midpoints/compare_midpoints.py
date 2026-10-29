#!/usr/bin/env python3
"""Compare the native midpoint profiles saved by benchmark_midpoints.sh.

Write a TSV row per group to stdout. Griffin always writes mean-one profiles, so
normalize cfDNAlab's counts over the same complete window before comparison.
Keep every 1 bp position, including zeros. Do not shift, smooth or rebin profiles.
"""

import argparse
import csv
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from cfdnalab import read_midpoints


def read_cfdnalab(path):
    """Read unsmoothed 1 bp counts, site counts and the window width.

    Require a single fragment length bin and a contiguous interval-relative
    position axis. Group metadata supplies the actual number of eligible sites.
    """
    result = read_midpoints(path)
    if len(result.length_bins()) != 1:
        raise ValueError("Expected a single fragment length bin in cfDNAlab output")
    positions = result.positions()
    width = len(positions)
    expected = np.arange(width)
    if width == 0 or width % 2:
        raise ValueError("Expected an even, nonzero profile width")
    if not np.array_equal(positions["position_bin_start_bp"], expected) or not np.array_equal(
        positions["position_bin_end_bp"], expected + 1
    ):
        raise ValueError("cfDNAlab positions must be contiguous 1 bp bins")
    counts = result.counts_array()[:, 0, :].astype(float)
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("Expected finite, nonnegative unsmoothed midpoint counts")
    profiles = {}
    for row, profile in zip(result.group_metadata().itertuples(index=False), counts):
        group = str(row.group_name)
        if group in profiles or row.eligible_intervals <= 0:
            raise ValueError(f"Duplicate group or no eligible sites: {group}")
        profiles[group] = (profile, int(row.eligible_intervals))
    return profiles, width


def read_griffin(path, width):
    """Read native mean-one profiles by their signed position labels.

    Griffin's other columns contain features and metadata. Select only the exact
    shared position range and check that the file contains no extra positions.
    """
    table = pd.read_csv(path, sep="\t", dtype={"site_name": str}, keep_default_na=False)
    position_columns = [str(position) for position in range(-width // 2, width // 2)]
    native_positions = {column for column in table.columns if column.lstrip("-").isdigit()}
    if native_positions != set(position_columns):
        raise ValueError("Griffin position labels do not match the cfDNAlab window")
    profiles = {}
    for row_index, row in table.iterrows():
        group = row["site_name"]
        site_count = float(row["number_of_sites"])
        if group in profiles or not np.isfinite(site_count) or site_count <= 0 or not site_count.is_integer():
            raise ValueError(f"Duplicate group or invalid site count: {group}")
        if str(row["smoothing"]).lower() != "false" or str(row["CNA_normalization"]).lower() != "false":
            raise ValueError("Expected unsmoothed Griffin output without site-wise CNA normalization")
        # Griffin writes an empty field for a NaN, including an empty group's
        # undefined mean-one profile. Preserve that distinction from a real zero
        profile = pd.to_numeric(row[position_columns].replace("", np.nan)).to_numpy(dtype=float)
        mean_per_site = float(row["mean_reads_per_bp_in_normalization_window"])
        if not np.isfinite(mean_per_site) or mean_per_site < 0:
            raise ValueError(f"Invalid Griffin mean-count metadata for {group}")
        if np.isfinite(profile).all() and not (profile < 0).any():
            # Its native TSV rounds values to five decimal places
            # A nonempty sparse profile can have mean-count metadata rounded
            # to zero. Determine emptiness from the profile, not that metadata
            if not np.isclose(profile.mean(), 1.0, rtol=0, atol=1e-5):
                raise ValueError(f"Griffin profile does not have mean one: {group}")
        elif mean_per_site != 0 or not np.isnan(profile).all():
            raise ValueError(f"Unexpected zero-coverage Griffin output for {group}")
        profiles[group] = (profile, int(site_count), mean_per_site)
    return profiles


def main():
    """Check site/window agreement and write concordance metrics per group."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", type=Path, help="The benchmark results/profiles directory")
    arguments = parser.parse_args()
    cfdnalab, width = read_cfdnalab(arguments.profiles / "cfdnalab" / "midpoint_profiles.zarr")
    griffin = read_griffin(arguments.profiles / "griffin" / "sample.uncorrected.coverage.tsv", width)
    if set(cfdnalab) != set(griffin):
        raise ValueError("The tools did not report the same group set")
    columns = [
        "group", "sites", "positions", "cfdnalab_mean_midpoints_per_site_bp",
        "griffin_mean_midpoints_per_site_bp", "pearson_r", "mae_mean_one",
        "rmse_mean_one", "max_abs_mean_one", "status",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for group in sorted(cfdnalab):
        counts, site_count = cfdnalab[group]
        griffin_profile, griffin_sites, griffin_mean = griffin[group]
        if site_count != griffin_sites:
            raise ValueError(f"Different site counts for {group}: {site_count} vs {griffin_sites}")
        count_mean = float(counts.mean())
        row = {
            "group": group, "sites": site_count, "positions": width,
            "cfdnalab_mean_midpoints_per_site_bp": count_mean / site_count,
            "griffin_mean_midpoints_per_site_bp": griffin_mean,
            "pearson_r": "NA", "mae_mean_one": "NA", "rmse_mean_one": "NA",
            "max_abs_mean_one": "NA", "status": "both-zero-coverage",
        }
        if count_mean == 0 and np.isfinite(griffin_profile).all():
            row["status"] = "cfdnalab-zero-coverage"
        elif count_mean > 0 and np.isnan(griffin_profile).all():
            row["status"] = "griffin-zero-coverage"
        if count_mean > 0 and np.isfinite(griffin_profile).all():
            normalized = counts / count_mean
            difference = normalized - griffin_profile
            row["status"] = "ok"
            if np.std(normalized) > 0 and np.std(griffin_profile) > 0:
                row["pearson_r"] = float(np.corrcoef(normalized, griffin_profile)[0, 1])
            row["mae_mean_one"] = float(np.mean(np.abs(difference)))
            row["rmse_mean_one"] = float(np.sqrt(np.mean(difference**2)))
            row["max_abs_mean_one"] = float(np.max(np.abs(difference)))
        writer.writerow(row)


if __name__ == "__main__":
    main()
