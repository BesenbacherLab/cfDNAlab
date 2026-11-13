#!/usr/bin/env python3
"""Compare the native global profiles saved by benchmark_end_motifs.sh.

Write one TSV row per pair of saved profiles to stdout. Compare the full labeled
motif set using frequencies, because FinaleToolkit's output has no raw counts.
cfDNAPro supplies native left-end profiles. Comparisons against both ends
use an explicitly labeled side-mirroring approximation, without doubling counts.
"""

import argparse
import csv
from itertools import combinations, product
from pathlib import Path
import sys

import numpy as np
import pandas as pd

# The shared checks live in benchmarks/, above this feature's script directory
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_metrics import checked_metric, pearson_correlation


def load_profile(directory, motifs):
    """Read a native tool output onto the same inward-oriented motif axis.

    cfDNAlab's one-sided labels have a leading underscore, which is removed for
    matching to FinaleToolkit labels. All formats must contain every ACGT motif.
    Reject invalid labels/values and return normalized frequencies,
    native counted-end totals when available, and the native end selection.
    """
    # The saved native file identifies the tool (see save_profile in ends/commands.sh)
    zarr_path = directory / "end_motifs.zarr"
    tsv_path = directory / "motifs.tsv"
    cfdnapro_path = directory / "cfdnapro-motifs.tsv"
    end_selection = "both"
    if zarr_path.is_dir():
        # cfDNAlab: counts of both fragment ends, read with the cfdnalab package
        from cfdnalab import read_end_motifs

        result = read_end_motifs(zarr_path)
        if result.row_mode != "global":
            raise ValueError(f"Expected a global cfDNAlab profile: {directory.name}")
        table = result.data_frame()
        entries = zip(table["motif"], table["count"])
        contains_counts = True
    elif tsv_path.is_file():
        # FinaleToolkit: a headerless motif/frequency table for both ends. It has
        # no raw counts, so its frequencies must already sum to 1
        with tsv_path.open() as source:
            entries = list(csv.reader(source, delimiter="\t"))
        contains_counts = False
    elif cfdnapro_path.is_file():
        # cfDNAPro: callMotif()'s table of left-end motifs, with counts in column n
        with cfdnapro_path.open() as source:
            table = list(csv.DictReader(source, delimiter="\t"))
        entries = [(row["motif"], row["n"]) for row in table]
        contains_counts = True
        end_selection = "left"
    else:
        raise ValueError(f"No supported native motif profile in {directory}")

    values = {}
    for label, value in entries:
        motif = str(label).removeprefix("_")
        if motif in values:
            raise ValueError(f"Duplicate motif {motif} in {directory.name}")
        values[motif] = float(value)
    if set(values) != set(motifs):
        raise ValueError(
            f"{directory.name} does not contain exactly the expected ACGT motif set"
        )
    ordered_values = np.array([values[motif] for motif in motifs], dtype=float)
    if not np.isfinite(ordered_values).all() or (ordered_values < 0).any():
        raise ValueError(f"Non-finite or negative values in {directory.name}")
    total = checked_metric(
        float(ordered_values.sum()), f"motif total for {directory.name}"
    )
    if total <= 0:
        raise ValueError(f"No counted motif mass in {directory.name}")
    if not contains_counts and not np.isclose(total, 1.0, rtol=1e-6, atol=1e-9):
        raise ValueError(
            f"FinaleToolkit frequencies do not sum to 1 in {directory.name}: {total}"
        )
    return ordered_values / total, total if contains_counts else "NA", end_selection


def motif_diversity(frequencies):
    """Calculate normalized Shannon entropy, with zero-frequency terms omitted."""
    positive = frequencies[frequencies > 0]
    value = float(-np.sum(positive * np.log(positive)) / np.log(len(frequencies)))
    return checked_metric(value, "normalized motif diversity", maximum=1.0)


def main():
    """Load the small native profiles and print all pairwise comparisons."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "profiles", type=Path, help="The benchmark results/profiles directory"
    )
    parser.add_argument(
        "--k",
        type=int,
        choices=(3, 4),
        default=4,
        help="Motif length used by the commands",
    )
    arguments = parser.parse_args()
    motifs = ["".join(bases) for bases in product("ACGT", repeat=arguments.k)]
    directories = sorted(path for path in arguments.profiles.iterdir() if path.is_dir())
    if len(directories) < 2:
        parser.error("At least two saved tool profiles are required")
    profiles = {path.name: load_profile(path, motifs) for path in directories}

    columns = [
        "left_profile",
        "right_profile",
        "k",
        "motifs",
        "left_counted_ends",
        "right_counted_ends",
        "left_end_selection",
        "right_end_selection",
        "side_comparison",
        "pearson_r",
        "spearman_r",
        "max_abs_log2_ratio_abundant",
        "mae_frequency",
        "rmse_frequency",
        "max_abs_frequency_difference",
        "total_variation",
        "left_mds",
        "right_mds",
        "abs_mds_difference",
    ]
    writer = csv.DictWriter(
        sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for left_name, right_name in combinations(profiles, 2):
        left_frequency, left_total, left_ends = profiles[left_name]
        right_frequency, right_total, right_ends = profiles[right_name]
        difference = left_frequency - right_frequency
        # Correlation is undefined for constant profiles, even if they are equal
        correlation = pearson_correlation(left_frequency, right_frequency)
        left_ranks = pd.Series(left_frequency).rank().to_numpy()
        right_ranks = pd.Series(right_frequency).rank().to_numpy()
        rank_correlation = pearson_correlation(left_ranks, right_ranks)
        # Include a motif when either profile reaches 1e-4. Zero on the other
        # side gives an infinite ratio, without introducing a pseudocount
        abundant = np.maximum(left_frequency, right_frequency) >= 1e-4
        log_ratio = "NA"
        if abundant.any():
            if (
                (left_frequency[abundant] == 0) | (right_frequency[abundant] == 0)
            ).any():
                log_ratio = "Inf"
            else:
                # Separate logarithms avoid overflow in a ratio of positive frequencies
                log_ratio = checked_metric(
                    float(np.abs(
                        np.log2(left_frequency[abundant]) - np.log2(right_frequency[abundant])
                    ).max()),
                    "absolute motif log2 ratio",
                )
        left_mds = motif_diversity(left_frequency)
        right_mds = motif_diversity(right_frequency)
        row = {
            "left_profile": left_name,
            "right_profile": right_name,
            "k": arguments.k,
            "motifs": len(motifs),
            "left_counted_ends": left_total,
            "right_counted_ends": right_total,
            "left_end_selection": left_ends,
            "right_end_selection": right_ends,
            "side_comparison": "same-sides"
            if left_ends == right_ends
            else "approximate-side-mirroring",
            "pearson_r": correlation,
            "spearman_r": rank_correlation,
            "max_abs_log2_ratio_abundant": log_ratio,
            "mae_frequency": float(np.mean(np.abs(difference))),
            "rmse_frequency": float(np.sqrt(np.mean(difference**2))),
            "max_abs_frequency_difference": float(np.max(np.abs(difference))),
            "total_variation": float(np.sum(np.abs(difference)) / 2),
            "left_mds": left_mds,
            "right_mds": right_mds,
            "abs_mds_difference": abs(left_mds - right_mds),
        }
        # These quantities are bounded because both profiles sum to one
        for metric in (
            "mae_frequency",
            "rmse_frequency",
            "max_abs_frequency_difference",
            "total_variation",
            "abs_mds_difference",
        ):
            row[metric] = checked_metric(row[metric], metric, maximum=1.0)
        writer.writerow(row)


if __name__ == "__main__":
    main()
