#!/usr/bin/env python3
"""Compare the native global profiles saved by benchmark_end_motifs.sh.

Write one TSV row per pair of saved profiles to stdout. Compare the full labeled
motif set using frequencies, because FinaleToolkit's output has no raw counts.
FrEIA's right-end labels are complemented before pooling its two ends.
"""

import argparse
import csv
from itertools import combinations, product
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from cfdnalab import read_end_motifs


def load_profile(directory, motifs):
    """Read a native tool output onto the same inward-oriented motif axis.

    cfDNAlab's one-sided labels have a leading underscore, which is removed for
    matching to FinaleToolkit labels. Require every expected ACGT motif, including
    zero-count motifs, and reject duplicate labels or non-finite/negative values.
    Return normalized frequencies and the raw counted-end total when available.
    """
    zarr_path = directory / "end_motifs.zarr"
    tsv_path = directory / "motifs.tsv"
    freia_path = directory / "T__sample.pq"
    if zarr_path.is_dir():
        result = read_end_motifs(zarr_path)
        if result.row_mode != "global":
            raise ValueError(f"Expected a global cfDNAlab profile: {directory.name}")
        table = result.data_frame()
        entries = zip(table["motif"], table["count"])
        contains_counts = True
    elif tsv_path.is_file():
        with tsv_path.open() as source:
            entries = list(csv.reader(source, delimiter="\t"))
        contains_counts = False
    elif freia_path.is_file():
        table = pd.read_parquet(freia_path)
        # CountLen is the count in each length bin. CountSamp is repeated across
        # length bins and must not be summed. The native --both call keeps each
        # end separate. Its right-end strings are reversed but not complemented
        complement = str.maketrans("ACGT", "TGCA")
        counts = dict.fromkeys(motifs, 0)
        for row in table.itertuples(index=False):
            if row.WhichEnd == "leftmost":
                motif = row.Base
            elif row.WhichEnd == "rightmost":
                motif = row.Base.translate(complement)
            else:
                raise ValueError(f"Unexpected FrEIA end label: {row.WhichEnd}")
            if motif not in counts:
                raise ValueError(f"Unexpected FrEIA motif: {motif}")
            counts[motif] += row.CountLen
        entries = counts.items()
        contains_counts = True
    else:
        raise ValueError(f"No supported native motif profile in {directory}")

    values = {}
    for label, value in entries:
        motif = str(label).removeprefix("_")
        if motif in values:
            raise ValueError(f"Duplicate motif {motif} in {directory.name}")
        values[motif] = float(value)
    if set(values) != set(motifs):
        raise ValueError(f"{directory.name} does not contain exactly the expected ACGT motif set")
    ordered_values = np.array([values[motif] for motif in motifs], dtype=float)
    if not np.isfinite(ordered_values).all() or (ordered_values < 0).any():
        raise ValueError(f"Non-finite or negative values in {directory.name}")
    total = float(ordered_values.sum())
    if total <= 0:
        raise ValueError(f"No counted motif mass in {directory.name}")
    if not contains_counts and not np.isclose(total, 1.0, rtol=1e-6, atol=1e-9):
        raise ValueError(f"FinaleToolkit frequencies do not sum to 1 in {directory.name}: {total}")
    return ordered_values / total, total if contains_counts else "NA"


def motif_diversity(frequencies):
    """Calculate normalized Shannon entropy, with zero-frequency terms omitted."""
    positive = frequencies[frequencies > 0]
    return float(-np.sum(positive * np.log(positive)) / np.log(len(frequencies)))


def main():
    """Load the small native profiles and print all pairwise comparisons."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", type=Path, help="The benchmark results/profiles directory")
    parser.add_argument("--k", type=int, choices=(3, 4), default=4, help="Motif length used by the commands")
    arguments = parser.parse_args()
    motifs = ["".join(bases) for bases in product("ACGT", repeat=arguments.k)]
    directories = sorted(path for path in arguments.profiles.iterdir() if path.is_dir())
    if len(directories) < 2:
        parser.error("At least two saved tool profiles are required")
    profiles = {path.name: load_profile(path, motifs) for path in directories}

    columns = [
        "left_profile", "right_profile", "k", "motifs", "left_counted_ends", "right_counted_ends",
        "pearson_r", "mae_frequency", "rmse_frequency", "max_abs_frequency_difference",
        "total_variation", "left_mds", "right_mds", "abs_mds_difference",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for left_name, right_name in combinations(profiles, 2):
        left_frequency, left_total = profiles[left_name]
        right_frequency, right_total = profiles[right_name]
        difference = left_frequency - right_frequency
        # Correlation is undefined for constant profiles, even if they are equal
        correlation = "NA"
        if np.std(left_frequency) > 0 and np.std(right_frequency) > 0:
            correlation = float(np.corrcoef(left_frequency, right_frequency)[0, 1])
        left_mds = motif_diversity(left_frequency)
        right_mds = motif_diversity(right_frequency)
        row = {
            "left_profile": left_name,
            "right_profile": right_name,
            "k": arguments.k,
            "motifs": len(motifs),
            "left_counted_ends": left_total,
            "right_counted_ends": right_total,
            "pearson_r": correlation,
            "mae_frequency": float(np.mean(np.abs(difference))),
            "rmse_frequency": float(np.sqrt(np.mean(difference**2))),
            "max_abs_frequency_difference": float(np.max(np.abs(difference))),
            "total_variation": float(np.sum(np.abs(difference)) / 2),
            "left_mds": left_mds,
            "right_mds": right_mds,
            "abs_mds_difference": abs(left_mds - right_mds),
        }
        writer.writerow(row)


if __name__ == "__main__":
    main()
