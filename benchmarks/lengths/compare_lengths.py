#!/usr/bin/env python3
"""Compare native histograms saved by benchmark_lengths.sh.

Write a TSV row per pair of profiles to stdout. Use every integer length in the
requested inclusive range, adding zero counts for absent bins. Compare raw counts
and frequencies normalized within that range. Do not smooth or shift profiles.
"""

import argparse
import csv
from io import StringIO
from itertools import combinations
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from cfdnalab import read_lengths


def read_picard(path):
    """Select the native all-reads FR histogram, excluding RF and tandem counts.

    Picard writes metrics and a histogram in the same file. Read only the section
    following its HISTOGRAM marker. A missing FR column means no FR pairs were
    reported, provided a valid histogram section is present.
    """
    lines = path.read_text().splitlines()
    section_start = next((index for index, line in enumerate(lines) if line.startswith("## HISTOGRAM")), None)
    if section_start is None:
        raise ValueError(f"Picard histogram section is missing: {path}")
    histogram_lines = []
    for line in lines[section_start + 1:]:
        if not line:
            break
        histogram_lines.append(line)
    table = pd.read_csv(StringIO("\n".join(histogram_lines)), sep="\t")
    if "insert_size" not in table:
        raise ValueError(f"Picard histogram lacks insert_size: {path}")
    if "All_Reads.fr_count" not in table:
        return pd.Series(dtype=float)
    return pd.Series(table["All_Reads.fr_count"].to_numpy(), index=table["insert_size"])


def read_profile(directory, lengths):
    """Load a native global count table and align its bins to integer lengths.

    cfDNAlab's axis must exactly match the requested range. Other formats can
    omit zero bins. samtools and Picard also report lengths outside our range,
    which are discarded before normalizing. No normalization uses those tails.
    """
    if directory.name.startswith("cfdnalab"):
        result = read_lengths(directory / "length_counts.tsv.zst")
        table = result.data_frame()
        if not np.array_equal(table["length_start_bp"], lengths) or not np.array_equal(
            table["length_end_bp"], lengths + 1
        ):
            raise ValueError("Expected a global cfDNAlab histogram with the exact requested 1 bp length bins")
        counts = pd.Series(table["count"].to_numpy(), index=table["length_start_bp"])
    elif directory.name.startswith("cfdnapro"):
        table = pd.read_csv(directory / "lengths.tsv", sep="\t")
        counts = pd.Series(table["All_Reads.fr_count"].to_numpy(), index=table["insert_size"])
    elif directory.name == "finaletoolkit":
        table = pd.read_csv(directory / "lengths.tsv", sep="\t")
        if not np.array_equal(table["min"], table["max"]):
            raise ValueError("Expected FinaleToolkit bins covering a single fragment length each")
        counts = pd.Series(table["count"].to_numpy(), index=table["min"])
    elif directory.name == "samtools":
        with (directory / "stats.txt").open() as source:
            # IS columns are label, length, all pairs, inward, outward and other
            entries = [line.rstrip("\n").split("\t") for line in source if line.startswith("IS\t")]
        if not entries:
            raise ValueError("samtools did not report an IS histogram")
        counts = pd.Series([float(row[3]) for row in entries], index=[int(row[1]) for row in entries])
    elif directory.name == "picard":
        counts = read_picard(directory / "metrics.txt")
    else:
        raise ValueError(f"Unknown profile directory: {directory}")

    if counts.index.has_duplicates:
        raise ValueError(f"Duplicate length bins in {directory.name}")
    native_lengths = np.asarray(counts.index, dtype=float)
    native_counts = counts.to_numpy(dtype=float)
    if not np.isfinite(native_lengths).all() or not np.equal(native_lengths, np.floor(native_lengths)).all():
        raise ValueError(f"Non-integer length labels in {directory.name}")
    if not np.isfinite(native_counts).all() or (native_counts < 0).any():
        raise ValueError(f"Non-finite or negative counts in {directory.name}")
    if directory.name not in {"samtools", "picard"} and not np.isin(native_lengths, lengths).all():
        raise ValueError(f"Unexpected lengths outside the requested range in {directory.name}")
    return counts.reindex(lengths, fill_value=0).to_numpy(dtype=float)


def summarize_profile(counts, lengths):
    """Return the retained total, frequencies and simple biological summaries.

    All summaries use the requested length range. Short fraction means lengths
    <=150 bp. A tied mode uses the shortest tied length. Empty profiles have
    undefined frequencies and summaries, rather than fabricated zero values.
    """
    total = float(counts.sum())
    if total == 0:
        return total, None, {"mean_bp": "NA", "mode_bp": "NA", "short_fraction": "NA"}
    frequencies = counts / total
    summaries = {
        "mean_bp": float(np.dot(lengths, frequencies)),
        "mode_bp": int(lengths[np.argmax(counts)]),
        "short_fraction": float(frequencies[lengths <= 150].sum()),
    }
    return total, frequencies, summaries


def main():
    """Read saved profiles and report all pairwise count/distribution differences."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", type=Path, help="The benchmark results/profiles directory")
    parser.add_argument("--min-length", type=int, default=30)
    parser.add_argument("--max-length", type=int, default=1000)
    arguments = parser.parse_args()
    if not 10 <= arguments.min_length <= arguments.max_length <= 50000:
        parser.error("Use an inclusive length range within 10-50000 bp")
    lengths = np.arange(arguments.min_length, arguments.max_length + 1)
    directories = sorted(path for path in arguments.profiles.iterdir() if path.is_dir())
    if len(directories) < 2:
        parser.error("At least two saved profiles are required")
    profiles = {directory.name: read_profile(directory, lengths) for directory in directories}
    summaries = {name: summarize_profile(counts, lengths) for name, counts in profiles.items()}

    columns = [
        "left_profile", "right_profile", "comparison", "min_length", "max_length",
        "left_counted_fragments", "right_counted_fragments", "net_count_difference",
        "sum_abs_count_difference", "mae_counts", "rmse_counts", "fraction_identical_bins",
        "pearson_r", "mae_frequency", "rmse_frequency", "max_abs_frequency_difference",
        "total_variation", "max_cdf_difference", "wasserstein_bp",
        "left_mean_bp", "right_mean_bp", "left_mode_bp", "right_mode_bp",
        "left_short_fraction", "right_short_fraction", "status",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for left_name, right_name in combinations(profiles, 2):
        left_counts, right_counts = profiles[left_name], profiles[right_name]
        left_total, left_frequency, left_summary = summaries[left_name]
        right_total, right_frequency, right_summary = summaries[right_name]
        comparison = "cross-tool-or-filter"
        if {left_name, right_name} == {"cfdnapro", "cfdnapro-outer-span"}:
            comparison = "span-definition-diagnostic"
        if (left_name == "cfdnalab-autosomes") != (right_name == "cfdnalab-autosomes"):
            comparison = "different-chromosome-scope"
        count_difference = left_counts - right_counts
        row = {column: "NA" for column in columns}
        row.update({
            "left_profile": left_name, "right_profile": right_name, "comparison": comparison,
            "min_length": arguments.min_length, "max_length": arguments.max_length,
            "left_counted_fragments": left_total, "right_counted_fragments": right_total,
            "net_count_difference": left_total - right_total,
            "sum_abs_count_difference": float(np.abs(count_difference).sum()),
            "mae_counts": float(np.abs(count_difference).mean()),
            "rmse_counts": float(np.sqrt(np.mean(count_difference**2))),
            "fraction_identical_bins": float(np.mean(count_difference == 0)),
            "status": "zero-count-profile",
        })
        for side, summary in (("left", left_summary), ("right", right_summary)):
            for label, value in summary.items():
                row[f"{side}_{label}"] = value
        if left_frequency is not None and right_frequency is not None:
            difference = left_frequency - right_frequency
            cumulative_difference = np.cumsum(difference)
            if np.std(left_frequency) > 0 and np.std(right_frequency) > 0:
                row["pearson_r"] = float(np.corrcoef(left_frequency, right_frequency)[0, 1])
            row.update({
                "mae_frequency": float(np.abs(difference).mean()),
                "rmse_frequency": float(np.sqrt(np.mean(difference**2))),
                "max_abs_frequency_difference": float(np.abs(difference).max()),
                "total_variation": float(np.abs(difference).sum() / 2),
                "max_cdf_difference": float(np.abs(cumulative_difference).max()),
                # Lengths are spaced 1 bp apart, so integrating the absolute CDF
                # difference gives the Wasserstein distance directly in bp
                "wasserstein_bp": float(np.abs(cumulative_difference[:-1]).sum()),
                "status": "ok",
            })
        writer.writerow(row)


if __name__ == "__main__":
    main()
