#!/usr/bin/env python3
"""Compare native histograms saved by benchmark_lengths.sh.

Each subdirectory of the profiles directory is named after its method and holds
that method's native histogram (see save_profile in lengths/commands.sh). READERS
lists the reader for each method.

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

# The shared checks live in benchmarks/, above this feature's script directory
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_metrics import checked_metric, pearson_correlation


def read_cfdnalab(directory):
    """Read cfDNAlab's histogram, a table with a row per length bin.

    The job requests 1 bp bins, so each bin's start (`length_start_bp`) is its
    length, and its exclusive end must be the next length.
    """
    table = read_lengths(directory / "length_counts.tsv.zst").data_frame()
    if not np.array_equal(table["length_end_bp"], table["length_start_bp"] + 1):
        raise ValueError(
            "Expected a global cfDNAlab histogram with the exact requested 1 bp length bins"
        )
    return pd.Series(table["count"].to_numpy(), index=table["length_start_bp"])


def read_cfdnapro(directory):
    """Read callLength()'s table: insert_size and All_Reads.fr_count columns.

    callLength() can add a `prop` column with each length's proportion. It is
    not a count and is ignored.
    """
    table = pd.read_csv(directory / "lengths.tsv", sep="\t")
    return pd.Series(table["All_Reads.fr_count"].to_numpy(), index=table["insert_size"])


def read_finaletoolkit(directory):
    """Read frag-length-bins' table: min, max and count columns per bin.

    The job requests 1 bp bins, so min and max must be the same length.
    """
    table = pd.read_csv(directory / "lengths.tsv", sep="\t")
    if not np.array_equal(table["min"], table["max"]):
        raise ValueError(
            "Expected FinaleToolkit bins covering a single fragment length each"
        )
    return pd.Series(table["count"].to_numpy(), index=table["min"])


def read_samtools(directory):
    """Read the insert size (IS) section of samtools stats' report.

    IS rows have the columns: label, length, all pairs, inward pairs, outward pairs
    and other pairs. The inward-pair count is used, which matches the other tools'
    forward-reverse fragments.
    """
    with (directory / "stats.txt").open() as source:
        entries = [
            line.rstrip("\n").split("\t")
            for line in source
            if line.startswith("IS\t")
        ]
    if not entries:
        raise ValueError("samtools did not report an IS histogram")
    return pd.Series(
        [float(row[3]) for row in entries], index=[int(row[1]) for row in entries]
    )


def read_picard(directory):
    """Read the all-reads FR column of Picard's filtered concordance histogram.

    This reads concordance.metrics.txt, the run made with the shared MAPQ and flag
    filters after all measurements. The timed metrics.txt is not used, because
    Picard has no MAPQ option. Picard writes metrics and a histogram in the same
    file, so only the section after its HISTOGRAM marker is read. RF and tandem
    counts are excluded. A missing FR column means no FR pairs were reported,
    provided a valid histogram section is present.
    """
    path = directory / "concordance.metrics.txt"
    lines = path.read_text().splitlines()
    section_start = next(
        (index for index, line in enumerate(lines) if line.startswith("## HISTOGRAM")),
        None,
    )
    if section_start is None:
        raise ValueError(f"Picard histogram section is missing: {path}")
    histogram_lines = []
    for line in lines[section_start + 1 :]:
        if not line:
            break
        histogram_lines.append(line)
    table = pd.read_csv(StringIO("\n".join(histogram_lines)), sep="\t")
    if "insert_size" not in table:
        raise ValueError(f"Picard histogram lacks insert_size: {path}")
    if "All_Reads.fr_count" not in table:
        return pd.Series(dtype=float)
    return pd.Series(table["All_Reads.fr_count"].to_numpy(), index=table["insert_size"])


# The reader for each method's saved histogram, and how its lengths relate to the
# requested range:
#   "exact"     every requested length and nothing else (cfDNAlab writes the requested bins)
#   "inside"    only requested lengths. Lengths without fragments may be left out
#   "anything"  may include lengths outside the range, which are dropped before any
#               normalization (samtools and Picard report their whole histogram)
READERS = {
    "cfdnalab": (read_cfdnalab, "exact"),
    "cfdnalab-autosomes": (read_cfdnalab, "exact"),
    "cfdnalab-proper": (read_cfdnalab, "exact"),
    "cfdnapro": (read_cfdnapro, "inside"),
    "cfdnapro-outer-span": (read_cfdnapro, "inside"),
    "finaletoolkit": (read_finaletoolkit, "inside"),
    "samtools": (read_samtools, "anything"),
    "picard": (read_picard, "anything"),
}


def read_profile(directory, lengths):
    """Load a method's saved histogram as counts for every requested length.

    The directory name is the method. Its reader is chosen from READERS. Lengths
    without a count are zero. Counts must be finite and nonnegative, and length
    labels must be unique integers.

    Returns a float array with a count for each value in `lengths`, in order.
    """
    if directory.name not in READERS:
        raise ValueError(f"Unknown profile directory: {directory}")
    reader, length_axis = READERS[directory.name]
    counts = reader(directory)

    if counts.index.has_duplicates:
        raise ValueError(f"Duplicate length bins in {directory.name}")
    native_lengths = np.asarray(counts.index, dtype=float)
    native_counts = counts.to_numpy(dtype=float)
    if (
        not np.isfinite(native_lengths).all()
        or not np.equal(native_lengths, np.floor(native_lengths)).all()
    ):
        raise ValueError(f"Non-integer length labels in {directory.name}")
    if not np.isfinite(native_counts).all() or (native_counts < 0).any():
        raise ValueError(f"Non-finite or negative counts in {directory.name}")
    if length_axis == "exact" and not np.array_equal(native_lengths, lengths):
        raise ValueError(
            "Expected a global cfDNAlab histogram with the exact requested 1 bp length bins"
        )
    if length_axis == "inside" and not np.isin(native_lengths, lengths).all():
        raise ValueError(
            f"Unexpected lengths outside the requested range in {directory.name}"
        )
    # Requested lengths without a count become zero. Other lengths are dropped
    return counts.reindex(lengths, fill_value=0).to_numpy(dtype=float)


def summarize_profile(counts, lengths):
    """Return the retained total, frequencies and simple biological summaries.

    All summaries use the requested length range. Short fraction means lengths
    <=150 bp. A tied mode uses the shortest tied length. Empty profiles have
    undefined frequencies and summaries, rather than fabricated zero values.
    """
    total = checked_metric(float(counts.sum()), "fragment total")
    if total == 0:
        return total, None, {"mean_bp": "NA", "mode_bp": "NA", "short_fraction": "NA"}
    frequencies = counts / total
    summaries = {
        "mean_bp": checked_metric(
            float(np.dot(lengths, frequencies)), "mean fragment length",
            minimum=float(lengths[0]), maximum=float(lengths[-1]),
        ),
        "mode_bp": int(lengths[np.argmax(counts)]),
        "short_fraction": checked_metric(
            float(frequencies[lengths <= 150].sum()), "short fraction", maximum=1.0
        ),
    }
    return total, frequencies, summaries


def main():
    """Read saved profiles and report all pairwise count/distribution differences."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "profiles", type=Path, help="The benchmark results/profiles directory"
    )
    parser.add_argument("--min-length", type=int, default=30)
    parser.add_argument("--max-length", type=int, default=1000)
    arguments = parser.parse_args()
    if not 10 <= arguments.min_length <= arguments.max_length <= 50000:
        parser.error("Use an inclusive length range within 10-50000 bp")
    lengths = np.arange(arguments.min_length, arguments.max_length + 1)
    directories = sorted(path for path in arguments.profiles.iterdir() if path.is_dir())
    if len(directories) < 2:
        parser.error("At least two saved profiles are required")
    profiles = {
        directory.name: read_profile(directory, lengths) for directory in directories
    }
    summaries = {
        name: summarize_profile(counts, lengths) for name, counts in profiles.items()
    }

    columns = [
        "left_profile",
        "right_profile",
        "comparison",
        "min_length",
        "max_length",
        "left_counted_fragments",
        "right_counted_fragments",
        "net_count_difference",
        "sum_abs_count_difference",
        "mae_counts",
        "rmse_counts",
        "fraction_identical_bins",
        "pearson_r",
        "mae_frequency",
        "rmse_frequency",
        "max_abs_frequency_difference",
        "total_variation",
        "max_cdf_difference",
        "wasserstein_bp",
        "left_mean_bp",
        "right_mean_bp",
        "left_mode_bp",
        "right_mode_bp",
        "left_short_fraction",
        "right_short_fraction",
        "status",
    ]
    writer = csv.DictWriter(
        sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    # Every pair of saved profiles is compared twice: as raw counts per length, which
    # also show how many fragments each method kept, and as frequencies within the
    # requested range, which compare the shape of the distributions
    for left_name, right_name in combinations(profiles, 2):
        left_counts, right_counts = profiles[left_name], profiles[right_name]
        left_total, left_frequency, left_summary = summaries[left_name]
        right_total, right_frequency, right_summary = summaries[right_name]
        # The label tells cross-tool comparisons apart from the two diagnostics:
        # cfDNAPro's outer-span definition and the autosome-only cfDNAlab run
        comparison = "cross-tool-or-filter"
        if {left_name, right_name} == {"cfdnapro", "cfdnapro-outer-span"}:
            comparison = "span-definition-diagnostic"
        if (left_name == "cfdnalab-autosomes") != (right_name == "cfdnalab-autosomes"):
            comparison = "different-chromosome-scope"
        count_difference = left_counts - right_counts
        row = {column: "NA" for column in columns}
        row.update(
            {
                "left_profile": left_name,
                "right_profile": right_name,
                "comparison": comparison,
                "min_length": arguments.min_length,
                "max_length": arguments.max_length,
                "left_counted_fragments": left_total,
                "right_counted_fragments": right_total,
                "net_count_difference": checked_metric(
                    left_total - right_total, "net count difference", minimum=None
                ),
                "sum_abs_count_difference": checked_metric(
                    float(np.abs(count_difference).sum()), "sum of absolute count differences"
                ),
                "mae_counts": checked_metric(
                    float(np.abs(count_difference).mean()), "count MAE"
                ),
                "rmse_counts": checked_metric(
                    float(np.sqrt(np.mean(count_difference**2))), "count RMSE"
                ),
                "fraction_identical_bins": checked_metric(
                    float(np.mean(count_difference == 0)), "identical bin fraction", maximum=1.0
                ),
                "status": "zero-count-profile",
            }
        )
        for side, summary in (("left", left_summary), ("right", right_summary)):
            for label, value in summary.items():
                row[f"{side}_{label}"] = value
        # Frequency metrics need both profiles to have fragments. Total variation is
        # half the summed absolute difference, and max_cdf_difference is the largest
        # gap between the cumulative distributions (the Kolmogorov-Smirnov statistic)
        if left_frequency is not None and right_frequency is not None:
            difference = left_frequency - right_frequency
            cumulative_difference = np.cumsum(difference)
            row["pearson_r"] = pearson_correlation(left_frequency, right_frequency)
            row.update(
                {
                    "mae_frequency": float(np.abs(difference).mean()),
                    "rmse_frequency": float(np.sqrt(np.mean(difference**2))),
                    "max_abs_frequency_difference": float(np.abs(difference).max()),
                    "total_variation": float(np.abs(difference).sum() / 2),
                    "max_cdf_difference": float(np.abs(cumulative_difference).max()),
                    # Lengths are spaced 1 bp apart, so integrating the absolute CDF
                    # difference gives the Wasserstein distance directly in bp
                    "wasserstein_bp": float(np.abs(cumulative_difference[:-1]).sum()),
                    "status": "ok",
                }
            )
            for metric in (
                "mae_frequency",
                "rmse_frequency",
                "max_abs_frequency_difference",
                "total_variation",
                "max_cdf_difference",
            ):
                row[metric] = checked_metric(row[metric], metric, maximum=1.0)
            row["wasserstein_bp"] = checked_metric(
                row["wasserstein_bp"], "Wasserstein distance",
                maximum=float(lengths[-1] - lengths[0]),
            )
        writer.writerow(row)


if __name__ == "__main__":
    main()
