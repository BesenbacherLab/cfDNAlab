#!/usr/bin/env python3
"""Compare saved coverage profiles that measure the same quantity.

benchmark_coverage.sh runs this after all timed runs, on the results' profiles/
directory. Each subdirectory is named after its method and holds that method's
native output (see save_profile in coverage/commands.sh).

Methods are grouped by the quantity they measure (QUANTITY_BY_METHOD below), and
only methods with the same quantity are compared:

- Per-base read depth and per-base fragment coverage are compared base by base on
  each pair's shared chromosomes. Tracks are never expanded into genomic arrays.
  They are read as runs of constant depth, and each run counts with its number of
  bases. Temporary per-chromosome run files reconcile different native chromosome
  orders and keep memory bounded. Base-resolution MAE is the primary measure.
- Regional read depth (covered bases per window) and regional fragment counts are
  compared window by window, at the exact coordinates of the shared regions.bed
  beside profiles/.

Outputs:

- stdout: a TSV row of summary metrics per compared pair.
- residuals.tsv, beside profiles/: for each per-base pair, the number of bases
  with each residual value (left - right).
- coverage_jaccard.tsv, beside profiles/: for each per-base pair and chromosome,
  the bases each side covers and their overlap.

Requires the zstd executable for cfDNAlab's compressed tables.
"""

import argparse
from collections import Counter
from contextlib import closing, contextmanager
import csv
import gzip
from itertools import combinations
import math
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
from typing import NamedTuple

# The shared checks live in benchmarks/, above this feature's script directory
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_metrics import checked_metric, correlation_value


# The quantity each coverage method measures. Only methods with the same quantity
# are compared with each other. The labels appear in the quantity column of the output
QUANTITY_BY_METHOD = {
    # Read depth at individual bases
    "cfdnalab": "read-depth",
    "cfdnalab-autosomes": "read-depth",
    "mosdepth": "read-depth",
    "samtools-depth": "read-depth",
    "pandepth": "read-depth",
    # Fragment coverage at individual bases, including the gap between mates
    "cfdnalab-fragment": "fragment-coverage",
    "cfdnalab-fragment-proper": "fragment-coverage",
    "mosdepth-fragment": "fragment-coverage",
    # Covered bases summed per 500 kb window
    "cfdnalab-regions": "regional-read-depth",
    "mosdepth-regions": "regional-read-depth",
    "samtools-regions": "regional-read-depth",
    "pandepth-regions": "regional-read-depth",
    # Fragments counted per 500 kb window
    "cfdnalab-region-counts": "regional-fragment-counts",
    "cfdnalab-region-counts-proper": "regional-fragment-counts",
    "finaletoolkit-coverage": "regional-fragment-counts",
}
PER_BASE_QUANTITIES = ("read-depth", "fragment-coverage")


class TrackFormat(NamedTuple):
    """How a per-base method's saved track is stored.

    Fields
    ------
    - `filename`:
        The saved file in profiles/METHOD.
    - `one_row_per_position`:
        True for chromosome/position/depth rows. False for chromosome/start/end/depth
        rows, which are 0-based half-open runs.
    - `first_position`:
        1 for 1-based positions, 0 for 0-based positions. Only used for position rows.
    - `sparse`:
        True when zero-depth positions are left out. Otherwise the track must cover
        every base of every processed chromosome, so a truncated table fails.
    """

    filename: str
    one_row_per_position: bool
    first_position: int
    sparse: bool


# The saved per-base tracks. See save_profile in coverage/commands.sh
CFDNALAB_TRACK = TrackFormat("fcoverage.per_position.bedgraph.zst", False, 0, True)
MOSDEPTH_TRACK = TrackFormat("coverage.per-base.bed.gz", False, 0, False)
TRACK_FORMAT_BY_METHOD = {
    # cfDNAlab's bedGraph: runs, with zero-depth runs left out
    "cfdnalab": CFDNALAB_TRACK,
    "cfdnalab-autosomes": CFDNALAB_TRACK,
    "cfdnalab-fragment": CFDNALAB_TRACK,
    "cfdnalab-fragment-proper": CFDNALAB_TRACK,
    # mosdepth's per-base BED: runs covering every base
    "mosdepth": MOSDEPTH_TRACK,
    "mosdepth-fragment": MOSDEPTH_TRACK,
    # samtools depth: a 1-based row per covered position, zero depths left out
    "samtools-depth": TrackFormat("positions.tsv.gz", True, 1, True),
    # PanDepth -a: a 0-based row for every position, zero depths included
    "pandepth": TrackFormat("positions.tsv.gz", True, 0, False),
}

RESIDUAL_COLUMNS = [
    "left_profile", "right_profile", "quantity", "same_native_scope", "chromosomes",
    "residual", "bases", "fraction_of_bases",
]
JACCARD_COLUMNS = [
    "left_profile", "right_profile", "quantity", "same_native_scope", "chromosomes",
    "chromosome", "chromosome_order", "bases", "left_nonzero_bases", "right_nonzero_bases",
    "nonzero_union_bases", "nonzero_intersection_bases", "coverage_jaccard",
    "possible_samtools_truncation",
]


# ---------------------------------------------------------------------------
# Reading native files
# ---------------------------------------------------------------------------


@contextmanager
def open_track(path):
    """Open a plain, gzip or zstd text file for reading.

    zstd files are read through the zstd executable, which reads every frame of a
    file made of concatenated frames. A decompression error raises an error.
    """
    if path.suffix == ".zst":
        with tempfile.TemporaryFile(mode="w+") as errors:
            process = subprocess.Popen(
                ["zstd", "-dc", str(path)],
                stdout=subprocess.PIPE,
                stderr=errors,
                text=True,
            )
            try:
                yield process.stdout
                status = process.wait()
                if status:
                    errors.seek(0)
                    raise ValueError(f"Could not decompress {path}: {errors.read()}")
            finally:
                process.stdout.close()
                if process.poll() is None:
                    process.terminate()
                process.wait()
    else:
        with gzip.open(path, "rt") if path.suffix == ".gz" else path.open() as source:
            yield source


def read_scope(directory):
    """Read the chromosomes a per-base method processed, with their lengths.

    chrom.sizes lists every processed chromosome, including those with zero
    coverage. A chromosome missing from a sparse track is therefore zero
    everywhere, while a chromosome missing from chrom.sizes was not processed.
    """
    scope = {}
    with (directory / "chrom.sizes").open() as source:
        for chromosome, length in csv.reader(source, delimiter="\t"):
            if chromosome in scope or int(length) <= 0:
                raise ValueError(
                    f"Invalid or duplicate chromosome in {directory}/chrom.sizes"
                )
            scope[chromosome] = int(length)
    if not scope:
        raise ValueError(f"Empty processed scope: {directory}")
    return scope


def spool_track(directory, scope, destination):
    """Check a native per-base track and save its runs, one file per chromosome.

    Each output file holds the chromosome's (start, end, depth) runs in 0-based
    half-open coordinates. Adjacent runs with equal depth are merged. Splitting by
    chromosome lets two tracks with different chromosome orders be read side by
    side later, without holding either in memory.

    The track must be sorted, with no overlapping runs, inside the declared scope,
    and with finite nonnegative depths. A dense track (see TrackFormat) must also
    cover every base of every chromosome in its scope.

    Returns a dictionary from chromosome to its run file. Run files of chromosomes
    absent from a sparse track are not created, and are read as zero coverage.
    """
    method = directory.name
    track_format = TRACK_FORMAT_BY_METHOD[method]
    sparse = track_format.sparse
    expected_columns = 3 if track_format.one_row_per_position else 4
    paths = {
        chromosome: destination / f"{index}.tsv"
        for index, chromosome in enumerate(scope)
    }
    previous_chromosome = None
    previous_end = 0
    visited = set()
    pending = None
    output = None
    try:
        with open_track(directory / track_format.filename) as source:
            for line_number, line in enumerate(source, start=1):
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split()
                if len(fields) != expected_columns:
                    raise ValueError(
                        f"Unexpected track columns: {method}, line {line_number}"
                    )
                chromosome = fields[0]
                if track_format.one_row_per_position:
                    # A position becomes a 1 bp run in 0-based coordinates
                    start = int(fields[1]) - track_format.first_position
                    end, depth = start + 1, float(fields[2])
                else:
                    start, end, depth = int(fields[1]), int(fields[2]), float(fields[3])
                if chromosome not in scope or not 0 <= start < end <= scope[chromosome]:
                    raise ValueError(
                        f"Out-of-scope interval: {method}, line {line_number}"
                    )
                if not math.isfinite(depth) or depth < 0:
                    raise ValueError(
                        f"Invalid coverage value: {method}, line {line_number}"
                    )
                if chromosome != previous_chromosome:
                    # A new chromosome starts. Each chromosome must appear in a single block
                    if chromosome in visited:
                        raise ValueError(
                            f"Repeated chromosome block: {method}, {chromosome}"
                        )
                    if (
                        previous_chromosome is not None
                        and not sparse
                        and previous_end != scope[previous_chromosome]
                    ):
                        raise ValueError(
                            f"Incomplete native track: {method}, {previous_chromosome}"
                        )
                    if output:
                        output.write("\t".join(map(str, pending)) + "\n")
                        output.close()
                    output = paths[chromosome].open("w")
                    visited.add(chromosome)
                    previous_chromosome, previous_end, pending = chromosome, 0, None
                if start < previous_end or (not sparse and start != previous_end):
                    raise ValueError(
                        f"Overlapping, unsorted or incomplete runs: {method}, line {line_number}"
                    )
                # Merge with the pending run when it continues it at the same depth
                if pending and pending[1] == start and pending[2] == depth:
                    pending = (pending[0], end, depth)
                else:
                    if pending:
                        output.write("\t".join(map(str, pending)) + "\n")
                    pending = (start, end, depth)
                previous_end = end
        if not sparse and (
            visited != set(scope) or previous_end != scope[previous_chromosome]
        ):
            raise ValueError(f"Incomplete native track: {method}")
        if output:
            output.write("\t".join(map(str, pending)) + "\n")
    finally:
        if output:
            output.close()
    return paths


def zero_filled_runs(path, chromosome_length):
    """Yield a chromosome's runs, with zero-depth runs filling the gaps.

    The runs then cover [0, chromosome_length) without gaps, in half-open
    coordinates. A missing run file means zero coverage on the whole chromosome.
    """
    position = 0
    if path.exists():
        with path.open() as source:
            for fields in csv.reader(source, delimiter="\t"):
                start, end, depth = int(fields[0]), int(fields[1]), float(fields[2])
                if start > position:
                    yield position, start, 0.0
                yield start, end, depth
                position = end
    if position < chromosome_length:
        yield position, chromosome_length, 0.0


def read_bed_rows(path):
    """Return the tab-separated fields of every row, skipping empty and '#' rows."""
    with open_track(path) as source:
        return [
            fields
            for fields in csv.reader(source, delimiter="\t")
            if fields and not fields[0].startswith("#")
        ]


def read_cfdnalab_windows(path, value_column):
    """Read cfDNAlab's per-window table, which has a header row."""
    with open_track(path) as source:
        return [
            (row["chromosome"], row["start"], row["end"], row[value_column])
            for row in csv.DictReader(source, delimiter="\t")
        ]


def read_regional_table(directory):
    """Read a regional method's native table as (chromosome, start, end, value) rows.

    The value is the covered-base total for regional read depth, and the fragment
    count for regional fragment counts. Coordinates are returned as 0-based
    half-open BED coordinates.
    """
    method = directory.name
    if method == "cfdnalab-regions":
        # total_coverage sums the depth of every base in the window, which is the
        # covered-base total
        return read_cfdnalab_windows(directory / "fcoverage.total.tsv.zst", "total_coverage")
    if method in ("cfdnalab-region-counts", "cfdnalab-region-counts-proper"):
        # With --normalize-by-length, each fragment adds a total of 1 across the
        # windows it overlaps, so total_fragment_mass estimates the fragment count
        return read_cfdnalab_windows(
            directory / "length_normalized.fcoverage.total.tsv.zst", "total_fragment_mass"
        )
    if method == "mosdepth-regions":
        # The last column is the mean depth in the window. Mean depth times the
        # window's width gives the covered-base total
        return [
            (fields[0], fields[1], fields[2], float(fields[-1]) * (int(fields[2]) - int(fields[1])))
            for fields in read_bed_rows(directory / "coverage.regions.bed.gz")
        ]
    if method == "samtools-regions":
        # bedcov repeats the BED columns and adds the covered-base total last
        return [fields[:3] + fields[-1:] for fields in read_bed_rows(directory / "samtools.bedcov.tsv")]
    if method == "pandepth-regions":
        # The v2.26 -w writer uses inclusive 1-based Start and End, then Length,
        # CoveredSite and TotalDepth. Its file has no BED RegionID column. Start - 1
        # gives the BED start, and the inclusive End equals the BED end
        return [
            (fields[0], int(fields[1]) - 1, fields[2], fields[5])
            for fields in read_bed_rows(directory / "pandepth.win.stat.gz")
        ]
    if method == "finaletoolkit-coverage":
        # FinaleToolkit repeats the BED columns and adds the fragment count last
        return [
            fields[:3] + fields[-1:]
            for fields in read_bed_rows(directory / "finaletoolkit.coverage.bed.gz")
        ]
    raise ValueError(f"Unknown regional coverage method: {method}")


def regional_counts(directory, expected_regions):
    """Read a regional method's values, keyed by (chromosome, start, end).

    Every window must be one of the shared BED's windows and appear once. Values
    must be finite and nonnegative, and FinaleToolkit's fragment counts must be
    whole numbers. Fractional values from other methods are kept as they are.

    PanDepth's window loop (j < chromosome_length) leaves out a final window of
    exactly 1 bp. That is the only missing window accepted, and the comparison
    reports it as a missing region. Any other missing or unexpected window fails.
    """
    counts = {}
    for chromosome, start, end, count in read_regional_table(directory):
        key = str(chromosome), int(start), int(end)
        count = float(count)
        if (
            key in counts
            or not math.isfinite(count)
            or count < 0
            or (directory.name == "finaletoolkit-coverage" and not count.is_integer())
        ):
            raise ValueError(f"Duplicate region or invalid count in {directory}")
        counts[key] = count
    missing = expected_regions - set(counts)
    allowed_missing = {
        key for key in expected_regions
        if key[2] - key[1] == 1
        and key[2] == max(region[2] for region in expected_regions if region[0] == key[0])
    } if directory.name == "pandepth-regions" else set()
    if set(counts) - expected_regions or not missing <= allowed_missing:
        raise ValueError(
            f"Regional row identity differs from the shared BED: {directory}"
        )
    return counts


def read_shared_regions(path):
    """Read the shared regions.bed as a list of unique (chromosome, start, end) windows."""
    with path.open() as source:
        regions = [
            (fields[0], int(fields[1]), int(fields[2]))
            for fields in csv.reader(source, delimiter="\t")
        ]
    if not regions or len(set(regions)) != len(regions):
        raise ValueError("Shared BED must contain unique regions")
    return regions


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


class Correlation:
    """Accumulate a weighted Pearson correlation one value pair at a time.

    The running means, variances and covariance are updated with Welford's method,
    which stays numerically stable for long tracks.
    """

    def __init__(self):
        self.count = 0
        self.left_mean = self.right_mean = 0.0
        self.left_variance = self.right_variance = self.covariance = 0.0

    def add(self, left, right, weight=1):
        """Add a value pair. Per-base runs are weighted by bases, windows by 1."""
        self.count += weight
        left_delta, right_delta = left - self.left_mean, right - self.right_mean
        relative_weight = weight / self.count
        self.left_mean += left_delta * relative_weight
        self.right_mean += right_delta * relative_weight
        self.left_variance += weight * left_delta * (left - self.left_mean)
        self.right_variance += weight * right_delta * (right - self.right_mean)
        self.covariance += weight * left_delta * (right - self.right_mean)

    def value(self):
        """Return Pearson r, or NA for a constant profile (including all zeros)."""
        # The shared check rejects nonfinite arithmetic and excessive overshoots
        return correlation_value(
            self.covariance, self.left_variance, self.right_variance
        )


class ResidualCounts:
    """Count the bases with each residual value (left - right), with bounded memory.

    Counts are collected in memory in batches of up to 4096 distinct residuals and
    then added to a temporary SQLite table on disk. Integer depths give few distinct
    residuals, but fractional depths could give one per run, so the disk table keeps
    memory bounded either way. The SQLite cache is limited as well.
    """

    def __init__(self, path):
        self.database = sqlite3.connect(path)
        self.database.execute("PRAGMA cache_size = -2048")
        self.database.execute(
            "CREATE TABLE residual_counts (residual REAL PRIMARY KEY, bases INTEGER NOT NULL)"
        )
        self.batch = Counter()

    def add(self, residual, bases):
        """Count `bases` more bases with this residual."""
        self.batch[residual] += bases
        if len(self.batch) >= 4096:
            self.flush()

    def flush(self):
        """Add the in-memory batch to the table on disk."""
        self.database.executemany(
            "INSERT INTO residual_counts (residual, bases) VALUES (?, ?) "
            "ON CONFLICT(residual) DO UPDATE SET bases = residual_counts.bases + excluded.bases",
            self.batch.items(),
        )
        self.batch.clear()

    def sorted_counts(self):
        """Yield (residual, bases) pairs in increasing residual order."""
        self.flush()
        yield from self.database.execute(
            "SELECT residual, bases FROM residual_counts ORDER BY residual"
        )

    def clear(self):
        """Remove all counts, before the next pair."""
        self.database.execute("DELETE FROM residual_counts")
        self.database.commit()

    def close(self):
        """Close the database."""
        self.database.close()


def compare_tracks(left_paths, right_paths, scope, residual_counts):
    """Compare two per-base tracks at base resolution over the shared scope.

    The two tracks' runs are walked side by side. Each step covers a stretch of
    bases where both depths are constant, and counts with its number of bases
    (`weight`). Every base in the scope is compared, including bases where both
    depths are zero. For example, MAE = sum(|left - right| * width) / bases.

    The "union" metrics only count bases where either track has nonzero depth.
    Coverage Jaccard is (bases where both are nonzero) / (bases where either is).
    The residual of every base is added to `residual_counts`.

    Returns the summary metrics, plus per-chromosome overlap in "chromosome_metrics".
    """
    bases = identical = union_bases = union_identical = intersection_bases = 0
    below = above = 0
    absolute_error = squared_error = signed_error = left_total = right_total = 0.0
    correlation = Correlation()
    chromosome_metrics = []
    for chromosome, length in scope.items():
        chromosome_union = chromosome_intersection = 0
        left_nonzero = right_nonzero = 0
        left_runs = iter(zero_filled_runs(left_paths[chromosome], length))
        right_runs = iter(zero_filled_runs(right_paths[chromosome], length))
        left_run, right_run = next(left_runs), next(right_runs)
        position = 0
        while position < length:
            # The step ends where the first of the two current runs ends
            end = min(left_run[1], right_run[1], length)
            weight = end - position
            left, right = left_run[2], right_run[2]
            difference = left - right
            bases += weight
            identical += weight if difference == 0 else 0
            below += weight if difference < 0 else 0
            above += weight if difference > 0 else 0
            if left != 0 or right != 0:
                union_bases += weight
                chromosome_union += weight
                union_identical += weight if difference == 0 else 0
            if left != 0 and right != 0:
                intersection_bases += weight
                chromosome_intersection += weight
            left_nonzero += weight if left != 0 else 0
            right_nonzero += weight if right != 0 else 0
            absolute_error += abs(difference) * weight
            squared_error += difference**2 * weight
            signed_error += difference * weight
            left_total += left * weight
            right_total += right * weight
            # Base weights give Pearson over individual bases without binning
            correlation.add(left, right, weight)
            residual_counts.add(difference, weight)
            position = end
            if end == left_run[1] and end < length:
                left_run = next(left_runs)
            if end == right_run[1] and end < length:
                right_run = next(right_runs)
        chromosome_metrics.append({
            "chromosome": chromosome,
            "bases": length,
            "left_nonzero_bases": left_nonzero,
            "right_nonzero_bases": right_nonzero,
            "nonzero_union_bases": chromosome_union,
            "nonzero_intersection_bases": chromosome_intersection,
            "coverage_jaccard": checked_metric(
                chromosome_intersection / chromosome_union, "chromosome coverage Jaccard", maximum=1.0
            ) if chromosome_union else "NA",
        })
    return {
        "bases": bases,
        "fraction_identical": checked_metric(
            identical / bases, "identical fraction", maximum=1.0
        ),
        "mae": checked_metric(absolute_error / bases, "mae"),
        "rmse": checked_metric(math.sqrt(squared_error / bases), "rmse"),
        "mean_signed_difference": checked_metric(
            signed_error / bases, "mean signed difference", minimum=None
        ),
        "fraction_left_below": checked_metric(
            below / bases, "left-below fraction", maximum=1.0
        ),
        "fraction_left_above": checked_metric(
            above / bases, "left-above fraction", maximum=1.0
        ),
        "nonzero_union_bases": union_bases,
        "nonzero_intersection_bases": intersection_bases,
        "coverage_jaccard": checked_metric(
            intersection_bases / union_bases, "coverage Jaccard", maximum=1.0
        ) if union_bases else "NA",
        "chromosome_metrics": chromosome_metrics,
        "union_fraction_identical": checked_metric(
            union_identical / union_bases, "union identical fraction", maximum=1.0
        )
        if union_bases
        else "NA",
        "union_mae": checked_metric(absolute_error / union_bases, "union MAE")
        if union_bases
        else "NA",
        "union_rmse": checked_metric(
            math.sqrt(squared_error / union_bases), "union RMSE"
        )
        if union_bases
        else "NA",
        "left_integrated_depth": checked_metric(left_total, "left integrated depth"),
        "right_integrated_depth": checked_metric(right_total, "right integrated depth"),
        "pearson_r_1bp": correlation.value(),
    }


def compare_regions(left, right):
    """Compare two methods' values in the same windows, without rounding fractions.

    Each window counts once. Returns the window count, both totals and their ratio,
    and the per-window errors and Pearson r.
    """
    differences = [left[key] - right[key] for key in left]
    correlation = Correlation()
    for key in left:
        correlation.add(left[key], right[key])
    left_total = checked_metric(sum(left.values()), "left summed region counts")
    right_total = checked_metric(sum(right.values()), "right summed region counts")
    return {
        "regions": len(left),
        "left_summed_region_counts": left_total,
        "right_summed_region_counts": right_total,
        "regional_total_ratio": checked_metric(
            left_total / right_total, "regional total ratio"
        )
        if right_total
        else "NA",
        "net_count_difference": checked_metric(
            sum(differences), "net count difference", minimum=None
        ),
        "mae": checked_metric(
            sum(abs(value) for value in differences) / len(left), "regional MAE"
        ),
        "rmse": checked_metric(
            math.sqrt(sum(value**2 for value in differences) / len(left)), "regional RMSE"
        ),
        "fraction_identical": checked_metric(
            sum(value == 0 for value in differences) / len(left),
            "identical region fraction", maximum=1.0,
        ),
        "pearson_r_regions": correlation.value(),
    }


def possible_samtools_truncation(left_name, right_name, chromosomes):
    """Flag partial overlap followed by a zero-coverage suffix in BAM order.

    samtools depth could stop early without an error, leaving later chromosomes
    empty. This looks for a chromosome with partial overlap, after which samtools
    reports nothing while the other method still has coverage. It is a diagnostic
    pattern, not a failure criterion. Chromosomes with no coverage from either tool
    cannot tell truncation apart from real zeros, so they are ignored.

    Returns the flagged chromosomes, or an empty set.
    """
    if "samtools-depth" not in (left_name, right_name):
        return set()
    samtools_column = "left_nonzero_bases" if left_name == "samtools-depth" else "right_nonzero_bases"
    other_column = "right_nonzero_bases" if left_name == "samtools-depth" else "left_nonzero_bases"
    informative = [row for row in chromosomes if row["nonzero_union_bases"]]
    for index, row in enumerate(informative[:-1]):
        later = informative[index + 1:]
        if (
            0 < row["coverage_jaccard"] < 1
            and row[samtools_column] > 0
            and all(following[samtools_column] == 0 and following[other_column] > 0 for following in later)
        ):
            return {entry["chromosome"] for entry in informative[index:]}
    return set()


# ---------------------------------------------------------------------------
# The comparison steps
# ---------------------------------------------------------------------------


def plan_per_base_comparisons(scopes, allow_common_scope):
    """Pair per-base profiles of the same quantity and choose each pair's chromosomes.

    `scopes` maps each per-base method to its processed chromosomes (read_scope).
    All pairs are checked here, before any large track is read. Two profiles that
    processed different chromosomes are compared on the chromosomes they share only
    with --common-scope. The job passes it, so the optional autosome configuration
    can be compared with the all-contig profiles. Shared chromosomes must have the
    same length in both profiles.

    Returns a list of (left method, right method, shared scope, whether the two
    native scopes are identical).
    """
    comparisons = []
    for left_name, right_name in combinations(scopes, 2):
        if QUANTITY_BY_METHOD[left_name] != QUANTITY_BY_METHOD[right_name]:
            continue
        left_scope, right_scope = scopes[left_name], scopes[right_name]
        same_scope = left_scope == right_scope
        if not same_scope and not allow_common_scope:
            raise ValueError(
                f"Scopes differ: {left_name}, {right_name}. Use --common-scope to compare their intersection"
            )
        scope = {
            chromosome: length
            for chromosome, length in left_scope.items()
            if chromosome in right_scope
        }
        if not scope or any(
            right_scope[chromosome] != length for chromosome, length in scope.items()
        ):
            raise ValueError("Empty shared scope or mismatched chromosome lengths")
        comparisons.append((left_name, right_name, scope, same_scope))
    return comparisons


def compare_per_base_pairs(comparisons, tracks, residual_counts, residual_writer):
    """Compare every planned per-base pair and write its residual histogram.

    Returns the summary rows and the per-chromosome overlap rows.
    """
    rows = []
    chromosome_rows = []
    for left_name, right_name, scope, same_scope in comparisons:
        comparison = {
            "left_profile": left_name,
            "right_profile": right_name,
            "quantity": QUANTITY_BY_METHOD[left_name],
            "same_native_scope": same_scope,
            "chromosomes": ",".join(scope),
        }
        metrics = compare_tracks(
            tracks[left_name], tracks[right_name], scope, residual_counts
        )
        chromosomes = metrics.pop("chromosome_metrics")
        flagged = possible_samtools_truncation(left_name, right_name, chromosomes)
        if flagged:
            print(
                f"Possible samtools-depth truncation in {left_name} versus {right_name}. "
                "Partial chromosome overlap is followed by zero samtools coverage on later "
                "covered chromosomes in BAM order. Inspect coverage_jaccard.tsv and the run logs.",
                file=sys.stderr,
            )
        chromosome_rows.extend(
            {**comparison, **row, "chromosome_order": index,
             "possible_samtools_truncation": row["chromosome"] in flagged}
            for index, row in enumerate(chromosomes, 1)
        )
        rows.append({**comparison, **metrics})
        for residual, bases in residual_counts.sorted_counts():
            residual_writer.writerow(
                {
                    **comparison,
                    "residual": residual,
                    "bases": bases,
                    "fraction_of_bases": checked_metric(
                        bases / metrics["bases"], "residual fraction", maximum=1.0
                    ),
                }
            )
        residual_counts.clear()
    return rows, chromosome_rows


def compare_regional_pairs(regions, shared_region_count):
    """Compare every pair of regional profiles with the same quantity.

    Pairs are compared on the windows both report. All methods report every
    shared window, except PanDepth's missing final 1 bp window (see
    regional_counts). The missing windows are counted for each side.
    """
    rows = []
    for left_name, right_name in combinations(regions, 2):
        quantity = QUANTITY_BY_METHOD[left_name]
        if quantity != QUANTITY_BY_METHOD[right_name]:
            continue
        common_regions = regions[left_name].keys() & regions[right_name].keys()
        if not common_regions:
            raise ValueError(f"No shared regional windows: {left_name}, {right_name}")
        left_counts = {key: value for key, value in regions[left_name].items() if key in common_regions}
        right_counts = {key: regions[right_name][key] for key in left_counts}
        metrics = compare_regions(left_counts, right_counts)
        metrics["left_missing_regions"] = shared_region_count - len(regions[left_name])
        metrics["right_missing_regions"] = shared_region_count - len(regions[right_name])
        if quantity == "regional-read-depth":
            # These totals count covered bases rather than fragments
            metrics["left_integrated_depth"] = metrics.pop("left_summed_region_counts")
            metrics["right_integrated_depth"] = metrics.pop("right_summed_region_counts")
            metrics["net_coverage_difference"] = metrics.pop("net_count_difference")
        rows.append(
            {
                "left_profile": left_name,
                "right_profile": right_name,
                "quantity": quantity,
                "same_native_scope": regions[left_name].keys() == regions[right_name].keys(),
                "chromosomes": ",".join(
                    sorted({key[0] for key in common_regions})
                ),
                **metrics,
            }
        )
    return rows


@contextmanager
def replace_when_complete(path):
    """Write to a temporary file beside `path`, and move it into place on success.

    If the block raises an error, the temporary file is removed and an earlier
    complete file at `path` stays unchanged.
    """
    pending = tempfile.NamedTemporaryFile(
        mode="w", dir=path.parent, prefix=f".{path.stem}-", suffix=path.suffix, delete=False
    )
    pending_path = Path(pending.name)
    try:
        with pending:
            yield pending
    except BaseException:
        pending_path.unlink(missing_ok=True)
        raise
    os.replace(pending_path, path)


def write_table(destination, rows, columns):
    """Write rows as a tab-separated table, with NA for missing values."""
    writer = csv.DictWriter(
        destination, fieldnames=columns, delimiter="\t", restval="NA", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)


def main():
    """Compare like quantities, recording the evaluated chromosome scope per row."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", type=Path)
    parser.add_argument(
        "--common-scope",
        action="store_true",
        help="Explicitly allow intersection of declared scopes",
    )
    arguments = parser.parse_args()
    profiles = arguments.profiles
    results = profiles.parent

    directories = sorted(path for path in profiles.iterdir() if path.is_dir())
    for directory in directories:
        if directory.name not in QUANTITY_BY_METHOD:
            raise ValueError(f"Unknown coverage method: {directory.name}")
    per_base_directories = [
        path for path in directories if QUANTITY_BY_METHOD[path.name] in PER_BASE_QUANTITIES
    ]
    regional_directories = [path for path in directories if path not in per_base_directories]
    # All pair scopes are checked before large tracks are read or results replaced
    scopes = {path.name: read_scope(path) for path in per_base_directories}
    comparisons = plan_per_base_comparisons(scopes, arguments.common_scope)

    with (
        tempfile.TemporaryDirectory(prefix="coverage-concordance-") as temporary,
        replace_when_complete(results / "residuals.tsv") as residual_output,
        replace_when_complete(results / "coverage_jaccard.tsv") as jaccard_output,
        closing(ResidualCounts(Path(temporary) / "residuals.sqlite")) as residual_counts,
    ):
        residual_writer = csv.DictWriter(
            residual_output, fieldnames=RESIDUAL_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        residual_writer.writeheader()

        # Regional profiles are small and are read completely
        regions = {}
        shared_region_count = 0
        if regional_directories:
            shared_regions = read_shared_regions(results / "regions.bed")
            shared_region_count = len(shared_regions)
            regions = {
                path.name: regional_counts(path, set(shared_regions))
                for path in regional_directories
            }

        # Per-base tracks are checked and split into per-chromosome run files
        tracks = {}
        for directory in per_base_directories:
            destination = Path(temporary) / directory.name
            destination.mkdir()
            tracks[directory.name] = spool_track(directory, scopes[directory.name], destination)

        rows, chromosome_rows = compare_per_base_pairs(
            comparisons, tracks, residual_counts, residual_writer
        )
        rows += compare_regional_pairs(regions, shared_region_count)
        if not rows:
            parser.error("No pair of saved profiles measures the same quantity")
        # This small table records overlap per chromosome in BAM order. It is
        # separate from the exact residual histogram, which stays base-weighted
        write_table(jaccard_output, chromosome_rows, JACCARD_COLUMNS)

    # Different quantities have different metrics, so the columns are the union of
    # all rows' metrics, in the order they first appear
    columns = list(dict.fromkeys(column for row in rows for column in row))
    write_table(sys.stdout, rows, columns)


if __name__ == "__main__":
    main()
