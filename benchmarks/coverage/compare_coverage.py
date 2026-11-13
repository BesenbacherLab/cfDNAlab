#!/usr/bin/env python3
"""Compare retained coverage tracks without expanding them into genomic arrays.

Run boundaries carry base weights. Temporary per-contig run files reconcile
different native contig orders and keep memory bounded. Regional read depth and
fragment counts are compared separately at exact shared BED coordinates. Requires zstd for
cfDNAlab tracks and length-normalized regional coverage tables.
Always save exact base-weighted residual counts to residuals.tsv beside profiles/.
Summary metrics go to stdout. Base-resolution MAE is the primary similarity measure.
"""

import argparse
from collections import Counter
from contextlib import ExitStack, closing, contextmanager
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

# The shared checks live in benchmarks/, above this feature's script directory
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_metrics import checked_metric, correlation_value


REGIONAL_COVERAGE_METHODS = {
    "cfdnalab-regions", "mosdepth-regions", "samtools-regions", "pandepth-regions"
}


@contextmanager
def open_track(path):
    """Read all compressed frames and propagate decompression errors."""
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
    """Require explicit processed contigs, including those with zero coverage."""
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
    """Validate native runs, then coalesce equal adjacent values on disk.

    cfDNAlab and samtools depth may omit zero positions. Other native formats
    must cover every base of their declared scope, so truncated tables fail.
    """
    method = directory.name
    positional = method in ("samtools-depth", "pandepth")
    sparse = method.startswith("cfdnalab") or method == "samtools-depth"
    filename = (
        "positions.tsv.gz"
        if positional
        else (
            "fcoverage.per_position.bedgraph.zst"
            if sparse
            else "coverage.per-base.bed.gz"
        )
    )
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
        with open_track(directory / filename) as source:
            for line_number, line in enumerate(source, start=1):
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split()
                if len(fields) != (3 if positional else 4):
                    raise ValueError(
                        f"Unexpected track columns: {method}, line {line_number}"
                    )
                chromosome = fields[0]
                if positional:
                    start = int(fields[1]) - (1 if method == "samtools-depth" else 0)
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
    """Fill only declared sparse gaps, preserving half-open coordinates."""
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


class Correlation:
    """Accumulate weighted paired values with stable covariance and variance."""

    def __init__(self):
        self.count = 0
        self.left_mean = self.right_mean = 0.0
        self.left_variance = self.right_variance = self.covariance = 0.0

    def add(self, left, right, weight=1):
        """Weight constant depth runs by bases, or regional counts equally."""
        self.count += weight
        left_delta, right_delta = left - self.left_mean, right - self.right_mean
        relative_weight = weight / self.count
        self.left_mean += left_delta * relative_weight
        self.right_mean += right_delta * relative_weight
        self.left_variance += weight * left_delta * (left - self.left_mean)
        self.right_variance += weight * right_delta * (right - self.right_mean)
        self.covariance += weight * left_delta * (right - self.right_mean)

    def value(self):
        """Constant profiles have undefined correlation, including equal zeros."""
        # The shared check rejects nonfinite arithmetic and excessive overshoots
        return correlation_value(
            self.covariance, self.left_variance, self.right_variance
        )


def flush_residual_counts(database, counts):
    """Merge a bounded batch into exact residual counts stored on disk."""
    database.executemany(
        "INSERT INTO residual_counts (residual, bases) VALUES (?, ?) "
        "ON CONFLICT(residual) DO UPDATE SET bases = residual_counts.bases + excluded.bases",
        counts.items(),
    )
    counts.clear()


def compare_tracks(left_paths, right_paths, scope, residual_database):
    """Integrate errors, Pearson and residual counts at base resolution."""
    bases = identical = union_bases = union_identical = intersection_bases = 0
    below = above = 0
    absolute_error = squared_error = signed_error = left_total = right_total = 0.0
    correlation = Correlation()
    residual_counts = Counter()
    chromosome_metrics = []
    for chromosome, length in scope.items():
        chromosome_union = chromosome_intersection = 0
        left_nonzero = right_nonzero = 0
        left_runs = iter(zero_filled_runs(left_paths[chromosome], length))
        right_runs = iter(zero_filled_runs(right_paths[chromosome], length))
        left_run, right_run = next(left_runs), next(right_runs)
        position = 0
        while position < length:
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
            # Length weights give Pearson over individual bases without binning
            correlation.add(left, right, weight)
            residual_counts[difference] += weight
            if len(residual_counts) >= 4096:
                flush_residual_counts(residual_database, residual_counts)
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
    flush_residual_counts(residual_database, residual_counts)
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


def regional_counts(directory, expected_regions):
    """Read fragment or covered-base counts at exact shared coordinates.

    mosdepth reports mean depth, converted to covered-base totals by window width.
    PanDepth's native inclusive coordinates convert back to BED coordinates.
    """
    is_cfdnalab = directory.name.startswith("cfdnalab")
    if is_cfdnalab:
        length_normalized = "region-counts" in directory.name
        filename = (
            "length_normalized.fcoverage.total.tsv.zst"
            if length_normalized else "fcoverage.total.tsv.zst"
        )
        value_column = "total_fragment_mass" if length_normalized else "total_coverage"
        with open_track(directory / filename) as source:
            entries = [
                (row["chromosome"], row["start"], row["end"], row[value_column])
                for row in csv.DictReader(source, delimiter="\t")
            ]
    elif directory.name == "mosdepth-regions":
        with open_track(directory / "coverage.regions.bed.gz") as source:
            entries = [
                (
                    fields[0], fields[1], fields[2],
                    float(fields[-1]) * (int(fields[2]) - int(fields[1])),
                )
                for fields in csv.reader(source, delimiter="\t")
                if fields and not fields[0].startswith("#")
            ]
    elif directory.name == "samtools-regions":
        with open_track(directory / "samtools.bedcov.tsv") as source:
            entries = [
                fields[:3] + fields[-1:]
                for fields in csv.reader(source, delimiter="\t")
                if fields and not fields[0].startswith("#")
            ]
    elif directory.name == "pandepth-regions":
        # The v2.26 -w writer uses inclusive 1-based Start/End and puts TotalDepth
        # after Length and CoveredSite. Its file has no BED RegionID column
        with open_track(directory / "pandepth.win.stat.gz") as source:
            entries = [
                (fields[0], int(fields[1]) - 1, fields[2], fields[5])
                for fields in csv.reader(source, delimiter="\t")
                if fields and not fields[0].startswith("#")
            ]
    else:
        with open_track(directory / "finaletoolkit.coverage.bed.gz") as source:
            entries = [
                fields[:3] + fields[-1:]
                for fields in csv.reader(source, delimiter="\t")
                if fields and not fields[0].startswith("#")
            ]
    counts = {}
    for chromosome, start, end, count in entries:
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
    # PanDepth's j < chromosome_length loop omits a final window of exactly 1 bp
    # No other missing windows or unexpected coordinates are silently accepted
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


def compare_regions(left, right):
    """Compare regional totals without rounding fractional values."""
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

    This is a diagnostic pattern, not a failure criterion. Chromosomes with no
    coverage from either tool cannot distinguish truncation from real zeros.
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
    directories = sorted(path for path in arguments.profiles.iterdir() if path.is_dir())
    regional_directories = [
        path
        for path in directories
        if (
            "region-counts" in path.name or path.name == "finaletoolkit-coverage"
            or path.name in REGIONAL_COVERAGE_METHODS
        )
    ]
    track_scopes = {
        path.name: read_scope(path)
        for path in directories
        if path not in regional_directories
    }
    comparisons = []
    # Check all pair scopes before loading large tracks or replacing any results
    for left_name, right_name in combinations(track_scopes, 2):
        if ("fragment" in left_name) != ("fragment" in right_name):
            continue
        left_scope, right_scope = track_scopes[left_name], track_scopes[right_name]
        same_scope = left_scope == right_scope
        if not same_scope and not arguments.common_scope:
            parser.error(
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
    rows = []
    chromosome_rows = []
    with (
        tempfile.TemporaryDirectory(prefix="coverage-concordance-") as temporary,
        ExitStack() as resources,
    ):
        tracks, regions = {}, {}
        # Publish only a completed histogram, keeping previous successful results on error
        pending_output = tempfile.NamedTemporaryFile(
            mode="w",
            dir=arguments.profiles.parent,
            prefix=".residuals-",
            suffix=".tsv",
            delete=False,
        )
        pending_path = Path(pending_output.name)
        resources.callback(pending_path.unlink, missing_ok=True)
        output = resources.enter_context(pending_output)
        residual_writer = csv.DictWriter(
            output,
            delimiter="\t",
            lineterminator="\n",
            fieldnames=[
                "left_profile",
                "right_profile",
                "quantity",
                "same_native_scope",
                "chromosomes",
                "residual",
                "bases",
                "fraction_of_bases",
            ],
        )
        residual_writer.writeheader()
        residual_database = resources.enter_context(
            closing(sqlite3.connect(Path(temporary) / "residuals.sqlite"))
        )
        # Bound the database cache as well as the Python residual-count batch
        residual_database.execute("PRAGMA cache_size = -2048")
        residual_database.execute(
            "CREATE TABLE residual_counts (residual REAL PRIMARY KEY, bases INTEGER NOT NULL)"
        )
        if regional_directories:
            with (arguments.profiles.parent / "regions.bed").open() as source:
                expected = [
                    (fields[0], int(fields[1]), int(fields[2]))
                    for fields in csv.reader(source, delimiter="\t")
                ]
            if not expected or len(set(expected)) != len(expected):
                raise ValueError("Shared BED must contain unique regions")
            regions = {
                path.name: regional_counts(path, set(expected))
                for path in regional_directories
            }
        for directory in directories:
            if directory in regional_directories:
                continue
            scope = track_scopes[directory.name]
            destination = Path(temporary) / directory.name
            destination.mkdir()
            tracks[directory.name] = spool_track(directory, scope, destination)
        for left_name, right_name, scope, same_scope in comparisons:
            comparison = {
                "left_profile": left_name,
                "right_profile": right_name,
                "quantity": "fragment-coverage"
                if "fragment" in left_name
                else "read-depth",
                "same_native_scope": same_scope,
                "chromosomes": ",".join(scope),
            }
            metrics = compare_tracks(
                tracks[left_name], tracks[right_name], scope, residual_database
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
            for residual, bases in residual_database.execute(
                "SELECT residual, bases FROM residual_counts ORDER BY residual"
            ):
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
            residual_database.execute("DELETE FROM residual_counts")
            residual_database.commit()
        for left_name, right_name in combinations(regions, 2):
            covered_bases = left_name in REGIONAL_COVERAGE_METHODS
            if covered_bases != (right_name in REGIONAL_COVERAGE_METHODS):
                continue
            common_regions = regions[left_name].keys() & regions[right_name].keys()
            if not common_regions:
                raise ValueError(f"No shared regional windows: {left_name}, {right_name}")
            left_counts = {key: value for key, value in regions[left_name].items() if key in common_regions}
            right_counts = {key: regions[right_name][key] for key in left_counts}
            metrics = compare_regions(left_counts, right_counts)
            metrics["left_missing_regions"] = len(expected) - len(regions[left_name])
            metrics["right_missing_regions"] = len(expected) - len(regions[right_name])
            if covered_bases:
                # These totals count covered bases rather than fragments
                metrics["left_integrated_depth"] = metrics.pop("left_summed_region_counts")
                metrics["right_integrated_depth"] = metrics.pop("right_summed_region_counts")
                metrics["net_coverage_difference"] = metrics.pop("net_count_difference")
            rows.append(
                {
                    "left_profile": left_name,
                    "right_profile": right_name,
                    "quantity": "regional-read-depth" if covered_bases else "regional-fragment-counts",
                    "same_native_scope": regions[left_name].keys() == regions[right_name].keys(),
                    "chromosomes": ",".join(
                        sorted({key[0] for key in common_regions})
                    ),
                    **metrics,
                }
            )
        if not rows:
            parser.error("No pair of saved profiles measures the same quantity")
        # This small table records overlap per chromosome in BAM order. It is
        # separate from the exact residual histogram, which remains base-weighted
        with tempfile.NamedTemporaryFile(
            mode="w", dir=arguments.profiles.parent, prefix=".coverage-jaccard-", delete=False
        ) as diagnostics:
            diagnostics_path = Path(diagnostics.name)
            resources.callback(diagnostics_path.unlink, missing_ok=True)
            diagnostic_writer = csv.DictWriter(
                diagnostics, delimiter="\t", lineterminator="\n",
                fieldnames=["left_profile", "right_profile", "quantity", "same_native_scope",
                            "chromosomes", "chromosome", "chromosome_order", "bases",
                            "left_nonzero_bases", "right_nonzero_bases", "nonzero_union_bases",
                            "nonzero_intersection_bases", "coverage_jaccard", "possible_samtools_truncation"],
            )
            diagnostic_writer.writeheader()
            diagnostic_writer.writerows(chromosome_rows)
        os.replace(diagnostics_path, arguments.profiles.parent / "coverage_jaccard.tsv")
        output.close()
        os.replace(pending_path, arguments.profiles.parent / "residuals.tsv")
    columns = list(dict.fromkeys(column for row in rows for column in row))
    writer = csv.DictWriter(
        sys.stdout,
        fieldnames=columns,
        delimiter="\t",
        restval="NA",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
