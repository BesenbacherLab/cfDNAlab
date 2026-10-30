#!/usr/bin/env python3
"""Compare retained coverage tracks without expanding them into genomic arrays.

Run boundaries carry base weights. Temporary per-contig run files reconcile
different native contig orders and keep memory bounded. Regional counts are
compared separately by their exact shared BED coordinates. Requires zstd for
cfDNAlab tracks and the comparison environment for regional length tables.
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
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile


@contextmanager
def open_track(path):
    """Read all compressed frames and propagate decompression errors."""
    if path.suffix == ".zst":
        with tempfile.TemporaryFile(mode="w+") as errors:
            process = subprocess.Popen(["zstd", "-dc", str(path)], stdout=subprocess.PIPE,
                                       stderr=errors, text=True)
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
        with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as source:
            yield source


def read_scope(directory):
    """Require explicit processed contigs, including those with zero coverage."""
    scope = {}
    with (directory / "chrom.sizes").open() as source:
        for chromosome, length in csv.reader(source, delimiter="\t"):
            if chromosome in scope or int(length) <= 0:
                raise ValueError(f"Invalid or duplicate chromosome in {directory}/chrom.sizes")
            scope[chromosome] = int(length)
    if not scope:
        raise ValueError(f"Empty processed scope: {directory}")
    return scope


def spool_track(directory, scope, destination):
    """Validate native runs, then coalesce equal adjacent values on disk.

    Only cfDNAlab's sparse format may omit zero positions. Other native formats
    must cover every base of their declared scope, so truncated tables fail.
    """
    method = directory.name
    positional = method in ("samtools-depth", "pandepth")
    sparse = method.startswith("cfdnalab")
    filename = "positions.tsv.gz" if positional else (
        "fcoverage.per_position.bedgraph.zst" if sparse else "coverage.per-base.bed.gz")
    paths = {chromosome: destination / f"{index}.tsv" for index, chromosome in enumerate(scope)}
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
                    raise ValueError(f"Unexpected track columns: {method}, line {line_number}")
                chromosome = fields[0]
                if positional:
                    start = int(fields[1]) - (1 if method == "samtools-depth" else 0)
                    end, depth = start + 1, float(fields[2])
                else:
                    start, end, depth = int(fields[1]), int(fields[2]), float(fields[3])
                if chromosome not in scope or not 0 <= start < end <= scope[chromosome]:
                    raise ValueError(f"Out-of-scope interval: {method}, line {line_number}")
                if not math.isfinite(depth) or depth < 0:
                    raise ValueError(f"Invalid coverage value: {method}, line {line_number}")
                if chromosome != previous_chromosome:
                    if chromosome in visited:
                        raise ValueError(f"Repeated chromosome block: {method}, {chromosome}")
                    if previous_chromosome is not None and not sparse and previous_end != scope[previous_chromosome]:
                        raise ValueError(f"Incomplete native track: {method}, {previous_chromosome}")
                    if output:
                        output.write("\t".join(map(str, pending)) + "\n")
                        output.close()
                    output = paths[chromosome].open("w")
                    visited.add(chromosome)
                    previous_chromosome, previous_end, pending = chromosome, 0, None
                if start < previous_end or (not sparse and start != previous_end):
                    raise ValueError(f"Overlapping, unsorted or incomplete runs: {method}, line {line_number}")
                if pending and pending[1] == start and pending[2] == depth:
                    pending = (pending[0], end, depth)
                else:
                    if pending:
                        output.write("\t".join(map(str, pending)) + "\n")
                    pending = (start, end, depth)
                previous_end = end
        if not sparse and (visited != set(scope) or previous_end != scope[previous_chromosome]):
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
        if self.left_variance <= 0 or self.right_variance <= 0:
            return "NA"
        return max(-1.0, min(1.0, self.covariance / math.sqrt(self.left_variance * self.right_variance)))


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
    bases = identical = union_bases = union_identical = 0
    below = above = 0
    absolute_error = squared_error = signed_error = left_total = right_total = 0.0
    correlation = Correlation()
    residual_counts = Counter()
    for chromosome, length in scope.items():
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
                union_identical += weight if difference == 0 else 0
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
    flush_residual_counts(residual_database, residual_counts)
    return {
        "bases": bases, "fraction_identical": identical / bases,
        "mae": absolute_error / bases, "rmse": math.sqrt(squared_error / bases),
        "mean_signed_difference": signed_error / bases,
        "fraction_left_below": below / bases, "fraction_left_above": above / bases,
        "nonzero_union_bases": union_bases,
        "union_fraction_identical": union_identical / union_bases if union_bases else "NA",
        "union_mae": absolute_error / union_bases if union_bases else "NA",
        "union_rmse": math.sqrt(squared_error / union_bases) if union_bases else "NA",
        "left_integrated_depth": left_total, "right_integrated_depth": right_total,
        "pearson_r_1bp": correlation.value(),
    }


def regional_counts(directory, expected_regions):
    """Require the exact same unique BED rows and a single length bin."""
    if directory.name.startswith("cfdnalab"):
        from cfdnalab import read_lengths, WindowedLengthCounts

        result = read_lengths(directory / "length_counts.tsv.zst")
        if not isinstance(result, WindowedLengthCounts) or len(result.length_bins()) != 1:
            raise ValueError("Expected a BED count table with a single length bin")
        table = result.data_frame()
        entries = table[["chrom", "start", "end", "count"]].itertuples(index=False, name=None)
    else:
        with open_track(directory / "finaletoolkit.coverage.bed.gz") as source:
            entries = [fields[:3] + fields[-1:] for fields in csv.reader(source, delimiter="\t")
                       if fields and not fields[0].startswith("#")]
    counts = {}
    for chromosome, start, end, count in entries:
        key = str(chromosome), int(start), int(end)
        count = float(count)
        if key in counts or not math.isfinite(count) or count < 0 or not count.is_integer():
            raise ValueError(f"Duplicate region or invalid raw count in {directory}")
        counts[key] = count
    if set(counts) != expected_regions:
        raise ValueError(f"Regional row identity differs from the shared BED: {directory}")
    return counts


def compare_regions(left, right):
    """Keep summed regional counts distinct from unique-fragment totals."""
    differences = [left[key] - right[key] for key in left]
    correlation = Correlation()
    for key in left:
        correlation.add(left[key], right[key])
    left_total, right_total = sum(left.values()), sum(right.values())
    return {"regions": len(left), "left_summed_region_counts": left_total,
            "right_summed_region_counts": right_total,
            "regional_total_ratio": left_total / right_total if right_total else "NA",
            "net_count_difference": sum(differences),
            "mae": sum(abs(value) for value in differences) / len(left),
            "rmse": math.sqrt(sum(value**2 for value in differences) / len(left)),
            "fraction_identical": sum(value == 0 for value in differences) / len(left),
            "pearson_r_regions": correlation.value()}


def main():
    """Compare like quantities, recording the evaluated chromosome scope per row."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profiles", type=Path)
    parser.add_argument("--common-scope", action="store_true", help="Explicitly allow intersection of declared scopes")
    arguments = parser.parse_args()
    directories = sorted(path for path in arguments.profiles.iterdir() if path.is_dir())
    rows = []
    with tempfile.TemporaryDirectory(prefix="coverage-concordance-") as temporary, ExitStack() as resources:
        tracks, regions = {}, {}
        output = resources.enter_context((arguments.profiles.parent / "residuals.tsv").open("w"))
        residual_writer = csv.DictWriter(output, delimiter="\t", lineterminator="\n", fieldnames=[
            "left_profile", "right_profile", "quantity", "same_native_scope", "chromosomes",
            "residual", "bases", "fraction_of_bases",
        ])
        residual_writer.writeheader()
        residual_database = resources.enter_context(closing(sqlite3.connect(Path(temporary) / "residuals.sqlite")))
        # Bound the database cache as well as the Python residual-count batch
        residual_database.execute("PRAGMA cache_size = -2048")
        residual_database.execute("CREATE TABLE residual_counts (residual REAL PRIMARY KEY, bases INTEGER NOT NULL)")
        regional_directories = [path for path in directories if "region-counts" in path.name or path.name == "finaletoolkit-coverage"]
        if regional_directories:
            with (arguments.profiles.parent / "regions.bed").open() as source:
                expected = [(fields[0], int(fields[1]), int(fields[2])) for fields in csv.reader(source, delimiter="\t")]
            if not expected or len(set(expected)) != len(expected):
                raise ValueError("Shared BED must contain unique regions")
            regions = {path.name: regional_counts(path, set(expected)) for path in regional_directories}
        for directory in directories:
            if directory in regional_directories:
                continue
            scope = read_scope(directory)
            destination = Path(temporary) / directory.name
            destination.mkdir()
            tracks[directory.name] = scope, spool_track(directory, scope, destination)
        for left_name, right_name in combinations(tracks, 2):
            # Read depth and fragment coverage answer different questions
            if ("fragment" in left_name) != ("fragment" in right_name):
                continue
            left_scope, left_paths = tracks[left_name]
            right_scope, right_paths = tracks[right_name]
            same_scope = left_scope == right_scope
            if not same_scope and not arguments.common_scope:
                parser.error(f"Scopes differ: {left_name}, {right_name}. Use --common-scope to compare their intersection")
            scope = {chromosome: length for chromosome, length in left_scope.items() if chromosome in right_scope}
            if not scope or any(right_scope[chromosome] != length for chromosome, length in scope.items()):
                raise ValueError("Empty shared scope or mismatched chromosome lengths")
            comparison = {"left_profile": left_name, "right_profile": right_name,
                          "quantity": "fragment-coverage" if "fragment" in left_name else "read-depth",
                          "same_native_scope": same_scope, "chromosomes": ",".join(scope)}
            metrics = compare_tracks(left_paths, right_paths, scope, residual_database)
            rows.append({**comparison, **metrics})
            for residual, bases in residual_database.execute("SELECT residual, bases FROM residual_counts ORDER BY residual"):
                residual_writer.writerow({**comparison, "residual": residual, "bases": bases,
                                          "fraction_of_bases": bases / metrics["bases"]})
            residual_database.execute("DELETE FROM residual_counts")
            residual_database.commit()
        for left_name, right_name in combinations(regions, 2):
            rows.append({"left_profile": left_name, "right_profile": right_name, "quantity": "regional-fragment-counts",
                         "same_native_scope": True, "chromosomes": ",".join(sorted({key[0] for key in regions[left_name]})),
                         **compare_regions(regions[left_name], regions[right_name])})
    if not rows:
        parser.error("No pair of saved profiles measures the same quantity")
    columns = list(dict.fromkeys(column for row in rows for column in row))
    writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter="\t", restval="NA", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
