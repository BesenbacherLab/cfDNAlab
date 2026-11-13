"""Regression tests for benchmark outputs, installation and input preparation.

The job scripts, commands and measurement are tested in test_runners.py.

Expected numerical values are manually derived from the synthetic inputs.

Run on Linux with Bash 4+, zstd and the comparison Python packages. The Snyder
pair reconstruction cases additionally need pysam. No benchmark tools or Slurm
allocation are needed for these small fixtures.
"""

import csv
import gzip
import importlib.util
from itertools import product
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


BENCHMARKS = Path(__file__).resolve().parents[1]


def load_script(name, relative_path):
    """Import a standalone script without running its CLI."""
    specification = importlib.util.spec_from_file_location(
        name, BENCHMARKS / relative_path
    )
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


timings = load_script("timing_summary", "summarize_timings.py")
coverage = load_script("coverage_comparison", "coverage/compare_coverage.py")


class CorrelationTests(unittest.TestCase):
    """Check correlation bounds directly, including nonfinite arithmetic results."""

    def test_accepts_valid_values_and_clips_overshoots_up_to_tolerance(self):
        # Unit variances give a denominator of 1, so covariance is the raw r
        # Overshoots through 1e-5 are accepted on either side of the valid range
        cases = (
            (-1.00001, -1.0),
            (-1.000005, -1.0),
            (-1.0, -1.0),
            (-0.25, -0.25),
            (0.0, 0.0),
            (0.25, 0.25),
            (1.0, 1.0),
            (1.000005, 1.0),
            (1.00001, 1.0),
        )
        for covariance, expected in cases:
            with self.subTest(covariance=covariance):
                correlation = coverage.Correlation()
                correlation.left_variance = correlation.right_variance = 1.0
                correlation.covariance = covariance

                self.assertEqual(correlation.value(), expected)

    def test_rejects_larger_overshoots_and_nonfinite_values(self):
        # The next float beyond each tolerated boundary must fail
        # NaN must fail too, rather than reaching the clamp and becoming +1
        outside_limit = math.nextafter(1.00001, math.inf)
        for covariance in (
            -outside_limit, outside_limit, -1.1, 1.1, math.nan, -math.inf, math.inf
        ):
            with self.subTest(covariance=covariance):
                correlation = coverage.Correlation()
                correlation.left_variance = correlation.right_variance = 1.0
                correlation.covariance = covariance

                with self.assertRaisesRegex(RuntimeError, "Invalid correlation value"):
                    correlation.value()


class CoverageTests(unittest.TestCase):
    """Exercise native formats and base weights at the CLI boundary."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.profiles = Path(self.temporary.name) / "profiles"
        self.profiles.mkdir()

    def profile(self, name, contents, scope="chr1\t100\n", second_frame=None):
        """Create a tiny native compressed track and declared scope."""
        directory = self.profiles / name
        directory.mkdir()
        (directory / "chrom.sizes").write_text(scope)
        if name.startswith("cfdnalab"):
            data = subprocess.check_output(
                ["zstd", "-q", "-c"], input=contents.encode()
            )
            if second_frame is not None:
                data += subprocess.check_output(
                    ["zstd", "-q", "-c"], input=second_frame.encode()
                )
            (directory / "fcoverage.per_position.bedgraph.zst").write_bytes(data)
        else:
            filename = (
                "positions.tsv.gz"
                if name in ("samtools-depth", "pandepth")
                else "coverage.per-base.bed.gz"
            )
            with gzip.open(directory / filename, "wt") as destination:
                destination.write(contents)

    def compare(self, *options):
        """Call the supported CLI, preserving error output for assertions."""
        return subprocess.run(
            [
                sys.executable,
                str(BENCHMARKS / "coverage/compare_coverage.py"),
                str(self.profiles),
                *options,
            ],
            text=True,
            capture_output=True,
        )

    def test_weights_errors_by_bases_and_reports_nonzero_union(self):
        # A difference of 1 on 10/100 bases gives MAE .1 and RMSE sqrt(.1)
        self.profile("cfdnalab", "chr1\t0\t10\t1\n")
        self.profile("mosdepth", "chr1\t0\t100\t0\n")
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(int(row["bases"]), 100)
        self.assertAlmostEqual(float(row["fraction_identical"]), 0.9)
        self.assertAlmostEqual(float(row["mae"]), 0.1)
        self.assertAlmostEqual(float(row["rmse"]), math.sqrt(0.1))
        self.assertEqual(int(row["nonzero_union_bases"]), 10)
        self.assertAlmostEqual(float(row["union_mae"]), 1)
        self.assertEqual(row["pearson_r_1bp"], "NA")

    def test_sparse_samtools_reports_jaccard_and_truncation_in_bam_order(self):
        # Each chromosome has 6 covered bases in cfDNAlab. samtools shares
        # 6, 3 and 0 bases respectively, so genome-wide Jaccard is 9/18 = 0.5
        # The intentionally non-lexicographic header order is the BAM order
        scope = "chr2\t10\nchr10\t10\nchr1\t10\n"
        self.profile("cfdnalab", "chr2\t2\t8\t1\nchr10\t2\t8\t1\nchr1\t2\t8\t1\n", scope)
        positions = "".join(f"chr2\t{position}\t1\n" for position in range(3, 9))
        positions += "".join(f"chr10\t{position}\t1\n" for position in range(3, 6))
        self.profile("samtools-depth", positions, scope)
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row, = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
        self.assertEqual(int(row["nonzero_intersection_bases"]), 9)
        self.assertEqual(int(row["nonzero_union_bases"]), 18)
        self.assertEqual(float(row["coverage_jaccard"]), 0.5)
        with (self.profiles.parent / "coverage_jaccard.tsv").open() as source:
            diagnostics = list(csv.DictReader(source, delimiter="\t"))
        self.assertEqual([row["chromosome"] for row in diagnostics], ["chr2", "chr10", "chr1"])
        self.assertEqual([float(row["coverage_jaccard"]) for row in diagnostics], [1.0, 0.5, 0.0])
        self.assertEqual([row["possible_samtools_truncation"] for row in diagnostics], ["False", "True", "True"])
        self.assertIn("Possible samtools-depth truncation", result.stderr)

    def test_internal_zero_chromosome_does_not_flag_truncation(self):
        # Coverage resumes after the zero-overlap chromosome, so there is no
        # zero-coverage suffix even though native chromosome Jaccards differ
        chromosomes = [
            {"chromosome": "chr1", "nonzero_union_bases": 6, "coverage_jaccard": 0.5,
             "left_nonzero_bases": 6, "right_nonzero_bases": 3},
            {"chromosome": "chr2", "nonzero_union_bases": 6, "coverage_jaccard": 0.0,
             "left_nonzero_bases": 6, "right_nonzero_bases": 0},
            {"chromosome": "chr3", "nonzero_union_bases": 6, "coverage_jaccard": 1.0,
             "left_nonzero_bases": 6, "right_nonzero_bases": 6},
        ]
        self.assertEqual(coverage.possible_samtools_truncation("cfdnalab", "samtools-depth", chromosomes), set())

    def test_empty_sparse_samtools_output_represents_zero_coverage(self):
        self.profile("cfdnalab", "")
        self.profile("samtools-depth", "")
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row, = csv.DictReader(result.stdout.splitlines(), delimiter="\t")
        self.assertEqual(float(row["mae"]), 0)
        self.assertEqual(row["coverage_jaccard"], "NA")
        self.assertNotIn("Possible samtools-depth truncation", result.stderr)

    def test_rejects_overflow_in_a_derived_coverage_metric(self):
        # Each depth and squared difference is finite, but the weighted sum overflows
        self.profile("cfdnalab", "chr1\t0\t100\t1e154\n")
        self.profile("mosdepth", "chr1\t0\t100\t0\n")
        histogram = self.profiles.parent / "residuals.tsv"
        histogram.write_text("Previous successful histogram\n")

        result = self.compare()

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Invalid rmse", result.stderr)
        self.assertEqual(histogram.read_text(), "Previous successful histogram\n")

    def test_residual_distribution_preserves_sign_and_counts_genomic_bases(self):
        # Residuals -2, 0, +1 cover 20, 70, 30 bases across both contigs
        # MAE=(40+30)/120=7/12, RMSE=sqrt(110/120), signed mean=-10/120
        scope = "chr1\t100\nchr2\t20\n"
        self.profile("cfdnalab", "chr1\t0\t10\t3\nchr2\t0\t20\t1\n", scope)
        self.profile(
            "mosdepth", "chr2\t0\t20\t0\nchr1\t0\t30\t2\nchr1\t30\t100\t0\n", scope
        )
        output = self.profiles.parent / "residuals.tsv"

        result = self.compare()

        self.assertEqual(result.returncode, 0, result.stderr)
        summary = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        with output.open() as source:
            residuals = list(csv.DictReader(source, delimiter="\t"))
        self.assertEqual([float(row["residual"]) for row in residuals], [-2, 0, 1])
        self.assertEqual([int(row["bases"]) for row in residuals], [20, 70, 30])
        for row, expected_fraction in zip(residuals, (1 / 6, 7 / 12, 1 / 4)):
            self.assertAlmostEqual(float(row["fraction_of_bases"]), expected_fraction)
            self.assertEqual(row["left_profile"], "cfdnalab")
            self.assertEqual(row["right_profile"], "mosdepth")
            self.assertEqual(row["quantity"], "read-depth")
            self.assertEqual(row["chromosomes"], "chr1,chr2")
            self.assertEqual(row["same_native_scope"], "True")
        self.assertEqual(int(summary["bases"]), 120)
        self.assertAlmostEqual(float(summary["mae"]), 7 / 12)
        self.assertAlmostEqual(float(summary["rmse"]), math.sqrt(11 / 12))
        self.assertAlmostEqual(float(summary["mean_signed_difference"]), -1 / 12)

    def test_residual_distribution_includes_all_zero_profiles(self):
        # An omitted sparse zero track and an explicit zero track cover 15,000 bases
        self.profile("cfdnalab", "", "chr1\t15000\n")
        self.profile("mosdepth", "chr1\t0\t15000\t0\n", "chr1\t15000\n")
        output = self.profiles.parent / "residuals.tsv"

        result = self.compare()

        self.assertEqual(result.returncode, 0, result.stderr)
        with output.open() as source:
            residuals = list(csv.DictReader(source, delimiter="\t"))
        self.assertEqual(len(residuals), 1)
        self.assertEqual(float(residuals[0]["residual"]), 0)
        self.assertEqual(int(residuals[0]["bases"]), 15000)
        self.assertEqual(float(residuals[0]["fraction_of_bases"]), 1)
        summary = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(float(summary["mae"]), 0)
        self.assertEqual(summary["union_mae"], "NA")
        self.assertEqual(summary["pearson_r_1bp"], "NA")
        # Rerunning regenerates the fixed output rather than appending or failing
        repeated = self.compare()
        self.assertEqual(repeated.returncode, 0, repeated.stderr)
        with output.open() as source:
            regenerated = list(csv.DictReader(source, delimiter="\t"))
        self.assertEqual(regenerated, residuals)

    def test_residual_distribution_preserves_fractional_differences(self):
        # Residual -.25 covers 75 bases and +.5 covers 25: MAE=.3125
        self.profile("cfdnalab", "chr1\t0\t25\t1.25\nchr1\t25\t100\t0.5\n")
        self.profile("mosdepth", "chr1\t0\t100\t0.75\n")
        output = self.profiles.parent / "residuals.tsv"

        result = self.compare()

        self.assertEqual(result.returncode, 0, result.stderr)
        with output.open() as source:
            residuals = list(csv.DictReader(source, delimiter="\t"))
        self.assertEqual([float(row["residual"]) for row in residuals], [-0.25, 0.5])
        self.assertEqual([int(row["bases"]) for row in residuals], [75, 25])
        summary = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertAlmostEqual(float(summary["mae"]), 0.3125)
        self.assertEqual(summary["pearson_r_1bp"], "NA")

    def test_residual_counts_survive_large_value_sets_and_multiple_pairs(self):
        # Depths 0..4096 each cover a base, then depth 1 covers ten more
        # The repeated residual must accumulate after a batch spills to disk
        scope = "chr1\t4107\n"
        self.profile("cfdnalab", "", scope)
        depths = "".join(
            f"chr1\t{position}\t{position + 1}\t{position}\n"
            for position in range(4097)
        )
        self.profile("mosdepth", depths + "chr1\t4097\t4107\t1\n", scope)
        self.profile(
            "pandepth",
            "".join(f"chr1\t{position}\t0\n" for position in range(4107)),
            scope,
        )
        output = self.profiles.parent / "residuals.tsv"

        result = self.compare()

        self.assertEqual(result.returncode, 0, result.stderr)
        comparisons = {}
        with output.open() as source:
            for row in csv.DictReader(source, delimiter="\t"):
                comparisons.setdefault(
                    (row["left_profile"], row["right_profile"]), []
                ).append(row)
        self.assertEqual(len(comparisons), 3)
        negative = comparisons[("cfdnalab", "mosdepth")]
        zero = comparisons[("cfdnalab", "pandepth")]
        positive = comparisons[("mosdepth", "pandepth")]
        self.assertEqual(
            [float(row["residual"]) for row in negative], list(range(-4096, 1))
        )
        self.assertEqual([int(row["bases"]) for row in negative], [1] * 4095 + [11, 1])
        self.assertEqual(
            [float(row["residual"]) for row in positive], list(range(4097))
        )
        self.assertEqual([int(row["bases"]) for row in positive], [1, 11] + [1] * 4095)
        self.assertEqual(len(zero), 1)
        self.assertEqual(float(zero[0]["residual"]), 0)
        self.assertEqual(int(zero[0]["bases"]), 4107)
        for residuals in comparisons.values():
            self.assertEqual(sum(int(row["bases"]) for row in residuals), 4107)
            self.assertAlmostEqual(
                sum(float(row["fraction_of_bases"]) for row in residuals), 1
            )

    def test_tile_boundaries_contig_order_and_concatenated_frames_are_equivalent(self):
        scope = "chr1\t100\nchr2\t80\n"
        self.profile(
            "cfdnalab",
            "chr1\t0\t30\t2\n",
            scope,
            second_frame="chr1\t30\t100\t2\nchr2\t0\t80\t3\n",
        )
        self.profile("mosdepth", "chr2\t0\t80\t3\nchr1\t0\t100\t2\n", scope)
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(float(row["fraction_identical"]), 1)
        self.assertEqual(float(row["mae"]), 0)
        self.assertEqual(float(row["left_integrated_depth"]), 440)

    def test_position_origins_agree(self):
        self.profile(
            "samtools-depth",
            "".join(f"chr1\t{position}\t2\n" for position in range(1, 101)),
        )
        self.profile(
            "pandepth", "".join(f"chr1\t{position}\t2\n" for position in range(100))
        )
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(float(row["fraction_identical"]), 1)
        self.assertEqual(float(row["left_integrated_depth"]), 200)

    def test_rejects_truncated_dense_track(self):
        self.profile("cfdnalab", "chr1\t0\t100\t1\n")
        self.profile("mosdepth", "chr1\t0\t99\t1\n")
        result = self.compare()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Incomplete native track", result.stderr)

    def test_rejects_overlapping_or_out_of_bounds_runs(self):
        self.profile("mosdepth", "chr1\t0\t100\t1\n")
        for contents in ("chr1\t0\t50\t1\nchr1\t49\t100\t1\n", "chr1\t0\t101\t1\n"):
            directory = self.profiles / "cfdnalab"
            if directory.exists():
                for path in directory.iterdir():
                    path.unlink()
                directory.rmdir()
            self.profile("cfdnalab", contents)
            self.assertNotEqual(self.compare().returncode, 0)

    def test_pearson_uses_individual_bases_instead_of_binned_means(self):
        # In 1000-base units, weights are (10,3,2) for pairs (1,2), (3,6), (3,0)
        # Sums: left=25, right=38, left squared=55, right squared=148, product=74
        # Pearson=(15*74-25*38)/sqrt((15*55-625)*(15*148-1444))=4/sqrt(97)
        # The 10 kb means would instead give r=1, hiding within-bin disagreement
        # Errors are 1 on 10,000 bases and 3 on 5,000 bases
        self.profile(
            "cfdnalab", "chr1\t0\t10000\t1\nchr1\t10000\t15000\t3\n", "chr1\t15000\n"
        )
        self.profile(
            "mosdepth",
            "chr1\t0\t10000\t2\nchr1\t10000\t13000\t6\nchr1\t13000\t15000\t0\n",
            "chr1\t15000\n",
        )
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertAlmostEqual(float(row["pearson_r_1bp"]), 4 / math.sqrt(97))
        self.assertAlmostEqual(float(row["mae"]), 5 / 3)
        self.assertAlmostEqual(float(row["rmse"]), math.sqrt(11 / 3))

    def test_scope_intersection_requires_explicit_selection(self):
        self.profile("cfdnalab", "chr1\t0\t100\t1\n", "chr1\t100\nchr2\t80\n")
        self.profile("mosdepth", "chr1\t0\t100\t1\n")
        self.assertNotEqual(self.compare().returncode, 0)
        output = self.profiles.parent / "residuals.tsv"
        result = self.compare("--common-scope")
        self.assertEqual(result.returncode, 0, result.stderr)
        row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(row["same_native_scope"], "False")
        self.assertEqual(row["chromosomes"], "chr1")
        with output.open() as source:
            residuals = list(csv.DictReader(source, delimiter="\t"))
        self.assertEqual(len(residuals), 1)
        self.assertEqual(float(residuals[0]["residual"]), 0)
        self.assertEqual(int(residuals[0]["bases"]), 100)
        self.assertEqual(residuals[0]["chromosomes"], "chr1")
        self.assertEqual(residuals[0]["same_native_scope"], "False")

    def test_regional_comparison_preserves_fractional_counts_and_coordinates(self):
        # Fragments [10,40), [40,80), [60,90) give normalized totals (1.25,1.75)
        # The 40 bp boundary fragment contributes 10/40 and 30/40 to the windows
        # Any-overlap counts are (2,2). Differences (-.75,-.25) give MAE .5,
        # RMSE sqrt(5)/4 and summed-count ratio 3/4
        (self.profiles.parent / "regions.bed").write_text(
            "chr1\t0\t50\t.\nchr1\t50\t100\t.\n"
        )
        left = self.profiles / "cfdnalab-region-counts"
        left.mkdir()
        table = (
            "chromosome\tstart\tend\ttotal_fragment_mass\tblacklisted_positions\n"
            "chr1\t0\t50\t1.25\t0\nchr1\t50\t100\t1.75\t0\n"
        )
        (left / "length_normalized.fcoverage.total.tsv.zst").write_bytes(
            subprocess.check_output(["zstd", "-q", "-c"], input=table.encode())
        )
        right = self.profiles / "finaletoolkit-coverage"
        right.mkdir()
        with gzip.open(right / "finaletoolkit.coverage.bed.gz", "wt") as destination:
            destination.write("chr1\t50\t100\t.\t2.0\nchr1\t0\t50\t.\t2.0\n")
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(float(row["mae"]), 0.5)
        self.assertAlmostEqual(float(row["rmse"]), math.sqrt(5) / 4)
        self.assertEqual(float(row["regional_total_ratio"]), 0.75)
        self.assertEqual(float(row["left_summed_region_counts"]), 3)
        self.assertEqual(float(row["right_summed_region_counts"]), 4)
        self.assertEqual(float(row["net_count_difference"]), -1)
        self.assertEqual(row["pearson_r_regions"], "NA")
        # Regional-only comparisons keep the automatic residual file header-only
        with (self.profiles.parent / "residuals.tsv").open() as source:
            self.assertEqual(list(csv.DictReader(source, delimiter="\t")), [])
        # A shifted endpoint must fail rather than align by row number
        with gzip.open(right / "finaletoolkit.coverage.bed.gz", "wt") as destination:
            destination.write("chr1\t0\t51\t.\t1\nchr1\t50\t100\t.\t3\n")
        self.assertNotEqual(self.compare().returncode, 0)

    def test_regional_reader_rejects_invalid_counts_and_fractional_finale_counts(self):
        expected_regions = {("chr1", 0, 50)}
        for method, invalid_values in (
            ("cfdnalab-region-counts", ("-0.25", "nan", "inf")),
            ("finaletoolkit-coverage", ("-1", "nan", "inf", "0.5")),
        ):
            directory = self.profiles / method
            directory.mkdir()
            for invalid_value in invalid_values:
                with self.subTest(method=method, value=invalid_value):
                    if method.startswith("cfdnalab"):
                        table = (
                            "chromosome\tstart\tend\ttotal_fragment_mass\tblacklisted_positions\n"
                            f"chr1\t0\t50\t{invalid_value}\t0\n"
                        )
                        (directory / "length_normalized.fcoverage.total.tsv.zst").write_bytes(
                            subprocess.check_output(["zstd", "-q", "-c"], input=table.encode())
                        )
                    else:
                        with gzip.open(directory / "finaletoolkit.coverage.bed.gz", "wt") as destination:
                            destination.write(f"chr1\t0\t50\t.\t{invalid_value}\n")
                    with self.assertRaises(ValueError):
                        coverage.regional_counts(directory, expected_regions)

    def test_native_regional_formats_report_pandepth_missing_one_base_tail(self):
        # A 50 bp window at mean depth 2 has 100 covered bases. The final
        # 1 bp window at depth 3 has 3. PanDepth omits that last 1 bp window
        # Its inclusive [1,50] becomes BED [0,50), with covered-base total 100
        (self.profiles.parent / "regions.bed").write_text(
            "chr1\t0\t50\t.\nchr1\t50\t51\t.\n"
        )
        for method in (
            "cfdnalab-regions", "mosdepth-regions", "samtools-regions", "pandepth-regions"
        ):
            directory = self.profiles / method
            directory.mkdir()
            if method == "cfdnalab-regions":
                table = (
                    "chromosome\tstart\tend\ttotal_coverage\tblacklisted_positions\n"
                    "chr1\t0\t50\t100\t0\nchr1\t50\t51\t3\t0\n"
                )
                (directory / "fcoverage.total.tsv.zst").write_bytes(
                    subprocess.check_output(["zstd", "-q", "-c"], input=table.encode())
                )
            elif method == "mosdepth-regions":
                with gzip.open(directory / "coverage.regions.bed.gz", "wt") as destination:
                    destination.write("chr1\t0\t50\t2.00\nchr1\t50\t51\t3.00\n")
            elif method == "samtools-regions":
                (directory / "samtools.bedcov.tsv").write_text(
                    "chr1\t50\t51\t.\t3\nchr1\t0\t50\t.\t100\n"
                )
            else:
                with gzip.open(directory / "pandepth.win.stat.gz", "wt") as destination:
                    destination.write(
                        "#Chr\tStart\tEnd\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n"
                        "chr1\t1\t50\t50\t50\t100\t100.00\t2.00\n"
                        "##RegionLength: 50\tCoveredSite: 50\n"
                    )
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(len(rows), 6)
        for row in rows:
            self.assertEqual(row["quantity"], "regional-read-depth")
            includes_pandepth = "pandepth-regions" in (row["left_profile"], row["right_profile"])
            self.assertEqual(int(row["regions"]), 1 if includes_pandepth else 2)
            self.assertEqual(float(row["mae"]), 0)
            self.assertEqual(float(row["left_integrated_depth"]), 100 if includes_pandepth else 103)
            self.assertEqual(float(row["right_integrated_depth"]), 100 if includes_pandepth else 103)
            self.assertEqual(int(row["left_missing_regions"]) + int(row["right_missing_regions"]), int(includes_pandepth))
            self.assertEqual(float(row["net_coverage_difference"]), 0)
        # A fragment-count table must not be paired with covered-base totals
        counts = self.profiles / "finaletoolkit-coverage"
        counts.mkdir()
        with gzip.open(counts / "finaletoolkit.coverage.bed.gz", "wt") as destination:
            destination.write("chr1\t0\t50\t.\t2\nchr1\t50\t51\t.\t1\n")
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))), 6)


class TimingTests(unittest.TestCase):
    """Use literal timing/accounting fixtures rather than invoking Slurm."""

    def test_join_preserves_correct_step_memory_and_missing_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "timings.tsv"
            path.write_text("run_order\tjob_id\tstep_id\n1\t123\t0\n2\t123\t1\n")
            rows = timings.read_runs(
                [path],
                {
                    "123.0": {
                        "NTasks": "1",
                        "MaxRSS": "2048K",
                        "State": "COMPLETED",
                        "ExitCode": "0:0",
                    }
                },
            )
            self.assertEqual(rows[0]["slurm_max_rss_kib"], 2048)
            self.assertEqual(rows[1]["slurm_max_rss_kib"], "NA")

    def test_sample_statistics_use_sample_sd_and_inclusive_quartiles(self):
        self.assertEqual(timings.distribution([1, 3, 5]), [3, 3, 2, 3, 2, 4])
        self.assertEqual(timings.distribution([4]), [1, 4, "NA", 4, 4, 4])

    def test_two_method_budget_includes_a_warmup_and_three_measurements(self):
        columns = "run_order repetition warmup pilot method tool_version cores node wall_seconds user_seconds system_seconds flush_seconds process_plus_flush_seconds output_bytes max_rss_kib tool_exit_code step_exit_code job_id step_id".split()
        values = [
            "1 1 false 1 cfdnalab checkout 1 node 720 700 20 0 720 100 200 0 0 123 0".split(),
            "2 1 false 1 slow release 1 node 7200 7000 200 0 7200 100 200 0 0 123 1".split(),
        ]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "pilot.tsv"
            with path.open("w") as destination:
                writer = csv.writer(destination, delimiter="\t")
                writer.writerow(columns)
                writer.writerows(values)
            result = subprocess.run(
                [
                    sys.executable,
                    str(BENCHMARKS / "summarize_timings.py"),
                    str(path),
                    "--estimate",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            # 4*(720+7200) = 31680 seconds = 8.8 hours, or 11 with 25% margin
            self.assertIn("8.80 h raw", result.stderr)
            self.assertIn("11.00 h with margin", result.stderr)

    def test_griffin_budget_never_adds_a_full_input_warmup(self):
        # At 30 h per run, three Griffin measurements take 90 h raw or 112.5 h
        # with 25% margin. A deliberately large normal warmup count adds nothing
        runs = [
            {
                "method": method,
                "cores": "4",
                "pilot": "1",
                "tool_exit_code": "0",
                "step_exit_code": "0",
                "slurm_state": "NA",
                "slurm_exit_code": "NA",
            }
            for method in ("griffin-gc",)
        ]
        summary = [
            {
                "method": method,
                "cores": "4",
                "process_plus_flush_seconds_mean": 108000,
            }
            for method in ("griffin-gc",)
        ]
        report = timings.pilot_budgets(runs, ["method", "cores"], summary, 10, 3, 0.25)
        self.assertEqual(len(report), 1)
        for budget in report:
            self.assertEqual(budget["planned_runs"], 3)
            self.assertEqual(budget["budget_minutes"], 6750)
            self.assertEqual(budget["budget_status"], "complete")

    def test_accounting_oom_excludes_an_apparently_successful_native_run(self):
        self.assertFalse(
            timings.successful(
                {
                    "tool_exit_code": "0",
                    "step_exit_code": "0",
                    "slurm_state": "OUT_OF_MEMORY",
                    "slurm_exit_code": "0:125",
                }
            )
        )


class RunControlTests(unittest.TestCase):
    """Check Picard's unusual version status and the installer without tool installs."""

    def test_picard_accepts_only_expected_status_and_version(self):
        with tempfile.TemporaryDirectory() as temporary:
            executable = Path(temporary) / "picard"
            for version, status, succeeds in (
                ("Version:3.5.0", 1, True),
                ("Version:3.5.0", 2, False),
                ("Version:3.4.0", 1, False),
                ("error opening picard-3.5.0.jar", 1, False),
                ("", 1, False),
            ):
                executable.write_text(
                    f"#!/usr/bin/env bash\nprintf '%s\\n' '{version}'\nexit {status}\n"
                )
                executable.chmod(0o755)
                result = subprocess.run(
                    [
                        "bash",
                        "-euc",
                        'source "$BENCH_COMMON"; benchmark_picard_version "$PICARD_FIXTURE"',
                    ],
                    env={
                        **os.environ,
                        "BENCH_COMMON": str(BENCHMARKS / "common.sh"),
                        "PICARD_FIXTURE": str(executable),
                    },
                    capture_output=True,
                )
                self.assertEqual(result.returncode == 0, succeeds)

    def test_installer_resumes_and_never_completes_a_failed_reinstall(self):
        # A fake Conda boundary permits info but fails installation commands
        # Each reinstall reason must reach that failure and remain outdated
        for reason in ("complete", "outdated", "version", "missing", "force"):
            with (
                self.subTest(reason=reason),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                installation, fake_bin, conda_base = (
                    root / "installation with spaces",
                    root / "bin",
                    root / "conda",
                )
                fake_bin.mkdir()
                (conda_base / "etc/profile.d").mkdir(parents=True)
                (conda_base / "etc/profile.d/conda.sh").write_text("# Conda fixture\n")
                executable = fake_bin / "conda"
                executable.write_text(
                    r'''#!/usr/bin/env bash
if [[ "$1" == info ]]; then
    for setting in CONDA_PKGS_DIRS CONDA_ENVS_PATH CONDA_REGISTER_ENVS \
        PIP_CACHE_DIR PIP_USER PIP_PREFIX PIP_TARGET CARGO_HOME \
        APPTAINER_CACHEDIR APPTAINER_TMPDIR APPTAINER_CONFIGDIR \
        XDG_CACHE_HOME XDG_CONFIG_HOME XDG_DATA_HOME MPLCONFIGDIR \
        TMPDIR TMP TEMP PYTHONDONTWRITEBYTECODE HOME; do
        printf '%s\t%s\n' "$setting" "${!setting}"
    done > "$BENCH_STORAGE_RECORD"
    printf '%s\n' "$BENCH_FAKE_CONDA_BASE"
else
    exit 9
fi
'''
                )
                executable.chmod(0o755)
                (installation / "envs/comparison/bin").mkdir(parents=True)
                if reason != "missing":
                    (installation / "envs/comparison/bin/python").touch()
                version = "old" if reason == "version" else "checkout"
                status = "outdated" if reason == "outdated" else "complete"
                manifest = installation / "manifest.tsv"
                manifest.write_text(
                    f"tool\tversion\tstatus\tartifacts\ncomparison\t{version}\t{status}\tenvs/comparison/bin/python\n"
                )
                command = [
                    "bash",
                    str(BENCHMARKS / "install_tools.sh"),
                    str(installation),
                    "--only",
                    "comparison",
                ]
                if reason == "force":
                    command += ["--force", "comparison"]
                storage_record = root / "storage.tsv"
                # Inherited cache locations must not override the supplied installation root
                inherited_paths = {
                    name: str(root / "outside-installation" / name)
                    for name in (
                        "CONDA_PKGS_DIRS", "CONDA_ENVS_PATH", "PIP_CACHE_DIR",
                        "PIP_PREFIX", "PIP_TARGET", "CARGO_HOME", "APPTAINER_CACHEDIR",
                        "APPTAINER_TMPDIR", "APPTAINER_CONFIGDIR", "XDG_CACHE_HOME",
                        "XDG_CONFIG_HOME", "XDG_DATA_HOME", "MPLCONFIGDIR", "TMPDIR",
                        "TMP", "TEMP",
                    )
                }
                result = subprocess.run(
                    command,
                    env={
                        **os.environ,
                        **inherited_paths,
                        "PATH": f"{fake_bin}:{os.environ['PATH']}",
                        "BENCH_FAKE_CONDA_BASE": str(conda_base),
                        "BENCH_STORAGE_RECORD": str(storage_record),
                    },
                    capture_output=True,
                    text=True,
                )
                if reason == "complete":
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn("Keeping comparison checkout", result.stdout)
                else:
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(
                        "comparison\tcheckout\toutdated\t", manifest.read_text()
                    )
                # Capture settings at the first Conda call, before any package work
                storage = dict(csv.reader(
                    storage_record.read_text().splitlines(), delimiter="\t"
                ))
                expected_paths = {
                    "CONDA_PKGS_DIRS": "cache/conda/pkgs",
                    "CONDA_ENVS_PATH": "envs",
                    "PIP_CACHE_DIR": "cache/pip",
                    "CARGO_HOME": "cache/cargo",
                    "APPTAINER_CACHEDIR": "cache/apptainer",
                    "APPTAINER_TMPDIR": "tmp",
                    "APPTAINER_CONFIGDIR": "config/apptainer",
                    "XDG_CACHE_HOME": "cache",
                    "XDG_CONFIG_HOME": "config",
                    "XDG_DATA_HOME": "share",
                    "MPLCONFIGDIR": "cache/matplotlib",
                    "TMPDIR": "tmp", "TMP": "tmp", "TEMP": "tmp",
                }
                for setting, relative_path in expected_paths.items():
                    expected_path = installation.resolve() / relative_path
                    self.assertEqual(storage[setting], str(expected_path), setting)
                    self.assertTrue(expected_path.is_dir(), setting)
                self.assertEqual(storage["CONDA_REGISTER_ENVS"], "false")
                self.assertEqual(storage["PIP_USER"], "false")
                self.assertEqual(storage["PYTHONDONTWRITEBYTECODE"], "1")
                self.assertEqual(storage["PIP_PREFIX"], "")
                self.assertEqual(storage["PIP_TARGET"], "")
                self.assertEqual(storage["HOME"], os.environ["HOME"])
                self.assertFalse((root / "outside-installation").exists())

    def test_freia_environment_is_conda_only_and_uses_the_supplied_root(self):
        # The upstream pip packages are not needed by the benchmark's native scripts
        # Every Conda dependency, including its existing version pin, must survive
        dependencies = (
            "biopython=1.78", "fastparquet=0.5", "logomaker=0.8", "mpmath=1.1.0",
            "numpy=1.18.5", "pandas=1.1.5", "patsy=0.5.1", "pip=21.0.1",
            "pysam=0.16.0.1", "python=3.8", "samtools=1.12",
            "scikit-learn=0.24.1", "scikit-bio=0.5.6", "seaborn", "umap-learn",
        )
        upstream = (
            "name: FrEIA\nprefix: ~/miniconda3/envs/FrEIA\n"
            "channels:\n  - conda-forge\n  - bioconda\ndependencies:\n"
            + "".join(f"  - {dependency}\n" for dependency in dependencies)
        ).replace("  - pip=21.0.1\n", "  - pip=21.0.1\n  - pip:\n    - statannot\n    - combat\n")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "installation with spaces"
            source = root / "FrEIA-source/workflow/envs/FrEIA_env.yaml"
            source.parent.mkdir(parents=True)
            source.write_text(upstream)
            arguments = root / "conda-arguments.txt"

            # Fake the source checkout and stop at Conda. No downloads or installation occur
            result = subprocess.run(
                ["bash", "-euc", r'''
source <(sed -n '/^install_freia() {$/,/^}$/p' "$BENCH_INSTALLER")
install_root=$BENCH_INSTALL_ROOT
git() { return 0; }
conda() {
    printf '%s\n' "$@" > "$BENCH_CONDA_ARGUMENTS"
    return 9
}
install_freia
'''],
                env={
                    **os.environ,
                    "BENCH_INSTALLER": str(BENCHMARKS / "install_tools.sh"),
                    "BENCH_INSTALL_ROOT": str(root),
                    "BENCH_CONDA_ARGUMENTS": str(arguments),
                },
                capture_output=True, text=True,
            )

            self.assertEqual(result.returncode, 9, result.stderr)
            generated = (root / "freia-environment.yaml").read_text()
            self.assertEqual(
                generated.split("dependencies:\n")[0],
                "channels:\n  - conda-forge\n  - bioconda\n  - nodefaults\n",
            )
            self.assertEqual(
                generated.split("dependencies:\n")[1].splitlines(),
                [f"  - {dependency}" for dependency in dependencies],
            )
            self.assertEqual(arguments.read_text().splitlines(), [
                "env", "create", "--yes", "--prefix", str(root / "envs/freia"),
                "--file", str(root / "freia-environment.yaml"),
            ])
            self.assertEqual(source.read_text(), upstream)

    def test_finaletoolkit_installs_loess_with_pip_before_packing(self):
        # Local shell functions replace Conda, Python and FinaleToolkit
        # No package manager or installed tool is invoked by this fixture
        revision = "d16632ec5ae8999e774234457f6b0f9405c8f590"
        source_requirement = f"FinaleToolkit @ git+https://github.com/epifluidlab/FinaleToolkit.git@{revision}"
        pip_command = f"python -m pip install loess==2.1.2 {source_requirement}"
        for pip_status in (0, 9):
            with self.subTest(pip_status=pip_status), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                (root / "records").mkdir()
                packaging = root / "envs/packaging/bin"
                packaging.mkdir(parents=True)
                pack = packaging / "conda-pack"
                pack.write_text('#!/usr/bin/env bash\nprintf "pack\\n" >> "$BENCH_CALLS"\n')
                pack.chmod(0o755)
                calls = root / "calls.txt"
                result = subprocess.run(
                    ["bash", "-euc", r'''
source <(sed -n '/^install_finaletoolkit() {$/,/^}$/p' "$BENCH_INSTALLER")
install_root=$BENCH_INSTALL_ROOT
# The function uses the installer's source pin as well as its installation root
declare -A versions=([finaletoolkit]="$BENCH_FINALE_REVISION")
conda() {
    printf 'conda %s\n' "$*" >> "$BENCH_CALLS"
    if [[ "$*" == *loess=* ]]; then return 7; fi
    if [[ "$1" == activate ]]; then export CONDA_PREFIX=$2; fi
}
python() {
    printf 'python %s\n' "$*" >> "$BENCH_CALLS"
    if [[ "$1" == -m && "$2" == pip && "$3" == install ]]; then
        # Keep the VCS requirement as a single argument, including its spaces
        [[ "$#" == 5 && "$4" == loess==2.1.2 && "$5" == "$BENCH_FINALE_REQUIREMENT" ]] || return 17
        return "$BENCH_PIP_STATUS"
    fi
    if [[ "$*" == '-m pip freeze' ]]; then printf 'loess==2.1.2\n'; fi
}
finaletoolkit() { printf 'version\n' >> "$BENCH_CALLS"; }
install_finaletoolkit
'''],
                    env={
                        **os.environ,
                        "BENCH_INSTALLER": str(BENCHMARKS / "install_tools.sh"),
                        "BENCH_INSTALL_ROOT": str(root),
                        "BENCH_CALLS": str(calls),
                        "BENCH_PIP_STATUS": str(pip_status),
                        "BENCH_FINALE_REVISION": revision,
                        "BENCH_FINALE_REQUIREMENT": source_requirement,
                    },
                    capture_output=True, text=True,
                )
                commands = calls.read_text().splitlines()
                self.assertIn(pip_command, commands)
                revision_file = root / "envs/finaletoolkit/share/finaletoolkit/revision.txt"
                if pip_status:
                    self.assertEqual(result.returncode, pip_status, result.stderr)
                    self.assertNotIn("pack", commands)
                    self.assertNotIn("version", commands)
                    self.assertFalse(revision_file.exists())
                else:
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertLess(commands.index(pip_command), commands.index("python -m pip check"))
                    self.assertLess(commands.index("python -m pip check"), commands.index("pack"))
                    self.assertEqual(revision_file.read_text(), revision + "\n")
                    self.assertEqual(
                        (root / "records/finaletoolkit.pip.txt").read_text(),
                        "loess==2.1.2\n",
                    )


class BiologicalInputTests(unittest.TestCase):
    """Keep reconstructed TLEN signs and filtering flags explicit."""

    @unittest.skipUnless(
        importlib.util.find_spec("pysam"), "Requires preparation dependency pysam"
    )
    def test_reconstructed_forward_mate_has_positive_tlen_and_preserves_filtering_flags(
        self,
    ):
        import pysam

        preparation = load_script("snyder_preparation", "data/prepare_snyder.py")
        for reverse in (False, True):
            record = pysam.AlignedSegment()
            record.query_name = "M_known_merged"
            record.query_sequence = "A" * 30
            record.reference_id, record.reference_start = 0, 100
            record.cigarstring = "30M"
            record.flag = 1024 | 512 | 256 | 2048 | (16 if reverse else 0)
            read1, read2 = preparation.reconstruct_pair(record)
            forward = read2 if reverse else read1
            backward = read1 if reverse else read2
            self.assertEqual(
                (forward.template_length, backward.template_length), (30, -30)
            )
            for mate in (read1, read2):
                self.assertEqual(
                    mate.flag & (1024 | 512 | 256 | 2048), 1024 | 512 | 256 | 2048
                )
            record.query_name = "arbitrary_singleton"
            with self.assertRaisesRegex(ValueError, "Unexpected unpaired name"):
                preparation.reconstruct_pair(record)


class MotifTests(unittest.TestCase):
    """Compare native left-end frequencies and reject damaged both-end output."""

    def test_absent_abundant_motif_retains_an_infinite_log_ratio(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, present in (("left", "AAA"), ("right", "CCC")):
                directory = root / name
                directory.mkdir()
                with (directory / "motifs.tsv").open("w") as destination:
                    for bases in product("ACGT", repeat=3):
                        motif = "".join(bases)
                        destination.write(f"{motif}\t{int(motif == present)}\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(BENCHMARKS / "ends/compare_end_motifs.py"),
                    str(root), "--k", "3",
                ],
                text=True, capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
            # Two bins differ by 1 among 64 motifs, giving TV=1 and MAE=2/64
            self.assertEqual(row["max_abs_log2_ratio_abundant"], "Inf")
            self.assertEqual(float(row["total_variation"]), 1.0)
            self.assertEqual(float(row["mae_frequency"]), 1 / 32)
            self.assertAlmostEqual(float(row["pearson_r"]), -1 / 63)
            self.assertAlmostEqual(float(row["spearman_r"]), -1 / 63)

            # Positive frequencies near the float minimum still have finite log ratios
            # log2(1 / 10^-320) is about 320 * log2(10). Binary rounding adds <1e-4
            for name, large, tiny in (("left", "AAA", "CCC"), ("right", "CCC", "AAA")):
                with (root / name / "motifs.tsv").open("w") as destination:
                    for bases in product("ACGT", repeat=3):
                        motif = "".join(bases)
                        frequency = 1.0 if motif == large else 1e-320 if motif == tiny else 0.0
                        destination.write(f"{motif}\t{frequency}\n")
            finite_result = subprocess.run(
                [sys.executable, str(BENCHMARKS / "ends/compare_end_motifs.py"),
                 str(root), "--k", "3"], text=True, capture_output=True,
            )
            self.assertEqual(finite_result.returncode, 0, finite_result.stderr)
            finite_row = next(csv.DictReader(finite_result.stdout.splitlines(), delimiter="\t"))
            self.assertAlmostEqual(
                float(finite_row["max_abs_log2_ratio_abundant"]),
                320 * math.log2(10), delta=1e-4,
            )

    @unittest.skipUnless(
        importlib.util.find_spec("pyarrow") is not None,
        "Requires native Parquet support",
    )
    def test_left_profiles_retain_native_counts_and_label_side_mirroring(self):
        import pandas as pd

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ("finaletoolkit", "cfdnapro", "freia"):
                (root / name).mkdir()
            motifs = ["".join(bases) for bases in product("ACGT", repeat=3)]
            with (root / "finaletoolkit/motifs.tsv").open("w") as destination:
                for motif in motifs:
                    destination.write(
                        f"{motif}\t{0.5 if motif in ('AAA', 'CCC') else 0}\n"
                    )
            with (root / "cfdnapro/cfdnapro-motifs.tsv").open("w") as destination:
                destination.write("motif\tn\tfraction\n")
                for motif in motifs:
                    count = {"AAA": 6, "CCC": 2}.get(motif, 0)
                    destination.write(f"{motif}\t{count}\t{count / 8}\n")
            # CountSamp repeats the motif total per length. Only CountLen is additive
            pd.DataFrame(
                {
                    "Base": ["AAA", "AAA", "CCC"],
                    "Length": [30, 40, 30],
                    "CountLen": [2, 3, 5],
                    "CountSamp": [5, 5, 5],
                    "WhichEnd": ["leftmost"] * 3,
                    "Unit": ["GRCh38"] * 3,
                }
            ).to_parquet(root / "freia/left-motifs.pq", index=False)
            result = subprocess.run(
                [
                    sys.executable,
                    str(BENCHMARKS / "ends/compare_end_motifs.py"),
                    str(root),
                    "--k",
                    "3",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = {
                (row["left_profile"], row["right_profile"]): row
                for row in csv.DictReader(result.stdout.splitlines(), delimiter="\t")
            }
            self.assertEqual(len(rows), 3)
            mirrored = rows[("finaletoolkit", "freia")]
            self.assertEqual(mirrored["right_end_selection"], "left")
            self.assertEqual(mirrored["side_comparison"], "approximate-side-mirroring")
            self.assertEqual(float(mirrored["right_counted_ends"]), 10)
            self.assertEqual(float(mirrored["total_variation"]), 0)
            native_left = rows[("cfdnapro", "freia")]
            self.assertEqual(native_left["side_comparison"], "same-sides")
            self.assertEqual(float(native_left["left_counted_ends"]), 8)
            self.assertEqual(float(native_left["right_counted_ends"]), 10)
            # .75/.25 versus .5/.5 gives TV .25 and mean error .5/64
            self.assertAlmostEqual(float(native_left["total_variation"]), 0.25)
            self.assertAlmostEqual(float(native_left["mae_frequency"]), 0.0078125)

    @unittest.skipUnless(
        importlib.util.find_spec("pyarrow") is not None,
        "Requires native Parquet support",
    )
    def test_a_right_end_table_cannot_be_labeled_as_intact_left_output(self):
        import pandas as pd

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pd.DataFrame(
                {
                    "Base": ["AAA"],
                    "Length": [30],
                    "CountLen": [4],
                    "WhichEnd": ["rightmost"],
                    "Unit": ["GRCh38"],
                }
            ).to_parquet(root / "left-motifs.pq")
            comparator = load_script("motif_comparator", "ends/compare_end_motifs.py")
            with self.assertRaisesRegex(ValueError, "leftmost"):
                comparator.load_profile(
                    root, ["".join(bases) for bases in product("ACGT", repeat=3)]
                )

    def test_existing_freia_profile_is_skipped_without_recovery(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ("finaletoolkit", "second-profile"):
                directory = root / name
                directory.mkdir()
                with (directory / "motifs.tsv").open("w") as destination:
                    for bases in product("ACGT", repeat=3):
                        motif = "".join(bases)
                        destination.write(f"{motif}\t{1 if motif == 'AAA' else 0}\n")
            (root / "freia").mkdir()
            (root / "freia/T__sample.pq").write_text(
                "Deliberately unreadable damaged aggregate"
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(BENCHMARKS / "ends/compare_end_motifs.py"),
                    str(root),
                    "--k",
                    "3",
                ],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Skipping FrEIA", result.stderr)
            rows = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(float(rows[0]["total_variation"]), 0)
            self.assertNotIn("freia", result.stdout)


if __name__ == "__main__":
    unittest.main()
