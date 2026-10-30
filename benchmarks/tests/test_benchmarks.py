"""Hand-derived benchmark regressions, left unexecuted for the repository owner.

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
    specification = importlib.util.spec_from_file_location(name, BENCHMARKS / relative_path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


timings = load_script("timing_summary", "summarize_timings.py")


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
            data = subprocess.check_output(["zstd", "-q", "-c"], input=contents.encode())
            if second_frame is not None:
                data += subprocess.check_output(["zstd", "-q", "-c"], input=second_frame.encode())
            (directory / "fcoverage.per_position.bedgraph.zst").write_bytes(data)
        else:
            filename = "positions.tsv.gz" if name in ("samtools-depth", "pandepth") else "coverage.per-base.bed.gz"
            with gzip.open(directory / filename, "wt") as destination:
                destination.write(contents)

    def compare(self, *options):
        """Call the supported CLI, preserving error output for assertions."""
        return subprocess.run([sys.executable, str(BENCHMARKS / "coverage/compare_coverage.py"),
                               str(self.profiles), *options], text=True, capture_output=True)

    def test_weights_errors_by_bases_and_reports_nonzero_union(self):
        # A difference of 1 on 10/100 bases gives MAE .1 and RMSE sqrt(.1)
        self.profile("cfdnalab", "chr1\t0\t10\t1\n")
        self.profile("mosdepth", "chr1\t0\t100\t0\n")
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(int(row["bases"]), 100)
        self.assertAlmostEqual(float(row["fraction_identical"]), .9)
        self.assertAlmostEqual(float(row["mae"]), .1)
        self.assertAlmostEqual(float(row["rmse"]), math.sqrt(.1))
        self.assertEqual(int(row["nonzero_union_bases"]), 10)
        self.assertAlmostEqual(float(row["union_mae"]), 1)
        self.assertEqual(row["pearson_r_1bp"], "NA")

    def test_residual_distribution_preserves_sign_and_counts_genomic_bases(self):
        # Residuals -2, 0, +1 cover 20, 70, 30 bases across both contigs
        # MAE=(40+30)/120=7/12, RMSE=sqrt(110/120), signed mean=-10/120
        scope = "chr1\t100\nchr2\t20\n"
        self.profile("cfdnalab", "chr1\t0\t10\t3\nchr2\t0\t20\t1\n", scope)
        self.profile("mosdepth", "chr2\t0\t20\t0\nchr1\t0\t30\t2\nchr1\t30\t100\t0\n", scope)
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
        self.assertEqual([float(row["residual"]) for row in residuals], [-.25, .5])
        self.assertEqual([int(row["bases"]) for row in residuals], [75, 25])
        summary = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertAlmostEqual(float(summary["mae"]), .3125)
        self.assertEqual(summary["pearson_r_1bp"], "NA")

    def test_residual_counts_survive_large_value_sets_and_multiple_pairs(self):
        # Depths 0..4096 each cover a base, then depth 1 covers ten more
        # The repeated residual must accumulate after a batch spills to disk
        scope = "chr1\t4107\n"
        self.profile("cfdnalab", "", scope)
        depths = "".join(f"chr1\t{position}\t{position + 1}\t{position}\n" for position in range(4097))
        self.profile("mosdepth", depths + "chr1\t4097\t4107\t1\n", scope)
        self.profile("pandepth", "".join(f"chr1\t{position}\t0\n" for position in range(4107)), scope)
        output = self.profiles.parent / "residuals.tsv"

        result = self.compare()

        self.assertEqual(result.returncode, 0, result.stderr)
        comparisons = {}
        with output.open() as source:
            for row in csv.DictReader(source, delimiter="\t"):
                comparisons.setdefault((row["left_profile"], row["right_profile"]), []).append(row)
        self.assertEqual(len(comparisons), 3)
        negative = comparisons[("cfdnalab", "mosdepth")]
        zero = comparisons[("cfdnalab", "pandepth")]
        positive = comparisons[("mosdepth", "pandepth")]
        self.assertEqual([float(row["residual"]) for row in negative], list(range(-4096, 1)))
        self.assertEqual([int(row["bases"]) for row in negative], [1] * 4095 + [11, 1])
        self.assertEqual([float(row["residual"]) for row in positive], list(range(4097)))
        self.assertEqual([int(row["bases"]) for row in positive], [1, 11] + [1] * 4095)
        self.assertEqual(len(zero), 1)
        self.assertEqual(float(zero[0]["residual"]), 0)
        self.assertEqual(int(zero[0]["bases"]), 4107)
        for residuals in comparisons.values():
            self.assertEqual(sum(int(row["bases"]) for row in residuals), 4107)
            self.assertAlmostEqual(sum(float(row["fraction_of_bases"]) for row in residuals), 1)

    def test_tile_boundaries_contig_order_and_concatenated_frames_are_equivalent(self):
        scope = "chr1\t100\nchr2\t80\n"
        self.profile("cfdnalab", "chr1\t0\t30\t2\n", scope,
                     second_frame="chr1\t30\t100\t2\nchr2\t0\t80\t3\n")
        self.profile("mosdepth", "chr2\t0\t80\t3\nchr1\t0\t100\t2\n", scope)
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(float(row["fraction_identical"]), 1)
        self.assertEqual(float(row["mae"]), 0)
        self.assertEqual(float(row["left_integrated_depth"]), 440)

    def test_position_origins_agree(self):
        self.profile("samtools-depth", "".join(f"chr1\t{position}\t2\n" for position in range(1, 101)))
        self.profile("pandepth", "".join(f"chr1\t{position}\t2\n" for position in range(100)))
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
        self.profile("cfdnalab", "chr1\t0\t10000\t1\nchr1\t10000\t15000\t3\n", "chr1\t15000\n")
        self.profile("mosdepth", "chr1\t0\t10000\t2\nchr1\t10000\t13000\t6\nchr1\t13000\t15000\t0\n", "chr1\t15000\n")
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

    def test_regional_comparison_checks_coordinates_and_raw_counts(self):
        # Counts (2, 1) versus (1, 3) give differences (1, -2), MAE 1.5
        (self.profiles.parent / "regions.bed").write_text("chr1\t0\t50\t.\nchr1\t50\t100\t.\n")
        left = self.profiles / "cfdnalab-region-counts"
        left.mkdir()
        table = "chrom\tstart\tend\tcount_30_1001\nchr1\t0\t50\t2\nchr1\t50\t100\t1\n"
        (left / "length_counts.tsv.zst").write_bytes(subprocess.check_output(["zstd", "-q", "-c"], input=table.encode()))
        right = self.profiles / "finaletoolkit-coverage"
        right.mkdir()
        with gzip.open(right / "finaletoolkit.coverage.bed.gz", "wt") as destination:
            destination.write("chr1\t50\t100\t.\t3.0\nchr1\t0\t50\t.\t1.0\n")
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        row = next(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
        self.assertEqual(float(row["mae"]), 1.5)
        self.assertAlmostEqual(float(row["rmse"]), math.sqrt(2.5))
        self.assertEqual(float(row["regional_total_ratio"]), .75)
        self.assertEqual(float(row["pearson_r_regions"]), -1)
        # Regional-only comparisons keep the automatic residual file header-only
        with (self.profiles.parent / "residuals.tsv").open() as source:
            self.assertEqual(list(csv.DictReader(source, delimiter="\t")), [])
        # A shifted endpoint must fail rather than align by row number
        with gzip.open(right / "finaletoolkit.coverage.bed.gz", "wt") as destination:
            destination.write("chr1\t0\t51\t.\t1\nchr1\t50\t100\t.\t3\n")
        self.assertNotEqual(self.compare().returncode, 0)


class TimingTests(unittest.TestCase):
    """Use literal timing/accounting fixtures rather than invoking Slurm."""

    def test_join_preserves_correct_step_memory_and_missing_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "timings.tsv"
            path.write_text("run_order\tjob_id\tstep_id\n1\t123\t0\n2\t123\t1\n")
            rows = timings.read_runs([path], {"123.0": {"NTasks": "1", "MaxRSS": "2048K", "State": "COMPLETED", "ExitCode": "0:0"}})
            self.assertEqual(rows[0]["slurm_max_rss_kib"], 2048)
            self.assertEqual(rows[1]["slurm_max_rss_kib"], "NA")

    def test_sample_statistics_use_sample_sd_and_inclusive_quartiles(self):
        self.assertEqual(timings.distribution([1, 3, 5]), [3, 3, 2, 3, 2, 4])
        self.assertEqual(timings.distribution([4]), [1, 4, "NA", 4, 4, 4])

    def test_two_method_budget_includes_all_seven_invocations(self):
        columns = "run_order repetition warmup pilot method tool_version cores node wall_seconds user_seconds system_seconds flush_seconds process_plus_flush_seconds output_bytes max_rss_kib tool_exit_code step_exit_code job_id step_id".split()
        values = ["1 1 false 1 cfdnalab checkout 1 node 720 700 20 0 720 100 200 0 0 123 0".split(),
                  "2 1 false 1 slow release 1 node 7200 7000 200 0 7200 100 200 0 0 123 1".split()]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "pilot.tsv"
            with path.open("w") as destination:
                writer = csv.writer(destination, delimiter="\t")
                writer.writerow(columns)
                writer.writerows(values)
            result = subprocess.run([sys.executable, str(BENCHMARKS / "summarize_timings.py"), str(path), "--estimate"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            # 7*(720+7200) = 55440 seconds = 15.4 hours
            self.assertIn("15.40 h raw", result.stderr)
            self.assertIn("19.25 h with margin", result.stderr)

    def test_accounting_oom_excludes_an_apparently_successful_native_run(self):
        self.assertFalse(timings.successful({"tool_exit_code": "0", "step_exit_code": "0",
                                           "slurm_state": "OUT_OF_MEMORY", "slurm_exit_code": "0:125"}))


class RunControlTests(unittest.TestCase):
    """Check scheduling and Picard's unusual version status without tool installs."""

    def test_each_job_runs_all_selected_methods_at_a_single_core_setting(self):
        script = '''source "$BENCH_COMMON"
benchmark_controls
BENCH_FEATURE=lengths
BENCH_WARMUPS=2
BENCH_REPETITIONS=5
methods_to_run=(cfdnalab cfdnapro finaletoolkit)
benchmark_schedule
'''
        for cores in (1, 4, 12, 32):
            with self.subTest(cores=cores):
                result = subprocess.run(["bash", "-euc", script], env={
                    **os.environ, "BENCH_COMMON": str(BENCHMARKS / "common.sh"),
                    "BENCH_CORES": str(cores), "BENCH_PILOT": "0", "BENCH_EXPECTED_NODE": "",
                    "BENCH_SETUP_PYTHON": sys.executable,
                }, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                rows = list(csv.reader(result.stdout.splitlines(), delimiter="\t"))
                # Three selected methods each get two warmups and five measurements
                self.assertEqual(len(rows), 21)
                self.assertEqual([row[0] for row in rows], ["cfdnalab"] * 7 + ["cfdnapro"] * 7 + ["finaletoolkit"] * 7)
                self.assertEqual(sum(row[3] == "true" for row in rows), 6)
                self.assertEqual({row[1] for row in rows}, {str(cores)})

    def test_multiple_core_settings_are_rejected_for_a_single_job(self):
        result = subprocess.run(["bash", "-euc", 'source "$BENCH_COMMON"; benchmark_controls'], env={
            **os.environ, "BENCH_COMMON": str(BENCHMARKS / "common.sh"), "BENCH_CORES": "1 4",
            "BENCH_PILOT": "0", "BENCH_EXPECTED_NODE": "", "BENCH_SETUP_PYTHON": sys.executable,
        }, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("single positive core count", result.stderr)

    def test_picard_accepts_only_expected_status_and_version(self):
        with tempfile.TemporaryDirectory() as temporary:
            executable = Path(temporary) / "picard"
            for version, status, succeeds in (("Version:3.5.0", 1, True), ("Version:3.5.0", 2, False), ("Version:3.4.0", 1, False), ("error opening picard-3.5.0.jar", 1, False), ("", 1, False)):
                executable.write_text(f"#!/usr/bin/env bash\nprintf '%s\\n' '{version}'\nexit {status}\n")
                executable.chmod(0o755)
                result = subprocess.run(["bash", "-euc", 'source "$BENCH_COMMON"; benchmark_picard_version "$PICARD_FIXTURE"'],
                                        env={**os.environ, "BENCH_COMMON": str(BENCHMARKS / "common.sh"), "PICARD_FIXTURE": str(executable)}, capture_output=True)
                self.assertEqual(result.returncode == 0, succeeds)

    def test_installer_resumes_and_never_completes_a_failed_reinstall(self):
        # A fake Conda boundary permits info but fails installation commands
        # Each reinstall reason must reach that failure and remain outdated
        for reason in ("complete", "outdated", "version", "missing", "force"):
            with self.subTest(reason=reason), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                installation, fake_bin, conda_base = root / "installation", root / "bin", root / "conda"
                fake_bin.mkdir()
                (conda_base / "etc/profile.d").mkdir(parents=True)
                (conda_base / "etc/profile.d/conda.sh").write_text("# Conda fixture\n")
                executable = fake_bin / "conda"
                executable.write_text('#!/usr/bin/env bash\nif [[ "$1" == info ]]; then printf "%s\\n" "$BENCH_FAKE_CONDA_BASE"; else exit 9; fi\n')
                executable.chmod(0o755)
                (installation / "envs/comparison/bin").mkdir(parents=True)
                if reason != "missing":
                    (installation / "envs/comparison/bin/python").touch()
                version = "old" if reason == "version" else "checkout"
                status = "outdated" if reason == "outdated" else "complete"
                manifest = installation / "manifest.tsv"
                manifest.write_text(f"tool\tversion\tstatus\tartifacts\ncomparison\t{version}\t{status}\tenvs/comparison/bin/python\n")
                command = ["bash", str(BENCHMARKS / "install_tools.sh"), str(installation), "--only", "comparison"]
                if reason == "force":
                    command += ["--force", "comparison"]
                result = subprocess.run(command, env={**os.environ, "PATH": f"{fake_bin}:{os.environ['PATH']}",
                                                       "BENCH_FAKE_CONDA_BASE": str(conda_base)}, capture_output=True, text=True)
                if reason == "complete":
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn("Keeping comparison checkout", result.stdout)
                else:
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("comparison\tcheckout\toutdated\t", manifest.read_text())


class BiologicalInputTests(unittest.TestCase):
    """Keep reconstructed TLEN signs and filtering flags explicit."""

    @unittest.skipUnless(importlib.util.find_spec("pysam"), "Requires preparation dependency pysam")
    def test_reconstructed_forward_mate_has_positive_tlen_and_preserves_filtering_flags(self):
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
            self.assertEqual((forward.template_length, backward.template_length), (30, -30))
            for mate in (read1, read2):
                self.assertEqual(mate.flag & (1024 | 512 | 256 | 2048), 1024 | 512 | 256 | 2048)
            record.query_name = "arbitrary_singleton"
            with self.assertRaisesRegex(ValueError, "Unexpected unpaired name"):
                preparation.reconstruct_pair(record)


class MotifTests(unittest.TestCase):
    """The damaged FrEIA aggregate must never enter concordance."""

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
            (root / "freia/T__sample.pq").write_text("Deliberately unreadable damaged aggregate")
            result = subprocess.run([sys.executable, str(BENCHMARKS / "ends/compare_end_motifs.py"), str(root), "--k", "3"], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Skipping FrEIA", result.stderr)
            rows = list(csv.DictReader(result.stdout.splitlines(), delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(float(rows[0]["total_variation"]), 0)
            self.assertNotIn("freia", result.stdout)


if __name__ == "__main__":
    unittest.main()
