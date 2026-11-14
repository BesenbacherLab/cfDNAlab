"""Regression tests for submission, setup checks, concordance, pilot budgets and native formats.

Expected numerical values are manually derived from the synthetic inputs.

Native histogram readers require the comparison environment. No Slurm allocation
or real BAM is used. Measurement and the job scripts are tested in test_runners.py.
"""

import csv
from contextlib import redirect_stdout
import gzip
import importlib.util
from io import StringIO
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


BENCHMARKS = Path(__file__).resolve().parents[1]
EXCLUSION_FILENAMES = (
    "encode_unified_GRCh38_exclusion_list.bed", "hg38_centromeres.bed", "hg38_gaps.bed",
    "hg38_fix_patches.bed", "hg38_alternative_haplotypes.bed",
)


class SubmissionTests(unittest.TestCase):
    """Check allocation requests without contacting Slurm or running a benchmark."""

    def test_only_production_submissions_request_the_whole_node(self):
        for feature, job_filename in (
            ("coverage", "benchmark_coverage.sh"),
            ("lengths", "benchmark_lengths.sh"),
            ("ends", "benchmark_end_motifs.sh"),
            ("midpoints", "benchmark_midpoints.sh"),
        ):
            # Unconditional header options would also apply to pilots
            header = (BENCHMARKS / feature / job_filename).read_text().splitlines()
            self.assertNotIn("#SBATCH --exclusive", header)
            self.assertNotIn("#SBATCH --mem=0", header)
            self.assertNotIn("#SBATCH --cpus-per-task=32", header)
            for pilot, memory_override in (
                (None, None), ("0", None), ("1", None), ("1", "80G"), ("invalid", None)
            ):
                with (
                    self.subTest(feature=feature, pilot=pilot, memory_override=memory_override),
                    tempfile.TemporaryDirectory() as temporary,
                ):
                    root = Path(temporary)
                    call_log = root / "sbatch-arguments"
                    pilot_log = root / "pilot"
                    sbatch = root / "sbatch"
                    sbatch.write_text(
                        '#!/usr/bin/env bash\n'
                        'printf "%s\\n" "$@" > "$SUBMISSION_LOG"\n'
                        'printf "%s\\n" "$BENCH_PILOT" > "$PILOT_LOG"\n'
                        'exit 7\n'
                    )
                    sbatch.chmod(0o755)
                    environment = {
                        **os.environ,
                        "PATH": str(root) + ":" + os.environ["PATH"],
                        "BENCH_CORES": "4",
                        "SUBMISSION_LOG": str(call_log),
                        "PILOT_LOG": str(pilot_log),
                    }
                    if feature == "midpoints":
                        references = root / "references"
                        references.mkdir()
                        for exclusion_filename in EXCLUSION_FILENAMES:
                            (references / exclusion_filename).touch()
                        # The launcher's format sample accepts a 2001 bp BED4
                        # without requiring a BAM or any installed native tool
                        (root / "sites.bed").write_text("chr1\t1000\t3001\tTFBS\n")
                        environment.update({
                            "BENCH_INTERVALS": str(root / "sites.bed"),
                            "BENCH_RESULTS": str(root / "results"),
                            "BENCH_SETUP_PYTHON": sys.executable,
                            "BENCH_METHODS": "cfdnalab",
                            "BENCH_GRIFFIN_REF_DIR": str(references),
                        })
                    environment.pop("BENCH_PILOT", None)
                    if pilot is not None:
                        environment["BENCH_PILOT"] = pilot
                    extra_options = ["--time=00:10:00", "--comment=benchmark debug"]
                    if memory_override is not None:
                        extra_options.append("--mem=" + memory_override)
                    result = subprocess.run(
                        [
                            "bash", str(BENCHMARKS / "submit_benchmark.sh"), feature,
                            *extra_options,
                        ],
                        cwd=BENCHMARKS.parent,
                        env=environment,
                        text=True,
                        capture_output=True,
                    )
                    if pilot == "invalid":
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn("BENCH_PILOT", result.stderr)
                        self.assertFalse(call_log.exists())
                        continue
                    # The launcher returns sbatch's status and preserves quoted arguments
                    self.assertEqual(result.returncode, 7, result.stderr)
                    expected = ["--export=ALL", "--cpus-per-task=4"]
                    if pilot != "1":
                        expected += ["--exclusive", "--mem=0"]
                    else:
                        expected += ["--mem=50G"]
                    expected += [
                        *extra_options,
                        f"benchmarks/{feature}/{job_filename}",
                    ]
                    self.assertEqual(call_log.read_text().splitlines(), expected)
                    self.assertEqual(pilot_log.read_text(), (pilot or "0") + "\n")


class AutomaticConcordanceTests(unittest.TestCase):
    """Compare saved profiles after timing without native tools or Slurm."""

    def test_coverage_concordance_and_histogram_are_saved_for_pilots_and_full_runs(self):
        for pilot in ("0", "1"):
            with self.subTest(pilot=pilot), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                results = root / "results"
                for method in ("cfdnalab", "mosdepth"):
                    directory = results / "profiles" / method
                    directory.mkdir(parents=True)
                    (directory / "chrom.sizes").write_text("chr1\t10\n")
                (results / "profiles/cfdnalab/fcoverage.per_position.bedgraph.zst").write_text(
                    "chr1\t0\t5\t2\n"
                )
                with gzip.open(results / "profiles/mosdepth/coverage.per-base.bed.gz", "wt") as output:
                    output.write("chr1\t0\t5\t3\nchr1\t5\t10\t0\n")
                # The comparison environment must supply its own zstd on PATH
                # Only decompression is faked. The real comparator calculates the metrics
                comparison_bin = root / "comparison environment" / "bin"
                comparison_bin.mkdir(parents=True)
                comparison_python = comparison_bin / "python"
                comparison_python.write_text('#!/usr/bin/env bash\nexec "$TEST_PYTHON" "$@"\n')
                comparison_python.chmod(0o755)
                zstd = comparison_bin / "zstd"
                zstd.write_text(
                    '#!/usr/bin/env bash\n'
                    'printf "comparison zstd used\\n" > "$ZSTD_LOG"\n'
                    'cat "$2"\n'
                )
                zstd.chmod(0o755)
                timings = results / "timings.tsv"
                timings.write_text("native timing records remain unchanged\n")
                result = subprocess.run(
                    ["bash", "-euc", 'source "$BENCH_COMMON"\ncompare_saved_profiles "$COMPARATOR" --common-scope'],
                    env={
                        **os.environ,
                        "BENCH_COMMON": str(BENCHMARKS / "common.sh"),
                        "BENCH_COMPARE_PYTHON": str(comparison_python),
                        "BENCH_RESULTS": str(results),
                        "BENCH_PILOT": pilot,
                        "COMPARATOR": str(BENCHMARKS / "coverage/compare_coverage.py"),
                        "TEST_PYTHON": sys.executable,
                        "ZSTD_LOG": str(root / "zstd-log"),
                        "TMPDIR": temporary,
                    },
                    text=True,
                    capture_output=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                with (results / "concordance.tsv").open() as source:
                    rows = list(csv.DictReader(source, delimiter="\t"))
                # Five bases differ by 1 and five match. MAE = 5 / 10 = 0.5
                self.assertEqual(len(rows), 1)
                self.assertEqual(float(rows[0]["mae"]), 0.5)
                with (results / "residuals.tsv").open() as source:
                    residuals = list(csv.DictReader(source, delimiter="\t"))
                self.assertEqual(
                    [(float(row["residual"]), int(row["bases"])) for row in residuals],
                    [(-1.0, 5), (0.0, 5)],
                )
                self.assertEqual(timings.read_text(), "native timing records remain unchanged\n")
                self.assertEqual((root / "zstd-log").read_text(), "comparison zstd used\n")

    def test_only_profiles_of_the_same_quantity_form_a_comparable_pair(self):
        # Per-base read depth, per-base fragment coverage and regional fragment
        # counts are three different quantities, so the coverage job skips
        # concordance. Adding pandepth gives a second read-depth profile
        for methods, expected in (
            (("cfdnalab", "mosdepth-fragment", "finaletoolkit-coverage"), "no pair"),
            (("cfdnalab", "mosdepth-fragment", "finaletoolkit-coverage", "pandepth"), "pair"),
        ):
            with self.subTest(methods=methods), tempfile.TemporaryDirectory() as temporary:
                results = Path(temporary)
                for method in methods:
                    (results / "profiles" / method).mkdir(parents=True)
                result = subprocess.run(
                    [
                        "bash", "-euc",
                        'source "$BENCH_COMMON"\nsource "$COVERAGE_COMMANDS"\n'
                        'if saved_profiles_share_a_quantity; then echo pair; else echo "no pair"; fi',
                    ],
                    env={
                        **os.environ,
                        "BENCH_COMMON": str(BENCHMARKS / "common.sh"),
                        "COVERAGE_COMMANDS": str(BENCHMARKS / "coverage/commands.sh"),
                        "BENCH_RESULTS": str(results),
                    },
                    text=True,
                    capture_output=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout, expected + "\n")

    def test_comparison_failure_preserves_previous_metrics_and_prints_the_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            for method in ("cfdnalab", "samtools"):
                (results / "profiles" / method).mkdir(parents=True)
            (results / "concordance.tsv").write_text("previous successful metrics\n")
            comparator = results / "broken-comparator.py"
            comparator.write_text('print("partial output")\nraise ValueError("invalid profile")\n')
            result = subprocess.run(
                ["bash", "-euc", 'source "$BENCH_COMMON"\ncompare_saved_profiles "$COMPARATOR"'],
                env={
                    **os.environ,
                    "BENCH_COMMON": str(BENCHMARKS / "common.sh"),
                    "BENCH_COMPARE_PYTHON": sys.executable,
                    "BENCH_RESULTS": str(results),
                    "COMPARATOR": str(comparator),
                },
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("invalid profile", result.stderr)
            self.assertIn("invalid profile", (results / "concordance-error.log").read_text())
            self.assertEqual((results / "concordance.tsv").read_text(), "previous successful metrics\n")
            self.assertFalse((results / ".concordance.tsv").exists())


class SetupAndScopeTests(unittest.TestCase):
    """Reject incomplete setup before staging or decompressing native data."""

    def test_each_runner_uses_the_compute_nodes_tmpdir_before_staging(self):
        for feature, filename, directory_prefix in (
            ("coverage", "benchmark_coverage.sh", "coverage"),
            ("lengths", "benchmark_lengths.sh", "lengths"),
            ("ends", "benchmark_end_motifs.sh", "end-motifs"),
            ("midpoints", "benchmark_midpoints.sh", "midpoints"),
        ):
            for scratch_available in (True, False):
                with (
                    self.subTest(feature=feature, scratch_available=scratch_available),
                    tempfile.TemporaryDirectory() as temporary,
                ):
                    root = Path(temporary)
                    scratch = root / "node scratch"
                    scratch.mkdir()
                    fake_bin = root / "bin"
                    fake_bin.mkdir()
                    call_log = root / "mktemp-arguments"
                    # Stop at directory creation so no BAM, tool or Slurm command runs
                    mktemp = fake_bin / "mktemp"
                    mktemp.write_text(
                        '#!/usr/bin/env bash\n'
                        'printf "%s\\n" "$@" > "$SCRATCH_CALL_LOG"\n'
                        'exit 9\n'
                    )
                    mktemp.chmod(0o755)
                    environment = {
                        **os.environ,
                        "PATH": str(fake_bin) + ":" + os.environ["PATH"],
                        "SLURM_SUBMIT_DIR": str(BENCHMARKS.parent),
                        "SLURM_JOB_ID": "123",
                        "SLURM_JOB_NUM_NODES": "1",
                        "SLURM_CPUS_PER_TASK": "4",
                        "BENCH_BAM": "/unused/input.bam",
                        "BENCH_INDEX": "/unused/input.bam.bai",
                        "BENCH_INTERVALS": "/unused/sites.bed",
                        "BENCH_GRIFFIN_REF_DIR": "/unused/griffin-references",
                        "BENCH_SAMTOOLS_ARCHIVE": "/unused/samtools.tar.gz",
                        "BENCH_CFDNA": "/unused/cfdna",
                        "BENCH_REF_2BIT": "/unused/reference.2bit",
                        "BENCH_RESULTS": str(root / "results"),
                        "BENCH_CORES": "4",
                        "BENCH_METHODS": "cfdnalab-reference" if feature == "ends" else "cfdnalab",
                        "BENCH_PILOT": "1",
                        "BENCH_EXPECTED_NODE": "",
                        "BENCH_SETUP_PYTHON": sys.executable,
                        "BENCH_RUN_TIMEOUT_SECONDS": "0",
                        "SCRATCH_CALL_LOG": str(call_log),
                    }
                    environment.pop("BENCH_LOCAL_ROOT", None)
                    environment.pop("TMPDIR", None)
                    if scratch_available:
                        environment["TMPDIR"] = str(scratch)
                    result = subprocess.run(
                        ["bash", str(BENCHMARKS / feature / filename)],
                        env=environment,
                        text=True,
                        capture_output=True,
                    )
                    if scratch_available:
                        self.assertEqual(result.returncode, 9, result.stderr)
                        self.assertEqual(
                            call_log.read_text().splitlines(),
                            ["-d", str(scratch / f"{directory_prefix}-123.XXXXXX")],
                        )
                    else:
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn("TMPDIR", result.stderr)
                        self.assertFalse(call_log.exists())
                    self.assertFalse((root / "results").exists())
                    self.assertEqual(list(scratch.iterdir()), [])

    def test_each_runner_requires_an_explicit_method_block_before_staging(self):
        for feature, filename in (
            ("coverage", "benchmark_coverage.sh"),
            ("lengths", "benchmark_lengths.sh"),
            ("ends", "benchmark_end_motifs.sh"),
            ("midpoints", "benchmark_midpoints.sh"),
        ):
            with (
                self.subTest(feature=feature),
                tempfile.TemporaryDirectory() as temporary,
            ):
                environment = {
                    **os.environ,
                    "SLURM_SUBMIT_DIR": str(BENCHMARKS.parent),
                    "SLURM_JOB_ID": "123",
                    "BENCH_CORES": "4",
                    "BENCH_BAM": "/unused/input.bam",
                    "BENCH_INDEX": "/unused/input.bam.bai",
                    "BENCH_INTERVALS": "/unused/sites.bed",
                    "BENCH_SAMTOOLS_ARCHIVE": "/unused/samtools.tar.gz",
                    "TMPDIR": temporary,
                    "BENCH_RESULTS": temporary + "/results",
                }
                environment.pop("BENCH_METHODS", None)
                result = subprocess.run(
                    ["bash", str(BENCHMARKS / feature / filename)],
                    env=environment,
                    text=True,
                    capture_output=True,
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("BENCH_METHODS", result.stderr)
                self.assertFalse((Path(temporary) / "results").exists())

    def test_griffin_requires_a_job_of_its_own_before_inputs_are_staged(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = subprocess.run(
                ["bash", str(BENCHMARKS / "midpoints/benchmark_midpoints.sh")],
                env={
                    **os.environ,
                    "SLURM_SUBMIT_DIR": str(BENCHMARKS.parent),
                    "SLURM_JOB_ID": "123",
                    "SLURM_JOB_NUM_NODES": "1",
                    "SLURM_CPUS_PER_TASK": "4",
                    "BENCH_METHODS": "cfdnalab griffin-gc",
                    # Keep the GC range valid so this reaches the job-isolation check
                    "BENCH_MIN_LENGTH": "100",
                    "BENCH_MAX_LENGTH": "200",
                    "BENCH_MIN_MAPQ": "20",
                    "BENCH_COMPARE_PYTHON": sys.executable,
                    "BENCH_CORES": "4",
                    "BENCH_PILOT": "0",
                    "BENCH_BAM": "/unused/input.bam",
                    "BENCH_INTERVALS": "/unused/sites.bed",
                    "TMPDIR": temporary,
                    "BENCH_RESULTS": temporary + "/results",
                    "BENCH_EXPECTED_NODE": "",
                    "BENCH_SETUP_PYTHON": sys.executable,
                    "BENCH_RUN_TIMEOUT_SECONDS": "0",
                },
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Run Griffin in its own job", result.stderr)
            self.assertEqual(list(Path(temporary).iterdir()), [])

    def test_scope_mismatch_fails_before_track_loading_and_keeps_previous_histogram(
        self,
    ):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles = root / "profiles"
            for name, scope in (
                ("cfdnalab", "chr1\t100\nchr2\t100\n"),
                ("mosdepth", "chr1\t100\n"),
            ):
                directory = profiles / name
                directory.mkdir(parents=True)
                (directory / "chrom.sizes").write_text(scope)
            # No track files exist: scope validation must precede decompression
            (root / "residuals.tsv").write_text("previous successful histogram\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(BENCHMARKS / "coverage/compare_coverage.py"),
                    str(profiles),
                ],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Scopes differ", result.stderr)
            self.assertEqual(
                (root / "residuals.tsv").read_text(), "previous successful histogram\n"
            )

    def test_invalid_track_after_scope_preflight_keeps_previous_histogram(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles = root / "profiles"
            for name in ("samtools-depth", "pandepth"):
                directory = profiles / name
                directory.mkdir(parents=True)
                (directory / "chrom.sizes").write_text("chr1\t100\n")
                with gzip.open(directory / "positions.tsv.gz", "wt") as destination:
                    # A single covered base cannot be a complete 100 bp dense track
                    destination.write(
                        "chr1\t1\t2\n" if name == "samtools-depth" else "chr1\t0\t2\n"
                    )
            (root / "residuals.tsv").write_text("previous successful histogram\n")
            result = subprocess.run(
                [
                    sys.executable,
                    str(BENCHMARKS / "coverage/compare_coverage.py"),
                    str(profiles),
                ],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Incomplete native track", result.stderr)
            self.assertEqual(
                (root / "residuals.tsv").read_text(), "previous successful histogram\n"
            )
            self.assertEqual(list(root.glob(".residuals-*")), [])


class PilotBudgetTests(unittest.TestCase):
    """Keep unavailable runtimes unknown while reporting useful successful budgets."""

    def test_partial_and_all_failed_pilots_report_unknown_budgets(self):
        for include_success in (True, False):
            with (
                self.subTest(include_success=include_success),
                tempfile.TemporaryDirectory() as temporary,
            ):
                path = Path(temporary) / "pilot.tsv"
                base = dict(
                    run_order="1",
                    repetition="1",
                    warmup="false",
                    pilot="1",
                    method="slow",
                    tool_version="pinned",
                    cores="4",
                    node="node",
                    job_id="123",
                    step_id="0",
                    wall_seconds="NA",
                    user_seconds="NA",
                    system_seconds="NA",
                    flush_seconds="NA",
                    process_plus_flush_seconds="NA",
                    max_rss_kib="NA",
                    output_bytes="NA",
                    tool_exit_code="124",
                    step_exit_code="124",
                )
                rows = [base]
                if include_success:
                    rows.append(
                        {
                            **base,
                            "run_order": "2",
                            "step_id": "1",
                            "method": "cfdnalab",
                            "wall_seconds": "60",
                            "user_seconds": "50",
                            "system_seconds": "10",
                            "flush_seconds": "0",
                            "process_plus_flush_seconds": "60",
                            "max_rss_kib": "10",
                            "output_bytes": "100",
                            "tool_exit_code": "0",
                            "step_exit_code": "0",
                        }
                    )
                with path.open("w") as destination:
                    writer = csv.DictWriter(
                        destination, fieldnames=list(base), delimiter="\t"
                    )
                    writer.writeheader()
                    writer.writerows(rows)
                result = subprocess.run(
                    [
                        sys.executable,
                        str(BENCHMARKS / "summarize_timings.py"),
                        str(path),
                        "--estimate",
                    ],
                    text=True,
                    capture_output=True,
                )
                self.assertNotEqual(result.returncode, 0)
                report = {
                    row["method"]: row
                    for row in csv.DictReader(
                        result.stdout.splitlines(), delimiter="\t"
                    )
                }
                self.assertEqual(report["slow"]["budget_minutes"], "NA")
                self.assertIn("incomplete", result.stderr.lower())
                self.assertNotIn("Traceback", result.stderr)
                if include_success:
                    # Four 60-second invocations with 25% margin need five minutes
                    self.assertAlmostEqual(
                        float(report["cfdnalab"]["budget_minutes"]), 5
                    )


@unittest.skipUnless(
    all(
        importlib.util.find_spec(name) is not None
        for name in ("cfdnalab", "numpy", "pandas")
    ),
    "Requires the benchmark comparison environment",
)
class NativeComparatorTests(unittest.TestCase):
    """Check hand-derived native-format boundaries without running native tools."""

    def comparator(self, relative_path):
        """Load the standalone comparator for a small native-format fixture."""
        specification = importlib.util.spec_from_file_location(
            "review_comparator", BENCHMARKS / relative_path
        )
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        return module

    def test_length_normalization_rejects_an_overflowed_total(self):
        import numpy as np

        comparator = self.comparator("lengths/compare_lengths.py")
        # Finite bins must not sum to infinity and then normalize into two zeros
        with self.assertRaisesRegex(RuntimeError, "Invalid fragment total"):
            comparator.summarize_profile(
                np.array([1e308, 1e308]), np.array([30, 31])
            )

    def test_motif_normalization_rejects_an_overflowed_total(self):
        comparator = self.comparator("ends/compare_end_motifs.py")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory / "cfdnapro-motifs.tsv").write_text(
                "motif\tn\nAAA\t1e308\nCCC\t1e308\n"
            )
            with self.assertRaisesRegex(RuntimeError, "Invalid motif total"):
                comparator.load_profile(directory, ["AAA", "CCC"])

    def test_midpoint_normalization_rejects_an_overflowed_mean(self):
        import numpy as np

        comparator = self.comparator("midpoints/compare_midpoints.py")
        # Every count is finite, but the sum used by the mean exceeds float range
        with (
            mock.patch.object(
                comparator, "read_cfdnalab",
                return_value=({"group": (np.full(666, 1e308), 2)}, 9990, 15),
            ),
            mock.patch.object(
                comparator, "read_griffin",
                return_value={"group": (np.ones(132), 2, 1.0)},
            ),
            mock.patch.object(sys, "argv", ["compare_midpoints.py", "profiles"]),
            redirect_stdout(StringIO()),
        ):
            with self.assertRaisesRegex(RuntimeError, "Invalid midpoint mean"):
                comparator.main()

    def test_corrected_midpoint_comparison_selects_corrected_profiles(self):
        import numpy as np

        comparator = self.comparator("midpoints/compare_midpoints.py")
        # Weighted counts average to 2, giving normalized halves 0.5 and 1.5
        # Griffin's halves are 0.75 and 1.25, so every absolute residual is 0.25
        output = StringIO()
        with (
            mock.patch.object(comparator, "read_cfdnalab", return_value=(
                {"TFBS": (np.array([1.0] * 333 + [3.0] * 333), 2)}, 9990, 15,
            )) as cfdnalab_reader,
            mock.patch.object(comparator, "read_griffin", return_value={
                "TFBS": (np.array([0.75] * 66 + [1.25] * 66), 2, 0.5),
            }) as griffin_reader,
            mock.patch.object(sys, "argv", ["compare_midpoints.py", "profiles", "--gc-corrected"]),
            redirect_stdout(output),
        ):
            comparator.main()
        cfdnalab_reader.assert_called_once_with(Path("profiles/cfdnalab-gc/midpoint_profiles.zarr"))
        griffin_reader.assert_called_once_with(Path("profiles/griffin-gc/sample.GC_corrected.coverage.tsv"), 9990, 15)
        row, = csv.DictReader(StringIO(output.getvalue()), delimiter="\t")
        self.assertEqual(row["gc_correction"], "corrected")
        self.assertEqual(row["status"], "ok")
        self.assertAlmostEqual(float(row["pearson_r"]), 1.0)
        for metric in ("mae_mean_one", "rmse_mean_one", "max_abs_mean_one"):
            self.assertAlmostEqual(float(row[metric]), 0.25)
        self.assertAlmostEqual(float(row["cfdnalab_mean_midpoints_per_site_bp"]), 1.0)

    def test_samtools_tails_do_not_enter_the_requested_histogram(self):
        import numpy as np

        comparator = self.comparator("lengths/compare_lengths.py")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "samtools"
            directory.mkdir()
            (directory / "stats.txt").write_text(
                "IS\t29\t999\t999\t0\t0\nIS\t30\t1\t1\t0\t0\n"
                "IS\t1000\t2\t2\t0\t0\nIS\t1001\t999\t999\t0\t0\n"
            )
            counts = comparator.read_profile(directory, np.arange(30, 1001))
            self.assertEqual(counts[0], 1)
            self.assertEqual(counts[-1], 2)
            self.assertEqual(counts.sum(), 3)

    def test_midpoint_comparison_keeps_smoothing_undershoots_and_native_site_counts(self):
        import numpy as np
        import pandas as pd

        comparator = self.comparator("midpoints/compare_midpoints.py")
        # Both profiles have mean 1. Their halves differ by exactly 0.25,
        # including negative smoothing values. Pearson r and r squared are 1
        counts = np.array([-0.5] * 333 + [2.5] * 333)
        native_result = SimpleNamespace(
            length_bins=lambda: [(100, 201)],
            positions=lambda: pd.DataFrame({
                "position_bin_start_bp": range(0, 9990, 15), "position_bin_end_bp": range(15, 9991, 15),
            }),
            counts_array=lambda: counts[None, None, :],
            group_metadata=lambda: pd.DataFrame({"group_name": ["TFBS"], "eligible_intervals": [1]}),
        )
        with tempfile.TemporaryDirectory() as temporary:
            profiles = Path(temporary)
            (profiles / "griffin-gc").mkdir()
            with (profiles / "griffin-gc/sample.uncorrected.coverage.tsv").open("w") as destination:
                writer = csv.writer(destination, delimiter="\t")
                writer.writerow([
                    "site_name", "number_of_sites", "smoothing", "CNA_normalization",
                    "mean_reads_per_bp_in_normalization_window",
                ] + list(map(str, range(-990, 990, 15))))
                writer.writerow(["TFBS", 2, True, False, 0.5] + [-0.25] * 66 + [2.25] * 66)
            output = StringIO()
            with (
                mock.patch.object(comparator, "read_midpoints", return_value=native_result),
                mock.patch.object(sys, "argv", ["compare_midpoints.py", str(profiles)]),
                redirect_stdout(output),
            ):
                comparator.main()
            row, = csv.DictReader(StringIO(output.getvalue()), delimiter="\t")
            self.assertEqual(row["cfdnalab_sites"], "1")
            self.assertEqual(row["griffin_sites"], "2")
            self.assertEqual(row["positions"], "132")
            self.assertEqual(row["window_bp"], "1980")
            self.assertEqual(row["resolution_bp"], "15")
            self.assertEqual(row["status"], "ok")
            self.assertAlmostEqual(float(row["pearson_r"]), 1)
            self.assertAlmostEqual(float(row["pearson_r_squared"]), 1)
            for metric in ("mae_mean_one", "rmse_mean_one", "max_abs_mean_one"):
                self.assertAlmostEqual(float(row[metric]), 0.25)

    def test_valid_picard_histogram_without_fr_is_an_empty_fr_profile(self):
        import numpy as np

        comparator = self.comparator("lengths/compare_lengths.py")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "picard"
            directory.mkdir()
            (directory / "concordance.metrics.txt").write_text(
                "## HISTOGRAM\tjava.lang.Integer\n"
                "insert_size\tAll_Reads.rf_count\n30\t5\n1000\t4\n\n"
            )
            counts = comparator.read_profile(directory, np.arange(30, 1001))
            self.assertEqual(counts.sum(), 0)

    def test_picard_concordance_reads_filtered_histogram_and_cfdnapro_ignores_prop(self):
        import numpy as np

        comparator = self.comparator("lengths/compare_lengths.py")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            picard = root / "picard"
            picard.mkdir()
            header = "## HISTOGRAM\tjava.lang.Integer\ninsert_size\tAll_Reads.fr_count\n"
            (picard / "metrics.txt").write_text(header + "30\t99\n31\t88\n32\t77\n")
            (picard / "concordance.metrics.txt").write_text(header + "30\t0\n31\t2\n32\t3\n")
            cfdnapro = root / "cfdnapro"
            cfdnapro.mkdir()
            (cfdnapro / "lengths.tsv").write_text(
                "insert_size\tAll_Reads.fr_count\tprop\n30\t0\t0\n31\t2\t0.4\n32\t3\t0.6\n"
            )
            # Concordance uses the filtered counts 0,2,3, not Picard's timed
            # 99,88,77. cfDNAPro's optional proportion column is not another count
            for directory in (picard, cfdnapro):
                np.testing.assert_array_equal(
                    comparator.read_profile(directory, np.array([30, 31, 32])), [0, 2, 3]
                )

    def test_griffin_rounded_density_zero_does_not_erase_a_positive_profile(self):
        import numpy as np

        comparator = self.comparator("midpoints/compare_midpoints.py")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "griffin.tsv"
            columns = [
                "site_name",
                "number_of_sites",
                "smoothing",
                "CNA_normalization",
                "mean_reads_per_bp_in_normalization_window",
            ] + list(map(str, range(-990, 990, 15)))
            with path.open("w") as destination:
                writer = csv.writer(destination, delimiter="\t")
                writer.writerow(columns)
                # The saved subset can have a mean different from one
                writer.writerow(["sparse", 2, False, False, 0] + [4] * 5 + [0] * 127)
                writer.writerow(["empty", 2, False, False, 0] + [""] * 132)
            profiles = comparator.read_griffin(path, 9990, 15)
            self.assertEqual(profiles["sparse"][1:], (2, 0))
            np.testing.assert_array_equal(profiles["sparse"][0], [4] * 5 + [0] * 127)
            self.assertTrue(np.isnan(profiles["empty"][0]).all())
            self.assertEqual(profiles["empty"][2], 0)

    def test_griffin_blank_density_is_undefined_rather_than_zero(self):
        import numpy as np

        comparator = self.comparator("midpoints/compare_midpoints.py")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "griffin.tsv"
            columns = [
                "site_name", "number_of_sites", "smoothing", "CNA_normalization",
                "mean_reads_per_bp_in_normalization_window",
            ] + list(map(str, range(-990, 990, 15)))
            with path.open("w") as destination:
                writer = csv.writer(destination, delimiter="\t")
                writer.writerow(columns)
                writer.writerow(["masked", 2, True, False, ""] + [""] * 132)
            profiles = comparator.read_griffin(path, 9990, 15)
            self.assertTrue(np.isnan(profiles["masked"][0]).all())
            self.assertTrue(np.isnan(profiles["masked"][2]))

            # Undefined density is only valid when the profile is undefined too
            with path.open("w") as destination:
                writer = csv.writer(destination, delimiter="\t")
                writer.writerow(columns)
                writer.writerow(["invalid", 2, True, False, ""] + [1] * 132)
            with self.assertRaisesRegex(ValueError, "Invalid Griffin mean-count metadata"):
                comparator.read_griffin(path, 9990, 15)

    def test_midpoint_normalizes_all_666_bins_before_comparing_132_saved_bins(self):
        import numpy as np
        import pandas as pd

        comparator = self.comparator("midpoints/compare_midpoints.py")
        # The 666 counts sum to 1332, giving full-window mean 2. The central
        # 132 counts are all 4 and must become 2, not be renormalized to 1
        # Dividing the full mean by two sites gives density 1 per site per bp
        counts = np.array([1.0] * 264 + [2.0] * 3 + [4.0] * 132 + [2.0] * 267)
        native_result = SimpleNamespace(
            length_bins=lambda: [(100, 201)],
            positions=lambda: pd.DataFrame({
                "position_bin_start_bp": range(0, 9990, 15),
                "position_bin_end_bp": range(15, 9991, 15),
            }),
            counts_array=lambda: counts[None, None, :],
            group_metadata=lambda: pd.DataFrame({"group_name": ["TFBS"], "eligible_intervals": [2]}),
        )
        with tempfile.TemporaryDirectory() as temporary:
            profiles = Path(temporary)
            (profiles / "griffin-gc").mkdir()
            with (profiles / "griffin-gc/sample.uncorrected.coverage.tsv").open("w") as destination:
                writer = csv.writer(destination, delimiter="\t")
                writer.writerow([
                    "site_name", "number_of_sites", "smoothing", "CNA_normalization",
                    "mean_reads_per_bp_in_normalization_window",
                ] + list(map(str, range(-990, 990, 15))))
                writer.writerow(["TFBS", 2, False, False, 1] + list(counts[267:399] / 2))
            output = StringIO()
            with (
                mock.patch.object(comparator, "read_midpoints", return_value=native_result),
                mock.patch.object(sys, "argv", ["compare_midpoints.py", str(profiles)]),
                redirect_stdout(output),
            ):
                comparator.main()
            row, = csv.DictReader(StringIO(output.getvalue()), delimiter="\t")
            self.assertEqual(row["positions"], "132")
            self.assertEqual(row["window_bp"], "1980")
            self.assertEqual(row["resolution_bp"], "15")
            self.assertAlmostEqual(float(row["cfdnalab_mean_midpoints_per_site_bp"]), 1.0)
            self.assertEqual(row["pearson_r"], "NA")
            self.assertEqual(row["normalization_window_bp"], "9990")
            for metric in ("mae_mean_one", "rmse_mean_one", "max_abs_mean_one"):
                self.assertEqual(float(row[metric]), 0.0)

    def test_midpoint_comparison_rejects_misaligned_native_bins(self):
        comparator = self.comparator("midpoints/compare_midpoints.py")
        # Griffin rounds both windows to [-4995, 4995), whereas cfDNAlab
        # starts its first bin at -5000. Equal step sizes do not align the bins
        for width in (10000, 10001):
            with self.subTest(width=width):
                with self.assertRaisesRegex(ValueError, "does not align"):
                    comparator.read_griffin(Path("unused.tsv"), width, 15)

    def test_midpoint_metrics_keep_zero_positions_and_distinguish_empty_groups(self):
        import numpy as np
        import pandas as pd

        comparator = self.comparator("midpoints/compare_midpoints.py")
        groups = [
            "reversed",
            "both-empty",
            "cfdnalab-empty",
            "griffin-empty",
            "constant",
        ]
        counts = np.array([[2] * 333 + [0] * 333, [0] * 666, [0] * 666, [2] * 666, [2] * 666])
        # Fake only the external reader result, exercising the real comparator CLI
        native_result = SimpleNamespace(
            length_bins=lambda: [(100, 201)],
            positions=lambda: pd.DataFrame(
                {
                    "position_bin_start_bp": range(0, 9990, 15),
                    "position_bin_end_bp": range(15, 9991, 15),
                }
            ),
            counts_array=lambda: counts[:, None, :],
            group_metadata=lambda: pd.DataFrame(
                {"group_name": groups, "eligible_intervals": [2] * 5}
            ),
        )
        with tempfile.TemporaryDirectory() as temporary:
            profiles = Path(temporary) / "profiles"
            (profiles / "griffin-gc").mkdir(parents=True)
            columns = [
                "site_name",
                "number_of_sites",
                "smoothing",
                "CNA_normalization",
                "mean_reads_per_bp_in_normalization_window",
            ] + list(map(str, range(-990, 990, 15)))
            with (profiles / "griffin-gc/sample.uncorrected.coverage.tsv").open(
                "w"
            ) as destination:
                writer = csv.writer(destination, delimiter="\t")
                writer.writerow(columns)
                writer.writerow(
                    ["reversed", 2, False, False, 0.5] + [0] * 66 + [2] * 66
                )
                writer.writerow(["both-empty", 2, False, False, 0] + [""] * 132)
                writer.writerow(["cfdnalab-empty", 2, False, False, 1] + [1] * 132)
                writer.writerow(["griffin-empty", 2, False, False, 0] + [""] * 132)
                writer.writerow(["constant", 2, False, False, 1] + [1] * 132)
            output = StringIO()
            with (
                mock.patch.object(
                    comparator, "read_midpoints", return_value=native_result
                ),
                mock.patch.object(sys, "argv", ["compare_midpoints.py", str(profiles)]),
                redirect_stdout(output),
            ):
                comparator.main()
            rows = {
                row["group"]: row
                for row in csv.DictReader(
                    output.getvalue().splitlines(), delimiter="\t"
                )
            }
            # Opposite half-window peaks: 132 absolute differences of two, r=-1
            self.assertAlmostEqual(float(rows["reversed"]["pearson_r"]), -1)
            self.assertEqual(float(rows["reversed"]["mae_mean_one"]), 2)
            self.assertEqual(float(rows["reversed"]["rmse_mean_one"]), 2)
            self.assertEqual(
                float(rows["reversed"]["cfdnalab_mean_midpoints_per_site_bp"]), 0.5
            )
            for group, status in (
                ("both-empty", "both-zero-coverage"),
                ("cfdnalab-empty", "cfdnalab-zero-coverage"),
                ("griffin-empty", "griffin-zero-coverage"),
            ):
                self.assertEqual(rows[group]["status"], status)
                self.assertEqual(rows[group]["mae_mean_one"], "NA")
                self.assertEqual(rows[group]["pearson_r"], "NA")
            self.assertEqual(rows["constant"]["pearson_r"], "NA")
            self.assertEqual(float(rows["constant"]["mae_mean_one"]), 0)


if __name__ == "__main__":
    unittest.main()
