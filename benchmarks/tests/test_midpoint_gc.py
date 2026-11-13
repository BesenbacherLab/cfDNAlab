"""Check the corrected midpoint workflows using small stand-ins for native tools.

No Slurm jobs, installations or real BAM processing are performed. Each fake
estimator writes a distinct bias value so profiling must use the current result.
"""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest


BENCHMARKS = Path(__file__).resolve().parents[1]


class CorrectedMidpointWorkflowTests(unittest.TestCase):
    """Exercise the native command sequence and its failure boundary."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        (self.work / "bin").mkdir()
        (self.work / "griffin-env/bin").mkdir(parents=True)
        shutil.copyfile(
            BENCHMARKS / "midpoints/griffin_profile.sh",
            self.work / "griffin_profile.sh",
        )
        self.environment = {
            **os.environ,
            "BENCH_WORK": str(self.work),
            "BENCH_MIN_LENGTH": "100",
            "BENCH_MAX_LENGTH": "200",
            "BENCH_MIN_MAPQ": "30",
            "BENCH_HALF_WIDTH": "1000",
            "BENCH_WINDOW_WIDTH": "2000",
            "BENCH_CHROMOSOMES": "chr1,chr2",
            "BENCH_GC_REFERENCE_NAME": "autosomal-regions",
            "FAKE_FAIL_STAGE": "",
        }
        # These stand-ins represent the external command boundary. A later step
        # can only succeed if it receives the estimate written by its predecessor
        native = f"#!{sys.executable}\n" + textwrap.dedent(r'''
            import json
            import os
            from pathlib import Path
            import sys

            work = Path(os.environ['BENCH_WORK'])
            arguments = sys.argv[1:]
            stage = arguments[0]
            if stage.endswith('.py'):
                stage = Path(stage).name
            log = work / 'calls.jsonl'
            with log.open('a') as destination:
                destination.write(json.dumps([stage, arguments]) + '\n')
            if os.environ['FAKE_FAIL_STAGE'] == stage:
                sys.exit(17)
            def argument(name):
                return arguments[arguments.index(name) + 1]
            bias = work / 'output/gc-bias'
            if stage == 'gc-bias':
                estimate = bias / 'gc_bias_correction.zarr'
                estimate.mkdir(exist_ok=True)
                generation = len(log.read_text().splitlines())
                (estimate / 'value').write_text(str(generation))
            elif stage == 'midpoints':
                value = (Path(argument('--gc-file')) / 'value').read_text()
                profile = work / 'output/midpoint_profiles.zarr'
                profile.mkdir(exist_ok=True)
                (profile / 'value').write_text(value)
            elif stage == 'griffin_GC_counts.py':
                (bias / 'GC_counts').mkdir(exist_ok=True)
                generation = len(log.read_text().splitlines())
                (bias / 'GC_counts/sample.GC_counts.txt').write_text(str(generation))
            elif stage == 'griffin_GC_bias.py':
                value = (bias / 'GC_counts/sample.GC_counts.txt').read_text()
                (bias / 'GC_bias').mkdir(exist_ok=True)
                (bias / 'GC_bias/sample.GC_bias.txt').write_text(value)
            elif stage == 'griffin_coverage.py':
                value = Path(argument('--GC_bias')).read_text()
                (work / 'used-bias').write_text(value)
            elif stage == 'griffin_merge_sites.py':
                profile = work / 'output/results/sample'
                profile.mkdir(exist_ok=True)
                (profile / 'sample.GC_corrected.coverage.tsv').write_text(
                    (work / 'used-bias').read_text()
                )
            else:
                raise ValueError('Unexpected native call: ' + stage)
        ''')
        for path in (self.work / "bin/cfdna", self.work / "griffin-env/bin/python"):
            path.write_text(native)
            path.chmod(0o755)

    def run_workflow(self, method):
        """Run only the workflow helper, with all native programs replaced."""
        return subprocess.run(
            ["bash", str(BENCHMARKS / "midpoints/gc_corrected_profile.sh"), method, "4"],
            env=self.environment, text=True, capture_output=True,
        )

    def calls(self):
        """Return the fake tools' actual argument lists in invocation order."""
        return [
            json.loads(line)
            for line in (self.work / "calls.jsonl").read_text().splitlines()
        ]

    def test_cfdnalab_reestimates_bias_before_each_corrected_profile(self):
        for expected_generation in (1, 3):
            result = self.run_workflow("cfdnalab-gc")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                (self.work / "output/midpoint_profiles.zarr/value").read_text(),
                str(expected_generation),
            )
        calls = self.calls()
        self.assertEqual([call[0] for call in calls], ["gc-bias", "midpoints"] * 2)
        for stage, arguments in calls:
            self.assertEqual(arguments[arguments.index("--n-threads") + 1], "4")
            self.assertEqual(arguments[arguments.index("--min-mapq") + 1], "30")
            if stage == "midpoints":
                length_index = arguments.index("--length-bins")
                self.assertEqual(arguments[length_index + 1:length_index + 3], ["100", "201"])
                self.assertEqual(arguments[arguments.index("--bin-size") + 1], "15")
                self.assertNotIn("--smoothing", arguments)
            else:
                # Estimation uses native defaults, independently of the TFBS
                # windows and Griffin's estimation BED
                for option in ("--by-bed", "--chromosomes", "--min-fragment-length", "--max-fragment-length"):
                    self.assertNotIn(option, arguments)

    def test_griffin_reestimates_bias_before_each_corrected_profile(self):
        for expected_generation in (1, 5):
            result = self.run_workflow("griffin-gc")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                (self.work / "output/results/sample/sample.GC_corrected.coverage.tsv").read_text(),
                str(expected_generation),
            )
        calls = self.calls()
        self.assertEqual(
            [call[0] for call in calls],
            ["griffin_GC_counts.py", "griffin_GC_bias.py", "griffin_coverage.py", "griffin_merge_sites.py"] * 2,
        )
        arguments = calls[0][1]
        self.assertEqual(arguments[arguments.index("--CPU") + 1], "4")
        self.assertEqual(arguments[arguments.index("--map_q") + 1], "30")
        # Both estimation stages use Griffin's release range. Coverage keeps the
        # narrower midpoint range, so estimation smoothing has neighbouring lengths
        for stage, arguments in calls:
            if stage in ("griffin_GC_counts.py", "griffin_GC_bias.py"):
                expected_range = ["15", "500"]
            elif stage == "griffin_coverage.py":
                expected_range = ["100", "200"]
            else:
                continue
            length_index = arguments.index("--size_range")
            self.assertEqual(arguments[length_index + 1:length_index + 3], expected_range)

    def test_failed_sample_estimation_stops_before_profiling(self):
        cases = (
            ("cfdnalab-gc", "gc-bias", ["gc-bias"]),
            ("griffin-gc", "griffin_GC_counts.py", ["griffin_GC_counts.py"]),
            ("griffin-gc", "griffin_GC_bias.py", ["griffin_GC_counts.py", "griffin_GC_bias.py"]),
        )
        for method, failure, expected_stages in cases:
            with self.subTest(method=method, failure=failure):
                (self.work / "calls.jsonl").unlink(missing_ok=True)
                self.environment["FAKE_FAIL_STAGE"] = failure
                result = self.run_workflow(method)
                self.assertEqual(result.returncode, 17, result.stderr)
                self.assertEqual([call[0] for call in self.calls()], expected_stages)
                self.assertFalse((self.work / "output/midpoint_profiles.zarr").exists())
                self.assertFalse((self.work / "used-bias").exists())

    def test_griffin_retains_wide_windows_and_native_profiling_steps(self):
        for width, half_width, end in ((9990, "4995", "4995"), (10001, "5000", "5001")):
            with self.subTest(width=width):
                (self.work / "calls.jsonl").unlink(missing_ok=True)
                self.environment["BENCH_WINDOW_WIDTH"] = str(width)
                self.environment["BENCH_HALF_WIDTH"] = half_width
                result = self.run_workflow("griffin-gc")
                self.assertEqual(result.returncode, 0, result.stderr)
                for stage, arguments in self.calls():
                    if stage in ("griffin_coverage.py", "griffin_merge_sites.py"):
                        index = arguments.index("--norm_window")
                        self.assertEqual(arguments[index + 1:index + 3], ["-" + half_width, end])
                        self.assertEqual(arguments[arguments.index("--mappability_bw") + 1], str(self.work / "mappability.bw"))
                        self.assertEqual(arguments[arguments.index("--mappability_correction") + 1], "False")
                    if stage == "griffin_merge_sites.py":
                        index = arguments.index("--save_window")
                        self.assertEqual(arguments[index + 1:index + 3], ["-1000", "1000"])
                        self.assertEqual(arguments[arguments.index("--step") + 1], "15")
                        for option in ("--exclude_outliers", "--exclude_zero_mappability", "--smoothing"):
                            self.assertEqual(arguments[arguments.index(option) + 1], "True")
                        self.assertEqual(arguments[arguments.index("--CNA_normalization") + 1], "False")
                        self.assertEqual(arguments[arguments.index("--fft_index") + 1], "10")
                        index = arguments.index("--fft_window")
                        self.assertEqual(arguments[index + 1:index + 3], ["-960", "960"])
                        index = arguments.index("--center_window")
                        self.assertEqual(arguments[index + 1:index + 3], ["-30", "30"])
                        index = arguments.index("--exclude_paths")
                        self.assertEqual(len(arguments[index + 1:arguments.index("--exclude_zero_mappability")]), 5)

    def test_griffin_rejects_profiles_outside_estimation_range_before_setup(self):
        for minimum, maximum in ((14, 200), (100, 501)):
            with self.subTest(minimum=minimum, maximum=maximum):
                results = self.work / "results"
                result = subprocess.run(
                    ["bash", str(BENCHMARKS / "midpoints/benchmark_midpoints.sh")],
                    env={
                        **self.environment,
                        "SLURM_SUBMIT_DIR": str(BENCHMARKS.parent),
                        "SLURM_JOB_ID": "123",
                        "SLURM_JOB_NUM_NODES": "1",
                        "SLURM_CPUS_PER_TASK": "4",
                        "BENCH_BAM": "/unused/input.bam",
                        "BENCH_INTERVALS": "/unused/sites.bed",
                        "BENCH_RESULTS": str(results),
                        "BENCH_CORES": "4",
                        "BENCH_METHODS": "griffin-gc",
                        "BENCH_PILOT": "0",
                        "BENCH_SETUP_PYTHON": sys.executable,
                        "BENCH_EXPECTED_NODE": "",
                        "BENCH_RUN_TIMEOUT_SECONDS": "0",
                        "BENCH_MIN_LENGTH": str(minimum),
                        "BENCH_MAX_LENGTH": str(maximum),
                        "TMPDIR": str(self.work),
                    },
                    text=True, capture_output=True,
                )
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn("Griffin GC estimation uses 15-500 bp", result.stderr)
                self.assertFalse(results.exists())


if __name__ == "__main__":
    unittest.main()
