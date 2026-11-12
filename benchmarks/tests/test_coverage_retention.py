"""Check that saved coverage profiles retain every processed chromosome.

These cases use synthetic depth tables and do not run Slurm or native tools.
Linux is required because the result-saving helper uses GNU sync -f.
"""

import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


BENCHMARKS = Path(__file__).resolve().parents[1]


@unittest.skipUnless(sys.platform.startswith("linux"), "Requires GNU sync -f")
class CoverageRetentionTests(unittest.TestCase):
    """Keep whole-genome scope when saving positional depth outputs."""

    def test_first_measured_depth_profile_keeps_all_chromosomes(self):
        for method in ("samtools-depth", "pandepth"):
            with (
                self.subTest(method=method),
                tempfile.TemporaryDirectory() as temporary,
            ):
                root = Path(temporary)
                work, results = root / "work", root / "results"
                coverage = work / "coverage"
                coverage.mkdir(parents=True)
                results.mkdir()
                # The old retention rule dropped chr1 despite processing it
                chromosome_sizes = "chr1\t2\nchr21\t1\nchr22\t1\n"
                (work / "chrom.sizes").write_text(chromosome_sizes)
                (work / "capture.chrom.sizes").write_text("chr21\t1\nchr22\t1\n")
                if method == "samtools-depth":
                    depths = "chr1\t1\t2\nchr1\t2\t0\nchr21\t1\t3\nchr22\t1\t4\n"
                    profile = coverage / "samtools.depth.tsv"
                    profile.write_text(depths)
                else:
                    depths = "chr1\t0\t2\nchr1\t1\t0\nchr21\t0\t3\nchr22\t0\t4\n"
                    profile = coverage / "pandepth.SiteDepth.gz"
                    native_compressed = gzip.compress(depths.encode())
                    profile.write_bytes(native_compressed)
                (work / "profile.path").write_text(str(profile) + "\n")
                (work / "time.tsv").write_text("fixture\t1\t1\t0\t10\t0\n")
                result = subprocess.run(
                    [
                        "bash", "-euo", "pipefail", "-c",
                        'source "$BENCH_COMMON"\n'
                        'run_order=2\n'
                        'benchmark_finish_run "$SELECTED_METHOD" 2 false 0\n',
                    ],
                    env={
                        **os.environ,
                        "BENCH_COMMON": str(BENCHMARKS / "common.sh"),
                        "BENCH_WORK": str(work),
                        "BENCH_RESULTS": str(results),
                        "BENCH_PILOT": "0",
                        "BENCH_WARMUPS": "1",
                        "SELECTED_METHOD": method,
                    },
                    text=True,
                    capture_output=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                retained = results / "profiles" / method
                with gzip.open(retained / "positions.tsv.gz", "rt") as source:
                    self.assertEqual(source.read(), depths)
                self.assertEqual(
                    (retained / "chrom.sizes").read_text(), chromosome_sizes
                )
                if method == "pandepth":
                    # Already-compressed native output should be copied intact
                    self.assertEqual(
                        (retained / "positions.tsv.gz").read_bytes(), native_compressed
                    )


if __name__ == "__main__":
    unittest.main()
