"""Tests for the benchmark job scripts, their native commands and the measurement step.

These tests call the real scripts and shell functions. Native tools, srun and the
comparison Python are replaced by small fakes at the command boundary, so no tool
installation, BAM or Slurm allocation is needed. Expected values are derived by
hand from the fixtures.

Requirements:
- Bash 4 or later as `bash` on PATH. The scripts use associative arrays.
- The measurement and whole-job tests also need Linux: GNU time at /usr/bin/time,
  GNU coreutils (sync -f, du -sb) and the CPU topology under /sys. They are skipped
  on other systems. The fake srun binds each step to a single CPU, so the real
  placement check in measure_command.sh passes with BENCH_CORES=1.
"""

import csv
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import textwrap
import unittest


BENCHMARKS = Path(__file__).resolve().parents[1]
COMMON = BENCHMARKS / "common.sh"
MEASURE_COMMAND = BENCHMARKS / "measure_command.sh"
EXCLUSION_FILENAMES = (
    "encode_unified_GRCh38_exclusion_list.bed", "hg38_centromeres.bed", "hg38_gaps.bed",
    "hg38_fix_patches.bed", "hg38_alternative_haplotypes.bed",
)

# A scratch path with a space, so the tests also check that arguments stay whole
WORK = "/scratch/job 7/work"

# Columns that common.sh fills for every run. prepare_command fills the others
IDENTITY_AND_MEASUREMENT_COLUMNS = {
    "run_order", "repetition", "warmup", "method", "cores", "job_id", "pilot",
    "step_exit_code", "step_id", "node", "cpu_model", "cpu_ids", "physical_core_ids",
    "started_utc", "wall_seconds", "user_seconds", "system_seconds", "max_rss_kib",
    "tool_exit_code", "flush_seconds", "process_plus_flush_seconds", "output_bytes",
}


def single_cpu():
    """Return a CPU this process may use, or None where affinity is unavailable."""
    if not hasattr(os, "sched_getaffinity"):
        return None
    return min(os.sched_getaffinity(0))


def measurement_supported():
    """Whether this system has what measure_command.sh needs."""
    cpu = single_cpu()
    return (
        sys.platform.startswith("linux")
        and Path("/usr/bin/time").exists()
        and cpu is not None
        and Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/core_id").exists()
    )


requires_measurement = unittest.skipUnless(
    measurement_supported(),
    "Requires Linux with GNU time, GNU coreutils and CPU topology in /sys",
)


def clean_environment(**variables):
    """Return the test environment without benchmark or Slurm settings from the caller."""
    environment = {
        name: value for name, value in os.environ.items()
        if not name.startswith(("BENCH_", "SLURM_"))
    }
    environment.update(variables)
    return environment


def write_executable(path, text):
    """Write a script and make it executable."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(0o755)


def hostname():
    """Return the name measure_command.sh compares with BENCH_NODE."""
    return subprocess.run(["hostname"], text=True, capture_output=True, check=True).stdout.strip()


def read_record(path):
    """Read measure_command.sh's "name<TAB>value" lines in order."""
    return [tuple(line.split("\t", 1)) for line in path.read_text().splitlines()]


def read_timings(path):
    """Read timings.tsv as a list of dictionaries."""
    with path.open() as source:
        return list(csv.DictReader(source, delimiter="\t"))


# ---------------------------------------------------------------------------
# Native commands in <feature>/commands.sh
# ---------------------------------------------------------------------------


def prepare(feature, method, variables, script_prefix=""):
    """Call prepare_command from FEATURE/commands.sh and return what it set.

    Returns a dictionary with the command arguments, the row settings, the
    expected output, the stdout file, whether an empty output is allowed and the
    feature's timing columns.
    """
    script = (
        'source "$COMMON"\nsource "$COMMANDS"\n' + script_prefix
        + 'prepare_command "$METHOD"\n'
        'printf "command\\t%s\\n" "${command[@]}"\n'
        'for key in "${!row[@]}"; do printf "row\\t%s\\t%s\\n" "$key" "${row[$key]}"; done\n'
        'printf "expected_output\\t%s\\n" "$expected_output"\n'
        'printf "stdout_file\\t%s\\n" "$stdout_file"\n'
        'printf "empty_output_allowed\\t%s\\n" "$empty_output_allowed"\n'
        'printf "timing_column\\t%s\\n" "${timing_columns[@]}"\n'
    )
    result = subprocess.run(
        ["bash", "-c", "set -euo pipefail\n" + script],
        env=clean_environment(
            COMMON=str(COMMON),
            COMMANDS=str(BENCHMARKS / feature / "commands.sh"),
            METHOD=method,
            **variables,
        ),
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        raise AssertionError(f"prepare_command {method} failed: {result.stderr}")
    prepared = {"command": [], "row": {}, "timing_columns": []}
    for line in result.stdout.splitlines():
        kind, value = line.split("\t", 1)
        if kind == "command":
            prepared["command"].append(value)
        elif kind == "row":
            key, setting = value.split("\t", 1)
            prepared["row"][key] = setting
        elif kind == "timing_column":
            prepared["timing_columns"].append(value)
        else:
            prepared[kind] = value
    return prepared


def argument_after(command, option):
    """Return the argument that follows an option."""
    return command[command.index(option) + 1]


LENGTH_VARIABLES = {
    "BENCH_WORK": WORK,
    "cores": "4",
    "bam": f"{WORK}/sample.bam",
    "output": f"{WORK}/output",
    "min_mapq": "30",
    "min_length": "30",
    "max_length": "1000",
    "length_end": "1001",
    "autosomes": "chr1,chr2",
    "genome": "hg38",
    "cfdna": f"{WORK}/bin/cfdna",
    "samtools": f"{WORK}/samtools-env/bin/samtools",
    "picard": f"{WORK}/picard-env/bin/picard",
    "finale_python": f"{WORK}/finale-env/bin/python",
    "finale": f"{WORK}/finale-env/bin/finaletoolkit",
    "apptainer": "/usr/bin/apptainer",
    "cfdnapro_image": f"{WORK}/cfdnapro.sif",
    "cfdna_version": "cfdna fixture",
    "samtools_version": "samtools fixture",
    "picard_version": "picard fixture",
    "finale_version": "finale fixture",
    "cfdnapro_version": "cfdnapro fixture",
}
LENGTH_METHODS = (
    "cfdnalab", "cfdnalab-autosomes", "cfdnalab-proper", "cfdnapro",
    "cfdnapro-outer-span", "finaletoolkit", "samtools", "picard",
)

COVERAGE_VARIABLES = {
    **LENGTH_VARIABLES,
    "regions": f"{WORK}/regions.bed",
    "mosdepth": f"{WORK}/bin/mosdepth",
    "pandepth": f"{WORK}/bin/pandepth",
    "mosdepth_version": "mosdepth fixture",
    "pandepth_version": "pandepth fixture",
}
COVERAGE_QUANTITIES = {
    "cfdnalab": "per-base-read-depth",
    "cfdnalab-autosomes": "per-base-read-depth",
    "mosdepth": "per-base-read-depth",
    "samtools-depth": "per-base-read-depth",
    "pandepth": "per-base-read-depth",
    "cfdnalab-fragment": "per-base-fragment-coverage",
    "cfdnalab-fragment-proper": "per-base-fragment-coverage",
    "mosdepth-fragment": "per-base-fragment-coverage",
    "cfdnalab-regions": "regional-read-depth-500kb",
    "mosdepth-regions": "regional-read-depth-500kb",
    "samtools-regions": "regional-read-depth-500kb",
    "pandepth-regions": "regional-read-depth-500kb",
    "cfdnalab-region-counts": "region-counts-500kb",
    "cfdnalab-region-counts-proper": "region-counts-500kb",
    "finaletoolkit-coverage": "region-counts-500kb",
}

END_MOTIF_VARIABLES = {
    **LENGTH_VARIABLES,
    "k": "4",
    "reference_2bit": f"{WORK}/reference.2bit",
    "freia_version": "freia fixture",
}
END_MOTIF_METHODS = (
    "cfdnalab-reference", "cfdnalab-reference-proper", "cfdnalab-reference-proper-wide",
    "cfdnalab-reference-autosomes", "cfdnalab-read-aligned", "cfdnalab-read",
    "cfdnalab-read-include", "finaletoolkit", "cfdnapro", "freia",
)

MIDPOINT_VARIABLES = {
    "BENCH_WORK": WORK,
    "cores": "4",
    "bam": f"{WORK}/sample.bam",
    "output": f"{WORK}/output",
    "min_mapq": "20",
    "min_length": "100",
    "max_length": "200",
    "length_end": "201",
    "cfdna": f"{WORK}/bin/cfdna",
    "cfdna_version": "cfdna fixture",
    "griffin_version": "griffin fixture",
    "BENCH_CHROMOSOMES": "chr1,chr2",
    "BENCH_SITES": "4",
    "BENCH_GROUPS": "2",
    "BENCH_WINDOW_WIDTH": "9990",
}
MIDPOINT_EXCLUSIONS = [f"{WORK}/exclusions/{filename}" for filename in EXCLUSION_FILENAMES]
# The exclusions array cannot pass through the environment, so it is set in the script
MIDPOINT_PREFIX = "exclusions=(" + " ".join(f'"{path}"' for path in MIDPOINT_EXCLUSIONS) + ")\n"


class CommandRowTests(unittest.TestCase):
    """Every method must fill exactly the settings columns of its feature's timings.tsv."""

    def test_every_method_fills_exactly_the_settings_columns(self):
        for feature, methods, variables, prefix in (
            ("lengths", LENGTH_METHODS, LENGTH_VARIABLES, ""),
            ("coverage", tuple(COVERAGE_QUANTITIES), COVERAGE_VARIABLES, ""),
            ("ends", END_MOTIF_METHODS, END_MOTIF_VARIABLES, ""),
            ("midpoints", ("cfdnalab", "cfdnalab-gc", "griffin-gc"), MIDPOINT_VARIABLES, MIDPOINT_PREFIX),
        ):
            for method in methods:
                with self.subTest(feature=feature, method=method):
                    prepared = prepare(feature, method, variables, prefix)
                    columns = prepared["timing_columns"]
                    settings = set(prepared["row"])

                    # write_timing_row needs a value for every column and no others.
                    # common.sh adds the identity and measurement columns itself
                    self.assertEqual(len(columns), len(set(columns)))
                    self.assertEqual(settings & IDENTITY_AND_MEASUREMENT_COLUMNS, set())
                    self.assertEqual(settings | IDENTITY_AND_MEASUREMENT_COLUMNS, set(columns))
                    self.assertTrue(prepared["command"])
                    self.assertTrue(prepared["expected_output"].startswith(f"{WORK}/output/"))

    def test_unknown_methods_are_rejected(self):
        for feature, variables in (
            ("lengths", LENGTH_VARIABLES), ("coverage", COVERAGE_VARIABLES),
            ("ends", END_MOTIF_VARIABLES), ("midpoints", MIDPOINT_VARIABLES),
        ):
            with self.subTest(feature=feature):
                with self.assertRaisesRegex(AssertionError, "Unknown"):
                    prepare(feature, "not-a-method", variables)


class LengthCommandTests(unittest.TestCase):
    """The literal fragment length commands and their recorded settings."""

    def test_cfdnalab_counts_every_length_in_1_bp_bins_on_all_contigs(self):
        prepared = prepare("lengths", "cfdnalab", LENGTH_VARIABLES)

        # --length-bins START:END:STEP has an exclusive end, so 1001 includes 1000 bp
        self.assertEqual(prepared["command"], [
            f"{WORK}/bin/cfdna", "lengths", "--bam", f"{WORK}/sample.bam",
            "--output-dir", f"{WORK}/output", "--length-bins", "30:1001:1",
            "--indel-mode", "ignore", "--clip-mode", "aligned", "--min-mapq", "30",
            "--n-threads", "4", "--chromosomes", "all",
        ])
        self.assertEqual(prepared["expected_output"], f"{WORK}/output/length_counts.tsv.zst")
        self.assertEqual(prepared["row"]["tool_threads"], "4")
        self.assertEqual(prepared["row"]["length_definition"], "directional-aligned-span")
        self.assertEqual(prepared["row"]["scope"], "all")
        self.assertEqual(prepared["empty_output_allowed"], "false")

    def test_cfdnalab_variants_change_only_scope_or_pairing(self):
        base = prepare("lengths", "cfdnalab", LENGTH_VARIABLES)["command"]
        autosomes = prepare("lengths", "cfdnalab-autosomes", LENGTH_VARIABLES)
        proper = prepare("lengths", "cfdnalab-proper", LENGTH_VARIABLES)

        # The autosome call replaces "all" with the BAM's autosome names
        self.assertEqual(autosomes["command"], base[:-1] + ["chr1,chr2"])
        self.assertEqual(autosomes["row"]["scope"], "autosomes")
        # The proper-pair call appends a single option
        self.assertEqual(proper["command"], base + ["--require-proper-pair"])
        self.assertEqual(proper["row"]["proper_pair"], "required")

    def test_samtools_prints_its_report_into_the_output_directory(self):
        prepared = prepare("lengths", "samtools", LENGTH_VARIABLES)

        # -i 1001 puts the overflow bin above the requested maximum of 1000 bp
        self.assertEqual(prepared["command"], [
            f"{WORK}/samtools-env/bin/samtools", "stats", "-@", "4",
            "--input-fmt-option", "filter=mapq >= 30",
            "-f", "1", "-F", "3852", "-m", "1", "-i", "1001", f"{WORK}/sample.bam",
        ])
        self.assertEqual(prepared["stdout_file"], f"{WORK}/output/stats.txt")
        self.assertEqual(prepared["expected_output"], f"{WORK}/output/stats.txt")
        self.assertEqual(prepared["row"]["tool_threads"], "4")

    def test_timed_picard_reads_the_bam_directly_without_a_mapq_filter(self):
        prepared = prepare("lengths", "picard", LENGTH_VARIABLES)

        self.assertEqual(prepared["command"], [
            f"{WORK}/picard-env/bin/picard", "-Xmx2g", "CollectInsertSizeMetrics",
            f"I={WORK}/sample.bam", f"O={WORK}/output/metrics.txt", "HISTOGRAM_WIDTH=1000",
            f"TMP_DIR={WORK}/tmp",
        ])
        # Picard has no MAPQ or thread option, which the row must say
        self.assertEqual(prepared["row"]["min_mapq"], "NA")
        self.assertEqual(prepared["row"]["tool_threads"], "NA")
        self.assertEqual(prepared["stdout_file"], f"{WORK}/tool-stdout.log")

    def test_cfdnapro_curated_and_outer_span_runs_differ_only_in_curation(self):
        curated = prepare("lengths", "cfdnapro", LENGTH_VARIABLES)
        outer = prepare("lengths", "cfdnapro-outer-span", LENGTH_VARIABLES)

        # The container sees $BENCH_WORK as /work, so /work/output is the output directory
        self.assertEqual(curated["command"], [
            "/usr/bin/apptainer", "exec", "--cleanenv", "--containall",
            "--bind", f"{WORK}:/work", "--pwd", "/work",
            f"{WORK}/cfdnapro.sif", "Rscript", "--vanilla", "/work/cfdnapro_lengths.R",
            "/work/sample.bam", "/work/output/lengths.tsv", "30", "1000", "30", "hg38", "TRUE",
        ])
        self.assertEqual(outer["command"], curated["command"][:-1] + ["FALSE"])
        self.assertEqual(curated["row"]["length_definition"], "directional-aligned-span")
        self.assertEqual(outer["row"]["length_definition"], "outer-aligned-span")
        self.assertEqual(curated["row"]["tool_threads"], "NA")

    def test_finaletoolkit_runs_its_launcher_with_the_packed_python(self):
        prepared = prepare("lengths", "finaletoolkit", LENGTH_VARIABLES)

        self.assertEqual(prepared["command"], [
            f"{WORK}/finale-env/bin/python", f"{WORK}/finale-env/bin/finaletoolkit",
            "frag-length-bins", f"{WORK}/sample.bam", "--bin-size", "1",
            "--min-length", "30", "--max-length", "1000", "--min-mapq", "30",
            "--output", f"{WORK}/output/lengths.tsv",
        ])
        self.assertEqual(prepared["row"]["proper_pair"], "required")
        self.assertEqual(prepared["row"]["tool_threads"], "NA")


class CoverageCommandTests(unittest.TestCase):
    """The literal coverage commands, thread options and compared quantities."""

    def test_every_method_reports_the_quantity_it_measures(self):
        for method, quantity in COVERAGE_QUANTITIES.items():
            with self.subTest(method=method):
                prepared = prepare("coverage", method, COVERAGE_VARIABLES)
                self.assertEqual(prepared["row"]["workload"], quantity)

    def test_every_thread_option_receives_the_allocated_core_count(self):
        def thread_option(method):
            """The native thread option of the tool behind a method."""
            if method.startswith("cfdnalab"):
                return "--n-threads"
            if method.startswith(("mosdepth", "finaletoolkit")):
                return "--threads"
            if method.startswith("pandepth"):
                return "-t"
            return "-@"

        for cores in ("1", "4", "12", "32"):
            for method in COVERAGE_QUANTITIES:
                with self.subTest(cores=cores, method=method):
                    prepared = prepare("coverage", method, {**COVERAGE_VARIABLES, "cores": cores})
                    command, row = prepared["command"], prepared["row"]
                    if method == "samtools-regions":
                        # bedcov has no thread option
                        self.assertNotIn("-@", command)
                        self.assertEqual(row["tool_threads"], "NA")
                    else:
                        self.assertEqual(argument_after(command, thread_option(method)), cores)
                        self.assertEqual(row["tool_threads"], cores)

    def test_per_base_read_depth_commands(self):
        cfdnalab = prepare("coverage", "cfdnalab", COVERAGE_VARIABLES)
        autosomes = prepare("coverage", "cfdnalab-autosomes", COVERAGE_VARIABLES)
        mosdepth = prepare("coverage", "mosdepth", COVERAGE_VARIABLES)
        samtools = prepare("coverage", "samtools-depth", COVERAGE_VARIABLES)
        pandepth = prepare("coverage", "pandepth", COVERAGE_VARIABLES)

        # Read depth leaves out the gap between mates
        self.assertEqual(cfdnalab["command"], [
            f"{WORK}/bin/cfdna", "fcoverage", "--bam", f"{WORK}/sample.bam",
            "--output-dir", f"{WORK}/output", "--min-mapq", "30",
            "--min-fragment-length", "30", "--max-fragment-length", "1000",
            "--n-threads", "4", "--chromosomes", "all", "--ignore-gap",
        ])
        self.assertEqual(argument_after(autosomes["command"], "--chromosomes"), "chr1,chr2")
        self.assertEqual(autosomes["row"]["scope"], "autosomes")
        self.assertEqual(mosdepth["command"], [
            f"{WORK}/bin/mosdepth", "--threads", "4", "--include-flag", "1", "--flag", "3852",
            "--mapq", "30", "--min-frag-len", "30", "--max-frag-len", "1000",
            f"{WORK}/output/coverage", f"{WORK}/sample.bam",
        ])
        self.assertEqual(samtools["command"], [
            f"{WORK}/samtools-env/bin/samtools", "depth", "-@", "4",
            "--input-fmt-option",
            "filter=(tlen >= 30 && tlen <= 1000) || (tlen <= -30 && tlen >= -1000)",
            "-s", "-q", "0", "-Q", "30", "-G", "3852", "--require-flags", "1",
            "-o", f"{WORK}/output/samtools.depth.tsv", f"{WORK}/sample.bam",
        ])
        # Only samtools depth's sparse output may be empty: it then means zero coverage
        self.assertEqual(samtools["empty_output_allowed"], "true")
        for prepared in (cfdnalab, autosomes, mosdepth, pandepth):
            self.assertEqual(prepared["empty_output_allowed"], "false")
        self.assertEqual(pandepth["command"], [
            f"{WORK}/bin/pandepth", "-i", f"{WORK}/sample.bam", "-o", f"{WORK}/output/pandepth",
            "-q", "30", "-x", "3852", "-t", "4", "-a",
        ])
        self.assertEqual(pandepth["row"]["length_filter"], "unsupported")

    def test_fragment_coverage_includes_the_gap_between_mates(self):
        cfdnalab = prepare("coverage", "cfdnalab-fragment", COVERAGE_VARIABLES)
        proper = prepare("coverage", "cfdnalab-fragment-proper", COVERAGE_VARIABLES)
        mosdepth = prepare("coverage", "mosdepth-fragment", COVERAGE_VARIABLES)

        self.assertNotIn("--ignore-gap", cfdnalab["command"])
        self.assertEqual(proper["command"], cfdnalab["command"] + ["--require-proper-pair"])
        self.assertEqual(proper["row"]["proper_pair"], "required")
        self.assertEqual(mosdepth["command"][-3:], ["--fragment-mode", f"{WORK}/output/coverage", f"{WORK}/sample.bam"])
        self.assertEqual(mosdepth["row"]["proper_pair"], "required")

    def test_regional_commands_write_window_tables_without_per_base_output(self):
        cfdnalab = prepare("coverage", "cfdnalab-regions", COVERAGE_VARIABLES)
        mosdepth = prepare("coverage", "mosdepth-regions", COVERAGE_VARIABLES)
        samtools = prepare("coverage", "samtools-regions", COVERAGE_VARIABLES)
        pandepth = prepare("coverage", "pandepth-regions", COVERAGE_VARIABLES)
        counts = prepare("coverage", "cfdnalab-region-counts", COVERAGE_VARIABLES)
        finale = prepare("coverage", "finaletoolkit-coverage", COVERAGE_VARIABLES)

        self.assertEqual(cfdnalab["command"][-7:], [
            "--chromosomes", "all", "--ignore-gap", "--by-size", "500000", "--per-window", "total",
        ])
        self.assertNotIn("--normalize-by-length", cfdnalab["command"])
        self.assertEqual(cfdnalab["expected_output"], f"{WORK}/output/fcoverage.total.tsv.zst")
        # The precision setting is part of the recorded command
        self.assertEqual(mosdepth["command"][:3], ["env", "MOSDEPTH_PRECISION=5", f"{WORK}/bin/mosdepth"])
        self.assertEqual(mosdepth["command"][-5:], [
            "--by", "500000", "--no-per-base", f"{WORK}/output/coverage", f"{WORK}/sample.bam",
        ])
        self.assertEqual(mosdepth["expected_output"], f"{WORK}/output/coverage.regions.bed.gz")
        # bedcov prints its table, which goes straight into the output directory
        self.assertEqual(samtools["command"], [
            f"{WORK}/samtools-env/bin/samtools", "bedcov", "--input-fmt-option",
            "filter=flag.paired && ((tlen >= 30 && tlen <= 1000) || (tlen <= -30 && tlen >= -1000))",
            "-j", "-Q", "30", "-G", "3852", f"{WORK}/regions.bed", f"{WORK}/sample.bam",
        ])
        self.assertEqual(samtools["stdout_file"], f"{WORK}/output/samtools.bedcov.tsv")
        self.assertEqual(samtools["expected_output"], samtools["stdout_file"])
        self.assertEqual(argument_after(pandepth["command"], "-w"), "500000")
        self.assertNotIn("-a", pandepth["command"])
        self.assertEqual(pandepth["expected_output"], f"{WORK}/output/pandepth.win.stat.gz")
        self.assertIn("--normalize-by-length", counts["command"])
        self.assertEqual(counts["expected_output"], f"{WORK}/output/length_normalized.fcoverage.total.tsv.zst")
        # Leaving these options unset keeps midpoint assignment, unit scale and no
        # genome-wide normalization
        self.assertEqual(finale["command"][2:5], ["coverage", f"{WORK}/sample.bam", f"{WORK}/regions.bed"])
        for option in ("--intersect-policy", "--scale-factor", "--normalize"):
            self.assertNotIn(option, finale["command"])
        self.assertEqual(finale["row"]["proper_pair"], "required")


class EndMotifCommandTests(unittest.TestCase):
    """The literal end-motif commands and their recorded settings."""

    def test_cfdnalab_configurations_choose_motif_source_clipping_and_scope(self):
        # method: (source, clip, proper_pair, scope, --chromosomes value)
        expected = {
            "cfdnalab-reference": ("reference", "aligned", "not-required", "all", "all"),
            "cfdnalab-reference-proper": ("reference", "aligned", "required", "all", "all"),
            "cfdnalab-reference-autosomes": ("reference", "aligned", "not-required", "autosomes", "chr1,chr2"),
            "cfdnalab-read-aligned": ("read", "aligned", "not-required", "all", "all"),
            "cfdnalab-read": ("read", "skip", "not-required", "all", "all"),
            "cfdnalab-read-include": ("read", "include-at-aligned-boundary", "not-required", "all", "all"),
        }
        for method, (source, clip, proper_pair, scope, chromosomes) in expected.items():
            with self.subTest(method=method):
                prepared = prepare("ends", method, END_MOTIF_VARIABLES)
                command, row = prepared["command"], prepared["row"]

                self.assertEqual(command[:12], [
                    f"{WORK}/bin/cfdna", "ends", "--bam", f"{WORK}/sample.bam",
                    "--output-dir", f"{WORK}/output", "--k-inside", "4", "--k-outside", "0",
                    "--all-motifs", "--min-mapq",
                ])
                self.assertEqual(argument_after(command, "--source-inside"), source)
                self.assertEqual(argument_after(command, "--clip-strategy"), clip)
                self.assertEqual(argument_after(command, "--chromosomes"), chromosomes)
                self.assertEqual(argument_after(command, "--n-threads"), "4")
                # Only reference motifs need the reference sequence
                self.assertEqual("--ref-2bit" in command, source == "reference")
                self.assertEqual("--require-proper-pair" in command, proper_pair == "required")
                self.assertEqual(
                    (row["source"], row["clip"], row["proper_pair"], row["scope"]),
                    (source, clip, proper_pair, scope),
                )
                self.assertEqual(prepared["expected_output"], f"{WORK}/output/end_motifs.zarr")

    def test_wide_configuration_records_requested_and_effective_bounds_separately(self):
        prepared = prepare("ends", "cfdnalab-reference-proper-wide", END_MOTIF_VARIABLES)
        command, row = prepared["command"], prepared["row"]

        self.assertEqual(argument_after(command, "--min-fragment-length"), "10")
        self.assertEqual(argument_after(command, "--max-fragment-length"), "50000")
        self.assertEqual((row["requested_min_length"], row["requested_max_length"]), ("30", "1000"))
        self.assertEqual((row["effective_min_length"], row["effective_max_length"]), ("10", "50000"))

    def test_other_tools_use_their_native_calls(self):
        finale = prepare("ends", "finaletoolkit", END_MOTIF_VARIABLES)
        cfdnapro = prepare("ends", "cfdnapro", END_MOTIF_VARIABLES)
        freia = prepare("ends", "freia", END_MOTIF_VARIABLES)

        self.assertEqual(finale["command"], [
            f"{WORK}/finale-env/bin/python", f"{WORK}/finale-env/bin/finaletoolkit", "end-motifs",
            f"{WORK}/sample.bam", f"{WORK}/reference.2bit", "-k", "4",
            "--min-length", "30", "--max-length", "1000", "--min-mapq", "30",
            "--threads", "4", "--strand", "both", "--output", f"{WORK}/output/motifs.tsv",
        ])
        self.assertEqual((finale["row"]["proper_pair"], finale["row"]["clip"]), ("required", "unrestricted"))
        self.assertEqual(cfdnapro["command"][-7:], [
            "/work/sample.bam", "/work/output", "30", "1000", "30", "hg38", "4",
        ])
        self.assertEqual(cfdnapro["row"]["tool_threads"], "NA")
        self.assertEqual(cfdnapro["row"]["length_filter"], "applied-in-r")
        self.assertEqual(cfdnapro["row"]["scope"], "bam-contigs-in-cfdnapro-reference")
        self.assertEqual(freia["command"], ["bash", f"{WORK}/freia_motifs.sh", "4"])
        self.assertEqual(freia["row"]["source"], "read")
        self.assertEqual(
            freia["expected_output"],
            f"{WORK}/output/4_FrEIA/3_Abundances/benchmark/Genome_Lvl/T__sample.pq",
        )


class MidpointCommandTests(unittest.TestCase):
    """The literal midpoint commands and their recorded settings."""

    def test_uncorrected_cfdnalab_uses_15_bp_bins_and_the_shared_exclusions(self):
        prepared = prepare("midpoints", "cfdnalab", MIDPOINT_VARIABLES, MIDPOINT_PREFIX)

        # --length-bins has an exclusive end, so 201 includes 200 bp fragments
        self.assertEqual(prepared["command"], [
            f"{WORK}/bin/cfdna", "midpoints", "--bam", f"{WORK}/sample.bam",
            "--output-dir", f"{WORK}/output", "--intervals", f"{WORK}/sites.bed",
            "--chromosomes", "chr1,chr2", "--length-bins", "100", "201",
            "--min-mapq", "20", "--bin-size", "15", "--n-threads", "4",
            "--blacklist", *MIDPOINT_EXCLUSIONS,
        ])
        row = prepared["row"]
        self.assertEqual(
            (row["sites"], row["groups"], row["window_bp"], row["resolution_bp"], row["gc_correction"]),
            ("4", "2", "9990", "15", "none"),
        )
        self.assertEqual(row["midpoint_policy"], "coordinate-random-even")

    def test_gc_corrected_methods_time_their_whole_workflow(self):
        for method, expected_output, policy in (
            ("cfdnalab-gc", f"{WORK}/output/midpoint_profiles.zarr", "coordinate-random-even"),
            ("griffin-gc", f"{WORK}/output/results/sample/sample.GC_corrected.coverage.tsv", "right-central-even"),
        ):
            with self.subTest(method=method):
                prepared = prepare("midpoints", method, MIDPOINT_VARIABLES, MIDPOINT_PREFIX)

                self.assertEqual(prepared["command"], ["bash", f"{WORK}/gc_corrected_profile.sh", method, "4"])
                self.assertEqual(prepared["expected_output"], expected_output)
                self.assertEqual(prepared["row"]["gc_correction"], "corrected")
                self.assertEqual(prepared["row"]["midpoint_policy"], policy)
                self.assertEqual(prepared["row"]["tool_threads"], "4")


# ---------------------------------------------------------------------------
# Saved profiles (save_profile in <feature>/commands.sh)
# ---------------------------------------------------------------------------


class SaveProfileTests(unittest.TestCase):
    """The first measured output is copied into profiles/METHOD with the files concordance reads."""

    def save(self, feature, method, work):
        """Run save_profile with the output in WORK/output. Returns the saved directory."""
        results = work / "results"
        result = subprocess.run(
            ["bash", "-c", 'set -euo pipefail\nsource "$COMMON"\nsource "$COMMANDS"\nsave_profile "$METHOD"'],
            env=clean_environment(
                COMMON=str(COMMON),
                COMMANDS=str(BENCHMARKS / feature / "commands.sh"),
                METHOD=method,
                BENCH_WORK=str(work),
                BENCH_RESULTS=str(results),
                output=str(work / "output"),
            ),
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return results / "profiles" / method

    def test_positional_depth_profiles_keep_every_processed_chromosome(self):
        # chr1 has zero depth at its second base. The saved scope must still list
        # chr1, chr21 and chr22, so the comparator reads chr1 as processed
        chromosome_sizes = "chr1\t2\nchr21\t1\nchr22\t1\n"
        for method, native_name, depths in (
            ("samtools-depth", "samtools.depth.tsv", "chr1\t1\t2\nchr21\t1\t3\nchr22\t1\t4\n"),
            ("pandepth", "pandepth.SiteDepth.gz", "chr1\t0\t2\nchr1\t1\t0\nchr21\t0\t3\nchr22\t0\t4\n"),
        ):
            with self.subTest(method=method), tempfile.TemporaryDirectory() as temporary:
                work = Path(temporary)
                (work / "output").mkdir()
                (work / "chrom.sizes").write_text(chromosome_sizes)
                native = work / "output" / native_name
                if method == "pandepth":
                    native.write_bytes(gzip.compress(depths.encode()))
                else:
                    native.write_text(depths)

                saved = self.save("coverage", method, work)

                # samtools' text is compressed after timing. PanDepth's file is copied unchanged
                with gzip.open(saved / "positions.tsv.gz", "rt") as source:
                    self.assertEqual(source.read(), depths)
                if method == "pandepth":
                    self.assertEqual((saved / "positions.tsv.gz").read_bytes(), native.read_bytes())
                self.assertEqual((saved / "chrom.sizes").read_text(), chromosome_sizes)

    def test_autosome_profile_lists_only_the_autosomes_it_processed(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            (work / "output").mkdir()
            (work / "output/fcoverage.per_position.bedgraph.zst").write_text("native track\n")
            (work / "chrom.sizes").write_text("chr1\t10\nchrX\t5\nchr22\t3\nchrEBV\t2\n")

            saved = self.save("coverage", "cfdnalab-autosomes", work)

            self.assertEqual((saved / "chrom.sizes").read_text(), "chr1\t10\nchr22\t3\n")
            self.assertEqual((saved / "fcoverage.per_position.bedgraph.zst").read_text(), "native track\n")

    def test_regional_profiles_keep_only_their_window_table(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            (work / "output").mkdir()
            (work / "output/coverage.regions.bed.gz").write_text("native windows\n")
            (work / "chrom.sizes").write_text("chr1\t10\n")

            saved = self.save("coverage", "mosdepth-regions", work)

            # Regional profiles are compared on the shared regions.bed instead
            self.assertEqual(sorted(path.name for path in saved.iterdir()), ["coverage.regions.bed.gz"])

    def test_left_motif_profiles_keep_counts_and_reference_together(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            freia_table = work / "output/4_FrEIA/3_Abundances/benchmark/Genome_Lvl/T__sample.pq"
            freia_table.parent.mkdir(parents=True)
            freia_table.write_text("native left counts\n")
            (work / "output/cfdnapro-motifs.tsv").write_text("motif counts\n")
            (work / "output/reference.tsv").write_text("native reference manifest\n")
            (work / "output/contig_scope.tsv").write_text("chromosome\tretained\nchr1\tTRUE\n")

            freia = self.save("ends", "freia", work)
            cfdnapro = self.save("ends", "cfdnapro", work)

            # The rename says these are FrEIA's intact left-end counts
            self.assertEqual(sorted(path.name for path in freia.iterdir()), ["left-motifs.pq"])
            self.assertEqual((freia / "left-motifs.pq").read_text(), "native left counts\n")
            self.assertEqual((cfdnapro / "cfdnapro-motifs.tsv").read_text(), "motif counts\n")
            self.assertEqual((cfdnapro / "reference.tsv").read_text(), "native reference manifest\n")
            self.assertEqual((cfdnapro / "contig_scope.tsv").read_text(), "chromosome\tretained\nchr1\tTRUE\n")

    def test_cfdnalab_motif_store_is_copied_as_a_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            (work / "output/end_motifs.zarr/counts").mkdir(parents=True)
            (work / "output/end_motifs.zarr/counts/0").write_text("chunk\n")

            saved = self.save("ends", "cfdnalab-reference", work)

            self.assertEqual((saved / "end_motifs.zarr/counts/0").read_text(), "chunk\n")

    def test_griffin_keeps_both_tables_and_the_sample_gc_estimate(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            tables = work / "output/results/sample"
            tables.mkdir(parents=True)
            (tables / "sample.GC_corrected.coverage.tsv").write_text("corrected\n")
            (tables / "sample.uncorrected.coverage.tsv").write_text("uncorrected\n")
            (work / "output/gc-bias/GC_bias").mkdir(parents=True)
            (work / "output/gc-bias/GC_bias/sample.GC_bias.txt").write_text("bias\n")

            saved = self.save("midpoints", "griffin-gc", work)

            self.assertEqual((saved / "sample.GC_corrected.coverage.tsv").read_text(), "corrected\n")
            self.assertEqual((saved / "sample.uncorrected.coverage.tsv").read_text(), "uncorrected\n")
            self.assertEqual((saved / "sample-gc-bias/GC_bias/sample.GC_bias.txt").read_text(), "bias\n")

    def test_timed_picard_histogram_is_saved_for_later_comparison(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            (work / "output").mkdir()
            (work / "output/metrics.txt").write_text("timed histogram\n")

            saved = self.save("lengths", "picard", work)

            self.assertEqual((saved / "metrics.txt").read_text(), "timed histogram\n")


# ---------------------------------------------------------------------------
# Helper scripts for multi-command workflows
# ---------------------------------------------------------------------------


class FreiaWorkflowTests(unittest.TestCase):
    """FrEIA's filter, extraction and counting, with fakes for its samtools and Python."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name) / "scratch with spaces"
        (self.work / "output").mkdir(parents=True)
        (self.work / "sample.bam").write_bytes(b"filtered fixture\n")
        # The fake filter copies the BAM to stdout. The fake indexer only accepts
        # stdin and saves those bytes, so reading the finished file again cannot pass
        write_executable(self.work / "freia/bin/samtools", f"#!{sys.executable}\n" + textwrap.dedent(r"""
            import os
            from pathlib import Path
            import sys

            stage, *arguments = sys.argv[1:]
            with open(os.environ["CALL_LOG"], "a") as log:
                log.write(stage + "\n")
            assert arguments[arguments.index("-@") + 1] == "4"
            if stage == "view":
                assert arguments[arguments.index("-q") + 1] == "30"
                assert arguments[:-2] == ["-@", "4", "-h", "-q", "30", "-F", "4", "-F", "256", "-F", "1024", "-F", "2048"]
                sys.stdout.buffer.write(Path(arguments[-1]).read_bytes())
            elif stage == "index":
                assert arguments[-2] == "-"
                Path(arguments[-1]).write_bytes(sys.stdin.buffer.read())
            else:
                raise ValueError(stage)
        """))
        write_executable(self.work / "freia/bin/python", f"#!{sys.executable}\n" + textwrap.dedent(r"""
            import json
            import os
            from pathlib import Path
            import sys

            script = Path(sys.argv[1]).name
            with open(os.environ["CALL_LOG"], "a") as log:
                log.write(json.dumps([script, sys.argv[2:]]) + "\n")
            if script == os.environ.get("FAIL_SCRIPT"):
                raise SystemExit(5)
        """))
        self.environment = clean_environment(
            BENCH_WORK=str(self.work), BENCH_MIN_MAPQ="30", BENCH_MIN_LENGTH="30",
            BENCH_MAX_LENGTH="1000", BENCH_CONTIGS="chr1,chr2,chrEBV",
            CALL_LOG=str(self.work / "calls"),
        )

    def run_workflow(self, **variables):
        return subprocess.run(
            ["bash", str(BENCHMARKS / "ends/freia_motifs.sh"), "4"],
            env={**self.environment, **variables}, text=True, capture_output=True,
        )

    def test_filtered_bytes_are_saved_and_streamed_to_the_indexer_before_extraction(self):
        result = self.run_workflow()

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.work / "output/filtered.bam").read_bytes(), b"filtered fixture\n")
        self.assertEqual((self.work / "output/filtered.bam.bai").read_bytes(), b"filtered fixture\n")
        calls = (self.work / "calls").read_text().splitlines()
        # The filter and the indexer run together in a pipeline, so their order varies
        self.assertEqual(sorted(calls[:2]), ["index", "view"])
        output = f"{self.work}/output"
        self.assertEqual(json.loads(calls[2]), [
            "1_extract_fragment_ends.py",
            ["-b", f"{output}/filtered.bam", "-o", f"{output}/sample.pq",
             "-c", "chr1,chr2,chrEBV", "-t", "4", "-n", "3"],
        ])
        # FrEIA's upper length bound is exclusive, so 1000 bp is passed as 1001
        self.assertEqual(json.loads(calls[3]), [
            "2_data_transformation.py",
            ["-i", f"{output}/sample.pq", "-o", f"{output}/", "-st", f"{self.work}/freia-samples.tsv",
             "-fsmin", "30", "-fsmax", "1001", "-subs", "1", "--bootstrap_sample", "0"],
        ])
        # The trace shows each native call with its variables filled in
        trace = (self.work / "command-trace.log").read_text()
        self.assertIn("1_extract_fragment_ends.py", trace)
        self.assertIn("-c chr1,chr2,chrEBV -t 4 -n 3", trace)
        self.assertIn("-fsmin 30 -fsmax 1001", trace)

    def test_failed_extraction_stops_before_counting(self):
        result = self.run_workflow(FAIL_SCRIPT="1_extract_fragment_ends.py")

        self.assertEqual(result.returncode, 5)
        calls = (self.work / "calls").read_text().splitlines()
        self.assertEqual(len(calls), 3)
        self.assertNotIn("2_data_transformation.py", calls[-1])


class PicardFilteredHistogramTests(unittest.TestCase):
    """The untimed filtered Picard run reads filtered BAM from samtools."""

    def test_filtered_pipeline_writes_a_separate_histogram(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            profile = work / "results/profiles/picard"
            profile.mkdir(parents=True)
            write_executable(
                work / "samtools-env/bin/samtools",
                '#!/usr/bin/env bash\nset -euo pipefail\n'
                '[[ "$*" == "view -@ 4 -u -q 30 -f 1 -F 3852 "* ]]\n'
                'printf "filtered BAM fixture"\n',
            )
            write_executable(
                work / "picard-env/bin/picard",
                f"#!{sys.executable}\n"
                "from pathlib import Path\nimport sys\n"
                'assert sys.argv[1:4] == ["-Xmx2g", "CollectInsertSizeMetrics", "I=/dev/stdin"]\n'
                'assert "HISTOGRAM_WIDTH=1000" in sys.argv\n'
                'assert sys.stdin.read() == "filtered BAM fixture"\n'
                'output = next(argument[2:] for argument in sys.argv if argument.startswith("O="))\n'
                'Path(output).write_text("filtered histogram fixture\\n")\n',
            )
            timed = profile / "metrics.txt"
            timed.write_text("timed native histogram\n")

            result = subprocess.run(
                ["bash", str(BENCHMARKS / "lengths/picard_filtered_histogram.sh"), "4",
                 str(profile / "concordance.metrics.txt")],
                env=clean_environment(BENCH_WORK=str(work), BENCH_MIN_MAPQ="30", BENCH_MAX_LENGTH="1000"),
                text=True, capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(timed.read_text(), "timed native histogram\n")
            self.assertEqual((profile / "concordance.metrics.txt").read_text(), "filtered histogram fixture\n")


class SitePreparationTests(unittest.TestCase):
    """prepare_sites.py gives both tools the same windows in their own formats."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        (self.work / "bam.chrom.sizes").write_text("chr1\t6000\nchr2\t9000\nchrX\t100\n")
        # The first row sets the width: 3001 - 1000 = 2001, so half the width is 1000
        # BED4 rows are unstranded. A BED5 row's fifth column and a BED6 row's
        # sixth column are strands. The repeated TFBS window is kept
        (self.work / "input.bed").write_text(
            "chr2\t1000\t3001\tTFBS\n"
            "chr1\t1000\t3001\tTFBS\t-\n"
            "# a comment line\n"
            "chr1\t500\t2501\tCTCF\t0\t+\n"
            "chr1\t1000\t3001\tTFBS\t.\n"
        )

    def prepare_sites(self, *options):
        return subprocess.run(
            [sys.executable, str(BENCHMARKS / "midpoints/prepare_sites.py"), str(self.work),
             "--min-length", "100", "--max-length", "102", *options],
            text=True, capture_output=True,
        )

    def test_windows_become_sorted_cfdnalab_rows_and_griffin_anchors(self):
        result = self.prepare_sites()

        self.assertEqual(result.returncode, 0, result.stderr)
        # Width, half width, 4 sites, 2 groups and the chromosomes with sites
        self.assertEqual(result.stdout, "2001\t1000\t4\t2\tchr1,chr2\n")
        # Sorted by chromosome, start and end. The two equal chr1 TFBS windows keep
        # their input order ("-" before ".")
        self.assertEqual((self.work / "sites.bed").read_text(), (
            "chr1\t500\t2501\tCTCF\t0\t+\n"
            "chr1\t1000\t3001\tTFBS\t0\t-\n"
            "chr1\t1000\t3001\tTFBS\t0\t.\n"
            "chr2\t1000\t3001\tTFBS\t0\t.\n"
        ))
        site_lists = json.loads((self.work / "sites.yaml").read_text())["site_lists"]
        self.assertEqual(site_lists, {
            "CTCF": str(self.work / "sites/0.tsv"), "TFBS": str(self.work / "sites/1.tsv"),
        })
        # "+" and "." anchors are start + 1000. Griffin's reverse window needs the
        # other central boundary, end - 1000 = 2001, to cover the same [1000, 3001)
        self.assertEqual((self.work / "sites/0.tsv").read_text(), "Chrom\tposition\tStrand\nchr1\t1500\t+\n")
        self.assertEqual(
            (self.work / "sites/1.tsv").read_text(),
            "Chrom\tposition\tStrand\nchr1\t2001\t-\nchr1\t2000\t.\nchr2\t2000\t.\n",
        )
        # Only chromosomes with sites, in sorted order, with the BAM's lengths
        self.assertEqual((self.work / "chrom.sizes").read_text(), "chr1\t6000\nchr2\t9000\n")
        self.assertFalse((self.work / "warmup-sites.tsv").exists())
        self.assertFalse((self.work / "unit_gc_bias.tsv").exists())

    def test_griffin_jobs_also_get_the_demo_sites_and_unit_gc_weights(self):
        demo = self.work / "griffin-env/share/griffin/CTCF.hg38.1000.txt"
        demo.parent.mkdir(parents=True)
        demo.write_text("Chrom\tposition\n" + "chr2\t3000\n" * 1000)

        result = self.prepare_sites("--griffin")

        self.assertEqual(result.returncode, 0, result.stderr)
        demo_rows = (self.work / "warmup-sites.tsv").read_text().splitlines()
        self.assertEqual(demo_rows[0], "Chrom\tposition\tStrand")
        self.assertEqual(demo_rows[1:], ["chr2\t3000\t."] * 1000)
        self.assertEqual((self.work / "warmup.chrom.sizes").read_text(), "chr2\t9000\n")
        self.assertEqual(
            json.loads((self.work / "warmup-sites.yaml").read_text()),
            {"site_lists": {"CTCF-demo": str(self.work / "warmup-sites.tsv")}},
        )
        # Lengths 100-102 with 0..length GC bases: 101 + 102 + 103 = 306 rows
        bias_rows = (self.work / "unit_gc_bias.tsv").read_text().splitlines()
        self.assertEqual(bias_rows[0], "length\tnum_GC\tsmoothed_GC_bias")
        self.assertEqual(len(bias_rows) - 1, 306)
        self.assertEqual((bias_rows[1], bias_rows[-1]), ("100\t0\t1.0", "102\t102\t1.0"))

    def test_incomplete_demo_or_empty_bed_fails(self):
        demo = self.work / "griffin-env/share/griffin/CTCF.hg38.1000.txt"
        demo.parent.mkdir(parents=True)
        demo.write_text("Chrom\tposition\n" + "chr2\t3000\n" * 999)
        self.assertNotEqual(self.prepare_sites("--griffin").returncode, 0)

        (self.work / "input.bed").write_text("# only a comment\n")
        result = self.prepare_sites()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No midpoint sites", result.stderr)


# ---------------------------------------------------------------------------
# Shared steps in common.sh
# ---------------------------------------------------------------------------


class JobSettingsTests(unittest.TestCase):
    """check_job_settings rejects invalid settings before anything is copied."""

    def check(self, **overrides):
        """Run check_job_settings with valid settings, changed by `overrides`.

        An override of None removes that variable. BENCH_RESULTS="existing" points
        it at a directory that already exists.
        """
        with tempfile.TemporaryDirectory() as temporary:
            environment = clean_environment(
                COMMON=str(COMMON),
                SLURM_JOB_ID="123", SLURM_JOB_NUM_NODES="1", SLURM_CPUS_PER_TASK="4",
                BENCH_BAM="/unused/input.bam", TMPDIR=temporary,
                BENCH_RESULTS=temporary + "/results", BENCH_CORES="4",
                BENCH_METHODS="cfdnalab", BENCH_WARMUPS="1", BENCH_REPETITIONS="3",
                BENCH_MIN_MAPQ="30", BENCH_MIN_LENGTH="30", BENCH_MAX_LENGTH="1000",
                BENCH_SETUP_PYTHON=sys.executable,
            )
            for name, value in overrides.items():
                if value is None:
                    environment.pop(name, None)
                elif name == "BENCH_RESULTS" and value == "existing":
                    environment[name] = temporary
                else:
                    environment[name] = value
            return subprocess.run(
                ["bash", "-c",
                 'set -euo pipefail\nsource "$COMMON"\ncheck_job_settings lengths\n'
                 'printf "%s\\n" "${methods[@]}" "warmups=$BENCH_WARMUPS" "repetitions=$BENCH_REPETITIONS"'],
                env=environment,
                text=True, capture_output=True,
            )

    def test_valid_settings_list_the_methods(self):
        result = self.check(BENCH_METHODS="cfdnalab samtools", BENCH_COMPARE_PYTHON=sys.executable)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "cfdnalab\nsamtools\nwarmups=1\nrepetitions=3\n")

    def test_a_pilot_runs_each_method_once_without_warmup(self):
        result = self.check(BENCH_PILOT="1")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "cfdnalab\nwarmups=0\nrepetitions=1\n")

    def test_invalid_settings_are_rejected_with_an_explanation(self):
        for overrides, message in (
            # A repeated method would double its repetitions
            ({"BENCH_METHODS": "cfdnalab cfdnalab"}, "Duplicate method"),
            # Each job measures a single core count
            ({"BENCH_CORES": "1 4"}, "single positive core count"),
            ({"BENCH_CORES": "8"}, "exceeds allocated CPUs"),
            ({"BENCH_METHODS": None}, "BENCH_METHODS"),
            ({"TMPDIR": None}, "TMPDIR"),
            ({"BENCH_RESULTS": "existing"}, "already exists"),
            ({"SLURM_JOB_NUM_NODES": "2"}, "one allocated node"),
            ({"BENCH_PILOT": "yes"}, "BENCH_PILOT"),
            ({"BENCH_REPETITIONS": "0"}, "at least one repetition"),
            ({"BENCH_MIN_LENGTH": "9"}, "10-50000"),
            ({"BENCH_MIN_MAPQ": "-1"}, "BENCH_MIN_MAPQ"),
            # Two methods are compared at the end, which needs the comparison Python
            ({"BENCH_METHODS": "cfdnalab samtools"}, "BENCH_COMPARE_PYTHON"),
        ):
            with self.subTest(overrides=overrides):
                result = self.check(**overrides)

                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                self.assertEqual(result.stdout, "")


class TimingRowTests(unittest.TestCase):
    """write_timing_row writes the row by column name and refuses incomplete rows."""

    def write(self, row_assignment):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            result = subprocess.run(
                ["bash", "-c",
                 'set -euo pipefail\nsource "$COMMON"\n'
                 'timing_columns=(method wall_seconds step_exit_code)\n'
                 'create_results_directory\n' + row_assignment + '\nwrite_timing_row\n'],
                env=clean_environment(COMMON=str(COMMON), BENCH_RESULTS=str(results / "results")),
                text=True, capture_output=True,
            )
            timings = (results / "results/timings.tsv").read_text()
            return result, timings

    def test_values_are_written_in_column_order(self):
        result, timings = self.write("row=([step_exit_code]=0 [method]=samtools [wall_seconds]=1.5)")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(timings, "method\twall_seconds\tstep_exit_code\nsamtools\t1.5\t0\n")

    def test_missing_or_unknown_columns_stop_without_writing(self):
        for assignment, message in (
            ("row=([method]=samtools [wall_seconds]=1.5)", "No value for timing column 'step_exit_code'"),
            ("row=([method]=samtools [wall_seconds]=1.5 [step_exit_code]=0 [scop]=all)", "'scop' is not a timing column"),
        ):
            with self.subTest(assignment=assignment):
                result, timings = self.write(assignment)

                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                self.assertEqual(timings, "method\twall_seconds\tstep_exit_code\n")


# ---------------------------------------------------------------------------
# measure_command.sh
# ---------------------------------------------------------------------------


@requires_measurement
class MeasureCommandTests(unittest.TestCase):
    """The measurement step, run on a single CPU as the fake srun would bind it."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        (self.work / "output").mkdir()
        self.cpu = single_cpu()
        self.environment = clean_environment(
            BENCH_WORK=str(self.work), BENCH_CORES="1", BENCH_NODE=hostname(),
            BENCH_RUN_TIMEOUT_SECONDS="0", SLURM_STEP_ID="7",
        )

    def measure(self, command, expected_output="output/result.txt", empty_output_allowed="false", **variables):
        """Run measure_command.sh with `command` in bash, bound to a single CPU."""
        cpu = self.cpu
        return subprocess.run(
            ["bash", str(MEASURE_COMMAND), str(self.work / "tool-stdout.log"),
             str(self.work / expected_output), empty_output_allowed, "--", "bash", "-c", command],
            env={**self.environment, **variables},
            text=True, capture_output=True,
            preexec_fn=lambda: os.sched_setaffinity(0, {cpu}),
        )

    def record(self):
        return dict(read_record(self.work / "measurement.tsv"))

    def test_successful_run_records_placement_time_flush_and_size(self):
        result = self.measure(
            'printf "native result" > "$BENCH_WORK/output/result.txt"; echo printed; echo warning >&2'
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        names = [name for name, _value in read_record(self.work / "measurement.tsv")]
        self.assertEqual(names, [
            "step_id", "node", "cpu_model", "cpu_ids", "physical_core_ids", "started_utc",
            "wall_seconds", "user_seconds", "system_seconds", "max_rss_kib", "tool_exit_code",
            "flush_seconds", "process_plus_flush_seconds", "output_bytes",
        ])
        record = self.record()
        self.assertEqual(record["step_id"], "7")
        self.assertEqual(record["node"], hostname())
        self.assertEqual(record["cpu_ids"], str(self.cpu))
        self.assertEqual(record["tool_exit_code"], "0")
        self.assertGreater(int(record["max_rss_kib"]), 0)
        # The combined time is wall + flush, rounded to 3 decimals
        self.assertAlmostEqual(
            float(record["process_plus_flush_seconds"]),
            float(record["wall_seconds"]) + float(record["flush_seconds"]), places=3,
        )
        # The output directory holds at least the 13 bytes of "native result"
        self.assertGreaterEqual(int(record["output_bytes"]), 13)
        self.assertEqual((self.work / "tool-stdout.log").read_text(), "printed\n")
        self.assertEqual((self.work / "tool-stderr.log").read_text(), "warning\n")

    def test_exit_and_signal_statuses_are_recorded_without_flush_or_size(self):
        for command, expected_status in (("exit 3", 3), ("kill -9 $$", 137)):
            with self.subTest(command=command):
                (self.work / "measurement.tsv").unlink(missing_ok=True)

                result = self.measure(command)

                self.assertEqual(result.returncode, expected_status)
                record = self.record()
                self.assertEqual(record["tool_exit_code"], str(expected_status))
                self.assertIn("wall_seconds", record)
                self.assertNotIn("flush_seconds", record)
                self.assertNotIn("output_bytes", record)

    def test_time_limit_wraps_the_command_inside_gnu_time(self):
        # Only the outbound timeout program is faked, avoiding clock-dependent waits
        fake_bin = self.work / "bin"
        write_executable(
            fake_bin / "timeout",
            '#!/usr/bin/env bash\nprintf "%s\\n" "$@" > "$BENCH_WORK/timeout.args"\nexit 124\n',
        )

        result = self.measure(
            "exit 0", BENCH_RUN_TIMEOUT_SECONDS="43",
            PATH=str(fake_bin) + ":" + self.environment["PATH"],
        )

        self.assertEqual(result.returncode, 124)
        self.assertEqual(
            (self.work / "timeout.args").read_text().splitlines(),
            ["--verbose", "--kill-after=60", "43", "bash", "-c", "exit 0"],
        )
        self.assertEqual(self.record()["tool_exit_code"], "124")

    def test_missing_or_empty_output_fails_after_measuring(self):
        missing = self.measure("true", expected_output="output/missing.txt")
        self.assertEqual(missing.returncode, 1)
        self.assertIn("expected output is missing", missing.stderr)
        # The command itself succeeded and was fully measured
        self.assertEqual(self.record()["tool_exit_code"], "0")
        self.assertIn("output_bytes", self.record())

        for allowed, expected_status in (("false", 1), ("true", 0)):
            with self.subTest(empty_output_allowed=allowed):
                (self.work / "measurement.tsv").unlink()
                result = self.measure(': > "$BENCH_WORK/output/empty.tsv"', "output/empty.tsv", allowed)
                self.assertEqual(result.returncode, expected_status, result.stderr)

    def test_wrong_placement_stops_before_running_the_command(self):
        # The step can use a single CPU, but the job expects two cores or another node
        for variables in ({"BENCH_CORES": "2"}, {"BENCH_NODE": "another-node"}):
            with self.subTest(variables=variables):
                (self.work / "measurement.tsv").unlink(missing_ok=True)

                result = self.measure('touch "$BENCH_WORK/command-ran"', **variables)

                self.assertEqual(result.returncode, 1)
                self.assertIn("Unexpected CPU placement", result.stderr)
                self.assertFalse((self.work / "command-ran").exists())
                record = self.record()
                self.assertEqual(record["step_id"], "7")
                self.assertNotIn("started_utc", record)
                self.assertNotIn("tool_exit_code", record)

    def test_missing_separator_is_a_usage_error(self):
        result = subprocess.run(
            ["bash", str(MEASURE_COMMAND), "stdout.log", "output", "false", "true"],
            env=self.environment, text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("Usage", result.stderr)


# ---------------------------------------------------------------------------
# Whole jobs with fake srun and tools
# ---------------------------------------------------------------------------


# Stands in for srun. It logs its options, gives each step a new SLURM_STEP_ID and
# runs the step's program bound to a single CPU, like an exclusive single-core step
FAKE_SRUN = textwrap.dedent(r"""
    import os
    from pathlib import Path
    import subprocess
    import sys

    arguments = sys.argv[1:]
    log = Path(os.environ["SRUN_LOG"])
    with log.open("a") as destination:
        destination.write("\t".join(arguments) + "\n")
    step_id = len(log.read_text().splitlines()) - 1
    # srun's own options come before the program it starts
    program = arguments[arguments.index("bash"):]
    cpu = min(os.sched_getaffinity(0))
    result = subprocess.run(
        program,
        env={**os.environ, "SLURM_STEP_ID": str(step_id)},
        preexec_fn=lambda: os.sched_setaffinity(0, {cpu}),
    )
    raise SystemExit(result.returncode)
""")


def write_archive(path, environment_directory):
    """Pack a fake Conda environment the way install_tools.sh's archives are packed."""
    # conda-unpack is run by the setup Python. A comment is a valid Python script
    (environment_directory / "bin").mkdir(parents=True, exist_ok=True)
    (environment_directory / "bin/conda-unpack").write_text("# Relocation fixture\n")
    with tarfile.open(path, "w:gz") as destination:
        destination.add(environment_directory, arcname=".")


@requires_measurement
class LengthsJobTests(unittest.TestCase):
    """Run benchmark_lengths.sh end to end with fake tools and a fake srun."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        fake_bin = self.root / "bin"
        (self.root / "scratch").mkdir()
        write_executable(fake_bin / "srun", f"#!{sys.executable}\n" + FAKE_SRUN)
        write_executable(fake_bin / "ldd", "#!/usr/bin/env bash\nprintf '\\tstatically linked\\n'\n")
        write_executable(fake_bin / "scontrol", "#!/usr/bin/env bash\necho fixture configuration\n")
        # The fake cfDNAlab logs which configuration ran and writes a labeled table
        write_executable(fake_bin / "cfdna", textwrap.dedent(r"""
            #!/usr/bin/env bash
            if [[ "$1" == --version ]]; then
                if [[ "${FAIL_CFDNA_VERSION:-0}" == 1 ]]; then echo "cfdna version fixture failed" >&2; exit 6; fi
                echo "cfdna 0.0.0-fixture"
                exit 0
            fi
            method=cfdnalab
            if [[ " $* " == *" --require-proper-pair "* ]]; then method=cfdnalab-proper; fi
            printf '%s\n' "$method" >> "$CALL_LOG"
            if [[ "$method" == cfdnalab && "${FAIL_CFDNALAB:-0}" == 1 ]]; then
                echo "cfdnalab fixture failed" >&2
                exit 3
            fi
            while (($#)); do
                if [[ "$1" == --output-dir ]]; then output=$2; fi
                shift
            done
            printf 'length counts from %s\n' "$method" > "$output/length_counts.tsv.zst"
        """).lstrip())
        # The fake samtools reads a two-contig header. Its stats report counts its
        # own calls, so the saved report shows which repetition it came from
        samtools_environment = self.root / "samtools-env"
        write_executable(samtools_environment / "bin/samtools", textwrap.dedent(r"""
            #!/usr/bin/env bash
            case "$1" in
                idxstats) printf 'chr1\t100\t0\t0\nchr2\t50\t0\t0\n*\t0\t0\t0\n' ;;
                --version-only) echo 1.25-fixture ;;
                stats)
                    printf 'samtools\n' >> "$CALL_LOG"
                    if [[ "${FAIL_SAMTOOLS_STATS:-0}" == 1 ]]; then echo "samtools stats fixture failed" >&2; exit 4; fi
                    calls=$(grep -c '^samtools$' "$CALL_LOG")
                    printf 'IS\t30\t%s\t%s\t0\t0\n' "$calls" "$calls"
                    ;;
                *) exit 9 ;;
            esac
        """).lstrip())
        write_archive(self.root / "samtools.tar.gz", samtools_environment)
        write_executable(self.root / "comparison-python", f"#!{sys.executable}\n" + textwrap.dedent(r"""
            import os
            from pathlib import Path
            import sys

            with open(os.environ["CALL_LOG"], "a") as log:
                log.write("concordance\n")
            Path(os.environ["COMPARISON_ARGUMENTS"]).write_text("\n".join(sys.argv[1:]) + "\n")
            print("mae_frequency\n0")
        """))
        self.bam = self.root / "input.bam"
        self.bam.write_text("BAM contents are only read by the fake tools\n")
        Path(str(self.bam) + ".bai").touch()
        self.environment = clean_environment(
            PATH=str(fake_bin) + ":" + os.environ["PATH"],
            SLURM_SUBMIT_DIR=str(BENCHMARKS.parent),
            SLURM_JOB_ID="123",
            SLURM_JOB_NUM_NODES="1",
            SLURM_CPUS_PER_TASK="1",
            TMPDIR=str(self.root / "scratch"),
            BENCH_BAM=str(self.bam),
            BENCH_CFDNA=str(fake_bin / "cfdna"),
            BENCH_SAMTOOLS_ARCHIVE=str(self.root / "samtools.tar.gz"),
            BENCH_COMPARE_PYTHON=str(self.root / "comparison-python"),
            BENCH_SETUP_PYTHON=sys.executable,
            BENCH_CORES="1",
            BENCH_PILOT="0",
            SRUN_LOG=str(self.root / "srun.log"),
            CALL_LOG=str(self.root / "calls"),
            COMPARISON_ARGUMENTS=str(self.root / "comparison-arguments"),
        )

    def run_job(self, results_name, **settings):
        results = self.root / results_name
        result = subprocess.run(
            ["bash", str(BENCHMARKS / "lengths/benchmark_lengths.sh")],
            env={**self.environment, "BENCH_RESULTS": str(results), **settings},
            text=True, capture_output=True,
        )
        return result, results

    def calls(self):
        return (self.root / "calls").read_text().splitlines()

    def test_failed_method_is_recorded_and_the_other_methods_finish(self):
        result, results = self.run_job(
            "results", BENCH_METHODS="cfdnalab cfdnalab-proper samtools",
            BENCH_WARMUPS="0", BENCH_REPETITIONS="2", FAIL_CFDNALAB="1",
        )

        # The job reports the failure only after the other methods and concordance
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertEqual(self.calls(), [
            "cfdnalab", "cfdnalab-proper", "cfdnalab-proper", "samtools", "samtools", "concordance",
        ])
        rows = read_timings(results / "timings.tsv")
        self.assertEqual([row["method"] for row in rows], [
            "cfdnalab", "cfdnalab-proper", "cfdnalab-proper", "samtools", "samtools",
        ])
        self.assertEqual([row["run_order"] for row in rows], ["1", "2", "3", "4", "5"])
        self.assertEqual([row["repetition"] for row in rows], ["1", "1", "2", "1", "2"])
        self.assertEqual([row["step_exit_code"] for row in rows], ["3", "0", "0", "0", "0"])
        self.assertEqual([row["tool_exit_code"] for row in rows], ["3", "0", "0", "0", "0"])
        # Each run was its own step, numbered 0-4 by the fake srun
        self.assertEqual([row["step_id"] for row in rows], ["0", "1", "2", "3", "4"])
        for row in rows:
            self.assertEqual((row["warmup"], row["cores"], row["job_id"], row["pilot"]), ("false", "1", "123", "0"))
            self.assertEqual(row["tool_threads"], "1")
            self.assertEqual(row["node"], hostname())
            self.assertGreaterEqual(float(row["wall_seconds"]), 0)
        # The failed run keeps its command and error, and saved no profile
        self.assertEqual((results / "runs/1/stderr.log").read_text(), "cfdnalab fixture failed\n")
        self.assertIn(" lengths --bam ", (results / "runs/1/command.sh").read_text())
        self.assertIn(" --n-threads 1 --chromosomes all", (results / "runs/1/command.sh").read_text())
        self.assertFalse((results / "profiles/cfdnalab").exists())
        # The first measured repetition of each successful method is saved
        self.assertEqual(
            (results / "profiles/cfdnalab-proper/length_counts.tsv.zst").read_text(),
            "length counts from cfdnalab-proper\n",
        )
        self.assertEqual((results / "profiles/samtools/stats.txt").read_text(), "IS\t30\t1\t1\t0\t0\n")
        # samtools prints its report, which is also kept as that run's stdout
        self.assertEqual((results / "runs/5/stdout.log").read_text(), "IS\t30\t2\t2\t0\t0\n")
        self.assertEqual((results / "concordance.tsv").read_text(), "mae_frequency\n0\n")
        # Every run was an exact single-core step of this job on the reserved node
        steps = (self.root / "srun.log").read_text().splitlines()
        self.assertEqual(len(steps), 5)
        for step in steps:
            options = step.split("\t")
            self.assertEqual(options[:8], [
                "--jobid=123", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
                "--threads-per-core=1", "--cpu-bind=threads", "--export=ALL",
            ])
            self.assertEqual(options[8], "bash")
            self.assertTrue(options[9].endswith("/measure_command.sh"))
        self.assertIn("Saved run 5", (results / "progress.log").read_text())

    def test_warmups_are_recorded_and_the_first_measured_output_is_saved(self):
        result, results = self.run_job(
            "results", BENCH_METHODS="cfdnalab samtools",
            BENCH_WARMUPS="1", BENCH_REPETITIONS="2", BENCH_MAX_LENGTH="200",
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls(), ["cfdnalab"] * 3 + ["samtools"] * 3 + ["concordance"])
        rows = read_timings(results / "timings.tsv")
        self.assertEqual([row["warmup"] for row in rows], ["true", "false", "false"] * 2)
        self.assertEqual([row["repetition"] for row in rows], ["1", "2", "3"] * 2)
        self.assertEqual({row["step_exit_code"] for row in rows}, {"0"})
        # The warmup was samtools call 1, so the first measured report is call 2
        self.assertEqual((results / "profiles/samtools/stats.txt").read_text(), "IS\t30\t2\t2\t0\t0\n")
        # --length-bins has an exclusive end: 200 bp becomes 201
        self.assertIn("--length-bins 30:201:1", (results / "runs/1/command.sh").read_text())
        self.assertEqual((self.root / "comparison-arguments").read_text().splitlines(), [
            str(BENCHMARKS / "lengths/compare_lengths.py"), str(results / "profiles"),
            "--min-length", "30", "--max-length", "200",
        ])
        self.assertIn("warmups\t1\n", (results / "run.tsv").read_text())

    def test_failure_of_the_last_method_keeps_earlier_results(self):
        result, results = self.run_job(
            "results", BENCH_METHODS="cfdnalab cfdnalab-proper samtools",
            BENCH_WARMUPS="0", BENCH_REPETITIONS="2", FAIL_SAMTOOLS_STATS="1",
        )

        self.assertEqual(result.returncode, 1, result.stderr)
        rows = read_timings(results / "timings.tsv")
        self.assertEqual([row["step_exit_code"] for row in rows], ["0", "0", "0", "0", "4"])
        self.assertEqual((results / "runs/5/stderr.log").read_text(), "samtools stats fixture failed\n")
        self.assertTrue((results / "profiles/cfdnalab/length_counts.tsv.zst").exists())
        self.assertTrue((results / "profiles/cfdnalab-proper/length_counts.tsv.zst").exists())
        self.assertFalse((results / "profiles/samtools").exists())
        self.assertTrue((results / "concordance.tsv").exists())
        self.assertIn("Saved run 5", (results / "progress.log").read_text())

    def test_a_broken_tool_fails_before_the_bam_is_copied(self):
        result, results = self.run_job(
            "results", BENCH_METHODS="cfdnalab samtools", FAIL_CFDNA_VERSION="1",
        )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("cfdna version fixture failed", result.stderr)
        # Setup records exist, but no run has started
        self.assertTrue((results / "run.tsv").exists())
        self.assertEqual(read_timings(results / "timings.tsv"), [])
        progress = (results / "progress.log").read_text()
        self.assertIn("Checking tool versions", progress)
        self.assertNotIn("Copying BAM", progress)
        self.assertFalse((self.root / "srun.log").exists())


@requires_measurement
class GriffinJobTests(unittest.TestCase):
    """Run benchmark_midpoints.sh for griffin-gc with fake Griffin scripts and srun."""

    def test_demo_warmup_is_a_separate_step_before_three_measured_repetitions(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / "bin"
            (root / "scratch").mkdir()
            write_executable(fake_bin / "srun", f"#!{sys.executable}\n" + FAKE_SRUN)
            write_executable(fake_bin / "scontrol", "#!/usr/bin/env bash\necho fixture\n")
            samtools_environment = root / "samtools-env"
            write_executable(samtools_environment / "bin/samtools", '#!/usr/bin/env bash\nprintf "chr1\\t6000\\t0\\t0\\n"\n')
            write_archive(root / "samtools.tar.gz", samtools_environment)
            griffin_environment = root / "griffin-env"
            scripts = griffin_environment / "share/griffin"
            scripts.mkdir(parents=True)
            (scripts / "version.txt").write_text("v0.2.0 fixture\n")
            for script in ("griffin_GC_counts.py", "griffin_GC_bias.py"):
                (scripts / script).touch()
            # Repeated rows are intentional. The demo must keep exactly 1000 rows
            (scripts / "CTCF.hg38.1000.txt").write_text("Chrom\tposition\n" + "chr1\t3000\n" * 1000)
            # The fake Griffin Python logs every native call. Its merge step writes
            # both tables, labeled with the number of the merge call
            write_executable(griffin_environment / "bin/python", f"#!{sys.executable}\n" + textwrap.dedent(r"""
                import os
                from pathlib import Path
                import sys

                if sys.argv[1] == "-c":
                    # The setup import check does not process any data
                    raise SystemExit(0)
                script = Path(sys.argv[1]).name
                arguments = sys.argv[2:]
                log = Path(os.environ["NATIVE_CALL_LOG"])
                with log.open("a") as destination:
                    destination.write("\t".join([script, *arguments]) + "\n")
                if script == "griffin_merge_sites.py":
                    merge_calls = sum(line.startswith(script) for line in log.read_text().splitlines())
                    tables = Path(arguments[arguments.index("--results_dir") + 1]) / "sample"
                    tables.mkdir(parents=True, exist_ok=True)
                    (tables / "sample.GC_corrected.coverage.tsv").write_text(f"corrected from merge {merge_calls}\n")
                    (tables / "sample.uncorrected.coverage.tsv").write_text(f"uncorrected from merge {merge_calls}\n")
            """))
            write_archive(root / "griffin.tar.gz", griffin_environment)
            bam = root / "input.bam"
            bam.touch()
            Path(str(bam) + ".bai").touch()
            fasta = root / "reference.fa"
            fasta.touch()
            Path(str(fasta) + ".fai").write_text("chr1\t6000\t0\t6000\t6001\n")
            intervals = root / "production.bed"
            production_sites = (
                "chr1\t1000\t3001\tproduction\t+\nchr1\t1000\t3001\tproduction\t+\n"
                "chr1\t1000\t3001\treverse\t-\nchr1\t1000\t3001\tunstranded\t.\n"
            )
            intervals.write_text(production_sites)
            references = root / "references"
            references.mkdir()
            for filename in EXCLUSION_FILENAMES:
                (references / filename).touch()
            (references / "mappability.bw").touch()
            (references / "gc-regions.bed").write_text("chr1\t0\t6000\n")
            gc_frequencies = references / "genome_GC_frequency"
            gc_frequencies.mkdir()
            for fragment_length in range(15, 501):
                (gc_frequencies / f"gc-regions.{fragment_length}bp.GC_frequency.txt").write_text("fixture\n")
            results = root / "results"

            result = subprocess.run(
                ["bash", str(BENCHMARKS / "midpoints/benchmark_midpoints.sh")],
                env=clean_environment(
                    PATH=str(fake_bin) + ":" + os.environ["PATH"],
                    SLURM_SUBMIT_DIR=str(BENCHMARKS.parent),
                    SLURM_JOB_ID="123",
                    SLURM_JOB_NUM_NODES="1",
                    SLURM_CPUS_PER_TASK="1",
                    TMPDIR=str(root / "scratch"),
                    BENCH_BAM=str(bam),
                    BENCH_REF_FASTA=str(fasta),
                    BENCH_INTERVALS=str(intervals),
                    BENCH_RESULTS=str(results),
                    BENCH_SAMTOOLS_ARCHIVE=str(root / "samtools.tar.gz"),
                    BENCH_GRIFFIN_ARCHIVE=str(root / "griffin.tar.gz"),
                    BENCH_GRIFFIN_REF_DIR=str(references),
                    BENCH_GRIFFIN_MAPPABILITY_BW=str(references / "mappability.bw"),
                    BENCH_GRIFFIN_GC_REGIONS=str(references / "gc-regions.bed"),
                    BENCH_GRIFFIN_REF_GC=str(gc_frequencies),
                    BENCH_METHODS="griffin-gc",
                    BENCH_CORES="1",
                    # Griffin's full-input warmups are replaced by the demo
                    BENCH_WARMUPS="10",
                    BENCH_REPETITIONS="3",
                    BENCH_PILOT="0",
                    BENCH_SETUP_PYTHON=sys.executable,
                    SRUN_LOG=str(root / "srun.log"),
                    NATIVE_CALL_LOG=str(root / "native-calls"),
                ),
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            # The demo is the first of four steps. It runs Griffin's profiling helper
            # directly, while the three measured runs go through measure_command.sh
            steps = [line.split("\t") for line in (root / "srun.log").read_text().splitlines()]
            self.assertEqual(len(steps), 4)
            self.assertTrue(steps[0][9].endswith("/griffin_profile.sh"))
            self.assertIn("--cpus-per-task=1", steps[0])
            self.assertTrue(all(step[9].endswith("/measure_command.sh") for step in steps[1:]))
            calls = [line.split("\t") for line in (root / "native-calls").read_text().splitlines()]
            demo_calls = [call for call in calls if any(argument.endswith("warmup-sites.yaml") for argument in call)]
            self.assertEqual([call[0] for call in demo_calls], ["griffin_coverage.py", "griffin_merge_sites.py"])
            for call in demo_calls:
                # A 2001 bp window is [-1000, 1001) around Griffin's anchor
                self.assertEqual(call[call.index("--norm_window") + 1:call.index("--norm_window") + 3], ["-1000", "1001"])
                self.assertEqual(argument_after(call, "--CPU"), "1")
            self.assertEqual(argument_after(demo_calls[0], "--size_range"), "100")
            self.assertEqual(argument_after(demo_calls[0], "--map_quality"), "20")
            # Every measured repetition estimates GC bias from scratch, then profiles
            measured_scripts = [call[0] for call in calls if call not in demo_calls]
            self.assertEqual(measured_scripts, [
                "griffin_GC_counts.py", "griffin_GC_bias.py", "griffin_coverage.py", "griffin_merge_sites.py",
            ] * 3)
            rows = read_timings(results / "timings.tsv")
            self.assertEqual([row["repetition"] for row in rows], ["1", "2", "3"])
            for row in rows:
                self.assertEqual(
                    (row["warmup"], row["window_bp"], row["gc_correction"], row["tool_threads"], row["step_exit_code"]),
                    ("false", "2001", "corrected", "1", "0"),
                )
            # The demo's merge was call 1, so the first measured repetition is call 2
            saved = results / "profiles/griffin-gc"
            self.assertEqual((saved / "sample.GC_corrected.coverage.tsv").read_text(), "corrected from merge 2\n")
            self.assertEqual((saved / "sample.uncorrected.coverage.tsv").read_text(), "uncorrected from merge 2\n")
            self.assertTrue((saved / "sample-gc-bias").is_dir())
            # The first measured run's trace has its GC estimation and both profiling
            # calls, with values filled in, but not the demo's calls
            trace = (results / "runs/1/command-trace.log").read_text()
            self.assertIn("griffin_GC_counts.py", trace)
            self.assertIn("--size_range 15 500 --CPU 1", trace)
            self.assertIn("griffin_merge_sites.py", trace)
            self.assertNotIn("warmup-sites.yaml", trace)
            self.assertEqual(intervals.read_text(), production_sites)
            self.assertEqual(len((results / "sites.bed").read_text().splitlines()), 4)
            self.assertIn("warmups\t0\n", (results / "run.tsv").read_text())


if __name__ == "__main__":
    unittest.main()
