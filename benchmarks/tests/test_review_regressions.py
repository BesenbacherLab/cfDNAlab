"""Hand-derived review regressions, prepared for the owner and not executed.

Timing integration cases require Linux GNU time. Native histogram readers also
require the comparison environment. No Slurm allocation or real BAM is used.
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
import tarfile
import tempfile
import textwrap
from types import SimpleNamespace
import unittest
from unittest import mock


BENCHMARKS = Path(__file__).resolve().parents[1]


class SetupAndScopeTests(unittest.TestCase):
    """Reject incomplete setup before staging or decompressing native data."""

    def test_each_runner_requires_an_explicit_method_block_before_staging(self):
        for feature, filename in (('coverage', 'benchmark_coverage.sh'), ('lengths', 'benchmark_lengths.sh'),
                                  ('ends', 'benchmark_end_motifs.sh'), ('midpoints', 'benchmark_midpoints.sh')):
            with self.subTest(feature=feature), tempfile.TemporaryDirectory() as temporary:
                environment = {**os.environ, 'SLURM_SUBMIT_DIR': str(BENCHMARKS.parent),
                               'SLURM_JOB_ID': '123', 'BENCH_CORES': '4', 'BENCH_BAM': '/unused/input.bam',
                               'BENCH_INDEX': '/unused/input.bam.bai', 'BENCH_INTERVALS': '/unused/sites.bed',
                               'BENCH_SAMTOOLS_ARCHIVE': '/unused/samtools.tar.gz',
                               'BENCH_LOCAL_ROOT': temporary, 'BENCH_RESULTS': temporary + '/results'}
                environment.pop('BENCH_METHODS', None)
                result = subprocess.run(['bash', str(BENCHMARKS / feature / filename)],
                                        env=environment, text=True, capture_output=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('BENCH_METHODS', result.stderr)
                self.assertFalse((Path(temporary) / 'results').exists())

    def test_scope_mismatch_fails_before_track_loading_and_keeps_previous_histogram(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles = root / 'profiles'
            for name, scope in (('cfdnalab', 'chr1\t100\nchr2\t100\n'), ('mosdepth', 'chr1\t100\n')):
                directory = profiles / name
                directory.mkdir(parents=True)
                (directory / 'chrom.sizes').write_text(scope)
            # No track files exist: scope validation must precede decompression
            (root / 'residuals.tsv').write_text('previous successful histogram\n')
            result = subprocess.run([sys.executable, str(BENCHMARKS / 'coverage/compare_coverage.py'),
                                     str(profiles)], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Scopes differ', result.stderr)
            self.assertEqual((root / 'residuals.tsv').read_text(), 'previous successful histogram\n')

    def test_duplicate_method_blocks_are_rejected_instead_of_doubling_repetitions(self):
        result = subprocess.run(['bash', '-euc', 'source "$BENCH_COMMON"; benchmark_controls'], env={
            **os.environ, 'BENCH_COMMON': str(BENCHMARKS / 'common.sh'),
            'BENCH_CORES': '4', 'BENCH_METHODS': 'cfdnalab cfdnalab',
            'BENCH_PILOT': '0', 'BENCH_EXPECTED_NODE': '', 'BENCH_RUN_TIMEOUT_SECONDS': '0',
            'BENCH_SETUP_PYTHON': sys.executable,
        }, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Duplicate method', result.stderr)

    def test_invalid_track_after_scope_preflight_keeps_previous_histogram(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profiles = root / 'profiles'
            for name in ('samtools-depth', 'pandepth'):
                directory = profiles / name
                directory.mkdir(parents=True)
                (directory / 'chrom.sizes').write_text('chr1\t100\n')
                with gzip.open(directory / 'positions.tsv.gz', 'wt') as destination:
                    # A single covered base cannot be a complete 100 bp dense track
                    destination.write('chr1\t1\t2\n' if name == 'samtools-depth' else 'chr1\t0\t2\n')
            (root / 'residuals.tsv').write_text('previous successful histogram\n')
            result = subprocess.run([sys.executable, str(BENCHMARKS / 'coverage/compare_coverage.py'),
                                     str(profiles)], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Incomplete native track', result.stderr)
            self.assertEqual((root / 'residuals.tsv').read_text(), 'previous successful histogram\n')
            self.assertEqual(list(root.glob('.residuals-*')), [])


@unittest.skipUnless(sys.platform.startswith('linux'), 'Requires Linux GNU time and sync -f')
class RunLifecycleTests(unittest.TestCase):
    """Exercise native statuses and retained failure diagnostics through the helpers."""

    def test_exit_and_signal_statuses_are_recorded_without_fabricating_success(self):
        for command, expected_status in (('exit 3', 3), ('kill -9 $$', 137)):
            with self.subTest(command=command), tempfile.TemporaryDirectory() as temporary:
                work = Path(temporary)
                (work / 'output').mkdir()
                script = '''source "$BENCH_COMMON"
benchmark_time_command fixture "$BENCH_WORK/output" "$BENCH_WORK/tool.log" "$BENCH_WORK/tool-error.log" bash -c "$NATIVE_COMMAND"
'''
                result = subprocess.run(['bash', '-euo', 'pipefail', '-c', script], env={
                    **os.environ, 'BENCH_COMMON': str(BENCHMARKS / 'common.sh'), 'BENCH_WORK': str(work),
                    'BENCH_RUN_TIMEOUT_SECONDS': '0', 'NATIVE_COMMAND': command,
                }, text=True, capture_output=True)
                self.assertEqual(result.returncode, expected_status, result.stderr)
                fields = (work / 'time.tsv').read_text().strip().split('\t')
                self.assertEqual(fields[0], 'fixture')
                self.assertEqual(len(fields), 6)
                self.assertEqual(int(fields[-1]), expected_status)
                self.assertFalse((work / 'io.tsv').exists())

    def test_timeout_wrapper_status_is_recorded(self):
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            (work / 'output').mkdir()
            (work / 'bin').mkdir()
            # Fake only the outbound timeout command, avoiding clock-dependent waits
            timeout = work / 'bin/timeout'
            timeout.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$@" > "$BENCH_WORK/timeout.args"\nexit 124\n')
            timeout.chmod(0o755)
            script = '''source "$BENCH_COMMON"
benchmark_time_command fixture "$BENCH_WORK/output" "$BENCH_WORK/tool.log" "$BENCH_WORK/tool-error.log" bash -c 'exit 0'
'''
            result = subprocess.run(['bash', '-euo', 'pipefail', '-c', script], env={
                **os.environ, 'BENCH_COMMON': str(BENCHMARKS / 'common.sh'), 'BENCH_WORK': str(work),
                'BENCH_RUN_TIMEOUT_SECONDS': '43', 'PATH': str(work / 'bin') + ':' + os.environ['PATH'],
            }, text=True, capture_output=True)
            self.assertEqual(result.returncode, 124, result.stderr)
            self.assertEqual((work / 'timeout.args').read_text().splitlines()[:3],
                             ['--verbose', '--kill-after=60', '43'])
            self.assertEqual((work / 'time.tsv').read_text().strip().split('\t')[-1], '124')

    def test_recorded_native_failure_preserves_logs_and_allows_the_next_method(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            work, results = root / 'work', root / 'results'
            work.mkdir()
            results.mkdir()
            (work / 'output').mkdir()
            (work / 'time.tsv').write_text('fixture\t1\t0\t0\t10\t3\n')
            (work / 'tool-error.log').write_text('native failure\n')
            (work / 'output/stats.txt').write_text('partial native stdout\n')
            (work / 'stdout.path').write_text(str(work / 'output/stats.txt') + '\n')
            script = '''source "$BENCH_COMMON"
run_order=1
benchmark_finish_run samtools 4 1 true 3
printf 'next method\n' > "$BENCH_RESULTS/continued"
'''
            result = subprocess.run(['bash', '-euo', 'pipefail', '-c', script], env={
                **os.environ, 'BENCH_COMMON': str(BENCHMARKS / 'common.sh'), 'BENCH_WORK': str(work),
                'BENCH_RESULTS': str(results), 'BENCH_PILOT': '0', 'BENCH_WARMUPS': '2',
            }, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((results / 'continued').read_text(), 'next method\n')
            self.assertEqual((results / 'failed-run-1/tool-error.log').read_text(), 'native failure\n')
            self.assertEqual((results / 'failed-run-1/tool-stdout.log').read_text(), 'partial native stdout\n')
            self.assertEqual((results / 'timings.tsv').read_text().strip().split('\t')[-1], '3')
            self.assertFalse((work / 'output').exists())

    def test_missing_timing_record_is_a_fatal_bookkeeping_error(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            work, results = root / 'work', root / 'results'
            work.mkdir()
            results.mkdir()
            script = '''source "$BENCH_COMMON"
run_order=1
benchmark_finish_run samtools 4 1 true 3
printf 'must not continue\n' > "$BENCH_RESULTS/continued"
'''
            result = subprocess.run(['bash', '-euo', 'pipefail', '-c', script], env={
                **os.environ, 'BENCH_COMMON': str(BENCHMARKS / 'common.sh'), 'BENCH_WORK': str(work),
                'BENCH_RESULTS': str(results), 'BENCH_PILOT': '0', 'BENCH_WARMUPS': '2',
            }, text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((results / 'continued').exists())
            self.assertFalse((results / 'timings.tsv').exists())

    def test_lengths_runner_skips_failed_method_repetitions_and_finishes_other_methods(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fake_bin = root / 'bin'
            fake_bin.mkdir()
            scratch = root / 'scratch'
            scratch.mkdir()
            environment_root = root / 'samtools-env'
            (environment_root / 'bin').mkdir(parents=True)
            samtools = environment_root / 'bin/samtools'
            samtools.write_text('#!/usr/bin/env bash\ncase "$1" in idxstats) printf "chr1\\t100\\t0\\t0\\n" ;; --version-only) echo 1.25 ;; *) exit 9 ;; esac\n')
            samtools.chmod(0o755)
            (environment_root / 'bin/conda-unpack').write_text('# Relocation boundary fixture\n')
            archive = root / 'samtools.tar.gz'
            with tarfile.open(archive, 'w:gz') as destination:
                destination.add(environment_root, arcname='.')
            for name, body in (('cfdna', 'echo "cfdna fixture"'), ('ldd', 'echo "statically linked"'),
                               ('scontrol', 'echo "fixture configuration"')):
                executable = fake_bin / name
                executable.write_text('#!/usr/bin/env bash\n' + body + '\n')
                executable.chmod(0o755)
            srun = fake_bin / 'srun'
            srun.write_text(f'#!{sys.executable}\n' + textwrap.dedent(r'''
                import os
                from pathlib import Path
                import sys

                order, repetition, warmup, method, cores = sys.argv[sys.argv.index('--') + 1:]
                work = Path(os.environ['BENCH_WORK'])
                results = Path(os.environ['BENCH_RESULTS'])
                with Path(os.environ['CALL_LOG']).open('a') as destination:
                    destination.write(method + '\n')
                if method == 'cfdnalab':
                    (work / 'tool-error.log').write_text('native boundary failed\n')
                    raise SystemExit(3)
                output = work / 'output'
                output.mkdir()
                profile = output / 'stats.txt'
                profile.write_text('IS\t30\t1\t1\t0\t0\n')
                (work / 'profile.path').write_text(str(profile) + '\n')
                columns = (results / 'timings.tsv').read_text().splitlines()[0].split('\t')
                columns = columns[:columns.index('flush_seconds')]
                values = {'run_order': order, 'repetition': repetition, 'warmup': warmup,
                          'method': method, 'cores': cores, 'job_id': '123', 'step_id': order,
                          'wall_seconds': '1', 'user_seconds': '1', 'system_seconds': '0',
                          'max_rss_kib': '10', 'tool_exit_code': '0'}
                (work / 'time.tsv').write_text('\t'.join(values.get(column, 'NA') for column in columns) + '\n')
                (work / 'io.tsv').write_text('0\t1\t20\n')
            '''))
            srun.chmod(0o755)
            bam = root / 'input.bam'
            bam.write_text('BAM data is handled by the fake samtools boundary\n')
            Path(str(bam) + '.bai').touch()
            results = root / 'results'
            result = subprocess.run(['bash', str(BENCHMARKS / 'lengths/benchmark_lengths.sh')], env={
                **os.environ, 'PATH': str(fake_bin) + ':' + os.environ['PATH'],
                'SLURM_SUBMIT_DIR': str(BENCHMARKS.parent), 'SLURM_JOB_ID': '123',
                'SLURM_JOB_NUM_NODES': '1', 'SLURM_CPUS_PER_TASK': '32',
                'BENCH_BAM': str(bam), 'BENCH_INDEX': str(bam) + '.bai',
                'BENCH_LOCAL_ROOT': str(scratch), 'BENCH_RESULTS': str(results),
                'BENCH_CFDNA': str(fake_bin / 'cfdna'), 'BENCH_SAMTOOLS_ARCHIVE': str(archive),
                'BENCH_CORES': '4', 'BENCH_METHODS': 'cfdnalab samtools',
                'BENCH_WARMUPS': '0', 'BENCH_REPETITIONS': '2', 'BENCH_PILOT': '0',
                'BENCH_SETUP_PYTHON': sys.executable, 'BENCH_EXPECTED_NODE': '',
                'BENCH_RUN_TIMEOUT_SECONDS': '0', 'CALL_LOG': str(root / 'calls'),
            }, text=True, capture_output=True)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertEqual((root / 'calls').read_text().splitlines(), ['cfdnalab', 'samtools', 'samtools'])
            with (results / 'timings.tsv').open() as source:
                rows = list(csv.DictReader(source, delimiter='\t'))
            self.assertEqual([row['step_exit_code'] for row in rows], ['3', '0', '0'])
            for row in rows:
                self.assertNotIn(None, row)
                self.assertNotIn(None, row.values())
            self.assertTrue((results / 'profiles/samtools/stats.txt').exists())
            self.assertEqual((results / 'failed-run-1/tool-error.log').read_text(), 'native boundary failed\n')


class PilotBudgetTests(unittest.TestCase):
    """Keep unavailable runtimes unknown while reporting useful successful budgets."""

    def test_partial_and_all_failed_pilots_report_unknown_budgets(self):
        for include_success in (True, False):
            with self.subTest(include_success=include_success), tempfile.TemporaryDirectory() as temporary:
                path = Path(temporary) / 'pilot.tsv'
                base = dict(run_order='1', repetition='1', warmup='false', pilot='1', method='slow',
                            tool_version='pinned', cores='4', node='node', job_id='123', step_id='0',
                            wall_seconds='NA', user_seconds='NA', system_seconds='NA', flush_seconds='NA',
                            process_plus_flush_seconds='NA', max_rss_kib='NA', output_bytes='NA',
                            tool_exit_code='124', step_exit_code='124')
                rows = [base]
                if include_success:
                    rows.append({**base, 'run_order': '2', 'step_id': '1', 'method': 'cfdnalab',
                                 'wall_seconds': '60', 'user_seconds': '50', 'system_seconds': '10',
                                 'flush_seconds': '0', 'process_plus_flush_seconds': '60',
                                 'max_rss_kib': '10', 'output_bytes': '100',
                                 'tool_exit_code': '0', 'step_exit_code': '0'})
                with path.open('w') as destination:
                    writer = csv.DictWriter(destination, fieldnames=list(base), delimiter='\t')
                    writer.writeheader()
                    writer.writerows(rows)
                result = subprocess.run([sys.executable, str(BENCHMARKS / 'summarize_timings.py'),
                                         str(path), '--estimate'], text=True, capture_output=True)
                self.assertNotEqual(result.returncode, 0)
                report = {row['method']: row for row in csv.DictReader(result.stdout.splitlines(), delimiter='\t')}
                self.assertEqual(report['slow']['budget_minutes'], 'NA')
                self.assertIn('incomplete', result.stderr.lower())
                self.assertNotIn('Traceback', result.stderr)
                if include_success:
                    # Seven 60-second invocations with 25% margin need 8.75 minutes
                    self.assertAlmostEqual(float(report['cfdnalab']['budget_minutes']), 8.75)


@unittest.skipUnless(all(importlib.util.find_spec(name) is not None for name in ('cfdnalab', 'numpy', 'pandas')),
                     'Requires the benchmark comparison environment')
class NativeComparatorTests(unittest.TestCase):
    """Check hand-derived native-format boundaries without running native tools."""

    def comparator(self, relative_path):
        """Load the standalone comparator for a small native-format fixture."""
        specification = importlib.util.spec_from_file_location('review_comparator', BENCHMARKS / relative_path)
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        return module

    def test_samtools_tails_do_not_enter_the_requested_histogram(self):
        import numpy as np

        comparator = self.comparator('lengths/compare_lengths.py')
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / 'samtools'
            directory.mkdir()
            (directory / 'stats.txt').write_text('IS\t29\t999\t999\t0\t0\nIS\t30\t1\t1\t0\t0\n'
                                                'IS\t1000\t2\t2\t0\t0\nIS\t1001\t999\t999\t0\t0\n')
            counts = comparator.read_profile(directory, np.arange(30, 1001))
            self.assertEqual(counts[0], 1)
            self.assertEqual(counts[-1], 2)
            self.assertEqual(counts.sum(), 3)

    def test_valid_picard_histogram_without_fr_is_an_empty_fr_profile(self):
        import numpy as np

        comparator = self.comparator('lengths/compare_lengths.py')
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / 'picard'
            directory.mkdir()
            (directory / 'metrics.txt').write_text('## HISTOGRAM\tjava.lang.Integer\n'
                                                  'insert_size\tAll_Reads.rf_count\n30\t5\n1000\t4\n\n')
            counts = comparator.read_profile(directory, np.arange(30, 1001))
            self.assertEqual(counts.sum(), 0)

    def test_griffin_rounded_density_zero_does_not_erase_a_positive_profile(self):
        import numpy as np

        comparator = self.comparator('midpoints/compare_midpoints.py')
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'griffin.tsv'
            columns = ['site_name', 'number_of_sites', 'smoothing', 'CNA_normalization',
                       'mean_reads_per_bp_in_normalization_window'] + list(map(str, range(-10, 10)))
            with path.open('w') as destination:
                writer = csv.writer(destination, delimiter='\t')
                writer.writerow(columns)
                # Five positions at four and fifteen zeros give mean one
                writer.writerow(['sparse', 2, False, False, 0] + [4] * 5 + [0] * 15)
                writer.writerow(['empty', 2, False, False, 0] + [''] * 20)
            profiles = comparator.read_griffin(path, 20)
            self.assertEqual(profiles['sparse'][1:], (2, 0))
            np.testing.assert_array_equal(profiles['sparse'][0], [4] * 5 + [0] * 15)
            self.assertTrue(np.isnan(profiles['empty'][0]).all())

    def test_midpoint_metrics_keep_zero_positions_and_distinguish_empty_groups(self):
        import numpy as np
        import pandas as pd

        comparator = self.comparator('midpoints/compare_midpoints.py')
        groups = ['reversed', 'both-empty', 'cfdnalab-empty', 'griffin-empty', 'constant']
        counts = np.array([[2] * 10 + [0] * 10, [0] * 20, [0] * 20, [2] * 20, [2] * 20])
        # Fake only the external reader result, exercising the real comparator CLI
        native_result = SimpleNamespace(
            length_bins=lambda: [(100, 201)],
            positions=lambda: pd.DataFrame({'position_bin_start_bp': range(20),
                                             'position_bin_end_bp': range(1, 21)}),
            counts_array=lambda: counts[:, None, :],
            group_metadata=lambda: pd.DataFrame({'group_name': groups, 'eligible_intervals': [2] * 5}),
        )
        with tempfile.TemporaryDirectory() as temporary:
            profiles = Path(temporary) / 'profiles'
            (profiles / 'griffin').mkdir(parents=True)
            columns = ['site_name', 'number_of_sites', 'smoothing', 'CNA_normalization',
                       'mean_reads_per_bp_in_normalization_window'] + list(map(str, range(-10, 10)))
            with (profiles / 'griffin/sample.uncorrected.coverage.tsv').open('w') as destination:
                writer = csv.writer(destination, delimiter='\t')
                writer.writerow(columns)
                writer.writerow(['reversed', 2, False, False, .5] + [0] * 10 + [2] * 10)
                writer.writerow(['both-empty', 2, False, False, 0] + [''] * 20)
                writer.writerow(['cfdnalab-empty', 2, False, False, 1] + [1] * 20)
                writer.writerow(['griffin-empty', 2, False, False, 0] + [''] * 20)
                writer.writerow(['constant', 2, False, False, 1] + [1] * 20)
            output = StringIO()
            with mock.patch.object(comparator, 'read_midpoints', return_value=native_result), \
                    mock.patch.object(sys, 'argv', ['compare_midpoints.py', str(profiles)]), redirect_stdout(output):
                comparator.main()
            rows = {row['group']: row for row in csv.DictReader(output.getvalue().splitlines(), delimiter='\t')}
            # Opposite half-window peaks: 20 absolute differences of two, r=-1
            self.assertAlmostEqual(float(rows['reversed']['pearson_r']), -1)
            self.assertEqual(float(rows['reversed']['mae_mean_one']), 2)
            self.assertEqual(float(rows['reversed']['rmse_mean_one']), 2)
            self.assertEqual(float(rows['reversed']['cfdnalab_mean_midpoints_per_site_bp']), .5)
            for group, status in (('both-empty', 'both-zero-coverage'), ('cfdnalab-empty', 'cfdnalab-zero-coverage'),
                                  ('griffin-empty', 'griffin-zero-coverage')):
                self.assertEqual(rows[group]['status'], status)
                self.assertEqual(rows[group]['mae_mean_one'], 'NA')
                self.assertEqual(rows[group]['pearson_r'], 'NA')
            self.assertEqual(rows['constant']['pearson_r'], 'NA')
            self.assertEqual(float(rows['constant']['mae_mean_one']), 0)


if __name__ == '__main__':
    unittest.main()
