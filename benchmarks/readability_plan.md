# Benchmark readability plan

Status (2026-11-13): implemented. The tests are written but have not been run.

Goal: a reviewer can open a feature folder and see (1) the exact native command for each method, (2) what is timed, (3) which output is kept and how it is compared. The experimental settings, commands, filters and metrics stay as they are. Ideas outside this plan go in `readability_ideas.md`, not into the code.

## What makes the current code hard to read

1. Each of the four runners repeats about 250 lines of the same machinery: setting checks, the run schedule, a `srun ... <<'STEP'` heredoc that runs a second Bash program, a 25-line CPU-topology check, and the result bookkeeping.
2. The heredoc program passes results back through five scratch files (`time.tsv`, `io.tsv`, `profile.path`, `stdout.path`, `step.id`).
3. Timing rows are built positionally: `printf -v prefix '%s\t%s\t...'` with 20-25 fields, plus `fallback[17]=...` arrays whose indexes must match the header by hand.
4. `common.sh` holds tool-specific rules: FrEIA's file rename, samtools/PanDepth gzip handling, cfDNAPro's extra files, Griffin's uncorrected table, autosome chromosome lists, and which coverage methods measure the same quantity (from name fragments such as `*fragment*`).
5. Commands are assembled from `command+=` fragments across `if` branches, so the full command for a method is never written in one place.
6. Comparators partly choose readers from name fragments (`startswith("cfdnalab")`, `"region-counts" in name`, `"fragment" in name`).
7. Tests cut code out of the runners by splitting on comments and heredoc markers, so they break on any edit and do not test what actually runs.

## New structure

```
benchmarks/
  submit_benchmark.sh      unchanged
  common.sh                shared mechanical steps only (checks, staging, run loop helpers, concordance call)
  measure_command.sh       NEW. Runs inside each srun step: checks CPUs, times the command, measures flush and size
  <feature>/commands.sh    NEW. The native command for every method, the output kept, and its timing-table settings
  <feature>/benchmark_*.sh the Slurm job: settings, staging, warmups/measurements, concordance
  ends/freia_motifs.sh                  NEW. FrEIA's three-step workflow, moved out of an inline bash -c string
  lengths/picard_filtered_histogram.sh  NEW. The untimed filtered Picard run, moved out of a heredoc
  midpoints/prepare_sites.py            NEW. The site/Griffin-table preparation, moved out of an 80-line Python heredoc
```

### commands.sh, one per feature

One `case` branch per method, with the **full literal command** (no `command+=` assembly), using the variable names already used in `experiment_commands.md` (`$cfdna`, `$bam`, `$output`, `$cores`, `$min_mapq`, `$min_length`, `$max_length`, `$length_end`, ...). The code and the document can then be compared line by line. Example:

```bash
picard)
    # Picard's documented standalone call with its recommended 2 GB heap. It has no
    # MAPQ option, so the timed run keeps all MAPQs. The filtered histogram used for
    # concordance is made after all measurements (picard_filtered_histogram.sh)
    command=("$picard" -Xmx2g CollectInsertSizeMetrics
        "I=$bam" "O=$output/metrics.txt" "HISTOGRAM_WIDTH=$max_length" "TMP_DIR=$BENCH_WORK/tmp")
    expected_output="$output/metrics.txt"
    row[tool_version]=$picard_version
    row[min_mapq]=NA
    ;;
```

A second function, `save_profile METHOD`, lists exactly which files each method keeps for concordance (for example FrEIA's `T__sample.pq` saved as `left-motifs.pq`). These rules move out of `common.sh`.

### The job script

Reads top to bottom: check settings, stage inputs, measure, compare. The measurement loop is written out in each job script instead of being hidden in `common.sh`:

```bash
for method in "${methods[@]}"; do
    for ((repetition = 1; repetition <= BENCH_WARMUPS + BENCH_REPETITIONS; repetition++)); do
        prepare_command "$method"               # commands.sh
        run_in_slurm_step "$method" "$repetition"   # common.sh: srun measure_command.sh
        save_run_record                         # common.sh: timings.tsv row, command and logs
        if ((step_exit_code == 0 && repetition == BENCH_WARMUPS + 1)); then
            save_profile "$method"              # commands.sh
        fi
        clear_run_output                        # common.sh
        if ((step_exit_code != 0)); then job_exit_code=1; break; fi
    done
done
```

The only feature-specific line inside a loop is FrEIA's worker-count correction in the ends job, placed between the step and the row it corrects.

### Every run stays a separate Slurm step on the same exclusive node (unchanged)

- **One allocation.** A full run is a single-node job with `--exclusive --mem=0`, requested by `submit_benchmark.sh`. The job script checks that it got one node.
- **One `srun` per repetition.** Each warmup and measured repetition is its own `srun` step inside that allocation, with the same options as today:

  ```bash
  srun --jobid="$SLURM_JOB_ID" --exact --nodes=1 --ntasks=1 \
      --cpus-per-task="$cores" --threads-per-core=1 --cpu-bind=threads --export=ALL \
      bash "$BENCH_WORK/measure_command.sh" ... -- "${command[@]}"
  ```

  Slurm binds the step to exactly `cores` physical cores, so the tool cannot use more. Each step keeps its own step ID, which is used to look up its memory in `sacct` afterwards.
- **The CPU check stays inside each step.** `measure_command.sh` checks it before timing anything: the step must be on the job's node, see exactly `cores` logical CPUs, and those CPUs must be `cores` distinct physical cores. A wrong placement fails the run instead of timing it.
- **Untimed work also gets its own step.** The Griffin demo warmup and the filtered Picard histogram stay separate `srun` steps with the same binding, so their memory is not counted in a measured step.
- **What changes is only how the step's code is passed.** Today it is a heredoc. After the refactor it is a script file plus the command as arguments.

### measure_command.sh

Called as `srun ... bash measure_command.sh STDOUT_FILE EXPECTED_OUTPUT EMPTY_OUTPUT_ALLOWED -- "${command[@]}"`. The command arrives as arguments, so the timed command is the exact array that is also written to the run's record. It writes `name<TAB>value` lines (`step_id`, `node`, `cpu_ids`, `wall_seconds`, `tool_exit_code`, `flush_seconds`, ...). A field is simply missing if the step stopped before measuring it, and becomes `NA` in the row.

### Timing rows by name

Each feature lists its timing columns once, in the current order. As implemented, the list is in `commands.sh`, next to the settings that fill it, so a test can check every method's row against it. The row is an associative array (`row[scope]=autosomes`) filled by `commands.sh` and by the measurement, and written column by column. A missing value is an error, not a silent empty field. The column names and order in `timings.tsv` stay as they are, so `summarize_timings.py` is unaffected.

## Small bookkeeping differences (please confirm)

These are not scientific changes, but they change what ends up in a results directory:

1. Every run gets `runs/<run_order>/` with `command.sh` (the exact quoted command), `stdout.log` and `stderr.log`. This replaces `failed-run-<run_order>/`, which only existed for failures.
2. A failed CPU-placement check records `tool_exit_code=NA` instead of `1`, because the tool never ran. `step_exit_code` stays 1.
3. One output check for all features: the expected output must exist, and a file must be non-empty. The exception is samtools depth, where an empty file means zero coverage. Today ends and midpoints only check existence.
4. `chrom.sizes` is saved only in the positional coverage profiles, where the comparator reads it. Today it is also copied into length, motif, midpoint and regional coverage profiles, where nothing reads it.
5. The coverage scratch output directory becomes `$BENCH_WORK/output` like the other features (`$coverage_output` -> `$output` in `experiment_commands.md`).
6. `mosdepth-regions` runs as `env MOSDEPTH_PRECISION=5 mosdepth ...`, so the setting is visible in the recorded command. Today the step exports it before timing. Same effect.
7. FrEIA's two-line `samples.tsv` is written once during setup, outside the output directory, instead of before every repetition. `output_bytes` drops by about 30 bytes.
8. Midpoint steps no longer pre-create empty `tmp/` and `results/` folders. Griffin's helper already creates them inside timing, and cfDNAlab does not use them.
9. Midpoint jobs no longer call the comparator at the end. With one method per job it was always skipped.
10. On failure, the Slurm log shows the tool's stderr and points to `runs/<run_order>/`. Today it also prints stdout, which for samtools stats is the whole report.

## Comparators

Keep the command-line interface, the output files, their columns and values, and the function names the tests import (`read_profile`, `summarize_profile`, `load_profile`, `read_cfdnalab`, `read_griffin`, `Correlation`, `main`). Within that:

- Lengths and coverage: replace the name-fragment checks with explicit tables from method name to reader, and in coverage to compared quantity. An unknown method name is an error. The coverage job's check for a comparable pair is `saved_profiles_share_a_quantity` in `coverage/commands.sh`, which uses the same explicit quantity list as the `workload` column.
- Every reader gets a docstring saying which columns it reads, their units and the coordinate conversion.
- `compare_coverage.py`: split the 180-line `main` into named steps, and put the SQLite residual counter in a small class whose comment explains why it is used (bounded memory). Its behavior stays as it is.
- End motifs already choose the reader from the native file name, and midpoints compare two fixed methods. These only get comment and structure changes.

## Tests

I will not run them, per AGENTS.md. I write the expectations by hand, run only syntax checks (`bash -n`, `py_compile`), and you run the suite.

- Unchanged: comparator, metric, timing-summary, installer, launcher, Snyder reconstruction and midpoint GC-workflow tests. Fixture folders that are not real method names get renamed if the explicit reader tables require it.
- Removed: tests that cut code out of the runners, and tests of the removed `common.sh` functions.
- New `tests/test_runners.py`:
  - `measure_command.sh`: success record, exit and signal statuses, the timeout wrapper, missing or empty output, and wrong placement (Linux, run under `taskset`).
  - Each `commands.sh`: the command and timing settings for every method, by calling `prepare_command` directly.
  - Each `save_profile`: the saved file names and contents.
  - `freia_motifs.sh` and `picard_filtered_histogram.sh` with fake tools (ported from the current tests).
  - `prepare_sites.py`: anchors for both strands and odd widths, groups, chromosome order and demo files.
  - The settings checks and the concordance call: skip, failure and success.
  - A lengths job end to end with fake `srun` and tools: a failing method skips its remaining repetitions, other methods finish, rows and the first measured profile are saved, and the job exits 1.
  - The Griffin demo running as its own step before the measurements.

## Documentation

Update only statements that become false: paths and variable names in `experiment_commands.md` (plus links to each `commands.sh`), the implementation, failure and results sections of `benchmarks_spec.md`, and the failure and test paragraphs of `README.md` and the feature READMEs.

## Differences from Codex's auditability plan

- A single shared `measure_command.sh` that receives the command as arguments, instead of a `run_method.sh` worker per feature. The commands stay in `commands.sh` and are passed in literally, so there is one measurement implementation.
- No new "experiment contract" tables. `experiment_commands.md` already is that document, and the code will use its variable names.
- No before/after traces or metric reproduction from real native outputs. There are none in this checkout, so that belongs to real-tool validation.
- Agreed: residual storage stays in SQLite, and the installer, summary script and data preparation are out of scope.

## Order of work

1. `common.sh`, `measure_command.sh`, and the lengths `commands.sh` and job script.
2. Ends, midpoints and coverage in the same pattern. The old `common.sh` functions go once all four are switched. `benchmark_picard_version` stays because the installer uses it.
3. Comparators.
4. Tests.
5. Documentation and `readability_ideas.md`.
