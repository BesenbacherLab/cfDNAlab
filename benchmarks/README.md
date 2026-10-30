# Benchmarks

Specifications cover [installation and execution](benchmarks_spec.md) and [methods and concordance](benchmarks_concordance_spec.md). The [remaining plan](benchmarks_review_plan.md) lists unfinished implementation and verification. All benchmark plans and specs live in this directory.

The [new review plan](new_review_plan.md) separates implemented repairs from experimental decisions. The [independent review and Claude annotations](codex_review.md) record their evidence and verification limits.

## Installation

Run `bash benchmarks/install_tools.sh /path/to/benchmark-tools` on Linux x86_64 with Bash 4 or later, using a dedicated directory. Conda environments go under `envs/`. Standalone binaries, source checkouts, packed archives and the container image stay at the installation root. The script prints the paths needed by the benchmark runners. Conda, Apptainer, Git, wget and tar must already be available. Building cfDNAlab also needs GNU make and C/C++ compilers.

The benchmark nodes also need Slurm, GNU time at `/usr/bin/time`, GNU coreutils, `ldd`, and Python 3 for untimed setup. These are system requirements for running the jobs and measuring resources.

The installation routes below were checked against the linked guides and package recipes on **2026-10-29**. Software versions are pinned. This is an audit of the instructions and published packages, not a completed installation on the Linux cluster.

| Tool or helper | Installation | Alignment with its guide and reason for this choice |
| --- | --- | --- |
| cfDNAlab | Current checkout, Cargo, `envs/cfdnalab` | Uses the dependencies in the [repository guide](../README.md#installation), including Rust 1.94.0, LLVM/Clang 21 and Fontconfig. `--locked` uses Cargo.lock. `--path` selects the checkout being benchmarked, and `--root` puts the executable inside its environment. |
| mosdepth 0.3.14 | Authors' standalone Linux binary | The [guide](https://github.com/brentp/mosdepth#installation) recommends the release binary and states that it is static. Conda is also supported, but would add an environment without removing any build work. |
| samtools 1.25 | Bioconda, `envs/samtools` | Uses the [Bioconda package](https://anaconda.org/bioconda/samtools), including HTSlib and its shared libraries. No source build or packaging-only Python is needed in this environment. |
| PanDepth 2.26 | Authors' Linux x86_64 release archive | Follows the prebuilt-binary route in the [guide](https://github.com/HuiyangYu/PanDepth#1-install), updating its example's version to 2.26. A [Conda package](https://anaconda.org/bioconda/pandepth) also exists, but this download already avoids compilation. |
| FinaleToolkit 1.1.0 | Bioconda, Python 3.11, `envs/finaletoolkit` | The [authors' guide](https://epifluidlab.github.io/FinaleToolkit/documentation/user_guide/installation.html) uses pip in a clean Python 3.11 environment. We use the same release through its [Bioconda recipe](https://bioconda.github.io/recipes/finaletoolkit/README.html) to obtain compiled dependencies as packages. Extra CLI dependencies for `filter-file` are unnecessary for the benchmarked commands. |
| cfDNAPro 1.7.4 | Authors' Apptainer image, pinned by digest | Follows the [recommended container installation](https://github.com/hw538/cfDNAPro#installation). The image contains R and the required reference packages. Its labels identify source commit `ab9eeb1` and version 1.7.4, built on October 26. Bioconductor and [Bioconda](https://anaconda.org/bioconda/bioconductor-cfdnapro) use a separate version sequence. We retain the image named in the authors' current guide. |
| Picard 3.5.0 | Bioconda `picard-slim`, Java 17, `envs/picard` | The [recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/picard-slim/meta.yaml) installs the official release JAR. It omits R. [CollectInsertSizeMetrics](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/CollectInsertSizeMetrics.java) only checks for R when a PDF histogram is requested, which our call does not do. |
| FrEIA | Commit `2c32fb9`, authors' dependency pins with `nodefaults`, `envs/freia` | The [guide](https://github.com/mouliere-lab/FrEIA#installation-guide) delegates tool dependencies to its environment files. We preserve its dependency pins and add `nodefaults` to a local YAML copy to exclude inherited channels. It includes downstream plotting and analysis packages that the benchmark does not use. Keeping the supplied file avoids maintaining our own dependency selection. Snakemake, mamba and remapping environments are omitted because we call native scripts on an existing BAM. |
| Griffin v0.2.0 | Release scripts and Conda dependencies, `envs/griffin` | Keeps the relevant library versions in the [wiki](https://github.com/adoebley/Griffin/wiki), including `nomkl`. Python 3.7.12 replaces its unavailable 3.7.4 pin. Bioconda supplies the specified pybedtools 0.8.0 and pysam 0.19.0 without local compilation. PyYAML 5.3.1 is added because the scripts import it. The full revised Linux solve remains to be checked. Snakemake, matplotlib and samtools are omitted because these direct profiling calls do not use them. |
| Output comparison | Local `py-cfdnalab` plus zstd, `envs/comparison` | Uses the [Python package](../py-cfdnalab/README.md) from this checkout. zstd reads concatenated coverage frames. FrEIA concordance is disabled. This environment is only used after timing. |
| Snyder preparation | Bioconda pysam 0.24.1, Python 3.11, `envs/snyder-preparation` | Follows [pysam's recommended Conda/Bioconda installation](https://pysam.readthedocs.io/en/latest/installation.html), avoiding a possible pip source build. Not part of benchmark timing. |
| Environment packaging | conda-pack 0.9.2, `envs/packaging` | A single helper environment packs the other environments by their paths, following the [conda-pack guide](https://conda.github.io/conda-pack/). It is not a benchmarked tool and is not copied into every tool environment. |

FrEIA has no published release, so its source commit is pinned. Griffin's newest published tag is marked as a prerelease upstream. Its old dependency versions belong to that release's guide, rather than an older Griffin comparison. The scripts copied into the FrEIA and Griffin archives remain unmodified.

The environment settings have these specific purposes:

| Setting | Purpose |
| --- | --- |
| `PYTHONNOUSERSITE=1` | Python normally searches a per-user package directory such as `~/.local/lib/python3.11/site-packages`, even inside a Conda environment. This disables that additional directory while retaining the environment's packages. It prevents a package installed with `pip --user` from overriding the benchmark version. See [Python's definition](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONNOUSERSITE) and [FinaleToolkit's troubleshooting advice](https://epifluidlab.github.io/FinaleToolkit/documentation/user_guide/installation.html#troubleshooting). |
| `unset PYTHONPATH PYTHONHOME` | Removes inherited extra module directories and an inherited Python installation prefix. Both could make the selected Python use files outside its environment. See the [Python environment-variable documentation](https://docs.python.org/3/using/cmdline.html#environment-variables). |
| `CONDA_CHANNEL_PRIORITY=strict` | Applies the channel priority recommended by [Bioconda](https://bioconda.github.io/#with-conda) and FrEIA, including when creating FrEIA's YAML environment. The commands list conda-forge before bioconda and use `--override-channels`. FrEIA uses `nodefaults` in its derived YAML. This only changes this script's environment, without editing the user's `.condarc`. |
| `LIBCLANG_PATH="$CONDA_PREFIX/lib"` | Lets cfDNAlab's build find the Conda libclang used to generate Rust bindings. This comes from cfDNAlab's installation guide and is confined to its build section. |
| `JAVA_HOME` | Picard's packaged launcher honors this variable. It selects the installed Java during the version check and the staged Java during benchmarks. The runner also clears inherited JVM option variables so outside settings do not override the recorded call. |

The Python settings appear once in the installer and again in each benchmark job because those are separate shell processes. Running the installer with `bash` does not change the submitting shell's environment.

Packing is needed for the agreed node-local benchmark setup. Copying a Conda executable alone can leave it reading its interpreter or libraries from the original installation. The runners extract the complete archive and run `conda-unpack` before timing. That relocation script needs only Python's standard library. `BENCH_SETUP_PYTHON` (default `python3`) executes it outside timing, including its sibling progress module. Only selected competitors are unpacked. Every runner also uses the samtools archive to read the BAM header. Build and pack on the same Linux platform used for benchmarking, as required by conda-pack.

The runners copy the standalone cfDNAlab and PanDepth executables, so their runtime shared libraries must be available from the compute node's system installation. The runners check their Linux linkage with `ldd` before benchmarking. A small pilot must also exercise default QC plotting and fonts. It cannot be verified from this macOS workspace. The static mosdepth binary has no such package dependencies.

## Resuming installations

The installer records `tool`, `version`, `status` and expected artifacts in `manifest.tsv`. A matching completed version with existing artifacts skips. A changed version, missing artifact or outdated status reinstalls that tool, with a separate log in `logs/`. Cleanup is limited to that tool's known files. cfDNAlab keeps an incremental Cargo target directory under `build/cfdnalab/`.

```bash
bash benchmarks/install_tools.sh /path/to/benchmark-tools --only cfdnalab,samtools
bash benchmarks/install_tools.sh /path/to/benchmark-tools --force finaletoolkit
bash benchmarks/install_tools.sh /path/to/benchmark-tools --mark-outdated cfdnalab
```

Use `--force all` for a full reinstall. After source or dependency changes without a version change, explicitly force or mark the relevant tool outdated. `--force` applies to the tools selected by `--only`, which otherwise defaults to all. There is no automatic source-change detection. Environment exports, pip records and the image identity stay under `records/`. `cfdnalab-build.tsv` records the installed checkout and dirty flag, rather than a later checkout state. The installer prints `BENCH_INSTALL_MANIFEST` and `BENCH_CFDNA_BUILD` so runners can copy these small records beside results. Retain the full installation records with the experiment.

Picard's version probe accepts its documented status 1 only with a recognizable 3.5.0 version. FinaleToolkit's added runtime pins and conda-pack 0.9.2 need Linux solve/relocation verification. Availability and source inspection are not evidence of a completed installation.

## Runs and budgets

Follow the [input preparation](data/README.md), then the [coverage](coverage/README.md), [lengths](lengths/README.md), [end-motif](ends/README.md) or [midpoint](midpoints/README.md) instructions. Submit **from the repository root** so `SLURM_SUBMIT_DIR` locates helpers. Jobs reserve an exclusive node and all its memory. Node generation, partition and wall time must match the cluster.

All runners require `BENCH_CORES` to specify a single core count for that job. The experiment uses **1, 4, 12 and 32 in separate jobs**. `BENCH_METHODS` is required and has no default. Supply each method label once. Set `BENCH_WARMUPS` (2) and `BENCH_REPETITIONS` (5 measured). `BENCH_PILOT=1` selects zero warmups and one measurement at the job's core count. Use a new `BENCH_RESULTS` directory for every job. Only selected tools are required, plus the samtools archive and setup Python. References are required only by selected methods that use them.

The four budgets retain the serial baseline, 4 cores for mosdepth's documented limited scaling, the suggested 12-core configuration, and 32 cores for the upper scaling point. Mosdepth's [guide](https://github.com/brentp/mosdepth/blob/v0.3.14/README.md#speed-and-memory-comparison) describes a plateau near four decompression threads. The runner passes `cores - 1` additional decompression threads, so its 4-core allocation uses three.

Each selected method runs as a contiguous block at the job's core count, with its warmups before measurements. No method is automatically excluded. For cfDNAPro and FinaleToolkit lengths, first compare **1 versus 4 cores in separate jobs**. If wall time and effective CPU use show no benefit beyond run-to-run variation, stop higher-core measurements for that workload and report the tested 1/4-core result. If they improve, add 12- and 32-core jobs. This is an experimental decision, not an automatic scheduler rule.

First check wiring/relocation/fonts with a small input. Then pilot the selected methods at each planned core count in separate jobs on the actual input to budget all invocations. Pilot cfDNAPro's memory progressively, starting near 0.1x, before choosing the common node. Even a single-core-setting pilot can exceed a short allocation when many slow methods are selected. Split methods too when needed.

```bash
python3 benchmarks/summarize_timings.py /path/to/pilot/timings.tsv \
    --estimate --warmups 2 --repetitions 5 > pilot-summary.tsv
```

The estimate writes `budget_minutes`, planned invocation count, failed pilot count and completeness for each observed configuration, adding 25% by default. Failed configurations without a successful pilot have an `NA` budget. It still writes the report, labels any successful subtotal as incomplete and returns nonzero. Add staging/cleanup separately. It cannot budget configurations absent from its inputs. Combining TSVs reports a matrix subtotal, not permission to combine core settings in a job. At 12 minutes per cfDNAlab run and 120 minutes per slow-tool run, those two tools at **one core alone** need 15 h 24 min for two warmups plus five measurements.

For a 24-hour allocation limit, start with a separate job per feature/core budget, then split by method when needed. In minutes, budget `H + 1.25 * 7 * sum(method_runtime_at_this_core_count)`, where `H` covers staging, profile preservation and cleanup. Seven invocations include two warmups and five measurements. With no allowance for `H`, selected methods must total at most about 164.6 minutes per pass to fit 24 hours.

The following working scaling scenario starts at **12 minutes per cfDNAlab configuration and 120 minutes per competitor at 1 core**. It assumes linear speedup for parallel paths, with mosdepth's speedup capped at 4x, and adds 25%. The two apparently serial length tools retain 120 minutes per run for their initial 1/4-core budgets. Linear scaling is optimistic for partially parallel or I/O-limited commands, including samtools and Picard. Replace these assumptions with pilot measurements. Times are rounded up to the next minute and exclude `H`.

| Feature/core job | 1 core | 4 cores | 12 cores | 32 cores |
| --- | --- | --- | --- | --- |
| Coverage, 3 cfDNAlab + 5 competitors | 92 h 45 min | 23 h 12 min | 13 h 34 min | 10 h 34 min |
| Lengths: cfDNAlab/native, autosomes, samtools, Picard | 38 h 30 min | 9 h 38 min | 3 h 13 min | 1 h 13 min |
| cfDNAPro lengths, its own job | 17 h 30 min | 17 h 30 min | Decide after 1/4 | Decide after 1/4 |
| FinaleToolkit lengths, its own job | 17 h 30 min | 17 h 30 min | Decide after 1/4 | Decide after 1/4 |
| End motifs, 2 cfDNAlab + 1 competitor | 21 h | 5 h 15 min | 1 h 45 min | 40 min |
| Midpoints, 1 cfDNAlab + 1 competitor | 19 h 15 min | 4 h 49 min | 1 h 37 min | 37 min |

Every table cell describes a separate feature/method-block job at exactly that column's core count. The two serial length tools are included as separate 1- and 4-core jobs. Their 12/32-core measurements are deferred pending the comparison. The other lengths block uses `BENCH_METHODS="cfdnalab cfdnalab-autosomes samtools picard"`. The script runs any methods explicitly selected, at any requested single core count.

The 1-core coverage and other-lengths jobs still need method splits. The 4-core coverage estimate leaves less than 49 minutes for `H`, so split it further if staging and preservation need longer. At 1 core, a 120-minute method alone needs 14 hours raw or 17 h 30 min with margin. Pairing it with a 12-minute cfDNAlab configuration needs 15 h 24 min raw or 19 h 15 min with margin. Two 120-minute methods already need 28 hours raw. Use the actual pilot sum to choose blocks within each core count.

Prefer the same physical node across jobs when available. Otherwise use matched hardware and retain the node identity with every timing. Cross-core speedups then also reflect variation between nodes. Each new job includes its own warmups. To pin a specific node, after setting input/tool paths:

```bash
export BENCH_METHODS=cfdnapro
export BENCH_CORES=1
export BENCH_EXPECTED_NODE=chosen-node
export BENCH_RESULTS=/path/to/lengths-cfdnapro-1core
sbatch --nodelist="$BENCH_EXPECTED_NODE" --time=YOUR_PILOT_BASED_BUDGET \
    benchmarks/lengths/benchmark_lengths.sh
```

## Timing and Slurm memory

Each TSV keeps warmups, pilots and failures. Summaries exclude them from production results. Actual run order, version, filters, scope, hostname, logical/physical CPU IDs, job/step IDs and UTC start remain inspectable. Native worker counts are distinct from the allocated physical core budget. Placement is checked inside every single-task step. Hard confinement also depends on the site's Slurm cgroup configuration.

| Columns | Meaning |
| --- | --- |
| `wall_seconds`, `user_seconds`, `system_seconds` | GNU time for the native command or documented pipeline, including native output generation and compression. |
| `flush_seconds` | A separately timed `sync -f` on the output filesystem immediately after command completion. This can include other filesystem writes. |
| `process_plus_flush_seconds` | Process wall time plus filesystem flush time, the main proposed end-to-end result. It is not exact per-tool durable-write latency. |
| `output_bytes` | `du -sb` of the native output directory, including required intermediates. |
| `max_rss_kib` | GNU time's maximum RSS, kept separately from Slurm step memory. |
| `tool_exit_code`, `step_exit_code`, `pilot` | Native/step status and whether the row belongs to a pilot. Unavailable measurements are `NA`. |

Collect completed Slurm accounting **once after the job** and keep the raw export:

```bash
sacct -j JOB_ID --parsable2 --units=K \
    --format=JobIDRaw,NTasks,MaxRSS,State,ExitCode > accounting.tsv
python3 benchmarks/summarize_timings.py /path/to/results/timings.tsv \
    --accounting accounting.tsv --runs-output runs-with-memory.tsv > summary.tsv
```

Supply multiple timing TSVs only for blocks of the same feature/experiment with the same biological input. For multiple jobs, collect their IDs together using `sacct -j ID1,ID2`. The summary joins exact job/step IDs and reports sample size, mean/SD, median and quartiles per configuration. `runs-with-memory.tsv` retains every measured point and failure for plotting. Missing/delayed accounting stays `NA`, so recollect and repeat the join later.

Use `slurm_max_rss_kib` as the primary memory result after checking the cluster's accounting plugin, frequency, descendant accounting and file-cache policy once. `slurm-config.txt` records the available configuration. Slurm cgroup v2 MaxRSS can include file cache, so label that policy in results. The removed process-tree sampler is no longer part of the workload. See [Slurm sacct](https://slurm.schedmd.com/sacct.html) and [cgroup v2 limitations](https://slurm.schedmd.com/cgroup_v2.html#limitations).

Effective cores are `(user_seconds + system_seconds) / wall_seconds`. Divide by `cores` for allocated utilization. Use process wall time for this denominator, rather than the flush-inclusive duration.

Inputs and environments are copied to local scratch before warmups. Scheduler startup, profile preservation, cleanup and concordance stay outside measurement. Profiles are copied **after srun** from the first measured repetition at the job's single core count, so copying does not inflate the tool step's accounting peak. FrEIA profiles are excluded. No cache eviction is attempted. Staging does not guarantee that the full input remains in page cache.

A recorded native, placement, missing-output or flush failure skips the remaining repetitions of that method and continues other selected methods. Its timing row and available tool logs remain in results, and the job exits nonzero after finishing those other methods. Setup and results-bookkeeping errors remain fatal. Slurm or kernel termination of the batch shell can prevent final bookkeeping and continuation.

`BENCH_RUN_TIMEOUT_SECONDS=0` disables per-invocation timeouts by default. Set a positive integer to run the native command under GNU `timeout --verbose --kill-after=60`, inside GNU time. The limit is recorded in `run.tsv`. Leave room for the extra 60-second termination grace, staging, flush and bookkeeping below the allocation limit. Status 124 normally identifies an expired limit. Status 137 alone cannot distinguish a forced timeout kill from another SIGKILL, so retain timeout diagnostics and Slurm accounting. Choose actual limits from pilots. Timeouts are unsuccessful observations and never enter successful production averages.

## Prepared regression checks

The hand-derived fixtures in [tests/test_benchmarks.py](tests/test_benchmarks.py) cover coverage formats/weights/scope, residual distributions, base-resolution Pearson, regional row identity, scheduling, installer resumption, Picard status, Slurm joins, reconstructed TLEN signs and FrEIA exclusion. [Review regressions](tests/test_review_regressions.py) add failure continuation, signal/timeout statuses, partial pilot reports, scope preflight, native length boundaries and midpoint empty-profile behavior. They were not run during implementation. On Linux, use the comparison Python with zstd on PATH to run `python -m unittest discover -s benchmarks/tests`. The reconstruction cases require pysam and explicitly skip when it is unavailable. Linux environment solves, relocation, fonts and real Slurm accounting still need the small wiring pilot.
