# Benchmarks

These benchmarks compare cfDNAlab with native tools for coverage, fragment lengths, end motifs and midpoint profiles. Production runs measure repeated commands on exclusive Slurm nodes and preserve outputs for numerical comparison. Pilots can share nodes for faster debugging.

Start with installation below and the feature-specific run instructions. Detailed references describe [installation and execution](benchmarks_spec.md) and [methods and concordance](benchmarks_concordance_spec.md).

## Installation

Run `bash benchmarks/install_tools.sh /path/to/benchmark-tools` on Linux x86_64 with Bash 4 or later, using a dedicated directory. Conda environments go under `envs/`. Standalone binaries, source checkouts, packed archives and the container image stay at the installation root. The script prints the paths needed by the benchmark runners. Conda, Apptainer, Git, wget and tar must already be available. Building cfDNAlab also needs GNU make and C/C++ compilers.

The installation root also holds `cache/` for downloaded/extracted Conda packages, pip wheels, Cargo dependencies, Apptainer layers and font caches, plus `tmp/` for unpacking and package builds. Installer configuration and application metadata use `config/` and `share/`. These locations are set before the first Conda call. Selecting an environment with `--prefix` alone does not move Conda's package cache. The installer overrides inherited cache paths for this process, disables environment registration in `~/.conda/environments.txt`, and leaves the submitting shell's settings unchanged. Existing files in earlier cache locations are not moved or deleted.

The benchmark nodes also need Slurm, GNU time at `/usr/bin/time`, GNU coreutils, `ldd`, and Python 3 for untimed setup. These are system requirements for running the jobs and measuring resources.

The installation routes below were checked against the linked guides and package recipes on **2026-10-29**. Software versions are pinned. This is an audit of the instructions and published packages, not a completed installation on the Linux cluster.

| Tool or helper | Installation | Alignment with its guide and reason for this choice |
| --- | --- | --- |
| cfDNAlab | Current checkout, Cargo, `envs/cfdnalab` | Uses the dependencies in the [repository guide](../README.md#installation), including Rust 1.94.0, LLVM/Clang 21 and Fontconfig. `--locked` uses Cargo.lock. `--path` selects the checkout being benchmarked, and `--root` puts the executable inside its environment. |
| mosdepth 0.3.14 | Authors' standalone Linux binary | The [guide](https://github.com/brentp/mosdepth#installation) recommends the release binary and states that it is static. Conda is also supported, but would add an environment without removing any build work. |
| samtools 1.25 | Bioconda, `envs/samtools` | Uses the [Bioconda package](https://anaconda.org/bioconda/samtools), including HTSlib and its shared libraries. No source build or packaging-only Python is needed in this environment. |
| PanDepth 2.26 | Authors' Linux x86_64 release archive | Follows the prebuilt-binary route in the [guide](https://github.com/HuiyangYu/PanDepth#1-install), updating its example's version to 2.26. A [Conda package](https://anaconda.org/bioconda/pandepth) also exists, but this download already avoids compilation. |
| FinaleToolkit `d16632ec` | Upstream Git commit via pip, Conda compiled dependencies, Python 3.11 | Pins merge commit `d16632ec5ae8999e774234457f6b0f9405c8f590` containing [PR #194](https://github.com/epifluidlab/FinaleToolkit/pull/194). Both motif length filtering and midpoint window membership are fixed. The packed environment carries its source commit into timing rows. `pip check`, the Conda export and pip freeze record the dependencies. Reinstall this tool with `bash benchmarks/install_tools.sh ROOT --only finaletoolkit` when replacing the old 1.1.0 archive. |
| cfDNAPro 1.7.4 | Authors' Apptainer image, pinned by digest | Follows the [recommended container installation](https://github.com/hw538/cfDNAPro#installation). The image contains R and the required reference packages. Its labels identify source commit `ab9eeb1` and version 1.7.4, built on October 26. Bioconductor and [Bioconda](https://anaconda.org/bioconda/bioconductor-cfdnapro) use a separate version sequence. We retain the image named in the authors' current guide. |
| Picard 3.5.0 | Bioconda `picard-slim`, Java 17, `envs/picard` | The [recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/picard-slim/meta.yaml) installs the official release JAR. It omits R. [CollectInsertSizeMetrics](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/CollectInsertSizeMetrics.java) only checks for R when a PDF histogram is requested, which our call does not do. |
| FrEIA | Commit `2c32fb9`, authors' Conda pins with `nodefaults`, `envs/freia` | The [guide](https://github.com/mouliere-lab/FrEIA#installation-guide) delegates tool dependencies to its environment files. The derived YAML preserves all Conda dependencies and pins, adds `nodefaults`, and removes the upstream environment name/prefix. Its pip-only `statannot` and `combat` packages are omitted because neither the benchmark's extraction/abundance scripts nor their shared helper imports them. This avoids a PyPI installation step. Both native scripts must successfully print argument help before packaging, checking imports without processing a BAM. Other downstream Conda packages remain. Snakemake, mamba and remapping environments are omitted because we call native scripts on an existing BAM. |
| Griffin v0.2.0 | Release scripts and Conda dependencies, `envs/griffin` | Keeps the relevant library versions in the [wiki](https://github.com/adoebley/Griffin/wiki), including `nomkl`. Python 3.7.12 replaces its unavailable 3.7.4 pin. Bioconda supplies the specified pybedtools 0.8.0 and pysam 0.19.0 without local compilation. PyYAML 5.3.1 is added because the scripts import it. Matplotlib-base 3.4.3 supports native sample GC-bias plots. The archive includes profiling, sample GC estimation and reference GC-generation scripts. The full revised Linux solve remains to be checked. Snakemake and samtools are omitted because these direct calls do not use them. |
| Output comparison | Local `py-cfdnalab`, zstd and pyarrow, `envs/comparison` | Uses the [Python package](../py-cfdnalab/README.md) from this checkout. zstd reads concatenated coverage frames and pyarrow reads FrEIA's native left-end Parquet table. This environment is used after timing for output comparisons. |
| Snyder preparation | Bioconda pysam 0.24.1, Python 3.11, `envs/snyder-preparation` | Follows [pysam's recommended Conda/Bioconda installation](https://pysam.readthedocs.io/en/latest/installation.html), avoiding a possible pip source build. Not part of benchmark timing. |
| Environment packaging | conda-pack 0.9.2, `envs/packaging` | A single helper environment packs the other environments by their paths, following the [conda-pack guide](https://conda.github.io/conda-pack/). It is not a benchmarked tool and is not copied into every tool environment. |

FrEIA has no published release, so its source commit is pinned. Griffin's newest published tag is marked as a prerelease upstream. Its old dependency versions belong to that release's guide, rather than an older Griffin comparison. The scripts copied into the FrEIA and Griffin archives remain unmodified.

Griffin's archive also includes the release's 1000-site CTCF demo for its short warmup. Existing installations need a refreshed Griffin archive and pyarrow in the comparison environment. Use `--only griffin,comparison --force griffin,comparison` to reinstall these components.

The environment settings have these specific purposes:

| Setting | Purpose |
| --- | --- |
| `PYTHONNOUSERSITE=1` | Python normally searches a per-user package directory such as `~/.local/lib/python3.11/site-packages`, even inside a Conda environment. This disables that additional directory while retaining the environment's packages. It prevents a package installed with `pip --user` from overriding the benchmark version. See [Python's definition](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONNOUSERSITE) and [FinaleToolkit's troubleshooting advice](https://epifluidlab.github.io/FinaleToolkit/documentation/user_guide/installation.html#troubleshooting). |
| `unset PYTHONPATH PYTHONHOME` | Removes inherited extra module directories and an inherited Python installation prefix. Both could make the selected Python use files outside its environment. See the [Python environment-variable documentation](https://docs.python.org/3/using/cmdline.html#environment-variables). |
| `CONDA_CHANNEL_PRIORITY=strict` | Applies the channel priority recommended by [Bioconda](https://bioconda.github.io/#with-conda) and FrEIA, including when creating FrEIA's YAML environment. The commands list conda-forge before bioconda and use `--override-channels`. FrEIA uses `nodefaults` in its derived YAML. This only changes this script's environment, without editing the user's `.condarc`. |
| `CONDA_PKGS_DIRS`, `CONDA_ENVS_PATH` | Select `<root>/cache/conda/pkgs` for downloaded and extracted packages and `<root>/envs` for environments. The shared cache beside the existing Conda installation is not used. See [Conda's cache-location guide](https://docs.conda.io/projects/conda/en/stable/user-guide/configuration/custom-env-and-pkg-locations.html). |
| `PIP_CACHE_DIR`, `CARGO_HOME`, `APPTAINER_CACHEDIR` | Keep pip downloads/wheels, Rust dependency downloads and container layers under `<root>/cache`. |
| `TMPDIR`, `TMP`, `TEMP`, `APPTAINER_TMPDIR` | Select `<root>/tmp` for package/build temporary files and container image assembly. The installer's container version checks bind this directory as `/tmp`. |
| `LIBCLANG_PATH="$CONDA_PREFIX/lib"` | Lets cfDNAlab's build find the Conda libclang used to generate Rust bindings. This comes from cfDNAlab's installation guide and is confined to its build section. |
| `JAVA_HOME` | Picard's packaged launcher honors this variable. It selects the installed Java during the version check and the staged Java during benchmarks. The runner also clears inherited JVM option variables so outside settings do not override the recorded call. |

The Python settings appear once in the installer and again in each benchmark job because those are separate shell processes. Running the installer with `bash` does not change the submitting shell's environment.

Packing lets complete Conda environments run from node-local scratch. Copying a Conda executable alone can leave it reading its interpreter or libraries from the original installation. The runners extract the complete archive and run `conda-unpack` before timing. That relocation script needs only Python's standard library. `BENCH_SETUP_PYTHON` (default `python3`) executes it outside timing, including its sibling progress module. Only selected competitors are unpacked. FinaleToolkit's original CLI entrypoint is executed with its packed environment's Python for both version checks and timed runs. This avoids depending on a versioned Python name such as `python3.11` being on the compute node's PATH. Every runner also uses the samtools archive to read the BAM header. Build and pack on the same Linux platform used for benchmarking, as required by conda-pack.

The runners copy the standalone cfDNAlab and PanDepth executables, so their runtime shared libraries must be available from the compute node's system installation. The runners check their Linux linkage with `ldd` before benchmarking. Run a small input on the compute node to verify default QC plotting and fonts as well. The static mosdepth binary has no such package dependencies.

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

Follow the [input preparation](data/README.md), then the [coverage](coverage/README.md), [lengths](lengths/README.md), [end-motif](ends/README.md) or [midpoint](midpoints/README.md) instructions. Submit **from the repository root on the login node** using `bash benchmarks/submit_benchmark.sh FEATURE [sbatch options]`, where `FEATURE` is `coverage`, `lengths`, `ends` or `midpoints`. The launcher requests `BENCH_CORES` CPUs. Full runs reserve an exclusive node and all its memory. With `BENCH_PILOT=1`, the launcher allows shared nodes and requests 50 GB by default, overridable with `--mem`. Node generation, partition and wall time must match the cluster. Jobs use `$TMPDIR` supplied on the compute node, so no scratch path is needed in the submission settings. Each job creates its own working directory below it, directs tool temporary files there and removes only that working directory on exit.

All runners require `BENCH_CORES` to specify a single core count for that job. The experiment uses **1, 4, 12 and 32 in separate jobs**. `BENCH_METHODS` is required and has no default. Supply each method label once. Defaults are `BENCH_WARMUPS=1` and `BENCH_REPETITIONS=3` measured. `BENCH_PILOT=1` selects zero full-input warmups and one measurement at the job's core count. Use a new `BENCH_RESULTS` directory for every job. Only selected tools are required, plus the samtools archive and setup Python. Jobs selecting multiple methods also require `BENCH_COMPARE_PYTHON`, printed by the installer. References are required only by selected methods that use them.

Pilots and full runs automatically write `concordance.tsv` after all tool attempts finish, outside tool timings and step memory measurements. Coverage also writes `residuals.tsv`, containing only the residual histogram. The comparison uses the saved profiles and the job's biological settings. Coverage, length and motif jobs without a comparable pair explain why concordance is skipped. A comparison failure fails the job and preserves its error in `concordance-error.log`, along with the completed timings and profiles. Midpoint jobs run a single method and do not call the comparison. Single-tool jobs in separate allocations still need their profiles combined for cross-job comparisons.

Griffin runs in its own job with no full-input warmups. Before measurements, including a pilot, it runs the release's 1000-site CTCF demo untimed on the same staged BAM/reference and node. This gives Griffin a practical warmup without adding another full production run. The demo has its own Slurm step so its memory is excluded from measured step peaks. Griffin jobs can exceed 24 hours. Choose their wall time from pilots, including all measured repetitions and demo/setup overhead.

The four budgets retain the serial baseline, 4 cores for mosdepth's documented limited scaling, the suggested 12-core configuration, and 32 cores for the upper scaling point. Mosdepth's [guide](https://github.com/brentp/mosdepth/blob/v0.3.14/README.md#speed-and-memory-comparison) describes a plateau near four decompression threads. The runner passes the allocated core count directly, so its 4-core allocation uses `--threads 4`.

Each selected method runs as a contiguous block at the job's core count, with its warmups before measurements. No method is automatically excluded. For cfDNAPro and FinaleToolkit lengths, first compare **1 versus 4 cores in separate jobs**. If wall time and effective CPU use show no benefit beyond run-to-run variation, stop higher-core measurements for that workload and report the tested 1/4-core result. If they improve, add 12- and 32-core jobs. This is an experimental decision, not an automatic scheduler rule.

Allocate N cores and pass N directly to every available native thread/core/worker argument, including options that count additional decoding threads. Reduce that argument only if the tool's documentation explicitly instructs it, with any proposed exception discussed before changing the benchmark. Tools without a native option still run at their selected core settings. The runners do not override numerical-library thread settings by feature.

First check environment relocation, native outputs and fonts with a small input. Then pilot the selected methods at each planned core count in separate jobs on the actual input to budget all invocations. Shared-node pilot timings are rough estimates and are excluded from production summaries. Pilot cfDNAPro's memory progressively, starting near 0.1x, before choosing the common node. Even a single-core-setting pilot can exceed a short allocation when many slow methods are selected. Split methods too when needed.

```bash
python3 benchmarks/summarize_timings.py /path/to/pilot/timings.tsv \
    --estimate --warmups 1 --repetitions 3 > pilot-summary.tsv
```

The estimate writes `budget_minutes`, planned invocation count, failed pilot count and completeness for each observed configuration, adding 25% by default. It budgets four full-input invocations normally and three for Griffin. Failed configurations without a successful pilot have an `NA` budget. It still writes the report, labels any successful subtotal as incomplete and returns nonzero. Add staging, concordance, cleanup and the Griffin demo separately. It cannot budget configurations absent from its inputs. When combining TSVs from several core settings, budget each job separately from the reported matrix subtotal. At 12 minutes per cfDNAlab run and 120 minutes per slow-tool run, those two tools at **one core alone** need 8 h 48 min for a warmup plus three measurements, or 11 hours with margin.

Where allocations are limited to 24 hours, start with a separate job per feature/core budget, then split by method when needed. In minutes, budget `H + 1.25 * 4 * sum(method_runtime_at_this_core_count)`, where `H` covers staging, profile preservation, concordance and cleanup. Four invocations include a warmup and three measurements. With no allowance for `H`, selected methods must total at most 288 minutes per pass to fit 24 hours. Griffin is scheduled separately with a sufficiently long allocation.

The following working scaling scenario starts at **12 minutes per cfDNAlab configuration and 120 minutes per competitor at 1 core**. It assumes linear speedup for parallel paths, with mosdepth's speedup capped at 4x, and adds 25%. The two apparently serial length tools retain 120 minutes per run for their initial 1/4-core budgets. Linear scaling is optimistic for partially parallel or I/O-limited commands, including samtools and Picard. Replace these assumptions with pilot measurements. Times are rounded up to the next minute and exclude `H`.

| Feature/core job | 1 core | 4 cores | 12 cores | 32 cores |
| --- | --- | --- | --- | --- |
| Coverage, 3 cfDNAlab + 5 competitors | 53 h | 13 h 15 min | 7 h 45 min | 6 h 2 min |
| Lengths: cfDNAlab/native, autosomes, samtools, Picard | 22 h | 5 h 30 min | 1 h 50 min | 42 min |
| cfDNAPro lengths, its own job | 10 h | 10 h | Decide after 1/4 | Decide after 1/4 |
| FinaleToolkit lengths, its own job | 10 h | 10 h | Decide after 1/4 | Decide after 1/4 |
| End motifs, 2 cfDNAlab + FinaleToolkit | 12 h | 3 h | 1 h | 23 min |
| Midpoints, cfDNAlab only | 1 h | 15 min | 5 min | 2 min |

Every table cell describes a separate feature/method-block job at exactly that column's core count. The two serial length tools are included as separate 1- and 4-core jobs. Their 12/32-core measurements are deferred pending the comparison. The other lengths block uses `BENCH_METHODS="cfdnalab cfdnalab-autosomes samtools picard"`. The script runs any methods explicitly selected, at any requested single core count.

The 1-core coverage job still needs method splits. The illustrated 1-core lengths block leaves two hours for `H` within a 24-hour allocation. At 1 core, a 120-minute method alone needs 8 hours raw or 10 hours with margin. Two such methods need 16 hours raw or 20 hours with margin. Use the actual pilot sum to choose blocks within each core count. FrEIA and cfDNAPro motifs need their own full-input pilot budgets.

Griffin has one `griffin-gc` job at each core setting, timing sample bias estimation plus profiling. That job retains both raw and corrected profiles. No additional uncorrected Griffin jobs are needed. If a full run takes 30 hours, three measurements need 90 hours raw or 112 h 30 min with margin, plus its short demo and setup overhead. This is a runtime scenario, not a measured estimate. Splitting measured repetitions across jobs remains an experimental choice.

Prefer the same physical node across jobs when available. Otherwise use matched hardware and retain the node identity with every timing. Cross-core speedups then also reflect variation between nodes. Each new job includes its own warmups. To pin a specific node, after setting input/tool paths:

```bash
export BENCH_METHODS=cfdnapro
export BENCH_CORES=1
export BENCH_PILOT=0
export BENCH_EXPECTED_NODE=chosen-node
export BENCH_RESULTS=/path/to/lengths-cfdnapro-1core
bash benchmarks/submit_benchmark.sh lengths \
    --nodelist="$BENCH_EXPECTED_NODE" --time=YOUR_PILOT_BASED_BUDGET
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

Use `slurm_max_rss_kib` as the primary memory result after checking the cluster's accounting plugin, frequency, descendant accounting and file-cache policy once. `slurm-config.txt` records the available configuration. Slurm cgroup v2 MaxRSS can include file cache, so label that policy in results. Memory is collected from Slurm accounting after the job, without an additional sampler in the timed workload. See [Slurm sacct](https://slurm.schedmd.com/sacct.html) and [cgroup v2 limitations](https://slurm.schedmd.com/cgroup_v2.html#limitations).

Effective cores are `(user_seconds + system_seconds) / wall_seconds`. Divide by `cores` for allocated utilization. Use process wall time for this denominator, rather than the flush-inclusive duration.

Inputs and environments are copied to local scratch before warmups. Scheduler startup, profile preservation, cleanup and concordance stay outside measurement. Profiles are copied **after srun** from the first measured repetition at the job's single core count, so copying does not inflate the tool step's accounting peak. FrEIA's intact left-only profile is retained with that end selection identified. No cache eviction is attempted. Staging does not guarantee that the full input remains in page cache.

A recorded native, placement, missing-output or flush failure skips the remaining repetitions of that method and continues other selected methods. Its timing row remains in results, and `runs/<run_order>/` keeps the exact command and the tool's stdout and stderr, as for every run, and the job exits nonzero after finishing those other methods. Setup and results-bookkeeping errors remain fatal. Slurm or kernel termination of the batch shell can prevent final bookkeeping and continuation. Previously saved timing rows and measured profiles remain available when a later run fails or the allocation ends. The active run may not finish saving its measurements or output on a hard termination.

The results directory, timing-table header and `run.tsv` are created before input staging. `progress.log` records initialization, environment extraction, input copying, run starts and completed saves, with the same messages printed to Slurm's output. A header-only `timings.tsv` means no completed run has been recorded yet. Check `progress.log` to distinguish setup from an active run. Coverage and lengths check tool versions before copying the large BAM, so a broken launcher fails before that input transfer.

`BENCH_RUN_TIMEOUT_SECONDS=0` disables per-invocation timeouts by default. Set a positive integer to run the native command under GNU `timeout --verbose --kill-after=60`, inside GNU time. The limit is recorded in `run.tsv`. Leave room for the extra 60-second termination grace, staging, flush and bookkeeping below the allocation limit. Status 124 normally identifies an expired limit. Status 137 alone cannot distinguish a forced timeout kill from another SIGKILL, so retain timeout diagnostics and Slurm accounting. Choose actual limits from pilots. Timeouts are unsuccessful observations and never enter successful production averages.

## Regression tests and validation

Run the automated tests directly in a normal Linux shell. They do not submit Slurm jobs or require an interactive allocation. Slurm commands and benchmark executables are replaced with local fakes. Inputs are small synthetic files created in temporary directories. No production BAM, reference files or benchmark tool installations are needed.

The tests in [test_runners.py](tests/test_runners.py) cover every method's native command and timing settings, the saved profile files, FrEIA's and Picard's helper pipelines, midpoint site preparation, setting checks, the measurement step (placement, signal/timeout statuses, missing outputs) and whole lengths and Griffin jobs with a fake srun, including failure continuation. The tests in [test_benchmarks.py](tests/test_benchmarks.py) cover coverage formats/weights/scope, residual distributions, base-resolution Pearson, fractional regional counts and row identity, installer resumption, Picard status, Slurm joins, reconstructed TLEN signs and native left-motif counts with exclusion of damaged FrEIA aggregates. [Additional regression tests](tests/test_review_regressions.py) cover submission, automatic concordance, partial pilot reports, scope validation before loading tracks, native length boundaries and midpoint empty-profile behavior. Expected numerical values are manually derived from small synthetic inputs.

The full suite needs Python >=3.11, Bash >=4, GNU time at `/usr/bin/time`, GNU coreutils and the `zstd` executable on PATH. The local Python package supplies the comparator dependencies. pyarrow enables the Parquet cases and pysam enables the reconstruction cases. From the repository root, with Conda available:

```bash
conda create -n benchmark-tests -c conda-forge python=3.11 pip zstd pyarrow pysam
conda activate benchmark-tests
python -m pip install ./py-cfdnalab
python -m unittest discover -s benchmarks/tests
```

An existing comparison environment can also run the tests. Missing pyarrow or pysam skips their respective cases. Linux-specific lifecycle tests skip on macOS, so use Linux for the complete suite. The setup above assumes the required Bash and GNU system utilities are already installed.

Real-tool validation is separate and does require a compute-node allocation. After installing the benchmark tools, run each selected native tool and its comparator on a small input to check environment relocation, native output, plotting/fonts and actual Slurm placement/accounting. Then use full-input pilots to establish resource requirements before production measurements. The automated tests do not perform these cluster checks. No measured benchmark results are included in these documents. The result tables are templates.


## User's notes

```bash
source benchmarks/local_settings.sh
BENCH_RESULTS="$benchmark_root/results/coverage-pilot-4cores" BENCH_METHODS="cfdnalab mosdepth" BENCH_CORES=4 BENCH_PILOT=1 bash benchmarks/submit_benchmark.sh coverage
```
