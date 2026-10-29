# Benchmarks

## Installation

Run `bash benchmarks/install_tools.sh /path/to/benchmark-tools` on Linux x86_64, using a new directory. Conda environments go under `envs/`. Standalone binaries, source checkouts, packed archives and the container image stay at the installation root. The script prints the paths needed by the benchmark runners. Conda, Apptainer, Git, wget and tar must already be available. Building cfDNAlab also needs GNU make and C/C++ compilers.

The benchmark nodes also need Slurm, GNU time at `/usr/bin/time`, Linux `ps`, and GNU coreutils. These are system requirements for running the jobs and measuring resources.

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
| FrEIA | Commit `2c32fb9`, authors' unchanged environment YAML, `envs/freia` | The [guide](https://github.com/mouliere-lab/FrEIA#installation-guide) delegates tool dependencies to its environment files. We create `FrEIA_env.yaml` directly and preserve its pins. It includes downstream plotting and analysis packages that the benchmark does not use. Keeping the supplied file avoids maintaining our own dependency selection. Snakemake, mamba and remapping environments are omitted because we call native scripts on an existing BAM. |
| Griffin v0.2.0 | Release scripts and Conda dependencies, `envs/griffin` | Keeps the Python and relevant library versions in the [wiki](https://github.com/adoebley/Griffin/wiki), including its `nomkl` workaround. Bioconda supplies the specified pybedtools 0.8.0 and pysam 0.19.0 without local compilation. PyYAML is added because the scripts import it. Snakemake, matplotlib and samtools are omitted because these direct profiling calls do not use them. |
| Output comparison | Local `py-cfdnalab` plus Conda fastparquet, `envs/comparison` | Uses the [Python package](../py-cfdnalab/README.md) from this checkout. fastparquet reads FrEIA's saved tables. This environment is only used after timing. |
| Snyder preparation | Bioconda pysam 0.24.1, Python 3.11, `envs/snyder-preparation` | Follows [pysam's recommended Conda/Bioconda installation](https://pysam.readthedocs.io/en/latest/installation.html), avoiding a possible pip source build. Not part of benchmark timing. |
| Environment packaging | conda-pack 0.8.1, `envs/packaging` | A single helper environment packs the other environments by their paths, following the [conda-pack guide](https://conda.github.io/conda-pack/). It is not a benchmarked tool and is not copied into every tool environment. |

FrEIA has no published release, so its source commit is pinned. Griffin's newest published tag is marked as a prerelease upstream. Its old dependency versions belong to that release's guide, rather than an older Griffin comparison. The scripts copied into the FrEIA and Griffin archives remain unmodified.

The environment settings have these specific purposes:

| Setting | Purpose |
| --- | --- |
| `PYTHONNOUSERSITE=1` | Python normally searches a per-user package directory such as `~/.local/lib/python3.11/site-packages`, even inside a Conda environment. This disables that additional directory while retaining the environment's packages. It prevents a package installed with `pip --user` from overriding the benchmark version. See [Python's definition](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONNOUSERSITE) and [FinaleToolkit's troubleshooting advice](https://epifluidlab.github.io/FinaleToolkit/documentation/user_guide/installation.html#troubleshooting). |
| `unset PYTHONPATH PYTHONHOME` | Removes inherited extra module directories and an inherited Python installation prefix. Both could make the selected Python use files outside its environment. See the [Python environment-variable documentation](https://docs.python.org/3/using/cmdline.html#environment-variables). |
| `CONDA_CHANNEL_PRIORITY=strict` | Applies the channel priority recommended by [Bioconda](https://bioconda.github.io/#with-conda) and FrEIA, including when creating FrEIA's YAML environment. The commands list conda-forge before bioconda. This only changes this script's environment, without editing the user's `.condarc`. |
| `LIBCLANG_PATH="$CONDA_PREFIX/lib"` | Lets cfDNAlab's build find the Conda libclang used to generate Rust bindings. This comes from cfDNAlab's installation guide and is confined to its build section. |
| `JAVA_HOME` | Picard's packaged launcher honors this variable. It selects the installed Java during the version check and the staged Java during benchmarks. The runner also clears inherited JVM option variables so outside settings do not override the recorded call. |

The Python settings appear once in the installer and again in each benchmark job because those are separate shell processes. Running the installer with `bash` does not change the submitting shell's environment.

Packing is needed for the agreed node-local benchmark setup. Copying a Conda executable alone can leave it reading its interpreter or libraries from the original installation. The runners extract the complete archive and run `conda-unpack` before timing. That relocation script needs only Python's standard library, so FinaleToolkit's staged Python also relocates the samtools and Picard environments. Griffin's Python can run the generated relocation script without installing conda-pack into its old environment. Build and pack on the same Linux platform used for benchmarking, as required by conda-pack.

The runners copy the standalone cfDNAlab and PanDepth executables, so their runtime shared libraries must be available from the compute node's system installation. Check their Linux linkage before benchmarking. It cannot be verified from this macOS workspace. The static mosdepth binary has no such package dependencies.

## Runs and timing TSVs

Install tools with [install_tools.sh](install_tools.sh), then follow the command-specific README or the coverage script's submission example. Submit jobs **from the repository root**, so the runners can copy their helper scripts using `SLURM_SUBMIT_DIR`. Each job reserves an exclusive node and all its memory. Every tool receives the same core budgets, **1, 4, 8, 12 and 32**, through sequential steps within that allocation. Node generation, wall time and partition must match the cluster's settings.

Each timing TSV contains a row per attempted run, including warmups and failures. It is suitable for supplementary data: columns have a header, units are explicit, versions are retained, and there are no input/output paths in the timing rows. Filter `warmup=false`, `tool_exit_code=0` and `step_exit_code=0` for measured successful runs. `run_order` is the actual execution order. `repetition` numbers the warmup and measured blocks together, with three warmups and ten measured blocks by default.

| Columns | Meaning |
| --- | --- |
| `method`, `tool_version` | Native tool/configuration and its recorded version. The length benchmark also identifies cfDNAPro's actual container image by SHA-256. |
| `cores` | Allocated physical core budget, checked against the step's actual affinity. Native worker/thread settings are recorded separately. |
| `node`, `cpu_model`, `cpu_ids`, `physical_core_ids` | Actual hostname, CPU model, allowed Linux logical CPU IDs and unique socket:core IDs. These make placement inspectable. |
| `job_id`, `step_id`, `started_utc` | Scheduler identifiers and UTC start time. |
| `wall_seconds`, `user_seconds`, `system_seconds` | GNU time measurements for the native command or documented complete pipeline. |
| `max_rss_kib` | GNU time's kernel-reported maximum RSS. This does not measure the simultaneous sum of a worker process tree. |
| `sampled_peak_rss_kib` | Largest observed sum of resident memory across the native command and all descendant processes, including pipeline stages and workers. |
| `sampled_mean_rss_kib` | Time-weighted average of those process-tree RSS samples over the command's observed lifetime. |
| `memory_samples`, `memory_interval_seconds` | Number of snapshots and nominal interval, currently one second. |
| `tool_exit_code`, `step_exit_code` | Native timed command status and enclosing scheduler-step status. Missing measurements are `NA`. |

Calculate effective cores as `(user_seconds + system_seconds) / wall_seconds`. Divide by `cores` for allocated CPU utilization. Plotting time against allocated cores shows scaling even for tools without a native thread option. No utilization outcome is assumed in advance.

[measure_memory.py](measure_memory.py) uses Linux `ps` and the Python standard library. It sums process RSS, so shared resident pages can be counted in multiple processes. It excludes the sampler and GNU time itself, and does not include filesystem cache outside those resident mappings. The mean uses elapsed-time weights, including the final partial interval. Peaks between snapshots can be missed, so the sampled peak is explicitly labeled. The separate GNU time maximum remains available. The sampler shares the same core allocation and adds a small monitoring workload whose overhead has not been measured.

Inputs, tool environments and helper scripts are copied to node-local scratch before warmups. Scheduler startup, input staging, output preservation, cleanup and concordance analysis are outside GNU time's timer. Native startup and output writing remain timed. No cache eviction is attempted. Native result files are retained only where needed for output concordance. They are separate from the timing TSV.
