# Benchmarks

These benchmarks compare cfDNAlab with native tools for coverage, fragment lengths, end motifs and midpoint profiles. Production runs measure repeated commands on exclusive Slurm nodes and preserve outputs for numerical comparison.

Start with installation below and the commented [production run procedure](#runs-and-budgets). Detailed references describe [installation and execution](benchmarks_spec.md) and [methods and concordance](benchmarks_concordance_spec.md).

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
| FinaleToolkit `d16632ec` | Upstream Git commit via pip, Conda compiled dependencies, Python 3.11 | Pins merge commit `d16632ec5ae8999e774234457f6b0f9405c8f590` containing [PR #194](https://github.com/epifluidlab/FinaleToolkit/pull/194). That PR fixes the 1.1.0 release's motif length filtering and midpoint window membership. The packed environment carries its source commit into timing rows. `pip check`, the Conda export and pip freeze record the dependencies. |
| cfDNAPro 1.7.4 | Authors' Apptainer image, pinned by digest | Follows the [recommended container installation](https://github.com/hw538/cfDNAPro#installation). The image contains R and the required reference packages. Its labels identify source commit `ab9eeb1` and version 1.7.4, built on October 26. Bioconductor and [Bioconda](https://anaconda.org/bioconda/bioconductor-cfdnapro) use a separate version sequence. We retain the image named in the authors' current guide. |
| Picard 3.5.0 | Bioconda `picard-slim`, Java 17, `envs/picard` | The [recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/picard-slim/meta.yaml) installs the official release JAR. It omits R. [CollectInsertSizeMetrics](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/CollectInsertSizeMetrics.java) only checks for R when a PDF histogram is requested, which our call does not do. |
| Griffin v0.2.0 | Release scripts and Conda dependencies, `envs/griffin` | Keeps the relevant library versions in the [wiki](https://github.com/adoebley/Griffin/wiki), including `nomkl`. Python 3.7.12 replaces its unavailable 3.7.4 pin. Bioconda supplies the specified pybedtools 0.8.0 and pysam 0.19.0 without local compilation. PyYAML 5.3.1 is added because the scripts import it. Matplotlib-base 3.4.3 supports native sample GC-bias plots. The archive includes profiling, sample GC estimation and reference GC-generation scripts. The full revised Linux solve remains to be checked. Snakemake and samtools are omitted because these direct calls do not use them. |
| Output comparison | Local `py-cfdnalab` and zstd, `envs/comparison` | Uses the [Python package](../py-cfdnalab/README.md) from this checkout. zstd reads concatenated coverage frames. This environment is used after timing for output comparisons. |
| Snyder preparation | Bioconda pysam 0.24.1, Python 3.11, `envs/snyder-preparation` | Follows [pysam's recommended Conda/Bioconda installation](https://pysam.readthedocs.io/en/latest/installation.html), avoiding a possible pip source build. Not part of benchmark timing. |
| Environment packaging | conda-pack 0.9.2, `envs/packaging` | A single helper environment packs the other environments by their paths, following the [conda-pack guide](https://conda.github.io/conda-pack/). It is not a benchmarked tool and is not copied into every tool environment. |

Griffin's newest published tag is marked as a prerelease upstream. Its old dependency versions belong to that release's guide, rather than an older Griffin comparison. The scripts copied into the Griffin archive remain unmodified.

Griffin's archive also includes the release's 1000-site CTCF demo for its short warmup.

The environment settings have these specific purposes:

| Setting | Purpose |
| --- | --- |
| `PYTHONNOUSERSITE=1` | Python normally searches a per-user package directory such as `~/.local/lib/python3.11/site-packages`, even inside a Conda environment. This disables that additional directory while retaining the environment's packages. It prevents a package installed with `pip --user` from overriding the benchmark version. See [Python's definition](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONNOUSERSITE) and [FinaleToolkit's troubleshooting advice](https://epifluidlab.github.io/FinaleToolkit/documentation/user_guide/installation.html#troubleshooting). |
| `unset PYTHONPATH PYTHONHOME` | Removes inherited extra module directories and an inherited Python installation prefix. Both could make the selected Python use files outside its environment. See the [Python environment-variable documentation](https://docs.python.org/3/using/cmdline.html#environment-variables). |
| `CONDA_CHANNEL_PRIORITY=strict` | Applies the channel priority recommended by [Bioconda](https://bioconda.github.io/#with-conda). The commands list conda-forge before bioconda and use `--override-channels`. This only changes this script's environment, without editing the user's `.condarc`. |
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

This is the production run procedure for all four features. The feature READMEs explain native counting differences. The executable method definitions are in [coverage/commands.sh](coverage/commands.sh), [lengths/commands.sh](lengths/commands.sh), [ends/commands.sh](ends/commands.sh) and [midpoints/commands.sh](midpoints/commands.sh). [experiment_commands.md](experiment_commands.md) expands those native commands for inspection.

### Production matrix

Submit the normal timing rows below at **1, 4, 12 and 32 cores**, giving **24 exclusive-node timing jobs**. Methods within a job run sequentially. Each method receives a full-input warmup and three measured repetitions, except Griffin, which receives its short untimed release demo followed by three full-input measurements. Every native thread/worker argument receives the allocated core count directly. Commands without a native thread option still run at all four core settings, with their actual CPU use recorded. The proper-pair diagnostics run once each at 32 cores in three separate concordance jobs, giving **27 submissions in total**. Their rows are outside the normal timing summaries.

| Feature / job | Methods per job | Filters | Reason |
| --- | --- | --- | --- |
| Coverage | 13 | MAPQ >=30, inclusive lengths 30-1000 bp | Compare read depth, fragment coverage, regional read depth and regional fragment counts separately, timing the normal cfDNAlab pairing policy. The autosome-only read-depth run also measures cfDNAlab's usual chromosome scope. |
| Lengths | 5 | MAPQ >=30, inclusive lengths 30-1000 bp | Compare normal cfDNAlab directional spans against curated cfDNAPro and the TLEN-based tools. |
| End motifs | 3 | MAPQ >=30, inclusive lengths 30-1000 bp, 4-mers | Compare normal cfDNAlab reference-derived motifs against FinaleToolkit and native cfDNAPro left ends. Left-versus-both comparisons retain their approximation label. |
| Midpoints / cfDNAlab uncorrected | 1 | MAPQ >=20, inclusive profiling lengths 100-200 bp | Compare raw cfDNAlab profiles with the raw table already produced by the Griffin job. |
| Midpoints / cfDNAlab GC-corrected | 1 | MAPQ >=20, inclusive profiling lengths 100-200 bp | Time cfDNAlab sample GC estimation plus corrected profiling. |
| Midpoints / Griffin GC-corrected | 1 | MAPQ >=20, inclusive profiling lengths 100-200 bp | Time Griffin sample GC estimation plus profiling. This job retains both raw and corrected profiles for the two midpoint comparisons. |

Each proper-pair diagnostic adds flag 0x2 selection to the corresponding normal cfDNAlab configuration. Each runs once without a warmup at 32 cores, with the same MAPQ and length bounds as the normal comparison. The two coverage diagnostics share a job. Lengths and end motifs each have a diagnostic job. These checks do not replace the normal cfDNAlab timing runs.

| 32-core diagnostic | Comparator |
| --- | --- |
| Coverage `cfdnalab-fragment-proper` | `mosdepth-fragment` |
| Coverage `cfdnalab-region-counts-proper` | `finaletoolkit-coverage` |
| Lengths `cfdnalab-proper` | `finaletoolkit` |
| End motifs `cfdnalab-reference-proper` | `finaletoolkit` |

The main coverage, length and motif cfDNAlab calls explicitly use all BAM contigs. The additional `cfdnalab-autosomes` coverage run measures the usual cfDNAlab chromosome scope and the effect of processing fewer chromosomes. The motif series uses 4-mers. Midpoint profiles use the same prepared 9,990 bp grouped site windows, 15 bp bins and native 165 bp smoothing.

The tools apply these filters through different native selection rules. PanDepth has no fragment length filter. Timed Picard keeps its native MAPQ/QC selection, and the runner makes a separately filtered histogram for concordance after timing. cfDNAPro motif lengths are subset in R inside timing, and its contigs are limited to those present in its bundled reference. Midpoint sample GC estimation retains the tools' native estimation ranges: cfDNAlab inherits the precomputed reference package's range, normally 30-1000 bp, while Griffin uses 15-500 bp. Both midpoint profilers use 100-200 bp. These differences are recorded and explained in the feature READMEs.

Each successful cfDNAPro motif invocation also retains a [stage timing breakdown](ends/README.md#cfdnapro-motifs) at `runs/<run_order>/cfdnapro-stage-timings.tsv`. It separates native BAM reading/curation and motif counting from the benchmark's R length filter, while the complete process remains the end-to-end timing result.

### Local paths and preparation

Run the commands in **Bash 4 or later from the repository root on the login node**, after [input preparation](data/README.md) and installation. Keep absolute paths in the untracked `benchmarks/local_settings.sh`. That file sets `benchmark_root` and exports the following input and tool settings. It should only set variables, without downloads, reference calculations or submissions.

| Settings | Required contents |
| --- | --- |
| `benchmark_root` | Parent directory for the new results series |
| `BENCH_BAM` | Prepared coordinate-sorted hg38 BAM. `BENCH_INDEX` can override its default `<BAM>.bai` path. Coverage, lengths and motifs require BAI |
| `BENCH_CFDNA`, `BENCH_MOSDEPTH`, `BENCH_PANDEPTH` | Executable paths printed by the installer |
| `BENCH_SAMTOOLS_ARCHIVE`, `BENCH_FINALE_ARCHIVE`, `BENCH_PICARD_ARCHIVE`, `BENCH_GRIFFIN_ARCHIVE`, `BENCH_CFDNAPRO_IMAGE`, `BENCH_COMPARE_PYTHON` | Installed archive, image and comparison Python paths printed by the installer |
| `BENCH_REF_2BIT`, `BENCH_REF_FASTA` | Matching alignment reference in 2bit and uncompressed FASTA formats. Griffin requires an existing `<FASTA>.fai`, overridable with `BENCH_FASTA_INDEX` |
| `BENCH_INTERVALS` | Prepared grouped 9,990 bp site BED described in the [midpoint input instructions](midpoints/README.md#inputs-and-run) |
| `BENCH_CFDNA_REF_GC` | Reusable cfDNAlab reference GC Zarr package. Prepare it once as described in the [GC reference instructions](midpoints/README.md#gc-corrected-comparison) |
| `BENCH_GRIFFIN_REF_DIR`, `BENCH_GRIFFIN_GC_REGIONS`, `BENCH_GRIFFIN_REF_GC`, `BENCH_GRIFFIN_MAPPABILITY_BW` | Griffin exclusion BED directory, bundled estimation BED/reference GC tables and the separately downloaded hg38 k100 Umap bigWig, described in the [midpoint instructions](midpoints/README.md) |

`BENCH_SETUP_PYTHON` defaults to `python3` and runs untimed setup only. Optional `BENCH_INSTALL_MANIFEST` and `BENCH_CFDNA_BUILD` preserve installer/build records with each job. `$TMPDIR` is supplied by the compute node, so the settings file does not need a scratch path.

### Submit the production jobs

The six wall-time variables below name **complete jobs**, including setup, every selected method's repetitions and any concordance work. The three midpoint names explicitly include their feature and correction workflow. They do not set cfDNAlab's budget in other features. The configured limits below are allocation limits rather than measured runtime estimates. A single value per timing job type must accommodate that type's slowest core setting. The 32-core diagnostic jobs reuse their feature's limit but execute each diagnostic only once. Griffin's limit must cover three complete estimation-plus-profiling runs and its demo/setup overhead. If a complete timing job exceeds the partition limit, split its method list into named blocks with new result directories and retain all selected method/core combinations.

Adjust `--constraint=gen6` and add a partition in `slurm_options` if required by the cluster. Prefer the same physical node when available. Setting `BENCH_EXPECTED_NODE` in the local settings both adds `--nodelist` below and enables the runner's hostname check. Otherwise use matched hardware and retain the recorded node identity when interpreting scaling.

```bash
set -euo pipefail
source benchmarks/local_settings.sh
: "${benchmark_root:?Set benchmark_root in benchmarks/local_settings.sh}"

# Limits for the complete multi-method jobs, including cfDNAlab and competitors
coverage_wall_time="11:00:00"
lengths_wall_time="11:00:00"
end_motifs_wall_time="11:00:00"

# Limits for the three separate midpoint jobs
midpoints_uncorrected_wall_time="11:00:00"
midpoints_cfdnalab_gc_wall_time="11:00:00"
midpoints_griffin_gc_wall_time="1-23:00:00"

# Use exclusive production allocations with a warmup and three measurements
# Griffin replaces the full-input warmup with its short untimed demo internally
export BENCH_PILOT=0
export BENCH_WARMUPS=1 BENCH_REPETITIONS=3
# Preserve an invocation timeout from local settings, defaulting to disabled
export BENCH_RUN_TIMEOUT_SECONDS=${BENCH_RUN_TIMEOUT_SECONDS:-0}
export BENCH_MIN_MAPQ=30 BENCH_MIN_LENGTH=30 BENCH_MAX_LENGTH=1000
export BENCH_K=4 BENCH_GENOME=hg38

# Create a new series and feature parents, leaving each job's directory uncreated
# The runners reject existing BENCH_RESULTS directories and create their own
mkdir -p "$benchmark_root/results"
results_root=$(mktemp -d "$benchmark_root/results/production-$(date -u +%Y%m%dT%H%M%SZ)-XXXXXX")
results_root=$(cd "$results_root" && pwd)
mkdir -p "$results_root"/{coverage,lengths,ends,midpoints,concordance,slurm}
mkdir -p "$results_root/midpoints"/{cfdnalab,cfdnalab-gc,griffin-gc}
# Save this assignment to restore the results path in a later shell
printf 'results_root=%q\n' "$results_root"

# --parsable returns job IDs, saved below for the final sacct export
# Each job writes its Slurm logs outside the directory created by the runner
slurm_options=(
    --parsable
    --constraint=gen6
    --output="$results_root/slurm/%x-%j.out"
    --error="$results_root/slurm/%x-%j.err"
)
if [[ -n ${BENCH_EXPECTED_NODE:-} ]]; then
    slurm_options+=(--nodelist="$BENCH_EXPECTED_NODE")
fi

coverage_methods=(
    # Read depth with mate overlap merged where the native tool supports it
    cfdnalab mosdepth samtools-depth pandepth
    # Measure read-depth runtime with cfDNAlab's usual autosome-only scope
    cfdnalab-autosomes
    # Fragment coverage with the normal cfDNAlab pairing policy
    cfdnalab-fragment mosdepth-fragment
    # Regional read depth computed directly in shared 500 kb windows
    cfdnalab-regions mosdepth-regions samtools-regions pandepth-regions
    # Regional fragment counts with the normal cfDNAlab pairing policy
    cfdnalab-region-counts finaletoolkit-coverage
)
length_methods=(
    # Normal directional spans compared with curated cfDNAPro
    cfdnalab cfdnapro
    # Dedicated cfDNA histogram and established insert-size comparators
    finaletoolkit samtools picard
)
end_methods=(
    # Reference-derived 4-mers with the normal cfDNAlab pairing policy
    cfdnalab-reference
    # Both-end FinaleToolkit and native left-end cfDNAPro, both with k=4
    finaletoolkit cfdnapro
)

for cores in 1 4 12 32; do
    export BENCH_CORES="$cores"

    # Coverage: all four quantities, compared separately by the runner
    BENCH_METHODS="${coverage_methods[*]}" \
    BENCH_RESULTS="$results_root/coverage/${cores}cores" \
        bash benchmarks/submit_benchmark.sh coverage \
        "${slurm_options[@]}" --time="$coverage_wall_time" \
        | tee -a "$results_root/job-ids.txt"

    # Lengths: normal configurations used for timing
    BENCH_METHODS="${length_methods[*]}" \
    BENCH_RESULTS="$results_root/lengths/${cores}cores" \
        bash benchmarks/submit_benchmark.sh lengths \
        "${slurm_options[@]}" --time="$lengths_wall_time" \
        | tee -a "$results_root/job-ids.txt"

    # End motifs: the 4-mer comparison
    BENCH_METHODS="${end_methods[*]}" \
    BENCH_RESULTS="$results_root/ends/${cores}cores" \
        bash benchmarks/submit_benchmark.sh ends \
        "${slurm_options[@]}" --time="$end_motifs_wall_time" \
        | tee -a "$results_root/job-ids.txt"

    # Raw midpoint profiles, compared with Griffin's retained raw table
    # Midpoint-specific filters affect only this submission
    BENCH_MIN_MAPQ=20 BENCH_MIN_LENGTH=100 BENCH_MAX_LENGTH=200 \
    BENCH_METHODS=cfdnalab \
    BENCH_RESULTS="$results_root/midpoints/cfdnalab/${cores}cores" \
        bash benchmarks/submit_benchmark.sh midpoints \
        "${slurm_options[@]}" --job-name="midpoints-cfdnalab-${cores}cores" \
        --time="$midpoints_uncorrected_wall_time" \
        | tee -a "$results_root/job-ids.txt"

    # cfDNAlab sample GC estimation plus corrected midpoint profiling
    BENCH_MIN_MAPQ=20 BENCH_MIN_LENGTH=100 BENCH_MAX_LENGTH=200 \
    BENCH_METHODS=cfdnalab-gc \
    BENCH_RESULTS="$results_root/midpoints/cfdnalab-gc/${cores}cores" \
        bash benchmarks/submit_benchmark.sh midpoints \
        "${slurm_options[@]}" --job-name="midpoints-cfdnalab-gc-${cores}cores" \
        --time="$midpoints_cfdnalab_gc_wall_time" \
        | tee -a "$results_root/job-ids.txt"

    # Griffin sample GC estimation plus profiling, retaining raw and corrected tables
    BENCH_MIN_MAPQ=20 BENCH_MIN_LENGTH=100 BENCH_MAX_LENGTH=200 \
    BENCH_METHODS=griffin-gc \
    BENCH_RESULTS="$results_root/midpoints/griffin-gc/${cores}cores" \
        bash benchmarks/submit_benchmark.sh midpoints \
        "${slurm_options[@]}" --job-name="midpoints-griffin-gc-${cores}cores" \
        --time="$midpoints_griffin_gc_wall_time" \
        | tee -a "$results_root/job-ids.txt"
done

# Concordance-only checks at 32 cores, each run once without a warmup
# Separate directories keep their rows outside the normal timing summaries
BENCH_CORES=32 BENCH_WARMUPS=0 BENCH_REPETITIONS=1 \
BENCH_METHODS="cfdnalab-fragment-proper cfdnalab-region-counts-proper" \
BENCH_RESULTS="$results_root/concordance/coverage" \
    bash benchmarks/submit_benchmark.sh coverage \
    "${slurm_options[@]}" --job-name=coverage-proper-concordance-32cores \
    --time="$coverage_wall_time" | tee -a "$results_root/job-ids.txt"

BENCH_CORES=32 BENCH_WARMUPS=0 BENCH_REPETITIONS=1 \
BENCH_METHODS=cfdnalab-proper \
BENCH_RESULTS="$results_root/concordance/lengths" \
    bash benchmarks/submit_benchmark.sh lengths \
    "${slurm_options[@]}" --job-name=lengths-proper-concordance-32cores \
    --time="$lengths_wall_time" | tee -a "$results_root/job-ids.txt"

BENCH_CORES=32 BENCH_WARMUPS=0 BENCH_REPETITIONS=1 \
BENCH_METHODS=cfdnalab-reference-proper \
BENCH_RESULTS="$results_root/concordance/ends" \
    bash benchmarks/submit_benchmark.sh ends \
    "${slurm_options[@]}" --job-name=ends-proper-concordance-32cores \
    --time="$end_motifs_wall_time" | tee -a "$results_root/job-ids.txt"
```

Each timing core directory and concordance-job directory contains `timings.tsv`, `run.tsv`, `runs/`, `profiles/` and setup/progress records. Coverage, lengths and motifs calculate concordance automatically after all tool attempts when a comparable pair is present. Coverage also writes residual and Jaccard tables. The proper-pair checks are compared with their normal 32-core competitor profiles below. Each midpoint job has a single method, so its raw/corrected concordance is also calculated across jobs below. Jobs stage inputs and environments on node-local scratch, leave source inputs unchanged and do not evict caches. All attempted runs, including warmups and failures, remain inspectable.

### Completed-job accounting and concordance

Run this block **after all submitted jobs have finished**. In a new Bash session, first source `benchmarks/local_settings.sh` and restore the `results_root=...` assignment printed above. Completed Slurm accounting supplies the step memory values and final statuses. This collection stays outside benchmark timing. The job ID file comes from `sbatch --parsable`, which can append a cluster name after a semicolon, so `cut` below retains only job IDs.

```bash
set -euo pipefail
# A new analysis directory keeps each raw export with the summaries that use it
# The summarizer also requires new --runs-output files
analysis_root=$(mktemp -d "$results_root/analysis-XXXXXX")
printf 'analysis_root=%q\n' "$analysis_root"

# Collect accounting once for the whole series, including the measured srun steps
job_ids=$(cut -d ';' -f 1 "$results_root/job-ids.txt" | paste -sd, -)
sacct -j "$job_ids" --parsable2 --units=K \
    --format=JobIDRaw,NTasks,MaxRSS,State,ExitCode \
    > "$analysis_root/accounting.tsv"

# Keep each feature's timing summaries and every original row together
for feature in coverage lengths ends; do
    mkdir -p "$analysis_root/$feature"
    python3 benchmarks/summarize_timings.py \
        "$results_root/$feature/"*cores/timings.tsv \
        --accounting "$analysis_root/accounting.tsv" \
        --runs-output "$analysis_root/$feature/runs-with-memory.tsv" \
        > "$analysis_root/$feature/summary.tsv"
done

# Pool midpoint timings only with other midpoint jobs on the same inputs
mkdir -p "$analysis_root/midpoints"
python3 benchmarks/summarize_timings.py \
    "$results_root/midpoints/"*/*cores/timings.tsv \
    --accounting "$analysis_root/accounting.tsv" \
    --runs-output "$analysis_root/midpoints/runs-with-memory.tsv" \
    > "$analysis_root/midpoints/summary.tsv"

# Compare the 32-core proper-pair checks with their matching competitors
# These directories are outside all timing-summary globs above
for feature in coverage lengths ends; do
    comparison="$analysis_root/concordance/$feature"
    mkdir -p "$comparison/profiles"
    for profile in "$results_root/concordance/$feature/profiles/"*; do
        ln -s "$profile" "$comparison/profiles/"
    done
done

# Coverage has a fragment-coverage pair and a regional fragment-count pair
for method in mosdepth-fragment finaletoolkit-coverage; do
    ln -s "$results_root/coverage/32cores/profiles/$method" \
        "$analysis_root/concordance/coverage/profiles/$method"
done
ln -s "$results_root/coverage/32cores/regions.bed" \
    "$analysis_root/concordance/coverage/regions.bed"
for feature in lengths ends; do
    ln -s "$results_root/$feature/32cores/profiles/finaletoolkit" \
        "$analysis_root/concordance/$feature/profiles/finaletoolkit"
done

# The comparison environment supplies zstd for the compressed native tables
PATH="$(dirname "$BENCH_COMPARE_PYTHON"):$PATH" \
    "$BENCH_COMPARE_PYTHON" benchmarks/coverage/compare_coverage.py \
    "$analysis_root/concordance/coverage/profiles" \
    > "$analysis_root/concordance/coverage/concordance.tsv"
PATH="$(dirname "$BENCH_COMPARE_PYTHON"):$PATH" \
    "$BENCH_COMPARE_PYTHON" benchmarks/lengths/compare_lengths.py \
    "$analysis_root/concordance/lengths/profiles" --min-length 30 --max-length 1000 \
    > "$analysis_root/concordance/lengths/concordance.tsv"
"$BENCH_COMPARE_PYTHON" benchmarks/ends/compare_end_motifs.py \
    "$analysis_root/concordance/ends/profiles" --k 4 \
    > "$analysis_root/concordance/ends/concordance.tsv"

# Compare matching midpoint core settings, reusing Griffin's two native tables
# Absolute symlinks avoid copying the saved profiles
for cores in 1 4 12 32; do
    comparison="$analysis_root/midpoints/${cores}cores"
    mkdir -p "$comparison/profiles"
    for method in cfdnalab cfdnalab-gc griffin-gc; do
        ln -s "$results_root/midpoints/$method/${cores}cores/profiles/$method" \
            "$comparison/profiles/$method"
    done
    "$BENCH_COMPARE_PYTHON" benchmarks/midpoints/compare_midpoints.py \
        "$comparison/profiles" \
        > "$comparison/concordance-uncorrected.tsv"
    "$BENCH_COMPARE_PYTHON" benchmarks/midpoints/compare_midpoints.py \
        "$comparison/profiles" --gc-corrected \
        > "$comparison/concordance-corrected.tsv"
done
```

Accounting can be delayed. Missing memory remains `NA`, so rerun the completed-job block later if needed. Failed native/Slurm runs stay in the joined rows and are excluded from successful production averages. Preserve each generated analysis directory with its results.

The runners still support `BENCH_PILOT=1` for debugging with zero full-input warmups and one measurement on a shared-node allocation. The production procedure above uses `BENCH_PILOT=0` throughout.

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

The [completed-job procedure](#completed-job-accounting-and-concordance) collects accounting once for the full production series and retains the raw export with each analysis.

Supply multiple timing TSVs only for blocks of the same feature/experiment with the same biological input. For multiple jobs, collect their IDs together using `sacct -j ID1,ID2`. The summary joins exact job/step IDs and reports sample size, mean/SD, median and quartiles per configuration. `runs-with-memory.tsv` retains every measured point and failure for plotting. Missing/delayed accounting stays `NA`, so recollect and repeat the join later.

Use `slurm_max_rss_kib` as the primary memory result after checking the cluster's accounting plugin, frequency, descendant accounting and file-cache policy once. `slurm-config.txt` records the available configuration. Slurm cgroup v2 MaxRSS can include file cache, so label that policy in results. Memory is collected from Slurm accounting after the job, without an additional sampler in the timed workload. See [Slurm sacct](https://slurm.schedmd.com/sacct.html) and [cgroup v2 limitations](https://slurm.schedmd.com/cgroup_v2.html#limitations).

Effective cores are `(user_seconds + system_seconds) / wall_seconds`. Divide by `cores` for allocated utilization. Use process wall time for this denominator, rather than the flush-inclusive duration.

Inputs and environments are copied to local scratch before warmups. Scheduler startup, profile preservation, cleanup and concordance stay outside measurement. Profiles are copied **after srun** from the first measured repetition at the job's single core count, so copying does not inflate the tool step's accounting peak. No cache eviction is attempted. Staging does not guarantee that the full input remains in page cache.

A recorded native, placement, missing-output or flush failure skips the remaining repetitions of that method and continues other selected methods. Its timing row remains in results, and the job exits nonzero after finishing those other methods. As for every run, `runs/<run_order>/` keeps the exact command and the tool's stdout and stderr. Setup and results-bookkeeping errors remain fatal. Slurm or kernel termination of the batch shell can prevent final bookkeeping and continuation. Previously saved timing rows and measured profiles remain available when a later run fails or the allocation ends. The active run may not finish saving its measurements or output on a hard termination.

The results directory, timing-table header and `run.tsv` are created before input staging. `progress.log` records initialization, environment extraction, input copying, run starts and completed saves, with the same messages printed to Slurm's output. A header-only `timings.tsv` means no completed run has been recorded yet. Check `progress.log` to distinguish setup from an active run. Coverage and lengths check tool versions before copying the large BAM, so a broken launcher fails before that input transfer.

`BENCH_RUN_TIMEOUT_SECONDS=0` disables per-invocation timeouts by default. Set a positive integer to run the native command under GNU `timeout --verbose --kill-after=60`, inside GNU time. The limit is recorded in `run.tsv`. Leave room for the extra 60-second termination grace, staging, flush and bookkeeping below the allocation limit. Status 124 normally identifies an expired limit. Status 137 alone cannot distinguish a forced timeout kill from another SIGKILL, so retain timeout diagnostics and Slurm accounting. Choose an invocation limit that leaves room for flush and bookkeeping within the job allocation. Timeouts are unsuccessful observations and never enter successful production averages.

## Regression tests and validation

Run the automated tests directly in a normal Linux shell. They do not submit Slurm jobs or require an interactive allocation. Slurm commands and benchmark executables are replaced with local fakes. Inputs are small synthetic files created in temporary directories. No production BAM, reference files or benchmark tool installations are needed.

The tests in [test_runners.py](tests/test_runners.py) cover every method's native command and timing settings, the saved profile files, Picard's helper pipeline, midpoint site preparation, setting checks, the measurement step (placement, signal/timeout statuses, missing outputs) and whole lengths and Griffin jobs with a fake srun, including failure continuation. The tests in [test_benchmarks.py](tests/test_benchmarks.py) cover coverage formats/weights/scope, residual distributions, base-resolution Pearson, fractional regional counts and row identity, installer resumption, Picard status, Slurm joins, reconstructed TLEN signs and native left-motif counts. [Additional regression tests](tests/test_review_regressions.py) cover submission, automatic concordance, partial pilot reports, scope validation before loading tracks, native length boundaries and midpoint empty-profile behavior. Expected numerical values are manually derived from small synthetic inputs.

The full suite needs Python >=3.11, Bash >=4, GNU time at `/usr/bin/time`, GNU coreutils and the `zstd` executable on PATH. The local Python package supplies the comparator dependencies. pysam enables the reconstruction cases. From the repository root, with Conda available:

```bash
conda create -n benchmark-tests -c conda-forge python=3.11 pip zstd pysam
conda activate benchmark-tests
python -m pip install ./py-cfdnalab
python -m unittest discover -s benchmarks/tests
```

An existing comparison environment can also run the tests. Missing pysam skips the reconstruction cases. Linux-specific lifecycle tests skip on macOS, so use Linux for the complete suite. The setup above assumes the required Bash and GNU system utilities are already installed.

Real-tool validation is separate and does require a compute-node allocation. After installing the benchmark tools, run each selected native tool and its comparator on a small input to check environment relocation, native output, plotting/fonts and actual Slurm placement/accounting. The automated tests do not perform these cluster checks. No measured benchmark results are included in these documents. The result tables are templates.
