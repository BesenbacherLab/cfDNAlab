# Benchmark installation and execution specification

This describes implemented benchmark installation and execution behavior. Scientific comparisons are specified in [benchmarks_concordance_spec.md](benchmarks_concordance_spec.md). Unimplemented work and outstanding verification are listed in [benchmarks_review_plan.md](benchmarks_review_plan.md).

The benchmark runs selected native commands repeatedly, records timing and Slurm step identifiers, and preserves measured outputs for concordance. It does not change cfDNAlab's counting definitions or patch competitors. There is no automatic job planner, source fingerprinting system or custom memory sampler.

## Implementation

| Responsibility | Source |
| --- | --- |
| Installer and explicit reinstall controls | [install_tools.sh](install_tools.sh) |
| Shared setup, scheduling, measurement and preservation | [common.sh](common.sh) |
| Coverage | [benchmark_coverage.sh](coverage/benchmark_coverage.sh) |
| Lengths | [benchmark_lengths.sh](lengths/benchmark_lengths.sh) |
| End motifs | [benchmark_end_motifs.sh](ends/benchmark_end_motifs.sh) |
| Midpoints | [benchmark_midpoints.sh](midpoints/benchmark_midpoints.sh) |
| Summaries, accounting join and estimates | [summarize_timings.py](summarize_timings.py) |
| BAM preparation utility | [prepare_snyder.py](data/prepare_snyder.py) |
| Prepared regression fixtures | [test_benchmarks.py](tests/test_benchmarks.py) |

Submission examples and operating instructions are in the [benchmark README](README.md).

## Installation

### Platform and versions

Installation targets Linux x86_64 with Bash 4 or later. Conda, Apptainer, Git, wget, tar, GNU make and C/C++ build tools are external prerequisites. A dedicated installation root contains `envs/`, current attempt logs under `logs/`, dependency records under `records/`, and standalone binaries, source checkouts, images and archives at the root.

| Installation ID | Requested version/source | Main retained artifact |
| --- | --- | --- |
| `packaging` | conda-pack 0.9.2, Python 3.11 | `envs/packaging/bin/conda-pack` |
| `cfdnalab` | Current checkout, version from Cargo.toml | `envs/cfdnalab/bin/cfdna`, `cfdnalab-build.tsv` |
| `mosdepth` | 0.3.14 | Authors' standalone executable |
| `samtools` | 1.25 | Conda executable and `samtools.tar.gz` |
| `pandepth` | 2.26 | Authors' Linux release executable |
| `finaletoolkit` | 1.1.0, Python 3.11 | Conda executable and `finaletoolkit.tar.gz` |
| `cfdnapro` | 1.7.4, source `ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c` | Apptainer image pinned by digest |
| `picard` | 3.5.0, Java 17 | Conda launcher and `picard.tar.gz` |
| `freia` | `2c32fb9d9fab774f4ae57931df7fd7d75654d343` | Unmodified scripts and `freia.tar.gz` |
| `griffin` | v0.2.0 | Unmodified scripts and `griffin.tar.gz` |
| `comparison` | Local py-cfdnalab checkout, Python 3.11 and zstd | `envs/comparison/bin/python` |
| `snyder-preparation` | pysam 0.24.1, Python 3.11 | `envs/snyder-preparation/bin/python` |

The cfDNAPro image digest is `sha256:26a28821a272856e77fdf3d13d902ce8ca5c1e60475b1d4f2f076ef2ac3e2b45`.

cfDNAlab builds with Rust 1.94.0, LLVM/Clang 21, Fontconfig, zstandard and Perl. Its Cargo target directory is `<root>/build/cfdnalab`, surviving explicit reinstalls for incremental builds. The build record contains the installed checkout commit, dirty flag and compiler version.

FinaleToolkit additionally pins pysam 0.24.0, pybigwig 0.3.25, py2bit 0.3.3, numba 0.63.1, scipy 1.17.0, statsmodels 0.14.6, tqdm 4.67.3, loess 2.1.2 and matplotlib 3.10.8. These recipe pins do not independently fix every transitive dependency. The solved environment export records the complete result.

Griffin uses Python 3.7.12, `nomkl`, numpy 1.21.2, pandas 1.3.2, scipy 1.7.1, pybigwig 0.3.17, PyYAML 5.3.1, bedtools 2.29.2, pybedtools 0.8.0 and pysam 0.19.0. Python 3.7.12 replaces the unavailable 3.7.4 guide pin.

Conda creates use explicit channels, `--override-channels` and strict channel priority. FrEIA preserves its supplied dependency pins in a derived YAML with `nodefaults`. Its source checkout and derived YAML remain at the installation root. Python user packages and inherited Python import paths are disabled during installation and job setup.

Picard's version probe accepts status 0 or 1 only with a `Version:3.5.0` field. Missing/wrong versions or other statuses fail. The installer and lengths setup share this behavior.

### Resumption and explicit controls

The interface is `ROOT [--only TOOLS] [--force TOOLS] [--mark-outdated TOOLS]`. Selections are comma-separated installation IDs or `all`. The default installation set is all tools. Force applies within the selected set. Packed-tool subsets automatically include the packaging helper first.

`manifest.tsv` contains `tool`, `version`, `status`, `artifacts`. A tool skips only when its record is `complete`, its version matches and all expected artifacts exist. Version mismatch, missing artifacts, outdated status or explicit force causes reinstallation. Marking outdated updates records and exits without installing.

A reinstall records `outdated` before work and `complete` only after successful installation and artifact checks. Steps run under active shell error handling, with failures reported through per-tool logs. Cleanup removes only known paths for that tool and refuses cleanup through symlinked environment paths. Unrelated completed tools remain available.

Same-version source/dependency changes require explicitly forcing or marking outdated. The recorded commit supplies provenance and does not automatically trigger rebuilding. There are no recipe hashes or source-tree checksum scans on resume.

Each Conda step saves an explicit environment export. FrEIA and comparison environments also save pip records. Image metadata and cfDNAPro/GenomicAlignments/Rsamtools versions are retained. Path exports are printed only for matching completed installations with existing artifacts.

## Biological inputs

### Existing BAM and reference

Runners consume an existing coordinate-sorted BAM and index. They copy inputs into scratch and do not modify source files. Coverage, lengths and ends require BAI. Midpoints copies the supplied index with its extension.

The selected BH01 route is the published hs37d5 input, deterministic preparation, then the [pinned lionheart hg38 remap](https://github.com/BesenbacherLab/lionheart/blob/085324a/remap/workflow.py), including duplicate marking. The user's already remapped BAM is the experimental input. Extra preparation statistics do not require rebuilding it.

Reference-derived end motifs require `BENCH_REF_2BIT`. Griffin requires uncompressed `BENCH_REF_FASTA` and its existing FAI. These must represent the alignment reference. Equal contig names/lengths do not prove equal reference bases. Actual checksums, remapping versions/commands and final-input QC are experimental records supplied outside timing, described in [data/README.md](data/README.md).

### Preparation utility

The current utility hashes `seed:query_name` using 64-bit BLAKE2b and compares it to the requested retention threshold. Seed defaults to 42, threads to 4, and retention probability must be in `(0,1]`. Both mates and additional alignments with the same name share the sampling decision. Retained existing paired records pass through unchanged.

Unpaired records must have an `M_` name. Validation occurs before sampling, including records that would otherwise be discarded. It stops at the first offending name and reports its input record number.

Each retained recognized merged record becomes two mates at the original aligned coordinates, preserving stored sequence, qualities, CIGAR and filtering flags. Read1 keeps the original strand, read2 takes the opposite strand. Mapped forward mates have positive TLEN and reverse mates negative TLEN. Unmapped pairs keep TLEN zero and are not marked proper. Mapped synthetic pairs are marked proper as a preparation convention. SA/XA tags are removed and mate CIGAR/MAPQ tags updated.

The input header is preserved. Existing output BAM/BAI paths are rejected, output is indexed and record/reconstruction counts are logged. Separate original mate sequences/qualities are not recovered. Probability 0.15625 targets approximately 15x from the published approximately 96x depth and does not measure post-remap depth. BAM-to-FASTQ conversion discards preparation TLEN and remapping recalculates alignment fields.

### Supplied midpoint sites

`BENCH_INTERVALS` is the user's BED. No site-generation helper, selection, deduplication, artificial strands, 2 kb width, site-count limit, TF source or TSS group is imposed.

Rows require chromosome, start, end and group. BED5 uses field 5 as strand, BED6 uses field 6, BED4 is unstranded. Strands are `+`, `-` or `.`. Groups must be nonempty without whitespace. Intervals must be nonempty, in bounds, of equal even width and on BAM contigs. When Griffin is selected, those contigs must have matching FASTA lengths, and each group must not have entirely numeric chromosome names.

Setup retains all rows, including repeats and overlaps, sorts chromosome/start/end and writes normalized BED6 with score zero. Griffin receives equivalent per-group anchors centered at `start + width/2`. This converts the supplied windows to native formats. The normalized BED and source SHA-256 are retained with results.

## Run controls

Submit from the repository root so `SLURM_SUBMIT_DIR` locates helpers. Nodes need Slurm, GNU time at `/usr/bin/time`, GNU coreutils, `ldd` and setup Python 3. Selected tool paths/archives come from installer exports.

| Variable | Default and behavior |
| --- | --- |
| `BENCH_BAM` | Required existing BAM |
| `BENCH_INDEX` | `<BENCH_BAM>.bai` |
| `BENCH_LOCAL_ROOT` | Required existing node-local scratch root |
| `BENCH_RESULTS` | Required new results directory |
| `BENCH_METHODS` | Feature-specific space-separated method labels |
| `BENCH_CORES` | `1 4 8 12 32` |
| `BENCH_WARMUPS` | **2** per scheduled method/core |
| `BENCH_REPETITIONS` | **5 measured** per scheduled method/core |
| `BENCH_PILOT` | `0`, or `1` to use zero warmups and one measurement |
| `BENCH_SETUP_PYTHON` | `python3`, untimed setup/relocation only |
| `BENCH_EXPECTED_NODE` | Optional exact hostname check before staging |
| `BENCH_SERIAL_ENDPOINTS` | `1` restricts cfDNAPro/FinaleToolkit lengths to selected 1/32-core endpoints. `0` uses all selected budgets |
| `BENCH_MIN_MAPQ` | 30 |
| `BENCH_MIN_LENGTH`, `BENCH_MAX_LENGTH` | Inclusive 30-1000 bp, midpoint defaults 100-200 bp |
| `BENCH_K` | Ends only, 4 with 3 also accepted |
| `BENCH_GENOME` | cfDNAPro label, `hg38` |
| `BENCH_FASTA_INDEX` | `<BENCH_REF_FASTA>.fai` for Griffin |
| `BENCH_CAPTURE_CHROMOSOMES` | Override for saved coverage subset |
| `BENCH_INSTALL_MANIFEST`, `BENCH_CFDNA_BUILD` | Optional records copied beside results |

Methods/cores must be nonempty. Core budgets must fit the allocation. Integer filters and repetitions are validated, measured repetitions must be positive, MAPQ is within 0-255, and fragment bounds satisfy `10 <= min <= max <= 50000`.

All runners use the samtools archive for header preparation and setup Python for relocation. Competitor environments and references are required only by selected methods using them. Apptainer is required only for cfDNAPro. Staged native executables are checked with `ldd`, allowing static binaries and rejecting missing shared libraries. Font availability still requires a real small-input run.

## Scheduling and placement

Jobs request one exclusive node, one task, 32 CPUs per task, one hardware thread per core, all node memory and no requeue. Partition, node selector and wall time are submission choices. There is no universal production wall-time allocation.

Each method runs as a contiguous block. Its warmup repetitions precede measurements. Core budgets shuffle within each repetition, with actual order retained in the timing TSV. Production repetitions 1-2 are warmups and 3-7 are measured.

cfDNAPro and FinaleToolkit lengths use requested 1/32-core endpoints unless `BENCH_SERIAL_ENDPOINTS=0`. Selecting neither endpoint fails. Picard keeps the full requested matrix because its timed pipeline includes samtools and the JVM.

Sequential single-task `srun` calls attach to the existing allocation with exact resources and thread binding. Inside each step, hostname and logical/physical CPU counts must match the expected node and core budget. Affinity and socket:core IDs are recorded. Hard confinement and descendant accounting depend on cluster configuration.

Inputs/environments are staged once before warmups. Packed environments are relocated with `conda-unpack` through setup Python. Host caches are not explicitly evicted. Coverage, ends and midpoints limit numerical libraries to a thread, while lengths allows them the allocated budget. Native worker/thread settings remain distinct from allocated cores.

Oversized series are split explicitly by method/core and submitted with `--nodelist` and `BENCH_EXPECTED_NODE` to retain the same physical node. Each block includes its own warmups. Scripts do not coordinate multiple jobs.

## Measurement and failures

GNU time measures the native command or documented complete pipeline, including startup, output writing/compression, compulsory calculations and enabled default plots. Setup, scheduler startup, placement checks, profile copying, cleanup and concordance are outside the process timer.

After native success, a second timer measures `sync -f` on the output filesystem. This can include other filesystem writes and is a filesystem flush boundary, not exact per-tool durable-write latency. `du -sb` measures apparent native output-directory size including required intermediates. Staging writes are flushed before the first invocation.

| Measurement | Meaning |
| --- | --- |
| `wall_seconds`, `user_seconds`, `system_seconds` | Native process/pipeline GNU time |
| `max_rss_kib` | GNU time maximum RSS, separate from Slurm memory |
| `tool_exit_code` | GNU time native exit status |
| `flush_seconds` | Post-command filesystem flush duration |
| `process_plus_flush_seconds` | Process wall plus flush, the selected end-to-end result |
| `output_bytes` | Apparent output-directory bytes from `du -sb` |
| `step_exit_code` | Scheduler-step status including placement/output/flush failures |

Rows additionally record run order, repetition, warmup/pilot labels, method/version, core/native workers, feature filters/scope, node/CPU placement, job/step IDs and UTC start. Schemas have 31 coverage, 30 lengths, 35 ends and 32 midpoint columns.

Native, placement, expected-output or flush failure stops the job after saving available timing information. Missing measurements are `NA`. Pre-timer failures retain known run/method/core/job identifiers and any recorded step ID. Failed tool logs are copied to `failed-run-<run_order>/` before scratch removal when parent bookkeeping can complete. Whole-allocation termination can prevent that final bookkeeping.

## Result artifacts

`timings.tsv` is written incrementally and keeps warmups, pilots and failures. `run.tsv` records feature, node, requested methods/cores, repetitions, benchmark commit/dirty flag and input paths. Optional installation/build TSVs are copied beside it. Available Slurm configuration is captured once per job, with diagnostics and a warning if unavailable.

One native output per method is preserved from the first measured repetition at the first requested capture core budget. Endpoint-restricted lengths uses the first selected endpoint. Copying happens after `srun`, outside tool-step memory accounting. Other temporary outputs are removed. FrEIA profiles are never preserved for concordance.

Coverage keeps full compressed cfDNAlab/mosdepth tracks. Saved samtools/PanDepth positions cover chr21/chr22 or 21/22, falling back to the first BAM contig when absent. `BENCH_CAPTURE_CHROMOSOMES` overrides the comma-separated subset and unknown requested names fail. Full native output is still generated inside timing. Subset extraction/recompression occurs afterward. Each coverage profile has its processed `chrom.sizes`, and regional outputs share saved `regions.bed`.

Midpoints also retains normalized `sites.bed` and source `sites.sha256`. The original input BED remains unchanged.

## Slurm accounting and summaries

Collect accounting once after completion using `sacct --parsable2 --units=K --format=JobIDRaw,NTasks,MaxRSS,State,ExitCode` and retain the raw export. The summary joins `JobIDRaw` to exact `job_id.step_id`. Conflicting records or a known task count other than 1 fail. Missing records/blank memory remain `NA`.

The primary memory result is `slurm_max_rss_kib` under the site's documented accounting policy. Its step includes shell/time wrappers, the native command and children, flush and size measurement. Staging and profile copying are outside it. Plugin, frequency, descendant coverage and file-cache treatment must be established on the cluster. See [sacct](https://slurm.schedmd.com/sacct.html) and [cgroup v2 limitations](https://slurm.schedmd.com/cgroup_v2.html#limitations).

Production summaries select successful `warmup=false`, `pilot=0` rows. Native and step statuses must be zero. Known accounting must be `COMPLETED` and `0:0`. Missing accounting does not invalidate process timings. Versions, nodes and recorded command settings remain separate groups. TSVs supplied together are expected to belong to the same feature/experiment and biological input.

Per configuration, process-plus-flush time, process/flush separately, Slurm/GNU memory, output bytes and effective cores each report `n`, mean, sample SD, median and inclusive lower/upper quartiles. A singleton has undefined sample SD. Missing metrics have a smaller `n`. Effective cores use `(user_seconds + system_seconds) / wall_seconds`, with zero wall time yielding `NA`. Allocated utilization divides effective cores by the core budget.

`--runs-output` writes original rows with joined memory/status to a new file, retaining points and failures. The utility does not create plots.

### Pilot estimates

`--estimate` summarizes successful measured pilot rows and defaults to two planned warmups, five planned measurements and a 25% margin. It sums observed configuration mean process-plus-flush times and multiplies by seven. Staging/cleanup are added separately. An attempted method/core without a successful measured pilot row causes an incomplete-budget error. Absent configurations cannot be estimated.

At working estimates of 12 minutes for cfDNAlab and 120 minutes for a comparator, those two methods at one core alone require **15 h 24 min** for seven invocations each, or **19 h 15 min** with margin. Other methods/core budgets add work. A serial 120-minute comparator at both 1 and 32 cores requires 28 hours before margin.

## Verification boundary

Source inspection and static/Cargo checks do not establish successful Linux installation or correct cluster measurements. Prepared fixtures have not been run during implementation. Linux solves/relocation, controlled native cases, accounting, real-input resource use and publication results remain in the [plan](benchmarks_review_plan.md).
