# Benchmark folder review (Claude, 2026-10-30)

Scope: everything under `benchmarks/`. Competitor behavior was checked against the pinned upstream source (links inline). Nothing was executed. *(estimate)* marks numbers I have not measured.

Decisions this review follows:
- Each job gets an exclusive node. Features can run in parallel on separate nodes of the same CPU generation. Within a feature, the tools are compared in the same job. Griffin is the exception: it gets its own node because it is so slow.
- 1 warmup + 3 measured runs for every tool.
- The pilot runs at 4 cores first; 32 cores only if 1-core time is hard to estimate.
- Serial tools are run at every core count to show that they don't scale.
- Griffin uses the paper's full site set.
- FrEIA runs its default workflow (left ends), with counts doubled for comparison.
- cfDNAPro is added for end motifs.

## Summary

- **Within a job, the measurement design is sound.** Each job takes an exclusive node, and every run is bound to its core budget and checked. Inputs are staged to local scratch, so every tool reads from a warm page cache. Native output is written inside the timer, and the filesystem flush is timed separately. Keep this design.
- **Fix three things before the pilot.**
  - A failure or a hang loses the rest of a 24 h job.
  - The default method lists create infeasible jobs.
  - Numerical-library thread settings differ between features.
- **Griffin on its own published site set** (377 TF lists × 10,000 sites) is likely to take tens of hours at 1 core *(estimate)*, while cfDNAlab barely depends on the number of sites. That is the clearest speed result in the benchmark, so present it as such (section 4).
- **End motifs get three competitors:** FinaleToolkit, FrEIA's default left-end workflow and cfDNAPro's `callMotif()`. Both FrEIA and cfDNAPro count only the left end of each fragment (section 5).
- **Much of the harness guards against things that won't happen.** Section 6 lists what to remove.

---

## 1. Fix before the pilot

1. **One failure ends the job.** `benchmark_finish_run` returns the failure status ([common.sh:140-147](common.sh#L140-L147)) and `set -e` stops the runner. In a job holding several methods, an out-of-memory failure in an early method discards all the methods after it. Fix: record the row, skip the remaining runs of *that* method, continue with the next method, and exit non-zero at the end.
2. **There is no per-run time limit.** This affects FinaleToolkit, Griffin and any run that is longer than the job.
   - FinaleToolkit and Griffin use `multiprocessing.Pool`. If the kernel kills a worker, `Pool.map` waits forever (CPython bpo-22393), and the job runs until Slurm stops it.
   - When Slurm kills the job, the row and logs of the run in progress are never written. A run that is simply longer than the job (Griffin at 1 core, section 4) is lost the same way.
   - Fix: add an optional `BENCH_RUN_TIMEOUT_SECONDS` and apply it as `timeout --kill-after=60 "$limit" ...` inside GNU time. Exit status 124 then becomes a recorded "did not finish within X h".
3. **The default method lists create infeasible jobs.**
   - Coverage defaults to 8 methods ([coverage runner:66](coverage/benchmark_coverage.sh#L66)), which the README itself estimates at 92 h for 1 core.
   - Lengths puts cfDNAPro and FinaleToolkit in the same job as everything else ([lengths runner:29](lengths/benchmark_lengths.sh#L29)).
   - Fix: remove the defaults and require `BENCH_METHODS`.
4. **Numerical-library thread settings differ between features.** Lengths sets the OMP/BLAS variables to `cores` ([lengths runner:167-168](lengths/benchmark_lengths.sh#L167-L168)); the other three runners set them to 1. FinaleToolkit runs under both settings. None of the benchmarked code paths use BLAS, OpenMP or parallel numba, so the measured effect should be nil. Fix: set a single policy in `common.sh`. `1` is simplest.
5. **Change the repetition defaults to 1 warmup and 3 measured runs.** The current values are 2 and 5, set in all four runners, the docs and `summarize_timings.py --estimate`.
6. **Schedule the serial length tools at every core count.** [README.md:68](README.md#L68) and [:87-88](README.md#L87-L88) defer cfDNAPro and FinaleToolkit lengths at 12 and 32 cores, pending a 1 vs 4 comparison. Schedule them at all four core counts.
7. **`--estimate` aborts on any failed pilot configuration** ([summarize_timings.py:151-155](summarize_timings.py#L151-L155)), and it prints only a grand total. Fix: print `budget_minutes = mean × 4 × 1.25` per configuration and list the failed configurations.
8. **The coverage comparator stops partway through.** The saved samtools-depth and PanDepth profiles cover chr21/22, while cfDNAlab and mosdepth cover all contigs. Without `--common-scope`, the comparator stops at the first mismatched pair, after computing every full-genome pair before it ([compare_coverage.py:316-317](coverage/compare_coverage.py#L316-L317)). Fix: always compare on the shared chromosomes and remove the flag. `same_native_scope` is already recorded on each row.

## 2. Competitor commands against the pinned source

**No command fails, and none silently computes something other than what the docs say.** Confirmed:
- Picard does not need R ([CollectInsertSizeMetrics.java:92-93](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/CollectInsertSizeMetrics.java#L92-L93)), and `INPUT=/dev/stdin` works.
- samtools stats `-i 1001` puts all inserts ≥1001 into bin 1001, which the comparator discards ([stats.c:1325-1326](https://github.com/samtools/samtools/blob/1.25/stats.c#L1325-L1326), [:2450](https://github.com/samtools/samtools/blob/1.25/stats.c#L2450)).
- FinaleToolkit `frag-length-bins` has no thread option and writes rows with `min == max` ([_frag_length.py:333-349](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/frag/_frag_length.py#L333-L349), [:483-485](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/frag/_frag_length.py#L483-L485)).
- FinaleToolkit `end-motifs` ignores the length bounds ([_end_motifs.py:107-151](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/frag/_end_motifs.py#L107-L151)).
- The cfDNAPro `bamfile=` shortcut note is true ([read_bam_insert_metrics.R:64-71](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/read_bam_insert_metrics.R#L64-L71)).
- All output formats match what the comparators parse.

What affects how the timings should be described:

| Tool | Source finding | What to write |
| --- | --- | --- |
| samtools depth `-aa` | `-@` only adds decompression threads. Counting and the `fputs` text output are single-threaded ([bam2depth.c:857-866](https://github.com/samtools/samtools/blob/1.25/bam2depth.c#L857-L866)). | Mostly measures writing about 3.1×10⁹ lines (about 55 GB *(estimate)*). Expect a flat curve over cores. Check scratch space in the wiring run. |
| PanDepth `-a` | `-t` assigns a worker per contig group, capped at 24 for hg38. `SiteDepth.gz` is written after the workers finish, through single-threaded gzip ([PanDepth.cpp:4264-4284](https://github.com/HuiyangYu/PanDepth/blob/v2.26/src/PanDepth.cpp#L4264-L4284)). | Dominated by output, with a serial tail. |
| mosdepth | `--threads` only adds decompression threads. Per-base output is single-threaded BGZF ([mosdepth.nim:642-643](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L642-L643)). | Plateaus after a few threads, as its own README says. |
| samtools stats | Does mandatory extra work per read: checksums, per-cycle composition and coverage distribution ([stats.c:806-821](https://github.com/samtools/samtools/blob/1.25/stats.c#L806-L821)). | Describe it as "full QC", not a histogram kernel. The README already does. |
| Picard | Counts read2 only, with no MAPQ, proper-pair or QC-fail check ([InsertSizeMetricsCollector.java:71-84](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/directed/InsertSizeMetricsCollector.java#L71-L84)). Default STRICT validation checks every record. | Read2's MAPQ in the upstream `samtools view -q` decides inclusion. Keep the native default and say so. |
| FinaleToolkit coverage | Each worker gets a static chunk of `len/workers` windows ([_coverage.py:244-248](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/frag/_coverage.py#L244-L248)). | Expect imperfect scaling. This is native behavior. |
| cfDNAPro | No parallelism on this path. Reads the whole BAM into memory ([readBam.R:247-251](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/readBam.R#L247-L251)). | Runs at 4/12/32 cores show this directly. Report peak memory as well as time. |

Keep the flush time. On a large-memory node, samtools' 55 GB can still be sitting in the page cache when the process exits. Process time alone would then understate its cost.

## 3. Pilot and job plan

**Wiring run** (tiny BAM, as in [data/README.md:41-48](data/README.md#L41-L48)). Run every method once at 4 cores, then **run the four comparators on its outputs**. This is the only check of native-output parsing before you spend days of node time. Also check free scratch space.

**Pilot: full input, 4 cores, `BENCH_PILOT=1`.** Run a job per feature with every method, and Griffin on site subsets (section 4). This gives T(4), effective cores and peak memory. For CPU-bound tools on a warm cache, `user + system` seconds at 4 cores should be an upper bound on 1-core wall time *(inference)*, because parallel overhead only adds CPU time. If that isn't enough to judge the 1-core runs, add 32-core pilots.

**Packing jobs.**
- **Budget:** each configuration needs `4 × T × 1.25`. A job holds up to about 22 h, which leaves room for staging.
- **Features in parallel:** coverage, lengths, ends and midpoints run on separate exclusive nodes with the same `--constraint` (CPU generation). Only staging and result copying touch the shared filesystem, and both are outside the timer. So jobs on different nodes don't affect each other's timings.
- **Within a feature:** all of a feature's tools at a given core count go in the same job. Griffin is the exception (section 4). If a block exceeds about 22 h, split it into as few jobs as needed. The 1-core blocks will need this. Each configuration, cfDNAlab included, runs in exactly one job, so every configuration × core count has exactly 1+3 runs.
- **Runs longer than a job:** if a single run can't fit in 24 h, record it as a timeout (1.2).

**Node time** *(estimate)*. This uses the README's placeholder runtimes: cfDNAlab 12 min and competitors 120 min at 1 core, linear scaling, mosdepth capped at 4×, serial length tools, 1/4/12/32 cores, 1+3 runs and a 25% margin.

| Part | Node-hours |
| --- | --- |
| Current matrix without Griffin (10 competitor and 8 cfDNAlab configurations; serial tools at all four core counts) | about 208 h |
| FrEIA left-end workflow, serial at every core count | about 40 h |
| cfDNAPro end motifs, serial at every core count | about 40 h |
| Griffin on the full TFBS set, plus cfDNAlab midpoints (section 4) | about 70–140 h |
| **Total node time** | **about 360–430 h** |

Calendar time depends on how many same-generation nodes you get at once, because no job exceeds about 24 h. The 1-core competitor runs and Griffin account for most of the node time. Replace these numbers with the pilot measurements.

## 4. Griffin: use the paper's full site set and make the result prominent

**Which sites.** The paper's TFBS lists are in [adoebley/Griffin_analyses `sites/TFBS/10000_unfiltered_sites_CIS_BP_v2`](https://github.com/adoebley/Griffin_analyses/tree/75467e438a88029206d8d778ddce05eb0ab6efa8/sites/TFBS/10000_unfiltered_sites_CIS_BP_v2). It contains 377 files, `<TF>.hg38.10000.txt`, about **3.77M sites** in total. Cite commit `75467e4`.
- Columns are `Chrom Start End uniprotId Gene exp.count peak.count position`. There is no `Strand` column, so the sites are unstranded.
- Convert them once to BED4: `Chrom, position - HALF, position + HALF, <TF>`. The harness then gives Griffin the same windows.
- **Use ±1000 bp (2 kb).** That is Griffin's own `save_window`. The harness uses the BED width for every Griffin window, including normalization. Griffin's default `norm_window` is ±5000, which would make it slower still, so ±1000 favors Griffin. Say so in the methods.
- Griffin v0.2.0 itself ships only a 1,000-site CTCF demo list. That is too small to show anything.

**Why it is slow** (confirmed in [griffin_coverage.py v0.2.0](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_coverage.py)):
- For every forward read near a site, it fetches the reference sequence and computes GC in Python. Per fragment, that means `np.array(list(seq))`, four `np.isin` calls and a new `np.random.default_rng`. All of this happens before checking whether the midpoint is in the window ([L366-L393](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_coverage.py#L366-L393)).
- It reopens the BAM and FASTA for every merged interval ([L345-L348](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_coverage.py#L345-L348)).
- It writes all intervals into two bigWigs serially, in the parent process.

**Expected cost** *(estimate)*: roughly 10 µs per window bp plus 10 ms per site, per core.
- At 2 kb that is about 30 ms per site, or **up to about 30 h for a single 1-core run** on 3.77M sites. Overlap between the 377 lists lowers this, because Griffin fetches merged intervals once.
- cfDNAlab's midpoint cost barely depends on the number of sites. The gap should be the most striking speed result in the benchmark. Put the absolute times side by side, together with site count, window and group count.
- With 377 groups, the merge step parallelizes well ([griffin_merge_sites.py:609-612](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py#L609-L612)). The serial bigWig export of 3.77M intervals limits scaling at 32 cores.
- The parent process holds every interval's data frame until export, which is tens of GB *(estimate)*. That is fine with your RAM.

**Running it.**
- **Pilot:** at 4 cores, run random 1% and 4% subsets of every list. Fit time against merged window bp and extrapolate to the full set. Griffin's own log timestamps ("Done with fetch", "Done with export") separate the parallel and serial parts.
- **Jobs:** Griffin gets its own node of the same CPU generation. cfDNAlab midpoints runs once per core count, in a single one of the Griffin jobs, not in every job. A Griffin configuration whose four runs don't fit in a job gets a job per run, and those jobs can run in parallel. At this scale, node-to-node variation and warmup effects are noise: staging already warms the page cache, and startup takes seconds against hours. Keep 3 measured runs.
- **Summaries:** `summarize_timings.py` groups rows by `node` ([summarize_timings.py:93-98](summarize_timings.py#L93-L98)), so Griffin runs from separate nodes show up as separate rows. Pool those Griffin rows by `cpu_model`.
- **1 core:** if a 1-core run exceeds the 24 h limit, it cannot complete. Use a longer partition if one exists; otherwise record it once as "> 24 h" with the timeout from 1.2. Repeating a censored run adds nothing.

## 5. End motifs: FinaleToolkit, FrEIA (left ends) and cfDNAPro

### FrEIA, using its default workflow
- **Drop `--both`.** The authors' Snakemake rule doesn't pass it ([2_FrEIA.smk:101-111](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/rules/2_FrEIA/2_FrEIA.smk#L101-L111)), and the aggregation bug occurs only with two ends.
- **What it counts.** The left end only: the forward mate's first n aligned read bases, 5'→3' ([1_extract_fragment_ends.py:45-51](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/scripts/1_FrEIA/1_extract_fragment_ends.py#L45-L51)).
  - Its trinucleotide table is fixed at K=3. Compare it with `cfdnalab-read-aligned` at `BENCH_K=3`: read bases with soft clips excluded match `query_alignment_sequence`.
- **Comparator.**
  - Read `T__sample.pq` (add `pyarrow` to the comparison environment).
  - Take a per-motif count from the `leftmost` rows: `CountSamp` is already the per-motif total on each length row ([2_data_transformation.py:130-134](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/scripts/1_FrEIA/2_data_transformation.py#L130-L134)).
  - Multiply by 2 before comparing counts. Frequencies are unaffected by the factor. Doubling assumes the right ends have the same motif distribution as the left ends (strand symmetry), which should hold closely genome-wide *(inference)*.
- **Remove the `--both` special cases.** These are the comparator skip, the capture exclusion ([common.sh:148](common.sh#L148)) and the worker-count rewrite ([ends:327-334](ends/benchmark_end_motifs.sh#L327-L334)). For this BAM, extraction is always a single worker: the filtered BAM exceeds the 8000 MiB threshold ([1_extract_fragment_ends.py:137-140](https://github.com/mouliere-lab/FrEIA/blob/2c32fb9d9fab774f4ae57931df7fd7d75654d343/workflow/scripts/1_FrEIA/1_extract_fragment_ends.py#L137-L140)). State that in the docs instead.

Text for the methods (replace the bracket with your issue link):

> FrEIA was run with its default aggregation, which counts the 5' motif of the left (forward-strand) end of each fragment. Its `--both` option, which adds right ends, is not part of the authors' workflow. In the pinned version (`2c32fb9`), it assigns the right-end table to the left-end table's row index. Right-end rows are therefore dropped or padded with missing values whenever the two ends have different numbers of (length, motif) combinations [issue]. For comparison with both-end profiles, FrEIA's left-end counts were doubled. This does not affect motif frequencies.

**On the "second bug".** I overstated it. With `--both`, FrEIA reads the right end as reversed reference-strand bases, without complementing them. The code comments this as intended. It is a convention, and it only matters if right-end rows are pooled with conventional 5' motifs. It is no longer relevant now that `--both` is not used.

### cfDNAPro `callMotif()`
It is already in the installed image (from the [callMotif source](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/callMotif.R) and [get_motif](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/internalFunctions.R#L676-L733)).

**What it computes.**
- `callMotif(fragments, motif_type = "s", motif_length = K)` takes reference bases `[start, start + K)` of the unstranded fragment. That is the left end only, like FrEIA. Double its counts the same way.
- `"e"` returns the right end's reference bases without reverse complement, so don't use it as a 5' motif.
- Output is `motif`, `n`, `fraction`, with all 4^K motifs (missing ones as 0) and motifs containing N removed.
- Bases come from the image's `BSgenome.Hsapiens.UCSC.hg38`. Compare against `cfdnalab-reference`.

**Implementation.**
- Add a `cfdnapro_end_motifs.R` with the same `readBam()` call as [cfdnapro_lengths.R](lengths/cfdnapro_lengths.R).
- Subset to 30–1000 bp: `callMotif()` has no length bounds, and a one-line `width()` filter matches cfDNAlab's bounds.
- Then call `callMotif` and write `motif` and `n`.
- Add a `cfdnapro` method to the ends runner, reusing the lengths runner's Apptainer command.

**What it shows.** It is timed from BAM to profile, like every other method. It loads the whole BAM into R and has no parallelism, which is exactly the contrast you want to show.

## 6. Overengineering to remove

Each item guards against something you control or that the wiring run shows immediately. Each line names the code to delete.

**Runners and `common.sh`**
- **Hostname checks** inside every step ([coverage:373](coverage/benchmark_coverage.sh#L373), and the same line in each runner) and `BENCH_EXPECTED_NODE` ([common.sh:28-31](common.sh#L28-L31)). An `srun --jobid` step in a single-node allocation cannot run elsewhere, and `run.tsv` records the node for each job.
- **CPU-topology block:** 18 lines copied four times ([coverage:248-262](coverage/benchmark_coverage.sh#L248-L262)). Keep the `nproc == cores` check and record `Cpus_allowed_list`. Save `lscpu -p` once per job if you want the socket and core map.
- **Fallback-row reconstruction for failures before GNU time** (`step.id`, column indices per runner): [coverage:395-407](coverage/benchmark_coverage.sh#L395-L407), [lengths:304-316](lengths/benchmark_lengths.sh#L304-L316), [ends:344-356](ends/benchmark_end_motifs.sh#L344-L356), [midpoints:330-342](midpoints/benchmark_midpoints.sh#L330-L342). The Slurm log already shows what failed. A row with method, repetition and status is enough.
- **Leftovers from multi-core jobs:** `cores_to_run` arrays, per-core validation loops, `BENCH_CAPTURE_CORES`, the `cores` column in `order.tsv` and the `cores == capture_cores` test ([coverage:71-90](coverage/benchmark_coverage.sh#L71-L90), [:188](coverage/benchmark_coverage.sh#L188), [common.sh:137-148](common.sh#L137-L148)). `BENCH_CORES` is a single value now.
- **Argument validation:** integer and range checks for MAPQ, lengths and repetitions, the BAI-vs-CSI check and index-extension handling ([coverage:77-96](coverage/benchmark_coverage.sh#L77-L96) and the equivalents). The tools validate their own arguments.
- **Autosome and capture-contig detection for both `chr1` and `1` naming,** the `BENCH_CAPTURE_CHROMOSOMES` override and the first-contig fallback ([coverage:154-161](coverage/benchmark_coverage.sh#L154-L161), [:172-185](coverage/benchmark_coverage.sh#L172-L185), [common.sh:167-168](common.sh#L167-L168)). Hard-code `chr1..chr22` and `chr21,chr22`.
- **`ldd` parsing in `benchmark_stage_binary`** ([common.sh:59-68](common.sh#L59-L68)). A missing library fails the wiring run anyway.
- **Picard version probe** that accepts status 1 with a regex ([common.sh:6-15](common.sh#L6-L15)), and the PanDepth version scraping that fails on mismatch ([coverage:168-171](coverage/benchmark_coverage.sh#L168-L171)). Just record the version output.
- **Midpoint BED checks** for ambiguous strand columns, all-numeric chromosome names, BAM/FASTA length match and `.gz` FASTA ([midpoints:70-73](midpoints/benchmark_midpoints.sh#L70-L73), [:148-149](midpoints/benchmark_midpoints.sh#L148-L149), [:156-159](midpoints/benchmark_midpoints.sh#L156-L159), [:176-177](midpoints/benchmark_midpoints.sh#L176-L177)). Keep the equal/even-width check and the strand parsing.
- **End-motif header check:** the BAM/2bit contig comparison and the 1 Mb-multiple guard ([ends:139-157](ends/benchmark_end_motifs.sh#L139-L157)). FinaleToolkit validates contigs itself at startup, and the 1 Mb bug doesn't affect hg38's main chromosomes. A sentence in the docs is enough.

**Installer**
- **The `snyder-preparation` environment** ([install_tools.sh:318-329](install_tools.sh#L318-L329)). The BAM exists. Keep `prepare_snyder.py` as a record of how it was made, and drop its test.
- Keep the manifest and resume logic.

**Summarizer and comparators**
- **`summarize_timings.py`:** the conflicting-accounting and NTasks checks and the MaxRSS unit-suffix parser ([summarize_timings.py:18-39](summarize_timings.py#L18-L39), [:49-59](summarize_timings.py#L49-L59)). With `--units=K` the unit is always K.
- **`compare_coverage.py`:** the SQLite residual store with 4096-value batching ([:168-176](coverage/compare_coverage.py#L168-L176), [:291-294](coverage/compare_coverage.py#L291-L294)). Residuals of unweighted integer depths take a few hundred distinct values, so a `Counter` holds them.
  - A larger, optional change: per-chromosome NumPy arrays instead of the pure-Python run merge. That is less code and far faster; each full-genome pair currently takes on the order of an hour *(estimate)*.
- **Comparator value checks:** the non-finite and negative checks in every comparator, and Griffin's mean-one tolerance and empty-profile status categories ([compare_midpoints.py:64-81](midpoints/compare_midpoints.py#L64-L81)). Keep the structural checks (axis, labels, group/site agreement); they catch real format mistakes.

**Docs.** README.md, two specs, the plan and four feature READMEs add up to about 1,000 lines. Run controls, Slurm memory, staging, failure rows and budgets are each described three to five times, and the copies have started to disagree. For example, [lengths/README.md:150-155](lengths/README.md#L150-L155) and [midpoints/README.md:90-91](midpoints/README.md#L90-L91) list cores 1/4/8/12/32. Proposal:
- `README.md` covers how to run, measurement definitions and the job plan.
- `METHODS.md` covers per-feature configurations and native differences, with source links. It includes the FrEIA paragraph above.
- Each feature README keeps only its run example.
- The specs are deleted, and the plan moves to `.AI/docs/future/`.

## 7. Tests

- **Missing:** there are no tests for `compare_lengths.py`, `compare_midpoints.py` or `benchmark_time_command`. Worth adding:
  - samtools IS rows at 29/30/1000/1001 (expect 29 and 1001 to be dropped);
  - a Picard histogram without an FR column;
  - a command that exits 3 and one killed with `kill -9 $$`, passed through `benchmark_time_command`;
  - a small FrEIA left-end Parquet fixture and a cfDNAPro motif table, both doubled.
- **Remove together with the code above:** the Picard status test, the FrEIA skip test, the Snyder reconstruction test and the SQLite batch test.
