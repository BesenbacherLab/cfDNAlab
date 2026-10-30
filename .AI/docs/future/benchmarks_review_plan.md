# Benchmark review and proposed changes

Status: benchmark repairs implemented. Linux installation, accounting and experimental runs remain to be verified. Reviewed on 2026-10-29 against checkout `d366d65bc8ea485e3d5a8f48cc7526def36eba83`.

## Assessment and scope

The benchmark needs changes before publication runs. The strongest findings in `benchmarks/claude_review.md` are the Picard version-check failure, Griffin's unavailable Python pin, FrEIA's incorrect right-end aggregation, missing coverage concordance, and misleading process-tree memory comparisons. Several proposed remedies need correction before implementation.

I reread every existing file in `benchmarks/`, the relevant cfDNAlab counting/filtering and Python-reader code, and pinned upstream source for the disputed behavior. I independently queried official package metadata for Griffin's Python availability, FinaleToolkit's dependencies and conda-pack's current version. No tests, benchmarks, installations or Conda solves were run. Linux binary linkage, cluster behavior, full-input runtime/memory, and the final remapped BAM were not verified here. Claude's sampled BH01 percentages and installation solves are reported evidence, not measurements reproduced by this review.

The decisions D1-D5 recorded in Claude's document are the working preferences for this plan: retain the hg38 remap route, use five measured repetitions, run on the same physical node in contiguous method blocks, retain native per-position output costs, pilot cfDNAPro's memory use, and add mosdepth fragment mode. The implementation follows these preferences and the later decision to exclude broken FrEIA output from concordance. No external issue submission is authorized here.

The plan changes benchmark scripts and documentation. It does not change cfDNAlab's biological definitions or patch competing tools. In particular, fragment spans remain directional: `forward.pos` to `reverse.reference_end`.

Implementation should stay small: fix the known failures, run each selected method for the requested repetitions, save its timing and Slurm step memory, and compare a retained output. Use explicit reinstall controls and a short installation record. Automatic source-change detection, a custom memory monitor, a new helper environment and compulsory drift-monitoring runs are outside this plan.

## Findings checked against the review

The references below use Claude's finding numbers. Local line references describe the reviewed checkout.

| Finding | Assessment | Proposed disposition |
| --- | --- | --- |
| P0-1: hs37d5 input versus hg38 assumptions | The published-input concern and absence of a documented local remap recipe are valid. The linked lionheart workflow does perform the claimed trimming, bwa mem alignment, sorting and duplicate marking. Its existence does not verify the final benchmark BAM. | Document the selected route and require final-input provenance/QC before calling this resolved. |
| P0-2: FrEIA `--both` | Confirmed in pinned `BaseCompAbCalc`. A data frame initialized outside the end loop retains the left end's index. The comparator currently sums potentially truncated `CountLen` values. | Exclude FrEIA concordance until the native aggregation is fixed. Optional timing is labeled as producing an incorrect aggregate. |
| P0-3: schedules | The schedules need substantially longer allocations. At the user's working estimates of 12 minutes for cfDNAlab and 120 minutes for a slow tool, two warmups plus five measurements at one core already consumes 15 h 24 min for those two methods alone. | Budget the complete method/core/repetition matrix, then split into method blocks when it exceeds the cluster limit. |
| P0-4: coverage concordance | Confirmed. The coverage runner deletes all outputs and there is neither a comparator nor a command README. | Add retained outputs, a streaming comparator and a settings/differences README. |
| P0-5: Picard version status | Confirmed by [Barclay's version return](https://github.com/broadinstitute/barclay/blob/5.0.0/src/main/java/org/broadinstitute/barclay/argparser/CommandLineArgumentParser.java#L151-L158) and [Picard's exit path](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/cmdline/CommandLineProgram.java#L191-L193). `--version` returns 1, and `install_tools.sh:127` and `lengths/benchmark_lengths.sh:122` treat that as fatal. | Handle that specific probe status and validate the reported version. |
| P0-6: Griffin environment | Confirmed. The conda-forge release API has [no Python 3.7.4 release](https://api.anaconda.org/release/conda-forge/python/3.7.4) and does have [linux-64 builds of 3.7.12](https://api.anaconda.org/release/conda-forge/python/3.7.12). | Pin Python 3.7.12, preserve the other relevant versions, and verify the full solve on Linux later. |
| P1-1: proper pairs | FinaleToolkit's proper-pair/read1 selection is confirmed. The quoted 12-16% is from sampled original alignments and cannot describe the remap. | Retain native and proper-pair diagnostics, measure their effect on the actual input, and distinguish their timing claims. |
| P1-2: fragment coverage | Confirmed that current coverage calls only use `--ignore-gap` and mosdepth has fragment mode. The proposed real-BAM equality claim omits CIGAR differences. | Add fragment coverage with exact controlled checks and real-input concordance, without promising equality. |
| P1-3: output costs and ordering | Native formats, interleaving and untimed `sync` are confirmed. A post-command filesystem flush is not an exact measurement of that tool's own disk I/O. | Keep native outputs, change block scheduling as requested, and report process time plus a separately labeled filesystem-flush measurement. |
| P1-4: duplicate marking | The selected remap includes MarkDuplicates. No final BAM/metrics were supplied to independently verify its execution. | Save actual duplicate metrics and keep FrEIA's additional coordinate deduplication visible. |
| P1-5: reconstructed mates | The current preparation code gives read1 positive TLEN regardless of strand. Remapping ordinary inward pairs recomputes it, so that original sign convention should not be reported as a Griffin limitation on the remap. | Correct the preparation convention and restrict synthetic pairing to recognized merged reads. Document that the final alignments are authoritative. |
| P1-6: memory | Confirmed. Summed process RSS can count shared pages repeatedly. Each invocation already has its own single-task Slurm step. | Use that step's Slurm memory accounting, check the site's accounting settings once, and retain GNU time separately. Remove the custom sampler from the timed path. |
| P1-7: cfDNAPro memory | Whole-file loading and lack of an explicit parallel reading path are confirmed. The numerical memory estimate and failure at 15x are unverified. | Pilot progressively, choose a common suitable node, and record resource failures as outcomes. |
| P1-8: FrEIA `-n 10` | Confirmed. The shipped config uses three extracted bases. Ten changes both N filtering and stored data. | Use `-n 3` and correct the associated workload description. |
| P1-9: FinaleToolkit end-motif bounds | Confirmed in [the end-motif reader call](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/frag/_end_motifs.py#L107-L120) and [alignment iterator](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/io/alignment.py#L242-L268). The requested 30-1000 bp filter is ignored. | Use a single 10-50,000 bp cfDNAlab diagnostic as a practical way to remove its usual length restriction for this dataset. No separate investigation of fragments above 50,000 bp is needed. |
| P1-10: length/orientation explanation | The distinction between start-order and 5'-end orientation is real. The review's Picard strict inequality is inaccurate at the boundary. samtools' integer halving truncation is confirmed. | Correct the exact rules and scope the bwa mem TLEN equivalence to ordinary inward alignments. |
| P1-11: midpoint sites | The user provides the site BED. The claim that no current group exercises reverse strands cannot be established without that input. | Use the supplied BED and record its provenance. Do not generate or select sites. |
| P1-12: provenance | Confirmed. Version strings alone do not record the checkout, resolved dependencies and biological inputs. | Save a short installation record, environment exports and input identifiers with the results. |
| P1-13: coverage interface | Confirmed. Methods, cores and repetitions are hard-coded, and timing settings/output capture lag behind the other runners. | Bring coverage under the same explicit run controls. |

### Corrections that materially change the proposed fixes

1. **mosdepth fragment mode does not reproduce cfDNAlab's treatment of deletions.** `src/commands/fcoverage/fcoverage.rs:1008` passes a D/N trigger of 1 to the segmented fragment iterator. `src/shared/fragment/segment_fragment.rs:230` constructs mapped segments, adds the inter-mate gap when requested, merges them, and clips them to `forward.pos` through `reverse.reference_end`. A D/N position unsupported by the other mate remains uncovered. [mosdepth v0.3.14](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L324-L334) instead fills its complete TLEN interval in fragment mode. Matching proper-pair flags cannot remove this difference.

2. **mosdepth's ordinary overlap correction also depends on proper pairs.** Its [overlap branch](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim#L269-L319) is conditional on `proper_pair`. Non-proper overlapping mates can therefore contribute twice there while cfDNAlab pairs and merges them. Add this to the read-depth comparison's expected differences.

3. **FrEIA concordance is disabled.** The pinned native aggregate can lose complete right-end motifs. The user explicitly rejected recovering counts from this damaged output. Keep optional native-workflow timing separate and omit its profile from concordance until the upstream aggregation is fixed.

4. **The wide end-motif diagnostic addresses the ignored 30-1000 bp filter.** Setting cfDNAlab to its supported 10-50,000 bp range effectively removes that usual restriction for this cfDNA dataset. Fragments above 50,000 bp are not a meaningful expected population here. My earlier emphasis on that formal upper bound was unnecessary. Record the settings and retain the real remaining differences: both-mate filtering, CIGAR rules, orientation and FinaleToolkit's window queries.

5. **Picard and cfDNAPro have slightly different exact orientation boundaries.** [htsjdk 5.0.0](https://github.com/samtools/htsjdk/blob/5.0.0/src/main/java/htsjdk/samtools/SamPairUtil.java#L101-L136) classifies FR with forward 5' position **<=** reverse 5' position, using the mate CIGAR when available. In cfDNAlab's coordinate terms this admits `forward.pos < reverse.reference_end`. [cfDNAPro](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/internalFunctions.R#L578-L605) uses a strict comparison in 1-based inclusive coordinates, equivalent to `forward.pos < reverse.reference_end - 1`. Their difference is outside the allowed >=10 bp benchmark range. cfDNAlab and samtools inward selection instead require forward/reverse start order, allowing coincident starts.

6. **Use Slurm step memory and keep flush timing explicit.** Each existing `srun --ntasks=1` invocation supplies an accounting unit for a tool and its children. Collect its completed `MaxRSS` through [sacct](https://slurm.schedmd.com/sacct.html). Check the accounting plugin, frequency and file-cache treatment once on the cluster. [Slurm's cgroup v2 accounting](https://slurm.schedmd.com/cgroup_v2.html#limitations) can include file cache and is not the same quantity as process RSS. There is no need to add a custom PSS monitor. Separately, GNU `sync -f` flushes the filesystem containing the argument, so describe its duration as filesystem flush time. [Coreutils implementation](https://github.com/coreutils/coreutils/blob/v9.8/src/sync.c#L127-L149).

7. **Keep installation state explicit.** Skip a tool when its completed version matches the requested version and its expected files exist. Reinstall it when the version changes, it is marked outdated, or `--force TOOL` is given. After same-version source changes, explicitly force or mark that tool outdated. Record the source commit for provenance, without building an automatic source-fingerprinting system.

## Implementation order

### 1. Repair installation and make setup depend on selected methods

Files: `benchmarks/install_tools.sh`, the four runners, and `benchmarks/README.md`.

- Capture Picard's probe output and exit status explicitly. Accept the documented version-probe status only when a recognizable version field equals 3.5.0. Reject missing/mismatched output and unrelated failures. Do not use an unrestricted `|| true` as the entire validation. Apply the same behavior to installer verification and runner setup. Preserve Java 17 and Picard 3.5.0.
- Replace Griffin's Python 3.7.4 pin with 3.7.12 and explicitly pin PyYAML after a Linux solve. Record the departure from the wiki. Package availability confirms the original blocker, not the solvability of the entire revised environment.
- Add `--override-channels` to `conda create` calls. For FrEIA's YAML, make a derived local copy adding `nodefaults` to its channels while retaining dependency pins and recording both files. This isolates the installation without rejecting a user's unrelated global Conda settings. [Conda documents this equivalent of channel override](https://docs.conda.io/projects/conda/en/stable/user-guide/tasks/manage-environments.html#creating-an-environment-file-manually).
- Update conda-pack to 0.9.2 after relocation verification. Its generated `conda-unpack` imports a sibling `conda_unpack_progress.py`. Preserve both files and ordinary script-directory imports. Keep packaging separate from timed tool environments. [Pinned source](https://github.com/conda/conda-pack/blob/0.9.2/conda_pack/core.py#L1094-L1097).
- Add FinaleToolkit's pinned runtime requirements from its [1.1.0 package metadata](https://pypi.org/pypi/FinaleToolkit/1.1.0/json), including the Python 3.11 SciPy requirement. Export the complete solved environment. The suggested six extra pins do not alone freeze numpy, pandas or every dependency. Do not describe version requirements as proof of upstream test coverage.
- Validate method selection before requiring paths or staging competitors. A cfDNAlab-only lengths run must not require Picard, cfDNAPro or Apptainer. Apply this principle across all four runners, including reference inputs only needed by selected methods.
- Reuse the existing samtools/Python tooling for header preparation where needed and document those setup dependencies. Removing the memory sampler removes its Python dependency. Do not add another packed helper environment just for benchmark monitoring.
- Replace coverage's silent `command -v` exits with clear missing-tool diagnostics for every selected executable, not just cfDNAlab.
- Keep the current standalone mosdepth/PanDepth routes unless node verification shows a concrete reason to change them. Do not reinstall or update unrelated competitors merely to enlarge this change.
- Add a node-side preflight using the staged binaries and environments. Check cfDNAlab/PanDepth linkage with `ldd` where relevant, and exercise cfDNAlab's default lengths/midpoints plotting on the small input. Fontconfig is loaded at runtime, so linkage inspection alone cannot verify fonts. Keep the same enabled outputs during verification and measurement.

Acceptance: each method's setup has only its documented tool/support dependencies, version failures are actionable, and a failed setup cannot be mistaken for a completed installation.

### 2. Add simple installation resumption and reinstall controls

Files: `benchmarks/install_tools.sh` and `benchmarks/README.md`.

- Put each tool/helper installation in a named step with a log. After verification and archive creation succeed, record the tool, pinned version/source revision, artifact paths and completion status in a small TSV.
- Resume by skipping completed entries whose version matches and whose expected files exist. A version change or `outdated` status causes that step to rerun. Add `--force TOOL[,TOOL]` and `--force all`, plus `--only TOOL[,TOOL]` for subsets. These controls are sufficient for this benchmark.
- Document that source or dependency changes without a version change require explicitly forcing or marking the tool outdated. Do not add automatic recipe hashes, recursive source fingerprints or checksum scans on every resume.
- Remove only that step's known partial files under the installation root. Preserve completed tools and logs. Write completion only after the entire step succeeds, keeping shell failures visible through functions and logging pipelines.
- Keep cfDNAlab's Cargo target directory under the installation root for incremental rebuilding. Record its actual build commit and dirty flag. The record describes the installed binary even if the checkout later changes.
- Save environment exports and pip package records once after installation. Keep the existing container digest/image identity. Copy this small provenance record alongside results, using the existing version, node, job and step fields in timing rows.
- Print `BENCH_*` exports only for completed tools. No new installation-state abstraction is needed beyond the TSV and per-tool functions.

Acceptance: an interrupted installation resumes, version changes and explicit outdated/force settings reinstall only the intended tools, and a failed step is never recorded as complete.

### 3. Document and validate the biological inputs

Files: new `benchmarks/data/README.md`, move `benchmarks/coverage/prepare_snyder.py` to `benchmarks/data/prepare_snyder.py`, and update references in existing benchmark documentation.

- Document the published BH01 source, original assembly, preparation arguments, seed, retention probability, merged-read reconstruction counts, and remapping route. Link the pinned [lionheart README](https://github.com/BesenbacherLab/lionheart/blob/085324a/remap/README.md) and [workflow](https://github.com/BesenbacherLab/lionheart/blob/085324a/remap/workflow.py). Capture the actual versions/commands used, since that workflow's installation example does not pin every dependency.
- Identify the exact GRCh38 no-alt reference and derived FASTA/FAI/2bit files by checksums. Equal contig names and lengths alone do not verify equal reference bases.
- Save final BAM/BAI checksums, header, flagstat, idxstats and MarkDuplicates metrics. Do not invent values or infer completion from the example filename. Describe 15x as a target until measured on the remapped input, and specify whether depth means sequenced-base depth or fragment coverage.
- Measure proper-pair effects using the actual native/proper cfDNAlab outputs under the benchmark filters. flagstat gives useful record-level QC but not the conditional pair inclusion fraction after all benchmark filters. Do not reuse the published-BAM 12-16% figure as a remap result.
- Restrict reconstruction to recognized `M_` merged reads. Fail with a count and example for unexpected retained unpaired names instead of turning an arbitrary singleton into a complete fragment. That is a biological correctness issue, not just defensive input handling.
- Give the reconstructed forward mate positive TLEN and the reverse mate negative TLEN for mapped spans in the allowed domain. Preserve unmapped behavior and existing filtering flags. Explain that synthetic mates reuse the merged sequence and do not recover original separate reads/qualities. The later alignment recomputes TLEN and proper-pair status.
- Remove the midpoint README's reconstruction-sign explanation as an expected difference for the remapped input. Retain accurate native Griffin filtering differences.
- Fail early when Griffin's selected chromosome names are entirely numeric, since its pandas type inference conflicts with string chromosome arguments. Validate the selected BAM/reference contigs and window bounds before timing.
- Use the user's supplied midpoint BED through `BENCH_INTERVALS`. Preserve its coordinates, groups, strands and repeated/overlapping rows. Record its checksum. Do not generate sites, select subsets or prescribe biological sources, site counts, strands or a window size.
- Reject invalid explicit BED5/BED6 strand tokens with a row number. Do not silently convert them to `.`. Preserve explicitly unstranded rows and the existing duplicate/overlap policy.

Acceptance: the same biological inputs can be identified and reconstructed from recorded commands, and no final-input property is claimed merely because the workflow is capable of producing it.

### 4. Repair end-motif comparison and define comparable configurations

Files: `benchmarks/ends/benchmark_end_motifs.sh`, `compare_end_motifs.py`, `README.md`, and the corresponding length/coverage configuration descriptions.

FrEIA exclusion, updated by the user's implementation instruction:

- Use `-n 3` for optional native-workflow timing, retaining its native filtering and exclusive maximum.
- Do not retain or recover the damaged both-end aggregate for concordance. Explicitly skip any old FrEIA profile in the comparator.
- Label the known incorrect aggregate and exclude it from equivalent-output speed claims until fixed upstream. Keep the issue draft without submitting it.

Other configurations:

- Keep normal cfDNAlab configurations and the existing proper-pair diagnostics. Add proper-pair regional and fragment coverage configurations. Describe them as closer filter alignment, not identical populations.
- Add the 10-50,000 bp proper-pair/reference diagnostic to remove the usual cfDNAlab length restriction for concordance. Record its effective bounds and FinaleToolkit's ignored bounds. No extra out-of-range census or >50,000 bp fixture is needed.
- Keep optional 3-mer FrEIA timing separate from 4-mer FinaleToolkit concordance and the read/reference/clipping diagnostics. Avoid conflating autosome-only outputs with all-contig results.
- Report default and closer-filter concordance separately. If a speed ratio is described as a matched-filter comparison, time that configuration too. Timings of cfDNAlab's normal configuration remain useful but must not inherit the matched label solely because a different configuration was checked for concordance.
- Correct the lengths README with the exact orientation rules, samtools' floor-after-halving behavior, and bwa mem's 5'-end TLEN construction for ordinary inward pairs. Do not promise equivalence for arbitrary aligners, abnormal pairs or stale TLEN values.
- Cite the actual bwa version's TLEN calculation. The independently inspected [bwa 0.7.19 implementation](https://github.com/lh3/bwa/blob/v0.7.19/bwamem.c#L887-L892) uses the forward `pos` and reverse `reference_end - 1` for ordinary inward pairs. The preparation manifest must identify which version produced the benchmark BAM.
- Add Spearman correlation and a clearly defined high-frequency log-ratio diagnostic. Specify whether the >=1e-4 cutoff applies to either or both profiles. A qualifying motif with zero frequency on the other side should have an explicit infinite/undefined result, rather than an undocumented pseudocount.
- If FinaleToolkit is selected, guard against contigs whose lengths are exact multiples of 1,000,000. Its pinned window builder omits the last full window in that case. Name the upstream limitation in the error rather than presenting the reference as biologically invalid. Preserve the known boundary-crossing double-count limitation.

Acceptance: FrEIA is excluded, damaged tables cannot pass as complete, actual effective length ranges are visible, and configuration labels support only the comparisons actually performed.

### 5. Add fragment coverage and a bounded-memory concordance analysis

Files: `benchmarks/coverage/benchmark_coverage.sh`, new `compare_coverage.py`, and new `README.md`.

- Add `cfdnalab-fragment`, `cfdnalab-fragment-proper`, `mosdepth-fragment` and `cfdnalab-region-counts-proper`. Retain clearly named read-depth, regional count and autosome diagnostics. mosdepth fragment mode uses `--fragment-mode` with the selected flag/MAPQ/length settings.
- Retain a native output for each concordance configuration from the first successful measured run at the capture core budget. Do the same in other runners, replacing their present capture from repetition 1, which is a warmup under production defaults. Preserve the shared BED and scope metadata beside outputs.
- Stream comparisons over run boundaries, not arrays with a value per genomic base. Normalize contig order using a saved header/manifest, since outputs can order contigs differently. Validate sorted, non-overlapping, in-bounds intervals and supported finite depths.
- Interpret omitted runs as zero only within a declared processed scope. Distinguish zero coverage from an absent/incomplete track. Split or merge equivalent adjacent runs without treating tile boundaries as differences. Read all concatenated zstd frames.
- Weight identical-base fractions, signed difference distributions, MAE and RMSE by segment length. Also report the nonzero union separately so uncovered reference bases do not dominate apparent agreement. Compute 10 kb-bin correlation from integrated depths, with defined treatment of shorter final bins and constant profiles.
- Compare regional fragment counts by chromosome/start/end, with a verified row identity. Report count totals, differences, MAE/RMSE and correlation. The ratio of regional totals is not a unique-fragment ratio, since a boundary-crossing fragment can contribute to adjacent windows.
- Retain compressed cfDNAlab/mosdepth tracks for full-scope concordance. Extract a declared chr21/chr22 subset of samtools/PanDepth output after timing, unless pilot storage permits preserving the full tables. Record the subset scope in every comparison. Do not equate different scopes in summary tables.
- Document the remaining pair/MAPQ/flags/orientation/TLEN/CIGAR differences. In particular, include D/N treatment for fragment mode, proper-pair-dependent overlap removal in ordinary mosdepth, and FinaleToolkit regional boundary misses from fetching read1 alignments rather than every overlapping fragment.
- Make mosdepth the required added fragment comparator. deepTools fragment coverage can be a later optional addition using read1+proper-pair inclusion (`66`) and verified 4.0 CLI settings. It also retains selection/CIGAR differences. Do not expand this repair into adding Qualimap, bedtools, megadepth or Triton without a specific scientific comparison to support.

Acceptance: controlled equivalent cases give exact counts, native biological differences remain visible on real data, and whole-genome comparison does not require whole-genome dense arrays.

### 6. Make timing, memory and pilot scheduling consistent

Files: all runners, `benchmarks/README.md`, and removal of the unused `benchmarks/measure_memory.py` once Slurm accounting replaces it.

With the user's working runtime estimates, the minimum production budget for just two methods at one core is:

| Method | Time per invocation | Two warmups plus five measured runs |
| --- | --- | --- |
| cfDNAlab | 12 min | 1 h 24 min |
| Slow comparator | 120 min | 14 h |
| Total, before other core budgets or methods | | 15 h 24 min |

Under the old three-warmup/ten-measurement schedule, those same two one-core blocks would take 28 h 36 min. A serial comparator repeated at both 1 and 32 cores would take about 28 h under the proposed seven-invocation schedule. Therefore 12 h cannot be the allocation for the full experiment. Pilots refine the larger budget and tell us how to split it.

- Add validated methods/cores/warmup/repetition controls to coverage. Use two warmups and five measured repetitions as the production default. Keep the current `1 4 8 12 32` budgets.
- Add `BENCH_PILOT=1` consistently, selecting zero warmups and one measured repetition while preserving requested methods/core budgets. Mark pilot results explicitly. Distinguish a small wiring check, a full-input resource pilot, and a concordance/configuration pilot.
- Document a reproducible small input made by read-name-hash sampling, for example `samtools view -b -s 42.01` followed by indexing, plus a small midpoint site file. Both mates receive the same sampling decision. Reduce the fraction further if 1% is still too large for a quick wiring check. Never infer final-input runtime from this small input alone.
- Run contiguous method blocks on the same physical node, with each method's warmups immediately before its measurements. Shuffle its requested core budgets within each repetition and retain the actual order in the existing TSV. Extra anchor runs are optional follow-up only if the pilot shows unexplained drift.
- Run cfDNAPro and FinaleToolkit's length histogram at the requested 1- and 32-core endpoints, since their inspected counting paths are serial. Keep Picard on the selected matrix because its samtools stage and JVM can use more than one core.
- Sum the pilot duration of each requested method/core combination, multiply by the scheduled warmup-plus-measurement count, and allow about 20-30% plus staging/cleanup time. Use that calculation for each feature and for the entire experiment. Do not apply a universal 12 h limit, including to a full-matrix pilot that itself contains many long runs.
- Split oversized series by method, or by core budget if a method block is still too long, using the existing method/core environment variables. Pin these sequential jobs to the same physical node with `--nodelist` and check the hostname. Each separately submitted block needs its own warmup budget. No automatic job-planning framework is needed.
- Use the existing single-task `srun` step as the memory-measurement unit. Keep staging, profile copying and concordance outside that step so they do not dominate its accounting peak. Save the existing `job_id` and `step_id` for every invocation.
- After a job finishes, collect completed step accounting once, for example `sacct -j JOB_ID --parsable2 --units=K --format=JobIDRaw,NTasks,MaxRSS,State,ExitCode`, and join by job/step ID during summary generation. Keep the raw accounting export. Do not poll Slurm inside timed commands. If database records are not yet available, leave values `NA` and repeat this collection later.
- Use Slurm's reported step peak as the primary memory result, retaining GNU time's `max_rss_kib` separately. In the pilot, check the site's accounting plugin, frequency, descendant accounting and file-cache policy once. With cgroup v2, Slurm's `MaxRSS` can include file cache. Label it accurately and use the same policy for every tool. If the site's accounting cannot supply a useful measurement, report that limitation before adding another measurement system.
- Remove the custom process-tree sampler from the timed path. There is no requirement for sampled mean RSS/PSS or for a new monitoring runtime.
- Measure flush duration immediately after command completion, before output preservation/cleanup. Flush staging writes before the first series. Record process wall time, flush time and their sum, plus output bytes from `du`. Check scratch-space requirements during the pilot.
- Use the combined process-plus-filesystem-flush time as the proposed main end-to-end number under D3, while retaining process time. Call it a filesystem-flush boundary, not exact per-tool durable-write latency. Treat a flush failure as a failed run. Keep CPU utilization based on process CPU time divided by process wall time, not the flush-inclusive denominator.
- Keep the existing node/core-placement checks and record the relevant Slurm configuration once. Per-run cache-residency and CPU-state monitoring are not required. Investigate them only if the pilot shows a concrete problem.
- Retain native per-position output generation and compression costs. Label samtools full QC, the samtools/Picard pipeline, FrEIA's intermediates and Griffin's compulsory reference/GC/bigWig work. Show observed scaling, but do not interpret output-dominated curves as pure counting-kernel scaling.
- Pilot cfDNAPro first on about 0.1x, then increase input size if needed to choose a node with demonstrated headroom. Use that node for every compared method.
- Preserve immediate failure rows, correct TSV column counts and actionable messages. Retain failed logs in results before scratch cleanup, including interrupted/OOM outcomes where possible. Report failed resource configurations and node memory instead of silently omitting them. Missing measurements stay `NA`. Continue or resume remaining methods explicitly rather than interpreting absent rows as success.
- Keep the existing TSV format straightforward. Add the needed measurements and preserve correct columns on success and failure. A new generic TSV-writing layer is unnecessary.

Acceptance: every published timing has a documented schedule, comparable measurement units and status, and pilot estimates account for the work actually scheduled.

### 7. Finish reporting and targeted documentation

- Summarize the timing TSV and Slurm accounting with sample size, mean/SD, median/IQR and the measured points. Exclude warmups and failed runs. Keep bootstrap intervals and additional plotting/statistics optional until the basic results are available.
- Update all command examples and result placeholders to the same production defaults. Remove arbitrary production wall-time promises and point to pilot estimates.
- Fix the Griffin padding statement to the sufficient condition `floor(max_length / 2) <= min_length`. Fetch losses become possible beyond that condition, from `max_length >= 2 * min_length + 2`, depending on read alignment extent.
- Correct the inconsistent verification date when documenting the fresh audit. Keep a distinction between source inspection, package availability, successful installation and actual run verification.
- Preserve the existing valid native flags, output readers, inward motif orientation, exclusive FrEIA maximum, length-axis boundaries, native midpoint normalization and unmodified competitor sources unless a specific later check contradicts them.
- Keep the current cfDNAFE exclusion. Its additional CLI keyword issue is not needed to justify the already identified motif-definition mismatch, and this review does not endorse the claim that every upstream branch/fork has been audited.
- Defer a 1x/5x/15x depth series until the corrected primary comparisons and resource pilots succeed. It is useful for a scaling-with-input-size claim but is additional experimental scope, not a prerequisite for a clearly scoped fixed-input benchmark. Choose any fastest competitor from measured comparable results.

## Verification to prepare during implementation

Do not run tests as part of this plan or its implementation without a new instruction overriding the repository's prohibition. Prepare deterministic regression fixtures with expectations derived by hand, and leave execution to the user. Static shell/Python syntax checks can be reported separately. Run the required Cargo checks if implementation later changes Rust code or Rust tests. The implementation also runs the repository-required Cargo checks. Tests remain unexecuted.

Meaningful regression cases to prepare:

| Area | Hand-derived expectation |
| --- | --- |
| Picard probe | Recognizable 3.5.0 with the known status 1 passes. Wrong version, empty output or an unrelated failure fails. A cfDNAlab-only selection never probes Picard. |
| Installation state | A completed matching version with its expected files skips. A changed version, missing archive, outdated status or explicit force reinstalls that tool. Same-version source changes are handled by marking outdated or forcing. A failed subcommand never creates a completion record. |
| FrEIA exclusion | Existing FrEIA directories are explicitly skipped. Native incorrect profiles are not retained or recovered. |
| Coverage gap and deletion | Forward `pos=100`, `30M`, and reverse `pos=170`, `30M` give a 100 bp directional span, including a 40 bp gap. `--ignore-gap` covers 60 bases. Replacing forward CIGAR with `15M2D15M` retains that span but leaves positions `[115,117)` uncovered in cfDNAlab fragment coverage. mosdepth fragment mode covers all 100 bases. |
| Coverage comparison | Adjacent/tile-split runs and omitted declared zero runs compare identically. A depth difference of 1 over 10 of 100 evaluated bases gives MAE 0.1, RMSE `sqrt(0.1)` and identical fraction 0.9. Reject overlapping/out-of-bounds runs and unintended scope differences. |
| Length/orientation | Ordinary inward spans match bwa mem TLEN. Reverse-start-first overlapping pairs expose the known selection difference. An odd number of passing samtools alignment records in a length bin is halved and truncated. All fragment fixtures are >=10 bp and normal benchmark fixtures respect the 30 bp inclusion minimum. |
| Motif range/window | A valid 20 bp pair demonstrates that the normal 30 bp cfDNAlab bound differs from the wide diagnostic. A contig exactly 1 Mb long triggers the selected FinaleToolkit limitation guard. |
| Midpoint input | Equal even windows on valid `+`, `-` and `.` rows retain their intended orientation and row counts. Invalid explicit strand and all-numeric Griffin chromosome input fail before timing. |
| Timing/memory | Pilot and production matrices have the intended invocation counts and estimated total hours. Success, placement failure, missing output and flush failure produce valid status rows. Saved Slurm records join to the correct job/step, missing accounting stays `NA`, and the existing GNU time metric remains separately labeled. |
| Input preparation | Only retained recognized merged names are reconstructed. Forward/reverse TLEN signs are opposite in the intended direction. Unknown singletons fail clearly, and unmapped/secondary/supplementary/duplicate/QC-fail flags remain intact. |

User-run Linux verification after the code changes should cover revised Conda solves, interrupted installer resumption, relocated environments, runtime linkage/fonts and actual CLI output generation on a tiny input. Then run a concordance pilot, a full-input resource pilot covering the production matrix, and production measurements on the recorded node. No publication equality, memory-capacity or wall-time conclusion should be asserted before the corresponding measurements exist.

## Implementation record (2026-10-30)

- Installer resumption, selected-tool setup, Picard/Griffin fixes and explicit dependency records are implemented.
- All runners use two warmups/five measurements, contiguous method blocks, pilot controls, separately timed filesystem flushes, output sizes and first-measured profile capture outside the tool step.
- The process-tree sampler is removed. A standalone summary joins sacct records and budgets the observed pilot matrix.
- Coverage has fragment/proper-pair configurations, explicit retained scope and a streaming comparator. End motifs include the wide diagnostic, added frequency metrics and FrEIA exclusion.
- BAM preparation is under `benchmarks/data/`, with corrected TLEN signs, merged-read validation and provenance instructions. Midpoint sites remain user-provided input with no preparation helper.
- Deterministic regression fixtures are prepared for user execution. Static checks do not establish Linux runtime correctness.

## Remaining inputs before publication runs

- Actual BAM/reference/site files and their provenance/QC records.
- The user's midpoint site BED and its provenance.
- A cluster node, local filesystem and wall-time limit validated by the resource pilot.
- Whether optional deepTools coverage and a depth series support claims the paper actually intends to make.

These do not block implementing the benchmark repairs. They do block calling the final experimental design and measurements complete.
