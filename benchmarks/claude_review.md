# Review of `benchmarks/`

Reviewed on 2026-10-29 at commit `922f717c` (branch `benchmarking`). I read every file in `benchmarks/` in full. Tool behaviour was checked against upstream source at the pinned versions:

| Tool | Version read |
| --- | --- |
| Griffin | v0.2.0 |
| FinaleToolkit | v1.1.0 |
| FrEIA | `2c32fb9` |
| cfDNAFE | `90c6979` |
| cfDNAPro | `ab9eeb1` |
| samtools/HTSlib | 1.25 |
| Picard / htsjdk | 3.5.0 / 5.0.0 |
| mosdepth | v0.3.14 |
| PanDepth | v2.26 |
| conda-pack | main |

cfDNAlab and `py-cfdnalab` were checked against this repository's code. **No tool was installed or run.** Three exceptions:
- Conda solvability was checked with `conda create --dry-run` solves, which only run the solver.
- Release binaries were inspected with `file` and `objdump`.
- A 20 MB range sample of the public `BH01.bam` (the header, the index and three 3 MB chunks) was parsed to check record types.

`file:line` citations refer to the pinned upstream version unless the path starts with `benchmarks/`, `src/` or `py-cfdnalab/`. Claims marked *estimate* or *inferred* were not confirmed from source or data.

Severity:
- **P0**: the run fails, or a reported number would be wrong.
- **P1**: a reviewer could reasonably reject the comparison or its interpretation.
- **P2**: fix before the final runs.
- **P3**: wording or minor.

## Summary

All cfDNAlab calls are valid and compute what the READMEs say. The Python reader usage in the comparators is correct, as are the Griffin calls, the FinaleToolkit, samtools and Picard settings, and the conda-pack relocation. The section "Verified as correct" lists these so they are not changed by mistake.

The benchmark is still not ready to run. There are four problems:

1. **The installer and the lengths runner stop at once** (P0-5, P0-6).
   - `picard … --version` exits with status 1. Under `set -e`, this kills `install_tools.sh` before FrEIA, Griffin and the helper environments are created. It also kills every lengths job during setup, even when Picard is not selected.
   - Separately, `python=3.7.4` does not exist on conda-forge, so the Griffin environment cannot be created.
2. **The FrEIA comparison would be wrong** (P0-2). FrEIA's `--both` aggregation misaligns the right-end rows, and the comparator sums exactly those rows. This is an upstream bug; an issue text is in `benchmarks/freia_issue.md`.
3. **The schedules cannot finish** (P0-3). For example, coverage needs 455 runs within 12 h, an average of 95 s per run, while samtools alone writes about 50 GB of text per run.
4. **Coverage has no concordance analysis** (P0-4). It also never tests cfDNAlab's main quantity, fragment coverage, against `mosdepth --fragment-mode` (P1-2).

The hg38 input itself is settled (P0-1): `prepare_snyder.py`, then the lionheart remap workflow with bwa mem and Picard MarkDuplicates. This resolves P1-4 and P1-5. The benchmark only needs to document and link those steps.

The largest scientific issue is proper pairs (P1-1). In a sample of the published BH01 BAM, 12–16% of pairs are not flagged proper; remeasure on the remap. Every FinaleToolkit comparison therefore needs the matching cfDNAlab `--require-proper-pair` configuration as the reported comparison, alongside cfDNAlab's default.

## Decisions (answered by Ludvig on 2026-10-29)

| # | Decision | Outcome |
| --- | --- | --- |
| D1 | Input data route (P0-1) | **Done:** `prepare_snyder.py`, then the lionheart remap workflow to hg38, with duplicates marked. Codex: link it as described in P0-1. |
| D2 | Repetitions and core budgets (P0-3) | **5 measured repetitions.** Set wall times from a pilot that covers every core budget (P0-3). Run tools with no parallel code path (cfDNAPro, FinaleToolkit `frag-length-bins`) only at 1 and 32 cores. Run **everything on the same node**, in contiguous per-method blocks, with warmups immediately before each method's measured runs (P1-3). |
| D3 | samtools `depth -aa` and PanDepth `-a` in per-base timing (P1-3) | **Keep them**; their output cost is their own. Per-method blocks keep it from affecting other tools. Time output flushing for every tool (P1-3), and report bytes written. |
| D4 | cfDNAPro memory (P1-7) | Pilot on a 0.1x downsample. Then choose a node type with enough memory for every tool, for example ≥400 GB, and run **all** tools on that node type. |
| D5 | Extra comparators (P1-2) | **Add `mosdepth --fragment-mode`** with an equality check. Use deepTools only for fragment coverage (P1-2); it does not compute the length, end-motif or midpoint quantities exactly. |

## P0: blockers

### P0-1. The input BAM is GRCh37/hs37d5 with numeric contig names, but everything assumes hg38

**Status: resolved** (D1). The hg38 input was produced in three steps:

1. **`benchmarks/coverage/prepare_snyder.py`** on the published BAM. This downsamples by read name and reconstructs about 3.5 million pairs from merged `M_` reads. That count is consistent with merged reads being about 1.8% of fragments at 15x.
2. **The lionheart remap workflow** (`BesenbacherLab/lionheart`, `remap/`; latest commit touching it is `085324a`, 2026-03-24):
   - `samtools collate | samtools fastq` (reverse-flagged reads are written reverse-complemented, so each reconstructed pair becomes an ordinary fully overlapping read pair);
   - `guessadapt` + `cutadapt` (`-q 20`, `-m 20`);
   - `bwa mem -p -Y -K 100000000` against the GRCh38 no-alt analysis set with UCSC names;
   - `samtools sort`;
   - Picard `MarkDuplicates`.
3. **Indexing** with BAI.

**Consequences, all verified in source:**
- bwa mem computes TLEN from the mates' 5' ends (`bwamem.c:888-891`), so TLEN is recomputed during the remap.
  - For inward pairs, `|TLEN|` equals cfDNAlab's `reverse.reference_end - forward.pos`. The TLEN-based tools' length definition is therefore identical to cfDNAlab's (P1-10).
  - The forward mate always gets the positive sign, including for coincident mates, so P1-5 is moot.
- Duplicates are marked, so P1-4 is resolved.
- The `chr` names fix the Griffin problem described below.

**What remains for Codex**
- Write `benchmarks/data/README.md`. It should cite both steps with the lionheart commit, the `prepare_snyder.py` arguments (keep probability, seed) and the reference files.
- Move `prepare_snyder.py` there, since it is input preparation, not part of the coverage benchmark.
- Record the md5, `samtools flagstat` and `idxstats` of the final BAM.
- Add a guard that fails early if contig names are all numeric and Griffin is selected.

The original finding is kept below for the record.

**Evidence**
- The Snyder README (`shendurelab/cfDNA` README.md:34-36) says reads were aligned with BWA aln/sampe to GRCh37, 1000 Genomes phase 2 (hs37d5).
- The `BH01.bam` header has 86 references: `1..22, X, Y, MT, GL…, NC_007605, hs37d5`. Download: <https://kircherlab.bihealth.org/download/cfDNA/BH01.bam> (217 GB).
- The repository assumes hg38 with UCSC names in `benchmarks/ends/README.md:3,48`, `benchmarks/lengths/README.md:30,91,98,122`, `benchmarks/midpoints/README.md:40` and `benchmarks/lengths/benchmark_lengths.sh:36`.
- `benchmarks/coverage/prepare_snyder.py` neither realigns nor renames contigs. No script produces the `BH01.hg38.paired.15x.bam` the READMEs use.

**Consequences if the prepared BAM keeps hs37d5 names**
- **Griffin finds no sites.**
  - `griffin_functions.py:9` reads the site table with `pd.read_csv`, so all-numeric `Chrom` values become integers.
  - `griffin_functions.py:14` keeps `Chrom.isin(chroms)`, and the `--chroms` values are strings, so every site is removed.
  - `griffin_coverage.py:471-472` filters the header-less chrom sizes table the same way.
  - The midpoint benchmark would fail or produce empty profiles. This follows from the source and pandas dtype inference; it was not run.
- **cfDNAPro runs with the wrong genome label.**
  - `genome_label="hg38"` works only by coincidence: the only change is renaming `MT` to `chrM`, and both are 16569 bp (`R/readBam.R:312-335`, `R/genomeReference.R:4-42`).
  - The label is still wrong for the paper.
  - `hg19` would stop with a length conflict, since hg19 `chrM` is 16571 bp. This is inferred from the code.
- **Ends and midpoints need hs37d5 references.** They would need an hs37d5 `.2bit` and an hs37d5 FASTA/`.fai`, not hg38 ones.

**Provenance:** see the status above. State in the Methods that the alignments differ from the published BH01 BAM that FinaleToolkit's paper benchmarked.

### P0-2. FrEIA's `--both` output misaligns the right-end rows, and the comparator sums them

**Evidence (verified in `2_data_transformation.py` at `2c32fb9`)**
- `BaseCompAbCalc` creates `BCDf` once, outside the per-end loop (lines 121-124).
- For the left end, `BCDf.CountLen = hBCDf.strand` assigns a Series to an empty frame, which gives `BCDf` the left end's index.
- For the right end, the same assignments are aligned to that left-end index (pandas 1.1.5, `frame.py:3722-3731`):
  - When the right end has more (length, motif) rows, its last rows (the longest lengths) are silently dropped.
  - When it has fewer, the extra rows get `NaN` in `Base`, and `benchmarks/ends/compare_end_motifs.py:53` crashes on `.translate`.
- The authors' own rule never uses `--both` (`workflow/rules/2_FrEIA.smk:103-110`).
- **This is an upstream bug.** `--both` is a documented option of FrEIA's script, and using it is necessary to count both ends. The comparator's error is not guarding against the bug; `-n 10` (P1-8) is a separate, local misconfiguration.
- The upstream fix is to create `BCDf` inside the `for e in whichEnd` loop. File an upstream issue and cite it in the Methods next to the workaround. Do not rely on a fix arriving, since the repository has had no commits since 2024-04.

**Fix in `benchmarks/ends/compare_end_motifs.py:43-60`**
- Drop rows where `Base` is null.
- Take **`CountSamp` once per (`WhichEnd`, `Base`)**. It is the total per motif, computed on the full table before the misaligned assignment, so it is complete for every motif that keeps at least one row.
- Assert that the left and right totals are equal. Each retained fragment contributes one end of each.
- Assert that all 64 motifs are present for both ends.
- Describe the upstream bug in `benchmarks/ends/README.md` (FrEIA section).

### P0-3. The schedules cannot finish in the example wall times, and a failed run stops the job

| Runner | Runs (methods × cores × 13 repetitions) | Example `--time` | Average allowed per run |
| --- | --- | --- | --- |
| coverage | 7 × 5 × 13 = 455 | 12 h (`benchmark_coverage.sh:40`) | 95 s |
| lengths | 6 × 5 × 13 = 390 | 48 h (`lengths/README.md:119`) | 7.4 min |
| ends | 3 × 5 × 13 = 195 | 24 h (`ends/README.md:92`) | 7.4 min |
| midpoints | 2 × 5 × 13 = 130 | 24 h (`midpoints/README.md:62`) | 11 min |

These budgets are not realistic:
- `samtools depth -aa` writes about 3.1 × 10⁹ text lines (43–55 GB) per run (P1-3).
- The cfDNAPro, FinaleToolkit `frag-length-bins`, Picard collector, FrEIA extraction and Griffin aggregation steps are single-process Python, R or Java passes over a 15x BAM.
- The jobs use `--no-requeue` and exit at the first failed step.

**Fix**
1. **Pilot timing over every core budget.**
   - Run each method once at **every** core budget on the final BAM: `BENCH_WARMUPS=0 BENCH_REPETITIONS=1 BENCH_CORES="1 4 8 12 32"`. The README examples use `BENCH_CORES="1"`, which cannot predict the rest.
   - The pilot is exactly one repetition block, so a job with `W` warmups and `R` measured repetitions needs about `(W + R) × pilot wall time × 1.2`:
     - 1 + 5 → about 7.2 × the pilot;
     - 1 + 10 → about 13.2 × the pilot;
     - 3 + 10 → about 15.6 × the pilot.
   - Add `benchmarks/estimate_walltime.py`. It should read a pilot `timings.tsv` and print this per method, plus the per-method `--time` for one-job-per-method submission.
   - If a job exceeds the cluster's wall-time limit, split by method.
   - Also record scratch use and peak memory from the pilot.
2. **Same node for everything.** If the whole series exceeds the cluster's wall-time limit, split it into several jobs pinned to the **same physical node** (`--nodelist=<node>`), one or more methods per job, and verify the hostname in every row (already recorded). Never spread methods over different nodes.
3. **Warmups and repetitions.**
   - Use `BENCH_WARMUPS=1`. The page cache is already warmed by the staging copy; the warmup mainly absorbs first-run effects.
   - Use `BENCH_REPETITIONS=5`, or 3 for methods whose runs exceed 30 min.
   - If runs on an exclusive node vary by a few percent (*assumption*; check it in the pilot or the first block), the standard error with 5 repetitions is about 1%. That is far below the differences between tools that are likely to be reported.
   - Use 10 repetitions only for methods that are cheap to repeat. Check the variance in the pilot or the first block.
4. **Core budgets for tools without a parallel code path.** cfDNAPro (no BiocParallel, one `readGAlignmentPairs` call, `R/readBam.R:247`) and FinaleToolkit `frag-length-bins` (fixed single reader thread, `io/alignment.py:107`): run only at 1 core and at 32 cores. This is enough to show the absence of scaling, and the README should say so. Keep Picard on all budgets, because its samtools input stage uses `-@`.
5. Add the coverage runner's missing `BENCH_METHODS`, `BENCH_CORES`, `BENCH_WARMUPS` and `BENCH_REPETITIONS` overrides (P1-13) so this is possible there too.

### P0-4. The coverage benchmark has no concordance analysis and no README

**Evidence**
- `benchmarks/coverage/benchmark_coverage.sh:18-19,380` deletes every output.
- There is no `compare_coverage.py` and no `coverage/README.md`.
- Output concordance is one of the two main measures in the application note.

**Fix**
1. **Keep outputs from repetition 1 at the first core budget**, as the other runners do.
   - Skip retaining the samtools and PanDepth per-position tables because of their size, or keep only chr21 and chr22 for them.
   - mosdepth's per-base output and cfDNAlab's bedgraph are run-length encoded and small enough to keep.
2. **Write `benchmarks/coverage/compare_coverage.py`.**
   - **Per-base:** stream-merge the run-length tracks and report:
     - the fraction of bases with identical depth;
     - MAE and RMSE;
     - the distribution of differences;
     - Pearson r over 10 kb bins.
   - **Inputs to handle:** cfDNAlab omits zero runs unless `--keep-zero-runs` is set (`src/commands/fcoverage/writers.rs:1118-1122`), and splits rows at 10 Mb tile edges. Its `.zst` file is several concatenated frames, so read it with `zstd -dc`, or a reader that handles multiple frames.
   - **Regional:** compare `cfdnalab-region-counts` with `finaletoolkit-coverage`: r, MAE, RMSE, ratio of totals, and per-window differences.
3. **Write `benchmarks/coverage/README.md`** with a settings and expected-differences table like the other commands.

### P0-5. Picard's `--version` exits with status 1, which stops the installer and every lengths job

**Evidence (Picard 3.5.0 source; not run)**
- Barclay 5.0.0 `CommandLineArgumentParser.parseArguments` prints the version to stderr and returns `false` (lines 156-158).
- Picard's `parseArgs` passes `false` through (`CommandLineProgram.java:319-335`), and `instanceMain` returns 1 (lines 209-210).
- `PicardCommandLine.main` calls `System.exit` with that value (line 115). The Bioconda launcher `exec`s Java, so the status reaches the shell.

**Consequences**
- `benchmarks/install_tools.sh:127` runs inside a subshell that inherits `set -e`. The subshell exits, the parent `set -e` stops, and the Picard archive, FrEIA, Griffin, the comparison environment and the Snyder preparation environment are never created.
- `benchmarks/lengths/benchmark_lengths.sh:122` captures the version with `$(… --version 2>&1 | tr …)` under `pipefail`. Every lengths job therefore exits during setup, whichever methods are selected.

**Fix**
- Capture the output, tolerate the status, and verify the content, for example `version=$(picard CollectInsertSizeMetrics --version 2>&1 || true)`, then fail unless it contains `3.5.0`.
- Separately, only stage and version-check the tools whose methods are selected (the ends runner already does this for FrEIA). The lengths runner currently requires every archive (`benchmark_lengths.sh:21-25`).

### P0-6. The Griffin environment cannot be created

**Evidence**
- `benchmarks/install_tools.sh:163` pins `python=3.7.4`. conda-forge has no 3.7.4 build in any subdir; its 3.7 builds are 3.7.0–3.7.3, 3.7.5, 3.7.6, 3.7.8, 3.7.9, 3.7.10 and 3.7.12.
- A `conda create --dry-run` (linux-64, strict priority) fails with `PackagesNotFoundError: python=3.7.4`.
- With `defaults` present, the build there is rejected by strict priority. The wiki's 3.7.4 came from `defaults`.

**Fix**
- Use `python=3.7` (resolves to 3.7.12). Keep every other wiki pin; with `pyyaml=5.3.1` added, the dry-run solves.
- Document in `benchmarks/README.md` that the patch version differs from the wiki because conda-forge never built 3.7.4.

## P1: scientific validity

### P1-1. Proper-pair selection removes 12–16% of BH01 pairs in FinaleToolkit

**Evidence**
- FinaleToolkit's reader requires flag 0x2 and uses read1, its TLEN and read1's MAPQ only (`io/alignment.py:60-71,248-260`). This applies to `coverage`, `end-motifs` and `frag-length-bins`.
- In the BH01 sample, 12–16% of mapped, same-chromosome, opposite-strand paired records were not proper, mostly with |TLEN| 201–1000. BWA sampe flags pairs as proper from its estimated insert-size distribution.

**Consequences**
- FinaleToolkit's length histogram is depleted in the dinucleosome range.
- Its regional counts are lower.
- Its end motifs come from a shorter fragment population.
- None of this reflects how the tools count. It reflects the aligner's proper-pair flag.

**Fix** (agreed: run both with and without the proper-pair requirement)
- **Coverage:** add `cfdnalab-region-counts-proper` (the same call plus `--require-proper-pair`) and `cfdnalab-fragment-proper` (P1-2).
- **Lengths and ends:** keep `cfdnalab-proper` and `cfdnalab-reference-proper`.
- **No extra timing cost.** The proper-pair variants are needed for concordance only, so run them in the concordance pilot (one run each) and leave them out of `BENCH_METHODS` for the timing series. Time cfDNAlab's normal configuration.
- **Result tables:**
  - The *matched* FinaleToolkit and mosdepth-fragment comparisons use the cfDNAlab proper-pair configuration.
  - The *native* comparison uses cfDNAlab's default.
  - Also report the proper-pair effect itself (cfDNAlab against cfDNAlab-proper).
- **README:** state that the proper-pair flag depends on the aligner. The 12–16% was measured on the published BH01 BAM (BWA aln/sampe). bwa mem also infers proper pairs from its estimated insert-size distribution, so measure the rate again on the hg38 remap (`samtools flagstat`).

### P1-2. Coverage never compares cfDNAlab's fragment coverage

**Evidence**
- Without `--ignore-gap`, `cfdna fcoverage` counts the full fragment span, which is the quantity cfDNA analyses use.
- The benchmark only runs `--ignore-gap`, which counts read-aligned bases (`benchmarks/coverage/benchmark_coverage.sh:250-253`).
- mosdepth v0.3.14 has `-a/--fragment-mode`: it counts the full fragment from read1 using TLEN, for proper pairs only (`mosdepth.nim:346-353`; README line 60).

**Fix**
- Add `cfdnalab-fragment` (the default fcoverage call, without `--ignore-gap`) and `mosdepth-fragment` (`--fragment-mode` with the same flags, MAPQ and length bounds).
- Add `cfdnalab-fragment-proper` (`--require-proper-pair`), because mosdepth's fragment mode only counts proper pairs.
- Make these the headline per-base comparison. Keep `--ignore-gap` against mosdepth's default mode as the read-depth comparison.
- Add the equality check to `compare_coverage.py` (P0-4).

**What mosdepth's fragment mode counts** (`mosdepth.nim:346-353`)
- Only read1 of proper, non-supplementary pairs.
- The span is `[min(start, matepos), min(start, matepos) + |TLEN|)`.
- `--mapq` and `-F` are checked on read1 alone.
- On a bwa-mem BAM, `|TLEN|` is `reverse.reference_end - forward.pos` (P0-1), so for ordinary inward pairs this is exactly cfDNAlab's span. **Expect near-identical depth** against `cfdnalab-fragment-proper`. The residual differences should come only from:
  - mate MAPQ, since mosdepth does not check read2;
  - pairs whose reverse mate starts first, which cfDNAlab rejects;
  - read2-only flags.
- Report the fraction of identical bases, and categorize the remaining positions if feasible.

**deepTools 4.0.0 (optional, coverage only)**
- Use `bamCoverage --extendReads --samFlagInclude 66 --samFlagExclude 3852 --minMappingQuality 30 --minFragmentLength 30 --maxFragmentLength 1000 --binSize 1 --outFileFormat bedgraph -p N` against `cfdnalab-fragment-proper`.
- **`--samFlagInclude 66` (read1 + proper pair) is required.** `--samFlagInclude` requires *all* listed bits (`filtering.rs:217-220`).
  - **Bit 64 (read1):** the bamCoverage loop counts every record that passes the filters (`src/covcalc.rs:435-480`; there is no read2 skip). Each proper-pair record is extended to the full fragment (`src/filtering.rs:460-470`: forward `[start, start+|TLEN|)`, reverse `[mpos, reference_end)`). Without bit 64, every fragment is counted twice.
  - **Bit 2 (proper pair):** records that are not proper pairs are extended by a fixed length, the median of the proper-pair `|TLEN|` values (`filtering.rs:166-199, 448-480`), not by their own mate. Their coverage would not be a fragment span.
  - The matched cfDNAlab configuration is therefore `cfdnalab-fragment-proper`.
- Without `--region`, each chromosome is one task (`covcalc.rs:1630-1660`, `region_divider`), so `-p` parallelizes across chromosomes, bounded by chr1. The bamCoverage path fetches each chromosome without padding (`covcalc.rs:277`); because tasks are whole chromosomes, this loses nothing.
- Verify the help text of the rewritten 4.0 CLI (`bamCoverage2.py`) before relying on any other option.
- Do not use deepTools for lengths (`bamPEFragmentSize` samples regions and double counts at chunk boundaries) or for midpoints (`--MNase` counts 2–4 central bases). The READMEs already explain both.

### P1-3. Per-base timings are dominated by writing very different outputs

| Tool | What it writes |
| --- | --- |
| `samtools depth -aa` | One plain-text line per position, about 3.1 × 10⁹ lines (43–55 GB) per run, with no compression option (`bam2depth.c:857-866`). |
| PanDepth `-a` | Every position, gzip-compressed single-threaded after its workers finish (`PanDepth.cpp:4266-4282`, `gzstream.c:53-61`). Its workers are also capped by contig groups (25 for hs37d5), bounded by chr1 (`PanDepth.cpp:4161-4187`). |
| mosdepth | Run-length encoded BGZF at level 1, single-threaded (`mosdepth.nim:642-651`). |
| cfDNAlab | Run-length encoded zstd, without zero runs. |

Their output cost is part of their native workload, so keep it in the timing (D3). Three points still matter.

**1. The current runners interleave tools.**
- Each repetition block shuffles *all* method × core combinations (`shuf "$BENCH_WORK/combinations.tsv"` in every runner). Tools do not each run their own warmups and repetitions consecutively.
- Interleaving is a deliberate design: it spreads slow drift (thermal state, background load) evenly over the methods.
- The cost is that one tool's leftovers can affect the next, different tool. Examples are page-cache eviction after a large write or a large memory peak (cfDNAPro), and dirty pages still being written back.
- **Fix: run contiguous per-method blocks on one node** (required by D2).
  - Run the methods in a fixed, documented order. Shuffling the block order in a single job is only one random draw and protects against nothing. The anchor runs below are what detect drift.
  - For each method in turn, run its `BENCH_WARMUPS` warmups, then all of its measured repetitions. Shuffle only that method's core budgets within each repetition.
  - Warmups therefore always run immediately before the same tool's measured runs, and no other tool runs in between.
  - To check for drift over the job, run a short cfDNAlab **anchor** (for example 8 cores, 1 warmup and 2 measured runs, labelled `anchor`) before the first method block and after every method block. Report the anchor times.

**2. Page-cache eviction is probably small anyway.**
- On a node with at least 400 GB, the BAM plus 55 GB of samtools output fit in the page cache, so little or nothing should be evicted (*inference*).
- Record the BAM's page-cache residency before each run (`fincore`, P2-4) to confirm this, rather than adding a re-warm step.

**3. Write-back is outside the timer for every tool.**
- Each runner calls `sync` after the timed command. GNU time stops when the process exits, which can be before its dirty pages reach disk.
- The timings therefore measure producing the output into the page cache, not flushing it to disk. Tools whose large output is still dirty in memory at exit would get their disk writing for free. The README currently says output writing is timed (`benchmarks/README.md:67`).
- **Fix: time the flush for every tool and report both numbers.**
  - Immediately after the timed command, time `sync -f "$output_dir"` (GNU coreutils; it flushes the file system holding the outputs) and record it as `flush_seconds`.
  - Report `wall_seconds` (process only) and `wall_seconds + flush_seconds` (output on disk). Use the second as the headline, since every cfDNAlab run writes its outputs too.
  - Keep the existing untimed `sync` as a final safety step.

**Fix**
- Record the bytes written per run (`du -sb` of the output directory, outside the timer) as a TSV column.
- Confirm node-local scratch has room for the BAM, the environments and about 60 GB of samtools output.
- Describe samtools and PanDepth as "writes every position as text" workloads, and keep them out of scaling fits. A reader should see that their flat scaling comes from output, not counting.

### P1-4. BH01 has no duplicate flags

**Status: resolved.** The hg38 remap marks duplicates with Picard `MarkDuplicates` (P0-1). In the Methods, state that FrEIA additionally removes coordinate duplicates itself (`2_data_transformation.py:329-333`) while the other tools use the flag. The original finding follows.

**Evidence**
- No sampled record had flag 0x400. The published BAM was not deduplicated; `FilterUniqueBAM.py` was used only inside the authors' analysis pipes (`shendurelab/cfDNA` README.md:80-86).
- Every `-F 1024` and duplicate filter in the benchmark therefore does nothing.
- FrEIA is the only tool that deduplicates, by coordinates (`2_data_transformation.py:329-333`).

**Fix**
- Check whether the hg38 remap marked duplicates (`samtools flagstat`). If not, run `samtools fixmate -m` and `samtools markdup` after downsampling.
- Otherwise state in the Methods that duplicates are not marked, and that FrEIA alone removes coordinate duplicates.

### P1-5. `prepare_snyder.py` causes Griffin to drop reverse-strand merged fragments

**Status: moot for the hg38 input.** The remap converts the reconstructed pairs to FASTQ and realigns them with bwa mem, which recomputes TLEN and gives the forward mate the positive sign (P0-1). What remains:
- **README:** remove the reconstructed-mate paragraph at `benchmarks/midpoints/README.md:18`.
- **`prepare_snyder.py` (optional):** if it is rerun, restricting reconstruction to `M_` names (see the last fix bullet) is cheap insurance. The reported count of about 3.5 million reconstructed pairs is already consistent with merged reads only.

**Original assessment:** this is an artefact of our preparation, not a Griffin limitation to report.
- Griffin assumes the forward-strand mate has the positive TLEN.
- The SAM specification leaves the sign arbitrary for coincident mates, but bwa-mem's convention puts it on the forward mate (`bwamem.c:888-891`). So does every BAM a Griffin user is likely to have.
- Our script chose the opposite sign for half the merged reads. Fix the script rather than describe this as an expected Griffin difference.

**Evidence**
- The script gives read1, which keeps the original strand, a positive TLEN (`benchmarks/coverage/prepare_snyder.py:66-69`).
- For a merged read that was reverse-strand, the forward-strand mate therefore has a *negative* TLEN.
- Griffin only uses `is_reverse == False and template_length > 0` (`griffin_coverage.py:370`), so it drops about half of all merged fragments.
- The SAM specification allows either sign when the mates have the same leftmost position, as long as the signs differ (`SAMv1.tex:593-597`).
- `benchmarks/midpoints/README.md:18` describes this as Griffin behaviour, but it is caused by our preparation.
- Merged reads are only about 1.8% of fragments, and few fall in 100–200 bp at MAPQ ≥30, so the effect is small. It is systematic, however.

**Fix**
- Give the positive TLEN to the **forward-strand** mate. FinaleToolkit, samtools, Picard and mosdepth are unaffected by the sign:
  - FinaleToolkit handles both signs (`io/alignment.py:252-258`).
  - samtools classifies coincident mates as inward (`stats.c:1347-1350`).
  - Picard uses `MC` to classify them as FR (`SamPairUtil.java:114-118`). Keep setting `MC`.
  - mosdepth uses `abs(isize)`.
- Rewrite `benchmarks/midpoints/README.md:18`.
- Only reconstruct pairs for records whose names start with `M_`, the authors' convention for merged reads (`FilterUniqueBAM.py:364`). Count records by name prefix, and fail on unexpected unpaired names. A non-merged singleton would otherwise become a fake fragment of read length. None were seen in the sample, but the full 217 GB file was not scanned.

### P1-6. Summed RSS overstates the memory of multi-process tools

**Evidence**
- `benchmarks/measure_memory.py:16-33` sums RSS across the process tree.
- Forked worker pools share copy-on-write pages, and each worker counts them in full. This affects:
  - FinaleToolkit's `Pool` (`frag/_coverage.py:244-248`, the end-motif `Pool`);
  - Griffin's `Pool(CPU)` in both scripts;
  - FrEIA's extraction `Pool`.
- Pipelines of separate programs, such as samtools | Picard, only share library pages, so their sum is close to correct.
- The overstatement grows with the core budget. A reviewer comparing memory against cfDNAlab, which is a single process, will see an artefact that favours cfDNAlab.

**Fix**
- Sample PSS from `/proc/<pid>/smaps_rollup`; a sum of PSS counts shared pages once. Keep RSS as a second column.
- For single-process tools, report GNU time's maximum RSS as the primary number, since it is exact there.
- Read `/proc` directly instead of forking `ps` every second.
- The sampler shares the step's CPUs, so measure its overhead in the pilot: cfDNAlab at 1 core with and without the sampler.
- If `smaps_rollup` is slow for processes over 100 GB, sample PSS less often for those trees and document the interval.

### P1-7. cfDNAPro is memory-bound and single-threaded

**Evidence**
- `readBam` reads the whole BAM with a single `readGAlignmentPairs` call, with no chunking or BiocParallel (`R/readBam.R:247`, `GenomicAlignments R/readGAlignments.R:81-87,272-282`).
- The OMP/BLAS thread variables passed in `benchmarks/lengths/benchmark_lengths.sh:232` have no effect.
- Estimated peak at 15x: 120–250 GB (*estimate*).
- The FinaleToolkit paper reports that cfDNAPro end motifs exceeded 256 GB on BH01 (Bioinformatics Advances 2025, PMC12597888).
- `readGAlignmentPairs` stops the whole run (`stopifnot`) if two mates disagree on the proper-pair or secondary bits (`readGAlignments.R:224-232`).

**Why the estimate is so large** (*estimate*, not measured)
- At 15x there are roughly 6 × 10⁸ alignment records.
- `use_names=TRUE`, which the outward-pair check needs, keeps every read name as an R string, at about 50–100 bytes each.
- `readGAlignmentPairs` then makes full subset copies while pairing (`readGAlignments.R:197,214-215`), and `readBam` keeps further copies during curation.
- The FinaleToolkit paper's report did not state the depth in the part checked here. If it used the full-depth BH01 BAM, the 15x figure would be several times smaller. The pilot settles this.

**Fix**
- Pilot on a 0.1x downsample and extrapolate memory and time (D4).
- Choose a node type with enough memory for every tool, and run **all** tools on that node type, so timings remain comparable. `--mem=0` already gives the job all of that node's memory.
- Run at 1 and 32 cores only (P0-3).
- If the 15x run still fails with an out-of-memory error, record it as a result, with the node memory, rather than silently dropping cfDNAPro.
- Do not patch cfDNAPro. Reading per chromosome (`ScanBamParam(which=)`) would bound its memory, because pairs across chromosomes are discarded anyway. But the benchmark measures the released tool. The paper can explain the memory use as whole-file loading.

### P1-8. FrEIA `-n 10` is not the authors' setting

**Evidence**
- The authors' rule passes `-n {nBase}` (`workflow/rules/2_FrEIA/1_preprocessing.smk:73-78`), and the shipped config sets `nBase: 3` (`config/config.yaml:37`). 10 is only the script's own default.
- With 10, FrEIA's N filter (`2_data_transformation.py:336-337`) removes pairs that have an N anywhere in the 10 extracted bases, not only in the 3-mer. It also increases memory.

**Fix**
- Use `-n 3` at `benchmarks/ends/benchmark_end_motifs.sh:268`.
- Update the comment at lines 252-253 and `benchmarks/ends/README.md:12` ("10-base sequence").

### P1-9. The ends benchmark has no cfDNAlab configuration matching FinaleToolkit's ignored length bounds

**Evidence**
- FinaleToolkit 1.1.0 accepts `--min-length`/`--max-length` for `end-motifs` but never applies them (`frag/_end_motifs.py:107-120`, `io/alignment.py:242-268`). The README states this correctly.
- Every cfDNAlab run still filters 30–1000 bp, so length filtering is mixed into every FinaleToolkit difference.

**Fix**
- Add `cfdnalab-reference-proper-unbounded`: `--source-inside reference --require-proper-pair --min-fragment-length 10 --max-fragment-length 50000`, the widest range the runner accepts.
- Report it as the matched FinaleToolkit comparison. The remaining differences should then be:
  - TLEN against the actual span;
  - CIGAR filters;
  - FinaleToolkit's double counting of read1 alignments that cross its 1 Mb window boundaries (`_motif_common.py:541-576`, `io/alignment.py:245`; about 1e-4 of fragments, *estimate*).

### P1-10. The lengths README misstates the orientation rules

**Evidence**
- There are two rules, not four:
  - cfDNAlab keeps a pair when `forward.pos <= reverse.pos` (`src/shared/fragment/minimal_fragment.rs:179-180`). samtools' `IS` inward column uses the same rule (`stats.c:1329-1352`).
  - Picard's FR uses forward start < reverse end (`SamPairUtil.java:101-136`). cfDNAPro's check is `start(first) < end(second)` in 1-based inclusive coordinates (`R/internalFunctions.R:592-599`), which is `forward.pos < reverse.reference_end - 1` in 0-based terms. This only differs for 1 bp spans.
  - Overlapping pairs where the reverse mate starts first are therefore kept by Picard and cfDNAPro, and rejected by cfDNAlab and samtools.
- samtools stores halved bins as integers and rounds down (`stats.c:1576-1578`, `stats_isize.c:139-145`). A bin with an odd number of records, for example when mates are filtered separately by the MAPQ input filter, loses half a fragment.

**Length definition**
- The README calls the TLEN-based tools' length "absolute TLEN" and warns that it may differ from cfDNAlab's span.
- The hg38 remap uses bwa mem (P0-1), which computes `|TLEN|` from the 5' ends (`bwamem.c:888-891`). For inward pairs this equals `reverse.reference_end - forward.pos`, the same length as cfDNAlab's.
- Remaining length differences between cfDNAlab and samtools, Picard or FinaleToolkit then come from *selection*: orientation rules, MAPQ on one mate or both, proper pairs, and halving. They do not come from the length definition.
- Say this in the README. It strengthens the equality claim.

**Fix** — update the "Expected differences" column and the paragraph at `benchmarks/lengths/README.md:58`:
- State the rules precisely.
- Change `forward.pos < reverse.reference_end` to the 1-based form above.
- Add the truncation to the samtools row.
- Add the length-definition point above.

### P1-11. The midpoint site set is not specified

**Evidence**
- `benchmarks/midpoints/README.md:34` asks for "the biological site set" but never names one, nor the window width or the number of sites.
- Concordance and timing both depend on that choice.
- No current group exercises reverse-strand handling, although it is verified to match (Griffin flips the array per site, `griffin_merge_sites.py:372-373`; cfDNAlab mirrors per site, `src/commands/midpoints/strand.rs:15-16`).

**Fix**
- Fix the site set in the README and in a preparation script. For example:
  - a few transcription factor binding site lists in Griffin's own format (CTCF plus one or two others), each limited to a fixed top N by peak count. Griffin supports this through `--number_of_sites`/`--sort_by`, but the runner passes all sites, so preselect them in the preparation script;
  - plus a stranded housekeeping TSS group to exercise `-` strands.
- Fix the window, for example ±1000 bp. Record the source URLs and checksums.

### P1-12. Provenance is incomplete for a paper

**Fix**
- **cfDNAlab:**
  - `cfdna --version` prints only `cfdna 0.9.0` (`cli_app.rs:48`), while the benchmark builds an arbitrary checkout.
  - Record `git rev-parse HEAD` and a dirty flag in the installer output and in every TSV row.
  - Tag a release before the final runs and cite it.
- **Conda environments:** write `conda list --explicit` for each environment next to its archive.
- **cfDNAPro image:**
  - The pinned digest is exactly tag `1.7.4`, with `org.opencontainers.image.revision=ab9eeb1…`, so the installer comment is correct.
  - Still record the `apptainer inspect` labels and the `GenomicAlignments` and `Rsamtools` versions inside the image with each run, because `packageVersion()` alone cannot identify the commit.
  - In the paper, call it "cfDNAPro GitHub `ab9eeb1` (DESCRIPTION 1.7.4), container digest …". Bioconductor's release is 1.18.0, a different version sequence, and its `callMotif` is a stub.
- **Comparison environment:** record `pip freeze`.
- **Input data:** see the data script in P0-1.

### P1-13. The coverage runner is inconsistent with the other three

**Evidence**
- It hard-codes cores, warmups, repetitions and methods (`benchmarks/coverage/benchmark_coverage.sh:72-75,180`), so it has no pilot mode.
- Its TSV lacks the settings columns the others record (MAPQ, length bounds, scope, proper pair).

**Fix**
- Use the same pattern as the other runners: environment overrides with validation, settings columns, and profile capture (P0-4).

## P2: fix before the final runs

1. **Installation.** See the installation section: channel handling, FinaleToolkit pins and conda-pack 0.9.2.
2. **Smoke-test cfDNAlab on a compute node from the installer** (referred to as P2-2).
   - The runners copy `cfdna` out of its Conda environment.
   - `plotters` loads fontconfig at run time via `dlopen` (`Cargo.toml:86`), so `ldd` cannot show the dependency.
   - `lengths` and `midpoints` write a QC plot by default. If fontconfig or a sans-serif font is missing on the node, the plot error propagates after the main outputs are written, and the run exits non-zero (`src/commands/lengths/lengths.rs:737-747`, `src/commands/midpoints/plotting.rs:101-103`).
   - Run a tiny `cfdna lengths` call on a compute node, and `ldd` both `cfdna` and `pandepth`, failing on `not found`.
3. **Silent exit in the coverage runner.** `benchmarks/coverage/benchmark_coverage.sh:62-64` uses `${BENCH_CFDNA:-$(command -v cfdna)}`. If the tool is missing, `set -e` exits with no message. Use `: "${BENCH_CFDNA:?...}"` as the other runners do.
4. **Record state that can explain timing outliers**, outside the timer:
   - the BAM's page-cache residency before each run (`fincore sample.bam`);
   - the CPU frequency governor and turbo state, once per job;
   - bytes written (P1-3).
5. **Statistics.**
   - Report the median and IQR, or the minimum, alongside mean and SD.
   - Plot every measured point.
   - Give speed-up ratios with bootstrap intervals.
6. **Depth series.** Timing on a single 15x input cannot show how runtime scales with input size or the fixed start-up cost (R and container start-up, JVM, Python imports). Add 1x, 5x and 15x (and 30x if feasible) at a fixed core count, for cfDNAlab and the fastest competitor per feature.
7. **Ends comparator metrics.** Pearson r over 256 frequencies is dominated by abundant motifs. Also report Spearman r and the maximum |log2 ratio| among motifs with frequency ≥1e-4.
8. **Make the installer resumable** (requested). The installer is a single `set -e` script, so one failure (P0-5, P0-6) means reinstalling every tool that had already finished. Proposed design:
   - **Steps.** Put each tool in its own function (`install_cfdnalab`, `install_mosdepth`, …), run in the current order. Write each step's output to `$install_root/logs/<tool>.log`, and report which step failed.
   - **Completion record.** When a step finishes *and* its verification passes (version check, archive packed), append a row to `$install_root/manifest.tsv`:
     - `tool`, `version`, `route`;
     - the artefact paths (environment prefix, archive, binary or image);
     - a `sha256` of each archive, binary or image;
     - `completed_utc`;
     - for cfDNAlab, the repository commit and a dirty flag.

     Write the row atomically (temporary file, then `mv`). A tool counts as done only if it has a row *and* its artefacts still exist.
   - **Rerunning.** The default resumes: skip completed tools, and redo the rest. Before retrying, remove only that tool's own partial artefacts inside `$install_root` (its `envs/<tool>` prefix, archive and download directory), never anything outside.
   - **Flags.**
     - `--force TOOL[,TOOL]` or `--force all` redoes the named steps.
     - `--only TOOL[,TOOL]` runs a subset.
     - For cfDNAlab, rebuild automatically if the checkout's commit or dirty state differs from the manifest. Set `CARGO_TARGET_DIR=$install_root/build/cfdnalab` so rebuilds are incremental.
   - **Output.** Print the `BENCH_*` exports from the manifest (completed tools only) at the end. Replace "use a new directory" in `benchmarks/README.md` with the resume behaviour.
   - **Scope.** No `--continue` flag is needed if resuming is the default. If a stricter default is preferred, refuse a non-empty `$install_root` unless `--resume` is given.
9. **Pilot and smoke-test presets.**
   - There is no `--pilot` flag today. The ends, lengths and midpoints runners can already run a pilot through environment variables (`BENCH_WARMUPS=0 BENCH_REPETITIONS=1 BENCH_CORES=… BENCH_METHODS=…`, see each README). The coverage runner cannot (P1-13).
   - Add a uniform `BENCH_PILOT=1` to all four runners. It should set no warmups and one repetition, keep `BENCH_CORES` (so P0-3's wall-time pilot can cover every budget), and write to a results directory marked as a pilot.
   - Document a **smoke-test input** for checking the wiring in minutes before the real pilot. Use a small BAM made with `samtools view -b -s 42.01` (1% of fragments, chosen by read-name hash so both mates are kept), and for midpoints, a small site file.
10. **FinaleToolkit window edge case.** If a contig's length is an exact multiple of 1,000,000, its last 1 Mb window is never generated (`_motif_common.py:541-576`; reported by source reading, not run). Assert this cannot happen in the runner's header check (`benchmarks/ends/benchmark_end_motifs.sh:126-144`).

## Installation (`install_tools.sh`) by tool

Versions were checked on 2026-10-29 against GitHub releases, PyPI and `api.anaconda.org` for conda-forge and Bioconda. Solvability was checked with `conda create --dry-run` (linux-64, glibc 2.28, strict priority). Release binaries were inspected with `file` and `objdump` but not run.

The installation policy is "latest version, prefer Conda over building from source". Every route below follows it; only cfDNAlab itself is compiled.

| Tool | Pinned | Latest | Conda | Route in script | Verdict |
| --- | --- | --- | --- | --- | --- |
| cfDNAlab | checkout | — | — | Cargo build with Rust 1.94.0 and clang 21 from conda-forge, as in the repository README | OK. Rust 1.94.0 matches `rust-version`; the latest conda-forge Rust is 1.99.1. Record the commit (P1-12). |
| mosdepth | 0.3.14 | 0.3.14 (2026-04-24) | Bioconda 0.3.14 | Release binary, statically linked | OK, latest. This is the upstream README's first route. |
| samtools | 1.25 | 1.25 (2026-10-27, security fixes) | Bioconda 1.25 with HTSlib 1.25 | Bioconda | OK, latest. |
| PanDepth | 2.26 | 2.26 (2025-03-25) | Bioconda 2.26 | Release tarball; the layout `PanDepth-2.26-Linux-x86_64/pandepth` is correct | OK, latest. The binary is dynamically linked: libstdc++ (GLIBCXX ≥3.4.11), glibc ≥2.14, plus a stale RPATH. Either keep it with the `ldd` check (P2-2) or switch to Bioconda. |
| FinaleToolkit | 1.1.0 | 1.1.0 (2026-07-23) | Bioconda 1.1.0, noarch, `python >=3.10,<3.13`; includes pysam, py2bit, pybigwig | Bioconda, Python 3.11 | Version OK, **but see dependency drift below.** |
| cfDNAPro | image digest `26a28821…` | GitHub `master` = `ab9eeb1`, DESCRIPTION 1.7.4 | Bioconductor 3.23: 1.18.0; Bioconda: 1.16.1 (separate histories) | The authors' Apptainer image (their recommended "Option 1") | OK. The digest is exactly tag `1.7.4`, with the labels `revision=ab9eeb1…` and `version=1.7.4`, built 2026-10-26. This is the newest code. |
| Picard | 3.5.0 (picard-slim) + OpenJDK 17 | 3.5.0 (2026-07-31); requires Java 17 | Bioconda picard-slim 3.5.0 (`openjdk >=17`) | Bioconda | Version OK. The launcher honours `JAVA_HOME` and passes `-Xm*`, `-XX*` and `-D*` to the JVM. **The `--version` check fails (P0-5).** |
| FrEIA | `2c32fb9` | `2c32fb9` is HEAD; no releases | — | The authors' `FrEIA_env.yaml` | OK. It solves (Python 3.8.12). The YAML lists only conda-forge and Bioconda. The `pip:` section installs unpinned `statannot` and `combat` from PyPI; record `pip freeze`. |
| Griffin | v0.2.0 | v0.2.0 (2022-05-10, prerelease) | none | Wiki pins via conda-forge and Bioconda | **Unsolvable (P0-6).** The other pins exist as Bioconda py37 builds and are not shadowed by conda-forge. |
| conda-pack | 0.8.1 | **0.9.2** (2026-06-16, noarch, `python >=3.10`) | conda-forge 0.9.2 | conda-forge | **Update to 0.9.2.** The 0.9.x `conda-unpack` imports `bin/conda_unpack_progress.py` from its own directory. That works with the runners' `python …/bin/conda-unpack` calls, but would not work under `python -I` or `-P`. The code has no syntax newer than Python 3.6 (it was not run under 3.7). |
| pysam (Snyder preparation) | 0.24.1 | 0.24.1 (2026-09-07) | Bioconda | Bioconda | OK, latest. Only needed if the hg38 remap still contains merged single-end reads (P0-1). |

**Channel handling (P2)**
- Without `--override-channels`, channels in the user's `.condarc` are appended at the lowest priority.
- If one of them is `defaults`:
  - Strict priority still prefers conda-forge and Bioconda, but `defaults` can supply packages that are missing from both.
  - Recent Miniconda and Anaconda installers bundle the Anaconda terms-of-service plugin. It can block a non-interactive job until the terms are accepted.
- **Fix:**
  - Add `--override-channels` to every `conda create`.
  - `conda env create` (used for FrEIA) has no `--override-channels`. Before it, fail if `conda config --show channels` lists `defaults`.

**FinaleToolkit dependency drift (P2)**
- The PyPI 1.1.0 metadata pins the versions the release was tested with: `pysam==0.24.0`, `pybigwig==0.3.25`, `py2bit==0.3.3`, `numba==0.63.1`, `scipy==1.17.0`, `statsmodels==0.14.6`.
- The Bioconda recipe leaves them unpinned. Today it resolves pandas 3.1.0, numpy 2.4.6, py2bit 1.0.1 and numba 0.68.0. A pandas or numpy major version could change behaviour in an untested way.
- **Fix:** add the PyPI pins to the `conda create` line. A dry-run with them solves (numpy 2.3.5). This keeps the Conda route and installs FinaleToolkit's tested dependencies.

**Other installer points**
- `install_tools.sh` prints paths but does not check the copied binaries' linkage on a compute node; see P2-2.
- Write `conda list --explicit` for every environment next to its archive (P1-12).
- The installer requires every tool even when only one benchmark is run. `--only` in the resumable design (P2-8) covers this.

**Comparators a reviewer may expect** (not inspected in source unless noted):
- deepTools 4.0.0 (2026-09-05, Bioconda; new multi-threaded Rust core). This is the most likely expectation, for `bamCoverage` and `bamPEFragmentSize`; the README's deepTools claims were verified in source.
- Triton 2.0.0 (GavinHaLab, 2026), the same lab's successor to Griffin. It is not on Bioconda. Expect "why Griffin and not Triton?" and answer it in the README.
- Qualimap 2.3, for insert-size QC.
- bedtools genomecov and megadepth, as coverage baselines.

## P3: wording and minor

- `benchmarks/midpoints/README.md:20`: the padding threshold is off by one. A forward read is guaranteed inside Griffin's fetch window when `floor(max/2) <= min`. Misses become possible from `max >= 2*min + 2`, not `max > 2*min`.
- `benchmarks/midpoints/benchmark_midpoints.sh:139-142` silently turns an invalid strand token into `.`. This contradicts "fails on invalid input" (`benchmarks/midpoints/README.md:34`). Raise an error instead.
- `benchmarks/ends/README.md:43` says the pins were verified on 2026-10-28; `benchmarks/README.md:9` and `install_tools.sh:11` say 2026-10-29.
- Core budgets: 1, 2, 4, 8, 16, 32 gives evenly spaced points on a log axis for scaling plots, unlike 1, 4, 8, 12, 32.
- The GNU time `--format` string embeds tool versions and the CPU model, where `%` and `\` are interpreted. The current values are safe; pass identifiers some other way if this becomes a problem.
- FinaleToolkit's CLI default for `end-motifs` MAPQ is 20 (`cli/commands/__init__.py:493`), unlike its Python API (30). The runner sets 30 explicitly, so no change is needed. Mention it if defaults are discussed.
- cfDNAFE (excluded): besides the documented motif bug, its `bam` CLI silently ignores the chromosome and length options, because the keyword names don't match (`bam.py:27-46` against `Run_BamProcess.py:10-24`). Add this to `benchmarks/ends/README.md` to strengthen the exclusion.

## Verified as correct (do not change)

**cfDNAlab** (repository source):
- Every flag in all four runners exists and parses.
- `--chromosomes all` takes every `@SQ` contig.
- Default outputs are `fcoverage.per_position.bedgraph.zst`, `length_counts.tsv.zst` (with and without `--by-bed`), `end_motifs.zarr` (global) and `midpoint_profiles.zarr`.
- `--length-bins 30 1001` gives one half-open bin; `30:1001:1` gives 971 bins of 1 bp.
- `--ref-2bit` is not needed for `--source-inside read --k-outside 0`.
- `--n-threads` caps CPU use: BAM decompression runs inside the Rayon workers, and no htslib threads are started.
- The hard-clip, terminal-indel and 256 bp soft-clip claims in `benchmarks/ends/README.md:21` hold.
- The even-length midpoint is seeded from a coordinate hash.
- `py-cfdnalab` calls in all three comparators match the API, including the `_ACGT` labels, the length bin columns and the midpoint array shape and group order.

**Griffin v0.2.0:**
- Every argument passed exists, and no required argument is missing.
- The output paths and TSV columns match.
- Unit GC bias gives weight 1.0, and fragments with N in the reference are kept in the uncorrected output.
- Selection is inclusive `min <= |TLEN| <= max`, and the even-length midpoint is the right central base.
- Window coordinates match cfDNAlab exactly for `+`, `.` and `-` sites.
- The group profile is the site mean divided by its window mean, which matches the comparator's normalization.
- The lexicographic `chrom.sizes` is required by pyBigWig (`pyBigWig.c:1574-1579`).
- The scripts on main are identical to v0.2.0.

**FinaleToolkit v1.1.0:**
- All flags used by the three runners exist.
- The `end-motifs` TSV has no header and lists all 4^k motifs as frequencies.
- `frag-length-bins` writes `min max count` with `min == max` at bin size 1, and its bounds are inclusive.
- `coverage` normalization is off by default, and its bounds are inclusive.

**FrEIA:**
- The samtools filter matches the authors' rule; repeated `-F` options accumulate in samtools 1.12.
- The sample sheet is whitespace-separated.
- The output path and the exclusive `-fsmax` are as the runner expects.
- The 8000 MiB serial threshold matches exactly.
- The right-end labels are reversed but not complemented, so complementing them in the comparator is correct.

**samtools 1.25:**
- `stats` accepts `--input-fmt-option filter=` and applies it to every read; `-f` and `-F` exist.
- `-i 1001` places overflow in the last bin; `-m 1` truncates nothing.
- The `IS` inward count is column 4, matching `row[3]` in the comparator.
- `depth` has `--require-flags`. `-G` adds to the default exclusions, and `-s` handles overlapping mates.

**Picard 3.5.0:**
- `HISTOGRAM_FILE` is optional in 3.5.0 (it was required in 3.4.0), so **keep the pin**.
- `HISTOGRAM_WIDTH` overrides the deviation trim.
- The `All_Reads.fr_count` column name is correct.
- `/dev/stdin` input works.
- Java 17 is required.

**PanDepth v2.26:** the `sed` regex extracts `v2.26`.

**cfDNAPro:**
- The README's MAPQ 40 default and the `bamfile=` path not forwarding `...` are both correct, so supplying the fragment object is necessary.
- The output columns are `insert_size` and `All_Reads.fr_count`.
- Bounds are inclusive, and zero bins are filled.

**cfDNAFE:** the motif bug described in `benchmarks/ends/README.md:38` is real and is unfixed in every commit, branch and fork.

**conda-pack:**
- `conda-unpack` is generated even for environments without Python, and it uses only the standard library.
- Its code uses no syntax newer than Python 3.6 (read, not run), so Griffin's Python 3.7 can execute it.
- Prefix length is not a problem, because binary files are packed with the long build placeholder.
- Running it with another environment's Python, as the runners do, works: the target prefix is taken from the script's own location.
