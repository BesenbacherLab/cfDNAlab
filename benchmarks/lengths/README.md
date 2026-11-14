# Fragment length benchmark comparators

Compare cfDNAlab against **cfDNAPro, FinaleToolkit, samtools and Picard**. cfDNAPro checks the fragment definition, FinaleToolkit provides a direct cfDNA histogram command, and samtools and Picard provide established insert-size comparators.

Compare global histograms with one bin per fragment length from **30 through 1000 bp inclusive**, MAPQ >=30, and no GC correction, genomic scaling or blacklists. These are bins on the fragment length axis, not genomic intervals. Keep the prepared BAM unchanged and let each tool apply its filters.

| Tool | Native calculation and settings | Role | Expected differences |
| --- | --- | --- | --- |
| cfDNAPro 1.7.4 | `readBam(curate_start_and_end=TRUE, galp_mapqFilter=30)`, with explicit matching flag exclusions, then `callLength(frag_obj=..., isize_min=30L, isize_max=1000L)` | Primary definition check | Uses the same directional aligned span as cfDNAlab. Overlapping outward-pair handling and native pairing can differ, as described below. Neither tool needs proper-pair filtering for this comparison. |
| FinaleToolkit `d16632ec` | `frag-length-bins BAM --bin-size 1 --min-length 30 --max-length 1000 --min-mapq 30 --output lengths.tsv` | Direct cfDNA CLI comparator | Uses absolute TLEN from read1 and requires proper pairs. Does not verify both mates' MAPQ. This command has no worker-count option. |
| samtools stats (full QC report, insert-size section extracted) 1.25 | `stats --input-fmt-option 'filter=mapq >= 30' -f 1 -F 3852 -m 1 -i 1001`, using its `IS` inward-pair column and retaining lengths 30-1000 afterward | Established compiled baseline | Uses absolute TLEN. The last bin collects lengths at or above its maximum, hence 1001 rather than 1000. Counts passing alignments, halves each bin and truncates the integer result, so independently filtering mates can change counts. Computes many additional QC statistics. |
| Picard 3.5.0 | `picard -Xmx2g CollectInsertSizeMetrics I=sample.bam O=metrics.txt HISTOGRAM_WIDTH=1000` | Native insert-size timing | No MAPQ or QC-fail prefilter inside timing. Concordance uses the separately generated, untimed filtered histogram. No PDF. |

For cfDNAlab, use `lengths --length-bins 30:1001:1 --indel-mode ignore --clip-mode aligned --min-mapq 30`. Use the same explicit chromosome selection on both sides. The production comparison uses all BAM contigs. A separately labeled `--require-proper-pair` run diagnoses FinaleToolkit's selection difference once at 32 cores for concordance. Normal cfDNAlab runs supply the timing comparison at all four core settings.

## cfDNAPro settings and definition

The inspected author source reports version **1.7.4**, commit `ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c`, checked on 2026-10-29. Its curation constructs each fragment from the forward alignment's start to the reverse alignment's end. The inclusive R coordinates give the same width as cfDNAlab's half-open `forward.pos` to `reverse.reference_end` interval. Neither calculation extends soft clipping or adjusts the length for insertions and deletions in this configuration.

Read the BAM explicitly with the shared MAPQ and flag settings, then pass its fragment object to `callLength()`, the entry point used in cfDNAPro's tutorial and paper. `readBam()` defaults to MAPQ 40, so MAPQ 30 remains explicit. The resulting histogram may include a zero-filling `prop` column, which the comparator ignores.

```r
library(cfDNAPro)

# Select all BAM contigs to match the production comparison scope
# Choose the genome label matching the alignment reference
fragments <- readBam(
    bamfile = "sample.bam",
    genome_label = "hg38",
    chromosome_to_keep = FALSE,
    curate_start_and_end = TRUE,
    galp_mapqFilter = 30,
    galp_flag = Rsamtools::scanBamFlag(
        isPaired = TRUE,
        isUnmappedQuery = FALSE,
        hasUnmappedMate = FALSE,
        isSecondaryAlignment = FALSE,
        isSupplementaryAlignment = FALSE,
        isDuplicate = FALSE,
        isNotPassingQualityControls = FALSE
    ),
    outdir = NA
)

# Return the native histogram without saving a large per-fragment RDS file
counts <- callLength(
    frag_obj = fragments,
    isize_min = 30L,
    isize_max = 1000L
)
write.table(counts, "lengths.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

The runner times the complete R process, including package loading, BAM reading, pairing, curation and histogram writing. Timing only the second function would exclude most of the work. Record peak memory as well, since cfDNAPro materializes the paired alignments and fragments. A speed difference must be measured rather than inferred from R versus Rust.

Keep the native cfDNAPro filters visible. Its default `scanBamFlag` excludes unmapped, secondary, supplementary and duplicate records, but does not explicitly exclude QC-fail records. The benchmark supplies `isNotPassingQualityControls=FALSE` and `hasUnmappedMate=FALSE` through `readBam(galp_flag=...)` to match cfDNAlab's exclusions. Alignment fields and NM/MD tag loading retain the package defaults. Its strict 1-based 5'-end comparison is equivalent to `forward.pos < reverse.reference_end - 1`, whereas cfDNAlab additionally requires `forward.pos <= reverse.pos`. Overlapping pairs with the reverse alignment starting first can therefore remain a source of differences even with matched flags. GenomicAlignments pairing and reference-bound checks can also affect selection.

For ordinary inward bwa mem alignments, stored TLEN spans the forward `pos` through the reverse `reference_end`, matching the directional aligned length. This is not a guarantee for arbitrary aligners or abnormal pairs. [bwa 0.7.19 source](https://github.com/lh3/bwa/blob/v0.7.19/bwamem.c#L887-L892) supplies the inspected rule. Record the actual remapping version. samtools inward selection uses forward/reverse start order, allowing coincident starts. Picard's FR rule instead compares 5' positions with <=, equivalent to `forward.pos < reverse.reference_end`. Its one-base boundary difference from cfDNAPro is outside this benchmark's >=10 bp range.

## Timing scope

samtools `stats` has no dedicated MAPQ switch. Its `-q` controls base-quality trimming. The runner uses HTSlib's general `--input-fmt-option 'filter=mapq >= 30'` instead, applying MAPQ in the same BAM read pass. Select the `IS` inward-pair column after timing. This does not establish that both mates pass MAPQ or reproduce cfDNAlab pairing. Timing includes the full native command needed to obtain the histogram, including its other statistics.

Picard timing uses its documented standalone `picard -Xmx2g CollectInsertSizeMetrics I=sample.bam O=metrics.txt HISTOGRAM_WIDTH=1000` call. It has no MAPQ option, so timed output retains low-MAPQ and QC-fail pairs and timing rows record `min_mapq=NA`. Once after the timed runs, an untimed `samtools view -u -q 30 -f 1 -F 3852 | picard -Xmx2g CollectInsertSizeMetrics I=/dev/stdin ...` pipeline writes `concordance.metrics.txt`. Only that filtered histogram enters concordance. The extra pass adds job wall time but no reported tool time. Picard's optional PDF is omitted because picard-slim lacks R. cfDNAlab's default QC plot remains timed.

Use the commented [production run procedure](../README.md#runs-and-budgets) for the complete 1/4/12/32-core matrix. It times the five normal configurations at each core setting, including cfDNAPro and FinaleToolkit. The proper-pair cfDNAlab configuration runs once at 32 cores for concordance and is excluded from timing summaries. The runner executes every explicitly selected method and does not automatically skip tools. cfDNAlab's `--n-threads` and samtools `-@` receive the allocated core count directly, including the samtools filter in the separate untimed Picard concordance pipeline. Picard has no thread option and receives no JVM processor override. No numerical-library thread overrides are added. The timing TSV records `NA` for a command without an explicit native thread setting.

Plot wall time against allocated cores for every method. Calculate effective cores as `(user_seconds + system_seconds) / wall_seconds` and allocated CPU utilization as that value divided by `cores`. This shows actual utilization without assuming that a tool will or will not scale. Picard's utilization reflects its standalone native call, and cfDNAPro's includes container startup. Include native output generation in timing and perform concordance calculations afterward.

Compare raw totals and per-bin count differences as well as normalized frequencies. Report Pearson r and MAE/RMSE of the frequencies. High correlation alone can conceal a systematic count difference. Supply zero counts for absent length bins before comparing.

## Install and run

Use the shared [../install_tools.sh](../install_tools.sh). It installs the authors' cfDNAPro 1.7.4 Apptainer image by immutable digest, including its R dependencies and reference packages. Picard 3.5.0 uses Bioconda's `picard-slim` package with Java 17 under `envs/picard`. It contains the official JAR but omits R, which this command only needs for its optional PDF. samtools 1.25 uses a separate Bioconda environment under `envs/samtools`, packed as `samtools.tar.gz` with its shared libraries. The runner copies the image and packed environments locally. It records the actual cfDNAPro package version and image SHA-256 in each timing row. Apptainer must already be available on the cluster. See the [installation audit](../README.md#installation) for sources and deliberate adaptations.

The shared [production run procedure](../README.md#runs-and-budgets) contains the local path settings, commented submissions, new result directories and accounting commands. `BENCH_INDEX` can override the usual `.bam.bai` path. `BENCH_MIN_LENGTH`, `BENCH_MAX_LENGTH` and `BENCH_MIN_MAPQ` override the common filters. The automatic comparison uses the same bounds. `BENCH_GENOME=hg38` selects cfDNAPro's UCSC hg38 reference label and must match the alignment reference and chromosome naming. The input remains a coordinate-sorted BAM with its original non-autosomes and alignment flags.

For full runs, the submission launcher reserves an exclusive node and all its memory. Every sequential `srun --jobid` attaches to that allocation and checks the actual hostname and logical/physical CPU IDs. Hard confinement also depends on the cluster's Slurm cgroup configuration. Each method runs as a contiguous block at the job's core count, retaining actual order in `timings.tsv`. Normal timing jobs use a warmup followed by three measured repetitions. The proper-pair concordance job runs once without a warmup at 32 cores. No caches are evicted.

The timing result is a single `timings.tsv`, with a row for every attempted run, including warmups and failures. A native histogram per method is retained under `profiles/` from the first measured repetition at the job's single core count. These small files support concordance analysis. Other temporary outputs and logs are deleted. Failed diagnostics are retained, the remaining repetitions of that method are skipped and other selected methods continue. The job exits nonzero after those methods finish. See the [shared failure and timeout behavior](../README.md#timing-and-slurm-memory). Input staging, scheduler startup, profile copying, cleanup and concordance calculation are outside the timer. Native startup and output generation remain timed, including cfDNAlab's default QC plot and cfDNAPro's Apptainer launch. GNU time's maximum RSS is not the simultaneous total memory of a process tree.

Pilots and full runs automatically write `concordance.tsv` after all tool attempts finish, using this job's inclusive length bounds. This comparison is outside every tool's timing and step memory measurements. Single-method jobs retain their profile and explain that there is no pair to compare. A comparison failure fails the job, with details in `concordance-error.log`, while preserving completed timings and profiles.

## Supplementary result templates

The tables below contain no measured results. Fill them from the saved concordance metrics and completed timing/accounting records.

| cfDNAlab vs tool | Pearson r | MAE / RMSE | Count difference | Wasserstein, bp | Main expected differences |
| --- | --- | --- | --- | --- | --- |
| cfDNAPro 1.7.4 | [ ] | [ ] / [ ] | [ ] | [ ] | Pair selection and overlapping outward-pair checks. |
| FinaleToolkit `d16632ec` | [ ] | [ ] / [ ] | [ ] | [ ] | Proper-pair selection, read1 MAPQ and TLEN. |
| samtools stats (full QC report, insert-size section extracted) 1.25 | [ ] | [ ] / [ ] | [ ] | [ ] | TLEN, per-record filters and halved counts. |
| Picard 3.5.0 | [ ] | [ ] / [ ] | [ ] | [ ] | TLEN from read2 and orientation classification. |

All rows use whole-reference counts, MAPQ >=30 and lengths 30-1000 bp. MAE/RMSE use within-range frequencies. Count differences are cfDNAlab minus comparator.

| Tool / configuration | Cores | Warmups / measured | Mean process+flush, s | SD seconds | Effective cores |
| --- | --- | --- | --- | --- | --- |
| cfDNAlab, all contigs | 1 / 4 / 12 / 32, separate jobs | 1 / 3 each | [ ] | [ ] | [ ] |
| cfDNAPro 1.7.4 | 1 / 4 / 12 / 32, separate jobs | 1 / 3 each | [ ] | [ ] | [ ] |
| FinaleToolkit `d16632ec` | 1 / 4 / 12 / 32, separate jobs | 1 / 3 each | [ ] | [ ] | [ ] |
| samtools stats (full QC report, insert-size section extracted) 1.25 | 1 / 4 / 12 / 32, separate jobs | 1 / 3 each | [ ] | [ ] | [ ] |
| Picard 3.5.0 | 1 / 4 / 12 / 32, separate jobs | 1 / 3 each | [ ] | [ ] | [ ] |

Expand timing rows to a core budget per row in the final A4 supplement. Put workload notes below the table and report maximum RSS separately if needed.

## Source references

- cfDNAlab: [length construction](../../src/shared/fragment/indel_counting_fragment.rs), [orientation checks](../../src/shared/fragment/minimal_fragment.rs), [read filtering](../../src/shared/read.rs), [CLI settings](../../src/commands/lengths/config.rs).
- cfDNAPro: [readBam](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/readBam.R), [curation and orientation](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/internalFunctions.R), [histogram](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/callLength.R).
- FinaleToolkit: [histogram implementation](https://github.com/epifluidlab/FinaleToolkit/blob/d16632ec5ae8999e774234457f6b0f9405c8f590/src/finaletoolkit/frag/_frag_length.py), [alignment reader](https://github.com/epifluidlab/FinaleToolkit/blob/d16632ec5ae8999e774234457f6b0f9405c8f590/src/finaletoolkit/io/alignment.py).
- samtools: [stats manual](https://www.htslib.org/doc/samtools-stats.html), [1.25 counting and output](https://github.com/samtools/samtools/blob/1.25/stats.c), [HTSlib input filter parsing](https://github.com/samtools/htslib/blob/1.25/hts.c), [filtering during BAM reads](https://github.com/samtools/htslib/blob/1.25/sam.c).
- Picard: [3.5.0 release](https://github.com/broadinstitute/picard/releases/tag/3.5.0), [command](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/CollectInsertSizeMetrics.java), [collector](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/directed/InsertSizeMetricsCollector.java).

Join completed Slurm step memory after the job. Timing rows also record filesystem flush time, process-plus-flush time and output bytes. See the [shared timing TSV documentation](../README.md) for units and measurement limits. Submit from the repository root so helper scripts can be staged before warmups.
