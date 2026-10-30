# Fragment length benchmark comparators

Compare cfDNAlab against **cfDNAPro, FinaleToolkit, samtools and Picard**. cfDNAPro checks the fragment definition, FinaleToolkit provides a direct cfDNA histogram command, and samtools and Picard provide established insert-size comparators. No sample analyses or timings have been run while preparing these scripts.

Compare global histograms with one bin per fragment length from **30 through 1000 bp inclusive**, MAPQ >=30, and no GC correction, genomic scaling or blacklists. These are bins on the fragment length axis, not genomic intervals. Keep the prepared BAM unchanged and let each tool apply its filters.

| Tool | Native calculation and settings | Role | Expected differences |
| --- | --- | --- | --- |
| cfDNAPro 1.7.4 | `readBam(curate_start_and_end=TRUE, galp_mapqFilter=30)`, with explicit matching flag exclusions, then `read_bam_insert_metrics(fragment_obj=..., isize_min=30L, isize_max=1000L)` | Primary definition check | Uses the same directional aligned span as cfDNAlab. Overlapping outward-pair handling and native pairing can differ, as described below. Neither tool needs proper-pair filtering for this comparison. |
| FinaleToolkit 1.1.0 | `frag-length-bins BAM --bin-size 1 --min-length 30 --max-length 1000 --min-mapq 30 --output lengths.tsv` | Direct cfDNA CLI comparator | Uses absolute TLEN from read1 and requires proper pairs. Does not verify both mates' MAPQ. This command has no worker-count option. |
| samtools 1.25 | `stats --input-fmt-option 'filter=mapq >= 30' -f 1 -F 3852 -m 1 -i 1001`, using its `IS` inward-pair column and retaining lengths 30-1000 afterward | Established compiled baseline | Uses absolute TLEN. The last bin collects lengths at or above its maximum, hence 1001 rather than 1000. Counts passing alignments, halves each bin and truncates the integer result, so independently filtering mates can change counts. Computes many additional QC statistics. |
| Picard 3.5.0 | `samtools view -u -q 30 -f 1 -F 3852`, piped into `CollectInsertSizeMetrics HISTOGRAM_WIDTH=1000 MINIMUM_PCT=0 INCLUDE_DUPLICATES=false`, using the `All_Reads.fr_count` histogram column | Familiar insert-size comparator | Uses TLEN from the second mate. Explicit histogram width prevents automatic tail trimming. Input filtering remains inside timing. Compare only lengths 30-1000. |

For cfDNAlab, use `lengths --length-bins 30:1001:1 --indel-mode ignore --clip-mode aligned --min-mapq 30`. Use the same explicit chromosome selection on both sides. Keep a separate `cfdnalab-autosomes` run when reporting whole-reference tools alongside cfDNAlab's usual autosome selection. A separately labeled `--require-proper-pair` run can diagnose FinaleToolkit's selection difference. It should not replace the normal cfDNAlab run.

## cfDNAPro settings and definition

The inspected author source reports version **1.7.4**, commit `ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c`, checked on 2026-10-29. Its curation constructs each fragment from the forward alignment's start to the reverse alignment's end. The inclusive R coordinates give the same width as cfDNAlab's half-open `forward.pos` to `reverse.reference_end` interval. Neither calculation extends soft clipping or adjusts the length for insertions and deletions in this configuration.

There are two important configuration details. `readBam()` currently defaults to **MAPQ 40**, although the README says 30. Also, `read_bam_insert_metrics(bamfile=..., ...)` does not forward `...` to `readBam()`. Read the BAM explicitly, then pass its fragment object to the native histogram function:

```r
library(cfDNAPro)

# Select all BAM contigs. For the autosome comparison, supply their actual names
# instead, for example paste0("chr", 1:22) for a UCSC-style human BAM
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
    galp_tag = character(),
    outdir = NA
)

# Return the native histogram without saving a large per-fragment RDS file
counts <- read_bam_insert_metrics(
    fragment_obj = fragments,
    isize_min = 30L,
    isize_max = 1000L
)
write.table(counts, "lengths.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

The runner times the complete R process, including package loading, BAM reading, pairing, curation and histogram writing. Timing only the second function would exclude most of the work. Record peak memory as well, since cfDNAPro materializes the paired alignments and fragments. A speed difference must be measured rather than inferred from R versus Rust.

Keep the native cfDNAPro filters visible. Its default `scanBamFlag` excludes unmapped, secondary, supplementary and duplicate records, but does not explicitly exclude QC-fail records. The benchmark supplies `isNotPassingQualityControls=FALSE` and `hasUnmappedMate=FALSE` through `readBam(galp_flag=...)` to match cfDNAlab's exclusions. NM/MD tags are not loaded for this length-only calculation. Its strict 1-based 5'-end comparison is equivalent to `forward.pos < reverse.reference_end - 1`, whereas cfDNAlab additionally requires `forward.pos <= reverse.pos`. Overlapping pairs with the reverse alignment starting first can therefore remain a source of differences even with matched flags. GenomicAlignments pairing and reference-bound checks can also affect selection.

## How much does the fragment definition matter in this sample?

Include `cfdnapro-outer-span` in the concordance pilot. This uses the same cfDNAPro version, MAPQ, flag selection and pairing, changing only `curate_start_and_end` to `FALSE`. The native uncurated `granges(galp)` uses the outer alignment span, from the minimum mate start to the maximum mate end. This is a diagnostic configuration of the same tool, not another competing tool or version.

Compare `cfdnapro` against `cfdnapro-outer-span` to measure the effect on this sample's histogram. Report total retained counts, frequency errors and correlation, mean length, mode, the fraction <=150 bp, and Wasserstein distance in bp. A small Wasserstein distance and small changes in those summaries indicate a small distribution-level effect in this sample. The result does not establish that the definitions are interchangeable in other samples.

Length inclusion is applied after choosing the span, so a changed definition can move fragments into or out of the requested range. Native reference-bound checks also run after curation. The comparison therefore reports count changes as well as normalized distributions. It does not infer the fraction of individual fragments whose lengths changed from aggregate histograms.

The outer-span diagnostic does not test whether the aligner's stored TLEN equals that outer span. The samtools, Picard and FinaleToolkit comparisons assess the observed distributions from TLEN-based tools, but also include their selection differences. Do not attribute their entire discrepancy to the fragment definition.

For ordinary inward bwa mem alignments, stored TLEN spans the forward `pos` through the reverse `reference_end`, matching the directional aligned length. This is not a guarantee for arbitrary aligners or abnormal pairs. [bwa 0.7.19 source](https://github.com/lh3/bwa/blob/v0.7.19/bwamem.c#L887-L892) supplies the inspected rule. Record the actual remapping version. samtools inward selection uses forward/reverse start order, allowing coincident starts. Picard's FR rule instead compares 5' positions with <=, equivalent to `forward.pos < reverse.reference_end`. Its one-base boundary difference from cfDNAPro is outside this benchmark's >=10 bp range.

## Other tools and timing scope

samtools `stats` has no dedicated MAPQ switch. Its `-q` controls base-quality trimming. The runner uses HTSlib's general `--input-fmt-option 'filter=mapq >= 30'` instead, applying MAPQ in the same BAM read pass. Select the `IS` inward-pair column after timing. This does not establish that both mates pass MAPQ or reproduce cfDNAlab pairing. Label its runtime as full QC, since the histogram is part of a larger native calculation. Do not present it as a dedicated histogram kernel benchmark.

Picard has no MAPQ setting for this command, so the runner streams an uncompressed, filtered BAM from samtools into its native command. Both processes share the allocated CPU mask. The complete pipeline remains timed. The JVM sees the full core budget through `-XX:ActiveProcessorCount`, with a 4 GiB heap limit. Histogram width and orientation settings above are essential to preserving the requested length range. PDF plotting is omitted.

deepTools `bamPEFragmentSize` is lower priority. It samples genomic regions by default. `--distanceBetweenBins 0 --outRawFragmentLengths lengths.tsv` removes the sampling gaps, but its indexed queries can count alignments crossing chunk boundaries more than once. Its native selection also lacks the MAPQ/duplicate/QC-fail filtering needed here and falls back to other mapped reads in chunks without proper read1 records. These differences make it a less direct full-histogram comparison.

Use the established exclusive-node, local-scratch timing procedure. Use **1, 4, 12 and 32 cores in separate jobs** for parallel methods. First measure cfDNAPro and FinaleToolkit lengths at 1 versus 4 cores, also in separate jobs. If wall time and CPU use show no consistent benefit, stop higher-core measurements for this workload. If they improve, add 12/32-core jobs. The runner executes every explicitly selected method and does not automatically skip tools. cfDNAlab's `--n-threads` uses that budget. samtools `-@` receives `cores - 1`, since it counts additional decoding threads. Numerical libraries and the JVM may use the allocated budget, rather than being forced to a single thread. The timing TSV records `NA` for a command without an explicit native thread setting.

Plot wall time against allocated cores for every method. Calculate effective cores as `(user_seconds + system_seconds) / wall_seconds` and allocated CPU utilization as that value divided by `cores`. This shows actual utilization without assuming that a tool will or will not scale. Picard's utilization includes its samtools input filter, and cfDNAPro's includes container startup. Include native output generation in timing and perform concordance calculations afterward.

Compare raw totals and per-bin count differences as well as normalized frequencies. Report Pearson r and MAE/RMSE of the frequencies. High correlation alone can conceal a systematic count difference. Supply zero counts for absent length bins before comparing.

## Install and run

Use the shared [../install_tools.sh](../install_tools.sh). It installs the authors' cfDNAPro 1.7.4 Apptainer image by immutable digest, including its R dependencies and reference packages. Picard 3.5.0 uses Bioconda's `picard-slim` package with Java 17 under `envs/picard`. It contains the official JAR but omits R, which this command only needs for its optional PDF. samtools 1.25 uses a separate Bioconda environment under `envs/samtools`, packed as `samtools.tar.gz` with its shared libraries. The runner copies the image and packed environments locally. It records the actual cfDNAPro package version and image SHA-256 in each timing row. Apptainer must already be available on the cluster. See the [installation audit](../README.md#installation) for sources and deliberate adaptations.

```bash
bash benchmarks/install_tools.sh /path/to/benchmark-tools

export BENCH_BAM=/path/to/BH01.hg38.paired.15x.bam
export BENCH_CFDNA=/path/to/benchmark-tools/envs/cfdnalab/bin/cfdna
export BENCH_SAMTOOLS_ARCHIVE=/path/to/benchmark-tools/samtools.tar.gz
export BENCH_FINALE_ARCHIVE=/path/to/benchmark-tools/finaletoolkit.tar.gz
export BENCH_CFDNAPRO_IMAGE=/path/to/benchmark-tools/cfdnapro-1.7.4.sif
export BENCH_PICARD_ARCHIVE=/path/to/benchmark-tools/picard.tar.gz
export BENCH_LOCAL_ROOT=/path/to/node-local-scratch
export BENCH_GENOME=hg38

# Run the four comparisons and the outer-span/proper-pair diagnostics once
export BENCH_RESULTS=/path/to/lengths-pilot
export BENCH_CORES="1"
export BENCH_PILOT=1
export BENCH_METHODS="cfdnalab cfdnalab-autosomes cfdnalab-proper cfdnapro cfdnapro-outer-span finaletoolkit samtools picard"
sbatch --constraint=gen6 --time=YOUR_PILOT_BASED_BUDGET benchmarks/lengths/benchmark_lengths.sh

# After the pilot finishes, compare all saved native histograms
/path/to/benchmark-tools/envs/comparison/bin/python \
    benchmarks/lengths/compare_lengths.py "$BENCH_RESULTS/profiles" \
    > "$BENCH_RESULTS/concordance.tsv"

# Submit a separate job and results directory for each core setting and method block
export BENCH_RESULTS=/path/to/lengths-timing-4cores-parallel
export BENCH_CORES=4
export BENCH_PILOT=0
export BENCH_WARMUPS=2
export BENCH_REPETITIONS=5
export BENCH_METHODS="cfdnalab cfdnalab-autosomes samtools picard"
sbatch --constraint=gen6 --time=YOUR_PILOT_BASED_BUDGET benchmarks/lengths/benchmark_lengths.sh
```

Set wall time from pilots at each core setting. Each job uses a single core count. Split methods further as needed to fit the allocation limit, using matched hardware as described in [the shared run controls](../README.md#runs-and-budgets). The example times the other length methods at 4 cores. Submit cfDNAPro and FinaleToolkit separately at 1 and 4 cores, and decide on higher-core runs after inspecting their measurements. `BENCH_INDEX` can override the usual `.bam.bai` path. `BENCH_MIN_LENGTH`, `BENCH_MAX_LENGTH` and `BENCH_MIN_MAPQ` override the common filters. Supply the same nondefault length bounds to the comparator. `BENCH_GENOME` selects cfDNAPro's reference label and must match the alignment reference and chromosome naming. The default `hg38` means UCSC-style names. The input remains a coordinate-sorted BAM with its original non-autosomes and alignment flags.

The runner reserves an exclusive node and all its memory. Every sequential `srun --jobid` attaches to that allocation and checks the actual hostname and logical/physical CPU IDs. Hard confinement also depends on the cluster's Slurm cgroup configuration. Prefer the same physical node for split jobs when available, as described in the shared run controls. Each method runs as a contiguous block at the job's core count, retaining actual order in `timings.tsv`. There are two warmups followed by five measured repetitions by default. No caches are evicted.

The timing result is a single `timings.tsv`, with a row for every run, including warmups and failures. A native histogram per method is retained under `profiles/` from the first measured repetition at the capture core budget. These small files support concordance analysis. Other temporary outputs and logs are deleted. Failed commands are printed in the Slurm log, and a failed timing row is saved before stopping. Input staging, scheduler startup, profile copying, cleanup and concordance calculation are outside the timer. Native startup and output generation remain timed, including cfDNAlab's default QC plot and cfDNAPro's Apptainer launch. GNU time's maximum RSS is not the simultaneous total memory of a process tree.

## Supplementary result placeholders

| cfDNAlab vs tool | Pearson r | MAE / RMSE | Count difference | Wasserstein, bp | Main expected differences |
| --- | --- | --- | --- | --- | --- |
| cfDNAPro 1.7.4 | [ ] | [ ] / [ ] | [ ] | [ ] | Pair selection and overlapping outward-pair checks. |
| FinaleToolkit 1.1.0 | [ ] | [ ] / [ ] | [ ] | [ ] | Proper-pair selection, read1 MAPQ and TLEN. |
| samtools 1.25 | [ ] | [ ] / [ ] | [ ] | [ ] | TLEN, per-record filters and halved counts. |
| Picard 3.5.0 | [ ] | [ ] / [ ] | [ ] | [ ] | TLEN from read2 and orientation classification. |

All rows use whole-reference counts, MAPQ >=30 and lengths 30-1000 bp. MAE/RMSE use within-range frequencies. Count differences are cfDNAlab minus comparator. Interpret the separately labeled autosome run as a different scope.

| Span definition in cfDNAPro | Retained fragments | Mean, bp | Mode, bp | Fraction <=150 bp |
| --- | --- | --- | --- | --- |
| Directional, curated | [ ] | [ ] | [ ] | [ ] |
| Outer alignment span | [ ] | [ ] | [ ] | [ ] |

Report the definition diagnostic's Pearson r, frequency MAE/RMSE and Wasserstein distance alongside these rows. Both configurations use the same native pair/flag/MAPQ settings. Scope the conclusion to this sample and length range.

| Tool / configuration | Cores | Warmups / measured | Mean process+flush, s | SD seconds | Effective cores |
| --- | --- | --- | --- | --- | --- |
| cfDNAlab, all contigs | 1 / 4 / 8 / 12 / 32 | 2 / 5 each | [ ] | [ ] | [ ] |
| cfDNAlab, autosomes | 1 / 4 / 8 / 12 / 32 | 2 / 5 each | [ ] | [ ] | [ ] |
| cfDNAPro 1.7.4 | 1 / 32 | 2 / 5 each | [ ] | [ ] | [ ] |
| FinaleToolkit 1.1.0 | 1 / 32 | 2 / 5 each | [ ] | [ ] | [ ] |
| samtools 1.25 | 1 / 4 / 8 / 12 / 32 | 2 / 5 each | [ ] | [ ] | [ ] |
| Picard 3.5.0 | 1 / 4 / 8 / 12 / 32 | 2 / 5 each | [ ] | [ ] | [ ] |

Expand timing rows to a core budget per row in the final A4 supplement. Put workload notes below the table and report maximum RSS separately if needed.

## Source references

- cfDNAlab: [length construction](../../src/shared/fragment/indel_counting_fragment.rs), [orientation checks](../../src/shared/fragment/minimal_fragment.rs), [read filtering](../../src/shared/read.rs), [CLI settings](../../src/commands/lengths/config.rs).
- cfDNAPro: [readBam](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/readBam.R), [curation and orientation](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/internalFunctions.R), [histogram](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/read_bam_insert_metrics.R).
- FinaleToolkit: [histogram implementation](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/frag/_frag_length.py), [alignment reader](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/io/alignment.py).
- samtools: [stats manual](https://www.htslib.org/doc/samtools-stats.html), [1.25 counting and output](https://github.com/samtools/samtools/blob/1.25/stats.c), [HTSlib input filter parsing](https://github.com/samtools/htslib/blob/1.25/hts.c), [filtering during BAM reads](https://github.com/samtools/htslib/blob/1.25/sam.c).
- Picard: [3.5.0 release](https://github.com/broadinstitute/picard/releases/tag/3.5.0), [command](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/CollectInsertSizeMetrics.java), [collector](https://github.com/broadinstitute/picard/blob/3.5.0/src/main/java/picard/analysis/directed/InsertSizeMetricsCollector.java).
- deepTools: [4.0.0 fragment sampling and filtering](https://github.com/deeptools/deepTools/blob/4.0.0/pydeeptools/deeptools/getFragmentAndReadSize.py).
- Outer-span diagnostic: [cfDNAPro's curation toggle](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/readBam.R), [GenomicAlignments range union](https://github.com/Bioconductor/GenomicAlignments/blob/devel/R/GAlignmentPairs-class.R).

Join completed Slurm step memory after the job. Timing rows also record filesystem flush time, process-plus-flush time and output bytes. See the [shared timing TSV documentation](../README.md) for units and measurement limits. Submit from the repository root so helper scripts can be staged before warmups.
