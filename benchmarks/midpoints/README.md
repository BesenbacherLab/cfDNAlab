# Midpoint profile concordance and speed

Compare `cfdna midpoints` with Griffin's native BAM-to-profile workflow on the same sample and grouped sites. Use the shared [../install_tools.sh](../install_tools.sh). Griffin v0.2.0 is its newest published version, marked as a prerelease upstream. The coverage and aggregation scripts in the current main branch are identical to that release. No benchmark data have been run while preparing these scripts.

## Comparable settings

The baseline uses unsmoothed 1 bp profiles of fragments 100-200 bp inclusive, MAPQ >=30, and no blacklist, mappability mask, outlier removal, GC correction or genomic scaling. This length range follows Griffin's nucleosome-profile example. Both tools retain their native alignment filtering and pairing behavior. Proper-pair selection is not added to cfDNAlab.

| Tool | Settings | Expected differences |
| --- | --- | --- |
| cfDNAlab | `midpoints --length-bins 100 201 --min-mapq 30 --bin-size 1 --smoothing none`, shared intervals and their chromosomes | Uses actual paired spans from `forward.pos` to `reverse.reference_end` and checks both mates. Even-length fragments reproducibly choose either central base. Writes summed counts across each group. |
| Griffin v0.2.0 | `griffin_coverage.py --size_range 100 200 --map_quality 30`, then `griffin_merge_sites.py --step 1 --smoothing False --CNA_normalization False --individual False --exclude_paths none --exclude_outliers False --exclude_zero_mappability False` | Uses forward records with positive TLEN. Always uses the right central base for even lengths. Averages sites and always normalizes the group profile to mean one. |

The comparator divides each cfDNAlab group profile by its mean over the complete window. This equals Griffin's group normalization when site sets and counts agree. It does not normalize each site separately, change midpoint choices, smooth, shift or rebin the profiles. All zero positions remain. Check Pearson r together with MAE/RMSE and each tool's mean midpoints per site per bp, since normalization can conceal differences in total counts. Griffin rounds its saved values and mean-count metadata to five decimal places.

The source also shows filtering differences. Griffin checks the forward record's MAPQ, duplicate and QC-fail flags and paired status. It does not explicitly reject secondary/supplementary records or validate the actual mate's MAPQ, chromosome and orientation. cfDNAlab performs these checks. Neither baseline requires flag 0x2. The source-specific behavior stays in the comparison.

The input preparation now assigns positive TLEN to the reconstructed forward mate. The shared bwa mem remap recomputes TLEN and flags, so final alignments determine inclusion. A reverse record with positive TLEN would be omitted by Griffin, but that is not an expected reconstruction difference after this remap. See [input provenance](../data/README.md).

Griffin pads queried regions by the **minimum** fragment length. The 100-200 bp baseline provides at least half the maximum length as padding. The sufficient padding condition is `floor(max_length / 2) <= min_length`. Fetch losses become possible from `max_length >= 2 * min_length + 2`, depending on read alignment extent. Group aggregation parallelizes across site lists, so a single group cannot use all aggregation workers even when counting scales.

Pilot the actual site set before scheduling all repetitions. Griffin holds per-base dictionaries across merged fetch regions and dense site-by-position matrices during aggregation. Large or heavily overlapping site sets can therefore require substantial memory.

Sources: [cfDNAlab counting](../../src/commands/midpoints/midpoints.rs), [midpoint choices](../../src/shared/midpoint.rs), [window orientation](../../src/commands/midpoints/strand.rs), [Griffin coverage](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_coverage.py), [aggregation](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py), [site handling](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_functions.py).

## Griffin's compulsory GC work

Griffin requires a reference FASTA and a GC-bias table even when using its **uncorrected** output. Its native counting script always reads reference fragment sequences, calculates GC content and writes both raw and GC-weighted bigWigs. There is no CLI switch to disable that work.

The runner supplies a small GC table with every bias set to 1. This gives unit weights without estimating a sample-specific correction. It compares the native `uncorrected.coverage.tsv`. Griffin's unavoidable reference reading, unit-weight calculation, duplicate weighted output and native feature calculations remain inside the timer. They are disclosed workload differences. No Griffin source is patched. This baseline does **not** validate the biological effect of GC correction. A corrected-profile comparison must use independently prepared correction inputs and report that the two tools estimate bias differently.

## Inputs and run

Supply the prepared coordinate-sorted BAM and its index, the matching **uncompressed FASTA plus existing `.fai`**, and your BED-like file with `chromosome, start, end, group` and optional strand columns through `BENCH_INTERVALS`. BED6 uses column 6 for strand, as in cfDNAlab. Use tab-separated rows without a header. All intervals must have the same even width and lie entirely inside their chromosomes. Each interval becomes the identical symmetric window in Griffin, centered at `start + width/2`. Reverse-stranded windows are reversed by both tools. The runner retains all supplied site rows, including repeats and overlaps, and fails on invalid input instead of silently dropping sites. All-numeric chromosome names and invalid strand tokens fail before Griffin timing.

```bash
# The shared installer prints these tool paths
export BENCH_CFDNA=/path/to/benchmark-tools/envs/cfdnalab/bin/cfdna
export BENCH_GRIFFIN_ARCHIVE=/path/to/benchmark-tools/griffin.tar.gz
export BENCH_BAM=/path/to/BH01.hg38.paired.15x.bam
export BENCH_REF_FASTA=/path/to/the-alignment-reference.fa
export BENCH_INTERVALS=/path/to/grouped-midpoint-windows.bed
export BENCH_SAMTOOLS_ARCHIVE=/path/to/benchmark-tools/samtools.tar.gz
export BENCH_LOCAL_ROOT=/path/to/node-local-scratch

# Concordance pilot: run each tool once
export BENCH_RESULTS=/path/to/midpoint-pilot
export BENCH_CORES="1"
export BENCH_PILOT=1
sbatch --constraint=gen6 --time=YOUR_PILOT_BASED_BUDGET benchmarks/midpoints/benchmark_midpoints.sh

# After the job finishes, write a row per group with agreement metrics
/path/to/benchmark-tools/envs/comparison/bin/python \
    benchmarks/midpoints/compare_midpoints.py "$BENCH_RESULTS/profiles" \
    > "$BENCH_RESULTS/concordance.tsv"

# Submit a separate job and results directory for each core setting
export BENCH_RESULTS=/path/to/midpoint-timing-4cores
export BENCH_CORES=4
export BENCH_PILOT=0
export BENCH_WARMUPS=2
export BENCH_REPETITIONS=5
sbatch --constraint=gen6 --time=YOUR_PILOT_BASED_BUDGET benchmarks/midpoints/benchmark_midpoints.sh
```

Set the generation selector and wall time from the cluster settings and full-matrix pilot. Split by core budget and then by method as needed to fit the allocation limit, using matched hardware as described in [the shared run controls](../README.md#runs-and-budgets). `BENCH_INDEX` and `BENCH_FASTA_INDEX` can override the usual index paths. `BENCH_MIN_LENGTH`, `BENCH_MAX_LENGTH` and `BENCH_MIN_MAPQ` override the common filters. `BENCH_METHODS="cfdnalab griffin"` selects the two tools.

The runner reserves an exclusive node and all its memory. Sequential `srun --jobid` steps use that allocation and check their actual node and allowed logical/physical CPU IDs. Hard confinement also depends on Slurm's cgroup setup. All inputs and Griffin's environment are staged on local scratch before warmups. Methods run in contiguous blocks with warmups before measurements at the job's single core count. Submit separate jobs and results directories for the other core settings. No caches are evicted. Both tools use exactly the chromosomes present in the supplied site set, while the BAM retains its other chromosomes and records.

The only timing result is `timings.tsv`, with a row for every run, including warmups and failures. Native profiles are also retained under `profiles/cfdnalab/` and `profiles/griffin/` from the first measured repetition at the capture core budget. Other temporary files and logs are deleted. Errors are printed to the Slurm log. Source files are not modified or downloaded.

Timing includes native startup, counting, output writing and Griffin's bigWig-to-profile aggregation. cfDNAlab's native default QC plot also remains timed. Input staging, scheduler startup, profile preservation, cleanup and concordance calculations are outside the timer. GNU time's maximum RSS is not the simultaneous total memory of a multiprocessing tree.

## Other software checked

[deepTools 4.0.0](https://github.com/deeptools/deepTools/releases/tag/4.0.0) provides `bamCoverage --MNase`. Its current [native counting code](https://github.com/deeptools/deepTools/blob/4.0.0/src/filtering.rs) counts two central bases for even insert lengths and four for odd lengths. `--centerReads` centers a read-sized span rather than reducing a fragment to a 1 bp midpoint. These settings do not directly reproduce this single-midpoint definition. A wider-kernel comparison could be separately labeled, but it is not included in this baseline. See also the [current CLI](https://github.com/deeptools/deepTools/blob/4.0.0/pydeeptools/deeptools/bamCoverage2.py).

[LIQUORICE](https://www.nature.com/articles/s41467-021-23445-w) describes a fragment-coverage workflow rather than this 1 bp midpoint-count profile. [FinaleToolkit's regional coverage and WPS outputs](https://epifluidlab.github.io/FinaleToolkit/documentation/user_guide/features.html) are also different quantities. No 1 bp genomic intervals file is generated to force these tools into a midpoint comparison.

## Supplementary placeholders

| Group | Sites | Pearson r, 1 bp | MAE / RMSE, mean-one | Max. difference | Mean midpoint density, cfDNAlab / Griffin |
| --- | --- | --- | --- | --- | --- |
| [group] | [ ] | [ ] | [ ] / [ ] | [ ] | [ ] / [ ] |

Expected differences: even-length midpoint choices, actual paired spans versus TLEN, filtering, and Griffin's five-decimal output. Site counts and coordinate axes must agree before interpreting these metrics.

| Tool | Cores | Warmups / measured | Mean process+flush, s | SD seconds | Slurm peak, KiB |
| --- | --- | --- | --- | --- | --- |
| cfDNAlab | 1 / 4 / 8 / 12 / 32 | 2 / 5 each | [ ] | [ ] | [ ] |
| Griffin v0.2.0 | 1 / 4 / 8 / 12 / 32 | 2 / 5 each | [ ] | [ ] | [ ] |

Expand the timing rows to a core budget per row for the final A4 supplement. Report the shared site set, groups, window width, 1 bp resolution and 100-200 bp length range alongside the table, and place the workload notes below it.

Join completed Slurm step memory after the job. Timing rows also record filesystem flush time, process-plus-flush time and output bytes. See the [shared timing TSV documentation](../README.md) for units and measurement limits. Submit from the repository root so helper scripts can be staged before warmups.
