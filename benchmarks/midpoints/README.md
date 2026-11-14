# Midpoint profile concordance and speed

Compare `cfdna midpoints` with Griffin's native BAM-to-profile workflow on the same sample and grouped sites. Use the shared [../install_tools.sh](../install_tools.sh). The benchmark uses unmodified Griffin v0.2.0 scripts, a release marked as a prerelease upstream.

## Comparable settings

The benchmark uses the same grouped sites, inclusive fragment lengths 100-200 bp and MAPQ >=20 for both profilers. Both sample GC estimators use MAPQ >=20 and the native estimation ranges described below. Supply 9,990 bp BED windows (`center - 4995`, `center + 4995`). cfDNAlab saves all 666 bins of 15 bp. Griffin normalizes over that same window but saves its official central 132 bins, `[-990, 990)`. Both use native 165 bp order-3 Savitzky-Golay smoothing. No per-site CNA correction or fragment-count weights are added.

| Tool | Settings | Expected differences |
| --- | --- | --- |
| cfDNAlab | `midpoints --length-bins 100 201 --min-mapq 20 --bin-size 15`, default smoothing, shared exclusion BEDs | Uses paired spans from `forward.pos` to `reverse.reference_end`, checks both mates and applies native blacklist filtering. Even-length fragments reproducibly choose either central base. Writes group sums averaged across the bases in each bin. |
| Griffin v0.2.0 | `griffin_coverage.py`, then `griffin_merge_sites.py`, native exclusion BEDs, zero-mappability masking, outlier removal and 165 bp smoothing | Uses forward records with positive TLEN and the right central base for even lengths. Averages sites and normalizes the combined group profile to mean one. |

Griffin receives a centered anchor per site. Its normalization window is `[-4995, 4995)`, matching cfDNAlab's complete BED window. Its official `--save_window -1000 1000` rounds inward to `[-990, 990)` at the 15 bp step. Saving the central window excludes Griffin's polynomial smoothing edge bins and does not change its full-window normalization. The runner never resizes the supplied BED. Misaligned bin coordinates fail comparison.

Griffin keeps its official center feature window `[-30, 30)`, FFT window `[-960, 960)` and FFT index 10. FFT feature calculation is compulsory and remains timed. Mappability correction is off, while zero-mappability exclusion is on. Both tools use the guide's 15 bp step and MAPQ 20. Griffin's final plotting rule is omitted. Its GC-bias plots and cfDNAlab's native QC plots remain timed.

The comparator divides cfDNAlab's complete 666-bin group profile by its full-window mean, then selects the same 132 central bins as Griffin. It never renormalizes the saved subset. Griffin's cropped output therefore need not have mean one. cfDNAlab averages per-bp counts within bins, while Griffin sums them. This constant factor cancels during normalization. cfDNAlab smooths before binning and Griffin bins before smoothing over 11 bins. The comparator adds no smoothing, shifting or rebinning. It reports signed Pearson r, r squared, MAE/RMSE and maximum absolute difference.

The same five public exclusion BEDs go through each tool's native implementation. cfDNAlab can exclude whole intervals near blacklists, while Griffin masks individual positions and removes outliers. Griffin also uses its Umap zero-mappability mask. These are native algorithm differences, not forced-equivalent filters. cfDNAlab's eligible interval count and Griffin's input site-list count are reported separately. Griffin's count includes masked sites. Finite negative smoothing undershoots are retained. If Griffin leaves undefined positions, metrics for that group are unavailable instead of silently comparing a subset. Constant-profile correlations are also unavailable. Griffin rounds saved values to five decimal places.

The source also shows filtering differences. Griffin checks the forward record's MAPQ, duplicate and QC-fail flags and paired status. It does not explicitly reject secondary/supplementary records or validate the actual mate's MAPQ, chromosome and orientation. cfDNAlab performs these checks. Neither baseline requires flag 0x2. The source-specific behavior stays in the comparison.

The preparation utility assigns positive TLEN to the reconstructed forward mate. The shared bwa mem remap recomputes TLEN and flags, so final alignments determine inclusion. A reverse record with positive TLEN would be omitted by Griffin, but that is not an expected reconstruction difference after this remap. See [input provenance](../data/README.md).

Griffin pads queried regions by the **minimum** fragment length. The 100-200 bp baseline provides at least half the maximum length as padding. The sufficient padding condition is `floor(max_length / 2) <= min_length`. Fetch losses become possible from `max_length >= 2 * min_length + 2`, depending on read alignment extent. Group aggregation parallelizes across site lists, so a single group cannot use all aggregation workers even when counting scales.

Griffin holds per-base dictionaries across merged fetch regions and dense site-by-position matrices during aggregation. Large or heavily overlapping site sets can therefore require substantial memory.

The production source is Griffin's [public paper TFBS set](https://github.com/adoebley/Griffin_analyses/tree/75467e438a88029206d8d778ddce05eb0ab6efa8/sites/TFBS/10000_unfiltered_sites_CIS_BP_v2), pinned to the linked revision. Supply its prepared BED windows rather than selecting sites in the runner. Record the actual site/group counts, window width and derivation with the results.

Sources: [cfDNAlab counting](../../src/commands/midpoints/midpoints.rs), [midpoint choices](../../src/shared/midpoint.rs), [window orientation](../../src/commands/midpoints/strand.rs), [Griffin coverage](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_coverage.py), [aggregation](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_merge_sites.py), [site handling](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_functions.py).

## Griffin workflow and retained profiles

Griffin requires a reference FASTA and a GC-bias table even when using its **uncorrected** output. Its native counting script always reads reference fragment sequences, calculates GC content and writes both raw and GC-weighted bigWigs. There is no CLI switch to disable that work.

The timed Griffin workload is `griffin-gc`: sample GC counting, bias fitting and both profiling scripts. It writes both raw and corrected profile tables, which are retained from its first measured repetition. The raw table is compared with uncorrected `cfdnalab`; the corrected table is compared with `cfdnalab-gc`. There is no separate unit-weight production Griffin job. A unit GC table is used only by the short untimed demo warmup.

## GC-corrected comparison

Run `cfdnalab-gc` and `griffin-gc` in separate jobs for the timing comparison. Each timed repetition estimates sample GC bias and produces corrected profiles. Only reference GC distributions are precomputed. Run uncorrected `cfdnalab` separately if raw-profile concordance is also wanted, using the uncorrected table already saved by `griffin-gc`.

| Workload | Commands inside each timing |
| --- | --- |
| `cfdnalab-gc` | `cfdna gc-bias`, then `cfdna midpoints --gc-file ... --ref-2bit ...` |
| `griffin-gc` | Native `griffin_GC_counts.py --size_range 15 500`, `griffin_GC_bias.py --size_range 15 500`, then `griffin_coverage.py` with the profiling range and `griffin_merge_sites.py`. Retains `sample.GC_corrected.coverage.tsv` |

Both workflows keep native bias fitting, smoothing, exclusion rules and diagnostic plots. Correction algorithms are not forced to match. cfDNAlab uses directional paired spans, GC percentages and its native reference end offset (10 bp by default). Griffin uses record/TLEN selection, whole-fragment GC counts and its native smoothing. Its estimator counts qualifying forward and reverse records independently. Corrected-profile differences therefore reflect both bias estimates and midpoint counting. Mean per-site/per-bp values represent weighted midpoint counts in this comparison.

cfDNAlab uses its default reference and sample GC-estimation settings. Its reference package is calculated once per alignment reference and reused across all BAMs, repetitions and core settings. The default package covers 30-1000 bp, and sample GC estimation inherits that range. The midpoint profiles use the 100-200 bp range.

Griffin's downloaded source already includes its hg38 GC-estimation BED and precomputed reference GC-frequency tables under `Griffin-source/Ref/`. Use these bundled inputs for the hg38 benchmark. Its estimator uses `chr1` through `chr22` and handles the supplied estimation regions itself. This BED is only used by `griffin-gc`. It does not replace or modify `BENCH_INTERVALS` or restrict cfDNAlab's reference calculation.

Griffin's GC estimation uses its release's [15-500 bp range](https://github.com/adoebley/Griffin/blob/v0.2.0/snakemakes/griffin_GC_and_mappability_correction/config/config.yaml), independently of 100-200 bp profiling. Setup stages every reference table from 15 through 500 bp and records the bounds. Restricting estimation to the profiling range would remove records and neighbouring lengths used in native smoothing. Both estimation paths and both profilers use MAPQ >=20 for this comparison.

Add these paths to the local settings file, in addition to the existing BAM, site, FASTA and tool paths. Reuse `BENCH_REF_2BIT` if it is already exported for end motifs:

```bash
# cfDNAlab
export BENCH_REF_2BIT=/path/to/the-alignment-reference.2bit
export BENCH_CFDNA_REF_GC=/path/to/cfdnalab-reference-gc/ref_gc_package.zarr
# Griffin
export BENCH_GRIFFIN_GC_REGIONS=/path/to/benchmark-tools/Griffin-source/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed
export BENCH_GRIFFIN_REF_GC=/path/to/benchmark-tools/Griffin-source/Ref/genome_GC_frequency
```

The FASTA and 2bit must contain the same alignment reference. `BENCH_CFDNA_REF_GC` points to the reusable default reference package. Native loading checks reference compatibility, and midpoint correction checks that the requested length range is covered. Griffin expects a file named `<GC_BED_basename_without_extension>.<length>bp.GC_frequency.txt` for every length from 15 through 500 bp in `BENCH_GRIFFIN_REF_GC`. The pinned release includes all of these tables. They must correspond to the supplied reference, estimation BED and read length. A supplied reference table from another assembly or region set is not interchangeable.

Run cfDNAlab's reference calculation once in a compute allocation, outside benchmark measurements. The following command supplies only the required paths, leaving all calculation settings at their defaults. Keep this command outside the settings file, which should only export paths:

```bash
"$BENCH_CFDNA" ref-gc-bias \
    --ref-2bit "$BENCH_REF_2BIT" \
    --output-dir "$(dirname "$BENCH_CFDNA_REF_GC")"
```

Griffin's bundled hg38 reference tables are already downloaded by the installer, which prints both Griffin input exports. No separate Griffin reference-generation job is needed for these inputs. Its native [reference generator](https://github.com/adoebley/Griffin/blob/v0.2.0/scripts/griffin_calc_GC_frequency.py) is only needed for a different reference or estimation region set. It runs once per fragment length with an explicit read length. Its chromosome-boundary branches contain an undefined `chrom_sizes_dict` variable, which can make custom reference preparation fail upstream. The bundled tables avoid rerunning that generator.

The Griffin archive must contain its GC scripts and matplotlib for native bias plots. If those files are missing, rerun `bash benchmarks/install_tools.sh ROOT --only griffin`. All installation files remain under ROOT.

With those references prepared, use the commented [production run procedure](../README.md#runs-and-budgets). It submits uncorrected cfDNAlab, GC-corrected cfDNAlab and GC-corrected Griffin as three separate midpoint jobs at each of 1, 4, 12 and 32 cores. Each midpoint job has its own result directory and an explicitly named midpoint wall-time variable.

Griffin keeps its short, untimed demo with unit weights, without estimating sample bias during that demo. Every measured `griffin-gc` repetition performs the full sample estimation and production profiling. Its first measured estimation has no preceding estimation warmup in that job. Full `cfdnalab-gc` runs warm up the entire corrected workflow. The first measured bias estimate and its diagnostics are retained under `profiles/<method>/sample-gc-bias/`. Inputs are staged before timing, and both `run.tsv` and the timing rows distinguish corrected runs from the uncorrected baseline.

Budget corrected Griffin jobs for all three complete estimation-plus-profiling runs, plus the demo and setup. Its [release cluster configuration](https://github.com/adoebley/Griffin/blob/v0.2.0/snakemakes/griffin_GC_and_mappability_correction/config/cluster_slurm.yaml) requests 48 hours with eight CPUs for GC counting alone and warns that large BAMs can take days. Those are upstream resource requests, not measured runtimes for this input. No complete-job time limit is established here.

The shared [completed-job commands](../README.md#completed-job-accounting-and-concordance) compare both raw and corrected profiles at matching core settings. `--gc-corrected` selects `cfdnalab-gc` and Griffin's corrected table. The default comparator pairs uncorrected `cfdnalab` with `profiles/griffin-gc/sample.uncorrected.coverage.tsv`. In both comparisons, cfDNAlab is normalized over the full window before selecting the shared central bins. Griffin's saved profiles already use full-window normalization.

## Inputs and run

The same exclusion BEDs are used by both tools. Point `BENCH_GRIFFIN_REF_DIR` to the existing installation's `Griffin-source/Ref` directory. It contains the five BEDs used by Griffin's release config. Griffin additionally needs the hg38 k100 Umap bigWig. Download that reference once outside timed jobs:

```bash
export BENCH_GRIFFIN_REF_DIR=/path/to/benchmark-tools/Griffin-source/Ref
export BENCH_GRIFFIN_MAPPABILITY_BW="$BENCH_GRIFFIN_REF_DIR/k100.Umap.MultiTrackMappability.bw"
wget -O "$BENCH_GRIFFIN_MAPPABILITY_BW" \
    https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw
```

Keep the two exports in the local settings file. Downloading is a separate preparation command. No environment reinstall is required for an existing Griffin installation. The launcher checks reference-file existence before submission. Benchmark jobs copy these references to node-local scratch before timing.

Supply the prepared coordinate-sorted BAM and its index, the matching **uncompressed FASTA plus existing `.fai`** for Griffin, and your BED-like file with `chromosome, start, end, group` through `BENCH_INTERVALS`. Use tab-separated rows without a header. BED4 is unstranded. Optional strands follow cfDNAlab's BED5/BED6 convention. Supply equal-width windows within their chromosomes. For production, use BED `[center-4995, center+4995)`, 9,990 bp total, to align both tools' 15 bp bins. Reverse-stranded windows are reversed by both tools. The runner retains all supplied site rows, including repeats and overlaps, without selecting or resizing sites. Before submission, the launcher checks only the first five lines for the expected columns and integer coordinates. Input validation belongs to the native tools.

Use the shared [production run procedure](../README.md#runs-and-budgets) for the exported paths, all three workloads, fresh result directories and cross-job comparisons. Each midpoint job uses a single core setting and method. `BENCH_INDEX` and `BENCH_FASTA_INDEX` can override the usual index paths. `BENCH_MIN_LENGTH`, `BENCH_MAX_LENGTH` and `BENCH_MIN_MAPQ` override the profiling filters. Keep the same production BED and settings across the jobs.

Griffin receives a short, untimed warmup on the release's [1000-site CTCF demo](https://github.com/adoebley/Griffin/blob/v0.2.0/demo/griffin_nucleosome_profiling_demo_files/sites/CTCF.hg38.1000.txt), using the same staged production BAM/reference and core budget. It runs the same two native scripts and production window/filter settings on those demo coordinates in a separate Slurm step. This gives Griffin a practical warmup without another full production run. Demo timings and memory do not enter measured results. Logs and the demo checksum are retained, and demo outputs are removed before measurements. The supplied production sites are unchanged.

The production cfDNAlab schedule is a full-input warmup plus three measurements. Griffin has no full-input warmup and performs three measured runs after its demo. Every new Griffin job receives the demo.

For full runs, the submission launcher reserves an exclusive node and all its memory. Sequential `srun --jobid` steps use that allocation and check their actual node and allowed logical/physical CPU IDs. Hard confinement also depends on Slurm's cgroup setup. All inputs and Griffin's environment are staged on local scratch before warmups. No caches are evicted. Both tools profile the chromosomes present in the supplied site set, while the BAM retains its other chromosomes and records.

`timings.tsv` contains a row for every attempted full-input run, including full-input warmups and failures, with `gc_correction=none` or `corrected`. Griffin's demo is untimed, with logs and its checksum retained separately. A native profile is saved under `profiles/<method>/` from the first measured repetition at the job's single core count. Corrected workloads also retain that repetition's sample bias and diagnostics. Other temporary files and logs are deleted. Failed diagnostics are retained and the remaining repetitions of that method are skipped. See the [shared failure and timeout behavior](../README.md#timing-and-slurm-memory). Source files are not modified or downloaded by benchmark jobs.

Timing includes native startup, counting, output writing and Griffin's bigWig-to-profile aggregation. cfDNAlab's native default QC plot also remains timed. Input staging, scheduler startup, profile preservation, cleanup and concordance calculations are outside the timer. GNU time's maximum RSS is not the simultaneous total memory of a multiprocessing tree.

## Supplementary result templates

The tables below contain no measured results. Fill them from the saved concordance metrics and completed timing/accounting records.

| Group | cfDNAlab eligible sites / Griffin listed sites | Pearson r / r squared, 15 bp | MAE / RMSE, mean-one | Max. difference | Mean midpoint density, cfDNAlab / Griffin |
| --- | --- | --- | --- | --- | --- |
| [group] | [ ] | [ ] | [ ] / [ ] | [ ] | [ ] / [ ] |

Expected differences: even-length midpoint choices, actual paired spans versus TLEN, filtering, and Griffin's five-decimal output. Coordinate axes must agree. The cfDNAlab count is eligible intervals. The Griffin count is sites in its selected site list, including masked sites.

| Tool | Cores | Warmups / measured | Mean process+flush, s | SD seconds | Slurm peak, KiB |
| --- | --- | --- | --- | --- | --- |
| cfDNAlab | 1 / 4 / 12 / 32, separate jobs | 1 / 3 each | [ ] | [ ] | [ ] |
| Griffin v0.2.0 | 1 / 4 / 12 / 32, separate jobs | 1000-site demo / 3 full-input each | [ ] | [ ] | [ ] |

Expand the timing rows to a core budget per row for the final A4 supplement. Report the shared site set, groups, window width, 15 bp resolution and 100-200 bp length range alongside the table, and place the workload notes below it.

Join completed Slurm step memory after the job. Timing rows also record filesystem flush time, process-plus-flush time and output bytes. See the [shared timing TSV documentation](../README.md) for units and measurement limits. Submit from the repository root so helper scripts can be staged before warmups.
