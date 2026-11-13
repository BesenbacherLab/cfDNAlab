# Benchmark methods and concordance specification

This reference describes native tool configurations, counting differences, output formats and concordance metrics. Installation, execution, accounting and inputs are specified in [benchmarks_spec.md](benchmarks_spec.md).

Fragment coverage, read depth, regional fragment counts, global length histograms, end-motif frequencies and midpoint profiles are separate quantities. Normalized agreement does not establish equal total counts or individual fragment selection. Speed claims apply to the configurations actually timed, including their native compulsory work and filter differences.

## Shared scientific settings

cfDNAlab paired spans remain directional: `forward.pos` to `reverse.reference_end`. Competitors retain their native pairing, TLEN and CIGAR rules. Proper-pair diagnostics add flag 0x2 selection without making all inclusion rules identical.

Normal settings request MAPQ >=30 and inclusive fragment lengths 30-1000 bp. Midpoints uses 100-200 bp and MAPQ >=20 for profiling and sample GC estimation. cfDNAlab length bins end at `maximum + 1`. Global scope is all BAM contigs unless labeled autosomes or cfDNAPro motif reference compatibility. Midpoints uses chromosomes present in supplied sites.

Baseline workloads do not add GC correction, blacklists or genomic scaling. Griffin's mandatory GC work remains with unit weights and its uncorrected output selected. Additional `cfdnalab-gc` and `griffin-gc` midpoint workloads estimate and apply native sample GC corrections. Default enabled cfDNAlab plots remain timed. Competitor source is unmodified.

Calculated metrics must be finite and respect their mathematical bounds. Correlations use [-1, 1]. Fractions, total variation, maximum CDF difference and normalized motif diversity use [0, 1]. Absolute boundary overshoots up to 1e-5 are clipped, and larger errors fail the comparison. Correlations are checked before clipping, including when calculated from NumPy covariance. Constant profiles retain `NA` correlation. Existing undefined metrics retain `NA`, and motif log ratios retain `Inf` when an abundant motif has zero frequency in either profile. Normalization totals must be finite before frequencies are calculated.

## Coverage

Sources: [runner](coverage/benchmark_coverage.sh), [comparator](coverage/compare_coverage.py), [operating notes](coverage/README.md).

### Methods

| Method | Configuration and native output |
| --- | --- |
| `cfdnalab-fragment` | `fcoverage` including inter-mate gaps, normal pairs/all contigs. `fcoverage.per_position.bedgraph.zst` |
| `cfdnalab-fragment-proper` | Same plus `--require-proper-pair` |
| `cfdnalab` | `fcoverage --ignore-gap`, paired read depth with overlap merged |
| `cfdnalab-autosomes` | Read-depth configuration restricted to actual BAM autosomes |
| `mosdepth-fragment` | `--fragment-mode --include-flag 1 --flag 3852`, requested MAPQ/lengths. `coverage.per-base.bed.gz` |
| `mosdepth` | Same native filters without fragment mode |
| `samtools-depth` | `depth -s -q 0 -Q MAPQ -G 3852 --require-flags 1` with HTSlib inclusive absolute TLEN bounds, uncompressed positions |
| `pandepth` | `-a -q MAPQ -x 3852 -t CORES`, compressed SiteDepth and summaries |
| `cfdnalab-regions` | `fcoverage --ignore-gap --by-size 500000 --per-window total`, requested MAPQ/lengths/all contigs. `fcoverage.total.tsv.zst` |
| `mosdepth-regions` | Ordinary mosdepth filters plus `--by 500000 --no-per-base`, `MOSDEPTH_PRECISION=5`. `coverage.regions.bed.gz` |
| `samtools-regions` | `bedcov -j -Q MAPQ -G 3852`, HTSlib paired/TLEN filter without threads. `samtools.bedcov.tsv` |
| `pandepth-regions` | `-w 500000` native fixed windows, `-q MAPQ -x 3852 -t CORES`, no `-a`. `pandepth.win.stat.gz` |
| `cfdnalab-region-counts` | `fcoverage --by-size 500000 --per-window total --normalize-by-length`, requested MAPQ/lengths, all contigs. `length_normalized.fcoverage.total.tsv.zst` |
| `cfdnalab-region-counts-proper` | Same normalized totals plus proper-pair requirement |
| `finaletoolkit-coverage` | `coverage` with native midpoint assignment and unit scale, requested MAPQ/lengths, no normalization. Compressed BED |

Choose methods explicitly through `BENCH_METHODS`. A complete baseline contains `cfdnalab-fragment`, `mosdepth-fragment`, `cfdnalab`, `mosdepth`, `samtools-depth`, `pandepth`, `cfdnalab-regions`, `mosdepth-regions`, `samtools-regions`, `pandepth-regions`, `cfdnalab-region-counts`, `finaletoolkit-coverage`, split into feasible jobs at each core setting. Proper-pair/autosome diagnostics are optional. Every available native thread/worker option receives the allocated core count directly, including the decoding options in mosdepth and samtools.

Regional methods use the same disjoint 500 kb coordinates from the BAM index, including shorter final windows. cfDNAlab generates them with `--by-size 500000`, while FinaleToolkit reads the shared BED4. cfDNAlab's length-normalized coverage gives each fragment total weight 1 before output rounding, split between windows when it crosses a boundary. FinaleToolkit assigns that fragment to the window containing its midpoint. Their regional counts are therefore approximate comparisons with different boundary handling.

The four `*-regions` methods measure summed read depth in 500 kb windows. cfDNAlab uses `--ignore-gap` without length normalization. mosdepth writes means with `MOSDEPTH_PRECISION=5` and no per-base output. samtools bedcov uses the shared BED4 and `-j`, with no thread option (`tool_threads=NA`). Its paired/TLEN expression is samtools' global HTSlib input filter. PanDepth uses `-w 500000` with native filters and double-counted mate overlap. Native outputs and compulsory summaries remain timed.

### Native differences

- cfDNAlab checks both mates and their directional span. Fragment coverage includes inter-mate gaps, excluding D/N positions unsupported by the other mate. Read depth additionally omits the inter-mate gap.
- mosdepth fragment mode requires proper pairs and fills complete TLEN intervals including D/N positions. Ordinary overlap correction depends on proper pairs, allowing some non-proper overlapping mates to count twice.
- samtools depth applies inclusive absolute TLEN bounds through HTSlib's input-expression filter. It does not reconstruct directional fragment spans or check that both mates pass MAPQ. PanDepth lacks fragment length filtering, paired-read inclusion and overlap removal. Its unsupported length filtering is recorded.
- FinaleToolkit regional coverage uses proper read1/TLEN selection. Its indexed read1 fetch can miss fragments whose spans alone overlap window boundaries. The end-motif length-bound bug does not establish the same bug in regional coverage.

See [cfDNAlab coverage](../src/commands/fcoverage/fcoverage.rs), [segmented fragment construction](../src/shared/fragment/segment_fragment.rs) and [pinned mosdepth source](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim). These native differences can produce unequal coverage on the same BAM.

### Loading, scope and memory

Supported retained formats are cfDNAlab zstd bedGraph, mosdepth gzip bedGraph and samtools/PanDepth `positions.tsv.gz`. Every processed chromosome is retained. samtools text is compressed outside timing, and PanDepth's compressed table is copied intact. samtools 1-based positions convert to 0-based half-open intervals. PanDepth positions are already 0-based. zstd reads all concatenated frames and propagates decompression errors.

Each per-base profile requires nonempty `chrom.sizes` with unique positive lengths. Runs must be finite/nonnegative, in bounds, sorted and non-overlapping within contiguous contig blocks. Different tools can order those blocks differently. Repeated blocks fail.

cfDNAlab and samtools depth may omit zero runs within their declared scope. Dense native formats must cover every base of every declared contig, including zeros. Missing/truncated dense tracks fail. Temporary per-contig run files merge adjacent equal values and reconcile contig order, avoiding per-base genomic arrays. `TMPDIR` selects comparison scratch.

Scopes must match unless `--common-scope` explicitly selects their contig intersection. Shared lengths must match and overlap must be nonempty. All per-base pair scopes are checked before large tracks are loaded. Results record evaluated `chromosomes` and `same_native_scope`. Fragment coverage and read depth are paired separately. Regional counts form a separate comparison set.

### Per-base metrics

The primary coverage similarity measure is **base-resolution MAE**: the sum of absolute depth differences at individual bases divided by the number of evaluated bases. Include bases where both profiles are zero. It is the average absolute depth difference per base, without binning or normalization. Zero means identical coverage throughout the evaluated scope. Weighting constant-run differences by their genomic lengths gives the same result as evaluating every base separately.

RMSE, signed mean, identical fraction, nonzero-union metrics, base-resolution Pearson and the residual distribution are secondary diagnostics. All per-base metrics use individual bases without genomic binning.

Segment weights are genomic lengths and signed differences are left minus right:

| Fields | Meaning |
| --- | --- |
| `mae` | **Primary:** base-resolution mean absolute depth difference |
| `bases`, `fraction_identical` | Evaluated bases and exact equal-depth fraction |
| `rmse` | Base-resolution root mean squared depth difference |
| `mean_signed_difference` | Base-weighted mean left minus right |
| `fraction_left_below`, `fraction_left_above` | Base fractions with negative/positive differences |
| `left_integrated_depth`, `right_integrated_depth` | Sum of depth times segment length |
| `nonzero_union_bases` | Bases with either depth nonzero |
| `union_fraction_identical`, `union_mae`, `union_rmse` | Nonzero-union metrics, `NA` for empty union |
| `pearson_r_1bp` | Pearson correlation of individual-base depths |

Pearson uses online covariance/variance weighted by genomic run lengths. This equals correlation over individual bases without expanding genomic arrays. Constant profiles have `NA` correlation, including equal zero profiles. Tile-split equal runs do not create errors.

### Residual distribution

Every comparison writes `residuals.tsv` beside the supplied `profiles/` directory, replacing that generated file only after all comparisons succeed. An error preserves the previous successful histogram. It contains one row per observed residual `left - right` for every supported per-base pair, sorted numerically within each pair. It includes zero residuals and preserves fractional values. Columns are `left_profile`, `right_profile`, `quantity`, `same_native_scope`, `chromosomes`, `residual`, `bases`, `fraction_of_bases`.

`bases` counts genomic bases with that residual, not runs. `fraction_of_bases` divides by all evaluated bases for that pair, including shared zeros. Counts sum to the summary's `bases`, and the distribution's base-weighted mean absolute residual equals the primary MAE. Scope and quantity match the corresponding summary row. Regional fragment counts are not included in this distribution.

Aggregation buffers at most 4096 distinct residual values before merging them into a temporary SQLite table. The database has a 2 MiB cache budget, resides under comparison scratch selected by `TMPDIR`, and is removed on exit. This bounds memory without binning or discarding residuals. The TSV can be used to plot a base-weighted residual histogram. Plot generation is separate from comparison.

### Regional metrics

The shared BED must have unique chromosome/start/end rows. cfDNAlab supplies the native `total_fragment_mass` column from its length-normalized `fcoverage` total table. FinaleToolkit supplies raw BED counts. Counts must be finite and nonnegative without duplicate regions. cfDNAlab values may be fractional and are not rounded by the comparator. FinaleToolkit values must remain integers. Each table must exactly match the shared coordinate set, independently of row order.

Results include region count, left/right summed regional counts, left/right total ratio, net difference, MAE/RMSE, identical fraction and Pearson r. FinaleToolkit's midpoint rule assigns each selected fragment to a single disjoint window, while cfDNAlab splits boundary mass fractionally. Zero right total yields `NA` ratio, including two zero totals. Constant-count correlation is `NA`. Quantity is `regional-fragment-counts`.

Regional read-depth inputs use cfDNAlab total_coverage, samtools bedcov's final column, PanDepth TotalDepth from its fixed-window table, and mosdepth means multiplied by actual widths. PanDepth uses 1-based inclusive coordinates. Tables match the shared BED except that PanDepth may omit final 1 bp windows. Compare reported shared windows and record left_missing_regions/right_missing_regions. Other missing or unexpected coordinates fail. Five-decimal mosdepth precision limits rounding to 2.5 covered bases per 500 kb window. Metrics remain covered-base totals, ratio, MAE/RMSE, identical fraction and Pearson r. Regional metrics do not replace base-resolution concordance.

## Fragment lengths

Sources: [runner](lengths/benchmark_lengths.sh), [R call](lengths/cfdnapro_lengths.R), [comparator](lengths/compare_lengths.py), [definition notes](lengths/README.md).

| Method | Configuration |
| --- | --- |
| `cfdnalab` | Global `lengths --length-bins MIN:MAX+1:1 --indel-mode ignore --clip-mode aligned`, all contigs |
| `cfdnalab-autosomes` | Same on actual BAM autosomes |
| `cfdnalab-proper` | All-contig histogram plus proper-pair requirement |
| `cfdnapro` | Explicit `readBam(curate_start_and_end=TRUE)`, requested MAPQ and paired/mapped/mate-mapped/nonsecondary/nonsupplementary/nonduplicate/QC-pass flags, then `callLength(frag_obj=...)` |
| `cfdnapro-outer-span` | Same with curation FALSE |
| `finaletoolkit` | `frag-length-bins --bin-size 1` with requested inclusive bounds/MAPQ |
| `samtools` | Full `stats`, `filter=mapq >= MAPQ`, `-f 1 -F 3852 -m 1 -i MAX+1` |
| `picard` | Timed standalone `picard -Xmx2g CollectInsertSizeMetrics I=sample.bam O=metrics.txt HISTOGRAM_WIDTH=MAX`, no MAPQ/QC-fail prefilter or PDF. Concordance reads the untimed filtered `concordance.metrics.txt` |

Choose methods explicitly through `BENCH_METHODS`. The lengths baseline contains cfDNAlab, cfDNAPro, FinaleToolkit, samtools and Picard with the named scope/definition diagnostics. samtools is labeled "samtools stats (full QC report, insert-size section extracted)". Plain Picard and R/container startup remain timed. After timing, the common-filter `samtools view -@ CORES -u -q MAPQ -f 1 -F 3852 | picard -Xmx2g ...` pipeline generates Picard's concordance histogram once. The timed native histogram is retained separately. Picard's PDF is omitted, while cfDNAlab's native QC plot stays timed.

Curated cfDNAPro uses the directional aligned span and uncurated mode uses the outer alignment span. Native pairing, orientation and reference-bound checks can differ. FinaleToolkit selects proper read1/TLEN. samtools supplies inward IS counts and halves passing alignment counts with integer truncation. Picard supplies its all-reads FR histogram.

For ordinary inward bwa mem alignments, stored TLEN matches the directional aligned span, scoped to the [inspected bwa 0.7.19 rule](https://github.com/lh3/bwa/blob/v0.7.19/bwamem.c#L887-L892). samtools/cfDNAlab inward selection uses forward/reverse start order. Picard's [5'-position <= comparison](https://github.com/samtools/htsjdk/blob/5.0.0/src/main/java/htsjdk/samtools/SamPairUtil.java#L101-L136) admits `forward.pos < reverse.reference_end`. cfDNAPro's [strict 1-based comparison](https://github.com/hw538/cfDNAPro/blob/ab9eeb166b8614b9ec2ce4ceda190dc295c65a7c/R/internalFunctions.R#L578-L605) is equivalent to `forward.pos < reverse.reference_end - 1`. Their one-base difference is outside the >=10 bp domain here. Arbitrary aligners/abnormal pairs are not assumed equivalent.

The comparator aligns every integer length in the requested inclusive range. cfDNAlab's 1 bp axis must match exactly. Other absent bins become zero. Duplicate/noninteger labels, non-finite/negative counts and unexpected out-of-range bins fail. samtools/Picard tails are discarded before normalization. A valid Picard histogram without FR supplies zero FR counts, while missing histogram sections fail.

Profiles report total counts, within-range frequencies, mean, shortest tied mode and fraction <=150 bp. Pairwise results include net/summed absolute count differences, count MAE/RMSE, identical-bin fraction, frequency Pearson r, frequency MAE/RMSE/max error, total variation, maximum CDF difference and Wasserstein distance in bp. On the 1 bp axis, Wasserstein sums absolute CDF differences between adjacent lengths.

Zero totals retain raw-count results with normalized metrics `NA` and status `zero-count-profile`. Constant frequencies have `NA` correlation. Curated/outer cfDNAPro is labeled `span-definition-diagnostic`, autosome/whole-reference pairs `different-chromosome-scope`, and others `cross-tool-or-filter`. Histograms cannot identify the fraction of individual fragments with changed lengths.

## End motifs

Sources: [runner](ends/benchmark_end_motifs.sh), [comparator](ends/compare_end_motifs.py), [native limitations](ends/README.md).

cfDNAlab uses `ends --k-inside K --k-outside 0 --all-motifs` with requested filters:

| Method | Source, clipping, pairs and scope |
| --- | --- |
| `cfdnalab-reference` | Reference, aligned boundaries, normal pairs, all contigs |
| `cfdnalab-reference-proper` | Same with proper pairs |
| `cfdnalab-reference-proper-wide` | Proper/reference, fixed effective 10-50,000 bp |
| `cfdnalab-reference-autosomes` | Reference/aligned, actual autosomes |
| `cfdnalab-read` | Read source, skip soft-clipped ends, all contigs |
| `cfdnalab-read-aligned` | Read source at aligned boundaries |
| `cfdnalab-read-include` | Read source, `include-at-aligned-boundary` |
| `finaletoolkit` | Reference `end-motifs --strand both`, requested K/bounds/MAPQ and allocated workers |
| `freia` | Native three-base extraction and default left-only aggregation |
| `cfdnapro` | Curated `readBam()` followed by reference-derived left-end `callMotif(motif_type="s")` |

The 4-mer comparison includes reference/native, reference/proper, FinaleToolkit and cfDNAPro. `BENCH_METHODS` must be supplied explicitly. K=3 is also supported and required for FrEIA. Both-end cfDNAlab/FinaleToolkit profiles and native left-end FrEIA/cfDNAPro profiles retain their end selection labels.

FinaleToolkit is pinned to commit `d16632ec5ae8999e774234457f6b0f9405c8f590`, with inclusive motif length bounds and unique midpoint window ownership fixed. Rows record effective bounds (`length_filter=applied`) and the source commit. The optional wide cfDNAlab diagnostic is no longer a workaround. Fragments longer than twice their read length can still be lost when read1 does not overlap the window containing their midpoint.

FinaleToolkit selection requires identical BAM/2bit contig dictionaries. Lengths divisible by 1,000,000 fail with the pinned last-window omission named. Read-only cfDNAlab and FrEIA-only selections need no 2bit. FinaleToolkit's read1 window queries can still double-count boundaries, and proper read1/TLEN selection differs from cfDNAlab pair/CIGAR checks. See the pinned [end-motif implementation](https://github.com/epifluidlab/FinaleToolkit/blob/d16632ec5ae8999e774234457f6b0f9405c8f590/src/finaletoolkit/frag/_end_motifs.py).

FrEIA timing includes the authors' streaming `samtools view | tee | samtools index -` filtering/indexing pipeline, extraction `-n 3`, unchanged default left-only aggregation, exclusive maximum `MAX+1`, and no intentional downsampling/bootstrap. Native coordinate deduplication, filtering N at either end, intermediates and mononucleotide work remain. The samtools filtering and indexing commands each receive N directly through `-@`. Extraction receives N workers, with the native size safeguard still forcing one worker above 8000 MiB filtered BAM size. The actual extraction-worker count is recorded.

FrEIA's broken `--both` aggregate remains excluded. The runner uses intact left-only output and retains it as `left-motifs.pq`. The comparator verifies genome-level `WhichEnd=leftmost`, then sums additive `CountLen` over fragment lengths for each motif. `CountSamp` repeats a motif total across lengths and is not summed. Missing ACGT motifs receive zero. Old `T__sample.pq` files remain skipped without recovering or reconstructing their damaged counts.

cfDNAPro uses the same explicit `readBam` MAPQ/flags as its curated length comparison, all contigs and inclusive width filtering. Native `callMotif(motif_type="s", motif_length=K, integrate_mut=FALSE, downsample_ref=FALSE)` counts the leftmost K reference bases on the forward strand. It writes `motif`, `n` and `fraction`. BAM reading, curation, reference loading, motif extraction and writing stay timed. There is no native thread option for this call. The image's selected BSgenome label/package/version/sequence metadata are retained in `reference.tsv`. Reference identity with the alignment FASTA/2bit must be verified before attributing differences solely to tools. Native pair/orientation/reference-bound checks and cfDNAlab end/CIGAR checks remain different.

Profiles must contain exactly all `4^K` inward ACGT labels after filling FrEIA's sparse zeros, without duplicate/non-finite/negative values. cfDNAlab's leading underscore is removed, raw global counts are loaded natively, cfDNAPro native counts are read directly and FinaleToolkit TSV frequencies read directly. Total mass must be positive. FinaleToolkit sums must be close to 1 with `rtol=1e-6`, `atol=1e-9`. Profiles normalize by their native counted-end totals, which are `NA` for FinaleToolkit.

Pairwise rows report `left_end_selection`, `right_end_selection` and `side_comparison`. A left-end versus both-end comparison is labeled `approximate-side-mirroring`. This assumes the two sides have sufficiently similar inward-oriented motif frequencies for the intended approximate comparison, and does not establish that their observed frequencies are identical. Comparisons selecting the same sides are labeled `same-sides`. Native left counts remain native left counts. No doubled count estimate is produced, and these approximations do not support an equivalent-output timing claim.

All supported profiles are compared pairwise. Metrics are Pearson r, Spearman r using average ranks, frequency MAE/RMSE/max error, total variation, normalized Shannon diversity and its absolute difference. Diversity is `-sum(p * log(p)) / log(4^K)`, omitting zero terms. Constant correlations are `NA`.

`max_abs_log2_ratio_abundant` uses motifs reaching >=1e-4 in either profile. A qualifying zero on the other side gives `Inf`, without a pseudocount. No qualifying motif gives `NA`. Configuration labels and metadata retain scope/filter/source distinctions. cfDNAFE and original FRAGMA are excluded from this runnable set.

## Midpoint profiles

Sources: [runner](midpoints/benchmark_midpoints.sh), [comparator](midpoints/compare_midpoints.py), [workload notes](midpoints/README.md).

Select `cfdnalab`, `cfdnalab-gc` and `griffin-gc` in separate jobs. Only the latter two form a timed comparison. Griffin's full workflow also writes the raw table used for uncorrected concordance. All midpoint jobs use MAPQ >=20, inclusive 100-200 bp fragments, supplied 9,990 bp windows and 15 bp bins. No extra site selection or proper-pair requirement is added. Site derivation is documented in [the command report](experiment_commands.md#site-preparation).

Griffin gets one short untimed warmup using its release's 1000-site CTCF demo, the same production BAM/reference, length/MAPQ/window settings and core allocation. Its separate Slurm step and output are excluded from measured timings, memory and concordance. There are no full-input Griffin warmups. Production site coordinates are unchanged. This is a practical warmup intended to give Griffin a fair start without another full production run.

cfDNAlab uses `midpoints --length-bins MIN MAX+1 --bin-size 15`, default smoothing and the shared exclusion BEDs, actual-pair checks and deterministic choice between the two central bases of even fragments. Griffin retains forward-record/TLEN selection and the right central base. Fragment/filter populations can differ.

The sole timed Griffin workload is `griffin-gc`: GC estimation, both profiling scripts, raw/corrected bigWigs, aggregation, normalization and compulsory features. Both tables are retained from its first measured repetition. The official --save_window -1000 1000 rounds to [-990, 990), or 132 bins. Normalization remains [-4995, 4995), or 666 bins. Native exclusion BEDs, zero-mappability masking, outlier removal and smoothing stay on. Mappability correction, per-site CNA correction and individual output stay off. Griffin's final plotting rule is omitted.

Corrected `cfdnalab-gc` timing includes `gc-bias`, then `midpoints --gc-file ... --ref-2bit ...`. Corrected `griffin-gc` timing includes native GC counting and bias fitting with the release's 15-500 bp estimation range, followed by both profiling scripts with the agreed 100-200 bp range, selecting `sample.GC_corrected.coverage.tsv`. Both rebuild sample bias on every repetition, retaining native bias fitting/smoothing/exclusions and plots. cfDNAlab's reference is calculated once with native defaults and reused across samples and benchmark runs. Sample estimation keeps its native windows and chromosome selection, inheriting the default reference's 30-1000 bp range. Griffin uses its bundled hg38 region BED and precomputed 15-500 bp GC-frequency tables, already downloaded under `Griffin-source/Ref/`. Setup records the Griffin estimation bounds in `run.tsv`. These inputs are separate from the supplied profiling sites and from cfDNAlab's reference package. Its per-length files use that BED's original basename. Both tools' midpoint estimation and profiling use MAPQ >=20, matching Griffin's workflow. References must match the alignment assembly. Their native GC coordinates, fragment selection, end offsets, estimation ranges and fitting methods differ, so corrected-profile agreement is not an isolated test of midpoint counting.

Griffin uses its short untimed unit-weight demo without sample-bias estimation for warmup. Corrected cfDNAlab warmups include estimation. Corrected profiles and first measured estimates are saved in their *-gc directories. --gc-corrected pairs the corrected tables. The default comparator pairs cfdnalab with griffin-gc/sample.uncorrected.coverage.tsv. Mean densities represent weighted mass for corrected profiles.

Griffin fetch padding uses minimum length. `floor(maximum / 2) <= minimum` is sufficient, satisfied by 100-200 bp. Fetch losses become possible from `maximum >= 2 * minimum + 2`, depending on alignment extent. Aggregation parallelism depends on supplied group count.

The comparator requires one cfDNAlab length bin and contiguous native bins. It validates bin alignment and exactly the Griffin saved-position labels. Divide cfDNAlab's complete 666 bins by their full-window mean before selecting the 132 starts in [-990, 990). Do not renormalize the saved subset. Finite smoothing undershoots remain. Undefined Griffin positions make metrics unavailable rather than changing the compared subset.

cfDNAlab group sums are averaged across bases within bins. Griffin sums those bases, but the constant factor cancels during normalization. cfDNAlab smooths before binning; Griffin bins first and smooths over 11 bins. The comparator adds neither smoothing nor rebinning. Native full-window per-site/per-bp density metadata stays separate.

Griffin rounds profiles/mean metadata to five decimals. The saved central profile need not have mean one because normalization precedes cropping. Sparse positive profiles can have density metadata rounded to zero. All-NaN profiles with zero density indicate zero coverage. Undefined density or positive density with undefined saved positions receives an unavailable-profile status.

Rows report group, cfDNAlab eligible intervals (cfdnalab_sites), Griffin listed sites including masked sites (griffin_sites), compared bins (positions=132), compared width (window_bp=1980), normalization_window_bp=9990, resolution_bp=15, native densities, r, r squared, errors and status. Constant profiles have unavailable correlation. Empty and undefined-profile statuses stay distinct.

Positional metrics also report nonzero_intersection_bases and coverage_jaccard (intersection / union, or NA if both are empty). coverage_jaccard.tsv records these per chromosome in BAM order. A partial Jaccard followed by zero samtools coverage on every later informative chromosome flags possible truncation without failing on a similarity threshold. Missing samtools positions are zero.

cfDNAPro motifs use timed user-side R length subsetting and retain BAM contigs present in the selected BSgenome, recorded in contig_scope.tsv. For this BAM the excluded contig is chrEBV. Lengths keep whole-BAM scope. FrEIA runs its three per-sample extraction/abundance rules directly without Snakemake and stops before group comparison/scoring.
