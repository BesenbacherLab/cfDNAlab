# Benchmark methods and concordance specification

This specifies implemented native configurations and comparisons. Installation, execution, accounting and inputs are specified in [benchmarks_spec.md](benchmarks_spec.md). Remaining implementation and verification are in [benchmarks_review_plan.md](benchmarks_review_plan.md).

Fragment coverage, read depth, regional fragment counts, global length histograms, end-motif frequencies and midpoint profiles are separate quantities. Normalized agreement does not establish equal total counts or individual fragment selection. Speed claims apply to the configurations actually timed, including their native compulsory work and filter differences.

## Shared scientific settings

cfDNAlab paired spans remain directional: `forward.pos` to `reverse.reference_end`. Competitors retain their native pairing, TLEN and CIGAR rules. Proper-pair diagnostics add flag 0x2 selection without making all inclusion rules identical.

Normal settings request MAPQ >=30 and inclusive fragment lengths 30-1000 bp. Midpoints uses 100-200 bp. cfDNAlab length bins end at `maximum + 1`. Global scope is all BAM contigs unless labeled autosomes. Midpoints uses chromosomes present in supplied sites.

GC correction, blacklists and genomic scaling are not added. Griffin's mandatory GC work remains with unit weights and its uncorrected output selected. Default enabled cfDNAlab plots remain timed. Competitor source is unmodified.

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
| `samtools-depth` | `depth -aa -s -q 0 -Q MAPQ -G 3852 --require-flags 1`, uncompressed positions |
| `pandepth` | `-a -q MAPQ -x 3852 -t CORES`, compressed SiteDepth and summaries |
| `cfdnalab-region-counts` | `lengths --by-bed --assign-by any --length-bins MIN MAX+1 --indel-mode ignore --clip-mode aligned`, all contigs. `length_counts.tsv.zst` |
| `cfdnalab-region-counts-proper` | Same regional counts plus proper-pair requirement |
| `finaletoolkit-coverage` | `coverage --intersect-policy any --scale-factor 1`, requested MAPQ/lengths, no normalization. Compressed BED |

Choose methods explicitly through `BENCH_METHODS`. A complete baseline contains `cfdnalab-fragment`, `mosdepth-fragment`, `cfdnalab`, `mosdepth`, `samtools-depth`, `pandepth`, `cfdnalab-region-counts`, `finaletoolkit-coverage`, split into feasible jobs at each core setting. Proper-pair/autosome diagnostics are optional. cfDNAlab/PanDepth receive the allocated cores, mosdepth/samtools get `cores - 1` additional decompression threads, and FinaleToolkit gets the allocated worker count.

Regional methods share disjoint 500 kb BED4 windows generated from the BAM index, including shorter final windows. Each fragment counts in every overlapping window, so summed regional counts are not unique-fragment totals.

### Native differences

- cfDNAlab checks both mates and their directional span. Fragment coverage includes inter-mate gaps, excluding D/N positions unsupported by the other mate. Read depth additionally omits the inter-mate gap.
- mosdepth fragment mode requires proper pairs and fills complete TLEN intervals including D/N positions. Ordinary overlap correction depends on proper pairs, allowing some non-proper overlapping mates to count twice.
- samtools depth lacks a matching length filter or complete-pair/MAPQ check. PanDepth additionally lacks paired-read inclusion and overlap removal. Unsupported length filtering is recorded.
- FinaleToolkit regional coverage uses proper read1/TLEN selection. Its indexed read1 fetch can miss fragments whose spans alone overlap window boundaries. The end-motif length-bound bug does not establish the same bug in regional coverage.

See [cfDNAlab coverage](../src/commands/fcoverage/fcoverage.rs), [segmented fragment construction](../src/shared/fragment/segment_fragment.rs) and [pinned mosdepth source](https://github.com/brentp/mosdepth/blob/v0.3.14/mosdepth.nim). General real-BAM equality is not promised.

### Loading, scope and memory

Supported retained formats are cfDNAlab zstd bedGraph, mosdepth gzip bedGraph and samtools/PanDepth `positions.tsv.gz`. samtools 1-based positions convert to 0-based half-open intervals. PanDepth positions are already 0-based. zstd reads all concatenated frames and propagates decompression errors.

Each per-base profile requires nonempty `chrom.sizes` with unique positive lengths. Runs must be finite/nonnegative, in bounds, sorted and non-overlapping within contiguous contig blocks. Different tools can order those blocks differently. Repeated blocks fail.

Only cfDNAlab's sparse format may omit zero runs, within its declared scope. Dense native formats must cover every base of every declared contig, including zeros. Missing/truncated dense tracks fail. Temporary per-contig run files merge adjacent equal values and reconcile contig order, avoiding per-base genomic arrays. `TMPDIR` selects comparison scratch.

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

The shared BED must have unique chromosome/start/end rows. cfDNAlab input is a windowed count table with a single length bin. FinaleToolkit supplies raw BED counts. Counts must be finite nonnegative integers without duplicate regions. Each table must exactly match the shared coordinate set, independently of row order.

Results include region count, left/right summed regional counts, left/right total ratio, net difference, MAE/RMSE, identical fraction and Pearson r. Zero right total yields `NA` ratio, including two zero totals. Constant-count correlation is `NA`. Quantity is `regional-fragment-counts`.

## Fragment lengths

Sources: [runner](lengths/benchmark_lengths.sh), [R call](lengths/cfdnapro_lengths.R), [comparator](lengths/compare_lengths.py), [definition notes](lengths/README.md).

| Method | Configuration |
| --- | --- |
| `cfdnalab` | Global `lengths --length-bins MIN:MAX+1:1 --indel-mode ignore --clip-mode aligned`, all contigs |
| `cfdnalab-autosomes` | Same on actual BAM autosomes |
| `cfdnalab-proper` | All-contig histogram plus proper-pair requirement |
| `cfdnapro` | Explicit `readBam(curate_start_and_end=TRUE)`, requested MAPQ and paired/mapped/mate-mapped/nonsecondary/nonsupplementary/nonduplicate/QC-pass flags, then `read_bam_insert_metrics(fragment_obj=...)` |
| `cfdnapro-outer-span` | Same with curation FALSE |
| `finaletoolkit` | `frag-length-bins --bin-size 1` with requested inclusive bounds/MAPQ |
| `samtools` | Full `stats`, `filter=mapq >= MAPQ`, `-f 1 -F 3852 -m 1 -i MAX+1` |
| `picard` | Timed `samtools view -u -q MAPQ -f 1 -F 3852` pipeline into `CollectInsertSizeMetrics HISTOGRAM_WIDTH=MAX MINIMUM_PCT=0 INCLUDE_DUPLICATES=false`, no PDF |

Choose methods explicitly through `BENCH_METHODS`. The baseline contains `cfdnalab`, `cfdnalab-autosomes`, `cfdnapro`, `finaletoolkit`, `samtools`, `picard`, split into feasible jobs at each core setting. Proper-pair/outer-span diagnostics are optional. Full R/container startup and Picard's filter/JVM pipeline remain timed. Picard uses a 4 GiB heap and allocated `ActiveProcessorCount`. samtools gets `cores - 1` additional threads.

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
| `freia` | Optional native three-base extraction/aggregation timing only |

The baseline selection is reference/native, reference/proper and FinaleToolkit with K=4. `BENCH_METHODS` must be supplied explicitly. K=3 is supported. FrEIA requires K=3 and explicitly warns that its pinned aggregate is incorrect.

FinaleToolkit 1.1.0 ignores requested end-motif bounds. Rows retain requested bounds, set effective bounds `NA` and label `length_filter=ignored-by-release`. The wide diagnostic practically removes cfDNAlab's usual 30-1000 bp restriction for this dataset. No >50,000 bp census is added.

FinaleToolkit selection requires identical BAM/2bit contig dictionaries. Lengths divisible by 1,000,000 fail with the pinned last-window omission named. Read-only cfDNAlab and FrEIA-only selections need no 2bit. FinaleToolkit's read1 window queries can still double-count boundaries, and proper read1/TLEN selection differs from cfDNAlab pair/CIGAR checks. See the pinned [end-motif implementation](https://github.com/epifluidlab/FinaleToolkit/blob/v1.1.0/src/finaletoolkit/frag/_end_motifs.py).

FrEIA timing includes BAM filtering/indexing, extraction `-n 3`, unchanged `--both` aggregation, exclusive maximum `MAX+1`, and no intentional downsampling/bootstrap. Native coordinate deduplication, N filtering, intermediates and mononucleotide work remain. Extraction workers are recorded as 1 above 8000 MiB filtered BAM size.

**FrEIA concordance is disabled until native aggregation is fixed.** Its pinned assignment can drop/pad right-end rows. No CountSamp/CountLen recovery or replacement aggregation is used. Its profile is not retained, old profile directories are skipped, and optional timings cannot support equivalent-output claims.

Profiles must contain exactly all `4^K` inward ACGT labels including zeros, without duplicates/non-finite/negative values. cfDNAlab's leading underscore is removed, raw global counts are loaded natively and FinaleToolkit TSV frequencies read directly. Total mass must be positive. FinaleToolkit sums must be close to 1 with `rtol=1e-6`, `atol=1e-9`. Profiles normalize by counted-end totals, which are `NA` for FinaleToolkit.

All supported profiles are compared pairwise. Metrics are Pearson r, Spearman r using average ranks, frequency MAE/RMSE/max error, total variation, normalized Shannon diversity and its absolute difference. Diversity is `-sum(p * log(p)) / log(4^K)`, omitting zero terms. Constant correlations are `NA`.

`max_abs_log2_ratio_abundant` uses motifs reaching >=1e-4 in either profile. A qualifying zero on the other side gives `Inf`, without a pseudocount. No qualifying motif gives `NA`. Configuration labels and metadata retain scope/filter/source distinctions. cfDNAFE and original FRAGMA are excluded from this runnable set.

## Midpoint profiles

Sources: [runner](midpoints/benchmark_midpoints.sh), [comparator](midpoints/compare_midpoints.py), [workload notes](midpoints/README.md).

The complete set is selected explicitly with `BENCH_METHODS="cfdnalab griffin"`. Both use the user's supplied windows/chromosomes, inclusive 100-200 bp fragments, MAPQ >=30 and unsmoothed 1 bp profiles. Width comes from the BED. Proper-pair filtering is not added.

cfDNAlab uses `midpoints --length-bins MIN MAX+1 --smoothing none --bin-size 1`, actual-pair checks and deterministic choice between the two central bases of even fragments. Griffin retains forward-record/TLEN selection and the right central base. Fragment/filter populations can differ.

Griffin timing includes both native scripts, reference/GC work, raw/GC-weighted bigWigs, group aggregation, final mean-one normalization and features. Setup supplies unit GC bias over requested lengths and all possible GC counts. Comparison selects `sample.uncorrected.coverage.tsv`. Masking, mappability correction, outlier exclusion, site-wise CNA normalization, smoothing and individual-site output are disabled.

Griffin fetch padding uses minimum length. `floor(maximum / 2) <= minimum` is sufficient, satisfied by 100-200 bp. Fetch losses become possible from `maximum >= 2 * minimum + 2`, depending on alignment extent. Aggregation parallelism depends on supplied group count.

The comparator requires one cfDNAlab length bin, nonzero even width and contiguous 1 bp positions. Group identities, positive site counts and Griffin signed position labels must match. Unexpected smoothing/CNA settings and invalid profiles fail.

cfDNAlab summed counts divide by their complete-window mean to produce a mean-one profile. No site-wise normalization, shift, smoothing, rebinning or zero removal occurs. Raw mean per-site/per-bp densities remain separate.

Griffin rounds profiles/mean metadata to five decimals. Finite profiles need mean one within absolute tolerance 1e-5. Sparse positive profiles may have mean metadata rounded to zero, so that alone does not mark emptiness. Empty Griffin profiles contain all NaN values and zero mean metadata.

Rows report group, sites, positions, both raw densities, Pearson r, normalized MAE/RMSE/max difference and status. Empty statuses distinguish `both-zero-coverage`, `cfdnalab-zero-coverage`, `griffin-zero-coverage`, with normalized metrics `NA`. Constant positive profiles have undefined correlation even when equal. Density differences remain relevant because normalization can hide inclusion changes.
