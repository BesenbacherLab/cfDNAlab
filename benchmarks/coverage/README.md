# Coverage concordance and timing

Use the shared [installer](../README.md#installation) and [prepared hg38 input](../data/README.md). Run fragment coverage, read depth and regional fragment counts as separate comparisons. The baseline requests MAPQ >=30, lengths 30-1000 bp, all BAM contigs and no GC correction, blacklist or genomic scaling. Native differences remain visible rather than being patched out.

| Method | Quantity and native output | Remaining differences |
| --- | --- | --- |
| `cfdnalab-fragment` | Paired fragment coverage, sparse compressed bedGraph | Directional `forward.pos` to `reverse.reference_end`. Both mates pass filters. Includes inter-mate gaps but excludes D/N positions unsupported by the other mate. |
| `cfdnalab-fragment-proper` | Same, requiring flag 0x2 | Closer filter alignment to mosdepth fragment mode. Both-mate filtering, orientation, TLEN and CIGAR differences remain. |
| `mosdepth-fragment` | `--fragment-mode`, compressed bedGraph | Requires proper pairs and fills the complete TLEN interval, including D/N positions. Native read/TLEN selection differs from actual paired spans. |
| `cfdnalab`, `cfdnalab-autosomes` | Paired read depth with `--ignore-gap`, compressed bedGraph | Merges mate overlap. The autosome configuration has a smaller scope, explicitly recorded. |
| `mosdepth` | Ordinary overlap-corrected read depth, compressed bedGraph | Overlap correction depends on proper pairs, so non-proper overlapping mates can count twice. Native CIGAR and record filtering remain. |
| `samtools-depth` | `depth -aa -s -q 0 -Q 30 -G 3852 --require-flags 1` with HTSlib TLEN bounds, uncompressed positions | Inclusive absolute TLEN filtering, without reconstructing fragment spans or checking that both mates pass MAPQ. Its native position is 1-based. |
| `pandepth` | `-a -q 30 -x 3852`, compressed SiteDepth and summaries | No fragment length filter, paired-read inclusion requirement or mate-overlap removal. Its native position is 0-based. |
| `cfdnalab-regions` | `fcoverage --ignore-gap --by-size 500000 --per-window total`, compressed total coverage | Sums read depth with mate overlap merged, without length normalization. |
| `mosdepth-regions` | Ordinary mosdepth with `--by 500000 --no-per-base`, compressed regional means | Same native selection/overlap behavior as per-base mosdepth. Means are printed to its default precision. |
| `samtools-regions` | `bedcov -j -Q 30 -G 3852`, HTSlib paired/TLEN filters and decoding threads | Sums read depth without mate-overlap removal. No per-base output. |
| `pandepth-regions` | `-b` with the shared windows converted to native inclusive coordinates, MAPQ/flags/threads as above | Uses native `TotalDepth`, without per-base output. No fragment length filter or mate-overlap removal. |
| `cfdnalab-region-counts` | `fcoverage --by-size 500000 --per-window total --normalize-by-length` | Estimates fragment counts by summing unit-weight fragment mass. A boundary fragment's contribution is split between windows. |
| `cfdnalab-region-counts-proper` | Same, requiring flag 0x2 | Closer filter alignment to FinaleToolkit, retaining actual paired-span and both-mate checks. |
| `finaletoolkit-coverage` | Raw counts, `coverage --intersect-policy any --scale-factor 1`, compressed BED | Proper read1/TLEN selection. Indexed fetch requires read1 to overlap the window, which can miss fragments whose span alone overlaps a boundary. |

The common exclusion mask removes unmapped, mate-unmapped, secondary, supplementary, QC-fail and duplicate records. A flag mask alone does not give tools identical pair selection. The end-motif release's ignored length bounds are a separate bug and do not imply that FinaleToolkit's regional coverage ignores them.

samtools depth applies the requested inclusive length bounds to both positive and negative TLEN using `--input-fmt-option 'filter=(tlen >= MIN && tlen <= MAX) || (tlen <= -MIN && tlen >= -MAX)'`. This uses the [HTSlib input-expression filter](https://www.htslib.org/doc/samtools.html#FILTER_EXPRESSIONS) during the native BAM read pass, without a separate filtering command.

For example, forward `pos=100`, `30M`, reverse `pos=170`, `30M` spans 100 bp. Fragment coverage fills the 40 bp inter-mate gap, whereas read depth covers 60 bp. With forward `15M2D15M`, cfDNAlab fragment coverage leaves `[115,117)` uncovered unless supported by the other mate. mosdepth fragment mode fills all 100 bp. Native CIGAR and pairing differences can therefore produce unequal coverage on the same BAM.

## Run and compare

Set the input/tool paths printed by the installer and use a new results directory. Every run also needs `BENCH_SAMTOOLS_ARCHIVE` for untimed header preparation and Python 3 for environment relocation. Jobs selecting multiple methods require `BENCH_COMPARE_PYTHON`, also printed by the installer.

```bash
export BENCH_BAM=/path/to/BH01.hg38.paired.15x.bam
# The batch job uses TMPDIR supplied on the compute node
export BENCH_RESULTS=/path/to/coverage-concordance
export BENCH_CORES=1
export BENCH_PILOT=1
export BENCH_METHODS="cfdnalab-fragment cfdnalab-fragment-proper mosdepth-fragment cfdnalab mosdepth samtools-depth pandepth cfdnalab-regions mosdepth-regions samtools-regions pandepth-regions cfdnalab-region-counts cfdnalab-region-counts-proper finaletoolkit-coverage"
bash benchmarks/submit_benchmark.sh coverage --time=YOUR_PILOT_BASED_BUDGET

# The job writes concordance.tsv and residuals.tsv automatically after timing
```

Pilots and full runs calculate concordance after all tool attempts finish. Results are saved as `concordance.tsv`, including the primary base-resolution `mae`, and `residuals.tsv`, the residual histogram. Only profiles measuring the same quantity are paired. Different retained chromosome scopes are compared on their intersection, with the evaluated scope recorded in each row. A job without a comparable pair prints why concordance was skipped. A comparison failure fails the job and keeps its error in `concordance-error.log`, while preserving the completed timings and profiles.

Full coverage tracks are preserved after the first measured step at the job's single core count, with every processed chromosome retained. samtools and PanDepth generate their complete native outputs **inside timing**. PanDepth's compressed positions are copied intact, and samtools text is compressed after the Slurm tool step. Each profile has an explicit `chrom.sizes`, so all-contig comparisons evaluate the full BAM reference scope.

The comparator requires equal declared scope unless `--common-scope` explicitly allows intersection. It checks every per-base pair's scope before loading tracks. Every result records evaluated chromosomes and whether the original scopes match. A missing chromosome in a sparse cfDNAlab track means zero only when declared in its processed scope. Other formats must contain complete runs/positions for every declared contig, so truncated files fail. Runs must be sorted within each contiguous chromosome block, non-overlapping, in bounds and finite/nonnegative. Temporary coalesced run files reconcile arbitrary contig order without allocating per-base arrays. Comparison scratch uses `TMPDIR`, which the job sets below its own working directory. zstd reads all concatenated frames.

The primary similarity measure is **base-resolution MAE**: `sum(abs(left_depth - right_depth)) / evaluated_bases`, summing over individual bases in the evaluated scope, including shared zeros. This is the average absolute depth difference per base, without binning or normalization. Constant runs are weighted by their genomic lengths, so bedGraph segmentation does not change the result.

RMSE, identical fraction, mean signed difference and fractions with cfDNAlab/left values above or below the other profile are secondary diagnostics. Nonzero-union metrics exclude bases where both are zero. The secondary `pearson_r_1bp` correlates individual-base depths using genomic run lengths as weights, without binning. Constant-profile correlations and empty-union metrics are `NA`. Read-depth and fragment-coverage tracks are not pooled.

Every comparison saves the full residual distribution to `residuals.tsv` beside `profiles/`, replacing that generated file after successful completion. Failures preserve any previous histogram. It contains profile pair, quantity, scope, residual `left - right`, `bases` and `fraction_of_bases` for all supported per-base pairs. Rows are sorted by residual within each pair and include zero residuals. Counts are weighted by genomic bases and fractions divide by all evaluated bases, matching the primary MAE's scope. Bounded batches merge into temporary SQLite storage under `TMPDIR`. This provides the values for a residual histogram without expanding tracks into arrays. Regional counts remain separate and produce no residual-distribution rows.

Regional comparison requires the exact unique chromosome/start/end rows in the shared `regions.bed`. cfDNAlab uses `--by-size 500000` to generate those coordinates, including shorter final windows, and writes `length_normalized.fcoverage.total.tsv.zst`. Each fragment contributes total weight 1 before output rounding, split across windows when it crosses a boundary. FinaleToolkit counts a boundary fragment once in every overlapping window, so its summed regional counts can exceed its unique-fragment total. The comparator preserves cfDNAlab's fractional values and reports totals, total ratio, net difference, MAE/RMSE, identical fraction and correlation. These are approximate count comparisons with different boundary handling.

The four `*-regions` methods form a separate **regional read-depth** comparison on the same 500 kb windows. Values are summed per-base read depths, also called covered-base totals. cfDNAlab uses unnormalized `fcoverage --ignore-gap`, [mosdepth](https://github.com/brentp/mosdepth/blob/v0.3.14/README.md) uses native regional means, [samtools](https://www.htslib.org/doc/samtools-bedcov.html) uses `bedcov -j`, and [PanDepth](https://github.com/HuiyangYu/PanDepth/blob/main/src/PanDepth.cpp) uses native BED statistics. PanDepth 2.26 treats BED starts/ends as 1-based inclusive, so setup writes a scratch-only three-column copy with start increased by 1. The comparator converts its output back to BED coordinates and reads `TotalDepth`. mosdepth means are multiplied by each actual window width, including shorter final windows, to recover covered-base totals at the precision of its native output. Its default two-decimal means introduce rounding differences. Comparisons retain the tools' native mate-overlap and filtering differences. They are labeled `regional-read-depth` and never paired with normalized fragment counts or positional tracks.

For a regional-only job, select `BENCH_METHODS="cfdnalab-regions mosdepth-regions samtools-regions pandepth-regions"`. Submit separate jobs at 1, 4, 12 and 32 cores with the usual warmup/repetition settings. Each native command writes regional output directly, without generating a full per-base table. Native summaries/indexes remain timed. samtools bedcov receives `cores - 1` BAM decoding threads through `--input-fmt-option nthreads=...`, which [HTSlib applies to its BGZF reader](https://github.com/samtools/htslib/blob/1.25/hts.c).

## Production measurements

Use `BENCH_PILOT=1` on the actual production method/core matrix to estimate all work, splitting pilots when necessary. The submission launcher allows shared nodes for pilots and requests 50 GB unless `--mem` is passed. Then set `BENCH_PILOT=0`, a warmup and three measurements, and submit core/method blocks that fit the allocation limit on matched hardware. Full runs reserve an exclusive node and all its memory. Prefer the same physical node when available. See [shared run controls and Slurm memory](../README.md#runs-and-budgets). `BENCH_MIN_MAPQ`, `BENCH_MIN_LENGTH` and `BENCH_MAX_LENGTH` override the filters. The timing TSV records applied/unsupported length filtering, proper-pair selection and scope.

Keep process-plus-filesystem-flush time, output bytes and the joined Slurm peak alongside process time and GNU time RSS. Native positional output costs stay in the comparison. Include sample size, SD and median/quartiles and retain all measured points. The table below is a result template, not measured concordance data.

| Comparison | Scope | MAE, 1 bp (primary) | RMSE, 1 bp | Identical bases | Nonzero-union MAE | Pearson r, 1 bp |
| --- | --- | --- | --- | --- | --- | --- |
| Fragment/proper vs mosdepth fragment | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] |
| Native fragment vs mosdepth fragment | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] |
| Paired read depth vs mosdepth | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] |
| Paired read depth vs samtools / PanDepth | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] |

Report regional counts in a separate table and label native/proper diagnostics separately. Only configurations actually timed can support a corresponding speed comparison.
