# Coverage concordance and timing

Use the shared [installer](../README.md#installation) and [prepared hg38 input](../data/README.md). Run fragment coverage, read depth and regional fragment counts as separate comparisons. The baseline requests MAPQ >=30, lengths 30-1000 bp, all BAM contigs and no GC correction, blacklist or genomic scaling. Native differences remain visible rather than being patched out.

| Method | Quantity and native output | Remaining differences |
| --- | --- | --- |
| `cfdnalab-fragment` | Paired fragment coverage, sparse compressed bedGraph | Directional `forward.pos` to `reverse.reference_end`. Both mates pass filters. Includes inter-mate gaps but excludes D/N positions unsupported by the other mate. |
| `cfdnalab-fragment-proper` | Same, requiring flag 0x2 | Closer filter alignment to mosdepth fragment mode. Both-mate filtering, orientation, TLEN and CIGAR differences remain. |
| `mosdepth-fragment` | `--fragment-mode`, compressed bedGraph | Requires proper pairs and fills the complete TLEN interval, including D/N positions. Native read/TLEN selection differs from actual paired spans. |
| `cfdnalab`, `cfdnalab-autosomes` | Paired read depth with `--ignore-gap`, compressed bedGraph | Merges mate overlap. The autosome configuration has a smaller scope, explicitly recorded. |
| `mosdepth` | Ordinary overlap-corrected read depth, compressed bedGraph | Overlap correction depends on proper pairs, so non-proper overlapping mates can count twice. Native CIGAR and record filtering remain. |
| `samtools-depth` | `depth -aa -s -q 0 -Q 30 -G 3852 --require-flags 1`, uncompressed positions | No fragment length filter or complete-pair/MAPQ check. Its native position is 1-based. |
| `pandepth` | `-a -q 30 -x 3852`, compressed SiteDepth and summaries | No fragment length filter, paired-read inclusion requirement or mate-overlap removal. Its native position is 0-based. |
| `cfdnalab-region-counts` | `lengths --by-bed --assign-by any --length-bins 30 1001 --indel-mode ignore --clip-mode aligned` | Counts each included fragment in every overlapping shared 500 kb window. It is not mean base depth. |
| `cfdnalab-region-counts-proper` | Same, requiring flag 0x2 | Closer filter alignment to FinaleToolkit, retaining actual paired-span and both-mate checks. |
| `finaletoolkit-coverage` | Raw counts, `coverage --intersect-policy any --scale-factor 1`, compressed BED | Proper read1/TLEN selection. Indexed fetch requires read1 to overlap the window, which can miss fragments whose span alone overlaps a boundary. |

The common exclusion mask removes unmapped, mate-unmapped, secondary, supplementary, QC-fail and duplicate records. A flag mask alone does not give tools identical pair selection. The end-motif release's ignored length bounds are a separate bug and do not imply that FinaleToolkit's regional coverage ignores them.

For example, forward `pos=100`, `30M`, reverse `pos=170`, `30M` spans 100 bp. Fragment coverage fills the 40 bp inter-mate gap, whereas read depth covers 60 bp. With forward `15M2D15M`, cfDNAlab fragment coverage leaves `[115,117)` uncovered unless supported by the other mate. mosdepth fragment mode fills all 100 bp. General real-BAM equality is therefore not an acceptance condition.

## Run and compare

Set the input/tool paths printed by the installer and use a new results directory. Every run also needs `BENCH_SAMTOOLS_ARCHIVE` for untimed header preparation and Python 3 for environment relocation.

```bash
export BENCH_BAM=/path/to/BH01.hg38.paired.15x.bam
export BENCH_LOCAL_ROOT=/path/to/node-local-scratch
export BENCH_RESULTS=/path/to/coverage-concordance
export BENCH_CORES=1
export BENCH_PILOT=1
export BENCH_METHODS="cfdnalab-fragment cfdnalab-fragment-proper mosdepth-fragment cfdnalab mosdepth samtools-depth pandepth cfdnalab-region-counts cfdnalab-region-counts-proper finaletoolkit-coverage"
sbatch --time=YOUR_PILOT_BASED_BUDGET benchmarks/coverage/benchmark_coverage.sh

# After completion, compare like quantities in their declared common scope
"$BENCH_COMPARE_PYTHON" benchmarks/coverage/compare_coverage.py \
    "$BENCH_RESULTS/profiles" --common-scope > "$BENCH_RESULTS/concordance.tsv"
```

Full compressed cfDNAlab/mosdepth tracks are preserved after the first measured step at the first selected core budget. samtools and PanDepth still generate their complete native outputs **inside timing**. Their retained profiles contain only chr21/chr22 (or 21/22), falling back to the first contig when absent. Override with `BENCH_CAPTURE_CHROMOSOMES=chr21,chr22`. Each profile has an explicit `chrom.sizes`. Preservation/compression of the selected subset happens outside the Slurm tool step.

The comparator requires equal declared scope unless `--common-scope` explicitly allows intersection. Every result records evaluated chromosomes and whether the original scopes match. A missing chromosome in a sparse cfDNAlab track means zero only when declared in its processed scope. Other formats must contain complete runs/positions for every declared contig, so truncated files fail. Runs must be sorted within each contiguous chromosome block, non-overlapping, in bounds and finite/nonnegative. Temporary coalesced run files reconcile arbitrary contig order without allocating per-base arrays. Set `TMPDIR` to suitable scratch for this untimed comparison. zstd reads all concatenated frames.

Base metrics are weighted by genomic bases: identical fraction, MAE/RMSE, mean signed difference and fractions with cfDNAlab/left values above or below the other profile. Nonzero-union metrics exclude bases where both are zero. The 10 kb Pearson correlation uses integrated mean depth, with the actual width of each shorter final bin and equal weight per bin. Constant-profile correlations and empty-union metrics are `NA`. Read-depth and fragment-coverage tracks are not pooled.

Regional comparison requires the exact unique chromosome/start/end rows in the shared `regions.bed`. It reports totals, total ratio, net difference, MAE/RMSE, identical fraction and correlation. Summed regional counts are **not unique-fragment totals**, because a boundary-crossing fragment can count in adjacent windows.

## Production measurements

Use `BENCH_PILOT=1` on the actual production method/core matrix to estimate all work, splitting pilots when necessary. Then set `BENCH_PILOT=0`, two warmups and five measurements, and submit method/core blocks with pilot-based budgets on the same physical node. See [shared run controls and Slurm memory](../README.md#runs-and-budgets). `BENCH_MIN_MAPQ`, `BENCH_MIN_LENGTH` and `BENCH_MAX_LENGTH` override the filters. The timing TSV records applied/unsupported length filtering, proper-pair selection and scope.

Keep process-plus-filesystem-flush time, output bytes and the joined Slurm peak alongside process time and GNU time RSS. Native positional output costs stay in the comparison. Include sample size, SD and median/quartiles and retain all measured points. No timings or real-input concordance values have been measured here.

| Comparison | Scope | Identical bases | MAE / RMSE | Nonzero-union MAE | Pearson r, 10 kb |
| --- | --- | --- | --- | --- | --- |
| Fragment/proper vs mosdepth fragment | [ ] | [ ] | [ ] / [ ] | [ ] | [ ] |
| Native fragment vs mosdepth fragment | [ ] | [ ] | [ ] / [ ] | [ ] | [ ] |
| Paired read depth vs mosdepth | [ ] | [ ] | [ ] / [ ] | [ ] | [ ] |
| Paired read depth vs samtools / PanDepth | [ ] | [ ] | [ ] / [ ] | [ ] | [ ] |

Report regional counts in a separate table and label native/proper diagnostics separately. Only configurations actually timed can support a corresponding speed comparison.
