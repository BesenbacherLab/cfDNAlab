# Handle extreme fragment coverage

A short genomic region with extremely high fragment coverage can dominate a coverage profile.
`cfdna outliers` finds unusually high fragment coverage and calculates weights that reduce its
influence during counting.

Use your usual blacklist first. This command then looks for sample-specific pileups over short
genomic intervals outside known problematic regions. It supplements the blacklist.

High coverage can also be biological. The command leaves broad or clustered calls uncorrected
by default, but coverage alone cannot establish their cause. Inspect the results before adopting
the correction in your analysis.

## Run detection

```bash
cfdna outliers \
  --bam sample.bam \
  --output-dir results/outliers \
  --output-prefix sample \
  --n-threads 12 \
  --blacklist blacklist.bed \
  --plot-preserved-regions
```

The model estimates expected coverage from the surrounding region and calls bases that reach
an unusually high coverage cutoff. Adjacent called bases form a call. Positions outside the
supplied blacklists are eligible, including positions with zero coverage. Blacklisted positions
are excluded from detection and the coverage estimate.

The default model is zero-inflated negative binomial (ZINB), which allows extra variation in
coverage. `--model zip` selects the simpler zero-inflated Poisson model for comparison.
Detection uses raw fragment coverage, before GC correction or normalization.

`--plot-preserved-regions` adds local coverage plots and requires extra BAM reads around those
regions. Detection uses two BAM passes. Model fitting and preservation use their recorded
histograms and calls. Plotting adds local rereads, which can become substantial with many
groups. Progress logs identify each group and requested span. The close-up uses finer bins
from the same reread and includes nearby corrected calls as well as preserved calls.
Omit plotting if you only want the tables and the overall diagnostic plot.

Unmapped, secondary, supplementary, duplicate and QC-failed alignments are excluded. Mapping
quality and fragment length filters apply. Pairs need mapped mates on the same chromosome
and opposite strands. Requiring a proper-pair flag is optional. `--reads-are-fragments`
counts each accepted alignment separately. Alignment D/N gaps contribute no coverage, and
`--ignore-gap` also omits the inter-mate gap. These positions remain eligible unless masked.
The detector does not read reference sequence, so unmasked reference gaps remain eligible zeros.

Local contexts expand until their span and accepted midpoint fragment support are sufficient.
A chromosome that cannot provide both uses a pooled model from all selected chromosomes.
No pooled model is fitted when every chromosome has enough support. Automatic alpha depends
on the selected eligible span. Detect once over the intended chromosomes and consume subsets.
Using the same manual alpha in separate runs cannot reproduce a different pooled fallback.

## Inspect the results

The filenames below start with `sample.` when using the example's output prefix.

- **Start with `outliers.summary.tsv`.** It reports how many calls and bases were selected for
  correction or preserved. "Preserved" means left unchanged. Base percentages use the
  number of eligible positions on the analysed chromosomes as their denominator.
- **Open `outliers.global_fit.png`.** The upper panels compare the observed coverage distribution
  with the sample-wide model when a chromosome needs global fallback. Otherwise the plot
  states that the global fit was not needed, and its TSV keeps observed counts with `NA`
  for unavailable fitted counts. The bottom panel shows the local means and cutoffs used along
  chromosomes. Calls spread across much of the genome, or a broad mismatch between observed and
  fitted coverage, deserve investigation before applying the correction.
- **Read `outliers.audit.tsv` for individual calls.** It records their coverage, proposed weight,
  and the reason for correcting or preserving them. `outliers.models.tsv` gives the background
  estimate and cutoff used at each location.
- **Use `outliers.preserved_regions.tsv` to locate groups left uncorrected.** Nearby preserved
  calls are grouped using the windows that informed their decisions. Optional numbered PNGs
  show surrounding coverage and a closer view. These are approximate locations of preserved
  signal, not estimated biological breakpoints.

`outliers.README.txt` accompanies every run with the settings used and explanations of every
file and column. Coordinates follow BED conventions: starts are 0-based and ends are excluded.

### Why a short call can be preserved

The preservation check considers nearby calls together. By default, it leaves calls uncorrected
when more than 5% of a surrounding 100 kb window was called. This helps protect broader coverage
elevations whose called bases may be scattered across many short intervals. Calls connected by
their flanks also share the preservation decision.

The percentage excludes blacklisted bases. Special handling near chromosome ends and large
blacklisted regions prevents a tiny usable span from giving a misleadingly high percentage.
See `--min-protection-positions` in the CLI help for that rule. In a fully eligible default
100 kb window, an isolated contiguous call longer than 5 kb exceeds the cutoff. Exactly
5 kb does not. Smaller eligible spans can lower that boundary toward 2.5 kb at the denominator
floor. Nearby calls and merged flanks can preserve shorter calls. An event's biological width
is not necessarily the number of its bases that cross the calling threshold.

A preserved group is not necessarily an amplification. Broad technical pileups may also be
preserved, and small biological events may still be corrected. The window and percentage
settings control this tradeoff. Lowering `--max-outlier-percentage` preserves more calls.

## Apply correction weights

Use the keep-weight table with `fcoverage`, keeping the detection filters and blacklist:

```bash
cfdna fcoverage \
  --bam sample.bam \
  --output-dir results/coverage \
  --blacklist blacklist.bed \
  --outlier-weights results/outliers/sample.outliers.keep_weights.tsv
```

A keep weight is the fraction of a fragment's contribution to retain. For example, 0.1 retains
10%, while zero excludes the fragment. The default target is the estimated background mean.
Preserved calls are absent from the correction table.

Even a 1 bp overlap applies the full regional weight to the whole fragment, including its
coverage outside the call. If a fragment overlaps several calls, it receives the smallest
weight. The weights use exact called intervals, without flanks.

This can reduce ordinary neighboring signal. For stationary starts, fixed contiguous 30 bp
fragments, a single 1 bp call and weight 0.2, ordinary coverage 30 becomes 6 at the call,
6.8 a base outside it, 18 at distance 15, and 30 at distance 30. This calculation describes
ordinary background under those assumptions. Variable fragment lengths, alignment gaps,
artifact fragments, other calls and other corrections can change the observed shape.

GC correction can be used together with outlier weights. Their weights multiply. For the first
inspection, start without GC correction, genomic scaling, or length normalization so coverage
can be compared directly with the detector's mean and cutoff.

### Use a blacklist instead

The two BED outputs offer exclusion as an alternative to fractional weighting:

- `outliers.exact.bed` contains the exact called bases selected for correction.
- `outliers.flanked.bed` adds a buffer on each side and joins overlapping buffers. The buffer
  defaults to the maximum accepted fragment length and can be changed with `--blacklist-flank`.

Both omit preserved calls. Add either alongside your existing blacklist, for example:

```bash
cfdna fcoverage \
  --bam sample.bam \
  --output-dir results/coverage \
  --blacklist blacklist.bed \
  --blacklist results/outliers/sample.outliers.flanked.bed
```

In `fcoverage`, blacklists mask coverage positions. A fragment can still contribute outside
the masked interval. Outlier weights instead reduce the complete overlapping fragment. The
receiving command's rules determine what is excluded when using BED files elsewhere.

## Check the effect on coverage

Called-base percentages describe genomic extent. They do not measure how much coverage is
removed, because weights affect complete fragments. Compare matched raw and weighted coverage
around the largest pileups, preserved groups, and nearby background.

Create `inspection_loci.bed` with three columns naming the chromosome, start, and end of each
region to inspect. These commands write per-position coverage using the same inputs and filters:

```bash
cfdna fcoverage \
  --bam sample.bam --output-dir results/inspection --output-prefix raw \
  --blacklist blacklist.bed \
  --by-bed inspection_loci.bed --per-window unique-positions --decimals 6

cfdna fcoverage \
  --bam sample.bam --output-dir results/inspection --output-prefix weighted \
  --blacklist blacklist.bed \
  --by-bed inspection_loci.bed --per-window unique-positions --decimals 6 \
  --outlier-weights results/outliers/sample.outliers.keep_weights.tsv
```

If you changed mapping quality, fragment length, chromosome selection, or gap handling during
detection, use the same settings in both coverage commands.

From a repository checkout, `scripts/plot_outliers.R` can turn these tracks into PNG and PDF
comparisons showing raw, weighted, and removed coverage alongside calls and cutoffs:

```bash
Rscript scripts/plot_outliers.R \
  results/inspection/raw.fcoverage.per_position.bedgraph.zst \
  results/inspection/weighted.fcoverage.per_position.bedgraph.zst \
  results/outliers/sample.outliers.models.tsv \
  results/outliers/sample.outliers.audit.tsv \
  inspection_loci.bed results/inspection/locus
```

The script uses base R and the `zstd` executable for compressed inputs. It also accepts
decompressed bedGraphs.

For all detection and preservation settings, run `cfdna outliers --help` or read the
[CLI reference](../generated/cli/outliers.md).
