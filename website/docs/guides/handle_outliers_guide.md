# Detect Extreme Fragment Coverage Outliers

`cfdna outliers` detects positions with much higher raw fragment coverage than expected for the
sample and genomic region. These extreme pileups can dominate coverage-based analyses, even when
they cover only a small part of the genome.

Use an existing blacklist to exclude known problematic regions, then use `outliers` to detect
remaining sample-specific extremes. Apply the existing blacklist during detection and downstream
counting. The detected regions supplement that blacklist.

The default model is zero-inflated negative binomial (ZINB), which permits more ordinary coverage
variation than Poisson. Use `--model zip` to compare ZIP with the same preservation and weighting
rules. Both models target unusually high coverage, rather than zeros or low coverage.

The command does not determine why a pileup is present. High coverage can result from mismapping,
residual duplicates, difficult reference regions, or a real biological signal. The output is
therefore a sample-specific set of candidate regions. Broad or clustered calls are automatically
withheld from correction by the preservation rule below.

## Basic command

```bash
cfdna outliers \
  --bam <sample>.bam \
  --output-dir <sample_directory>/outliers \
  --output-prefix <sample_id> \
  --n-threads 12 \
  --blacklist <path>/hg38-blacklist.v2.bed
```

See the [`cfdna outliers` CLI reference](../generated/cli/outliers.md) for every option.

Positions covered by an input blacklist do not influence the fitted coverage expectations and are
not reported as outliers. A fragment that crosses a blacklisted position can still contribute at
its other covered positions. The eligible-position count used for the automatic tail probability
also excludes blacklisted bases. Unblacklisted positions with zero coverage remain eligible.

## Use the detected regions

The command writes two actionable BED files. Preserved candidates are omitted from both:

- `outliers.exact.bed` contains only positions whose coverage reached the outlier threshold. Use
  this when you want the least aggressive blacklist.

- `outliers.flanked.bed` adds a buffer around each exact region. Use this when fragments associated
  with the pileup may also affect nearby positions. The default buffer on each side is the maximum
  accepted fragment length. You can change it with `--blacklist-flank`.

The added flanks are a practical exclusion buffer. They are not additional outlier calls. The
command that receives either BED file determines whether the intervals exclude positions,
fragments, windows, or another command-specific unit.

For example, the flanked calls can be added alongside an existing blacklist:

```bash
cfdna fcoverage \
  ... \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --blacklist <sample_directory>/outliers/<sample_id>.outliers.flanked.bed
```

`outliers.keep_weights.tsv` provides fractional weights as an alternative to removing complete
regions. It is a sparse file, so omitted positions have a weight of `1.0`. `fcoverage` currently
consumes these weights. Any positive overlap of a fragment's complete `pos` to `reference_end`
span applies the full regional weight to the whole fragment, including contributions outside
the interval. If a fragment overlaps multiple regions, it receives the smallest keep weight.
Support in other commands is planned.

For fractional correction of the remaining extremes, retain the existing blacklist:

```bash
cfdna fcoverage \
  ... \
  --blacklist <path>/hg38-blacklist.v2.bed \
  --outlier-weights <sample_directory>/outliers/<sample_id>.outliers.keep_weights.tsv
```

## Automatic preservation

The first version targets smaller residual pileups. Every exact call starts a physical 100 kb
window. Each call receives the largest called-base score from evaluated forward windows it
overlaps, including earlier calls' windows. The score counts exact called bases before flanking,
divided by eligible bases, with a 50 kb denominator floor. Eligible sampling zeros are included.

When its forward window has fewer than 50 kb eligible bases and no supported earlier window
overlaps the call, the command checks a backward window ending at the exact call's end. If this
has sufficient eligible support, it replaces the undersupported forward evidence. Otherwise the
forward score retains its denominator floor. Windows never grow across blacklisted gaps to
collect eligible positions. The audit records the selected direction and actual support.

A score strictly above 5% flags a call. All candidates are flanked, then overlapping/touching
flanks are merged. If any member is flagged, the whole merged region is preserved and all its
member calls are withheld from actionable weights and BEDs. Its original evidence remains in
the audit. This does not change the baseline fit. Merging does not create pooled weights.

Change these settings with `--protection-window-size`, `--min-protection-positions`,
`--max-outlier-percentage`, and `--ambiguous-region-action`. Lower percentage cutoffs preserve
more regions. `correct` applies proposed corrections despite a flag. `error` stops before
publishing outputs. The default action is `preserve`. Changing `--blacklist-flank` can join calls
into a region with a shared action.

This is a density heuristic. Preserved elevations are not proven amplifications. A narrow
biological event can still be corrected, and a dense technical event can be preserved.

## Inspect and interpret the calls

Start with `outliers.global_fit.png`. It shows:

1. The observed sample-wide coverage distribution and the fitted reference distributions near the
   coverage values used for outlier detection.
2. The complete observed coverage histogram, including the most extreme coverage values.
3. The expected coverage and the two thresholds used across consecutive genomic regions.

The sample-wide thresholds in the first two panels are references. The thresholds in the bottom
panel are the ones that were actually used in each part of the genome. This matters when coverage
varies across chromosomes or along a chromosome.

Read `outliers.summary.tsv` for candidate, actionable, and preserved exact-base percentages on
analysed eligible positions. Inspect `outliers.audit.tsv` for proposed weights, selected windows,
and reasons. Then compare the plot with `outliers.exact.bed`. A small number of isolated regions above otherwise
stable thresholds is consistent with the intended use of the detector. Calls spread across large
parts of the genome, or a fitted distribution that differs broadly from the observed histogram,
suggest that the model may not describe that sample well enough for automatic blacklisting.

The output identifies unusually high coverage, not technical noise. Before adopting the calls as a
routine blacklist, inspect representative samples and consider whether the regions overlap known
repeats, problematic reference regions, copy-number changes, or plausible biological signals.

## First BAM inspection workflow

Run detection with the existing blacklist first. The default correction target is the active
component's fitted mean. `outliers.models.tsv` also reports the unconditional mean including
structural zeros. A regional weight applies to complete overlapping fragments, so exact called
percentages do not measure removed coverage or the effect on a feature set.

Create a three-column `inspection_loci.bed` with regions around actionable calls, preserved
candidates where available, and nearby background. From a repository checkout, these commands
produce matched raw and weighted coverage for those loci:

```bash
sample_bam="sample.bam"
known_blacklist="blacklist.bed"

cfdna outliers \
  --bam "$sample_bam" --output-dir results/outliers --output-prefix sample \
  --blacklist "$known_blacklist" \
  --min-mapq 30 --min-fragment-length 30 --max-fragment-length 1000

cfdna fcoverage \
  --bam "$sample_bam" --output-dir results/inspection --output-prefix raw \
  --blacklist "$known_blacklist" \
  --min-mapq 30 --min-fragment-length 30 --max-fragment-length 1000 \
  --by-bed inspection_loci.bed --per-window unique-positions --decimals 6

cfdna fcoverage \
  --bam "$sample_bam" --output-dir results/inspection --output-prefix weighted \
  --blacklist "$known_blacklist" \
  --min-mapq 30 --min-fragment-length 30 --max-fragment-length 1000 \
  --by-bed inspection_loci.bed --per-window unique-positions --decimals 6 \
  --outlier-weights results/outliers/sample.outliers.keep_weights.tsv

Rscript scripts/plot_outliers.R \
  results/inspection/raw.fcoverage.per_position.bedgraph.zst \
  results/inspection/weighted.fcoverage.per_position.bedgraph.zst \
  results/outliers/sample.outliers.models.tsv \
  results/outliers/sample.outliers.audit.tsv \
  inspection_loci.bed results/inspection/locus
```

The helper uses base R and the `zstd` executable for compressed inputs. Decompressed bedGraphs
are also accepted. It saves PNG and PDF panels showing raw/weighted coverage, active and
unconditional means, T2, removed coverage, and candidate decisions. The displayed removed
percentage is based on integrated coverage over the plotted locus. It is not a fraction of
unique molecules or a genome-wide estimate.

Keep chromosome selection, mapping/fragment filters, gap policy, and blacklists consistent.
Start without GC correction, genomic scaling, or length normalization so model thresholds and
coverage have the same units. GC weights can subsequently be multiplied with outlier weights.
Applying outlier scaling during GC/scaling estimation is a future integration.

For a model comparison, rerun detection in a separate directory with `--model zip`, keeping
other settings fixed. Zero-weight fragments are skipped before correction and counting.
Partial weights retain the existing normalization behavior.

## Output files

| Output | What it is useful for |
| --- | --- |
| `outliers.audit.tsv` | Every candidate, proposed weight, selected window support, and final decision |
| `outliers.summary.tsv` | Candidate/actionable/preserved genomic extents and model/fallback counts |
| `outliers.exact.bed` | Blacklisting or inspecting the exact detected regions |
| `outliers.flanked.bed` | Applying a more conservative blacklist around the detected regions |
| `outliers.keep_weights.tsv` | Fractional fragment correction with `fcoverage --outlier-weights` |
| `outliers.global_fit.png` | Checking the complete coverage distribution and how thresholds vary across the genome |
| `outliers.models.tsv` | Finding the fitted expectation, actual threshold, and fit diagnostics for a genomic region |
| `outliers.global_fit.tsv` | Reading the complete sample-wide observed and fitted distributions in another program |
| `outliers.histograms.tsv.zst` | Auditing raw coverage counts and fragment support across consecutive genomic regions |

All coordinates are 0-based and half-open. Text outputs begin with `# key=value` lines that record
the command settings and definitions needed to interpret the columns.

## Detailed model diagnostics

This section is only needed when you want to investigate why a region was called or assess whether
the fitted model is appropriate for your data.

Fragment coverage can contain structural zeros and extra variation in positive counts.
ZINB combines a structural-zero component with a negative-binomial component having mean `mu`
and variance `mu + dispersion * mu^2`. ZIP is available through `--model zip` and is the
zero-dispersion boundary. It first fits the complete coverage
distribution and obtains an initial threshold, `T1`. It then refits the expected coverage using
values below `T1` while accounting for the excluded upper tail. The final threshold, `T2`, is
calculated from this second estimate. Original positions with coverage at or above `T2` are called
as outliers.

The command estimates these values separately across the genome. In `outliers.models.tsv`, `start`
and `end` identify the region where a model and its `final_threshold` were used. `context_start` and
`context_end` identify the broader region whose coverage informed that model.

### Coverage expectations and thresholds

- `underlying_mean` is the expected raw positional coverage after accounting for both the active
  coverage component and excess zeros.

- `initial_threshold` is `T1`. It protects the expected coverage estimate from the extreme values
  that the command is trying to detect.

- `second_fit_retained_positions` is the number of positions below `T1` that were available for the
  second fit.

- `final_threshold` is `T2`, the coverage threshold actually used to call positions in that row's
  genomic region.

- `initial_component_mean` and `underlying_component_mean` are active-component means, not unconditional means
  across all positions. `initial_zero_inflation` and `underlying_zero_inflation` describe
  the probability assigned to the structural-zero component. The unconditional mean is
  `(1 - zero_inflation) * component_mean`. `local-mean` correction uses the component mean.

### Which data informed the model

- `model_source` is `local` when nearby coverage from the same chromosome was used. It is
  `global_fallback` when the complete chromosome could not provide the requested genomic span or
  number of fragments, so the sample-wide model was used instead. Fallback is an audit flag. It
  does not by itself mean that a call is wrong.

- `eligible_positions` is the number of unmasked positions in the recorded context. For fallback
  rows, the applied model was fitted globally, rather than from that chromosome context.

- `fragment_support` is the number of accepted fragments whose unblacklisted midpoint falls between
  `context_start` and `context_end`. Together with `eligible_positions`, it shows how much position
  and fragment data informed the row.

### Positive-coverage variance

The variance columns describe how much nonzero coverage varies. The expected variance uses
the selected model conditional on coverage being positive. ZINB includes its fitted dispersion.
`initial_dispersion` and `underlying_dispersion` record both stages, with corresponding
`*_poisson_boundary` flags. Dispersion zero is also a valid result when ZINB was requested.

- `observed_positive_coverage_mean` and `observed_positive_coverage_variance` summarize the original
  positive coverage values, including the extreme values.

- `underlying_expected_positive_coverage_variance` is the positive-coverage variance expected
  from the fitted model.

- `positive_coverage_variance_ratio` is the observed variance divided by the fitted expectation. A
  value near `1` is consistent with model-expected spread among positive values. A value above `1`
  means the observed positive coverage is more variable than the model expects. A value below `1`
  means it is less variable.

The variance ratio is a diagnostic, not a calling rule or a parameter uncertainty estimate.
A high value can reflect either model mismatch or the pileups being detected, since its observed
side includes the extreme tail.

The observed positive-coverage metrics are `NaN` only when the fitted region contains no positive
coverage values. In that case, the row uses the sample-wide fallback model.

### Observed and expected tail counts

The `observed_initial_tail_positions` and `expected_initial_tail_positions` columns compare the
number of positions at or above `T1` with the fitted expectation. The corresponding `final` columns
make the same comparison at `T2`.

An observed excess at `T2` shows the high-coverage tail that produced the calls. If the observed and
expected tails differ broadly across most genomic regions or across many representative samples,
the selected model may be a poor description of ordinary coverage variation in that dataset.

These diagnostics are not p-values. Adjacent positions share fragments and are therefore not
independent observations.

## Statistical references

- [Negative-binomial parameterization and survival probabilities](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html)


- [Lambert's definition of zero-inflated Poisson regression](https://www.stat.cmu.edu/technometrics/90-00/vol-34-01/v3401001.pdf)
- [Stan's definitions for truncated discrete distributions](https://mc-stan.org/docs/reference-manual/statements.html#truncated-distributions)
- [NIST's definition of the inclusive discrete survival function](https://www.itl.nist.gov/div898/software/dataplot/refman2/ch8/intro.pdf)
