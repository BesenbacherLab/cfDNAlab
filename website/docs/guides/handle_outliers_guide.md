# Detect Extreme Fragment Coverage Outliers

`cfdna outliers` detects positions with much higher raw fragment coverage than expected for the
sample and genomic region. These extreme pileups can dominate coverage-based analyses, even when
they cover only a small part of the genome.

The command does not determine why a pileup is present. High coverage can result from mismapping,
residual duplicates, difficult reference regions, or a real biological signal. The output is
therefore a sample-specific set of candidate regions to inspect and, when appropriate, blacklist.

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
its other covered positions.

## Use the detected regions

The command writes two BED files so that you can use the calls as blacklists:

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
regions. It is a sparse file, so omitted positions have a weight of `1.0`. These weights are
intended for planned downstream correction and are not currently consumed by other cfDNAlab
commands. Use a BED output when you need a blacklist now.

## Inspect and interpret the calls

Start with `outliers.global_fit.png`. It shows:

1. The observed sample-wide coverage distribution and the fitted reference distributions near the
   coverage values used for outlier detection.
2. The complete observed coverage histogram, including the most extreme coverage values.
3. The expected coverage and the two thresholds used across consecutive genomic regions.

The sample-wide thresholds in the first two panels are references. The thresholds in the bottom
panel are the ones that were actually used in each part of the genome. This matters when coverage
varies across chromosomes or along a chromosome.

Compare the plot with `outliers.exact.bed`. A small number of isolated regions above otherwise
stable thresholds is consistent with the intended use of the detector. Calls spread across large
parts of the genome, or a fitted distribution that differs broadly from the observed histogram,
suggest that the model may not describe that sample well enough for automatic blacklisting.

The output identifies unusually high coverage, not technical noise. Before adopting the calls as a
routine blacklist, inspect representative samples and consider whether the regions overlap known
repeats, problematic reference regions, copy-number changes, or plausible biological signals.

## Output files

| Output | What it is useful for |
| --- | --- |
| `outliers.exact.bed` | Blacklisting or inspecting the exact detected regions |
| `outliers.flanked.bed` | Applying a more conservative blacklist around the detected regions |
| `outliers.keep_weights.tsv` | Fractional correction in a future workflow, currently not consumed by other commands |
| `outliers.global_fit.png` | Checking the complete coverage distribution and how thresholds vary across the genome |
| `outliers.models.tsv` | Finding the fitted expectation, actual threshold, and fit diagnostics for a genomic region |
| `outliers.global_fit.tsv` | Reading the complete sample-wide observed and fitted distributions in another program |
| `outliers.histograms.tsv.zst` | Auditing raw coverage counts and fragment support across consecutive genomic regions |

All coordinates are 0-based and half-open. Text outputs begin with `# key=value` lines that record
the command settings and definitions needed to interpret the columns.

## Detailed model diagnostics

This section is only needed when you want to investigate why a region was called or assess whether
the fitted model is appropriate for your data.

Fragment coverage can contain more zero-coverage positions than an ordinary Poisson model allows,
so the command uses a zero-inflated Poisson (ZIP) model. It first fits the complete coverage
distribution and obtains an initial threshold, `T1`. It then refits the expected coverage using
values below `T1` while accounting for the excluded upper tail. The final threshold, `T2`, is
calculated from this second estimate. Original positions with coverage at or above `T2` are called
as outliers.

The command estimates these values separately across the genome. In `outliers.models.tsv`, `start`
and `end` identify the region where a model and its `final_threshold` were used. `context_start` and
`context_end` identify the broader region whose coverage informed that model.

### Coverage expectations and thresholds

- `underlying_mean` is the expected raw positional coverage after accounting for both the Poisson
  coverage component and excess zeros.

- `initial_threshold` is `T1`. It protects the expected coverage estimate from the extreme values
  that the command is trying to detect.

- `second_fit_retained_positions` is the number of positions below `T1` that were available for the
  second fit.

- `final_threshold` is `T2`, the coverage threshold actually used to call positions in that row's
  genomic region.

- `initial_lambda` and `underlying_lambda` are means of the Poisson component, not the expected
  coverage across all positions. `initial_zero_inflation` and `underlying_zero_inflation` describe
  the probability assigned to the structural-zero component. For a ZIP model, the expected
  coverage is `(1 - zero_inflation) * lambda`.

### Which data informed the model

- `model_source` is `local` when nearby coverage from the same chromosome was used. It is
  `global_fallback` when the complete chromosome could not provide the requested genomic span or
  number of fragments, so the sample-wide model was used instead. Fallback is an audit flag. It
  does not by itself mean that a call is wrong.

- `eligible_positions` is the number of unmasked positions used to estimate the model.

- `fragment_support` is the number of accepted fragments whose midpoint falls between
  `context_start` and `context_end`. Together with `eligible_positions`, it shows how much position
  and fragment data informed the row.

### Positive-coverage variance

The variance columns describe how much the nonzero coverage values vary. For an ordinary Poisson
distribution, the mean and variance are both `lambda`. Because these diagnostics consider only
positions with coverage greater than zero, the expected variance is adjusted for that condition.

- `observed_positive_coverage_mean` and `observed_positive_coverage_variance` summarize the original
  positive coverage values, including the extreme values.

- `underlying_zip_expected_positive_coverage_variance` is the positive-coverage variance expected
  from the fitted model.

- `positive_coverage_variance_ratio` is the observed variance divided by the fitted expectation. A
  value near `1` is consistent with Poisson-like spread among positive values. A value above `1`
  means the observed positive coverage is more variable than the model expects. A value below `1`
  means it is less variable.

The variance ratio is included because ZIP can account for excess zeros but cannot otherwise make
positive coverage more variable than a Poisson distribution. It is a warning about possible model
mismatch, not a calling rule. A high value can also be caused by the pileups being detected because
the observed variance includes the extreme tail.

The observed positive-coverage metrics are `NaN` only when the fitted region contains no positive
coverage values. In that case, the row uses the sample-wide fallback model.

### Observed and expected tail counts

The `observed_initial_tail_positions` and `expected_initial_tail_positions` columns compare the
number of positions at or above `T1` with the fitted expectation. The corresponding `final` columns
make the same comparison at `T2`.

An observed excess at `T2` shows the high-coverage tail that produced the calls. If the observed and
expected tails differ broadly across most genomic regions or across many representative samples,
the ZIP model may be a poor description of ordinary coverage variation in that dataset.

These diagnostics are not p-values. Adjacent positions share fragments and are therefore not
independent observations.

## Statistical references

- [Lambert's definition of zero-inflated Poisson regression](https://www.stat.cmu.edu/technometrics/90-00/vol-34-01/v3401001.pdf)
- [Stan's definitions for truncated discrete distributions](https://mc-stan.org/docs/reference-manual/statements.html#truncated-distributions)
- [NIST's definition of the inclusive discrete survival function](https://www.itl.nist.gov/div898/software/dataplot/refman2/ch8/intro.pdf)
