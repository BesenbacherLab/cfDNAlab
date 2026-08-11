# Outlier command specification

`cfdna outliers` detects sample-specific regions with extreme uncorrected positional fragment
coverage. It writes correction weights and diagnostic artifacts but does not currently apply those
weights in another command.

## Coordinates, fragments, and ordering

All genomic coordinates are 0-based and half-open. Paired fragment spans are defined directionally
from `forward.pos` to `reverse.reference_end`. Unpaired reads used as fragments span `read.pos` to
`read.reference_end`.

Positional coverage follows mapped reference segments. Reference `D` and `N` gaps do not contribute
coverage. The `--ignore-gap` option controls only whether the inter-mate gap contributes for paired
fragments. Input-blacklisted positions are excluded from fitting and calling, although a fragment
crossing them can still contribute at its other mapped positions.

Rows follow the user-selected chromosome order, then increasing genomic coordinates. Intervals in
the correction and BED outputs are non-overlapping within a chromosome.

## Detection model

Each selected chromosome is divided into consecutive model cores of `--stride` bases, with a
shorter final core when necessary. Each model core receives a two-stage zero-inflated Poisson
model. A ZIP variable is a structural zero with probability `pi` and a `Poisson(lambda)` value
otherwise:

```text
P(X = 0) = pi + (1 - pi) exp(-lambda)
P(X = k) = (1 - pi) exp(-lambda) lambda^k / k!, k >= 1
E[X] = (1 - pi) lambda
```

The first fit uses the complete context histogram. `T1` is the smallest positive integer satisfying
the inclusive tail condition `P(X >= T1) <= alpha`. Coverage values at or above `T1` are excluded
from the second likelihood. The second fit is conditioned on the retained support `X < T1`, but its
parameters describe the underlying untruncated ZIP. `T2` is calculated from that underlying model
and is the smallest positive integer satisfying `P(X >= T2) <= alpha`. Original coverage, not the
tail-excluded histogram, is called at `X >= T2`.

Local contexts expand in complete model cores until both the minimum span and midpoint-assigned
fragment support are present. A core uses the global model only if its complete chromosome cannot
meet those requirements. A complete local context that cannot be fitted is an error rather than an
unrecorded fallback. The global model and its thresholds are diagnostic references. Calls use the
model recorded for each core.

The automatic `alpha = multiplier / eligible_positions` limits the fitted sum of marginal tail
probabilities to the multiplier. This is an expected-count interpretation, not a family-wise error
guarantee. Positional coverages are spatially dependent, fitted parameters are uncertain, and ZIP
cannot model positive-count overdispersion. The model and global-fit outputs provide positive-count
variance and initial- and final-tail calibration summaries that must be assessed on representative
data.

Raw coverage is stored as `f32` by the shared coverage implementation. The command rejects coverage
at or above 16,777,216 because larger consecutive integers are not exactly representable in that
type.

Processing tiles and model cores are different. `--tile-size` controls the large BAM intervals used
only for parallel processing and memory limits. `--stride` controls the smaller model cores that
receive histograms, local fits, and thresholds. The first BAM pass builds core histograms within
each processing tile. The second pass applies the fitted core models within each processing tile.
Results are joined across both kinds of boundary.

## Interpreting model metrics

The model TSV is intended to answer three questions: which model was used, what threshold did it
produce, and how well did its positive distribution and upper tail match the observations.

### Model and threshold columns

- `model_source` is `local` when the core's chromosome supplied the required fitting context and
  `global_fallback` otherwise. Frequent fallback means the configured span or fragment-support
  requirement is unsuitable for those chromosomes.
- `eligible_positions` and `fragment_support` describe the context used for that core. They help
  distinguish a stable fit from a small or heavily masked context.
- `initial_lambda` is the mean of the Poisson component in the complete first fit.
  `initial_zero_inflation` is its structural-zero probability.
- `initial_threshold` is `T1`. Original observations at or above it are omitted from the second
  likelihood so that the extreme tail does not determine the baseline model.
- `second_fit_retained_positions` is the number of original observations below `T1` available to
  that likelihood.
- `underlying_lambda` and `underlying_zero_inflation` describe the untruncated ZIP estimated by the
  second fit. `underlying_mean = (1 - underlying_zero_inflation) * underlying_lambda` is its expected
  positional coverage, including zeros.
- `final_threshold` is `T2`, the actual threshold used to call original positions in that core.

### Positive-coverage spread

Variance measures the average squared spread of coverage values around their mean. It is not the
uncertainty of a fitted parameter. For an ordinary `Poisson(lambda)` variable, both the mean and
variance equal `lambda`. The diagnostic conditions on coverage being positive, however. Positive
ZIP observations follow a zero-truncated Poisson with:

```text
positive_mean = lambda / (1 - exp(-lambda))
positive_variance = positive_mean * (1 + lambda - positive_mean)
```

- `observed_positive_coverage_mean` and `observed_positive_coverage_variance` are population moments
  calculated from all original positive coverage values in the fitting context, including values
  in the extreme tail.
- `underlying_zip_expected_positive_coverage_variance` is the positive-coverage variance expected from
  `underlying_lambda` under the second-fit ZIP.
- `positive_coverage_variance_ratio` is observed variance divided by the underlying ZIP expectation.
  A value near `1` is Poisson-like spread, a value above `1` is more variable, and a value below `1`
  is less variable.

The ratio is a coarse model diagnostic, not a hypothesis test, and it does not affect fitting or
calling. Because its observed side includes the extreme values being detected, a high ratio can
reflect a sparse pileup, general positive-count overdispersion, or both. It cannot distinguish those
explanations by itself.

### Tail calibration

- `observed_initial_tail_positions` and `expected_initial_tail_positions` compare the original count
  at or above `T1` with the expectation from the initial ZIP.
- `observed_final_tail_positions` and `expected_final_tail_positions` compare the original count at
  or above `T2` with the expectation from the underlying second-fit ZIP.

These pairs are more direct checks of the tail used for exclusion and calling. Large observed excess
at `T2` is expected in a context containing a genuine pileup, but widespread or systematic excess
across samples can indicate that ZIP underestimates ordinary positive-tail variation. Adjacent
coverage values share fragments, so neither the variance ratio nor the tail comparisons are
independent-observation tests or p-values.

A practical audit is to check fallback use and context support first, then inspect the variance
ratio and observed-versus-expected tail counts, and finally use the global plot and complete
histogram to determine whether discrepancies are localized pileups or broad lack of fit.

## Regional correction weight

Contiguous called positions are joined even across model-core and tile boundaries. For a called
region:

```text
observed_mass = sum of original positional coverage
target_mass = sum of the selected target at each position
keep_weight = min(1, target_mass / observed_mass)
```

The target is the core's underlying ZIP mean, its final threshold `T2`, or zero according to
`--target`. Because models can change within a region, target mass is accumulated position by
position with the model actually used there.

The keep-weight TSV is sparse. Omitted positions have weight `1.0`. The planned downstream rule is
to assign a fragment the minimum weight among regions overlapping its complete `pos` to
`reference_end` span. No downstream command implements that rule yet.

## Output files

An optional output prefix is followed by a dot before each suffix below.

### Which output should be used?

The BED files are the immediately usable exclusion outputs. They can be passed to cfDNAlab
commands, or other tools, that accept BED blacklists. The receiving command determines whether a
blacklisted interval masks positions, fragments, windows, or another command-specific unit.

- Use `outliers.exact.bed` to blacklist only positions that actually crossed their local calling
  threshold. It preserves the detector's exact evidence and is the less aggressive choice.
- Use `outliers.flanked.bed` when the blacklist should conservatively include nearby positions that
  may be affected by fragments associated with the pileup or when interval-boundary uncertainty is
  undesirable.
- Use `outliers.keep_weights.tsv` for fractional correction once downstream fragment-weight loading
  is implemented. Unlike a BED blacklist, it is intended to reduce contributions rather than remove
  them completely.

The remaining outputs are diagnostics. The histogram and model TSVs support detailed auditing, the
global-fit TSV supports programmatic assessment, and the PNG provides quick visual quality control.

### `outliers.keep_weights.tsv`

The primary sparse correction file has columns:

```text
chromosome  start  end  keep_weight
```

Weights are finite and lie in `[0, 1]`. The metadata records the implicit omitted weight and planned
minimum-over-overlaps fragment rule. This file is the planned input for fractional fragment
correction, but no command consumes it yet.

### `outliers.exact.bed`

Three BED columns contain the exact connected called regions:

```text
chromosome  start  end
```

This file preserves the exact threshold-crossing positions. It can be used as a blacklist when the
user does not want the detector to infer an exclusion distance beyond the observed call.

### `outliers.flanked.bed`

The exact regions are expanded on both sides by `--blacklist-flank`, clipped to chromosome bounds,
and merged when expanded intervals overlap or touch. If the option is omitted, the maximum accepted
fragment length is used as the flank. This is the more conservative blacklist: it is useful when
nearby positions may also be influenced by fragments associated with the extreme pileup. Flanking
by the maximum accepted fragment length covers the greatest possible extension of an accepted
fragment beyond an exact call. The flank is an operational exclusion buffer, not a claim that the
added positions independently crossed the outlier threshold. Flanking does not change the exact
BED or the keep-weight intervals.

### `outliers.histograms.tsv.zst`

The sparse per-core histograms have columns:

```text
chromosome  start  end  fragment_support  coverage  observed_positions
```

`fragment_support` is the number of accepted fragments with an unblacklisted midpoint, assigned
once to the model core containing that midpoint. `observed_positions` is the number of eligible
positions at that exact integer raw coverage. Missing coverage values have zero observed positions.
The sum of `observed_positions` within a core is its eligible-position count. A fully blacklisted
core is retained as an explicit `coverage=0, observed_positions=0` row so its interval and zero
fragment support remain visible.

Use this file to verify the raw coverage distribution and support underlying a core, or to perform
independent model checks without rereading the BAM. It is not a blacklist or correction input.

### `outliers.models.tsv`

Each model core has a row with these columns:

```text
chromosome
start
end
context_start
context_end
model_source
eligible_positions
fragment_support
initial_lambda
initial_zero_inflation
initial_threshold
second_fit_retained_positions
underlying_lambda
underlying_zero_inflation
underlying_mean
final_threshold
observed_positive_coverage_mean
observed_positive_coverage_variance
underlying_zip_expected_positive_coverage_variance
positive_coverage_variance_ratio
observed_initial_tail_positions
expected_initial_tail_positions
observed_final_tail_positions
expected_final_tail_positions
```

`model_source` is `local` or `global_fallback`. `eligible_positions` and `fragment_support` describe
the recorded fitting context. `initial_threshold` is `T1`, `final_threshold` is `T2`, and
`second_fit_retained_positions` counts observations with `X < T1`. The positive-coverage metrics and
their interpretation are defined above. The initial and final tail columns compare observed and
fitted counts at `X >= T1` and `X >= T2`, respectively.
Observed positive moments and their ratio are `NaN` when an all-zero incomplete context uses the
global fallback model.

Use this file to see the actual model and thresholds applied in each part of the genome, identify
global fallback use, and assess positive-count dispersion and tail calibration.

### `outliers.global_fit.tsv`

Metadata comments record the global model parameters, both thresholds, positive-count dispersion,
and initial- and final-tail calibration. Data columns are:

```text
coverage  observed_positions  initial_fitted_positions  underlying_fitted_positions
```

Rows include every integer coverage from zero through the maximum observed coverage, including
zero-count gaps. The fitted columns are expected position counts under the global initial and
underlying models. This file never truncates the observed extreme tail at a model threshold.

Use this file for scripts that compare the complete global observed distribution with both fitted
distributions and their tail calibration.

### `outliers.global_fit.png`

The plot contains the detailed global fitted distribution, the complete original histogram through
the maximum observed coverage, and a chromosome-position track of every core's actual `T1`, `T2`,
and underlying fitted mean. The global thresholds are reference annotations. Separate minimum,
arithmetic mean, and maximum annotations summarize the local thresholds used across cores.

Use this plot for rapid visual quality control. It is a summary, while the TSV outputs retain the
exact values needed for detailed or automated assessment.

## Shared text metadata

Every text output begins with `# key=value` metadata. It records the package and command, input BAM,
selected chromosomes, coordinate convention, model and target, actual tail probability and its
mode, automatic multiplier and denominator, core and context settings, fragment-support definition,
tile size, fragment filters, raw coverage definition, inter-mate gap behavior, blacklist paths, and
blacklist flank. The BED files use legal comment lines before their data rows.

## Current scope

The command implements detection, diagnostics, regional weights, and exact and flanked BED outputs.
It does not implement hysteresis, a downstream weight loader, correction in other commands, or
cohort-level merging of sample-specific calls.
