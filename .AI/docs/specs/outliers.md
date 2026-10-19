# Outlier command specification

`cfdna outliers` detects sample-specific regions with extreme uncorrected positional fragment
coverage. It writes correction weights and diagnostic artifacts. `fcoverage --outlier-weights`
currently consumes the keep weights. Support in other commands is planned.

The intended workflow excludes known problematic regions with an existing blacklist before
detecting the remaining sample-specific extremes. The calls and correction weights supplement
that blacklist. The automatic tail-probability denominator counts only unblacklisted positions
on selected chromosomes, including positions with zero coverage.

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

Each selected chromosome is divided into consecutive `--stride` cores. The default model is
zero-inflated negative binomial (ZINB). `--model zip` selects the retained analytic ZIP fitter.
Both use a structural-zero mixture probability `pi`. The active component has mean `mu` and
variance `mu + dispersion * mu^2`. Dispersion zero gives Poisson and therefore ZIP.

ZINB profiles `pi` analytically and optimizes log mean and log dispersion with three deterministic
Nelder-Mead starts. It compares converged candidates with the explicitly fitted ZIP boundary.
Numerical failure is an error, not an automatic model switch. Fits reaching the upper numerical
bounds also fail. Log-space likelihood calculations preserve small retained probabilities.

The complete context histogram supplies the initial fit. T1 is the smallest positive integer with
`P(X >= T1) <= alpha`. Exclude observations at or above T1, then maximize the likelihood conditioned
on `X < T1`. Its parameters describe the underlying untruncated distribution. T2 is calculated
from that distribution using the same inclusive tail condition. Call original raw coverage at
`X >= T2`. Preservation never adds observations back or starts another baseline fit.

Local contexts expand in complete model cores until both the minimum span and midpoint-assigned
fragment support are present. A core uses the global model only if its complete chromosome cannot
meet those requirements. A complete local context that cannot be fitted is an error rather than an
unrecorded fallback. The global model and its thresholds are diagnostic references. Calls use the
model recorded for each core.

The automatic `alpha = multiplier / eligible_positions` limits the fitted sum of marginal tail
probabilities to the multiplier. This is an expected-count interpretation, not a family-wise error
guarantee. Positional coverages are spatially dependent and fitted parameters are uncertain.
ZINB permits positive-count overdispersion but still assumes a particular extreme-tail shape. The model and global-fit outputs provide positive-count
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
- `eligible_positions` and `fragments_by_midpoint` describe the context used for that core. They help
  distinguish a stable fit from a small or heavily masked context.
- `initial_component_mean` is the active-component mean in the complete first fit.
  `initial_zero_inflation` is its structural-zero probability.
- `initial_threshold` is `T1`. Original observations at or above it are omitted from the second
  likelihood so that the extreme tail does not determine the baseline model.
- `second_fit_retained_positions` is the number of original observations below `T1` available to
  that likelihood.
- `underlying_component_mean` and `underlying_zero_inflation` describe the untruncated model estimated by the
  second fit. `underlying_mean = (1 - underlying_zero_inflation) * underlying_component_mean` is its expected
  positional coverage, including zeros.
- `final_threshold` is `T2`, the actual threshold used to call original positions in that core.

### Positive-coverage spread

Variance measures the average squared spread of coverage values. It is not fitted-parameter
uncertainty. Conditional on being positive, both models exclude the structural-zero component:

```text
positive_mean = mu / (1 - P(component = 0))
positive_variance = positive_mean * (1 + dispersion * mu + mu - positive_mean)
```

`observed_positive_coverage_mean` and `observed_positive_coverage_variance` describe the complete
original positive histogram, including extremes. `underlying_expected_positive_coverage_variance`
is the second-fit model's expectation. Their `positive_coverage_variance_ratio` is near 1 for
comparable spread. It is diagnostic only. Sparse pileups and general model mismatch can both
increase it. All-zero incomplete contexts using global fallback have unavailable observed
positive metrics (NaN).

`initial_dispersion`, `underlying_dispersion`, and the corresponding `*_poisson_boundary` flags
record both fits. A selected ZINB model can fit its Poisson boundary without changing the requested
model family.

### Tail calibration

- `observed_initial_tail_positions` and `expected_initial_tail_positions` compare the original count
  at or above `T1` with the expectation from the initial model.
- `observed_final_tail_positions` and `expected_final_tail_positions` compare the original count at
  or above `T2` with the expectation from the underlying second-fit model.

These pairs are more direct checks of the tail used for exclusion and calling. Large observed excess
at `T2` is expected in a context containing a genuine pileup, but widespread or systematic excess
across samples can indicate that the chosen model underestimates ordinary positive-tail variation. Adjacent
coverage values share fragments, so neither the variance ratio nor the tail comparisons are
independent-observation tests or p-values.

A practical audit is to check fallback use and context support first, then inspect the variance
ratio and observed-versus-expected tail counts, and finally use the global plot and complete
histogram to determine whether discrepancies are localized pileups or broad lack of fit.

## Preservation decisions

Defaults are a 100 kb physical reference window, 50 kb minimum eligible support/denominator floor,
5% cutoff, and `preserve` action. All are configurable. Every exact call starts a forward window.
Its score is `100 * called_eligible_bases / max(eligible_bases, minimum_support)`. If any overlapping
forward window has at least the minimum eligible support, each call uses the maximum only among
those supported windows, including earlier windows. Undersupported windows cannot override them.
The numerator is the union of exact called bases. Eligible zeros count in the denominator.

If no overlapping forward window has minimum support, evaluate
`[max(0, call_end - window_size), call_end)`. Called bases are already known from exact
intervals. If this backward window has enough eligible support, its score replaces undersupported
forward evidence for that call. Otherwise retain the original forward maximum with the denominator
floor. Record the selected window's direction, span, called bases, and actual eligible support.
Windows are clipped at chromosome ends and never extended to collect eligible bases across gaps.

Scores strictly above the cutoff flag a call. Flank all candidates, merge touching/overlapping
flanks, and propagate a shared decision through each merged region. `preserve` withholds all member
calls from correction and actionable BED files if any member is flagged. `correct` applies the
individual proposed weights despite the flag. `error` stops before publishing outputs. Preserved
candidates remain in the audit. A 100% cutoff disables density flags.

The rule targets smaller residual pileups. It is a coverage-density heuristic, not a biological
classification. Narrow biological events can be corrected and dense technical events preserved.

## Regional correction weight

Contiguous called positions are joined even across model-core and tile boundaries. For a called
region:

```text
observed_mass = sum of original positional coverage
target_mass = sum of the selected target at each position
keep_weight = min(1, target_mass / observed_mass)
```

The local-mean target is the core's active-component mean. Other `--target` choices retain
its final threshold T2 or use zero. Because models can change within a region, target mass is accumulated position by
position with the model actually used there.

Selected regions retain individual exact weights. Merging does not pool masses or weight flanks.
The keep-weight TSV is sparse. Omitted positions have weight `1.0`. `fcoverage` assigns a fragment
the minimum weight among regions overlapping its complete `pos` to `reference_end` span.
Any positive overlap applies the full regional weight to all counted contributions from that
fragment, without dilution by overlap fraction.

## Output files

An optional output prefix is followed by a dot before each suffix below.

`outliers.README.txt` accompanies each run with a file guide, column definitions, and the actual
preservation settings. Human-readable diagnostic explanations belong there rather than in repeated
setting-like metadata comments.

### Which output should be used?

The BED files contain only actionable regions. Preserved candidates are retained in the audit.
The BED files are the immediately usable exclusion outputs. They can be passed to cfDNAlab
commands, or other tools, that accept BED blacklists. The receiving command determines whether a
blacklisted interval masks positions, fragments, windows, or another command-specific unit.

- Use `outliers.exact.bed` to blacklist only positions that actually crossed their local calling
  threshold. It preserves the detector's exact evidence and is the less aggressive choice.
- Use `outliers.flanked.bed` when the blacklist should conservatively include nearby positions that
  may be affected by fragments associated with the pileup or when interval-boundary uncertainty is
  undesirable.
- Use `outliers.keep_weights.tsv` with `fcoverage --outlier-weights` for fractional fragment
  correction. A zero keep weight excludes overlapping fragments completely.

The remaining outputs are diagnostics. The histogram and model TSVs support detailed auditing, the
global-fit TSV supports programmatic assessment, and the PNG provides quick visual quality control.

### `outliers.keep_weights.tsv`

The primary sparse correction file has columns:

```text
chromosome  start  end  keep_weight
```

Weights are finite and lie in `[0, 1]`. The metadata records the implicit omitted weight and
minimum-over-overlaps fragment rule. This file is consumed by `fcoverage --outlier-weights`.

### `outliers.exact.bed`

Three BED columns contain the actionable exact connected called regions:

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
added positions independently crossed the outlier threshold. Exact coordinates and individual weights remain unchanged for selected calls. Flanking can change
which calls share a preservation decision and are therefore withheld.

### `outliers.audit.tsv` and `outliers.summary.tsv`

The audit contains every exact candidate, observed/target masses, mean/max coverage,
`proposed_keep_weight`, selected window evidence, merged bounds/maximum score, action, and reason.
First/last model-core starts identify relevant rows in the model table. Connected calls can span
multiple fits. The audit deliberately lacks the operational `keep_weight` column.

The summary contains per-chromosome and total eligible positions, candidate/actionable/preserved
call counts and exact bases/percentages, longest calls, merged-region counts, actionable flanked
span and eligible extent, and model/fallback/Poisson-boundary core counts. Every percentage uses
that row's eligible positions. Percentages are NaN for zero eligible support. Flanked span includes
masked bases, while flanked eligible extent excludes them. These are genomic extents, not
fragment counts or removed coverage mass.
The `longest_candidate_call_bp`, `longest_actionable_call_bp`, and `longest_preserved_call_bp`
columns are maximum single exact call lengths before flanking. The total row takes their maxima
over chromosomes. Short calls can be preserved through a shared window score or merged flanks.

### `outliers.preserved_regions.tsv` and optional PNGs

The preserved-region table is always written, including a header when empty. Columns are
`region_id`, `chromosome`, `start`, `end`, `calls`, `called_bases`, `maximum_window_percentage`,
`window_start`, `window_end`, `plot_start`, `plot_end`, and `plot_file`.
Each existing merged preserved region remains intact. Join groups transitively when the spans
containing their selected windows overlap. Merely touching window spans remain separate.
Group bounds contain exact member calls, including gaps but without flanks. `called_bases` sums
only exact lengths. These are inspection groups, not inferred biological breakpoints.

Plot contexts are centered on grouped call bounds with configurable total-width multiplier 10
and minimum span 100 kb. Include the selected window span and clip to chromosome bounds.
The context is recorded even without plots. `plot_file` is `NA` unless plots were requested.
`--plot-preserved-regions` writes a numbered PNG per group. Indexed BAM reads cover only plot
contexts, processed in chunks bounded by `--tile-size` with the same filters and fragment halos.
Plot reads do not affect detection counters. At most 2,000 reference bins retain mean and peak
coverage over unblacklisted bases, including true zeros. Masks are unavailable, not zeros.
Context and close-up panels show the active-component mean, T2, preserved calls, and masked gaps.
Table and PNGs participate in final-output registration. Reporting does not change corrections.

### `outliers.histograms.tsv.zst`

The sparse per-core histograms have columns:

```text
chromosome  start  end  fragments_by_midpoint  coverage  observed_positions
```

`fragments_by_midpoint` is the number of accepted fragments with an unblacklisted midpoint, assigned
once to the model core containing that midpoint. `observed_positions` is the number of eligible
positions at that exact integer raw coverage. Missing coverage values have zero observed positions.
The fragment count repeats across coverage rows within a core and must not be summed across them.
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
fragments_by_midpoint
initial_component_mean
initial_zero_inflation
initial_threshold
second_fit_retained_positions
underlying_component_mean
underlying_zero_inflation
underlying_mean
final_threshold
observed_positive_coverage_mean
observed_positive_coverage_variance
underlying_expected_positive_coverage_variance
positive_coverage_variance_ratio
observed_initial_tail_positions
expected_initial_tail_positions
observed_final_tail_positions
expected_final_tail_positions
initial_dispersion
underlying_dispersion
initial_poisson_boundary
underlying_poisson_boundary
```

`model_source` is `local` or `global_fallback`. `eligible_positions` and `fragments_by_midpoint` describe
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
and both active-component and unconditional fitted means. The global thresholds are reference annotations. Separate minimum,
arithmetic mean, and maximum annotations summarize the local thresholds used across cores.

Use this plot for rapid visual quality control. It is a summary, while the TSV outputs retain the
exact values needed for detailed or automated assessment.

## Shared text metadata

Every text table begins with `# key=value` metadata. It records the package and command, input BAM,
coordinate convention, model and target, actual tail probability and its mode, automatic multiplier
and denominator, core and context settings, tile size, fragment filters, raw coverage definition,
inter-mate gap behavior, blacklist flank, and preservation settings.
`analysed_chromosomes` contains the sorted analysed names and lengths as JSON objects with `name`
and `size` fields. This replaces separate name-only and reference-length records. Counting can use
a subset of these chromosomes. Each requested chromosome must be recorded with a matching length.
`blacklists` records the supplied blacklist paths. Diagnostic explanations are in the output README.
The BED files use legal comment lines before their data rows.

## Current scope

The command implements detection, diagnostics, regional weights, and exact and flanked BED outputs.
`fcoverage` implements downstream weight loading and application. Zero weights skip correction,
length-normalization work, and counting. Partial weights retain existing normalization semantics.
Recorded reference/chromosome mismatches fail, filter differences warn, and custom tables without
recorded settings remain accepted with a warning. No BAM or blacklist identity is compared. Hysteresis, correction in
additional commands, and cohort-level merging of sample-specific calls remain unimplemented.
