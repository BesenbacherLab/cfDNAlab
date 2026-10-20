# Outlier command specification

`cfdna outliers` detects smaller sample-specific pileups in raw positional fragment coverage,
after ordinary blacklisting. It writes sparse correction weights, exact and flanked BED files,
and inspection outputs. `fcoverage --outlier-weights` is the current downstream integration.

CLI help is the user-facing option reference. Each run's `outliers.README.txt` explains files
and columns. This spec records the implementation invariants. Unresolved decisions and future
work are in the [remaining-work plan](../future/outliers_remaining_work.md).

## Cargo features

`cmd_outliers` enables the command and its public runner, independently of `cli` and `plotters`.
It requires `cmd_fcoverage` for shared fragment coverage construction and the optional `argmin`,
`argmin-math` and `statrs` dependencies for ZINB fitting. `cli` adds Clap parsing and dispatch.
Equivalent CLI rendering remains available without `cli`.

`plotters` compiles both PNG writers and preserved-context BAM rereads. Without it, all numerical
analysis and text outputs remain available. The normal output manifest contains ten text files,
`OutliersRunResult.output_plot` is `None`, and preserved-region `plot_file` values are `NA`.
With it, `output_plot` is `Some(path)` and the global PNG is included in the manifest.
`--plot-preserved-regions` requires `plotters` and otherwise fails validation before BAM access.
Plot-context coordinates are always recorded. Text metadata records `plotting_available`.

## Coordinates and fragment coverage

- Intervals are 0-based and half-open. Rows follow selected chromosome order, then genomic order.
  Correction and BED intervals do not overlap within a chromosome.
- Paired fragment spans are directional, from `forward.pos` to `reverse.reference_end`.
  Unpaired reads used as fragments span `read.pos` to `read.reference_end`.
- Coverage follows mapped reference segments. Reference `D` and `N` gaps do not contribute.
  `--ignore-gap` controls the inter-mate gap for paired fragments.
- Blacklisted bases are excluded from histograms, fitting, calls, and eligible-position counts.
  A fragment crossing a blacklist can still contribute coverage outside it.
  Positions outside the blacklists remain eligible even when their coverage is zero.
- Shared raw coverage is stored as `f32`. Reject stored coverage at or above 16,777,216.
  That boundary is representable, but the next integer is not and can round down to it.

## Processing and fitting contexts

`--tile-size` controls BAM processing tiles and memory use. `--stride` controls consecutive
model cores, each with a histogram, model, and calling threshold. The first BAM pass builds
core histograms. The second applies fitted models to raw coverage. Connected calls are joined
across both tile and core boundaries, preserving observed and target masses.

Assign accepted fragments once to the core containing their midpoint, excluding blacklisted
midpoints. Partial cores are combined across processing tiles. Missing support indices are errors,
not zero support.

A local fitting context starts with its core and expands by complete cores on both available
sides until it meets both `--bin-size` and `--min-context-fragments`. With 500 kb cores, centered
spans grow as 0.5, 1.5, ..., 5.5 Mb for a 5 Mb minimum. Ends expand on the available side.
A core uses the sample-wide model only when its whole chromosome cannot meet both requirements.
A sufficiently supported local context that cannot be fitted must fail explicitly. Determine
global fallback need from each full chromosome's physical span and summed midpoint support.
Fit the pooled model once only if needed. Keep the pooled histogram for alpha and observations
regardless. Local fitting receives the model family independently of the optional global model.
Collect independent local fit failures across selected chromosomes before returning diagnostics.

`start` and `end` in the model table bound the core where that model is applied.
`context_start`, `context_end`, `eligible_positions`, and `fragments_by_midpoint` describe the
recorded chromosome context, including fallback rows. Fallback parameters come from the
sample-wide histogram, not that recorded context. Overlapping context counts are not additive.

## Statistical model

ZINB is the default. `--model zip` selects ZIP. Both combine a structural-zero probability `pi`
with an active coverage component of mean `mu`. The active variance is
`mu + dispersion * mu^2`. Zero dispersion gives Poisson and therefore ZIP.

ZINB profiles `pi` analytically and fits log mean and log dispersion with three deterministic
Nelder-Mead starts. Compare converged valid candidates with an explicitly fitted ZIP boundary.
A selected boundary is part of ZINB fitting, not a silent fallback after optimization failure.
Numerical failures and fits reaching upper numerical bounds are errors. Preserve stable
log-space likelihoods. ZINB failures identify the fitting stage. Local fitting errors include
the chromosome, core/context, and position/fragment support.

Both models follow the same two-stage procedure:

1. Fit the complete eligible coverage histogram.
2. Find T1, the smallest positive integer satisfying `P(X >= T1) <= alpha`.
3. Exclude observations at or above T1. Fit the likelihood conditioned on `X < T1`, including
   its normalization. Its parameters describe the full underlying distribution.
4. Find T2 using that underlying distribution and the same inclusive survival condition.
5. Call original raw positions at coverage `>= T2`.

Fail explicitly if T1 leaves insufficient positive support. Do not raise T1 silently or bypass
refitting. Preservation is decided after calling and never adds observations back to the fit.

Reject T1=2 for both free-zero mixtures. With retained levels 0 and 1, their underlying
component mean and zero-inflation probability are not separately determined. For ZINB at T1=3
with observed levels 1 and 2, reject a truncated ZIP boundary candidate with positive pi above
`1e-12`. It identifies an interior dispersion/zero-inflation ridge with different full tails.
Keep equality, numerical near-equality and valid constrained boundary fits. These guards do
not establish identification in every larger-support or near-boundary case.

Preflight all initial optimizer vertices before argmin setup. Likelihood errors remain fatal,
including trial errors. Re-evaluate converged and selected candidates for finite likelihood.
ZINB optimizer diagnostics record converged/unfinished counts and parameter/threshold ranges
across valid near-optimal candidates, including the ZIP boundary when comparable. The average
negative-log-likelihood tolerance is recorded (`1e-10`). These are numerical comparisons,
not confidence intervals. An unfinished start with materially better cost remains an error.

Second-fit errors include T1, alpha and first-fit parameters. Smaller alpha may retain more
positive support, but does not fix an initial-fit failure. Never clamp T2 to T1. Persist
`final_threshold_exceeds_initial` as an inspection flag.

Keeping a called elevation uncorrected does not make it the background for detecting other
pileups. Only values at or above T1 are excluded from refitting. Lower values still contribute,
including within regions that are later preserved.

ZIP caches log factorials from 0 through 255 and uses the existing Stirling approximation above
that range. Preserve inclusive tail and predecessor checks at numerical boundaries.

Automatic `alpha = multiplier / eligible_positions` uses eligible positions on selected
chromosomes, including zeros. This bounds the fitted sum of marginal tail probabilities by the
multiplier, not a family-wise error guarantee. Coverage shares fragments between positions,
and fitted parameters are uncertain. The required 0.1x to 200x operating range remains a
validation requirement, not an established calibration result.

### Diagnostic definitions

Model and global-fit outputs retain both fits' means, zero inflation, dispersion, thresholds,
and Poisson-boundary flags. `underlying_mean = (1 - pi) * mu` includes structural zeros.
It differs from the active-component mean used by the `local-mean` target. Neither is the
positive conditional mean `mu / (1 - P(component = 0))`. Ordinary sampling zeros remain part
of the active component at all depths. Fitted structural zeros do not identify specific gaps,
unmappable positions, or deletions. For ZIP specifically, `E[X | X > 0]` is
`mu / (1 - exp(-mu))`, not `mu`.

Observed positive moments describe the complete original positive histogram, including extremes.
Variance is the population mean squared deviation. For either active component:

```text
positive_mean = mu / (1 - P(component = 0))
positive_variance = positive_mean * (1 + dispersion * mu + mu - positive_mean)
positive_coverage_variance_ratio = observed_positive_variance / expected_positive_variance
```

Expected positive variance comes from the second fit. All-zero incomplete contexts using global
fallback have unavailable observed positive moments and ratio (`NaN`).

Observed and expected initial-tail counts compare original bases at `X >= T1` with the first-fit
expectation. Final-tail counts compare original bases at `X >= T2` with the second-fit expectation.
These and the variance ratio are diagnostics, not extra calling rules or independent tests.
Localized excess can reflect pileups. Widespread excess can indicate ordinary-tail model mismatch.

### Rate variation and model limits

For Poisson counts conditional on a varying expected rate, the marginal variance is
`mu + Var(rate) = mu + mu^2 * rate_CV^2`. Here `rate_CV` is the standard deviation of underlying
expected rates divided by their mean, not the CV of observed counts. Sampling variation remains
at every depth. Rate heterogeneity grows in relative importance as depth increases.

NB represents a gamma mixture of Poisson rates. It allows extra variation but does not guarantee
correct extreme tails for every rate distribution. Equal means and variances do not imply equal
tails. Multiplying the tail-probability cutoff by squared depth does not model this variation.
Simulated CV choices are sensitivity assumptions, not estimates of actual library bias.

## Preservation

Defaults are a physical 100 kb window, 50 kb minimum eligible support and denominator, a strictly
greater than 5% cutoff, and `preserve` action. All are configurable.

Start a forward window at every exact call's start and clip to chromosome bounds. Its score is:

```text
100 * called_eligible_bases / max(eligible_bases, minimum_support)
```

The numerator counts the union of all exact candidate bases before flanking or preservation.
The denominator includes eligible zeros. Windows keep physical genomic spans and never extend
across blacklisted gaps to collect more positions.

For each call:

1. If any overlapping forward window meets minimum eligible support, use the maximum only
   among those windows, including windows started by earlier calls.
2. Otherwise evaluate `[max(0, call_end - window_size), call_end)`. If this backward window meets
   minimum support, use its percentage in place of undersupported forward evidence.
3. Otherwise retain the original forward maximum with the denominator floor.

Record the selected window's direction, interval, called bases, and actual eligible support.
Scores strictly above the cutoff flag a call. Equality does not. A 100% cutoff disables flags.

Flank all candidate calls, clip to chromosome bounds, and merge touching or overlapping flanks.
Preserve membership through transitive merges. If any member is flagged, apply the configured
action to every member: `preserve` omits their weights and actionable BED regions, `correct`
retains individual proposed weights, and `error` reports every failing merged region across
selected chromosomes after the second BAM pass before stopping without normal results.
All candidates and decisions remain in the audit when outputs are written.

Do not subtract preserved bases from flanks, pool masses into merged-region weights, or add
fragment-level protection overrides. Preservation leaves all fit parameters and thresholds
unchanged. It is a density heuristic, not biological classification. Narrow biological events
can be corrected and dense technical events preserved.

## Correction weights and application

For each connected exact call:

```text
observed_mass = sum of original positional coverage
target_mass = sum of the selected target at each position
keep_weight = min(1, target_mass / observed_mass)
```

The default `local-mean` target is the second-fit active-component mean. `threshold` uses T2,
and `zero` uses zero. The active-component target represents the intended background at an
active locus. Extra structural zeros elsewhere in the context do not lower that target.
This is a correction choice, not an implication that a positive or threshold-crossing count has
conditional expectation `mu`. Keep the same target definition at shallow and high depth.

Accumulate target mass using the model applied at each position, including when a connected call
crosses cores. Calls selected for correction retain individual weights.

The sparse table has columns `chromosome`, `start`, `end`, and `keep_weight`. Weights are finite
and in `[0, 1]`. Missing overlaps have weight 1. No flanks are added to correction intervals.
`fcoverage` takes the minimum weight over any positive overlap with the original fragment's
`pos` to `reference_end` span, before trimming. Apply that full weight to every counted
contribution from the fragment, without dilution by overlap fraction.

Zero outlier weights skip GC correction, length-normalization work, trimming, coverage insertion,
and counted-fragment accounting. Partial weights multiply GC weights and retain existing
normalization formulas and ownership rules. Exclusion counters count tile encounters.
Downweighted counters follow counted tile insertions, not unique molecules. Zero-weight
exclusion encounters can differ between GC modes because of padding-only fragments. This
reporting discrepancy is deferred and does not establish a coverage difference.

Called-base percentages measure genomic extent. Whole-fragment weighting also affects bases
outside calls, so audit masses and genomic percentages do not measure total coverage reduction.

## Output contracts

An optional output prefix precedes each filename with a dot. Text tables begin with `# key=value`
comments. Genomic intervals use the coordinate convention above.

| File | Required content |
| --- | --- |
| `outliers.keep_weights.tsv` | Actionable individual exact intervals and keep weights. No preserved calls or flanks. |
| `outliers.exact.bed` | The same actionable exact calls, as three BED columns. |
| `outliers.flanked.bed` | Actionable merged flanks, with chromosome clipping. Default flank is the maximum accepted fragment length. |
| `outliers.audit.tsv` | Every candidate's interval, observed/target masses, mean/max coverage, proposed weight, selected window evidence, merged bounds/score, action/reason, and first/last model-core starts. Use `proposed_keep_weight`, so this table cannot be used accidentally as correction input. |
| `outliers.summary.tsv` | Chromosome and total eligible bases, candidate/actionable/preserved exact counts and percentages, longest single calls, merged-region counts, actionable flanked extents, and model/fallback/boundary counts. |
| `outliers.preserved_regions.tsv` | Always-written preserved inspection groups and plot coordinates, including a header when empty. |
| `outliers.preserved_regions.region_000001.png` | Optional local coverage figure for each inspection group. |
| `outliers.models.tsv` | One row per core with context, support, both fits, applied threshold, and diagnostics defined above. |
| `outliers.histograms.tsv.zst` | Sparse per-core raw coverage histograms and midpoint-assigned fragment counts. |
| `outliers.global_fit.tsv` | Complete pooled observed histogram. Expected counts and parameter metadata only when a global fallback fit was needed. Otherwise expected columns are `NA` and `global_fit_status=not_needed`. |
| `outliers.global_fit.png` | Global observed/fitted views and genomic tracks of local T1/T2 and both means. Requires the `plotters` feature. |
| `outliers.README.txt` | File and column explanations with the actual preservation settings. |

Summary percentages use the same row's `eligible_positions` and are `NaN` when it is zero.
Exact bases exclude flanks and gaps. Flanked spans include blacklisted bases and joined gaps,
while eligible flanked extent excludes blacklists. Longest-call columns are maximum individual
exact lengths before flanking. The total row takes maxima for those columns, not sums.

Histograms have columns `chromosome`, `start`, `end`, `fragments_by_midpoint`, `coverage`, and
`observed_positions`. Fragment counts repeat on each coverage row and must not be summed across
those repetitions. Missing coverage levels have no positions. A fully blacklisted core has an
explicit `coverage=0, observed_positions=0` row. Global-fit rows instead include every integer
coverage through the observed maximum, even zero-count gaps. Never truncate the observed tail.

Analysis and report-writing failures publish only distinct `outliers.error.txt` and available
`outliers.error_histogram.tsv` or `outliers.error_audit.tsv`. The human report includes settings,
stage, original errors and companion paths. Histograms are uncompressed sparse counts for failing
contexts. Pooled rows use `GLOBAL` and `NA` coordinates. Ambiguity audits cover every failed merged
region and its member evidence, with window contributors and actions marked `withheld_error`. They omit
unrelated successful decisions. Use already collected evidence without extra BAM scans.

Normal files are staged until all analysis and reporting finish, then published with the shared
`FinalOutputFiles::move_into_place` helper. Publication errors can leave some normal outputs
replaced or newly present. Report any diagnostic-writing failure alongside the original error.
Existing successful files are retained after analysis failures. Error files use fixed names per
prefix and attempt identifiers. Reports list only companions from their attempt. Successful retries
may leave old error files or optional PNGs. Consumers must check the current exit status.

### Preserved inspection groups and plots

Keep each merged preserved region intact. Join groups when the spans containing their selected
preservation windows overlap, including transitive overlaps. Merely touching spans do not join.
Bounds contain exact member calls, including gaps but without flanks. `called_bases` sums exact
lengths. `maximum_window_percentage` is the highest selected member score. That window's numerator
includes all candidates in it, not only group members. These bounds are not biological breakpoints.

`plot_start` and `plot_end` are centered on group bounds with total-width multiplier 10 and
minimum 100 kb by default, then include the selected window span and clip to chromosome bounds.
Plot contexts do not affect grouping. `plot_file` is `NA` unless plots were requested.

`--plot-preserved-regions` rereads only these contexts using unchanged detection filters, in
chunks bounded by `--tile-size` with complete fragment fetch halos. Reads do not enter detection
statistics. Overview and close-up each use at most 2,000 independently chosen reference bins,
accumulated from the same chunks. Retain mean, peak and peak position over eligible positions,
including genuine zeros. Blacklisted spans are unavailable, not zeros. Both panels show active
mean, T2, all visible candidate calls with actual decisions, and masks. Progress identifies each
group and requested reread span. Old optional PNGs are not removed when reusing a prefix.
Reporting changes no detection or correction decisions. Table and optional PNGs use final-output registration.

## Correction combinations and current scope

`fcoverage` multiplies outlier weights with GC-file or GC-tag weights and existing normalization
weights. This does not require reconstructing raw coverage in the consumer or assuming that
the corrections are statistically independent. Producing GC/scaling files from the same
outlier-weighted fragments remains future integration. Residual extremes may affect those
calibration estimates. Their practical influence must be measured after ordinary blacklisting.
When producer integration is added, estimate and apply corrections from a consistent weighted
fragment signal. Preserve current normalization formulas and ownership in the meantime.

Keep additional consumers, GC/scaling producer integration, normalization changes, smarter
preservation evidence, hysteresis, and cohort aggregation outside this inspection milestone.
Move overlap-resolver state out of `fcoverage` only when another consumer needs it.
The user must explicitly request any investigation or change to restore-mean weighting.

## Recorded settings

Common metadata records model/target, actual tail probability and denominator, core/context/tile
settings, fragment filters, raw coverage and gap definitions, blacklists, and preservation settings.
`analysed_chromosomes` contains sorted JSON objects with `name` and `size` for analysed chromosomes
only, including those without calls. Downstream counting may select a subset. Known requested
chromosome or reference-length mismatches fail. Filter differences warn. Custom tables without
recorded settings remain usable with a warning. Do not compare BAM identities or blacklists.
`automatic_probability_denominator` and `tail_probability_multiplier` are `NA` with manual alpha.
Obsolete hysteresis metadata is omitted. `detector` retains the two-stage algorithm identity. Document existing CLI/metadata
mappings until the owner reviews renaming both together.

## Mathematical references

- [Negative-binomial parameterization and survival probabilities](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html)
- [Lambert's zero-inflated Poisson model](https://www.stat.cmu.edu/technometrics/90-00/vol-34-01/v3401001.pdf)
- [Normalization of truncated distributions](https://mc-stan.org/docs/reference-manual/statements.html#truncated-distributions)
