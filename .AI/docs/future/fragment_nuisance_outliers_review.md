# Review of the fragment nuisance outlier implementation

Date: 2026-08-11

This review covers the current implementation against
[`fragment_nuisance_outliers_plan.md`](fragment_nuisance_outliers_plan.md). The command now performs
detection and writes diagnostics and correction weights. Applying those weights in other commands
remains intentionally out of scope for this stage.

## Assessment

No blocking correctness defect remains in the reviewed implementation. The statistical model,
two-pass coverage flow, adaptive local contexts, recorded fallback behavior, regional weights, and
output schemas agree with the plan. The implementation is also consistent with the surrounding
commands: BAM reading, fragment filtering, coverage construction, and tile processing remain in the
command module, while output writing is separated into an output module.

The command should still be treated as scientifically provisional until its positive-tail fit has
been assessed on representative samples and sequencing depths. That is a model-validation
requirement, not an identified implementation error. The command is not yet wired into downstream
commands, as intended.

## Correctness

### ZIP model and the two fitted distributions

[Lambert's original ZIP paper](https://www.stat.cmu.edu/technometrics/90-00/vol-34-01/v3401001.pdf)
defines a mixture of a structural zero with probability `pi` and a `Poisson(lambda)` value with
probability `1 - pi`:

```text
P(X = 0) = pi + (1 - pi) exp(-lambda)
P(X = k) = (1 - pi) exp(-lambda) lambda^k / k!, k >= 1
E[X] = (1 - pi) lambda
```

The implementation preserves the distinction between the Poisson-component mean `lambda` and the
overall ZIP mean. Conditional on `X > 0`, the factor `1 - pi` cancels, so the positive counts have a
zero-truncated Poisson distribution. The initial fit correctly estimates `lambda` from that
conditional distribution and then estimates `pi` from the zero frequency. When the unconstrained
solution would require negative zero inflation, it correctly uses the ordinary Poisson boundary
`pi = 0`.

The second likelihood retains observations with `X < T1`. The
[Stan truncation definition](https://mc-stan.org/docs/reference-manual/statements.html#truncated-distributions)
requires dividing each retained mass by the probability of the retained support. With
`F = P_Poisson(X < T1)`, `p0 = exp(-lambda)`, and `r` equal to the retained zero fraction, the
right-truncated ZIP solution is:

```text
pi = (r F - p0) / ((1 - p0) - r (1 - F))
```

[`model.rs`](../../../src/commands/outliers/model.rs) uses this likelihood. Its resulting parameters
describe the underlying untruncated distribution, not the retained histogram. Calculating `T2` and
the default regional target from that underlying distribution is therefore correct.

### Inclusive discrete tail

The [NIST probability reference](https://www.itl.nist.gov/div898/software/dataplot/refman2/ch8/intro.pdf)
defines the discrete survival probability at `k` as `P(X >= k) = 1 - CDF(k) + PMF(k)`. Stan likewise
distinguishes an inclusive discrete CDF, `P(X <= k)`, from the strict complementary CDF,
`P(X > k)`. The implementation correctly finds the smallest positive integer satisfying:

```text
P_ZIP(X >= k) = (1 - pi) P_Poisson(X >= k) <= alpha
```

It compares log probabilities, avoiding loss of the extreme tail through subtraction from one. A
regression test fixes the previously exposed case `lambda = 1`, `pi = 0`, `alpha = 1e-100` at the
correct threshold `70`. Other tests check the inclusive boundary, minimality, zero-inflation scaling,
monotonicity, and invalid probabilities.

### Coverage, contexts, and calling

The first BAM pass constructs sparse per-core coverage histograms and midpoint-assigned fragment
support. Positional coverage follows mapped reference segments, so reference `D` and `N` gaps do not
contribute; `--ignore-gap` controls only the inter-mate gap. Blacklisted positions are excluded from
both fitting and calling.

Each core expands its context in complete cores until the configured span and fragment support are
met. The global model is used only when the complete chromosome cannot meet those requirements, and
that fallback is recorded. A complete local context that is degenerate fails explicitly instead of
silently changing models. Calls in the second BAM pass use each core's actual local or recorded
fallback model. The global model and global thresholds are diagnostic references only.

The global diagnostic TSV includes the complete original histogram through the maximum observed
coverage. The plot also shows that complete histogram and a genomic track of every core's actual
`T1`, `T2`, and underlying ZIP mean. Its minimum, arithmetic mean, and maximum annotations summarize
the local thresholds actually used rather than minimum or maximum observed coverage.

Regional mass is accumulated position by position, so a region crossing a model boundary uses the
model assigned at each position. The three targets implement the specified formulas, and sparse
omitted weights mean `1.0`.

### Numerical and accounting safeguards

Coverage histograms are sparse, avoiding allocation proportional to the largest outlier. The
command rejects raw coverage at or above `16,777,216`, where the shared `f32` coverage storage can no
longer represent every integer exactly. Histogram counts, fragment support, and regional mass use
checked accumulation. Conservation checks cover eligible positions, fragment support, and the
global histogram.

The automatic rule `alpha = multiplier / eligible_positions` limits the fitted sum of marginal tail
probabilities to the multiplier by linearity of expectation. It is not a family-wise error guarantee:
adjacent coverages are dependent, the model is fitted to the same observations, and parameter
uncertainty is not included.

## Consistency and simplicity

The command follows the existing CLI, configuration, chromosome selection, fragment filtering,
blacklist, tiling, temporary-output, and run-result conventions. It is included in the default
command feature set and has programmatic configuration and CLI round-trip coverage.

The module boundary is deliberately modest:

- [`outliers.rs`](../../../src/commands/outliers/outliers.rs) owns orchestration, BAM reading,
  filtering, fragment iteration, tile coverage, histogram merging, calling, and tile-result merging.
- [`outputs.rs`](../../../src/commands/outliers/outputs.rs) owns output paths, schemas, metadata, and
  serialization.
- [`model.rs`](../../../src/commands/outliers/model.rs),
  [`striding.rs`](../../../src/commands/outliers/striding.rs), and
  [`plotting.rs`](../../../src/commands/outliers/plotting.rs) contain their focused numerical work.

This matches the neighboring `fcoverage`, `ends`, and `midpoints` commands, whose main command modules
also retain BAM filtering and tile processing. Moving raw tile-coverage construction into a dedicated
helper would reduce consistency without isolating a reusable abstraction. No trait hierarchy or
generic detector framework is justified here.

Names now describe the statistical roles directly: `underlying_fit`,
`second_fit_retained_positions`, `underlying_mean`, and `fragment_support`. The CLI retains `stride`
and `bin-size` because those terms are established across commands, while the specification explains
their exact roles for outlier fitting.

## Tests

The focused tests cover both ZIP fits, ordinary-Poisson boundaries, exact `T1` exclusion, inclusive
tail thresholds, the extreme-probability regression, diagnostic moments, all-zero fallback,
monotonicity, invalid probabilities, sparse histograms, coverage validation, context assignment,
fallback recording, and histogram/support conservation.

The public command test constructs hand-derived BAM fixtures and exercises both passes. It checks:

- all seven output artifacts and their metadata and schemas
- exact model thresholds and diagnostics
- the complete global histogram, including an observed extreme tail
- core, blacklist, tile, and chromosome boundaries
- a fully blacklisted core and a short-chromosome global fallback
- exact and flanked region merging
- underlying-mean, threshold, and zero target weights
- execution with two worker threads

The principal remaining automated coverage gaps are a public end-to-end fixture for unpaired reads
and the inter-mate `--ignore-gap` setting, including mapped reference `D` or `N` gaps, and an explicit
comparison of single-threaded and multithreaded output ordering. Those paths use shared machinery
and have lower risk than the command-specific model and boundary logic, but integration tests would
protect the documented semantics.

Per repository policy, the tests were compiled but not executed during this review. The user should
run them before accepting the implementation.

## Documentation

The public configuration API, command help, statistical functions, output writers, and non-obvious
numerical invariants are documented in the code. In particular, the code documentation records the
ZIP source, truncated-likelihood normalization, inclusive-tail convention, expected-count
interpretation, positive-tail limitation, spatial dependence, and exact-integer coverage limit.

The public [`fragment outliers guide`](../../../website/docs/guides/handle_outliers_guide.md) and
generated CLI reference explain output selection, fragment and coordinate semantics, local versus
global models, fallback behavior, regional weights, and how to interpret every model diagnostic.
The model and global-fit TSVs also carry concise metric definitions in their metadata so detached
outputs remain interpretable. The current [`outliers` specification](../specs/outliers.md) remains
the internal implementation contract rather than the user documentation.

## Remaining scientific validation

ZIP permits excess zeros but leaves the positive counts zero-truncated Poisson. It cannot absorb
positive-count overdispersion merely by changing `pi`. Before the detector becomes an accepted
scientific default or its weights are consumed downstream, representative samples and downsampled
depths should be inspected for:

- observed-to-fitted positive-count variance ratios
- observed and expected counts at `T1` and `T2`
- stability of local thresholds and calls across depth
- persistent local positive-tail lack of fit that would justify a different count model

These quantities are now present in the model TSV and global diagnostics. The remaining decision is
empirical model adequacy, not missing observability in the implementation.

## Verification performed

- `cargo check`
- `cargo check --tests --features "testing cmd_outliers"`
- `cargo check --all-features`
- `cargo fmt --all -- --check`

The all-feature check reports only pre-existing dead-code warnings in unrelated commands.
