//! Zero-inflated Poisson fitting and upper-tail thresholds for raw positional coverage.
//!
//! # ZIP model
//!
//! A zero-inflated Poisson (ZIP) variable is a mixture. With probability `pi` it is a
//! structural zero. With probability `1 - pi` it is drawn from `Poisson(lambda)`. Therefore
//!
//! ```text
//! P(X = 0) = pi + (1 - pi) exp(-lambda)
//! P(X = k) = (1 - pi) exp(-lambda) lambda^k / k!, k >= 1
//! E[X] = (1 - pi) lambda
//! ```
//!
//! This is the ZIP definition in [Lambert (1992), pages 3-4][lambert]. It gives two fitting
//! invariants used here. Positive ZIP observations follow a zero-truncated Poisson because the
//! factor `1 - pi` cancels after conditioning on `X >= 1`. Zero inflation changes the probability
//! of zero, but it cannot change the relative shape of the positive Poisson counts.
//!
//! # Two-stage fit
//!
//! The initial model is fitted to the complete histogram. Its inclusive upper-tail threshold
//! `T1` excludes observations with `X >= T1` from the second fit. The second likelihood is
//! conditioned on the retained support `X < T1`, following the standard definition of a
//! [right-truncated discrete distribution][stan-truncation]. Its parameters still describe the
//! underlying untruncated ZIP. The final threshold `T2` is therefore calculated from that
//! untruncated distribution, not from the retained conditional distribution.
//!
//! Both thresholds use the inclusive discrete survival probability `P(X >= k)`. This equals
//! `1 - CDF(k) + PMF(k)`, as defined by the [NIST probability reference][nist-survival]. It is not
//! the point probability `P(X = k)` or the strict complementary CDF `P(X > k)`.
//!
//! # Interpretation limits
//!
//! Lambert's regression model assumes independent responses. Positional fragment coverages are
//! spatially dependent because a fragment contributes to multiple positions. This module fits a
//! shared marginal distribution from histogram counts. It does not model spatial dependence,
//! parameter uncertainty, or positive-count overdispersion. The observed positive tail must
//! therefore be checked against the fitted Poisson tail before treating the thresholds as
//! calibrated.
//!
//! [lambert]: https://www.stat.cmu.edu/technometrics/90-00/vol-34-01/v3401001.pdf
//! [stan-truncation]: https://mc-stan.org/docs/reference-manual/statements.html#truncated-distributions
//! [nist-survival]: https://www.itl.nist.gov/div898/software/dataplot/refman2/ch8/intro.pdf

use std::collections::BTreeMap;

use anyhow::{Context, Result, bail, ensure};

const LAMBDA_BISECTION_ITERATIONS: usize = 96;
const MIN_POSITIVE_LAMBDA: f64 = 1e-12;
const CONDITIONAL_MEAN_TOLERANCE: f64 = 1e-12;
const LOG_TWO_PI: f64 = 1.837_877_066_409_345_3;

/// Sparse eligible-position counts keyed by integer raw coverage.
pub(crate) type CoverageCounts = BTreeMap<u32, u64>;

/// Parameters of an underlying, untruncated zero-inflated Poisson distribution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ZipParameters {
    /// Mean of the Poisson component, not the overall ZIP mean.
    pub(crate) lambda: f64,
    /// Probability of selecting the structural-zero component.
    pub(crate) zero_inflation: f64,
}

impl ZipParameters {
    /// Expected positional coverage, including structural zeros.
    #[inline]
    pub(crate) fn mean(self) -> f64 {
        (1.0 - self.zero_inflation) * self.lambda
    }

    /// Return the probability mass `P(X = coverage)` under this ZIP distribution.
    ///
    /// At zero, this includes both structural zeros and zeros from the Poisson component. At
    /// positive coverage, only the non-structural Poisson component contributes.
    pub(crate) fn probability_mass(self, coverage: u32) -> f64 {
        let poisson_probability = poisson_probability(self.lambda, coverage);
        if coverage == 0 {
            self.zero_inflation + (1.0 - self.zero_inflation) * poisson_probability
        } else {
            (1.0 - self.zero_inflation) * poisson_probability
        }
    }

    /// Return the inclusive positive-tail probability `P(X >= minimum_coverage)`.
    pub(crate) fn survival_at_or_above(self, minimum_coverage: u32) -> Result<f64> {
        if minimum_coverage == 0 {
            return Ok(1.0);
        }
        let poisson_survival = if minimum_coverage as f64 >= self.lambda {
            poisson_log_survival_at_or_above(self.lambda, minimum_coverage)?.exp()
        } else {
            let log_cdf_before_threshold = poisson_log_cdf(self.lambda, minimum_coverage - 1);
            -log_cdf_before_threshold.exp_m1()
        };
        Ok((1.0 - self.zero_inflation) * poisson_survival)
    }
}

/// Observed positive-count dispersion and initial- and final-tail calibration for a fitted ZIP.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ZipFitDiagnostics {
    /// Mean among original observed positions with positive coverage.
    pub(crate) observed_positive_coverage_mean: f64,
    /// Population variance among original observed positions with positive coverage.
    pub(crate) observed_positive_coverage_variance: f64,
    /// Positive-coverage variance expected from the fitted underlying ZIP.
    pub(crate) underlying_zip_expected_positive_coverage_variance: f64,
    /// Observed positive-coverage variance divided by the underlying ZIP expectation.
    pub(crate) positive_coverage_variance_ratio: f64,
    /// Observed positions with coverage greater than or equal to `T1`.
    pub(crate) observed_initial_tail_positions: u64,
    /// Positions expected at or above `T1` under the initial ZIP.
    pub(crate) expected_initial_tail_positions: f64,
    /// Observed positions with coverage greater than or equal to `T2`.
    pub(crate) observed_final_tail_positions: u64,
    /// Positions expected at or above `T2` under the fitted underlying ZIP.
    pub(crate) expected_final_tail_positions: f64,
}

/// Initial and right-truncation-corrected ZIP fits for a coverage histogram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TwoStageZipModel {
    /// Underlying ZIP fitted to the complete, original histogram.
    pub(crate) initial: ZipParameters,
    /// Smallest positive `T1` with `P_initial(X >= T1) <= tail_probability`.
    pub(crate) initial_threshold: u32,
    /// Underlying untruncated ZIP estimated from observations with `X < T1`.
    pub(crate) underlying_fit: ZipParameters,
    /// Smallest positive `T2` with `P_final(X >= T2) <= tail_probability`.
    pub(crate) final_threshold: u32,
    /// Number of histogram observations with `X < T1` used by the second fit.
    pub(crate) second_fit_retained_positions: u64,
}

/// Fit a complete histogram, exclude its called upper tail, and refit the underlying ZIP.
///
/// The second fit excludes every observation with coverage greater than or equal to `T1`. It
/// conditions its likelihood on the remaining support instead of treating the retained histogram
/// as an untruncated sample. No excluded probability mass is moved to `T1`, and the procedure does
/// not iterate.
pub(crate) fn fit_two_stage_zip(
    histogram: &CoverageCounts,
    tail_probability: f64,
) -> Result<TwoStageZipModel> {
    validate_tail_probability(tail_probability)?;
    let initial = fit_zip(histogram)?;
    let initial_threshold = zip_upper_tail_threshold(initial, tail_probability)?;
    let (underlying_fit, second_fit_retained_positions) =
        fit_right_truncated_zip(histogram, initial_threshold)?;
    let final_threshold = zip_upper_tail_threshold(underlying_fit, tail_probability)?;

    Ok(TwoStageZipModel {
        initial,
        initial_threshold,
        underlying_fit,
        final_threshold,
        second_fit_retained_positions,
    })
}

/// Compare observed positive counts and both fitted tails with a two-stage ZIP model.
pub(crate) fn diagnose_zip_fit(
    histogram: &CoverageCounts,
    model: TwoStageZipModel,
) -> Result<ZipFitDiagnostics> {
    let (positive_positions, positive_coverage_sum) = positive_summary(histogram, None)?;
    let underlying_zip_positive_coverage_mean =
        zero_truncated_poisson_mean(model.underlying_fit.lambda);
    let underlying_zip_expected_positive_coverage_variance = underlying_zip_positive_coverage_mean
        * (1.0 + model.underlying_fit.lambda - underlying_zip_positive_coverage_mean);
    ensure!(
        underlying_zip_expected_positive_coverage_variance.is_finite()
            && underlying_zip_expected_positive_coverage_variance > 0.0,
        "ZIP fit produced invalid positive-count variance {}",
        underlying_zip_expected_positive_coverage_variance
    );
    let (
        observed_positive_coverage_mean,
        observed_positive_coverage_variance,
        positive_coverage_variance_ratio,
    ) = if positive_positions == 0 {
        // An incomplete all-zero local context can validly use the global fallback model.
        // Its positive-count diagnostics are undefined and are written as NaN
        (f64::NAN, f64::NAN, f64::NAN)
    } else {
        let observed_positive_coverage_mean = positive_coverage_sum / positive_positions as f64;
        let observed_positive_squared_deviation_sum = histogram
            .iter()
            .filter(|(coverage, _)| **coverage > 0)
            .map(|(&coverage, &positions)| {
                let deviation = coverage as f64 - observed_positive_coverage_mean;
                deviation * deviation * positions as f64
            })
            .sum::<f64>();
        let observed_positive_coverage_variance =
            observed_positive_squared_deviation_sum / positive_positions as f64;
        (
            observed_positive_coverage_mean,
            observed_positive_coverage_variance,
            observed_positive_coverage_variance
                / underlying_zip_expected_positive_coverage_variance,
        )
    };
    let eligible_positions = histogram.values().try_fold(0_u64, |total, &positions| {
        total
            .checked_add(positions)
            .context("eligible-position count overflow")
    })?;
    let observed_initial_tail_positions = observed_tail_positions(
        histogram,
        model.initial_threshold,
        "initial-tail position count overflow",
    )?;
    let expected_initial_tail_positions = eligible_positions as f64
        * model
            .initial
            .survival_at_or_above(model.initial_threshold)?;
    let observed_final_tail_positions = observed_tail_positions(
        histogram,
        model.final_threshold,
        "final-tail position count overflow",
    )?;
    let expected_final_tail_positions = eligible_positions as f64
        * model
            .underlying_fit
            .survival_at_or_above(model.final_threshold)?;

    Ok(ZipFitDiagnostics {
        observed_positive_coverage_mean,
        observed_positive_coverage_variance,
        underlying_zip_expected_positive_coverage_variance,
        positive_coverage_variance_ratio,
        observed_initial_tail_positions,
        expected_initial_tail_positions,
        observed_final_tail_positions,
        expected_final_tail_positions,
    })
}

fn observed_tail_positions(
    histogram: &CoverageCounts,
    threshold: u32,
    overflow_message: &'static str,
) -> Result<u64> {
    histogram
        .range(threshold..)
        .try_fold(0_u64, |total, (_, &positions)| {
            total.checked_add(positions).context(overflow_message)
        })
}

fn validate_tail_probability(tail_probability: f64) -> Result<()> {
    ensure!(
        tail_probability.is_finite() && tail_probability > 0.0 && tail_probability <= 0.5,
        "tail probability must be finite and in (0, 0.5], got {}",
        tail_probability
    );
    Ok(())
}

/// Fit a constant ZIP model to an untruncated histogram.
///
/// The interior maximum-likelihood estimate first matches the positive-count mean to
/// `E[X | X >= 1] = lambda / (1 - exp(-lambda))`, then solves the observed zero fraction for
/// `zero_inflation`. If that solution would require negative zero inflation, the constrained
/// maximum is the ordinary Poisson boundary with `zero_inflation = 0`.
fn fit_zip(histogram: &CoverageCounts) -> Result<ZipParameters> {
    let total_positions = histogram.values().try_fold(0_u64, |total, &positions| {
        total
            .checked_add(positions)
            .context("eligible-position count overflow")
    })?;
    let zero_positions = histogram.get(&0).copied().unwrap_or(0);
    let (positive_positions, positive_coverage_sum) = positive_summary(histogram, None)?;

    ensure!(
        total_positions > 0,
        "cannot fit ZIP to a histogram with no eligible positions"
    );
    ensure!(
        positive_positions > 0,
        "cannot fit ZIP because every eligible position has zero coverage"
    );

    // Conditional on being positive, ZIP values follow a zero-truncated Poisson because the
    // structural-zero probability cancels. This lets us fit lambda before zero inflation
    let positive_coverage_mean = positive_coverage_sum / positive_positions as f64;
    let lambda =
        fit_lambda_from_conditional_mean(positive_coverage_mean, zero_truncated_poisson_mean)?;
    let observed_zero_fraction = zero_positions as f64 / total_positions as f64;
    let poisson_zero_probability = (-lambda).exp();
    let zero_inflation_denominator = 1.0 - poisson_zero_probability;

    // Solve P(X = 0) = pi + (1 - pi) * Poisson(0 | lambda) for pi
    if zero_inflation_denominator > 0.0 {
        let zero_inflation =
            (observed_zero_fraction - poisson_zero_probability) / zero_inflation_denominator;
        if zero_inflation >= 0.0 && zero_inflation < 1.0 {
            return Ok(ZipParameters {
                lambda,
                zero_inflation,
            });
        }
    }

    // ZIP reduces to an ordinary Poisson when the data contain no excess zeros
    let poisson_lambda = positive_coverage_sum / total_positions as f64;
    ensure_valid_parameters(ZipParameters {
        lambda: poisson_lambda,
        zero_inflation: 0.0,
    })
}

/// Estimate an underlying ZIP from observations retained on `0 <= X < exclusion_threshold`.
///
/// Conditioning retained positive values on `1 <= X < exclusion_threshold` removes
/// `zero_inflation`, so their mean identifies `lambda`. The retained zero fraction then identifies
/// `zero_inflation` only after accounting for the Poisson probability of retention. If the implied
/// value would be negative, the constrained fit uses the right-truncated Poisson boundary.
fn fit_right_truncated_zip(
    histogram: &CoverageCounts,
    exclusion_threshold: u32,
) -> Result<(ZipParameters, u64)> {
    ensure!(
        exclusion_threshold > 1,
        "initial ZIP threshold {} leaves no positive coverage values for the second fit",
        exclusion_threshold
    );
    let retained_positions =
        histogram
            .range(..exclusion_threshold)
            .try_fold(0_u64, |total, (_, &positions)| {
                total
                    .checked_add(positions)
                    .context("second-fit retained-position count overflow")
            })?;
    let zero_positions = histogram.get(&0).copied().unwrap_or(0);
    let (positive_positions, positive_coverage_sum) =
        positive_summary(histogram, Some(exclusion_threshold))?;

    ensure!(
        positive_positions > 0,
        "initial ZIP threshold {} removed every positive coverage value from the second fit",
        exclusion_threshold
    );

    let maximum_retained_coverage = exclusion_threshold - 1;
    // The retained positive values are conditioned on 1 <= coverage < exclusion_threshold
    let positive_coverage_mean = positive_coverage_sum / positive_positions as f64;
    let lambda = fit_lambda_from_conditional_mean(positive_coverage_mean, |lambda| {
        truncated_poisson_mean(lambda, 1, maximum_retained_coverage)
    })?;
    let poisson_retention_probability = poisson_log_cdf(lambda, maximum_retained_coverage).exp();
    let poisson_zero_probability = (-lambda).exp();
    let retained_zero_fraction = zero_positions as f64 / retained_positions as f64;

    // Solve the retained-sample zero fraction for zero inflation. Retention changes the
    // denominator because values at or above the exclusion threshold are absent
    let zero_inflation_denominator = (1.0 - poisson_zero_probability)
        - retained_zero_fraction * (1.0 - poisson_retention_probability);

    if zero_inflation_denominator > 0.0 {
        let zero_inflation = (retained_zero_fraction * poisson_retention_probability
            - poisson_zero_probability)
            / zero_inflation_denominator;
        if zero_inflation >= 0.0 && zero_inflation < 1.0 {
            return Ok((
                ZipParameters {
                    lambda,
                    zero_inflation,
                },
                retained_positions,
            ));
        }
    }

    // Use the boundary model when the retained data do not contain excess zeros
    let retained_coverage_mean = positive_coverage_sum / retained_positions as f64;
    let poisson_lambda = fit_lambda_from_conditional_mean(retained_coverage_mean, |lambda| {
        truncated_poisson_mean(lambda, 0, maximum_retained_coverage)
    })?;
    Ok((
        ensure_valid_parameters(ZipParameters {
            lambda: poisson_lambda,
            zero_inflation: 0.0,
        })?,
        retained_positions,
    ))
}

fn positive_summary(
    histogram: &CoverageCounts,
    exclusive_upper_bound: Option<u32>,
) -> Result<(u64, f64)> {
    let mut positive_positions = 0_u64;
    let mut positive_coverage_sum = 0.0_f64;
    for (&coverage, &positions) in histogram.iter().filter(|(coverage, _)| {
        **coverage > 0 && exclusive_upper_bound.is_none_or(|bound| **coverage < bound)
    }) {
        positive_positions = positive_positions
            .checked_add(positions)
            .ok_or_else(|| anyhow::anyhow!("positive-position count overflow"))?;
        positive_coverage_sum += coverage as f64 * positions as f64;
    }
    Ok((positive_positions, positive_coverage_sum))
}

fn ensure_valid_parameters(parameters: ZipParameters) -> Result<ZipParameters> {
    ensure!(
        parameters.lambda.is_finite() && parameters.lambda > 0.0,
        "ZIP fit produced invalid lambda {}",
        parameters.lambda
    );
    ensure!(
        parameters.zero_inflation.is_finite()
            && parameters.zero_inflation >= 0.0
            && parameters.zero_inflation < 1.0,
        "ZIP fit produced invalid zero inflation {}",
        parameters.zero_inflation
    );
    Ok(parameters)
}

/// Fit the Poisson rate from a sample mean observed after truncation.
///
/// For a Poisson distribution retained within a fixed coverage range, the maximum-likelihood
/// estimate of `lambda` makes the model's conditional mean equal the observed conditional mean.
/// For an ordinary untruncated Poisson this is simply `lambda = observed_mean`. Under truncation,
/// the conditional mean is a nonlinear function of `lambda`, so that equality must be solved
/// numerically.
///
/// `conditional_mean_for_lambda` calculates the model mean under the same condition used to
/// calculate `observed_conditional_mean`. For example, the initial ZIP fit compares the observed
/// mean above zero with `E[X | X >= 1]`. The second fit compares the observed positive mean below
/// the exclusion threshold with `E[X | 1 <= X < exclusion_threshold]`.
///
/// The model's conditional mean increases with `lambda`, so the function can use binary search
/// over a continuous range of lambda values. This numeric form of binary search is usually called
/// the bisection method. It first finds a lower lambda whose model mean is too small and an upper
/// lambda whose model mean is large enough. The desired lambda must lie between them. It then
/// repeatedly tests the midpoint and keeps the half that can still contain the solution.
fn fit_lambda_from_conditional_mean(
    observed_conditional_mean: f64,
    conditional_mean_for_lambda: impl Fn(f64) -> f64,
) -> Result<f64> {
    ensure!(
        observed_conditional_mean.is_finite() && observed_conditional_mean >= 0.0,
        "invalid observed conditional mean {} during Poisson fitting",
        observed_conditional_mean
    );

    // Lambda must remain positive. Treat a mean at the numerical lower boundary as lambda ~= 0
    let minimum_conditional_mean = conditional_mean_for_lambda(MIN_POSITIVE_LAMBDA);
    if observed_conditional_mean <= minimum_conditional_mean + CONDITIONAL_MEAN_TOLERANCE {
        return Ok(MIN_POSITIVE_LAMBDA);
    }

    // The minimum positive lambda is a known lower bound. Start the upper-bound search at the
    // observed mean, or at 1 for means below 1
    let mut lower_lambda_bound = MIN_POSITIVE_LAMBDA;
    let mut upper_lambda_bound = observed_conditional_mean.max(1.0);

    // Binary search needs an upper bound known to be above the solution. If the current value
    // still predicts too small a mean, double it and check again. This preliminary exponential
    // search reaches a sufficiently large upper bound quickly without assuming a fixed maximum
    while conditional_mean_for_lambda(upper_lambda_bound) < observed_conditional_mean {
        upper_lambda_bound *= 2.0;
        if !upper_lambda_bound.is_finite() || upper_lambda_bound > u32::MAX as f64 {
            bail!(
                "could not find an upper Poisson lambda bound for observed conditional mean {}",
                observed_conditional_mean
            );
        }
    }

    // The solution is now between the bounds. Binary search tests their midpoint. Because the
    // model mean increases with lambda, a candidate mean that is too small moves the lower bound
    // up. Otherwise the upper bound moves down. Each iteration halves the remaining range
    for _ in 0..LAMBDA_BISECTION_ITERATIONS {
        let candidate_lambda = lower_lambda_bound + (upper_lambda_bound - lower_lambda_bound) / 2.0;
        let candidate_conditional_mean = conditional_mean_for_lambda(candidate_lambda);
        if candidate_conditional_mean < observed_conditional_mean {
            lower_lambda_bound = candidate_lambda;
        } else {
            upper_lambda_bound = candidate_lambda;
        }
    }

    Ok(lower_lambda_bound + (upper_lambda_bound - lower_lambda_bound) / 2.0)
}

/// Return `E[X | X >= 1]` for a Poisson distribution.
#[inline]
fn zero_truncated_poisson_mean(lambda: f64) -> f64 {
    lambda / -(-lambda).exp_m1()
}

/// Return `E[X | minimum_coverage <= X <= maximum_coverage]` for a Poisson distribution.
fn truncated_poisson_mean(lambda: f64, minimum_coverage: u32, maximum_coverage: u32) -> f64 {
    let log_lambda = lambda.ln();

    // A Poisson probability is exp(-lambda) * lambda^coverage / coverage!. The exp(-lambda)
    // factor is identical at every coverage and therefore cancels from the conditional mean
    //
    // Summing lambda^coverage / coverage! directly can overflow even when the final mean is
    // well behaved. Calculate its logarithm instead, then subtract the largest logarithm before
    // converting back with exp(). Every resulting scaled mass is at most 1. The same scale factor
    // appears in the numerator and denominator of the mean, so it also cancels
    let mut maximum_log_unnormalized_mass = f64::NEG_INFINITY;
    for coverage in minimum_coverage..=maximum_coverage {
        let log_unnormalized_mass = coverage as f64 * log_lambda - log_factorial(coverage);
        maximum_log_unnormalized_mass = maximum_log_unnormalized_mass.max(log_unnormalized_mass);
    }

    let mut scaled_mass_sum = 0.0_f64;
    let mut scaled_coverage_sum = 0.0_f64;
    for coverage in minimum_coverage..=maximum_coverage {
        let log_unnormalized_mass = coverage as f64 * log_lambda - log_factorial(coverage);
        let scaled_mass = (log_unnormalized_mass - maximum_log_unnormalized_mass).exp();
        scaled_mass_sum += scaled_mass;
        scaled_coverage_sum += coverage as f64 * scaled_mass;
    }
    scaled_coverage_sum / scaled_mass_sum
}

/// Return `ln(P(X <= maximum_coverage))` using a stable log-sum-exp calculation.
fn poisson_log_cdf(lambda: f64, maximum_coverage: u32) -> f64 {
    let log_lambda = lambda.ln();
    let mut maximum_log_weight = f64::NEG_INFINITY;
    for coverage in 0..=maximum_coverage {
        maximum_log_weight =
            maximum_log_weight.max(coverage as f64 * log_lambda - log_factorial(coverage));
    }
    let scaled_probability_sum = (0..=maximum_coverage)
        .map(|coverage| {
            let log_weight = coverage as f64 * log_lambda - log_factorial(coverage);
            (log_weight - maximum_log_weight).exp()
        })
        .sum::<f64>();
    -lambda + maximum_log_weight + scaled_probability_sum.ln()
}

/// Return the smallest positive `k` with inclusive ZIP survival `P(X >= k) <= tail_probability`.
///
/// Structural zeros cannot contribute when `k >= 1`, so the ZIP survival is the Poisson survival
/// multiplied by `1 - zero_inflation`.
pub(crate) fn zip_upper_tail_threshold(
    parameters: ZipParameters,
    tail_probability: f64,
) -> Result<u32> {
    validate_tail_probability(tail_probability)?;
    ensure_valid_parameters(parameters)?;

    // Structural zeros do not contribute to P(X >= k) for k >= 1. The ZIP upper tail is therefore
    // the Poisson upper tail multiplied by the non-structural fraction
    let equivalent_poisson_tail_probability = tail_probability / (1.0 - parameters.zero_inflation);
    if equivalent_poisson_tail_probability >= 1.0 {
        return Ok(1);
    }
    poisson_upper_tail_threshold(parameters.lambda, equivalent_poisson_tail_probability)
}

/// Return the smallest positive integer `k` for which `P(X >= k) <= tail_probability`.
fn poisson_upper_tail_threshold(lambda: f64, tail_probability: f64) -> Result<u32> {
    let initial_coverage_f64 = lambda.ceil().max(1.0);
    ensure!(
        initial_coverage_f64 < u32::MAX as f64,
        "Poisson lambda {} is too large for u32 coverage thresholds",
        lambda
    );
    let mut threshold = initial_coverage_f64 as u32;
    let log_tail_probability = tail_probability.ln();
    let mut log_survival_at_threshold = poisson_log_survival_at_or_above(lambda, threshold)?;

    // Starting at or above lambda permits stable direct upper-tail summation. If this initial
    // threshold already passes, walk left by adding the preceding point mass in log space
    if log_survival_at_threshold <= log_tail_probability {
        while threshold > 1 {
            let previous_coverage = threshold - 1;
            let log_survival_at_previous_coverage = log_add_exp(
                log_survival_at_threshold,
                poisson_log_probability(lambda, previous_coverage),
            );
            if log_survival_at_previous_coverage > log_tail_probability {
                break;
            }
            threshold = previous_coverage;
            log_survival_at_threshold = log_survival_at_previous_coverage;
        }
        return Ok(threshold);
    }

    // Find a passing upper bound exponentially, then binary-search the monotone survival
    // function. Recomputing each candidate tail in log space avoids subtracting a tiny remaining
    // tail from a value near one
    let mut failing_lower_bound = threshold;
    let mut search_step = 1_u32;
    let mut passing_upper_bound = loop {
        let candidate = failing_lower_bound
            .checked_add(search_step)
            .unwrap_or(u32::MAX);
        let candidate_log_survival = poisson_log_survival_at_or_above(lambda, candidate)?;
        if candidate_log_survival <= log_tail_probability {
            break candidate;
        }
        ensure!(
            candidate < u32::MAX,
            "Poisson coverage threshold exceeds the supported u32 range"
        );
        failing_lower_bound = candidate;
        search_step = search_step.saturating_mul(2);
    };

    while passing_upper_bound - failing_lower_bound > 1 {
        let candidate = failing_lower_bound + (passing_upper_bound - failing_lower_bound) / 2;
        let candidate_log_survival = poisson_log_survival_at_or_above(lambda, candidate)?;
        if candidate_log_survival <= log_tail_probability {
            passing_upper_bound = candidate;
        } else {
            failing_lower_bound = candidate;
        }
    }
    Ok(passing_upper_bound)
}

/// Calculate `ln(P(X >= minimum_coverage))` relative to the first upper-tail term.
///
/// The minimum must be at or above `lambda` so every following relative probability decreases.
/// Keeping the first probability in log space prevents underflow at extreme tail probabilities.
fn poisson_log_survival_at_or_above(lambda: f64, minimum_coverage: u32) -> Result<f64> {
    ensure!(
        minimum_coverage as f64 >= lambda,
        "internal Poisson survival start {} is below lambda {}",
        minimum_coverage,
        lambda
    );
    let mut relative_probability = 1.0_f64;
    let mut relative_tail_sum = 1.0_f64;
    let mut coverage = u64::from(minimum_coverage);

    loop {
        coverage = coverage
            .checked_add(1)
            .ok_or_else(|| anyhow::anyhow!("Poisson survival summation overflow"))?;
        relative_probability *= lambda / coverage as f64;
        relative_tail_sum += relative_probability;
        if relative_probability <= relative_tail_sum * f64::EPSILON {
            break;
        }
    }
    Ok(poisson_log_probability(lambda, minimum_coverage) + relative_tail_sum.ln())
}

#[inline]
fn log_add_exp(first_log_value: f64, second_log_value: f64) -> f64 {
    let larger_log_value = first_log_value.max(second_log_value);
    let smaller_log_value = first_log_value.min(second_log_value);
    larger_log_value + (smaller_log_value - larger_log_value).exp().ln_1p()
}

#[inline]
fn poisson_log_probability(lambda: f64, coverage: u32) -> f64 {
    -lambda + coverage as f64 * lambda.ln() - log_factorial(coverage)
}

fn poisson_probability(lambda: f64, coverage: u32) -> f64 {
    poisson_log_probability(lambda, coverage).exp()
}

/// Calculate `ln(value!)`, using Stirling's series when direct summation would be expensive.
fn log_factorial(value: u32) -> f64 {
    if value < 256 {
        return (2..=value).map(|number| (number as f64).ln()).sum();
    }

    let value = value as f64;
    let inverse = 1.0 / value;
    (value + 0.5) * value.ln() - value + 0.5 * LOG_TWO_PI + inverse / 12.0 - inverse.powi(3) / 360.0
        + inverse.powi(5) / 1260.0
}

#[cfg(test)]
mod tests {
    include!("model_tests.rs");
}
