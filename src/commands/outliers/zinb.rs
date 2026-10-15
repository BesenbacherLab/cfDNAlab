//! Maximum-likelihood ZINB fits, including correctly normalized right truncation.
//!
//! Zero inflation is profiled analytically. Only log component mean and log dispersion require
//! numerical optimization. The separately fitted ZIP boundary is an actual likelihood candidate,
//! not a fallback for unsuccessful positive-dispersion optimization.

use anyhow::{Context, Result, ensure};
use argmin::{
    core::{CostFunction, Executor, State, TerminationReason, TerminationStatus},
    solver::neldermead::NelderMead,
};
use statrs::function::gamma::ln_gamma;

use super::{
    distribution::CoverageParameters,
    model::{
        CoverageCounts, fit_right_truncated_zip, fit_zip, log_factorial, poisson_log_probability,
    },
};

const MINIMUM_LOG_PARAMETER: f64 = -20.0;
const MAXIMUM_DISPERSION: f64 = 100.0;
const LIKELIHOOD_TOLERANCE: f64 = 1e-12;

/// Sparse retained observations, normalized so optimizer tolerances are independent of span.
#[derive(Clone)]
struct ProfiledHistogram {
    /// Raw coverage levels and their fractions among retained eligible positions in this fit.
    observations: Vec<(u32, f64)>,
    /// Observed zero fraction among those retained positions, including sampling zeros.
    zero_fraction: f64,
    /// Exclusive T1 bound for the second fit, or no coverage cutoff for the complete first fit.
    truncation: Option<u32>,
    /// Numerical bound on log component mean for optimization, not a genomic window or T2.
    maximum_log_mean: f64,
}

impl ProfiledHistogram {
    /// Prepare one supplied fitting histogram, which can describe a local context or the genome.
    fn new(histogram: &CoverageCounts, truncation: Option<u32>) -> Result<Self> {
        let retained = histogram.iter().filter(|(coverage, positions)| {
            **positions > 0 && truncation.is_none_or(|threshold| **coverage < threshold)
        });
        let positions = retained.clone().try_fold(0_u64, |total, (_, &count)| {
            total
                .checked_add(count)
                .context("ZINB retained-position count overflow")
        })?;
        ensure!(
            positions > 0,
            "cannot fit ZINB without retained eligible positions"
        );
        let observations: Vec<_> = retained
            .map(|(&coverage, &count)| (coverage, count as f64 / positions as f64))
            .collect();
        let maximum = observations
            .last()
            .context("ZINB has no retained observations")?
            .0;
        ensure!(
            maximum > 0,
            "cannot fit ZINB because every retained position has zero coverage"
        );
        Ok(Self {
            zero_fraction: histogram.get(&0).copied().unwrap_or(0) as f64 / positions as f64,
            observations,
            truncation,
            maximum_log_mean: (2.0 * (maximum as f64 + 1.0)).max(20.0).ln(),
        })
    }

    /// Negative log likelihood per retained position, with the optimal constrained pi.
    ///
    /// Mean and dispersion describe the underlying active component before T1 exclusion.
    /// Truncation changes the likelihood normalization, rather than the meaning of the parameters.
    fn evaluate(&self, mean: f64, dispersion: f64) -> Result<(f64, f64)> {
        // These probabilities refer to the active NB/Poisson component before zero inflation
        // In truncated fits, log_retained covers 0 <= X < T1 and log_positive covers 1 <= X < T1
        let log_zero = component_log_probability(mean, dispersion, 0);
        let (log_retained, log_positive) = if let Some(threshold) = self.truncation {
            let mut log_positive = f64::NEG_INFINITY;
            // Sum the retained support directly. Subtracting a near-one tail can lose all
            // retained probability while the optimizer explores a mean well above T1
            for coverage in 1..threshold {
                log_positive = log_add(
                    log_positive,
                    component_log_probability(mean, dispersion, coverage),
                );
            }
            (log_add(log_zero, log_positive), log_positive)
        } else {
            (0.0, (-log_zero.exp_m1()).ln())
        };
        ensure!(
            log_retained.is_finite() && log_positive.is_finite(),
            "nonfinite ZINB retained probability"
        );
        let conditional_zero = (log_zero - log_retained).exp();
        if self.zero_fraction <= conditional_zero {
            // This component already predicts at least the observed zero fraction, so adding
            // structural zeros cannot improve the fit. Use the constrained boundary pi = 0
            let log_likelihood: f64 = self
                .observations
                .iter()
                .map(|&(coverage, weight)| {
                    weight * (component_log_probability(mean, dispersion, coverage) - log_retained)
                })
                .sum();
            return Ok((-log_likelihood, 0.0));
        }

        // With interior pi, the retained zero fraction is fitted exactly. Positive values
        // follow the component conditional on 1 <= X < T1 (or X >= 1 in the complete fit)
        let denominator = (1.0 - self.zero_fraction) * -log_zero.exp_m1()
            + self.zero_fraction * log_positive.exp();
        let zero_inflation =
            log_retained.exp() * (self.zero_fraction - conditional_zero) / denominator;
        ensure!(
            zero_inflation.is_finite() && (0.0..1.0).contains(&zero_inflation),
            "invalid profiled structural-zero probability {zero_inflation}"
        );
        let positive_likelihood: f64 = self
            .observations
            .iter()
            .filter(|(coverage, _)| *coverage > 0)
            .map(|&(coverage, weight)| {
                weight * (component_log_probability(mean, dispersion, coverage) - log_positive)
            })
            .sum();
        let log_likelihood = self.zero_fraction * self.zero_fraction.ln()
            + (1.0 - self.zero_fraction) * (-self.zero_fraction).ln_1p()
            + positive_likelihood;
        Ok((-log_likelihood, zero_inflation))
    }
}

impl CostFunction for ProfiledHistogram {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, parameters: &Self::Param) -> Result<f64, argmin::core::Error> {
        if parameters.len() != 2
            || !parameters.iter().all(|value| value.is_finite())
            || !(MINIMUM_LOG_PARAMETER..=self.maximum_log_mean).contains(&parameters[0])
            || !(MINIMUM_LOG_PARAMETER..=MAXIMUM_DISPERSION.ln()).contains(&parameters[1])
        {
            // Out-of-domain trial vertices are rejected, never clamped into fake optima
            return Ok(f64::INFINITY);
        }
        self.evaluate(parameters[0].exp(), parameters[1].exp())
            .map(|(cost, _)| cost)
    }
}

/// Fit the complete or right-truncated histogram with three deterministic dispersion starts.
pub(super) fn fit_zinb(
    histogram: &CoverageCounts,
    truncation: Option<u32>,
) -> Result<CoverageParameters> {
    let problem = ProfiledHistogram::new(histogram, truncation)?;
    // Fit dispersion = 0 on exactly the same retained observations, so the Poisson boundary
    // competes with positive-dispersion fits through their average log likelihoods
    let boundary: CoverageParameters = match truncation {
        Some(threshold) => fit_right_truncated_zip(histogram, threshold)?.0,
        None => fit_zip(histogram)?,
    }
    .into();
    let (mut best_cost, _) = problem.evaluate(boundary.component_mean, 0.0)?;
    let mut best = boundary;
    let initial_log_mean = boundary
        .component_mean
        .ln()
        .clamp(MINIMUM_LOG_PARAMETER + 0.1, problem.maximum_log_mean - 0.1);
    let mut converged_starts = 0;
    // An unfinished start with a materially better likelihood would make accepting another fit
    // misleading. Keep its best cost so the final checks can reject that unresolved comparison
    let mut failed_best_cost = f64::INFINITY;
    let mut outcomes = Vec::new();
    for initial_dispersion in [0.0005_f64, 0.02, 0.3] {
        let initial_log_dispersion = initial_dispersion.ln();
        let solver = NelderMead::new(vec![
            vec![initial_log_mean, initial_log_dispersion],
            vec![initial_log_mean + 0.05, initial_log_dispersion],
            vec![initial_log_mean, initial_log_dispersion + 0.2],
        ])
        .with_sd_tolerance(1e-13)?;
        let result = Executor::new(problem.clone(), solver)
            .configure(|state| state.max_iters(2_000))
            .ctrlc(false)
            .run()
            .context("executing ZINB Nelder-Mead optimization")?;
        let state = result.state();
        outcomes.push(state.get_termination_status().to_string());
        let cost = state.get_best_cost();
        if state.get_termination_status()
            != &TerminationStatus::Terminated(TerminationReason::SolverConverged)
            || !cost.is_finite()
        {
            failed_best_cost = failed_best_cost.min(cost);
            continue;
        }
        converged_starts += 1;
        let parameters = state
            .get_best_param()
            .context("converged ZINB optimizer has no parameters")?;
        if cost < best_cost - LIKELIHOOD_TOLERANCE {
            let component_mean = parameters[0].exp();
            let dispersion = parameters[1].exp();
            let (_, zero_inflation) = problem.evaluate(component_mean, dispersion)?;
            best = CoverageParameters {
                component_mean,
                dispersion,
                zero_inflation,
            };
            best_cost = cost;
        }
    }
    ensure!(
        converged_starts > 0,
        "all positive-dispersion ZINB starts failed: {}",
        outcomes.join(", ")
    );
    ensure!(
        failed_best_cost >= best_cost - 1e-10,
        "an unconverged ZINB start improved materially on every accepted fit: {}",
        outcomes.join(", ")
    );
    ensure!(
        best.dispersion < MAXIMUM_DISPERSION * 0.999,
        "fitted ZINB dispersion {} reached its numerical upper bound",
        best.dispersion
    );
    ensure!(
        best.dispersion == 0.0 || best.component_mean < problem.maximum_log_mean.exp() * 0.999,
        "fitted ZINB mean reached its numerical upper bound"
    );
    Ok(best)
}

/// Stable NB2 log mass, including the exact Poisson boundary.
///
/// A Stirling difference evaluates the rising factorial for large shape without subtracting
/// almost equal log-gamma values. Small shapes use log-gamma directly. This needs no dense array
/// extending to the largest observed pileup.
pub(super) fn component_log_probability(mean: f64, dispersion: f64, coverage: u32) -> f64 {
    if dispersion == 0.0 {
        return poisson_log_probability(mean, coverage);
    }
    let product = dispersion * mean;
    let log_product = product.ln_1p();
    let log_zero = -mean * (log_product / product);
    if coverage == 0 {
        return log_zero;
    }
    let count = coverage as f64;
    let shape = 1.0 / dispersion;
    let log_rising = if shape >= 16.0 {
        let ratio = dispersion * count;
        shape * log_one_plus_minus_value(ratio)
            + (count - 0.5) * ratio.ln_1p()
            + stirling_correction(shape + count)
            - stirling_correction(shape)
    } else {
        ln_gamma(shape + count) - ln_gamma(shape) - count * shape.ln()
    };
    log_zero + log_rising - log_factorial(coverage) + count * (mean.ln() - log_product)
}

/// Evaluate log(1 + value) - value without cancellation near zero.
fn log_one_plus_minus_value(value: f64) -> f64 {
    if value >= 0.01 {
        return value.ln_1p() - value;
    }
    let mut power = -value * value;
    let mut sum = power / 2.0;
    for denominator in 3..=12 {
        power *= -value;
        sum += power / denominator as f64;
    }
    sum
}

/// Stirling correction through inverse ninth power, accurate for arguments >= 16.
fn stirling_correction(value: f64) -> f64 {
    let inverse = 1.0 / value;
    inverse / 12.0 - inverse.powi(3) / 360.0 + inverse.powi(5) / 1260.0 - inverse.powi(7) / 1680.0
        + inverse.powi(9) / 1188.0
}

fn log_add(first: f64, second: f64) -> f64 {
    if first == f64::NEG_INFINITY {
        return second;
    }
    if second == f64::NEG_INFINITY {
        return first;
    }
    let maximum = first.max(second);
    maximum + (first.min(second) - maximum).exp().ln_1p()
}

#[cfg(test)]
mod tests {
    include!("zinb_tests.rs");
}
