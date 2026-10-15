//! Shared parameters, two-stage fitting, and diagnostics for ZIP and ZINB coverage models.

use anyhow::{Context, Result, ensure};
use statrs::function::beta::checked_beta_reg;

use super::{
    config::OutlierModel,
    model::{
        CoverageCounts, TwoStageZipModel, ZipParameters, fit_two_stage_zip,
        observed_tail_positions, positive_summary,
    },
    zinb::{component_log_probability, fit_zinb},
};

/// Parameters describe the underlying distribution, including after truncated fitting.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CoverageParameters {
    /// Mean of the active Poisson or negative-binomial component, including its sampling zeros.
    pub(crate) component_mean: f64,
    /// Extra component variance divided by squared component mean. Zero is Poisson.
    pub(crate) dispersion: f64,
    /// Probability of selecting a structural zero rather than the active component.
    pub(crate) zero_inflation: f64,
}

impl From<ZipParameters> for CoverageParameters {
    fn from(parameters: ZipParameters) -> Self {
        Self {
            component_mean: parameters.lambda,
            dispersion: 0.0,
            zero_inflation: parameters.zero_inflation,
        }
    }
}

impl CoverageParameters {
    /// Unconditional mean, including structural zeros.
    pub(crate) fn mean(self) -> f64 {
        if self.dispersion == 0.0 {
            return ZipParameters {
                lambda: self.component_mean,
                zero_inflation: self.zero_inflation,
            }
            .mean();
        }
        (1.0 - self.zero_inflation) * self.component_mean
    }

    /// Mean conditional on positive coverage, independent of structural-zero inflation.
    pub(crate) fn positive_mean(self) -> f64 {
        self.component_mean
            / -component_log_probability(self.component_mean, self.dispersion, 0).exp_m1()
    }

    /// Probability mass at an integer raw coverage value.
    pub(crate) fn probability_mass(self, coverage: u32) -> f64 {
        if self.dispersion == 0.0 {
            return ZipParameters {
                lambda: self.component_mean,
                zero_inflation: self.zero_inflation,
            }
            .probability_mass(coverage);
        }
        let component =
            component_log_probability(self.component_mean, self.dispersion, coverage).exp();
        (1.0 - self.zero_inflation) * component
            + if coverage == 0 {
                self.zero_inflation
            } else {
                0.0
            }
    }

    /// Inclusive upper tail. The NB beta argument is calculated without subtracting from 1.
    pub(crate) fn survival_at_or_above(self, minimum_coverage: u32) -> Result<f64> {
        if minimum_coverage == 0 {
            return Ok(1.0);
        }
        if self.dispersion == 0.0 {
            return ZipParameters {
                lambda: self.component_mean,
                zero_inflation: self.zero_inflation,
            }
            .survival_at_or_above(minimum_coverage);
        }
        let product = self.dispersion * self.component_mean;
        let failure_probability = product / (1.0 + product);
        let component_tail = checked_beta_reg(
            minimum_coverage as f64,
            1.0 / self.dispersion,
            failure_probability,
        )
        .context("evaluating the inclusive negative-binomial upper tail")?;
        let tail = (1.0 - self.zero_inflation) * component_tail;
        ensure!(
            tail.is_finite() && (0.0..=1.0).contains(&tail),
            "invalid coverage survival probability {tail}"
        );
        Ok(tail)
    }
}

/// Complete and right-truncated fits used for one global or local coverage context.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TwoStageCoverageModel {
    pub(crate) family: OutlierModel,
    /// First fit to all eligible observations in the supplied local context or global histogram.
    pub(crate) initial: CoverageParameters,
    /// Smallest count with an initial inclusive tail at or below the requested probability.
    pub(crate) initial_threshold: u32,
    /// Underlying distribution inferred from observations below T1 with truncation accounted for.
    pub(crate) underlying_fit: CoverageParameters,
    /// Inclusive T2 cutoff applied to original positional coverage during the calling pass.
    pub(crate) final_threshold: u32,
    /// Eligible positions below T1 in this fitting histogram, including zero-coverage positions.
    pub(crate) second_fit_retained_positions: u64,
}

impl From<TwoStageZipModel> for TwoStageCoverageModel {
    fn from(model: TwoStageZipModel) -> Self {
        Self {
            family: OutlierModel::Zip,
            initial: model.initial.into(),
            initial_threshold: model.initial_threshold,
            underlying_fit: model.underlying_fit.into(),
            final_threshold: model.final_threshold,
            second_fit_retained_positions: model.second_fit_retained_positions,
        }
    }
}

/// Fit the selected model coherently through both stages without adding preserved values back.
pub(crate) fn fit_two_stage(
    histogram: &CoverageCounts,
    probability: f64,
    family: OutlierModel,
) -> Result<TwoStageCoverageModel> {
    if family == OutlierModel::Zip {
        return fit_two_stage_zip(histogram, probability).map(Into::into);
    }
    let initial = fit_zinb(histogram, None).context("initial complete ZINB fit")?;
    let initial_threshold = upper_tail_threshold(initial, probability)?;
    ensure!(
        initial_threshold > 1,
        "initial ZINB threshold {initial_threshold} at probability {probability} leaves no positive second-fit support. Choose a smaller tail probability"
    );
    // Excluding the initial tail changes the fitted background only. Preservation is decided
    // after calling, and preserved candidates are never added back for another baseline fit
    let underlying_fit = fit_zinb(histogram, Some(initial_threshold))
        .with_context(|| format!("right-truncated ZINB refit below T1={initial_threshold}"))?;
    let final_threshold = upper_tail_threshold(underlying_fit, probability)?;
    let second_fit_retained_positions =
        histogram
            .range(..initial_threshold)
            .try_fold(0_u64, |total, (_, &positions)| {
                total
                    .checked_add(positions)
                    .context("retained-position count overflow")
            })?;
    Ok(TwoStageCoverageModel {
        family,
        initial,
        initial_threshold,
        underlying_fit,
        final_threshold,
        second_fit_retained_positions,
    })
}

/// Find the smallest positive integer whose inclusive tail is no greater than the cutoff.
pub(crate) fn upper_tail_threshold(
    parameters: CoverageParameters,
    probability: f64,
) -> Result<u32> {
    ensure!(
        probability.is_finite() && probability > 0.0 && probability <= 0.5,
        "tail probability must be finite and in (0, 0.5], got {probability}"
    );
    ensure!(
        parameters.component_mean.is_finite()
            && parameters.component_mean > 0.0
            && parameters.component_mean < u32::MAX as f64,
        "invalid component mean for threshold search"
    );
    let mut lower = 1_u32;
    let mut upper = parameters.component_mean.ceil().max(1.0) as u32;
    while parameters.survival_at_or_above(upper)? > probability {
        upper = upper
            .checked_mul(2)
            .context("coverage threshold exceeds u32 representation")?;
    }
    while lower < upper {
        let midpoint = lower + (upper - lower) / 2;
        if parameters.survival_at_or_above(midpoint)? <= probability {
            upper = midpoint;
        } else {
            lower = midpoint + 1;
        }
    }
    Ok(lower)
}

/// Observed positive-count dispersion and initial- and final-tail calibration for a fitted model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FitDiagnostics {
    /// Mean among original observed positions with positive coverage.
    pub(crate) observed_positive_coverage_mean: f64,
    /// Population variance among original observed positions with positive coverage.
    pub(crate) observed_positive_coverage_variance: f64,
    /// Positive-coverage variance expected from the fitted underlying model.
    pub(crate) underlying_expected_positive_coverage_variance: f64,
    /// Observed positive-coverage variance divided by the underlying model expectation.
    pub(crate) positive_coverage_variance_ratio: f64,
    /// Observed positions with coverage greater than or equal to `T1`.
    pub(crate) observed_initial_tail_positions: u64,
    /// Positions expected at or above `T1` under the initial model.
    pub(crate) expected_initial_tail_positions: f64,
    /// Observed positions with coverage greater than or equal to `T2`.
    pub(crate) observed_final_tail_positions: u64,
    /// Positions expected at or above `T2` under the fitted underlying model.
    pub(crate) expected_final_tail_positions: f64,
}

/// Compare observed positive counts and both fitted tails with a two-stage coverage model.
pub(crate) fn diagnose_fit(
    histogram: &CoverageCounts,
    model: TwoStageCoverageModel,
) -> Result<FitDiagnostics> {
    let (positive_positions, positive_coverage_sum) = positive_summary(histogram, None)?;
    let underlying_positive_coverage_mean = model.underlying_fit.positive_mean();
    let underlying_expected_positive_coverage_variance = underlying_positive_coverage_mean
        * (1.0
            + model.underlying_fit.dispersion * model.underlying_fit.component_mean
            + model.underlying_fit.component_mean
            - underlying_positive_coverage_mean);
    ensure!(
        underlying_expected_positive_coverage_variance.is_finite()
            && underlying_expected_positive_coverage_variance > 0.0,
        "coverage fit produced invalid positive-count variance {}",
        underlying_expected_positive_coverage_variance
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
            observed_positive_coverage_variance / underlying_expected_positive_coverage_variance,
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

    Ok(FitDiagnostics {
        observed_positive_coverage_mean,
        observed_positive_coverage_variance,
        underlying_expected_positive_coverage_variance,
        positive_coverage_variance_ratio,
        observed_initial_tail_positions,
        expected_initial_tail_positions,
        observed_final_tail_positions,
        expected_final_tail_positions,
    })
}

#[cfg(test)]
mod tests {
    include!("distribution_tests.rs");
}
