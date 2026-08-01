//! LIONHEART-compatible mixture distributions, correction sequence, and BFGS optimizer.

use std::collections::BTreeMap;

use anyhow::{Context, Result, ensure};
use statrs::distribution::{Continuous, ContinuousCDF, StudentsT};
use statrs::function::erf::erf;

use super::package::MixtureFitParameters;
use super::reducer::OverlappingLengthStatistics;
use super::scipy_bfgs::minimize_bfgs;

/// Empirical single-fragment spread used by LIONHEART before depth scaling.
pub(crate) const BASE_SIGMA: f64 = 8.026_649_608_460_776;
/// Degrees of freedom of every Student-t component.
pub(crate) const STUDENT_T_DEGREES_OF_FREEDOM: f64 = 5.0;
/// Mean fragment length of the final symmetric target distribution.
pub(crate) const TARGET_MEAN_FRAGMENT_LENGTH: f64 = 166.0;
/// Quadratic skewness penalty used by the initial noise-and-skew fit.
const FIRST_SKEWNESS_PENALTY: f64 = 0.005;
/// Stronger quadratic skewness penalty used after noise and skew correction.
const REFIT_SKEWNESS_PENALTY: f64 = 0.1;

/// Complete in-memory result of LIONHEART's two-fit correction sequence.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct OverlappingLengthModel {
    /// Midpoint of every configured average-length bin.
    pub bin_midpoints: Vec<f64>,
    /// Mean observed signal by length bin, normalized to mean one.
    pub observed_bias: Vec<f64>,
    /// Smoothed mixture curve from the first fit.
    pub first_fitted_bias: Vec<f64>,
    /// Observed curve after noise and linear-skew division.
    pub first_corrected_bias: Vec<f64>,
    /// Smoothed mixture curve fitted to `first_corrected_bias`.
    pub second_fitted_bias: Vec<f64>,
    /// Symmetric refit-scale mixture centered at 166 bp.
    pub target_bias: Vec<f64>,
    /// First-stage observed-to-fitted noise division factors.
    pub noise_division_factors: Vec<f64>,
    /// Normalized linear division factors derived from the first fitted skewness.
    pub skew_division_factors: Vec<f64>,
    /// Division factors shifting the normalized curve to the target mean.
    pub mean_shift_division_factors: Vec<f64>,
    /// Final lookup weights equal to the reciprocal product of all three division factors.
    pub combined_weights: Vec<f64>,
    /// Optimized parameters from the initial fit.
    pub initial_fit: MixtureFitParameters,
    /// Optimized parameters from the fit after noise and skew correction.
    pub refit: MixtureFitParameters,
}

/// Fit the two LIONHEART mixture models from a single sweep's additive statistics.
///
/// The fitting sequence intentionally preserves LIONHEART's normalization order and constants:
///
/// 1. Normalize observed per-bin signal to mean one and fit the raw mixture.
/// 2. Smooth the fitted curve and derive the noise and linear-skew division factors.
/// 3. Normalize the stored per-bin means by those two factors and fit the mixture again.
/// 4. Construct a symmetric target at 166 bp using the refitted scale.
/// 5. Derive the mean-shift division factor and multiply the reciprocal factors.
///
/// Both mixture fits use the same raw integer depth-frequency table. Corrected or scaled coverage
/// never controls the `1 / sqrt(depth)` sampling spread.
pub(crate) fn fit_overlapping_length_model(
    statistics: &OverlappingLengthStatistics,
) -> Result<OverlappingLengthModel> {
    ensure!(
        statistics.eligible_covered_bases > 0,
        "no eligible covered bases were available for overlapping fragment length fitting"
    );
    // Construct the observed curve retained from the single genomic sweep
    let bin_midpoints = statistics
        .length_bin_edges
        .windows(2)
        .map(|pair| (pair[0] + pair[1]) / 2.0)
        .collect::<Vec<_>>();
    let observed_bias = scale_to_mean_one(&statistics.observed_means()?)?;
    let start_mean =
        mean_midpoint_of_bins_above_half_normalized_signal(&bin_midpoints, &observed_bias)?;
    let lower_bound = *statistics
        .length_bin_edges
        .first()
        .context("missing first overlapping fragment length bin edge")?;
    let upper_bound = *statistics
        .length_bin_edges
        .last()
        .context("missing final overlapping fragment length bin edge")?;

    // Fit raw observed bias with LIONHEART's initial parameters and light skew penalty
    let initial_fit = optimize_mixture(
        &bin_midpoints,
        &statistics.raw_depth_frequencies,
        &observed_bias,
        lower_bound,
        upper_bound,
        MixtureFitParameters {
            scale_multiplier: 8.0,
            skewness: -0.5,
            mean_fragment_length: start_mean,
        },
        FIRST_SKEWNESS_PENALTY,
    )?;
    // Smooth only after optimization, matching the Python implementation's fitting sequence
    let first_fitted_unsmoothed = mixture_distribution(
        &bin_midpoints,
        &statistics.raw_depth_frequencies,
        lower_bound,
        upper_bound,
        initial_fit,
    )?;
    let first_fitted_bias = smooth_with_five_point_gaussian_kernel(&first_fitted_unsmoothed);

    // Separate short-scale residual noise from the parametric fitted distribution
    let noise_division_factors = scale_to_mean_one(&divide_elementwise(
        &observed_bias,
        &first_fitted_bias,
        "noise correction",
    )?)?;
    let mean_midpoint = bin_midpoints.iter().sum::<f64>() / bin_midpoints.len() as f64;
    // Preserve LIONHEART's normalized linear skew correction exactly
    let raw_skew_division_factors = bin_midpoints
        .iter()
        .map(|midpoint| {
            midpoint * initial_fit.skewness + mean_midpoint * (1.0 - initial_fit.skewness) + 1.0
        })
        .collect::<Vec<_>>();
    let skew_division_factors = scale_to_mean_one(&raw_skew_division_factors)?;
    // Bin-wise constant division factors make this identical to normalizing every position and
    // re-averaging
    let first_corrected_bias = scale_to_mean_one(&divide_elementwise(
        &divide_elementwise(&observed_bias, &noise_division_factors, "noise correction")?,
        &skew_division_factors,
        "skew correction",
    )?)?;

    // Refit the intermediate curve with the stronger skewness penalty
    let refit_start_mean =
        mean_midpoint_of_bins_above_half_normalized_signal(&bin_midpoints, &first_corrected_bias)?;
    let refit = optimize_mixture(
        &bin_midpoints,
        &statistics.raw_depth_frequencies,
        &first_corrected_bias,
        lower_bound,
        upper_bound,
        MixtureFitParameters {
            scale_multiplier: 8.0,
            skewness: -0.5,
            mean_fragment_length: refit_start_mean,
        },
        REFIT_SKEWNESS_PENALTY,
    )?;
    let second_fitted_bias = smooth_with_five_point_gaussian_kernel(&mixture_distribution(
        &bin_midpoints,
        &statistics.raw_depth_frequencies,
        lower_bound,
        upper_bound,
        refit,
    )?);
    // Keep the refitted scale but force the target to mean 166 and zero skew
    let target_bias = scale_to_mean_one(&mixture_distribution(
        &bin_midpoints,
        &statistics.raw_depth_frequencies,
        lower_bound,
        upper_bound,
        MixtureFitParameters {
            scale_multiplier: refit.scale_multiplier,
            skewness: 0.0,
            mean_fragment_length: TARGET_MEAN_FRAGMENT_LENGTH,
        },
    )?)?;
    // The final division shifts the intermediate curve onto the symmetric target
    let mean_shift_division_factors = scale_to_mean_one(&divide_elementwise(
        &first_corrected_bias,
        &target_bias,
        "mean correction",
    )?)?;
    // fcoverage applies one scalar lookup rather than repeating three divisions
    let combined_weights = noise_division_factors
        .iter()
        .zip(&skew_division_factors)
        .zip(&mean_shift_division_factors)
        .enumerate()
        .map(|(bin_index, ((noise, skew), mean))| {
            let combined_division_factor = noise * skew * mean;
            ensure!(
                combined_division_factor.is_finite() && combined_division_factor > 0.0,
                "combined overlapping fragment length division factor is invalid in bin {}: {}",
                bin_index,
                combined_division_factor
            );
            Ok(1.0 / combined_division_factor)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(OverlappingLengthModel {
        bin_midpoints,
        observed_bias,
        first_fitted_bias,
        first_corrected_bias,
        second_fitted_bias,
        target_bias,
        noise_division_factors,
        skew_division_factors,
        mean_shift_division_factors,
        combined_weights,
        initial_fit,
        refit,
    })
}

/// Build the depth-weighted skewed Student-t mixture for a parameter set.
///
/// Each positive raw depth contributes a component whose spread is
/// `BASE_SIGMA / sqrt(depth) * scale_multiplier`. Component weights are the genome-wide raw-depth
/// frequencies. The clipping scale follows LIONHEART and is calculated once with the unscaled
/// base sigma times the fitted multiplier.
fn mixture_distribution(
    bin_midpoints: &[f64],
    raw_depth_frequencies: &BTreeMap<u32, u64>,
    lower_bound: f64,
    upper_bound: f64,
    parameters: MixtureFitParameters,
) -> Result<Vec<f64>> {
    ensure!(
        parameters.scale_multiplier.is_finite() && parameters.scale_multiplier > 0.0,
        "mixture scale multiplier must be finite and positive"
    );
    ensure!(
        parameters.skewness.is_finite() && parameters.mean_fragment_length.is_finite(),
        "mixture skewness and mean must be finite"
    );
    let total_depth_observations = raw_depth_frequencies.values().sum::<u64>();
    ensure!(
        total_depth_observations > 0,
        "raw depth frequency table is empty"
    );

    // Calculate the probability mass inside the configured fragment length interval
    let clipping_distribution = StudentsT::new(
        parameters.mean_fragment_length,
        BASE_SIGMA * parameters.scale_multiplier,
        STUDENT_T_DEGREES_OF_FREEDOM,
    )
    .context("create clipping Student-t distribution")?;
    let t_scale_factor =
        clipping_distribution.cdf(upper_bound) - clipping_distribution.cdf(lower_bound);
    ensure!(
        t_scale_factor.is_finite() && t_scale_factor > 0.0,
        "Student-t clipping scale is invalid"
    );

    // Accumulate only observed positive depths, weighted by their genomic frequencies
    let mut composite = vec![0.0; bin_midpoints.len()];
    for (&depth, &count) in raw_depth_frequencies {
        if depth == 0 || count == 0 {
            continue;
        }
        let probability = count as f64 / total_depth_observations as f64;
        let depth_sigma = BASE_SIGMA / (depth as f64).sqrt() * parameters.scale_multiplier;
        let distribution = skewed_t_pdf(
            bin_midpoints,
            parameters.mean_fragment_length,
            depth_sigma,
            parameters.skewness,
        )?;
        for (target, value) in composite.iter_mut().zip(distribution) {
            *target += value / (t_scale_factor * depth_sigma) * probability;
        }
    }
    scale_to_mean_one(&composite)
}

/// Evaluate LIONHEART's skewed Student-t density and normalize its sampled trapezoid area.
///
/// This is the same pragmatic skew function used by LIONHEART. The Student-t PDF is multiplied by
/// `1 + erf(skewness * (x - location) / (scale * sqrt(2)))`.
fn skewed_t_pdf(x: &[f64], location: f64, scale: f64, skewness: f64) -> Result<Vec<f64>> {
    let distribution = StudentsT::new(location, scale, STUDENT_T_DEGREES_OF_FREEDOM)
        .context("create depth Student-t distribution")?;
    let root_two = 2.0_f64.sqrt();
    let mut values = x
        .iter()
        .map(|value| {
            distribution.pdf(*value)
                * (1.0 + erf(skewness * (*value - location) / (scale * root_two)))
        })
        .collect::<Vec<_>>();
    let area = trapezoid_integral(x, &values)?;
    ensure!(
        area.is_finite() && area > 0.0,
        "skewed Student-t density has invalid area"
    );
    for value in &mut values {
        *value /= area;
    }
    Ok(values)
}

/// Integrate sampled values over possibly nonuniform bin midpoints with the trapezoid rule.
fn trapezoid_integral(x: &[f64], values: &[f64]) -> Result<f64> {
    ensure!(
        x.len() == values.len() && x.len() >= 2,
        "trapezoid integration requires matching arrays with at least two values"
    );
    Ok(x.windows(2)
        .zip(values.windows(2))
        .map(|(x_pair, value_pair)| (x_pair[1] - x_pair[0]) * (value_pair[0] + value_pair[1]) / 2.0)
        .sum())
}

/// Apply LIONHEART's five-point Gaussian kernel with standard deviation one.
///
/// Values beyond the ends are treated as zero, matching NumPy's `convolve(..., mode="same")` for
/// a five-value kernel when the fitted curve is at least as long as the kernel.
fn smooth_with_five_point_gaussian_kernel(values: &[f64]) -> Vec<f64> {
    let mut kernel =
        [-2.0_f64, -1.0, 0.0, 1.0, 2.0].map(|position| (-0.5 * position.powi(2)).exp());
    let kernel_sum = kernel.iter().sum::<f64>();
    for weight in &mut kernel {
        *weight /= kernel_sum;
    }
    let mut output = vec![0.0; values.len()];
    for (output_index, output_value) in output.iter_mut().enumerate() {
        for (kernel_index, weight) in kernel.iter().enumerate() {
            let input_index = output_index as isize + kernel_index as isize - 2;
            if input_index >= 0 && (input_index as usize) < values.len() {
                *output_value += values[input_index as usize] * weight;
            }
        }
    }
    output
}

/// Scale a finite vector so the arithmetic mean of its values is one.
///
/// This changes only the vertical scale of a sampled curve. It does not calculate or change the
/// fitted fragment length location parameter, which is also called a mean in the model.
fn scale_to_mean_one(values: &[f64]) -> Result<Vec<f64>> {
    ensure!(!values.is_empty(), "cannot scale an empty distribution");
    let arithmetic_mean = values.iter().sum::<f64>() / values.len() as f64;
    ensure!(
        arithmetic_mean.is_finite() && arithmetic_mean != 0.0,
        "distribution values have invalid arithmetic mean {}",
        arithmetic_mean
    );
    let scaled_values = values
        .iter()
        .map(|value| value / arithmetic_mean)
        .collect::<Vec<_>>();
    ensure!(
        scaled_values.iter().all(|value| value.is_finite()),
        "scaling distribution values to arithmetic mean one produced non-finite values"
    );
    Ok(scaled_values)
}

/// Divide two equal-length curves with contextual errors for invalid division factors.
fn divide_elementwise(numerator: &[f64], denominator: &[f64], label: &str) -> Result<Vec<f64>> {
    ensure!(
        numerator.len() == denominator.len(),
        "{} arrays have different lengths",
        label
    );
    numerator
        .iter()
        .zip(denominator)
        .enumerate()
        .map(|(bin_index, (numerator, denominator))| {
            ensure!(
                denominator.is_finite() && *denominator != 0.0,
                "{} division factor is invalid in bin {}",
                label,
                bin_index
            );
            Ok(numerator / denominator)
        })
        .collect()
}

/// Choose an initial fitted location from bins above half the curve's arithmetic mean.
///
/// The curve has already been scaled so its arithmetic mean is one. This helper therefore averages
/// the fragment length bin midpoints whose relative signal is greater than 0.5. The numerical
/// threshold originates in LIONHEART, but LIONHEART applies it to per-position coverage before
/// averaging the corresponding per-position overlap lengths. Applying it to the retained binned
/// curve is the single-sweep approximation specified for this implementation. It is not
/// mathematically identical to LIONHEART's positional initialization.
fn mean_midpoint_of_bins_above_half_normalized_signal(
    midpoints: &[f64],
    normalized_signal: &[f64],
) -> Result<f64> {
    ensure!(
        midpoints.len() == normalized_signal.len(),
        "fragment length midpoints and normalized signal have different lengths"
    );
    let mut midpoint_sum = 0.0;
    let mut selected_bin_count = 0_usize;
    for (&midpoint, &value) in midpoints.iter().zip(normalized_signal) {
        if value > 0.5 {
            midpoint_sum += midpoint;
            selected_bin_count += 1;
        }
    }
    ensure!(
        selected_bin_count > 0,
        "no overlapping fragment length bins have normalized coverage above 0.5"
    );
    Ok(midpoint_sum / selected_bin_count as f64)
}

/// Calculate mean squared curve error plus the fit-specific quadratic skewness penalty.
fn objective(
    parameters: [f64; 3],
    midpoints: &[f64],
    depths: &BTreeMap<u32, u64>,
    observed: &[f64],
    lower: f64,
    upper: f64,
    skewness_penalty: f64,
) -> Result<f64> {
    let fitted = mixture_distribution(
        midpoints,
        depths,
        lower,
        upper,
        MixtureFitParameters {
            scale_multiplier: parameters[0],
            skewness: parameters[1],
            mean_fragment_length: parameters[2],
        },
    )?;
    let mse = fitted
        .iter()
        .zip(observed)
        .map(|(fit, observation)| (fit - observation).powi(2))
        .sum::<f64>()
        / fitted.len() as f64;
    Ok(mse + skewness_penalty * parameters[1].powi(2))
}

#[allow(clippy::too_many_arguments)]
/// Optimize mixture scale, skewness, and mean from a LIONHEART-compatible starting point.
fn optimize_mixture(
    midpoints: &[f64],
    depths: &BTreeMap<u32, u64>,
    observed: &[f64],
    lower: f64,
    upper: f64,
    initial: MixtureFitParameters,
    skewness_penalty: f64,
) -> Result<MixtureFitParameters> {
    let evaluate = |point: [f64; 3]| {
        objective(
            point,
            midpoints,
            depths,
            observed,
            lower,
            upper,
            skewness_penalty,
        )
    };
    let point = minimize_bfgs(
        [
            initial.scale_multiplier,
            initial.skewness,
            initial.mean_fragment_length,
        ],
        &evaluate,
    )?;
    Ok(MixtureFitParameters {
        scale_multiplier: point[0],
        skewness: point[1],
        mean_fragment_length: point[2],
    })
}

#[cfg(test)]
mod tests {
    include!("model_tests.rs");
}
