use super::*;
use crate::commands::overlapping_lengths_correction::reducer::{
    LengthBinDepthStatistics, build_length_bin_edges,
};

fn assert_curves_close(left: &[f64], right: &[f64], tolerance: f64) {
    assert_eq!(left.len(), right.len());
    for (bin_index, (&left_value, &right_value)) in left.iter().zip(right).enumerate() {
        assert!(
            (left_value - right_value).abs() <= tolerance,
            "curves differ in bin {bin_index}: {left_value} versus {right_value}"
        );
    }
}

/// Evaluate the LIONHEART mixture equation independently of `mixture_distribution`.
///
/// This reference keeps the depth loop, skewed Student-t calculation, clipping factor, and final
/// normalization in one test helper. It deliberately does not call the production distribution
/// helpers, so changes to their normalization order are detectable.
fn reference_mixture_distribution(
    bin_midpoints: &[f64],
    raw_depth_frequencies: &BTreeMap<u32, u64>,
    lower_bound: f64,
    upper_bound: f64,
    parameters: MixtureFitParameters,
) -> Vec<f64> {
    let total_frequency = raw_depth_frequencies.values().sum::<u64>() as f64;
    let clipping_distribution = StudentsT::new(
        parameters.mean_fragment_length,
        BASE_SIGMA * parameters.scale_multiplier,
        STUDENT_T_DEGREES_OF_FREEDOM,
    )
    .expect("reference clipping distribution should be valid");
    let clipping_factor =
        clipping_distribution.cdf(upper_bound) - clipping_distribution.cdf(lower_bound);
    let root_two = 2.0_f64.sqrt();
    let mut mixture = vec![0.0; bin_midpoints.len()];

    for (&depth, &frequency) in raw_depth_frequencies {
        let depth_scale = BASE_SIGMA / (depth as f64).sqrt() * parameters.scale_multiplier;
        let distribution = StudentsT::new(
            parameters.mean_fragment_length,
            depth_scale,
            STUDENT_T_DEGREES_OF_FREEDOM,
        )
        .expect("reference depth distribution should be valid");
        let mut component = bin_midpoints
            .iter()
            .map(|midpoint| {
                distribution.pdf(*midpoint)
                    * (1.0
                        + erf(
                            parameters.skewness
                                * (*midpoint - parameters.mean_fragment_length)
                                / (depth_scale * root_two),
                        ))
            })
            .collect::<Vec<_>>();
        let component_area = bin_midpoints
            .windows(2)
            .zip(component.windows(2))
            .map(|(midpoints, values)| {
                (midpoints[1] - midpoints[0]) * (values[0] + values[1]) / 2.0
            })
            .sum::<f64>();
        for value in &mut component {
            *value /= component_area;
        }
        let depth_probability = frequency as f64 / total_frequency;
        for (mixture_value, component_value) in mixture.iter_mut().zip(component) {
            *mixture_value +=
                component_value / (clipping_factor * depth_scale) * depth_probability;
        }
    }

    let mixture_mean = mixture.iter().sum::<f64>() / mixture.len() as f64;
    mixture
        .into_iter()
        .map(|value| value / mixture_mean)
        .collect()
}

fn synthetic_statistics() -> Result<(OverlappingLengthStatistics, MixtureFitParameters)> {
    let length_bin_edges = build_length_bin_edges(100, 220, 3)?;
    let bin_midpoints = length_bin_edges
        .windows(2)
        .map(|edges| (edges[0] + edges[1]) / 2.0)
        .collect::<Vec<_>>();
    let bases_per_bin = 100_000_u64;
    let eligible_covered_bases = bases_per_bin * bin_midpoints.len() as u64;
    let raw_depth_frequencies = BTreeMap::from([
        (1_u32, 400_000_u64),
        (2, 1_200_000),
        (3, 1_600_000),
        (4, 800_000),
    ]);
    assert_eq!(
        raw_depth_frequencies.values().sum::<u64>(),
        eligible_covered_bases
    );
    let generating_parameters = MixtureFitParameters {
        scale_multiplier: 8.0,
        skewness: 0.0,
        mean_fragment_length: 160.0,
    };
    let observed_bias = mixture_distribution(
        &bin_midpoints,
        &raw_depth_frequencies,
        100.0,
        220.0,
        generating_parameters,
    )?;
    let observed_signal_sums = observed_bias
        .iter()
        .map(|value| value * bases_per_bin as f64)
        .collect::<Vec<_>>();
    let length_bin_depth_statistics = bin_midpoints
        .iter()
        .map(|&midpoint| {
            BTreeMap::from([
                (
                    1,
                    LengthBinDepthStatistics {
                        position_count: 10_000,
                        average_length_sum: midpoint * 10_000.0,
                    },
                ),
                (
                    2,
                    LengthBinDepthStatistics {
                        position_count: 30_000,
                        average_length_sum: midpoint * 30_000.0,
                    },
                ),
                (
                    3,
                    LengthBinDepthStatistics {
                        position_count: 40_000,
                        average_length_sum: midpoint * 40_000.0,
                    },
                ),
                (
                    4,
                    LengthBinDepthStatistics {
                        position_count: 20_000,
                        average_length_sum: midpoint * 20_000.0,
                    },
                ),
            ])
        })
        .collect();

    Ok((
        OverlappingLengthStatistics {
            length_bin_edges,
            length_bin_base_counts: vec![bases_per_bin; bin_midpoints.len()],
            observed_signal_sums,
            raw_depth_frequencies,
            length_bin_depth_statistics,
            eligible_covered_bases,
        },
        generating_parameters,
    ))
}

#[test]
fn gaussian_smoothing_matches_zero_padded_same_convolution() {
    let smoothed = smooth_with_five_point_gaussian_kernel(&[1.0, 1.0, 1.0, 1.0, 1.0]);
    assert!(smoothed[0] < smoothed[2]);
    assert!((smoothed[2] - 1.0).abs() < 1.0e-12);
}

#[test]
fn gaussian_smoothing_uses_five_point_standard_deviation_one_kernel() {
    let smoothed = smooth_with_five_point_gaussian_kernel(&[0.0, 0.0, 1.0, 0.0, 0.0]);

    assert!((smoothed.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
    assert!((smoothed[0] - smoothed[4]).abs() < 1.0e-12);
    assert!((smoothed[1] - smoothed[3]).abs() < 1.0e-12);
    assert!((smoothed[1] / smoothed[2] - (-0.5_f64).exp()).abs() < 1.0e-12);
    assert!((smoothed[0] / smoothed[2] - (-2.0_f64).exp()).abs() < 1.0e-12);
}

#[test]
fn mixture_matches_lionheart_equation_for_fixed_depth_frequencies() -> Result<()> {
    let bin_midpoints = build_length_bin_edges(100, 220, 3)?
        .windows(2)
        .map(|edges| (edges[0] + edges[1]) / 2.0)
        .collect::<Vec<_>>();
    let raw_depth_frequencies = BTreeMap::from([(1_u32, 300_u64), (2, 500), (4, 200)]);
    let parameters = MixtureFitParameters {
        scale_multiplier: 2.5,
        skewness: -0.3,
        mean_fragment_length: 164.0,
    };

    let fitted = mixture_distribution(
        &bin_midpoints,
        &raw_depth_frequencies,
        100.0,
        220.0,
        parameters,
    )?;
    let reference = reference_mixture_distribution(
        &bin_midpoints,
        &raw_depth_frequencies,
        100.0,
        220.0,
        parameters,
    );

    assert_curves_close(&fitted, &reference, 1.0e-12);
    assert!((fitted.iter().sum::<f64>() / fitted.len() as f64 - 1.0).abs() < 1.0e-12);
    Ok(())
}

#[test]
fn mixture_is_invariant_to_common_depth_frequency_scaling() -> Result<()> {
    let bin_midpoints = build_length_bin_edges(100, 220, 3)?
        .windows(2)
        .map(|edges| (edges[0] + edges[1]) / 2.0)
        .collect::<Vec<_>>();
    let small_frequencies = BTreeMap::from([(1_u32, 3_u64), (2, 5), (4, 2)]);
    let large_frequencies = BTreeMap::from([
        (1_u32, 3_000_000_u64),
        (2, 5_000_000),
        (4, 2_000_000),
    ]);
    let parameters = MixtureFitParameters {
        scale_multiplier: 2.5,
        skewness: -0.3,
        mean_fragment_length: 164.0,
    };

    let small = mixture_distribution(
        &bin_midpoints,
        &small_frequencies,
        100.0,
        220.0,
        parameters,
    )?;
    let large = mixture_distribution(
        &bin_midpoints,
        &large_frequencies,
        100.0,
        220.0,
        parameters,
    )?;

    assert_curves_close(&small, &large, 1.0e-12);
    Ok(())
}

#[test]
fn three_division_factors_combine_multiplicatively() {
    let noise = [2.0, 0.5];
    let skew = [0.5, 2.0];
    let mean = [4.0, 0.25];
    let combined = noise
        .iter()
        .zip(skew)
        .zip(mean)
        .map(|((noise, skew), mean)| 1.0 / (noise * skew * mean))
        .collect::<Vec<_>>();
    assert_eq!(combined, vec![0.25, 4.0]);
}

#[test]
fn initial_location_is_weighted_by_covered_positions_and_uses_actual_average_lengths() -> Result<()>
{
    let mut statistics = OverlappingLengthStatistics::new(vec![100.0, 110.0, 120.0])?;
    statistics.add_position(101.25, 0.1, 1)?;
    statistics.add_position(114.0, 10.0, 2)?;
    statistics.add_position(115.5, 10.0, 2)?;

    // Signal magnitude and the unweighted midpoint mean are irrelevant. LIONHEART initializes
    // from positional average lengths with positive raw coverage: (101.25 + 114 + 115.5) / 3.
    let initial_location = statistics.mean_average_overlapping_length()?;

    assert!((initial_location - 110.25).abs() < 1.0e-12);
    Ok(())
}

#[test]
fn refit_sampling_filters_corrected_depth_and_uses_numpy_ties_to_even_rounding() -> Result<()> {
    let mut statistics = OverlappingLengthStatistics::new(vec![100.0, 110.0, 120.0])?;
    statistics.add_position(101.0, 1.0, 1)?;
    statistics.add_position(102.0, 1.0, 1)?;
    statistics.add_position(103.0, 3.0, 3)?;
    statistics.add_position(117.0, 5.0, 5)?;

    let sampling = refit_sampling_statistics(&statistics, &[2.0, 2.0], &[1.0, 1.0])?;

    // Depth 1 becomes exactly 0.5 and is excluded. NumPy rounds both 1.5 and 2.5 to the even 2;
    // ordinary half-away-from-zero rounding would incorrectly place the latter at depth 3.
    assert_eq!(sampling.depth_frequencies, BTreeMap::from([(2, 2)]));
    assert!((sampling.mean_average_length - 110.0).abs() < 1.0e-12);
    Ok(())
}

#[test]
fn complete_two_fit_model_recovers_synthetic_mixture_and_target_correction() -> Result<()> {
    // Four million bases are represented by sufficient statistics rather than individual values
    let (statistics, generating_parameters) = synthetic_statistics()?;

    let model = fit_overlapping_length_model(&statistics)?;

    let fitted_objective = objective(
        [
            model.initial_fit.scale_multiplier,
            model.initial_fit.skewness,
            model.initial_fit.mean_fragment_length,
        ],
        &model.bin_midpoints,
        &statistics.raw_depth_frequencies,
        &model.observed_bias,
        100.0,
        220.0,
        FIRST_SKEWNESS_PENALTY,
    )?;
    assert!(fitted_objective < 1.0e-8);
    assert!(
        (model.initial_fit.scale_multiplier - generating_parameters.scale_multiplier).abs() < 0.1
    );
    assert!((model.initial_fit.skewness - generating_parameters.skewness).abs() < 0.01);
    assert!(
        (model.initial_fit.mean_fragment_length - generating_parameters.mean_fragment_length).abs()
            < 0.1
    );

    // Both persisted fitted curves are the five-point Gaussian smoothing of their fitted mixtures
    let first_fitted_unsmoothed = mixture_distribution(
        &model.bin_midpoints,
        &statistics.raw_depth_frequencies,
        100.0,
        220.0,
        model.initial_fit,
    )?;
    assert_curves_close(
        &model.first_fitted_bias,
        &smooth_with_five_point_gaussian_kernel(&first_fitted_unsmoothed),
        1.0e-12,
    );
    let expected_noise_division_factors = scale_to_mean_one(&divide_elementwise(
        &model.observed_bias,
        &model.first_fitted_bias,
        "test noise correction",
    )?)?;
    assert_curves_close(
        &model.noise_division_factors,
        &expected_noise_division_factors,
        1.0e-12,
    );
    let refit_sampling = refit_sampling_statistics(
        &statistics,
        &model.noise_division_factors,
        &model.skew_division_factors,
    )?;
    let second_fitted_unsmoothed = mixture_distribution(
        &model.bin_midpoints,
        &refit_sampling.depth_frequencies,
        100.0,
        220.0,
        model.refit,
    )?;
    assert_curves_close(
        &model.second_fitted_bias,
        &smooth_with_five_point_gaussian_kernel(&second_fitted_unsmoothed),
        1.0e-12,
    );
    let refit_start_objective = objective(
        [8.0, -0.5, refit_sampling.mean_average_length],
        &model.bin_midpoints,
        &refit_sampling.depth_frequencies,
        &model.first_corrected_bias,
        100.0,
        220.0,
        REFIT_SKEWNESS_PENALTY,
    )?;
    let refit_objective = objective(
        [
            model.refit.scale_multiplier,
            model.refit.skewness,
            model.refit.mean_fragment_length,
        ],
        &model.bin_midpoints,
        &refit_sampling.depth_frequencies,
        &model.first_corrected_bias,
        100.0,
        220.0,
        REFIT_SKEWNESS_PENALTY,
    )?;
    assert!(refit_objective <= refit_start_objective);

    for curve in [
        &model.observed_bias,
        &model.first_fitted_bias,
        &model.first_corrected_bias,
        &model.second_fitted_bias,
        &model.target_bias,
        &model.noise_division_factors,
        &model.skew_division_factors,
        &model.mean_shift_division_factors,
        &model.combined_weights,
    ] {
        assert_eq!(curve.len(), model.bin_midpoints.len());
        assert!(curve.iter().all(|value| value.is_finite() && *value > 0.0));
    }
    for normalized_curve in [
        &model.observed_bias,
        &model.first_corrected_bias,
        &model.target_bias,
        &model.noise_division_factors,
        &model.skew_division_factors,
        &model.mean_shift_division_factors,
    ] {
        assert!(
            (normalized_curve.iter().sum::<f64>() / normalized_curve.len() as f64 - 1.0).abs()
                < 1.0e-10
        );
    }

    let corrected_observed = model
        .observed_bias
        .iter()
        .zip(&model.combined_weights)
        .map(|(observed, weight)| observed * weight)
        .collect::<Vec<_>>();
    let corrected_observed = scale_to_mean_one(&corrected_observed)?;
    assert_curves_close(&corrected_observed, &model.target_bias, 1.0e-10);
    for bin_index in 0..model.bin_midpoints.len() {
        let expected_weight = 1.0
            / (model.noise_division_factors[bin_index]
                * model.skew_division_factors[bin_index]
                * model.mean_shift_division_factors[bin_index]);
        assert!((model.combined_weights[bin_index] - expected_weight).abs() < 1.0e-12);
    }
    Ok(())
}
