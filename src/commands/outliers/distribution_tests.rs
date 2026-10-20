use super::*;

#[test]
fn two_stage_zinb_recovers_a_positive_dispersion_geometric_mixture() {
    // P(0)=2/3 and P(k)=2^(-k)/3 for k>0. The tiny tail above 40 is omitted
    // P(X>=19)>1e-6 and P(X>=20)<1e-6, so both thresholds are 20
    let mut histogram: CoverageCounts = (1..=40)
        .map(|coverage| (coverage, 1_u64 << (40 - coverage)))
        .collect();
    histogram.insert(0, 1_u64 << 41);
    let model = fit_two_stage(&histogram, 1e-6, OutlierModel::Zinb).unwrap();
    assert_eq!(model.initial_threshold, 20);
    assert_eq!(model.final_threshold, 20);
    assert_eq!(model.second_fit_retained_positions, 3_298_532_786_176);
    for parameters in [model.initial, model.underlying_fit] {
        assert!(
            (parameters.component_mean - 1.0).abs() < 0.01,
            "{parameters:?}"
        );
        assert!((parameters.dispersion - 1.0).abs() < 0.01, "{parameters:?}");
        assert!(
            (parameters.zero_inflation - 1.0 / 3.0).abs() < 0.01,
            "{parameters:?}"
        );
    }
}

#[test]
fn final_threshold_can_exceed_the_initial_threshold() {
    // Complete Poisson mean=.2 gives T1=3 at alpha=.001175. Correct truncation
    // gives mean=(sqrt(1.36)-.8)/1.8, whose tail at 3 exceeds alpha and at 4 passes
    let histogram = [(0, 810), (1, 180), (2, 10)].into_iter().collect();
    let model = fit_two_stage(&histogram, 0.001175, OutlierModel::Zip).unwrap();
    assert_eq!(model.initial_threshold, 3);
    assert_eq!(model.final_threshold, 4);
}

#[test]
fn geometric_component_has_hand_derived_masses_mean_and_inclusive_threshold() {
    // NB(mean=1, dispersion=1) is geometric with P(Y=k)=2^(-k-1).
    // Mixing 25% structural zeros gives P(0)=5/8 and P(X>=k)=3/4 * 2^(-k).
    let parameters = CoverageParameters {
        component_mean: 1.0,
        dispersion: 1.0,
        zero_inflation: 0.25,
    };
    for (coverage, expected) in [(0, 0.625), (1, 0.1875), (2, 0.09375), (3, 0.046875)] {
        assert!((parameters.probability_mass(coverage) - expected).abs() < 1e-12);
    }
    assert_eq!(parameters.mean(), 0.75);
    assert!((parameters.positive_mean() - 2.0).abs() < 1e-12);
    let threshold = upper_tail_threshold(parameters, 0.1).expect("geometric threshold");
    assert_eq!(threshold, 3);
    assert!((parameters.survival_at_or_above(3).unwrap() - 0.09375).abs() < 1e-12);
    assert!(parameters.survival_at_or_above(2).unwrap() > 0.1);
}

#[test]
fn zero_dispersion_has_hand_derived_zip_masses_and_threshold() {
    let parameters = CoverageParameters {
        component_mean: 0.5,
        dispersion: 0.0,
        zero_inflation: 0.2,
    };
    let poisson_zero = (-0.5_f64).exp();
    assert!((parameters.probability_mass(0) - (0.2 + 0.8 * poisson_zero)).abs() < 1e-14);
    assert!((parameters.probability_mass(1) - 0.4 * poisson_zero).abs() < 1e-14);
    assert!(
        (parameters.survival_at_or_above(1).unwrap() - 0.8 * (1.0 - poisson_zero)).abs() < 1e-14
    );
    // The mixture tails at 2 and 3 are about 0.0722 and 0.0115, respectively
    assert_eq!(upper_tail_threshold(parameters, 0.05).unwrap(), 3);
    assert_eq!(parameters.mean(), 0.4);
}

#[test]
fn positive_variance_uses_the_negative_binomial_component() {
    let parameters = CoverageParameters {
        component_mean: 1.0,
        dispersion: 1.0,
        zero_inflation: 0.25,
    };
    let model = TwoStageCoverageModel {
        family: OutlierModel::Zinb,
        initial: parameters,
        initial_threshold: 3,
        underlying_fit: parameters,
        final_threshold: 3,
        second_fit_retained_positions: 7,
        initial_optimizer: None,
        underlying_optimizer: None,
    };
    let histogram = [(0, 4), (1, 2), (3, 2)].into_iter().collect();
    let diagnostics = diagnose_fit(&histogram, model).unwrap();
    // A positive geometric variable has mean 2 and variance 2. Observed [1,1,3,3] has variance 1
    assert!((diagnostics.underlying_expected_positive_coverage_variance - 2.0).abs() < 1e-12);
    assert!((diagnostics.positive_coverage_variance_ratio - 0.5).abs() < 1e-12);
    assert!((diagnostics.expected_final_tail_positions - 0.75).abs() < 1e-12);
}

#[test]
fn diagnostics_keep_stage_specific_tails_and_use_population_positive_variance() {
    let underlying = CoverageParameters {
        component_mean: 1.0,
        dispersion: 1.0,
        zero_inflation: 0.25,
    };
    let model = TwoStageCoverageModel {
        family: OutlierModel::Zinb,
        initial: CoverageParameters {
            zero_inflation: 0.5,
            ..underlying
        },
        initial_threshold: 2,
        underlying_fit: underlying,
        final_threshold: 3,
        second_fit_retained_positions: 6,
        initial_optimizer: None,
        underlying_optimizer: None,
    };
    let histogram = [(0, 4), (1, 2), (2, 1), (3, 1), (100, 0)]
        .into_iter()
        .collect();

    let diagnostics = diagnose_fit(&histogram, model).unwrap();

    // Positive observations [1,1,2,3] have mean 7/4 and variance 15/4 - (7/4)^2 = 11/16
    assert!((diagnostics.observed_positive_coverage_mean - 7.0 / 4.0).abs() < 1e-12);
    assert!((diagnostics.observed_positive_coverage_variance - 11.0 / 16.0).abs() < 1e-12);
    assert!((diagnostics.underlying_expected_positive_coverage_variance - 2.0).abs() < 1e-12);
    assert!((diagnostics.positive_coverage_variance_ratio - 11.0 / 32.0).abs() < 1e-12);
    // Inclusive geometric tails give 8 * (1/2) * 2^-2 = 1 initially
    // The underlying final expectation is 8 * (3/4) * 2^-3 = 3/4
    assert_eq!(diagnostics.observed_initial_tail_positions, 2);
    assert_eq!(diagnostics.observed_final_tail_positions, 1);
    assert!((diagnostics.expected_initial_tail_positions - 1.0).abs() < 1e-12);
    assert!((diagnostics.expected_final_tail_positions - 0.75).abs() < 1e-12);
}

#[test]
fn zero_only_fallback_context_has_undefined_observed_positive_moments() {
    let parameters = CoverageParameters {
        component_mean: 1.0,
        dispersion: 1.0,
        zero_inflation: 0.25,
    };
    let model = TwoStageCoverageModel {
        family: OutlierModel::Zinb,
        initial: parameters,
        initial_threshold: 2,
        underlying_fit: parameters,
        final_threshold: 3,
        second_fit_retained_positions: 8,
        initial_optimizer: None,
        underlying_optimizer: None,
    };
    let histogram = [(0, 8), (100, 0)].into_iter().collect();

    let diagnostics = diagnose_fit(&histogram, model).unwrap();

    // No observed positives means undefined moments, rather than a measured variance of zero
    assert!(diagnostics.observed_positive_coverage_mean.is_nan());
    assert!(diagnostics.observed_positive_coverage_variance.is_nan());
    assert!(diagnostics.positive_coverage_variance_ratio.is_nan());
    assert!((diagnostics.underlying_expected_positive_coverage_variance - 2.0).abs() < 1e-12);
    assert_eq!(diagnostics.observed_initial_tail_positions, 0);
    assert_eq!(diagnostics.observed_final_tail_positions, 0);
    assert!((diagnostics.expected_initial_tail_positions - 1.5).abs() < 1e-12);
    assert!((diagnostics.expected_final_tail_positions - 0.75).abs() < 1e-12);
}
