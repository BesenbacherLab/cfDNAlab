use super::*;

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
    };
    let histogram = [(0, 4), (1, 2), (3, 2)].into_iter().collect();
    let diagnostics = diagnose_fit(&histogram, model).unwrap();
    // A positive geometric variable has mean 2 and variance 2. Observed [1,1,3,3] has variance 1
    assert!((diagnostics.underlying_expected_positive_coverage_variance - 2.0).abs() < 1e-12);
    assert!((diagnostics.positive_coverage_variance_ratio - 0.5).abs() < 1e-12);
    assert!((diagnostics.expected_final_tail_positions - 0.75).abs() < 1e-12);
}
