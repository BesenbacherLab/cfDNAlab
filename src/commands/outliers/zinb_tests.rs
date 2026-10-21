use super::*;

#[test]
fn two_retained_levels_cannot_identify_a_free_zero_mixture() {
    let histogram = [(0, 90), (1, 10)].into_iter().collect();
    assert!(
        fit_zinb(&histogram, Some(2), 1e-6)
            .unwrap_err()
            .to_string()
            .contains("unidentified")
    );
}

#[test]
fn three_level_ridge_is_rejected_but_its_poisson_boundary_is_retained() {
    // Counts 12,2,1 have the same retained proportions under both
    // NB(mean=1,dispersion=1,pi=1/2) and ZIP(mean=1,pi=5/(e+5))
    let histogram = [(0, 12), (1, 2), (2, 1)].into_iter().collect();
    assert!(
        fit_zinb(&histogram, Some(3), 1e-6)
            .unwrap_err()
            .to_string()
            .contains("ridge")
    );
    // At ratio n2/n1=1/2, the Poisson conditional zero fraction is 2/5
    // Equality and a roundoff-sized excess are allowed, a clear excess is rejected
    for zeros in [4_000_000_000_000_u64, 4_000_000_000_001] {
        let counts = [(0, zeros), (1, 4_000_000_000_000), (2, 2_000_000_000_000)]
            .into_iter()
            .collect();
        let boundary = fit_right_truncated_zip(&counts, 3).unwrap().0;
        check_truncated_identification(&counts, Some(3), boundary.into()).unwrap();
    }
    let counts = [
        (0, 4_000_000_001_000),
        (1, 4_000_000_000_000),
        (2, 2_000_000_000_000),
    ]
    .into_iter()
    .collect();
    let boundary = fit_right_truncated_zip(&counts, 3).unwrap().0;
    assert!(check_truncated_identification(&counts, Some(3), boundary.into()).is_err());
    let sparse = [(0, 150), (1, 150)].into_iter().collect();
    let parameters = fit_zinb(&sparse, Some(3), 1e-6).unwrap().parameters;
    assert_eq!(parameters.dispersion, 0.0);
    assert_eq!(parameters.zero_inflation, 0.0);
}

#[test]
fn stable_log_mass_matches_geometric_and_integer_shape_probabilities() {
    for coverage in 0..20 {
        assert!(
            (component_log_probability(1.0, 1.0, coverage)
                + (coverage as f64 + 1.0) * 2.0_f64.ln())
            .abs()
                < 1e-11
        );
    }
    // At shape=16 and mean=16, p=1/2. P(0)=2^-16 and P(1)=8*2^-16
    assert!((component_log_probability(16.0, 1.0 / 16.0, 0).exp() - 1.0 / 65_536.0).abs() < 1e-14);
    assert!((component_log_probability(16.0, 1.0 / 16.0, 1).exp() - 8.0 / 65_536.0).abs() < 1e-14);
    for coverage in [0, 1, 17, 30, 200] {
        assert!(
            (component_log_probability(17.0, 1e-10, coverage)
                - poisson_log_probability(17.0, coverage))
            .abs()
                < 2e-6
        );
    }
}

#[test]
fn truncated_profile_accounts_for_retention_when_solving_zero_inflation() {
    // At mean=1, dispersion=1, pi=1/2, masses at 0,1,2 are 3/4,1/8,1/16.
    // Retention below 3 gives counts 12,2,1 out of 15 and recovers underlying pi=1/2
    let histogram = [(0, 12), (1, 2), (2, 1), (100, 20)].into_iter().collect();
    let problem = ProfiledHistogram::new(&histogram, Some(3)).unwrap();
    let (cost, zero_inflation) = problem.evaluate(1.0, 1.0).unwrap();
    let expected =
        -(12.0 * (12.0_f64 / 15.0).ln() + 2.0 * (2.0_f64 / 15.0).ln() + (1.0_f64 / 15.0).ln())
            / 15.0;
    assert!((zero_inflation - 0.5).abs() < 1e-12);
    assert!((cost - expected).abs() < 1e-12);
}

#[test]
fn truncated_geometric_histogram_recovers_its_underlying_component() {
    // Counts 2^(15-k), k=0..15, are exactly proportional to a geometric(1/2)
    // retained below 16. Correct truncation estimates underlying mean=1 and dispersion=1
    let mut histogram: CoverageCounts = (0..16)
        .map(|coverage| (coverage, 1_u64 << (15 - coverage)))
        .collect();
    histogram.insert(100, 100);
    let parameters = fit_zinb(&histogram, Some(16), 1e-6).unwrap().parameters;
    assert!(
        (parameters.component_mean - 1.0).abs() < 0.01,
        "{parameters:?}"
    );
    assert!((parameters.dispersion - 1.0).abs() < 0.01, "{parameters:?}");
    assert!(parameters.zero_inflation < 0.01);
}

#[test]
fn underdispersed_counts_select_the_explicit_poisson_boundary() {
    // Zeros and ones alone have less spread than a Poisson. Extra NB dispersion cannot help
    let histogram = [(0, 100), (1, 100)].into_iter().collect();
    let parameters = fit_zinb(&histogram, None, 1e-6).unwrap().parameters;
    assert_eq!(parameters.dispersion, 0.0);
    assert_eq!(parameters.component_mean, 0.5);
    assert_eq!(parameters.zero_inflation, 0.0);
}

#[test]
fn truncated_profile_constrains_zero_inflation_when_observed_zeros_are_too_few() {
    let histogram = [(0, 1), (1, 1), (2, 1)].into_iter().collect();
    let problem = ProfiledHistogram::new(&histogram, Some(3)).unwrap();

    let (cost, zero_inflation) = problem.evaluate(1.0, 1.0).unwrap();

    // Geometric masses below 3 normalize to [4/7,2/7,1/7]
    // Observed zeros are only 1/3, so adding structural zeros cannot improve this component
    let expected = -((4.0_f64 / 7.0).ln() + (2.0_f64 / 7.0).ln() + (1.0_f64 / 7.0).ln()) / 3.0;
    assert_eq!(zero_inflation, 0.0);
    assert!((cost - expected).abs() < 1e-12);
}

#[test]
fn complete_profile_without_observed_zeros_retains_the_component_zero_probability() {
    let histogram = [(0, 0), (1, 2), (2, 1)].into_iter().collect();
    let problem = ProfiledHistogram::new(&histogram, None).unwrap();

    let (cost, zero_inflation) = problem.evaluate(1.0, 1.0).unwrap();

    // P(1)=1/4 and P(2)=1/8. Absence of observed zeros does not condition the model on X>0
    assert_eq!(zero_inflation, 0.0);
    assert!((cost - 7.0 * 2.0_f64.ln() / 3.0).abs() < 1e-12);
}

#[test]
fn truncated_profile_ignores_excluded_counts_and_normalizes_sample_size() {
    let histogram = [(0, 12), (1, 2), (2, 1)].into_iter().collect();
    // Scale all retained counts by seven. Even u64::MAX excluded observations must not
    // enter the retained count sum, objective normalization, or component-mean search bound
    let scaled_with_tail = [(0, 84), (1, 14), (2, 7), (3, u64::MAX), (100, u64::MAX)]
        .into_iter()
        .collect();
    let original = ProfiledHistogram::new(&histogram, Some(3)).unwrap();
    let scaled = ProfiledHistogram::new(&scaled_with_tail, Some(3)).unwrap();
    let expected_cost =
        -(12.0 * (12.0_f64 / 15.0).ln() + 2.0 * (2.0_f64 / 15.0).ln() + (1.0_f64 / 15.0).ln())
            / 15.0;

    for problem in [original, scaled] {
        let (cost, zero_inflation) = problem.evaluate(1.0, 1.0).unwrap();
        assert!((cost - expected_cost).abs() < 1e-12);
        assert!((zero_inflation - 0.5).abs() < 1e-12);
        assert!((problem.maximum_log_mean - 20.0_f64.ln()).abs() < 1e-12);
    }
}

#[test]
fn complete_profile_normalizes_integer_fitting_frequencies_without_changing_the_likelihood() {
    // At mean=1 and dispersion=1 the NB component is geometric, with masses
    // [1/2,1/4,1/8] at [0,1,2]. Both histograms have fractions [1/2,1/3,1/6]
    // Their zero fraction matches the component, giving pi=0 and cost (5/3)*ln(2)
    for histogram in [
        [(0, 6), (1, 4), (2, 2)].into_iter().collect(),
        [(0, 42), (1, 28), (2, 14)].into_iter().collect(),
    ] {
        let problem = ProfiledHistogram::new(&histogram, None).unwrap();
        let (cost, zero_inflation) = problem.evaluate(1.0, 1.0).unwrap();

        assert!((cost - (5.0 / 3.0) * 2.0_f64.ln()).abs() < 1e-12);
        assert!(zero_inflation.abs() < 1e-12);
    }
}

#[test]
fn profile_rejects_empty_zero_only_and_overflowing_retained_samples() {
    for (histogram, expected_reason) in [
        (CoverageCounts::new(), "without retained eligible positions"),
        (
            [(0, 0), (1, 0)].into_iter().collect(),
            "without retained eligible positions",
        ),
        (
            [(0, 10), (1, 0), (3, 100)].into_iter().collect(),
            "every retained position has zero coverage",
        ),
    ] {
        let error = ProfiledHistogram::new(&histogram, Some(3)).err().unwrap();
        assert!(error.to_string().contains(expected_reason), "{error:#}");
    }
    let overflowing = [(0, u64::MAX), (1, 1)].into_iter().collect();
    let error = ProfiledHistogram::new(&overflowing, None).err().unwrap();
    assert!(
        error
            .to_string()
            .contains("retained-position count overflow")
    );
}

#[test]
fn three_level_sample_without_excess_zeros_is_not_rejected_as_an_interior_ridge() {
    let histogram = [(0, 0), (1, 2), (2, 1)].into_iter().collect();

    let (boundary, retained_positions) = fit_right_truncated_zip(&histogram, 3).unwrap();

    // With pi=0 and X<3, E[X]=(lambda+lambda^2)/(1+lambda+lambda^2/2)
    // Setting this to the observed mean 4/3 gives lambda^2-lambda-4=0
    let expected_mean = (1.0 + 17.0_f64.sqrt()) / 2.0;
    assert_eq!(retained_positions, 3);
    assert_eq!(boundary.zero_inflation, 0.0);
    assert!((boundary.lambda - expected_mean).abs() < 1e-10);
    check_truncated_identification(&histogram, Some(3), boundary.into()).unwrap();
}
