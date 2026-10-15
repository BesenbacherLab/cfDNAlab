use super::*;

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
    let parameters = fit_zinb(&histogram, Some(16)).unwrap();
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
    let parameters = fit_zinb(&histogram, None).unwrap();
    assert_eq!(parameters.dispersion, 0.0);
    assert_eq!(parameters.component_mean, 0.5);
    assert_eq!(parameters.zero_inflation, 0.0);
}
