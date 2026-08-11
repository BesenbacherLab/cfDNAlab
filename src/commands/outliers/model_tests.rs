use super::*;

fn coverage_counts(dense_counts: &[u64]) -> CoverageCounts {
    dense_counts
        .iter()
        .enumerate()
        .filter(|(_, positions)| **positions > 0)
        .map(|(coverage, &positions)| (coverage as u32, positions))
        .collect()
}

fn zip_one_histogram_with_extreme_tail() -> CoverageCounts {
    // These rounded counts represent 1,000,000 positions from ZIP(lambda=1, pi=0.5).
    // Their weighted coverage sum is exactly 500,000 before the added extreme tail.
    let mut histogram = vec![0_u64; 21];
    histogram[0] = 683_939;
    histogram[1] = 183_940;
    histogram[2] = 91_970;
    histogram[3] = 30_657;
    histogram[4] = 7_664;
    histogram[5] = 1_533;
    histogram[6] = 256;
    histogram[7] = 36;
    histogram[8] = 5;
    histogram[20] = 100;
    coverage_counts(&histogram)
}

#[test]
fn returns_four_as_poisson_one_five_percent_upper_tail_threshold() {
    // P(X >= 3) is about 0.0803, while P(X >= 4) is about 0.0190.
    let threshold = zip_upper_tail_threshold(
        ZipParameters {
            lambda: 1.0,
            zero_inflation: 0.0,
        },
        0.05,
    )
    .expect("valid Poisson threshold");

    assert_eq!(threshold, 4);
}

#[test]
fn returns_seventy_for_poisson_one_at_extreme_tail_probability() {
    // For lambda = 1, P(X = 69) = exp(-1) / 69! is about 2.15e-99, so P(X >= 69) is
    // greater than 1e-100. Starting at 70, each successive probability is at most
    // 1/71 of the previous value. Therefore P(X >= 70) is less than
    // (exp(-1) / 70!) / (1 - 1/71), which is about 3.12e-101.
    let threshold = zip_upper_tail_threshold(
        ZipParameters {
            lambda: 1.0,
            zero_inflation: 0.0,
        },
        1e-100,
    )
    .expect("valid extreme Poisson threshold");

    assert_eq!(threshold, 70);
}

#[test]
fn right_truncated_fit_recovers_underlying_zip_after_extreme_tail_removal() {
    let histogram = zip_one_histogram_with_extreme_tail();

    // Excluding coverage >= 7 removes the synthetic extreme tail and conditions the likelihood
    // on the retained support. The underlying parameters should remain lambda=1 and pi=0.5.
    let (fit, retained_positions) =
        fit_right_truncated_zip(&histogram, 7).expect("right-truncated ZIP fit");

    assert_eq!(retained_positions, 999_959);
    assert!((fit.lambda - 1.0).abs() < 0.01, "lambda={}", fit.lambda);
    assert!(
        (fit.zero_inflation - 0.5).abs() < 0.01,
        "zero_inflation={}",
        fit.zero_inflation
    );
    assert!((fit.mean() - 0.5).abs() < 0.01, "mean={}", fit.mean());
}

#[test]
fn two_stage_fit_calls_original_histogram_but_targets_tail_excluded_mean() {
    let histogram = zip_one_histogram_with_extreme_tail();

    let model = fit_two_stage_zip(&histogram, 0.0001).expect("two-stage ZIP fit");

    // The complete fit calls coverage >= 7. Removing those bins leaves exactly 999,959
    // observations and recovers the underlying ZIP(lambda = 1, pi = 0.5), whose 0.0001 upper-tail
    // threshold is also 7
    assert_eq!(model.initial_threshold, 7);
    assert_eq!(model.second_fit_retained_positions, 999_959);
    assert!((model.underlying_fit.mean() - 0.5).abs() < 0.01);
    assert_eq!(model.final_threshold, 7);
    assert!(
        model
            .underlying_fit
            .survival_at_or_above(model.final_threshold)
            .expect("final threshold survival")
            <= 0.0001
    );
    assert!(
        model
            .underlying_fit
            .survival_at_or_above(model.final_threshold - 1)
            .expect("predecessor survival")
            > 0.0001
    );
}

#[test]
fn rejects_histogram_without_positive_coverage() {
    let error = fit_two_stage_zip(&coverage_counts(&[1_000]), 0.01)
        .expect_err("all-zero histogram should not produce a ZIP fit");

    assert!(error.to_string().contains("zero coverage"));
}

#[test]
fn untruncated_fit_uses_the_poisson_boundary_without_excess_zeros() {
    let histogram = [(0, 100), (1, 100)].into_iter().collect();

    let fit = fit_zip(&histogram).expect("Poisson-boundary ZIP fit");

    assert_eq!(fit.zero_inflation, 0.0);
    assert!((fit.lambda - 0.5).abs() < 1e-12);
}

#[test]
fn right_truncated_fit_uses_the_poisson_boundary_and_matches_the_retained_mean() {
    let histogram = [(0, 500), (1, 400), (2, 100)].into_iter().collect();

    let (fit, retained_positions) =
        fit_right_truncated_zip(&histogram, 4).expect("right-truncated Poisson fit");

    assert_eq!(retained_positions, 1_000);
    assert_eq!(fit.zero_inflation, 0.0);
    // The complete retained sample has mean (400 + 2 * 100) / 1,000 = 0.6. The fitted
    // underlying lambda is slightly above 0.6 because its likelihood conditions on X < 4.
    assert!((truncated_poisson_mean(fit.lambda, 0, 3) - 0.6).abs() < 1e-12);
    assert!(fit.lambda > 0.6);
}

#[test]
fn rejects_initial_threshold_that_leaves_no_positive_second_fit_support() {
    let histogram = [(0, 100), (1, 10)].into_iter().collect();

    let error = fit_right_truncated_zip(&histogram, 1)
        .expect_err("T1=1 must leave no retained positive support");

    assert!(error.to_string().contains("leaves no positive coverage"));
}

#[test]
fn zero_inflation_can_make_one_the_smallest_valid_positive_threshold() {
    let threshold = zip_upper_tail_threshold(
        ZipParameters {
            lambda: 5.0,
            zero_inflation: 0.99,
        },
        0.05,
    )
    .expect("valid zero-inflated threshold");

    // Only one percent of values enter the Poisson component, so P(X >= 1) is below 0.01.
    assert_eq!(threshold, 1);
}

#[test]
fn thresholds_satisfy_the_inclusive_tail_boundary_and_probability_monotonicity() {
    let parameters = ZipParameters {
        lambda: 3.0,
        zero_inflation: 0.25,
    };
    let probabilities = [0.1, 0.01, 0.000_1];
    let mut previous_threshold = 0;

    for probability in probabilities {
        let threshold =
            zip_upper_tail_threshold(parameters, probability).expect("valid ZIP threshold");
        assert!(threshold >= previous_threshold);
        assert!(
            parameters
                .survival_at_or_above(threshold)
                .expect("threshold survival")
                <= probability
        );
        if threshold > 1 {
            assert!(
                parameters
                    .survival_at_or_above(threshold - 1)
                    .expect("predecessor survival")
                    > probability
            );
        }
        previous_threshold = threshold;
    }
}

#[test]
fn rejects_nonfinite_and_out_of_range_tail_probabilities() {
    let parameters = ZipParameters {
        lambda: 1.0,
        zero_inflation: 0.0,
    };

    for probability in [f64::NAN, f64::INFINITY, 0.0, -0.1, 0.500_001] {
        let error = zip_upper_tail_threshold(parameters, probability)
            .expect_err("invalid probability must be rejected");
        assert!(error.to_string().contains("tail probability"));
    }
}

#[test]
fn inclusive_survival_at_one_excludes_only_poisson_zero() {
    let parameters = ZipParameters {
        lambda: 1.0,
        zero_inflation: 0.0,
    };

    let survival = parameters
        .survival_at_or_above(1)
        .expect("valid Poisson survival");

    assert!((survival - (1.0 - (-1.0_f64).exp())).abs() < 1e-15);
}

#[test]
fn diagnostics_report_hand_derived_positive_variance_and_final_tail() {
    let histogram = [(0, 4), (1, 2), (3, 2), (5, 1)].into_iter().collect();
    let parameters = ZipParameters {
        lambda: 1.0,
        zero_inflation: 0.0,
    };
    let model = TwoStageZipModel {
        initial: parameters,
        initial_threshold: 4,
        underlying_fit: parameters,
        final_threshold: 4,
        second_fit_retained_positions: 8,
    };

    let diagnostics = diagnose_zip_fit(&histogram, model).expect("valid ZIP diagnostics");

    // Positive values are [1, 1, 3, 3, 5], with mean 2.6 and population variance 2.24.
    assert!((diagnostics.observed_positive_coverage_mean - 2.6).abs() < 1e-12);
    assert!((diagnostics.observed_positive_coverage_variance - 2.24).abs() < 1e-12);
    // For Poisson(1), Var(X | X > 0) is approximately 0.661303.
    assert!(
        (diagnostics.underlying_zip_expected_positive_coverage_variance
            - 0.661_303_112_661_534)
            .abs()
            < 1e-12
    );
    assert_eq!(diagnostics.observed_initial_tail_positions, 1);
    assert!((diagnostics.expected_initial_tail_positions - 0.170_893_411_885).abs() < 1e-12);
    assert_eq!(diagnostics.observed_final_tail_positions, 1);
    // Nine positions times P(Poisson(1) >= 4) is approximately 0.170893.
    assert!((diagnostics.expected_final_tail_positions - 0.170_893_411_885).abs() < 1e-12);
}

#[test]
fn diagnostics_mark_positive_statistics_unavailable_for_an_all_zero_context() {
    let histogram = [(0, 100)].into_iter().collect();
    let parameters = ZipParameters {
        lambda: 1.0,
        zero_inflation: 0.5,
    };
    let model = TwoStageZipModel {
        initial: parameters,
        initial_threshold: 4,
        underlying_fit: parameters,
        final_threshold: 4,
        second_fit_retained_positions: 100,
    };

    let diagnostics = diagnose_zip_fit(&histogram, model).expect("fallback diagnostics");

    assert!(diagnostics.observed_positive_coverage_mean.is_nan());
    assert!(diagnostics.observed_positive_coverage_variance.is_nan());
    assert!(diagnostics.positive_coverage_variance_ratio.is_nan());
    assert_eq!(diagnostics.observed_initial_tail_positions, 0);
    assert!(diagnostics.expected_initial_tail_positions > 0.0);
    assert_eq!(diagnostics.observed_final_tail_positions, 0);
    assert!(diagnostics.expected_final_tail_positions > 0.0);
}
