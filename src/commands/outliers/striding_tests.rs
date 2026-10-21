use super::*;

fn cores(number_of_cores: usize) -> Vec<CoreHistogram> {
    (0..number_of_cores)
        .map(|index| {
            let start = index as u32 * 500_000;
            let mut histogram =
                CoreHistogram::new(Interval::new(start, start + 500_000).expect("valid core"));
            histogram.counts = [(0, 250_000), (1, 250_000)].into_iter().collect();
            histogram.fragment_support = 100;
            histogram
        })
        .collect()
}

#[test]
fn centered_context_expands_to_five_point_five_megabases() {
    let histograms = cores(13);

    let (combined, context, eligible_positions, complete) =
        build_context(&histograms, &[500_000; 13], 6, 5_000_000, 1_000).expect("valid context");

    assert!(complete);
    assert_eq!(context.as_tuple(), (500_000, 6_000_000));
    assert_eq!(context.len(), 5_500_000);
    assert_eq!(combined.eligible_positions().unwrap(), 5_500_000);
    assert_eq!(eligible_positions, 5_500_000);
    assert_eq!(combined.fragment_support, 1_100);
}

#[test]
fn seven_megabase_minimum_keeps_the_symmetric_whole_core_expansion_rule() {
    let histograms = cores(17);

    let (combined, context, eligible_positions, complete) =
        build_context(&histograms, &[500_000; 17], 8, 7_000_000, 1_000).unwrap();

    // Thirteen centered cores supply only 6.5 Mb. Fifteen supply 7.5 Mb
    assert!(complete);
    assert_eq!(context.as_tuple(), (500_000, 8_000_000));
    assert_eq!(eligible_positions, 7_500_000);
    assert_eq!(combined.fragment_support, 1_500);
    assert_eq!(combined.eligible_positions().unwrap(), 7_500_000);
}

#[test]
fn triangular_frequencies_follow_distance_without_changing_coverage_levels() {
    let mut histograms = cores(5);
    for (index, histogram) in histograms.iter_mut().enumerate() {
        histogram.counts = [(index as u32, index as u64 + 1)].into_iter().collect();
    }
    let original = histograms.clone();

    let weighted = triangular_context_counts(
        &histograms,
        2,
        Interval::new(0, 2_500_000).unwrap(),
        &TriangularWeightCache::default(),
    )
    .unwrap();

    // Multiply frequencies [1,2,3,4,5] by the triangle [1,2,3,2,1]
    // The single observation at coverage zero remains present at the outer edge
    assert_eq!(
        weighted,
        [(0, 1), (1, 4), (2, 9), (3, 8), (4, 5)]
            .into_iter()
            .collect()
    );
    assert_eq!(histograms, original);
}

#[test]
fn histogram_pooling_excludes_cores_outside_the_recorded_context() {
    let mut histograms = cores(5);
    histograms[0].counts = [(100, 1)].into_iter().collect();
    histograms[1].counts = [(1, 2)].into_iter().collect();
    histograms[2].counts = [(1, 3)].into_iter().collect();
    histograms[3].counts = [(1, 5)].into_iter().collect();
    histograms[4].counts = [(200, 1)].into_iter().collect();

    let weighted = triangular_context_counts(
        &histograms,
        2,
        Interval::new(500_000, 2_000_000).unwrap(),
        &TriangularWeightCache::default(),
    )
    .unwrap();

    // Pool the three included cores at coverage 1: 2*1 + 3*2 + 5*1 = 13
    // The high-coverage observations outside the recorded context must not enter the fit
    assert_eq!(weighted, [(1, 13)].into_iter().collect());
}

#[test]
fn triangle_keeps_masked_distances_and_uses_actual_support_of_sparse_and_partial_cores() {
    let mut histograms = cores(5);
    histograms[0].counts = [(0, 10)].into_iter().collect();
    histograms[1].counts.clear();
    histograms[2].counts = [(2, 1)].into_iter().collect();
    histograms[3].counts = [(3, 2)].into_iter().collect();
    histograms[4].interval = Interval::new(2_000_000, 2_000_030).unwrap();
    histograms[4].counts = [(4, 5)].into_iter().collect();

    let weighted = triangular_context_counts(
        &histograms,
        2,
        Interval::new(0, 2_000_030).unwrap(),
        &TriangularWeightCache::default(),
    )
    .unwrap();

    // The empty core still occupies a distance step. Do not multiply frequencies by an
    // eligible fraction again, or give the sparse central core the mass of a complete core
    assert_eq!(
        weighted,
        [(0, 10), (2, 3), (3, 4), (4, 5)].into_iter().collect()
    );
}

#[test]
fn blend_fraction_normalizes_actual_eligible_support_before_mixing() {
    // A left core supplies 100 zeros, the center supplies 50 ones, and the right
    // core is masked. The [1,2,1] profile therefore has 200 fitting observations
    let uniform = [(0, 100), (1, 50)].into_iter().collect();
    let triangular = [(0, 100), (1, 100)].into_iter().collect();

    let blended = blend_context_counts(&uniform, triangular, 150, 0.2).unwrap();

    // Equal-total mixing gives 0.8*{100,50} + 0.2*(150/200)*{100,100} = {95,55}
    // A common million-fold scale preserves fractional frequencies before integer fitting
    assert_eq!(
        blended,
        [(0, 95_000_000), (1, 55_000_000)].into_iter().collect()
    );
    assert_eq!(blended.values().sum::<u64>(), 150_000_000);
    assert_eq!(uniform, [(0, 100), (1, 50)].into_iter().collect());
}

#[test]
fn blend_endpoints_preserve_the_original_integer_frequencies_exactly() {
    let uniform = [(0, 100), (1, 50)].into_iter().collect();
    let triangular: CoverageCounts = [(0, 100), (1, 100)].into_iter().collect();

    assert_eq!(
        blend_context_counts(&uniform, triangular.clone(), 150, 0.0).unwrap(),
        uniform
    );
    assert_eq!(
        blend_context_counts(&uniform, triangular.clone(), 150, 1.0).unwrap(),
        triangular
    );
}

#[test]
fn blending_identical_distributions_preserves_their_relative_frequencies() {
    let uniform = [(1, 99), (20, 1)].into_iter().collect();
    let triangular = [(1, 198), (20, 2)].into_iter().collect();

    let blended = blend_context_counts(&uniform, triangular, 100, 0.2).unwrap();

    assert_eq!(
        blended,
        [(1, 99_000_000), (20, 1_000_000)].into_iter().collect()
    );
}

#[test]
fn blend_discretization_keeps_a_rare_coverage_level_with_subunit_frequency() {
    let uniform = [(0, 999_999_999), (100, 1)].into_iter().collect();
    let triangular = [(0, 999_999_999_000_000_000), (100, 1)]
        .into_iter()
        .collect();

    let blended =
        blend_context_counts(&uniform, triangular, 1_000_000_000, 1.0 - f64::EPSILON).unwrap();

    // The rare level has approximately 0.001 scaled units. Rounding it to zero would
    // erase an observed coverage level. Keep the smallest positive fitting frequency
    assert_eq!(blended[&100], 1);
    assert_eq!(blended.len(), 2);
}

#[test]
fn both_fitting_stages_use_the_selected_weighting_but_report_actual_position_counts() {
    let mut histograms = cores(3);
    for (core_index, histogram) in histograms.iter_mut().enumerate() {
        let start = core_index as u32 * 100;
        histogram.interval = Interval::new(start, start + 100).unwrap();
        histogram.counts = [(0, 50), (1, 25), (2, 25)].into_iter().collect();
        histogram.fragment_support = 1;
    }
    histograms[1].counts = [(1, 50), (2, 49), (20, 1)].into_iter().collect();
    let chromosomes = vec!["chr1".to_string()];
    let histograms = [("chr1".to_string(), histograms)].into_iter().collect();
    let pooled = sum_global_histogram(&histograms).unwrap();

    // Uniform frequencies have initial mean 318/300 and retained mean 298/299
    // Triangular frequencies have initial mean 486/400 and retained mean 446/398
    // A 20% blend gives {0:95,1:102.5,2:101.4,20:1.1}. Its initial mean is
    // 327.3/300 = 1.091 and its retained mean is 305.3/298.9 = 3053/2989
    for (weighting, expected_initial_mean, expected_retained_mean) in [
        (ContextWeighting::Uniform, 1.06, 298.0 / 299.0),
        (ContextWeighting::Triangular, 1.215, 223.0 / 199.0),
        (
            ContextWeighting::TriangularBlend(0.2),
            1.091,
            3053.0 / 2989.0,
        ),
        (ContextWeighting::TriangularBlend(0.0), 1.06, 298.0 / 299.0),
        (ContextWeighting::TriangularBlend(1.0), 1.215, 223.0 / 199.0),
    ] {
        let fitted = fit_core_models(
            &chromosomes,
            &histograms,
            &pooled,
            OutlierModel::Zip,
            0.05,
            300,
            3,
            weighting,
        )
        .unwrap();

        let model = fitted.by_chromosome["chr1"][1].as_ref().unwrap();
        // The center's [1,2,1] triangle gives frequencies {0:100,1:150,2:148,20:2}
        // All choices have too few zeros for positive pi. Their Poisson fits have
        // P(X>=3) > 0.05 and P(X>=4) < 0.05, so T1 is exactly 4 for every choice
        assert_eq!(model.source, ModelSource::Local);
        assert_eq!(model.core.as_tuple(), (100, 200));
        assert_eq!(model.context.as_tuple(), (0, 300));
        assert!((model.model.initial.component_mean - expected_initial_mean).abs() < 1e-12);
        assert_eq!(model.model.initial.zero_inflation, 0.0);
        assert_eq!(model.model.initial_threshold, 4);
        assert_eq!(model.eligible_positions, 300);
        assert_eq!(model.fragment_support, 3);
        assert_eq!(model.model.second_fit_retained_positions, 299);
        assert_eq!(model.diagnostics.observed_initial_tail_positions, 1);
        // Observed positive counts remain unweighted: 100 ones, 99 twos and a single 20
        // Their mean is 318/200 = 1.59 and variance is 896/200 - 1.59^2 = 1.9519
        assert!((model.diagnostics.observed_positive_coverage_mean - 1.59).abs() < 1e-12);
        assert!((model.diagnostics.observed_positive_coverage_variance - 1.9519).abs() < 1e-12);
        // The second fit has pi=0 too. Its conditional mean on X<4 must equal the
        // retained mean under the selected weighting. This is the closed-form
        // truncated Poisson mean, independent of the fitter's conditional-mean helper
        let mean = model.model.underlying_fit.component_mean;
        let conditional_mean = (mean + mean.powi(2) + mean.powi(3) / 2.0)
            / (1.0 + mean + mean.powi(2) / 2.0 + mean.powi(3) / 6.0);
        assert_eq!(model.model.underlying_fit.zero_inflation, 0.0);
        assert!((conditional_mean - expected_retained_mean).abs() < 1e-10);
        assert!(fitted.global.is_none());
    }
}

#[test]
fn second_fit_excludes_a_neighboring_cores_tail_even_when_the_center_has_no_outliers() {
    let mut histograms = cores(2);
    for (core_index, histogram) in histograms.iter_mut().enumerate() {
        let start = core_index as u32 * 100;
        histogram.interval = Interval::new(start, start + 100).unwrap();
        histogram.fragment_support = 1;
    }
    histograms[0].counts = [(0, 50), (1, 49), (20, 1)].into_iter().collect();
    histograms[1].counts = [(0, 50), (1, 50)].into_iter().collect();
    let chromosomes = vec!["chr1".to_string()];
    let histograms = [("chr1".to_string(), histograms)].into_iter().collect();
    let pooled = sum_global_histogram(&histograms).unwrap();

    let fitted = fit_core_models(
        &chromosomes,
        &histograms,
        &pooled,
        OutlierModel::Zip,
        0.05,
        200,
        1,
        ContextWeighting::Uniform,
    )
    .unwrap();

    let model = fitted.by_chromosome["chr1"][1].as_ref().unwrap();
    // Full-context mean is (99+20)/200 = 0.595, with too few zeros for pi>0
    // Poisson(0.595) has P(X>=2)>0.05 and P(X>=3)<0.05, so T1=3
    assert!((model.model.initial.component_mean - 0.595).abs() < 1e-12);
    assert_eq!(model.model.initial_threshold, 3);
    assert_eq!(model.model.second_fit_retained_positions, 199);
    // Only the neighbor has an extreme base. The second fit must still exclude it
    // The remaining {0:100,1:99} has conditional mean 99/199 on X<3
    let mean = model.model.underlying_fit.component_mean;
    let conditional_mean = (mean + mean.powi(2)) / (1.0 + mean + mean.powi(2) / 2.0);
    assert_eq!(model.model.underlying_fit.zero_inflation, 0.0);
    assert!((conditional_mean - 99.0 / 199.0).abs() < 1e-10);
}

#[test]
fn triangular_multipliers_do_not_supply_missing_position_or_fragment_support() {
    let mut histograms = cores(3);
    histograms[1].counts.clear();
    histograms[1].fragment_support = 0;
    let chromosomes = vec!["chr1".to_string()];
    let histograms = [("chr1".to_string(), histograms)].into_iter().collect();
    let pooled = sum_global_histogram(&histograms).unwrap();

    for (minimum_positions, minimum_fragments) in [(1_000_001, 1), (1_000_000, 201)] {
        let fitted = fit_core_models(
            &chromosomes,
            &histograms,
            &pooled,
            OutlierModel::Zip,
            0.01,
            minimum_positions,
            minimum_fragments,
            ContextWeighting::Triangular,
        )
        .unwrap();

        // The two eligible cores supply exactly 1 Mb and 200 fragments. Multipliers cannot
        // supply an additional required base or fragment, or change the unweighted pooled fit
        let pooled_model = fitted.global.unwrap();
        assert!((pooled_model.initial.component_mean - 0.5).abs() < 1e-12);
        assert_eq!(pooled_model.initial.zero_inflation, 0.0);
        assert!(fitted.by_chromosome["chr1"][1].is_none());
        for core_index in [0, 2] {
            let model = fitted.by_chromosome["chr1"][core_index].as_ref().unwrap();
            assert_eq!(model.source, ModelSource::GlobalFallback);
            assert_eq!(model.eligible_positions, 1_000_000);
            assert_eq!(model.fragment_support, 200);
            assert_eq!(model.model, pooled_model);
        }
    }
}

#[test]
fn failed_triangular_contexts_record_raw_counts_and_the_fitting_frequencies_separately() {
    let mut histograms = cores(3);
    for histogram in &mut histograms {
        histogram.counts = [(0, 500_000)].into_iter().collect();
    }
    let chromosomes = vec!["chr1".to_string()];
    let histograms = [("chr1".to_string(), histograms)].into_iter().collect();
    let pooled = sum_global_histogram(&histograms).unwrap();

    let error = fit_core_models(
        &chromosomes,
        &histograms,
        &pooled,
        OutlierModel::Zip,
        0.01,
        1_500_000,
        300,
        ContextWeighting::Triangular,
    )
    .unwrap_err();

    let failures = error.downcast_ref::<AnalysisFailures>().unwrap();
    assert_eq!(failures.messages.len(), 3);
    // Each context has 1.5 Mb actual zeros. Edge multipliers [3,2,1] sum to 6
    // The central triangle [1,2,1] sums to 4, giving different fitting frequencies
    let rows: Vec<_> = failures
        .histogram
        .as_ref()
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .collect();
    assert_eq!(
        rows,
        vec![
            "chr1\t0\t500000\t0\t1500000\t300\t0\t1500000\t3000000",
            "chr1\t500000\t1000000\t0\t1500000\t300\t0\t1500000\t2000000",
            "chr1\t1000000\t1500000\t0\t1500000\t300\t0\t1500000\t3000000",
        ]
    );
}

#[test]
fn triangular_frequencies_reject_multiplication_and_sum_overflow() {
    let mut histograms = cores(2);
    let weight_cache = TriangularWeightCache::default();
    let context = Interval::new(0, 1_000_000).unwrap();
    histograms[0].counts = [(1, u64::MAX / 2 + 1)].into_iter().collect();
    histograms[1].counts = [(1, 1)].into_iter().collect();

    // The center's multiplier is 2, overflowing before the neighboring frequency is added
    let error = triangular_context_counts(&histograms, 0, context, &weight_cache).unwrap_err();
    assert!(error.to_string().contains("multiplication overflow"));

    // Multiplication now gives MAX-1, but adding the neighbor's two positions still overflows
    histograms[0].counts = [(1, u64::MAX / 2)].into_iter().collect();
    histograms[1].counts = [(1, 2)].into_iter().collect();
    let error = triangular_context_counts(&histograms, 0, context, &weight_cache).unwrap_err();
    assert!(error.to_string().contains("sum overflow"));
}

#[test]
fn edge_context_continues_on_available_side() {
    let histograms = cores(13);

    // Either edge needs ten whole cores. Each contributes 500 kb eligible and 100 fragments
    for (core_index, expected_context) in [(0, (0, 5_000_000)), (12, (1_500_000, 6_500_000))] {
        let (combined, context, eligible_positions, complete) =
            build_context(&histograms, &[500_000; 13], core_index, 5_000_000, 1_000)
                .expect("valid edge context");

        assert!(complete);
        assert_eq!(context.as_tuple(), expected_context);
        assert_eq!(eligible_positions, 5_000_000);
        assert_eq!(combined.eligible_positions().unwrap(), 5_000_000);
        assert_eq!(combined.fragment_support, 1_000);
    }
}

#[test]
fn a_single_core_at_both_support_minima_needs_no_expansion_or_fallback() {
    let histograms = cores(1);

    for (minimum_positions, minimum_fragments, expected_complete) in [
        (500_000, 100, true),
        (500_001, 100, false),
        (500_000, 101, false),
    ] {
        let (combined, context, eligible_positions, complete) = build_context(
            &histograms,
            &[500_000],
            0,
            minimum_positions,
            minimum_fragments,
        )
        .unwrap();

        assert_eq!(complete, expected_complete);
        assert_eq!(context.as_tuple(), (0, 500_000));
        assert_eq!(eligible_positions, 500_000);
        assert_eq!(
            combined.counts,
            [(0, 250_000), (1, 250_000)].into_iter().collect()
        );
        assert_eq!(combined.fragment_support, 100);
        assert_eq!(
            needs_global_fallback(
                &histograms,
                &[500_000],
                minimum_positions,
                minimum_fragments,
            )
            .unwrap(),
            !expected_complete
        );
    }
}

#[test]
fn incomplete_contig_is_marked_for_global_fallback() {
    let histograms = cores(4);

    let (combined, context, eligible_positions, complete) =
        build_context(&histograms, &[500_000; 4], 1, 5_000_000, 1_000)
            .expect("valid short context");

    assert!(!complete);
    assert_eq!(context.as_tuple(), (0, 2_000_000));
    assert_eq!(eligible_positions, 2_000_000);
    assert_eq!(combined.fragment_support, 400);
}

#[test]
fn fallback_need_checks_full_chromosome_eligible_positions_and_summed_support() {
    let histograms = cores(4);
    assert!(!needs_global_fallback(&histograms, &[500_000; 4], 2_000_000, 400).unwrap());
    assert!(needs_global_fallback(&histograms, &[500_000; 4], 2_000_001, 400).unwrap());
    assert!(needs_global_fallback(&histograms, &[500_000; 4], 2_000_000, 401).unwrap());
}

#[test]
fn masked_positions_require_wider_centered_contexts() {
    let mut histograms = cores(13);
    for (index, histogram) in histograms.iter_mut().enumerate() {
        if index != 6 {
            histogram.counts = [(0, 125_000), (1, 125_000)].into_iter().collect();
        }
    }

    let mut eligible_positions_by_core = [250_000; 13];
    eligible_positions_by_core[6] = 500_000;
    let (combined, context, eligible_positions, complete) =
        build_context(&histograms, &eligible_positions_by_core, 6, 2_000_000, 100).unwrap();

    // The central 500 kb plus six half-eligible neighbors supply exactly 2 Mb eligible
    // Five cores span 2.5 Mb physically but contain only 1.5 Mb eligible positions
    assert!(complete);
    assert_eq!(context.as_tuple(), (1_500_000, 5_000_000));
    assert_eq!(combined.eligible_positions().unwrap(), 2_000_000);
    assert_eq!(eligible_positions, 2_000_000);
    assert_eq!(
        combined.counts,
        [(0, 1_000_000), (1, 1_000_000)].into_iter().collect()
    );
    assert_eq!(combined.fragment_support, 700);
}

#[test]
fn edge_context_crosses_masked_cores_without_counting_them_as_support() {
    let mut histograms = cores(7);
    for histogram in &mut histograms[1..5] {
        histogram.counts.clear();
        histogram.fragment_support = 0;
    }

    let (combined, context, eligible_positions, complete) = build_context(
        &histograms,
        &[500_000, 0, 0, 0, 0, 500_000, 500_000],
        0,
        1_000_000,
        100,
    )
    .unwrap();

    // Cores 0 and 5 supply 1 Mb eligible. The intervening 2 Mb mask stays in the coordinates
    assert!(complete);
    assert_eq!(context.as_tuple(), (0, 3_000_000));
    assert_eq!(combined.eligible_positions().unwrap(), 1_000_000);
    assert_eq!(eligible_positions, 1_000_000);
    assert_eq!(combined.fragment_support, 200);
}

#[test]
fn eligible_zeros_count_toward_context_size_independently_of_fragment_support() {
    let mut histograms = cores(3);
    histograms[1].counts = [(0, 500_000)].into_iter().collect();

    let (combined, context, eligible_positions, complete) =
        build_context(&histograms, &[500_000; 3], 0, 1_000_000, 200).unwrap();

    assert!(complete);
    assert_eq!(context.as_tuple(), (0, 1_000_000));
    assert_eq!(
        combined.counts,
        [(0, 750_000), (1, 250_000)].into_iter().collect()
    );
    assert_eq!(combined.eligible_positions().unwrap(), 1_000_000);
    assert_eq!(eligible_positions, 1_000_000);
    assert_eq!(combined.fragment_support, 200);
}

#[test]
fn enough_eligible_positions_does_not_bypass_fragment_support() {
    let histograms = cores(5);

    for (minimum_fragments, expected_complete) in [(500, true), (501, false)] {
        let (combined, context, eligible_positions, complete) =
            build_context(&histograms, &[500_000; 5], 2, 500_000, minimum_fragments).unwrap();

        assert_eq!(complete, expected_complete);
        assert_eq!(context.as_tuple(), (0, 2_500_000));
        assert_eq!(combined.eligible_positions().unwrap(), 2_500_000);
        assert_eq!(eligible_positions, 2_500_000);
        assert_eq!(combined.fragment_support, 500);
    }
}

#[test]
fn physically_long_chromosome_requires_fallback_when_eligible_bases_are_insufficient() {
    let mut histograms = cores(4);
    histograms[1].counts = [(0, 250_000)].into_iter().collect();

    // A 2 Mb chromosome has 1.75 Mb eligible. One additional required base needs fallback
    let eligible_positions_by_core = [500_000, 250_000, 500_000, 500_000];
    assert!(
        !needs_global_fallback(&histograms, &eligible_positions_by_core, 1_750_000, 400).unwrap()
    );
    assert!(
        needs_global_fallback(&histograms, &eligible_positions_by_core, 1_750_001, 400).unwrap()
    );
    let (combined, context, eligible_positions, complete) =
        build_context(&histograms, &eligible_positions_by_core, 1, 1_750_001, 400).unwrap();
    assert!(!complete);
    assert_eq!(context.as_tuple(), (0, 2_000_000));
    assert_eq!(combined.eligible_positions().unwrap(), 1_750_000);
    assert_eq!(eligible_positions, 1_750_000);
    assert_eq!(combined.fragment_support, 400);
}

#[test]
fn partially_masked_terminal_core_contributes_only_its_observed_eligible_bases() {
    let mut histograms = cores(2);
    let mut terminal = CoreHistogram::new(Interval::new(1_000_000, 1_000_030).unwrap());
    terminal.counts = [(0, 10), (1, 10)].into_iter().collect();
    terminal.fragment_support = 1;
    histograms.push(terminal);

    // The last 30 bp core contributes 20 eligible bases, rather than its physical length or stride
    for (minimum_positions, expected_complete) in [(1_000_020, true), (1_000_021, false)] {
        let (combined, context, eligible_positions, complete) = build_context(
            &histograms,
            &[500_000, 500_000, 20],
            2,
            minimum_positions,
            201,
        )
        .unwrap();

        assert_eq!(complete, expected_complete);
        assert_eq!(context.as_tuple(), (0, 1_000_030));
        assert_eq!(combined.eligible_positions().unwrap(), 1_000_020);
        assert_eq!(eligible_positions, 1_000_020);
        assert_eq!(combined.fragment_support, 201);
        assert_eq!(
            needs_global_fallback(&histograms, &[500_000, 500_000, 20], minimum_positions, 201)
                .unwrap(),
            !expected_complete
        );
    }
}

#[test]
fn supported_local_models_skip_the_pooled_fit() {
    let chromosomes = vec!["chr1".to_string()];
    let histograms = [("chr1".to_string(), cores(2))].into_iter().collect();
    // The deliberately unusable pooled histogram cannot matter when local contexts suffice
    let mut pooled = CoreHistogram::new(Interval::new(0, 1).unwrap());
    pooled.counts.insert(0, 1_000);
    let FittedModels {
        global,
        by_chromosome: models,
    } = fit_core_models(
        &chromosomes,
        &histograms,
        &pooled,
        OutlierModel::Zip,
        0.01,
        500_000,
        100,
        ContextWeighting::Uniform,
    )
    .unwrap();
    assert!(global.is_none());
    assert_eq!(models["chr1"].len(), 2);
    assert!(
        models["chr1"]
            .iter()
            .all(|model| model.as_ref().unwrap().source == ModelSource::Local)
    );
}

#[test]
fn fully_blacklisted_cores_keep_empty_grid_slots_between_eligible_models() {
    let chromosomes = vec!["chr1".to_string()];
    let mut histograms = cores(5);
    for index in [0, 2, 4] {
        histograms[index].counts.clear();
        histograms[index].fragment_support = 0;
    }
    let histograms = [("chr1".to_string(), histograms)].into_iter().collect();
    let pooled = sum_global_histogram(&histograms).unwrap();

    let fitted = fit_core_models(
        &chromosomes,
        &histograms,
        &pooled,
        OutlierModel::Zip,
        0.01,
        500_000,
        100,
        ContextWeighting::Uniform,
    )
    .unwrap();

    // Masking the first, middle and last cores must not move either remaining grid slot
    let models = &fitted.by_chromosome["chr1"];
    assert!(fitted.global.is_none());
    assert_eq!(models.len(), 5);
    for index in [0, 2, 4] {
        assert!(models[index].is_none());
    }
    for (index, bounds) in [(1, (500_000, 1_000_000)), (3, (1_500_000, 2_000_000))] {
        let model = models[index].as_ref().unwrap();
        assert_eq!(model.core.as_tuple(), bounds);
        assert_eq!(model.context.as_tuple(), bounds);
        assert_eq!(model.eligible_positions, 500_000);
    }
}

#[test]
fn fully_blacklisted_chromosome_does_not_request_an_unused_global_fit() {
    let chromosomes = vec!["chr1".to_string(), "masked".to_string()];
    let mut masked = cores(1);
    // A stored zero-count bin is also empty support, rather than observed zero coverage
    masked[0].counts = [(0, 0)].into_iter().collect();
    masked[0].fragment_support = 0;
    let histograms = [
        ("chr1".to_string(), cores(2)),
        ("masked".to_string(), masked),
    ]
    .into_iter()
    .collect();
    let mut unusable_pooled = CoreHistogram::new(Interval::new(0, 1).unwrap());
    unusable_pooled.counts.insert(0, 1_000);

    let fitted = fit_core_models(
        &chromosomes,
        &histograms,
        &unusable_pooled,
        OutlierModel::Zip,
        0.01,
        500_000,
        100,
        ContextWeighting::Uniform,
    )
    .unwrap();

    assert!(fitted.global.is_none());
    assert!(fitted.by_chromosome["chr1"].iter().all(Option::is_some));
    assert_eq!(fitted.by_chromosome["masked"].len(), 1);
    assert!(fitted.by_chromosome["masked"][0].is_none());
}

#[test]
fn global_fallback_is_assigned_only_to_cores_with_eligible_positions() {
    let chromosomes = vec!["chr1".to_string()];
    let mut histograms = cores(3);
    for index in [0, 2] {
        histograms[index].counts.clear();
        histograms[index].fragment_support = 0;
    }
    let histograms = [("chr1".to_string(), histograms)].into_iter().collect();
    let pooled = sum_global_histogram(&histograms).unwrap();

    let fitted = fit_core_models(
        &chromosomes,
        &histograms,
        &pooled,
        OutlierModel::Zip,
        0.01,
        1_000_000,
        100,
        ContextWeighting::Uniform,
    )
    .unwrap();

    // The 1.5 Mb chromosome supplies only 500 kb eligible bases toward the requested 1 Mb
    // Its sole eligible core needs the pooled model. Neither masked core acquires a model
    let models = &fitted.by_chromosome["chr1"];
    assert_eq!(models.len(), 3);
    assert!(models[0].is_none());
    assert!(models[2].is_none());
    let model = models[1].as_ref().unwrap();
    assert_eq!(model.core.as_tuple(), (500_000, 1_000_000));
    assert_eq!(model.context.as_tuple(), (0, 1_500_000));
    assert_eq!(model.eligible_positions, 500_000);
    assert_eq!(model.fragment_support, 100);
    assert_eq!(model.source, ModelSource::GlobalFallback);
    assert_eq!(Some(model.model), fitted.global);
}

#[test]
fn unused_context_does_not_add_a_failure_for_an_eligible_zero_coverage_core() {
    let chromosomes = vec!["chr1".to_string()];
    let masked = CoreHistogram::new(Interval::new(0, 500_000).unwrap());
    let mut observed_zero = CoreHistogram::new(Interval::new(500_000, 1_000_000).unwrap());
    observed_zero.counts.insert(0, 500_000);
    observed_zero.fragment_support = 1_000;
    let histograms = [("chr1".to_string(), vec![masked, observed_zero])]
        .into_iter()
        .collect();
    let pooled = sum_global_histogram(&histograms).unwrap();

    let error = fit_core_models(
        &chromosomes,
        &histograms,
        &pooled,
        OutlierModel::Zip,
        0.01,
        500_000,
        1_000,
        ContextWeighting::Uniform,
    )
    .unwrap_err();

    // Real eligible zeros still fail their required fit. The masked core needs no fit at all
    let failures = error.downcast_ref::<AnalysisFailures>().unwrap();
    assert_eq!(failures.messages.len(), 1);
    assert!(failures.messages[0].contains("core (500000, 1000000)"));
    assert!(failures.messages[0].contains("zero coverage"));
}

#[test]
fn a_required_pooled_fit_failure_is_reported_with_its_histogram() {
    let chromosomes = vec!["chr1".to_string()];
    let mut core = CoreHistogram::new(Interval::new(0, 500_000).unwrap());
    core.counts.insert(0, 500_000);
    let histograms = [("chr1".to_string(), vec![core])].into_iter().collect();
    let pooled = sum_global_histogram(&histograms).unwrap();
    let error = fit_core_models(
        &chromosomes,
        &histograms,
        &pooled,
        OutlierModel::Zip,
        0.01,
        1_000_000,
        1,
        ContextWeighting::Uniform,
    )
    .unwrap_err();
    let failures = error.downcast_ref::<AnalysisFailures>().unwrap();
    assert_eq!(failures.messages.len(), 1);
    assert!(failures.messages[0].contains("required global"));
    assert!(
        failures
            .histogram
            .as_ref()
            .unwrap()
            .contains("GLOBAL\tNA\tNA\tNA\tNA\t0\t0\t500000")
    );
}

#[test]
fn short_contig_assigns_the_global_fallback_model_to_every_core() {
    let chromosomes = vec!["chr1".to_string()];
    let histograms = cores(4);
    let histograms_by_chromosome = [("chr1".to_string(), histograms)].into_iter().collect();
    let global_histogram = sum_global_histogram(&histograms_by_chromosome).unwrap();

    let FittedModels {
        global: global_model,
        by_chromosome: models,
    } = fit_core_models(
        &chromosomes,
        &histograms_by_chromosome,
        &global_histogram,
        OutlierModel::Zip,
        0.01,
        5_000_000,
        1_000,
        ContextWeighting::Uniform,
    )
    .expect("global fallback models");

    let chromosome_models = models.get("chr1").expect("chromosome models");
    assert_eq!(chromosome_models.len(), 4);
    for model in chromosome_models {
        let model = model.as_ref().unwrap();
        assert_eq!(model.source, ModelSource::GlobalFallback);
        assert_eq!(Some(model.model), global_model);
        assert_eq!(model.context.as_tuple(), (0, 2_000_000));
    }
}

#[test]
fn complete_degenerate_context_fails_instead_of_silently_falling_back() {
    let chromosomes = vec!["chr1".to_string()];
    let mut histogram = CoreHistogram::new(Interval::new(0, 500_000).expect("valid core interval"));
    histogram.counts.insert(0, 500_000);
    histogram.fragment_support = 1_000;
    let histograms_by_chromosome = [("chr1".to_string(), vec![histogram])]
        .into_iter()
        .collect();

    let error = fit_core_models(
        &chromosomes,
        &histograms_by_chromosome,
        &sum_global_histogram(&histograms_by_chromosome).unwrap(),
        OutlierModel::Zip,
        0.01,
        500_000,
        1_000,
        ContextWeighting::Uniform,
    )
    .expect_err("a complete all-zero context must fail its local fit");

    assert!(format!("{error:#}").contains("zero coverage"));
}

#[test]
fn core_histogram_excludes_masked_positions() {
    let mut histogram = CoreHistogram::new(Interval::new(0, 4).expect("valid histogram interval"));

    histogram
        .add_coverage(&[0.0, 1.0, 2.0, 2.0], Some(&[0, 1, 0, 0]))
        .expect("valid integer coverage");

    assert_eq!(histogram.counts, [(0, 1), (2, 2)].into_iter().collect());
    assert_eq!(
        histogram.eligible_positions().expect("eligible positions"),
        3
    );
}

#[test]
fn rejects_coverage_at_the_f32_exact_integer_boundary() {
    let mut histogram = CoreHistogram::new(Interval::new(0, 1).expect("valid histogram interval"));

    let error = histogram
        .add_coverage(&[16_777_216.0], None)
        .expect_err("inexact integer coverage must be rejected");

    assert!(error.to_string().contains("f32 exact-integer limit"));
}

#[test]
fn rejects_invalid_raw_coverage_values_and_mask_lengths() {
    for invalid_coverage in [f32::NAN, f32::INFINITY, -1.0, 1.5] {
        let mut histogram =
            CoreHistogram::new(Interval::new(0, 1).expect("valid histogram interval"));
        let error = histogram
            .add_coverage(&[invalid_coverage], None)
            .expect_err("invalid raw coverage must be rejected");
        assert!(error.to_string().contains("coverage"));
    }

    let mut histogram = CoreHistogram::new(Interval::new(0, 2).expect("valid histogram interval"));
    let error = histogram
        .add_coverage(&[0.0, 1.0], Some(&[0]))
        .expect_err("mask and coverage lengths must agree");
    assert!(error.to_string().contains("blacklist-mask length"));
}

#[test]
fn global_histogram_conserves_counts_and_fragment_support_across_chromosomes() {
    let mut chr1_first = CoreHistogram::new(Interval::new(0, 100).expect("valid first core"));
    chr1_first.counts = [(0, 70), (1, 30)].into_iter().collect();
    chr1_first.fragment_support = 2;
    let mut chr1_second = CoreHistogram::new(Interval::new(100, 200).expect("valid second core"));
    chr1_second.counts = [(0, 80), (2, 20)].into_iter().collect();
    chr1_second.fragment_support = 3;
    let mut chr2 = CoreHistogram::new(Interval::new(0, 50).expect("valid chr2 core"));
    chr2.counts = [(0, 40), (1, 10)].into_iter().collect();
    chr2.fragment_support = 1;
    let histograms = [
        ("chr1".to_string(), vec![chr1_first, chr1_second]),
        ("chr2".to_string(), vec![chr2]),
    ]
    .into_iter()
    .collect();

    let global = sum_global_histogram(&histograms).expect("global histogram");

    assert_eq!(
        global.counts,
        [(0, 190), (1, 40), (2, 20)].into_iter().collect()
    );
    assert_eq!(
        global.eligible_positions().expect("eligible positions"),
        250
    );
    assert_eq!(global.fragment_support, 6);
}

#[test]
fn required_global_failure_does_not_hide_independent_local_failures() {
    let mut local_cores = cores(2);
    for core in &mut local_cores {
        core.counts = [(0, 500_000)].into_iter().collect();
    }
    let mut short_core = CoreHistogram::new(Interval::new(0, 250_000).unwrap());
    short_core.counts.insert(0, 250_000);
    let chromosomes = vec!["chr_short".to_string(), "chr_local".to_string()];
    let histograms = [
        ("chr_local".to_string(), local_cores),
        ("chr_short".to_string(), vec![short_core]),
    ]
    .into_iter()
    .collect();
    let pooled = sum_global_histogram(&histograms).unwrap();

    let error = fit_core_models(
        &chromosomes,
        &histograms,
        &pooled,
        OutlierModel::Zip,
        0.01,
        500_000,
        100,
        ContextWeighting::Uniform,
    )
    .unwrap_err();

    // The short chromosome needs the failed pooled model. Its dependency is reported once
    // Both supported 500 kb local contexts still require their own failed-fit evidence
    let failures = error.downcast_ref::<AnalysisFailures>().unwrap();
    assert_eq!(failures.messages.len(), 3);
    assert!(failures.messages[0].contains("required global"));
    assert!(failures.messages[1].contains("chr_local core (0, 500000)"));
    assert!(failures.messages[2].contains("chr_local core (500000, 1000000)"));
    let rows: Vec<_> = failures
        .histogram
        .as_ref()
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .collect();
    assert_eq!(
        rows,
        vec![
            "GLOBAL\tNA\tNA\tNA\tNA\t200\t0\t1250000\t1250000",
            "chr_local\t0\t500000\t0\t500000\t100\t0\t500000\t500000",
            "chr_local\t500000\t1000000\t500000\t1000000\t100\t0\t500000\t500000",
        ]
    );
    assert!(failures.audit.is_none());
}

#[test]
fn failed_context_histograms_exclude_successful_local_contexts() {
    let good_core = cores(1);
    let mut bad_cores = cores(2);
    for core in &mut bad_cores {
        core.counts = [(0, 500_000)].into_iter().collect();
    }
    let chromosomes = vec!["chr_good".to_string(), "chr_bad".to_string()];
    let histograms = [
        ("chr_good".to_string(), good_core),
        ("chr_bad".to_string(), bad_cores),
    ]
    .into_iter()
    .collect();
    let pooled = sum_global_histogram(&histograms).unwrap();

    let error = fit_core_models(
        &chromosomes,
        &histograms,
        &pooled,
        OutlierModel::Zip,
        0.01,
        500_000,
        100,
        ContextWeighting::Uniform,
    )
    .unwrap_err();

    let failures = error.downcast_ref::<AnalysisFailures>().unwrap();
    assert_eq!(failures.messages.len(), 2);
    assert!(
        failures
            .messages
            .iter()
            .all(|message| message.starts_with("chr_bad core"))
    );
    let histogram = failures.histogram.as_ref().unwrap();
    assert!(!histogram.contains("chr_good"));
    assert!(!histogram.contains("GLOBAL\tNA"));
    assert!(histogram.contains("chr_bad\t0\t500000\t0\t500000\t100\t0\t500000"));
    assert!(histogram.contains("chr_bad\t500000\t1000000\t500000\t1000000\t100\t0\t500000"));
}

#[test]
fn histogram_and_support_overflow_are_errors_instead_of_wrapped_counts() {
    let interval = Interval::new(0, 30).unwrap();
    let mut histogram = CoreHistogram::new(interval);
    histogram.counts.insert(1, u64::MAX);
    let error = histogram.add_coverage(&[1.0], None).unwrap_err();
    assert!(error.to_string().contains("histogram count overflow"));

    let mut additional = CoreHistogram::new(interval);
    additional.counts.insert(1, 1);
    let error = histogram.add_histogram(&additional).unwrap_err();
    assert!(error.to_string().contains("histogram count overflow"));

    histogram.counts = [(0, u64::MAX), (1, 1)].into_iter().collect();
    let error = histogram.eligible_positions().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("eligible-position count overflow")
    );

    let mut support = CoreHistogram::new(interval);
    support.fragment_support = u64::MAX;
    let mut extra_support = CoreHistogram::new(interval);
    extra_support.fragment_support = 1;
    let error = support.add_histogram(&extra_support).unwrap_err();
    assert!(error.to_string().contains("fragment support overflow"));
}

#[test]
fn fallback_support_checks_reject_overflow_and_empty_core_lists() {
    assert!(
        needs_global_fallback(&[], &[], 1, 1)
            .unwrap_err()
            .to_string()
            .contains("no model-core histograms")
    );

    let mut histograms = cores(2);
    histograms[0].counts = [(0, u64::MAX)].into_iter().collect();
    histograms[1].counts = [(1, 1)].into_iter().collect();
    let error = needs_global_fallback(&histograms, &[u64::MAX, 1], 1, 1).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("eligible-position count overflow")
    );

    histograms[0].counts = [(0, 500_000)].into_iter().collect();
    histograms[0].fragment_support = u64::MAX;
    histograms[1].fragment_support = 1;
    let error = needs_global_fallback(&histograms, &[500_000, 1], 1, 1).unwrap_err();
    assert!(error.to_string().contains("fragment support overflow"));
}

#[test]
fn incremental_context_eligibility_rejects_overflow_across_different_coverage_bins() {
    let mut histograms = cores(2);
    histograms[0].counts = [(0, u64::MAX)].into_iter().collect();
    histograms[1].counts = [(1, 1)].into_iter().collect();
    for histogram in &mut histograms {
        histogram.fragment_support = 0;
    }

    // Missing fragment support forces expansion. Neither individual frequency overflows
    // Adding the second core would overflow only the total eligible-position count
    let error = build_context(&histograms, &[u64::MAX, 1], 0, 1, 1).unwrap_err();

    assert!(
        error
            .to_string()
            .contains("context eligible-position count overflow")
    );
}
