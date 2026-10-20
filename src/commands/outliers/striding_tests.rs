use super::*;

fn cores(number_of_cores: usize) -> Vec<CoreHistogram> {
    (0..number_of_cores)
        .map(|index| {
            let start = index as u32 * 500_000;
            let mut histogram =
                CoreHistogram::new(Interval::new(start, start + 500_000).expect("valid core"));
            histogram.counts = [(0, 500), (1, 500)].into_iter().collect();
            histogram.fragment_support = 100;
            histogram
        })
        .collect()
}

#[test]
fn centered_context_expands_to_five_point_five_megabases() {
    let histograms = cores(13);

    let (combined, context, complete) =
        build_context(&histograms, 6, 5_000_000, 1_000).expect("valid context");

    assert!(complete);
    assert_eq!(context.as_tuple(), (500_000, 6_000_000));
    assert_eq!(context.len(), 5_500_000);
    assert_eq!(combined.fragment_support, 1_100);
}

#[test]
fn edge_context_continues_on_available_side() {
    let histograms = cores(13);

    let (combined, context, complete) =
        build_context(&histograms, 0, 5_000_000, 1_000).expect("valid edge context");

    assert!(complete);
    assert_eq!(context.as_tuple(), (0, 5_000_000));
    assert_eq!(combined.fragment_support, 1_000);
}

#[test]
fn incomplete_contig_is_marked_for_global_fallback() {
    let histograms = cores(4);

    let (combined, context, complete) =
        build_context(&histograms, 1, 5_000_000, 1_000).expect("valid short context");

    assert!(!complete);
    assert_eq!(context.as_tuple(), (0, 2_000_000));
    assert_eq!(combined.fragment_support, 400);
}

#[test]
fn fallback_need_checks_full_chromosome_span_and_summed_support() {
    let histograms = cores(4);
    assert!(!needs_global_fallback(&histograms, 2_000_000, 400).unwrap());
    assert!(needs_global_fallback(&histograms, 2_000_001, 400).unwrap());
    assert!(needs_global_fallback(&histograms, 2_000_000, 401).unwrap());
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
        assert_eq!(model.eligible_positions, 1_000);
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
        2_000_000,
        100,
    )
    .unwrap();

    // The chromosome is only 1.5 Mb long, so its sole eligible core needs the pooled model
    let models = &fitted.by_chromosome["chr1"];
    assert_eq!(models.len(), 3);
    assert!(models[0].is_none());
    assert!(models[2].is_none());
    let model = models[1].as_ref().unwrap();
    assert_eq!(model.core.as_tuple(), (500_000, 1_000_000));
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
        1_000_000,
        1_000,
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
            "GLOBAL\tNA\tNA\tNA\tNA\t200\t0\t1250000",
            "chr_local\t0\t500000\t0\t500000\t100\t0\t500000",
            "chr_local\t500000\t1000000\t500000\t1000000\t100\t0\t500000",
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
    assert!(histogram.eligible_positions().is_err());

    let mut support = CoreHistogram::new(interval);
    support.fragment_support = u64::MAX;
    let mut extra_support = CoreHistogram::new(interval);
    extra_support.fragment_support = 1;
    let error = support.add_histogram(&extra_support).unwrap_err();
    assert!(error.to_string().contains("fragment support overflow"));
}
