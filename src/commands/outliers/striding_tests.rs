use super::*;

fn cores(number_of_cores: usize) -> Vec<CoreHistogram> {
    (0..number_of_cores)
        .map(|index| {
            let start = index as u32 * 500_000;
            let mut histogram = CoreHistogram::new(
                Interval::new(start, start + 500_000).expect("valid core"),
            );
            histogram.counts = [(0, 500), (1, 500)].into_iter().collect();
            histogram.fragment_support = 100;
            histogram
        })
        .collect()
}

fn example_global_model() -> TwoStageZipModel {
    let parameters = crate::commands::outliers::model::ZipParameters {
        lambda: 1.0,
        zero_inflation: 0.5,
    };
    TwoStageZipModel {
        initial: parameters,
        initial_threshold: 5,
        underlying_fit: parameters,
        final_threshold: 5,
        second_fit_retained_positions: 1_000,
    }
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
fn short_contig_assigns_the_global_fallback_model_to_every_core() {
    let chromosomes = vec!["chr1".to_string()];
    let histograms = cores(4);
    let histograms_by_chromosome = [("chr1".to_string(), histograms)].into_iter().collect();
    let global_model = example_global_model();

    let models = fit_core_models(
        &chromosomes,
        &histograms_by_chromosome,
        global_model,
        0.01,
        5_000_000,
        1_000,
    )
    .expect("global fallback models");

    let chromosome_models = models.get("chr1").expect("chromosome models");
    assert_eq!(chromosome_models.len(), 4);
    for model in chromosome_models {
        assert_eq!(model.source, ModelSource::GlobalFallback);
        assert_eq!(model.zip, global_model);
        assert_eq!(model.context.as_tuple(), (0, 2_000_000));
    }
}

#[test]
fn complete_degenerate_context_fails_instead_of_silently_falling_back() {
    let chromosomes = vec!["chr1".to_string()];
    let mut histogram =
        CoreHistogram::new(Interval::new(0, 500_000).expect("valid core interval"));
    histogram.counts.insert(0, 500_000);
    histogram.fragment_support = 1_000;
    let histograms_by_chromosome =
        [("chr1".to_string(), vec![histogram])].into_iter().collect();

    let error = fit_core_models(
        &chromosomes,
        &histograms_by_chromosome,
        example_global_model(),
        0.01,
        500_000,
        1_000,
    )
    .expect_err("a complete all-zero context must fail its local fit");

    assert!(error.to_string().contains("zero coverage"));
}

#[test]
fn core_histogram_excludes_masked_positions() {
    let mut histogram =
        CoreHistogram::new(Interval::new(0, 4).expect("valid histogram interval"));

    histogram
        .add_coverage(&[0.0, 1.0, 2.0, 2.0], Some(&[0, 1, 0, 0]))
        .expect("valid integer coverage");

    assert_eq!(histogram.counts, [(0, 1), (2, 2)].into_iter().collect());
    assert_eq!(histogram.eligible_positions().expect("eligible positions"), 3);
}

#[test]
fn rejects_coverage_at_the_f32_exact_integer_boundary() {
    let mut histogram =
        CoreHistogram::new(Interval::new(0, 1).expect("valid histogram interval"));

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

    let mut histogram =
        CoreHistogram::new(Interval::new(0, 2).expect("valid histogram interval"));
    let error = histogram
        .add_coverage(&[0.0, 1.0], Some(&[0]))
        .expect_err("mask and coverage lengths must agree");
    assert!(error.to_string().contains("blacklist-mask length"));
}

#[test]
fn global_histogram_conserves_counts_and_fragment_support_across_chromosomes() {
    let mut chr1_first =
        CoreHistogram::new(Interval::new(0, 100).expect("valid first core"));
    chr1_first.counts = [(0, 70), (1, 30)].into_iter().collect();
    chr1_first.fragment_support = 2;
    let mut chr1_second =
        CoreHistogram::new(Interval::new(100, 200).expect("valid second core"));
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
    assert_eq!(global.eligible_positions().expect("eligible positions"), 250);
    assert_eq!(global.fragment_support, 6);
}
