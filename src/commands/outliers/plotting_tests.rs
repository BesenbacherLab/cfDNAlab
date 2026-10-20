use super::*;

#[test]
fn summarizes_minimum_mean_and_maximum_used_thresholds() {
    let summary = summarize_thresholds(&[Some(7.0), Some(10.0), Some(13.0)]).unwrap();

    assert_eq!(
        summary,
        ThresholdSummary {
            minimum: 7,
            mean_across_cores: 10.0,
            maximum: 13,
        }
    );
}

#[test]
fn threshold_summary_excludes_fully_blacklisted_grid_slots() {
    let summary = summarize_thresholds(&[None, Some(7.0), None, Some(13.0), None]).unwrap();

    // Only two actual cutoffs contribute. Empty slots are neither zeros nor extra observations
    assert_eq!(summary.minimum, 7);
    assert_eq!(summary.mean_across_cores, 10.0);
    assert_eq!(summary.maximum, 13);
}

#[test]
fn threshold_summary_requires_an_actual_fitted_model() {
    for thresholds in [&[][..], &[None, None][..]] {
        assert!(summarize_thresholds(thresholds).is_err());
    }
}

/// Assign fixed model values so diagnostic ordering is independent of numerical fitting.
fn diagnostic_model(start: u32, threshold: u32, mean: f64) -> CoreModel {
    use crate::commands::outliers::{
        config::OutlierModel,
        distribution::{CoverageParameters, diagnose_fit},
        striding::ModelSource,
    };
    use crate::shared::interval::Interval;
    let parameters = CoverageParameters {
        component_mean: mean,
        dispersion: 0.0,
        zero_inflation: 0.5,
    };
    let model = TwoStageCoverageModel {
        family: OutlierModel::Zip,
        initial: parameters,
        initial_threshold: threshold,
        underlying_fit: parameters,
        final_threshold: threshold,
        second_fit_retained_positions: 100,
        initial_optimizer: None,
        underlying_optimizer: None,
    };
    let core = Interval::new(start, start + 500_000).unwrap();
    CoreModel {
        core,
        context: core,
        eligible_positions: 100,
        fragment_support: 10,
        source: ModelSource::Local,
        model,
        diagnostics: diagnose_fit(&[(0, 50), (1, 50)].into_iter().collect(), model).unwrap(),
    }
}

#[test]
fn diagnostic_tracks_keep_blacklisted_gaps_and_chromosome_boundaries_at_original_indices() {
    let chromosomes = vec!["masked".to_string(), "chr1".to_string(), "chr2".to_string()];
    let models = [
        ("masked".to_string(), vec![None, None]),
        (
            "chr1".to_string(),
            vec![Some(diagnostic_model(0, 7, 2.0)), None],
        ),
        (
            "chr2".to_string(),
            vec![None, Some(diagnostic_model(500_000, 13, 4.0))],
        ),
    ]
    .into_iter()
    .collect();

    let diagnostics = collect_local_model_diagnostics(&chromosomes, &models).unwrap();

    assert_eq!(
        diagnostics.initial_thresholds,
        vec![None, None, Some(7.0), None, None, Some(13.0)]
    );
    assert_eq!(diagnostics.final_thresholds, diagnostics.initial_thresholds);
    assert_eq!(
        diagnostics.underlying_means,
        vec![None, None, Some(1.0), None, None, Some(2.0)]
    );
    assert_eq!(
        diagnostics.active_component_means,
        vec![None, None, Some(2.0), None, None, Some(4.0)]
    );
    assert_eq!(diagnostics.chromosome_boundaries, vec![2, 4]);
    assert_eq!(
        diagnostics.chromosome_labels,
        vec![(0, "masked".into()), (2, "chr1".into()), (4, "chr2".into())]
    );
    assert_eq!(diagnostics.final_summary.mean_across_cores, 10.0);
}
