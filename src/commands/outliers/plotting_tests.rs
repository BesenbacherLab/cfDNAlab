use super::*;

#[test]
fn summarizes_minimum_mean_and_maximum_used_thresholds() {
    let summary = summarize_thresholds(&[7.0, 10.0, 13.0]);

    assert_eq!(
        summary,
        ThresholdSummary {
            minimum: 7,
            mean_across_cores: 10.0,
            maximum: 13,
        }
    );
}
