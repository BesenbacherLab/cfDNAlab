use super::*;
use crate::commands::{
    cli_common::{ChromosomeArgs, IOCArgs},
    outliers::config::OutliersConfig,
};

#[test]
fn summaries_separate_exact_bases_flanked_spans_and_masked_bases() {
    let mut config = OutliersConfig::new(
        IOCArgs {
            bam: "unused.bam".into(),
            output_dir: "unused".into(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".into()]),
            chromosomes_file: None,
        },
    );
    config.set_blacklist_flank(Some(1_000));
    let regions = vec![CalledRegion {
        interval: Interval::new(10_000, 10_030).unwrap(),
        observed_mass: 600.0,
        target_mass: 30.0,
        maximum_coverage: 20,
    }];
    let blacklist = Interval::from_tuples(&[(9_000_u64, 9_500)]).unwrap();
    let decisions = evaluate_regions(&config, &regions, &blacklist, 200_000).unwrap();
    let summary = summarize_calls(
        &regions,
        &decisions,
        199_500,
        &IntervalLengths::new(blacklist),
        &[],
    );
    assert_eq!(summary.candidate_bases, 30);
    assert_eq!(summary.actionable_bases, 30);
    assert_eq!(summary.preserved_bases, 0);
    assert_eq!(summary.actionable_flanked_span, 2_030);
    assert_eq!(summary.actionable_flanked_eligible_bases, 1_530);
    assert_eq!(summary.percentage(30), 100.0 * 30.0 / 199_500.0);
    let mut total = CallSummary::default();
    total.add(&summary);
    assert_eq!(
        total.candidate_bases,
        summary.actionable_bases + summary.preserved_bases
    );
}

#[test]
fn no_calls_still_report_eligible_denominator_and_zero_percentages() {
    let summary = CallSummary {
        eligible_positions: 100_000,
        ..CallSummary::default()
    };
    let mut output = Vec::new();
    summary.write(&mut output, "chromosome", "chr1").unwrap();
    let text = String::from_utf8(output).unwrap();
    let fields: Vec<_> = text.trim_end().split('\t').collect();
    assert_eq!(fields.len(), 23);
    assert_eq!(fields[2], "100000");
    assert_eq!(fields[7], "0");
    assert_eq!(fields[9], "0");
    assert_eq!(fields[11], "0");
}
