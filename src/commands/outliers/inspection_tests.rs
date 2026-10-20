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

#[test]
fn longest_preserved_is_a_single_exact_call_not_cluster_extent_or_total_bases() {
    let config = OutliersConfig::new(
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
    // A 651 bp call and nineteen 300 bp calls occupy 6,351 of the first call's
    // 100 kb window. All twenty share its 6.351% score and are preserved at 5%
    let mut regions = vec![CalledRegion {
        interval: Interval::new(10_000, 10_651).unwrap(),
        observed_mass: 13_020.0,
        target_mass: 651.0,
        maximum_coverage: 20,
    }];
    for call_index in 0..19_u32 {
        let start = 11_000 + call_index * 1_000;
        regions.push(CalledRegion {
            interval: Interval::new(start, start + 300).unwrap(),
            observed_mass: 6_000.0,
            target_mass: 300.0,
            maximum_coverage: 20,
        });
    }
    let decisions = evaluate_regions(&config, &regions, &[], 200_000).unwrap();

    let summary = summarize_calls(
        &regions,
        &decisions,
        200_000,
        &IntervalLengths::new(vec![]),
        &[],
    );

    assert_eq!(summary.preserved_calls, 20);
    assert_eq!(summary.preserved_bases, 6_351);
    assert_eq!(summary.longest_preserved, 651);
    assert_eq!(summary.preserved_merged_regions, 1);
    let mut total = CallSummary::default();
    total.add(&summary);
    total.add(&summary);
    assert_eq!(total.preserved_bases, 12_702);
    assert_eq!(total.longest_preserved, 651);
}

/// Build exact candidates with a known proposed weight of 1/20, without fitting a model.
fn audit_call(start: u32, end: u32) -> CalledRegion {
    CalledRegion {
        interval: Interval::new(start, end).unwrap(),
        observed_mass: f64::from(end - start) * 20.0,
        target_mass: f64::from(end - start),
        maximum_coverage: 20,
    }
}

#[test]
fn error_audit_collects_all_chromosomes_and_clips_only_the_contribution_count() {
    let mut config = OutliersConfig::new(
        IOCArgs {
            bam: "unused.bam".into(),
            output_dir: "unused".into(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr2".into(), "chr1".into()]),
            chromosomes_file: None,
        },
    );
    config.set_blacklist_flank(Some(0));
    config.set_ambiguous_region_action(AmbiguousRegionAction::Error);
    let chromosomes = vec!["chr2".into(), "chr1".into()];
    let contigs = Contigs {
        contigs: [("chr1".into(), (0, 500_000)), ("chr2".into(), (1, 500_000))]
            .into_iter()
            .collect(),
    };
    let metadata = OutputMetadata {
        config: &config,
        chromosomes: &chromosomes,
        contigs: &contigs,
        tail_probability: 1e-6,
        eligible_positions: 1_000_000,
        blacklist_flank: 0,
    };
    let calls = [
        (
            "chr1".into(),
            vec![
                audit_call(10_000, 16_000),
                audit_call(109_990, 110_030),
                audit_call(300_000, 300_030),
            ],
        ),
        ("chr2".into(), vec![audit_call(50_000, 56_000)]),
    ]
    .into_iter()
    .collect();

    let error = collect_decisions(&metadata, &calls, &FxHashMap::default()).unwrap_err();

    let failures = error.downcast_ref::<AnalysisFailures>().unwrap();
    assert_eq!(failures.messages.len(), 3);
    assert!(failures.messages[0].starts_with("chr2:50000-56000:"));
    assert!(failures.messages[1].starts_with("chr1:10000-16000:"));
    assert!(failures.messages[2].starts_with("chr1:109990-110030:"));
    assert!(failures.histogram.is_none());
    let rows: Vec<Vec<_>> = failures
        .audit
        .as_ref()
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 3);
    assert!(
        rows.iter()
            .all(|row| row.len() == 15 && row[13] == "withheld_error")
    );
    assert_eq!(rows[0][0], "chr2");
    // The first window includes 6,000 bp from the first call and only 10 bp from the second
    // Both use that 6.01% window, even though the complete second call spans 40 bp
    for row in &rows[1..] {
        assert_eq!(row[0], "chr1");
        assert_eq!(row[5], "true");
        assert_eq!(
            &row[6..11],
            &["10000", "110000", "forward", "100000", "6010"]
        );
        assert!((row[11].parse::<f64>().unwrap() - 6.01).abs() < 1e-12);
        assert!((row[12].parse::<f64>().unwrap() - 0.05).abs() < 1e-12);
        assert_eq!(row[14], "10000-16000,109990-110030");
    }
    assert_eq!(&rows[1][1..5], &["10000", "16000", "10000", "16000"]);
    assert_eq!(&rows[2][1..5], &["109990", "110030", "109990", "110030"]);
    assert!(!failures.audit.as_ref().unwrap().contains("300000"));
}

#[test]
fn error_audit_includes_unflagged_members_of_a_failed_transitive_flank_group() {
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
    config.set_blacklist_flank(Some(60_000));
    config.set_ambiguous_region_action(AmbiguousRegionAction::Error);
    let chromosomes = vec!["chr1".into()];
    let contigs = Contigs {
        contigs: [("chr1".into(), (0, 700_000))].into_iter().collect(),
    };
    let metadata = OutputMetadata {
        config: &config,
        chromosomes: &chromosomes,
        contigs: &contigs,
        tail_probability: 1e-6,
        eligible_positions: 700_000,
        blacklist_flank: 60_000,
    };
    let calls = [(
        "chr1".into(),
        vec![
            audit_call(10_000, 16_000),
            audit_call(136_000, 136_030),
            audit_call(256_030, 256_060),
            audit_call(450_000, 450_030),
        ],
    )]
    .into_iter()
    .collect();

    let error = collect_decisions(&metadata, &calls, &FxHashMap::default()).unwrap_err();

    let failures = error.downcast_ref::<AnalysisFailures>().unwrap();
    assert_eq!(failures.messages.len(), 1);
    assert!(failures.messages[0].contains("3 member calls"));
    let rows: Vec<Vec<_>> = failures
        .audit
        .as_ref()
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 3);
    // Flanks [0,76000), [76000,196030), [196030,316060) touch transitively
    // Only the first call has a 6% window. The two sparse neighbors each have a 0.03% window
    for row in &rows {
        assert_eq!(&row[1..3], &["0", "316060"]);
        assert_eq!(row[13], "withheld_error");
    }
    assert_eq!(
        rows.iter().map(|row| row[5]).collect::<Vec<_>>(),
        vec!["true", "false", "false"]
    );
    assert_eq!(
        rows.iter().map(|row| row[3]).collect::<Vec<_>>(),
        vec!["10000", "136000", "256030"]
    );
    assert!(!failures.audit.as_ref().unwrap().contains("450000"));
}
