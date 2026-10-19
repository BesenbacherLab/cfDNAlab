use super::*;
use crate::commands::cli_common::{ChromosomeArgs, IOCArgs};
use crate::commands::outliers::{config::AmbiguousRegionAction, preservation::evaluate_regions};

fn config() -> OutliersConfig {
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
    config.set_blacklist_flank(Some(0));
    config
}

fn call(start: u32, end: u32) -> CalledRegion {
    CalledRegion {
        interval: Interval::new(start, end).unwrap(),
        observed_mass: f64::from(end - start) * 20.0,
        target_mass: f64::from(end - start),
        maximum_coverage: 20,
    }
}

#[test]
fn shared_windows_group_disconnected_calls_with_exact_counts_and_wider_context() {
    let config = config();
    let calls = [
        call(200_000, 203_000),
        call(210_000, 213_000),
        call(1_000_000, 1_006_000),
    ];
    let decisions = evaluate_regions(&config, &calls, &[], 2_000_000).unwrap();

    let groups =
        collect_preserved_regions("chr1", &calls, &decisions, &config, 2_000_000, 7).unwrap();

    // Two calls share a 6% window despite their 7 kb gap. Their reported span is 13 kb,
    // but only 6 kb is called. The 130 kb centered plot also expands to include that window
    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0].region_id, "region_000007");
    assert_eq!(groups[0].interval.as_tuple(), (200_000, 213_000));
    assert_eq!(groups[0].call_indices, vec![0, 1]);
    assert_eq!(groups[0].called_bases, 6_000);
    assert_eq!(groups[0].maximum_window_percentage, 6.0);
    assert_eq!(groups[0].window_span.as_tuple(), (200_000, 300_000));
    assert_eq!(groups[0].plot_interval.as_tuple(), (141_500, 300_000));
    assert_eq!(groups[1].region_id, "region_000008");
    assert_eq!(groups[1].call_indices, vec![2]);
    assert_eq!(groups[1].called_bases, 6_000);
    assert_eq!(groups[1].plot_interval.as_tuple(), (953_000, 1_100_000));
}

#[test]
fn overlapping_selected_windows_join_transitively() {
    let config = config();
    let calls = [
        call(100_000, 106_000),
        call(195_000, 201_000),
        call(290_000, 296_000),
    ];
    let decisions = evaluate_regions(&config, &calls, &[], 500_000).unwrap();

    let groups =
        collect_preserved_regions("chr1", &calls, &decisions, &config, 500_000, 1).unwrap();

    // The first selected window is [100k, 200k), with 11k called bases. The third call
    // selects [195k, 295k), also 11%, connecting all three originally separate merged regions
    assert_eq!(decisions.merged_regions.len(), 3);
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].call_indices, vec![0, 1, 2]);
    assert_eq!(groups[0].interval.as_tuple(), (100_000, 296_000));
    assert_eq!(groups[0].window_span.as_tuple(), (100_000, 295_000));
    assert_eq!(groups[0].called_bases, 18_000);
    assert_eq!(groups[0].maximum_window_percentage, 11.0);
}

#[test]
fn touching_evidence_windows_remain_separate_even_when_plot_contexts_overlap() {
    let config = config();
    let calls = [call(200_000, 206_000), call(300_000, 306_000)];
    let decisions = evaluate_regions(&config, &calls, &[], 1_000_000).unwrap();

    let groups =
        collect_preserved_regions("chr1", &calls, &decisions, &config, 1_000_000, 1).unwrap();

    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0].window_span.as_tuple(), (200_000, 300_000));
    assert_eq!(groups[1].window_span.as_tuple(), (300_000, 400_000));
    assert!(groups[0].plot_interval.intersects(groups[1].plot_interval));
}

#[test]
fn reporting_keeps_calls_preserved_through_flank_merging_in_their_original_group() {
    let mut config = config();
    config.set_blacklist_flank(Some(60_000));
    let calls = [
        call(10_000, 16_000),
        call(136_000, 136_030),
        call(256_030, 256_060),
    ];
    let decisions = evaluate_regions(&config, &calls, &[], 500_000).unwrap();

    let groups =
        collect_preserved_regions("chr1", &calls, &decisions, &config, 500_000, 1).unwrap();

    assert!(!decisions.calls[1].flagged);
    assert!(!decisions.calls[2].flagged);
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].call_indices, vec![0, 1, 2]);
    assert_eq!(groups[0].called_bases, 6_060);
    assert_eq!(groups[0].interval.as_tuple(), (10_000, 256_060));
    assert!(
        groups[0]
            .plot_interval
            .contains_interval(groups[0].window_span)
    );
}

#[test]
fn empty_and_explicitly_corrected_candidates_produce_no_preserved_groups() {
    let mut config = config();
    let empty = evaluate_regions(&config, &[], &[], 500_000).unwrap();
    assert!(
        collect_preserved_regions("chr1", &[], &empty, &config, 500_000, 1)
            .unwrap()
            .is_empty()
    );

    config.set_ambiguous_region_action(AmbiguousRegionAction::Correct);
    let calls = [call(100_000, 106_000)];
    let decisions = evaluate_regions(&config, &calls, &[], 500_000).unwrap();
    assert!(decisions.calls[0].flagged);
    assert!(
        collect_preserved_regions("chr1", &calls, &decisions, &config, 500_000, 1)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn plot_contexts_clip_chromosome_ends_and_respect_custom_settings() {
    let mut config = config();
    // A 40 kb group requests 400 kb around its midpoint. Clipping that centered span
    // must cover this complete 200 kb chromosome, even though the group is near its end
    assert_eq!(
        plot_context(
            Interval::new(150_000, 190_000).unwrap(),
            Interval::new(150_000, 200_000).unwrap(),
            &config,
            200_000
        )
        .unwrap()
        .as_tuple(),
        (0, 200_000)
    );
    assert_eq!(
        plot_context(
            Interval::new(0, 1_000).unwrap(),
            Interval::new(0, 100_000).unwrap(),
            &config,
            200_000
        )
        .unwrap()
        .as_tuple(),
        (0, 100_000)
    );
    assert_eq!(
        plot_context(
            Interval::new(196_000, 200_000).unwrap(),
            Interval::new(100_000, 200_000).unwrap(),
            &config,
            200_000
        )
        .unwrap()
        .as_tuple(),
        (100_000, 200_000)
    );

    config.set_preserved_plot_context_multiplier(2.5);
    config.set_preserved_plot_min_span(1);
    // A 10 kb group requests 25 kb centered coverage, then includes its 20 kb evidence window
    let context = plot_context(
        Interval::new(200_000, 210_000).unwrap(),
        Interval::new(200_000, 220_000).unwrap(),
        &config,
        500_000,
    )
    .unwrap();
    assert_eq!(context.as_tuple(), (192_500, 220_000));
}
