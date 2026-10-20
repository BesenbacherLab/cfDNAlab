use super::*;
use crate::commands::cli_common::{ChromosomeArgs, IOCArgs};

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
        observed_mass: (end - start) as f64 * 20.0,
        target_mass: (end - start) as f64,
        maximum_coverage: 20,
    }
}

#[test]
fn earlier_window_maximum_updates_later_calls_and_cutoff_equality_is_not_flagged() {
    let regions = [call(10_000, 12_500), call(20_000, 22_500)];
    let decisions = evaluate_regions(&config(), &regions, &[], 500_000).unwrap();
    // The first window contains 5,000 called bases in 100,000 eligible bases, exactly 5%
    for decision in &decisions.calls {
        assert_eq!(decision.window.score, 5.0);
        assert_eq!(decision.window.interval.start(), 10_000);
        assert!(!decision.flagged);
    }
}

#[test]
fn isolated_boundary_call_uses_supported_backward_context() {
    let decisions = evaluate_regions(&config(), &[call(196_000, 200_000)], &[], 200_000).unwrap();
    let window = decisions.calls[0].window;
    assert!(window.backward);
    assert_eq!(window.interval.as_tuple(), (100_000, 200_000));
    assert_eq!(window.called_positions, 4_000);
    assert_eq!(window.eligible_positions, 100_000);
    assert_eq!(window.score, 4.0);
    assert!(!decisions.merged_regions[0].preserved);
}

#[test]
fn supported_earlier_window_avoids_backward_fallback() {
    let regions = [call(110_000, 110_030), call(196_000, 200_000)];
    let decisions = evaluate_regions(&config(), &regions, &[], 200_000).unwrap();
    let window = decisions.calls[1].window;
    // The earlier window has 4,030 called bases in 90 kb. Prefer its supported 4.478%
    // over the later call's short window, whose 4,000 bases would give a floored 8%
    assert!(!window.backward);
    assert_eq!(window.interval.as_tuple(), (110_000, 200_000));
    assert_eq!(window.eligible_positions, 90_000);
    assert_eq!(window.called_positions, 4_030);
    assert!((window.score - 100.0 * 4_030.0 / 90_000.0).abs() < 1e-12);
    assert!(!decisions.calls[1].flagged);
}

#[test]
fn supported_current_window_beats_an_earlier_undersupported_maximum() {
    let regions = [call(100_000, 103_000), call(190_000, 190_100)];
    let blacklist = Interval::from_tuples(&[(110_000_u64, 185_000)]).unwrap();

    let decisions = evaluate_regions(&config(), &regions, &blacklist, 400_000).unwrap();

    // The first forward window has only 25 kb eligible and 3,100 called bases, giving
    // 6.2% with the 50 kb floor. Its originating call instead uses a supported backward
    // window containing 3,000 called bases in 100 kb, giving 3%
    assert!(decisions.calls[0].window.backward);
    assert_eq!(
        decisions.calls[0].window.interval.as_tuple(),
        (3_000, 103_000)
    );
    assert_eq!(decisions.calls[0].window.score, 3.0);
    // The later call's own supported window contains only its 100 called bases.
    // The earlier unsupported 6.2% window must not override this 0.1% evidence
    let window = decisions.calls[1].window;
    assert!(!window.backward);
    assert_eq!(window.interval.as_tuple(), (190_000, 290_000));
    assert_eq!(window.eligible_positions, 100_000);
    assert_eq!(window.called_positions, 100);
    assert!((window.score - 0.1).abs() < 1e-12);
    assert!(decisions.calls.iter().all(|decision| !decision.flagged));
    assert!(
        decisions
            .merged_regions
            .iter()
            .all(|region| !region.preserved)
    );
}

#[test]
fn a_window_at_minimum_support_participates_in_the_maximum() {
    let regions = [call(100_000, 103_000), call(190_000, 190_100)];
    let blacklist = Interval::from_tuples(&[(110_000_u64, 160_000)]).unwrap();

    let decisions = evaluate_regions(&config(), &regions, &blacklist, 400_000).unwrap();

    // Exactly 50 kb remains in the first window. It qualifies, so both calls use
    // its 3,100 / 50,000 = 6.2% rather than the later window's supported 0.1%
    for decision in &decisions.calls {
        assert_eq!(decision.window.interval.as_tuple(), (100_000, 200_000));
        assert_eq!(decision.window.eligible_positions, 50_000);
        assert_eq!(decision.window.called_positions, 3_100);
        assert!((decision.window.score - 6.2).abs() < 1e-12);
        assert!(!decision.window.backward);
        assert!(decision.flagged);
    }
}

#[test]
fn unsupported_windows_retain_their_maximum_when_no_supported_alternative_exists() {
    let regions = [call(100_000, 101_000), call(110_000, 112_000)];
    let blacklist = Interval::from_tuples(&[(0_u64, 90_000), (120_000, 200_000)]).unwrap();

    let decisions = evaluate_regions(&config(), &regions, &blacklist, 200_000).unwrap();

    // Only a 30 kb eligible island exists, so no forward or backward window can qualify.
    // The first forward window contains 3,000 called bases and scores 6% with the floor.
    // The second contains 2,000 and scores 4%. Both calls retain the first window's 6%
    for decision in &decisions.calls {
        assert_eq!(decision.window.interval.as_tuple(), (100_000, 200_000));
        assert_eq!(decision.window.eligible_positions, 20_000);
        assert_eq!(decision.window.called_positions, 3_000);
        assert_eq!(decision.window.score, 6.0);
        assert!(!decision.window.backward);
        assert!(decision.flagged);
    }
}

#[test]
fn a_supported_window_ending_at_the_call_start_is_not_an_alternative() {
    let regions = [call(10_000, 16_000), call(110_000, 111_000)];

    let decisions = evaluate_regions(&config(), &regions, &[], 130_000).unwrap();

    // The first forward window ends exactly at the second call's start. It cannot
    // prevent fallback. The second call's backward window contains 5,000 bases from
    // the first call and its own 1,000 bases in 100 kb, giving 6%
    let window = decisions.calls[1].window;
    assert!(window.backward);
    assert_eq!(window.interval.as_tuple(), (11_000, 111_000));
    assert_eq!(window.eligible_positions, 100_000);
    assert_eq!(window.called_positions, 6_000);
    assert_eq!(window.score, 6.0);
}

#[test]
fn masking_counts_sampling_zeros_and_never_collects_bases_across_a_large_gap() {
    let blacklist = Interval::from_tuples(&[(20_000_u64, 5_020_000)]).unwrap();
    let decisions =
        evaluate_regions(&config(), &[call(16_000, 20_000)], &blacklist, 6_000_000).unwrap();
    let window = decisions.calls[0].window;
    // Only 4,000 eligible positions remain forward. Backward also lacks 50 kb support,
    // so the 50 kb floor gives 8%, without reaching beyond the 5 Mb mask
    assert!(!window.backward);
    assert_eq!(window.interval.as_tuple(), (16_000, 116_000));
    assert_eq!(window.eligible_positions, 4_000);
    assert_eq!(window.called_positions, 4_000);
    assert_eq!(window.score, 8.0);
    assert!(decisions.merged_regions[0].preserved);
}

#[test]
fn preservation_propagates_through_all_touching_flanks() {
    let mut config = config();
    // Large configurable flanks connect three calls whose forward windows do not overlap.
    // Only the first call exceeds 5%, so the other two exercise shared preservation
    config.set_blacklist_flank(Some(60_000));
    let regions = [
        call(10_000, 16_000),
        call(136_000, 136_030),
        call(256_030, 256_060),
    ];
    let decisions = evaluate_regions(&config, &regions, &[], 500_000).unwrap();
    assert_eq!(decisions.merged_regions.len(), 1);
    assert_eq!(
        decisions.merged_regions[0].interval.as_tuple(),
        (0, 316_060)
    );
    assert!(decisions.calls[0].flagged);
    assert!(!decisions.calls[1].flagged);
    assert!(!decisions.calls[2].flagged);
    assert!(decisions.merged_regions[0].preserved);
    assert_eq!(decisions.calls[2].merged_region_index, 0);
    config.set_ambiguous_region_action(AmbiguousRegionAction::Correct);
    assert!(
        !evaluate_regions(&config, &regions, &[], 500_000)
            .unwrap()
            .merged_regions[0]
            .preserved
    );
    config.set_ambiguous_region_action(AmbiguousRegionAction::Error);
    let evidence = evaluate_regions(&config, &regions, &[], 500_000).unwrap();
    assert!(evidence.merged_regions[0].flagged);
    assert!(!evidence.merged_regions[0].preserved);
}

#[test]
fn empty_call_set_and_interval_prefix_boundaries_are_valid() {
    assert!(
        evaluate_regions(&config(), &[], &[], 200_000)
            .unwrap()
            .calls
            .is_empty()
    );
    let lengths =
        IntervalLengths::new(Interval::from_tuples(&[(10_u64, 20), (15, 30), (40, 50)]).unwrap());
    assert_eq!(lengths.within(Interval::new(20, 40).unwrap()), 10);
    assert_eq!(lengths.within(Interval::new(30, 40).unwrap()), 0);
    assert_eq!(lengths.within(Interval::new(0, 100).unwrap()), 30);
}
