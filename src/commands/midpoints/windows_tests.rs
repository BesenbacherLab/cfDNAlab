use super::{interval_with_margin_overlaps_blacklist, prepare_count_windows};
use crate::shared::{
    bam::Contigs,
    interval::{IndexedInterval, Interval},
};
use fxhash::FxHashMap;

fn contigs(entries: &[(&str, u32)]) -> Contigs {
    let mut contigs = Contigs {
        contigs: FxHashMap::default(),
    };
    for (tid, &(name, len)) in entries.iter().enumerate() {
        contigs.contigs.insert(name.to_string(), (tid as i32, len));
    }
    contigs
}

#[test]
fn prepare_count_windows_expands_by_smoothing_flank() {
    let mut windows = FxHashMap::default();
    windows.insert(
        "chr1".to_string(),
        vec![IndexedInterval::new(100, 110, 3_u64).expect("valid window")],
    );

    let (prepared, stats) = prepare_count_windows(
        windows,
        &contigs(&[("chr1", 500)]),
        &FxHashMap::default(),
        4,
        0,
        false,
    )
    .expect("expanded interval should fit chromosome bounds");

    let chr_windows = prepared.get("chr1").expect("chr1 should be retained");
    assert_eq!(chr_windows[0].as_tuple(), (96, 114, 3));
    assert_eq!(stats.loaded_after_chromosome_filtering, 1);
    assert_eq!(stats.retained_for_counting, 1);
}

#[test]
fn prepare_count_windows_drops_blacklist_margin_overlaps() {
    let mut windows = FxHashMap::default();
    windows.insert(
        "chr1".to_string(),
        vec![
            IndexedInterval::new(100, 110, 0_u64).expect("valid first window"),
            IndexedInterval::new(200, 210, 1_u64).expect("valid second window"),
        ],
    );
    let mut blacklist = FxHashMap::default();
    blacklist.insert(
        "chr1".to_string(),
        vec![Interval::new(85, 90).expect("valid blacklist")],
    );

    let (prepared, stats) =
        prepare_count_windows(windows, &contigs(&[("chr1", 500)]), &blacklist, 0, 15, true)
            .expect("blacklist prefiltering should succeed");

    let chr_windows = prepared.get("chr1").expect("chr1 should be retained");
    assert_eq!(chr_windows.len(), 1);
    assert_eq!(chr_windows[0].idx(), 1);
    assert_eq!(stats.dropped_by_blacklist_prefilter, 1);
    assert_eq!(stats.retained_for_counting, 1);
}

#[test]
fn prepare_count_windows_margin_can_include_fragment_radius_and_smoothing_flank() {
    let mut windows = FxHashMap::default();
    windows.insert(
        "chr1".to_string(),
        vec![IndexedInterval::new(100, 110, 0_u64).expect("valid window")],
    );
    let mut blacklist = FxHashMap::default();
    blacklist.insert(
        "chr1".to_string(),
        vec![Interval::new(91, 92).expect("valid blacklist")],
    );

    // This mirrors `ceil(max_fragment_length / 2) + smoothing_flank` with values 6 + 3.
    // The blacklist touches the margin-expanded interval [91,119), but it would not touch the
    // fragment-radius-only interval [94,116). That makes the smoothing-flank contribution visible.
    let (prepared, stats) =
        prepare_count_windows(windows, &contigs(&[("chr1", 500)]), &blacklist, 3, 9, true)
            .expect("blacklist prefiltering should succeed");

    assert!(prepared.get("chr1").expect("chr1 should exist").is_empty());
    assert_eq!(stats.dropped_by_blacklist_prefilter, 1);
    assert_eq!(stats.retained_for_counting, 0);
}

#[test]
fn keep_blacklisted_intervals_disables_interval_prefilter() {
    let mut windows = FxHashMap::default();
    windows.insert(
        "chr1".to_string(),
        vec![IndexedInterval::new(100, 110, 0_u64).expect("valid window")],
    );
    let mut blacklist = FxHashMap::default();
    blacklist.insert(
        "chr1".to_string(),
        vec![Interval::new(85, 90).expect("valid blacklist")],
    );

    let (prepared, stats) = prepare_count_windows(
        windows,
        &contigs(&[("chr1", 500)]),
        &blacklist,
        0,
        15,
        false,
    )
    .expect("disabled prefilter should keep the interval");

    assert_eq!(prepared.get("chr1").expect("chr1 should exist").len(), 1);
    assert_eq!(stats.dropped_by_blacklist_prefilter, 0);
    assert_eq!(stats.retained_for_counting, 1);
}

#[test]
fn prepare_count_windows_reports_raw_interval_outside_chromosome_bounds_without_smoothing_advice() {
    let mut windows = FxHashMap::default();
    windows.insert(
        "chr1".to_string(),
        vec![IndexedInterval::new(490, 510, 0_u64).expect("valid unchecked window")],
    );

    let error = prepare_count_windows(
        windows,
        &contigs(&[("chr1", 500)]),
        &FxHashMap::default(),
        0,
        0,
        false,
    )
    .expect_err("raw interval should fail when it extends past the chromosome end");
    let message = error.to_string();

    assert!(
        message.contains("Invalid midpoint interval chr1:490-510 extends beyond chromosome length 500"),
        "unexpected raw boundary error: {message}"
    );
    assert!(
        !message.contains("smooth"),
        "raw boundary error should not suggest smoothing changes: {message}"
    );
}

#[test]
fn prepare_count_windows_reports_outside_chromosome_before_blacklist_prefiltering() {
    let mut windows = FxHashMap::default();
    windows.insert(
        "chr1".to_string(),
        vec![IndexedInterval::new(490, 510, 0_u64).expect("valid unchecked window")],
    );
    let mut blacklist = FxHashMap::default();
    blacklist.insert(
        "chr1".to_string(),
        vec![Interval::new(489, 491).expect("valid blacklist")],
    );

    let error = prepare_count_windows(
        windows,
        &contigs(&[("chr1", 500)]),
        &blacklist,
        0,
        15,
        true,
    )
    .expect_err("out-of-chromosome intervals should fail before blacklist prefiltering");
    let message = error.to_string();

    assert!(
        message.contains("Invalid midpoint interval chr1:490-510 extends beyond chromosome length 500"),
        "unexpected out-of-chromosome error: {message}"
    );
}

#[test]
fn prepare_count_windows_requires_flank_to_fit_chromosome_bounds() {
    let mut windows = FxHashMap::default();
    windows.insert(
        "chr1".to_string(),
        vec![IndexedInterval::new(2, 12, 0_u64).expect("valid window")],
    );

    let error = prepare_count_windows(
        windows,
        &contigs(&[("chr1", 500)]),
        &FxHashMap::default(),
        4,
        0,
        false,
    )
    .expect_err("flanked interval should fail when it crosses the chromosome start");
    let message = error.to_string();

    assert!(
        message.contains("within 4 bp of the chromosome start"),
        "unexpected flank-boundary error: {message}"
    );
}

#[test]
fn interval_margin_overlap_uses_half_open_coordinates() {
    let blacklist = vec![Interval::new(90, 95).expect("valid blacklist")];
    let mut bl_ptr = 0usize;

    assert!(!interval_with_margin_overlaps_blacklist(
        Interval::new(100, 110).expect("valid interval"),
        5,
        &blacklist,
        &mut bl_ptr,
    ));

    let mut bl_ptr = 0usize;
    assert!(interval_with_margin_overlaps_blacklist(
        Interval::new(100, 110).expect("valid interval"),
        6,
        &blacklist,
        &mut bl_ptr,
    ));
}
