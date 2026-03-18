use anyhow::Result;

use cfdnalab::{
    commands::gc_bias::counting::{
        GCCounts, build_gc_prefixes, count_reference_gc_and_length_by_window,
    },
    shared::{
        blacklist::apply_blacklist_mask_to_seq,
        interval::{IndexedInterval, Interval},
    },
};

#[test]
fn counts_gc_for_each_window_with_end_offset() -> Result<()> {
    // Arrange: Two windows of equal size with start positions seeded at each window start
    // End offset trims one base on each side so GC is counted on the inner span.
    let seq = b"ACGTACGTACGT".to_vec();
    let prefixes = build_gc_prefixes(&seq);
    let windows = IndexedInterval::from_tuples(&[(0, 6, 0), (6, 12, 1)])?;
    let starts = vec![0usize, 6usize];
    let mut counts_by_bin = vec![
        GCCounts::new(4, 6, 1, (0, 0))?,
        GCCounts::new(4, 6, 1, (0, 0))?,
    ];

    // Act
    count_reference_gc_and_length_by_window(
        &mut counts_by_bin,
        &prefixes,
        (4, 7),
        windows.as_slice(),
        starts.as_slice(),
        seq.len() as u64,
        1.0,
        1,
        1,
    );

    // Assert: Window 0 sees GC=2 for lengths 4, 5, and 6; window 1 sees an increasing GC profile
    let window0 = &counts_by_bin[0];
    assert_eq!(window0.get(4, 2).unwrap(), 1.0);
    assert_eq!(window0.get(5, 2).unwrap(), 1.0);
    assert_eq!(window0.get(6, 2).unwrap(), 1.0);

    let window1 = &counts_by_bin[1];
    assert_eq!(window1.get(4, 0).unwrap(), 1.0);
    assert_eq!(window1.get(5, 1).unwrap(), 1.0);
    assert_eq!(window1.get(6, 2).unwrap(), 1.0);
    Ok(())
}

#[test]
fn skips_counts_after_blacklist_removes_acgt_support() -> Result<()> {
    // Arrange: Blacklist the middle of the fragment so only half the bases remain ACGT
    let mut seq = b"ACGT".to_vec();
    let blacklist_intervals = Interval::from_tuples(&[(1, 3)])?;
    apply_blacklist_mask_to_seq(&mut seq, &blacklist_intervals, 0);
    let prefixes = build_gc_prefixes(&seq);
    let windows = IndexedInterval::from_tuples(&[(0, 4, 0)])?;
    let starts = vec![0usize];
    let mut counts_by_bin = vec![GCCounts::new(4, 4, 0, (0, 0))?];

    // Act
    count_reference_gc_and_length_by_window(
        &mut counts_by_bin,
        &prefixes,
        (4, 5),
        windows.as_slice(),
        starts.as_slice(),
        seq.len() as u64,
        1.0,
        1,
        0,
    );

    // Assert: Masking drops the ACGT fraction below 1.0 so no counts are emitted
    assert_eq!(counts_by_bin[0].sum(), 0.0);
    Ok(())
}
