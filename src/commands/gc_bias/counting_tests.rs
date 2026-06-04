use super::{
    GCCounts, build_gc_prefixes, count_reference_gc_and_length_by_window,
    get_gc_integer_percentage_for_window,
};
use crate::shared::{
    blacklist::apply_blacklist_mask_to_seq,
    interval::{IndexedInterval, Interval},
};
use anyhow::Result;

#[test]
fn gc_prefix_helpers_return_prefix_differences_for_checked_intervals() -> Result<()> {
    // Sequence A C N G T A:
    // - [1,5) = C N G T -> GC=2, ACGT=3
    // - [2,4) = N G     -> GC=1, ACGT=1
    let seq = b"ACNGTA".to_vec();
    let prefixes = build_gc_prefixes(&seq);
    let full_interval = Interval::new(1usize, 5usize)?;
    let inner_interval = Interval::new(2usize, 4usize)?;

    assert_eq!(prefixes.gc_count(full_interval)?, 2);
    assert_eq!(prefixes.acgt_count(full_interval)?, 3);

    assert_eq!(prefixes.gc_count(inner_interval)?, 1);
    assert_eq!(prefixes.acgt_count(inner_interval)?, 1);
    Ok(())
}

#[test]
fn gc_prefix_helpers_error_when_interval_exceeds_prefix_bounds() -> Result<()> {
    // Sequence A C G T has prefix length 5, so the largest valid half-open interval end is 4.
    // Asking for [1,5) should therefore fail before any subtraction is attempted.
    let seq = b"ACGT".to_vec();
    let prefixes = build_gc_prefixes(&seq);
    let invalid_interval = Interval::new(1usize, 5usize)?;

    let gc_err = prefixes
        .gc_count(invalid_interval)
        .expect_err("expected GC prefix bounds error");
    let acgt_err = prefixes
        .acgt_count(invalid_interval)
        .expect_err("expected ACGT prefix bounds error");

    assert!(
        gc_err
            .to_string()
            .contains("GC interval [1, 5) out of bounds"),
        "unexpected GC error: {gc_err}"
    );
    assert!(
        acgt_err
            .to_string()
            .contains("ACGT interval [1, 5) out of bounds"),
        "unexpected ACGT error: {acgt_err}"
    );
    Ok(())
}

#[test]
fn gc_integer_percentage_window_returns_some_none_and_error() -> Result<()> {
    // Sequence A C N G T:
    // - [0,5) = A C N G T -> GC=2, ACGT=4, so GC% = round(200/4) = 50
    // - [1,4) = C N G     -> GC=2, ACGT=2, so this stays valid when min_acgt_count=2
    // - [2,4) = N G       -> GC=1, ACGT=1, so min_acgt_count=2 should return Ok(None)
    // - [3,6) extends past the prefix arrays for a 5 bp sequence, so it should error
    let seq = b"ACNGT".to_vec();
    let prefixes = build_gc_prefixes(&seq);

    let full_window = Interval::new(0usize, 5usize)?;
    let low_support_window = Interval::new(2usize, 4usize)?;
    let invalid_window = Interval::new(3usize, 6usize)?;

    assert_eq!(
        get_gc_integer_percentage_for_window(&prefixes, full_window, 0.0, 1)?,
        Some(50)
    );
    assert_eq!(
        get_gc_integer_percentage_for_window(&prefixes, low_support_window, 0.0, 2)?,
        None
    );

    let err = get_gc_integer_percentage_for_window(&prefixes, invalid_window, 0.0, 1)
        .expect_err("expected out-of-bounds GC window error");
    assert!(
        err.to_string().contains("GC interval [3, 6) out of bounds"),
        "unexpected GC window error: {err}"
    );
    Ok(())
}

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
    )?;

    // Assert:
    // Window 0 counts fragments starting at 0 inside `ACGTAC`:
    // - len4 trimmed to `CG`   -> gc=2
    // - len5 trimmed to `CGT`  -> gc=2
    // - len6 trimmed to `CGTA` -> gc=2
    //
    // Window 1 counts fragments starting at 6 inside `GTACGT`:
    // - len4 trimmed to `TA`   -> gc=0
    // - len5 trimmed to `TAC`  -> gc=1
    // - len6 trimmed to `TACG` -> gc=2
    //
    // No other `(length, gc)` cells should receive any counts.
    let expected_window0 = &[(4_usize, 2_usize), (5, 2), (6, 2)];
    let expected_window1 = &[(4_usize, 0_usize), (5, 1), (6, 2)];

    for (window_index, (window_counts, expected_non_zero_cells)) in counts_by_bin
        .iter()
        .zip([expected_window0, expected_window1].iter())
        .enumerate()
    {
        assert_eq!(
            window_counts.sum() as f64,
            expected_non_zero_cells.len() as f64,
            "window {window_index} should contain exactly one count for each tested length"
        );

        for length in 4..=6 {
            let effective_length = length - 2;
            for gc_count in 0..=effective_length {
                let expected_value = if expected_non_zero_cells.contains(&(length, gc_count)) {
                    1.0
                } else {
                    0.0
                };
                assert_eq!(
                    window_counts.get(length, gc_count).unwrap(),
                    expected_value,
                    "window {window_index} expected count {expected_value} at length {length}, gc {gc_count}"
                );
            }
        }
    }
    Ok(())
}

#[test]
fn reference_gc_counts_use_tile_local_prefix_coordinates_after_late_sequence_load() -> Result<()> {
    // This mirrors ref-gc-bias after loading a late reference slice, e.g. absolute [900,964).
    // The count helper receives tile-local coordinates, so absolute window [930,941) is [30,41)
    // and absolute start 930 is local start 30.
    //
    // In the ACGT repeat, local [30,41) is:
    //   G T A C G T A C G T A
    // This has 5 GC bases and 11 ACGT bases, so length 11 / GC count 5 gets one count.
    let seq = b"ACGT".repeat(16);
    let prefixes = build_gc_prefixes(&seq);
    let windows = IndexedInterval::from_tuples(&[(30, 41, 0)])?;
    let starts = vec![30usize];
    let mut counts_by_bin = vec![GCCounts::new(11, 11, 0, (0, 0))?];

    // Act
    count_reference_gc_and_length_by_window(
        &mut counts_by_bin,
        &prefixes,
        (11, 12),
        windows.as_slice(),
        starts.as_slice(),
        seq.len() as u64,
        1.0,
        1,
        0,
    )?;

    // Assert
    assert_eq!(counts_by_bin[0].sum(), 1.0);
    assert_eq!(
        counts_by_bin[0]
            .get(11, 5)
            .expect("length 11 / GC count 5 should be in range"),
        1.0
    );
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
    )?;

    // Assert: Masking drops the ACGT fraction below 1.0 so no counts are emitted
    assert_eq!(counts_by_bin[0].sum(), 0.0);
    Ok(())
}
