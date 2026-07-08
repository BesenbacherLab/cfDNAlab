use crate::{
    Result,
    shared::interval::{IndexedInterval, Interval},
};

const MIXED_COVERING_BROAD_NARROW_QUERY_STARTS: [u64; 8] = [
    1_000_000, 1_125_000, 1_250_000, 1_500_000, 1_875_000, 2_000_000, 2_250_000, 2_750_000,
];

// Shared mixed-size BED fixture for comparing the generic and tiered BED overlap finders.
fn mixed_covering_broad_narrow_windows() -> Result<Vec<IndexedInterval<u64>>> {
    IndexedInterval::from_tuples(&[
        (0, 4_000_000, 900_u64),
        (999_000, 2_001_000, 901_u64),
        (1_000_500, 1_000_800, 902_u64),
        (1_124_500, 1_125_500, 903_u64),
        (1_200_000, 1_300_000, 904_u64),
        (1_200_001, 1_300_000, 905_u64),
        (1_250_750, 1_350_750, 906_u64),
        (1_250_800, 1_350_799, 907_u64),
        (1_499_000, 1_501_000, 908_u64),
        (1_500_250, 1_500_750, 909_u64),
        (1_874_000, 2_001_000, 910_u64),
        (1_875_900, 1_975_900, 911_u64),
        (1_875_901, 1_975_900, 912_u64),
        (1_999_000, 3_001_000, 913_u64),
        (2_000_000, 2_001_000, 914_u64),
        (2_250_500, 2_250_700, 915_u64),
        (2_749_500, 2_750_500, 916_u64),
        (2_750_500, 2_850_500, 917_u64),
        (2_751_000, 2_751_200, 918_u64),
        (3_100_000, 3_200_000, 919_u64),
    ])
}

mod fixed_width_overlap_cursor_tests {
    use super::super::{FixedWidthOverlapCursor, OverlappingWindows, find_overlapping_windows};
    use super::{
        Interval, MIXED_COVERING_BROAD_NARROW_QUERY_STARTS, Result,
        mixed_covering_broad_narrow_windows,
    };
    use crate::shared::interval::IndexedInterval;
    use std::collections::BTreeMap;

    fn kmer_windows() -> Result<Vec<IndexedInterval<u64>>> {
        IndexedInterval::from_tuples(&[
            (2, 6, 10_u64),
            (6, 9, 11_u64),
            (8, 12, 12_u64),
        ])
    }

    fn overlap_signature(overlaps: Option<&OverlappingWindows>) -> Vec<(usize, u64, u64, f64)> {
        overlaps
            .map(|overlaps| {
                overlaps
                    .windows
                    .iter()
                    .map(|window| {
                        (
                            window.idx,
                            window.start(),
                            window.end(),
                            window.overlap_fraction,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default()
    }

    fn row_overlap_signature(
        windows: &[IndexedInterval<u64>],
        overlaps: Option<&OverlappingWindows>,
    ) -> Vec<(u64, u64, u64, f64)> {
        overlaps
            .map(|overlaps| {
                overlaps
                    .windows
                    .iter()
                    .map(|window| {
                        (
                            windows[window.idx].idx(),
                            window.start(),
                            window.end(),
                            window.overlap_fraction,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default()
    }

    fn assert_signature_close(
        observed: Option<&OverlappingWindows>,
        expected: &[(usize, u64, u64, f64)],
    ) {
        let observed = overlap_signature(observed);
        assert_eq!(observed.len(), expected.len());
        for (observed_window, expected_window) in observed.iter().zip(expected.iter()) {
            assert_eq!(observed_window.0, expected_window.0);
            assert_eq!(observed_window.1, expected_window.1);
            assert_eq!(observed_window.2, expected_window.2);
            assert!(
                (observed_window.3 - expected_window.3).abs() < 1e-12,
                "observed fraction {} != expected fraction {}",
                observed_window.3,
                expected_window.3
            );
        }
    }

    fn assert_signature_close_unordered(
        observed: Option<&OverlappingWindows>,
        expected: &[(usize, u64, u64, f64)],
    ) {
        let mut observed = overlap_signature(observed);
        let mut expected = expected.to_vec();
        observed.sort_by_key(|window| (window.0, window.1, window.2));
        expected.sort_by_key(|window| (window.0, window.1, window.2));
        assert_eq!(observed.len(), expected.len());
        for (observed_window, expected_window) in observed.iter().zip(expected.iter()) {
            assert_eq!(observed_window.0, expected_window.0);
            assert_eq!(observed_window.1, expected_window.1);
            assert_eq!(observed_window.2, expected_window.2);
            assert!(
                (observed_window.3 - expected_window.3).abs() < 1e-12,
                "observed fraction {} != expected fraction {}",
                observed_window.3,
                expected_window.3
            );
        }
    }

    fn assert_row_signature_close(
        windows: &[IndexedInterval<u64>],
        observed: Option<&OverlappingWindows>,
        expected: &[(u64, u64, u64, f64)],
    ) {
        let observed = row_overlap_signature(windows, observed);
        assert_eq!(observed.len(), expected.len());
        for (observed_window, expected_window) in observed.iter().zip(expected.iter()) {
            assert_eq!(observed_window.0, expected_window.0);
            assert_eq!(observed_window.1, expected_window.1);
            assert_eq!(observed_window.2, expected_window.2);
            assert!(
                (observed_window.3 - expected_window.3).abs() < 1e-12,
                "observed fraction {} != expected fraction {}",
                observed_window.3,
                expected_window.3
            );
        }
    }

    fn assert_same_overlaps(
        cached: Option<&OverlappingWindows>,
        baseline: Option<&OverlappingWindows>,
    ) {
        let cached = overlap_signature(cached);
        let baseline = overlap_signature(baseline);
        assert_eq!(cached.len(), baseline.len());
        for (cached_window, baseline_window) in cached.iter().zip(baseline.iter()) {
            assert_eq!(cached_window.0, baseline_window.0);
            assert_eq!(cached_window.1, baseline_window.1);
            assert_eq!(cached_window.2, baseline_window.2);
            assert!(
                (cached_window.3 - baseline_window.3).abs() < 1e-12,
                "cached fraction {} != baseline fraction {}",
                cached_window.3,
                baseline_window.3
            );
        }
    }

    fn add_count_overlap_weights(
        counts: &mut BTreeMap<u64, f64>,
        windows: &[IndexedInterval<u64>],
        overlaps: Option<&OverlappingWindows>,
    ) {
        if let Some(overlaps) = overlaps {
            for overlapped_window in &overlaps.windows {
                let row_idx = windows[overlapped_window.idx].idx();
                *counts.entry(row_idx).or_default() += overlapped_window.overlap_fraction;
            }
        }
    }

    fn assert_count_close(counts: &BTreeMap<u64, f64>, row_idx: u64, expected: f64) {
        let observed = counts.get(&row_idx).copied().unwrap_or_default();
        assert!(
            (observed - expected).abs() < 1e-12,
            "row {row_idx}: observed {observed}, expected {expected}"
        );
    }

    fn assert_counts_close(counts: &BTreeMap<u64, f64>, expected: &[(u64, f64)]) {
        let expected_by_row: BTreeMap<u64, f64> = expected.iter().copied().collect();
        for row_idx in counts.keys() {
            assert!(
                expected_by_row.contains_key(row_idx),
                "unexpected observed row {row_idx}"
            );
        }
        for &(row_idx, expected_count) in expected {
            assert_count_close(counts, row_idx, expected_count);
        }
    }

    fn assert_overlap_fraction_error(error: crate::Error, expected_fraction: f64) {
        match error {
            crate::Error::OverlapFractionOutOfBounds { overlap_fraction } => {
                if expected_fraction.is_nan() {
                    assert!(overlap_fraction.is_nan());
                } else {
                    assert_eq!(overlap_fraction, expected_fraction);
                }
            }
            other => panic!("unexpected error: {other}"),
        }
    }

    #[test]
    fn finder_accumulates_mixed_covering_broad_and_narrow_window_hits() -> Result<()> {
        // Eight query intervals stand in for 1000 bp paired-fragment spans:
        // [1000000,1001000), [1125000,1126000), [1250000,1251000),
        // [1500000,1501000), [1875000,1876000), [2000000,2001000),
        // [2250000,2251000), and [2750000,2751000).
        //
        // The expected counts below are hand-derived as overlap_bp / 1000. They intentionally
        // include:
        // - covering rows that span multiple query intervals and tile+reach regions
        // - broad rows with length exactly 100 kb and larger
        // - narrow rows, including a 99,999 bp row
        // - rows that touch a query boundary or sit after all queries and must remain zero
        let chrom_len = 4_000_000;
        let windows = mixed_covering_broad_narrow_windows()?;
        let mut window_ptr = 0;
        let mut counts = BTreeMap::new();

        for query_start in MIXED_COVERING_BROAD_NARROW_QUERY_STARTS {
            let observed = find_overlapping_windows(
                chrom_len,
                &mut window_ptr,
                Some(windows.as_slice()),
                None,
                Interval::new(query_start, query_start + 1_000)?,
                0.0,
                1_000,
            )?;
            add_count_overlap_weights(&mut counts, windows.as_slice(), observed.as_ref());
        }

        assert_counts_close(
            &counts,
            &[
                (900, 8.0),
                (901, 6.0),
                (902, 0.3),
                (903, 0.5),
                (904, 1.0),
                (905, 1.0),
                (906, 0.25),
                (907, 0.2),
                (908, 1.0),
                (909, 0.5),
                (910, 2.0),
                (911, 0.1),
                (912, 0.099),
                (913, 3.0),
                (914, 1.0),
                (915, 0.2),
                (916, 0.5),
                (917, 0.5),
                (918, 0.0),
                (919, 0.0),
            ],
        );
        assert!((counts.values().copied().sum::<f64>() - 26.149).abs() < 1e-12);

        Ok(())
    }

    #[test]
    fn finder_keeps_covering_broad_window_available_for_later_queries() -> Result<()> {
        // The 2 Mb [0,2000000) row covers both queries. Later 200 kb, 2 kb, and 200 bp rows must
        // not make the streaming pointer forget that covering row while moving between queries.
        let chrom_len = 3_000_000;
        let windows = IndexedInterval::from_tuples(&[
            (0, 2_000_000, 500_u64),
            (900_000, 1_100_000, 501_u64),
            (999_500, 1_001_500, 502_u64),
            (1_000_100, 1_000_300, 503_u64),
            (1_400_000, 1_600_000, 504_u64),
            (1_499_500, 1_501_500, 505_u64),
            (1_500_100, 1_500_300, 506_u64),
        ])?;
        let mut window_ptr = 0;

        let first_query = find_overlapping_windows(
            chrom_len,
            &mut window_ptr,
            Some(windows.as_slice()),
            None,
            Interval::new(1_000_000, 1_001_000)?,
            0.0,
            0,
        )?;
        let second_query = find_overlapping_windows(
            chrom_len,
            &mut window_ptr,
            Some(windows.as_slice()),
            None,
            Interval::new(1_500_000, 1_501_000)?,
            0.0,
            0,
        )?;

        assert_signature_close_unordered(
            first_query.as_ref(),
            &[
                (0, 0, 2_000_000, 1.0),
                (1, 900_000, 1_100_000, 1.0),
                (2, 999_500, 1_001_500, 1.0),
                (3, 1_000_100, 1_000_300, 0.2),
            ],
        );
        assert_signature_close_unordered(
            second_query.as_ref(),
            &[
                (0, 0, 2_000_000, 1.0),
                (4, 1_400_000, 1_600_000, 1.0),
                (5, 1_499_500, 1_501_500, 1.0),
                (6, 1_500_100, 1_500_300, 0.2),
            ],
        );

        Ok(())
    }

    #[test]
    fn fixed_width_cursor_returns_manual_bed_overlap_signatures() -> Result<()> {
        // Each expected fraction is derived from a four-base query divided by the overlap length
        // with the BED window. Row ids come from `IndexedInterval.idx()`, not scan positions.
        let chrom_len = 14;
        let kmer_size = 4;
        let min_overlap_fraction = 1.0 / (kmer_size as f64 + 1.0);
        let windows = kmer_windows()?;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(windows.as_slice()),
            None,
            kmer_size,
            min_overlap_fraction,
            0,
        )?;
        let expected_by_start: [Vec<(u64, u64, u64, f64)>; 11] = [
            vec![(10, 2, 6, 2.0 / 4.0)],
            vec![(10, 2, 6, 3.0 / 4.0)],
            vec![(10, 2, 6, 4.0 / 4.0)],
            vec![(10, 2, 6, 3.0 / 4.0), (11, 6, 9, 1.0 / 4.0)],
            vec![(10, 2, 6, 2.0 / 4.0), (11, 6, 9, 2.0 / 4.0)],
            vec![
                (10, 2, 6, 1.0 / 4.0),
                (11, 6, 9, 3.0 / 4.0),
                (12, 8, 12, 1.0 / 4.0),
            ],
            vec![(11, 6, 9, 3.0 / 4.0), (12, 8, 12, 2.0 / 4.0)],
            vec![(11, 6, 9, 2.0 / 4.0), (12, 8, 12, 3.0 / 4.0)],
            vec![(11, 6, 9, 1.0 / 4.0), (12, 8, 12, 4.0 / 4.0)],
            vec![(12, 8, 12, 3.0 / 4.0)],
            vec![(12, 8, 12, 2.0 / 4.0)],
        ];
        let mut counts = BTreeMap::new();

        for (query_start, expected_windows) in expected_by_start.iter().enumerate() {
            let observed = cursor.find_overlaps(query_start as u64)?;

            assert_row_signature_close(windows.as_slice(), observed.as_ref(), expected_windows);
            add_count_overlap_weights(&mut counts, windows.as_slice(), observed.as_ref());
        }

        assert_count_close(&counts, 10, 3.75);
        assert_count_close(&counts, 11, 3.00);
        assert_count_close(&counts, 12, 3.75);
        assert_eq!(counts.len(), 3);

        Ok(())
    }

    #[test]
    fn threshold_0_51_requires_three_of_four_bases() -> Result<()> {
        // With k = 4, threshold 0.51 requires at least ceil(2.04) = 3 shared bases.
        let chrom_len = 8;
        let kmer_size = 4;
        let windows = IndexedInterval::from_tuples(&[(2, 6, 10_u64)])?;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(windows.as_slice()),
            None,
            kmer_size,
            0.51,
            0,
        )?;
        let expected_by_start: [Vec<(u64, u64, u64, f64)>; 5] = [
            vec![],
            vec![(10, 2, 6, 3.0 / 4.0)],
            vec![(10, 2, 6, 4.0 / 4.0)],
            vec![(10, 2, 6, 3.0 / 4.0)],
            vec![],
        ];

        for (query_start, expected_windows) in expected_by_start.iter().enumerate() {
            let observed = cursor.find_overlaps(query_start as u64)?;

            assert_row_signature_close(windows.as_slice(), observed.as_ref(), expected_windows);
        }

        Ok(())
    }

    #[test]
    fn zero_threshold_still_requires_positive_overlap() -> Result<()> {
        // Touching half-open intervals share zero bases and must not be returned.
        let chrom_len = 12;
        let kmer_size = 4;
        let windows = IndexedInterval::from_tuples(&[(4, 8, 10_u64)])?;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(windows.as_slice()),
            None,
            kmer_size,
            0.0,
            0,
        )?;

        assert!(cursor.find_overlaps(0)?.is_none());
        let one_base_overlap = cursor.find_overlaps(1)?;
        assert_row_signature_close(
            windows.as_slice(),
            one_base_overlap.as_ref(),
            &[(10, 4, 8, 1.0 / 4.0)],
        );
        assert!(cursor.find_overlaps(8)?.is_none());

        Ok(())
    }

    #[test]
    fn bed_window_end_is_clipped_to_chromosome() -> Result<()> {
        // The first BED row extends past the chromosome and is counted as [8, 10). The second row
        // starts at the chromosome end and never has a positive-width clipped span.
        let chrom_len = 10;
        let kmer_size = 4;
        let windows = IndexedInterval::from_tuples(&[(8, 15, 10_u64), (10, 12, 11_u64)])?;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(windows.as_slice()),
            None,
            kmer_size,
            0.0,
            0,
        )?;

        let one_base_overlap = cursor.find_overlaps(5)?;
        assert_row_signature_close(
            windows.as_slice(),
            one_base_overlap.as_ref(),
            &[(10, 8, 10, 1.0 / 4.0)],
        );
        let two_base_overlap = cursor.find_overlaps(6)?;
        assert_row_signature_close(
            windows.as_slice(),
            two_base_overlap.as_ref(),
            &[(10, 8, 10, 2.0 / 4.0)],
        );
        let clipped_query_overlap = cursor.find_overlaps(7)?;
        assert_row_signature_close(
            windows.as_slice(),
            clipped_query_overlap.as_ref(),
            &[(10, 8, 10, 2.0 / 3.0)],
        );

        Ok(())
    }

    #[test]
    fn clipped_query_uses_clipped_denominator() -> Result<()> {
        // Query start 8 with k = 4 is clipped to [8, 10), so fractions use length 2.
        let chrom_len = 10;
        let kmer_size = 4;
        let windows = IndexedInterval::from_tuples(&[(8, 10, 10_u64), (9, 10, 11_u64)])?;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(windows.as_slice()),
            None,
            kmer_size,
            0.0,
            0,
        )?;

        let observed = cursor.find_overlaps(8)?;
        assert_row_signature_close(
            windows.as_slice(),
            observed.as_ref(),
            &[(10, 8, 10, 1.0), (11, 9, 10, 0.5)],
        );

        Ok(())
    }

    #[test]
    fn nonzero_starting_pointer_keeps_original_scan_indices() -> Result<()> {
        // Starting at pointer 1 skips the first BED row, but returned indices still refer to the
        // original scan positions in the full slice.
        let chrom_len = 10;
        let kmer_size = 4;
        let windows =
            IndexedInterval::from_tuples(&[(0, 2, 20_u64), (2, 6, 21_u64), (5, 9, 22_u64)])?;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(windows.as_slice()),
            None,
            kmer_size,
            0.0,
            1,
        )?;

        let observed = cursor.find_overlaps(3)?;
        assert_signature_close(
            observed.as_ref(),
            &[(1, 2, 6, 3.0 / 4.0), (2, 5, 9, 2.0 / 4.0)],
        );

        Ok(())
    }

    #[test]
    fn fixed_size_clipped_query_uses_clipped_denominator() -> Result<()> {
        // Query start 7 with k = 4 is clipped to [7, 10), so fractions use length 3.
        let mut cursor = FixedWidthOverlapCursor::new(10, None, Some(4), 4, 0.0, 0)?;

        let observed = cursor.find_overlaps(7)?;
        assert_signature_close(
            observed.as_ref(),
            &[(1, 4, 8, 1.0 / 3.0), (2, 8, 10, 2.0 / 3.0)],
        );

        Ok(())
    }

    #[test]
    fn fixed_width_cursor_matches_finder_and_count_overlap_totals_in_kmer_loop() -> Result<()> {
        // Windows:
        //   row 10: [2, 6)
        //   row 11: [6, 9)
        //   row 12: [8, 12)
        //
        // Four-base k-mers scan starts 0..=10 on chr length 14.
        // Hand-derived count-overlap sums:
        //   row 10: 2/4 + 3/4 + 4/4 + 3/4 + 2/4 + 1/4 = 3.75
        //   row 11: 1/4 + 2/4 + 3/4 + 3/4 + 2/4 + 1/4 = 3.00
        //   row 12: 1/4 + 2/4 + 3/4 + 4/4 + 3/4 + 2/4 = 3.75
        let chrom_len = 14;
        let kmer_size = 4;
        let min_overlap_fraction = 1.0 / (kmer_size as f64 + 1.0);
        let windows = kmer_windows()?;
        let mut baseline_wd_ptr = 0;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(windows.as_slice()),
            None,
            kmer_size,
            min_overlap_fraction,
            0,
        )?;
        let mut counts = BTreeMap::new();

        for kmer_start in 0..=chrom_len - kmer_size {
            let query_interval = Interval::new(kmer_start, kmer_start + kmer_size)?;
            let baseline = find_overlapping_windows(
                chrom_len,
                &mut baseline_wd_ptr,
                Some(windows.as_slice()),
                None,
                query_interval,
                min_overlap_fraction,
                0,
            )?;
            let cached = cursor.find_overlaps(kmer_start)?;

            assert_same_overlaps(cached.as_ref(), baseline.as_ref());
            add_count_overlap_weights(&mut counts, windows.as_slice(), cached.as_ref());
        }

        assert_count_close(&counts, 10, 3.75);
        assert_count_close(&counts, 11, 3.00);
        assert_count_close(&counts, 12, 3.75);
        assert_eq!(counts.len(), 3);
        // The cursor should reuse cached BED candidates instead of rebuilding for every start
        assert!(cursor.refresh_count() < (chrom_len - kmer_size + 1) as usize);

        Ok(())
    }

    #[test]
    fn fixed_width_cursor_aggregates_grouped_windows_by_stored_group_idx() -> Result<()> {
        // Grouped BED windows carry `group_idx` in `IndexedInterval.idx()`. The overlap finder
        // returns scan positions, so grouped counting maps `overlap.idx` back through
        // `windows[overlap.idx].idx()` before incrementing a group row.
        //
        // Group 20 has two separate intervals:
        //   [2, 6): 2/4 + 3/4 + 4/4 + 3/4 + 2/4 + 1/4 = 3.75
        //   [8, 12): 1/4 + 2/4 + 3/4 + 4/4 + 3/4 + 2/4 = 3.75
        //
        // Group 21 has one interval:
        //   [6, 9): 1/4 + 2/4 + 3/4 + 3/4 + 2/4 + 1/4 = 3.00
        let chrom_len = 14;
        let kmer_size = 4;
        let min_overlap_fraction = 1.0 / (kmer_size as f64 + 1.0);
        let grouped_windows = IndexedInterval::from_tuples(&[
            (2, 6, 20_u64),
            (6, 9, 21_u64),
            (8, 12, 20_u64),
        ])?;
        let mut baseline_wd_ptr = 0;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(grouped_windows.as_slice()),
            None,
            kmer_size,
            min_overlap_fraction,
            0,
        )?;
        let mut group_counts = BTreeMap::new();

        for kmer_start in 0..=chrom_len - kmer_size {
            let query_interval = Interval::new(kmer_start, kmer_start + kmer_size)?;
            let baseline = find_overlapping_windows(
                chrom_len,
                &mut baseline_wd_ptr,
                Some(grouped_windows.as_slice()),
                None,
                query_interval,
                min_overlap_fraction,
                0,
            )?;
            let cached = cursor.find_overlaps(kmer_start)?;

            assert_same_overlaps(cached.as_ref(), baseline.as_ref());
            add_count_overlap_weights(
                &mut group_counts,
                grouped_windows.as_slice(),
                cached.as_ref(),
            );
        }

        assert_count_close(&group_counts, 20, 7.50);
        assert_count_close(&group_counts, 21, 3.00);
        assert_eq!(group_counts.len(), 2);

        Ok(())
    }

    #[test]
    fn fixed_width_cursor_applies_threshold_after_cache_reuse() -> Result<()> {
        // With k = 4 and threshold 0.75, cached candidate windows cannot be returned blindly.
        // The same candidate window must be re-evaluated at each start because its overlap
        // fraction can cross the threshold before that window enters or exits the candidate slice.
        let chrom_len = 14;
        let kmer_size = 4;
        let windows = kmer_windows()?;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(windows.as_slice()),
            None,
            kmer_size,
            0.75,
            0,
        )?;
        let expected_hits_by_start: [Vec<usize>; 11] = [
            vec![],
            vec![0],
            vec![0],
            vec![0],
            vec![],
            vec![1],
            vec![1],
            vec![2],
            vec![2],
            vec![2],
            vec![],
        ];

        for (kmer_start, expected_hits) in expected_hits_by_start.iter().enumerate() {
            let observed_hits: Vec<usize> = cursor
                .find_overlaps(kmer_start as u64)?
                .map(|overlaps| {
                    overlaps
                        .windows
                        .iter()
                        .map(|window| window.idx)
                        .collect()
                })
                .unwrap_or_default();

            assert_eq!(observed_hits.as_slice(), expected_hits.as_slice());
        }

        Ok(())
    }

    #[test]
    fn fixed_width_cursor_excludes_short_windows_for_full_overlap_threshold() -> Result<()> {
        // The middle window is [6, 9), so its maximum overlap with a 4 bp k-mer is 3/4.
        // A full-overlap threshold must therefore never emit it, even though it has a
        // constant-overlap span where the query fully covers the shorter window.
        let chrom_len = 14;
        let kmer_size = 4;
        let windows = kmer_windows()?;
        let mut cursor = FixedWidthOverlapCursor::new(
            chrom_len,
            Some(windows.as_slice()),
            None,
            kmer_size,
            1.0,
            0,
        )?;
        let expected_hits_by_start: [Vec<usize>; 11] = [
            vec![],
            vec![],
            vec![0],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![2],
            vec![],
            vec![],
        ];

        for (kmer_start, expected_hits) in expected_hits_by_start.iter().enumerate() {
            let observed_hits: Vec<usize> = cursor
                .find_overlaps(kmer_start as u64)?
                .map(|overlaps| {
                    overlaps
                        .windows
                        .iter()
                        .map(|window| window.idx)
                        .collect()
                })
                .unwrap_or_default();

            assert_eq!(observed_hits.as_slice(), expected_hits.as_slice());
        }

        Ok(())
    }

    #[test]
    fn fixed_width_cursor_handles_size_global_and_empty_bed_edges() -> Result<()> {
        let mut by_size = FixedWidthOverlapCursor::new(10, None, Some(5), 4, 0.0, 0)?;
        let split_window = by_size
            .find_overlaps(3)?
            .expect("query [3, 7) should hit both fixed windows");
        let split_signature = overlap_signature(Some(&split_window));
        assert_eq!(split_signature.len(), 2);
        assert_eq!(split_signature[0], (0, 0, 5, 0.5));
        assert_eq!(split_signature[1], (1, 5, 10, 0.5));

        let right_window = by_size
            .find_overlaps(6)?
            .expect("query [6, 10) should hit the right fixed window");
        assert_eq!(overlap_signature(Some(&right_window)), vec![(1, 5, 10, 1.0)]);
        assert!(by_size.find_overlaps(10)?.is_none());

        let mut global = FixedWidthOverlapCursor::new(10, None, None, 4, 0.0, 0)?;
        let clipped_global = global
            .find_overlaps(8)?
            .expect("query [8, 10) should hit the chromosome-wide window");
        assert_eq!(
            overlap_signature(Some(&clipped_global)),
            vec![(0, 0, 10, 1.0)]
        );

        let empty_windows = Vec::new();
        let mut empty_bed = FixedWidthOverlapCursor::new(
            10,
            Some(empty_windows.as_slice()),
            None,
            4,
            0.0,
            0,
        )?;
        assert!(empty_bed.find_overlaps(0)?.is_none());
        assert_eq!(empty_bed.refresh_count(), 1);

        Ok(())
    }

    #[test]
    fn overlap_finder_rejects_invalid_min_overlap_fractions() -> Result<()> {
        let chrom_len = 14;
        let kmer_size = 4;
        let windows = kmer_windows()?;

        for threshold in [-0.1, f64::NAN, 1.1] {
            let mut wd_ptr = 0;
            let error = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                Some(windows.as_slice()),
                None,
                Interval::new(2, 2 + kmer_size)?,
                threshold,
                0,
            )
            .expect_err("invalid overlap thresholds should be rejected");

            assert_overlap_fraction_error(error, threshold);
        }

        Ok(())
    }

    #[test]
    fn fixed_width_cursor_rejects_invalid_min_overlap_fractions() -> Result<()> {
        let chrom_len = 14;
        let kmer_size = 4;
        let windows = kmer_windows()?;

        for threshold in [-0.1, f64::NAN, 1.1] {
            let error = FixedWidthOverlapCursor::new(
                chrom_len,
                Some(windows.as_slice()),
                None,
                kmer_size,
                threshold,
                0,
            )
            .expect_err("invalid overlap thresholds should be rejected");

            assert_overlap_fraction_error(error, threshold);
        }

        Ok(())
    }
}

mod tile_bed_overlap_context_tests {
    use super::super::{
        ChromosomeBedWindows, OverlappingWindows, TileBedOverlapContext, TileBedWindowSpans,
        find_overlapping_windows,
    };
    use super::{
        MIXED_COVERING_BROAD_NARROW_QUERY_STARTS, mixed_covering_broad_narrow_windows,
    };
    use crate::{
        Result,
        shared::interval::{IndexedInterval, Interval},
        shared::tiled_run::TileWindowSpan,
    };

    fn full_span<T>(windows: &[T]) -> Option<TileWindowSpan> {
        (!windows.is_empty()).then_some(TileWindowSpan {
            first_idx: 0,
            last_idx_exclusive: windows.len(),
        })
    }

    fn full_bed_spans(chromosome_windows: &ChromosomeBedWindows) -> TileBedWindowSpans {
        TileBedWindowSpans {
            all_windows_span: full_span(chromosome_windows.all_windows.as_slice()),
            tier_spans: chromosome_windows
                .tiers
                .iter()
                .map(|tier| full_span(tier.windows.as_slice()))
                .collect(),
        }
    }

    fn sorted_overlap_signature(
        overlaps: Option<&OverlappingWindows>,
    ) -> Vec<(usize, u64, u64, f64)> {
        let mut signature: Vec<_> = overlaps
            .map(|overlaps| {
                overlaps
                    .windows
                    .iter()
                    .map(|window| {
                        (
                            window.idx,
                            window.start(),
                            window.end(),
                            window.overlap_fraction,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default();
        signature.sort_by(|left_hit, right_hit| {
            left_hit
                .0
                .cmp(&right_hit.0)
                .then(left_hit.1.cmp(&right_hit.1))
                .then(left_hit.2.cmp(&right_hit.2))
        });
        signature
    }

    fn sorted_context_bed_overlap_signature(
        overlaps: Option<&OverlappingWindows>,
    ) -> Vec<(usize, Option<u64>, u64, u64, f64)> {
        let mut signature: Vec<_> = overlaps
            .map(|overlaps| {
                overlaps
                    .windows
                    .iter()
                    .map(|window| {
                        (
                            window.idx,
                            window.output_idx,
                            window.start(),
                            window.end(),
                            window.overlap_fraction,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default();
        signature.sort_by(|left_hit, right_hit| {
            left_hit
                .0
                .cmp(&right_hit.0)
                .then(left_hit.1.cmp(&right_hit.1))
                .then(left_hit.2.cmp(&right_hit.2))
                .then(left_hit.3.cmp(&right_hit.3))
        });
        signature
    }

    fn sorted_generic_bed_overlap_signature(
        all_windows: &[IndexedInterval<u64>],
        overlaps: Option<&OverlappingWindows>,
    ) -> Vec<(usize, Option<u64>, u64, u64, f64)> {
        let mut signature: Vec<_> = overlaps
            .map(|overlaps| {
                overlaps
                    .windows
                    .iter()
                    .map(|window| {
                        (
                            window.idx,
                            Some(all_windows[window.idx].idx()),
                            window.start(),
                            window.end(),
                            window.overlap_fraction,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default();
        signature.sort_by(|left_hit, right_hit| {
            left_hit
                .0
                .cmp(&right_hit.0)
                .then(left_hit.1.cmp(&right_hit.1))
                .then(left_hit.2.cmp(&right_hit.2))
                .then(left_hit.3.cmp(&right_hit.3))
        });
        signature
    }

    fn assert_overlap_sets_close(
        observed: Option<&OverlappingWindows>,
        expected: Option<&OverlappingWindows>,
    ) {
        let observed = sorted_overlap_signature(observed);
        let expected = sorted_overlap_signature(expected);
        assert_eq!(observed.len(), expected.len());
        for (observed_hit, expected_hit) in observed.iter().zip(expected.iter()) {
            assert_eq!(observed_hit.0, expected_hit.0);
            assert_eq!(observed_hit.1, expected_hit.1);
            assert_eq!(observed_hit.2, expected_hit.2);
            assert!(
                (observed_hit.3 - expected_hit.3).abs() < 1e-12,
                "observed fraction {} != expected fraction {}",
                observed_hit.3,
                expected_hit.3
            );
        }
    }

    fn assert_bed_overlap_sets_close(
        all_windows: &[IndexedInterval<u64>],
        observed: Option<&OverlappingWindows>,
        expected: Option<&OverlappingWindows>,
    ) {
        let observed = sorted_context_bed_overlap_signature(observed);
        let expected = sorted_generic_bed_overlap_signature(all_windows, expected);
        assert_eq!(observed.len(), expected.len());
        for (observed_hit, expected_hit) in observed.iter().zip(expected.iter()) {
            assert_eq!(observed_hit.0, expected_hit.0);
            assert_eq!(observed_hit.1, expected_hit.1);
            assert_eq!(observed_hit.2, expected_hit.2);
            assert_eq!(observed_hit.3, expected_hit.3);
            assert!(
                (observed_hit.4 - expected_hit.4).abs() < 1e-12,
                "observed fraction {} != expected fraction {}",
                observed_hit.4,
                expected_hit.4
            );
        }
    }

    fn assert_manual_signature_close(
        observed: Option<&OverlappingWindows>,
        expected: &[(usize, u64, u64, f64)],
    ) {
        let observed = sorted_overlap_signature(observed);
        assert_eq!(observed.len(), expected.len());
        for (observed_hit, expected_hit) in observed.iter().zip(expected.iter()) {
            assert_eq!(observed_hit.0, expected_hit.0);
            assert_eq!(observed_hit.1, expected_hit.1);
            assert_eq!(observed_hit.2, expected_hit.2);
            assert!(
                (observed_hit.3 - expected_hit.3).abs() < 1e-12,
                "observed fraction {} != expected fraction {}",
                observed_hit.3,
                expected_hit.3
            );
        }
    }

    #[test]
    fn tiered_context_matches_generic_finder_for_mixed_covering_broad_and_narrow_windows(
    ) -> Result<()> {
        let chrom_len = 4_000_000;
        let broad_window_min_bp = 100_000;
        let windows = mixed_covering_broad_narrow_windows()?;
        let chromosome_windows = ChromosomeBedWindows::from_indexed_windows(
            windows.as_slice(),
            broad_window_min_bp,
        );
        let spans = full_bed_spans(&chromosome_windows);
        let min_overlap_fraction = 0.0;
        let look_back = 1_000;
        // The query groups emulate two 1 Mb tiles with 1 kb right reach. Each tile has a
        // tile+reach covering window that should move into the context's always-hit set.
        let tile_query_groups: [(Interval<u64>, &[u64]); 2] = [
            (
                Interval::new(999_000, 2_001_000)?,
                &MIXED_COVERING_BROAD_NARROW_QUERY_STARTS[..5],
            ),
            (
                Interval::new(1_999_000, 3_001_000)?,
                &MIXED_COVERING_BROAD_NARROW_QUERY_STARTS[5..],
            ),
        ];

        for (tile_assignment_envelope, query_starts) in tile_query_groups {
            let mut context = TileBedOverlapContext::new(
                chrom_len,
                &chromosome_windows,
                &spans,
                tile_assignment_envelope,
            )?;
            let mut baseline_wd_ptr = 0;

            for query_start in query_starts {
                let query_interval = Interval::new(*query_start, *query_start + 1_000)?;
                let baseline = find_overlapping_windows(
                    chrom_len,
                    &mut baseline_wd_ptr,
                    Some(windows.as_slice()),
                    None,
                    query_interval,
                    min_overlap_fraction,
                    look_back,
                )?;
                let observed = context.find_overlapping_windows(
                    query_interval,
                    min_overlap_fraction,
                    look_back,
                )?;

                assert_bed_overlap_sets_close(
                    windows.as_slice(),
                    observed.as_ref(),
                    baseline.as_ref(),
                );
            }
        }

        Ok(())
    }

    #[test]
    fn tiered_context_matches_generic_finder_for_nested_mixed_windows() -> Result<()> {
        let chrom_len = 10_000_000;
        let broad_window_min_bp = 100_000;
        let windows = IndexedInterval::from_tuples(&[
            // These stored idx values mimic grouped-BED group ids. The tiered context must still
            // return all_windows positions 0, 1, ... like the generic finder.
            (0, 10_000_000, 900_u64),
            (1_000_000, 1_100_005, 901_u64),
            (1_150_000, 1_150_020, 902_u64),
            (1_250_000, 2_000_000, 903_u64),
            (1_250_010, 1_250_030, 904_u64),
        ])?;
        let chromosome_windows = ChromosomeBedWindows::from_indexed_windows(
            windows.as_slice(),
            broad_window_min_bp,
        );
        let spans = full_bed_spans(&chromosome_windows);
        let tile_assignment_envelope = Interval::new(1_100_000, 1_300_000)?;
        let mut context = TileBedOverlapContext::new(
            chrom_len,
            &chromosome_windows,
            &spans,
            tile_assignment_envelope,
        )?;
        let mut baseline_wd_ptr = 0;
        let min_overlap_fraction = 1.0 / 201.0;
        let look_back = 200_000;
        let queries = [
            Interval::new(1_100_000, 1_100_010)?,
            Interval::new(1_150_000, 1_150_010)?,
            Interval::new(1_250_000, 1_250_020)?,
        ];

        for query_interval in queries {
            let baseline = find_overlapping_windows(
                chrom_len,
                &mut baseline_wd_ptr,
                Some(windows.as_slice()),
                None,
                query_interval,
                min_overlap_fraction,
                look_back,
            )?;
            let observed =
                context.find_overlapping_windows(query_interval, min_overlap_fraction, look_back)?;

            assert_overlap_sets_close(observed.as_ref(), baseline.as_ref());
        }

        let last_observed = context.find_overlapping_windows(
            Interval::new(1_250_010, 1_250_020)?,
            min_overlap_fraction,
            look_back,
        )?;
        assert_manual_signature_close(
            last_observed.as_ref(),
            &[
                (0, 0, 10_000_000, 1.0),
                (3, 1_250_000, 2_000_000, 1.0),
                (4, 1_250_010, 1_250_030, 1.0),
            ],
        );
        let mut output_signature: Vec<_> = last_observed
            .as_ref()
            .expect("last query should overlap windows")
            .windows
            .iter()
            .map(|window| (window.idx, window.output_idx))
            .collect();
        output_signature.sort_by_key(|(idx, _)| *idx);
        assert_eq!(
            output_signature,
            vec![(0, Some(900)), (3, Some(903)), (4, Some(904))]
        );

        Ok(())
    }

    #[test]
    fn tiered_context_uses_look_back_for_broad_pointer_retirement() -> Result<()> {
        let chrom_len = 300;
        let broad_window_min_bp = 20;
        let windows =
            IndexedInterval::from_tuples(&[(100, 135, 40_u64), (160, 180, 41_u64)])?;
        let chromosome_windows = ChromosomeBedWindows::from_indexed_windows(
            windows.as_slice(),
            broad_window_min_bp,
        );
        let spans = full_bed_spans(&chromosome_windows);
        let mut context = TileBedOverlapContext::new(
            chrom_len,
            &chromosome_windows,
            &spans,
            Interval::new(0, 300)?,
        )?;
        let mut baseline_wd_ptr = 0;
        let min_overlap_fraction = 1.0 / 11.0;
        let look_back = 100;

        let first_query = Interval::new(200, 210)?;
        let baseline_first = find_overlapping_windows(
            chrom_len,
            &mut baseline_wd_ptr,
            Some(windows.as_slice()),
            None,
            first_query,
            min_overlap_fraction,
            look_back,
        )?;
        let observed_first =
            context.find_overlapping_windows(first_query, min_overlap_fraction, look_back)?;
        assert_overlap_sets_close(observed_first.as_ref(), baseline_first.as_ref());

        let second_query = Interval::new(130, 140)?;
        let baseline_second = find_overlapping_windows(
            chrom_len,
            &mut baseline_wd_ptr,
            Some(windows.as_slice()),
            None,
            second_query,
            min_overlap_fraction,
            look_back,
        )?;
        let observed_second =
            context.find_overlapping_windows(second_query, min_overlap_fraction, look_back)?;
        assert_overlap_sets_close(observed_second.as_ref(), baseline_second.as_ref());
        assert_manual_signature_close(observed_second.as_ref(), &[(0, 100, 135, 0.5)]);

        Ok(())
    }

    #[test]
    fn narrow_window_containing_tile_envelope_still_matches_generic_finder() -> Result<()> {
        let chrom_len = 1_000;
        let broad_window_min_bp = 10_000;
        let windows = IndexedInterval::from_tuples(&[(0, 1_000, 70_u64)])?;
        let chromosome_windows = ChromosomeBedWindows::from_indexed_windows(
            windows.as_slice(),
            broad_window_min_bp,
        );
        assert!(chromosome_windows.tiers[0].windows.is_empty());
        assert_eq!(chromosome_windows.tiers[1].windows.len(), 1);
        let spans = full_bed_spans(&chromosome_windows);
        let mut context = TileBedOverlapContext::new(
            chrom_len,
            &chromosome_windows,
            &spans,
            Interval::new(100, 200)?,
        )?;
        let query_interval = Interval::new(120, 130)?;
        let mut baseline_wd_ptr = 0;

        let baseline = find_overlapping_windows(
            chrom_len,
            &mut baseline_wd_ptr,
            Some(windows.as_slice()),
            None,
            query_interval,
            1.0,
            0,
        )?;
        let observed = context.find_overlapping_windows(query_interval, 1.0, 0)?;

        assert_overlap_sets_close(observed.as_ref(), baseline.as_ref());
        assert_manual_signature_close(observed.as_ref(), &[(0, 0, 1_000, 1.0)]);

        Ok(())
    }
}
