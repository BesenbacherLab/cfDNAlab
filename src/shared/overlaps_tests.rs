mod fixed_width_overlap_cursor_tests {
    use super::super::{FixedWidthOverlapCursor, OverlappingWindows, find_overlapping_windows};
    use crate::{
        Result,
        shared::interval::{IndexedInterval, Interval},
    };
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
