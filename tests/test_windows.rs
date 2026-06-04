#[cfg(test)]
mod tests_flattening {
    use cfdnalab::shared::bed::*;
    use cfdnalab::shared::interval::{IndexedInterval, ScoredInterval, Span};

    // Helper: build a start-sorted Windows from (s,e) pairs (original_idx is dummy)
    fn mk_sorted(pairs: &[(u64, u64)]) -> Windows {
        let windows = pairs
            .iter()
            .enumerate()
            .map(|(i, &(start, end))| {
                IndexedInterval::new(start, end, i as u64).expect("test windows should be valid")
            })
            .collect();
        Windows::from_sorted(windows)
    }

    // Helper: assert strictly sorted and non-overlapping (touching should have been merged away)
    fn assert_sorted_non_overlapping(ws: &[IndexedInterval<u64>]) {
        for index in 1..ws.len() {
            assert!(
                ws[index - 1].start() <= ws[index].start(),
                "not sorted: prev.start={} > cur.start={}",
                ws[index - 1].start(),
                ws[index].start()
            );
            assert!(
                ws[index - 1].end() < ws[index].start(),
                "intervals overlap or still touch: prev={:?}, cur={:?}",
                ws[index - 1],
                ws[index]
            );
            assert!(
                ws[index - 1].start() < ws[index - 1].end(),
                "invalid interval with zero/negative length"
            );
        }
        if let Some(last) = ws.last() {
            assert!(
                last.start() < last.end(),
                "invalid interval with zero/negative length"
            );
        }
    }

    // Helper: assert indices are sequential starting at start_idx
    fn assert_sequential_indices(ws: &[IndexedInterval<u64>], start_idx: u64) {
        for (index, window) in ws.iter().enumerate() {
            assert_eq!(
                window.idx(),
                start_idx + index as u64,
                "non-sequential index at k={}",
                index
            );
        }
    }

    #[test]
    fn flatten_empty() {
        let w = Windows::from_sorted(Vec::new());
        let (flat, next) = w.into_flattened_reindexed(0);
        assert!(flat.as_slice().is_empty());
        assert_eq!(next, 0);
        assert_eq!(flat.span_start(), 0);
        assert_eq!(flat.span_end(), 0);
    }

    #[test]
    fn flatten_singleton() {
        let w = mk_sorted(&[(10, 20)]);
        let (flat, next) = w.into_flattened_reindexed(7);
        let a = flat.as_slice();
        assert_eq!(a[0].into_tuple(), (10, 20, 7));
        assert_eq!(next, 8);
        assert_eq!(flat.span_start(), 10);
        assert_eq!(flat.span_end(), 20);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 7);
    }

    #[test]
    fn flatten_no_merges() {
        // Non-overlapping and non-touching
        let w = mk_sorted(&[(10, 15), (20, 25), (30, 35)]);
        let (flat, next) = w.into_flattened_reindexed(0);
        let a = flat.as_slice();
        assert_eq!(a.len(), 3);
        assert_eq!(next, 3);
        // Starts/ends preserved
        assert_eq!(a[0].start(), 10);
        assert_eq!(a[0].end(), 15);
        assert_eq!(a[1].start(), 20);
        assert_eq!(a[1].end(), 25);
        assert_eq!(a[2].start(), 30);
        assert_eq!(a[2].end(), 35);
        assert_eq!(flat.span_start(), 10);
        assert_eq!(flat.span_end(), 35);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 0);
    }

    #[test]
    fn flatten_touching_merges() {
        // Touching intervals must merge (half-open semantics)
        let w = mk_sorted(&[(10, 15), (15, 20), (20, 30)]);
        let (flat, next) = w.into_flattened_reindexed(100);
        let a = flat.as_slice();
        assert_eq!(a[0].into_tuple(), (10, 30, 100));
        assert_eq!(next, 101);
        assert_eq!(flat.span_start(), 10);
        assert_eq!(flat.span_end(), 30);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 100);
    }

    #[test]
    fn flatten_overlapping_chain() {
        // Mixed: one disjoint small block and a chain that overlaps/touches
        let w = mk_sorted(&[(5, 7), (10, 14), (12, 16), (16, 19)]);
        let (flat, next) = w.into_flattened_reindexed(50);
        let a = flat.as_slice();
        // Expect (5,7) and (10,19)
        assert_eq!(a.len(), 2);
        assert_eq!(a[0].into_tuple(), (5, 7, 50));
        assert_eq!(a[1].into_tuple(), (10, 19, 51));
        assert_eq!(next, 52);
        assert_eq!(flat.span_start(), 5);
        assert_eq!(flat.span_end(), 19);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 50);
    }

    #[test]
    fn flatten_large_start_idx() {
        // Sanity: indices carry forward correctly from large start
        let w = mk_sorted(&[(0, 1), (2, 3), (4, 5)]);
        let (flat, next) = w.into_flattened_reindexed(1_000_000);
        let a = flat.as_slice();
        assert_eq!(a.len(), 3);
        assert_eq!(a[0].idx(), 1_000_000);
        assert_eq!(a[1].idx(), 1_000_001);
        assert_eq!(a[2].idx(), 1_000_002);
        assert_eq!(next, 1_000_003);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 1_000_000);
    }

    #[test]
    fn grouped_windows_sort_and_preserve_group_indices() {
        // Arrange:
        // - Inputs are unsorted by start, but group indices are row data and must survive sorting.
        // - After sorting by start we expect [10,15) idx 3, [15,18) idx 5, [20,30) idx 7.
        // - Span is therefore min start 10 and max end 30.
        let grouped = GroupedWindows::from_tuples(&[(20, 30, 7), (10, 15, 3), (15, 18, 5)], None)
            .expect("grouped test windows should be valid");

        let windows = grouped.windows_as_slice();

        assert_eq!(windows.len(), 3);
        assert_eq!(windows[0].into_tuple(), (10, 15, 3));
        assert_eq!(windows[1].into_tuple(), (15, 18, 5));
        assert_eq!(windows[2].into_tuple(), (20, 30, 7));
        assert_eq!(grouped.span_start(), 10);
        assert_eq!(grouped.span_end(), 30);
    }

    #[test]
    fn grouped_windows_sort_and_preserve_optional_strands() {
        // Arrange:
        // - Inputs are unsorted by start.
        // - Strand values are row data, so they must move with their windows.
        // - After sorting by start, [10,15) keeps Forward, [15,18) keeps Unstranded,
        //   and [20,30) keeps Reverse.
        let grouped = GroupedWindows::new(
            vec![
                IndexedInterval::new(20, 30, 7).expect("grouped interval should be valid"),
                IndexedInterval::new(10, 15, 3).expect("grouped interval should be valid"),
                IndexedInterval::new(15, 18, 5).expect("grouped interval should be valid"),
            ],
            Some(vec![Strand::Reverse, Strand::Forward, Strand::Unstranded]),
        );

        let windows = grouped.windows_as_slice();
        let strands = grouped
            .strands
            .as_ref()
            .expect("strand metadata should be retained");

        assert_eq!(windows[0].into_tuple(), (10, 15, 3));
        assert_eq!(strands[0], Strand::Forward);
        assert_eq!(windows[1].into_tuple(), (15, 18, 5));
        assert_eq!(strands[1], Strand::Unstranded);
        assert_eq!(windows[2].into_tuple(), (20, 30, 7));
        assert_eq!(strands[2], Strand::Reverse);
    }

    #[test]
    fn grouped_windows_span_uses_max_end_not_last_sorted_end() {
        // Sorting by start yields [10,40), [20,25), [30,32). The last sorted window ends at 32,
        // but the collection span must use the true maximum end 40.
        let grouped = GroupedWindows::new(
            vec![
                IndexedInterval::new(20, 25, 0).expect("grouped interval should be valid"),
                IndexedInterval::new(10, 40, 1).expect("grouped interval should be valid"),
                IndexedInterval::new(30, 32, 2).expect("grouped interval should be valid"),
            ],
            None,
        );

        assert_eq!(grouped.span(), Span::new(10, 40).unwrap());
    }

    #[test]
    fn grouped_windows_empty_has_zero_span() {
        let grouped = GroupedWindows::from_sorted(Vec::new(), None);

        assert!(grouped.windows_as_slice().is_empty());
        assert_eq!(grouped.span_start(), 0);
        assert_eq!(grouped.span_end(), 0);
    }

    #[test]
    fn grouped_windows_from_tuples_rejects_invalid_interval() {
        let error = GroupedWindows::from_tuples(&[(10, 10, 3)], None)
            .expect_err("invalid grouped interval should fail");

        assert_eq!(
            error.to_string(),
            "interval end (10) must be greater than start (10)"
        );
    }

    #[test]
    fn scored_windows_sort_and_preserve_scores() {
        // Arrange:
        // - Sorting by start should reorder [20,30) and [10,15) into [10,15), [20,30).
        // - Score and original index are interval data, so they must stay attached to their intervals.
        // - Span is the overall covered range [10,30).
        let scored = ScoredWindows::from_tuples(&[(20, 30, 7, 1.5), (10, 15, 3, 2.5)])
            .expect("scored test windows should be valid");

        let windows = scored.as_slice();

        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].into_tuple(), (10, 15, 3, 2.5));
        assert_eq!(windows[1].into_tuple(), (20, 30, 7, 1.5));
        assert_eq!(scored.span_start(), 10);
        assert_eq!(scored.span_end(), 30);
    }

    #[test]
    fn scored_windows_span_uses_max_end_not_last_sorted_end() {
        // Sorting by start yields [10,45), [20,25), [30,33). As with grouped windows, span_end
        // must be the global maximum end 45 rather than the last sorted end 33.
        let scored = ScoredWindows::new(vec![
            ScoredInterval::new(20, 25, 0, 0.5).expect("scored interval should be valid"),
            ScoredInterval::new(10, 45, 1, 1.5).expect("scored interval should be valid"),
            ScoredInterval::new(30, 33, 2, 2.5).expect("scored interval should be valid"),
        ]);

        assert_eq!(scored.span(), Span::new(10, 45).unwrap());
    }

    #[test]
    fn scored_windows_to_windows_drops_score_but_keeps_interval_and_index() {
        // Converting scored windows to plain windows should discard only the score field.
        // Interval bounds, original indices, and the collection span must remain unchanged.
        let scored = ScoredWindows::new(vec![
            ScoredInterval::new(5, 9, 11, 0.5).expect("scored interval should be valid"),
            ScoredInterval::new(10, 15, 12, 1.0).expect("scored interval should be valid"),
        ]);

        let plain = scored.to_windows();
        let windows = plain.as_slice();

        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].into_tuple(), (5, 9, 11));
        assert_eq!(windows[1].into_tuple(), (10, 15, 12));
        assert_eq!(plain.span(), Span::new(5, 15).unwrap());
    }

    #[test]
    fn windows_from_tuples_rejects_invalid_interval() {
        let error =
            Windows::from_tuples(&[(12, 12, 0)]).expect_err("invalid plain interval should fail");

        assert_eq!(
            error.to_string(),
            "interval end (12) must be greater than start (12)"
        );
    }

    #[test]
    fn scored_windows_from_tuples_rejects_invalid_interval() {
        let error = ScoredWindows::from_tuples(&[(20, 19, 7, 1.5)])
            .expect_err("invalid scored interval should fail");

        assert_eq!(
            error.to_string(),
            "interval end (19) must be greater than start (20)"
        );
    }
}
