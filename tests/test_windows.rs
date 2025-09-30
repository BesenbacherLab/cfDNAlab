#[cfg(test)]
mod tests_flattening {
    use cfdnalab::shared::bed::*;

    // Helper: build a start-sorted Windows from (s,e) pairs (original_idx is dummy)
    fn mk_sorted(pairs: &[(u64, u64)]) -> Windows {
        let v: Vec<(u64, u64, u64)> = pairs
            .iter()
            .enumerate()
            .map(|(i, &(s, e))| (s, e, i as u64))
            .collect();
        Windows::from_sorted(v)
    }

    // Helper: assert strictly sorted and non-overlapping (touching should have been merged away)
    fn assert_sorted_non_overlapping(ws: &[(u64, u64, u64)]) {
        for i in 1..ws.len() {
            assert!(
                ws[i - 1].0 <= ws[i].0,
                "not sorted: prev.start={} > cur.start={}",
                ws[i - 1].0,
                ws[i].0
            );
            assert!(
                ws[i - 1].1 < ws[i].0,
                "intervals overlap or still touch: prev={:?}, cur={:?}",
                ws[i - 1],
                ws[i]
            );
            assert!(
                ws[i - 1].0 < ws[i - 1].1,
                "invalid interval with zero/negative length"
            );
        }
        if let Some(last) = ws.last() {
            assert!(
                last.0 < last.1,
                "invalid interval with zero/negative length"
            );
        }
    }

    // Helper: assert indices are sequential starting at start_idx
    fn assert_sequential_indices(ws: &[(u64, u64, u64)], start_idx: u64) {
        for (k, &(_, _, idx)) in ws.iter().enumerate() {
            assert_eq!(idx, start_idx + k as u64, "non-sequential index at k={}", k);
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
        assert_eq!(a, &[(10, 20, 7)]);
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
        assert_eq!(a[0].0, 10);
        assert_eq!(a[0].1, 15);
        assert_eq!(a[1].0, 20);
        assert_eq!(a[1].1, 25);
        assert_eq!(a[2].0, 30);
        assert_eq!(a[2].1, 35);
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
        assert_eq!(a, &[(10, 30, 100)]);
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
        assert_eq!(a[0], (5, 7, 50));
        assert_eq!(a[1], (10, 19, 51));
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
        assert_eq!(a[0].2, 1_000_000);
        assert_eq!(a[1].2, 1_000_001);
        assert_eq!(a[2].2, 1_000_002);
        assert_eq!(next, 1_000_003);
        assert_sorted_non_overlapping(a);
        assert_sequential_indices(a, 1_000_000);
    }
}
