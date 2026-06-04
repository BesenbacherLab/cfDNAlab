mod tests_postprocess {

    use super::super::build_coverage_index;
    use crate::commands::prepare_windows::{
        config::{CoordinateSet, DedupKeep, DistancePolicy, MergeScope},
        labels::LabelTuple,
        postprocess::{
            deduplicate_identical, enforce_min_distance_within_group, partition_safe_and_tail,
        },
        prepare_windows::Window,
    };
    use std::sync::Arc;

    fn win(chrom: &str, start: u32, end: u32, group: &str, score: Option<f32>) -> Window {
        Window::from_bounds(
            Arc::<str>::from(chrom.to_string()),
            start,
            end,
            start,
            end,
            vec![LabelTuple::new(group.to_string())],
            group.to_string(),
            score,
        )
        .expect("test window should be valid")
    }

    fn snapshot(windows: &[Window]) -> Vec<(String, u32, u32, String, Option<f32>)> {
        windows
            .iter()
            .map(|w| {
                (
                    w.chrom.as_ref().to_string(),
                    w.resized_start(),
                    w.resized_end(),
                    w.group_key.clone(),
                    w.score,
                )
            })
            .collect()
    }

    #[test]
    fn build_coverage_index_combines_same_position_deltas_without_losing_coverage() {
        // Arrange
        // The boundaries at 15 and 20 contain both starts and ends. In particular, the net delta at
        // 15 is zero, so the segment split must be preserved without changing the running depth.
        let windows = vec![
            win("chr1", 10, 20, "group", None),
            win("chr1", 10, 15, "group", None),
            win("chr1", 15, 20, "group", None),
            win("chr1", 20, 30, "group", None),
        ];

        // Act
        let (boundaries, coverage_by_segment, coverage_prefix) =
            build_coverage_index(&windows, CoordinateSet::Resized);

        // Assert
        assert_eq!(boundaries, vec![10, 15, 20, 30]);
        assert_eq!(coverage_by_segment, vec![2, 2, 1]);
        assert_eq!(coverage_prefix, vec![0, 10, 20, 30]);
    }

    #[test]
    fn dedup_none_keeps_all_windows() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 10, 20, "g1", Some(2.0)),
        ];
        let result = deduplicate_identical(
            windows.clone(),
            DedupKeep::None,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(snapshot(&result), snapshot(&windows));
    }

    #[test]
    fn dedup_keep_first_prefers_first_duplicate() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 10, 20, "g1", Some(5.0)),
        ];
        let result =
            deduplicate_identical(windows, DedupKeep::KeepFirst, true, CoordinateSet::Resized);
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", Some(1.0))])
        );
    }

    #[test]
    fn dedup_keep_highest_score_uses_scores_when_available() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 10, 20, "g1", Some(5.0)),
            win("chr1", 10, 20, "g1", Some(2.5)),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", Some(5.0))])
        );
    }

    #[test]
    fn dedup_keep_highest_score_falls_back_without_scores() {
        let windows = vec![
            win("chr1", 10, 20, "g1", None),
            win("chr1", 10, 20, "g1", None),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepHighestScore,
            false,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", None)])
        );
    }

    #[test]
    fn dedup_keep_lowest_score_picks_smallest_score() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(3.0)),
            win("chr1", 10, 20, "g1", Some(1.5)),
            win("chr1", 10, 20, "g1", Some(4.0)),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepLowestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", Some(1.5))])
        );
    }

    #[test]
    fn dedup_keep_lowest_score_falls_back_without_scores() {
        let windows = vec![
            win("chr1", 10, 20, "g1", None),
            win("chr1", 10, 20, "g1", None),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepLowestScore,
            false,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 10, 20, "g1", None)])
        );
    }

    #[test]
    fn dedup_does_not_touch_unique_windows() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 30, 40, "g1", Some(2.0)),
            win("chr2", 5, 15, "", None),
        ];
        let result = deduplicate_identical(
            windows.clone(),
            DedupKeep::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(snapshot(&result), snapshot(&windows));
    }

    #[test]
    fn dedup_handles_multiple_duplicate_groups() {
        let windows = vec![
            win("chr1", 10, 20, "g1", Some(1.0)),
            win("chr1", 10, 20, "g1", Some(3.0)),
            win("chr1", 30, 40, "g2", Some(5.0)),
            win("chr1", 30, 40, "g2", Some(2.0)),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[
                win("chr1", 10, 20, "g1", Some(3.0)),
                win("chr1", 30, 40, "g2", Some(5.0)),
            ])
        );
    }

    #[test]
    fn dedup_keep_highest_score_prefers_non_none_scores() {
        let windows = vec![
            win("chr1", 0, 5, "g", None),
            win("chr1", 0, 5, "g", Some(1.0)),
            win("chr1", 0, 5, "g", Some(2.0)),
        ];
        let result = deduplicate_identical(
            windows,
            DedupKeep::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 0, 5, "g", Some(2.0))])
        );
    }

    #[test]
    fn enforce_min_distance_within_group_keep_first() {
        let windows = vec![
            win("chr1", 0, 10, "g", Some(1.0)),
            win("chr1", 4, 12, "g", Some(2.0)),
            win("chr1", 20, 30, "g", Some(3.0)),
        ];
        let result = enforce_min_distance_within_group(
            windows,
            Some(5),
            DistancePolicy::KeepFirst,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[
                win("chr1", 0, 10, "g", Some(1.0)),
                win("chr1", 20, 30, "g", Some(3.0))
            ])
        );
    }

    #[test]
    fn enforce_min_distance_within_group_keep_highest_score() {
        let windows = vec![
            win("chr1", 0, 10, "g", Some(1.0)),
            win("chr1", 4, 12, "g", Some(5.0)),
            win("chr1", 40, 50, "g", Some(2.0)),
        ];
        let result = enforce_min_distance_within_group(
            windows,
            Some(8),
            DistancePolicy::KeepHighestScore,
            true,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[
                win("chr1", 4, 12, "g", Some(5.0)),
                win("chr1", 40, 50, "g", Some(2.0))
            ])
        );
    }

    #[test]
    fn enforce_min_distance_within_group_keep_lowest_score_without_scores() {
        let windows = vec![
            win("chr1", 0, 5, "g", None),
            win("chr1", 3, 9, "g", None),
            win("chr1", 20, 25, "g", None),
        ];
        let result = enforce_min_distance_within_group(
            windows,
            Some(4),
            DistancePolicy::KeepLowestScore,
            false,
            CoordinateSet::Resized,
        );
        assert_eq!(
            snapshot(&result),
            snapshot(&[win("chr1", 0, 5, "g", None), win("chr1", 20, 25, "g", None)])
        );
    }

    #[test]
    fn partition_safe_and_tail_without_margin_writes_all() {
        let windows = vec![
            win("chr1", 0, 10, "g", None),
            win("chr1", 20, 30, "g", None),
        ];
        let (safe, tail) = partition_safe_and_tail(
            windows,
            None,
            MergeScope::Within,
            None,
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );
        assert_eq!(
            snapshot(&safe),
            snapshot(&[
                win("chr1", 0, 10, "g", None),
                win("chr1", 20, 30, "g", None)
            ])
        );
        assert!(tail.is_empty());
    }

    #[test]
    fn partition_safe_and_tail_retains_last_window_when_min_distance_crosses_chunk() {
        let windows = vec![
            win("chr1", 0, 5, "g1", None),
            win("chr1", 10, 15, "g1", None),
        ];
        let (safe, tail) = partition_safe_and_tail(
            windows,
            Some(4),
            MergeScope::Within,
            Some(0),
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );
        // The first window cannot reach the boundary. Only the suffix remains in tail
        assert_eq!(snapshot(&safe), snapshot(&[win("chr1", 0, 5, "g1", None)]));
        assert_eq!(
            snapshot(&tail),
            snapshot(&[win("chr1", 10, 15, "g1", None)])
        );
    }

    #[test]
    fn partition_safe_and_tail_across_scope_keeps_boundary_suffix_only() {
        let windows = vec![
            win("chr1", 0, 4, "g1", None),
            win("chr1", 5, 7, "g2", None),
            win("chr1", 20, 25, "g1", None),
        ];
        let (safe, tail) = partition_safe_and_tail(
            windows,
            Some(3),
            MergeScope::Across,
            Some(2),
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );
        assert_eq!(safe.len(), 2);
        assert_eq!(snapshot(&tail).len(), 1);
    }

    #[test]
    fn partition_safe_and_tail_across_scope_keeps_overlap_chain_in_tail() {
        let windows = vec![
            // Early windows that should finalize
            win("chr1", 0, 4, "g1", None),
            win("chr1", 15, 18, "g2", None),
            // Overlap/merge chain near the boundary - each link extends reach
            win("chr1", 40, 45, "g3", None),
            win("chr1", 47, 49, "g4", None), // within margin of previous end
            win("chr1", 52, 54, "g5", None), // extends chain to boundary so all three stay in tail
        ];
        let (safe, tail) = partition_safe_and_tail(
            windows,
            None,
            MergeScope::Across,
            Some(3),
            CoordinateSet::Resized,
            CoordinateSet::Resized,
            None,
            CoordinateSet::Resized,
            None,
        );
        assert_eq!(snapshot(&safe).len(), 2);
        assert_eq!(snapshot(&tail).len(), 3);
    }
}
