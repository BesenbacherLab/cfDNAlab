mod tests_mergers {
    use crate::commands::prepare_windows::{
        config::{CoordinateSet, MergeLabel, MergeScope},
        labels::LabelTuple,
        mergers::{merge_across_groups, merge_windows, merge_within_groups},
        prepare_windows::Window,
    };
    use std::sync::Arc;

    fn win(chrom: &str, start: u32, end: u32, group: &str) -> Window {
        Window::from_bounds(
            Arc::<str>::from(chrom.to_string()),
            start,
            end,
            start,
            end,
            vec![LabelTuple::new(group.to_string())],
            group.to_string(),
            None,
        )
        .expect("test window should be valid")
    }

    fn snapshot(windows: &[Window]) -> Vec<(String, u32, u32, Vec<String>)> {
        windows
            .iter()
            .map(|w| {
                (
                    w.chrom.as_ref().to_string(),
                    w.resized_start(),
                    w.resized_end(),
                    w.label_tuples.iter().map(|t| t.input.clone()).collect(),
                )
            })
            .collect()
    }

    #[test]
    fn merge_within_groups_merges_overlaps() {
        let windows = vec![
            win("chr1", 0, 5, "A"),
            win("chr1", 4, 8, "A"),
            win("chr1", 20, 25, "B"),
        ];
        let merged =
            merge_within_groups(windows, 2, MergeLabel::Join, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![
                ("chr1".into(), 0, 8, vec!["A".into()]),
                ("chr1".into(), 20, 25, vec!["B".into()]),
            ]
        );
    }

    #[test]
    fn merge_within_groups_respects_gap_threshold() {
        let windows = vec![win("chr1", 0, 4, "A"), win("chr1", 7, 10, "A")];
        let merged = merge_within_groups(
            windows.clone(),
            2,
            MergeLabel::Join,
            CoordinateSet::Resized,
            false,
        );
        assert_eq!(snapshot(&merged), snapshot(&windows));
    }

    #[test]
    fn merge_within_groups_bridges_gap_within_threshold() {
        let windows = vec![win("chr1", 0, 4, "A"), win("chr1", 6, 9, "A")];
        let merged =
            merge_within_groups(windows, 2, MergeLabel::Join, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![("chr1".into(), 0, 9, vec!["A".into()])]
        );
    }

    #[test]
    fn merge_across_groups_joins_labels() {
        let windows = vec![win("chr1", 0, 4, "G1"), win("chr1", 3, 6, "G2")];
        let merged =
            merge_across_groups(windows, 1, MergeLabel::Join, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![("chr1".into(), 0, 6, vec!["G1".into(), "G2".into()])]
        );
    }

    #[test]
    fn merge_across_groups_sorts_unsorted_input() {
        let windows = vec![win("chr1", 5, 7, "B"), win("chr1", 2, 6, "A")];
        let merged =
            merge_across_groups(windows, 1, MergeLabel::First, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![("chr1".into(), 2, 7, vec!["A".into()])]
        );
    }

    #[test]
    fn merge_across_groups_honors_first_label_policy() {
        let windows = vec![win("chr1", 0, 4, "G1"), win("chr1", 3, 6, "G2")];
        let merged =
            merge_across_groups(windows, 1, MergeLabel::First, CoordinateSet::Resized, false);
        assert_eq!(
            snapshot(&merged),
            vec![("chr1".into(), 0, 6, vec!["G1".into()])]
        );
    }

    #[test]
    fn merge_windows_respects_scope_none() {
        let windows = vec![win("chr1", 0, 5, "A")];
        let merged = merge_windows(
            windows.clone(),
            MergeScope::None,
            Some(3),
            MergeLabel::Join,
            CoordinateSet::Resized,
            false,
        );
        assert_eq!(snapshot(&merged), snapshot(&windows));
    }

    #[test]
    fn merge_windows_returns_original_when_gap_missing() {
        let windows = vec![win("chr1", 0, 5, "A"), win("chr1", 10, 12, "A")];
        let merged = merge_windows(
            windows.clone(),
            MergeScope::Within,
            None,
            MergeLabel::Join,
            CoordinateSet::Resized,
            false,
        );
        assert_eq!(snapshot(&merged), snapshot(&windows));
    }
}
