mod tests_visualize_positions {
    use std::num::NonZeroUsize;

    use cfdnalab::commands::visualize_positions::{
        PositionsSpec, RangeParseError, ReadClamp, ReferenceFrame, build_kmer_start_overlays,
        build_tracks_for_length, parse_positions,
    };

    fn default_step() -> NonZeroUsize {
        NonZeroUsize::new(1).unwrap()
    }

    fn take_linear_indices(
        length: u32,
        frame: ReferenceFrame,
        positions: &PositionsSpec,
        step: NonZeroUsize,
    ) -> Vec<Vec<i32>> {
        take_linear_indices_with_clamp(length, frame, positions, step, ReadClamp::None)
    }

    fn take_linear_indices_with_clamp(
        length: u32,
        frame: ReferenceFrame,
        positions: &PositionsSpec,
        step: NonZeroUsize,
        clamp: ReadClamp,
    ) -> Vec<Vec<i32>> {
        let viz = build_tracks_for_length(length, frame, positions, step, clamp);
        viz.tracks
            .iter()
            .map(|track| track.selected_indices.clone())
            .collect()
    }

    #[test]
    fn nearest_open_to_half_small_l() {
        let spec = parse_positions(ReferenceFrame::Nearest, "10..").unwrap();
        let tracks = take_linear_indices(18, ReferenceFrame::Nearest, &spec, default_step());
        assert_eq!(tracks.len(), 2);
        assert!(tracks[1].is_empty());
    }

    #[test]
    fn nearest_half_minus_k() {
        let spec = parse_positions(ReferenceFrame::Nearest, "5..half-3").unwrap();
        let tracks = take_linear_indices(151, ReferenceFrame::Nearest, &spec, default_step());
        let expected: Vec<i32> = (5..=72).collect();
        assert_eq!(tracks[1], expected);
        assert!(tracks[0].contains(&5));
        assert!(tracks[0].contains(&(151 - 5 + 1)));
    }

    #[test]
    fn left_opposite_end_bound() {
        let spec = parse_positions(ReferenceFrame::Left, "10..-10").unwrap();
        let tracks = take_linear_indices(100, ReferenceFrame::Left, &spec, default_step());
        let expected: Vec<i32> = (10..=90).collect();
        assert_eq!(tracks[0], expected);
    }

    #[test]
    fn right_opposite_end_bound() {
        let spec = parse_positions(ReferenceFrame::Right, "10..-10").unwrap();
        let tracks = take_linear_indices(101, ReferenceFrame::Right, &spec, default_step());
        let expected: Vec<i32> = (10..=91).collect();
        assert_eq!(tracks[0], expected);
    }

    #[test]
    fn per_end_two_tracks() {
        let spec = parse_positions(ReferenceFrame::PerEnd, "..5").unwrap();
        let tracks = take_linear_indices(120, ReferenceFrame::PerEnd, &spec, default_step());
        let expected: Vec<i32> = (1..=5).collect();
        assert_eq!(tracks.len(), 2);
        assert_eq!(tracks[0], expected);
        assert_eq!(tracks[1], expected);
    }

    #[test]
    fn per_end_stride_applies_independently() {
        let spec = parse_positions(ReferenceFrame::PerEnd, "1..10").unwrap();
        let viz = build_tracks_for_length(
            12,
            ReferenceFrame::PerEnd,
            &spec,
            NonZeroUsize::new(3).unwrap(),
            ReadClamp::None,
        );
        let left_track = viz
            .tracks
            .iter()
            .find(|track| track.name == "left")
            .expect("missing left track");
        let right_track = viz
            .tracks
            .iter()
            .find(|track| track.name == "right")
            .expect("missing right track");
        assert_eq!(left_track.selected_indices, vec![1, 4, 7, 10]);
        assert_eq!(right_track.selected_indices, vec![1, 4, 7, 10]);
    }

    #[test]
    fn left_trim_both_ends_extended() {
        let spec = parse_positions(ReferenceFrame::Left, "15..-15").unwrap();
        let tracks = take_linear_indices(80, ReferenceFrame::Left, &spec, default_step());
        let expected: Vec<i32> = (15..=65).collect();
        assert_eq!(tracks[0], expected);
    }

    #[test]
    fn selects_full_axis_when_positions_is_all_for_left_frame() {
        let spec = parse_positions(ReferenceFrame::Left, "..").unwrap();
        let length = 25;
        let tracks = take_linear_indices(length, ReferenceFrame::Left, &spec, default_step());
        let expected: Vec<i32> = (1..=length as i32).collect();
        assert_eq!(tracks[0], expected);
    }

    #[test]
    fn left_half_range_includes_first_half() {
        let spec = parse_positions(ReferenceFrame::Left, "..half").unwrap();
        let tracks = take_linear_indices(100, ReferenceFrame::Left, &spec, default_step());
        let expected: Vec<i32> = (1..=50).collect();
        assert_eq!(tracks[0], expected);
    }

    #[test]
    fn left_half_minus_offset() {
        let spec = parse_positions(ReferenceFrame::Left, "10..half-5").unwrap();
        let tracks = take_linear_indices(120, ReferenceFrame::Left, &spec, default_step());
        let expected_end = 120 / 2 - 5;
        let expected: Vec<i32> = (10..=expected_end).collect();
        assert_eq!(tracks[0], expected);
    }

    #[test]
    fn mid_symmetric_closed() {
        let spec = parse_positions(ReferenceFrame::Mid, "-10..10").unwrap();
        let tracks = take_linear_indices(99, ReferenceFrame::Mid, &spec, default_step());
        let expected: Vec<i32> = (-10..=10).collect();
        assert_eq!(tracks[0], expected);
    }

    #[test]
    fn mid_open_right() {
        let spec = parse_positions(ReferenceFrame::Mid, "..5").unwrap();
        let tracks = take_linear_indices(150, ReferenceFrame::Mid, &spec, default_step());
        let expected: Vec<i32> = (0..=5).collect();
        assert_eq!(tracks[0], expected);

        let legacy = parse_positions(ReferenceFrame::Mid, "..+5").unwrap();
        let legacy_tracks = take_linear_indices(150, ReferenceFrame::Mid, &legacy, default_step());
        assert_eq!(legacy_tracks[0], expected);
    }

    #[test]
    fn should_keep_origin_when_mid_stride_applied() {
        let spec = parse_positions(ReferenceFrame::Mid, "-6..6").unwrap();
        let step = NonZeroUsize::new(3).unwrap();
        let tracks = take_linear_indices(101, ReferenceFrame::Mid, &spec, step);
        assert_eq!(tracks[0], vec![-6, -3, 0, 3, 6]);
    }

    #[test]
    fn selects_full_axis_when_positions_is_all_for_mid_frame() {
        let length = 18;
        let spec = parse_positions(ReferenceFrame::Mid, "..").unwrap();
        let tracks = take_linear_indices(length, ReferenceFrame::Mid, &spec, default_step());
        let half = (length / 2) as i32;
        let axis_end = if length % 2 == 0 { half - 1 } else { half };
        let expected: Vec<i32> = (-half..=axis_end).collect();
        assert_eq!(tracks[0], expected);
    }

    #[test]
    fn stride_application() {
        let spec = parse_positions(ReferenceFrame::Left, "1..10").unwrap();
        let step = NonZeroUsize::new(3).unwrap();
        let tracks = take_linear_indices(20, ReferenceFrame::Left, &spec, step);
        assert_eq!(tracks[0], vec![1, 4, 7, 10]);
    }

    #[test]
    fn left_overlay_trims_tail_for_kmer_length() {
        let length = 30;
        let spec = parse_positions(ReferenceFrame::Left, "..").unwrap();
        let viz = build_tracks_for_length(
            length,
            ReferenceFrame::Left,
            &spec,
            default_step(),
            ReadClamp::None,
        );
        let k = 3u8;
        let overlays = build_kmer_start_overlays(ReferenceFrame::Left, length, &viz.tracks, &[k]);
        let overlay = overlays
            .iter()
            .find(|track| track.name == "left k-mer starts (k=3)")
            .expect("missing overlay for left frame");
        assert_eq!(overlay.selected_indices.first().copied(), Some(1));
        let expected_last = (length - u32::from(k) + 1) as i32;
        assert_eq!(
            overlay.selected_indices.last().copied(),
            Some(expected_last)
        );
    }

    #[test]
    fn right_overlay_shifts_to_start_coordinates() {
        let length = 30;
        let spec = parse_positions(ReferenceFrame::Right, "..").unwrap();
        let viz = build_tracks_for_length(
            length,
            ReferenceFrame::Right,
            &spec,
            default_step(),
            ReadClamp::None,
        );
        let k = 3u8;
        let overlays = build_kmer_start_overlays(ReferenceFrame::Right, length, &viz.tracks, &[k]);
        let overlay = overlays
            .iter()
            .find(|track| track.name == "right k-mer starts (k=3)")
            .expect("missing overlay for right frame");
        assert_eq!(overlay.selected_indices.first().copied(), Some(1));
        let expected_last = (length - u32::from(k) + 1) as i32;
        assert_eq!(
            overlay.selected_indices.last().copied(),
            Some(expected_last)
        );
    }

    #[test]
    fn nearest_center_double_count_guard() {
        let spec = parse_positions(ReferenceFrame::Nearest, "..half").unwrap();
        let tracks = take_linear_indices(100, ReferenceFrame::Nearest, &spec, default_step());
        assert_eq!(tracks[1].last().copied(), Some(50));
        assert_eq!(tracks[1].iter().filter(|&&v| v == 50).count(), 1);
        assert!(tracks[0].contains(&1));
        assert!(tracks[0].contains(&100));
    }

    #[test]
    fn nearest_guard_overlays_obey_midpoint_limits() {
        let length = 101;
        let spec = parse_positions(ReferenceFrame::Nearest, "..").unwrap();
        let viz = build_tracks_for_length(
            length,
            ReferenceFrame::Nearest,
            &spec,
            default_step(),
            ReadClamp::None,
        );

        let overlays =
            build_kmer_start_overlays(ReferenceFrame::Nearest, length, &viz.tracks, &[3]);
        assert_eq!(overlays.len(), 2);

        let fragment_overlay = overlays
            .iter()
            .find(|track| track.name == "fragment k-mer starts (k=3)")
            .expect("missing fragment overlay");
        let nearest_overlay = overlays
            .iter()
            .find(|track| track.name == "nearest k-mer starts (k=3)")
            .expect("missing nearest overlay");

        let len = length as u64;
        let k_span = 3u64;
        let half = len / 2;
        let left_max_start = half.saturating_sub(k_span);
        let right_min_anchor = half.saturating_add(k_span);

        assert!(
            fragment_overlay.selected_indices.iter().all(|&pos| {
                if pos <= 0 {
                    return false;
                }
                let pos_u64 = pos as u64;
                if pos_u64 <= half {
                    (pos_u64 - 1) <= left_max_start
                } else {
                    let anchor_offset = pos_u64 + k_span - 2;
                    anchor_offset >= right_min_anchor
                }
            }),
            "fragment overlay contains positions that cross the midpoint guard"
        );
        let expected_max_distance = {
            let left_max_distance = (half.saturating_sub(k_span) + 1).min(half);
            let right_min_anchor = half.saturating_add(k_span);
            let right_max_distance = len.saturating_sub(right_min_anchor);
            left_max_distance.max(right_max_distance)
        } as i32;
        assert_eq!(
            nearest_overlay.selected_indices.last().copied(),
            Some(expected_max_distance)
        );
        assert!(
            nearest_overlay
                .selected_indices
                .iter()
                .all(|&distance| distance > 0 && distance <= expected_max_distance),
            "folded overlay should stop at the guarded distance"
        );
    }

    #[test]
    fn selects_full_axis_when_positions_is_all_for_nearest_frame() {
        let length = 21;
        let spec = parse_positions(ReferenceFrame::Nearest, "..").unwrap();
        let viz = build_tracks_for_length(
            length,
            ReferenceFrame::Nearest,
            &spec,
            default_step(),
            ReadClamp::None,
        );
        let fragment = viz
            .tracks
            .iter()
            .find(|track| track.name == "fragment")
            .expect("missing fragment track");
        let nearest = viz
            .tracks
            .iter()
            .find(|track| track.name == "nearest")
            .expect("missing nearest track");

        let half = (length / 2) as i32;
        let expected_nearest: Vec<i32> = (1..=half).collect();
        assert_eq!(nearest.selected_indices, expected_nearest);

        let mut expected_fragment = Vec::new();
        for distance in 1..=half {
            expected_fragment.push(distance);
            expected_fragment.push(length as i32 - distance + 1);
        }
        expected_fragment.sort_unstable();
        expected_fragment.dedup();
        assert_eq!(fragment.selected_indices, expected_fragment);
    }

    #[test]
    fn left_clamp_nearest_read_truncates_second_half() {
        let spec = parse_positions(ReferenceFrame::Left, "1..100").unwrap();
        let tracks = take_linear_indices_with_clamp(
            100,
            ReferenceFrame::Left,
            &spec,
            default_step(),
            ReadClamp::Nearest,
        );
        assert_eq!(tracks[0].last().copied(), Some(50));
        assert!(!tracks[0].contains(&51));
    }

    #[test]
    fn right_clamp_nearest_read_truncates_first_half() {
        let spec = parse_positions(ReferenceFrame::Right, "1..100").unwrap();
        let tracks = take_linear_indices_with_clamp(
            100,
            ReferenceFrame::Right,
            &spec,
            default_step(),
            ReadClamp::Nearest,
        );
        assert_eq!(tracks[0].first().copied(), Some(51));
        assert!(!tracks[0].contains(&50));
    }

    #[test]
    fn per_end_clamp_both_reads_splits_tracks() {
        let spec = parse_positions(ReferenceFrame::PerEnd, "1..100").unwrap();
        let tracks = take_linear_indices_with_clamp(
            100,
            ReferenceFrame::PerEnd,
            &spec,
            default_step(),
            ReadClamp::Both,
        );
        assert_eq!(tracks.len(), 2);
        assert_eq!(tracks[0].last().copied(), Some(50));
        assert_eq!(tracks[1].first().copied(), Some(51));
    }

    #[test]
    fn bad_grammar_left_hyphen_range() {
        let err: RangeParseError = parse_positions(ReferenceFrame::Left, "1-10").unwrap_err();
        assert!(
            err.to_string().contains("Example: --positions 1..10"),
            "unexpected error: {}",
            err
        );
    }

    #[test]
    fn bad_negative_on_nearest() {
        let err = parse_positions(ReferenceFrame::Nearest, "10..-10").unwrap_err();
        assert!(
            err.to_string().contains("Example: --positions ..half"),
            "unexpected error: {}",
            err
        );
    }
}

#[cfg(test)]
mod tests_ticks {
    use cfdnalab::commands::visualize_positions::{
        Track,
        model::AxisBounds,
        render_ascii::{build_tick_lines, value_to_column},
    };

    #[test]
    fn overlapping_ticks_prefer_endpoint_label() {
        let track = Track {
            name: "test".to_string(),
            axis: AxisBounds::new(1, 121),
            selected_indices: Vec::new(),
        };
        let width = 8;
        let (ticks, labels) = build_tick_lines(&track, width);

        let start = track.axis.start as f64;
        let end = track.axis.end as f64;
        let end_column = value_to_column(end, start, end, width);
        assert_eq!(ticks.chars().nth(end_column), Some('|'));

        let end_label_len = track.axis.end.to_string().len();
        let end_start = end_column.saturating_sub(end_label_len.saturating_sub(1));
        let end_label: String = labels.chars().skip(end_start).take(end_label_len).collect();
        assert_eq!(end_label, track.axis.end.to_string());
        assert!(!labels.contains("1121"));
    }
}

mod tests_visualize_positions_config {
    use cfdnalab::commands::cli_common::FragmentPositionSelectionArgs;
    use cfdnalab::commands::visualize_positions::config::VisualizePositionsConfig;
    use cfdnalab::commands::visualize_positions::{
        BasesFrom, MismatchBasesFrom, ReferenceFrame, Style,
    };

    #[test]
    fn build_uses_expected_defaults() {
        let cfg = VisualizePositionsConfig {
            position_selection: FragmentPositionSelectionArgs {
                frame: ReferenceFrame::Left,
                positions: "1..5".to_string(),
                step: 1,
                bases_from: BasesFrom::PreferReads,
                mismatch_bases_from: MismatchBasesFrom::NearestRead,
            },
            lengths: Some(vec![120]),
            length_range: None,
            kmer_sizes: None,
            style: Style::Ascii,
            width: None,
            height: None,
            output: None,
            label: None,
            hide_index: false,
            show_half: false,
            hide_mid: false,
        };

        let viz = cfg.build().expect("config should build");
        assert_eq!(viz.frame, ReferenceFrame::Left);
        assert_eq!(viz.positions_input, "1..5");
        assert_eq!(viz.step.get(), 1);
        assert_eq!(viz.bases, BasesFrom::PreferReads);
        assert_eq!(viz.mismatch_bases_from, MismatchBasesFrom::NearestRead);
        assert!(viz.show_mid);
        assert_eq!(viz.fragment_lengths, vec![120]);
        assert_eq!(viz.style, Style::Ascii);
        assert_eq!(viz.width, 100);
        assert_eq!(viz.height, 120);
    }

    #[test]
    fn build_applies_overrides() {
        let cfg = VisualizePositionsConfig {
            position_selection: FragmentPositionSelectionArgs {
                frame: ReferenceFrame::Right,
                positions: "..half".to_string(),
                step: 3,
                bases_from: BasesFrom::NearestRead,
                mismatch_bases_from: MismatchBasesFrom::BaseQuality,
            },
            lengths: Some(vec![90, 120]),
            length_range: None,
            kmer_sizes: None,
            style: Style::Svg,
            width: Some(140),
            height: Some(200),
            output: None,
            label: Some("test".to_string()),
            hide_index: false,
            show_half: true,
            hide_mid: true,
        };

        let viz = cfg.build().expect("config should build");
        assert_eq!(viz.frame, ReferenceFrame::Right);
        assert_eq!(viz.step.get(), 3);
        assert_eq!(viz.bases, BasesFrom::NearestRead);
        assert_eq!(viz.mismatch_bases_from, MismatchBasesFrom::BaseQuality);
        assert!(!viz.show_mid);
        assert_eq!(viz.fragment_lengths, vec![90, 120]);
        assert_eq!(viz.style, Style::Svg);
        assert_eq!(viz.width, 140);
        assert_eq!(viz.height, 200);
        assert_eq!(viz.label.as_deref(), Some("test"));
        assert!(viz.show_index);
        assert!(viz.show_half);
    }

    #[test]
    fn svg_default_width_is_wider() {
        let cfg = VisualizePositionsConfig {
            position_selection: FragmentPositionSelectionArgs {
                frame: ReferenceFrame::Left,
                positions: "1..5".to_string(),
                step: 1,
                bases_from: BasesFrom::Reference,
                mismatch_bases_from: MismatchBasesFrom::NearestRead,
            },
            lengths: Some(vec![150]),
            length_range: None,
            kmer_sizes: None,
            style: Style::Svg,
            width: None,
            height: None,
            output: None,
            label: None,
            hide_index: true,
            show_half: false,
            hide_mid: false,
        };

        let viz = cfg.build().expect("config should build");
        assert_eq!(viz.width, 650);
    }

    #[test]
    fn build_rejects_zero_step() {
        let cfg = VisualizePositionsConfig {
            position_selection: FragmentPositionSelectionArgs {
                frame: ReferenceFrame::Left,
                positions: "1..5".to_string(),
                step: 0,
                bases_from: BasesFrom::PreferReads,
                mismatch_bases_from: MismatchBasesFrom::NearestRead,
            },
            lengths: Some(vec![100]),
            length_range: None,
            kmer_sizes: None,
            style: Style::Ascii,
            width: None,
            height: None,
            output: None,
            label: None,
            hide_index: true,
            show_half: false,
            hide_mid: false,
        };

        let err = cfg.build().expect_err("step of zero must error");
        assert!(err.to_string().contains("--step must be at least 1"));
    }

    #[test]
    fn build_rejects_fragments_shorter_than_minimum() {
        let cfg = VisualizePositionsConfig {
            position_selection: FragmentPositionSelectionArgs {
                frame: ReferenceFrame::Left,
                positions: "1..5".to_string(),
                step: 1,
                bases_from: BasesFrom::PreferReads,
                mismatch_bases_from: MismatchBasesFrom::NearestRead,
            },
            lengths: Some(vec![9]),
            length_range: None,
            kmer_sizes: None,
            style: Style::Ascii,
            width: None,
            height: None,
            output: None,
            label: None,
            hide_index: true,
            show_half: false,
            hide_mid: false,
        };

        let err = cfg.build().expect_err("length < 10 should fail");
        assert!(err.to_string().contains("10"));
    }

    #[test]
    fn build_rejects_mid_with_nearest_read() {
        let cfg = VisualizePositionsConfig {
            position_selection: FragmentPositionSelectionArgs {
                frame: ReferenceFrame::Mid,
                positions: "-5..5".to_string(),
                step: 1,
                bases_from: BasesFrom::NearestRead,
                mismatch_bases_from: MismatchBasesFrom::NearestRead,
            },
            lengths: Some(vec![100]),
            length_range: None,
            kmer_sizes: None,
            style: Style::Ascii,
            width: None,
            height: None,
            output: None,
            label: None,
            hide_index: true,
            show_half: false,
            hide_mid: false,
        };

        let err = cfg
            .build()
            .expect_err("mid frame with nearest-read should fail");
        assert!(
            err.to_string()
                .contains("`--bases-from nearest-read` is incompatible")
        );
    }
}

mod tests_ascii_render {
    use std::num::NonZeroUsize;

    use cfdnalab::commands::visualize_positions::model::AxisBounds;
    use cfdnalab::commands::visualize_positions::{
        BasesFrom, LengthVisualization, LinearRange, MismatchBasesFrom, PositionsSpec,
        ReferenceFrame, Style, Track, VizConfig, render_ascii,
    };

    fn base_config(width: usize) -> VizConfig {
        VizConfig {
            frame: ReferenceFrame::Left,
            positions: PositionsSpec::Linear(LinearRange::All),
            positions_input: "..".to_string(),
            step: NonZeroUsize::new(1).unwrap(),
            bases: BasesFrom::Reference,
            mismatch_bases_from: MismatchBasesFrom::NearestRead,
            kmer_sizes: None,
            fragment_lengths: vec![100],
            style: Style::Ascii,
            width,
            height: 120,
            output: None,
            label: None,
            show_index: false,
            show_half: false,
            show_mid: true,
        }
    }

    #[test]
    fn fills_contiguous_columns_when_width_small() {
        let track = Track {
            name: "fragment".to_string(),
            axis: AxisBounds::new(1, 100),
            selected_indices: (1..=100).collect(),
        };
        let viz = LengthVisualization {
            fragment_length: 100,
            tracks: vec![track],
        };
        let config = base_config(12);

        let ascii = render_ascii(&[viz], &config);
        let fragment_line = ascii
            .lines()
            .find(|line| line.starts_with("fragment"))
            .expect("missing fragment row");
        let bar = fragment_line
            .split(": ")
            .nth(1)
            .expect("missing bar segment");
        assert!(
            bar.chars().all(|ch| ch == '#'),
            "expected full coverage, got {}",
            bar
        );
    }

    #[test]
    fn preserves_gaps_for_sparse_selection() {
        let track = Track {
            name: "fragment".to_string(),
            axis: AxisBounds::new(1, 100),
            selected_indices: (1..=100).step_by(3).collect(),
        };
        let viz = LengthVisualization {
            fragment_length: 100,
            tracks: vec![track],
        };
        let config = base_config(12);

        let ascii = render_ascii(&[viz], &config);
        let fragment_line = ascii
            .lines()
            .find(|line| line.starts_with("fragment"))
            .expect("missing fragment row");
        let bar = fragment_line
            .split(": ")
            .nth(1)
            .expect("missing bar segment");
        assert!(
            bar.chars().any(|ch| ch == '.'),
            "sparse selection should leave gaps, got {}",
            bar
        );
        assert!(bar.chars().any(|ch| ch == '#'));
    }

    #[test]
    fn nearest_row_includes_max_distance_annotation() {
        let fragment = Track {
            name: "fragment".to_string(),
            axis: AxisBounds::new(1, 100),
            selected_indices: (1..=100).collect(),
        };
        let nearest = Track {
            name: "nearest".to_string(),
            axis: AxisBounds::new(1, 50),
            selected_indices: (1..=50).collect(),
        };
        let viz = LengthVisualization {
            fragment_length: 100,
            tracks: vec![fragment, nearest],
        };
        let mut config = base_config(100);
        config.frame = ReferenceFrame::Nearest;

        let ascii = render_ascii(&[viz], &config);
        assert!(ascii.contains("max distance 50"));
        assert!(ascii.contains("axis(nearest max=50)"));
    }
}
