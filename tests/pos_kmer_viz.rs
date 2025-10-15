use std::num::NonZeroUsize;

use cfdnalab::pos_kmer_viz::{
    PositionsSpec, RangeParseError, ReadClamp, ReferenceFrame, build_tracks_for_length,
    parse_positions,
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
fn stride_application() {
    let spec = parse_positions(ReferenceFrame::Left, "1..10").unwrap();
    let step = NonZeroUsize::new(3).unwrap();
    let tracks = take_linear_indices(20, ReferenceFrame::Left, &spec, step);
    assert_eq!(tracks[0], vec![1, 4, 7, 10]);
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
