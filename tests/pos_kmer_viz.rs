use std::num::NonZeroUsize;

use cfdnalab::pos_kmer_viz::{
    Anchor, PositionsSpec, RangeParseError, build_tracks_for_length, parse_positions,
};

fn default_step() -> NonZeroUsize {
    NonZeroUsize::new(1).unwrap()
}

fn take_linear_indices(
    length: u32,
    anchor: Anchor,
    positions: &PositionsSpec,
    step: NonZeroUsize,
) -> Vec<Vec<i32>> {
    let viz = build_tracks_for_length(length, anchor, positions, step);
    viz.tracks
        .iter()
        .map(|track| track.selected_indices.clone())
        .collect()
}

#[test]
fn nearest_open_to_half_small_l() {
    let spec = parse_positions(Anchor::Nearest, "10..").unwrap();
    let tracks = take_linear_indices(18, Anchor::Nearest, &spec, default_step());
    assert!(tracks[0].is_empty());
}

#[test]
fn nearest_half_minus_k() {
    let spec = parse_positions(Anchor::Nearest, "5..half-3").unwrap();
    let tracks = take_linear_indices(151, Anchor::Nearest, &spec, default_step());
    let expected: Vec<i32> = (5..=72).collect();
    assert_eq!(tracks[0], expected);
}

#[test]
fn left_opposite_end_bound() {
    let spec = parse_positions(Anchor::Left, "10..-10").unwrap();
    let tracks = take_linear_indices(100, Anchor::Left, &spec, default_step());
    let expected: Vec<i32> = (10..=90).collect();
    assert_eq!(tracks[0], expected);
}

#[test]
fn right_opposite_end_bound() {
    let spec = parse_positions(Anchor::Right, "10..-10").unwrap();
    let tracks = take_linear_indices(101, Anchor::Right, &spec, default_step());
    let expected: Vec<i32> = (10..=91).collect();
    assert_eq!(tracks[0], expected);
}

#[test]
fn per_end_two_tracks() {
    let spec = parse_positions(Anchor::PerEnd, "..5").unwrap();
    let tracks = take_linear_indices(120, Anchor::PerEnd, &spec, default_step());
    let expected: Vec<i32> = (1..=5).collect();
    assert_eq!(tracks.len(), 2);
    assert_eq!(tracks[0], expected);
    assert_eq!(tracks[1], expected);
}

#[test]
fn span_trim_both_ends() {
    let spec = parse_positions(Anchor::Span, "15..-15").unwrap();
    let tracks = take_linear_indices(80, Anchor::Span, &spec, default_step());
    let expected: Vec<i32> = (15..=65).collect();
    assert_eq!(tracks[0], expected);
}

#[test]
fn mid_symmetric_closed() {
    let spec = parse_positions(Anchor::Mid, "-10..10").unwrap();
    let tracks = take_linear_indices(99, Anchor::Mid, &spec, default_step());
    let expected: Vec<i32> = (-10..=10).collect();
    assert_eq!(tracks[0], expected);
}

#[test]
fn mid_open_right() {
    let spec = parse_positions(Anchor::Mid, "..5").unwrap();
    let tracks = take_linear_indices(150, Anchor::Mid, &spec, default_step());
    let expected: Vec<i32> = (0..=5).collect();
    assert_eq!(tracks[0], expected);

    let legacy = parse_positions(Anchor::Mid, "..+5").unwrap();
    let legacy_tracks = take_linear_indices(150, Anchor::Mid, &legacy, default_step());
    assert_eq!(legacy_tracks[0], expected);
}

#[test]
fn should_keep_origin_when_mid_stride_applied() {
    let spec = parse_positions(Anchor::Mid, "-6..6").unwrap();
    let step = NonZeroUsize::new(3).unwrap();
    let tracks = take_linear_indices(101, Anchor::Mid, &spec, step);
    assert_eq!(tracks[0], vec![-6, -3, 0, 3, 6]);
}

#[test]
fn stride_application() {
    let spec = parse_positions(Anchor::Left, "1..10").unwrap();
    let step = NonZeroUsize::new(3).unwrap();
    let tracks = take_linear_indices(20, Anchor::Left, &spec, step);
    assert_eq!(tracks[0], vec![1, 4, 7, 10]);
}

#[test]
fn nearest_center_double_count_guard() {
    let spec = parse_positions(Anchor::Nearest, "..half").unwrap();
    let tracks = take_linear_indices(100, Anchor::Nearest, &spec, default_step());
    assert_eq!(tracks[0].last().copied(), Some(50));
    assert_eq!(tracks[0].iter().filter(|&&v| v == 50).count(), 1);
}

#[test]
fn bad_grammar_left_hyphen_range() {
    let err: RangeParseError = parse_positions(Anchor::Left, "1-10").unwrap_err();
    assert!(
        err.to_string().contains("Example: --positions 1..10"),
        "unexpected error: {}",
        err
    );
}

#[test]
fn bad_negative_on_nearest() {
    let err = parse_positions(Anchor::Nearest, "10..-10").unwrap_err();
    assert!(
        err.to_string().contains("Example: --positions ..half"),
        "unexpected error: {}",
        err
    );
}
