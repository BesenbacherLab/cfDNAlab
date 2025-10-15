use std::num::NonZeroUsize;

use super::model::{
    AxisBounds, LengthVisualization, LinearRange, MidRange, NearestRange, PositionsSpec,
    ReferenceFrame, Track,
};

/// How aggressively the visualization should clamp selections to read coverage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadClamp {
    /// Keep every position described by the frame.
    None,
    /// Keep only the positions covered by the frame's nearest read.
    Nearest,
    /// Keep only the positions covered by either read.
    Both,
}

/// Build the set of tracks for a single fragment length.
pub fn build_tracks_for_length(
    length: u32,
    frame: ReferenceFrame,
    positions: &PositionsSpec,
    step: NonZeroUsize,
    read_clamp: ReadClamp,
) -> LengthVisualization {
    let tracks = match frame {
        ReferenceFrame::Left => vec![build_linear_track(
            "left",
            length,
            expect_linear(positions),
            step,
        )],
        ReferenceFrame::Right => vec![build_linear_track(
            "right",
            length,
            expect_linear(positions),
            step,
        )],
        ReferenceFrame::PerEnd => {
            let range = expect_linear(positions);
            vec![
                build_linear_track("left", length, range, step),
                build_linear_track("right", length, range, step),
            ]
        }
        ReferenceFrame::Nearest => build_nearest_tracks(length, expect_nearest(positions), step),
        ReferenceFrame::Mid => vec![build_mid_track(length, expect_mid(positions), step)],
    };

    let mut tracks = tracks;
    apply_read_clamp(&mut tracks, frame, length, read_clamp);

    LengthVisualization {
        fragment_length: length,
        tracks,
    }
}

fn build_linear_track(label: &str, length: u32, range: &LinearRange, step: NonZeroUsize) -> Track {
    let axis_end = length as i32;
    let axis = AxisBounds::new(1, axis_end);
    let indices = collect_linear_indices(length, range);
    Track {
        name: label.to_string(),
        axis,
        selected_indices: apply_stride(indices, step),
    }
}

fn collect_linear_indices(length: u32, range: &LinearRange) -> Vec<i32> {
    let axis_start = 1i32;
    let axis_end = length as i32;
    match *range {
        LinearRange::All => inclusive_range(axis_start, axis_end),
        LinearRange::Closed { start, end } => {
            clamp_range_to_domain(start as i32, end as i32, axis_start, axis_end)
                .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
        }
        LinearRange::From { start } => {
            clamp_range_to_domain(start as i32, axis_end, axis_start, axis_end)
                .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
        }
        LinearRange::To { end } => {
            clamp_range_to_domain(axis_start, end as i32, axis_start, axis_end)
                .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
        }
        LinearRange::TrimOtherEnd {
            start,
            other_end_trim,
        } => {
            let raw_end = length as i64 - other_end_trim as i64;
            if raw_end < 1 {
                Vec::new()
            } else {
                let end = raw_end.min(length as i64) as i32;
                clamp_range_to_domain(start as i32, end, axis_start, axis_end)
                    .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
            }
        }
        LinearRange::ToHalf { minus } => {
            let half = length / 2;
            if half == 0 {
                Vec::new()
            } else {
                let high = half.saturating_sub(minus) as i32;
                if high <= 0 {
                    Vec::new()
                } else {
                    clamp_range_to_domain(axis_start, high, axis_start, axis_end)
                        .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
                }
            }
        }
        LinearRange::FromToHalf { start, minus } => {
            let half = length / 2;
            if half == 0 {
                Vec::new()
            } else {
                let high = half.saturating_sub(minus) as i32;
                if high <= 0 {
                    Vec::new()
                } else {
                    clamp_range_to_domain(start as i32, high, axis_start, axis_end)
                        .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
                }
            }
        }
    }
}

fn build_nearest_tracks(length: u32, range: &NearestRange, step: NonZeroUsize) -> Vec<Track> {
    let half = (length / 2) as u32;
    let axis_end = if half == 0 { 1 } else { half as i32 };
    let folded = collect_nearest_indices(half, range);
    let folded = apply_stride(folded, step);
    let nearest_track = Track {
        name: "nearest".to_string(),
        axis: AxisBounds::new(1, axis_end),
        selected_indices: folded.clone(),
    };

    let fragment_track = Track {
        name: "fragment".to_string(),
        axis: AxisBounds::new(1, length as i32),
        selected_indices: unfold_nearest_indices(length, &folded),
    };

    vec![fragment_track, nearest_track]
}

fn collect_nearest_indices(half: u32, range: &NearestRange) -> Vec<i32> {
    if half == 0 {
        return Vec::new();
    }
    let axis_start = 1i32;
    let axis_end = half as i32;
    match *range {
        NearestRange::All => inclusive_range(axis_start, axis_end),
        NearestRange::Closed { start, end } => {
            clamp_range_to_domain(start as i32, end as i32, axis_start, axis_end)
                .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
        }
        NearestRange::From { start } => {
            clamp_range_to_domain(start as i32, axis_end, axis_start, axis_end)
                .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
        }
        NearestRange::ToHalf { minus } => {
            let high = half.saturating_sub(minus);
            if high == 0 {
                Vec::new()
            } else {
                clamp_range_to_domain(axis_start, high as i32, axis_start, axis_end)
                    .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
            }
        }
        NearestRange::FromToHalf { start, minus } => {
            let high = half.saturating_sub(minus);
            if high == 0 {
                Vec::new()
            } else {
                clamp_range_to_domain(start as i32, high as i32, axis_start, axis_end)
                    .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
            }
        }
    }
}

fn build_mid_track(length: u32, range: &MidRange, step: NonZeroUsize) -> Track {
    let axis = mid_axis_bounds(length);
    let indices = collect_mid_indices(length, range);
    Track {
        name: "mid".to_string(),
        axis,
        selected_indices: apply_mid_stride(indices, step),
    }
}

fn mid_axis_bounds(length: u32) -> AxisBounds {
    let half = (length / 2) as i32;
    if length % 2 == 0 {
        AxisBounds::new(-half, half - 1)
    } else {
        AxisBounds::new(-half, half)
    }
}

fn collect_mid_indices(length: u32, range: &MidRange) -> Vec<i32> {
    let axis = mid_axis_bounds(length);
    match *range {
        MidRange::All => inclusive_range(axis.start, axis.end),
        MidRange::Closed { neg, pos } => {
            let start = -(neg as i32);
            let end = pos as i32;
            clamp_range_to_domain(start, end, axis.start, axis.end)
                .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
        }
        MidRange::LeftOpen { neg } => {
            let start = -(neg as i32);
            clamp_range_to_domain(start, 0, axis.start, axis.end)
                .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
        }
        MidRange::RightOpen { pos } => {
            let end = pos as i32;
            clamp_range_to_domain(0, end, axis.start, axis.end)
                .map_or_else(Vec::new, |(s, e)| inclusive_range(s, e))
        }
    }
}

fn apply_mid_stride(indices: Vec<i32>, step: NonZeroUsize) -> Vec<i32> {
    if indices.len() <= 1 || step.get() == 1 {
        return indices;
    }

    if let Some(origin_idx) = indices.iter().position(|&value| value == 0) {
        let origin_idx = origin_idx as i64;
        let step_span = step.get() as i64;
        return indices
            .into_iter()
            .enumerate()
            .filter_map(|(idx, value)| {
                let idx = idx as i64;
                if (idx - origin_idx).rem_euclid(step_span) == 0 {
                    Some(value)
                } else {
                    None
                }
            })
            .collect();
    }

    apply_stride(indices, step)
}

fn apply_stride(indices: Vec<i32>, step: NonZeroUsize) -> Vec<i32> {
    if indices.len() <= 1 || step.get() == 1 {
        return indices;
    }
    indices
        .into_iter()
        .enumerate()
        .filter_map(|(idx, value)| {
            if idx % step.get() == 0 {
                Some(value)
            } else {
                None
            }
        })
        .collect()
}

fn inclusive_range(start: i32, end: i32) -> Vec<i32> {
    if start > end {
        Vec::new()
    } else {
        (start..=end).collect()
    }
}

fn clamp_range_to_domain(
    mut start: i32,
    mut end: i32,
    domain_start: i32,
    domain_end: i32,
) -> Option<(i32, i32)> {
    if domain_start > domain_end {
        return None;
    }
    if end < domain_start || start > domain_end {
        return None;
    }
    start = start.max(domain_start);
    end = end.min(domain_end);
    if start > end {
        None
    } else {
        Some((start, end))
    }
}

fn expect_linear(positions: &PositionsSpec) -> &LinearRange {
    match positions {
        PositionsSpec::Linear(range) => range,
        _ => panic!("expected linear range for linear frame"),
    }
}

fn expect_nearest(positions: &PositionsSpec) -> &NearestRange {
    match positions {
        PositionsSpec::Nearest(range) => range,
        _ => panic!("expected nearest range for nearest frame"),
    }
}

fn expect_mid(positions: &PositionsSpec) -> &MidRange {
    match positions {
        PositionsSpec::Mid(range) => range,
        _ => panic!("expected mid range for mid frame"),
    }
}

fn unfold_nearest_indices(length: u32, folded: &[i32]) -> Vec<i32> {
    if length == 0 || folded.is_empty() {
        return Vec::new();
    }
    let mut positions = Vec::with_capacity(folded.len() * 2);
    let max_pos = length as i32;
    for &distance in folded {
        if distance <= 0 {
            continue;
        }
        let left = distance;
        let right = max_pos - distance + 1;
        if (1..=max_pos).contains(&left) {
            positions.push(left);
        }
        if (1..=max_pos).contains(&right) {
            positions.push(right);
        }
    }
    positions.sort_unstable();
    positions.dedup();
    positions
}

fn apply_read_clamp(
    tracks: &mut [Track],
    frame: ReferenceFrame,
    length: u32,
    read_clamp: ReadClamp,
) {
    if matches!(read_clamp, ReadClamp::None) || length == 0 {
        return;
    }

    let half = ((length + 1) / 2) as i32;
    let right_start = (length as i32 + 1) - half;

    for track in tracks {
        match read_clamp {
            ReadClamp::None => {}
            ReadClamp::Nearest => clamp_track_nearest(track, frame, half, right_start),
            ReadClamp::Both => clamp_track_both_reads(track, frame, half, right_start),
        }
    }
}

fn clamp_track_nearest(track: &mut Track, frame: ReferenceFrame, half: i32, right_start: i32) {
    match frame {
        ReferenceFrame::Left => track.selected_indices.retain(|&idx| idx <= half),
        ReferenceFrame::Right => track.selected_indices.retain(|&idx| idx >= right_start),
        ReferenceFrame::PerEnd => {
            if track.name == "left" {
                track.selected_indices.retain(|&idx| idx <= half);
            } else if track.name == "right" {
                track.selected_indices.retain(|&idx| idx >= right_start);
            }
        }
        ReferenceFrame::Nearest => {
            if track.name == "fragment" {
                track
                    .selected_indices
                    .retain(|&idx| idx <= half || idx >= right_start);
            }
        }
        ReferenceFrame::Mid => {
            track.selected_indices.retain(|&idx| idx.abs() <= half);
        }
    }
}

fn clamp_track_both_reads(track: &mut Track, frame: ReferenceFrame, half: i32, right_start: i32) {
    match frame {
        ReferenceFrame::Left | ReferenceFrame::Right => {
            track
                .selected_indices
                .retain(|&idx| idx <= half || idx >= right_start);
        }
        ReferenceFrame::PerEnd => {
            if track.name == "left" {
                track.selected_indices.retain(|&idx| idx <= half);
            } else if track.name == "right" {
                track.selected_indices.retain(|&idx| idx >= right_start);
            }
        }
        ReferenceFrame::Nearest => {
            if track.name == "fragment" {
                track
                    .selected_indices
                    .retain(|&idx| idx <= half || idx >= right_start);
            }
        }
        ReferenceFrame::Mid => {
            track.selected_indices.retain(|&idx| idx.abs() <= half);
        }
    }
}
