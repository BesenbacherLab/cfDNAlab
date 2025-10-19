use std::num::NonZeroUsize;

use crate::commands::fragment_kmers::nearest_frame_guard::NearestFrameGuard;
use crate::commands::fragment_kmers::positions::{
    PositionGroup, PositionSelection, PositionSelectionCache,
};
use crate::commands::fragment_kmers::selection::{SelectionDecision, evaluate_selection};

use super::model::{AxisBounds, LengthVisualization, PositionsSpec, ReferenceFrame, Track};

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
    let cache = PositionSelectionCache::new(frame, positions, step, length, length)
        .expect("failed to build position cache");
    let selections = cache.offsets(length).unwrap_or(&[]);

    let mut tracks = match frame {
        ReferenceFrame::Left => {
            let indices = selections
                .iter()
                .filter(|sel| sel.group() == PositionGroup::Left)
                .map(|sel| (sel.offset() + 1) as i32)
                .collect();
            vec![Track {
                name: "left".to_string(),
                axis: AxisBounds::new(1, length as i32),
                selected_indices: dedup_sorted(indices),
            }]
        }
        ReferenceFrame::Right => {
            let indices = selections
                .iter()
                .filter(|sel| sel.group() == PositionGroup::Right)
                .map(|sel| (length as i32) - sel.offset() as i32)
                .collect();
            vec![Track {
                name: "right".to_string(),
                axis: AxisBounds::new(1, length as i32),
                selected_indices: dedup_sorted(indices),
            }]
        }
        ReferenceFrame::PerEnd => {
            let left = selections
                .iter()
                .filter(|sel| sel.group() == PositionGroup::Left)
                .map(|sel| (sel.offset() + 1) as i32)
                .collect();
            let right = selections
                .iter()
                .filter(|sel| sel.group() == PositionGroup::Right)
                .map(|sel| (length as i32) - sel.offset() as i32)
                .collect();
            vec![
                Track {
                    name: "left".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: dedup_sorted(left),
                },
                Track {
                    name: "right".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: dedup_sorted(right),
                },
            ]
        }
        ReferenceFrame::Nearest => {
            let fragment_positions = dedup_sorted(
                selections
                    .iter()
                    .map(|sel| (sel.offset() + 1) as i32)
                    .collect(),
            );
            let distances = fold_fragment_positions(length, &fragment_positions);
            vec![
                Track {
                    name: "fragment".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: fragment_positions,
                },
                Track {
                    name: "nearest".to_string(),
                    axis: AxisBounds::new(1, (length / 2).max(1) as i32),
                    selected_indices: distances,
                },
            ]
        }
        ReferenceFrame::Mid => {
            let center = (length as i64) / 2;
            let indices = selections
                .iter()
                .filter(|sel| sel.group() == PositionGroup::Mid)
                .map(|sel| sel.offset() as i64 - center)
                .filter(|val| *val >= i64::from(i32::MIN) && *val <= i64::from(i32::MAX))
                .map(|val| val as i32)
                .collect();
            vec![Track {
                name: "mid".to_string(),
                axis: mid_axis_bounds(length),
                selected_indices: dedup_sorted(indices),
            }]
        }
    };

    apply_read_clamp(&mut tracks, frame, length, read_clamp);

    LengthVisualization {
        fragment_length: length,
        tracks,
    }
}

/// Build helper tracks that illustrate the valid k-mer start bases for the requested kmer_sizes.
pub fn build_kmer_start_overlays(
    frame: ReferenceFrame,
    length: u32,
    positions: &PositionsSpec,
    step: NonZeroUsize,
    base_tracks: &[Track],
    kmer_sizes: &[u8],
) -> Vec<Track> {
    if length == 0 || kmer_sizes.is_empty() || base_tracks.is_empty() {
        return Vec::new();
    }

    let cache = match PositionSelectionCache::new(frame, positions, step, length, length) {
        Ok(cache) => cache,
        Err(_) => return Vec::new(),
    };
    let selections = cache.offsets(length).unwrap_or(&[]);

    match frame {
        ReferenceFrame::Left => build_left_overlays(length, selections, base_tracks, kmer_sizes),
        ReferenceFrame::Right => build_right_overlays(length, selections, base_tracks, kmer_sizes),
        ReferenceFrame::PerEnd => {
            build_per_end_overlays(length, selections, base_tracks, kmer_sizes)
        }
        ReferenceFrame::Nearest => {
            build_nearest_overlays(length, selections, base_tracks, kmer_sizes)
        }
        ReferenceFrame::Mid => build_mid_overlays(length, selections, base_tracks, kmer_sizes),
    }
}
fn build_left_overlays(
    length: u32,
    selections: &[PositionSelection],
    base_tracks: &[Track],
    kmer_sizes: &[u8],
) -> Vec<Track> {
    let base = match base_tracks.iter().find(|track| track.name == "left") {
        Some(track) => track,
        None => return Vec::new(),
    };
    let mut overlays = Vec::new();
    for &k in kmer_sizes {
        let mut overlay = base.clone();
        overlay.name = format!("{} k-mer starts (k={})", base.name, k);
        overlay.selected_indices = dedup_sorted(left_kmer_starts(length, selections, k));
        clamp_overlay_axis(&mut overlay, length, k);
        overlays.push(overlay);
    }
    overlays
}

fn build_right_overlays(
    length: u32,
    selections: &[PositionSelection],
    base_tracks: &[Track],
    kmer_sizes: &[u8],
) -> Vec<Track> {
    let base = match base_tracks.iter().find(|track| track.name == "right") {
        Some(track) => track,
        None => return Vec::new(),
    };
    let mut overlays = Vec::new();
    for &k in kmer_sizes {
        let mut overlay = base.clone();
        overlay.name = format!("{} k-mer starts (k={})", base.name, k);
        overlay.selected_indices = dedup_sorted(right_kmer_starts(length, selections, k));
        clamp_overlay_axis(&mut overlay, length, k);
        overlays.push(overlay);
    }
    overlays
}

fn build_per_end_overlays(
    length: u32,
    selections: &[PositionSelection],
    base_tracks: &[Track],
    kmer_sizes: &[u8],
) -> Vec<Track> {
    let left_base = match base_tracks.iter().find(|track| track.name == "left") {
        Some(track) => track,
        None => return Vec::new(),
    };
    let right_base = match base_tracks.iter().find(|track| track.name == "right") {
        Some(track) => track,
        None => return Vec::new(),
    };
    let mut overlays = Vec::new();
    for &k in kmer_sizes {
        let mut left_overlay = left_base.clone();
        left_overlay.name = format!("{} k-mer starts (k={})", left_base.name, k);
        left_overlay.selected_indices = dedup_sorted(left_kmer_starts(length, selections, k));
        clamp_overlay_axis(&mut left_overlay, length, k);
        overlays.push(left_overlay);

        let mut right_overlay = right_base.clone();
        right_overlay.name = format!("{} k-mer starts (k={})", right_base.name, k);
        right_overlay.selected_indices = dedup_sorted(right_kmer_starts(length, selections, k));
        clamp_overlay_axis(&mut right_overlay, length, k);
        overlays.push(right_overlay);
    }
    overlays
}

fn build_nearest_overlays(
    length: u32,
    selections: &[PositionSelection],
    base_tracks: &[Track],
    kmer_sizes: &[u8],
) -> Vec<Track> {
    let fragment_base = match base_tracks.iter().find(|track| track.name == "fragment") {
        Some(track) => track,
        None => return Vec::new(),
    };
    let nearest_base = match base_tracks.iter().find(|track| track.name == "nearest") {
        Some(track) => track,
        None => return Vec::new(),
    };

    let mut overlays = Vec::new();
    for &k in kmer_sizes {
        let fragment_starts = dedup_sorted(nearest_fragment_kmer_starts(length, selections, k));

        let mut fragment_overlay = fragment_base.clone();
        fragment_overlay.name = format!("{} k-mer starts (k={})", fragment_base.name, k);
        fragment_overlay.selected_indices = fragment_starts.clone();
        clamp_overlay_axis(&mut fragment_overlay, length, k);
        overlays.push(fragment_overlay);

        let mut nearest_overlay = nearest_base.clone();
        nearest_overlay.name = format!("{} k-mer starts (k={})", nearest_base.name, k);
        nearest_overlay.selected_indices = fold_fragment_positions(length, &fragment_starts);
        overlays.push(nearest_overlay);
    }
    overlays
}

fn build_mid_overlays(
    length: u32,
    selections: &[PositionSelection],
    base_tracks: &[Track],
    kmer_sizes: &[u8],
) -> Vec<Track> {
    let base = match base_tracks.iter().find(|track| track.name == "mid") {
        Some(track) => track,
        None => return Vec::new(),
    };
    let mut overlays = Vec::new();
    for &k in kmer_sizes {
        let mut overlay = base.clone();
        overlay.name = format!("{} k-mer starts (k={})", base.name, k);
        overlay.selected_indices = dedup_sorted(mid_kmer_starts(length, selections, k));
        overlays.push(overlay);
    }
    overlays
}

fn clamp_overlay_axis(overlay: &mut Track, length: u32, k: u8) {
    if let Some(max_start) = length
        .checked_sub(u32::from(k))
        .and_then(|value| value.checked_add(1))
    {
        let max_start_i32 = if max_start > i32::MAX as u32 {
            i32::MAX
        } else {
            max_start as i32
        };
        overlay.axis.end = overlay.axis.end.min(max_start_i32);
    }
}

fn default_ranges(length: u32, k_len: u32) -> (Option<(u64, u64)>, Option<(u64, u64)>) {
    if length == 0 {
        return (None, None);
    }
    let forward_max = length.saturating_sub(k_len) as u64;
    let forward = Some((0, forward_max));
    let reverse_min = k_len.saturating_sub(1) as u64;
    let reverse_max = length as u64 - 1;
    let reverse = if reverse_min > reverse_max {
        None
    } else {
        Some((reverse_min, reverse_max))
    };
    (forward, reverse)
}

fn left_kmer_starts(length: u32, selections: &[PositionSelection], k: u8) -> Vec<i32> {
    let k_len = u32::from(k);
    let guard = NearestFrameGuard::for_frame(ReferenceFrame::Left, length, k_len);
    let span = k_len as u64;
    let (forward_range, reverse_range) = default_ranges(length, k_len);
    selections
        .iter()
        .filter(|sel| sel.group() == PositionGroup::Left)
        .filter_map(|sel| {
            match evaluate_selection(
                *sel,
                guard.as_ref(),
                span,
                sel.offset() as u64,
                forward_range,
                reverse_range,
            ) {
                SelectionDecision::IncludeForward { start_offset_0 } => {
                    Some((start_offset_0 + 1) as i32)
                }
                _ => None,
            }
        })
        .collect()
}

fn right_kmer_starts(length: u32, selections: &[PositionSelection], k: u8) -> Vec<i32> {
    let k_len = u32::from(k);
    let guard = NearestFrameGuard::for_frame(ReferenceFrame::Right, length, k_len);
    let span = k_len as u64;
    let (forward_range, reverse_range) = default_ranges(length, k_len);
    selections
        .iter()
        .filter(|sel| sel.group() == PositionGroup::Right)
        .filter_map(|sel| {
            match evaluate_selection(
                *sel,
                guard.as_ref(),
                span,
                sel.offset() as u64,
                forward_range,
                reverse_range,
            ) {
                SelectionDecision::IncludeReverse { start_offset_0, .. } => {
                    let value = length as i64 - start_offset_0 as i64 - k_len as i64 + 1;
                    if value >= i64::from(i32::MIN) && value <= i64::from(i32::MAX) {
                        Some(value as i32)
                    } else {
                        None
                    }
                }
                _ => None,
            }
        })
        .collect()
}

fn nearest_fragment_kmer_starts(length: u32, selections: &[PositionSelection], k: u8) -> Vec<i32> {
    let k_len = u32::from(k);
    let guard = NearestFrameGuard::for_frame(ReferenceFrame::Nearest, length, k_len);
    let span = k_len as u64;
    let (forward_range, reverse_range) = default_ranges(length, k_len);
    selections
        .iter()
        .filter_map(|sel| {
            match evaluate_selection(
                *sel,
                guard.as_ref(),
                span,
                sel.offset() as u64,
                forward_range,
                reverse_range,
            ) {
                SelectionDecision::IncludeForward { start_offset_0 }
                | SelectionDecision::IncludeReverse { start_offset_0, .. } => {
                    Some((start_offset_0 + 1) as i32)
                }
                _ => None,
            }
        })
        .collect()
}

fn mid_kmer_starts(length: u32, selections: &[PositionSelection], k: u8) -> Vec<i32> {
    let k_len = u32::from(k);
    let guard = NearestFrameGuard::for_frame(ReferenceFrame::Mid, length, k_len);
    let span = k_len as u64;
    let (forward_range, reverse_range) = default_ranges(length, k_len);
    let center = (length as i64) / 2;
    selections
        .iter()
        .filter(|sel| sel.group() == PositionGroup::Mid)
        .filter_map(|sel| {
            match evaluate_selection(
                *sel,
                guard.as_ref(),
                span,
                sel.offset() as u64,
                forward_range,
                reverse_range,
            ) {
                SelectionDecision::IncludeForward { start_offset_0 } => {
                    let value = start_offset_0 as i64 - center;
                    if value >= i64::from(i32::MIN) && value <= i64::from(i32::MAX) {
                        Some(value as i32)
                    } else {
                        None
                    }
                }
                _ => None,
            }
        })
        .collect()
}

fn fold_fragment_positions(length: u32, starts: &[i32]) -> Vec<i32> {
    if length == 0 {
        return Vec::new();
    }
    let half = length / 2;
    let mut distances = Vec::with_capacity(starts.len());
    for &start in starts {
        if start <= 0 {
            continue;
        }
        let start_u32 = start as u32;
        let distance = if start_u32 <= half {
            start_u32
        } else {
            length - start_u32 + 1
        };
        if distance > 0 {
            distances.push(distance as i32);
        }
    }
    distances.sort_unstable();
    distances.dedup();
    distances
}

fn dedup_sorted(mut values: Vec<i32>) -> Vec<i32> {
    if values.len() <= 1 {
        return values;
    }
    values.sort_unstable();
    values.dedup();
    values
}

fn mid_axis_bounds(length: u32) -> AxisBounds {
    let half = (length / 2) as i32;
    if length % 2 == 0 {
        AxisBounds::new(-half, half - 1)
    } else {
        AxisBounds::new(-half, half)
    }
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
