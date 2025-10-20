use fxhash::FxHashMap;

use super::model::{AxisBounds, LengthVisualization, Track};
use crate::commands::fragment_kmers::nearest_frame_guard::NearestFrameGuard;
use crate::commands::fragment_kmers::parse::PositionalSelectionSpec;
use crate::commands::fragment_kmers::positions::{
    PositionGroup, PositionSelection, PositionSelectionCache, ReferenceFrame,
};
use crate::commands::fragment_kmers::selection::{SelectionDecision, evaluate_selection};

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
    position_selection_spec: &PositionalSelectionSpec,
    read_clamp: ReadClamp,
) -> LengthVisualization {
    let cache = PositionSelectionCache::new(
        vec![position_selection_spec.clone()],
        &[0u8],
        length,
        length,
    )
    .expect("failed to build position cache");
    let selections = cache.offsets(length, 1u8).unwrap_or(&[]);

    let mut tracks = match position_selection_spec.frame {
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
            let left_positions = dedup_sorted(
                selections
                    .iter()
                    .filter(|sel| sel.group() == PositionGroup::Left)
                    .filter_map(|sel| {
                        map_linear_position(length, sel.offset(), PositionGroup::Left)
                    })
                    .collect(),
            );
            let right_positions = dedup_sorted(
                selections
                    .iter()
                    .filter(|sel| sel.group() == PositionGroup::Right)
                    .filter_map(|sel| {
                        map_linear_position(length, sel.offset(), PositionGroup::Right)
                    })
                    .collect(),
            );
            let fragment_positions = dedup_sorted({
                let mut values = left_positions.clone();
                values.extend(right_positions.iter().copied());
                values
            });
            let distances = fold_fragment_positions(length, &fragment_positions);
            let mut tracks = vec![
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
            ];
            if !left_positions.is_empty() {
                tracks.push(Track {
                    name: "left".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: left_positions,
                });
            }
            if !right_positions.is_empty() {
                tracks.push(Track {
                    name: "right".to_string(),
                    axis: AxisBounds::new(1, length as i32),
                    selected_indices: right_positions,
                });
            }
            tracks
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

    apply_read_clamp(
        &mut tracks,
        position_selection_spec.frame,
        length,
        read_clamp,
    );

    LengthVisualization {
        fragment_length: length,
        tracks,
    }
}

/// Build helper tracks that illustrate the valid k-mer start bases for the requested kmer_sizes.
pub fn build_kmer_start_overlays(
    position_selection_spec: &PositionalSelectionSpec,
    length: u32,
    base_tracks: &[Track],
    kmer_sizes: &[u8],
) -> FxHashMap<u8, Vec<Track>> {
    let mut overlay_per_k: FxHashMap<u8, Vec<Track>> = FxHashMap::default();
    if length == 0 || kmer_sizes.is_empty() || base_tracks.is_empty() {
        return overlay_per_k;
    }

    let cache = match PositionSelectionCache::new(
        vec![position_selection_spec.clone()],
        kmer_sizes,
        length,
        length,
    ) {
        Ok(cache) => cache,
        Err(_) => return overlay_per_k,
    };

    for k in kmer_sizes {
        let selections = cache.offsets(length, *k).unwrap_or(&[]);

        let overlay = match position_selection_spec.frame {
            ReferenceFrame::Left => {
                build_left_overlays(length, selections, base_tracks, kmer_sizes)
            }
            ReferenceFrame::Right => {
                build_right_overlays(length, selections, base_tracks, kmer_sizes)
            }
            ReferenceFrame::PerEnd => {
                build_per_end_overlays(length, selections, base_tracks, kmer_sizes)
            }
            ReferenceFrame::Nearest => {
                build_nearest_overlays(length, selections, base_tracks, kmer_sizes)
            }
            ReferenceFrame::Mid => build_mid_overlays(length, selections, base_tracks, kmer_sizes),
        };
        overlay_per_k.insert(*k, overlay);
    }

    overlay_per_k
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
    let left_base = base_tracks.iter().find(|track| track.name == "left");
    let right_base = base_tracks.iter().find(|track| track.name == "right");

    let mut overlays = Vec::new();
    for &k in kmer_sizes {
        let (fragment_starts, left_starts, right_starts) =
            nearest_fragment_kmer_starts(length, selections, k);
        if fragment_starts.is_empty() && left_starts.is_empty() && right_starts.is_empty() {
            continue;
        }

        let mut fragment_overlay = fragment_base.clone();
        fragment_overlay.name = format!("{} k-mer starts (k={})", fragment_base.name, k);
        fragment_overlay.selected_indices = fragment_starts.clone();
        overlays.push(fragment_overlay);

        let mut nearest_overlay = nearest_base.clone();
        nearest_overlay.name = format!("{} k-mer starts (k={})", nearest_base.name, k);
        nearest_overlay.selected_indices = fold_fragment_positions(length, &fragment_starts);
        overlays.push(nearest_overlay);

        if let Some(base) = left_base {
            if !left_starts.is_empty() {
                let mut overlay = base.clone();
                overlay.name = format!("{} k-mer starts (k={})", base.name, k);
                overlay.selected_indices = left_starts;
                overlays.push(overlay);
            }
        }
        if let Some(base) = right_base {
            if !right_starts.is_empty() {
                let mut overlay = base.clone();
                overlay.name = format!("{} k-mer starts (k={})", base.name, k);
                overlay.selected_indices = right_starts;
                overlays.push(overlay);
            }
        }
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
    let span = k_len as u64;
    let (forward_range, reverse_range) = default_ranges(length, k_len);
    selections
        .iter()
        .filter(|sel| sel.group() == PositionGroup::Left)
        .filter_map(|sel| {
            match evaluate_selection(
                *sel,
                None,
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
    let span = k_len as u64;
    let (forward_range, reverse_range) = default_ranges(length, k_len);
    selections
        .iter()
        .filter(|sel| sel.group() == PositionGroup::Right)
        .filter_map(|sel| {
            match evaluate_selection(
                *sel,
                None,
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

fn nearest_fragment_kmer_starts(
    length: u32,
    selections: &[PositionSelection],
    k: u8,
) -> (Vec<i32>, Vec<i32>, Vec<i32>) {
    let k_len = u32::from(k);
    let guard = NearestFrameGuard::by_flag(true, length, k_len);
    let span = k_len as u64;
    let (forward_range, reverse_range) = default_ranges(length, k_len);
    let mut left = Vec::new();
    let mut right = Vec::new();
    for sel in selections {
        let start = match evaluate_selection(
            *sel,
            guard.as_ref(),
            span,
            sel.offset() as u64,
            forward_range,
            reverse_range,
        ) {
            SelectionDecision::IncludeForward { start_offset_0 } => {
                map_linear_position(length, start_offset_0 as u32, PositionGroup::Left)
            }
            SelectionDecision::IncludeReverse { start_offset_0, .. } => {
                let value = start_offset_0.saturating_add(1);
                if value == 0 || value > u64::from(length) {
                    None
                } else {
                    Some(value as i32)
                }
            }
            _ => None,
        };
        if let Some(value) = start {
            if sel.group() == PositionGroup::Right {
                right.push(value);
            } else {
                left.push(value);
            }
        }
    }
    let left = dedup_sorted(left);
    let right = dedup_sorted(right);
    let mut fragment = left.clone();
    fragment.extend(right.iter().copied());
    let fragment = dedup_sorted(fragment);
    (fragment, left, right)
}

fn mid_kmer_starts(length: u32, selections: &[PositionSelection], k: u8) -> Vec<i32> {
    let k_len = u32::from(k);
    let span = k_len as u64;
    let (forward_range, reverse_range) = default_ranges(length, k_len);
    let center = (length as i64) / 2;
    selections
        .iter()
        .filter(|sel| sel.group() == PositionGroup::Mid)
        .filter_map(|sel| {
            match evaluate_selection(
                *sel,
                None,
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

fn map_linear_position(length: u32, offset: u32, group: PositionGroup) -> Option<i32> {
    match group {
        PositionGroup::Left | PositionGroup::Right => {
            let value = offset.saturating_add(1);
            if value == 0 || value > length {
                None
            } else {
                Some(value as i32)
            }
        }
        PositionGroup::Mid => None,
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
        ReferenceFrame::Nearest => {
            if track.name == "fragment" {
                if let Some(&idx) = track
                    .selected_indices
                    .iter()
                    .find(|&&idx| !(idx <= half || idx >= right_start))
                {
                    panic!(
                        "Nearest-read clamp detected nearest fragment track index {} outside <= {} or >= {}.",
                        idx, half, right_start
                    );
                }
                track
                    .selected_indices
                    .retain(|&idx| idx <= half || idx >= right_start);
            } else if track.name == "left" {
                if let Some(&idx) = track.selected_indices.iter().find(|&&idx| idx > half) {
                    panic!(
                        "Nearest-read clamp detected nearest left track index {} outside <= {}.",
                        idx, half
                    );
                }
                track.selected_indices.retain(|&idx| idx <= half);
            } else if track.name == "right" {
                if let Some(&idx) = track
                    .selected_indices
                    .iter()
                    .find(|&&idx| idx < right_start)
                {
                    panic!(
                        "Nearest-read clamp detected nearest right track index {} outside >= {}.",
                        idx, right_start
                    );
                }
                track.selected_indices.retain(|&idx| idx >= right_start);
            }
        }
        ReferenceFrame::Mid => {
            track.selected_indices.retain(|&idx| idx.abs() <= half);
        }
        ReferenceFrame::Left => {
            track.selected_indices.retain(|&idx| idx <= half);
        }
        ReferenceFrame::Right => {
            track.selected_indices.retain(|&idx| idx >= right_start);
        }
        ReferenceFrame::PerEnd => {
            if track.name == "left" {
                track.selected_indices.retain(|&idx| idx <= half);
            } else if track.name == "right" {
                track.selected_indices.retain(|&idx| idx >= right_start);
            }
        }
    }
}

fn clamp_track_both_reads(track: &mut Track, frame: ReferenceFrame, half: i32, right_start: i32) {
    match frame {
        ReferenceFrame::Nearest => {
            if track.name == "fragment" {
                if let Some(&idx) = track
                    .selected_indices
                    .iter()
                    .find(|&&idx| !(idx <= half || idx >= right_start))
                {
                    panic!(
                        "Both-read clamp detected nearest fragment track index {} outside <= {} or >= {}.",
                        idx, half, right_start
                    );
                }
                track
                    .selected_indices
                    .retain(|&idx| idx <= half || idx >= right_start);
            } else if track.name == "left" {
                if let Some(&idx) = track.selected_indices.iter().find(|&&idx| idx > half) {
                    panic!(
                        "Both-read clamp detected nearest left track index {} outside <= {}.",
                        idx, half
                    );
                }
                track.selected_indices.retain(|&idx| idx <= half);
            } else if track.name == "right" {
                if let Some(&idx) = track
                    .selected_indices
                    .iter()
                    .find(|&&idx| idx < right_start)
                {
                    panic!(
                        "Both-read clamp detected nearest right track index {} outside >= {}.",
                        idx, right_start
                    );
                }
                track.selected_indices.retain(|&idx| idx >= right_start);
            }
        }
        ReferenceFrame::Mid => {
            track.selected_indices.retain(|&idx| idx.abs() <= half);
        }
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
    }
}
