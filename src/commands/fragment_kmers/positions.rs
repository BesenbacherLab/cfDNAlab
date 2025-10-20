use anyhow::{Result, bail};
use fxhash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::num::NonZeroUsize;
use strum_macros::{EnumCount as EnumCountMacro, EnumIter};

#[cfg(feature = "cli")]
use clap::ValueEnum;

use crate::commands::fragment_kmers::parse::PositionalSelectionSpec;
use crate::commands::visualize_positions::Track;
use crate::commands::visualize_positions::model::AxisBounds;
use crate::shared::kmers::kmer_codec::KmerSpec;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PositionOrientation {
    Forward,
    Reverse,
}

impl PositionOrientation {
    /// Get PositionOrientation from PositionGroup
    ///
    /// Left/Mid => forward, right => reverse
    pub fn from_position_group(group: PositionGroup) -> PositionOrientation {
        match group {
            PositionGroup::Left | PositionGroup::Mid => PositionOrientation::Forward,
            PositionGroup::Right => PositionOrientation::Reverse,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum PositionGroup {
    Left,
    Right,
    Mid,
}

/// The position specification tagged to allow frame-specific dispatch.
/// Each variant has range and step size
#[derive(Debug, Clone)]
pub enum PositionsSpec {
    Linear(LinearRange),
    Nearest(NearestRange),
    Mid(MidRange),
}

/// Range grammar used with the `nearest` frame.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NearestRange {
    /// Entire folded axis.
    All,
    Closed {
        start: u32,
        end: u32,
    },
    From {
        start: u32,
    },
    ToHalf {
        minus: u32,
    },
    FromToHalf {
        start: u32,
        minus: u32,
    },
}

/// Range grammar used with the `mid` frame.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MidRange {
    /// Entire symmetric axis.
    All,
    Closed {
        neg: u32,
        pos: u32,
    },
    LeftOpen {
        neg: u32,
    },
    RightOpen {
        pos: u32,
    },
}

/// Range grammar for frames that index strictly from one end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LinearRange {
    /// Entire fragment axis.
    All,
    /// Closed inclusive range `A..B`.
    Closed { start: u32, end: u32 },
    /// Open-right range `A:`.
    From { start: u32 },
    /// Open-left range `:B`.
    To { end: u32 },
    /// Opposite-end trimmed range `A..-B`.
    TrimOtherEnd { start: u32, other_end_trim: u32 },
    /// Range reaching up to the fragment half (`..half[-K]`).
    ToHalf { minus: u32 },
    /// Range starting at a fixed offset up to the fragment half (`A..half[-K]`).
    FromToHalf { start: u32, minus: u32 },
}

/// How to resolve overlapping read mismatches when choosing read-backed bases.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MismatchBasesFrom {
    #[default]
    NearestRead,
    BaseQuality,
    Reference,
}

impl MismatchBasesFrom {
    pub fn as_str(self) -> &'static str {
        match self {
            MismatchBasesFrom::NearestRead => "nearest-read",
            MismatchBasesFrom::BaseQuality => "base-quality",
            MismatchBasesFrom::Reference => "reference",
        }
    }
}

/// Whether the user wants to reason about read or reference coordinates.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BasesFrom {
    /// Always use reference positions regardless of read coverage.
    #[default]
    Reference,
    /// Prefer observed read coordinates, but fall back to the reference span when a read is missing.
    PreferReads,
    /// Only include positions covered by either read. Inferred mate gaps are skipped.
    Reads,
    /// Clamp to the read nearest to the frame origin.
    NearestRead,
}

impl BasesFrom {
    pub fn as_str(self) -> &'static str {
        match self {
            BasesFrom::Reference => "reference",
            BasesFrom::PreferReads => "prefer-reads",
            BasesFrom::Reads => "reads",
            BasesFrom::NearestRead => "nearest-read",
        }
    }
}

/// Enumeration of the available coordinate frames used to interpret positional selections.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Hash, EnumCountMacro, EnumIter)]
pub enum ReferenceFrame {
    #[default]
    Left,
    Right,
    PerEnd,
    Nearest,
    Mid,
}

impl ReferenceFrame {
    pub fn as_str(self) -> &'static str {
        match self {
            ReferenceFrame::Left => "left",
            ReferenceFrame::Right => "right",
            ReferenceFrame::PerEnd => "per-end",
            ReferenceFrame::Nearest => "nearest",
            ReferenceFrame::Mid => "mid",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PositionSelection {
    offset: u32,
    group: PositionGroup,
}

impl PositionSelection {
    #[inline]
    pub fn offset(&self) -> u32 {
        self.offset
    }

    #[inline]
    pub fn orientation(&self) -> PositionOrientation {
        PositionOrientation::from_position_group(self.group)
    }

    #[inline]
    pub fn group(&self) -> PositionGroup {
        self.group
    }
}

/// Cache of positional selections per fragment length.
pub struct PositionSelectionCache {
    pub min_length: u32,
    pub offsets: FxHashMap<u8, Vec<Vec<PositionSelection>>>,
    pub present_frames: FxHashSet<ReferenceFrame>,
    pub main_frame: ReferenceFrame,
}

impl PositionSelectionCache {
    /// Build a cache of selected offsets (0-based) for every fragment length in `[min_len, max_len]`.
    pub fn new(
        position_specs: Vec<PositionalSelectionSpec>,
        kmer_sizes: &[u8],
        min_len: u32,
        max_len: u32,
    ) -> Result<Self> {
        if min_len > max_len {
            bail!("min fragment length {min_len} exceeds max fragment length {max_len}");
        }

        let mut offsets: FxHashMap<u8, Vec<Vec<PositionSelection>>> =
            FxHashMap::with_capacity_and_hasher(kmer_sizes.len(), Default::default());
        for &k in kmer_sizes {
            offsets.insert(k, Vec::with_capacity((max_len - min_len + 1) as usize));
        }

        for length in min_len..=max_len {
            let per_spec_offsets = position_specs
                .iter()
                .map(|ps| offsets_for_length(length, ps))
                .collect::<Result<Vec<_>, _>>()?;

            for &k in kmer_sizes {
                let first_pairs: Vec<(u32, PositionSelection)> = per_spec_offsets[0]
                    .iter()
                    .filter_map(|ps| start_index_for_selection(k, ps).map(|s| (s, *ps)))
                    .collect();

                let mut offset_sets: Vec<FxHashSet<u32>> = Vec::with_capacity(position_specs.len());
                offset_sets.push(first_pairs.iter().map(|(s, _)| *s).collect());

                // Remaining specs: push their offset sets
                for spec_offsets in &per_spec_offsets[1..] {
                    let set_start_values: FxHashSet<u32> = spec_offsets
                        .iter()
                        .filter_map(|ps| start_index_for_selection(k, ps))
                        .collect();
                    offset_sets.push(set_start_values);
                }

                // Intersect offsets and keep first spec’s PositionSelection for survivors
                let surviving_offsets: FxHashSet<u32> = intersect_all(&offset_sets);

                let mut intersected_values: Vec<PositionSelection> = first_pairs
                    .into_iter()
                    .filter(|(start, _)| surviving_offsets.contains(start))
                    .map(|(_, ps)| ps)
                    .collect();

                // Sort by offset, orientation, group
                intersected_values.sort_unstable_by(|a, b| match a.offset().cmp(&b.offset()) {
                    std::cmp::Ordering::Equal => match orientation_order(a.orientation())
                        .cmp(&orientation_order(b.orientation()))
                    {
                        std::cmp::Ordering::Equal => {
                            group_order(a.group()).cmp(&group_order(b.group()))
                        }
                        other => other,
                    },
                    other => other,
                });

                let k_offsets = offsets.get_mut(&k).unwrap();
                k_offsets.push(intersected_values);
            }
        }

        let main_frame = position_specs[0].frame;
        let present_frames: FxHashSet<ReferenceFrame> =
            position_specs.iter().map(|ps| ps.frame).collect();

        Ok(Self {
            min_length: min_len,
            offsets,
            present_frames,
            main_frame,
        })
    }

    #[inline]
    pub fn offsets(&self, length: u32, k: u8) -> Option<&[PositionSelection]> {
        if length < self.min_length {
            return None;
        }
        let idx = (length - self.min_length) as usize;
        self.offsets.get(&k)?.get(idx).map(|v| v.as_slice())
    }

    #[inline]
    pub fn bounds(&self, length: u32, k: u8) -> Option<(u32, u32)> {
        let selections = self.offsets(length, k)?;
        if selections.is_empty() {
            None
        } else {
            let first = selections.first()?.offset();
            let last = selections.last()?.offset();
            Some((first, last))
        }
    }
}

#[inline]
fn start_index_for_selection(k: u8, position_selection: &PositionSelection) -> Option<u32> {
    match position_selection.orientation() {
        // Forward: start is the offset itself
        PositionOrientation::Forward => Some(position_selection.offset()),
        // Reverse: offset is the anchor (last base). Start = anchor - (k-1).
        PositionOrientation::Reverse => {
            let anchor = position_selection.offset();
            let back = k.saturating_sub(1) as u32;
            if anchor < back {
                None
            } else {
                Some(anchor - back)
            }
        }
    }
}

fn offsets_for_length(
    length: u32,
    position_spec: &PositionalSelectionSpec,
) -> Result<Vec<PositionSelection>> {
    let tracks = match position_spec.frame {
        ReferenceFrame::Left => vec![build_linear_track(
            "left",
            length,
            expect_linear(&position_spec.positions),
            position_spec.step,
        )],
        ReferenceFrame::Right => vec![build_linear_track(
            "right",
            length,
            expect_linear(&position_spec.positions),
            position_spec.step,
        )],
        ReferenceFrame::PerEnd => {
            let range = expect_linear(&position_spec.positions);
            vec![
                build_linear_track("left", length, range, position_spec.step),
                build_linear_track("right", length, range, position_spec.step),
            ]
        }
        ReferenceFrame::Nearest => build_nearest_tracks(
            length,
            expect_nearest(&position_spec.positions),
            position_spec.step,
        ),
        ReferenceFrame::Mid => vec![build_mid_track(
            length,
            expect_mid(&position_spec.positions),
            position_spec.step,
        )],
    };

    let mut offsets = Vec::new();
    for track in &tracks {
        match position_spec.frame {
            ReferenceFrame::Left => {
                offsets.extend(linear_offsets(track, length, false, PositionGroup::Left)?);
            }
            ReferenceFrame::Right => {
                offsets.extend(linear_offsets(track, length, true, PositionGroup::Right)?);
            }
            ReferenceFrame::PerEnd => {
                let is_right = track.name.eq_ignore_ascii_case("right");
                let group = if is_right {
                    PositionGroup::Right
                } else {
                    PositionGroup::Left
                };
                offsets.extend(linear_offsets(track, length, is_right, group)?);
            }
            ReferenceFrame::Nearest => {
                if track.name == "fragment" {
                    offsets.extend(nearest_offsets(track, length)?);
                }
            }
            ReferenceFrame::Mid => {
                offsets.extend(mid_offsets(track, length)?);
            }
        }
    }

    Ok(offsets)
}

fn linear_offsets(
    track: &Track,
    length: u32,
    from_right: bool,
    group: PositionGroup,
) -> Result<Vec<PositionSelection>> {
    let mut out = Vec::with_capacity(track.selected_indices.len());
    for &idx in &track.selected_indices {
        if idx <= 0 {
            continue;
        }
        let idx = idx as u32;
        if idx == 0 || idx > length {
            continue;
        }
        let offset = if from_right {
            length
                .checked_sub(idx)
                .ok_or_else(|| anyhow::anyhow!("invalid index {idx} for length {length}"))?
        } else {
            idx - 1
        };
        out.push(PositionSelection { offset, group });
    }
    Ok(out)
}

fn mid_offsets(track: &Track, length: u32) -> Result<Vec<PositionSelection>> {
    let mut out = Vec::with_capacity(track.selected_indices.len());
    let center = (length as i64) / 2;
    for &idx in &track.selected_indices {
        let offset = center + idx as i64;
        if offset < 0 || offset >= length as i64 {
            continue;
        }
        out.push(PositionSelection {
            offset: offset as u32,
            group: PositionGroup::Mid,
        });
    }
    Ok(out)
}

fn nearest_offsets(track: &Track, length: u32) -> Result<Vec<PositionSelection>> {
    let mut out = Vec::with_capacity(track.selected_indices.len());
    let half = length / 2;
    for &idx in &track.selected_indices {
        if idx <= 0 {
            continue;
        }
        let idx = idx as u32;
        if idx == 0 || idx > length {
            continue;
        }
        let orientation = if idx <= half {
            PositionOrientation::Forward
        } else {
            PositionOrientation::Reverse
        };
        out.push(PositionSelection {
            offset: idx - 1,
            group: if orientation == PositionOrientation::Forward {
                PositionGroup::Left
            } else {
                PositionGroup::Right
            },
        });
    }
    Ok(out)
}

#[inline]

fn build_linear_track(label: &str, length: u32, range: &LinearRange, step: NonZeroUsize) -> Track {
    let axis = AxisBounds::new(1, length as i32);
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

fn orientation_order(orientation: PositionOrientation) -> u8 {
    match orientation {
        PositionOrientation::Forward => 0,
        PositionOrientation::Reverse => 1,
    }
}

#[inline]
fn group_order(group: PositionGroup) -> u8 {
    match group {
        PositionGroup::Left => 0,
        PositionGroup::Right => 1,
        PositionGroup::Mid => 2,
    }
}

fn intersect_all<T>(sets: &[FxHashSet<T>]) -> FxHashSet<T>
where
    T: Eq + std::hash::Hash + Clone,
{
    if sets.is_empty() {
        return FxHashSet::default();
    }
    let mut idxs: Vec<usize> = (0..sets.len()).collect();
    idxs.sort_unstable_by_key(|&i| sets[i].len());
    let mut intersection: FxHashSet<T> = sets[idxs[0]].clone();
    for &i in &idxs[1..] {
        let s = &sets[i];
        intersection.retain(|x| s.contains(x));
        if intersection.is_empty() {
            break;
        }
    }
    intersection
}
