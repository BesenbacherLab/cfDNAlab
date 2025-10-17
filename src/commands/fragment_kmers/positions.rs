use std::num::NonZeroUsize;

use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};

use crate::commands::visualize_positions::select::{ReadClamp, build_tracks_for_length};
use crate::commands::visualize_positions::{PositionsSpec, ReferenceFrame, Track};

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
    min_length: u32,
    offsets: Vec<Vec<PositionSelection>>,
}

impl PositionSelectionCache {
    /// Build a cache of selected offsets (0-based) for every fragment length in `[min_len, max_len]`.
    pub fn new(
        frame: ReferenceFrame,
        positions: &PositionsSpec,
        step: NonZeroUsize,
        min_len: u32,
        max_len: u32,
    ) -> Result<Self> {
        if min_len > max_len {
            bail!("min fragment length {min_len} exceeds max fragment length {max_len}");
        }

        let mut offsets: Vec<Vec<PositionSelection>> =
            Vec::with_capacity((max_len - min_len + 1) as usize);
        for length in min_len..=max_len {
            let mut values = offsets_for_length(length, frame, positions, step)?;
            values.sort_unstable_by(|a, b| match a.offset.cmp(&b.offset) {
                std::cmp::Ordering::Equal => match orientation_order(a.orientation())
                    .cmp(&orientation_order(b.orientation()))
                {
                    std::cmp::Ordering::Equal => group_order(a.group).cmp(&group_order(b.group)),
                    other => other,
                },
                other => other,
            });
            values.dedup();
            offsets.push(values);
        }

        Ok(Self {
            min_length: min_len,
            offsets,
        })
    }

    #[inline]
    pub fn offsets(&self, length: u32) -> Option<&[PositionSelection]> {
        if length < self.min_length {
            return None;
        }
        let idx = (length - self.min_length) as usize;
        self.offsets.get(idx).map(|v| v.as_slice())
    }

    #[inline]
    pub fn bounds(&self, length: u32) -> Option<(u32, u32)> {
        let selections = self.offsets(length)?;
        if selections.is_empty() {
            None
        } else {
            let first = selections.first().unwrap().offset;
            let last = selections.last().unwrap().offset;
            Some((first, last))
        }
    }
}

fn offsets_for_length(
    length: u32,
    frame: ReferenceFrame,
    positions: &PositionsSpec,
    step: NonZeroUsize,
) -> Result<Vec<PositionSelection>> {
    let viz = build_tracks_for_length(length, frame, positions, step, ReadClamp::None);
    let mut offsets = Vec::new();
    match frame {
        ReferenceFrame::Left => {
            if let Some(track) = viz.tracks.first() {
                offsets.extend(linear_offsets(track, length, false, PositionGroup::Left)?);
            }
        }
        ReferenceFrame::Right => {
            if let Some(track) = viz.tracks.first() {
                offsets.extend(linear_offsets(track, length, true, PositionGroup::Right)?);
            }
        }
        ReferenceFrame::PerEnd => {
            for track in viz.tracks.iter() {
                let is_right = track.name.eq_ignore_ascii_case("right");
                let group = if is_right {
                    PositionGroup::Right
                } else {
                    PositionGroup::Left
                };
                offsets.extend(linear_offsets(track, length, is_right, group)?);
            }
        }
        ReferenceFrame::Nearest => {
            for track in viz.tracks.iter() {
                if track.name == "fragment" {
                    offsets.extend(nearest_offsets(track, length)?);
                }
            }
        }
        ReferenceFrame::Mid => {
            for track in viz.tracks.iter() {
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
