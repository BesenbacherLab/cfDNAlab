//! Tile-core difference-array helpers for overlapping fragment length model fitting.

use crate::shared::interval::Interval;

/// Add one covered segment to tile-core difference arrays.
///
/// The segment is clipped to the core before raw depth and fragment length sums are updated.
///
/// Returning `true` means at least one base of this segment belongs to the tile core.
pub(crate) fn add_segment_to_core_deltas(
    segment: Interval<u32>,
    core: Interval<u32>,
    fragment_length: i64,
    raw_depth_delta: &mut [i64],
    fragment_length_sum_delta: &mut [i64],
) -> bool {
    let Some(clipped) = segment.clip_to(core) else {
        return false;
    };
    let start = (clipped.start() - core.start()) as usize;
    let end = (clipped.end() - core.start()) as usize;
    raw_depth_delta[start] += 1;
    raw_depth_delta[end] -= 1;
    fragment_length_sum_delta[start] += fragment_length;
    fragment_length_sum_delta[end] -= fragment_length;
    true
}

#[cfg(test)]
mod tests {
    include!("tiling_tests.rs");
}
