//! Smoothing and baseline normalization helpers for Snyder-style WPS peaks.
//!
//! The functions in this module follow the narrative described in
//! `peak_calling_logic.md` and prioritize readability for newcomers.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::ops::Range;

use staged_sg_filter::sav_gol_f32;

/// Savitzky-Golay half window size (10 -> 21 bp total window).
const SG_HALF_WINDOW: usize = 10;
/// Full Savitzky-Golay window size used for smoothing (21 bp).
const SG_WINDOW_SIZE: usize = 2 * SG_HALF_WINDOW + 1;
/// Polynomial order preserved by the Savitzky-Golay kernel.
const SG_ORDER: usize = 2;

/// Normalize WPS values by subtracting a rolling median baseline.
///
/// Centers per-base WPS values around a local median so short-range spikes are easy to spot.
///
/// Maintains a sliding median over `baseline_reference`, honors optional mask
/// entries, and emits `f32::NAN` when the window lacks enough usable bases.
/// Feeding the raw WPS as `baseline_reference` and the smoothed values as
/// `numerator` reproduces Snyder's behavior.
///
/// Parameters
/// ----------
/// - `numerator`:
///   Per-position values after optional smoothing. Masked positions are ignored
///   and produce `NaN` in the output.
///
/// - `baseline_reference`:
///   Values used to compute the rolling median. Defaults to the raw WPS so the
///   baseline is not affected by smoothing.
///
/// - `mask`:
///   Optional mask that uses `0` for usable bases and `1` for excluded bases.
///   Any non-zero value is treated as masked. Masked entries are removed from
///   the median and force the output at the same index to `NaN`.
///
/// - `window_size`:
///   Number of bases inside the median window. Use an odd number to keep the
///   center aligned.
///
/// - `stride`:
///   Distance between centers. Use `1` to mirror Snyder's evaluation at every base.
///
/// - `min_unmasked`:
///   Minimum count of usable bases required to publish a value. Windows below
///   the threshold yield `NaN`.
///
/// Returns
/// -------
/// - `Vec<f32>`:
///   Baseline-adjusted signal ready for peak detection.
pub fn normalize_wps(
    numerator: &[f32],
    baseline_reference: &[f32],
    mask: Option<&[u8]>,
    window_size: usize,
    stride: usize,
    min_unmasked: usize,
) -> Vec<f32> {
    assert_eq!(
        numerator.len(),
        baseline_reference.len(),
        "numerator and baseline reference must have identical lengths"
    );
    if let Some(mask_slice) = mask {
        assert_eq!(
            numerator.len(),
            mask_slice.len(),
            "mask length must match the WPS series length"
        );
    }
    assert!(window_size > 0, "window size must be strictly positive");
    assert!(
        stride == 1,
        "normalization currently only supports stride == 1"
    );

    let len = numerator.len();
    if len == 0 {
        return Vec::new();
    }

    // Masking may use different sentinels (blacklist, edges)
    // so treat any non-zero byte as masked
    let mask_slice = mask.unwrap_or(&[]);
    let use_mask = !mask_slice.is_empty();

    let left_span = window_size / 2;
    let right_span = window_size.saturating_sub(left_span);
    let required = min_unmasked.min(window_size);

    let mut result = vec![f32::NAN; len];
    let mut median = SlidingMedian::new(len);

    // Tracks which indices are currently represented in the sliding median
    let mut active_range = Range {
        start: 0usize,
        end: 0usize,
    };

    for center in 0..len {
        let window_start = center.saturating_sub(left_span);
        let window_end = (center + right_span).min(len);

        // Extend window to include newly covered positions on the right
        while active_range.end < window_end {
            let idx = active_range.end;
            if value_usable(idx, baseline_reference, use_mask.then_some(mask_slice)) {
                median.insert(idx, baseline_reference[idx]);
            }
            active_range.end += 1;
        }

        // Shrink window to drop positions that are no longer covered
        while active_range.start < window_start {
            let idx = active_range.start;
            if value_usable(idx, baseline_reference, use_mask.then_some(mask_slice)) {
                median.remove(idx);
            }
            active_range.start += 1;
        }

        // Ignore masked centers or non-finite input values right away
        let center_masked = use_mask && mask_slice[center] != 0;
        if center_masked || !numerator[center].is_finite() {
            result[center] = f32::NAN;
            continue;
        }

        // Require enough unmasked bases before trusting the median
        let window_population = median.count();
        if window_population < required {
            result[center] = f32::NAN;
            continue;
        }

        let window_median = median
            .median()
            .expect("median unavailable despite population check");
        result[center] = numerator[center] - window_median;
    }

    result
}

/// Smooth the WPS signal with a 21 bp, order-2 Savitzky-Golay filter.
///
/// Mirrors the behavior from Snyder et al. by padding each contiguous unmasked
/// segment with reflected values. Masked bases act as hard edges and remain
/// `NaN` in the returned slice.
///
/// Parameters
/// ----------
/// - `wps_values`:
///   Raw WPS signal for a tile. Values can be finite floats or `NaN`.
///
/// - `mask`:
///   Optional mask matching the input length. Values are expected to be `0`
///   for usable bases and `1` for masked positions (blacklist or dilated
///   edges). Any non-zero value is treated as masked, so future sentinels are
///   handled as well. Masked entries partition the signal into independent
///   segments.
///
/// Returns
/// -------
/// - `Vec<f32>`:
///   Smoothed values aligned to the input positions.
pub fn smoothe_wps(wps_values: &[f32], mask: Option<&[u8]>) -> Vec<f32> {
    let len = wps_values.len();
    if let Some(mask_slice) = mask {
        assert_eq!(
            mask_slice.len(),
            len,
            "mask length must match the input WPS length"
        );
    }

    if len == 0 {
        return Vec::new();
    }

    let mask_slice = mask.unwrap_or(&[]);
    let use_mask = !mask_slice.is_empty();

    let mut smoothed = vec![f32::NAN; len];
    let mut idx = 0usize;

    while idx < len {
        if use_mask && mask_slice[idx] != 0 {
            // Skip masked bases entirely, output stays NaN for them
            idx += 1;
            continue;
        }
        // Remember where the current unmasked run begins
        let segment_start = idx;
        while idx < len && (!use_mask || mask_slice[idx] == 0) {
            idx += 1;
        }
        let segment_end = idx;
        let segment = &wps_values[segment_start..segment_end];
        // Smooth the current unmasked stretch in isolation
        apply_snyder_smoothing(segment, segment_start, &mut smoothed);
    }

    smoothed
}

/// Apply the Snyder-style smoothing to a contiguous unmasked segment.
///
/// Smooth each continuous stretch and write the results back into the
/// full-length buffer.
///
/// Mirrors the segment around its edges so the Savitzky-Golay window receives
/// enough context. Short segments therefore fall back to reflected padding
/// instead of producing `NaN`. Writes the result directly into `output` at the
/// provided offset.
///
/// Parameters
/// ----------
/// - `segment`:
///   Portion of the WPS signal without masked bases.
///
/// - `offset`:
///   Starting index in the original sequence where the segment begins.
///
/// - `output`:
///   Buffer receiving the smoothed values for the full tile.
fn apply_snyder_smoothing(segment: &[f32], offset: usize, output: &mut [f32]) {
    let seg_len = segment.len();
    if seg_len == 0 {
        return;
    }
    // Three scenarios to handle below:
    // 1) The filter window fits entirely inside the segment (no padding)
    // 2) Close to the left edge, where we mirror upcoming values
    // 3) Close to the right edge, where we mirror preceding values

    for local_idx in 0..seg_len {
        let absolute_idx = offset + local_idx;
        let smoothed = if seg_len >= SG_WINDOW_SIZE
            && local_idx >= SG_HALF_WINDOW
            && local_idx + SG_HALF_WINDOW < seg_len
        {
            // Case 1: full window fits inside the segment; use values as-is
            let window = &segment[local_idx - SG_HALF_WINDOW..=local_idx + SG_HALF_WINDOW];
            apply_savgol(window)
        } else if local_idx < SG_HALF_WINDOW {
            // Case 2: window would read past the left edge, mirror upcoming values
            let right_limit = (local_idx + SG_HALF_WINDOW + 1).min(seg_len);
            let edge_slice = &segment[local_idx..right_limit];
            let padded = build_left_edge_window(edge_slice);
            apply_savgol(&padded)
        } else {
            // Case 3: window would read past the right edge, mirror previous values
            let left_start = local_idx.saturating_sub(SG_HALF_WINDOW);
            let edge_slice = &segment[left_start..];
            let padded = build_right_edge_window(edge_slice);
            apply_savgol(&padded)
        };
        output[absolute_idx] = smoothed;
    }
}

/// Construct a mirrored window for the left edge of a segment.
pub fn build_left_edge_window(edge_slice: &[f32]) -> Vec<f32> {
    debug_assert!(!edge_slice.is_empty());
    let needed = SG_WINDOW_SIZE.saturating_sub(edge_slice.len());
    let base = edge_slice[0];
    let available = edge_slice.len().saturating_sub(1);
    let take = needed.min(available);

    let mut window = Vec::with_capacity(SG_WINDOW_SIZE);
    for &value in edge_slice[1..1 + take].iter().rev() {
        // Mirror the slice around the first base so the filter sees a smooth ramp
        window.push(base - (value - base).abs());
    }
    if take < needed {
        window.extend(std::iter::repeat(base).take(needed - take));
    }
    window.extend_from_slice(edge_slice);
    trim_or_pad(window, base)
}

/// Construct a mirrored window for the right edge of a segment.
pub fn build_right_edge_window(edge_slice: &[f32]) -> Vec<f32> {
    debug_assert!(!edge_slice.is_empty());
    let needed = SG_WINDOW_SIZE.saturating_sub(edge_slice.len());
    let base = edge_slice[edge_slice.len() - 1];
    let available = edge_slice.len().saturating_sub(1);
    let take = needed.min(available);
    let start = edge_slice.len().saturating_sub(take + 1);

    let mut window = edge_slice.to_vec();
    for &value in edge_slice[start..edge_slice.len() - 1].iter().rev() {
        // Mirror the slice around the last base in ascending order
        window.push(base + (value - base).abs());
    }
    if take < needed {
        window.extend(std::iter::repeat(base).take(needed - take));
    }
    trim_or_pad(window, base)
}

/// Ensure a mirrored edge window matches the Savitzky-Golay length.
fn trim_or_pad(mut window: Vec<f32>, fill: f32) -> Vec<f32> {
    if window.len() > SG_WINDOW_SIZE {
        window.truncate(SG_WINDOW_SIZE);
    } else if window.len() < SG_WINDOW_SIZE {
        window.extend(std::iter::repeat(fill).take(SG_WINDOW_SIZE - window.len()));
    }
    window
}

/// Run the Savitzky-Golay filter and return the centered value.
#[inline]
fn apply_savgol(window: &[f32]) -> f32 {
    debug_assert_eq!(window.len(), SG_WINDOW_SIZE);
    let mut buf = window.to_vec();
    sav_gol_f32::<SG_HALF_WINDOW, SG_ORDER>(&mut buf, window);
    buf[SG_HALF_WINDOW]
}

/// Check whether a position contributes to the rolling median.
///
/// Treats any non-zero mask entry as excluded so both blacklist and edge
/// markers are handled uniformly.
fn value_usable(idx: usize, values: &[f32], mask: Option<&[u8]>) -> bool {
    values[idx].is_finite() && mask.map_or(true, |m| m[idx] == 0)
}

/// Identify which heap currently owns a median entry.
#[derive(Clone, Copy)]
enum HeapSide {
    Lower,
    Upper,
}

/// Node stored in one of the median heaps (value plus its original index).
#[derive(Clone, Copy)]
struct MedianHeapEntry {
    value: OrderedF32,
    index: usize,
}

/// Float wrapper that provides total ordering for heap storage.
#[derive(Clone, Copy, Debug, PartialEq)]
struct OrderedF32(f32);

impl Eq for OrderedF32 {}

impl PartialEq for MedianHeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index && self.value == other.value
    }
}

impl Eq for MedianHeapEntry {}

impl PartialOrd for MedianHeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MedianHeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.value
            .cmp(&other.value)
            .then_with(|| self.index.cmp(&other.index))
    }
}

impl PartialOrd for OrderedF32 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(
            self.0
                .partial_cmp(&other.0)
                .unwrap_or_else(|| panic!("encountered NaN in rolling median")),
        )
    }
}

impl Ord for OrderedF32 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other)
            .expect("NaN encountered in rolling median")
    }
}

/// Sliding median with lazy deletion for Snyder-style baseline subtraction.
///
/// Tracks the median of a moving window without resorting the entire slice at
/// every step. Conceptually the window is split into two piles: the lower half
/// (a max-heap) and the upper half (a min-heap). Each update keeps the piles
/// balanced and drops stale entries just-in-time.
pub struct SlidingMedian {
    lower: BinaryHeap<MedianHeapEntry>,
    upper: BinaryHeap<std::cmp::Reverse<MedianHeapEntry>>,
    active: Vec<bool>,
    side: Vec<Option<HeapSide>>,
    lower_size: usize,
    upper_size: usize,
}

impl SlidingMedian {
    /// Create an empty sliding median.
    ///
    /// Allocates the internal arrays up to `capacity` so later inserts and
    /// removals remain `O(log n)`.
    ///
    /// Parameters
    /// ----------
    /// - `capacity`:
    ///   Maximum index that may appear in the sliding window.
    ///
    /// Returns
    /// -------
    /// - `Self`:
    ///   Ready-to-use median structure with empty heaps.
    pub fn new(capacity: usize) -> Self {
        Self {
            lower: BinaryHeap::new(),
            upper: BinaryHeap::new(),
            active: vec![false; capacity],
            side: vec![None; capacity],
            lower_size: 0,
            upper_size: 0,
        }
    }

    /// Add a new value to the sliding window.
    ///
    /// Chooses the heap based on the current maximum of the lower half, tags
    /// the index as active, and calls `rebalance` to maintain the size
    /// invariant.
    ///
    /// Parameters
    /// ----------
    /// - `index`:
    ///   Absolute position of the sample within the global series.
    ///
    /// - `value`:
    ///   Measured WPS value to insert.
    pub fn insert(&mut self, index: usize, value: f32) {
        let entry = MedianHeapEntry {
            value: OrderedF32(value),
            index,
        };
        if self.lower_size == 0 || self.lower.peek().map_or(true, |top| entry <= *top) {
            // Either the heaps are empty or this value belongs to the lower pile
            self.lower.push(entry);
            self.side[index] = Some(HeapSide::Lower);
            self.lower_size += 1;
        } else {
            // Otherwise place it into the upper pile (min-heap)
            self.upper.push(std::cmp::Reverse(entry));
            self.side[index] = Some(HeapSide::Upper);
            self.upper_size += 1;
        }
        self.active[index] = true;
        self.rebalance();
    }

    /// Remove a value from the sliding window if it is still active.
    ///
    /// Marks the index as inactive, adjusts the cached heap sizes, then
    /// rebalances the heaps to flush stale entries.
    ///
    /// Parameters
    /// ----------
    /// - `index`:
    ///   Sample index previously inserted via `insert`.
    pub fn remove(&mut self, index: usize) {
        if !self.active.get(index).copied().unwrap_or(false) {
            return;
        }
        self.active[index] = false;
        match self.side[index].take() {
            Some(HeapSide::Lower) => {
                self.lower_size = self.lower_size.saturating_sub(1);
            }
            Some(HeapSide::Upper) => {
                self.upper_size = self.upper_size.saturating_sub(1);
            }
            None => {}
        }
        self.rebalance();
    }

    /// Return the current median if the window is non-empty.
    ///
    /// Prunes stale heap heads, then either averages the
    /// two middle values or returns the max of the lower heap depending on the
    /// parity of the window size.
    ///
    /// Returns
    /// -------
    /// - `Option<f32>`:
    ///   `Some(median)` when at least one active value exists, otherwise `None`.
    pub fn median(&mut self) -> Option<f32> {
        self.rebalance();
        if self.lower_size + self.upper_size == 0 {
            return None;
        }
        if self.lower_size == self.upper_size {
            let left = self.lower.peek()?.value.0;
            let right = self.upper.peek()?.0.value.0;
            // Even population: average the middle pair
            Some((left + right) * 0.5)
        } else {
            // Odd population: the top of the lower heap is the median
            Some(self.lower.peek()?.value.0)
        }
    }

    /// Return the number of active values contributing to the median.
    ///
    /// Returns
    /// -------
    /// - `usize`:
    ///   Count of values currently represented in the heaps.
    pub fn count(&self) -> usize {
        self.lower_size + self.upper_size
        // Both heaps exclude inactive entries from the size counters
    }

    /// Restore the heap size invariant and prune stale nodes.
    ///
    /// Moves entries between heaps until their sizes differ by at most one
    /// while discarding inactive entries.
    pub fn rebalance(&mut self) {
        self.prune();
        while self.lower_size > self.upper_size + 1 {
            if let Some(entry) = self.lower.pop() {
                if !self.active[entry.index] {
                    continue;
                }
                // Move the largest value from the lower pile to the upper pile
                self.upper.push(std::cmp::Reverse(entry));
                self.side[entry.index] = Some(HeapSide::Upper);
                self.lower_size -= 1;
                self.upper_size += 1;
            } else {
                break;
            }
        }
        while self.lower_size < self.upper_size {
            if let Some(std::cmp::Reverse(entry)) = self.upper.pop() {
                if !self.active[entry.index] {
                    continue;
                }
                // Move the smallest value from the upper pile back to the lower pile
                self.lower.push(entry);
                self.side[entry.index] = Some(HeapSide::Lower);
                self.upper_size -= 1;
                self.lower_size += 1;
            } else {
                break;
            }
        }
    }

    /// Drop entries that were lazily marked as inactive.
    ///
    /// Pops from each heap until the front entry belongs to
    /// an active index.
    pub fn prune(&mut self) {
        while let Some(entry) = self.lower.peek() {
            if self.active[entry.index] {
                break;
            }
            // The head belongs to a stale index; drop it
            self.lower.pop();
        }
        while let Some(std::cmp::Reverse(entry)) = self.upper.peek() {
            if self.active[entry.index] {
                break;
            }
            // Same idea for the upper heap
            self.upper.pop();
        }
    }
}
