use std::cmp::Ordering;

const MIN_LENGTH: usize = 50;
const MAX_LENGTH: usize = 150;
const MAX_RUN_LENGTH: usize = MAX_LENGTH * 3;
const MAX_GAP: u64 = 5;

/// Detected peak with Snyder-style statistics.
#[derive(Clone, Debug)]
pub struct PeakCall {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub peak_position: u64,
    pub height: f32,
}

/// Call peaks on a normalized Snyder residual trace.
///
/// Peaks are detected by following the original Snyder et al. pipeline:
/// positive residual runs are grouped, short gaps are bridged, windows are
/// filtered by length, and the densest sub-region is emitted. The returned peaks
/// carry genomic coordinates (start inclusive, end exclusive) as well as the
/// the position of the maximum height inside the peak region. The `min_peak_height`
/// parameter controls the minimum residual height required to keep a peak.
pub fn call_peaks(
    chr: &str,
    start_offset: u64,
    normalized_wps_values: &[f32],
    mask: &[u8],
    min_peak_height: f32,
) -> Vec<PeakCall> {
    let mut peaks: Vec<PeakCall> = Vec::new();
    let mut positions: Vec<u64> = Vec::new();
    let mut values: Vec<f32> = Vec::new();
    let mut last_positive: Option<u64> = None;

    for (idx, &value) in normalized_wps_values.iter().enumerate() {
        let masked = mask.get(idx).copied().unwrap_or(0) != 0;
        let absolute_pos = start_offset + idx as u64;

        if masked || !value.is_finite() {
            finalize_run(
                chr,
                &mut peaks,
                &mut positions,
                &mut values,
                min_peak_height,
            );
            last_positive = None;
            continue;
        }

        if value > 0.0 {
            if let Some(last) = last_positive {
                let gap = absolute_pos.saturating_sub(last).saturating_sub(1);
                if gap > MAX_GAP {
                    finalize_run(
                        chr,
                        &mut peaks,
                        &mut positions,
                        &mut values,
                        min_peak_height,
                    );
                } else {
                    for step in 1..=gap {
                        positions.push(last + step);
                        values.push(0.0);
                    }
                }
            } else if !positions.is_empty() {
                finalize_run(
                    chr,
                    &mut peaks,
                    &mut positions,
                    &mut values,
                    min_peak_height,
                );
            }

            positions.push(absolute_pos);
            values.push(value);
            last_positive = Some(absolute_pos);
        }
    }

    finalize_run(
        chr,
        &mut peaks,
        &mut positions,
        &mut values,
        min_peak_height,
    );
    peaks
}

/// Finalize the current positive run and emit a peak if one is found.
fn finalize_run(
    chr: &str,
    peaks: &mut Vec<PeakCall>,
    positions: &mut Vec<u64>,
    values: &mut Vec<f32>,
    min_peak_height: f32,
) {
    if positions.is_empty() {
        return;
    }
    let len = positions.len();
    if len < MIN_LENGTH {
        positions.clear();
        values.clear();
        return;
    }
    let peaks_from_run = evaluate_run(
        chr,
        positions.as_slice(),
        values.as_slice(),
        min_peak_height,
    );
    peaks.extend(peaks_from_run);
    positions.clear();
    values.clear();
}

/// Evaluate one positive run and return the corresponding peaks.
fn evaluate_run(
    chr: &str,
    positions: &[u64],
    values: &[f32],
    min_peak_height: f32,
) -> Vec<PeakCall> {
    let len = positions.len();
    if len < MIN_LENGTH {
        return Vec::new();
    }

    if len > MAX_RUN_LENGTH {
        return Vec::new();
    }

    let med = median(values);
    let filtered: Vec<(u64, f32)> = positions
        .iter()
        .zip(values.iter())
        .filter(|(_, v)| **v >= med)
        .map(|(&pos, &val)| (pos, val))
        .collect();

    if filtered.is_empty() {
        return Vec::new();
    }

    let windows = continuous_windows(&filtered);
    let mut peaks = Vec::new();

    if len <= MAX_LENGTH {
        if let Some(best) = windows.into_iter().max_by(|a, b| compare_sum(a.sum, b.sum)) {
            if best.max_value > min_peak_height && best.length() >= MIN_LENGTH {
                peaks.push(PeakCall {
                    chromosome: chr.to_string(),
                    start: best.start,
                    end: best.end + 1,
                    peak_position: best.peak_position,
                    height: best.max_value,
                });
            }
        }
    } else {
        for window in windows {
            let length = window.length();
            if length >= MIN_LENGTH && length <= MAX_LENGTH && window.max_value > min_peak_height {
                peaks.push(PeakCall {
                    chromosome: chr.to_string(),
                    start: window.start,
                    end: window.end + 1,
                    peak_position: window.peak_position,
                    height: window.max_value,
                });
            }
        }
    }

    peaks
}

/// Compute the median of a slice (copying into a scratch buffer).
fn median(values: &[f32]) -> f32 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) * 0.5
    } else {
        sorted[mid]
    }
}

#[derive(Clone, Copy, Debug)]
struct RunWindow {
    sum: f32,
    start: u64,
    end: u64,
    peak_position: u64,
    max_value: f32,
}

impl RunWindow {
    fn new(pos: u64, value: f32) -> Self {
        Self {
            sum: value,
            start: pos,
            end: pos,
            peak_position: pos,
            max_value: value,
        }
    }

    fn extend(&mut self, pos: u64, value: f32) {
        self.end = pos;
        self.sum += value;
        if value > self.max_value {
            self.max_value = value;
            self.peak_position = pos;
        }
    }

    fn length(&self) -> usize {
        (self.end - self.start + 1) as usize
    }
}

/// Collapse contiguous positive values into windows tracking sum and peak location.
fn continuous_windows(filtered: &[(u64, f32)]) -> Vec<RunWindow> {
    let mut windows = Vec::new();
    let mut iter = filtered.iter();
    if let Some(&(pos, value)) = iter.next() {
        let mut current = RunWindow::new(pos, value);
        for &(p, v) in iter {
            if p == current.end + 1 {
                current.extend(p, v);
            } else {
                windows.push(current);
                current = RunWindow::new(p, v);
            }
        }
        windows.push(current);
    }
    windows
}

#[inline]
fn compare_sum(a: f32, b: f32) -> Ordering {
    a.partial_cmp(&b).unwrap_or(Ordering::Equal)
}
