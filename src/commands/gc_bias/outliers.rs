use ndarray::{Array2, Axis};
use std::cmp::Ordering;

/// How to detect outliers in a 1D vector (e.g., one length row).
#[derive(Debug, Clone, Copy)]
pub(crate) enum OutlierRule {
    /// No outlier handling.
    None,
    /// Clamp by quantiles, e.g. (0.03, 0.97) for GCfix.
    Quantile { lower: f32, upper: f32 },
    /// Tukey/IQR rule: [Q1 - k*IQR, Q3 + k*IQR]. GCParagon used large k (e.g., 8.0).
    TukeyIqr { k: f32 },
    /// Mean ± k·SD (classical, not robust).
    StdDev { k: f32 },
    /// Median ± k·(1.4826·MAD) (robust to outliers).
    Mad { k: f32 },
}

/// Whether to detect outliers per fragment length row or across the full matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum OutlierScope {
    /// Handle outliers independently per fragment length.
    PerLength,
    /// Handle outliers across the entire matrix.
    Global,
}

/// What to do with detected outliers.
#[derive(Debug, Clone, Copy)]
pub(crate) enum OutlierAction {
    /// Cap to the computed [lower, upper] bounds (winsorization).
    Winsorize,
    /// Replace with NaN (so a NaN-aware smoother can ignore them).
    MaskNaN,
}

/// Summary of outlier handling results.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(crate) struct OutlierStats {
    /// Number of finite, supported values checked against bounds.
    pub(crate) total_examined: usize,
    /// Number of values that were winsorized or masked as outliers.
    pub(crate) total_outliers_handled: usize,
    /// Number of finite values outside the support mask that were checked.
    pub(crate) unsupported_examined: usize,
    /// Number of unsupported values that were winsorized or masked.
    pub(crate) unsupported_outliers_handled: usize,
    /// Number of values adjusted by the hard safety clamp after outlier handling.
    pub(crate) hard_clamped: usize,
}

impl OutlierStats {
    #[inline]
    fn add(&mut self, other: OutlierStats) {
        self.total_examined += other.total_examined;
        self.total_outliers_handled += other.total_outliers_handled;
        self.unsupported_examined += other.unsupported_examined;
        self.unsupported_outliers_handled += other.unsupported_outliers_handled;
        self.hard_clamped += other.hard_clamped;
    }
}

/// Linear quantile interpolation for a sorted slice.
///
/// The percentile is treated as a zero-based position within the sorted values, which can land
/// between two indices.
///
/// Clamps the percentile to [0, 1], computes the position as a float (e.g., 2.4 means 40% of
/// the way from element 2 to 3), and interpolates to avoid stepwise jumps on short inputs.
///
/// Parameters
/// ----------
/// - `sorted`:
///     Values sorted in ascending order to sample from.
/// - `p`:
///     Target percentile expressed as a fraction; values outside [0, 1] are clamped.
///
/// Returns
/// -------
/// - Quantile value with linear interpolation, or NaN when the slice is empty.
#[inline]
pub(crate) fn interpolated_quantile(sorted: &[f32], p: f32) -> f32 {
    let len = sorted.len();
    if len == 0 {
        // Empty slice returns NaN to signal no data
        return f32::NAN;
    }
    // Clamp percentile so callers cannot request past the edges
    let clamped = p.clamp(0.0, 1.0);
    // Fractional position between indices (e.g., 2.4 is 40% of the way from 2 to 3)
    let position = clamped * (len - 1) as f32;
    let lower_idx = position.floor() as usize;
    let upper_idx = position.ceil() as usize;
    if lower_idx == upper_idx {
        return sorted[lower_idx];
    }
    // Linear blend between neighbors
    let blend = position - lower_idx as f32;
    let lower = sorted[lower_idx];
    let upper = sorted[upper_idx];
    lower + (upper - lower) * blend
}

/// Compute outlier bounds for `vals` per the rule. Ignores NaNs.
/// Returns `None` if rule is `None` or there's not enough data.
pub(crate) fn outlier_bounds(vals: &[f32], rule: OutlierRule) -> Option<(f32, f32)> {
    let mut sorted_values: Vec<f32> = vals.iter().copied().filter(|x| x.is_finite()).collect();
    if sorted_values.len() < 2 || matches!(rule, OutlierRule::None) {
        return None;
    }
    sorted_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    match rule {
        OutlierRule::None => None,
        OutlierRule::Quantile { lower, upper } => {
            if !(0.0..=1.0).contains(&lower) || !(0.0..=1.0).contains(&upper) || lower >= upper {
                return None;
            }
            Some((
                interpolated_quantile(&sorted_values, lower),
                interpolated_quantile(&sorted_values, upper),
            ))
        }
        OutlierRule::TukeyIqr { k } => {
            let q1 = interpolated_quantile(&sorted_values, 0.25);
            let q3 = interpolated_quantile(&sorted_values, 0.75);
            let iqr = (q3 - q1).max(0.0);
            Some((q1 - k * iqr, q3 + k * iqr))
        }
        OutlierRule::StdDev { k } => {
            let n = sorted_values.len() as f32;
            let mean = sorted_values.iter().sum::<f32>() / n;
            let var = sorted_values
                .iter()
                .map(|x| (x - mean).powi(2))
                .sum::<f32>()
                / n.max(1.0);
            let sd = var.sqrt();
            Some((mean - k * sd, mean + k * sd))
        }
        OutlierRule::Mad { k } => {
            // Median via interpolated quantile at 0.5 (midpoint when even-sized)
            let median = interpolated_quantile(&sorted_values, 0.5);
            let mut absolute_deviations: Vec<f32> =
                sorted_values.iter().map(|x| (x - median).abs()).collect();
            absolute_deviations.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
            // MAD scaled to be comparable to SD for normal data
            let mad = 1.4826_f32 * interpolated_quantile(&absolute_deviations, 0.5);
            Some((median - k * mad, median + k * mad))
        }
    }
}

/// Apply outlier handling to a correction matrix.
///
/// Supports per-length or global detection and respects an optional support mask (e.g., to
/// exclude extreme GC bins already handled elsewhere). Masked cells do not contribute to bound
/// estimation but are still clamped/masked using the bounds from supported cells. Non-finite
/// values are ignored during bound calculation, and MaskNaN writes NaN back into the matrix for
/// downstream handling. Stats separate supported vs unsupported cells; a later hard clamp is
/// tracked in `hard_clamped`.
pub(crate) fn apply_outliers_to_matrix(
    matrix: &mut Array2<f64>,
    support_mask: Option<&Array2<bool>>,
    scope: OutlierScope,
    rule: OutlierRule,
    action: OutlierAction,
) -> OutlierStats {
    if matches!(rule, OutlierRule::None) {
        return OutlierStats::default();
    }

    let apply_row = |row: &mut [f64], mask: Option<&[bool]>| -> OutlierStats {
        let mut vals: Vec<f32> = Vec::with_capacity(row.len());
        for (idx, v) in row.iter().enumerate() {
            if !mask.is_none_or(|m| m[idx]) {
                continue;
            }
            let f = *v as f32;
            if f.is_finite() {
                vals.push(f);
            }
        }

        if vals.len() < 2 {
            return OutlierStats::default();
        }

        let Some((lo, hi)) = outlier_bounds(&vals, rule) else {
            return OutlierStats::default();
        };

        let mut stats = OutlierStats::default();

        for (idx, v) in row.iter_mut().enumerate() {
            if !v.is_finite() {
                continue;
            }
            let in_support = mask.is_none_or(|m| m[idx]);
            if in_support {
                stats.total_examined += 1;
            } else {
                stats.unsupported_examined += 1;
            }

            match action {
                OutlierAction::Winsorize => {
                    if *v < lo as f64 {
                        *v = lo as f64;
                        if in_support {
                            stats.total_outliers_handled += 1;
                        } else {
                            stats.unsupported_outliers_handled += 1;
                        }
                    } else if *v > hi as f64 {
                        *v = hi as f64;
                        if in_support {
                            stats.total_outliers_handled += 1;
                        } else {
                            stats.unsupported_outliers_handled += 1;
                        }
                    }
                }
                OutlierAction::MaskNaN => {
                    if *v < lo as f64 || *v > hi as f64 {
                        *v = f64::NAN;
                        if in_support {
                            stats.total_outliers_handled += 1;
                        } else {
                            stats.unsupported_outliers_handled += 1;
                        }
                    }
                }
            }
        }
        stats
    };

    let mut total_stats = OutlierStats::default();

    match scope {
        OutlierScope::PerLength => {
            for (row_idx, mut row) in matrix.axis_iter_mut(Axis(0)).enumerate() {
                let mask_row = support_mask.map(|m| m.row(row_idx));
                let mask_values =
                    mask_row.map(|row_mask| row_mask.iter().copied().collect::<Vec<_>>());
                let mask_slice = mask_values.as_deref();

                // Fall back to a temporary buffer if ndarray gives us a non-contiguous row view
                let row_stats = if let Some(row_slice) = row.as_slice_mut() {
                    apply_row(row_slice, mask_slice)
                } else {
                    let mut row_values = row.iter().copied().collect::<Vec<_>>();
                    let row_stats = apply_row(&mut row_values, mask_slice);
                    for (slot, value) in row.iter_mut().zip(row_values.into_iter()) {
                        *slot = value;
                    }
                    row_stats
                };
                total_stats.add(row_stats);
            }
        }
        OutlierScope::Global => {
            let mask = support_mask.map(|m| m.iter().copied().collect::<Vec<_>>());
            let mask_slice = mask.as_deref();
            let stats = if let Some(matrix_slice) = matrix.as_slice_mut() {
                apply_row(matrix_slice, mask_slice)
            } else {
                let mut values = matrix.iter().copied().collect::<Vec<_>>();
                let stats = apply_row(&mut values, mask_slice);
                for (slot, value) in matrix.iter_mut().zip(values.into_iter()) {
                    *slot = value;
                }
                stats
            };
            total_stats.add(stats);
        }
    }
    total_stats
}
