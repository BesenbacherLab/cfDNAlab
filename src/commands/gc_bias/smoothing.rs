// src/gc/smoothing.rs

use ndarray::{Array2, ArrayViewMut2, s};

// TODO: Check what GCParagon does

/// Smoothing along the GC axis (per fragment-length row).
///
/// Apply this **before** normalization. Typically you outlier-handle first
/// (e.g., winsorize/mask), then smooth, then normalize.
///
/// All methods are **NaN-aware**: NaNs are ignored when computing local means.
/// If a window has no finite values, the output is set to `NaN` at that bin.
#[derive(Debug, Clone, Copy)]
pub enum GCSmoother1D {
    /// No smoothing.
    None,
    /// Moving-average (boxcar) with half-window `radius` (window = 2*radius+1).
    Boxcar { radius: usize },
    // Future:
    // Gaussian { sigma: f32, radius: usize },
    // Lowess { frac: f32, iters: usize },
}

impl GCSmoother1D {
    /// In-place smoothing of a `(n_lengths, n_gc)` weight matrix along the GC axis.
    ///
    /// - Operates row-wise (each fragment-length independently).
    /// - NaN-aware moving average: ignores NaNs in the window when averaging.
    /// - `radius = 0` is a no-op.
    pub fn apply(&self, w: &mut Array2<f32>) {
        match *self {
            GCSmoother1D::None => {}
            GCSmoother1D::Boxcar { radius } => {
                if radius == 0 || w.ncols() == 0 || w.nrows() == 0 {
                    return;
                }
                boxcar_rows_nan_aware(w.view_mut(), radius);
            }
        }
    }
}

// TODO: We're not smoothing counts but corrections
/// NaN-aware boxcar smoothing along columns, applied to each row independently.
/// Uses prefix sums of (values with NaNs->0) and counts of finite values to get O(1) window means.
fn boxcar_rows_nan_aware(mut w: ArrayViewMut2<f32>, radius: usize) {
    let (n_len, n_gc) = (w.nrows(), w.ncols());
    let win = 2 * radius + 1;
    if win <= 1 || n_gc == 0 || n_len == 0 {
        return;
    }

    // Temporary buffer for the smoothed row.
    let mut out = vec![0.0f32; n_gc];
    let mut pref_sum = vec![0.0f32; n_gc + 1];
    let mut pref_cnt = vec![0usize; n_gc + 1];

    for r in 0..n_len {
        // Build prefix sums for this row
        pref_sum.fill(0.0);
        pref_cnt.fill(0);
        for c in 0..n_gc {
            let x = w[(r, c)];
            let (val, cnt) = if x.is_finite() {
                (x, 1usize)
            } else {
                (0.0, 0usize)
            };
            pref_sum[c + 1] = pref_sum[c] + val;
            pref_cnt[c + 1] = pref_cnt[c] + cnt;
        }

        // Compute smoothed values
        for c in 0..n_gc {
            let lo = c.saturating_sub(radius);
            let hi = (c + radius + 1).min(n_gc);
            let sum = pref_sum[hi] - pref_sum[lo];
            let cnt = pref_cnt[hi] - pref_cnt[lo];
            out[c] = if cnt > 0 {
                sum / (cnt as f32)
            } else {
                f32::NAN
            };
        }

        // Write back
        w.slice_mut(s![r, ..])
            .assign(&ndarray::ArrayView1::from(&out));
    }
}
