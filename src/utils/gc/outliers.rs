use std::cmp::Ordering;

/// How to detect outliers in a 1D vector (e.g., one GC row).
#[derive(Debug, Clone, Copy)]
pub enum OutlierRule {
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

/// What to do with detected outliers.
#[derive(Debug, Clone, Copy)]
pub enum OutlierAction {
    /// Cap to the computed [lower, upper] bounds (winsorization).
    Winsorize,
    /// Replace with NaN (so a NaN-aware smoother can ignore them).
    MaskNaN,
}

/// Compute outlier bounds for `vals` per the rule. Ignores NaNs.
/// Returns `None` if rule is `None` or there’s not enough data.
pub fn outlier_bounds(vals: &[f32], rule: OutlierRule) -> Option<(f32, f32)> {
    let mut v: Vec<f32> = vals.iter().copied().filter(|x| x.is_finite()).collect();
    if v.len() < 2 || matches!(rule, OutlierRule::None) {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

    // tiny quantile helper (linear interpolation)
    let q = |p: f32| -> f32 {
        if v.is_empty() {
            return f32::NAN;
        }
        let p = p.clamp(0.0, 1.0);
        let x = p * ((v.len() - 1) as f32);
        let lo = x.floor() as usize;
        let hi = x.ceil() as usize;
        if lo == hi {
            v[lo]
        } else {
            let t = x - lo as f32;
            v[lo] * (1.0 - t) + v[hi] * t
        }
    };

    match rule {
        OutlierRule::None => None,
        OutlierRule::Quantile { lower, upper } => {
            if !(0.0..=1.0).contains(&lower) || !(0.0..=1.0).contains(&upper) || lower >= upper {
                return None;
            }
            Some((q(lower), q(upper)))
        }
        OutlierRule::TukeyIqr { k } => {
            let q1 = q(0.25);
            let q3 = q(0.75);
            let iqr = (q3 - q1).max(0.0);
            Some((q1 - k * iqr, q3 + k * iqr))
        }
        OutlierRule::StdDev { k } => {
            let n = v.len() as f32;
            let mean = v.iter().sum::<f32>() / n;
            let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f32>() / n.max(1.0);
            let sd = var.sqrt();
            Some((mean - k * sd, mean + k * sd))
        }
        OutlierRule::Mad { k } => {
            let med = q(0.5);
            let mut dev: Vec<f32> = v.into_iter().map(|x| (x - med).abs()).collect();
            dev.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
            // MAD scaled to be comparable to SD for normal data
            let mad = {
                let m = |p: f32| {
                    let n = dev.len();
                    if n == 0 {
                        return f32::NAN;
                    }
                    let x = p.clamp(0.0, 1.0) * ((n - 1) as f32);
                    let lo = x.floor() as usize;
                    let hi = x.ceil() as usize;
                    if lo == hi {
                        dev[lo]
                    } else {
                        let t = x - lo as f32;
                        dev[lo] * (1.0 - t) + dev[hi] * t
                    }
                };
                1.4826_f32 * m(0.5)
            };
            Some((med - k * mad, med + k * mad))
        }
    }
}

/// Apply `action` to values outside bounds (inclusive). Returns number of changes.
/// Call this **before smoothing**.
pub fn apply_outliers_inplace(row: &mut [f32], rule: OutlierRule, action: OutlierAction) -> usize {
    let Some((lo, hi)) = outlier_bounds(row, rule) else {
        return 0;
    };
    let mut changed = 0usize;
    match action {
        OutlierAction::Winsorize => {
            for x in row.iter_mut() {
                if !x.is_finite() {
                    continue;
                }
                if *x < lo {
                    *x = lo;
                    changed += 1;
                } else if *x > hi {
                    *x = hi;
                    changed += 1;
                }
            }
        }
        OutlierAction::MaskNaN => {
            for x in row.iter_mut() {
                if !x.is_finite() {
                    continue;
                }
                if *x < lo || *x > hi {
                    *x = f32::NAN;
                    changed += 1;
                }
            }
        }
    }
    changed
}
