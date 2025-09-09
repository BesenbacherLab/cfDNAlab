use ndarray::{Array2, s};

use crate::utils::gc::counting::GCCounts;

/// How to normalize weights within each fragment-length row.
///
/// - `RowMean`: divide the row by its **arithmetic** mean so the average weight is 1.
/// - `CountWeighted`: scale the row so the **count-weighted** sum is preserved:
///   `sum(w[c]*count[c])` equals the original row sum after scaling.
#[derive(Debug, Clone, Copy)]
pub enum WeightNormalizer {
    /// Divide by arithmetic mean → row’s average weight = 1 (totals may change).
    RowMean,
    /// Scale so `sum(w[c] * count[c]) == row_sum` (totals preserved; counts-weighted).
    CountWeighted,
    /// Scale all weights the same across lengths to get the original total sum.
    CountWeightedGlobal,
    /// Leave as-is.
    None,
}

impl WeightNormalizer {
    pub fn apply(&self, w: &mut Array2<f32>, gccounts: &GCCounts) {
        match *self {
            WeightNormalizer::RowMean => {
                let (n_len, n_gc) = (w.nrows(), w.ncols());
                for r in 0..n_len {
                    let mean_w = (0..n_gc).map(|c| w[(r, c)]).sum::<f32>() / n_gc as f32;
                    if mean_w.is_finite() && mean_w > 0.0 {
                        for c in 0..n_gc {
                            w[(r, c)] /= mean_w;
                        }
                    }
                }
            }
            WeightNormalizer::CountWeighted => {
                let (n_len, n_gc) = (w.nrows(), w.ncols());
                for r in 0..n_len {
                    let row_counts = &gccounts.counts[r];
                    let row_sum: f32 = row_counts.iter().map(|&x| x as f32).sum();
                    if row_sum > 0.0 {
                        let corrected_sum: f32 =
                            (0..n_gc).map(|c| w[(r, c)] * row_counts[c] as f32).sum();
                        if corrected_sum.is_finite() && corrected_sum > 0.0 {
                            let s = row_sum / corrected_sum;
                            for c in 0..n_gc {
                                w[(r, c)] *= s;
                            }
                        }
                    } else {
                        w.slice_mut(s![r, ..]).fill(1.0);
                    }
                }
            }
            WeightNormalizer::CountWeightedGlobal => {
                let (n_len, n_gc) = (w.nrows(), w.ncols());
                let total: f32 = (0..n_len)
                    .map(|r| gccounts.counts[r].iter().map(|&x| x as f32).sum::<f32>())
                    .sum();

                let corrected: f32 = (0..n_len)
                    .map(|r| {
                        (0..n_gc)
                            .map(|c| w[(r, c)] * gccounts.counts[r][c] as f32)
                            .sum::<f32>()
                    })
                    .sum();

                if corrected.is_finite() && corrected > 0.0 {
                    let s = total / corrected;
                    for r in 0..n_len {
                        for c in 0..n_gc {
                            w[(r, c)] *= s;
                        }
                    }
                }
            }
            WeightNormalizer::None => {}
        }
    }
}
