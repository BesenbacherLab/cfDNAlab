use crate::commands::gc_bias::counting::GCCounts;
use crate::commands::gc_bias::normalization::WeightNormalizer;
use crate::commands::gc_bias::reference::ExpectedSpec;
use crate::commands::gc_bias::smoothing::GCSmoother1D;
use ndarray::Array2;

// TODO: Consider reusing WindowSpec?
#[repr(u8)]
#[derive(Debug, Clone, Copy)]
pub enum CorrectionMode {
    Global = 0,
    BySize = 1,
    ByBed = 2,
}

/// End-to-end builder for GC-bias correction weights.
#[derive(Debug, Clone)]
pub struct CorrectionBuilder {
    pub expected: ExpectedSpec,
    pub pseudocount: f32,
    pub clamp_min: f32,
    pub clamp_max: f32,
    pub normalizer: WeightNormalizer,
    // Smoothing along GC axis
    pub smooth_gc: GCSmoother1D,
    // Future hooks:
    // pub interpolate_missing: Interp1D,
    // pub across_lengths: Option<SmoothAcrossLengths>,
}

impl Default for CorrectionBuilder {
    fn default() -> Self {
        Self {
            expected: ExpectedSpec::UniformPerLength,
            pseudocount: 1.0,
            clamp_min: 0.1,
            clamp_max: 10.0,
            normalizer: WeightNormalizer::None, // Consider which to make default!
            smooth_gc: GCSmoother1D::None,
        }
    }
}

impl CorrectionBuilder {
    /// Build a `Correction` (weights) from counts with the configured strategy.
    ///
    /// Pipeline (per length row):
    ///  - Compute expected counts E[r,c] from `expected` (uniform or reference).
    ///  - Raw weight: w = E / (count + pseudocount).
    ///  - Clamp to [clamp_min, clamp_max].
    ///  - Optional smoothing along GC axis.
    ///  - Normalize per `normalizer` (RowMean / CountWeighted / None).
    ///  - Final clamp (guard against blow-ups after normalization).
    pub fn build(&self, gccounts: &GCCounts) -> Correction {
        let n_len = gccounts.n_lengths();
        let n_gc = gccounts.n_gc_bins();

        // Expected counts matrix (same shape as counts)
        let exp = self.expected_counts(gccounts);

        // Raw weights
        let mut w = Array2::<f32>::zeros((n_len, n_gc));
        for r in 0..n_len {
            let row_counts = &gccounts.counts[r];
            for c in 0..n_gc {
                let denom = row_counts[c] as f32 + self.pseudocount;
                let ww = exp[(r, c)] / denom;
                w[(r, c)] = ww.clamp(self.clamp_min, self.clamp_max);
            }
        }

        // Smooth along GC axis (row-wise)
        self.smooth_gc.apply(&mut w);

        // Normalize rows
        self.normalizer.apply(&mut w, gccounts);

        // Clamp weights
        for r in 0..n_len {
            for c in 0..n_gc {
                w[(r, c)] = w[(r, c)].clamp(self.clamp_min, self.clamp_max);
            }
        }

        Correction {
            weights: w,
            gc_min: gccounts.gc_min,
            gc_max: gccounts.gc_max,
            len_min: gccounts.length_min,
            len_max: gccounts.length_max,
            mode: CorrectionMode::Global,
            bin_size: None,
            windows: None,
        }
    }

    /// Build expected counts `E[r,c]` from the chosen `expected` spec.
    fn expected_counts(&self, gccounts: &GCCounts) -> Array2<f32> {
        let n_len = gccounts.n_lengths();
        let n_gc = gccounts.n_gc_bins();
        let mut e = Array2::<f32>::zeros((n_len, n_gc));

        match &self.expected {
            ExpectedSpec::UniformPerLength => {
                for r in 0..n_len {
                    let row_sum: f32 = gccounts.counts[r].iter().map(|&x| x as f32).sum();
                    if row_sum == 0.0 {
                        // neutral row -> expected equals counts (weights->1 later)
                        for c in 0..n_gc {
                            e[(r, c)] = 1.0;
                        }
                    } else {
                        let tgt = row_sum / n_gc as f32;
                        for c in 0..n_gc {
                            e[(r, c)] = tgt;
                        }
                    }
                }
            }
            ExpectedSpec::Reference1D(p_gc) => {
                assert_eq!(p_gc.len(), n_gc, "Reference1D length != n_gc");
                // normalize ref distribution to sum=1 (defensive)
                let sum: f32 = p_gc.iter().sum();
                let norm = if sum.is_finite() && sum > 0.0 {
                    sum
                } else {
                    1.0
                };
                for r in 0..n_len {
                    let row_sum: f32 = gccounts.counts[r].iter().map(|&x| x as f32).sum();
                    for c in 0..n_gc {
                        e[(r, c)] = row_sum * (p_gc[c] / norm);
                    }
                }
            }
            ExpectedSpec::Reference2D(p_rg) => {
                assert_eq!(p_rg.nrows(), n_len, "Reference2D nrows != n_len");
                assert_eq!(p_rg.ncols(), n_gc, "Reference2D ncols != n_gc");
                for r in 0..n_len {
                    let row_sum: f32 = gccounts.counts[r].iter().map(|&x| x as f32).sum();
                    let row = p_rg.row(r);
                    let sum: f32 = row.iter().copied().sum();
                    let norm = if sum.is_finite() && sum > 0.0 {
                        sum
                    } else {
                        1.0
                    };
                    for c in 0..n_gc {
                        e[(r, c)] = row_sum * (row[c] / norm);
                    }
                }
            }
        }
        e
    }
}

#[derive(Debug, Clone)]
pub struct Correction {
    pub weights: Array2<f32>, // shape: (n_lengths, n_gc)
    pub gc_min: usize,
    pub gc_max: usize,
    pub len_min: usize,
    pub len_max: usize,
    pub mode: CorrectionMode,
    pub bin_size: Option<usize>,          // when mode == BySize
    pub windows: Option<Vec<(u64, u64)>>, // when mode == ByBed
}

// impl Correction {
//     /// Build flattening weights from a single `GCCounts`.
//     ///
//     /// Method:
//     /// - For each length row `r`, compute target per-bin:
//     ///       `target = (row_sum + alpha * n_gc) / n_gc`
//     /// - Raw weight: `w = target / (count + alpha)` (Laplace smoothing).
//     /// - Clamp to `[clamp_min, clamp_max]`.
//     /// - Normalize per `norm`:
//     ///   - `RowMean`: divide by arithmetic mean so row’s **average weight = 1**.
//     ///   - `CountWeighted`: scale so **sum(w * count) == row_sum** (totals preserved).
//     ///   - `None`: leave raw clamped weights.
//     pub fn from_counts(
//         gccounts: &GCCounts,
//         alpha: f32,
//         clamp_min: f32,
//         clamp_max: f32,
//         norm: WeightNormalizer,
//     ) -> Self {
//         let n_len = gccounts.n_lengths();
//         let n_gc = gccounts.n_gc_bins();
//         let mut w = Array2::<f32>::zeros((n_len, n_gc));

//         for r in 0..n_len {
//             let row = &gccounts.counts[r]; // same module ⇒ ok to read private field
//             let row_sum: f32 = row.iter().map(|&x| x as f32).sum();

//             if row_sum == 0.0 {
//                 // no data for this length: neutral weights
//                 w.slice_mut(s![r, ..]).fill(1.0);
//                 continue;
//             }

//             let target = (row_sum + alpha * (n_gc as f32)) / (n_gc as f32);

//             // raw, clamped weights
//             for c in 0..n_gc {
//                 let denom = row[c] as f32 + alpha;
//                 let mut ww = target / denom;
//                 if ww < clamp_min {
//                     ww = clamp_min;
//                 }
//                 if ww > clamp_max {
//                     ww = clamp_max;
//                 }
//                 w[(r, c)] = ww;
//             }

//             // Normalize the row
//             match norm {
//                 WeightNormalizer::RowMean => {
//                     let mean_w = (0..n_gc).map(|c| w[(r, c)]).sum::<f32>() / (n_gc as f32);
//                     if mean_w.is_finite() && mean_w > 0.0 {
//                         for c in 0..n_gc {
//                             w[(r, c)] /= mean_w;
//                         }
//                     }
//                 }
//                 WeightNormalizer::CountWeighted => {
//                     let corrected_sum: f32 = (0..n_gc).map(|c| w[(r, c)] * (row[c] as f32)).sum();
//                     if corrected_sum.is_finite() && corrected_sum > 0.0 {
//                         let scale = row_sum / corrected_sum;
//                         for c in 0..n_gc {
//                             w[(r, c)] *= scale;
//                         }
//                     }
//                 }
//                 WeightNormalizer::None => {}
//             }
//         }

//         Self {
//             weights: w,
//             gc_min: gccounts.gc_min,
//             gc_max: gccounts.gc_max,
//             len_min: gccounts.length_min,
//             len_max: gccounts.length_max,
//             mode: CorrectionMode::Global,
//             bin_size: None,
//             windows: None,
//         }
//     }

// /// Build a GC-only correction whose weights depend **only on GC bin** (shared across lengths).
// ///
// /// Let C_g be the total counts in GC bin g and p_g = (C_g + pseudocount) / (sum C + pseudocount * n_gc).
// /// Given a target distribution q_g (will be renormalized to sum=1; if None -> uniform),
// /// define w_g = (q_g / p_g), clamp to [clamp_min, clamp_max], and apply a single
// /// global scale so the grand total is preserved exactly:
// ///     s = (sum_{l,g} c_{l,g}) / (sum_{l,g} w_g * c_{l,g})
// /// Final weights are W_{l,g} = s * w_g (identical across l).
// pub fn gc_only_from_counts(
//     counts: &GCCounts,
//     target_q: Option<&[f32]>, // len = n_gc, will be renormalized; None => uniform
//     pseudocount: f32,
//     clamp_min: f32,
//     clamp_max: f32,
// ) -> Self {
//     let n_len = counts.n_lengths();
//     let n_gc = counts.n_gc_bins();

//     // 1) Column totals C_g and grand total
//     let mut col_tot = vec![0.0_f32; n_gc];
//     let mut grand: f32 = 0.0;
//     for r in 0..n_len {
//         for c in 0..n_gc {
//             let v = counts.counts[r][c] as f32;
//             col_tot[c] += v;
//             grand += v;
//         }
//     }

//     // 2) Observed p_g with pseudocount, and target q_g (normalized)
//     let denom = grand + pseudocount * (n_gc as f32);
//     let mut p = vec![0.0_f32; n_gc];
//     for c in 0..n_gc {
//         p[c] = (col_tot[c] + pseudocount) / denom;
//     }
//     let mut q = vec![0.0_f32; n_gc];
//     if let Some(tq) = target_q {
//         assert_eq!(tq.len(), n_gc, "target_q length != n_gc");
//         let s: f32 = tq.iter().copied().sum();
//         let s = if s.is_finite() && s > 0.0 { s } else { 1.0 };
//         for c in 0..n_gc {
//             q[c] = tq[c] / s;
//         }
//     } else {
//         let u = 1.0 / (n_gc as f32);
//         for c in 0..n_gc {
//             q[c] = u;
//         }
//     }

//     // 3) Per-GC weights w_g (clamped)
//     let mut w_gc = vec![1.0_f32; n_gc];
//     for c in 0..n_gc {
//         let mut w = if p[c] > 0.0 { q[c] / p[c] } else { 1.0 };
//         if w < clamp_min {
//             w = clamp_min;
//         }
//         if w > clamp_max {
//             w = clamp_max;
//         }
//         w_gc[c] = w;
//     }

//     // 4) Global scale to preserve the grand total exactly
//     let corrected_total: f32 = (0..n_gc).map(|c| w_gc[c] * col_tot[c]).sum();
//     let scale = if corrected_total.is_finite() && corrected_total > 0.0 {
//         grand / corrected_total
//     } else {
//         1.0
//     };

//     // 5) Expand to a full matrix: same weights for every row
//     let mut w = Array2::<f32>::zeros((n_len, n_gc));
//     for r in 0..n_len {
//         for c in 0..n_gc {
//             w[(r, c)] = scale * w_gc[c];
//         }
//     }

//     Self {
//         weights: w,
//         gc_min: counts.gc_min,
//         gc_max: counts.gc_max,
//         len_min: counts.length_min,
//         len_max: counts.length_max,
//         mode: CorrectionMode::Global, // metadata; set as you prefer
//         bin_size: None,
//         windows: None,
//     }
// }

// /// After applying this correction to `counts`, return the **length distribution**:
// /// a vector of per-length totals (summing across GC).
// pub fn length_distribution_after(&self, counts: &GCCounts) -> Vec<f64> {
//     assert!(
//         self.is_compatible_with(counts),
//         "Correction/GCCounts mismatch"
//     );
//     let (nr, nc) = (self.weights.nrows(), self.weights.ncols());
//     let mut out = vec![0.0_f64; nr];
//     for r in 0..nr {
//         let mut row_sum = 0.0_f64;
//         for c in 0..nc {
//             row_sum += self.weights[(r, c)] as f64 * counts.counts[r][c] as f64;
//         }
//         out[r] = row_sum;
//     }
//     out
// }

//     /// Convenience defaults for typical GC-bias ranges.
//     ///
//     /// Uses `alpha=1.0`, clamps to `[0.1, 10.0]`, and preserves per-length totals
//     /// (`WeightNormalizer::CountWeighted`). Change as needed.
//     #[inline]
//     pub fn from_counts_default(gccounts: &GCCounts) -> Self {
//         Self::from_counts(gccounts, 1.0, 0.1, 10.0, WeightNormalizer::CountWeighted)
//     }

//     /// Apply this correction to a `GCCounts`, returning a dense array of corrected values.
//     ///
//     /// Panics if the ranges/shapes don’t match (call `is_compatible_with` to check first).
//     pub fn apply_to_counts(&self, gccounts: &GCCounts) -> Array2<f32> {
//         assert!(
//             self.is_compatible_with(gccounts),
//             "Correction/GCCounts mismatch: correction(gc:[{}..={}], len:[{}..={}], shape=({},{}) ) \
//              vs counts(gc:[{}..={}], len:[{}..={}], shape=({},{}) )",
//             self.gc_min,
//             self.gc_max,
//             self.len_min,
//             self.len_max,
//             self.weights.nrows(),
//             self.weights.ncols(),
//             // TODO: Use the fmt of GCCounts here?? And make the same for the Correction obj?
//             gccounts.gc_min,
//             gccounts.gc_max,
//             gccounts.length_min,
//             gccounts.length_max,
//             gccounts.n_lengths(),
//             gccounts.n_gc_bins()
//         );

//         let mut out = Array2::<f32>::zeros((self.weights.nrows(), self.weights.ncols()));
//         for r in 0..self.weights.nrows() {
//             for c in 0..self.weights.ncols() {
//                 out[(r, c)] = self.weights[(r, c)] * gccounts.counts[r][c] as f32;
//             }
//         }
//         out
//     }

//     /// Check if this correction’s ranges/shape match a `GCCounts`.
//     #[inline]
//     pub fn is_compatible_with(&self, gccounts: &GCCounts) -> bool {
//         self.gc_min == gccounts.gc_min
//             && self.gc_max == gccounts.gc_max
//             && self.len_min == gccounts.length_min
//             && self.len_max == gccounts.length_max
//             && self.weights.nrows() == gccounts.n_lengths()
//             && self.weights.ncols() == gccounts.n_gc_bins()
//     }
// }
