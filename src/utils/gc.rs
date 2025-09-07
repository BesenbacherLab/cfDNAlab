use anyhow::Result;
use ndarray::{Array1, Array2, ArrayViewMut2, arr0, s};
use ndarray_npy::NpzWriter;
use std::fs::File;
use std::path::Path;

/// Prefix sums (cumsum) to compute GC and AT fractions while excluding Ns.
/// `gc[i]`   = # of G/C in seq[0..i)
/// `acgt[i]`= # of A/T/G/C in seq[0..i)
pub struct GCPrefixes {
    pub gc: Vec<u32>,
    pub acgt: Vec<u32>,
}

/// Build prefix-sums (cumsum) for GC and ACGT (non-N) counts over a byte slice.
/// This lets you compute GC% on A/T/G/C only, so AT% = 1 - GC%.
///
/// Ignores Ns.
pub fn build_gc_prefixes(seq: &[u8]) -> GCPrefixes {
    let mut gc = Vec::with_capacity(seq.len() + 1);
    let mut acgt = Vec::with_capacity(seq.len() + 1);
    gc.push(0);
    acgt.push(0);

    for &b in seq {
        let is_gc = matches!(b, b'G' | b'g' | b'C' | b'c') as u32;
        let is_acgt = matches!(b, b'A' | b'a' | b'T' | b't' | b'G' | b'g' | b'C' | b'c') as u32;

        gc.push(gc.last().copied().unwrap() + is_gc);
        acgt.push(acgt.last().copied().unwrap() + is_acgt);
    }

    GCPrefixes { gc, acgt }
}

/// Compute the GC fraction for a window [start, end), excluding 'N's.
///
/// `min_acgt_count`: Minimum number of actual ACGT bases counted in the window.
///   E.g. if most of the window is blacklisted or Ns.
///
/// Returns None if the window has no A/T/G/C bases.
#[inline]
pub fn get_gc_fraction_in_window(
    prefixes: &GCPrefixes,
    start: usize,
    end: usize,
    min_acgt_fraction: f32,
    min_acgt_count: u32,
) -> Option<f32> {
    debug_assert!(
        start < end && end <= prefixes.gc.len() - 1,
        "GC window [{}, {}) out of bounds (len={})",
        start,
        end,
        prefixes.gc.len() - 1
    );
    let gc = prefixes.gc[end] - prefixes.gc[start];
    let acgt = prefixes.acgt[end] - prefixes.acgt[start];
    let length = end - start;

    if acgt == 0 || acgt < min_acgt_count || (acgt as f32 / length as f32) < min_acgt_fraction {
        return None;
    }
    let gc_frac = gc as f32 / acgt as f32;
    Some(gc_frac)
}

/// Count matrix for fragment coverage across GC fraction bins and fragment lengths.
///
/// The matrix is two-dimensional:
/// - Rows correspond to fragment lengths.
/// - Columns correspond to GC fraction bins (0–100).
#[derive(Debug, Clone)]
pub struct GCCounts {
    counts: Vec<Vec<u32>>,
    gc_min: usize,
    gc_max: usize,
    length_min: usize,
    length_max: usize,
}

impl GCCounts {
    /// Create a new `GCCounts` with specified ranges and binning.
    ///
    /// Parameters
    /// ----------
    /// gc_min: usize
    ///     Minimum GC bin (inclusive).
    /// gc_max: usize
    ///     Maximum GC bin (inclusive).
    /// length_min: usize
    ///     Minimum fragment length (inclusive).
    /// length_max: usize
    ///     Maximum fragment length (inclusive).
    ///
    /// Returns
    /// -------
    /// counts: GCCounts
    ///     A `GCCounts` object with all counts initialized to zero.
    pub fn new(gc_min: usize, gc_max: usize, length_min: usize, length_max: usize) -> Self {
        let num_gc_bins = gc_max - gc_min + 1;
        let num_lengths = length_max - length_min + 1;
        let counts = vec![vec![0u32; num_gc_bins]; num_lengths];
        Self {
            counts,
            gc_min,
            gc_max,
            length_min,
            length_max,
        }
    }

    /// Check whether `(length, gc)` is within configured ranges.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length to test.
    /// gc: usize
    ///     GC bin to test.
    ///
    /// Returns
    /// -------
    /// ok: bool
    ///     True if both indices are in range.
    #[inline]
    fn in_bounds(&self, length: usize, gc: usize) -> bool {
        (self.length_min..=self.length_max).contains(&length)
            && (self.gc_min..=self.gc_max).contains(&gc)
    }

    /// Compute row/column indices for `(length, gc)` if in bounds.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute, not zero-based).
    /// gc: usize
    ///     GC bin (absolute, not zero-based).
    ///
    /// Returns
    /// -------
    /// idx: Option<(usize, usize)>
    ///     `(row, col)` zero-based indices if in range, otherwise `None`.
    #[inline]
    pub fn index_of(&self, length: usize, gc: usize) -> Option<(usize, usize)> {
        if self.in_bounds(length, gc) {
            Some((length - self.length_min, gc - self.gc_min))
        } else {
            None
        }
    }

    /// Increment the counter for a given fragment length and GC bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
    pub fn incr(&mut self, length: usize, gc: usize) {
        if let Some((r, c)) = self.index_of(length, gc) {
            self.counts[r][c] = self.counts[r][c].saturating_add(1);
        }
    }

    // Get the count at a given fragment length and GC bin.
    ///
    /// Parameters
    /// ----------
    /// length: usize
    ///     Fragment length (absolute).
    /// gc: usize
    ///     GC bin (absolute).
    ///
    /// Returns
    /// -------
    /// count: Option<u32>
    ///     The count if indices are in range, otherwise `None`.
    pub fn get(&self, length: usize, gc: usize) -> Option<u32> {
        self.index_of(length, gc).map(|(r, c)| self.counts[r][c])
    }

    /// Number of length rows.
    ///
    /// Returns
    /// -------
    /// n: usize
    ///     Count of length bins.
    #[inline]
    pub fn n_lengths(&self) -> usize {
        self.length_max - self.length_min + 1
    }

    /// Number of GC columns.
    ///
    /// Returns
    /// -------
    /// n: usize
    ///     Count of GC bins.
    #[inline]
    pub fn n_gc_bins(&self) -> usize {
        self.gc_max - self.gc_min + 1
    }

    /// Create a zero-initialized `GCCounts` with the same ranges and shape as `self`.
    ///
    /// Returns
    /// -------
    /// counts: GCCounts
    ///     New object with all counts set to zero and identical configuration.
    #[inline]
    pub fn zeroed_like(&self) -> Self {
        let n_len = self.n_lengths();
        let n_gc = self.n_gc_bins();
        Self {
            counts: vec![vec![0u32; n_gc]; n_len],
            gc_min: self.gc_min,
            gc_max: self.gc_max,
            length_min: self.length_min,
            length_max: self.length_max,
        }
    }

    /// Check if two `GCCounts` are compatible for merging (same ranges and shape).
    ///
    /// Parameters
    /// ----------
    /// other: &GCCounts
    ///     The other counts matrix to compare with.
    ///
    /// Returns
    /// -------
    /// ok: bool
    ///     `true` if the two objects have identical ranges and matrix shape.
    #[inline]
    pub fn is_compatible_with(&self, other: &Self) -> bool {
        self.gc_min == other.gc_min
            && self.gc_max == other.gc_max
            && self.length_min == other.length_min
            && self.length_max == other.length_max
            && self.n_lengths() == other.n_lengths()
            && self.n_gc_bins() == other.n_gc_bins()
    }

    /// Merge (sum) counts from `other` into `self` using saturating addition.
    ///
    /// Parameters
    /// ----------
    /// other: &GCCounts
    ///     Counts to add into `self`. Must be compatible (same ranges/shape).
    ///
    /// Returns
    /// -------
    /// result: Result<(), anyhow::Error>
    ///     An error is returned if the two objects are incompatible.
    pub fn merge_from(&mut self, other: &Self) -> anyhow::Result<()> {
        if !self.is_compatible_with(other) {
            return Err(anyhow::anyhow!(
                "incompatible GCCounts: self={} vs other={}",
                self,
                other
            ));
        }
        for (r, row_other) in other.counts.iter().enumerate() {
            let row_self = &mut self.counts[r];
            for (c, &v) in row_other.iter().enumerate() {
                row_self[c] = row_self[c].saturating_add(v);
            }
        }
        Ok(())
    }

    /// Collapse (sum) an iterator of `GCCounts` into a single object.
    ///
    /// All inputs must be compatible (same ranges/shape). Uses saturating addition.
    ///
    /// Parameters
    /// ----------
    /// iter: IntoIterator<Item = &GCCounts>
    ///     Collection of references to `GCCounts` to be summed.
    ///
    /// Returns
    /// -------
    /// total: GCCounts
    ///     The element-wise sum across all inputs.
    ///
    /// Examples
    /// --------
    /// ```ignore
    /// // Sum across chromosomes:
    /// // let by_chr: HashMap<String, GCCounts> = ...;
    /// let total = GCCounts::collapse(by_chr.values())?;
    /// ```
    pub fn collapse<'a, I>(iter: I) -> anyhow::Result<Self>
    where
        I: IntoIterator<Item = &'a GCCounts>,
    {
        let mut it = iter.into_iter();
        let first = it
            .next()
            .ok_or_else(|| anyhow::anyhow!("collapse requires at least one GCCounts"))?;
        let mut acc = first.clone(); // Start from the first
        for g in it {
            acc.merge_from(g)?; // Then merge the rest
        }
        Ok(acc)
    }
}

impl Default for GCCounts {
    /// Create an empty default `GCCounts` (0–100 GC, 20–600 length).
    fn default() -> Self {
        Self::new(0, 100, 20, 600)
    }
}

impl std::fmt::Display for GCCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GCCounts(gc:[{}..={}], len:[{}..={}], dims:({},{}) )",
            self.gc_min,
            self.gc_max,
            self.length_min,
            self.length_max,
            self.n_lengths(),
            self.n_gc_bins()
        )
    }
}

// TODO: Consider reusing WindowSpec?
#[repr(u8)]
#[derive(Debug, Clone, Copy)]
enum CorrectionMode {
    Global = 0,
    BySize = 1,
    ByBed = 2,
}

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

/// What expected/target distribution to correct **to**.
#[derive(Debug, Clone)]
pub enum ExpectedSpec {
    /// Uniform per length: each GC bin gets equal expected share of that row’s total.
    UniformPerLength,
    /// Reference across GC only (same for all lengths). Must be length n_gc and non-negative.
    Reference1D(Array1<f32>),
    /// Reference per (length, GC). Must be shape (n_len, n_gc) and non-negative.
    Reference2D(Array2<f32>),
}

/// Optional 1D smoothing of weights along the GC axis (per length row).
#[derive(Debug, Clone, Copy)]
pub enum Smooth1D {
    /// No smoothing.
    None,
    /// Boxcar / moving-average with half-window `radius` (total window = 2*radius+1).
    Boxcar { radius: usize },
    // Future:
    // Gaussian { sigma: f32, radius: usize },
    // Lowess { frac: f32, iters: usize },
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
    pub smooth_gc: Smooth1D,
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
            normalizer: WeightNormalizer::CountWeighted, // Consider which to make default!
            smooth_gc: Smooth1D::None,
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
        match self.smooth_gc {
            Smooth1D::None => {}
            Smooth1D::Boxcar { radius } if radius > 0 => {
                boxcar_smooth_rows(&mut w.view_mut(), radius);
            }
            _ => {}
        }

        // Normalize rows
        match self.normalizer {
            WeightNormalizer::RowMean => {
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
                        // neutral row → expected equals counts (weights→1 later)
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

// TODO: Check this, haven't looked at it yet!!!
/// In-place boxcar smoothing along GC axis for each row.
fn boxcar_smooth_rows(w: &mut ArrayViewMut2<f32>, radius: usize) {
    let (n_len, n_gc) = (w.nrows(), w.ncols());
    let win = 2 * radius + 1;
    if win <= 1 {
        return;
    }
    let mut tmp = Array2::<f32>::zeros((n_len, n_gc));
    for r in 0..n_len {
        // prefix sum for O(1) window sums
        let mut ps = vec![0.0_f32; n_gc + 1];
        for c in 0..n_gc {
            ps[c + 1] = ps[c] + w[(r, c)];
        }
        for c in 0..n_gc {
            let lo = c.saturating_sub(radius);
            let hi = (c + radius + 1).min(n_gc);
            let sum = ps[hi] - ps[lo];
            tmp[(r, c)] = sum / (hi - lo) as f32;
        }
    }
    w.assign(&tmp);
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
// /// Given a target distribution q_g (will be renormalized to sum=1; if None → uniform),
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

pub fn save_correction_npz(path: impl AsRef<Path>, corr: &Correction) -> Result<()> {
    let mut npz = NpzWriter::new(File::create(path)?);

    npz.add_array("weights", &corr.weights)?;
    npz.add_array("gc_min", &arr0(corr.gc_min as i64))?;
    npz.add_array("gc_max", &arr0(corr.gc_max as i64))?;
    npz.add_array("len_min", &arr0(corr.len_min as i64))?;
    npz.add_array("len_max", &arr0(corr.len_max as i64))?;
    npz.add_array("mode", &arr0(corr.mode as u8))?;

    match corr.mode {
        CorrectionMode::BySize => {
            let bs = corr.bin_size.expect("bin_size must be set for BySize");
            npz.add_array("bin_size", &arr0(bs as i64))?;
        }
        CorrectionMode::ByBed => {
            let wins = corr
                .windows
                .as_ref()
                .expect("windows must be set for ByBed");
            // Flatten Vec<(u64,u64)> → Array2<u64> (N,2)
            let mut flat = Vec::with_capacity(wins.len() * 2);
            for (s, e) in wins {
                flat.push(*s);
                flat.push(*e);
            }
            let arr = Array2::from_shape_vec((wins.len(), 2), flat)?;
            npz.add_array("windows", &arr)?;
        }
        CorrectionMode::Global => {}
    }

    npz.finish()?;
    Ok(())
}

// Multiple objects in one .npz
pub fn save_corrections_npz(path: impl AsRef<Path>, list: &[Correction]) -> Result<()> {
    let mut npz = NpzWriter::new(File::create(path)?);
    for (i, corr) in list.iter().enumerate() {
        let p = |k: &str| format!("obj{}_{}", i, k);
        npz.add_array(&p("weights"), &corr.weights)?;
        npz.add_array(&p("gc_min"), &arr0(corr.gc_min as i64))?;
        npz.add_array(&p("gc_max"), &arr0(corr.gc_max as i64))?;
        npz.add_array(&p("len_min"), &arr0(corr.len_min as i64))?;
        npz.add_array(&p("len_max"), &arr0(corr.len_max as i64))?;
        npz.add_array(&p("mode"), &arr0(corr.mode as u8))?;
        match corr.mode {
            CorrectionMode::BySize => {
                let bs = corr.bin_size.expect("bin_size must be set for BySize");
                npz.add_array(&p("bin_size"), &arr0(bs as i64))?;
            }
            CorrectionMode::ByBed => {
                let wins = corr
                    .windows
                    .as_ref()
                    .expect("windows must be set for ByBed");
                let mut flat = Vec::with_capacity(wins.len() * 2);
                for (s, e) in wins {
                    flat.push(*s);
                    flat.push(*e);
                }
                let arr = Array2::from_shape_vec((wins.len(), 2), flat)?;
                npz.add_array(&p("windows"), &arr)?;
            }
            CorrectionMode::Global => {}
        }
    }
    npz.finish()?;
    Ok(())
}
