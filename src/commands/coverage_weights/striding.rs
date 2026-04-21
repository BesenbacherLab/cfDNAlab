use anyhow::Result;
use fxhash::FxHashMap;

use crate::shared::interval::Interval;

/// Treat support smaller than this as zero so tiny numerical residues do not participate
/// in the global mean or blow up during inversion
const SUPPORT_FLOOR: f64 = 1e-10;

/// Represents a single stride bin with coverage and its overlapping mega-bin average coverage.
#[derive(Debug, Clone, Copy)]
pub struct StrideBin {
    /// Checked genomic span of the stride-bin
    pub interval: Interval<u32>,
    /// Average fragment coverage for the stride-bin
    pub average_coverage: f32,
    /// Average coverage of overlapping mega-bins
    pub average_overlap_coverage: f32,
    /// Scaling factor for normalizing coverage
    /// across the genome. Normalized across
    /// all stride-bins for a mean of 1.0.
    pub scaling_factor: f32,
}

impl StrideBin {
    /// Return the inclusive start coordinate of the stride-bin.
    #[inline]
    pub fn start(&self) -> u32 {
        self.interval.start()
    }

    /// Return the exclusive end coordinate of the stride-bin.
    #[inline]
    pub fn end(&self) -> u32 {
        self.interval.end()
    }

    /// Calculates the size (length) of the bin.
    pub fn size(&self) -> u32 {
        self.interval.len()
    }
}

#[inline]
fn is_effectively_zero(value: f64) -> bool {
    value.abs() < SUPPORT_FLOOR
}

/// Build the triangular integer weights, centered in the kernel.
///
/// Suppose `half_window = (bin_size / stride) - 1`, then `total_window = 2 * half_window + 1`.
///
/// E.g., `half_window = 2  ->  [1, 2, 3, 2, 1]`
fn triangular_weights(half_window: usize) -> Vec<usize> {
    // Total number of weights in the symmetric kernel
    let w_len = 2 * half_window + 1;
    let mut weights = Vec::with_capacity(w_len);

    // Build the triangle by walking from left to right.
    // For indices <= center, weights increase; then they decrease symmetrically.
    for k in 0..w_len {
        if k <= half_window {
            // Left side (including center): 1, 2, ..., half_window+1
            weights.push(k + 1);
        } else {
            // Right side mirrors the left: ..., 2, 1
            weights.push((2 * half_window - k) + 1);
        }
    }
    weights
}

/// Fill `bins[i].average_overlap_coverage` with a *triangularly weighted* average
/// of neighboring stride-bin averages around `i`.
///
/// Goal: For each stride-bin `i`, approximate the average coverage of all overlapping
/// "megabins" (large windows of size `bin_size`) *without* explicitly enumerating them.
/// Instead, use a fixed triangular kernel whose integer weights encode:
///   “How many megabins include this neighbor when centered at i?”
///
/// Key quantities:
/// - `half_window = (bin_size / stride) - 1`
///
///     The kernel radius, *measured in stride-bins*. The full kernel length is
///     `2*half_window + 1`. Example: if bin_size=5Mb and stride=0.5Mb,
///     then bin_size/stride = 10, half_window = 9, kernel length = 19.
///
/// - `weights = [1, 2, ..., half_window+1, ..., 2, 1]`.
///
///     The triangular profile. The center gets the largest weight.
///
/// What happens at chromosome edges?
/// - Near the left edge, there are fewer neighbors to the left of `i`. We *truncate*
///   the kernel on the left and start using it from some offset `w_start`.
/// - Near the right edge, we similarly truncate on the right.
/// - In both cases we still align the kernel's center to the current `i`.
///
/// Implementation outline per `i`:
/// 1) Pick the slice of bins we can actually use: `[start_i .. end_i)`
/// 2) Compute `w_start`, the index into `weights` that aligns the first usable bin with
///    the correct kernel position (so the kernel's center still targets `i`).
/// 3) Accumulate `sum( average_coverage[j] * weight[j] )` and `sum(weights)`
/// 4) Normalize: `average_overlap_coverage[i] = weighted_sum / sum_weights`.
///
/// Parameters
/// ----------
/// - bins:
///     Stride bins with `average_coverage` set to per-base averages (mask-adjusted if desired).
/// - bin_size:
///     Large window size; used only to derive the kernel radius.
/// - stride:
///     Stride size; used only to derive the kernel radius.
pub fn fill_triangular_overlap(bins: &mut Vec<StrideBin>, bin_size: u32, stride: u32) {
    // Kernel radius in *stride-bins*
    // If radius is 0, no neighbors -> identity
    let half_window = (bin_size / stride).saturating_sub(1) as usize;
    if half_window == 0 {
        // No overlap region: each bin's average = its coverage
        for b in bins.iter_mut() {
            b.average_overlap_coverage = b.average_coverage;
        }
        return;
    }

    // Precompute the triangular weights-kernel once: `[1,2,...,half_window+1,...,2,1]`
    let weights = triangular_weights(half_window);
    let n = bins.len();

    // Slide the kernel across all centers i
    // It may get truncated at the edges
    for i in 0..n {
        // Choose the usable neighborhood around i, clipped to genome edges
        // Target full window is [i-half_window, i+half_window], but we clamp to [0, n-1]
        let start_i = i.saturating_sub(half_window); // left bound (inclusive)
        let end_i = (i + half_window + 1).min(n); // right bound (exclusive)
        let bin_slice = &bins[start_i..end_i];
        let slice_len = bin_slice.len();

        // Compute where to start reading from `weights` so that its center aligns with `i`
        //
        //    Intuition:
        //    - In the interior (far from edges), we can use the whole kernel, so w_start = 0.
        //    - Near the left edge (small i), we are missing `i - start_i` neighbors on the left.
        //      We therefore skip that many weights from the *left side* of the kernel.
        //
        //    Formal:
        //    - `i - start_i` == number of bins available on the left side.
        //    - The kernel expects `half_window` bins on the left.
        //    - missing_left = half_window - (i - start_i)  (clamped at 0)
        //    - So start at weight index `w_start = missing_left`.
        //
        //    Example: half_window=2, i=0 -> start_i=0 -> i - start_i = 0 -> w_start = 2.
        //             We use weights[2..] = [3,2,1].
        let w_start = half_window.saturating_sub(i - start_i);

        let mut sum_cov = 0.0_f32; // Weighted sum of coverage densities
        let mut sum_w = 0usize; // Sum of integer weights actually used

        // Sum average_coverage * weight, and the weights
        for j in 0..slice_len {
            let w = weights[w_start + j];
            // Last stride-bin may be shorter, so weight by length (most == 1.0)
            let len_ratio = (bin_slice[j].size() as f32) / (stride as f32);
            sum_cov += bin_slice[j].average_coverage * (w as f32) * len_ratio;
            sum_w += w;
        }
        bins[i].average_overlap_coverage = if sum_w > 0 {
            sum_cov / (sum_w as f32)
        } else {
            0.0
        };
    }
}

/// Calculate the `StrideBin::scaling_factor` by dividing `StrideBin::average_overlap_coverage`
/// across all chromosomes.
///
/// Computes a global mean of `average_overlap_coverage` across all supported bins in `bins_by_chr`
/// and divides every bin's `average_overlap_coverage` by that mean so the new global mean is ~1.0.
/// Optionally weight the mean by bin length to better approximate a base-weighted genome mean.
///
/// Parameters
/// ----------
/// - bins_by_chr:
///     Map from chromosome to its stride bins
/// - length_weighted:
///     If true, weight each bin by its length; if false, weight all bins equally
/// - invert:
///     Invert the final scaling factor (1/x).
///     **NOTE**: Zero-values remain zero.
///
/// Returns
/// -------
/// - mean_before:
///     The global mean used for normalization (before scaling)
pub fn normalize_average_overlap_by_global_mean(
    bins_by_chr: &mut FxHashMap<String, Vec<StrideBin>>,
    length_weighted: bool,
    invert: bool,
) -> Result<f32> {
    let mut sum = 0.0_f64;
    let mut wsum = 0.0_f64;

    // Compute global mean over all chromosomes
    for bins in bins_by_chr.values() {
        for b in bins {
            let v = b.average_overlap_coverage as f64;
            if !v.is_finite() || is_effectively_zero(v) {
                continue; // Skip NaN/inf and bins without real support
            }
            let w = if length_weighted {
                b.size() as f64
            } else {
                1.0
            };
            // Skip zero-length bins just in case
            if w == 0.0 {
                continue;
            }
            sum += v * w;
            wsum += w;
        }
    }

    if wsum == 0.0 {
        anyhow::bail!("no bins to normalize or all had length 0");
    }

    let mean = sum / wsum;
    if !mean.is_finite() || mean <= 0.0 {
        anyhow::bail!("invalid global mean {}", mean);
    }

    // Calculate the scaling factors
    let inv_mean = 1.0_f64 / mean;
    for bins in bins_by_chr.values_mut() {
        for b in bins.iter_mut() {
            let v = b.average_overlap_coverage as f64;
            b.scaling_factor = if !v.is_finite() || is_effectively_zero(v) {
                0.0
            } else {
                let normalized = v * inv_mean;
                if invert {
                    (1.0 / normalized) as f32
                } else {
                    normalized as f32
                }
            };
        }
    }

    Ok(mean as f32)
}
