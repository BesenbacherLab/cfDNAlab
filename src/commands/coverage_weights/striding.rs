use anyhow::Result;
use fxhash::FxHashMap;

use crate::shared::interval::Interval;

/// Treat support smaller than this as zero so tiny numerical residues do not participate
/// in the global mean or blow up during inversion
const SUPPORT_FLOOR: f64 = 1e-10;

/// Represents a single stride bin and its triangularly smoothed value.
///
/// The smoothed value is the weighted average of the centered stride value and its
/// surrounding stride values under the triangular kernel used by `fill_triangular_overlap`.
/// The same field can represent smoothed average coverage or a smoothed
/// fragment count, depending on which command produced the raw stride value.
#[derive(Debug, Clone, Copy)]
pub(crate) struct StrideBin {
    /// Checked genomic span of the stride-bin
    pub(crate) interval: Interval<u32>,
    /// Number of non-blacklisted bases that support the stride-bin value
    pub(crate) eligible_positions: u32,
    /// Eligible support as a fraction of the configured stride.
    ///
    /// This is precomputed when stride bins are loaded so smoothing can weight
    /// short final bins and partly blacklisted bins without recalculating it
    /// for every overlapping kernel position.
    pub(crate) support_ratio: f64,
    /// Raw command-specific value for this stride bin.
    ///
    /// This is average coverage for `coverage-weights` and a fractional  
    /// fragment count for `fragment-count-weights`.
    pub(crate) stride_value: f32,
    /// Triangular weighted average of the stride values from the center and
    /// surrounding stride bins
    pub(crate) smoothed_value: f32,
    /// Multiplicative scaling factor for normalizing across the genome
    pub(crate) scaling_factor: f32,
}

impl StrideBin {
    /// Return the inclusive start coordinate of the stride-bin.
    #[inline]
    pub(crate) fn start(&self) -> u32 {
        self.interval.start()
    }

    /// Return the exclusive end coordinate of the stride-bin.
    #[inline]
    pub(crate) fn end(&self) -> u32 {
        self.interval.end()
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

/// Fill `bins[i].smoothed_value` with a triangularly weighted average
/// of the center stride value and surrounding stride values.
///
/// Each triangular kernel weight is multiplied by the bin's `support_ratio`, so short final bins
/// and partly blacklisted bins contribute less than full, unmasked stride bins.
///
/// Goal: For each stride-bin `i`, approximate the average value across all overlapping
/// "megabins" (large windows of size `bin_size`) **without** explicitly enumerating them.
///
/// The triangular kernel is aligned so its center weight lands on stride bin `i`.
/// Each kernel weight for the surrounding positions counts how many megabins, that overlap
/// `i`, also overlap the stride bin at that relative position.
///
/// Key quantities:
/// - `half_window = (bin_size / stride) - 1`
///
///     The kernel radius, *measured in stride-bins*. The full kernel length is
///   `2*half_window + 1`. Example: if bin_size=5Mb and stride=0.5Mb,
///   then bin_size/stride = 10, half_window = 9, kernel length = 19.
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
/// 3) Accumulate `sum(stride_value[j] * weight[j])` and `sum(weights)`.
///    Non-finite values are missing measurements and do not contribute to either sum.
///    The weight includes the triangular kernel weight and eligible-base support.
/// 4) Normalize: `smoothed_value[i] = weighted_sum / sum_weights`.
///    If no finite neighbor is available, the smoothed value is `NaN`.
///
/// Parameters
/// ----------
/// - bins:
///   Stride bins with `stride_value` set to the command-specific value:
///   average coverage for `coverage-weights` and fragment count for
///   `fragment-count-weights`.
///   `support_ratio` controls how much each stride contributes to smoothing.
/// - bin_size:
///   Large window size; used only to derive the kernel radius.
/// - stride:
///   Stride size; used only to derive the kernel radius.
pub(crate) fn fill_triangular_overlap(bins: &mut Vec<StrideBin>, bin_size: u32, stride: u32) {
    // Kernel radius in *stride-bins*
    // If radius is 0, no neighbors -> identity
    let half_window = (bin_size / stride).saturating_sub(1) as usize;
    if half_window == 0 {
        // No overlap region: each bin keeps its raw stride value
        for b in bins.iter_mut() {
            b.smoothed_value = b.stride_value;
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

        let mut weighted_value_sum = 0.0_f64;
        let mut weight_sum = 0.0_f64;

        // Sum stride value * weight, and the weights
        for j in 0..slice_len {
            let mut w = weights[w_start + j] as f64;
            let stride_value = bin_slice[j].stride_value as f64;
            if !stride_value.is_finite() {
                continue;
            }
            // Eligible support handles both short final bins and blacklist-masked bases
            if bin_slice[j].support_ratio == 0.0 {
                continue;
            }
            w *= bin_slice[j].support_ratio;
            weighted_value_sum += stride_value * w;
            weight_sum += w;
        }
        bins[i].smoothed_value = if weight_sum > 0.0 {
            (weighted_value_sum / weight_sum) as f32
        } else {
            f32::NAN
        };
    }
}

/// Calculate `StrideBin::scaling_factor` from each bin's smoothed value.
///
/// Here, "smoothed value" means the triangular weighted average already stored in
/// `StrideBin::smoothed_value` by `fill_triangular_overlap`.
///
/// This function computes a global mean across all supported bins in `bins_by_chr` and divides every
/// bin's smoothed value by that mean so the new global mean is ~1.0. Optionally weight the
/// mean by eligible bases to better approximate a base-weighted genome mean.
///
/// Parameters
/// ----------
/// - bins_by_chr:
///   Map from chromosome to its stride bins
/// - length_weighted:
///   If true, weight each bin by eligible positions; if false, weight all bins equally
/// - invert:
///   Invert the final scaling factor (1/x).
///   **NOTE**: Zero-values remain zero.
///
/// Returns
/// -------
/// - mean_before:
///   The global mean used for normalization (before scaling)
pub(crate) fn normalize_weighted_average_overlap_by_global_mean(
    bins_by_chr: &mut FxHashMap<String, Vec<StrideBin>>,
    length_weighted: bool,
    invert: bool,
) -> Result<f32> {
    let mut sum = 0.0_f64;
    let mut wsum = 0.0_f64;
    let mut total_bins = 0usize;
    let mut usable_bins = 0usize;

    // Compute global mean over all chromosomes.
    // A smoothed value can be finite even when the raw stride value is missing because smoothing
    // can interpolate across masked bins. Those rows are useful in the output, but they should not
    // define the global mean or receive a usable scaling factor.
    for bins in bins_by_chr.values() {
        for b in bins {
            total_bins += 1;
            let raw_value = b.stride_value as f64;
            let smoothed_value = b.smoothed_value as f64;
            if !raw_value.is_finite()
                || !smoothed_value.is_finite()
                || is_effectively_zero(smoothed_value)
                || b.eligible_positions == 0
            {
                continue; // Skip NaN/inf and bins without real support
            }
            usable_bins += 1;
            let w = if length_weighted {
                b.eligible_positions as f64
            } else {
                1.0
            };
            // Skip zero-length bins just in case
            if w == 0.0 {
                continue;
            }
            sum += smoothed_value * w;
            wsum += w;
        }
    }

    if wsum == 0.0 {
        if total_bins == 0 {
            anyhow::bail!("no stride bins were available to normalize");
        }
        if usable_bins == 0 {
            anyhow::bail!(
                "no usable finite non-zero smoothed fragment mass after filtering across {} stride bins. Check --chromosomes, --min-mapq, fragment length filters, blacklist, and GC correction inputs",
                total_bins
            );
        }
        anyhow::bail!(
            "internal error: total sum of 0 but found {} usable bins. Should be impossible, please report",
            usable_bins
        );
    }

    let mean = sum / wsum;
    if !mean.is_finite() || mean <= 0.0 {
        anyhow::bail!("invalid global mean {}", mean);
    }

    // Calculate the scaling factors
    let inv_mean = 1.0_f64 / mean;
    for bins in bins_by_chr.values_mut() {
        for b in bins.iter_mut() {
            let raw_value = b.stride_value as f64;
            let smoothed_value = b.smoothed_value as f64;
            b.scaling_factor = if !raw_value.is_finite()
                || !smoothed_value.is_finite()
                || is_effectively_zero(smoothed_value)
                || b.eligible_positions == 0
            {
                0.0
            } else {
                let normalized = smoothed_value * inv_mean;
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
