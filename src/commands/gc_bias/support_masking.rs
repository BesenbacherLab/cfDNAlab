use crate::commands::gc_bias::counting::calculate_gc_bin;
use ndarray::{Array2, ArrayBase, Data, Ix2, Zip};

pub struct StatsBySupportMask {
    pub sum_for_supported: f64,
    pub sum_for_unsupported: f64,
    pub n_supported: u64,
    pub n_unsupported: u64,
}

/// Get count and value-sums for all supported/unsupported bins.
pub fn stats_by_support_mask<S, M>(
    matrix: &ArrayBase<S, Ix2>,
    support_mask: &ArrayBase<M, Ix2>,
) -> StatsBySupportMask
where
    S: Data<Elem = f64>,
    M: Data<Elem = bool>,
{
    assert_eq!(
        matrix.dim(),
        support_mask.dim(),
        "Mask shape {:?} must match matrix shape {:?}",
        support_mask.dim(),
        matrix.dim()
    );

    let mut total_supported = 0.0;
    let mut total_unsupported = 0.0;
    let mut count_supported = 0;
    let mut count_unsupported = 0;

    Zip::from(matrix)
        .and(support_mask)
        .for_each(|value, &is_supported| {
            if is_supported {
                total_supported += *value;
                count_supported += 1;
            } else {
                total_unsupported += *value;
                count_unsupported += 1;
            }
        });

    StatsBySupportMask {
        sum_for_supported: total_supported,
        sum_for_unsupported: total_unsupported,
        n_supported: count_supported,
        n_unsupported: count_unsupported,
    }
}

/// Build mask that supports all elements that are not part of:
///
///   - The N most extreme GC bins on each side.
///   - The M shortest-length bins (often very sparse below ~70-100bp).
pub fn build_extreme_bins_support_mask(
    shape: (usize, usize),
    extreme_gc_bins_per_side: usize,
    short_length_bins: usize,
) -> Array2<bool> {
    let (num_length_bins, num_gc_bins) = shape;
    let gc_bins_to_mask = extreme_gc_bins_per_side.min(num_gc_bins);
    let length_bins_to_mask = short_length_bins.min(num_length_bins);

    if gc_bins_to_mask == 0 && length_bins_to_mask == 0 {
        return Array2::from_elem(shape, true);
    }

    let mut mask = Array2::from_elem(shape, true);

    if gc_bins_to_mask > 0 && num_gc_bins > 0 {
        let right_start = num_gc_bins - gc_bins_to_mask;
        for mut row in mask.outer_iter_mut() {
            for col_idx in 0..gc_bins_to_mask {
                row[col_idx] = false;
            }
            for col_idx in right_start..num_gc_bins {
                row[col_idx] = false;
            }
        }
    }

    if length_bins_to_mask > 0 && num_length_bins > 0 {
        for row_idx in 0..length_bins_to_mask {
            mask.row_mut(row_idx).fill(false);
        }
    }

    mask
}

pub fn set_masked_entries_to_value(matrix: &mut Array2<f64>, mask: &Array2<bool>, fill_value: f64) {
    Zip::from(matrix).and(mask).for_each(|value, &is_valid| {
        if !is_valid {
            *value = fill_value;
        }
    });
}

/* Reference-based masks */

/// Create mask of supported elements. Elements are usable
/// when they have a count of at least `threshold_per_mb`
/// per 1Mb of valid ACGT positions in the selected regions
/// of the genome.
///
/// **NOTE**: This does not consider the number of sampled starts.
/// The idea is that some elements are almost non-existent
/// (e.g. 100% GC in an 800bp fragment interval), so no matter
/// the number of sampled starts they will have almost no counts.
pub fn create_support_mask_threshold_per_mb(
    counts: &[Array2<f64>],
    num_acgt_positions: u64,
    threshold_per_mb: f64,
) -> Option<Array2<bool>> {
    let global_counts = sum_arrays(counts)?;

    // Need at least a count of `threshold_per_mb` per 1Mb valid positions
    let threshold = num_acgt_positions as f64 / 1000000 as f64 * threshold_per_mb;

    // Create mask of usable elements
    let mut mask = Array2::from_elem(global_counts.dim(), true);
    for ((row, col), &value) in global_counts.indexed_iter() {
        mask[(row, col)] = value >= threshold;
    }

    Some(mask)
}

/// Create mask of usable elements. Elements are usable
/// when they have a non-zero count in any of the windows.
pub fn create_support_mask(counts: &[Array2<f64>]) -> Option<Array2<bool>> {
    let global_counts = sum_arrays(counts)?;

    // Create mask of usable elements
    let mut mask = Array2::from_elem(global_counts.dim(), true);
    for ((row, col), &value) in global_counts.indexed_iter() {
        mask[(row, col)] = value > 0.;
    }

    Some(mask)
}

/// Sum a list of matrices.
fn sum_arrays(arrays: &[Array2<f64>]) -> Option<Array2<f64>> {
    let mut iter = arrays.iter();

    let mut sum = iter.next().cloned()?;

    for arr in iter {
        debug_assert_eq!(
            sum.dim(),
            arr.dim(),
            "All array components must share shape"
        );

        Zip::from(&mut sum).and(arr).for_each(|s, &v| *s += v);
    }
    Some(sum)
}

pub fn build_theoretical_support_mask(
    length_min: usize,
    length_max: usize,
    gc_min: usize,
    gc_max: usize,
    end_offset: usize,
) -> Array2<bool> {
    assert!(
        length_max >= length_min,
        "length range must be non-empty ({}..={})",
        length_min,
        length_max
    );
    assert!(
        gc_max >= gc_min,
        "GC bin range must be non-empty ({}..={})",
        gc_min,
        gc_max
    );

    let num_lengths = length_max - length_min + 1;
    let num_gc_bins = gc_max - gc_min + 1;
    let mut mask = Array2::from_elem((num_lengths, num_gc_bins), false);

    for length in length_min..=length_max {
        if length == 0 {
            continue;
        }
        let row_idx = length - length_min;
        let effective_length = length.saturating_sub(2 * end_offset);
        if effective_length == 0 {
            continue;
        }
        let acgt_count = effective_length as u64;
        for gc_count in 0..=effective_length {
            // Use the same integer rounding as the ref-gc-bias tool!
            let gc_bin = calculate_gc_bin(gc_count as u64, acgt_count) as u64;
            if gc_bin < gc_min as u64 {
                continue;
            }
            let col_idx = (gc_bin as usize) - gc_min;
            mask[(row_idx, col_idx)] = true;
        }
    }

    mask
}
