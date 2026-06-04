use anyhow::{Result, ensure};
use ndarray::Array2;

/// Smoothe a length-by-GC matrix using a separable Gaussian kernel.
///
/// Applies a Gaussian blur by first convolving the rows and then columns.
/// Using a separable (i.e., 2-step) kernel keeps the operation linear in `(rows + cols)`
/// instead of `(rows * cols)` for the same kernel size.
///
/// Parameters
/// ----------
/// - `counts`:
///     Input matrix where rows map to fragment lengths and columns map to GC bins.
/// - `sigma`:
///     Standard deviation for the Gaussian kernel (> 0). Larger values spread
///     mass over wider neighborhoods.
/// - `radius`:
///     Kernel radius in bins (> 0). The kernel spans `2 * radius + 1` bins.
///
/// Returns
/// -------
/// - `Array2<f64>`:
///     Smoothed matrix with the same shape as `counts`.
pub(crate) fn smoothe_counts_gaussian(
    counts: &Array2<f64>,
    sigma: f64,
    radius: usize,
) -> Array2<f64> {
    assert!(sigma > 0.0, "sigma must be positive");
    assert!(radius > 0, "radius must be > 0");

    let (n_rows, n_cols) = counts.dim();
    let kernel = gaussian_kernel(radius, sigma);

    // First pass: Smooth horizontally over each row
    let mut row_pass = Array2::<f64>::zeros((n_rows, n_cols));

    for row in 0..n_rows {
        for col in 0..n_cols {
            let mut acc = 0.0;
            for (k, &weight) in kernel.iter().enumerate() {
                let offset = k as isize - radius as isize;
                let src_col = (col as isize + offset).clamp(0, (n_cols - 1) as isize) as usize;
                acc += (counts[(row, src_col)]) * weight;
            }
            row_pass[(row, col)] = acc;
        }
    }

    // Second pass: Smooth vertically over the intermediate result
    let mut smoothed = Array2::<f64>::zeros((n_rows, n_cols));
    for col in 0..n_cols {
        for row in 0..n_rows {
            let mut acc = 0.0;
            for (k, &weight) in kernel.iter().enumerate() {
                let offset = k as isize - radius as isize;
                let src_row = (row as isize + offset).clamp(0, (n_rows - 1) as isize) as usize;
                acc += row_pass[(src_row, col)] * weight;
            }
            smoothed[(row, col)] = acc.max(0.0);
        }
    }

    smoothed
}

/// Build a normalized 1-D Gaussian kernel for the separable convolution.
///
/// Parameters
/// ----------
/// - `radius`:
///     Half-width of the kernel in bins. The full kernel length becomes `2 *
///     radius + 1`.
/// - `sigma`:
///     Standard deviation of the Gaussian. Controls how quickly the weights
///     decay from the center.
///
/// Returns
/// -------
/// - `Vec<f64>`:
///     Symmetric kernel whose entries sum to one.
fn gaussian_kernel(radius: usize, sigma: f64) -> Vec<f64> {
    let mut kernel = Vec::with_capacity(2 * radius + 1);
    let denom = 2.0 * sigma * sigma;

    for idx in 0..=2 * radius {
        let dist = idx as isize - radius as isize;
        let weight = (-((dist * dist) as f64) / denom).exp();
        kernel.push(weight);
    }

    let sum: f64 = kernel.iter().sum();
    for weight in kernel.iter_mut() {
        *weight /= sum;
    }

    kernel
}

/// Smooth a 1-D row in place with a clamped Gaussian kernel.
pub(crate) fn smooth_row_in_place(row: &mut [f64], sigma: f64, radius: usize) -> Result<()> {
    ensure!(sigma > 0.0, "sigma must be positive");
    ensure!(radius > 0, "radius must be positive");
    if row.is_empty() {
        return Ok(());
    }

    let kernel = gaussian_kernel(radius, sigma);
    let mut tmp = vec![0.0; row.len()];

    let n = row.len();
    for i in 0..n {
        let mut acc = 0.0;
        for (k, &w) in kernel.iter().enumerate() {
            let offset = k as isize - radius as isize; // Are sometimes negative
            let src = reflect_index(i as isize + offset, n);
            acc += row[src] * w;
        }
        tmp[i] = acc;
    }

    row.copy_from_slice(&tmp);
    Ok(())
}

/// Slice the row for `length` from a flat counts buffer and smooth it in place.
pub(crate) fn smooth_length_row_in_place(
    counts: &mut [f64],
    offsets: &[usize],
    length_min: usize,
    length: usize,
    sigma: f64,
    radius: usize,
) -> Result<()> {
    ensure!(length >= length_min);

    let row_idx = length - length_min;
    ensure!(row_idx + 1 < offsets.len());

    let start = offsets[row_idx];
    let end = offsets[row_idx + 1];
    ensure!(!(start >= end || end > counts.len()));

    // Split to get a mutable slice for this row
    let (_, tail) = counts.split_at_mut(start);
    let (row, _) = tail.split_at_mut(end - start);
    smooth_row_in_place(row, sigma, radius)?;
    Ok(())
}

#[inline]
fn reflect_index(i: isize, len: usize) -> usize {
    debug_assert!(len > 0);
    if i < 0 {
        (-i) as usize
    } else if i as usize >= len {
        (2 * len - 2).saturating_sub(i as usize)
    } else {
        i as usize
    }
}

#[cfg(test)]
mod tests {
    include!("smoothing_tests.rs");
}
