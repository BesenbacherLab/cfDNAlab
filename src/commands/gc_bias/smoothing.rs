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
pub fn smoothe_counts_gaussian(counts: &Array2<f64>, sigma: f64, radius: usize) -> Array2<f64> {
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

/// Given intuitive mass targets (e.g. 0.5 center, 0.4 for ±1, 0.1 for ±2),
/// fit a sigma that produces a discrete Gaussian matching those weights.
///
/// `targets` should be a slice of length `radius + 1`, where index 0 is the
/// desired mass on the center bin and index `k` is the combined mass for
/// both ±k neighbors.
pub fn fit_sigma_for_targets(radius: usize, targets: &[f64]) -> f64 {
    assert!(
        targets.len() == radius + 1,
        "targets must have entries for center + each offset"
    );
    let total: f64 = targets.iter().sum();
    assert!(
        (total - 1.0).abs() < 1e-6,
        "target weights must sum to 1.0 (got {total})"
    );

    let mut sigma = 1.0;
    for _ in 0..30 {
        let kernel = gaussian_kernel(radius, sigma);
        let mut error = 0.0;
        for offset in 0..=radius {
            let weight = if offset == 0 {
                kernel[radius]
            } else {
                kernel[radius + offset] * 2.0 // combine ±offset
            };
            error += (weight - targets[offset]).powi(2);
        }

        // crude gradient-free adjustment: compare farthest weight to target
        let far_weight = kernel.last().copied().unwrap_or(0.0) * 2.0;
        let far_target = *targets.last().unwrap();
        if far_target > 0.0 {
            let ratio = (far_weight / far_target).max(1e-6);
            sigma /= ratio.sqrt(); // shrink sigma if too broad, expand if too sharp
        }

        if error < 1e-8 {
            break;
        }
    }

    sigma
}
