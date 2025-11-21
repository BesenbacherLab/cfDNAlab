///! NOTE: This code was generated. TODO: Validate that it's correct.

/// Fill zero-valued histogram bins by fitting local weighted polynomials.
///
/// Technical details:
/// - Operates entirely in-place and only uses the original non-zero bins as anchors.
/// - Each zero run is interpolated independently to avoid bleeding across cliffs.
/// - Pads edges with zero-value pseudo anchors when all valid neighbours sit on one side.
/// - Panics if the histogram is empty or if either tuning parameter is zero.
///
/// Parameters
/// ----------
/// - `histogram`:
///     Dense 1D histogram to mutate.
/// - `polynomial_degree`:
///     Degree of the fitting polynomial (1 = linear, 2 = quadratic, etc.).
/// - `min_neighbours`:
///     Minimum total anchors (real or padded) required before fitting.
/// - `max_neighbours_per_side`:
///     Cap on how many anchors to take from each side of the zero run (keeps interpolation local).
pub fn fill_zero_bins_with_polynomial(
    histogram: &mut [f64],
    polynomial_degree: usize,
    min_neighbours: usize,
    max_neighbours_per_side: usize,
) {
    assert!(
        !histogram.is_empty(),
        "histogram must contain at least one bin before interpolation"
    );
    assert!(polynomial_degree > 0, "polynomial_degree must be >= 1");
    assert!(min_neighbours > 0, "min_neighbours must be >= 1");
    assert!(
        max_neighbours_per_side > 0,
        "max_neighbours_per_side must be >= 1"
    );

    // Snapshot the original non-zero bins so interpolation never feeds on newly synthesized values
    let anchors: Vec<(usize, f64)> = histogram
        .iter()
        .enumerate()
        .filter_map(|(bin_idx, &count)| {
            if count > 0.0 {
                Some((bin_idx, count))
            } else {
                None
            }
        })
        .collect();

    if anchors.len() < polynomial_degree + 1 {
        // Not enough information to fit anything meaningful
        return;
    }

    let mut cursor_idx = 0;
    while cursor_idx < histogram.len() {
        if histogram[cursor_idx] != 0.0 {
            cursor_idx += 1;
            continue;
        }

        // Identify contiguous zero run [run_start_idx, run_end_idx)
        let run_start_idx = cursor_idx;
        while cursor_idx < histogram.len() && histogram[cursor_idx] == 0.0 {
            cursor_idx += 1;
        }
        let run_end_idx = cursor_idx;

        // Fit once for the entire run and evaluate across the gap
        if let Some(coefficients) = fit_run_polynomial(
            run_start_idx,
            run_end_idx,
            polynomial_degree,
            min_neighbours,
            max_neighbours_per_side,
            &anchors,
        ) {
            for target_idx in run_start_idx..run_end_idx {
                let interpolated_value = evaluate_polynomial(&coefficients, target_idx as f64);
                histogram[target_idx] = interpolated_value.max(0.0); // Guard against tiny negatives
            }
        }
        // Leave untouched when we lack enough neighbours or the fit fails
    }
}

/// Build a polynomial for a zero run using nearby anchors (and optional zero padding).
fn fit_run_polynomial(
    run_start_idx: usize,
    run_end_idx: usize,
    polynomial_degree: usize,
    min_neighbours: usize,
    max_neighbours_per_side: usize,
    anchors: &[(usize, f64)],
) -> Option<Vec<f64>> {
    let required_points = polynomial_degree + 1;
    // Measure anchor distances from the midpoint of the zero run so both sides are comparable
    let run_center = if run_end_idx == run_start_idx {
        run_start_idx as f64
    } else {
        (run_start_idx + run_end_idx - 1) as f64 / 2.0
    };

    // Gather candidate anchors on each side together with their distance from the run center
    let mut left_bins: Vec<(f64, f64, f64)> = anchors
        .iter()
        .filter(|(bin_idx, _)| *bin_idx < run_start_idx)
        .map(|&(bin_idx, count_value)| {
            let distance = run_center - bin_idx as f64;
            (bin_idx as f64, count_value, distance)
        })
        .collect();
    let mut right_bins: Vec<(f64, f64, f64)> = anchors
        .iter()
        .filter(|(bin_idx, _)| *bin_idx >= run_end_idx)
        .map(|&(bin_idx, count_value)| {
            let distance = bin_idx as f64 - run_center;
            (bin_idx as f64, count_value, distance)
        })
        .collect();

    // Prefer the closest true neighbours first to keep interpolation local
    left_bins.sort_by(|lhs, rhs| lhs.2.partial_cmp(&rhs.2).unwrap());
    right_bins.sort_by(|lhs, rhs| lhs.2.partial_cmp(&rhs.2).unwrap());

    let left_take = max_neighbours_per_side.min(left_bins.len());
    let right_take = max_neighbours_per_side.min(right_bins.len());

    // Seed the weighted sample set with bounded neighbour lists from both sides
    let mut weighted_samples: Vec<(f64, f64, f64)> = Vec::with_capacity(left_take + right_take + 2);
    weighted_samples.extend(left_bins.into_iter().take(left_take).map(
        |(gc_idx, count_value, distance)| {
            let weight = 1.0 / (1.0 + distance);
            (gc_idx, count_value, weight)
        },
    ));
    weighted_samples.extend(right_bins.into_iter().take(right_take).map(
        |(gc_idx, count_value, distance)| {
            let weight = 1.0 / (1.0 + distance);
            (gc_idx, count_value, weight)
        },
    ));

    let total_required = required_points.max(min_neighbours);
    if weighted_samples.len() < total_required {
        // Zero-pad the missing side so edge gaps still get a gentle slope
        if left_take == 0 && right_take > 0 {
            weighted_samples.push((run_start_idx as f64 - 1.0, 0.0, 0.5));
        } else if right_take == 0 && left_take > 0 {
            weighted_samples.push((run_end_idx as f64 + 1.0, 0.0, 0.5));
        }
    }

    if weighted_samples.len() < total_required {
        return None;
    }

    fit_weighted_polynomial(weighted_samples.as_slice(), polynomial_degree)
}

#[inline]
fn evaluate_polynomial(coefficients: &[f64], x: f64) -> f64 {
    coefficients
        .iter()
        .enumerate()
        .fold(0.0, |acc, (idx, coefficient)| {
            acc + coefficient * x.powi(idx as i32)
        })
}

fn fit_weighted_polynomial(
    samples: &[(f64, f64, f64)], // (x, y, weight)
    polynomial_degree: usize,
) -> Option<Vec<f64>> {
    let num_coefficients = polynomial_degree + 1;
    let mut normal_matrix = vec![vec![0.0; num_coefficients]; num_coefficients]; // A^T W A
    let mut rhs_vector = vec![0.0; num_coefficients]; // A^T W y

    for &(x_coord, y_value, weight) in samples {
        let mut monomials = vec![1.0; num_coefficients];
        for degree_idx in 1..num_coefficients {
            monomials[degree_idx] = monomials[degree_idx - 1] * x_coord;
        }
        for row_idx in 0..num_coefficients {
            for col_idx in 0..num_coefficients {
                normal_matrix[row_idx][col_idx] += weight * monomials[row_idx] * monomials[col_idx];
            }
            rhs_vector[row_idx] += weight * monomials[row_idx] * y_value;
        }
    }

    solve_sym_posdef(&mut normal_matrix, &mut rhs_vector)
}

/// Solve a small symmetric positive-definite system via Gauss-Jordan elimination.
///
/// This treats `normal_matrix` as the left-hand side (A) and `rhs_vector` as the right-hand side (b),
/// and performs in-place elimination to produce the coefficient vector `x` such that `A * x = b`.
///
/// https://en.wikipedia.org/wiki/Gaussian_elimination
fn solve_sym_posdef(normal_matrix: &mut [Vec<f64>], rhs_vector: &mut [f64]) -> Option<Vec<f64>> {
    let matrix_size = rhs_vector.len();
    for pivot_idx in 0..matrix_size {
        // Pivot
        let pivot_value = normal_matrix[pivot_idx][pivot_idx];
        if pivot_value.abs() < 1e-12 {
            return None;
        }
        for col_idx in pivot_idx..matrix_size {
            normal_matrix[pivot_idx][col_idx] /= pivot_value;
        }
        rhs_vector[pivot_idx] /= pivot_value;

        // Eliminate
        for row_idx in 0..matrix_size {
            if row_idx == pivot_idx {
                continue;
            }
            let elimination_factor = normal_matrix[row_idx][pivot_idx];
            for col_idx in pivot_idx..matrix_size {
                normal_matrix[row_idx][col_idx] -=
                    elimination_factor * normal_matrix[pivot_idx][col_idx];
            }
            rhs_vector[row_idx] -= elimination_factor * rhs_vector[pivot_idx];
        }
    }
    Some(rhs_vector.to_vec())
}
