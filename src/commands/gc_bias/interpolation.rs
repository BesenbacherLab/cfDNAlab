use anyhow::{Result, ensure};

/// Interpolate unsupported bins using nearby supported anchors.
///
/// Treats every `false` entry in `support_mask` as a gap and fits a weighted
/// polynomial when at least `polynomial_degree + 1` genuine supported anchors
/// exist in the original row or column. Interpolated values are clamped to the
/// anchor range and can optionally mark the mask as supported. Edge runs can be
/// interpolated from one-sided real support by adding mirrored or boundary-valued
/// pseudo anchors. Runs are skipped only when the real anchor set is too small or
/// the local fit fails.
///
/// Parameters
/// ----------
/// - `histogram`:
///   Dense 1D histogram to mutate.
/// - `support_mask`:
///   Boolean mask flagging bins backed by real data (`true`) or unsupported gaps (`false`).
/// - `polynomial_degree`:
///   Degree of the fitting polynomial (1 = linear, 2 = quadratic, etc.).
/// - `min_neighbours`:
///   Minimum number of supported neighbouring bins required to fit an unsupported run.
///   A neighbouring bin is a bin outside the unsupported run with
///   `support_mask == true`. The closest neighbouring bins on the left and
///   right are used first. If one side has too few real neighbours, the fit
///   first adds artificial bins at distances mirrored from real neighbours on
///   the opposite side. If more bins are still needed, it adds artificial bins
///   farther outward using the nearest available boundary value.
/// - `max_neighbours_per_side`:
///   Cap on how many anchors to take from each side of the unsupported run.
///
/// Returns
/// -------
/// - `Ok(())`:
///   All unsupported runs were either interpolated or skipped due to insufficient anchors.
/// - `Err`:
///   Validation failed (mask length mismatch or invalid tuning parameters).
pub(crate) fn fill_unsupported_bins_with_polynomial(
    histogram: &mut [f64],
    support_mask: &mut [bool],
    polynomial_degree: usize,
    min_neighbours: usize,
    max_neighbours_per_side: usize,
    update_mask: bool,
) -> Result<()> {
    ensure!(
        !histogram.is_empty(),
        "histogram must contain at least one bin before interpolation"
    );
    ensure!(
        histogram.len() == support_mask.len(),
        "support mask must match histogram length"
    );
    ensure!(polynomial_degree > 0, "polynomial_degree must be >= 1");
    ensure!(min_neighbours > 0, "min_neighbours must be >= 1");
    ensure!(
        max_neighbours_per_side > 0,
        "max_neighbours_per_side must be >= 1"
    );

    // Capture the original support set once
    let anchors = collect_supported_anchors(histogram, support_mask);
    if anchors.len() < polynomial_degree + 1 {
        return Ok(());
    }

    let mut cursor_idx = 0;
    while cursor_idx < histogram.len() {
        // Skip spans that already have real support so we never overwrite measured values
        if support_mask[cursor_idx] {
            cursor_idx += 1;
            continue;
        }

        // Treat the unsupported stretch as a single interpolation run
        let run_start_idx = cursor_idx;
        while cursor_idx < histogram.len() && !support_mask[cursor_idx] {
            cursor_idx += 1;
        }
        let run_end_idx = cursor_idx;

        let left_anchor = if run_start_idx > 0 {
            histogram[run_start_idx - 1]
        } else {
            histogram.get(run_end_idx).copied().unwrap_or(0.0)
        };
        let right_anchor = histogram.get(run_end_idx).copied().unwrap_or(left_anchor);
        let lower_bound = left_anchor.min(right_anchor);
        let upper_bound = left_anchor.max(right_anchor);

        if let Some(polynomial_fit) = fit_run_polynomial(
            run_start_idx,
            run_end_idx,
            polynomial_degree,
            min_neighbours,
            max_neighbours_per_side,
            &anchors,
        ) {
            let mut any_value_updates = false;
            for target_idx in run_start_idx..run_end_idx {
                let mut interpolated_value =
                    evaluate_polynomial_fit(&polynomial_fit, target_idx as f64);
                if lower_bound != upper_bound {
                    interpolated_value = interpolated_value.max(lower_bound).min(upper_bound);
                } else {
                    interpolated_value = lower_bound;
                }
                let new_value = interpolated_value.max(0.0);
                if (histogram[target_idx] - new_value).abs() > f64::EPSILON {
                    histogram[target_idx] = new_value;
                    any_value_updates = true;
                }
                if update_mask {
                    support_mask[target_idx] = true;
                }
            }

            if any_value_updates {
                enforce_monotonic_segment(
                    &mut histogram[run_start_idx..run_end_idx],
                    left_anchor,
                    right_anchor,
                );
            }
        }
    }
    Ok(())
}

struct PolynomialFit {
    coefficients: Vec<f64>,
    center: f64,
    scale: f64,
}

/// Build a polynomial for a zero run using nearby anchors (and optional zero padding).
fn fit_run_polynomial(
    run_start_idx: usize,
    run_end_idx: usize,
    polynomial_degree: usize,
    min_neighbours: usize,
    max_neighbours_per_side: usize,
    anchors: &[(usize, f64)],
) -> Option<PolynomialFit> {
    let required_points = polynomial_degree + 1;
    // Measure anchor distances from the midpoint of the zero run so both sides are comparable
    let run_center = if run_end_idx == run_start_idx {
        run_start_idx as f64
    } else {
        (run_start_idx + run_end_idx - 1) as f64 / 2.0
    };

    // Gather candidate anchors on each side together with their distance from the run center
    let mut left_real: Vec<(f64, f64, f64)> = anchors
        .iter()
        .filter(|(bin_idx, _)| *bin_idx < run_start_idx)
        .map(|&(bin_idx, count_value)| {
            let distance = run_center - bin_idx as f64;
            (bin_idx as f64, count_value, distance)
        })
        .collect();
    let mut right_real: Vec<(f64, f64, f64)> = anchors
        .iter()
        .filter(|(bin_idx, _)| *bin_idx >= run_end_idx)
        .map(|&(bin_idx, count_value)| {
            let distance = bin_idx as f64 - run_center;
            (bin_idx as f64, count_value, distance)
        })
        .collect();

    // Prefer the closest true neighbours first to keep interpolation local
    left_real.sort_by(|lhs, rhs| lhs.2.total_cmp(&rhs.2));
    right_real.sort_by(|lhs, rhs| lhs.2.total_cmp(&rhs.2));
    let left_real_available = !left_real.is_empty();
    let right_real_available = !right_real.is_empty();

    // Take up to the requested number of neighbours on each side
    let mut left_selected: Vec<(f64, f64, f64)> = left_real
        .iter()
        .take(max_neighbours_per_side)
        .cloned()
        .collect();
    let mut right_selected: Vec<(f64, f64, f64)> = right_real
        .iter()
        .take(max_neighbours_per_side)
        .cloned()
        .collect();

    // Mirror zero-valued pseudo anchors to fill missing slots from the opposite side
    // When one side is missing, use the left/right boundary values to mirror
    let left_boundary = anchors
        .iter()
        .rev()
        .find(|(idx, _)| *idx < run_start_idx)
        .map(|(_, v)| *v)
        .unwrap_or(0.0);
    let right_boundary = anchors
        .iter()
        .find(|(idx, _)| *idx >= run_end_idx)
        .map(|(_, v)| *v)
        .unwrap_or(0.0);

    if left_selected.len() < max_neighbours_per_side {
        let needed = max_neighbours_per_side - left_selected.len();
        for (_, _, distance) in right_real.iter().take(needed) {
            left_selected.push((run_center - *distance, left_boundary, *distance));
        }
    }
    if right_selected.len() < max_neighbours_per_side {
        let needed = max_neighbours_per_side - right_selected.len();
        for (_, _, distance) in left_real.iter().take(needed) {
            right_selected.push((run_center + *distance, right_boundary, *distance));
        }
    }

    // If we still have gaps (not enough anchors overall),
    // extend outward with evenly spaced anchors (likely zeros)
    while left_selected.len() < max_neighbours_per_side {
        let dist = left_selected
            .last()
            .map(|(_, _, distance)| *distance + 1.0)
            .unwrap_or(1.0);
        left_selected.push((run_center - dist, left_boundary, dist));
    }
    while right_selected.len() < max_neighbours_per_side {
        let dist = right_selected
            .last()
            .map(|(_, _, distance)| *distance + 1.0)
            .unwrap_or(1.0);
        right_selected.push((run_center + dist, right_boundary, dist));
    }

    if left_selected.is_empty() && right_selected.is_empty() {
        return None;
    }

    // Seed the weighted sample set with bounded neighbour lists from both sides
    let mut weighted_samples: Vec<(f64, f64, f64)> =
        Vec::with_capacity(left_selected.len() + right_selected.len() + 2);
    for (gc_idx, count_value, distance) in &left_selected {
        let weight = 1.0 / (1.0 + (*distance * *distance));
        weighted_samples.push((*gc_idx, *count_value, weight));
    }
    for (gc_idx, count_value, distance) in &right_selected {
        let weight = 1.0 / (1.0 + (*distance * *distance));
        weighted_samples.push((*gc_idx, *count_value, weight));
    }

    // Normalize weights to sum to 1.0 so the solver retains scale sensitivity.
    let weight_sum: f64 = weighted_samples.iter().map(|(_, _, w)| *w).sum();
    if weight_sum > 0.0 {
        for (_, _, w) in &mut weighted_samples {
            *w /= weight_sum;
        }
    }

    let total_required = required_points.max(min_neighbours);
    if weighted_samples.len() < total_required {
        // When both sides run out of real/mirrored anchors (e.g., isolated runs),
        // synthesize evenly spaced zero points moving outward so the polynomial
        // system remains solvable. The absolute spacing is arbitrary; a linear
        // step keeps the code simple and biases the fit minimally.
        let mut extra = 1.0;
        while weighted_samples.len() < total_required {
            extra += 1.0;
            let weight = 1.0 / (1.0 + extra);
            if right_real_available {
                weighted_samples.push((run_center - extra, 0.0, weight));
            }
            if left_real_available {
                weighted_samples.push((run_center + extra, 0.0, weight));
            }
            if !left_real_available && !right_real_available {
                weighted_samples.push((run_center - extra, 0.0, weight));
                if weighted_samples.len() >= total_required {
                    break;
                }
                weighted_samples.push((run_center + extra, 0.0, weight));
            }
        }
    }

    // Center and scale x before building monomial normal equations. This keeps
    // powers like x^4 small without changing the fitted curve being evaluated.
    let coordinate_scale = weighted_samples
        .iter()
        .map(|(gc_idx, _, _)| (gc_idx - run_center).abs())
        .fold(0.0, f64::max)
        .max(1.0);
    let centered_samples: Vec<(f64, f64, f64)> = weighted_samples
        .iter()
        .map(|(gc_idx, count_value, weight)| {
            (
                (gc_idx - run_center) / coordinate_scale,
                *count_value,
                *weight,
            )
        })
        .collect();

    fit_weighted_polynomial(centered_samples.as_slice(), polynomial_degree).map(|coefficients| {
        PolynomialFit {
            coefficients,
            center: run_center,
            scale: coordinate_scale,
        }
    })
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

#[inline]
fn evaluate_polynomial_fit(polynomial_fit: &PolynomialFit, bin_idx: f64) -> f64 {
    let local_x = (bin_idx - polynomial_fit.center) / polynomial_fit.scale;
    evaluate_polynomial(&polynomial_fit.coefficients, local_x)
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

pub(crate) fn enforce_monotonic_segment(segment: &mut [f64], left_anchor: f64, right_anchor: f64) {
    if segment.is_empty() || (left_anchor - right_anchor).abs() < f64::EPSILON {
        return;
    }

    if left_anchor < right_anchor {
        let mut prev = left_anchor;
        for value in segment.iter_mut() {
            if *value < prev {
                *value = prev;
            }
            prev = *value;
        }
    } else if left_anchor > right_anchor {
        let mut prev = left_anchor;
        for value in segment.iter_mut() {
            if *value > prev {
                *value = prev;
            }
            prev = *value;
        }
    }
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

/// Extract `(index, value)` pairs for every bin backed by real data.
///
/// Used by both interpolation routines to keep the anchor list limited to genuine counts.
fn collect_supported_anchors(histogram: &[f64], support_mask: &[bool]) -> Vec<(usize, f64)> {
    histogram
        .iter()
        .zip(support_mask.iter())
        .enumerate()
        .filter_map(
            |(idx, (&value, &supported))| {
                if supported { Some((idx, value)) } else { None }
            },
        )
        .collect()
}

#[cfg(test)]
mod tests {
    include!("interpolation_tests.rs");
}
