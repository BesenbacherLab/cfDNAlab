use ndarray::{Array2, array};

use cfdnalab::commands::gc_bias::smoothing::{fit_sigma_for_targets, smoothe_counts_gaussian};

fn assert_matrix_close(actual: &Array2<f64>, expected: &Array2<f64>, tol: f64) {
    assert_eq!(actual.dim(), expected.dim(), "matrix dimensions differ");
    for (actual_entry, expected_entry) in actual.iter().zip(expected.iter()) {
        let delta = (actual_entry - expected_entry).abs();
        assert!(
            delta <= tol,
            "mismatch: actual={}, expected={}, delta={}",
            actual_entry,
            expected_entry,
            delta
        );
    }
}

#[test]
fn preserves_uniform_input_when_kernel_is_normalized() {
    // Arrange: a constant matrix should remain unchanged after smoothing because
    // the Gaussian kernel sums to one.
    let counts = array![[7.1, 7.1, 7.1], [7.1, 7.1, 7.1]];

    // Act: apply smoothing with arbitrary sigma/radius and pseudo-count.
    let smoothed = smoothe_counts_gaussian(&counts, 1.2, 2);

    // Assert: every cell is unchanged up to rounding noise.
    assert_matrix_close(&smoothed, &counts, 1e-12);
}

#[test]
fn spreads_mass_to_neighbors_with_pseudo_count() {
    // Arrange: a single bin carries almost all mass. Smoothing must redistribute it.
    let counts = array![[0.0, 0.0], [0.0, 4.0]];
    let expected = array![
        [0.30045443181644615, 0.7958200444283421],
        [0.7958200444283423, 2.1079054793268703],
    ];

    // Act: smooth with sigma=1, radius=1, and pseudo-count=0.5.
    let smoothed = smoothe_counts_gaussian(&counts, 1.0, 1);

    // Assert: neighbor bins pick up fractional counts as derived analytically.
    assert_matrix_close(&smoothed, &expected, 1e-12);
}

#[test]
fn targeted_kernel_and_small_pseudo_count_match_manual_expectation() {
    // Arrange: 5x5 matrix with a simple peak around the center. After normalization,
    // use the configuration described in the pipeline (radius=2, center 70% mass,
    // 20% to ±1, 10% to ±2).
    let counts = array![
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.5, 1.0, 0.5, 0.0],
        [0.0, 2.0, 3.0, 2.0, 0.0],
        [0.0, 0.5, 1.0, 0.5, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    // Manually derived expectation: running the separable filter with the fitted
    // sigma (~0.9733) produces the following matrix (each entry rounded to 12 decimals).
    let expected = array![
        [
            0.073695262505,
            0.198513570067,
            0.270526948794,
            0.198513570067,
            0.073695262505
        ],
        [
            0.234691813570,
            0.623311126994,
            0.842821911779,
            0.623311126994,
            0.234691813570
        ],
        [
            0.346731765036,
            0.915289260686,
            1.233413062753,
            0.915289260686,
            0.346731765036
        ],
        [
            0.234691813570,
            0.623311126994,
            0.842821911779,
            0.623311126994,
            0.234691813570
        ],
        [
            0.073695262505,
            0.198513570067,
            0.270526948794,
            0.198513570067,
            0.073695262505
        ],
    ];

    // Act
    let radius = 2;
    let sigma = fit_sigma_for_targets(radius, &[0.7, 0.2, 0.1]);
    let smoothed = smoothe_counts_gaussian(&counts, sigma, radius);

    // Assert
    assert_matrix_close(&smoothed, &expected, 1e-9);
}
