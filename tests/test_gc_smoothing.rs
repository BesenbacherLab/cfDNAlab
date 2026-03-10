use fxhash::FxHashMap;
use ndarray::{Array2, array};

use cfdnalab::commands::gc_bias::{
    binning::{BinnedAxis, CollapseAggregation, collapse_counts_by_bins},
    smoothing::smoothe_counts_gaussian,
};

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

/// Build a `BinnedAxis` that collapses every index into a single bin.
/// Used for small matrices where we want to merge all rows (or columns) together.
fn build_single_bin_axis(size: usize) -> BinnedAxis {
    let mut index_to_bin = FxHashMap::default();
    let mut bin_to_indices = FxHashMap::default();
    let indices: Vec<usize> = (0..size).collect();
    bin_to_indices.insert(0, indices.clone());
    for idx in indices {
        index_to_bin.insert(idx, 0);
    }
    BinnedAxis {
        index_to_bin,
        bin_to_indices,
        num_bins: 1,
    }
}

/// Build a `BinnedAxis` from explicit bin-to-index mappings (no greedy behavior).
/// Each tuple is `(bin_idx, indices)`; indices are grouped exactly as provided.
fn build_explicit_bins(mapping: &[(usize, Vec<usize>)]) -> BinnedAxis {
    let mut index_to_bin = FxHashMap::default();
    let mut bin_to_indices = FxHashMap::default();
    for (bin_idx, indices) in mapping {
        bin_to_indices.insert(*bin_idx, indices.clone());
        for &idx in indices {
            index_to_bin.insert(idx, *bin_idx);
        }
    }
    BinnedAxis {
        index_to_bin,
        bin_to_indices,
        num_bins: mapping.len(),
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
    assert_matrix_close(&smoothed, &expected, 1e-3);
}

#[test]
fn targeted_kernel_and_small_pseudo_count_match_manual_expectation() {
    // Arrange: 5x5 matrix with a simple peak around the center.
    let counts = array![
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.5, 1.0, 0.5, 0.0],
        [0.0, 2.0, 3.0, 2.0, 0.0],
        [0.0, 0.5, 1.0, 0.5, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
    ];

    // Gaussian settings
    let radius = 2;
    let sigma = 1.0;
    assert_eq!(
        radius, 2,
        "derivation below assumes radius 2 for explicit terms"
    );

    // Re-derive weights once so the algebra stays transparent.
    // Unnormalized 1-D Gaussian samples for offsets [-radius, ..., radius].
    let weight = |d: f64| (-(d * d) / (2.0 * sigma * sigma)).exp();
    let g2 = weight(radius as f64); // offset +/-2
    let g1 = weight(1.0); // offset +/-1
    let g0 = weight(0.0); // center
    let norm = 2.0 * g2 + 2.0 * g1 + g0;
    let w2 = g2 / norm;
    let w1 = g1 / norm;
    let w0 = g0 / norm;

    // Horizontal pass worked out by hand for each column and row (clamped edges).
    // let row0_h = [0.0; 5];
    let row1_h = [
        w1 * 0.5 + w2 * 1.0,
        w0 * 0.5 + w1 * 1.0 + w2 * 0.5,
        w0 * 1.0 + w1 * 0.5 + w1 * 0.5,
        w0 * 0.5 + w1 * 1.0 + w2 * 0.5,
        w1 * 0.5 + w2 * 1.0,
    ];
    let row2_h = [
        w1 * 2.0 + w2 * 3.0,
        w0 * 2.0 + w1 * 3.0 + w2 * 2.0,
        w1 * 2.0 + w0 * 3.0 + w1 * 2.0,
        w0 * 2.0 + w1 * 3.0 + w2 * 2.0,
        w1 * 2.0 + w2 * 3.0,
    ];
    let row3_h = row1_h;
    // let row4_h = row0_h;

    // Vertical pass unrolled with the same clamping rule.
    let expected_calc = array![
        [
            w1 * row1_h[0] + w2 * row2_h[0],
            w1 * row1_h[1] + w2 * row2_h[1],
            w1 * row1_h[2] + w2 * row2_h[2],
            w1 * row1_h[3] + w2 * row2_h[3],
            w1 * row1_h[4] + w2 * row2_h[4],
        ],
        [
            w0 * row1_h[0] + w1 * row2_h[0] + w2 * row3_h[0],
            w0 * row1_h[1] + w1 * row2_h[1] + w2 * row3_h[1],
            w0 * row1_h[2] + w1 * row2_h[2] + w2 * row3_h[2],
            w0 * row1_h[3] + w1 * row2_h[3] + w2 * row3_h[3],
            w0 * row1_h[4] + w1 * row2_h[4] + w2 * row3_h[4],
        ],
        [
            w1 * row1_h[0] + w0 * row2_h[0] + w1 * row3_h[0],
            w1 * row1_h[1] + w0 * row2_h[1] + w1 * row3_h[1],
            w1 * row1_h[2] + w0 * row2_h[2] + w1 * row3_h[2],
            w1 * row1_h[3] + w0 * row2_h[3] + w1 * row3_h[3],
            w1 * row1_h[4] + w0 * row2_h[4] + w1 * row3_h[4],
        ],
        [
            w0 * row3_h[0] + w1 * row2_h[0] + w2 * row1_h[0],
            w0 * row3_h[1] + w1 * row2_h[1] + w2 * row1_h[1],
            w0 * row3_h[2] + w1 * row2_h[2] + w2 * row1_h[2],
            w0 * row3_h[3] + w1 * row2_h[3] + w2 * row1_h[3],
            w0 * row3_h[4] + w1 * row2_h[4] + w2 * row1_h[4],
        ],
        [
            w1 * row3_h[0] + w2 * row2_h[0],
            w1 * row3_h[1] + w2 * row2_h[1],
            w1 * row3_h[2] + w2 * row2_h[2],
            w1 * row3_h[3] + w2 * row2_h[3],
            w1 * row3_h[4] + w2 * row2_h[4],
        ],
    ];

    // Same expectations, but frozen as literals for quick regression checks.
    let expected_numeric = array![
        [
            0.0786428276229825,
            0.205180120296464,
            0.2769933604483313,
            0.205180120296464,
            0.0786428276229825
        ],
        [
            0.23990695869970506,
            0.6182546706784065,
            0.8291631306895229,
            0.6182546706784065,
            0.23990695869970506
        ],
        [
            0.348701679478469,
            0.8939371765641258,
            1.195497596322409,
            0.8939371765641258,
            0.348701679478469
        ],
        [
            0.23990695869970506,
            0.6182546706784065,
            0.8291631306895229,
            0.6182546706784065,
            0.23990695869970506
        ],
        [
            0.0786428276229825,
            0.205180120296464,
            0.2769933604483313,
            0.205180120296464,
            0.0786428276229825
        ],
    ];

    // Act
    let smoothed = smoothe_counts_gaussian(&counts, sigma, radius);

    // Assert against the literal snapshot (allows small rounding slack).
    assert_matrix_close(&smoothed, &expected_numeric, 1e-3);
    // Assert the explicit hand calculation with tight tolerance to prove the algebra matches.
    assert_matrix_close(&smoothed, &expected_calc, 1e-12);
}

mod collapse_bins_tests {
    use super::*;

    fn simple_counts() -> Array2<f64> {
        array![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    }

    #[test]
    fn collapses_rows_by_sum_without_weights() {
        let counts = simple_counts();
        let bins = build_single_bin_axis(2);
        let collapsed = collapse_counts_by_bins(&counts, 0, &bins, CollapseAggregation::Sum, None)
            .expect("collapse should succeed");
        let expected = array![[5.0, 7.0, 9.0]];
        assert_matrix_close(&collapsed, &expected, 1e-12);
    }

    #[test]
    fn collapses_rows_by_mean_without_weights() {
        let counts = simple_counts();
        let bins = build_single_bin_axis(2);
        let collapsed = collapse_counts_by_bins(&counts, 0, &bins, CollapseAggregation::Mean, None)
            .expect("collapse should succeed");
        let expected = array![[(1.0 + 4.0) / 2.0, (2.0 + 5.0) / 2.0, (3.0 + 6.0) / 2.0]];
        assert_matrix_close(&collapsed, &expected, 1e-12);
    }

    #[test]
    fn collapses_rows_by_weighted_mean() {
        let counts = simple_counts();
        let bins = build_single_bin_axis(2);
        let mass = array![[1.0, 1.0, 1.0], [3.0, 3.0, 3.0]];
        let collapsed = collapse_counts_by_bins(
            &counts,
            0,
            &bins,
            CollapseAggregation::Mean,
            Some(mass.view()),
        )
        .expect("collapse should succeed");
        // Weights per row are mass sums: row0=1, row1=3 => denominator = 1 + 3 = 4 for every column.
        let expected = array![[
            (1.0 * 1.0 + 4.0 * 3.0) / 4.0,
            (2.0 * 1.0 + 5.0 * 3.0) / 4.0,
            (3.0 * 1.0 + 6.0 * 3.0) / 4.0,
        ]];
        assert_matrix_close(&collapsed, &expected, 1e-12);
    }

    #[test]
    fn collapses_columns_by_sum() {
        let counts = simple_counts();
        let bins = build_explicit_bins(&[(0, vec![0, 1]), (1, vec![2])]);
        let collapsed = collapse_counts_by_bins(&counts, 1, &bins, CollapseAggregation::Sum, None)
            .expect("collapse should succeed");
        let expected = array![[1.0 + 2.0, 3.0], [4.0 + 5.0, 6.0]];
        assert_matrix_close(&collapsed, &expected, 1e-12);
    }

    #[test]
    fn collapses_columns_by_mean_without_weights() {
        let counts = simple_counts();
        let bins = build_explicit_bins(&[(0, vec![0, 1]), (1, vec![2])]);
        let collapsed = collapse_counts_by_bins(&counts, 1, &bins, CollapseAggregation::Mean, None)
            .expect("collapse should succeed");
        let expected = array![[(1.0 + 2.0) / 2.0, 3.0], [(4.0 + 5.0) / 2.0, 6.0]];
        assert_matrix_close(&collapsed, &expected, 1e-12);
    }

    #[test]
    fn collapses_columns_by_weighted_mean() {
        let counts = simple_counts();
        let bins = build_explicit_bins(&[(0, vec![0, 1]), (1, vec![2])]);
        let mass = array![[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]];
        let collapsed = collapse_counts_by_bins(
            &counts,
            1,
            &bins,
            CollapseAggregation::Mean,
            Some(mass.view()),
        )
        .expect("collapse should succeed");
        // Column weights are mass sums: col0=2, col1=4, col2=6. Bin 0 covers col0+col1, so denominator = 2 + 4 = 6.
        let expected = array![
            [(1.0 * 2.0 + 2.0 * 4.0) / 6.0, 3.0,],
            [(4.0 * 2.0 + 5.0 * 4.0) / 6.0, 6.0,],
        ];
        assert_matrix_close(&collapsed, &expected, 1e-12);
    }

    #[test]
    fn errors_when_weights_given_for_sum() {
        let counts = simple_counts();
        let bins = build_single_bin_axis(2);
        let mass = array![[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]];
        let result = collapse_counts_by_bins(
            &counts,
            0,
            &bins,
            CollapseAggregation::Sum,
            Some(mass.view()),
        );
        assert!(result.is_err());
    }
}

mod binning_tests {
    use super::*;
    use cfdnalab::commands::gc_bias::binning::bin_greedily_by_mass;

    #[test]
    fn bins_all_mass_into_single_bin_when_threshold_is_high() {
        // Two rows each contribute exactly half the total mass (1 of 2).
        // Setting min_mass_pct to 99% forces all indices into one bin.
        let counts = array![[1.0, 0.0], [0.0, 1.0]];
        let bins = bin_greedily_by_mass(&counts, 0, 99.0, 1).expect("binning should succeed");
        assert_eq!(bins.num_bins, 1, "All rows must merge into a single bin");
        assert_eq!(
            bins.bin_to_indices.get(&0),
            Some(&vec![0, 1]),
            "Both indices should map to bin 0"
        );
    }

    #[test]
    fn splits_bins_when_mass_threshold_is_met() {
        // First three rows contribute exactly 30 units of mass (out of 70 total),
        // so once row 2 is included the threshold is met and a new bin starts.
        // The remaining row becomes the second bin.
        let counts = array![[10.0, 0.0], [0.0, 10.0], [10.0, 0.0], [40.0, 0.0]];
        let bins = bin_greedily_by_mass(&counts, 0, 30.0, 1).expect("binning should succeed");
        assert_eq!(bins.num_bins, 2, "Expected two bins");
        assert_eq!(
            bins.bin_to_indices.get(&0),
            Some(&vec![0, 1, 2]),
            "Indices before crossing the threshold stay in the first bin"
        );
        assert_eq!(
            bins.bin_to_indices.get(&1),
            Some(&vec![3]),
            "The remaining index forms the next bin after the threshold is met"
        );
    }

    #[test]
    fn merges_partial_tail_into_previous_bin_when_threshold_not_met() {
        // First three rows contribute 30 units of mass (out of 71 total ≈42%),
        // so a new bin befins for the fourht element, which is big enough
        // to make it's own bin. BUT the last bin is too small so it gets
        // merged into the second bin instead of being alone.
        let counts = array![
            [10.0, 0.0],
            [0.0, 10.0],
            [10.0, 0.0],
            [40.0, 0.0],
            [10.0, 0.0]
        ];
        let bins = bin_greedily_by_mass(&counts, 0, 30.0, 1).expect("binning should succeed");
        assert_eq!(bins.num_bins, 2, "Expected two bins");
        assert_eq!(
            bins.bin_to_indices.get(&0),
            Some(&vec![0, 1, 2]),
            "Indices before crossing the threshold stay in the first bin"
        );
        assert_eq!(
            bins.bin_to_indices.get(&1),
            Some(&vec![3, 4]),
            "window idx=3 is big enough on its own but idx=4 is too small to be alone so it merged into this window"
        );
    }

    #[test]
    fn handles_zero_total_mass() {
        let counts = array![[0.0, 0.0], [0.0, 0.0]];
        let bins = bin_greedily_by_mass(&counts, 0, 10.0, 1).expect("binning should succeed");
        assert_eq!(bins.num_bins, 0, "No bins should be created for zero mass");
        assert!(bins.index_to_bin.is_empty());
        assert!(bins.bin_to_indices.is_empty());
    }
}
