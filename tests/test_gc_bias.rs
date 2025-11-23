use fxhash::FxHashMap;
use ndarray::array;

use cfdnalab::commands::gc_bias::{
    binning::{BinnedAxis, bins_from_edges, compute_bin_edges},
    support_masking::build_extreme_gc_support_mask,
};

#[test]
fn masks_extreme_gc_bins_per_side_in_square_matrix() {
    // Arrange: 6x6 matrix with two extreme GC bins on each side.
    let expected = array![
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
        [false, false, true, true, false, false],
    ];

    // Act: build the support mask after binning.
    let mask = build_extreme_gc_support_mask((6, 6), 2);

    // Assert: the central two GC bins remain supported across all lengths.
    assert_eq!(mask, expected);
}

#[test]
fn round_trips_bins_to_edges_and_back() {
    // Arrange: build a simple BinnedAxis where bins group indices as [0-1], [2-4], and [5-7].
    let mut index_to_bin = FxHashMap::default();
    let mut bin_to_indices = FxHashMap::default();
    let bins: [Vec<usize>; 3] = [vec![0, 1], vec![2, 3, 4], vec![5, 6, 7]];
    for (bin_idx, indices) in bins.iter().enumerate() {
        bin_to_indices.insert(bin_idx, indices.clone());
        for &idx in indices {
            index_to_bin.insert(idx, bin_idx);
        }
    }
    let axis = BinnedAxis {
        index_to_bin,
        bin_to_indices,
        num_bins: 3,
    };

    // Act: compute edges then reconstruct the bins.
    let edges = compute_bin_edges(&axis, 0, 7).expect("edges should be computed");
    let reconstructed_axis = bins_from_edges(edges.as_slice()).expect("rebuild should work");

    // Assert: the derived edges match the expected bin boundaries, and the reconstructed
    // axis matches the original bin layout.
    assert_eq!(edges, vec![0, 2, 5, 7]);
    assert_eq!(reconstructed_axis.num_bins, axis.num_bins);
    assert_eq!(reconstructed_axis.bin_to_indices, axis.bin_to_indices);
    assert_eq!(reconstructed_axis.index_to_bin, axis.index_to_bin);
}
