use ndarray::array;

use cfdnalab::commands::gc_bias::gc_bias::build_extreme_gc_support_mask;

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
