use anyhow::Result;
use fxhash::FxHashMap;
use ndarray::array;
use tempfile::tempdir;

use cfdnalab::commands::gc_bias::{
    GC_CORRECTION_SCHEMA_VERSION,
    binning::{BinnedAxis, bins_from_edges, compute_bin_edges},
    correct::{GCCorrector, LengthAgnosticGCCorrector, MarginalizeLengthsWeightingScheme},
    gc_bias::interpolate_masked_corrections,
    package::GCCorrectionPackage,
    support_masking::build_extreme_bins_support_mask,
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
    let mask = build_extreme_bins_support_mask((6, 6), 2, 0);

    // Assert: the central two GC bins remain supported across all lengths.
    assert_eq!(mask, expected);
}

#[test]
fn masks_shortest_length_bins_in_matrix() {
    // Arrange: 5x4 matrix with one shortest length bin masked.
    let expected = array![
        [false, false, false, false],
        [true, true, true, true],
        [true, true, true, true],
        [true, true, true, true],
        [true, true, true, true],
    ];

    // Act: build the support mask after binning.
    let mask = build_extreme_bins_support_mask((5, 4), 0, 1);

    // Assert: the central three length bins remain supported across all GC bins.
    assert_eq!(mask, expected);
}

#[test]
fn interpolates_masked_short_length_row() -> Result<()> {
    // Arrange: first length row is masked; other rows are supported.
    let mut matrix = array![
        [0.0_f64, 0.0_f64],
        [2.0_f64, 2.0_f64],
        [4.0_f64, 4.0_f64],
        [6.0_f64, 6.0_f64],
    ];
    let mask = build_extreme_bins_support_mask((4, 2), 0, 1);

    // Act: interpolate masked bins.
    interpolate_masked_corrections(&mut matrix, &mask)?;

    // Assert: the masked first row is filled using neighbouring lengths.
    assert!((matrix[(0, 0)] - 2.0).abs() < 1e-6);
    assert!((matrix[(0, 1)] - 2.0).abs() < 1e-6);
    Ok(())
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

#[test]
fn provides_expected_weights_after_roundtrip() -> Result<()> {
    // Arrange: a package whose edges start at non-zero values so offset logic is exercised.
    let length_edges = vec![30, 34, 40];
    let gc_edges = vec![10, 60, 90];
    let correction_matrix = array![[1.5_f64, 2.0_f64], [0.5_f64, 0.75_f64]];
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 3,
        length_edges: length_edges.clone(),
        gc_edges: gc_edges.clone(),
        correction_matrix,
        length_bin_frequencies: array![1.0_f64, 1.0_f64],
    };
    let tmp_dir = tempdir()?;
    let pkg_path = tmp_dir.path().join("gc_package.npz");
    package.write_npz(&pkg_path)?;

    // Act: load the package and build a corrector.
    let loaded = GCCorrectionPackage::from_file(&pkg_path)?;
    let corrector = GCCorrector::from_package(&loaded)?;

    // Assert: fragments landing in each bin retrieve the expected weights.
    let weight_len31_gc20 = corrector.get_correction_weight(31, 20)?;
    assert!(
        (weight_len31_gc20 - 1.5).abs() < f64::EPSILON,
        "length 31 / GC 20 should map to 1.5"
    );

    let weight_len32_gc70 = corrector.get_correction_weight(32, 70)?;
    assert!(
        (weight_len32_gc70 - 2.0).abs() < f64::EPSILON,
        "length 32 / GC 70 should map to 2.0"
    );

    let weight_len38_gc55 = corrector.get_correction_weight(38, 55)?;
    assert!(
        (weight_len38_gc55 - 0.5).abs() < f64::EPSILON,
        "length 38 / GC 55 should map to 0.5"
    );

    let weight_len39_gc80 = corrector.get_correction_weight(39, 80)?;
    assert!(
        (weight_len39_gc80 - 0.75).abs() < f64::EPSILON,
        "length 39 / GC 80 should map to 0.75"
    );

    Ok(())
}

fn make_length_agnostic_package() -> GCCorrectionPackage {
    let correction_matrix = array![[1.0_f64, 2.0_f64], [3.0_f64, 5.0_f64]];
    GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![20, 30, 40],
        gc_edges: vec![0, 50, 100],
        correction_matrix,
        length_bin_frequencies: array![0.2_f64, 0.8_f64],
    }
}

#[test]
fn length_agnostic_equal_weighting_means_rows() -> Result<()> {
    let package = make_length_agnostic_package();
    let corrector = GCCorrector::from_package(&package)?;
    let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
        &corrector,
        &MarginalizeLengthsWeightingScheme::Equal,
    )?;

    assert!((agnostic.get_correction_weight(0)? - 2.0).abs() < 1e-12);
    assert!((agnostic.get_correction_weight(50)? - 3.5).abs() < 1e-12);
    Ok(())
}

#[test]
fn length_agnostic_coverage_weighting_uses_frequencies() -> Result<()> {
    let package = make_length_agnostic_package();
    let corrector = GCCorrector::from_package(&package)?;
    let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
        &corrector,
        &MarginalizeLengthsWeightingScheme::Coverage,
    )?;

    // Weighted average with frequencies [2, 8]
    assert!((agnostic.get_correction_weight(0)? - 2.6).abs() < 1e-12);
    assert!((agnostic.get_correction_weight(50)? - 4.4).abs() < 1e-12);
    Ok(())
}

#[test]
fn length_agnostic_max_coverage_picks_most_frequent_row() -> Result<()> {
    let package = make_length_agnostic_package();
    let corrector = GCCorrector::from_package(&package)?;
    let agnostic = LengthAgnosticGCCorrector::from_gc_corrector(
        &corrector,
        &MarginalizeLengthsWeightingScheme::MaxCoverage,
    )?;

    // Row with highest frequency is [3.0, 5.0]
    assert!((agnostic.get_correction_weight(0)? - 3.0).abs() < 1e-12);
    assert!((agnostic.get_correction_weight(50)? - 5.0).abs() < 1e-12);
    Ok(())
}
