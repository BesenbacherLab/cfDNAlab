use super::*;
use crate::{
    commands::gc_bias::{counting::build_gc_prefixes, package::GCCorrectionPackage},
    shared::constants::GC_CORRECTION_SCHEMA_VERSION,
};
use ndarray::{Array1, Array2, array};

fn one_bin_corrector() -> GCCorrector {
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![30, 40],
        gc_edges: vec![0, 100],
        correction_matrix: array![[2.0_f64]],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_signature: [0, 0],
    };
    GCCorrector::from_package(&package).expect("valid one-bin correction package")
}

fn length_trim_corrector(
    correction_values: Vec<f64>,
    length_bin_frequencies: Vec<f64>,
) -> GCCorrector {
    let num_bins = correction_values.len();
    let mut length_edges = Vec::with_capacity(num_bins + 1);
    for index in 0..num_bins {
        length_edges.push(10 + index as u32 * 10);
    }
    length_edges.push(10 + num_bins as u32 * 10);

    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges,
        gc_edges: vec![0, 100],
        correction_matrix: Array2::from_shape_vec((num_bins, 1), correction_values)
            .expect("test correction matrix shape should be valid"),
        length_bin_frequencies: Array1::from(length_bin_frequencies),
        reference_contig_signature: [0, 0],
    };
    GCCorrector::from_package(&package).expect("valid length-trim correction package")
}

fn collapsed_gc_value(
    corrector: &GCCorrector,
    gc_length_range: GCLengthRange,
    min_fragment_length: u32,
    max_fragment_length: u32,
    trim_rare: f64,
) -> f64 {
    collapsed_gc_value_with_weighting(
        corrector,
        &MarginalizeLengthsWeightingScheme::Equal,
        gc_length_range,
        min_fragment_length,
        max_fragment_length,
        trim_rare,
    )
}

fn collapsed_gc_value_with_weighting(
    corrector: &GCCorrector,
    weighting_scheme: &MarginalizeLengthsWeightingScheme,
    gc_length_range: GCLengthRange,
    min_fragment_length: u32,
    max_fragment_length: u32,
    trim_rare: f64,
) -> f64 {
    LengthAgnosticGCCorrector::from_gc_corrector(
        corrector,
        weighting_scheme,
        gc_length_range,
        trim_rare,
        min_fragment_length,
        max_fragment_length,
    )
    .expect("length-agnostic corrector should build")
    .correction_vector[0]
}

#[test]
fn length_agnostic_corrector_applies_trim_before_frequency_based_weighting() {
    // Arrange:
    // An 11% trim has budget 11 for frequencies [70, 20, 10].
    // The 10-frequency row is trimmed before weighting. The retained rows have
    // corrections [1, 10] and frequencies [70, 20].
    let corrector = length_trim_corrector(vec![1.0, 10.0, 100.0], vec![70.0, 20.0, 10.0]);

    // Act
    let frequency_weighted = collapsed_gc_value_with_weighting(
        &corrector,
        &MarginalizeLengthsWeightingScheme::Frequency,
        GCLengthRange::Package,
        10,
        40,
        0.11,
    );
    let max_frequency = collapsed_gc_value_with_weighting(
        &corrector,
        &MarginalizeLengthsWeightingScheme::MaxFrequency,
        GCLengthRange::Package,
        10,
        40,
        0.11,
    );

    // Assert
    // Frequency weighting after trimming gives (1 * 70 + 10 * 20) / 90 = 3.
    // Max-frequency weighting after trimming picks the retained 70-frequency row.
    assert!((frequency_weighted - 3.0).abs() < 1e-12);
    assert!((max_frequency - 1.0).abs() < 1e-12);
}

#[test]
fn length_agnostic_corrector_trims_rarest_bins_before_equal_weighting() {
    // Arrange:
    // Frequencies are [70, 20, 10], so an 11% trim has budget 11.
    // The row with frequency 10 is removed. The next rarest row has frequency
    // 20 and would exceed the budget, so rows with corrections 2 and 8 remain.
    let corrector = length_trim_corrector(vec![2.0, 8.0, 100.0], vec![70.0, 20.0, 10.0]);

    // Act
    let observed = collapsed_gc_value(&corrector, GCLengthRange::Package, 10, 40, 0.11);

    // Assert
    // Equal weighting after trimming gives (2 + 8) / 2 = 5.
    assert!((observed - 5.0).abs() < 1e-12);
}

#[test]
fn length_agnostic_corrector_keeps_rarest_bin_when_it_exceeds_trim_budget() {
    // Arrange:
    // Frequencies are [60, 40], so a 5% trim has budget 5.
    // Removing the rarest row would remove 40%, so both rows must remain.
    let corrector = length_trim_corrector(vec![2.0, 8.0], vec![60.0, 40.0]);

    // Act
    let observed = collapsed_gc_value(&corrector, GCLengthRange::Package, 10, 30, 0.05);

    // Assert
    // Equal weighting with both rows gives (2 + 8) / 2 = 5.
    assert!((observed - 5.0).abs() < 1e-12);
}

#[test]
fn length_agnostic_corrector_keeps_nearly_equal_frequency_group_when_only_part_fits_budget() {
    // Arrange:
    // Frequencies are [0.8, 0.1, 0.10000005], so a 10% trim has a budget just
    // above 0.1. The two rarest rows differ only below the tie tolerance and
    // together exceed the budget. Removing only the shorter near-tied row would
    // fit the budget but would introduce a length-order tie bias, so the whole
    // near-tied group is retained.
    let corrector = length_trim_corrector(vec![2.0, 8.0, 100.0], vec![0.8, 0.1, 0.10000005]);

    // Act
    let observed = collapsed_gc_value(&corrector, GCLengthRange::Package, 10, 40, 0.10);

    // Assert
    // Equal weighting with all rows retained gives (2 + 8 + 100) / 3.
    assert!((observed - (110.0 / 3.0)).abs() < 1e-12);
}

#[test]
fn length_agnostic_corrector_splits_frequency_group_when_difference_exceeds_tolerance() {
    // Arrange:
    // Frequencies are [0.8, 0.1, 0.1000002]. The selected total is 1.0000002,
    // so the tie tolerance is about 1.0000002e-7. The two rarest rows differ
    // by 2e-7, which is above tolerance. A 10% trim budget is 0.10000002, so
    // only the row with frequency 0.1 is removed.
    let corrector = length_trim_corrector(vec![2.0, 8.0, 100.0], vec![0.8, 0.1, 0.1000002]);

    // Act
    let observed = collapsed_gc_value(&corrector, GCLengthRange::Package, 10, 40, 0.10);

    // Assert
    // Equal weighting over retained rows gives (2 + 100) / 2 = 51.
    assert!((observed - 51.0).abs() < 1e-12);
}

#[test]
fn length_agnostic_corrector_trims_after_requested_range_selection() {
    // Arrange:
    // Package rows have corrections [2, 100, 10, 8] and frequencies [70, 4, 6, 20].
    // The requested range [20,39] selects only rows 1 and 2, with frequencies [4, 6].
    // A 50% trim has budget 5 within the selected rows, so row 1 is removed and
    // row 2 remains. Trimming before requested-range selection would produce a
    // different retained row, so this test fixes the operation order.
    let corrector =
        length_trim_corrector(vec![2.0, 100.0, 10.0, 8.0], vec![70.0, 4.0, 6.0, 20.0]);

    // Act
    let requested_observed =
        collapsed_gc_value(&corrector, GCLengthRange::Requested, 20, 39, 0.50);

    // Assert
    assert!((requested_observed - 10.0).abs() < 1e-12);

    // Act
    let package_observed = collapsed_gc_value(&corrector, GCLengthRange::Package, 20, 39, 0.50);

    // Assert:
    // With `package`, all rows are selected before trimming. The total frequency is
    // 70 + 4 + 6 + 20 = 100, so a 50% trim has budget 50. Rows with frequencies
    // 4, 6, and 20 are removed, leaving only the correction 2 row.
    assert!((package_observed - 2.0).abs() < 1e-12);
}

#[test]
fn length_agnostic_corrector_errors_when_trim_requested_with_zero_selected_frequency() {
    // Arrange
    let corrector = length_trim_corrector(vec![2.0, 8.0], vec![0.0, 0.0]);

    // Act
    let error = LengthAgnosticGCCorrector::from_gc_corrector(
        &corrector,
        &MarginalizeLengthsWeightingScheme::Equal,
        GCLengthRange::Package,
        0.05,
        10,
        30,
    )
    .expect_err("positive rare-bin trim should require nonzero selected frequencies");

    // Assert
    assert!(
        error
            .to_string()
            .contains("Cannot trim rare GC length bins because selected length-bin frequencies sum to zero"),
        "unexpected error: {error}"
    );
}

#[test]
fn correct_fragment_returns_none_when_fragment_length_is_below_package_range() {
    // Arrange: the package covers aligned lengths 30..=40, but the fragment's aligned reference
    // span is only 28 bp. This can happen in `ends --clip-strategy raw-shifted-boundary`, where
    // assignment length includes soft clips but GC correction still uses aligned reference bases.
    let corrector = one_bin_corrector();
    let gc_prefixes = build_gc_prefixes(&vec![b'G'; 28]);
    let fragment_interval = Interval::new(0_u64, 28_u64).expect("valid fragment interval");

    // Act
    let weight = corrector
        .correct_fragment(fragment_interval, &gc_prefixes)
        .expect("out-of-package length should be handled as an unusable GC weight");

    // Assert
    assert_eq!(weight, None);
}

#[test]
fn get_correction_weight_errors_without_underflow_when_fragment_length_is_below_package_range() {
    // Arrange: direct weight lookup is still a strict API, but it must fail with a useful error
    // rather than subtracting below the package minimum length.
    let corrector = one_bin_corrector();

    // Act
    let error = corrector
        .get_correction_weight(29, 50)
        .expect_err("length below package range should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("GC correction: unexpected fragment length 29")
    );
}
