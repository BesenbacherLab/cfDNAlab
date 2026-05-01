use super::*;
use crate::commands::gc_bias::{
    GC_CORRECTION_SCHEMA_VERSION, counting::build_gc_prefixes, package::GCCorrectionPackage,
};
use ndarray::array;

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
