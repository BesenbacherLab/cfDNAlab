use super::*;
use fxhash::FxHashMap;
use crate::shared::kmers::kmer_codec::{KmerSpec, build_kmer_specs};

fn spec_for_k(k: u8) -> KmerSpec {
    let specs = build_kmer_specs(&[k]).expect("valid k-mer spec");
    specs[&k].clone()
}

#[test]
fn ensure_dense_end_motif_output_size_accepts_small_matrix() {
    // Arrange / Act / Assert
    ensure_dense_end_motif_output_size(10, 20).expect("small dense matrix should be allowed");
}

#[test]
fn ensure_dense_end_motif_output_size_rejects_large_matrix() {
    // Arrange: just over the configured 5 GiB guard.
    let n_values = (MAX_DENSE_END_MOTIF_OUTPUT_BYTES / 8) + 1;

    // Act
    let err = ensure_dense_end_motif_output_size(1, n_values as usize)
        .expect_err("oversized dense matrix should be rejected");

    // Assert
    assert!(err.to_string().contains("Dense end-motif output would require"));
}

#[test]
fn ensure_all_motifs_enumeration_size_accepts_small_universe() {
    // Arrange / Act / Assert: 4^(1+1) = 16 motifs, which is trivially small.
    ensure_all_motifs_enumeration_size(1, 1, 10)
        .expect("small all-motifs universe should be allowed");
}

#[test]
fn ensure_all_motifs_enumeration_size_rejects_large_universe_before_enumeration() {
    // Arrange: 4^20 motifs is already far beyond the dense output guard for one window.
    let err = ensure_all_motifs_enumeration_size(10, 10, 1)
        .expect_err("large all-motifs universe should be rejected");

    // Assert
    assert!(err.to_string().contains("refusing to enumerate all motifs"));
}

#[test]
fn build_all_end_motif_order_returns_full_single_base_universe_without_collapse() {
    // Arrange
    let within_spec = spec_for_k(1);

    // Act
    let motifs = build_all_end_motif_order(Some(&within_spec), None, false).expect("motif order");

    // Assert
    assert_eq!(motifs, vec!["_A", "_C", "_G", "_T"]);
}

#[test]
fn build_all_end_motif_order_collapses_single_base_reverse_complements() {
    // Arrange
    let within_spec = spec_for_k(1);

    // Act
    let motifs =
        build_all_end_motif_order(Some(&within_spec), None, true).expect("canonical motif order");

    // Assert
    assert_eq!(motifs, vec!["_A", "_C"]);
}

#[test]
fn build_all_end_motif_order_enumerates_the_full_combined_outside_within_universe() {
    // Arrange
    let within_spec = spec_for_k(1);
    let outside_spec = spec_for_k(1);

    // Act
    let motifs = build_all_end_motif_order(Some(&within_spec), Some(&outside_spec), false)
        .expect("combined motif order");

    // Assert: labels are sorted lexicographically after formatting as `<outside>_<within>`.
    assert_eq!(
        motifs,
        vec![
            "A_A", "A_C", "A_G", "A_T", "C_A", "C_C", "C_G", "C_T", "G_A", "G_C", "G_G",
            "G_T", "T_A", "T_C", "T_G", "T_T",
        ]
    );
}

#[test]
fn collect_end_motif_order_returns_sorted_union_of_observed_sparse_motifs() {
    // Arrange: the sparse path should keep only observed motifs, but it must still produce one
    // stable global order across all windows.
    //
    // Mental derivation:
    // - observed union = {"_G", "AT_C", "_A"}
    // - the implementation uses a `BTreeSet`, so the final order is plain lexicographic order
    let bins = vec![
        FxHashMap::from_iter([
            ("_G".to_string(), 1.0),
            ("AT_C".to_string(), 2.0),
        ]),
        FxHashMap::from_iter([
            ("_A".to_string(), 3.0),
            ("AT_C".to_string(), 4.0),
        ]),
    ];

    // Act
    let motifs = collect_end_motif_order(&bins);

    // Assert: `collect_end_motif_order` uses a `BTreeSet`, so labels are sorted
    // lexicographically in plain string order.
    assert_eq!(motifs, vec!["AT_C", "_A", "_G"]);
}
