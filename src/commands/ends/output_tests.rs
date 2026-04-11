use super::*;
use crate::shared::kmers::kmer_codec::{KmerSpec, build_kmer_specs};
use fxhash::FxHashMap;
use std::sync::{Mutex, OnceLock};

fn output_size_env_lock() -> &'static Mutex<()> {
    static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
    LOCK.get_or_init(|| Mutex::new(()))
}

fn spec_for_k(k: u8) -> KmerSpec {
    let specs = build_kmer_specs(&[k]).expect("valid k-mer spec");
    specs[&k].clone()
}

fn expected_combined_1_plus_1_order_without_collapse() -> Vec<String> {
    let bases = ["A", "C", "G", "T"];
    let mut motifs = Vec::new();
    for outside in bases {
        for inside in bases {
            motifs.push(format!("{outside}_{inside}"));
        }
    }
    motifs
}

fn expected_combined_1_plus_1_order_with_collapse() -> Vec<String> {
    let bases = ["A", "C", "G", "T"];
    let mut motifs = Vec::new();
    for outside in ["A", "C"] {
        for inside in bases {
            motifs.push(format!("{outside}_{inside}"));
        }
    }
    motifs
}

fn expected_collapsed_combined_2_plus_2_order() -> Vec<String> {
    let bases = ["A", "C", "G", "T"];
    let mut motifs = Vec::new();
    for first_outside in ["A", "C"] {
        for second_outside in bases {
            for first_inside in bases {
                for second_inside in bases {
                    motifs.push(format!(
                        "{first_outside}{second_outside}_{first_inside}{second_inside}"
                    ));
                }
            }
        }
    }
    motifs
}

#[test]
fn ensure_dense_end_motif_output_size_accepts_small_matrix() {
    // Arrange / Act / Assert
    ensure_dense_end_motif_output_size(10, 20).expect("small dense matrix should be allowed");
}

#[test]
fn ensure_dense_end_motif_output_size_rejects_large_matrix() {
    // Arrange: just over the configured 5 GiB guard.
    let n_values = (DEFAULT_MAX_DENSE_END_MOTIF_OUTPUT_BYTES / 8) + 1;

    // Act
    let err = ensure_dense_end_motif_output_size(1, n_values as usize)
        .expect_err("oversized dense matrix should be rejected");

    // Assert
    assert!(err.to_string().contains("Dense end-motif output would require"));
}

#[test]
fn ensure_dense_end_motif_output_size_respects_env_override() {
    let _guard = output_size_env_lock()
        .lock()
        .expect("env-var test lock should not be poisoned");
    let overridden_limit = DEFAULT_MAX_DENSE_END_MOTIF_OUTPUT_BYTES + 8;

    // Safety: this test serializes environment mutation with a process-local mutex and restores
    // the variable before releasing the lock.
    unsafe {
        std::env::set_var(
            MAX_DENSE_END_MOTIF_OUTPUT_BYTES_ENV,
            overridden_limit.to_string(),
        );
    }
    let result = ensure_dense_end_motif_output_size(1, (overridden_limit / 8) as usize);
    // Safety: same reasoning as above; the mutation is serialized and immediately reverted.
    unsafe {
        std::env::remove_var(MAX_DENSE_END_MOTIF_OUTPUT_BYTES_ENV);
    }

    assert!(
        result.is_ok(),
        "env override should allow a matrix up to the configured byte limit: {result:?}"
    );
}

#[test]
fn ensure_dense_end_motif_output_size_rejects_invalid_env_override() {
    let _guard = output_size_env_lock()
        .lock()
        .expect("env-var test lock should not be poisoned");

    // Safety: this test serializes environment mutation with a process-local mutex and restores
    // the variable before releasing the lock.
    unsafe {
        std::env::set_var(MAX_DENSE_END_MOTIF_OUTPUT_BYTES_ENV, "not-a-number");
    }
    let err = ensure_dense_end_motif_output_size(1, 1)
        .expect_err("invalid env override should fail loudly");
    // Safety: same reasoning as above; the mutation is serialized and immediately reverted.
    unsafe {
        std::env::remove_var(MAX_DENSE_END_MOTIF_OUTPUT_BYTES_ENV);
    }

    assert!(
        err.to_string()
            .contains("CFDNALAB_ENDS_MAX_DENSE_OUTPUT_BYTES must be a positive integer byte count"),
        "unexpected error message: {err:#}"
    );
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
    let inside_spec = spec_for_k(1);

    // Act
    let motifs = build_all_end_motif_order(Some(&inside_spec), None, false).expect("motif order");

    // Assert
    assert_eq!(motifs, vec!["_A", "_C", "_G", "_T"]);
}

#[test]
fn build_all_end_motif_order_collapses_single_base_complements() {
    // Arrange: for k=1, complement and reverse-complement are identical, so this test only
    // checks the single-base universe reduction itself. The transform distinction is covered by
    // the combined-motif tests below.
    let inside_spec = spec_for_k(1);

    // Act
    let motifs =
        build_all_end_motif_order(Some(&inside_spec), None, true).expect("canonical motif order");

    // Assert
    assert_eq!(motifs, vec!["_A", "_C"]);
}

#[test]
fn build_all_end_motif_order_enumerates_the_full_combined_outside_inside_universe() {
    // Arrange
    let inside_spec = spec_for_k(1);
    let outside_spec = spec_for_k(1);

    // Act
    let motifs = build_all_end_motif_order(Some(&inside_spec), Some(&outside_spec), false)
        .expect("combined motif order");

    // Assert: labels are sorted lexicographically after formatting as `<outside>_<inside>`.
    assert_eq!(motifs, expected_combined_1_plus_1_order_without_collapse());
}

#[test]
fn build_all_end_motif_order_collapses_combined_even_length_motifs_by_same_orientation_complement()
{
    // Arrange: with k_outside=1 and k_inside=1, canonicalization is applied to the full
    // `outside || inside` string before formatting.
    //
    // First-principles derivation:
    // - compare the first base of the full motif to its complement
    // - A < T and C < G, so motifs starting with A or C stay as-is
    // - motifs starting with G or T collapse to complements starting with C or A
    //
    // Therefore the full canonical universe is exactly the eight labels whose outside base is
    // A or C.
    let inside_spec = spec_for_k(1);
    let outside_spec = spec_for_k(1);

    // Act
    let motifs = build_all_end_motif_order(Some(&inside_spec), Some(&outside_spec), true)
        .expect("collapsed combined motif order");

    // Assert
    assert_eq!(motifs, expected_combined_1_plus_1_order_with_collapse());
}

#[test]
fn build_all_end_motif_order_collapses_combined_odd_length_motifs_without_swapping_components() {
    // Arrange: k_outside=1 and k_inside=2 gives a 3-base full motif. This is the case that would
    // drift if collapse were done against revcomp(full_motif) instead of complement(full_motif).
    //
    // First-principles derivation:
    // - the full motif is compared as one `outside || inside` string
    // - the first base always decides the lexicographic winner against its complement
    // - canonical full motifs must therefore start with A or C, never G or T
    // - after splitting at k_outside=1, the exact dense universe is:
    //   A_<any 2-mer> and C_<any 2-mer>
    let inside_spec = spec_for_k(2);
    let outside_spec = spec_for_k(1);

    // Act
    let motifs = build_all_end_motif_order(Some(&inside_spec), Some(&outside_spec), true)
        .expect("collapsed odd-length combined motif order");

    // Assert: exact dense universe in sorted order.
    let bases = ["A", "C", "G", "T"];
    let mut expected = Vec::new();
    for outside in ["A", "C"] {
        for first_inside in bases {
            for second_inside in bases {
                expected.push(format!("{outside}_{first_inside}{second_inside}"));
            }
        }
    }
    assert_eq!(motifs, expected);

    // Spot-check the specific pair that distinguishes the intended contract from revcomp-based
    // collapsing on the full decoded motif.
    assert!(motifs.contains(&"C_AT".to_string()));
    assert!(!motifs.contains(&"G_TA".to_string()));
}

#[test]
fn build_all_end_motif_order_collapses_readable_combined_2_plus_2_examples() {
    // Arrange: k_outside=2 and k_inside=2 keeps both halves multi-base, but we only assert a few
    // hand-derived pairs instead of an opaque generated universe.
    //
    // First-principles examples:
    // - "GTAC" complements to "CATG", so the canonical label must be "CA_TG", not "GT_AC"
    // - "TGCA" complements to "ACGT", so the canonical label must be "AC_GT", not "TG_CA"
    // - "ACGT" is already canonical and must remain "AC_GT"
    //
    // These examples are enough to catch:
    // - using revcomp(full_motif) instead of complement(full_motif)
    // - splitting before canonicalization
    // - swapping `outside` and `inside` after collapse
    let inside_spec = spec_for_k(2);
    let outside_spec = spec_for_k(2);

    // Act
    let motifs = build_all_end_motif_order(Some(&inside_spec), Some(&outside_spec), true)
        .expect("collapsed 2+2 dense universe");

    // Assert: exact dense universe plus a few hand-derived labels that make the contract obvious.
    assert_eq!(motifs, expected_collapsed_combined_2_plus_2_order());

    assert!(motifs.contains(&"CA_TG".to_string()));
    assert!(!motifs.contains(&"GT_AC".to_string()));

    assert!(motifs.contains(&"AC_GT".to_string()));
    assert!(!motifs.contains(&"TG_CA".to_string()));

    // `AC_GT` is produced both directly and as the complement of `TG_CA`, so it must appear only
    // once in the dense universe.
    assert_eq!(
        motifs.iter().filter(|motif| motif.as_str() == "AC_GT").count(),
        1
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
