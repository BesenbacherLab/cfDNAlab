use super::*;
use crate::shared::kmers::kmer_codec::build_kmer_specs;
use anyhow::{bail, ensure, Result};

fn motif_labels(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| value.to_string()).collect()
}

fn store_positive_weight(weight: f64) -> Result<bool> {
    ensure!(weight >= 0.0, "negative weight {weight}");
    Ok(weight > 0.0)
}

fn allow_dense_output_size(_n_windows: usize, _n_motifs: usize) -> Result<()> {
    Ok(())
}

fn reject_dense_output_size(_n_windows: usize, _n_motifs: usize) -> Result<()> {
    bail!("dense output guard was called")
}

fn legacy_all_motifs(k: usize, spec: &KmerSpec) -> Vec<String> {
    let max_code = 5u64.pow(k as u32) - 1;
    (0..=max_code)
        .map(|code| spec.decode_kmer(code))
        .filter(|motif| !motif.contains('N'))
        .collect()
}

#[test]
fn all_motifs_returns_the_full_acgt_space() -> Result<()> {
    // For k = 2 there are exactly 4^2 = 16 A/C/G/T motifs. The helper must not include any motif
    // containing N because those are sentinel states, not countable reference k-mers.
    let specs = build_kmer_specs(&[2])?;

    let motifs = all_motifs(2, &specs);

    assert_eq!(motifs.len(), 16);
    assert!(motifs.contains(&"AA".to_string()));
    assert!(motifs.contains(&"AC".to_string()));
    assert!(motifs.contains(&"TT".to_string()));
    assert!(!motifs.iter().any(|motif| motif.contains('N')));
    Ok(())
}

#[test]
fn all_motifs_preserves_legacy_radix5_order_for_k3() -> Result<()> {
    // Arrange: this is the exact A/C/G/T order produced by filtering the radix-5 A/C/G/T/N
    // universe. The optimized radix-4 path must preserve both the labels and their order.
    let specs = build_kmer_specs(&[3])?;
    let spec = &specs[&3];
    let expected = vec![
        "AAA".to_string(),
        "AAC".to_string(),
        "AAG".to_string(),
        "AAT".to_string(),
        "ACA".to_string(),
        "ACC".to_string(),
        "ACG".to_string(),
        "ACT".to_string(),
        "AGA".to_string(),
        "AGC".to_string(),
        "AGG".to_string(),
        "AGT".to_string(),
        "ATA".to_string(),
        "ATC".to_string(),
        "ATG".to_string(),
        "ATT".to_string(),
        "CAA".to_string(),
        "CAC".to_string(),
        "CAG".to_string(),
        "CAT".to_string(),
        "CCA".to_string(),
        "CCC".to_string(),
        "CCG".to_string(),
        "CCT".to_string(),
        "CGA".to_string(),
        "CGC".to_string(),
        "CGG".to_string(),
        "CGT".to_string(),
        "CTA".to_string(),
        "CTC".to_string(),
        "CTG".to_string(),
        "CTT".to_string(),
        "GAA".to_string(),
        "GAC".to_string(),
        "GAG".to_string(),
        "GAT".to_string(),
        "GCA".to_string(),
        "GCC".to_string(),
        "GCG".to_string(),
        "GCT".to_string(),
        "GGA".to_string(),
        "GGC".to_string(),
        "GGG".to_string(),
        "GGT".to_string(),
        "GTA".to_string(),
        "GTC".to_string(),
        "GTG".to_string(),
        "GTT".to_string(),
        "TAA".to_string(),
        "TAC".to_string(),
        "TAG".to_string(),
        "TAT".to_string(),
        "TCA".to_string(),
        "TCC".to_string(),
        "TCG".to_string(),
        "TCT".to_string(),
        "TGA".to_string(),
        "TGC".to_string(),
        "TGG".to_string(),
        "TGT".to_string(),
        "TTA".to_string(),
        "TTC".to_string(),
        "TTG".to_string(),
        "TTT".to_string(),
    ];

    // Act
    let radix4_motifs = all_motifs(3, &specs);
    let legacy_motifs = legacy_all_motifs(3, spec);

    // Assert
    assert_eq!(radix4_motifs, expected);
    assert_eq!(legacy_motifs, expected);
    assert_eq!(radix4_motifs, legacy_motifs);
    Ok(())
}

#[test]
fn compacts_selected_motif_counts_in_label_order() -> Result<()> {
    // Arrange: target 2 is observed before target 0, but the output axis should still follow
    // parser-assigned label order.
    let labels = motif_labels(&["first", "second", "third"]);
    let counts_by_window = vec![(1_u64, vec![(2_u32, 2.0)]), (0, vec![(0, 1.5)])];

    // Act
    let (bins, motif_order) = postprocess_selected_motif_counts(
        counts_by_window,
        3,
        &labels,
        false,
        store_positive_weight,
        allow_dense_output_size,
    )?;

    // Assert
    assert_eq!(motif_order, motif_labels(&["first", "third"]));
    assert_eq!(bins.len(), 3);
    assert_eq!(bins[0].get(&0), Some(&1.5));
    assert_eq!(bins[1].get(&1), Some(&2.0));
    assert!(bins[2].is_empty());

    Ok(())
}

#[test]
fn retains_unobserved_selected_motif_targets_when_requested() -> Result<()> {
    // Arrange: only the middle target has counts, but include-all mode keeps the full target axis.
    let labels = motif_labels(&["first", "second", "third"]);
    let counts_by_window = vec![(0_u64, vec![(1_u32, 2.0)])];

    // Act
    let (bins, motif_order) = postprocess_selected_motif_counts(
        counts_by_window,
        2,
        &labels,
        true,
        store_positive_weight,
        allow_dense_output_size,
    )?;

    // Assert
    assert_eq!(motif_order, labels);
    assert_eq!(bins[0].get(&1), Some(&2.0));
    assert!(!bins[0].contains_key(&0));
    assert!(!bins[0].contains_key(&2));
    assert!(bins[1].is_empty());

    Ok(())
}

#[test]
fn filters_weights_rejected_by_selected_motif_weight_checker() -> Result<()> {
    // Arrange: zero weight is valid but should not create an observed output column.
    let labels = motif_labels(&["first", "second"]);
    let counts_by_window = vec![(0_u64, vec![(0_u32, 0.0), (1, 2.0)])];

    // Act
    let (bins, motif_order) = postprocess_selected_motif_counts(
        counts_by_window,
        1,
        &labels,
        false,
        store_positive_weight,
        allow_dense_output_size,
    )?;

    // Assert
    assert_eq!(motif_order, motif_labels(&["second"]));
    assert_eq!(bins[0].get(&0), Some(&2.0));
    assert!(!bins[0].contains_key(&1));

    Ok(())
}

#[test]
fn rejects_out_of_bounds_selected_motif_indices() {
    // Arrange
    let labels = motif_labels(&["first"]);
    let bad_window = vec![(1_u64, vec![(0_u32, 1.0)])];
    let bad_target = vec![(0_u64, vec![(1_u32, 1.0)])];

    // Act
    let window_error = postprocess_selected_motif_counts(
        bad_window,
        1,
        &labels,
        false,
        store_positive_weight,
        allow_dense_output_size,
    )
    .expect_err("out-of-bounds window should fail");
    let target_error = postprocess_selected_motif_counts(
        bad_target,
        1,
        &labels,
        false,
        store_positive_weight,
        allow_dense_output_size,
    )
    .expect_err("out-of-bounds target should fail");

    // Assert
    assert!(
        window_error.to_string().contains("window index 1"),
        "unexpected error: {window_error:#}"
    );
    assert!(
        target_error.to_string().contains("target index 1"),
        "unexpected error: {target_error:#}"
    );
}

#[test]
fn applies_dense_guard_only_when_all_selected_motif_targets_are_requested() {
    // Arrange
    let labels = motif_labels(&["first"]);
    let counts_by_window = vec![(0_u64, vec![(0_u32, 1.0)])];

    // Act
    postprocess_selected_motif_counts(
        counts_by_window.clone(),
        1,
        &labels,
        false,
        store_positive_weight,
        reject_dense_output_size,
    )
    .expect("sparse selected output should not call dense guard");
    let error = postprocess_selected_motif_counts(
        counts_by_window,
        1,
        &labels,
        true,
        store_positive_weight,
        reject_dense_output_size,
    )
    .expect_err("include-all selected output should call dense guard");

    // Assert
    assert!(
        error.to_string().contains("dense output guard"),
        "unexpected error: {error:#}"
    );
}
