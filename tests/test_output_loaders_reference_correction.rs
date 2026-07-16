#![cfg(all(feature = "cmd_ends", feature = "cmd_ref_kmers"))]
//! Public API tests for Rust reference correction in end-motif output loaders.

use cfdnalab::output_loaders::{
    EndMotifSparseEntry, EndMotifStorageMode, TwoSidedCorrectionMode, UnsupportedReferencePolicy,
    load_ends_output, load_ref_kmers_output,
};
use serde_json::{Map, Value, json};
use std::{fs, path::Path, sync::Arc};
use tempfile::TempDir;
use zarrs::{
    array::{ArrayBuilder, DataType, Element, builder::ArrayBuilderFillValue, data_type},
    filesystem::FilesystemStore,
    group::GroupBuilder,
};

#[derive(Debug, Clone, Copy)]
enum FixtureStorage {
    Dense,
    Sparse,
}

/// Verify the correction scale is stable across dense and sparse storage modes.
#[test]
fn corrected_counts_keep_reference_row_scale_across_storage_modes() -> anyhow::Result<()> {
    for end_storage in [FixtureStorage::Dense, FixtureStorage::Sparse] {
        for ref_storage in [FixtureStorage::Dense, FixtureStorage::Sparse] {
            let temp = TempDir::new()?;
            let ends_path = temp.path().join(format!("ends_{end_storage:?}.zarr"));
            let ref_path = temp.path().join(format!("ref_{ref_storage:?}.zarr"));

            write_windowed_end_motif_store(
                &ends_path,
                end_storage,
                &["_AA", "_CC", "_GG"],
                &[1.0, 0.0, 2.5, 0.5, 4.0, 0.0],
            )?;
            write_windowed_ref_kmer_store(
                &ref_path,
                ref_storage,
                &["AA", "CC", "GG"],
                &[1.0 / 3.0, 1.0 / 6.0, 0.5, 0.5, 0.25, 0.25],
            )?;

            let ends = load_ends_output(&ends_path)?;
            let ref_kmers = load_ref_kmers_output(&ref_path)?;

            let full = ends.select_corrected_counts(&ref_kmers).read()?;
            let selected = ends
                .select_corrected_counts(&ref_kmers)
                .windows(&[1, 0])
                .motifs_by_label(&["_CC", "_AA"])
                .read()?;

            assert_eq!(selected.row_indices(), &[1, 0]);
            assert_eq!(
                selected.motif_labels(),
                &["_CC".to_string(), "_AA".to_string()]
            );
            assert_eq!(
                selected.storage_mode(),
                match end_storage {
                    FixtureStorage::Dense => EndMotifStorageMode::Dense,
                    FixtureStorage::Sparse => EndMotifStorageMode::SparseCoo,
                }
            );
            assert_close(selected.count(0, 0).unwrap(), 16.0 / 3.0);
            assert_close(selected.count(0, 1).unwrap(), 1.0 / 3.0);
            assert_close(selected.count(1, 0).unwrap(), 0.0);
            assert_close(selected.count(1, 1).unwrap(), 1.0);

            // Motif selection happens after defining the reference scale. The
            // selected value must equal the corresponding value from the full
            // corrected matrix.
            assert_close(selected.count(0, 0).unwrap(), full.count(1, 1).unwrap());
        }
    }
    Ok(())
}

/// Verify sparse reference support counts use positive stored frequencies.
#[test]
fn sparse_corrected_counts_use_positive_reference_support_per_row() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("reference.ref_kmers.zarr");

    write_windowed_end_motif_store(
        &ends_path,
        FixtureStorage::Sparse,
        &["_AA", "_CC", "_GG"],
        &[1.0, 0.0, 2.5, 0.0, 4.0, 0.0],
    )?;
    write_windowed_ref_kmer_sparse_store(
        &ref_path,
        &["AA", "CC", "GG"],
        &[0, 0, 1, 1],
        &[0, 2, 0, 1],
        &[0.5, 0.5, 0.5, 0.5],
    )?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    let corrected = ends.select_corrected_counts(&ref_kmers).read()?;

    assert_eq!(corrected.storage_mode(), EndMotifStorageMode::SparseCoo);
    assert_eq!(
        corrected.sparse_counts()?.entries().collect::<Vec<_>>(),
        vec![
            EndMotifSparseEntry {
                row_index: 0,
                motif_index: 0,
                count: 1.0,
            },
            EndMotifSparseEntry {
                row_index: 0,
                motif_index: 2,
                count: 2.5,
            },
            EndMotifSparseEntry {
                row_index: 1,
                motif_index: 1,
                count: 4.0,
            },
        ]
    );
    Ok(())
}

/// Verify grouped selectors match reference rows by group name, not row index.
#[test]
fn corrected_counts_map_selected_group_rows_by_key() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("reference.ref_kmers.zarr");

    let ends_store = create_store(
        &ends_path,
        end_root_attributes(FixtureStorage::Dense, "grouped_bed"),
    )?;
    write_motif_axis(&ends_store, &["_A", "_C"])?;
    write_group_rows(&ends_store, &["alpha", "beta"])?;
    write_end_counts(
        &ends_store,
        FixtureStorage::Dense,
        2,
        2,
        &[2.0, 0.0, 0.0, 4.0],
    )?;

    let ref_store = create_store(
        &ref_path,
        ref_root_attributes(FixtureStorage::Dense, "grouped_bed", &["A", "C"]),
    )?;
    write_motif_axis(&ref_store, &["A", "C"])?;
    write_group_rows(&ref_store, &["beta", "alpha"])?;
    write_row_scaling_factors(&ref_store, 2)?;
    write_reference_footprint(&ref_store)?;
    write_ref_frequencies(
        &ref_store,
        FixtureStorage::Dense,
        2,
        2,
        &[0.5, 0.5, 0.25, 0.75],
    )?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;

    let corrected = ends
        .select_corrected_counts(&ref_kmers)
        .groups(&[1])
        .motifs_by_label(&["_C"])
        .read()?;

    assert_eq!(corrected.row_indices(), &[1]);
    assert_eq!(corrected.motif_labels(), &["_C".to_string()]);
    assert_close(corrected.count(0, 0).unwrap(), 4.0);
    Ok(())
}

/// Verify global reference correction is explicit for non-global end motifs.
#[test]
fn corrected_counts_require_global_bias_opt_in() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("global.ref_kmers.zarr");

    write_windowed_end_motif_store(
        &ends_path,
        FixtureStorage::Dense,
        &["_AA", "_CC", "_GG"],
        &[1.0, 0.0, 2.5, 0.5, 4.0, 0.0],
    )?;
    write_global_ref_kmer_store(
        &ref_path,
        FixtureStorage::Dense,
        &["AA", "CC", "GG"],
        &[1.0 / 3.0, 1.0 / 6.0, 0.5],
    )?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    let error = ends
        .select_corrected_counts(&ref_kmers)
        .read()
        .expect_err("global reference correction should require opt-in");
    assert!(error.to_string().contains("use_global_bias(true)"));

    let corrected = ends
        .select_corrected_counts(&ref_kmers)
        .use_global_bias(true)
        .read()?;
    assert_row_close(
        corrected.to_dense_matrix()?.values_row_major(),
        &[1.0, 0.0, 5.0 / 3.0, 0.5, 8.0, 0.0],
    );
    Ok(())
}

/// Verify global reference broadcasting works for each two-sided correction mode.
#[test]
fn two_sided_correction_modes_allow_global_bias_opt_in() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("global.ref_kmers.zarr");

    write_windowed_end_motif_store(
        &ends_path,
        FixtureStorage::Dense,
        &["A_C", "A_G", "T_C", "T_G"],
        &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
    )?;
    write_global_ref_kmer_store(
        &ref_path,
        FixtureStorage::Dense,
        &["AC", "AG", "TC", "TG"],
        &[0.25, 0.25, 0.25, 0.25],
    )?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    for (mode, expected) in [
        (
            TwoSidedCorrectionMode::Joint,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
        ),
        (
            TwoSidedCorrectionMode::Split,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
        ),
        (TwoSidedCorrectionMode::Outside, vec![3.0, 7.0, 11.0, 15.0]),
        (TwoSidedCorrectionMode::Inside, vec![4.0, 6.0, 12.0, 14.0]),
    ] {
        let corrected = ends
            .select_corrected_counts(&ref_kmers)
            .two_sided_correction(mode)
            .use_global_bias(true)
            .read()?;
        assert_row_close(corrected.to_dense_matrix()?.values_row_major(), &expected);
    }
    Ok(())
}

/// Verify the global-bias option is rejected when the reference is not global.
#[test]
fn corrected_counts_reject_global_bias_option_for_non_global_reference() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("reference.ref_kmers.zarr");

    write_windowed_end_motif_store(
        &ends_path,
        FixtureStorage::Dense,
        &["_AA", "_CC", "_GG"],
        &[1.0, 0.0, 2.5, 0.5, 4.0, 0.0],
    )?;
    write_windowed_ref_kmer_store(
        &ref_path,
        FixtureStorage::Dense,
        &["AA", "CC", "GG"],
        &[1.0 / 3.0, 1.0 / 6.0, 0.5, 0.5, 0.25, 0.25],
    )?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    let error = ends
        .select_corrected_counts(&ref_kmers)
        .use_global_bias(true)
        .read()
        .expect_err("global-bias option should require a global reference");

    assert!(
        error
            .to_string()
            .contains("use_global_bias(true) requires a global reference k-mer output")
    );
    Ok(())
}

/// Verify the global-bias option is harmless when both outputs are global.
#[test]
fn corrected_counts_allow_global_bias_option_for_global_pair() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("global.ref_kmers.zarr");

    write_global_end_motif_store(
        &ends_path,
        FixtureStorage::Dense,
        &["_AA", "_CC"],
        &[2.0, 4.0],
    )?;
    write_global_ref_kmer_store(&ref_path, FixtureStorage::Dense, &["AA", "CC"], &[0.5, 0.5])?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    let corrected = ends
        .select_corrected_counts(&ref_kmers)
        .use_global_bias(true)
        .read()?;

    assert_row_close(corrected.dense_counts()?.values_row_major(), &[2.0, 4.0]);
    Ok(())
}

/// Verify missing reference motifs error by default and can be kept as NaN.
#[test]
fn corrected_counts_can_keep_unsupported_positive_motifs_as_nan() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("reference.ref_kmers.zarr");

    write_global_end_motif_store(
        &ends_path,
        FixtureStorage::Dense,
        &["_AA", "_CC", "_GG"],
        &[0.0, 4.0, 2.0],
    )?;
    write_global_ref_kmer_store(&ref_path, FixtureStorage::Dense, &["AA", "GG"], &[0.5, 0.5])?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    let error = ends
        .select_corrected_counts(&ref_kmers)
        .read()
        .expect_err("positive missing reference motif should fail by default");
    let message = error.to_string();
    assert!(message.contains("positive-count end motifs have no positive correction denominator"));
    assert!(message.contains("_CC"));

    let corrected = ends
        .select_corrected_counts(&ref_kmers)
        .unsupported_reference_policy(UnsupportedReferencePolicy::KeepNaN)
        .read()?;
    let values = corrected.dense_counts()?.values_row_major();
    assert_eq!(values[0], 0.0);
    assert!(values[1].is_nan());
    assert_close(values[2], 2.0);
    Ok(())
}

/// Verify finite inputs cannot silently produce infinite corrected counts.
#[test]
fn corrected_counts_reject_non_finite_results() -> anyhow::Result<()> {
    for end_storage in [FixtureStorage::Dense, FixtureStorage::Sparse] {
        for ref_storage in [FixtureStorage::Dense, FixtureStorage::Sparse] {
            let temp = TempDir::new()?;
            let ends_path = temp.path().join(format!("ends_{end_storage:?}.zarr"));
            let ref_path = temp.path().join(format!("ref_{ref_storage:?}.zarr"));

            write_global_end_motif_store(&ends_path, end_storage, &["_AA"], &[1.0])?;
            write_global_ref_kmer_store(&ref_path, ref_storage, &["AA"], &[f64::from_bits(1)])?;

            let ends = load_ends_output(&ends_path)?;
            let ref_kmers = load_ref_kmers_output(&ref_path)?;
            let error = ends
                .select_corrected_counts(&ref_kmers)
                .read()
                .expect_err("an infinite corrected count should fail");
            let message = error.to_string();
            assert!(message.contains("non-finite corrected counts"));
            assert!(message.contains("_AA"));
        }
    }
    Ok(())
}

/// Verify global row broadcasting does not bypass motif-axis compatibility.
#[test]
fn global_correction_rejects_grouped_motifs_with_concrete_reference() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("global.ref_kmers.zarr");

    let ends_store = create_store(
        &ends_path,
        end_root_attributes_with_axis(FixtureStorage::Dense, "bed", "motif_group"),
    )?;
    write_motif_group_axis(&ends_store, &["left", "right"])?;
    write_window_rows(&ends_store)?;
    write_f64_array(
        &ends_store,
        "counts",
        &[2, 2],
        &["row", "motif"],
        &[1.0, 2.0, 3.0, 4.0],
        json!({}),
    )?;
    write_global_ref_kmer_store(&ref_path, FixtureStorage::Dense, &["AA", "CC"], &[0.5, 0.5])?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    for mode in [
        TwoSidedCorrectionMode::Joint,
        TwoSidedCorrectionMode::Split,
        TwoSidedCorrectionMode::Outside,
        TwoSidedCorrectionMode::Inside,
    ] {
        let mode_error = ends
            .select_corrected_counts(&ref_kmers)
            .two_sided_correction(mode)
            .use_global_bias(true)
            .read()
            .expect_err("motif-group outputs should reject explicit two-sided modes");
        assert!(mode_error.to_string().contains("motif-group"));
    }

    let error = ends
        .select_corrected_counts(&ref_kmers)
        .use_global_bias(true)
        .read()
        .expect_err("motif-axis mismatch should fail before global row broadcasting");

    assert!(
        error
            .to_string()
            .contains("grouped end-motif output requires grouped reference k-mer output")
    );
    Ok(())
}

/// Verify canonical reference k-mer outputs are rejected for correction.
#[test]
fn corrected_counts_reject_canonical_reference() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("canonical.ref_kmers.zarr");

    write_windowed_end_motif_store(
        &ends_path,
        FixtureStorage::Dense,
        &["_AA", "_CC", "_GG"],
        &[1.0, 0.0, 2.5, 0.5, 4.0, 0.0],
    )?;
    let mut reference_attributes = ref_root_attributes(FixtureStorage::Dense, "bed", &["AA", "CC"]);
    reference_attributes["canonical"] = json!(true);
    let ref_store = create_store(&ref_path, reference_attributes)?;
    write_motif_axis(&ref_store, &["AA", "CC"])?;
    write_window_rows(&ref_store)?;
    write_row_scaling_factors(&ref_store, 2)?;
    write_reference_footprint(&ref_store)?;
    write_ref_frequencies(
        &ref_store,
        FixtureStorage::Dense,
        2,
        2,
        &[0.5, 0.5, 0.5, 0.5],
    )?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    let error = ends
        .select_corrected_counts(&ref_kmers)
        .read()
        .expect_err("canonical reference correction should fail");

    assert!(error.to_string().contains("non-canonical"));
    Ok(())
}

/// Verify two-sided outputs require a mode and joint mode uses full reference motifs.
#[test]
fn two_sided_correction_requires_mode_and_joint_uses_full_reference_motifs() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("reference.ref_kmers.zarr");

    write_global_end_motif_store(
        &ends_path,
        FixtureStorage::Dense,
        &["A_C", "A_G", "T_C", "T_G"],
        &[1.0, 2.0, 3.0, 4.0],
    )?;
    write_global_ref_kmer_store(
        &ref_path,
        FixtureStorage::Dense,
        &["AC", "AG", "TC", "TG"],
        &[0.25, 0.25, 0.25, 0.25],
    )?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    let error = ends
        .select_corrected_counts(&ref_kmers)
        .read()
        .expect_err("two-sided correction should require an explicit mode");
    assert!(error.to_string().contains("two-sided"));

    let corrected = ends
        .select_corrected_counts(&ref_kmers)
        .two_sided_correction(TwoSidedCorrectionMode::Joint)
        .read()?;

    // Four positive joint motifs make the uniform baseline frequency 1/4.
    // Every reference frequency is also 1/4, so every correction factor is 1.
    assert_eq!(
        corrected.motif_labels(),
        &[
            "A_C".to_string(),
            "A_G".to_string(),
            "T_C".to_string(),
            "T_G".to_string(),
        ]
    );
    assert_row_close(
        corrected.dense_counts()?.values_row_major(),
        &[1.0, 2.0, 3.0, 4.0],
    );
    Ok(())
}

/// Verify an empty loaded motif axis does not require side-width inference.
#[test]
fn empty_motif_axis_returns_empty_correction_for_any_explicit_mode() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("empty.end_motifs.zarr");
    let ref_path = temp.path().join("reference.ref_kmers.zarr");

    write_global_end_motif_store(&ends_path, FixtureStorage::Dense, &[], &[])?;
    write_global_ref_kmer_store(&ref_path, FixtureStorage::Dense, &["AA"], &[1.0])?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    for mode in [
        TwoSidedCorrectionMode::Joint,
        TwoSidedCorrectionMode::Split,
        TwoSidedCorrectionMode::Outside,
        TwoSidedCorrectionMode::Inside,
    ] {
        let corrected = ends
            .select_corrected_counts(&ref_kmers)
            .two_sided_correction(mode)
            .read()?;
        assert_eq!(corrected.shape(), (1, 0));
        assert!(corrected.motif_labels().is_empty());
    }
    Ok(())
}

/// Verify every correction mode against the shared cross-language math fixture.
#[test]
fn correction_modes_match_shared_cross_language_expectations() -> anyhow::Result<()> {
    for end_storage in [FixtureStorage::Dense, FixtureStorage::Sparse] {
        for ref_storage in [FixtureStorage::Dense, FixtureStorage::Sparse] {
            let temp = TempDir::new()?;
            let ends_path = temp.path().join(format!("ends_{end_storage:?}.zarr"));
            let ref_path = temp.path().join(format!("ref_{ref_storage:?}.zarr"));

            write_global_end_motif_store(
                &ends_path,
                end_storage,
                &["A_C", "A_G", "T_C", "T_G"],
                &[2.0, 4.0, 6.0, 8.0],
            )?;
            write_global_ref_kmer_store(
                &ref_path,
                ref_storage,
                &["AC", "AG", "TC", "TG"],
                &[1.0 / 8.0, 1.0 / 8.0, 1.0 / 4.0, 1.0 / 2.0],
            )?;

            let ends = load_ends_output(&ends_path)?;
            let ref_kmers = load_ref_kmers_output(&ref_path)?;
            let joint = ends
                .select_corrected_counts(&ref_kmers)
                .two_sided_correction(TwoSidedCorrectionMode::Joint)
                .read()?;
            let split = ends
                .select_corrected_counts(&ref_kmers)
                .two_sided_correction(TwoSidedCorrectionMode::Split)
                .read()?;
            let outside = ends
                .select_corrected_counts(&ref_kmers)
                .two_sided_correction(TwoSidedCorrectionMode::Outside)
                .read()?;
            let inside = ends
                .select_corrected_counts(&ref_kmers)
                .two_sided_correction(TwoSidedCorrectionMode::Inside)
                .read()?;

            // Four positive reference motifs make the uniform frequency 1/4.
            // Relative to uniform, frequencies [1/8, 1/8, 1/4, 1/2] give
            // factors [1/2, 1/2, 1, 2] for [AC, AG, TC, TG]. Dividing counts
            // [2, 4, 6, 8] by them gives corrected counts [4, 8, 6, 4].
            assert_row_close(
                joint.to_dense_matrix()?.values_row_major(),
                &[4.0, 8.0, 6.0, 4.0],
            );
            // Each side has uniform frequency 1/2. Outside frequencies A=1/4
            // and T=3/4 give factors 1/2 and 3/2. Inside frequencies C=3/8
            // and G=5/8 give factors 3/4 and 5/4. For
            // [A_C, A_G, T_C, T_G], split multiplies matching side factors to
            // get [3/8, 5/8, 9/8, 15/8]. Dividing counts [2, 4, 6, 8] by them
            // gives [16/3, 32/5, 16/3, 64/15].
            assert_row_close(
                split.to_dense_matrix()?.values_row_major(),
                &[16.0 / 3.0, 32.0 / 5.0, 16.0 / 3.0, 64.0 / 15.0],
            );
            // Outside aggregates counts to [6, 14]. Relative to uniform 1/2,
            // reference frequencies [1/4, 3/4] give factors [1/2, 3/2].
            // Dividing the aggregated counts by them gives [12, 28/3].
            assert_row_close(
                outside.to_dense_matrix()?.values_row_major(),
                &[12.0, 28.0 / 3.0],
            );
            // Inside aggregates counts to [8, 12]. Relative to uniform 1/2,
            // reference frequencies [3/8, 5/8] give factors [3/4, 5/4].
            // Dividing the aggregated counts by them gives [32/3, 48/5].
            assert_row_close(
                inside.to_dense_matrix()?.values_row_major(),
                &[32.0 / 3.0, 48.0 / 5.0],
            );
            assert_eq!(
                split.storage_mode(),
                match end_storage {
                    FixtureStorage::Dense => EndMotifStorageMode::Dense,
                    FixtureStorage::Sparse => EndMotifStorageMode::SparseCoo,
                }
            );
        }
    }
    Ok(())
}

/// Verify side-mode aggregation cannot silently overflow finite input counts.
#[test]
fn side_correction_rejects_non_finite_aggregated_counts() -> anyhow::Result<()> {
    for end_storage in [FixtureStorage::Dense, FixtureStorage::Sparse] {
        let temp = TempDir::new()?;
        let ends_path = temp.path().join(format!("ends_{end_storage:?}.zarr"));
        let ref_path = temp.path().join("reference.ref_kmers.zarr");

        write_global_end_motif_store(
            &ends_path,
            end_storage,
            &["A_C", "A_G"],
            &[f64::MAX, f64::MAX],
        )?;
        write_global_ref_kmer_store(&ref_path, FixtureStorage::Dense, &["AC", "AG"], &[0.5, 0.5])?;

        let ends = load_ends_output(&ends_path)?;
        let ref_kmers = load_ref_kmers_output(&ref_path)?;
        let error = ends
            .select_corrected_counts(&ref_kmers)
            .two_sided_correction(TwoSidedCorrectionMode::Outside)
            .read()
            .expect_err("an infinite aggregated side count should fail");
        assert!(
            error
                .to_string()
                .contains("side-mode aggregation produced non-finite count for motif 'A_'")
        );
    }
    Ok(())
}

/// Verify unsupported-reference handling is applied after side aggregation.
#[test]
fn side_correction_applies_unsupported_policy_to_aggregated_counts() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("reference.ref_kmers.zarr");

    write_global_end_motif_store(
        &ends_path,
        FixtureStorage::Dense,
        &["A_C", "A_G"],
        &[2.0, 4.0],
    )?;
    write_global_ref_kmer_store(&ref_path, FixtureStorage::Dense, &["AC"], &[1.0])?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;
    let error = ends
        .select_corrected_counts(&ref_kmers)
        .two_sided_correction(TwoSidedCorrectionMode::Inside)
        .read()
        .expect_err("positive aggregated side count without reference support should fail");
    assert!(error.to_string().contains("_G"));

    let corrected = ends
        .select_corrected_counts(&ref_kmers)
        .two_sided_correction(TwoSidedCorrectionMode::Inside)
        .unsupported_reference_policy(UnsupportedReferencePolicy::KeepNaN)
        .read()?;
    assert_eq!(corrected.motif_labels(), &["_C", "_G"]);
    let values = corrected.dense_counts()?.values_row_major();
    assert_close(values[0], 2.0);
    assert!(values[1].is_nan());
    Ok(())
}

/// Verify side modes aggregate the loaded joint axis before side-label selection.
#[test]
fn side_correction_aggregates_joint_counts_in_loaded_axis_order() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let ends_path = temp.path().join("sample.end_motifs.zarr");
    let ref_path = temp.path().join("reference.ref_kmers.zarr");

    write_global_end_motif_store(
        &ends_path,
        FixtureStorage::Dense,
        &["T_G", "A_C", "T_C", "A_G"],
        &[8.0, 2.0, 6.0, 4.0],
    )?;
    write_global_ref_kmer_store(
        &ref_path,
        FixtureStorage::Dense,
        &["TG", "AC", "TC", "AG"],
        &[0.25, 0.5, 0.25, 0.0],
    )?;

    let ends = load_ends_output(&ends_path)?;
    let ref_kmers = load_ref_kmers_output(&ref_path)?;

    let outside = ends
        .select_corrected_counts(&ref_kmers)
        .two_sided_correction(TwoSidedCorrectionMode::Outside)
        .read()?;
    assert_eq!(
        outside.motif_labels(),
        &["T_".to_string(), "A_".to_string()]
    );
    assert_eq!(outside.motif_indices(), &[0, 1]);
    // The side axis follows first loaded occurrence: T_ before A_. The uniform
    // outside baseline and both reference frequencies are 0.5, so both
    // correction factors are 1 and corrected counts equal aggregated counts.
    assert_row_close(outside.dense_counts()?.values_row_major(), &[14.0, 6.0]);

    let selected_outside = ends
        .select_corrected_counts(&ref_kmers)
        .two_sided_correction(TwoSidedCorrectionMode::Outside)
        .motifs_by_label(&["A_"])
        .read()?;
    assert_eq!(selected_outside.motif_labels(), &["A_".to_string()]);
    assert_eq!(selected_outside.motif_indices(), &[1]);
    assert_row_close(selected_outside.dense_counts()?.values_row_major(), &[6.0]);

    let inside = ends
        .select_corrected_counts(&ref_kmers)
        .two_sided_correction(TwoSidedCorrectionMode::Inside)
        .read()?;
    assert_eq!(inside.motif_labels(), &["_G".to_string(), "_C".to_string()]);
    // The selected side order is [_G, _C]. Relative to uniform frequency 0.5,
    // reference frequencies G=0.25 and C=0.75 give factors 0.5 and 1.5.
    // Dividing aggregated counts G=12 and C=8 gives 24 and 16/3.
    assert_row_close(
        inside.dense_counts()?.values_row_major(),
        &[24.0, 16.0 / 3.0],
    );
    Ok(())
}

/// Verify mode validation keeps one-sided correction and side-mode selectors unambiguous.
#[test]
fn mode_validation_rejects_one_sided_modes_and_side_index_selectors() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let one_sided_ends_path = temp.path().join("one_sided.end_motifs.zarr");
    let one_sided_ref_path = temp.path().join("one_sided.ref_kmers.zarr");
    let two_sided_ends_path = temp.path().join("two_sided.end_motifs.zarr");
    let two_sided_ref_path = temp.path().join("two_sided.ref_kmers.zarr");

    write_global_end_motif_store(
        &one_sided_ends_path,
        FixtureStorage::Dense,
        &["_AA", "_CC"],
        &[2.0, 4.0],
    )?;
    write_global_ref_kmer_store(
        &one_sided_ref_path,
        FixtureStorage::Dense,
        &["AA", "CC"],
        &[0.5, 0.5],
    )?;
    let one_sided_ends = load_ends_output(&one_sided_ends_path)?;
    let one_sided_ref_kmers = load_ref_kmers_output(&one_sided_ref_path)?;
    let exact = one_sided_ends
        .select_corrected_counts(&one_sided_ref_kmers)
        .read()?;
    assert_row_close(exact.dense_counts()?.values_row_major(), &[2.0, 4.0]);
    for mode in [
        TwoSidedCorrectionMode::Joint,
        TwoSidedCorrectionMode::Split,
        TwoSidedCorrectionMode::Outside,
        TwoSidedCorrectionMode::Inside,
    ] {
        let error = one_sided_ends
            .select_corrected_counts(&one_sided_ref_kmers)
            .two_sided_correction(mode)
            .read()
            .expect_err("one-sided output should reject explicit two-sided modes");
        assert!(error.to_string().contains("one-sided"));
    }

    write_global_end_motif_store(
        &two_sided_ends_path,
        FixtureStorage::Dense,
        &["A_C", "A_G"],
        &[2.0, 4.0],
    )?;
    write_global_ref_kmer_store(
        &two_sided_ref_path,
        FixtureStorage::Dense,
        &["AC", "AG"],
        &[0.5, 0.5],
    )?;
    let two_sided_ends = load_ends_output(&two_sided_ends_path)?;
    let two_sided_ref_kmers = load_ref_kmers_output(&two_sided_ref_path)?;
    let index_error = two_sided_ends
        .select_corrected_counts(&two_sided_ref_kmers)
        .two_sided_correction(TwoSidedCorrectionMode::Outside)
        .motifs(&[0])
        .read()
        .expect_err("side modes should reject motif index selectors");
    assert!(index_error.to_string().contains("motif index selectors"));

    let wrong_axis_error = two_sided_ends
        .select_corrected_counts(&two_sided_ref_kmers)
        .two_sided_correction(TwoSidedCorrectionMode::Outside)
        .motifs_by_label(&["A_C"])
        .read()
        .expect_err("side modes should reject joint labels");
    assert!(
        wrong_axis_error
            .to_string()
            .contains("side-mode motif axis")
    );
    Ok(())
}

fn write_windowed_end_motif_store(
    path: &Path,
    storage: FixtureStorage,
    motif_labels: &[&str],
    dense_counts: &[f64],
) -> anyhow::Result<()> {
    let row_count = 2;
    let motif_count = motif_labels.len();
    let store = create_store(path, end_root_attributes(storage, "bed"))?;
    write_motif_axis(&store, motif_labels)?;
    write_window_rows(&store)?;
    write_end_counts(&store, storage, row_count, motif_count, dense_counts)
}

fn write_global_end_motif_store(
    path: &Path,
    storage: FixtureStorage,
    motif_labels: &[&str],
    dense_counts: &[f64],
) -> anyhow::Result<()> {
    let store = create_store(path, end_root_attributes(storage, "global"))?;
    write_motif_axis(&store, motif_labels)?;
    write_global_row_axis(&store)?;
    write_end_counts(&store, storage, 1, motif_labels.len(), dense_counts)
}

fn write_windowed_ref_kmer_store(
    path: &Path,
    storage: FixtureStorage,
    motif_labels: &[&str],
    dense_frequencies: &[f64],
) -> anyhow::Result<()> {
    let row_count = 2;
    let motif_count = motif_labels.len();
    let store = create_store(path, ref_root_attributes(storage, "bed", motif_labels))?;
    write_motif_axis(&store, motif_labels)?;
    write_window_rows(&store)?;
    write_row_scaling_factors(&store, row_count)?;
    write_reference_footprint(&store)?;
    write_ref_frequencies(&store, storage, row_count, motif_count, dense_frequencies)
}

fn write_windowed_ref_kmer_sparse_store(
    path: &Path,
    motif_labels: &[&str],
    row_indices: &[i32],
    motif_indices: &[i32],
    frequencies: &[f64],
) -> anyhow::Result<()> {
    let row_count = 2;
    let store = create_store(
        path,
        ref_root_attributes(FixtureStorage::Sparse, "bed", motif_labels),
    )?;
    write_motif_axis(&store, motif_labels)?;
    write_window_rows(&store)?;
    write_row_scaling_factors(&store, row_count)?;
    write_reference_footprint(&store)?;
    write_sparse_ref_frequencies(
        &store,
        row_count,
        motif_labels.len(),
        row_indices,
        motif_indices,
        frequencies,
    )
}

fn write_global_ref_kmer_store(
    path: &Path,
    storage: FixtureStorage,
    motif_labels: &[&str],
    dense_frequencies: &[f64],
) -> anyhow::Result<()> {
    let store = create_store(path, ref_root_attributes(storage, "global", motif_labels))?;
    write_motif_axis(&store, motif_labels)?;
    write_global_row_axis(&store)?;
    write_row_scaling_factors(&store, 1)?;
    write_reference_footprint(&store)?;
    write_ref_frequencies(&store, storage, 1, motif_labels.len(), dense_frequencies)
}

fn end_root_attributes(storage: FixtureStorage, row_mode: &str) -> Value {
    end_root_attributes_with_axis(storage, row_mode, "motif")
}

fn end_root_attributes_with_axis(
    storage: FixtureStorage,
    row_mode: &str,
    motif_axis_kind: &str,
) -> Value {
    let (storage_mode, primary_array, primary_group) = storage_metadata(storage, "counts");
    let mut attributes = json!({
        "cfdnalab_schema": "end_motif_counts",
        "cfdnalab_schema_version": 2,
        "storage_mode": storage_mode,
        "row_mode": row_mode,
        "motif_axis_kind": motif_axis_kind,
        "count_units": "weighted_end_motif_count",
        "primary_array": primary_array,
        "primary_group": primary_group,
    });
    if matches!(storage, FixtureStorage::Sparse) {
        attributes["sparse_format"] = json!("coo");
        attributes["sparse_indices_base"] = json!(0);
    }
    attributes
}

fn ref_root_attributes(storage: FixtureStorage, row_mode: &str, motif_labels: &[&str]) -> Value {
    let (storage_mode, primary_array, primary_group) = storage_metadata(storage, "frequencies");
    let mut attributes = json!({
        "cfdnalab_schema": "ref_kmer_frequencies",
        "cfdnalab_schema_version": 1,
        "storage_mode": storage_mode,
        "row_mode": row_mode,
        "motif_axis_kind": "motif",
        "value_units": "reference_kmer_frequency",
        "count_units": "reference_kmer_count",
        "row_scaling_factor_array": "row_scaling_factor",
        "count_reconstruction": "reference_kmer_count = frequency * row_scaling_factor[row]",
        "kmer_size": motif_labels.first().map_or(0, |label| label.len()),
        "canonical": false,
        "all_motifs": false,
        "assign_by": "count-overlap",
        "primary_array": primary_array,
        "primary_group": primary_group,
    });
    if matches!(storage, FixtureStorage::Sparse) {
        attributes["sparse_format"] = json!("coo");
        attributes["sparse_indices_base"] = json!(0);
    }
    attributes
}

fn storage_metadata(
    storage: FixtureStorage,
    primary_array_name: &str,
) -> (&'static str, Value, Value) {
    match storage {
        FixtureStorage::Dense => ("dense", json!(primary_array_name), Value::Null),
        FixtureStorage::Sparse => ("sparse_coo", Value::Null, json!("sparse")),
    }
}

fn write_end_counts(
    store: &Arc<FilesystemStore>,
    storage: FixtureStorage,
    row_count: usize,
    motif_count: usize,
    dense_counts: &[f64],
) -> anyhow::Result<()> {
    match storage {
        FixtureStorage::Dense => write_f64_array(
            store,
            "counts",
            &[row_count, motif_count],
            &["row", "motif"],
            dense_counts,
            json!({}),
        ),
        FixtureStorage::Sparse => {
            let (row_indices, motif_indices, counts) =
                dense_to_sparse_entries(dense_counts, motif_count);
            write_sparse_end_counts(
                store,
                row_count,
                motif_count,
                &row_indices,
                &motif_indices,
                &counts,
            )
        }
    }
}

fn write_ref_frequencies(
    store: &Arc<FilesystemStore>,
    storage: FixtureStorage,
    row_count: usize,
    motif_count: usize,
    dense_frequencies: &[f64],
) -> anyhow::Result<()> {
    match storage {
        FixtureStorage::Dense => write_f64_array(
            store,
            "frequencies",
            &[row_count, motif_count],
            &["row", "motif"],
            dense_frequencies,
            json!({}),
        ),
        FixtureStorage::Sparse => {
            let (row_indices, motif_indices, frequencies) =
                dense_to_sparse_entries(dense_frequencies, motif_count);
            write_sparse_ref_frequencies(
                store,
                row_count,
                motif_count,
                &row_indices,
                &motif_indices,
                &frequencies,
            )
        }
    }
}

fn dense_to_sparse_entries(values: &[f64], motif_count: usize) -> (Vec<i32>, Vec<i32>, Vec<f64>) {
    let mut row_indices = Vec::new();
    let mut motif_indices = Vec::new();
    let mut sparse_values = Vec::new();
    for (value_index, value) in values.iter().copied().enumerate() {
        if value == 0.0 {
            continue;
        }
        row_indices.push(i32::try_from(value_index / motif_count).unwrap());
        motif_indices.push(i32::try_from(value_index % motif_count).unwrap());
        sparse_values.push(value);
    }
    (row_indices, motif_indices, sparse_values)
}

fn write_motif_axis(store: &Arc<FilesystemStore>, labels: &[&str]) -> anyhow::Result<()> {
    let motif_count = labels.len();
    let motif_index = (0..motif_count)
        .map(|motif_index| i32::try_from(motif_index).expect("test motif index should fit i32"))
        .collect::<Vec<_>>();
    write_i32_array(
        store,
        "motif_index",
        &[motif_count],
        &["motif"],
        &motif_index,
        json!({}),
    )?;

    let motif_width = labels.first().map_or(0, |label| label.len());
    let motif_byte = (0..motif_width)
        .map(|motif_byte| i32::try_from(motif_byte).expect("test motif byte should fit i32"))
        .collect::<Vec<_>>();
    write_i32_array(
        store,
        "motif_byte",
        &[motif_width],
        &["motif_byte"],
        &motif_byte,
        json!({}),
    )?;

    let mut motif_ascii = Vec::with_capacity(motif_count.saturating_mul(motif_width));
    for label in labels {
        assert_eq!(label.len(), motif_width);
        motif_ascii.extend_from_slice(label.as_bytes());
    }
    write_u8_array(
        store,
        "motif_ascii",
        &[motif_count, motif_width],
        &["motif", "motif_byte"],
        &motif_ascii,
        json!({}),
    )
}

fn write_motif_group_axis(store: &Arc<FilesystemStore>, labels: &[&str]) -> anyhow::Result<()> {
    let motif_count = labels.len();
    let motif_index = (0..motif_count)
        .map(|motif_index| i32::try_from(motif_index).expect("test motif index should fit i32"))
        .collect::<Vec<_>>();
    write_i32_array(
        store,
        "motif_index",
        &[motif_count],
        &["motif"],
        &motif_index,
        json!({
            "label_field": "motif_group",
            "labels": labels,
        }),
    )
}

fn write_global_row_axis(store: &Arc<FilesystemStore>) -> anyhow::Result<()> {
    write_i32_array(
        store,
        "row",
        &[1],
        &["row"],
        &[0],
        json!({
            "label_field": "row_label",
            "labels": ["global"],
        }),
    )
}

fn write_window_rows(store: &Arc<FilesystemStore>) -> anyhow::Result<()> {
    write_i32_array(store, "row", &[2], &["row"], &[0, 1], json!({}))?;
    write_i32_array(
        store,
        "chromosome",
        &[1],
        &["chromosome"],
        &[0],
        json!({
            "label_field": "chromosome_name",
            "labels": ["chr1"],
        }),
    )?;
    write_i32_array(store, "row_chromosome", &[2], &["row"], &[0, 0], json!({}))?;
    write_i64_array(store, "row_start_bp", &[2], &["row"], &[0, 100], json!({}))?;
    write_i64_array(store, "row_end_bp", &[2], &["row"], &[50, 150], json!({}))?;
    write_f64_array(
        store,
        "blacklisted_fraction",
        &[2],
        &["row"],
        &[0.0, 0.0],
        json!({}),
    )
}

fn write_group_rows(store: &Arc<FilesystemStore>, group_labels: &[&str]) -> anyhow::Result<()> {
    let row_count = group_labels.len();
    let row = (0..row_count)
        .map(|row_index| i32::try_from(row_index).expect("test row index should fit i32"))
        .collect::<Vec<_>>();
    write_i32_array(store, "row", &[row_count], &["row"], &row, json!({}))?;
    write_i32_array(
        store,
        "group",
        &[row_count],
        &["row"],
        &row,
        json!({
            "label_field": "group_name",
            "labels": group_labels,
        }),
    )?;
    write_i32_array(
        store,
        "eligible_windows",
        &[row_count],
        &["row"],
        &vec![1; row_count],
        json!({}),
    )?;
    write_f64_array(
        store,
        "blacklisted_fraction",
        &[row_count],
        &["row"],
        &vec![0.0; row_count],
        json!({}),
    )
}

fn write_row_scaling_factors(store: &Arc<FilesystemStore>, row_count: usize) -> anyhow::Result<()> {
    write_f64_array(
        store,
        "row_scaling_factor",
        &[row_count],
        &["row"],
        &vec![1.0; row_count],
        json!({}),
    )
}

fn write_reference_footprint(store: &Arc<FilesystemStore>) -> anyhow::Result<()> {
    write_u8_array(
        store,
        "reference_contig_footprint_json",
        &[2],
        &["json_byte"],
        b"[]",
        json!({}),
    )
}

fn write_sparse_end_counts(
    store: &Arc<FilesystemStore>,
    row_count: usize,
    motif_count: usize,
    row_indices: &[i32],
    motif_indices: &[i32],
    counts: &[f64],
) -> anyhow::Result<()> {
    write_group(store, "/sparse", json!({}))?;
    write_i32_array(
        store,
        "sparse/row",
        &[row_indices.len()],
        &["nnz"],
        row_indices,
        json!({}),
    )?;
    write_i32_array(
        store,
        "sparse/motif",
        &[motif_indices.len()],
        &["nnz"],
        motif_indices,
        json!({}),
    )?;
    write_f64_array(
        store,
        "sparse/count",
        &[counts.len()],
        &["nnz"],
        counts,
        json!({}),
    )?;
    write_sparse_shape(store, row_count, motif_count)
}

fn write_sparse_ref_frequencies(
    store: &Arc<FilesystemStore>,
    row_count: usize,
    motif_count: usize,
    row_indices: &[i32],
    motif_indices: &[i32],
    frequencies: &[f64],
) -> anyhow::Result<()> {
    write_group(store, "/sparse", json!({}))?;
    write_i32_array(
        store,
        "sparse/row",
        &[row_indices.len()],
        &["nnz"],
        row_indices,
        json!({}),
    )?;
    write_i32_array(
        store,
        "sparse/motif",
        &[motif_indices.len()],
        &["nnz"],
        motif_indices,
        json!({}),
    )?;
    write_f64_array(
        store,
        "sparse/frequency",
        &[frequencies.len()],
        &["nnz"],
        frequencies,
        json!({}),
    )?;
    write_sparse_shape(store, row_count, motif_count)
}

fn write_sparse_shape(
    store: &Arc<FilesystemStore>,
    row_count: usize,
    motif_count: usize,
) -> anyhow::Result<()> {
    write_i32_array(
        store,
        "sparse/shape",
        &[2],
        &["sparse_dimension"],
        &[i32::try_from(row_count)?, i32::try_from(motif_count)?],
        json!({}),
    )?;
    write_i32_array(
        store,
        "sparse/sparse_dimension",
        &[2],
        &["sparse_dimension"],
        &[0, 1],
        json!({
            "label_field": "sparse_dimension_name",
            "labels": ["row", "motif"],
        }),
    )
}

fn create_store(path: &Path, attributes: Value) -> anyhow::Result<Arc<FilesystemStore>> {
    fs::create_dir_all(path)?;
    let store = Arc::new(FilesystemStore::new(path)?);
    write_group(&store, "/", attributes)?;
    Ok(store)
}

fn write_group(
    store: &Arc<FilesystemStore>,
    group_path: &str,
    attributes: Value,
) -> anyhow::Result<()> {
    let group = GroupBuilder::new()
        .attributes(json_object(attributes)?)
        .build(store.clone(), group_path)?;
    group.store_metadata()?;
    Ok(())
}

fn write_i32_array(
    store: &Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    dimension_names: &[&str],
    values: &[i32],
    attributes: Value,
) -> anyhow::Result<()> {
    write_array(
        store,
        name,
        shape,
        dimension_names,
        values,
        data_type::int32(),
        -1,
        attributes,
    )
}

fn write_i64_array(
    store: &Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    dimension_names: &[&str],
    values: &[i64],
    attributes: Value,
) -> anyhow::Result<()> {
    write_array(
        store,
        name,
        shape,
        dimension_names,
        values,
        data_type::int64(),
        -1,
        attributes,
    )
}

fn write_u8_array(
    store: &Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    dimension_names: &[&str],
    values: &[u8],
    attributes: Value,
) -> anyhow::Result<()> {
    write_array(
        store,
        name,
        shape,
        dimension_names,
        values,
        data_type::uint8(),
        u8::MAX,
        attributes,
    )
}

fn write_f64_array(
    store: &Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    dimension_names: &[&str],
    values: &[f64],
    attributes: Value,
) -> anyhow::Result<()> {
    write_array(
        store,
        name,
        shape,
        dimension_names,
        values,
        data_type::float64(),
        -1.0,
        attributes,
    )
}

#[allow(clippy::too_many_arguments)]
fn write_array<T>(
    store: &Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    dimension_names: &[&str],
    values: &[T],
    data_type: DataType,
    fill_value: T,
    attributes: Value,
) -> anyhow::Result<()>
where
    T: Element + Copy,
    ArrayBuilderFillValue: From<T>,
{
    let expected_len = shape
        .iter()
        .try_fold(1usize, |size, dimension| size.checked_mul(*dimension))
        .expect("test Zarr shape should not overflow");
    assert_eq!(values.len(), expected_len);
    let chunk_shape = shape
        .iter()
        .map(|dimension| (*dimension).max(1))
        .collect::<Vec<_>>();
    let array = ArrayBuilder::new(
        shape
            .iter()
            .map(|dimension| *dimension as u64)
            .collect::<Vec<_>>(),
        chunk_shape
            .iter()
            .map(|dimension| *dimension as u64)
            .collect::<Vec<_>>(),
        data_type,
        fill_value,
    )
    .dimension_names(Some(dimension_names.iter().copied()))
    .attributes(json_object(attributes)?)
    .build(store.clone(), format!("/{name}").as_str())?;
    array.store_metadata()?;
    if !values.is_empty() {
        array.store_chunk(&vec![0; shape.len()], values)?;
    }
    Ok(())
}

fn json_object(value: Value) -> anyhow::Result<Map<String, Value>> {
    match value {
        Value::Object(map) => Ok(map),
        _ => anyhow::bail!("expected JSON object"),
    }
}

fn assert_row_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (actual_value, expected_value) in actual.iter().zip(expected) {
        assert_close(*actual_value, *expected_value);
    }
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1.0e-12,
        "expected {actual} to equal {expected}"
    );
}
