#![cfg(feature = "cmd_ends")]
//! Public API tests for Rust output loaders for `cfdna ends`.
//!
//! These tests build tiny Zarr V3 stores with the public end-motif schema and
//! assert that the loader returns typed Rust metadata and native vector-backed
//! count containers.

use cfdnalab::output_loaders::{
    EndMotifAxisKind, EndMotifRowMetadata, EndMotifRowMode, EndMotifSparseEntry,
    EndMotifStorageMode, EndMotifWindowMode, load_ends_output,
};
use serde_json::{Map, Value, json};
use std::{fs, path::Path, sync::Arc};
use tempfile::TempDir;
use zarrs::{
    array::{ArrayBuilder, DataType, Element, builder::ArrayBuilderFillValue, data_type},
    filesystem::FilesystemStore,
    group::GroupBuilder,
};

/// Verify dense windowed end-motif stores load metadata, counts, and selections.
#[test]
fn load_ends_output_reads_dense_windowed_motif_store() -> anyhow::Result<()> {
    // Arrange:
    // Dense windowed output stores the full row-by-motif matrix and dictionary
    // encodes chromosome names through row_chromosome.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "dense",
            "row_mode": "bed",
            "motif_axis_kind": "motif",
            "count_units": "weighted_end_motif_count",
            "primary_array": "counts",
            "primary_group": null,
        }),
    )?;
    write_i32_array(&store, "motif_index", &[2], &["motif"], &[0, 1], json!({}))?;
    write_i32_array(
        &store,
        "motif_byte",
        &[3],
        &["motif_byte"],
        &[0, 1, 2],
        json!({}),
    )?;
    write_u8_array(
        &store,
        "motif_ascii",
        &[2, 3],
        &["motif", "motif_byte"],
        b"_AA_GG",
        json!({}),
    )?;
    write_i32_array(&store, "row", &[2], &["row"], &[0, 1], json!({}))?;
    write_i32_array(
        &store,
        "chromosome",
        &[2],
        &["chromosome"],
        &[0, 1],
        json!({
            "label_field": "chromosome_name",
            "labels": ["chr1", "chr2"],
        }),
    )?;
    write_i32_array(&store, "row_chromosome", &[2], &["row"], &[0, 1], json!({}))?;
    write_i64_array(&store, "row_start_bp", &[2], &["row"], &[10, 20], json!({}))?;
    write_i64_array(&store, "row_end_bp", &[2], &["row"], &[15, 25], json!({}))?;
    write_f64_array(
        &store,
        "blacklisted_fraction",
        &[2],
        &["row"],
        &[0.25, 0.0],
        json!({}),
    )?;
    write_f64_array(
        &store,
        "counts",
        &[2, 2],
        &["row", "motif"],
        &[1.0, 2.0, 0.0, 3.5],
        json!({}),
    )?;

    // Act
    let loaded = load_ends_output(&path)?;

    // Assert
    assert_eq!(loaded.storage_mode(), EndMotifStorageMode::Dense);
    assert_eq!(loaded.row_mode(), EndMotifRowMode::BedWindows);
    assert_eq!(loaded.motif_axis_kind(), EndMotifAxisKind::Motif);
    let output_metadata = loaded.output_metadata();
    assert_eq!(output_metadata.storage_mode, EndMotifStorageMode::Dense);
    assert_eq!(output_metadata.row_mode, EndMotifRowMode::BedWindows);
    assert_eq!(output_metadata.motif_axis_kind, EndMotifAxisKind::Motif);
    assert_eq!(output_metadata.row_count, 2);
    assert_eq!(output_metadata.motif_count, 2);
    assert_eq!(
        output_metadata.to_string(),
        "storage_mode=dense, row_mode=BED windows, motif_axis=motifs, row_count=2, motif_count=2"
    );
    assert_eq!(
        loaded.motif_labels(),
        &["_AA".to_string(), "_GG".to_string()]
    );
    assert_eq!(loaded.motif_index("_GG")?, 1);
    assert!(loaded.has_motif("_AA"));
    assert_eq!(loaded.dense_counts()?.shape(), (2, 2));
    assert_eq!(loaded.dense_counts()?.row_count(), 2);
    assert_eq!(loaded.dense_counts()?.column_count(), 2);
    assert_eq!(
        loaded.dense_counts()?.values_row_major(),
        &[1.0, 2.0, 0.0, 3.5]
    );
    assert_eq!(loaded.count(0, 1), Some(2.0));
    assert_eq!(loaded.count_for_motif(1, "_GG")?, Some(3.5));
    assert!(loaded.sparse_counts().is_err());

    let windows = loaded.window_metadata()?;
    assert_eq!(windows.len(), 2);
    assert_eq!(windows[0].index, 0);
    assert_eq!(windows[0].chrom, "chr1");
    assert_eq!(windows[0].interval.as_tuple(), (10, 15));
    assert_eq!(windows[0].blacklisted_fraction, Some(0.25));
    assert_eq!(windows[1].chrom, "chr2");
    assert_eq!(windows[1].interval.as_tuple(), (20, 25));
    assert_eq!(loaded.row_metadata().mode(), EndMotifRowMode::BedWindows);
    let cfdnalab::output_loaders::EndMotifRowMetadata::Windows { window_mode, .. } =
        loaded.row_metadata()
    else {
        panic!("expected window metadata");
    };
    assert_eq!(*window_mode, EndMotifWindowMode::Bed);

    let second_window = loaded.window(1)?.expect("second window should exist");
    assert_eq!(second_window.chrom, "chr2");

    let selected = loaded.select().windows(&[1, 0]).motifs(&[1]).read()?;
    assert_eq!(selected.storage_mode(), EndMotifStorageMode::Dense);
    assert_eq!(selected.row_indices(), &[1, 0]);
    assert_eq!(selected.motif_indices(), &[1]);
    assert_eq!(selected.motif_labels(), &["_GG".to_string()]);
    assert_eq!(selected.shape(), (2, 1));
    assert_eq!(selected.row_count(), 2);
    assert_eq!(selected.motif_count(), 1);
    assert_eq!(selected.count(0, 0), Some(3.5));
    assert_eq!(selected.dense_counts()?.values_row_major(), &[3.5, 2.0]);
    assert_eq!(
        selected
            .window_metadata()?
            .iter()
            .map(|window| (window.chrom.as_str(), window.interval.as_tuple()))
            .collect::<Vec<_>>(),
        vec![("chr2", (20, 25)), ("chr1", (10, 15))]
    );
    assert!(selected.group_metadata().is_err());

    let selected_by_label = loaded
        .select()
        .rows(&[0])
        .motifs_by_label(&["_AA"])
        .read()?;
    assert_eq!(selected_by_label.motif_indices(), &[0]);
    assert_eq!(selected_by_label.dense_counts()?.values_row_major(), &[1.0]);

    let duplicate_window_error = loaded
        .select()
        .windows(&[0, 0])
        .read()
        .expect_err("duplicate window selectors should fail");
    assert!(
        duplicate_window_error
            .to_string()
            .contains("duplicate value 0")
    );

    let duplicate_motif_error = loaded
        .select()
        .motifs_by_label(&["_AA", "_AA"])
        .read()
        .expect_err("duplicate motif-label selectors should fail");
    let wrong_mode_error = loaded
        .select()
        .groups(&[0])
        .read()
        .expect_err("windowed output should not provide group rows");
    let missing_motif_error = loaded
        .select()
        .motifs_by_label(&["_TT"])
        .read()
        .expect_err("missing motif label should fail");
    let bad_motif_index_error = loaded
        .select()
        .motifs(&[2])
        .read()
        .expect_err("motif index should be validated");
    let conflicting_row_selector_error = loaded
        .select()
        .rows(&[0])
        .windows(&[0])
        .read()
        .expect_err("conflicting row selectors should fail");
    let conflicting_motif_selector_error = loaded
        .select()
        .motifs(&[0])
        .motifs_by_label(&["_AA"])
        .read()
        .expect_err("conflicting motif selectors should fail");
    assert!(
        duplicate_motif_error
            .to_string()
            .contains("duplicate value '_AA'")
    );
    assert!(wrong_mode_error.to_string().contains("not grouped"));
    assert!(
        missing_motif_error
            .to_string()
            .contains("no motif label '_TT'")
    );
    assert!(
        bad_motif_index_error
            .to_string()
            .contains("motif index 2 is outside")
    );
    assert!(
        conflicting_row_selector_error
            .to_string()
            .contains("cannot combine rows() and windows() on the row axis")
    );
    assert!(
        conflicting_motif_selector_error
            .to_string()
            .contains("cannot combine motifs() and motifs_by_label() on the motif axis")
    );
    Ok(())
}

/// Verify dense global end-motif selections preserve global row metadata.
#[test]
fn load_ends_output_selects_dense_global_motifs_with_global_metadata() -> anyhow::Result<()> {
    // Arrange:
    // Global end-motif output has no window or group rows. A selection should
    // keep `Global` row metadata and still allow motif reordering by label.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "dense",
            "row_mode": "global",
            "motif_axis_kind": "motif",
            "count_units": "weighted_end_motif_count",
            "primary_array": "counts",
            "primary_group": null,
        }),
    )?;
    write_i32_array(&store, "motif_index", &[2], &["motif"], &[0, 1], json!({}))?;
    write_i32_array(
        &store,
        "motif_byte",
        &[3],
        &["motif_byte"],
        &[0, 1, 2],
        json!({}),
    )?;
    write_u8_array(
        &store,
        "motif_ascii",
        &[2, 3],
        &["motif", "motif_byte"],
        b"_AA_TT",
        json!({}),
    )?;
    write_i32_array(
        &store,
        "row",
        &[1],
        &["row"],
        &[0],
        json!({
            "label_field": "row_label",
            "labels": ["global"],
        }),
    )?;
    write_f64_array(
        &store,
        "counts",
        &[1, 2],
        &["row", "motif"],
        &[2.0, 5.0],
        json!({}),
    )?;
    let loaded = load_ends_output(&path)?;

    // Act
    let selected = loaded.select().motifs_by_label(&["_TT", "_AA"]).read()?;
    let row_error = loaded
        .select()
        .rows(&[0])
        .read()
        .expect_err("global output should not expose a selectable row axis");

    // Assert
    assert_eq!(selected.row_metadata(), &EndMotifRowMetadata::Global);
    assert_eq!(selected.row_indices(), &[0]);
    assert_eq!(selected.motif_indices(), &[1, 0]);
    assert_eq!(
        selected.motif_labels(),
        &["_TT".to_string(), "_AA".to_string()]
    );
    assert_eq!(selected.dense_counts()?.values_row_major(), &[5.0, 2.0]);
    assert!(
        row_error
            .to_string()
            .contains("global end-motif output has no selectable row axis")
    );
    assert!(selected.window_metadata().is_err());
    assert!(selected.group_metadata().is_err());
    Ok(())
}

/// Verify end-motif loaders reject malformed public labels from Zarr metadata.
#[test]
fn load_ends_output_rejects_invalid_public_labels() -> anyhow::Result<()> {
    // Arrange:
    // Concrete motifs are stored as fixed-width ASCII bytes. Motif groups use
    // JSON labels. Both become public selector strings and must be stable text.
    let temp = TempDir::new()?;
    let non_ascii_path = temp.path().join("non_ascii.end_motifs.zarr");
    let non_ascii_store = create_store(
        &non_ascii_path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "dense",
            "row_mode": "global",
            "motif_axis_kind": "motif",
            "count_units": "weighted_end_motif_count",
            "primary_array": "counts",
            "primary_group": null,
        }),
    )?;
    write_i32_array(
        &non_ascii_store,
        "motif_index",
        &[1],
        &["motif"],
        &[0],
        json!({}),
    )?;
    write_i32_array(
        &non_ascii_store,
        "motif_byte",
        &[2],
        &["motif_byte"],
        &[0, 1],
        json!({}),
    )?;
    write_u8_array(
        &non_ascii_store,
        "motif_ascii",
        &[1, 2],
        &["motif", "motif_byte"],
        b"\xC3\x85",
        json!({}),
    )?;
    write_i32_array(&non_ascii_store, "row", &[1], &["row"], &[0], json!({}))?;
    write_f64_array(
        &non_ascii_store,
        "counts",
        &[1, 1],
        &["row", "motif"],
        &[1.0],
        json!({}),
    )?;

    let control_label_path = temp.path().join("control_label.end_motifs.zarr");
    let control_label_store = create_store(
        &control_label_path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "dense",
            "row_mode": "global",
            "motif_axis_kind": "motif_group",
            "count_units": "weighted_end_motif_count",
            "primary_array": "counts",
            "primary_group": null,
        }),
    )?;
    write_i32_array(
        &control_label_store,
        "motif_index",
        &[1],
        &["motif"],
        &[0],
        json!({
            "label_field": "motif_group",
            "labels": ["bad\nlabel"],
        }),
    )?;
    write_i32_array(
        &control_label_store,
        "row",
        &[1],
        &["row"],
        &[0],
        json!({
            "label_field": "row_label",
            "labels": ["global"],
        }),
    )?;
    write_f64_array(
        &control_label_store,
        "counts",
        &[1, 1],
        &["row", "motif"],
        &[1.0],
        json!({}),
    )?;

    // Act
    let non_ascii_error =
        load_ends_output(&non_ascii_path).expect_err("non-ASCII motif bytes should fail");
    let control_label_error =
        load_ends_output(&control_label_path).expect_err("control-character label should fail");

    // Assert
    assert!(
        non_ascii_error
            .to_string()
            .contains("motif_ascii row 0 contains non-ASCII motif bytes")
    );
    assert!(
        control_label_error
            .to_string()
            .contains("Zarr label motif_group contains a control character")
    );
    Ok(())
}

/// Verify dense grouped motif-label selection keeps group metadata with counts.
#[test]
fn load_ends_output_selects_dense_grouped_labels_with_group_metadata() -> anyhow::Result<()> {
    // Arrange:
    // This follows the high-level grouped example path: select groups by name,
    // motifs by label, then pair selected group metadata with dense row values.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "dense",
            "row_mode": "grouped_bed",
            "motif_axis_kind": "motif_group",
            "count_units": "weighted_end_motif_count",
            "primary_array": "counts",
            "primary_group": null,
        }),
    )?;
    write_i32_array(
        &store,
        "motif_index",
        &[2],
        &["motif"],
        &[0, 1],
        json!({
            "label_field": "motif_group",
            "labels": ["left", "right"],
        }),
    )?;
    write_i32_array(&store, "row", &[2], &["row"], &[0, 1], json!({}))?;
    write_i32_array(
        &store,
        "group",
        &[2],
        &["row"],
        &[0, 1],
        json!({
            "label_field": "group_name",
            "labels": ["alpha", "beta"],
        }),
    )?;
    write_i32_array(
        &store,
        "eligible_windows",
        &[2],
        &["row"],
        &[1, 2],
        json!({}),
    )?;
    write_f64_array(
        &store,
        "blacklisted_fraction",
        &[2],
        &["row"],
        &[0.1, 0.2],
        json!({}),
    )?;
    write_f64_array(
        &store,
        "counts",
        &[2, 2],
        &["row", "motif"],
        &[1.0, 2.0, 3.0, 4.0],
        json!({}),
    )?;
    let loaded = load_ends_output(&path)?;

    // Act
    let selected = loaded
        .select()
        .groups_by_name(&["beta", "alpha"])
        .motifs_by_label(&["right", "left"])
        .read()?;
    let dense_counts = selected.to_dense_matrix()?;

    // Assert
    assert_eq!(selected.storage_mode(), EndMotifStorageMode::Dense);
    assert_eq!(selected.row_indices(), &[1, 0]);
    assert_eq!(selected.motif_indices(), &[1, 0]);
    assert_eq!(
        selected
            .group_metadata()?
            .iter()
            .map(|group| group.name.as_str())
            .collect::<Vec<_>>(),
        vec!["beta", "alpha"]
    );
    assert_eq!(dense_counts.values_row_major(), &[4.0, 3.0, 2.0, 1.0]);
    assert_eq!(
        selected
            .group_metadata()?
            .iter()
            .zip(dense_counts.rows())
            .map(|(group, motif_counts)| {
                (
                    group.name.as_str(),
                    motif_counts.iter().copied().sum::<f64>(),
                )
            })
            .collect::<Vec<_>>(),
        vec![("beta", 7.0), ("alpha", 3.0)]
    );
    Ok(())
}

/// Verify grouped end-motif outputs reject duplicate group names.
#[test]
fn load_ends_output_rejects_duplicate_group_names() -> anyhow::Result<()> {
    // Arrange:
    // Group names are public selectors. Duplicate labels would make
    // name-based lookup ambiguous, so the loader rejects them during load.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "dense",
            "row_mode": "grouped_bed",
            "motif_axis_kind": "motif_group",
            "count_units": "weighted_end_motif_count",
            "primary_array": "counts",
            "primary_group": null,
        }),
    )?;
    write_i32_array(
        &store,
        "motif_index",
        &[1],
        &["motif"],
        &[0],
        json!({
            "label_field": "motif_group",
            "labels": ["left"],
        }),
    )?;
    write_i32_array(&store, "row", &[2], &["row"], &[0, 1], json!({}))?;
    write_i32_array(
        &store,
        "group",
        &[2],
        &["row"],
        &[0, 1],
        json!({
            "label_field": "group_name",
            "labels": ["alpha", "alpha"],
        }),
    )?;
    write_i32_array(
        &store,
        "eligible_windows",
        &[2],
        &["row"],
        &[1, 1],
        json!({}),
    )?;
    write_f64_array(
        &store,
        "blacklisted_fraction",
        &[2],
        &["row"],
        &[0.0, 0.0],
        json!({}),
    )?;
    write_f64_array(
        &store,
        "counts",
        &[2, 1],
        &["row", "motif"],
        &[1.0, 2.0],
        json!({}),
    )?;

    // Act
    let error = load_ends_output(&path).expect_err("duplicate group names should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("group_names contains duplicate value 'alpha'")
    );
    Ok(())
}

/// Verify sparse grouped motif-group stores load COO counts and selections.
#[test]
fn load_ends_output_reads_sparse_grouped_motif_group_store() -> anyhow::Result<()> {
    // Arrange:
    // Sparse grouped output stores only observed row/motif pairs. Motif labels
    // are motif-group names in JSON attributes rather than motif_ascii rows.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "sparse_coo",
            "row_mode": "grouped_bed",
            "motif_axis_kind": "motif_group",
            "count_units": "weighted_end_motif_count",
            "primary_array": null,
            "primary_group": "sparse",
            "sparse_format": "coo",
            "sparse_indices_base": 0,
        }),
    )?;
    write_i32_array(
        &store,
        "motif_index",
        &[3],
        &["motif"],
        &[0, 1, 2],
        json!({
            "label_field": "motif_group",
            "labels": ["left", "right", "both"],
        }),
    )?;
    write_i32_array(&store, "row", &[2], &["row"], &[0, 1], json!({}))?;
    write_i32_array(
        &store,
        "group",
        &[2],
        &["row"],
        &[0, 1],
        json!({
            "label_field": "group_name",
            "labels": ["alpha", "beta"],
        }),
    )?;
    write_i32_array(
        &store,
        "eligible_windows",
        &[2],
        &["row"],
        &[2, 0],
        json!({}),
    )?;
    write_f64_array(
        &store,
        "blacklisted_fraction",
        &[2],
        &["row"],
        &[0.2, 0.0],
        json!({}),
    )?;
    write_group(&store, "/sparse", json!({}))?;
    write_i32_array(&store, "sparse/row", &[3], &["nnz"], &[0, 1, 1], json!({}))?;
    write_i32_array(
        &store,
        "sparse/motif",
        &[3],
        &["nnz"],
        &[2, 0, 2],
        json!({}),
    )?;
    write_f64_array(
        &store,
        "sparse/count",
        &[3],
        &["nnz"],
        &[4.0, 5.0, 6.0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "sparse/shape",
        &[2],
        &["sparse_dimension"],
        &[2, 3],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "sparse/sparse_dimension",
        &[2],
        &["sparse_dimension"],
        &[0, 1],
        json!({
            "label_field": "sparse_dimension_name",
            "labels": ["row", "motif"],
        }),
    )?;

    // Act
    let loaded = load_ends_output(&path)?;

    // Assert
    assert_eq!(loaded.storage_mode(), EndMotifStorageMode::SparseCoo);
    assert_eq!(loaded.row_mode(), EndMotifRowMode::Groups);
    assert_eq!(loaded.motif_axis_kind(), EndMotifAxisKind::MotifGroup);
    assert_eq!(
        loaded.motif_labels(),
        &["left".to_string(), "right".to_string(), "both".to_string()]
    );
    assert_eq!(loaded.count(0, 0), Some(0.0));
    assert_eq!(loaded.count(1, 2), Some(6.0));
    assert!(loaded.dense_counts().is_err());

    let groups = loaded.group_metadata()?;
    assert_eq!(groups[0].index, 0);
    assert_eq!(groups[0].name, "alpha");
    assert_eq!(groups[0].eligible_windows, 2);
    assert_eq!(groups[0].blacklisted_fraction, 0.2);
    assert_eq!(groups[1].name, "beta");
    assert_eq!(groups[1].eligible_windows, 0);
    assert_eq!(loaded.group_index("beta")?, 1);
    assert!(loaded.has_group("alpha"));
    assert_eq!(
        loaded.group(0)?.expect("first group should exist").name,
        "alpha"
    );

    let sparse = loaded.sparse_counts()?;
    assert_eq!(sparse.shape(), (2, 3));
    assert_eq!(sparse.nnz(), 3);
    assert_eq!(sparse.row_indices(), &[0, 1, 1]);
    assert_eq!(sparse.motif_indices(), &[2, 0, 2]);
    assert_eq!(sparse.counts(), &[4.0, 5.0, 6.0]);
    assert_eq!(
        sparse.entries().collect::<Vec<_>>(),
        vec![
            EndMotifSparseEntry {
                row_index: 0,
                motif_index: 2,
                count: 4.0,
            },
            EndMotifSparseEntry {
                row_index: 1,
                motif_index: 0,
                count: 5.0,
            },
            EndMotifSparseEntry {
                row_index: 1,
                motif_index: 2,
                count: 6.0,
            },
        ]
    );
    let sparse_lookup = sparse.to_lookup_index();
    assert_eq!(sparse_lookup.shape(), (2, 3));
    assert_eq!(sparse_lookup.count(1, 2), Some(6.0));
    assert_eq!(sparse_lookup.count(0, 1), Some(0.0));
    assert_eq!(sparse_lookup.count(2, 0), None);
    assert_eq!(
        sparse.to_dense_matrix()?.values_row_major(),
        &[0.0, 0.0, 4.0, 5.0, 0.0, 6.0]
    );

    let selected = loaded
        .select()
        .groups_by_name(&["beta", "alpha"])
        .motifs(&[2, 0])
        .read()?;
    assert_eq!(selected.storage_mode(), EndMotifStorageMode::SparseCoo);
    assert_eq!(selected.row_indices(), &[1, 0]);
    assert_eq!(selected.motif_indices(), &[2, 0]);
    assert_eq!(selected.row_count(), 2);
    assert_eq!(selected.motif_count(), 2);
    assert_eq!(
        selected.motif_labels(),
        &["both".to_string(), "left".to_string()]
    );
    assert_eq!(selected.count(0, 0), Some(6.0));
    assert_eq!(selected.count(1, 1), Some(0.0));
    assert_eq!(
        selected
            .group_metadata()?
            .iter()
            .map(|group| group.name.as_str())
            .collect::<Vec<_>>(),
        vec!["beta", "alpha"]
    );

    let selected_sparse = selected.sparse_counts()?;
    assert_eq!(selected_sparse.shape(), (2, 2));
    assert_eq!(selected_sparse.row_indices(), &[0, 0, 1]);
    assert_eq!(selected_sparse.motif_indices(), &[0, 1, 0]);
    assert_eq!(selected_sparse.counts(), &[6.0, 5.0, 4.0]);
    let selected_sparse_lookup = selected_sparse.to_lookup_index();
    assert_eq!(selected_sparse_lookup.count(0, 1), Some(5.0));
    assert_eq!(selected_sparse_lookup.count(1, 1), Some(0.0));
    assert_eq!(selected_sparse_lookup.count(2, 0), None);
    assert_eq!(
        selected.to_dense_matrix()?.values_row_major(),
        &[6.0, 5.0, 4.0, 0.0]
    );
    let dense_selection = selected.to_dense_matrix()?;
    assert_eq!(
        selected
            .group_metadata()?
            .iter()
            .zip(dense_selection.rows())
            .map(|(group, motif_counts)| {
                (
                    group.name.as_str(),
                    motif_counts.iter().copied().sum::<f64>(),
                )
            })
            .collect::<Vec<_>>(),
        vec![("beta", 11.0), ("alpha", 4.0)]
    );
    Ok(())
}

/// Verify sparse stores with no observed motifs load as empty-width matrices.
#[test]
fn load_ends_output_reads_sparse_global_store_with_no_motifs() -> anyhow::Result<()> {
    // Arrange:
    // The writer allows sparse global output with a valid row axis and an
    // empty motif axis when no motifs were observed.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "sparse_coo",
            "row_mode": "global",
            "motif_axis_kind": "motif",
            "count_units": "weighted_end_motif_count",
            "primary_array": null,
            "primary_group": "sparse",
            "sparse_format": "coo",
            "sparse_indices_base": 0,
        }),
    )?;
    write_i32_array(&store, "motif_index", &[0], &["motif"], &[], json!({}))?;
    write_i32_array(&store, "motif_byte", &[0], &["motif_byte"], &[], json!({}))?;
    write_u8_array(
        &store,
        "motif_ascii",
        &[0, 0],
        &["motif", "motif_byte"],
        b"",
        json!({}),
    )?;
    write_i32_array(
        &store,
        "row",
        &[1],
        &["row"],
        &[0],
        json!({
            "label_field": "row_label",
            "labels": ["global"],
        }),
    )?;
    write_group(&store, "/sparse", json!({}))?;
    write_i32_array(&store, "sparse/row", &[0], &["nnz"], &[], json!({}))?;
    write_i32_array(&store, "sparse/motif", &[0], &["nnz"], &[], json!({}))?;
    write_f64_array(&store, "sparse/count", &[0], &["nnz"], &[], json!({}))?;
    write_i32_array(
        &store,
        "sparse/shape",
        &[2],
        &["sparse_dimension"],
        &[1, 0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "sparse/sparse_dimension",
        &[2],
        &["sparse_dimension"],
        &[0, 1],
        json!({
            "label_field": "sparse_dimension_name",
            "labels": ["row", "motif"],
        }),
    )?;

    // Act
    let loaded = load_ends_output(&path)?;
    let selected = loaded.select().read()?;

    // Assert
    assert_eq!(loaded.storage_mode(), EndMotifStorageMode::SparseCoo);
    assert_eq!(loaded.row_mode(), EndMotifRowMode::Global);
    assert_eq!(loaded.motif_labels(), &Vec::<String>::new());
    assert_eq!(loaded.sparse_counts()?.shape(), (1, 0));
    assert_eq!(loaded.sparse_counts()?.nnz(), 0);
    assert_eq!(loaded.sparse_counts()?.to_dense_matrix()?.shape(), (1, 0));
    assert_eq!(selected.shape(), (1, 0));
    assert_eq!(selected.motif_indices(), &[] as &[usize]);
    assert_eq!(selected.sparse_counts()?.nnz(), 0);
    assert_eq!(
        selected.to_dense_matrix()?.values_row_major(),
        &[] as &[f64]
    );
    assert!(matches!(
        selected.row_metadata(),
        EndMotifRowMetadata::Global
    ));
    Ok(())
}

/// Verify unsupported end-motif schema versions are rejected.
#[test]
fn load_ends_output_rejects_wrong_schema_version() -> anyhow::Result<()> {
    // Arrange:
    // The Rust loader targets the current writer schema and reports a clear
    // error when the root metadata advertises another version.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 1,
            "storage_mode": "dense",
            "row_mode": "global",
            "motif_axis_kind": "motif",
            "count_units": "weighted_end_motif_count",
            "primary_array": "counts",
            "primary_group": null,
        }),
    )?;

    // Act
    let error = load_ends_output(&path).expect_err("schema version mismatch should fail");

    // Assert
    assert!(error.to_string().contains("schema version mismatch"));
    Ok(())
}

/// Verify dense count arrays must match row and motif axis sizes.
#[test]
fn load_ends_output_rejects_dense_count_shape_mismatch() -> anyhow::Result<()> {
    // Arrange:
    // Dense count arrays must match the row and motif axes. A mismatched shape
    // would make every row/motif lookup ambiguous.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "dense",
            "row_mode": "global",
            "motif_axis_kind": "motif",
            "count_units": "weighted_end_motif_count",
            "primary_array": "counts",
            "primary_group": null,
        }),
    )?;
    write_i32_array(&store, "motif_index", &[1], &["motif"], &[0], json!({}))?;
    write_i32_array(
        &store,
        "motif_byte",
        &[3],
        &["motif_byte"],
        &[0, 1, 2],
        json!({}),
    )?;
    write_u8_array(
        &store,
        "motif_ascii",
        &[1, 3],
        &["motif", "motif_byte"],
        b"_AA",
        json!({}),
    )?;
    write_i32_array(
        &store,
        "row",
        &[1],
        &["row"],
        &[0],
        json!({
            "label_field": "row_label",
            "labels": ["global"],
        }),
    )?;
    write_f64_array(
        &store,
        "counts",
        &[1, 2],
        &["row", "motif"],
        &[1.0, 2.0],
        json!({}),
    )?;

    // Act
    let error = load_ends_output(&path).expect_err("dense shape mismatch should fail");

    // Assert
    assert!(error.to_string().contains("dense end-motif counts shape"));
    Ok(())
}

/// Verify dense count arrays reject malformed numeric values.
#[test]
fn load_ends_output_rejects_non_finite_dense_counts() -> anyhow::Result<()> {
    // Arrange:
    // Dense end-motif counts are non-negative weighted counts. A non-finite
    // value should fail during load, before users can select or sum it.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "dense",
            "row_mode": "global",
            "motif_axis_kind": "motif",
            "count_units": "weighted_end_motif_count",
            "primary_array": "counts",
            "primary_group": null,
        }),
    )?;
    write_i32_array(&store, "motif_index", &[1], &["motif"], &[0], json!({}))?;
    write_i32_array(
        &store,
        "motif_byte",
        &[3],
        &["motif_byte"],
        &[0, 1, 2],
        json!({}),
    )?;
    write_u8_array(
        &store,
        "motif_ascii",
        &[1, 3],
        &["motif", "motif_byte"],
        b"_AA",
        json!({}),
    )?;
    write_i32_array(
        &store,
        "row",
        &[1],
        &["row"],
        &[0],
        json!({
            "label_field": "row_label",
            "labels": ["global"],
        }),
    )?;
    write_f64_array(
        &store,
        "counts",
        &[1, 1],
        &["row", "motif"],
        &[f64::NAN],
        json!({}),
    )?;

    // Act
    let error = load_ends_output(&path).expect_err("non-finite dense count should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("dense end-motif counts contain value outside finite and non-negative range")
    );
    Ok(())
}

/// Verify sparse COO count arrays reject malformed numeric values.
#[test]
fn load_ends_output_rejects_negative_sparse_counts() -> anyhow::Result<()> {
    // Arrange:
    // Sparse COO stores omit zeroes, but stored values are still non-negative
    // weighted counts and should be validated before lookup structures use
    // them.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "sparse_coo",
            "row_mode": "global",
            "motif_axis_kind": "motif",
            "count_units": "weighted_end_motif_count",
            "primary_array": null,
            "primary_group": "sparse",
            "sparse_format": "coo",
            "sparse_indices_base": 0,
        }),
    )?;
    write_i32_array(&store, "motif_index", &[1], &["motif"], &[0], json!({}))?;
    write_i32_array(
        &store,
        "motif_byte",
        &[3],
        &["motif_byte"],
        &[0, 1, 2],
        json!({}),
    )?;
    write_u8_array(
        &store,
        "motif_ascii",
        &[1, 3],
        &["motif", "motif_byte"],
        b"_AA",
        json!({}),
    )?;
    write_i32_array(
        &store,
        "row",
        &[1],
        &["row"],
        &[0],
        json!({
            "label_field": "row_label",
            "labels": ["global"],
        }),
    )?;
    write_group(&store, "/sparse", json!({}))?;
    write_i32_array(&store, "sparse/row", &[1], &["nnz"], &[0], json!({}))?;
    write_i32_array(&store, "sparse/motif", &[1], &["nnz"], &[0], json!({}))?;
    write_f64_array(&store, "sparse/count", &[1], &["nnz"], &[-1.0], json!({}))?;
    write_i32_array(
        &store,
        "sparse/shape",
        &[2],
        &["sparse_dimension"],
        &[1, 1],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "sparse/sparse_dimension",
        &[2],
        &["sparse_dimension"],
        &[0, 1],
        json!({
            "label_field": "sparse_dimension_name",
            "labels": ["row", "motif"],
        }),
    )?;

    // Act
    let error = load_ends_output(&path).expect_err("negative sparse count should fail");

    // Assert
    assert!(
        error.to_string().contains(
            "sparse end-motif counts contain value outside finite and non-negative range"
        )
    );
    Ok(())
}

/// Verify sparse COO coordinates must be sorted and unique.
#[test]
fn load_ends_output_rejects_unsorted_sparse_coordinates() -> anyhow::Result<()> {
    // Arrange:
    // Sparse COO coordinates are binary-searched by the public count accessors,
    // so they must be sorted and unique in source storage.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.end_motifs.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "end_motif_counts",
            "cfdnalab_schema_version": 2,
            "storage_mode": "sparse_coo",
            "row_mode": "grouped_bed",
            "motif_axis_kind": "motif_group",
            "count_units": "weighted_end_motif_count",
            "primary_array": null,
            "primary_group": "sparse",
            "sparse_format": "coo",
            "sparse_indices_base": 0,
        }),
    )?;
    write_i32_array(
        &store,
        "motif_index",
        &[2],
        &["motif"],
        &[0, 1],
        json!({
            "label_field": "motif_group",
            "labels": ["left", "right"],
        }),
    )?;
    write_i32_array(&store, "row", &[2], &["row"], &[0, 1], json!({}))?;
    write_i32_array(
        &store,
        "group",
        &[2],
        &["row"],
        &[0, 1],
        json!({
            "label_field": "group_name",
            "labels": ["alpha", "beta"],
        }),
    )?;
    write_i32_array(
        &store,
        "eligible_windows",
        &[2],
        &["row"],
        &[1, 1],
        json!({}),
    )?;
    write_f64_array(
        &store,
        "blacklisted_fraction",
        &[2],
        &["row"],
        &[0.0, 0.0],
        json!({}),
    )?;
    write_group(&store, "/sparse", json!({}))?;
    write_i32_array(&store, "sparse/row", &[2], &["nnz"], &[1, 0], json!({}))?;
    write_i32_array(&store, "sparse/motif", &[2], &["nnz"], &[0, 1], json!({}))?;
    write_f64_array(
        &store,
        "sparse/count",
        &[2],
        &["nnz"],
        &[2.0, 3.0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "sparse/shape",
        &[2],
        &["sparse_dimension"],
        &[2, 2],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "sparse/sparse_dimension",
        &[2],
        &["sparse_dimension"],
        &[0, 1],
        json!({
            "label_field": "sparse_dimension_name",
            "labels": ["row", "motif"],
        }),
    )?;

    // Act
    let error = load_ends_output(&path).expect_err("unsorted sparse coordinates should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("COO entries must be sorted and unique")
    );
    Ok(())
}

/// Create a temporary Zarr store with root attributes.
fn create_store(path: &Path, attributes: Value) -> anyhow::Result<Arc<FilesystemStore>> {
    fs::create_dir_all(path)?;
    let store = Arc::new(FilesystemStore::new(path)?);
    write_group(&store, "/", attributes)?;
    Ok(store)
}

/// Write one Zarr group metadata object.
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

/// Write an `i32` Zarr array fixture.
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

/// Write an `i64` Zarr array fixture.
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

/// Write a `u8` Zarr array fixture.
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

/// Write an `f64` Zarr array fixture.
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

/// Write a typed Zarr array fixture with one chunk.
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

/// Convert a JSON value into an object for Zarr attributes.
fn json_object(value: Value) -> anyhow::Result<Map<String, Value>> {
    match value {
        Value::Object(map) => Ok(map),
        other => anyhow::bail!("expected object attributes, got {other}"),
    }
}
