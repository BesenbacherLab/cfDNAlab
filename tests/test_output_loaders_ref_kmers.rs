#![cfg(feature = "cmd_ref_kmers")]
//! Public API tests for Rust output loaders for `cfdna ref-kmers`.
//!
//! These tests build tiny Zarr V3 stores with the public reference k-mer schema and assert that
//! the loader returns typed Rust metadata, native frequency containers, and reconstructed counts.

use cfdnalab::output_loaders::{
    RefKmerMotifAxisKind, RefKmerRowMetadata, RefKmerRowMode, RefKmerSparseCountEntry,
    RefKmerSparseFrequencyEntry, RefKmerStorageMode, load_ref_kmers_output,
};
use cfdnalab::run_like_cli::ref_kmers::RefKmerOrientation;
use serde_json::{Map, Value, json};
use std::{fs, path::Path, sync::Arc};
use tempfile::TempDir;
use zarrs::{
    array::{ArrayBuilder, DataType, Element, builder::ArrayBuilderFillValue, data_type},
    filesystem::FilesystemStore,
    group::GroupBuilder,
};

/// Verify dense global reference k-mer stores load metadata, frequencies, and counts.
#[test]
fn load_ref_kmers_output_reads_dense_global_store() -> anyhow::Result<()> {
    // Arrange:
    // Dense global output stores one row of frequencies. The row scaling factor
    // reconstructs counts from the frequency values.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.ref_kmers.zarr");
    let store = create_store(&path, dense_root_attributes(1, false))?;
    write_concrete_motif_axis(&store, &["A", "C"])?;
    write_global_row_axis(&store)?;
    write_row_scaling_factors(&store, &[4.0])?;
    write_reference_footprint(&store)?;
    write_f64_array(
        &store,
        "frequencies",
        &[1, 2],
        &["row", "motif"],
        &[0.75, 0.25],
        json!({}),
    )?;

    // Act
    let loaded = load_ref_kmers_output(&path)?;

    // Assert
    assert_eq!(loaded.storage_mode(), RefKmerStorageMode::Dense);
    assert_eq!(loaded.row_mode(), RefKmerRowMode::Global);
    assert_eq!(loaded.row_metadata(), &RefKmerRowMetadata::Global);
    assert_eq!(loaded.motif_axis_kind(), RefKmerMotifAxisKind::Motif);
    assert_eq!(loaded.kmer_size(), 1);
    assert!(!loaded.canonical());
    assert_eq!(loaded.orientation(), RefKmerOrientation::Both);
    assert!(!loaded.all_motifs());
    assert_eq!(loaded.assign_by(), "count-overlap");
    assert_eq!(loaded.motif_labels(), &["A".to_string(), "C".to_string()]);
    assert_eq!(loaded.motif_index("C")?, 1);
    assert_eq!(loaded.row_scaling_factors(), &[4.0]);
    assert_eq!(loaded.dense_frequencies()?.shape(), (1, 2));
    assert_eq!(
        loaded.dense_frequencies()?.values_row_major(),
        &[0.75, 0.25]
    );
    assert_eq!(loaded.frequency(0, 0), Some(0.75));
    assert_eq!(loaded.frequency_for_motif(0, "C")?, Some(0.25));
    assert_eq!(loaded.count(0, 0), Some(3.0));
    assert_eq!(loaded.count_for_motif(0, "C")?, Some(1.0));
    assert!(loaded.sparse_frequencies().is_err());
    assert_eq!(
        loaded.output_metadata().to_string(),
        "storage_mode=dense, row_mode=global, motif_axis=motifs, row_count=1, motif_count=2, kmer_size=1, canonical=false, orientation=both, all_motifs=false, assign_by=count-overlap"
    );
    Ok(())
}

/// Verify sparse stores expose sorted COO frequencies and a reusable lookup index.
#[test]
fn load_ref_kmers_output_exposes_sparse_lookup_index() -> anyhow::Result<()> {
    // Arrange:
    // Sparse output stores only non-zero frequency coordinates. Missing cells
    // remain implicit zeroes for point lookup and dense reconstruction.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.ref_kmers.zarr");
    let store = create_store(&path, sparse_root_attributes(1, false))?;
    write_concrete_motif_axis(&store, &["A", "C", "G"])?;
    write_global_row_axis(&store)?;
    write_row_scaling_factors(&store, &[2.0])?;
    write_reference_footprint(&store)?;
    write_sparse_metadata(&store, 1, 3, &[0, 0], &[0, 2], &[0.25, 0.75])?;

    // Act
    let loaded = load_ref_kmers_output(&path)?;

    // Assert
    assert_eq!(loaded.storage_mode(), RefKmerStorageMode::SparseCoo);
    assert_eq!(loaded.row_mode(), RefKmerRowMode::Global);
    assert_eq!(loaded.motif_axis_kind(), RefKmerMotifAxisKind::Motif);
    assert_eq!(loaded.frequency(0, 1), Some(0.0));
    assert_eq!(loaded.count(0, 2), Some(1.5));
    assert!(loaded.dense_frequencies().is_err());

    let sparse = loaded.sparse_frequencies()?;
    assert_eq!(sparse.shape(), (1, 3));
    assert_eq!(sparse.nnz(), 2);
    assert_eq!(sparse.row_indices(), &[0, 0]);
    assert_eq!(sparse.motif_indices(), &[0, 2]);
    assert_eq!(sparse.frequencies(), &[0.25, 0.75]);
    assert_eq!(
        sparse.entries().collect::<Vec<_>>(),
        vec![
            RefKmerSparseFrequencyEntry {
                row_index: 0,
                motif_index: 0,
                frequency: 0.25,
            },
            RefKmerSparseFrequencyEntry {
                row_index: 0,
                motif_index: 2,
                frequency: 0.75,
            },
        ]
    );
    let sparse_lookup = sparse.to_lookup_index();
    assert_eq!(sparse_lookup.shape(), (1, 3));
    assert_eq!(sparse_lookup.frequency(0, 2), Some(0.75));
    assert_eq!(sparse_lookup.frequency(0, 1), Some(0.0));
    assert_eq!(sparse_lookup.frequency(1, 0), None);
    assert_eq!(
        loaded.sparse_count_entries()?,
        vec![
            RefKmerSparseCountEntry {
                row_index: 0,
                motif_index: 0,
                count: 0.5,
            },
            RefKmerSparseCountEntry {
                row_index: 0,
                motif_index: 2,
                count: 1.5,
            },
        ]
    );
    assert_eq!(
        loaded.to_dense_frequency_matrix()?.values_row_major(),
        &[0.25, 0.0, 0.75]
    );
    assert_eq!(
        loaded.to_dense_count_matrix()?.values_row_major(),
        &[0.5, 0.0, 1.5]
    );
    Ok(())
}

/// Verify unsupported reference k-mer schema versions are rejected.
#[test]
fn load_ref_kmers_output_rejects_wrong_schema_version() -> anyhow::Result<()> {
    // Arrange:
    // The Rust loader targets the current writer schema and reports a clear
    // error when the root metadata advertises another version.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.ref_kmers.zarr");
    let mut attributes = dense_root_attributes(1, false);
    attributes["cfdnalab_schema_version"] = json!(1);
    create_store(&path, attributes)?;

    // Act
    let error = load_ref_kmers_output(&path).expect_err("schema version mismatch should fail");

    // Assert
    assert!(error.to_string().contains("schema version mismatch"));
    Ok(())
}

/// Verify unknown orientation metadata is rejected.
#[test]
fn load_ref_kmers_output_rejects_unknown_orientation() -> anyhow::Result<()> {
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.ref_kmers.zarr");
    let mut attributes = dense_root_attributes(1, false);
    attributes["orientation"] = json!("unknown");
    create_store(&path, attributes)?;

    let error = load_ref_kmers_output(&path).expect_err("unknown orientation should fail");

    assert!(
        error
            .to_string()
            .contains("unsupported reference k-mer orientation")
    );
    Ok(())
}

/// Verify count reconstruction metadata must match the current writer contract.
#[test]
fn load_ref_kmers_output_rejects_count_reconstruction_mismatch() -> anyhow::Result<()> {
    // Arrange:
    // This metadata defines how stored frequencies reconstruct counts. A
    // changed formula would make loader-side count methods wrong.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.ref_kmers.zarr");
    let mut attributes = dense_root_attributes(1, false);
    attributes["count_reconstruction"] = json!("count = frequency");
    create_store(&path, attributes)?;

    // Act
    let error =
        load_ref_kmers_output(&path).expect_err("count reconstruction mismatch should fail");

    // Assert
    assert!(error.to_string().contains("count_reconstruction"));
    Ok(())
}

/// Verify concrete reference k-mer labels must match k-mer metadata.
#[test]
fn load_ref_kmers_output_rejects_invalid_concrete_motif_labels() -> anyhow::Result<()> {
    // Arrange:
    // Concrete labels are interpreted as reference k-mers, so they must use
    // A/C/G/T bases and match the canonical setting recorded in root metadata.
    let temp = TempDir::new()?;
    let invalid_base_path = temp.path().join("invalid_base.ref_kmers.zarr");
    let invalid_base_store = create_store(&invalid_base_path, dense_root_attributes(2, false))?;
    write_i32_array(
        &invalid_base_store,
        "motif_index",
        &[1],
        &["motif"],
        &[0],
        json!({}),
    )?;
    write_i32_array(
        &invalid_base_store,
        "motif_byte",
        &[2],
        &["motif_byte"],
        &[0, 1],
        json!({}),
    )?;
    write_u8_array(
        &invalid_base_store,
        "motif_ascii",
        &[1, 2],
        &["motif", "motif_byte"],
        b"AN",
        json!({}),
    )?;

    let noncanonical_path = temp.path().join("noncanonical.ref_kmers.zarr");
    let noncanonical_store = create_store(&noncanonical_path, dense_root_attributes(2, true))?;
    write_concrete_motif_axis(&noncanonical_store, &["TG"])?;

    // Act
    let invalid_base_error = load_ref_kmers_output(&invalid_base_path)
        .expect_err("invalid concrete motif base should fail");
    let noncanonical_error = load_ref_kmers_output(&noncanonical_path)
        .expect_err("noncanonical concrete motif should fail");

    // Assert
    assert!(
        invalid_base_error
            .to_string()
            .contains("contains invalid base `N`")
    );
    assert!(
        noncanonical_error
            .to_string()
            .contains("should be represented as `CA`")
    );
    Ok(())
}

/// Verify dense frequency arrays must match row and motif axis sizes.
#[test]
fn load_ref_kmers_output_rejects_dense_frequency_shape_mismatch() -> anyhow::Result<()> {
    // Arrange:
    // Dense frequency arrays must match the row and motif axes. A mismatched
    // shape would make row/motif lookup ambiguous.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.ref_kmers.zarr");
    let store = create_store(&path, dense_root_attributes(1, false))?;
    write_concrete_motif_axis(&store, &["A"])?;
    write_global_row_axis(&store)?;
    write_row_scaling_factors(&store, &[1.0])?;
    write_reference_footprint(&store)?;
    write_f64_array(
        &store,
        "frequencies",
        &[1, 2],
        &["row", "motif"],
        &[0.5, 0.5],
        json!({}),
    )?;

    // Act
    let error = load_ref_kmers_output(&path).expect_err("dense shape mismatch should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("dense reference k-mer frequencies shape")
    );
    Ok(())
}

/// Verify frequency values and row scaling factors reject malformed numeric values.
#[test]
fn load_ref_kmers_output_rejects_bad_frequency_values_and_scaling_factors() -> anyhow::Result<()> {
    // Arrange:
    // Frequencies must be finite fractions, and row scaling factors must be
    // finite non-negative count totals used for count reconstruction.
    let temp = TempDir::new()?;
    let bad_frequency_path = temp.path().join("bad_frequency.ref_kmers.zarr");
    let bad_frequency_store = create_store(&bad_frequency_path, dense_root_attributes(1, false))?;
    write_concrete_motif_axis(&bad_frequency_store, &["A"])?;
    write_global_row_axis(&bad_frequency_store)?;
    write_row_scaling_factors(&bad_frequency_store, &[1.0])?;
    write_reference_footprint(&bad_frequency_store)?;
    write_f64_array(
        &bad_frequency_store,
        "frequencies",
        &[1, 1],
        &["row", "motif"],
        &[f64::NAN],
        json!({}),
    )?;

    let bad_scaling_path = temp.path().join("bad_scaling.ref_kmers.zarr");
    let bad_scaling_store = create_store(&bad_scaling_path, dense_root_attributes(1, false))?;
    write_concrete_motif_axis(&bad_scaling_store, &["A"])?;
    write_global_row_axis(&bad_scaling_store)?;
    write_row_scaling_factors(&bad_scaling_store, &[-1.0])?;

    // Act
    let bad_frequency_error = load_ref_kmers_output(&bad_frequency_path)
        .expect_err("non-finite dense frequency should fail");
    let bad_scaling_error = load_ref_kmers_output(&bad_scaling_path)
        .expect_err("negative row scaling factor should fail");

    // Assert
    assert!(
        bad_frequency_error
            .to_string()
            .contains("dense reference k-mer frequencies contain value outside [0, 1]")
    );
    assert!(
        bad_scaling_error
            .to_string()
            .contains("row_scaling_factor contains value outside finite and non-negative range")
    );
    Ok(())
}

/// Verify sparse COO coordinates must be sorted and unique.
#[test]
fn load_ref_kmers_output_rejects_unsorted_sparse_coordinates() -> anyhow::Result<()> {
    // Arrange:
    // Sparse COO coordinates are binary-searched by the public frequency
    // accessors, so they must be sorted and unique in source storage.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.ref_kmers.zarr");
    let store = create_store(&path, sparse_root_attributes(1, false))?;
    write_concrete_motif_axis(&store, &["A", "C"])?;
    write_global_row_axis(&store)?;
    write_row_scaling_factors(&store, &[2.0])?;
    write_reference_footprint(&store)?;
    write_sparse_metadata(&store, 1, 2, &[0, 0], &[1, 0], &[0.5, 0.5])?;

    // Act
    let error = load_ref_kmers_output(&path).expect_err("unsorted sparse coordinates should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("COO entries must be sorted and unique")
    );
    Ok(())
}

/// Root metadata for a dense reference k-mer Zarr fixture.
fn dense_root_attributes(kmer_size: u8, canonical: bool) -> Value {
    let mut attributes = common_root_attributes(kmer_size, canonical);
    attributes["storage_mode"] = json!("dense");
    attributes["primary_array"] = json!("frequencies");
    attributes["primary_group"] = Value::Null;
    attributes
}

/// Root metadata for a sparse reference k-mer Zarr fixture.
fn sparse_root_attributes(kmer_size: u8, canonical: bool) -> Value {
    let mut attributes = common_root_attributes(kmer_size, canonical);
    attributes["storage_mode"] = json!("sparse_coo");
    attributes["primary_array"] = Value::Null;
    attributes["primary_group"] = json!("sparse");
    attributes["sparse_format"] = json!("coo");
    attributes["sparse_indices_base"] = json!(0);
    attributes
}

/// Root metadata shared by dense and sparse reference k-mer fixtures.
fn common_root_attributes(kmer_size: u8, canonical: bool) -> Value {
    json!({
        "cfdnalab_schema": "ref_kmer_frequencies",
        "cfdnalab_schema_version": 2,
        "row_mode": "global",
        "motif_axis_kind": "motif",
        "value_units": "reference_kmer_frequency",
        "count_units": "reference_kmer_count",
        "row_scaling_factor_array": "row_scaling_factor",
        "count_reconstruction": "reference_kmer_count = frequency * row_scaling_factor[row]",
        "kmer_size": kmer_size,
        "canonical": canonical,
        "orientation": "both",
        "all_motifs": false,
        "assign_by": "count-overlap",
    })
}

/// Create a temporary Zarr store with root attributes.
fn create_store(path: &Path, attributes: Value) -> anyhow::Result<Arc<FilesystemStore>> {
    fs::create_dir_all(path)?;
    let store = Arc::new(FilesystemStore::new(path)?);
    write_group(&store, "/", attributes)?;
    Ok(store)
}

/// Write a concrete reference k-mer motif axis.
fn write_concrete_motif_axis(store: &Arc<FilesystemStore>, labels: &[&str]) -> anyhow::Result<()> {
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

/// Write a one-row global row axis.
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

/// Write row scaling factors for frequency-to-count reconstruction.
fn write_row_scaling_factors(store: &Arc<FilesystemStore>, values: &[f64]) -> anyhow::Result<()> {
    write_f64_array(
        store,
        "row_scaling_factor",
        &[values.len()],
        &["row"],
        values,
        json!({}),
    )
}

/// Write an empty reference contig footprint JSON array.
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

/// Write sparse COO metadata and arrays.
fn write_sparse_metadata(
    store: &Arc<FilesystemStore>,
    row_count: i32,
    motif_count: i32,
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
    write_i32_array(
        store,
        "sparse/shape",
        &[2],
        &["sparse_dimension"],
        &[row_count, motif_count],
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
    )?;
    Ok(())
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
