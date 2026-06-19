//! Shared Zarr I/O utilities.
//!
//! These helpers are deliberately low-level. They centralize the mechanical parts that should be
//! identical across cfDNAlab Zarr outputs:
//!
//! - opening a filesystem-backed Zarr store
//! - creating V3 arrays with zstd compression and native dimension names
//! - writing small coordinate/metadata arrays as a single chunk
//! - reading complete coordinate/metadata arrays
//! - validating cfDNAlab root schema attributes
//! - validating JSON attributes and public integer dtype narrowing
//! - validating public labels stored in JSON attributes
//!
//! Command modules still own their public schemas. This module should not know what a midpoint
//! profile, end motif, group, window, or length bin means.

use anyhow::{Context, Result, ensure};
use serde_json::{Map, Value};
use std::{fs, path::Path, sync::Arc};
use zarrs::{
    array::{
        ArrayBuilder, DataType, Element, FillValueMetadata, builder::ArrayBuilderFillValue,
        codec::ZstdCodec,
    },
    filesystem::FilesystemStore,
    group::GroupBuilder,
};

/// Default zstd level used by cfDNAlab Zarr stores.
///
/// Keeping this in one place avoids midpoint and ends drifting to different compression settings
/// unless downstream compatibility or performance tests give a concrete reason.
pub(crate) const DEFAULT_ZARR_ZSTD_LEVEL: i32 = 3;

/// Fill value for public `int32` coordinate and metadata arrays.
///
/// Several downstream readers treat a Zarr array's fill value as missing data. Public cfDNAlab
/// coordinate axes are zero-based, so using zero would turn valid index 0 values into missing
/// values in those readers. `-1` is outside the valid domain for these arrays and should only be
/// seen as chunk padding or metadata for empty arrays, not as real cfDNAlab data.
#[cfg(any(feature = "cmd_ends", feature = "cmd_midpoints"))]
pub(crate) const ZARR_INT32_FILL_VALUE: i32 = -1;

/// Fill value for public `int64` genomic coordinate arrays.
#[cfg(any(feature = "cmd_ends"))]
pub(crate) const ZARR_INT64_FILL_VALUE: i64 = -1;

/// Fill value for non-negative `float32` count arrays.
#[cfg(any(feature = "cmd_midpoints"))]
pub(crate) const ZARR_FLOAT32_FILL_VALUE: f32 = -1.0;

/// Fill value for non-negative `float64` count and fraction arrays.
#[cfg(any(feature = "cmd_ends"))]
pub(crate) const ZARR_FLOAT64_FILL_VALUE: f64 = -1.0;

/// Fill value for fixed-width ASCII label arrays.
///
/// Valid ASCII labels only use byte values `0..=127`, so `255` cannot be confused with a real
/// label byte. Do not reuse this for arbitrary numeric `uint8` arrays, where `255` may be a valid
/// data value.
#[cfg(any(feature = "cmd_ends"))]
pub(crate) const ZARR_ASCII_FILL_VALUE: u8 = u8::MAX;

/// Open or create a filesystem-backed Zarr store directory.
///
/// This does not remove or replace existing final outputs. Commands write into a temporary store
/// and rely on `FinalOutputFiles` to move the completed directory into place.
pub(crate) fn create_zarr_store(
    store_path: &Path,
    output_description: &str,
) -> Result<Arc<FilesystemStore>> {
    fs::create_dir_all(store_path).with_context(|| {
        format!(
            "create {output_description} Zarr store {}",
            store_path.display()
        )
    })?;
    Ok(Arc::new(FilesystemStore::new(store_path).with_context(
        || {
            format!(
                "open {output_description} Zarr store {}",
                store_path.display()
            )
        },
    )?))
}

/// Write root-level Zarr group metadata.
///
/// Root attributes are command schema contracts, not general command settings. The caller supplies
/// the schema-specific JSON object and this helper lets `zarrs` handle the V3 group metadata file.
pub(crate) fn write_zarr_root_metadata(
    store: Arc<FilesystemStore>,
    output_description: &str,
    attributes: Value,
) -> Result<()> {
    write_zarr_group_metadata(store, "/", output_description, attributes)
}

/// Write Zarr V3 group metadata for a non-root group.
///
/// Nested arrays such as `sparse/row` need their parent group metadata for readers that discover
/// arrays through the Zarr hierarchy rather than by direct filesystem paths.
pub(crate) fn write_zarr_group_metadata(
    store: Arc<FilesystemStore>,
    group_path: &str,
    output_description: &str,
    attributes: Value,
) -> Result<()> {
    let mut builder = GroupBuilder::new();
    builder.attributes(json_object(attributes)?);
    let group = builder
        .build(store, group_path)
        .with_context(|| format!("build {output_description} Zarr group {group_path}"))?;
    group
        .store_metadata()
        .with_context(|| format!("write {output_description} Zarr group {group_path} metadata"))
}

/// Write a small coordinate or metadata array as one chunk.
///
/// Zarr arrays can have zero-length dimensions. Such arrays have metadata but no chunks, so this
/// helper skips `store_chunk` when `values` is empty. Non-empty arrays are stored as one chunk
/// because coordinate and metadata arrays are usually read as complete vectors.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_single_chunk_zarr_array<T>(
    store: Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    dimension_names: &[&str],
    values: &[T],
    data_type: DataType,
    fill_value: T,
    attributes: Value,
) -> Result<()>
where
    T: Element + Copy,
    zarrs::array::builder::ArrayBuilderFillValue: From<T>,
{
    let expected_len =
        element_count(shape).with_context(|| format!("Zarr array {name} shape overflow"))?;
    ensure!(
        values.len() == expected_len,
        "Zarr array {name} has {} values but shape {:?} requires {}",
        values.len(),
        shape,
        expected_len
    );
    let chunk_shape: Vec<usize> = shape.iter().map(|dimension| (*dimension).max(1)).collect();
    let array = create_zarr_array(
        store,
        name,
        shape,
        &chunk_shape,
        dimension_names,
        data_type,
        fill_value,
        attributes,
    )?;
    if !values.is_empty() {
        let chunk_indices = vec![0; shape.len()];
        array
            .store_chunk(&chunk_indices, values)
            .with_context(|| format!("write Zarr chunk for array {name}"))?;
    }
    Ok(())
}

/// Create one Zarr V3 array and write its metadata.
///
/// This is the common `zarrs::ArrayBuilder` setup for cfDNAlab stores. It applies zstd,
/// native V3 dimension names, and caller-supplied attributes. Xarray reads V3 dimension names from
/// array metadata, so this helper deliberately does not add the V2-only `_ARRAY_DIMENSIONS`
/// attribute.
#[allow(clippy::too_many_arguments)]
pub(crate) fn create_zarr_array<T>(
    store: Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    chunk_shape: &[usize],
    dimension_names: &[&str],
    data_type: DataType,
    fill_value: T,
    attributes: Value,
) -> Result<zarrs::array::Array<FilesystemStore>>
where
    T: Element + Copy,
    ArrayBuilderFillValue: From<T>,
{
    create_zarr_array_with_fill_value(
        store,
        name,
        shape,
        chunk_shape,
        dimension_names,
        data_type,
        fill_value.into(),
        attributes,
    )
}

/// Create one Zarr V3 array using an explicit Zarr metadata fill value.
///
/// Most numeric arrays can use `create_zarr_array` with a Rust scalar fill value. Boolean arrays
/// need this lower-level helper because `zarrs` accepts boolean fill values through Zarr metadata,
/// not through the scalar conversion used for numeric primitive types.
#[allow(clippy::too_many_arguments)]
pub(crate) fn create_zarr_array_with_fill_value(
    store: Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    chunk_shape: &[usize],
    dimension_names: &[&str],
    data_type: DataType,
    fill_value: ArrayBuilderFillValue,
    attributes: Value,
) -> Result<zarrs::array::Array<FilesystemStore>> {
    ensure!(
        shape.len() == dimension_names.len(),
        "Zarr array {name} shape rank {} did not match {} dimension names",
        shape.len(),
        dimension_names.len()
    );
    element_count(shape).with_context(|| format!("Zarr array {name} shape overflow"))?;
    ensure!(
        shape.len() == chunk_shape.len(),
        "Zarr array {name} shape rank {} did not match chunk rank {}",
        shape.len(),
        chunk_shape.len()
    );
    ensure!(
        chunk_shape.iter().all(|dimension| *dimension > 0),
        "Zarr array {name} chunk shape {:?} contains a zero dimension",
        chunk_shape
    );

    let path = format!("/{name}");
    let array = ArrayBuilder::new(
        usize_slice_to_u64_vec(shape)?,
        usize_slice_to_u64_vec(chunk_shape)?,
        data_type,
        fill_value,
    )
    .bytes_to_bytes_codecs(vec![Arc::new(ZstdCodec::new(
        DEFAULT_ZARR_ZSTD_LEVEL.into(),
        false,
    ))])
    .dimension_names(Some(dimension_names.iter().copied()))
    .attributes(json_object(attributes)?)
    .build(store, path.as_str())
    .with_context(|| format!("build Zarr array {name}"))?;
    // Metadata is written before chunks inside the temporary store. The final output directory is
    // only moved into place after the command has written all arrays.
    array
        .store_metadata()
        .with_context(|| format!("write Zarr metadata for array {name}"))?;
    Ok(array)
}

/// Return a Zarr metadata fill value for boolean arrays.
pub(crate) fn bool_fill_value(value: bool) -> ArrayBuilderFillValue {
    FillValueMetadata::from(value).into()
}

#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_gc_bias",
    feature = "cmd_midpoints"
))]
pub(crate) use root_attribute_reader::read_zarr_root_attributes;

#[cfg(any(feature = "cmd_gc_bias"))]
pub(crate) use package_readers::{ensure_zarr_schema, read_zarr_array1, read_zarr_array2};

#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_gc_bias",
    feature = "cmd_midpoints"
))]
mod root_attribute_reader {
    use anyhow::{Context, Result};
    use serde_json::Value;
    use std::path::Path;

    /// Read root-level Zarr group attributes from `zarr.json`.
    pub(crate) fn read_zarr_root_attributes(path: &Path) -> Result<Value> {
        let metadata: Value =
            serde_json::from_str(&std::fs::read_to_string(path.join("zarr.json"))?)?;
        Ok(metadata
            .get("attributes")
            .cloned()
            .context("Zarr root metadata is missing attributes")?)
    }
}

#[cfg(any(feature = "cmd_gc_bias"))]
mod package_readers {
    use anyhow::{Context, Result, ensure};
    use ndarray::Array2;
    use serde_json::Value;
    use std::sync::Arc;
    use zarrs::{
        array::{Array, ElementOwned},
        filesystem::FilesystemStore,
    };

    /// Ensure a Zarr store advertises the expected cfDNAlab schema and version.
    pub(crate) fn ensure_zarr_schema(
        root_attributes: &Value,
        expected_schema: &str,
        expected_version: u32,
        package_name: &str,
    ) -> Result<()> {
        let schema = root_attributes
            .get("cfdnalab_schema")
            .and_then(Value::as_str);
        ensure!(
            schema == Some(expected_schema),
            "{package_name} schema mismatch: file={:?}, expected={expected_schema}",
            schema
        );
        let version = root_attributes
            .get("cfdnalab_schema_version")
            .and_then(Value::as_u64)
            .with_context(|| format!("{package_name} is missing cfdnalab_schema_version"))?;
        ensure!(
            version == u64::from(expected_version),
            "{package_name} schema version mismatch: file={}, expected={}; Incompatible with this version of cfDNAlab.",
            version,
            expected_version
        );
        Ok(())
    }

    /// Read a complete rank-1 Zarr array into memory.
    pub(crate) fn read_zarr_array1<T>(
        store: Arc<FilesystemStore>,
        array_path: &str,
    ) -> Result<Vec<T>>
    where
        T: ElementOwned,
    {
        let array = Array::open(store, array_path)?;
        Ok(array.retrieve_array_subset(&array.subset_all())?)
    }

    /// Read a complete rank-2 Zarr array into memory.
    pub(crate) fn read_zarr_array2<T>(
        store: Arc<FilesystemStore>,
        array_path: &str,
    ) -> Result<Array2<T>>
    where
        T: ElementOwned,
    {
        let array = Array::open(store, array_path)?;
        let shape = array.shape();
        ensure!(
            shape.len() == 2,
            "{array_path} must be a rank-2 array, found rank {}",
            shape.len()
        );
        let values: Vec<T> = array.retrieve_array_subset(&array.subset_all())?;
        let rows = usize::try_from(shape[0]).context("array row count exceeds usize")?;
        let cols = usize::try_from(shape[1]).context("array column count exceeds usize")?;
        Ok(Array2::from_shape_vec((rows, cols), values)?)
    }
}

/// Reject labels that cannot remain one stable public value.
///
/// JSON attributes can encode control characters, but downstream users often move labels into TSVs,
/// data frames, plots, and command-line examples. Rejecting control characters keeps labels
/// identical across public representations instead of silently rewriting them.
pub(crate) fn validate_zarr_label(value: &str, field_name: &str) -> Result<()> {
    ensure!(
        !value.chars().any(char::is_control),
        "{field_name} contains a control character and cannot be written without changing its value"
    );
    Ok(())
}

/// Build a zero-based signed integer coordinate axis.
///
/// cfDNAlab coordinate labels use `i32` because R's native integer vectors are signed. The helper
/// catches impossible axis sizes before they can truncate.
pub(crate) fn checked_index_axis(len: usize, axis_name: &str) -> Result<Vec<i32>> {
    (0..len)
        .map(|value| checked_i32(value, axis_name))
        .collect::<Result<Vec<_>>>()
}

/// Convert a metadata value to the public `i32` Zarr dtype.
pub(crate) fn checked_i32<T>(value: T, field_name: &str) -> Result<i32>
where
    T: TryInto<i32> + Copy + std::fmt::Display,
    T::Error: std::fmt::Debug,
{
    value
        .try_into()
        .map_err(|_| anyhow::anyhow!("{field_name} value {value} exceeds i32"))
}

/// Convert a metadata value to the public `i64` Zarr dtype.
#[cfg(any(feature = "cmd_ends"))]
pub(crate) fn checked_i64<T>(value: T, field_name: &str) -> Result<i64>
where
    T: TryInto<i64> + Copy + std::fmt::Display,
    T::Error: std::fmt::Debug,
{
    value
        .try_into()
        .map_err(|_| anyhow::anyhow!("{field_name} value {value} exceeds i64"))
}

/// Return the number of elements in a shape, or `None` on overflow.
pub(crate) fn element_count(shape: &[usize]) -> Option<usize> {
    shape
        .iter()
        .try_fold(1usize, |size, dimension| size.checked_mul(*dimension))
}

/// Convert Rust `usize` dimensions to the `u64` shape values expected by `zarrs`.
pub(crate) fn usize_slice_to_u64_vec(values: &[usize]) -> Result<Vec<u64>> {
    values
        .iter()
        .map(|value| {
            u64::try_from(*value).with_context(|| format!("Zarr dimension {value} exceeds u64"))
        })
        .collect()
}

/// Return an object map for Zarr attributes.
///
/// Callers pass attributes with `json!({ ... })`. Non-object values are schema mistakes in command
/// writers, so this helper errors instead of silently replacing them.
pub(crate) fn json_object(value: Value) -> Result<Map<String, Value>> {
    match value {
        Value::Object(map) => Ok(map),
        other => Err(anyhow::anyhow!(
            "Zarr attributes must be a JSON object, got {other}"
        )),
    }
}

#[cfg(test)]
mod tests {
    include!("zarr_tests.rs");
}
