//! Shared Zarr writing utilities.
//!
//! These helpers are deliberately low-level. They centralize the mechanical parts that should be
//! identical across cfDNAlab Zarr outputs:
//!
//! - opening a filesystem-backed Zarr store
//! - creating V3 arrays with zstd compression and native dimension names
//! - writing small coordinate/metadata arrays as a single chunk
//! - validating JSON attributes and public integer dtype narrowing
//! - validating public labels stored in JSON attributes
//!
//! Command modules still own their public schemas. This module should not know what a midpoint
//! profile, end motif, group, window, or length bin means.

use anyhow::{Context, Result, ensure};
use serde_json::{Map, Value};
use std::{fs, path::Path, sync::Arc};
use zarrs::{
    array::{ArrayBuilder, DataType, Element, codec::ZstdCodec},
    filesystem::FilesystemStore,
    group::GroupBuilder,
};

/// Default zstd level used by cfDNAlab Zarr stores.
///
/// Keeping this in one place avoids midpoint and ends drifting to different compression settings
/// unless downstream compatibility or performance tests give a concrete reason.
pub(crate) const DEFAULT_ZARR_ZSTD_LEVEL: i32 = 3;

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
    let mut builder = GroupBuilder::new();
    builder.attributes(json_object(attributes)?);
    let group = builder
        .build(store, "/")
        .with_context(|| format!("build {output_description} Zarr root group"))?;
    group
        .store_metadata()
        .with_context(|| format!("write {output_description} Zarr root metadata"))
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
    zarrs::array::builder::ArrayBuilderFillValue: From<T>,
{
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

/// Convert a metadata value to the public `u32` Zarr dtype.
pub(crate) fn checked_u32<T>(value: T, field_name: &str) -> Result<u32>
where
    T: TryInto<u32> + Copy + std::fmt::Display,
    T::Error: std::fmt::Debug,
{
    value
        .try_into()
        .map_err(|_| anyhow::anyhow!("{field_name} value {value} exceeds u32"))
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
