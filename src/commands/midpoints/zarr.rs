//! Midpoint profile Zarr writer.
//!
//! `zarrs` owns the Zarr metadata format, codec pipeline, chunk serialization, and string
//! encoding. This module owns the cfDNAlab schema that sits on top of Zarr:
//!
//! - which arrays are present
//! - which axis names and attributes downstream readers should see
//! - how internal `usize` and `u64` values are narrowed into public Zarr integer dtypes
//! - how the large count tensor is split into chunks before handing each chunk to `zarrs`
//!
//! Small coordinate arrays are written as one chunk. The primary `counts` tensor is chunked by this
//! module so a large run does not create one huge Zarr chunk.
//!
//! Coordinate and label arrays use signed `i32` values for smoother R handling. Count-like
//! metadata and byte lengths use unsigned `u32` values because those fields are never labels.

use crate::{
    commands::midpoints::{
        group_index::{MidpointGroupSummary, ordered_midpoint_group_summaries},
        postprocess::ProfileLayout,
    },
    shared::length_axis::LengthAxis,
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use ndarray::ArrayView3;
use serde_json::{Map, Value, json};
use std::{fs, path::Path, sync::Arc};
use zarrs::{
    array::{ArrayBuilder, DataType, Element, codec::ZstdCodec, data_type},
    filesystem::FilesystemStore,
    group::GroupBuilder,
};

/// Zstd compression level used for all midpoint Zarr arrays.
///
/// This keeps all arrays in the store compressed the same way unless a future downstream
/// compatibility test shows that string arrays need different handling.
const ZSTD_LEVEL: i32 = 3;

/// Soft target for count-array chunk size.
///
/// The exact chunk shape depends on the output shape, but this target keeps each count chunk near
/// sixteen MiB for `f32` values. Coordinate arrays are tiny and are stored as one chunk.
const TARGET_COUNT_CHUNK_CELLS: usize = 4_000_000;

/// Current midpoint Zarr schema version.
const CFDNALAB_MIDPOINT_SCHEMA_VERSION: u32 = 1;

/// Write midpoint profiles as a Zarr V3 hierarchy.
///
/// The public tensor is `counts[group, length_bin, position]`. The `group` coordinate is the same
/// zero-based `group_idx` written to `<prefix>.group_index.tsv`, and the count tensor uses that
/// value directly as its row index. The writer validates this before writing so group metadata
/// cannot silently drift away from the rows in `counts`.
///
/// Axis metadata is stored as ordinary arrays so downstream users can build dataframes without
/// consulting a sidecar file. General command settings stay in `midpoint_settings.json`.
pub(super) fn write_midpoint_profiles_zarr(
    store_path: &Path,
    counts: ArrayView3<'_, f32>,
    group_idx_to_name: &FxHashMap<u64, String>,
    eligible_interval_counts: &FxHashMap<u64, usize>,
    length_axis: &LengthAxis,
    profile_layout: ProfileLayout,
) -> Result<()> {
    // Align group metadata to the count tensor before any arrays are written
    let ordered_groups =
        ordered_midpoint_group_summaries(group_idx_to_name, eligible_interval_counts)?;
    let num_groups = ordered_groups.len();
    let num_length_bins = length_axis.num_bins();
    let num_positions = profile_layout.output_positions;
    ensure!(
        num_groups > 0,
        "midpoint Zarr output requires at least one group"
    );
    ensure!(
        num_length_bins > 0,
        "midpoint Zarr output requires at least one length bin"
    );
    ensure!(
        num_positions > 0,
        "midpoint Zarr output requires at least one position bin"
    );
    ensure!(
        counts.shape() == [num_groups, num_length_bins, num_positions],
        "midpoint profile shape {:?} did not match expected Zarr shape [{}, {}, {}]",
        counts.shape(),
        num_groups,
        num_length_bins,
        num_positions
    );

    let store = create_store(store_path)?;
    write_group_metadata(store.clone())?;

    // Public coordinate arrays and small metadata arrays
    let group: Vec<i32> = ordered_groups
        .iter()
        .map(|group| checked_i32(group.group_idx, "group_idx"))
        .collect::<Result<Vec<_>>>()?;
    let length_bin: Vec<i32> = checked_index_axis(num_length_bins, "length_bin")?;
    let position: Vec<i32> = checked_index_axis(num_positions, "position")?;
    let eligible_intervals: Vec<u32> = ordered_groups
        .iter()
        .map(|group| checked_u32(group.eligible_intervals, "eligible_intervals"))
        .collect::<Result<Vec<_>>>()?;
    let group_names: Vec<&str> = ordered_groups
        .iter()
        .map(|group| group.group_name)
        .collect();

    // Keep the byte arrays as a temporary compatibility fallback until downstream tests decide
    // whether all target R/Python readers handle Zarr string arrays cleanly
    let (group_name_utf8, group_name_nbytes, group_name_width) =
        encode_group_names(&ordered_groups)?;
    let group_name_byte: Vec<i32> = checked_index_axis(group_name_width, "group_name_byte")?;
    let (length_start_bp, length_end_bp) = length_axis_coordinate_arrays(length_axis);
    let (position_bin_start_bp, position_bin_end_bp) =
        position_bin_coordinate_arrays(profile_layout)?;

    // Primary count tensor
    let count_shape = [num_groups, num_length_bins, num_positions];
    let count_chunk_shape = counts_chunk_shape(count_shape)?;
    write_count_tensor_array(
        store.clone(),
        "counts",
        &count_shape,
        &count_chunk_shape,
        &["group", "length_bin", "position"],
        counts,
        json!({
            "long_name": "weighted midpoint count",
            "units": "weighted_midpoint_count",
        }),
    )?;

    // Group axis and group metadata
    write_single_chunk_array(
        store.clone(),
        "group",
        &[num_groups],
        &["group"],
        &group,
        data_type::int32(),
        0i32,
        json!({
            "long_name": "zero-based group index",
            "description": "Matches group_idx in group_index.tsv and indexes axis 0 of counts",
        }),
    )?;
    write_single_chunk_array(
        store.clone(),
        "eligible_intervals",
        &[num_groups],
        &["group"],
        &eligible_intervals,
        data_type::uint32(),
        0u32,
        json!({
            "long_name": "profile-eligible input intervals retained per group",
        }),
    )?;
    write_single_chunk_array(
        store.clone(),
        "group_name",
        &[num_groups],
        &["group"],
        &group_names,
        data_type::string(),
        "",
        json!({
            "long_name": "group name",
            "description": "Primary human-readable group labels. The string codec is recorded in the array metadata.",
            "compatibility_fallback": ["group_name_utf8", "group_name_nbytes"],
        }),
    )?;
    write_single_chunk_array(
        store.clone(),
        "group_name_utf8",
        &[num_groups, group_name_width],
        &["group", "group_name_byte"],
        &group_name_utf8,
        data_type::uint8(),
        0u8,
        json!({
            "long_name": "UTF-8 group names padded with zero bytes",
            "paired_length_array": "group_name_nbytes",
            "description": "Decode group_name_utf8[group, group_name_byte] row-by-row with group_name_nbytes, ignoring trailing zero padding",
        }),
    )?;
    write_single_chunk_array(
        store.clone(),
        "group_name_nbytes",
        &[num_groups],
        &["group"],
        &group_name_nbytes,
        data_type::uint32(),
        0u32,
        json!({
            "long_name": "number of UTF-8 bytes used by each group name",
        }),
    )?;
    write_single_chunk_array(
        store.clone(),
        "group_name_byte",
        &[group_name_width],
        &["group_name_byte"],
        &group_name_byte,
        data_type::int32(),
        0i32,
        json!({
            "long_name": "byte offset within padded UTF-8 group names",
        }),
    )?;

    // Length axis
    write_single_chunk_array(
        store.clone(),
        "length_bin",
        &[num_length_bins],
        &["length_bin"],
        &length_bin,
        data_type::int32(),
        0i32,
        json!({}),
    )?;
    write_single_chunk_array(
        store.clone(),
        "length_start_bp",
        &[num_length_bins],
        &["length_bin"],
        &length_start_bp,
        data_type::uint32(),
        0u32,
        json!({
            "long_name": "inclusive fragment length bin start",
            "units": "bp",
        }),
    )?;
    write_single_chunk_array(
        store.clone(),
        "length_end_bp",
        &[num_length_bins],
        &["length_bin"],
        &length_end_bp,
        data_type::uint32(),
        0u32,
        json!({
            "long_name": "exclusive fragment length bin end",
            "units": "bp",
        }),
    )?;

    // Position axis
    write_single_chunk_array(
        store.clone(),
        "position",
        &[num_positions],
        &["position"],
        &position,
        data_type::int32(),
        0i32,
        json!({}),
    )?;
    write_single_chunk_array(
        store.clone(),
        "position_bin_start_bp",
        &[num_positions],
        &["position"],
        &position_bin_start_bp,
        data_type::int32(),
        0i32,
        json!({
            "long_name": "inclusive interval-relative position bin start",
            "units": "bp",
        }),
    )?;
    write_single_chunk_array(
        store,
        "position_bin_end_bp",
        &[num_positions],
        &["position"],
        &position_bin_end_bp,
        data_type::int32(),
        0i32,
        json!({
            "long_name": "exclusive interval-relative position bin end",
            "units": "bp",
        }),
    )?;

    Ok(())
}

/// Encode variable-length group names as a padded byte matrix plus per-row byte lengths.
///
/// `group_name_utf8[group, group_name_byte]` stores raw UTF-8 bytes and pads shorter names with
/// zero bytes. `group_name_nbytes[group]` stores the real byte length for each row. The separate
/// `group_name_byte` coordinate labels the second dimension for readers that expose Zarr arrays as
/// labeled data.
fn encode_group_names(groups: &[MidpointGroupSummary<'_>]) -> Result<(Vec<u8>, Vec<u32>, usize)> {
    let width = groups
        .iter()
        .map(|group| group.group_name.len())
        .max()
        .unwrap_or(0)
        .max(1);
    let encoded_len = groups
        .len()
        .checked_mul(width)
        .context("group name byte matrix shape overflow")?;
    let mut encoded = vec![0u8; encoded_len];
    let mut lengths = Vec::with_capacity(groups.len());
    for (group_index, group) in groups.iter().enumerate() {
        let bytes = group.group_name.as_bytes();
        lengths.push(checked_u32(bytes.len(), "group_name_nbytes")?);
        let start = group_index * width;
        encoded[start..start + bytes.len()].copy_from_slice(bytes);
    }
    Ok((encoded, lengths, width))
}

/// Build half-open fragment length-bin coordinate arrays.
///
/// `LengthAxis` stores bin edges as `[start0, start1, ..., final_end]`. Zarr readers should not
/// have to reconstruct the half-open bins from settings JSON, so the writer stores start and end
/// arrays directly.
fn length_axis_coordinate_arrays(length_axis: &LengthAxis) -> (Vec<u32>, Vec<u32>) {
    let starts = length_axis.edges()[..length_axis.num_bins()].to_vec();
    let ends = length_axis.edges()[1..].to_vec();
    (starts, ends)
}

/// Build half-open position-bin coordinate arrays.
///
/// These are interval-relative profile bins after smoothing flanks have been trimmed. With
/// `--bin-size 1`, every bin is one base wide. With larger bins, the last bin can be shorter than
/// `bin_size`, so the explicit end array is part of the public analysis metadata.
fn position_bin_coordinate_arrays(profile_layout: ProfileLayout) -> Result<(Vec<i32>, Vec<i32>)> {
    let mut starts = Vec::with_capacity(profile_layout.output_positions);
    let mut ends = Vec::with_capacity(profile_layout.output_positions);
    for position_index in 0..profile_layout.output_positions {
        let start = position_index
            .checked_mul(profile_layout.bin_size as usize)
            .context("position axis start overflow")?;
        let end = (start + profile_layout.bin_size as usize).min(profile_layout.output_len);
        starts.push(checked_i32(start, "position_bin_start_bp")?);
        ends.push(checked_i32(end, "position_bin_end_bp")?);
    }
    Ok((starts, ends))
}

/// Build a zero-based signed integer coordinate axis.
///
/// Zarr coordinate arrays use `i32` because R's native integer vectors are signed. Keeping axis
/// labels signed avoids unsigned-specific reader behavior while still catching impossible axis
/// sizes before they can truncate.
fn checked_index_axis(len: usize, axis_name: &str) -> Result<Vec<i32>> {
    (0..len)
        .map(|value| checked_i32(value, axis_name))
        .collect::<Result<Vec<_>>>()
}

/// Convert a metadata value to the public `i32` Zarr dtype.
///
/// Group indices originate as `u64` and other axes originate as `usize`, but public coordinate
/// arrays use `i32` for R/Python reader compatibility. This helper keeps every narrowing point
/// explicit and erroring.
fn checked_i32<T>(value: T, field_name: &str) -> Result<i32>
where
    T: TryInto<i32> + Copy + std::fmt::Display,
    T::Error: std::fmt::Debug,
{
    value
        .try_into()
        .map_err(|_| anyhow::anyhow!("{field_name} value {value} exceeds i32"))
}

/// Convert metadata counts to the public `u32` Zarr dtype.
///
/// The underlying Rust counts use platform-sized integers. The public Zarr metadata uses `u32` to
/// keep small coordinate arrays compact while still failing clearly on unrealistic overflow.
fn checked_u32<T>(value: T, field_name: &str) -> Result<u32>
where
    T: TryInto<u32> + Copy + std::fmt::Display,
    T::Error: std::fmt::Debug,
{
    value
        .try_into()
        .map_err(|_| anyhow::anyhow!("{field_name} value {value} exceeds u32"))
}

/// Open or create the filesystem-backed Zarr store directory.
///
/// Final-output replacement is handled before this writer runs. This function only ensures the
/// temporary output path exists and lets `zarrs` wrap it as a store.
fn create_store(store_path: &Path) -> Result<Arc<FilesystemStore>> {
    fs::create_dir_all(store_path)
        .with_context(|| format!("create midpoint Zarr store {}", store_path.display()))?;
    Ok(Arc::new(FilesystemStore::new(store_path).with_context(
        || format!("open midpoint Zarr store {}", store_path.display()),
    )?))
}

/// Write root-level metadata for the midpoint profile store.
///
/// These attributes are not command settings. They are the short machine-readable schema contract
/// needed by downstream readers to identify the primary array and its axis order.
fn write_group_metadata(store: Arc<FilesystemStore>) -> Result<()> {
    let mut builder = GroupBuilder::new();
    builder.attributes(json_object(json!({
        "cfdnalab_schema": "midpoint_profiles",
        "cfdnalab_schema_version": CFDNALAB_MIDPOINT_SCHEMA_VERSION,
        "primary_array": "counts",
        "dimension_names": ["group", "length_bin", "position"],
        "count_units": "weighted_midpoint_count",
    }))?);
    let group = builder
        .build(store, "/")
        .context("build midpoint Zarr root group")?;
    group
        .store_metadata()
        .context("write midpoint Zarr root metadata")
}

/// Write a small array as one Zarr chunk.
///
/// Coordinate arrays and small metadata arrays are tiny compared with `counts`. Storing each as a
/// single chunk keeps the writer simple and makes downstream reads straightforward. `zarrs` handles
/// element encoding and compression after this helper validates that the supplied value count
/// matches the declared shape.
fn write_single_chunk_array<T>(
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
    let array = create_zarr_array(
        store,
        name,
        shape,
        shape,
        dimension_names,
        data_type,
        fill_value,
        attributes,
    )?;
    let chunk_indices = vec![0; shape.len()];
    array
        .store_chunk(&chunk_indices, values)
        .with_context(|| format!("write Zarr chunk for array {name}"))?;
    Ok(())
}

/// Write the primary `counts[group, length_bin, position]` tensor in bounded chunks.
///
/// The chunk shape keeps stored chunks reasonably sized even when the full tensor is large.
/// `zarrs::store_array_subset` writes the full logical array and handles boundary chunks whose
/// stored chunk shape extends past the array shape.
fn write_count_tensor_array(
    store: Arc<FilesystemStore>,
    name: &str,
    shape: &[usize; 3],
    chunk_shape: &[usize; 3],
    dimension_names: &[&str; 3],
    counts: ArrayView3<'_, f32>,
    attributes: Value,
) -> Result<()> {
    ensure!(
        counts.shape() == shape.as_slice(),
        "Zarr array {name} has count shape {:?} but declared shape is {:?}",
        counts.shape(),
        shape
    );
    let array = create_zarr_array(
        store,
        name,
        shape,
        chunk_shape,
        dimension_names,
        data_type::float32(),
        0.0f32,
        attributes,
    )?;

    let owned_counts = counts.as_standard_layout();
    let count_values = owned_counts
        .as_slice()
        .context("standard-layout midpoint counts were not contiguous")?;
    array
        .store_array_subset(&array.subset_all(), count_values)
        .with_context(|| format!("write Zarr array {name}"))?;
    Ok(())
}

/// Create one Zarr array and store its metadata.
///
/// This is the common `zarrs::ArrayBuilder` setup for all arrays in this store. It applies zstd,
/// native Zarr V3 dimension names, and caller-supplied attributes consistently. Xarray reads V3
/// dimension metadata from the array-level `dimension_names` field, so the V2-only
/// `_ARRAY_DIMENSIONS` attribute is deliberately not written here.
#[allow(clippy::too_many_arguments)]
fn create_zarr_array<T>(
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
    .bytes_to_bytes_codecs(vec![Arc::new(ZstdCodec::new(ZSTD_LEVEL.into(), false))])
    .dimension_names(Some(dimension_names.iter().copied()))
    .attributes(json_object(attributes)?)
    .build(store, path.as_str())
    .with_context(|| format!("build Zarr array {name}"))?;
    // Metadata is written before chunks inside the temporary store
    // The final output directory is only moved into place after all arrays are written
    array
        .store_metadata()
        .with_context(|| format!("write Zarr metadata for array {name}"))?;
    Ok(array)
}

/// Return the number of elements in a shape, or `None` on overflow.
fn element_count(shape: &[usize]) -> Option<usize> {
    shape
        .iter()
        .try_fold(1usize, |size, dimension| size.checked_mul(*dimension))
}

/// Choose a chunk shape for the count tensor.
///
/// If the full tensor is already small, it is written as one chunk. Otherwise the position axis is
/// filled first because common downstream reads are expected to take a group subset over all length
/// bins and positions. Length and group chunk sizes are then reduced to keep the cell count near
/// `TARGET_COUNT_CHUNK_CELLS`.
fn counts_chunk_shape(shape: [usize; 3]) -> Result<[usize; 3]> {
    let total_cells = shape
        .iter()
        .try_fold(1usize, |cells, dimension| cells.checked_mul(*dimension))
        .context("midpoint Zarr count shape overflow")?;
    if total_cells <= TARGET_COUNT_CHUNK_CELLS {
        return Ok(shape);
    }

    let position_chunk = shape[2].min(TARGET_COUNT_CHUNK_CELLS).max(1);
    let length_chunk = shape[1]
        .min((TARGET_COUNT_CHUNK_CELLS / position_chunk).max(1))
        .max(1);
    let group_chunk = shape[0]
        .min((TARGET_COUNT_CHUNK_CELLS / (position_chunk * length_chunk)).max(1))
        .max(1);
    Ok([group_chunk, length_chunk, position_chunk])
}

/// Convert Rust `usize` dimensions to the `u64` shape values expected by `zarrs`.
fn usize_slice_to_u64_vec(values: &[usize]) -> Result<Vec<u64>> {
    values
        .iter()
        .map(|value| {
            u64::try_from(*value).with_context(|| format!("Zarr dimension {value} exceeds u64"))
        })
        .collect()
}

/// Return an object map for attributes.
///
/// Callers pass attributes with `json!({ ... })`. Non-object values are schema mistakes in this
/// module, so they return an error instead of being silently rewritten.
fn json_object(value: Value) -> Result<Map<String, Value>> {
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
