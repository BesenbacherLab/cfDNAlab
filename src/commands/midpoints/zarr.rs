//! Midpoint profile Zarr writer.
//!
//! `zarrs` owns the Zarr metadata format, codec pipeline, and chunk serialization. This module
//! owns the cfDNAlab schema that sits on top of Zarr:
//!
//! - which arrays are present
//! - which axis names and attributes downstream readers should see
//! - how internal `usize` and `u64` values are narrowed into public Zarr integer dtypes
//! - how the large count tensor is split into chunks before handing each chunk to `zarrs`
//!
//! Small coordinate arrays are written as one chunk. The primary `counts` tensor is chunked by this
//! module so a large run does not create one huge Zarr chunk.
//!
//! Coordinate arrays use signed `i32` values for smoother R handling. Count-like metadata uses
//! signed `i32` values as well when values are small and non-negative, so R readers can keep native
//! integer columns without extra unsigned handling.

use crate::{
    commands::midpoints::{
        group_index::ordered_midpoint_group_summaries, postprocess::ProfileLayout,
    },
    shared::{
        length_axis::LengthAxis,
        zarr::{
            ZARR_FLOAT32_FILL_VALUE, ZARR_INT32_FILL_VALUE, checked_i32, checked_index_axis,
            create_zarr_array, create_zarr_store, validate_zarr_label,
            write_single_chunk_zarr_array, write_zarr_root_metadata,
        },
    },
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use ndarray::ArrayView3;
use serde_json::{Value, json};
use std::{path::Path, sync::Arc};
use zarrs::{array::data_type, filesystem::FilesystemStore};

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

    let store = create_zarr_store(store_path, "midpoint")?;
    write_group_metadata(store.clone())?;

    // Public coordinate arrays and small metadata arrays
    let group: Vec<i32> = ordered_groups
        .iter()
        .map(|group| checked_i32(group.group_idx, "group_idx"))
        .collect::<Result<Vec<_>>>()?;
    let length_bin: Vec<i32> = checked_index_axis(num_length_bins, "length_bin")?;
    let position: Vec<i32> = checked_index_axis(num_positions, "position")?;
    let eligible_intervals: Vec<i32> = ordered_groups
        .iter()
        .map(|group| checked_i32(group.eligible_intervals, "eligible_intervals"))
        .collect::<Result<Vec<_>>>()?;
    let group_names: Vec<&str> = ordered_groups
        .iter()
        .map(|group| group.group_name)
        .collect();
    for group_name in &group_names {
        validate_zarr_label(group_name, "group_name")?;
    }
    let (length_start_bp, length_end_bp) = length_axis_coordinate_arrays(length_axis)?;
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
    write_single_chunk_zarr_array(
        store.clone(),
        "group",
        &[num_groups],
        &["group"],
        &group,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based group index",
            "description": "Matches group_idx in group_index.tsv and indexes axis 0 of counts",
            "label_field": "group_name",
            "labels": group_names,
        }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "eligible_intervals",
        &[num_groups],
        &["group"],
        &eligible_intervals,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "profile-eligible input intervals retained per group",
        }),
    )?;
    // Length axis
    write_single_chunk_zarr_array(
        store.clone(),
        "length_bin",
        &[num_length_bins],
        &["length_bin"],
        &length_bin,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({}),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "length_start_bp",
        &[num_length_bins],
        &["length_bin"],
        &length_start_bp,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "inclusive fragment length bin start",
            "units": "bp",
        }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "length_end_bp",
        &[num_length_bins],
        &["length_bin"],
        &length_end_bp,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "exclusive fragment length bin end",
            "units": "bp",
        }),
    )?;

    // Position axis
    write_single_chunk_zarr_array(
        store.clone(),
        "position",
        &[num_positions],
        &["position"],
        &position,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({}),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "position_bin_start_bp",
        &[num_positions],
        &["position"],
        &position_bin_start_bp,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "inclusive interval-relative position bin start",
            "units": "bp",
        }),
    )?;
    write_single_chunk_zarr_array(
        store,
        "position_bin_end_bp",
        &[num_positions],
        &["position"],
        &position_bin_end_bp,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "exclusive interval-relative position bin end",
            "units": "bp",
        }),
    )?;

    Ok(())
}

/// Build half-open fragment length-bin coordinate arrays.
///
/// `LengthAxis` stores bin edges as `[start0, start1, ..., final_end]`. Zarr readers should not
/// have to reconstruct the half-open bins from settings JSON, so the writer stores start and end
/// arrays directly.
fn length_axis_coordinate_arrays(length_axis: &LengthAxis) -> Result<(Vec<i32>, Vec<i32>)> {
    let starts = length_axis.edges()[..length_axis.num_bins()]
        .iter()
        .map(|value| checked_i32(*value, "length_start_bp"))
        .collect::<Result<Vec<_>>>()?;
    let ends = length_axis.edges()[1..]
        .iter()
        .map(|value| checked_i32(*value, "length_end_bp"))
        .collect::<Result<Vec<_>>>()?;
    Ok((starts, ends))
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

/// Write root-level metadata for the midpoint profile store.
///
/// These attributes are not command settings. They are the short machine-readable schema contract
/// needed by downstream readers to identify the primary array and its axis order.
fn write_group_metadata(store: Arc<FilesystemStore>) -> Result<()> {
    write_zarr_root_metadata(
        store,
        "midpoint",
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": CFDNALAB_MIDPOINT_SCHEMA_VERSION,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )
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
        ZARR_FLOAT32_FILL_VALUE,
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

#[cfg(test)]
mod tests {
    include!("zarr_tests.rs");
}
