//! End-motif Zarr writer.
//!
//! `zarrs` owns the Zarr metadata, chunk encoding, and compression. This module owns the cfDNAlab
//! schema:
//!
//! - `counts[row, motif]` for dense output
//! - `sparse/{row,motif,count,shape,sparse_dimension}` for sparse COO output
//! - fixed-width motif labels and row metadata needed to interpret the counts without TSV sidecars
//!
//! The code below should only do work that is specific to this schema or to clearer validation.
//! Low-level details such as writing `zarr.json`, applying zstd, V3 dimension names, and chunk
//! serialization are delegated to `zarrs`.

use crate::shared::{
    bed::GroupedWindows,
    blacklist::compute_blacklist_overlap,
    interval::Interval,
    io::dot_join,
    windowing::WindowBinInfo,
    zarr::{
        ZARR_ASCII_FILL_VALUE, ZARR_FLOAT64_FILL_VALUE, ZARR_INT32_FILL_VALUE,
        ZARR_INT64_FILL_VALUE, checked_i32, checked_i64, checked_index_axis, create_zarr_array,
        create_zarr_store, validate_zarr_label, write_single_chunk_zarr_array,
        write_zarr_group_metadata, write_zarr_root_metadata,
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use ndarray::{Array2, ArrayView2};
use serde_json::json;
use std::{
    path::{Path, PathBuf},
    sync::Arc,
};
use zarrs::{array::data_type, filesystem::FilesystemStore};

const CFDNALAB_END_MOTIF_SCHEMA_VERSION: u32 = 1;

/// Soft target for dense count chunks.
///
/// This keeps dense `f64` chunks near sixteen MiB. Small dense matrices still use one chunk.
const TARGET_DENSE_COUNT_CHUNK_CELLS: usize = 2_000_000;

/// Row metadata for the public `row` axis.
pub(crate) enum EndMotifRowMetadata<'a> {
    /// One row covering all selected chromosomes.
    Global,
    /// Ordinary genomic windows in the same row order as the count matrix.
    Windows {
        bin_info: &'a [WindowBinInfo],
        row_mode: EndWindowRowMode,
    },
    /// Grouped BED rows in the same row order as the count matrix.
    Groups(Vec<EndGroupSummary<'a>>),
}

/// Source of ordinary genomic count rows.
pub(crate) enum EndWindowRowMode {
    Size,
    Bed,
}

/// One grouped-BED output row.
#[derive(Debug)]
pub(crate) struct EndGroupSummary<'a> {
    pub(crate) group_idx: u64,
    pub(crate) group_name: &'a str,
    pub(crate) eligible_windows: usize,
    pub(crate) blacklisted_fraction: f64,
}

/// Write dense or sparse end-motif counts as a self-contained Zarr V3 store.
pub(crate) fn write_end_motif_zarr(
    output_dir: &Path,
    prefix: &str,
    bins: &[FxHashMap<String, f64>],
    motifs: &[String],
    row_metadata: EndMotifRowMetadata<'_>,
    write_dense_output: bool,
) -> Result<PathBuf> {
    let store_path = output_dir.join(dot_join(&[prefix, "end_motifs.zarr"]));
    let store = create_zarr_store(&store_path, "end-motif")?;
    let storage_mode = if write_dense_output {
        "dense"
    } else {
        "sparse_coo"
    };
    let row_mode = row_mode_name(&row_metadata);
    write_root_metadata(store.clone(), storage_mode, row_mode)?;
    write_motif_metadata(store.clone(), motifs)?;
    write_row_metadata(store.clone(), row_metadata, bins.len())?;

    if write_dense_output {
        let counts = stack_end_motif_counts(bins, motifs)?;
        write_dense_counts(store, counts.view())?;
    } else {
        write_sparse_counts(store, bins, motifs)?;
    }

    Ok(store_path)
}

/// Build grouped row metadata in count-row order.
///
/// The ends reducer addresses grouped output rows directly by `group_idx`. This helper therefore
/// walks row indices `0..n_groups` and rejects missing indices rather than sorting metadata into a
/// different order after counts have already been built.
pub(crate) fn grouped_end_row_metadata<'a>(
    group_idx_to_name: &'a FxHashMap<u64, String>,
    chromosomes: &[String],
    grouped_windows_map: &FxHashMap<String, GroupedWindows>,
    blacklist_map: &FxHashMap<String, Vec<Interval<u64>>>,
) -> Result<Vec<EndGroupSummary<'a>>> {
    let mut eligible_windows_by_group: FxHashMap<u64, usize> = group_idx_to_name
        .keys()
        .map(|&group_idx| (group_idx, 0usize))
        .collect();
    let mut total_bp_by_group: FxHashMap<u64, u64> = group_idx_to_name
        .keys()
        .map(|&group_idx| (group_idx, 0u64))
        .collect();
    let mut blacklisted_bp_by_group: FxHashMap<u64, f64> = group_idx_to_name
        .keys()
        .map(|&group_idx| (group_idx, 0.0f64))
        .collect();

    for chromosome in chromosomes {
        let windows = grouped_windows_map
            .get(chromosome)
            .map(|windows| windows.windows_as_slice())
            .unwrap_or(&[]);
        let blacklist_intervals = blacklist_map
            .get(chromosome)
            .map(|intervals| intervals.as_slice())
            .unwrap_or(&[]);
        let mut blacklist_ptr = 0usize;
        for window in windows {
            let (start, end, group_idx) = window.as_tuple();
            ensure!(
                group_idx_to_name.contains_key(&group_idx),
                "grouped end-motif window references group_idx {group_idx}, but no group name exists for it"
            );
            let window_bp = end
                .checked_sub(start)
                .context("grouped end-motif window end must be >= start")?;
            let blacklist_fraction = compute_blacklist_overlap(
                blacklist_intervals,
                Interval::new(start, end)?,
                0,
                &mut blacklist_ptr,
            );
            *eligible_windows_by_group
                .get_mut(&group_idx)
                .expect("group_idx was validated above") += 1;
            *total_bp_by_group
                .get_mut(&group_idx)
                .expect("group_idx was validated above") += window_bp;
            *blacklisted_bp_by_group
                .get_mut(&group_idx)
                .expect("group_idx was validated above") += blacklist_fraction * window_bp as f64;
        }
    }

    let mut summaries = Vec::with_capacity(group_idx_to_name.len());
    for row_index in 0..group_idx_to_name.len() {
        let group_idx = row_index as u64;
        let Some(group_name) = group_idx_to_name.get(&group_idx) else {
            let mut observed_indices: Vec<u64> = group_idx_to_name.keys().copied().collect();
            observed_indices.sort_unstable();
            bail!(
                "end-motif group indices must match count rows 0..{} but observed {:?}",
                group_idx_to_name.len().saturating_sub(1),
                observed_indices
            );
        };
        validate_zarr_label(group_name, "group_name")?;
        let total_bp = total_bp_by_group.get(&group_idx).copied().unwrap_or(0);
        let blacklisted_fraction = if total_bp == 0 {
            0.0
        } else {
            blacklisted_bp_by_group
                .get(&group_idx)
                .copied()
                .unwrap_or(0.0)
                / total_bp as f64
        };
        summaries.push(EndGroupSummary {
            group_idx,
            group_name,
            eligible_windows: eligible_windows_by_group
                .get(&group_idx)
                .copied()
                .unwrap_or(0),
            blacklisted_fraction,
        });
    }

    Ok(summaries)
}

/// Stack sparse per-window motif maps into a dense matrix with a fixed column order.
pub(crate) fn stack_end_motif_counts(
    bins: &[FxHashMap<String, f64>],
    motifs: &[String],
) -> Result<Array2<f64>> {
    let mut counts = Array2::<f64>::zeros((bins.len(), motifs.len()));
    let motif_columns: FxHashMap<&String, usize> = motifs
        .iter()
        .enumerate()
        .map(|(column, motif)| (motif, column))
        .collect();

    for (row, bin) in bins.iter().enumerate() {
        for (motif, &count) in bin {
            let column = motif_columns
                .get(motif)
                .copied()
                .with_context(|| format!("missing output column for end-motif label '{motif}'"))?;
            counts[(row, column)] = count;
        }
    }

    Ok(counts)
}

fn row_mode_name(row_metadata: &EndMotifRowMetadata<'_>) -> &'static str {
    match row_metadata {
        EndMotifRowMetadata::Global => "global",
        EndMotifRowMetadata::Windows { row_mode, .. } => match row_mode {
            EndWindowRowMode::Size => "size",
            EndWindowRowMode::Bed => "bed",
        },
        EndMotifRowMetadata::Groups(_) => "grouped_bed",
    }
}

/// Write the root schema contract for the end-motif store.
///
/// Dense and sparse outputs expose different primary data locations. The root attributes are the
/// cheap way for downstream readers to discover which representation was written before opening
/// count arrays.
fn write_root_metadata(
    store: Arc<FilesystemStore>,
    storage_mode: &str,
    row_mode: &str,
) -> Result<()> {
    let mut attributes = json!({
        "cfdnalab_schema": "end_motif_counts",
        "cfdnalab_schema_version": CFDNALAB_END_MOTIF_SCHEMA_VERSION,
        "storage_mode": storage_mode,
        "row_mode": row_mode,
        "count_units": "weighted_end_motif_count",
        "primary_array": null,
        "primary_group": null,
    });
    if storage_mode == "dense" {
        attributes["primary_array"] = json!("counts");
    } else if storage_mode == "sparse_coo" {
        attributes["primary_group"] = json!("sparse");
        attributes["sparse_format"] = json!("coo");
        attributes["sparse_indices_base"] = json!(0);
    }
    write_zarr_root_metadata(store, "end-motif", attributes)
}

/// Write motif-axis metadata.
///
/// `motif_index` is the numeric count-column coordinate. The labels are stored as
/// `motif_ascii[motif, motif_byte]`, one ASCII byte per character. End-motif labels are fixed-width
/// for one run, so this avoids variable-length string support and avoids repeating labels in JSON.
fn write_motif_metadata(store: Arc<FilesystemStore>, motifs: &[String]) -> Result<()> {
    let motif_width = validated_motif_width(motifs)?;
    let motif_axis = checked_index_axis(motifs.len(), "motif")?;
    let motif_byte_axis = checked_index_axis(motif_width, "motif_byte")?;
    let motif_ascii = encode_motif_ascii(motifs, motif_width);

    write_single_chunk_zarr_array(
        store.clone(),
        "motif_index",
        &[motifs.len()],
        &["motif"],
        &motif_axis,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based motif column index",
            "label_array": "motif_ascii",
        }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "motif_byte",
        &[motif_width],
        &["motif_byte"],
        &motif_byte_axis,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "zero-based byte offset within motif labels" }),
    )?;
    write_single_chunk_zarr_array(
        store,
        "motif_ascii",
        &[motifs.len(), motif_width],
        &["motif", "motif_byte"],
        &motif_ascii,
        data_type::uint8(),
        ZARR_ASCII_FILL_VALUE,
        json!({
            "long_name": "fixed-width ASCII motif labels",
            "description": "Decode each [motif, motif_byte] row as ASCII to recover the motif label.",
        }),
    )?;
    Ok(())
}

/// Return the fixed byte width shared by all motif labels.
fn validated_motif_width(motifs: &[String]) -> Result<usize> {
    let Some(first_motif) = motifs.first() else {
        return Ok(0);
    };
    validate_motif_ascii(first_motif)?;
    let motif_width = first_motif.len();
    for motif in &motifs[1..] {
        validate_motif_ascii(motif)?;
        ensure!(
            motif.len() == motif_width,
            "end-motif labels must have one fixed ASCII width, but '{}' has {} bytes and expected {}",
            motif,
            motif.len(),
            motif_width
        );
    }
    Ok(motif_width)
}

/// Validate that a motif label can be stored as fixed-width ASCII bytes.
fn validate_motif_ascii(motif: &str) -> Result<()> {
    validate_zarr_label(motif, "motif")?;
    ensure!(
        motif.is_ascii(),
        "motif labels must be ASCII to be stored in motif_ascii"
    );
    Ok(())
}

/// Encode motif labels in row-major `[motif, motif_byte]` order.
fn encode_motif_ascii(motifs: &[String], motif_width: usize) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(motifs.len() * motif_width);
    for motif in motifs {
        bytes.extend_from_slice(motif.as_bytes());
    }
    bytes
}

/// Write row-axis metadata for the selected output mode.
///
/// Every end-motif output has a public `row` coordinate. The additional arrays differ by mode:
/// global output has one row label, windowed output has genomic coordinates, and grouped output
/// has group names and group-level eligibility metadata.
fn write_row_metadata(
    store: Arc<FilesystemStore>,
    row_metadata: EndMotifRowMetadata<'_>,
    n_rows: usize,
) -> Result<()> {
    ensure!(
        n_rows > 0,
        "end-motif Zarr output requires at least one row"
    );
    let row = checked_index_axis(n_rows, "row")?;
    let mut row_attributes = json!({ "long_name": "zero-based count row index" });
    if matches!(&row_metadata, EndMotifRowMetadata::Global) {
        row_attributes["label_field"] = json!("row_label");
        row_attributes["labels"] = json!(["global"]);
    }
    write_single_chunk_zarr_array(
        store.clone(),
        "row",
        &[n_rows],
        &["row"],
        &row,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        row_attributes,
    )?;

    match row_metadata {
        EndMotifRowMetadata::Global => {
            ensure!(n_rows == 1, "global end-motif output must have one row");
        }
        EndMotifRowMetadata::Windows { bin_info, .. } => {
            write_window_row_metadata(store, bin_info, n_rows)?;
        }
        EndMotifRowMetadata::Groups(groups) => {
            write_group_row_metadata(store, &groups, n_rows)?;
        }
    }
    Ok(())
}

/// Write genomic row metadata for size/bed window outputs.
///
/// Chromosome names are dictionary-encoded as a small chromosome axis plus `row_chromosome`, so
/// repeated chromosome names are not duplicated for every row.
fn write_window_row_metadata(
    store: Arc<FilesystemStore>,
    bin_info: &[WindowBinInfo],
    n_rows: usize,
) -> Result<()> {
    ensure!(
        bin_info.len() == n_rows,
        "end-motif window metadata rows ({}) did not match count rows ({})",
        bin_info.len(),
        n_rows
    );

    let chromosome_names = chromosome_axis_from_windows(bin_info)?;
    let chromosome_name_refs: Vec<&str> = chromosome_names.iter().map(String::as_str).collect();
    let chromosome_axis = checked_index_axis(chromosome_names.len(), "chromosome")?;
    let chromosome_lookup: FxHashMap<&str, usize> = chromosome_names
        .iter()
        .enumerate()
        .map(|(index, name)| (name.as_str(), index))
        .collect();
    let row_chromosome = bin_info
        .iter()
        .map(|entry| {
            let index = chromosome_lookup
                .get(entry.chromosome.as_str())
                .copied()
                .context("window chromosome missing from chromosome axis")?;
            checked_i32(index, "row_chromosome")
        })
        .collect::<Result<Vec<_>>>()?;
    let row_start_bp: Vec<i64> = bin_info
        .iter()
        .map(|entry| checked_i64(entry.start, "row_start_bp"))
        .collect::<Result<_>>()?;
    let row_end_bp: Vec<i64> = bin_info
        .iter()
        .map(|entry| checked_i64(entry.end, "row_end_bp"))
        .collect::<Result<_>>()?;
    let blacklisted_fraction: Vec<f64> = bin_info
        .iter()
        .map(|entry| entry.blacklisted_fraction)
        .collect();
    write_single_chunk_zarr_array(
        store.clone(),
        "chromosome",
        &[chromosome_names.len()],
        &["chromosome"],
        &chromosome_axis,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based chromosome index",
            "label_field": "chromosome_name",
            "labels": chromosome_name_refs,
        }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "row_chromosome",
        &[n_rows],
        &["row"],
        &row_chromosome,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "chromosome index for each count row" }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "row_start_bp",
        &[n_rows],
        &["row"],
        &row_start_bp,
        data_type::int64(),
        ZARR_INT64_FILL_VALUE,
        json!({ "long_name": "inclusive row start coordinate", "units": "bp" }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "row_end_bp",
        &[n_rows],
        &["row"],
        &row_end_bp,
        data_type::int64(),
        ZARR_INT64_FILL_VALUE,
        json!({ "long_name": "exclusive row end coordinate", "units": "bp" }),
    )?;
    write_single_chunk_zarr_array(
        store,
        "blacklisted_fraction",
        &[n_rows],
        &["row"],
        &blacklisted_fraction,
        data_type::float64(),
        ZARR_FLOAT64_FILL_VALUE,
        json!({ "long_name": "fraction of row coordinates overlapping blacklist intervals" }),
    )?;

    Ok(())
}

/// Write grouped-BED row metadata.
///
/// Group rows are already expected to be in count-row order when this function is called. The
/// writer stores the numeric `group` coordinate and keeps the human-readable group names in its
/// JSON attributes so downstream readers can label rows without opening a string array.
fn write_group_row_metadata(
    store: Arc<FilesystemStore>,
    groups: &[EndGroupSummary<'_>],
    n_rows: usize,
) -> Result<()> {
    ensure!(
        groups.len() == n_rows,
        "end-motif group metadata rows ({}) did not match count rows ({})",
        groups.len(),
        n_rows
    );
    for group in groups {
        validate_zarr_label(group.group_name, "group_name")?;
    }
    let group: Vec<i32> = groups
        .iter()
        .map(|group| checked_i32(group.group_idx, "group_idx"))
        .collect::<Result<Vec<_>>>()?;
    let group_name_refs: Vec<&str> = groups.iter().map(|group| group.group_name).collect();
    let eligible_windows: Vec<i32> = groups
        .iter()
        .map(|group| checked_i32(group.eligible_windows, "eligible_windows"))
        .collect::<Result<Vec<_>>>()?;
    let blacklisted_fraction: Vec<f64> = groups
        .iter()
        .map(|group| group.blacklisted_fraction)
        .collect();

    write_single_chunk_zarr_array(
        store.clone(),
        "group",
        &[n_rows],
        &["row"],
        &group,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based group index",
            "description": "For grouped BED output, this matches the count row index.",
            "label_field": "group_name",
            "labels": group_name_refs,
        }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "eligible_windows",
        &[n_rows],
        &["row"],
        &eligible_windows,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "eligible grouped BED windows retained per group" }),
    )?;
    write_single_chunk_zarr_array(
        store,
        "blacklisted_fraction",
        &[n_rows],
        &["row"],
        &blacklisted_fraction,
        data_type::float64(),
        ZARR_FLOAT64_FILL_VALUE,
        json!({ "long_name": "length-weighted blacklisted fraction across group windows" }),
    )?;

    Ok(())
}

/// Write dense `counts[row, motif]`.
///
/// The full logical array is handed to `zarrs::store_array_subset`; zarrs then handles chunk
/// splitting and boundary chunks. This avoids manual padding logic in the command writer.
fn write_dense_counts(store: Arc<FilesystemStore>, counts: ArrayView2<'_, f64>) -> Result<()> {
    ensure!(
        counts.shape().iter().all(|dimension| *dimension > 0),
        "dense end-motif counts cannot have an empty axis"
    );
    let shape = [counts.shape()[0], counts.shape()[1]];
    let chunk_shape = dense_count_chunk_shape(shape)?;
    let array = create_zarr_array(
        store,
        "counts",
        &shape,
        &chunk_shape,
        &["row", "motif"],
        data_type::float64(),
        ZARR_FLOAT64_FILL_VALUE,
        json!({
            "long_name": "weighted end-motif count",
            "units": "weighted_end_motif_count",
        }),
    )?;
    let owned_counts = counts.as_standard_layout();
    let values = owned_counts
        .as_slice()
        .context("standard-layout end-motif counts were not contiguous")?;
    array
        .store_array_subset(&array.subset_all(), values)
        .context("write dense end-motif Zarr counts")?;
    Ok(())
}

/// Choose a dense count chunk shape close to the target cell count.
///
/// Motif columns are kept together first because common downstream reads are expected to load all
/// motifs for a subset of rows.
fn dense_count_chunk_shape(shape: [usize; 2]) -> Result<[usize; 2]> {
    let total_cells = shape
        .iter()
        .try_fold(1usize, |cells, dimension| cells.checked_mul(*dimension))
        .context("end-motif dense Zarr count shape overflow")?;
    if total_cells <= TARGET_DENSE_COUNT_CHUNK_CELLS {
        return Ok(shape);
    }

    let motif_chunk = shape[1].min(TARGET_DENSE_COUNT_CHUNK_CELLS).max(1);
    let row_chunk = shape[0]
        .min((TARGET_DENSE_COUNT_CHUNK_CELLS / motif_chunk).max(1))
        .max(1);
    Ok([row_chunk, motif_chunk])
}

/// Write sparse counts in sorted COO form.
///
/// The sparse group contains parallel `row`, `motif`, and `count` arrays plus a `shape` array. The
/// entries are sorted by `(row, motif)` so downstream readers can reconstruct deterministic dense
/// matrices or sparse objects.
fn write_sparse_counts(
    store: Arc<FilesystemStore>,
    bins: &[FxHashMap<String, f64>],
    motifs: &[String],
) -> Result<()> {
    write_zarr_group_metadata(
        store.clone(),
        "/sparse",
        "end-motif sparse counts",
        json!({
            "long_name": "Sparse COO end-motif count arrays",
            "sparse_format": "coo",
            "sparse_indices_base": 0,
        }),
    )?;

    let motif_columns: FxHashMap<&String, i32> = motifs
        .iter()
        .enumerate()
        .map(|(column, motif)| Ok((motif, checked_i32(column, "motif index")?)))
        .collect::<Result<_>>()?;
    let mut entries = Vec::new();
    for (row, bin) in bins.iter().enumerate() {
        let row_index = checked_i32(row, "row index")?;
        for (motif, &count) in bin {
            let motif_index = motif_columns.get(motif).copied().with_context(|| {
                format!("missing sparse output column for end-motif label '{motif}'")
            })?;
            if count != 0.0 {
                entries.push((row_index, motif_index, count));
            }
        }
    }
    entries.sort_unstable_by(|left, right| (left.0, left.1).cmp(&(right.0, right.1)));
    for pair in entries.windows(2) {
        ensure!(
            (pair[0].0, pair[0].1) != (pair[1].0, pair[1].1),
            "duplicate sparse end-motif COO entry at row {} motif {}",
            pair[0].0,
            pair[0].1
        );
    }

    let row: Vec<i32> = entries.iter().map(|entry| entry.0).collect();
    let motif: Vec<i32> = entries.iter().map(|entry| entry.1).collect();
    let count: Vec<f64> = entries.iter().map(|entry| entry.2).collect();
    let shape = vec![
        checked_i32(bins.len(), "sparse row count")?,
        checked_i32(motifs.len(), "sparse motif count")?,
    ];
    let sparse_dimension = checked_index_axis(2, "sparse_dimension")?;
    let nnz = entries.len();

    write_single_chunk_zarr_array(
        store.clone(),
        "sparse/row",
        &[nnz],
        &["nnz"],
        &row,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "COO row index" }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "sparse/motif",
        &[nnz],
        &["nnz"],
        &motif,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "COO motif index" }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "sparse/count",
        &[nnz],
        &["nnz"],
        &count,
        data_type::float64(),
        ZARR_FLOAT64_FILL_VALUE,
        json!({
            "long_name": "weighted end-motif count",
            "units": "weighted_end_motif_count",
        }),
    )?;
    write_single_chunk_zarr_array(
        store.clone(),
        "sparse/shape",
        &[2],
        &["sparse_dimension"],
        &shape,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({ "long_name": "dense shape represented by sparse COO arrays" }),
    )?;
    write_single_chunk_zarr_array(
        store,
        "sparse/sparse_dimension",
        &[2],
        &["sparse_dimension"],
        &sparse_dimension,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based sparse shape dimension index",
            "label_field": "sparse_dimension_name",
            "labels": ["row", "motif"],
        }),
    )?;
    Ok(())
}

/// Return chromosome labels in first-seen row order.
///
/// Window rows store chromosome references as integer coordinates into this axis. Preserving
/// first-seen order keeps the mapping stable with the input row order while avoiding a lexical sort
/// that would be surprising for chromosome names such as `chr10`.
fn chromosome_axis_from_windows(bin_info: &[WindowBinInfo]) -> Result<Vec<String>> {
    let mut chromosome_names = Vec::new();
    let mut seen = FxHashMap::default();
    for entry in bin_info {
        validate_zarr_label(&entry.chromosome, "chromosome_name")?;
        if !seen.contains_key(entry.chromosome.as_str()) {
            seen.insert(entry.chromosome.as_str(), chromosome_names.len());
            chromosome_names.push(entry.chromosome.clone());
        }
    }
    ensure!(
        !chromosome_names.is_empty(),
        "windowed end-motif output requires at least one chromosome"
    );
    Ok(chromosome_names)
}

#[cfg(test)]
mod tests {
    include!("zarr_tests.rs");
}
