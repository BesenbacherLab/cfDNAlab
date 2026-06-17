//! End-motif Zarr writer.
//!
//! `zarrs` owns the Zarr metadata, chunk encoding, and compression. This module owns the cfDNAlab
//! schema:
//!
//! - `counts[row, motif]` for dense output
//! - `sparse/{row,motif,count,shape,sparse_dimension}` for sparse COO output
//! - motif or motif-group labels and row metadata needed to interpret the counts without joining TSV files
//!
//! The code below should only do work that is specific to this schema or to clearer validation.
//! Low-level details such as writing `zarr.json`, applying zstd, V3 dimension names, and chunk
//! serialization are delegated to `zarrs`.

use crate::{
    commands::ends::counting::EndMotifColumnKind,
    shared::{
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

const CFDNALAB_END_MOTIF_SCHEMA_VERSION: u32 = 2;

/// Soft target for dense count chunks.
///
/// This keeps dense `f64` chunks near sixteen MiB. Small dense matrices still use one chunk.
const TARGET_DENSE_COUNT_CHUNK_CELLS: usize = 2_000_000;

/// Row metadata for the public `row` axis.
///
/// End-motif counts always have a row axis, but the meaning of a row depends on the windowing mode.
/// This enum carries the extra metadata needed to make that row axis self-describing in the Zarr
/// store.
pub(crate) enum EndMotifRowMetadata<'a> {
    /// One row covering all selected chromosomes.
    Global,
    /// Genomic window rows in the same row order as the count matrix.
    Windows {
        /// Per-row chromosome, coordinate, and blacklist metadata
        bin_info: &'a [WindowBinInfo],
        /// Whether these rows came from size windows or BED windows
        row_mode: EndWindowRowMode,
    },
    /// Grouped BED rows in the same row order as the count matrix.
    Groups(Vec<EndGroupSummary<'a>>),
}

/// Source of genomic window count rows.
///
/// This is written into root metadata as `row_mode`, so downstream readers can distinguish
/// generated size windows from user-provided BED windows without inspecting coordinate spacing.
pub(crate) enum EndWindowRowMode {
    /// Generated fixed-size windows
    Size,
    /// User-provided BED windows
    Bed,
}

/// One grouped-BED output row.
///
/// Grouped BED output aggregates many genomic windows into one count row. The summary keeps the
/// public group label and enough metadata for downstream code to explain how much sequence
/// contributed to that row.
#[derive(Debug)]
pub(crate) struct EndGroupSummary<'a> {
    /// Zero-based group index, expected to match the count row index
    pub(crate) group_idx: u64,
    /// Public group label from the grouped BED input
    pub(crate) group_name: &'a str,
    /// Number of grouped BED windows contributing to the group
    pub(crate) eligible_windows: usize,
    /// Length-weighted blacklist overlap across the group's windows
    pub(crate) blacklisted_fraction: f64,
}

/// Write dense or sparse end-motif counts as a self-contained Zarr V3 store.
///
/// `bins` are already keyed by numeric column indices. The labels and `column_kind` describe how
/// those indices should be exposed to downstream readers. Concrete motif outputs use fixed-width
/// ASCII motif labels, while motif-group outputs store group labels as JSON metadata on
/// `motif_index`.
///
/// Parameters
/// ----------
/// - `output_dir`:
///   Directory where the Zarr store should be created
/// - `prefix`:
///   Output filename prefix
/// - `bins`:
///   Sparse row maps keyed by final zero-based motif or motif-group column index
/// - `column_labels`:
///   Public labels for the column axis in final column order
/// - `column_kind`:
///   Meaning of `column_labels`
/// - `row_metadata`:
///   Metadata for the public row axis
/// - `write_dense_output`:
///   Whether to write `counts[row, motif]` instead of sparse COO arrays
///
/// Returns
/// -------
/// - `Result<PathBuf>`:
///   Path to the completed temporary Zarr store
pub(crate) fn write_end_motif_zarr(
    output_dir: &Path,
    prefix: &str,
    bins: &[FxHashMap<u32, f64>],
    column_labels: &[String],
    column_kind: EndMotifColumnKind,
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
    write_root_metadata(store.clone(), storage_mode, row_mode, column_kind)?;
    write_motif_metadata(store.clone(), column_labels, column_kind)?;
    write_row_metadata(store.clone(), row_metadata, bins.len())?;

    if write_dense_output {
        let counts = stack_end_motif_counts(bins, column_labels.len())?;
        write_dense_counts(store, counts.view())?;
    } else {
        write_sparse_counts(store, bins, column_labels.len())?;
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

/// Stack sparse per-window motif maps into a dense matrix.
///
/// Sparse bins use numeric column ids rather than labels, so this function only has to validate
/// bounds and place values. Label interpretation stays in the metadata writer.
///
/// Parameters
/// ----------
/// - `bins`:
///   Sparse row maps keyed by final output column index
/// - `n_columns`:
///   Width of the final motif or motif-group axis
///
/// Returns
/// -------
/// - `Result<Array2<f64>>`:
///   Dense row-major count matrix
fn stack_end_motif_counts(bins: &[FxHashMap<u32, f64>], n_columns: usize) -> Result<Array2<f64>> {
    let mut counts = Array2::<f64>::zeros((bins.len(), n_columns));

    for (row, bin) in bins.iter().enumerate() {
        for (&target_idx, &count) in bin {
            // The parser and postprocessor should only produce valid final columns
            let column = target_idx as usize;
            ensure!(
                column < n_columns,
                "end-motif column index {} is out of bounds for {} columns",
                column,
                n_columns
            );
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
    column_kind: EndMotifColumnKind,
) -> Result<()> {
    let mut attributes = json!({
        "cfdnalab_schema": "end_motif_counts",
        "cfdnalab_schema_version": CFDNALAB_END_MOTIF_SCHEMA_VERSION,
        "storage_mode": storage_mode,
        "row_mode": row_mode,
        "motif_axis_kind": end_motif_column_kind_name(column_kind),
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

/// Return the schema string for the motif-axis kind.
///
/// This attribute tells downstream readers whether `motif_index` labels should be exposed as
/// motifs or motif groups.
fn end_motif_column_kind_name(column_kind: EndMotifColumnKind) -> &'static str {
    match column_kind {
        EndMotifColumnKind::Motif => "motif",
        EndMotifColumnKind::MotifGroup => "motif_group",
    }
}

/// Write motif-axis metadata.
///
/// `motif_index` is the numeric count-column coordinate. Ungrouped motif labels are stored as
/// fixed-width ASCII bytes. Grouped motif-file outputs store group labels in JSON attributes.
///
/// Parameters
/// ----------
/// - `store`:
///   Open Zarr store
/// - `labels`:
///   Final column labels in count-column order
/// - `column_kind`:
///   Whether labels are motifs or motif groups
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after the motif axis metadata has been written
fn write_motif_metadata(
    store: Arc<FilesystemStore>,
    labels: &[String],
    column_kind: EndMotifColumnKind,
) -> Result<()> {
    match column_kind {
        EndMotifColumnKind::Motif => write_motif_label_metadata(store, labels),
        EndMotifColumnKind::MotifGroup => write_motif_group_metadata(store, labels),
    }
}

/// Write metadata for concrete motif columns.
///
/// Concrete motifs have one fixed byte width within a run, so they are stored as an ASCII matrix
/// instead of a variable-length string array. This keeps the schema simple for readers that do not
/// support string arrays well.
fn write_motif_label_metadata(store: Arc<FilesystemStore>, motifs: &[String]) -> Result<()> {
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

/// Write metadata for grouped motifs-file columns.
///
/// Group labels can have different lengths, so they are written as JSON labels on `motif_index`
/// rather than through `motif_ascii`. Downstream readers can branch on `motif_axis_kind` and expose
/// these as motif-group labels without pretending they are DNA motifs.
fn write_motif_group_metadata(store: Arc<FilesystemStore>, motif_groups: &[String]) -> Result<()> {
    for motif_group in motif_groups {
        // Group names are public labels and must be safe in JSON metadata and data frame columns
        validate_zarr_label(motif_group, "motif_group")?;
    }
    let motif_axis = checked_index_axis(motif_groups.len(), "motif")?;
    let motif_group_refs: Vec<&str> = motif_groups.iter().map(String::as_str).collect();

    write_single_chunk_zarr_array(
        store,
        "motif_index",
        &[motif_groups.len()],
        &["motif"],
        &motif_axis,
        data_type::int32(),
        ZARR_INT32_FILL_VALUE,
        json!({
            "long_name": "zero-based motif group column index",
            "label_field": "motif_group",
            "labels": motif_group_refs,
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
/// The full logical array is handed to `zarrs::store_array_subset`. zarrs then handles chunk
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
///
/// Parameters
/// ----------
/// - `store`:
///   Open Zarr store
/// - `bins`:
///   Sparse row maps keyed by final output column index
/// - `n_columns`:
///   Width of the final motif or motif-group axis
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after the sparse group has been written
fn write_sparse_counts(
    store: Arc<FilesystemStore>,
    bins: &[FxHashMap<u32, f64>],
    n_columns: usize,
) -> Result<()> {
    let mut entries = Vec::new();
    for (row, bin) in bins.iter().enumerate() {
        let row_index = checked_i32(row, "row index")?;
        for (&target_idx, &count) in bin {
            // Bounds are checked here before entries are converted to i32 COO coordinates
            let column = target_idx as usize;
            ensure!(
                column < n_columns,
                "sparse end-motif column index {} is out of bounds for {} columns",
                column,
                n_columns
            );
            if count != 0.0 {
                entries.push((row_index, checked_i32(target_idx, "motif index")?, count));
            }
        }
    }
    write_sparse_entries(store, bins.len(), n_columns, entries)
}

/// Write prevalidated COO entries to the sparse group.
///
/// This helper is deliberately label-agnostic. Its contract is only that the rows and columns are
/// already numeric coordinates into the public axes written elsewhere in the store.
///
/// Parameters
/// ----------
/// - `store`:
///   Open Zarr store
/// - `n_rows`:
///   Number of rows in the represented dense matrix
/// - `n_motifs`:
///   Number of motif-axis columns in the represented dense matrix
/// - `entries`:
///   Sparse COO entries as `(row, motif, count)`
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after all sparse arrays have been written
fn write_sparse_entries(
    store: Arc<FilesystemStore>,
    n_rows: usize,
    n_motifs: usize,
    mut entries: Vec<(i32, i32, f64)>,
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

    // Sorting creates deterministic arrays and lets us reject duplicate coordinates cheaply
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
        checked_i32(n_rows, "sparse row count")?,
        checked_i32(n_motifs, "sparse motif count")?,
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
