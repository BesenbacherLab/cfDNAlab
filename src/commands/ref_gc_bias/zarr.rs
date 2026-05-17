//! Zarr package I/O for `ref-gc-bias`.
//!
//! This module owns only the public package container. The scientific arrays are
//! produced by `ref_gc_bias.rs` before they get here, and the loader returns the
//! same in-memory structures used by `gc-bias`.

use crate::{
    commands::gc_bias::load_reference_bias::{ReferenceGCData, ReferenceGCMetadata},
    shared::{
        constants::{GC_CORRECTION_SCHEMA_VERSION, MIN_ACGT_BASES_FOR_GC_FRACTION},
        reference::ContigFootprintEntry,
        zarr::{
            bool_fill_value, checked_index_axis, create_zarr_array_with_fill_value,
            create_zarr_store, ensure_zarr_schema, read_zarr_array1, read_zarr_array2,
            read_zarr_root_attributes, validate_zarr_label, write_single_chunk_zarr_array,
            write_zarr_root_metadata,
        },
    },
};
use anyhow::{Context, Result, ensure};
use ndarray::Array2;
use serde_json::{Value, json};
use std::{path::Path, sync::Arc};
use zarrs::{array::data_type, filesystem::FilesystemStore};

const REFERENCE_GC_SCHEMA: &str = "reference_gc_package";

/// Inputs needed to write a reference GC Zarr package.
pub struct ReferenceGCZarrPackage<'a> {
    pub counts: &'a Array2<f64>,
    pub support_unobservables: &'a Array2<bool>,
    pub support_outliers: &'a Array2<bool>,
    pub gc_percent_widths: &'a Array2<u16>,
    pub length_min: usize,
    pub length_max: usize,
    pub end_offset: u8,
    pub skip_interpolation: bool,
    pub smoothing_radius: u8,
    pub smoothing_sigma: f64,
    pub skip_smoothing: bool,
    pub chromosomes: &'a [String],
    pub reference_contig_footprint: &'a [ContigFootprintEntry],
}

/// Write the public reference GC package as a Zarr V3 store.
pub fn write_reference_gc_package_zarr(
    store_path: &Path,
    package: ReferenceGCZarrPackage<'_>,
) -> Result<()> {
    validate_reference_package_for_writing(&package)?;
    let store = create_zarr_store(store_path, "reference GC package")?;
    write_zarr_root_metadata(
        store.clone(),
        "reference GC package",
        json!({
            "cfdnalab_schema": REFERENCE_GC_SCHEMA,
            "cfdnalab_schema_version": GC_CORRECTION_SCHEMA_VERSION,
            "package_role": "reference_gc",
            "value_units": "reference_fragment_mass",
            "gc_percent_rounding": "integer_half_up",
            "minimum_acgt_bases_for_gc_fraction": MIN_ACGT_BASES_FOR_GC_FRACTION,
            "end_offset": package.end_offset,
            "skip_interpolation": package.skip_interpolation,
            "smoothing_radius": package.smoothing_radius,
            "smoothing_sigma": package.smoothing_sigma,
            "skip_smoothing": package.skip_smoothing,
        }),
    )?;

    let shape = [package.counts.nrows(), package.counts.ncols()];
    write_matrix(
        store.clone(),
        "counts",
        &shape,
        package.counts,
        data_type::float64(),
        0.0,
        json!({"long_name": "reference GC fragment mass"}),
    )?;
    write_bool_matrix(
        store.clone(),
        "support_mask_unobservables",
        &shape,
        package.support_unobservables,
        data_type::bool(),
        json!({"long_name": "theoretical reference support mask"}),
    )?;
    write_bool_matrix(
        store.clone(),
        "support_mask_outliers",
        &shape,
        package.support_outliers,
        data_type::bool(),
        json!({"long_name": "empirical reference support mask"}),
    )?;
    write_matrix(
        store.clone(),
        "gc_percent_widths",
        &shape,
        package.gc_percent_widths,
        data_type::uint16(),
        0u16,
        json!({"long_name": "number of integer GC counts represented by each GC percent bin"}),
    )?;

    let lengths = (package.length_min..=package.length_max)
        .map(|length| {
            i32::try_from(length).with_context(|| format!("fragment length {length} exceeds i32"))
        })
        .collect::<Result<Vec<_>>>()?;
    write_single_chunk_zarr_array(
        store.clone(),
        "length",
        &[lengths.len()],
        &["length"],
        &lengths,
        data_type::int32(),
        0,
        json!({"long_name": "fragment length in bp"}),
    )?;

    let gc_percent = checked_index_axis(package.counts.ncols(), "gc_percent")?;
    write_single_chunk_zarr_array(
        store.clone(),
        "gc_percent",
        &[gc_percent.len()],
        &["gc_percent"],
        &gc_percent,
        data_type::int32(),
        0,
        json!({"long_name": "integer GC percent"}),
    )?;

    for chromosome in package.chromosomes {
        validate_zarr_label(chromosome, "chromosome_name")?;
    }
    let chromosome = checked_index_axis(package.chromosomes.len(), "chromosome")?;
    write_single_chunk_zarr_array(
        store.clone(),
        "chromosome",
        &[chromosome.len()],
        &["chromosome"],
        &chromosome,
        data_type::int32(),
        0,
        json!({
            "long_name": "selected chromosome index",
            "label_field": "chromosome_name",
            "labels": package.chromosomes,
        }),
    )?;

    let reference_contig_footprint_json = serde_json::to_vec(package.reference_contig_footprint)?;
    write_single_chunk_zarr_array(
        store,
        "reference_contig_footprint_json",
        &[reference_contig_footprint_json.len()],
        &["json_byte"],
        &reference_contig_footprint_json,
        data_type::uint8(),
        0u8,
        json!({"long_name": "JSON-encoded reference contig footprint"}),
    )?;

    Ok(())
}

/// Read a reference GC package from a Zarr V3 store.
pub(crate) fn read_reference_gc_package_zarr(path: &Path) -> Result<ReferenceGCData> {
    ensure!(
        path.is_dir(),
        "Reference GC package path must point to an existing .zarr directory: {}",
        path.display()
    );
    ensure!(
        path.extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| extension.eq_ignore_ascii_case("zarr")),
        "Reference GC package path must point to a .zarr directory: {}",
        path.display()
    );
    let root = read_zarr_root_attributes(path)?;
    ensure_zarr_schema(
        &root,
        REFERENCE_GC_SCHEMA,
        GC_CORRECTION_SCHEMA_VERSION,
        "Reference GC package",
    )?;

    let store = Arc::new(FilesystemStore::new(path)?);
    let counts = read_zarr_array2::<f64>(store.clone(), "/counts")?;
    let unobservables_support_mask =
        read_zarr_array2::<bool>(store.clone(), "/support_mask_unobservables")?;
    let outliers_support_mask = read_zarr_array2::<bool>(store.clone(), "/support_mask_outliers")?;
    let gc_percent_widths = read_zarr_array2::<u16>(store.clone(), "/gc_percent_widths")?;
    let lengths = read_zarr_array1::<i32>(store.clone(), "/length")?;
    let gc_percent = read_zarr_array1::<i32>(store.clone(), "/gc_percent")?;
    let chromosome_axis = read_zarr_array1::<i32>(store.clone(), "/chromosome")?;
    let reference_contig_footprint_json =
        read_zarr_array1::<u8>(store, "/reference_contig_footprint_json")?;

    ensure!(
        !lengths.is_empty(),
        "Reference GC package length axis must not be empty"
    );
    ensure_contiguous_values(&lengths, "length")?;
    ensure_index_axis(&gc_percent, "gc_percent")?;
    ensure_index_axis(&chromosome_axis, "chromosome")?;

    let chromosomes = read_labels(path, "chromosome", "chromosome_name")?;
    ensure!(
        chromosomes.len() == chromosome_axis.len(),
        "chromosome labels length ({}) does not match chromosome axis length ({})",
        chromosomes.len(),
        chromosome_axis.len()
    );
    ensure!(
        !chromosomes.is_empty(),
        "chromosomes must contain at least one chromosome"
    );

    let min_fragment_length = usize::try_from(lengths[0]).context("length must be non-negative")?;
    let max_fragment_length =
        usize::try_from(*lengths.last().expect("length axis checked non-empty"))
            .context("length must be non-negative")?;
    let end_offset = read_u8_attr(&root, "end_offset")?;
    let smoothing_radius = read_u8_attr(&root, "smoothing_radius")?;
    let smoothing_sigma = read_f64_attr(&root, "smoothing_sigma")?;
    let skip_smoothing = read_bool_attr(&root, "skip_smoothing")?;

    ensure_reference_metadata(
        &counts,
        &unobservables_support_mask,
        &outliers_support_mask,
        &gc_percent_widths,
        min_fragment_length,
        max_fragment_length,
        end_offset,
        skip_smoothing,
        smoothing_radius,
        smoothing_sigma,
    )?;

    let reference_contig_footprint: Vec<ContigFootprintEntry> =
        serde_json::from_slice(&reference_contig_footprint_json)
            .context("invalid reference_contig_footprint_json in reference GC package")?;

    Ok(ReferenceGCData {
        counts,
        unobservables_support_mask,
        outliers_support_mask,
        gc_percent_widths,
        metadata: ReferenceGCMetadata {
            min_fragment_length,
            max_fragment_length,
            end_offset,
            chromosomes,
            reference_contig_footprint,
            skip_interpolation: read_bool_attr(&root, "skip_interpolation")?,
            smoothing_sigma,
            smoothing_radius,
            skip_smoothing,
        },
    })
}

fn validate_reference_package_for_writing(package: &ReferenceGCZarrPackage<'_>) -> Result<()> {
    let expected_rows = package
        .length_max
        .checked_sub(package.length_min)
        .context("length_min must be <= length_max")?
        + 1;
    ensure!(
        package.counts.nrows() == expected_rows,
        "counts rows ({}) do not match length range [{}, {}] ({})",
        package.counts.nrows(),
        package.length_min,
        package.length_max,
        expected_rows
    );
    ensure!(
        package.counts.ncols() > 0,
        "counts must contain at least one GC percent column"
    );
    ensure!(
        package.support_unobservables.dim() == package.counts.dim(),
        "support_mask_unobservables shape {:?} does not match counts shape {:?}",
        package.support_unobservables.dim(),
        package.counts.dim()
    );
    ensure!(
        package.support_outliers.dim() == package.counts.dim(),
        "support_mask_outliers shape {:?} does not match counts shape {:?}",
        package.support_outliers.dim(),
        package.counts.dim()
    );
    ensure!(
        package.gc_percent_widths.dim() == package.counts.dim(),
        "gc_percent_widths shape {:?} does not match counts shape {:?}",
        package.gc_percent_widths.dim(),
        package.counts.dim()
    );
    Ok(())
}

fn ensure_reference_metadata(
    counts: &Array2<f64>,
    unobservables_support_mask: &Array2<bool>,
    outliers_support_mask: &Array2<bool>,
    gc_percent_widths: &Array2<u16>,
    min_fragment_length: usize,
    max_fragment_length: usize,
    end_offset: u8,
    skip_smoothing: bool,
    smoothing_radius: u8,
    smoothing_sigma: f64,
) -> Result<()> {
    ensure!(
        min_fragment_length <= max_fragment_length,
        "length axis must be ordered. Found [{}, {}]",
        min_fragment_length,
        max_fragment_length
    );
    let expected_length_rows = max_fragment_length - min_fragment_length + 1;
    ensure!(
        counts.nrows() == expected_length_rows,
        "Reference GC package row count {} does not match length axis [{}, {}] (expected {})",
        counts.nrows(),
        min_fragment_length,
        max_fragment_length,
        expected_length_rows
    );
    let minimum_effective_length = MIN_ACGT_BASES_FOR_GC_FRACTION as usize;
    ensure!(
        min_fragment_length >= 2 * end_offset as usize + minimum_effective_length,
        "Reference GC package has invalid effective minimum length: min_fragment_length ({}) - 2 * end_offset ({}) must be >= {}",
        min_fragment_length,
        end_offset,
        minimum_effective_length
    );
    if !skip_smoothing {
        ensure!(
            smoothing_sigma.is_finite() && smoothing_sigma > 0.0,
            "smoothing_sigma in reference GC package must be finite and > 0 when smoothing is enabled"
        );
        ensure!(
            smoothing_radius > 0,
            "smoothing_radius in reference GC package must be > 0 when smoothing is enabled"
        );
    }
    ensure!(
        unobservables_support_mask.dim() == outliers_support_mask.dim(),
        "The two support masks must have the same shape. Unobservables ({:?}) != outliers ({:?}).",
        unobservables_support_mask.dim(),
        outliers_support_mask.dim(),
    );
    ensure!(
        unobservables_support_mask.dim() == counts.dim(),
        "Reference counts ({:?}) and support masks ({:?}) had incompatible shapes",
        counts.dim(),
        unobservables_support_mask.dim()
    );
    ensure!(
        gc_percent_widths.dim() == counts.dim(),
        "GC percent widths shape {:?} must match per-window counts shape {:?}",
        gc_percent_widths.dim(),
        counts.dim(),
    );
    Ok(())
}

fn write_matrix<T>(
    store: Arc<FilesystemStore>,
    name: &str,
    shape: &[usize; 2],
    values: &Array2<T>,
    data_type: zarrs::array::DataType,
    fill_value: T,
    attributes: Value,
) -> Result<()>
where
    T: zarrs::array::Element + Copy,
    zarrs::array::builder::ArrayBuilderFillValue: From<T>,
{
    let values = values.as_standard_layout();
    let slice = values
        .as_slice()
        .with_context(|| format!("{name} matrix should be contiguous after layout conversion"))?;
    write_single_chunk_zarr_array(
        store,
        name,
        shape,
        &["length", "gc_percent"],
        slice,
        data_type,
        fill_value,
        attributes,
    )
}

fn write_bool_matrix(
    store: Arc<FilesystemStore>,
    name: &str,
    shape: &[usize; 2],
    values: &Array2<bool>,
    data_type: zarrs::array::DataType,
    attributes: Value,
) -> Result<()> {
    let values = values.as_standard_layout();
    let slice = values
        .as_slice()
        .with_context(|| format!("{name} matrix should be contiguous after layout conversion"))?;
    let array = create_zarr_array_with_fill_value(
        store,
        name,
        shape,
        shape,
        &["length", "gc_percent"],
        data_type,
        bool_fill_value(false),
        attributes,
    )?;
    if !slice.is_empty() {
        array
            .store_chunk(&[0, 0], slice)
            .with_context(|| format!("write Zarr chunk for array {name}"))?;
    }
    Ok(())
}

fn read_labels(path: &Path, array_name: &str, expected_field: &str) -> Result<Vec<String>> {
    let metadata: Value = serde_json::from_str(&std::fs::read_to_string(
        path.join(array_name).join("zarr.json"),
    )?)?;
    let attributes = metadata
        .get("attributes")
        .context("Zarr array metadata is missing attributes")?;
    ensure!(
        attributes
            .get("label_field")
            .and_then(Value::as_str)
            .is_some_and(|field| field == expected_field),
        "{array_name} metadata must declare label_field = {expected_field}"
    );
    let labels = attributes
        .get("labels")
        .and_then(Value::as_array)
        .with_context(|| format!("{array_name} metadata is missing labels"))?;
    labels
        .iter()
        .map(|label| {
            label
                .as_str()
                .map(str::to_string)
                .with_context(|| format!("{array_name} label should be a string"))
        })
        .collect()
}

fn ensure_contiguous_values(values: &[i32], axis_name: &str) -> Result<()> {
    for pair in values.windows(2) {
        ensure!(
            pair[1] == pair[0] + 1,
            "{axis_name} axis must contain contiguous integer values"
        );
    }
    Ok(())
}

fn ensure_index_axis(values: &[i32], axis_name: &str) -> Result<()> {
    for (idx, value) in values.iter().enumerate() {
        ensure!(
            *value == i32::try_from(idx).context("axis index exceeds i32")?,
            "{axis_name} axis must be contiguous 0-based indices"
        );
    }
    Ok(())
}

fn read_bool_attr(root: &Value, name: &str) -> Result<bool> {
    root.get(name)
        .and_then(Value::as_bool)
        .with_context(|| format!("Zarr root attribute {name} must be a bool"))
}

fn read_u8_attr(root: &Value, name: &str) -> Result<u8> {
    let value = root
        .get(name)
        .and_then(Value::as_u64)
        .with_context(|| format!("{name} in reference GC package must be an unsigned integer"))?;
    u8::try_from(value).with_context(|| format!("{name} in reference GC package must fit in u8"))
}

fn read_f64_attr(root: &Value, name: &str) -> Result<f64> {
    root.get(name)
        .and_then(Value::as_f64)
        .with_context(|| format!("Zarr root attribute {name} must be a float"))
}
