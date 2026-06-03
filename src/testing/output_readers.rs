//! Readers for public cfDNAlab output artifacts.
//!
//! These helpers are for assertions in Rust tests. They intentionally read the
//! public file formats that cfDNAlab commands write, rather than reaching into
//! command internals. Use them when a test needs to compare numeric arrays or
//! compressed text without duplicating Zarr and zstd boilerplate.

use crate::shared::{constants::GC_CORRECTION_SCHEMA_VERSION, reference::ContigFootprintEntry};
use anyhow::{Context, Result, ensure};
use ndarray::{Array2, Array3};
use serde_json::Value;
use std::{
    fs::{File, OpenOptions},
    io::Read,
    path::Path,
    sync::Arc,
};
use zarrs::{
    array::{Array, ElementOwned},
    filesystem::FilesystemStore,
};
use zstd::stream::read::Decoder as ZstdDecoder;

/// Reference GC package written by `cfdna ref-gc-bias`.
///
/// Use this in tests that need to inspect the public reference-GC Zarr package
/// produced by the command. The struct mirrors the stable artifact fields that
/// downstream code can reasonably assert on: the numeric counts, the two
/// support masks, the GC-percent width correction array, and the package
/// metadata stored in Zarr attributes and coordinate arrays.
///
/// The arrays use the public package layout. Rows correspond to fragment
/// lengths in the inclusive metadata range `min_fragment_length..=max_fragment_length`.
/// Columns correspond to integer GC percent bins `0..=100`.
#[derive(Clone, Debug)]
pub struct ReferenceGCPackageOutput {
    /// Reference fragment mass by fragment length and integer GC percent.
    pub counts: Array2<f64>,
    /// Theoretical support mask based on reachable GC-percent bins.
    ///
    /// This is read from the public Zarr array `support_mask_unobservables`.
    /// A `true` value means the GC-percent bin is theoretically reachable for
    /// that effective fragment length after end trimming.
    pub unobservables_support_mask: Array2<bool>,
    /// Empirical support mask based on observed reference counts.
    ///
    /// This is read from the public Zarr array `support_mask_outliers`. A
    /// `true` value means the bin had enough empirical support to be treated as
    /// usable before downstream correction.
    pub outliers_support_mask: Array2<bool>,
    /// Number of raw GC-count states represented by each integer GC-percent bin.
    ///
    /// Widths are useful when tests need to distinguish raw counts from the
    /// width-corrected representation. Unreachable bins have width `0`.
    pub gc_percent_widths: Array2<u16>,
    /// Metadata read from root attributes and coordinate arrays.
    pub metadata: ReferenceGCPackageMetadata,
}

/// Metadata stored with a reference GC package.
///
/// These fields are the parts of the public `ref-gc-bias` artifact contract
/// that tests commonly need to check. They come from a mix of root attributes,
/// coordinate arrays, and the JSON-encoded reference contig footprint array.
///
/// The fragment length range is inclusive. `chromosomes` preserves the order
/// written by the command, which is also the order used by command-level
/// selection. `reference_contig_footprint` is the serialized reference identity
/// used by downstream GC correction compatibility checks.
#[derive(Clone, Debug, PartialEq)]
pub struct ReferenceGCPackageMetadata {
    /// Minimum fragment length represented by the package.
    pub min_fragment_length: usize,
    /// Maximum fragment length represented by the package.
    pub max_fragment_length: usize,
    /// Number of bases trimmed from each fragment end before GC counting.
    pub end_offset: u8,
    /// Chromosomes included in the run, in package order.
    pub chromosomes: Vec<String>,
    /// Reference contig footprint serialized by the command.
    pub reference_contig_footprint: Vec<ContigFootprintEntry>,
    /// Whether interpolation was skipped when building the package.
    pub skip_interpolation: bool,
    /// Gaussian smoothing sigma recorded by the package.
    pub smoothing_sigma: f64,
    /// Gaussian smoothing radius recorded by the package.
    pub smoothing_radius: u8,
    /// Whether smoothing was skipped when building the package.
    pub skip_smoothing: bool,
}

/// Read the `counts` array from a midpoint Zarr output.
///
/// `store_path` is the root directory of the Zarr store written by the command.
/// The helper opens the array at `/counts` and reads the entire array into
/// memory.
///
/// The returned array has shape `(group, length_bin, position)`, matching the
/// public midpoint Zarr schema and the in-memory layout used by cfDNAlab count
/// arrays. The helper checks that `/counts` is rank 3 before converting it to
/// `Array3<f32>`. It does not inspect store metadata beyond the array shape.
pub fn read_midpoint_zarr_counts<P: AsRef<Path>>(store_path: P) -> Result<Array3<f32>> {
    let array = open_zarr_array(store_path.as_ref(), "/counts")?;
    let shape = array.shape();
    ensure!(
        shape.len() == 3,
        "expected midpoint Zarr counts to be rank 3 but found rank {}",
        shape.len()
    );
    let values: Vec<f32> = array
        .retrieve_array_subset(&array.subset_all())
        .context("reading midpoint Zarr counts")?;
    let shape = (
        usize::try_from(shape[0]).context("group dimension exceeds usize")?,
        usize::try_from(shape[1]).context("length_bin dimension exceeds usize")?,
        usize::try_from(shape[2]).context("position dimension exceeds usize")?,
    );
    Array3::from_shape_vec(shape, values).context("building midpoint count array from Zarr values")
}

/// Read a one-dimensional signed-integer array from a midpoint Zarr output.
///
/// `store_path` is the root directory of the Zarr store written by the command.
/// `array_path` is the Zarr array path, for example `/positions`. The helper
/// verifies that the array is rank 1 and returns all values in stored order.
///
/// Use this for signed coordinate-like arrays in the public midpoint output.
/// The helper does not reinterpret coordinates or apply one-based conversion.
pub fn read_midpoint_zarr_i32_1d<P: AsRef<Path>>(
    store_path: P,
    array_path: &str,
) -> Result<Vec<i32>> {
    let array = open_zarr_array(store_path.as_ref(), array_path)?;
    ensure!(
        array.shape().len() == 1,
        "expected midpoint Zarr array {array_path} to be rank 1"
    );
    array
        .retrieve_array_subset(&array.subset_all())
        .with_context(|| format!("reading midpoint Zarr array {array_path}"))
}

/// Read a one-dimensional unsigned-integer array from a midpoint Zarr output.
///
/// `store_path` is the root directory of the Zarr store written by the command.
/// `array_path` is the Zarr array path, for example `/length_edges` when a
/// command writes unsigned length metadata. The helper verifies that the array
/// is rank 1 and returns all values in stored order.
///
/// Use this for unsigned metadata arrays. The helper does not validate that
/// values form valid edges or sorted coordinates.
pub fn read_midpoint_zarr_u32_1d<P: AsRef<Path>>(
    store_path: P,
    array_path: &str,
) -> Result<Vec<u32>> {
    let array = open_zarr_array(store_path.as_ref(), array_path)?;
    ensure!(
        array.shape().len() == 1,
        "expected midpoint Zarr array {array_path} to be rank 1"
    );
    array
        .retrieve_array_subset(&array.subset_all())
        .with_context(|| format!("reading midpoint Zarr array {array_path}"))
}

/// Read a reference GC Zarr package written by `cfdna ref-gc-bias`.
///
/// Use this when a test needs to assert the contents of the public
/// `<prefix>.ref_gc_package.zarr` or `ref_gc_package.zarr` artifact. The helper
/// opens the Zarr store, verifies that the root metadata declares the reference
/// GC package schema supported by this crate, and reads the full public arrays
/// into memory.
///
/// The helper reads the artifact as an external file format. It does not call
/// the crate-private production loader and it does not reuse command internals
/// that are allowed to change during implementation cleanup. That makes it
/// suitable for integration tests that are meant to stay in `tests/`.
///
/// Technical details
/// -----------------
/// - `counts`, `support_mask_unobservables`, `support_mask_outliers`, and
///   `gc_percent_widths` must be rank-2 arrays.
/// - The `length` coordinate must be present and non-empty. The first and last
///   values define the inclusive fragment length range in the returned metadata.
/// - Chromosome names are read from the `chromosome` array's `labels`
///   attribute, using `chromosome_name` as the expected label field.
/// - The reference footprint is decoded from `reference_contig_footprint_json`.
/// - The helper validates schema name and schema version, but it deliberately
///   does not duplicate every invariant checked by the production loader.
///
/// Parameters
/// ----------
/// - `package_path`:
///     Path to the root directory of the reference GC Zarr package.
///
/// Returns
/// -------
/// - `ReferenceGCPackageOutput`:
///     The package arrays and metadata needed for command-level artifact tests.
pub fn read_reference_gc_package<P: AsRef<Path>>(
    package_path: P,
) -> Result<ReferenceGCPackageOutput> {
    let package_path = package_path.as_ref();
    let root_attributes = read_zarr_root_attributes(package_path)?;
    ensure_reference_gc_schema(&root_attributes)?;

    let store = Arc::new(
        FilesystemStore::new(package_path)
            .with_context(|| format!("opening Zarr store {}", package_path.display()))?,
    );
    let counts = read_zarr_array2::<f64>(store.clone(), "/counts")?;
    let unobservables_support_mask =
        read_zarr_array2::<bool>(store.clone(), "/support_mask_unobservables")?;
    let outliers_support_mask = read_zarr_array2::<bool>(store.clone(), "/support_mask_outliers")?;
    let gc_percent_widths = read_zarr_array2::<u16>(store.clone(), "/gc_percent_widths")?;
    let lengths = read_zarr_array1::<i32>(store.clone(), "/length")?;
    let reference_contig_footprint_json =
        read_zarr_array1::<u8>(store, "/reference_contig_footprint_json")?;

    ensure!(
        !lengths.is_empty(),
        "reference GC package length axis must not be empty"
    );
    let min_fragment_length = usize::try_from(lengths[0]).context("length must be non-negative")?;
    let max_fragment_length =
        usize::try_from(*lengths.last().expect("length axis checked non-empty"))
            .context("length must be non-negative")?;
    let chromosomes = read_zarr_labels(package_path, "chromosome", "chromosome_name")?;
    let reference_contig_footprint: Vec<ContigFootprintEntry> =
        serde_json::from_slice(&reference_contig_footprint_json)
            .context("invalid reference_contig_footprint_json in reference GC package")?;

    Ok(ReferenceGCPackageOutput {
        counts,
        unobservables_support_mask,
        outliers_support_mask,
        gc_percent_widths,
        metadata: ReferenceGCPackageMetadata {
            min_fragment_length,
            max_fragment_length,
            end_offset: read_u8_attr(&root_attributes, "end_offset")?,
            chromosomes,
            reference_contig_footprint,
            skip_interpolation: read_bool_attr(&root_attributes, "skip_interpolation")?,
            smoothing_sigma: read_f64_attr(&root_attributes, "smoothing_sigma")?,
            smoothing_radius: read_u8_attr(&root_attributes, "smoothing_radius")?,
            skip_smoothing: read_bool_attr(&root_attributes, "skip_smoothing")?,
        },
    })
}

/// Read a zstd-compressed UTF-8 text file into a string.
///
/// Use this for command outputs such as compressed TSV files where the test
/// wants to assert exact text or parse rows manually. The helper fully
/// decompresses the file and decodes it as UTF-8. It does not normalize line
/// endings or trim trailing newlines.
pub fn read_zst_to_string(path: &Path) -> Result<String> {
    let reader = File::open(path)?;
    let mut decoder = ZstdDecoder::new(reader)?;
    let mut buffer = String::new();
    decoder.read_to_string(&mut buffer)?;
    Ok(buffer)
}

/// Read a zstd-compressed length-count TSV as text.
///
/// This is a named wrapper around `read_zst_to_string` for tests that document
/// they are reading the public `cfdna lengths` output format. The returned text
/// is unchanged after zstd decompression.
pub fn read_length_counts_text<P: AsRef<Path>>(path: P) -> Result<String> {
    read_zst_to_string(path.as_ref())
}

/// Read a length-count TSV and return only numeric `count_*` columns.
///
/// Metadata columns are ignored. The returned array has one row per TSV data
/// row and one column per `count_*` header column. Use this when the test only
/// cares about the length-count matrix and not the surrounding window or group
/// metadata.
///
/// `path` must point to the zstd-compressed text output. The helper finds the
/// first header column whose name starts with `count_` and treats that column
/// and all following columns as numeric count columns. Earlier columns are
/// skipped. Each data row must have the same number of tab-separated fields as
/// the header, and each count field must parse as `f64`.
///
/// The returned array shape is `(data_row_count, count_column_count)`. Header
/// names and metadata fields are not returned.
pub fn read_length_counts_tsv<P: AsRef<Path>>(path: P) -> Result<Array2<f64>> {
    let text = read_length_counts_text(path)?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .context("length counts TSV must have a header")?;
    let headers: Vec<&str> = header.split('\t').collect();
    let first_count_column = headers
        .iter()
        .position(|column| column.starts_with("count_"))
        .context("length counts TSV must contain count columns")?;
    let count_column_count = headers.len() - first_count_column;
    ensure!(
        count_column_count > 0,
        "length counts TSV must contain at least one count column"
    );

    let mut values = Vec::new();
    let mut row_count = 0;
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        ensure!(
            fields.len() == headers.len(),
            "length counts row must match the header column count"
        );
        for value in &fields[first_count_column..] {
            values.push(value.parse::<f64>()?);
        }
        row_count += 1;
    }

    Ok(Array2::from_shape_vec(
        (row_count, count_column_count),
        values,
    )?)
}

/// Create an empty file.
///
/// This is useful for command tests that need an existing placeholder path.
/// Existing files are opened without truncation, so existing contents are left
/// unchanged. Missing parent directories are not created.
pub fn touch_file<P: AsRef<Path>>(path: P) -> Result<()> {
    OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(false)
        .open(path)?;
    Ok(())
}

fn open_zarr_array(store_path: &Path, array_path: &str) -> Result<Array<FilesystemStore>> {
    let store = Arc::new(
        FilesystemStore::new(store_path)
            .with_context(|| format!("opening Zarr store {}", store_path.display()))?,
    );
    Array::open(store, array_path).with_context(|| {
        format!(
            "opening Zarr array {array_path} in {}",
            store_path.display()
        )
    })
}

fn read_zarr_root_attributes(package_path: &Path) -> Result<Value> {
    let metadata: Value =
        serde_json::from_str(&std::fs::read_to_string(package_path.join("zarr.json"))?)?;
    metadata
        .get("attributes")
        .cloned()
        .context("Zarr root metadata is missing attributes")
}

fn ensure_reference_gc_schema(root_attributes: &Value) -> Result<()> {
    let schema = root_attributes
        .get("cfdnalab_schema")
        .and_then(Value::as_str);
    ensure!(
        schema == Some("reference_gc_package"),
        "Reference GC package schema mismatch: file={schema:?}, expected=reference_gc_package"
    );
    let version = root_attributes
        .get("cfdnalab_schema_version")
        .and_then(Value::as_u64)
        .context("Reference GC package is missing cfdnalab_schema_version")?;
    ensure!(
        version == u64::from(GC_CORRECTION_SCHEMA_VERSION),
        "Reference GC package schema version mismatch: file={}, expected={}",
        version,
        GC_CORRECTION_SCHEMA_VERSION
    );
    Ok(())
}

fn read_zarr_array1<T>(store: Arc<FilesystemStore>, array_path: &str) -> Result<Vec<T>>
where
    T: ElementOwned,
{
    let array = Array::open(store, array_path)?;
    array
        .retrieve_array_subset(&array.subset_all())
        .with_context(|| format!("reading Zarr array {array_path}"))
}

fn read_zarr_array2<T>(store: Arc<FilesystemStore>, array_path: &str) -> Result<Array2<T>>
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
    let values: Vec<T> = array
        .retrieve_array_subset(&array.subset_all())
        .with_context(|| format!("reading Zarr array {array_path}"))?;
    let rows = usize::try_from(shape[0]).context("array row count exceeds usize")?;
    let cols = usize::try_from(shape[1]).context("array column count exceeds usize")?;
    Ok(Array2::from_shape_vec((rows, cols), values)?)
}

fn read_zarr_labels(
    package_path: &Path,
    array_name: &str,
    expected_field: &str,
) -> Result<Vec<String>> {
    let metadata: Value = serde_json::from_str(&std::fs::read_to_string(
        package_path.join(array_name).join("zarr.json"),
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

fn read_bool_attr(root_attributes: &Value, name: &str) -> Result<bool> {
    root_attributes
        .get(name)
        .and_then(Value::as_bool)
        .with_context(|| format!("Zarr root attribute {name} must be a bool"))
}

fn read_u8_attr(root_attributes: &Value, name: &str) -> Result<u8> {
    let value = root_attributes
        .get(name)
        .and_then(Value::as_u64)
        .with_context(|| format!("{name} in reference GC package must be an unsigned integer"))?;
    u8::try_from(value).with_context(|| format!("{name} in reference GC package must fit in u8"))
}

fn read_f64_attr(root_attributes: &Value, name: &str) -> Result<f64> {
    root_attributes
        .get(name)
        .and_then(Value::as_f64)
        .with_context(|| format!("Zarr root attribute {name} must be a float"))
}
