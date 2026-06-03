//! Readers for public cfDNAlab output artifacts.
//!
//! These helpers are for assertions in Rust tests. They intentionally read the
//! public file formats that cfDNAlab commands write, rather than reaching into
//! command internals. Use them when a test needs to compare numeric arrays or
//! compressed text without duplicating Zarr and zstd boilerplate.

use anyhow::{Context, Result, ensure};
use ndarray::{Array2, Array3};
use std::{
    fs::{File, OpenOptions},
    io::Read,
    path::Path,
    sync::Arc,
};
use zarrs::{array::Array, filesystem::FilesystemStore};
use zstd::stream::read::Decoder as ZstdDecoder;

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
