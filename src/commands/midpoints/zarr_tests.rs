use super::{
    checked_i32, checked_index_axis, checked_u32, counts_chunk_shape, create_store, element_count,
    encode_group_names, json_object, length_axis_coordinate_arrays, position_bin_coordinate_arrays,
    usize_slice_to_u64_vec, write_count_tensor_array, write_midpoint_profiles_zarr,
};
use crate::{
    commands::midpoints::{
        group_index::MidpointGroupSummary, postprocess::ProfileLayout,
        smoothing::MidpointSmoothing,
    },
    shared::length_axis::LengthAxis,
};
use fxhash::FxHashMap;
use ndarray::Array3;
use serde_json::Value;
use std::sync::Arc;
use tempfile::TempDir;
use zarrs::{array::Array, filesystem::FilesystemStore};

#[test]
fn midpoint_zarr_writes_counts_axes_and_group_metadata() {
    // Arrange:
    // Two groups, two length bins, and three positions exercise every axis while keeping the
    // single test chunk small enough to inspect directly.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("sample.midpoint_profiles.zarr");
    let counts = Array3::from_shape_vec(
        (2, 2, 3),
        vec![0.0, 1.0, 2.0, 10.0, 11.0, 12.0, 100.0, 101.0, 102.0, 110.0, 111.0, 112.0],
    )
    .expect("test count shape should be valid");
    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(1, "beta".to_string());
    group_idx_to_name.insert(0, "alpha".to_string());
    let mut eligible_interval_counts = FxHashMap::default();
    eligible_interval_counts.insert(0, 4usize);
    eligible_interval_counts.insert(1, 7usize);
    let length_axis =
        LengthAxis::new(vec![30, 40, 55]).expect("test length axis should be valid");
    let profile_layout = ProfileLayout::resolve(5, 2, MidpointSmoothing::None)
        .expect("test profile layout should be valid");

    // Act
    write_midpoint_profiles_zarr(
        &store_path,
        counts.view(),
        &group_idx_to_name,
        &eligible_interval_counts,
        &length_axis,
        profile_layout,
    )
    .expect("Zarr store should write");

    // Assert
    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(root_metadata["zarr_format"], 3);
    assert_eq!(root_metadata["node_type"], "group");
    assert_eq!(
        root_metadata["attributes"]["cfdnalab_schema"],
        "midpoint_profiles"
    );
    // Keep this literal so schema-version bumps require an explicit test update
    assert_eq!(
        root_metadata["attributes"]["cfdnalab_schema_version"],
        serde_json::json!(1)
    );

    let counts_metadata = read_json(&store_path.join("counts/zarr.json"));
    assert_eq!(counts_metadata["shape"], serde_json::json!([2, 2, 3]));
    assert_eq!(counts_metadata["data_type"], "float32");
    assert_eq!(
        counts_metadata["dimension_names"],
        serde_json::json!(["group", "length_bin", "position"])
    );
    assert!(
        counts_metadata["attributes"]
            .get("_ARRAY_DIMENSIONS")
            .is_none(),
        "Zarr V3 dimension names should not be duplicated as V2 _ARRAY_DIMENSIONS attrs"
    );
    assert_eq!(counts_metadata["codecs"][1]["name"], "zstd");
    assert_eq!(read_f32_chunk(&store_path.join("counts/c/0/0/0")), counts.as_slice().unwrap());
    assert_eq!(
        read_i32_chunk(&store_path.join("group/c/0")),
        vec![0, 1],
        "Zarr group coordinate must match group_idx and count row order"
    );

    assert_eq!(
        read_u32_chunk(&store_path.join("length_start_bp/c/0")),
        vec![30, 40]
    );
    assert_eq!(
        read_u32_chunk(&store_path.join("length_end_bp/c/0")),
        vec![40, 55]
    );
    assert_eq!(
        read_i32_chunk(&store_path.join("position_bin_start_bp/c/0")),
        vec![0, 2, 4]
    );
    assert_eq!(
        read_i32_chunk(&store_path.join("position_bin_end_bp/c/0")),
        vec![2, 4, 5]
    );
    assert_eq!(
        read_u32_chunk(&store_path.join("eligible_intervals/c/0")),
        vec![4, 7]
    );
    let group_name_metadata = read_json(&store_path.join("group_name/zarr.json"));
    assert_eq!(group_name_metadata["data_type"], "string");
    assert_eq!(group_name_metadata["codecs"][0]["name"], "vlen-utf8");
    assert_eq!(group_name_metadata["codecs"][1]["name"], "zstd");
    assert_eq!(
        read_string_array(&store_path, "/group_name"),
        vec!["alpha".to_string(), "beta".to_string()]
    );
    assert_eq!(
        read_group_names_from_byte_fallback(&store_path),
        vec!["alpha".to_string(), "beta".to_string()]
    );
    assert_eq!(
        read_i32_chunk(&store_path.join("group_name_byte/c/0")),
        vec![0, 1, 2, 3, 4]
    );
}

#[test]
fn midpoint_zarr_rejects_group_indices_that_do_not_match_count_rows() {
    // Arrange:
    // Count rows are addressed directly by group_idx. A sparse group_idx set would make it
    // ambiguous which metadata row describes count row 1, so the writer must fail instead.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("sample.midpoint_profiles.zarr");
    let counts = Array3::zeros((2, 1, 5));
    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(0, "alpha".to_string());
    group_idx_to_name.insert(2, "beta".to_string());
    let eligible_interval_counts = FxHashMap::default();
    let length_axis = LengthAxis::new(vec![30, 40]).expect("test length axis should be valid");
    let profile_layout = ProfileLayout::resolve(5, 1, MidpointSmoothing::None)
        .expect("test profile layout should be valid");

    // Act
    let error = write_midpoint_profiles_zarr(
        &store_path,
        counts.view(),
        &group_idx_to_name,
        &eligible_interval_counts,
        &length_axis,
        profile_layout,
    )
    .expect_err("non-contiguous group_idx values should be rejected");

    // Assert
    assert!(
        error
            .to_string()
            .contains("group indices must match count rows 0..1 but observed [0, 2]"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn count_tensor_array_round_trips_when_chunk_grid_has_partial_boundary_chunk() {
    // Arrange:
    // A [3, 2, 2] tensor with [2, 2, 2] chunks creates a partial final group chunk
    // This is the smallest useful regression case for boundary chunks that extend past the array shape
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("partial.midpoint_profiles.zarr");
    let store = create_store(&store_path).expect("Zarr store should be created");
    let counts = Array3::from_shape_vec(
        (3, 2, 2),
        vec![0.0, 1.0, 10.0, 11.0, 100.0, 101.0, 110.0, 111.0, 200.0, 201.0, 210.0, 211.0],
    )
    .expect("test count shape should be valid");

    // Act
    write_count_tensor_array(
        store.clone(),
        "counts",
        &[3, 2, 2],
        &[2, 2, 2],
        &["group", "length_bin", "position"],
        counts.view(),
        serde_json::json!({"units": "weighted_midpoint_count"}),
    )
    .expect("partial boundary chunk should write");

    // Assert
    let array = Array::open(store, "/counts").expect("counts array should open");
    let round_trip: Vec<f32> = array
        .retrieve_array_subset(&array.subset_all())
        .expect("counts array should read");
    assert_eq!(round_trip, counts.as_slice().unwrap());
}

#[test]
fn counts_chunk_shape_keeps_small_count_tensor_in_one_chunk() {
    // Arrange
    // 2 * 2 * 3 = 12 cells, well below the 4,000,000 cell target
    let shape = [2, 2, 3];

    // Act
    let chunk_shape = counts_chunk_shape(shape).expect("small shape should be valid");

    // Assert
    assert_eq!(chunk_shape, shape);
}

#[test]
fn counts_chunk_shape_chunks_large_count_tensor_across_group_axis() {
    // Arrange
    // The full position and length axes use 1,000,000 cells per group, so four groups fit in the
    // 4,000,000 cell target
    let shape = [10, 1_000, 1_000];

    // Act
    let chunk_shape = counts_chunk_shape(shape).expect("large shape should be valid");

    // Assert
    assert_eq!(chunk_shape, [4, 1_000, 1_000]);
}

#[test]
fn counts_chunk_shape_limits_wide_position_axis_first() {
    // Arrange
    // The position axis alone exceeds the target, so the position chunk is capped at 4,000,000 and
    // the remaining axes must shrink to one cell each
    let shape = [2, 3, 5_000_000];

    // Act
    let chunk_shape = counts_chunk_shape(shape).expect("wide shape should be valid");

    // Assert
    assert_eq!(chunk_shape, [1, 1, 4_000_000]);
}

#[test]
fn length_axis_coordinate_arrays_return_half_open_bin_edges() {
    // Arrange
    let length_axis =
        LengthAxis::new(vec![30, 40, 55]).expect("test length axis should be valid");

    // Act
    let (starts, ends) = length_axis_coordinate_arrays(&length_axis);

    // Assert
    assert_eq!(starts, vec![30, 40]);
    assert_eq!(ends, vec![40, 55]);
}

#[test]
fn position_bin_coordinate_arrays_record_short_final_bin() {
    // Arrange
    // Five output bases with bin size two gives bins [0, 2), [2, 4), and [4, 5)
    let profile_layout = ProfileLayout::resolve(5, 2, MidpointSmoothing::None)
        .expect("test profile layout should be valid");

    // Act
    let (starts, ends) =
        position_bin_coordinate_arrays(profile_layout).expect("position bins should fit i32");

    // Assert
    assert_eq!(starts, vec![0, 2, 4]);
    assert_eq!(ends, vec![2, 4, 5]);
}

#[test]
fn encode_group_names_builds_padded_utf8_fallback_matrix() {
    // Arrange
    // The longest name has width four, so the one-byte name is padded with three zero bytes
    let groups = vec![
        MidpointGroupSummary {
            group_idx: 0,
            group_name: "A",
            eligible_intervals: 1,
        },
        MidpointGroupSummary {
            group_idx: 1,
            group_name: "long",
            eligible_intervals: 2,
        },
    ];

    // Act
    let (encoded, lengths, width) =
        encode_group_names(&groups).expect("group names should encode");

    // Assert
    assert_eq!(width, 4);
    assert_eq!(lengths, vec![1, 4]);
    assert_eq!(encoded, vec![b'A', 0, 0, 0, b'l', b'o', b'n', b'g']);
}

#[test]
fn checked_integer_helpers_report_overflow_before_public_dtype_truncation() {
    // Arrange
    let i32_overflow_usize = i32::MAX as usize + 1;
    let i32_overflow_u64 = i32::MAX as u64 + 1;
    let u32_overflow_u64 = u32::MAX as u64 + 1;

    // Act and assert
    assert_eq!(checked_i32(42, "position").unwrap(), 42);
    assert!(
        checked_i32(i32_overflow_usize, "position")
            .unwrap_err()
            .to_string()
            .contains("position value")
    );
    assert_eq!(checked_i32(7u64, "group_idx").unwrap(), 7);
    assert!(
        checked_i32(i32_overflow_u64, "group_idx")
            .unwrap_err()
            .to_string()
            .contains("group_idx value")
    );
    assert_eq!(checked_u32(9usize, "eligible_intervals").unwrap(), 9);
    assert!(
        checked_u32(u32_overflow_u64, "eligible_intervals")
            .unwrap_err()
            .to_string()
            .contains("eligible_intervals value")
    );
}

#[test]
fn checked_index_axis_returns_zero_based_i32_values() {
    // Arrange
    let len = 4;

    // Act
    let axis = checked_index_axis(len, "position").expect("small axis should fit i32");

    // Assert
    assert_eq!(axis, vec![0, 1, 2, 3]);
}

#[test]
fn shape_and_json_helpers_validate_schema_metadata() {
    // Arrange
    let attributes = serde_json::json!({"long_name": "weighted midpoint count"});

    // Act
    let element_count_without_overflow = element_count(&[2, 3, 4]);
    let element_count_with_overflow = element_count(&[usize::MAX, 2]);
    let shape_as_u64 =
        usize_slice_to_u64_vec(&[2, 3, 4]).expect("small dimensions should fit u64");
    let parsed_attributes = json_object(attributes).expect("attributes should parse");
    let non_object_error = json_object(serde_json::json!([1, 2]))
        .expect_err("non-object attributes should be rejected");

    // Assert
    assert_eq!(element_count_without_overflow, Some(24));
    assert_eq!(element_count_with_overflow, None);
    assert_eq!(shape_as_u64, vec![2, 3, 4]);
    assert_eq!(
        parsed_attributes.get("long_name"),
        Some(&serde_json::json!("weighted midpoint count"))
    );
    assert!(
        parsed_attributes.get("_ARRAY_DIMENSIONS").is_none(),
        "Zarr V3 dimension names live in array metadata, not attrs"
    );
    assert!(
        non_object_error
            .to_string()
            .contains("Zarr attributes must be a JSON object")
    );
}

fn read_json(path: &std::path::Path) -> Value {
    let text = std::fs::read_to_string(path).expect("JSON should be readable");
    serde_json::from_str(&text).expect("JSON should parse")
}

fn read_chunk_bytes(path: &std::path::Path) -> Vec<u8> {
    let compressed = std::fs::read(path).expect("chunk should be readable");
    zstd::bulk::decompress(&compressed, 1 << 20).expect("chunk should decompress")
}

fn read_f32_chunk(path: &std::path::Path) -> Vec<f32> {
    read_chunk_bytes(path)
        .chunks_exact(4)
        .map(|bytes| f32::from_le_bytes(bytes.try_into().unwrap()))
        .collect()
}

fn read_u32_chunk(path: &std::path::Path) -> Vec<u32> {
    read_chunk_bytes(path)
        .chunks_exact(4)
        .map(|bytes| u32::from_le_bytes(bytes.try_into().unwrap()))
        .collect()
}

fn read_i32_chunk(path: &std::path::Path) -> Vec<i32> {
    read_chunk_bytes(path)
        .chunks_exact(4)
        .map(|bytes| i32::from_le_bytes(bytes.try_into().unwrap()))
        .collect()
}

fn read_string_array(store_path: &std::path::Path, array_path: &str) -> Vec<String> {
    let store = Arc::new(FilesystemStore::new(store_path).expect("Zarr store should open"));
    let array = Array::open(store, array_path).expect("Zarr string array should open");
    array
        .retrieve_array_subset(&array.subset_all())
        .expect("Zarr string array should read")
}

fn read_group_names_from_byte_fallback(store_path: &std::path::Path) -> Vec<String> {
    let name_bytes = read_chunk_bytes(&store_path.join("group_name_utf8/c/0/0"));
    let name_lengths = read_u32_chunk(&store_path.join("group_name_nbytes/c/0"));
    let width = name_bytes.len() / name_lengths.len();
    name_lengths
        .iter()
        .enumerate()
        .map(|(group_index, name_len)| {
            let start = group_index * width;
            let end = start + *name_len as usize;
            String::from_utf8(name_bytes[start..end].to_vec()).expect("valid UTF-8 group name")
        })
        .collect()
}
