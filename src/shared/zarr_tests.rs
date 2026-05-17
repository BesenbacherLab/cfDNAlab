use super::*;
use serde_json::json;
use tempfile::TempDir;
use zarrs::array::{Array, data_type};

#[test]
fn write_single_chunk_zarr_array_round_trips_vector_values_and_metadata() {
    // Arrange: a tiny coordinate array should be one chunk with V3 dimension names.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("store.zarr");
    let store = create_zarr_store(&store_path, "test").expect("store should open");
    let values = vec![10i32, 20, 30];

    // Act
    write_single_chunk_zarr_array(
        store.clone(),
        "position",
        &[3],
        &["position"],
        &values,
        data_type::int32(),
        0i32,
        json!({"long_name": "position index"}),
    )
    .expect("array should write");

    // Assert
    let array = Array::open(store, "/position").expect("array should open");
    let round_trip: Vec<i32> = array
        .retrieve_array_subset(&array.subset_all())
        .expect("array should read");
    assert_eq!(round_trip, values);
    let metadata = read_json(&store_path.join("position/zarr.json"));
    assert_eq!(metadata["zarr_format"], 3);
    assert_eq!(metadata["dimension_names"], json!(["position"]));
    assert!(metadata["attributes"].get("_ARRAY_DIMENSIONS").is_none());
}

#[test]
fn write_single_chunk_zarr_array_allows_zero_length_axes() {
    // Arrange: sparse outputs can have zero non-zero entries or zero motif labels.
    let temp = TempDir::new().expect("temp dir should be created");
    let store_path = temp.path().join("store.zarr");
    let store = create_zarr_store(&store_path, "test").expect("store should open");
    let values: Vec<u64> = Vec::new();

    // Act
    write_single_chunk_zarr_array(
        store.clone(),
        "sparse/row",
        &[0],
        &["nnz"],
        &values,
        data_type::uint64(),
        0u64,
        json!({}),
    )
    .expect("empty array should write metadata");

    // Assert
    let array = Array::open(store, "/sparse/row").expect("empty array should open");
    let round_trip: Vec<u64> = array
        .retrieve_array_subset(&array.subset_all())
        .expect("empty array should read");
    assert!(round_trip.is_empty());
}

#[test]
fn validate_zarr_label_rejects_control_characters() {
    // Arrange
    let label = "bad\tlabel";

    // Act
    let error =
        validate_zarr_label(label, "group_name").expect_err("control characters should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("group_name contains a control character")
    );
}

#[test]
fn checked_integer_helpers_report_overflow_before_dtype_truncation() {
    // Arrange
    let i32_overflow = i32::MAX as usize + 1;

    // Act and assert
    assert_eq!(checked_index_axis(3, "row").unwrap(), vec![0, 1, 2]);
    assert!(checked_i32(i32_overflow, "row").is_err());
}

#[test]
fn shape_and_attribute_helpers_reject_schema_mistakes() {
    // Arrange / Act
    let valid_object = json_object(json!({"units": "count"})).expect("object attrs should parse");
    let non_object_error = json_object(json!(["not", "object"]))
        .expect_err("non-object attrs should be rejected");

    // Assert
    assert_eq!(element_count(&[2, 3, 4]), Some(24));
    assert_eq!(element_count(&[usize::MAX, 2]), None);
    assert_eq!(usize_slice_to_u64_vec(&[2, 3]).unwrap(), vec![2, 3]);
    assert_eq!(valid_object.get("units"), Some(&json!("count")));
    assert!(
        non_object_error
            .to_string()
            .contains("Zarr attributes must be a JSON object")
    );
}

fn read_json(path: &std::path::Path) -> serde_json::Value {
    let text = std::fs::read_to_string(path).expect("JSON should be readable");
    serde_json::from_str(&text).expect("JSON should parse")
}
