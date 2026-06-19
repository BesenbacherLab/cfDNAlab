#![cfg(feature = "cmd_midpoints")]
//! Public API tests for Rust output loaders for `cfdna midpoints`.
//!
//! These tests build tiny Zarr V3 midpoint profile stores and assert that the
//! loader reads metadata eagerly while count array reads happen through the
//! public selection methods.

use cfdnalab::{
    interval::Interval,
    output_loaders::{MidpointGroupRow, load_midpoints_output},
};
use serde_json::{Map, Value, json};
use std::{fs, path::Path, sync::Arc};
use tempfile::TempDir;
use zarrs::{
    array::{ArrayBuilder, DataType, Element, builder::ArrayBuilderFillValue, data_type},
    filesystem::FilesystemStore,
    group::GroupBuilder,
};

/// Verify midpoint axes and count selections from a valid Zarr store.
#[test]
fn load_midpoints_output_reads_axes_and_count_selections() -> anyhow::Result<()> {
    // Arrange:
    // The count array is row-major over group, length bin, then position. The
    // values are chosen so each selected coordinate has a distinct expected
    // count that can be derived directly from its source position.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.midpoint_profiles.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_f32_array(
        &store,
        "counts",
        &[2, 2, 4],
        &["group", "length_bin", "position"],
        &[
            0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
        ],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "group",
        &[2],
        &["group"],
        &[0, 1],
        json!({
            "label_field": "group_name",
            "labels": ["alpha", "beta"],
        }),
    )?;
    write_i32_array(
        &store,
        "eligible_intervals",
        &[2],
        &["group"],
        &[3, 0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "length_bin",
        &[2],
        &["length_bin"],
        &[0, 1],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "length_start_bp",
        &[2],
        &["length_bin"],
        &[30, 100],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "length_end_bp",
        &[2],
        &["length_bin"],
        &[100, 151],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "position",
        &[4],
        &["position"],
        &[0, 1, 2, 3],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "position_bin_start_bp",
        &[4],
        &["position"],
        &[0, 5, 10, 15],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "position_bin_end_bp",
        &[4],
        &["position"],
        &[5, 10, 15, 18],
        json!({}),
    )?;

    // Act
    let loaded = load_midpoints_output(&path)?;

    // Assert
    assert_eq!(loaded.counts_shape(), (2, 2, 4));
    let output_metadata = loaded.output_metadata();
    assert_eq!(output_metadata.group_count, 2);
    assert_eq!(output_metadata.length_bin_count, 2);
    assert_eq!(output_metadata.min_fragment_length, 30);
    assert_eq!(output_metadata.max_fragment_length_exclusive, 151);
    assert_eq!(output_metadata.position_bin_count, 4);
    assert_eq!(output_metadata.min_position, 0);
    assert_eq!(output_metadata.max_position_exclusive, 18);
    assert_eq!(
        output_metadata.to_string(),
        "group_count=2, length_bin_count=2, fragment_length_range=[30, 151) bp, position_bin_count=4, position_range=[0, 18) bp"
    );
    assert_eq!(loaded.group_count(), 2);
    assert_eq!(loaded.length_bin_count(), 2);
    assert_eq!(loaded.position_bin_count(), 4);
    assert_eq!(
        loaded.group_metadata(),
        &[
            MidpointGroupRow {
                index: 0,
                name: "alpha".to_string(),
                eligible_intervals: 3,
            },
            MidpointGroupRow {
                index: 1,
                name: "beta".to_string(),
                eligible_intervals: 0,
            },
        ]
    );
    assert_eq!(loaded.group_index("beta")?, 1);
    assert!(loaded.has_group("alpha"));
    assert_eq!(loaded.length_bin_for_length(120), Some(1));
    assert_eq!(loaded.position_bin_for_position(16), Some(3));
    assert_eq!(loaded.length_bins()[0].as_tuple(), (30, 100, 0));
    assert_eq!(loaded.position_bins()[3].as_tuple(), (15, 18, 3));

    let full_counts = loaded.read_all_counts()?;
    assert_eq!(full_counts.shape(), (2, 2, 4));
    assert_eq!(
        full_counts.counts_row_major(),
        &[
            0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
        ]
    );

    let selected = loaded
        .select()
        .groups_by_name(&["beta", "alpha"])
        .length_range(Interval::new(90, 130)?)
        .positions(&[3, 1])
        .read()?;
    assert_eq!(selected.group_indices(), &[1, 0]);
    assert_eq!(selected.groups()[0].name, "beta");
    assert_eq!(selected.length_bins()[0].as_tuple(), (30, 100, 0));
    assert_eq!(selected.length_bins()[1].as_tuple(), (100, 151, 1));
    assert_eq!(selected.position_bins()[0].as_tuple(), (15, 18, 3));
    assert_eq!(selected.position_bins()[1].as_tuple(), (5, 10, 1));
    assert_eq!(selected.shape(), (2, 2, 2));
    assert_eq!(selected.group_count(), 2);
    assert_eq!(selected.length_bin_count(), 2);
    assert_eq!(selected.position_bin_count(), 2);
    assert_eq!(
        selected.counts_row_major(),
        &[11.0, 9.0, 15.0, 13.0, 3.0, 1.0, 7.0, 5.0]
    );
    assert_eq!(selected.count(1, 0, 1), Some(1.0));
    assert_eq!(selected.profile(0, 1), Some(&[15.0, 13.0][..]));
    let mut profile_summaries = Vec::new();
    for (selected_group_index, group) in selected.groups().iter().enumerate() {
        for (selected_length_index, length_bin) in selected.length_bins().iter().enumerate() {
            let profile = selected
                .profile(selected_group_index, selected_length_index)
                .expect("selected indices should be in bounds");
            profile_summaries.push((
                group.name.as_str(),
                length_bin.as_tuple(),
                profile.iter().copied().sum::<f32>(),
            ));
        }
    }
    assert_eq!(
        profile_summaries,
        vec![
            ("beta", (30, 100, 0), 20.0),
            ("beta", (100, 151, 1), 28.0),
            ("alpha", (30, 100, 0), 4.0),
            ("alpha", (100, 151, 1), 12.0),
        ]
    );

    let selected_contiguous_group_length = loaded
        .select()
        .groups(&[0, 1])
        .length_range(Interval::new(90, 130)?)
        .positions(&[3, 1])
        .read()?;
    assert_eq!(selected_contiguous_group_length.group_indices(), &[0, 1]);
    assert_eq!(
        selected_contiguous_group_length.counts_row_major(),
        &[3.0, 1.0, 7.0, 5.0, 11.0, 9.0, 15.0, 13.0]
    );

    let duplicate_group_error = loaded
        .select()
        .groups_by_name(&["alpha", "alpha"])
        .read()
        .expect_err("duplicate group-name selectors should fail");
    let bad_group_index_error = loaded
        .select()
        .groups(&[2])
        .read()
        .expect_err("group index should be validated");
    let duplicate_length_error = loaded
        .select()
        .length_bins(&[0, 0])
        .read()
        .expect_err("duplicate length-bin selectors should fail");
    let duplicate_position_error = loaded
        .select()
        .positions(&[1, 1])
        .read()
        .expect_err("duplicate position selectors should fail");
    let empty_position_error = loaded
        .select()
        .positions(&[])
        .read()
        .expect_err("empty position selectors should fail");
    let missing_length_range_error = loaded
        .select()
        .length_range(Interval::new(200, 250)?)
        .read()
        .expect_err("non-overlapping length range should fail");
    let missing_position_range_error = loaded
        .select()
        .position_range(Interval::new(20, 30)?)
        .read()
        .expect_err("non-overlapping position range should fail");
    let conflicting_group_selector_error = loaded
        .select()
        .groups(&[0])
        .groups_by_name(&["alpha"])
        .read()
        .expect_err("conflicting group selectors should fail");
    let conflicting_length_selector_error = loaded
        .select()
        .length_bins(&[0])
        .length_range(Interval::new(30, 100)?)
        .read()
        .expect_err("conflicting length selectors should fail");
    let conflicting_position_selector_error = loaded
        .select()
        .positions(&[0])
        .position_range(Interval::new(0, 5)?)
        .read()
        .expect_err("conflicting position selectors should fail");
    assert!(
        duplicate_group_error
            .to_string()
            .contains("duplicate value 'alpha'")
    );
    assert!(
        bad_group_index_error
            .to_string()
            .contains("group index 2 is outside")
    );
    assert!(
        duplicate_length_error
            .to_string()
            .contains("length bin indices contain duplicate value 0")
    );
    assert!(
        duplicate_position_error
            .to_string()
            .contains("position indices contain duplicate value 1")
    );
    assert!(
        empty_position_error
            .to_string()
            .contains("cannot select zero midpoint positions")
    );
    assert!(
        missing_length_range_error
            .to_string()
            .contains("does not overlap any midpoint length bins")
    );
    assert!(
        missing_position_range_error
            .to_string()
            .contains("does not overlap any midpoint position bins")
    );
    assert!(
        conflicting_group_selector_error
            .to_string()
            .contains("cannot combine groups() and groups_by_name() on the group axis")
    );
    assert!(
        conflicting_length_selector_error.to_string().contains(
            "cannot combine length_bins() and length_range() on the fragment length axis"
        )
    );
    assert!(
        conflicting_position_selector_error
            .to_string()
            .contains("cannot combine positions() and position_range() on the position axis")
    );
    Ok(())
}

/// Verify unsupported midpoint schema versions are rejected.
#[test]
fn load_midpoints_output_rejects_wrong_schema_version() -> anyhow::Result<()> {
    // Arrange:
    // The loader targets the current midpoint profile schema and rejects other
    // advertised versions before trying to read axis arrays.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.midpoint_profiles.zarr");
    create_store(
        &path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 2,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;

    // Act
    let error = load_midpoints_output(&path).expect_err("schema version mismatch should fail");

    // Assert
    assert!(error.to_string().contains("schema version mismatch"));
    Ok(())
}

/// Verify midpoint count arrays must be rank three.
#[test]
fn load_midpoints_output_rejects_bad_count_rank() -> anyhow::Result<()> {
    // Arrange:
    // Midpoint stores must expose a rank-3 count array. A rank-2 array cannot
    // be interpreted as `(group, length_bin, position)`.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.midpoint_profiles.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_f32_array(
        &store,
        "counts",
        &[1, 1],
        &["group", "length_bin"],
        &[1.0],
        json!({}),
    )?;

    // Act
    let error = load_midpoints_output(&path).expect_err("rank-2 counts should fail");

    // Assert
    assert!(error.to_string().contains("rank 3, found rank 2"));
    Ok(())
}

/// Verify midpoint coordinate axes must be zero-based.
#[test]
fn load_midpoints_output_rejects_non_zero_based_axes() -> anyhow::Result<()> {
    // Arrange:
    // Axis coordinate arrays are public indices. They must be zero-based and
    // contiguous so selectors and returned metadata agree.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.midpoint_profiles.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_f32_array(
        &store,
        "counts",
        &[1, 1, 1],
        &["group", "length_bin", "position"],
        &[1.0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "group",
        &[1],
        &["group"],
        &[0],
        json!({
            "label_field": "group_name",
            "labels": ["alpha"],
        }),
    )?;
    write_i32_array(
        &store,
        "eligible_intervals",
        &[1],
        &["group"],
        &[1],
        json!({}),
    )?;
    write_i32_array(&store, "length_bin", &[1], &["length_bin"], &[1], json!({}))?;
    write_i32_array(
        &store,
        "length_start_bp",
        &[1],
        &["length_bin"],
        &[30],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "length_end_bp",
        &[1],
        &["length_bin"],
        &[40],
        json!({}),
    )?;
    write_i32_array(&store, "position", &[1], &["position"], &[0], json!({}))?;
    write_i32_array(
        &store,
        "position_bin_start_bp",
        &[1],
        &["position"],
        &[0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "position_bin_end_bp",
        &[1],
        &["position"],
        &[1],
        json!({}),
    )?;

    // Act
    let error = load_midpoints_output(&path).expect_err("non-zero-based axis should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("length_bin must be a zero-based coordinate axis")
    );
    Ok(())
}

/// Verify midpoint interval axes must match the writer's sorted axis contract.
#[test]
fn load_midpoints_output_rejects_invalid_axis_intervals() -> anyhow::Result<()> {
    // Arrange:
    // Length bins and position bins are indexed interval axes. Gaps or
    // out-of-range length bins would make range selection disagree with the
    // command's generated axis.
    let temp = TempDir::new()?;
    let length_gap_path = temp.path().join("length_gap.midpoint_profiles.zarr");
    let length_gap_store = create_store(
        &length_gap_path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_minimal_midpoint_store(
        &length_gap_store,
        &[1, 2, 1],
        &[1.0, 2.0],
        &[30, 50],
        &[40, 60],
        &[0],
        &[1],
    )?;

    let below_min_path = temp.path().join("below_min.midpoint_profiles.zarr");
    let below_min_store = create_store(
        &below_min_path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_minimal_midpoint_store(
        &below_min_store,
        &[1, 1, 1],
        &[1.0],
        &[9],
        &[10],
        &[0],
        &[1],
    )?;

    let position_gap_path = temp.path().join("position_gap.midpoint_profiles.zarr");
    let position_gap_store = create_store(
        &position_gap_path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_minimal_midpoint_store(
        &position_gap_store,
        &[1, 1, 2],
        &[1.0, 2.0],
        &[30],
        &[40],
        &[0, 3],
        &[2, 4],
    )?;

    // Act
    let length_gap_error =
        load_midpoints_output(&length_gap_path).expect_err("gapped length axis should fail");
    let below_min_error =
        load_midpoints_output(&below_min_path).expect_err("length bin below minimum should fail");
    let position_gap_error =
        load_midpoints_output(&position_gap_path).expect_err("gapped position axis should fail");

    // Assert
    assert!(
        length_gap_error
            .to_string()
            .contains("length_bin intervals must be contiguous and sorted")
    );
    assert!(
        below_min_error
            .to_string()
            .contains("starts below minimum supported fragment length")
    );
    assert!(
        position_gap_error
            .to_string()
            .contains("position intervals must be contiguous and sorted")
    );
    Ok(())
}

/// Verify midpoint count reads reject non-finite values but allow negatives.
#[test]
fn midpoint_counts_reject_nonfinite_and_allow_negative_values() -> anyhow::Result<()> {
    // Arrange:
    // Smoothed midpoint profiles can contain finite negative values because the
    // smoothing kernel has negative coefficients. Non-finite values are still
    // corrupt profile values and should fail when counts are read.
    let temp = TempDir::new()?;
    let negative_path = temp.path().join("negative.midpoint_profiles.zarr");
    let negative_store = create_store(
        &negative_path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_minimal_midpoint_store(
        &negative_store,
        &[1, 1, 1],
        &[-1.5],
        &[30],
        &[40],
        &[0],
        &[1],
    )?;

    let non_finite_path = temp.path().join("non_finite.midpoint_profiles.zarr");
    let non_finite_store = create_store(
        &non_finite_path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_minimal_midpoint_store(
        &non_finite_store,
        &[1, 1, 1],
        &[f32::NAN],
        &[30],
        &[40],
        &[0],
        &[1],
    )?;

    // Act
    let negative_loaded = load_midpoints_output(&negative_path)?;
    let negative_counts = negative_loaded.read_all_counts()?;
    let non_finite_loaded = load_midpoints_output(&non_finite_path)?;
    let non_finite_error = non_finite_loaded
        .read_all_counts()
        .expect_err("non-finite midpoint count should fail");

    // Assert
    assert_eq!(negative_counts.counts_row_major(), &[-1.5]);
    assert!(
        non_finite_error
            .to_string()
            .contains("midpoint counts contain non-finite value")
    );
    Ok(())
}

/// Verify midpoint loaders reject malformed public labels from Zarr metadata.
#[test]
fn load_midpoints_output_rejects_control_character_group_labels() -> anyhow::Result<()> {
    // Arrange:
    // Group names become public selector strings. Control characters would not
    // stay stable when users move labels into tables or command examples.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.midpoint_profiles.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_f32_array(
        &store,
        "counts",
        &[1, 1, 1],
        &["group", "length_bin", "position"],
        &[1.0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "group",
        &[1],
        &["group"],
        &[0],
        json!({
            "label_field": "group_name",
            "labels": ["bad\nlabel"],
        }),
    )?;
    write_i32_array(
        &store,
        "eligible_intervals",
        &[1],
        &["group"],
        &[1],
        json!({}),
    )?;
    write_i32_array(&store, "length_bin", &[1], &["length_bin"], &[0], json!({}))?;
    write_i32_array(
        &store,
        "length_start_bp",
        &[1],
        &["length_bin"],
        &[30],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "length_end_bp",
        &[1],
        &["length_bin"],
        &[40],
        json!({}),
    )?;
    write_i32_array(&store, "position", &[1], &["position"], &[0], json!({}))?;
    write_i32_array(
        &store,
        "position_bin_start_bp",
        &[1],
        &["position"],
        &[0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "position_bin_end_bp",
        &[1],
        &["position"],
        &[1],
        json!({}),
    )?;

    // Act
    let error = load_midpoints_output(&path).expect_err("control-character label should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("Zarr label group_name contains a control character")
    );
    Ok(())
}

/// Verify duplicate midpoint group names are rejected.
#[test]
fn load_midpoints_output_rejects_duplicate_group_names() -> anyhow::Result<()> {
    // Arrange:
    // Group names support public lookup and name-based selectors. Duplicate
    // names would make those operations ambiguous.
    let temp = TempDir::new()?;
    let path = temp.path().join("sample.midpoint_profiles.zarr");
    let store = create_store(
        &path,
        json!({
            "cfdnalab_schema": "midpoint_profiles",
            "cfdnalab_schema_version": 1,
            "primary_array": "counts",
            "count_units": "weighted_midpoint_count",
        }),
    )?;
    write_f32_array(
        &store,
        "counts",
        &[2, 1, 1],
        &["group", "length_bin", "position"],
        &[1.0, 2.0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "group",
        &[2],
        &["group"],
        &[0, 1],
        json!({
            "label_field": "group_name",
            "labels": ["alpha", "alpha"],
        }),
    )?;
    write_i32_array(
        &store,
        "eligible_intervals",
        &[2],
        &["group"],
        &[1, 1],
        json!({}),
    )?;
    write_i32_array(&store, "length_bin", &[1], &["length_bin"], &[0], json!({}))?;
    write_i32_array(
        &store,
        "length_start_bp",
        &[1],
        &["length_bin"],
        &[30],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "length_end_bp",
        &[1],
        &["length_bin"],
        &[40],
        json!({}),
    )?;
    write_i32_array(&store, "position", &[1], &["position"], &[0], json!({}))?;
    write_i32_array(
        &store,
        "position_bin_start_bp",
        &[1],
        &["position"],
        &[0],
        json!({}),
    )?;
    write_i32_array(
        &store,
        "position_bin_end_bp",
        &[1],
        &["position"],
        &[1],
        json!({}),
    )?;

    // Act
    let error = load_midpoints_output(&path).expect_err("duplicate group names should fail");

    // Assert
    assert!(
        error
            .to_string()
            .contains("midpoint group name is not unique: alpha")
    );
    Ok(())
}

/// Write a minimal valid midpoint profile store around caller-provided axes.
fn write_minimal_midpoint_store(
    store: &Arc<FilesystemStore>,
    counts_shape: &[usize; 3],
    count_values: &[f32],
    length_starts: &[i32],
    length_ends: &[i32],
    position_starts: &[i32],
    position_ends: &[i32],
) -> anyhow::Result<()> {
    assert_eq!(counts_shape[0], 1);
    assert_eq!(counts_shape[1], length_starts.len());
    assert_eq!(counts_shape[1], length_ends.len());
    assert_eq!(counts_shape[2], position_starts.len());
    assert_eq!(counts_shape[2], position_ends.len());
    let length_indices = (0..(length_starts.len() as i32)).collect::<Vec<_>>();
    let position_indices = (0..(position_starts.len() as i32)).collect::<Vec<_>>();

    write_f32_array(
        store,
        "counts",
        counts_shape,
        &["group", "length_bin", "position"],
        count_values,
        json!({}),
    )?;
    write_i32_array(
        store,
        "group",
        &[1],
        &["group"],
        &[0],
        json!({
            "label_field": "group_name",
            "labels": ["alpha"],
        }),
    )?;
    write_i32_array(
        store,
        "eligible_intervals",
        &[1],
        &["group"],
        &[1],
        json!({}),
    )?;
    write_i32_array(
        store,
        "length_bin",
        &[length_indices.len()],
        &["length_bin"],
        &length_indices,
        json!({}),
    )?;
    write_i32_array(
        store,
        "length_start_bp",
        &[length_starts.len()],
        &["length_bin"],
        length_starts,
        json!({}),
    )?;
    write_i32_array(
        store,
        "length_end_bp",
        &[length_ends.len()],
        &["length_bin"],
        length_ends,
        json!({}),
    )?;
    write_i32_array(
        store,
        "position",
        &[position_indices.len()],
        &["position"],
        &position_indices,
        json!({}),
    )?;
    write_i32_array(
        store,
        "position_bin_start_bp",
        &[position_starts.len()],
        &["position"],
        position_starts,
        json!({}),
    )?;
    write_i32_array(
        store,
        "position_bin_end_bp",
        &[position_ends.len()],
        &["position"],
        position_ends,
        json!({}),
    )?;
    Ok(())
}

/// Create a temporary Zarr store with root attributes.
fn create_store(path: &Path, attributes: Value) -> anyhow::Result<Arc<FilesystemStore>> {
    fs::create_dir_all(path)?;
    let store = Arc::new(FilesystemStore::new(path)?);
    write_group(&store, "/", attributes)?;
    Ok(store)
}

/// Write one Zarr group metadata object.
fn write_group(
    store: &Arc<FilesystemStore>,
    group_path: &str,
    attributes: Value,
) -> anyhow::Result<()> {
    let group = GroupBuilder::new()
        .attributes(json_object(attributes)?)
        .build(store.clone(), group_path)?;
    group.store_metadata()?;
    Ok(())
}

/// Write an `i32` Zarr array fixture.
fn write_i32_array(
    store: &Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    dimension_names: &[&str],
    values: &[i32],
    attributes: Value,
) -> anyhow::Result<()> {
    write_array(
        store,
        name,
        shape,
        dimension_names,
        values,
        data_type::int32(),
        -1,
        attributes,
    )
}

/// Write an `f32` Zarr array fixture.
fn write_f32_array(
    store: &Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    dimension_names: &[&str],
    values: &[f32],
    attributes: Value,
) -> anyhow::Result<()> {
    write_array(
        store,
        name,
        shape,
        dimension_names,
        values,
        data_type::float32(),
        -1.0,
        attributes,
    )
}

/// Write a typed Zarr array fixture with one chunk.
#[allow(clippy::too_many_arguments)]
fn write_array<T>(
    store: &Arc<FilesystemStore>,
    name: &str,
    shape: &[usize],
    dimension_names: &[&str],
    values: &[T],
    data_type: DataType,
    fill_value: T,
    attributes: Value,
) -> anyhow::Result<()>
where
    T: Element + Copy,
    ArrayBuilderFillValue: From<T>,
{
    let expected_len = shape
        .iter()
        .try_fold(1usize, |size, dimension| size.checked_mul(*dimension))
        .expect("test Zarr shape should not overflow");
    assert_eq!(values.len(), expected_len);
    let chunk_shape = shape
        .iter()
        .map(|dimension| (*dimension).max(1))
        .collect::<Vec<_>>();
    let array = ArrayBuilder::new(
        shape
            .iter()
            .map(|dimension| *dimension as u64)
            .collect::<Vec<_>>(),
        chunk_shape
            .iter()
            .map(|dimension| *dimension as u64)
            .collect::<Vec<_>>(),
        data_type,
        fill_value,
    )
    .dimension_names(Some(dimension_names.iter().copied()))
    .attributes(json_object(attributes)?)
    .build(store.clone(), format!("/{name}").as_str())?;
    array.store_metadata()?;
    if !values.is_empty() {
        array.store_chunk(&vec![0; shape.len()], values)?;
    }
    Ok(())
}

/// Convert a JSON value into an object for Zarr attributes.
fn json_object(value: Value) -> anyhow::Result<Map<String, Value>> {
    match value {
        Value::Object(map) => Ok(map),
        other => anyhow::bail!("expected object attributes, got {other}"),
    }
}
