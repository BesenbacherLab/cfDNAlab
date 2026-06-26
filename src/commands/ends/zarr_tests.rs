use super::*;
use crate::shared::kmers::motifs_file::SelectedMotifColumnKind;
use crate::shared::bed::GroupedWindows;
use serde_json::Value;
use std::sync::Arc;
use tempfile::TempDir;
use zarrs::{array::Array, filesystem::FilesystemStore};

#[test]
fn dense_end_motif_zarr_writes_counts_motifs_and_window_metadata() {
    let temp = TempDir::new().expect("temp dir should be created");
    let motifs = vec!["_AA".to_string(), "_GG".to_string()];
    let bins = index_label_bins(
        &motifs,
        vec![
            FxHashMap::from_iter([("_AA", 1.0), ("_GG", 2.5)]),
            FxHashMap::from_iter([("_GG", 3.0)]),
        ],
    );
    let bin_info = vec![
        WindowBinInfo {
            chromosome: "chr2".to_string(),
            start: 10,
            end: 20,
            output_index: 0,
            blacklisted_fraction: 0.25,
        },
        WindowBinInfo {
            chromosome: "chr10".to_string(),
            start: 40,
            end: 60,
            output_index: 1,
            blacklisted_fraction: 0.0,
        },
    ];

    let store_path = write_end_motif_zarr(
        temp.path(),
        "sample",
        &bins,
        &motifs,
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Windows {
            bin_info: &bin_info,
            row_mode: EndWindowRowMode::Bed,
        },
        true,
    )
    .expect("dense Zarr output should write");

    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(root_metadata["attributes"]["cfdnalab_schema"], "end_motif_counts");
    // Keep this literal so schema-version bumps require an explicit test update.
    assert_eq!(
        root_metadata["attributes"]["cfdnalab_schema_version"],
        serde_json::json!(2)
    );
    assert_eq!(
        root_metadata["attributes"]["motif_axis_kind"],
        serde_json::json!("motif")
    );
    assert_eq!(root_metadata["attributes"]["storage_mode"], "dense");
    assert_eq!(root_metadata["attributes"]["row_mode"], "bed");
    assert_eq!(root_metadata["attributes"]["primary_array"], "counts");
    assert!(root_metadata["attributes"]["primary_group"].is_null());
    assert_eq!(
        read_f64_array(&store_path, "/counts"),
        vec![1.0, 2.5, 0.0, 3.0]
    );
    assert_eq!(read_i32_array(&store_path, "/motif_index"), vec![0, 1]);
    let motif_metadata = read_json(&store_path.join("motif_index/zarr.json"));
    assert_eq!(motif_metadata["attributes"]["label_array"], "motif_ascii");
    assert_eq!(read_i32_array(&store_path, "/motif_byte"), vec![0, 1, 2]);
    assert_eq!(
        read_u8_array(&store_path, "/motif_ascii"),
        b"_AA_GG".to_vec()
    );
    assert!(!store_path.join("motif").exists());
    assert!(!store_path.join("motif_utf8").exists());
    assert!(!store_path.join("motif_nbytes").exists());
    let chromosome_metadata = read_json(&store_path.join("chromosome/zarr.json"));
    assert_eq!(
        chromosome_metadata["attributes"]["label_field"],
        "chromosome_name"
    );
    assert_eq!(
        chromosome_metadata["attributes"]["labels"],
        serde_json::json!(["chr2", "chr10"])
    );
    assert!(!store_path.join("chromosome_name").exists());
    assert!(!store_path.join("chromosome_name_utf8").exists());
    assert!(!store_path.join("chromosome_name_nbytes").exists());
    assert!(!store_path.join("chromosome_name_byte").exists());
    assert_eq!(read_i32_array(&store_path, "/row_chromosome"), vec![0, 1]);
    assert_eq!(read_i64_array(&store_path, "/row_start_bp"), vec![10, 40]);
    assert_eq!(read_i64_array(&store_path, "/row_end_bp"), vec![20, 60]);
    assert_eq!(
        read_f64_array(&store_path, "/blacklisted_fraction"),
        vec![0.25, 0.0]
    );
}

#[test]
fn sparse_end_motif_zarr_writes_sorted_coo_arrays() {
    let temp = TempDir::new().expect("temp dir should be created");
    let motifs = vec!["_AA".to_string(), "_GG".to_string()];
    let bins = index_label_bins(
        &motifs,
        vec![FxHashMap::from_iter([("_GG", 2.5), ("_AA", 1.0)])],
    );

    let store_path = write_end_motif_zarr(
        temp.path(),
        "sample",
        &bins,
        &motifs,
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Global,
        false,
    )
    .expect("sparse Zarr output should write");

    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(root_metadata["attributes"]["storage_mode"], "sparse_coo");
    assert!(root_metadata["attributes"]["primary_array"].is_null());
    assert_eq!(root_metadata["attributes"]["primary_group"], "sparse");
    assert_eq!(root_metadata["attributes"]["sparse_format"], "coo");
    let sparse_group_metadata = read_json(&store_path.join("sparse/zarr.json"));
    assert_eq!(sparse_group_metadata["node_type"], "group");
    assert_eq!(
        sparse_group_metadata["attributes"]["sparse_format"],
        serde_json::json!("coo")
    );
    assert_eq!(read_i32_array(&store_path, "/sparse/row"), vec![0, 0]);
    assert_eq!(read_i32_array(&store_path, "/sparse/motif"), vec![0, 1]);
    assert_eq!(read_f64_array(&store_path, "/sparse/count"), vec![1.0, 2.5]);
    assert_eq!(read_i32_array(&store_path, "/sparse/shape"), vec![1, 2]);
    let sparse_dimension_metadata = read_json(&store_path.join("sparse/sparse_dimension/zarr.json"));
    assert_eq!(
        sparse_dimension_metadata["attributes"]["labels"],
        serde_json::json!(["row", "motif"])
    );
    assert_eq!(
        read_i32_array(&store_path, "/sparse/sparse_dimension"),
        vec![0, 1]
    );
    assert!(!store_path.join("sparse/dimension").exists());
}

#[test]
fn motif_group_zarr_writes_json_labels_without_motif_ascii() {
    let temp = TempDir::new().expect("temp dir should be created");
    let motif_groups = vec!["short".to_string(), "group-two".to_string()];
    let bins = index_label_bins(
        &motif_groups,
        vec![FxHashMap::from_iter([("short", 1.0), ("group-two", 2.5)])],
    );

    let store_path = write_end_motif_zarr(
        temp.path(),
        "sample",
        &bins,
        &motif_groups,
        SelectedMotifColumnKind::MotifGroup,
        EndMotifRowMetadata::Global,
        false,
    )
    .expect("motif-group Zarr output should write");

    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(
        root_metadata["attributes"]["cfdnalab_schema_version"],
        serde_json::json!(2)
    );
    assert_eq!(
        root_metadata["attributes"]["motif_axis_kind"],
        serde_json::json!("motif_group")
    );
    assert_eq!(read_i32_array(&store_path, "/motif_index"), vec![0, 1]);
    let motif_metadata = read_json(&store_path.join("motif_index/zarr.json"));
    assert_eq!(motif_metadata["attributes"]["label_field"], "motif_group");
    assert_eq!(
        motif_metadata["attributes"]["labels"],
        serde_json::json!(["short", "group-two"])
    );
    assert!(!store_path.join("motif_byte").exists());
    assert!(!store_path.join("motif_ascii").exists());
    assert_eq!(read_i32_array(&store_path, "/sparse/motif"), vec![0, 1]);
}

#[test]
fn selected_motif_zarr_writes_motif_ascii_labels_in_supplied_order() {
    // Arrange: selected ungrouped motifs are already final output labels when they reach the
    // writer. The writer must expose them as concrete motif labels, not as motif groups.
    let temp = TempDir::new().expect("temp dir should be created");
    let motifs = vec!["G_T".to_string(), "A_C".to_string()];
    let bins = index_label_bins(
        &motifs,
        vec![FxHashMap::from_iter([("A_C", 2.0), ("G_T", 1.0)])],
    );

    // Act
    let store_path = write_end_motif_zarr(
        temp.path(),
        "sample",
        &bins,
        &motifs,
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Global,
        false,
    )
    .expect("selected motif Zarr output should write");

    // Assert
    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(
        root_metadata["attributes"]["motif_axis_kind"],
        serde_json::json!("motif")
    );
    assert_eq!(read_i32_array(&store_path, "/motif_index"), vec![0, 1]);
    let motif_metadata = read_json(&store_path.join("motif_index/zarr.json"));
    assert_eq!(motif_metadata["attributes"]["label_array"], "motif_ascii");
    assert!(motif_metadata["attributes"].get("label_field").is_none());
    assert_eq!(read_i32_array(&store_path, "/motif_byte"), vec![0, 1, 2]);
    assert_eq!(
        read_u8_array(&store_path, "/motif_ascii"),
        b"G_TA_C".to_vec()
    );
    assert_eq!(read_i32_array(&store_path, "/sparse/motif"), vec![0, 1]);
    assert_eq!(read_f64_array(&store_path, "/sparse/count"), vec![1.0, 2.0]);
}

#[test]
fn dense_motif_group_zarr_writes_json_labels_and_counts_without_motif_ascii() {
    // Arrange: --all-motifs with a grouped motifs file keeps unobserved groups as zero-count
    // dense columns. Group labels are variable-width strings, so they must stay in JSON metadata.
    let temp = TempDir::new().expect("temp dir should be created");
    let motif_groups = vec![
        "left.hit".to_string(),
        "right-hit".to_string(),
        "unused_group".to_string(),
    ];
    let bins = index_label_bins(
        &motif_groups,
        vec![FxHashMap::from_iter([("left.hit", 2.0), ("right-hit", 1.0)])],
    );

    // Act
    let store_path = write_end_motif_zarr(
        temp.path(),
        "sample",
        &bins,
        &motif_groups,
        SelectedMotifColumnKind::MotifGroup,
        EndMotifRowMetadata::Global,
        true,
    )
    .expect("dense motif-group Zarr output should write");

    // Assert
    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(root_metadata["attributes"]["storage_mode"], "dense");
    assert_eq!(
        root_metadata["attributes"]["motif_axis_kind"],
        serde_json::json!("motif_group")
    );
    let motif_metadata = read_json(&store_path.join("motif_index/zarr.json"));
    assert_eq!(motif_metadata["attributes"]["label_field"], "motif_group");
    assert_eq!(
        motif_metadata["attributes"]["labels"],
        serde_json::json!(["left.hit", "right-hit", "unused_group"])
    );
    assert_eq!(read_i32_array(&store_path, "/motif_index"), vec![0, 1, 2]);
    assert_eq!(read_f64_array(&store_path, "/counts"), vec![2.0, 1.0, 0.0]);
    assert!(!store_path.join("motif_byte").exists());
    assert!(!store_path.join("motif_ascii").exists());
}

#[test]
fn sparse_end_motif_zarr_round_trips_to_dense_matrix_and_metadata() {
    // Arrange:
    // Three rows and three fixed-width motifs exercise sparse row order, motif order, and a
    // missing middle-row motif. The expected dense matrix is derived directly from the bins below
    // in motif order [_AA, _CC, _GG]:
    //
    // row 0: _AA = 1.0,  _CC = 0.0,  _GG = 2.5
    // row 1: _AA = 0.0,  _CC = 4.25, _GG = 0.0
    // row 2: _AA = 0.75, _CC = 0.0,  _GG = 8.0
    let temp = TempDir::new().expect("temp dir should be created");
    let motifs = vec!["_AA".to_string(), "_CC".to_string(), "_GG".to_string()];
    let bins = index_label_bins(
        &motifs,
        vec![
            FxHashMap::from_iter([("_GG", 2.5), ("_AA", 1.0)]),
            FxHashMap::from_iter([("_CC", 4.25)]),
            FxHashMap::from_iter([("_AA", 0.75), ("_CC", 0.0), ("_GG", 8.0)]),
        ],
    );
    let bin_info = vec![
        WindowBinInfo {
            chromosome: "chr2".to_string(),
            start: 10,
            end: 20,
            output_index: 0,
            blacklisted_fraction: 0.25,
        },
        WindowBinInfo {
            chromosome: "chr10".to_string(),
            start: 40,
            end: 60,
            output_index: 1,
            blacklisted_fraction: 0.0,
        },
        WindowBinInfo {
            chromosome: "chr2".to_string(),
            start: 90,
            end: 130,
            output_index: 2,
            blacklisted_fraction: 0.5,
        },
    ];

    // Act
    let store_path = write_end_motif_zarr(
        temp.path(),
        "sample",
        &bins,
        &motifs,
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Windows {
            bin_info: &bin_info,
            row_mode: EndWindowRowMode::Size,
        },
        false,
    )
    .expect("sparse Zarr output should write");

    // Assert: sparse arrays are sorted COO and can reconstruct the original dense matrix.
    assert_eq!(read_i32_array(&store_path, "/sparse/row"), vec![0, 0, 1, 2, 2]);
    assert_eq!(
        read_i32_array(&store_path, "/sparse/motif"),
        vec![0, 2, 1, 0, 2]
    );
    assert_eq!(
        read_f64_array(&store_path, "/sparse/count"),
        vec![1.0, 2.5, 4.25, 0.75, 8.0]
    );
    assert_eq!(read_i32_array(&store_path, "/sparse/shape"), vec![3, 3]);
    assert_eq!(
        read_sparse_counts_as_dense(&store_path),
        vec![1.0, 0.0, 2.5, 0.0, 4.25, 0.0, 0.75, 0.0, 8.0]
    );

    // Assert: sparse metadata is enough to label the reconstructed matrix.
    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(root_metadata["attributes"]["storage_mode"], "sparse_coo");
    assert_eq!(root_metadata["attributes"]["row_mode"], "size");
    assert!(root_metadata["attributes"]["primary_array"].is_null());
    assert_eq!(root_metadata["attributes"]["primary_group"], "sparse");
    let sparse_dimension_metadata = read_json(&store_path.join("sparse/sparse_dimension/zarr.json"));
    assert_eq!(
        sparse_dimension_metadata["attributes"]["labels"],
        serde_json::json!(["row", "motif"])
    );
    assert_eq!(
        read_i32_array(&store_path, "/sparse/sparse_dimension"),
        vec![0, 1]
    );
    assert_eq!(read_i32_array(&store_path, "/motif_index"), vec![0, 1, 2]);
    let motif_metadata = read_json(&store_path.join("motif_index/zarr.json"));
    assert_eq!(motif_metadata["attributes"]["label_array"], "motif_ascii");
    assert_eq!(read_i32_array(&store_path, "/motif_byte"), vec![0, 1, 2]);
    assert_eq!(
        read_u8_array(&store_path, "/motif_ascii"),
        b"_AA_CC_GG".to_vec()
    );
    assert!(!store_path.join("motif").exists());
    assert!(!store_path.join("motif_utf8").exists());
    assert!(!store_path.join("sparse/dimension").exists());
    let chromosome_metadata = read_json(&store_path.join("chromosome/zarr.json"));
    assert_eq!(
        chromosome_metadata["attributes"]["labels"],
        serde_json::json!(["chr2", "chr10"])
    );
    assert!(!store_path.join("chromosome_name").exists());
    assert!(!store_path.join("chromosome_name_utf8").exists());
    assert_eq!(read_i32_array(&store_path, "/row_chromosome"), vec![0, 1, 0]);
    assert_eq!(
        read_i64_array(&store_path, "/row_start_bp"),
        vec![10, 40, 90]
    );
    assert_eq!(
        read_i64_array(&store_path, "/row_end_bp"),
        vec![20, 60, 130]
    );
    assert_eq!(
        read_f64_array(&store_path, "/blacklisted_fraction"),
        vec![0.25, 0.0, 0.5]
    );
}

#[test]
fn sparse_end_motif_zarr_allows_no_observed_motifs() {
    let temp = TempDir::new().expect("temp dir should be created");
    let bins: Vec<FxHashMap<u32, f64>> = vec![FxHashMap::default()];
    let motifs = Vec::new();

    let store_path = write_end_motif_zarr(
        temp.path(),
        "empty",
        &bins,
        &motifs,
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Global,
        false,
    )
    .expect("empty sparse Zarr output should write");

    assert_eq!(read_i32_array(&store_path, "/motif_index"), Vec::<i32>::new());
    assert_eq!(read_i32_array(&store_path, "/motif_byte"), Vec::<i32>::new());
    assert_eq!(read_u8_array(&store_path, "/motif_ascii"), Vec::<u8>::new());
    assert!(!store_path.join("motif").exists());
    assert_eq!(read_i32_array(&store_path, "/sparse/row"), Vec::<i32>::new());
    assert_eq!(read_i32_array(&store_path, "/sparse/shape"), vec![1, 0]);
}

#[test]
fn grouped_end_motif_zarr_writes_group_metadata_and_dense_counts() {
    // Arrange:
    // Group rows must stay in count-row order. Variable-width group names exercise label attrs,
    // and an all-zero row checks that empty groups still have metadata.
    let temp = TempDir::new().expect("temp dir should be created");
    let motifs = vec!["_AA".to_string(), "_CT".to_string()];
    let bins = index_label_bins(
        &motifs,
        vec![
            FxHashMap::from_iter([("_AA", 1.0)]),
            FxHashMap::from_iter([("_CT", 2.5)]),
            FxHashMap::default(),
        ],
    );
    let groups = vec![
        EndGroupSummary {
            group_idx: 0,
            group_name: "A",
            eligible_windows: 1,
            blacklisted_fraction: 0.0,
        },
        EndGroupSummary {
            group_idx: 1,
            group_name: "long_group",
            eligible_windows: 2,
            blacklisted_fraction: 0.125,
        },
        EndGroupSummary {
            group_idx: 2,
            group_name: "mid",
            eligible_windows: 0,
            blacklisted_fraction: 0.0,
        },
    ];

    // Act
    let store_path = write_end_motif_zarr(
        temp.path(),
        "sample",
        &bins,
        &motifs,
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Groups(groups),
        true,
    )
    .expect("grouped dense Zarr output should write");

    // Assert
    let root_metadata = read_json(&store_path.join("zarr.json"));
    assert_eq!(root_metadata["attributes"]["storage_mode"], "dense");
    assert_eq!(root_metadata["attributes"]["row_mode"], "grouped_bed");
    assert_eq!(root_metadata["attributes"]["primary_array"], "counts");
    assert!(root_metadata["attributes"]["primary_group"].is_null());
    assert_eq!(
        read_f64_array(&store_path, "/counts"),
        vec![1.0, 0.0, 0.0, 2.5, 0.0, 0.0]
    );
    assert_eq!(read_i32_array(&store_path, "/group"), vec![0, 1, 2]);
    let group_metadata = read_json(&store_path.join("group/zarr.json"));
    assert_eq!(group_metadata["attributes"]["label_field"], "group_name");
    assert_eq!(
        group_metadata["attributes"]["labels"],
        serde_json::json!(["A", "long_group", "mid"])
    );
    assert!(!store_path.join("group_name").exists());
    assert!(!store_path.join("group_name_utf8").exists());
    assert!(!store_path.join("group_name_nbytes").exists());
    assert!(!store_path.join("group_name_byte").exists());
    assert_eq!(read_i32_array(&store_path, "/eligible_windows"), vec![1, 2, 0]);
    assert_eq!(
        read_f64_array(&store_path, "/blacklisted_fraction"),
        vec![0.0, 0.125, 0.0]
    );
}

#[test]
fn end_motif_zarr_rejects_control_characters_in_public_labels() {
    // Arrange: public labels are used in JSON attrs and dataframe columns downstream, so the writer
    // rejects control characters instead of silently rewriting them.
    let temp = TempDir::new().expect("temp dir should be created");
    let bad_motifs = vec!["bad\tmotif".to_string()];
    let bins = index_label_bins(
        &bad_motifs,
        vec![FxHashMap::from_iter([("bad\tmotif", 1.0)])],
    );

    // Act
    let motif_error = write_end_motif_zarr(
        temp.path(),
        "bad_motif",
        &bins,
        &bad_motifs,
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Global,
        false,
    )
    .expect_err("control characters in motif labels should fail");

    // Assert
    assert!(
        motif_error
            .to_string()
            .contains("motif contains a control character")
    );

    // Arrange: chromosome names use the same public-label rule.
    let good_motifs = vec!["_AA".to_string()];
    let good_bins = index_label_bins(&good_motifs, vec![FxHashMap::from_iter([("_AA", 1.0)])]);
    let bin_info = vec![WindowBinInfo {
        chromosome: "chr1\nbad".to_string(),
        start: 10,
        end: 20,
        output_index: 0,
        blacklisted_fraction: 0.0,
    }];

    // Act
    let chromosome_error = write_end_motif_zarr(
        temp.path(),
        "bad_chromosome",
        &good_bins,
        &good_motifs,
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Windows {
            bin_info: &bin_info,
            row_mode: EndWindowRowMode::Bed,
        },
        false,
    )
    .expect_err("control characters in chromosome labels should fail");

    // Assert
    assert!(
        chromosome_error
            .to_string()
            .contains("chromosome_name contains a control character")
    );

    // Arrange: grouped row names follow the same public-label rule even when the helper that
    // normally builds grouped metadata is bypassed.
    let bad_groups = vec![EndGroupSummary {
        group_idx: 0,
        group_name: "bad\ngroup",
        eligible_windows: 1,
        blacklisted_fraction: 0.0,
    }];

    // Act
    let group_error = write_end_motif_zarr(
        temp.path(),
        "bad_group",
        &good_bins,
        &good_motifs,
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Groups(bad_groups),
        false,
    )
    .expect_err("control characters in group names should fail");

    // Assert
    assert!(
        group_error
            .to_string()
            .contains("group_name contains a control character")
    );
}

#[test]
fn grouped_end_row_metadata_keeps_count_row_order_and_blacklist_fractions() {
    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(1, "beta".to_string());
    group_idx_to_name.insert(0, "alpha".to_string());
    let grouped_windows = GroupedWindows::from_tuples(
        &[
            (10, 20, 1),
            (20, 30, 0),
            (30, 50, 1),
        ],
        None,
    )
    .expect("grouped windows should build");
    let grouped_windows_map = FxHashMap::from_iter([("chr1".to_string(), grouped_windows)]);
    let blacklist_map = FxHashMap::from_iter([(
        "chr1".to_string(),
        vec![Interval::new(10, 15).expect("blacklist interval")],
    )]);

    let summaries = grouped_end_row_metadata(
        &group_idx_to_name,
        &["chr1".to_string()],
        &grouped_windows_map,
        &blacklist_map,
    )
    .expect("group metadata should build");

    assert_eq!(summaries.len(), 2);
    assert_eq!(summaries[0].group_idx, 0);
    assert_eq!(summaries[0].group_name, "alpha");
    assert_eq!(summaries[0].eligible_windows, 1);
    assert_eq!(summaries[0].blacklisted_fraction, 0.0);
    assert_eq!(summaries[1].group_idx, 1);
    assert_eq!(summaries[1].group_name, "beta");
    assert_eq!(summaries[1].eligible_windows, 2);
    assert_eq!(summaries[1].blacklisted_fraction, 5.0 / 30.0);
}

#[test]
fn grouped_end_row_metadata_rejects_non_contiguous_group_indices() {
    let group_idx_to_name =
        FxHashMap::from_iter([(0, "alpha".to_string()), (2, "beta".to_string())]);
    let error = grouped_end_row_metadata(
        &group_idx_to_name,
        &[],
        &FxHashMap::default(),
        &FxHashMap::default(),
    )
    .expect_err("non-contiguous group_idx values should fail");

    assert!(
        error
            .to_string()
            .contains("end-motif group indices must match count rows 0..1 but observed [0, 2]")
    );
}

#[test]
fn stack_end_motif_counts_places_values_in_expected_rows_and_columns() {
    let motifs = vec!["_AA".to_string(), "_GG".to_string()];
    let bins = index_label_bins(
        &motifs,
        vec![
            FxHashMap::from_iter([("_AA", 1.0), ("_GG", 2.5)]),
            FxHashMap::from_iter([("_GG", 3.0)]),
        ],
    );

    let counts = stack_end_motif_counts(&bins, 2).expect("dense counts should build");

    assert_eq!(counts.shape(), &[2, 2]);
    assert_eq!(counts[(0, 0)], 1.0);
    assert_eq!(counts[(0, 1)], 2.5);
    assert_eq!(counts[(1, 0)], 0.0);
    assert_eq!(counts[(1, 1)], 3.0);
}

#[test]
fn motif_metadata_rejects_variable_width_or_non_ascii_labels() {
    let temp = TempDir::new().expect("temp dir should be created");
    let bins: Vec<FxHashMap<u32, f64>> = vec![FxHashMap::default()];

    let variable_width_error = write_end_motif_zarr(
        temp.path(),
        "variable_width",
        &bins,
        &["_AA".to_string(), "_AAAA".to_string()],
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Global,
        false,
    )
    .expect_err("variable-width motif labels should fail");
    assert!(
        variable_width_error
            .to_string()
            .contains("end-motif labels must have one fixed ASCII width")
    );

    let non_ascii_error = write_end_motif_zarr(
        temp.path(),
        "non_ascii",
        &bins,
        &["_ÅA".to_string()],
        SelectedMotifColumnKind::Motif,
        EndMotifRowMetadata::Global,
        false,
    )
    .expect_err("non-ASCII motif labels should fail");
    assert!(
        non_ascii_error
            .to_string()
            .contains("motif labels must be ASCII")
    );
}

#[test]
fn dense_count_chunk_shape_keeps_small_matrix_in_one_chunk() {
    let shape = [2, 4];

    let chunk_shape = dense_count_chunk_shape(shape).expect("small shape should be valid");

    assert_eq!(chunk_shape, shape);
}

#[test]
fn dense_count_chunk_shape_chunks_large_matrix_by_row() {
    let shape = [10_000, 1_000];

    let chunk_shape = dense_count_chunk_shape(shape).expect("large shape should be valid");

    assert_eq!(chunk_shape, [2_000, 1_000]);
}

#[test]
fn dense_count_chunk_shape_limits_wide_motif_axis() {
    let shape = [2, 3_000_000];

    let chunk_shape = dense_count_chunk_shape(shape).expect("wide shape should be valid");

    assert_eq!(chunk_shape, [1, 2_000_000]);
}

fn index_label_bins(
    column_labels: &[String],
    bins: Vec<FxHashMap<&str, f64>>,
) -> Vec<FxHashMap<u32, f64>> {
    let column_by_label: FxHashMap<&str, u32> = column_labels
        .iter()
        .enumerate()
        .map(|(column_idx, label)| {
            (
                label.as_str(),
                u32::try_from(column_idx).expect("test column index should fit in u32"),
            )
        })
        .collect();

    bins.into_iter()
        .map(|bin| {
            bin.into_iter()
                .map(|(label, count)| {
                    let column_idx = column_by_label.get(label).copied().unwrap_or_else(|| {
                        panic!("test bin label `{label}` is missing from {column_labels:?}")
                    });
                    (column_idx, count)
                })
                .collect()
        })
        .collect()
}

fn read_json(path: &std::path::Path) -> Value {
    let text = std::fs::read_to_string(path).expect("JSON should be readable");
    serde_json::from_str(&text).expect("JSON should parse")
}

fn read_f64_array(store_path: &std::path::Path, array_path: &str) -> Vec<f64> {
    read_array(store_path, array_path)
}

fn read_i32_array(store_path: &std::path::Path, array_path: &str) -> Vec<i32> {
    read_array(store_path, array_path)
}

fn read_i64_array(store_path: &std::path::Path, array_path: &str) -> Vec<i64> {
    read_array(store_path, array_path)
}

fn read_u8_array(store_path: &std::path::Path, array_path: &str) -> Vec<u8> {
    read_array(store_path, array_path)
}

fn read_sparse_counts_as_dense(store_path: &std::path::Path) -> Vec<f64> {
    let row_indices = read_i32_array(store_path, "/sparse/row");
    let motif_indices = read_i32_array(store_path, "/sparse/motif");
    let counts = read_f64_array(store_path, "/sparse/count");
    let shape = read_i32_array(store_path, "/sparse/shape");
    assert_eq!(shape.len(), 2, "sparse shape should have row and motif dimensions");
    let n_rows = shape[0] as usize;
    let n_motifs = shape[1] as usize;
    let mut dense = vec![0.0; n_rows * n_motifs];
    for ((row, motif), count) in row_indices.iter().zip(motif_indices.iter()).zip(counts.iter()) {
        let dense_index = (*row as usize) * n_motifs + (*motif as usize);
        dense[dense_index] = *count;
    }
    dense
}

fn read_array<T>(store_path: &std::path::Path, array_path: &str) -> Vec<T>
where
    T: zarrs::array::ElementOwned,
{
    let store = Arc::new(FilesystemStore::new(store_path).expect("Zarr store should open"));
    let array = Array::open(store, array_path).expect("Zarr array should open");
    array
        .retrieve_array_subset(&array.subset_all())
        .expect("Zarr array should read")
}
