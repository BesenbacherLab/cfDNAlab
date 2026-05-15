use super::*;
use crate::{
    commands::{
        cli_common::{AssignToWindowArgs, ChromosomeArgs, IOCArgs},
        gc_bias::correct::{GCLengthRange, MarginalizeLengthsWeightingScheme},
        lengths::config::DEFAULT_OUTPUT_DECIMALS,
    },
    shared::{
        bed::GroupedWindows,
        clip_mode::ClipMode,
        indel_mode::IndelMode,
        interval::{IndexedInterval, Interval},
        windowing::WindowBinInfo,
    },
};
use fxhash::FxHashMap;
use serde_json::{Value, json};
use std::path::PathBuf;
use std::sync::Arc;
use tempfile::TempDir;

fn minimal_config(output_dir: PathBuf) -> LengthsConfig {
    LengthsConfig::new(
        IOCArgs {
            bam: output_dir.join("input.bam"),
            output_dir,
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    )
}

fn read_settings(output_dir: &std::path::Path, prefix: &str) -> Value {
    let settings_path = settings_path(output_dir, prefix);
    let settings_text =
        std::fs::read_to_string(settings_path).expect("settings sidecar should be readable");
    serde_json::from_str(&settings_text).expect("settings sidecar should be valid JSON")
}

fn settings_path(output_dir: &std::path::Path, prefix: &str) -> PathBuf {
    output_dir.join(crate::shared::io::dot_join(&[
        prefix,
        "fragment_length_settings.json",
    ]))
}

fn length_counts(axis: &Arc<LengthAxis>, counts: &[f64]) -> LengthCounts {
    LengthCounts {
        counts: counts.to_vec(),
        axis: Arc::clone(axis),
    }
}

fn read_text(path: &std::path::Path) -> String {
    std::fs::read_to_string(path).expect("TSV output should be readable")
}

#[test]
fn length_counts_tsv_writes_global_count_columns_without_row_key() {
    // Arrange:
    // Global output has no domain key. The row should contain only count columns whose names encode
    // the half-open length-bin intervals.
    let temp = TempDir::new().expect("tempdir");
    let output_path = temp.path().join("sample.length_counts.tsv");
    let axis = Arc::new(LengthAxis::new(vec![30, 31, 40]).expect("valid length axis"));
    let counts = vec![length_counts(&axis, &[12.0, 3.5])];

    // Act
    write_length_counts_tsv(
        &output_path,
        &counts,
        &axis,
        DEFAULT_OUTPUT_DECIMALS,
        LengthCountRowMetadata::Global,
    )
    .expect("global TSV should write");

    // Assert
    assert_eq!(
        read_text(&output_path),
        "count_30\tcount_31_40\n12\t3.5\n"
    );
}

#[test]
fn length_counts_tsv_writes_window_coordinates_without_row_index() {
    // Arrange:
    // Windowed output is keyed by genomic coordinates. The internal output index is intentionally
    // absent from the public table.
    let temp = TempDir::new().expect("tempdir");
    let output_path = temp.path().join("sample.length_counts.tsv");
    let axis = Arc::new(LengthAxis::new(vec![30, 31, 40]).expect("valid length axis"));
    let counts = vec![
        length_counts(&axis, &[12.0, 3.5]),
        length_counts(&axis, &[0.25, 7.0]),
    ];
    let windows = vec![
        WindowBinInfo {
            chromosome: "chr1".to_string(),
            start: 0,
            end: 100,
            output_index: 10,
            blacklisted_fraction: 0.25,
        },
        WindowBinInfo {
            chromosome: "chr1".to_string(),
            start: 100,
            end: 200,
            output_index: 11,
            blacklisted_fraction: 0.0,
        },
    ];

    // Act
    write_length_counts_tsv(
        &output_path,
        &counts,
        &axis,
        DEFAULT_OUTPUT_DECIMALS,
        LengthCountRowMetadata::Windows {
            windows: &windows,
            include_blacklisted_fraction: true,
        },
    )
    .expect("window TSV should write");

    // Assert
    assert_eq!(
        read_text(&output_path),
        concat!(
            "chrom\tstart\tend\tblacklisted_fraction\tcount_30\tcount_31_40\n",
            "chr1\t0\t100\t0.25\t12\t3.5\n",
            "chr1\t100\t200\t0\t0.25\t7\n",
        )
    );
}

#[test]
fn length_counts_tsv_rounds_counts_to_requested_decimals_and_blacklist_fraction_to_three() {
    // Arrange:
    // The writer rounds only the text representation. `--decimals` controls count columns, while
    // blacklist fractions keep enough precision to remain useful when counts are written as
    // integers.
    let temp = TempDir::new().expect("tempdir");
    let output_path = temp.path().join("sample.length_counts.tsv");
    let axis = Arc::new(LengthAxis::new(vec![30, 31, 40]).expect("valid length axis"));
    let counts = vec![length_counts(&axis, &[12.345_678, 0.004])];
    let windows = vec![WindowBinInfo {
        chromosome: "chr1".to_string(),
        start: 0,
        end: 100,
        output_index: 10,
        blacklisted_fraction: 1.0 / 3.0,
    }];

    // Act
    write_length_counts_tsv(
        &output_path,
        &counts,
        &axis,
        0,
        LengthCountRowMetadata::Windows {
            windows: &windows,
            include_blacklisted_fraction: true,
        },
    )
    .expect("window TSV should write");

    // Assert
    assert_eq!(
        read_text(&output_path),
        concat!(
            "chrom\tstart\tend\tblacklisted_fraction\tcount_30\tcount_31_40\n",
            "chr1\t0\t100\t0.333\t12\t0\n",
        )
    );
}

#[test]
fn length_counts_tsv_writes_group_names_and_eligible_windows_without_group_index() {
    // Arrange:
    // Grouped output is keyed by group name. `eligible_windows` gives the denominator for users who
    // want per-window means, including groups that have no retained windows.
    let temp = TempDir::new().expect("tempdir");
    let output_path = temp.path().join("sample.length_counts.tsv");
    let axis = Arc::new(LengthAxis::new(vec![30, 31, 40]).expect("valid length axis"));
    let counts = vec![
        length_counts(&axis, &[12.0, 3.5]),
        length_counts(&axis, &[0.0, 0.0]),
        length_counts(&axis, &[0.25, 7.0]),
    ];

    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(2, "groupC".to_string());
    group_idx_to_name.insert(0, "groupA".to_string());
    group_idx_to_name.insert(1, "groupWithoutWindows".to_string());

    let mut grouped_windows_map = FxHashMap::default();
    grouped_windows_map.insert(
        "chr1".to_string(),
        GroupedWindows::new(
            vec![
                IndexedInterval::new(10, 20, 0_u64).expect("valid grouped window"),
                IndexedInterval::new(30, 40, 0_u64).expect("valid grouped window"),
                IndexedInterval::new(50, 60, 2_u64).expect("valid grouped window"),
            ],
            None,
        ),
    );
    let chromosomes = vec!["chr1".to_string()];
    let blacklist_map: FxHashMap<String, Vec<Interval<u64>>> = FxHashMap::default();

    // Act
    write_length_counts_tsv(
        &output_path,
        &counts,
        &axis,
        DEFAULT_OUTPUT_DECIMALS,
        LengthCountRowMetadata::Groups {
            group_idx_to_name: &group_idx_to_name,
            chromosomes: &chromosomes,
            grouped_windows_map: &grouped_windows_map,
            blacklist_map: &blacklist_map,
            include_blacklisted_fraction: false,
        },
    )
    .expect("grouped TSV should write");

    // Assert
    assert_eq!(
        read_text(&output_path),
        concat!(
            "group_name\teligible_windows\tcount_30\tcount_31_40\n",
            "groupA\t2\t12\t3.5\n",
            "groupWithoutWindows\t0\t0\t0\n",
            "groupC\t1\t0.25\t7\n",
        )
    );
}

#[test]
fn length_counts_tsv_errors_when_group_indices_are_not_count_row_indices() {
    // Arrange:
    // Grouped count rows are indexed by group_idx. A non-contiguous mapping would place
    // `groupC` metadata next to row 1 counts even though its internal index is 2.
    let temp = TempDir::new().expect("tempdir");
    let output_path = temp.path().join("sample.length_counts.tsv");
    let axis = Arc::new(LengthAxis::new(vec![30, 31]).expect("valid length axis"));
    let counts = vec![
        length_counts(&axis, &[12.0]),
        length_counts(&axis, &[0.25]),
    ];

    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(0, "groupA".to_string());
    group_idx_to_name.insert(2, "groupC".to_string());

    let mut grouped_windows_map = FxHashMap::default();
    grouped_windows_map.insert(
        "chr1".to_string(),
        GroupedWindows::new(
            vec![
                IndexedInterval::new(10, 20, 0_u64).expect("valid grouped window"),
                IndexedInterval::new(50, 60, 2_u64).expect("valid grouped window"),
            ],
            None,
        ),
    );
    let chromosomes = vec!["chr1".to_string()];
    let blacklist_map: FxHashMap<String, Vec<Interval<u64>>> = FxHashMap::default();

    // Act
    let err = write_length_counts_tsv(
        &output_path,
        &counts,
        &axis,
        DEFAULT_OUTPUT_DECIMALS,
        LengthCountRowMetadata::Groups {
            group_idx_to_name: &group_idx_to_name,
            chromosomes: &chromosomes,
            grouped_windows_map: &grouped_windows_map,
            blacklist_map: &blacklist_map,
            include_blacklisted_fraction: false,
        },
    )
    .expect_err("non-contiguous grouped indices should not be written silently");

    // Assert
    assert!(
        err.to_string()
            .contains("grouped length count row 1 corresponds to group_idx 2")
    );
}

#[test]
fn length_counts_tsv_errors_when_grouped_window_references_unknown_group() {
    // Arrange:
    // The grouped window map should not contain group indices missing from the public name map.
    let temp = TempDir::new().expect("tempdir");
    let output_path = temp.path().join("sample.length_counts.tsv");
    let axis = Arc::new(LengthAxis::new(vec![30, 31]).expect("valid length axis"));
    let counts = vec![length_counts(&axis, &[12.0])];

    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(0, "groupA".to_string());

    let mut grouped_windows_map = FxHashMap::default();
    grouped_windows_map.insert(
        "chr1".to_string(),
        GroupedWindows::new(
            vec![IndexedInterval::new(10, 20, 1_u64).expect("valid grouped window")],
            None,
        ),
    );
    let chromosomes = vec!["chr1".to_string()];
    let blacklist_map: FxHashMap<String, Vec<Interval<u64>>> = FxHashMap::default();

    // Act
    let err = write_length_counts_tsv(
        &output_path,
        &counts,
        &axis,
        DEFAULT_OUTPUT_DECIMALS,
        LengthCountRowMetadata::Groups {
            group_idx_to_name: &group_idx_to_name,
            chromosomes: &chromosomes,
            grouped_windows_map: &grouped_windows_map,
            blacklist_map: &blacklist_map,
            include_blacklisted_fraction: false,
        },
    )
    .expect_err("unknown grouped window indices should not be dropped silently");

    // Assert
    assert!(
        err.to_string()
            .contains("grouped window references group_idx 1")
    );
}

#[test]
fn length_counts_tsv_errors_when_group_name_cannot_be_written_as_one_tsv_cell() {
    // Arrange:
    // Rewriting tabs or newlines would make distinct group names indistinguishable.
    let temp = TempDir::new().expect("tempdir");
    let output_path = temp.path().join("sample.length_counts.tsv");
    let axis = Arc::new(LengthAxis::new(vec![30, 31]).expect("valid length axis"));
    let counts = vec![length_counts(&axis, &[12.0])];

    let mut group_idx_to_name = FxHashMap::default();
    group_idx_to_name.insert(0, "group\tA".to_string());

    let mut grouped_windows_map = FxHashMap::default();
    grouped_windows_map.insert(
        "chr1".to_string(),
        GroupedWindows::new(
            vec![IndexedInterval::new(10, 20, 0_u64).expect("valid grouped window")],
            None,
        ),
    );
    let chromosomes = vec!["chr1".to_string()];
    let blacklist_map: FxHashMap<String, Vec<Interval<u64>>> = FxHashMap::default();

    // Act
    let err = write_length_counts_tsv(
        &output_path,
        &counts,
        &axis,
        DEFAULT_OUTPUT_DECIMALS,
        LengthCountRowMetadata::Groups {
            group_idx_to_name: &group_idx_to_name,
            chromosomes: &chromosomes,
            grouped_windows_map: &grouped_windows_map,
            blacklist_map: &blacklist_map,
            include_blacklisted_fraction: false,
        },
    )
    .expect_err("invalid TSV group names should fail");

    // Assert
    assert!(
        err.to_string()
            .contains("group_name contains a control character")
    );
}

#[test]
fn settings_writer_records_non_default_interpretation_fields() {
    // Arrange:
    // The sidecar should keep output interpretation fields, not ordinary filters.
    // This configuration exercises the non-default names and boolean indicators.
    let out_dir = TempDir::new().expect("tempdir");
    let mut config = minimal_config(out_dir.path().to_path_buf());
    config.set_indel_mode(IndelMode::Adjust);
    config.clip_mode = ClipMode::Adjust;
    config.max_soft_clips = 7;
    config.max_deletion_bases = 11;
    config.set_decimals(4);
    config.set_window_assignment(AssignToWindowArgs {
        assign_by: WindowAssigner::Proportion(0.5),
    });
    config.set_gc_length_weighting(MarginalizeLengthsWeightingScheme::Frequency);
    config.set_gc_length_range(GCLengthRange::Package);
    config.set_gc_length_trim_rare(0.05);
    config.gc.gc_file = Some(PathBuf::from("gc_bias_correction.npz"));
    config.scale_genome.scaling_factors = Some(PathBuf::from("scaling.tsv"));
    config.blacklist = Some(vec![PathBuf::from("blacklist.bed")]);
    config.set_min_mapq(17);

    let length_axis =
        LengthAxis::new(vec![10, 20, 31]).expect("test length axis should be valid");
    let window_spec = DistributionWindowSpec::Size(100);

    // Act
    write_fragment_length_settings_json(
        &settings_path(out_dir.path(), "sample"),
        &config,
        &window_spec,
        &length_axis,
    )
        .expect("settings sidecar should write");
    let settings = read_settings(out_dir.path(), "sample");

    // Assert
    assert_eq!(settings["length_axis"]["column_intervals"], json!("half_open"));
    assert_eq!(settings["length_axis"]["min_fragment_length"], json!(10));
    assert_eq!(settings["length_axis"]["max_fragment_length"], json!(30));
    assert_eq!(settings["length_axis"]["n_bins"], json!(2));
    assert_eq!(settings["length_axis"]["single_bp_bins"], json!(false));
    assert_eq!(
        settings["length_axis"]["bin_definition"],
        json!({"kind": "explicit_edges", "edges": [10, 20, 31]})
    );
    assert!(settings["length_axis"].get("edges").is_none());
    assert_eq!(settings["aggregation_level"], json!("windows"));
    assert!(settings.get("row_semantics").is_none());
    assert_eq!(settings["window_mode"], json!("by-size"));
    assert_eq!(settings["indel_mode"], json!("adjust"));
    assert_eq!(settings["clip_mode"], json!("adjust"));
    assert_eq!(settings["max_soft_clips"], json!(7));
    assert_eq!(settings["max_deletion_bases"], json!(11));
    assert_eq!(settings["assign_by"], json!("proportion=0.5"));
    assert_eq!(settings["decimals"], json!(4));
    assert_eq!(settings["gc_length_weighting"], json!("frequency"));
    assert_eq!(settings["gc_length_range"], json!("package"));
    assert_eq!(settings["gc_length_trim_rare"], json!(0.05));
    assert_eq!(settings["gc_correction_used"], json!(true));
    assert_eq!(settings["scaling_factors_used"], json!(true));
    assert_eq!(settings["blacklist_used"], json!(true));
    assert!(settings.get("min_mapq").is_none());
}

#[test]
fn settings_writer_uses_compact_stepped_range_definition_for_dense_default() {
    // Arrange:
    // The default per-bp axis has many edges. The sidecar should describe it
    // as a compact stepped range instead of writing every edge value.
    let out_dir = TempDir::new().expect("tempdir");
    let config = minimal_config(out_dir.path().to_path_buf());
    let edges: Vec<u32> = (30..=1001).collect();
    let length_axis = LengthAxis::new(edges).expect("test length axis should be valid");

    // Act
    write_fragment_length_settings_json(
        &settings_path(out_dir.path(), "dense"),
        &config,
        &DistributionWindowSpec::Global,
        &length_axis,
    )
    .expect("dense settings sidecar should write");
    let settings = read_settings(out_dir.path(), "dense");

    // Assert
    assert_eq!(settings["length_axis"]["n_bins"], json!(971));
    assert_eq!(settings["length_axis"]["single_bp_bins"], json!(true));
    assert_eq!(
        settings["length_axis"]["bin_definition"],
        json!({"kind": "stepped_range", "start": 30, "end": 1001, "step": 1})
    );
    assert!(settings["length_axis"].get("edges").is_none());
    assert_eq!(settings["gc_length_trim_rare"], json!(0.0));
    assert_eq!(settings["blacklist_used"], json!(false));
    assert_eq!(settings["decimals"], json!(DEFAULT_OUTPUT_DECIMALS));
}

#[test]
fn settings_writer_records_bed_and_grouped_window_semantics() {
    // Arrange
    let out_dir = TempDir::new().expect("tempdir");
    let config = minimal_config(out_dir.path().to_path_buf());
    let length_axis =
        LengthAxis::new(vec![30, 31, 32]).expect("test length axis should be valid");

    // Act
    write_fragment_length_settings_json(
        &settings_path(out_dir.path(), "bed"),
        &config,
        &DistributionWindowSpec::Bed(PathBuf::from("windows.bed")),
        &length_axis,
    )
    .expect("BED settings sidecar should write");
    write_fragment_length_settings_json(
        &settings_path(out_dir.path(), "groups"),
        &config,
        &DistributionWindowSpec::GroupedBed(PathBuf::from("groups.bed")),
        &length_axis,
    )
    .expect("grouped settings sidecar should write");

    let bed_settings = read_settings(out_dir.path(), "bed");
    let grouped_settings = read_settings(out_dir.path(), "groups");

    // Assert
    assert_eq!(bed_settings["aggregation_level"], json!("windows"));
    assert_eq!(bed_settings["window_mode"], json!("by-bed"));
    assert_eq!(grouped_settings["aggregation_level"], json!("groups"));
    assert_eq!(grouped_settings["window_mode"], json!("by-grouped-bed"));
}

#[test]
fn settings_writer_formats_all_simple_assignment_names() {
    // Arrange / Act / Assert
    assert_eq!(
        window_assigner_name(WindowAssigner::CountOverlap),
        "count-overlap"
    );
    assert_eq!(window_assigner_name(WindowAssigner::Any), "any");
    assert_eq!(window_assigner_name(WindowAssigner::All), "all");
    assert_eq!(window_assigner_name(WindowAssigner::Midpoint), "midpoint");
    assert_eq!(
        window_assigner_name(WindowAssigner::Proportion(0.125)),
        "proportion=0.125"
    );
}
