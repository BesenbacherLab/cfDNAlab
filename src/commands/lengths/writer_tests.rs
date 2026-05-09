use super::*;
use crate::{
    commands::{
        cli_common::{AssignToWindowArgs, ChromosomeArgs, IOCArgs},
        gc_bias::correct::{GCLengthRange, MarginalizeLengthsWeightingScheme},
    },
    shared::{clip_mode::ClipMode, indel_mode::IndelMode},
};
use serde_json::{Value, json};
use std::path::PathBuf;
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
    let settings_path = output_dir.join(dot_join(&[prefix, "fragment_length_settings.json"]));
    let settings_text =
        std::fs::read_to_string(settings_path).expect("settings sidecar should be readable");
    serde_json::from_str(&settings_text).expect("settings sidecar should be valid JSON")
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
    config.set_window_assignment(AssignToWindowArgs {
        assign_by: WindowAssigner::Proportion(0.5),
    });
    config.set_gc_length_weighting(MarginalizeLengthsWeightingScheme::Frequency);
    config.set_gc_length_range(GCLengthRange::Package);
    config.set_gc_length_trim_rare(0.05);
    config.gc.gc_file = Some(PathBuf::from("gc_bias_correction.npz"));
    config.scale_genome.scaling_factors = Some(PathBuf::from("scaling.tsv"));
    config.set_min_mapq(17);

    let length_axis =
        LengthAxis::new(vec![10, 20, 31]).expect("test length axis should be valid");
    let window_spec = DistributionWindowSpec::Size(100);

    // Act
    write_fragment_length_settings_json(&config, &window_spec, &length_axis, "sample")
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
    assert_eq!(settings["gc_length_weighting"], json!("frequency"));
    assert_eq!(settings["gc_length_range"], json!("package"));
    assert_eq!(settings["gc_length_trim_rare"], json!(0.05));
    assert_eq!(settings["gc_correction_used"], json!(true));
    assert_eq!(settings["scaling_factors_used"], json!(true));
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
        &config,
        &DistributionWindowSpec::Global,
        &length_axis,
        "dense",
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
        &config,
        &DistributionWindowSpec::Bed(PathBuf::from("windows.bed")),
        &length_axis,
        "bed",
    )
    .expect("BED settings sidecar should write");
    write_fragment_length_settings_json(
        &config,
        &DistributionWindowSpec::GroupedBed(PathBuf::from("groups.bed")),
        &length_axis,
        "groups",
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
