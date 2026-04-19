use super::*;
use crate::{
    commands::{
        cli_common::{ChromosomeArgs, IOCArgs},
        ends::config_structs::BaseQualityFilter,
    },
};
use serde_json::{Map, Value, json};
use std::path::Path;
use tempfile::TempDir;

fn minimal_config(output_dir: &Path) -> EndsConfig {
    EndsConfig::new(
        IOCArgs {
            bam: output_dir.join("dummy.bam"),
            output_dir: output_dir.to_path_buf(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        1,
        0,
    )
}

fn parse_json(text: &str) -> Value {
    serde_json::from_str(text).expect("settings sidecar should be valid JSON")
}

fn expected_settings_json(
    source_inside: &str,
    clip_strategy: &str,
    window_assignment: &str,
) -> Value {
    let mut expected = Map::new();
    expected.insert("source_inside".to_string(), json!(source_inside));
    expected.insert("clip_strategy".to_string(), json!(clip_strategy));
    expected.insert("window_assignment".to_string(), json!(window_assignment));
    #[cfg(feature = "ends_experimental")]
    expected.insert("collapse_complement".to_string(), json!(false));
    Value::Object(expected)
}

#[test]
fn stack_end_motif_counts_places_values_in_the_expected_rows_and_columns() {
    // Arrange: two windows and two motif columns.
    //
    // Mental derivation:
    // - row 0 contains `_A = 1.0` and `_G = 2.5`
    // - row 1 contains only `_G = 3.0`
    // - any missing row/column pair must stay at the dense default 0.0
    let bins = vec![
        FxHashMap::from_iter([("_A".to_string(), 1.0), ("_G".to_string(), 2.5)]),
        FxHashMap::from_iter([("_G".to_string(), 3.0)]),
    ];
    let motifs = vec!["_A".to_string(), "_G".to_string()];

    let matrix = stack_end_motif_counts(&bins, &motifs).expect("dense matrix should build");

    assert_eq!(matrix.shape(), &[2, 2]);
    assert_eq!(matrix[(0, 0)], 1.0);
    assert_eq!(matrix[(0, 1)], 2.5);
    assert_eq!(matrix[(1, 0)], 0.0);
    assert_eq!(matrix[(1, 1)], 3.0);
}

#[test]
fn stack_end_motif_counts_errors_when_a_bin_contains_an_unknown_motif() {
    // Arrange: the sparse bin mentions `_C`, but the declared dense column order only contains
    // `_A`. Silently dropping `_C` would corrupt the output, so this must be an error.
    let bins = vec![FxHashMap::from_iter([("_C".to_string(), 1.0)])];
    let motifs = vec!["_A".to_string()];

    let err = stack_end_motif_counts(&bins, &motifs)
        .expect_err("unknown motif labels should fail loudly");

    assert!(
        err.to_string()
            .contains("missing dense output column for motif label '_C'")
    );
}

#[test]
fn write_end_settings_json_writes_the_minimal_interpretation_sidecar() {
    // Arrange: the minimal default config has
    // - inside source: read
    // - clip strategy: skip
    // - window assignment: endpoint
    // - collapse_complement: false
    // Those are the fields currently retained in the sidecar.
    let out_dir = TempDir::new().expect("tempdir");
    let cfg = minimal_config(out_dir.path());

    // Act
    write_end_settings_json(out_dir.path(), "ends", &cfg).expect("settings json should write");
    let settings = std::fs::read_to_string(out_dir.path().join("ends.end_motif_settings.json"))
        .expect("settings json should be readable");

    // Assert
    assert_eq!(
        parse_json(&settings),
        expected_settings_json("read", "skip", "endpoint")
    );
}

#[test]
fn write_end_settings_json_includes_base_quality_filters_when_present() {
    // Arrange: the sidecar should retain any configured base-quality filters because they affect
    // which ends and fragments contribute to the output counts.
    let out_dir = TempDir::new().expect("tempdir");
    let mut cfg = minimal_config(out_dir.path());
    cfg.bq_filter = vec![
        "min in end >= 30"
            .parse::<BaseQualityFilter>()
            .expect("valid end filter"),
        "max in fragment < 20"
            .parse::<BaseQualityFilter>()
            .expect("valid fragment filter"),
    ];

    // Act
    write_end_settings_json(out_dir.path(), "ends", &cfg).expect("settings json should write");
    let settings = std::fs::read_to_string(out_dir.path().join("ends.end_motif_settings.json"))
        .expect("settings json should be readable");
    let parsed = parse_json(&settings);

    // Assert
    assert_eq!(
        parsed.get("bq_filters"),
        Some(&json!(["min in end >= 30", "max in fragment < 20"]))
    );
}

#[test]
fn window_assigner_name_formats_proportion_without_noisy_trailing_precision() {
    // Arrange / Act
    let exact_eighth = window_assigner_name(WindowMotifAssigner::Proportion(0.125));
    let exact_half = window_assigner_name(WindowMotifAssigner::Proportion(0.5));

    // Assert: the sidecar contract should keep simple decimal inputs readable and stable.
    assert_eq!(exact_eighth, "proportion=0.125");
    assert_eq!(exact_half, "proportion=0.5");
}
