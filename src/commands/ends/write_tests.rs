use super::*;
use crate::{
    commands::{
        cli_common::{ChromosomeArgs, IOCArgs},
        ends::{config_structs::BaseQualityFilter, counting::EndMotifColumnKind},
    },
    shared::indel_mode::IndelMotifFilterPolicy,
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
    k_inside: usize,
    k_outside: usize,
    source_inside: &str,
    clip_strategy: &str,
    window_assignment: &str,
    indel_filter: &str,
    effective_indel_filter: &str,
) -> Value {
    let mut expected = Map::new();
    expected.insert("k_inside".to_string(), json!(k_inside));
    expected.insert("k_outside".to_string(), json!(k_outside));
    expected.insert("all_motifs".to_string(), json!(false));
    expected.insert("motifs_file".to_string(), Value::Null);
    expected.insert("motifs_file_mode".to_string(), Value::Null);
    expected.insert("source_inside".to_string(), json!(source_inside));
    expected.insert("clip_strategy".to_string(), json!(clip_strategy));
    expected.insert("window_assignment".to_string(), json!(window_assignment));
    expected.insert("indel_filter".to_string(), json!(indel_filter));
    expected.insert(
        "effective_indel_filter".to_string(),
        json!(effective_indel_filter),
    );
    #[cfg(feature = "ends_experimental")]
    expected.insert("collapse_complement".to_string(), json!(false));
    Value::Object(expected)
}

#[test]
fn write_end_settings_json_writes_the_minimal_interpretation_sidecar() {
    // Arrange: the minimal default config has
    // - k_inside: 1
    // - k_outside: 0
    // - all_motifs: false
    // - no motifs file
    // - inside source: read
    // - clip strategy: skip
    // - window assignment: endpoint
    // - indel filter: auto
    // - effective indel filter: allow
    // - collapse_complement: false
    // Those are the fields currently retained in the sidecar.
    let out_dir = TempDir::new().expect("tempdir");
    let cfg = minimal_config(out_dir.path());

    // Act
    let settings_path = write_end_settings_json(out_dir.path(), "ends", &cfg, None)
        .expect("settings json should write");
    assert_eq!(settings_path, out_dir.path().join("ends.end_settings.json"));
    let settings =
        std::fs::read_to_string(settings_path).expect("settings json should be readable");

    // Assert
    assert_eq!(
        parse_json(&settings),
        expected_settings_json(1, 0, "read", "skip", "endpoint", "auto", "allow")
    );
}

#[test]
fn write_end_settings_json_includes_non_default_indel_filter() {
    // Arrange: the indel filter changes which motifs survive when indels overlap an end, so the
    // sidecar should keep the configured CLI spelling.
    let out_dir = TempDir::new().expect("tempdir");
    let mut cfg = minimal_config(out_dir.path());
    cfg.indel_filter = IndelMotifFilterPolicy::SkipAffectedFragment;

    // Act
    write_end_settings_json(out_dir.path(), "ends", &cfg, None)
        .expect("settings json should write");
    let settings = std::fs::read_to_string(out_dir.path().join("ends.end_settings.json"))
        .expect("settings json should be readable");
    let parsed = parse_json(&settings);

    // Assert
    assert_eq!(
        parsed.get("indel_filter"),
        Some(&json!("skip-affected-fragment"))
    );
    assert_eq!(
        parsed.get("effective_indel_filter"),
        Some(&json!("skip-affected-fragment"))
    );
}

#[test]
fn write_end_settings_json_resolves_auto_indel_filter_for_reference_inside_bases() {
    // Arrange: auto mode skips only affected ends when inside bases come from the reference.
    let out_dir = TempDir::new().expect("tempdir");
    let mut cfg = minimal_config(out_dir.path());
    cfg.source_inside = KmerSource::Reference;

    // Act
    write_end_settings_json(out_dir.path(), "ends", &cfg, None)
        .expect("settings json should write");
    let settings = std::fs::read_to_string(out_dir.path().join("ends.end_settings.json"))
        .expect("settings json should be readable");
    let parsed = parse_json(&settings);

    // Assert
    assert_eq!(parsed.get("indel_filter"), Some(&json!("auto")));
    assert_eq!(
        parsed.get("effective_indel_filter"),
        Some(&json!("skip-affected-end"))
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
    write_end_settings_json(out_dir.path(), "ends", &cfg, None)
        .expect("settings json should write");
    let settings = std::fs::read_to_string(out_dir.path().join("ends.end_settings.json"))
        .expect("settings json should be readable");
    let parsed = parse_json(&settings);

    // Assert
    assert_eq!(
        parsed.get("bq_filters"),
        Some(&json!(["min in end >= 30", "max in fragment < 20"]))
    );
}

#[test]
fn write_end_settings_json_includes_motifs_file_path_and_mode() {
    // Arrange: motifs-file runs change the count-column meaning, so the sidecar should record both
    // the source file path and whether the parsed file targeted motifs or motif groups.
    let out_dir = TempDir::new().expect("tempdir");
    let mut cfg = minimal_config(out_dir.path());
    let motifs_file = out_dir.path().join("selected_groups.tsv");
    cfg.motifs_file = Some(motifs_file.clone());
    cfg.all_motifs = true;

    // Act
    write_end_settings_json(
        out_dir.path(),
        "ends",
        &cfg,
        Some(EndMotifColumnKind::MotifGroup),
    )
    .expect("settings json should write");
    let settings = std::fs::read_to_string(out_dir.path().join("ends.end_settings.json"))
        .expect("settings json should be readable");
    let parsed = parse_json(&settings);

    // Assert
    assert_eq!(parsed.get("all_motifs"), Some(&json!(true)));
    assert_eq!(
        parsed.get("motifs_file"),
        Some(&json!(motifs_file.to_string_lossy()))
    );
    assert_eq!(parsed.get("motifs_file_mode"), Some(&json!("grouped")));
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
