use super::*;
use crate::commands::{
    cli_common::{ChromosomeArgs, IOCArgs},
    ends::config_structs::AssignMotifToWindowArgs,
};
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

#[test]
fn stack_end_motif_counts_places_values_in_the_expected_rows_and_columns() {
    // Arrange: two windows and two motif columns.
    //
    // Mental derivation:
    // - row 0 contains `_A = 1.0` and `_G = 2.5`
    // - row 1 contains only `_G = 3.0`
    // - any missing row/column pair must stay at the dense default 0.0
    let bins = vec![
        FxHashMap::from_iter([
            ("_A".to_string(), 1.0),
            ("_G".to_string(), 2.5),
        ]),
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

    assert!(err.to_string().contains("missing dense output column for motif label '_C'"));
}

#[test]
fn write_end_settings_json_serializes_proportion_assignment_without_float_noise() {
    // Arrange: `0.1 + 0.2` is the classic floating-point trap and often becomes
    // `0.30000000000000004`. The sidecar should still serialize the human meaning, `0.3`.
    let out_dir = TempDir::new().expect("tempdir");
    let mut cfg = minimal_config(out_dir.path());
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Proportion(0.1 + 0.2),
    });

    write_end_settings_json(out_dir.path(), "ends", &cfg).expect("settings json should write");
    let settings = std::fs::read_to_string(out_dir.path().join("ends.end_motif_settings.json"))
        .expect("settings json should be readable");

    assert!(settings.contains("\"window_assignment\": \"proportion=0.3\""));
}

#[test]
fn format_proportion_threshold_keeps_a_trailing_decimal_for_whole_numbers() {
    // Arrange / Act / Assert:
    // We keep `1.0` and `0.0` rather than collapsing them to `1` and `0` so the value still
    // reads like a proportion threshold in the sidecar.
    assert_eq!(format_proportion_threshold(1.0), "1.0");
    assert_eq!(format_proportion_threshold(0.0), "0.0");
}
