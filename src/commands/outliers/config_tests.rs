use super::*;

fn config() -> OutliersConfig {
    OutliersConfig::new(
        IOCArgs {
            bam: PathBuf::from("input.bam"),
            output_dir: PathBuf::from("out"),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    )
}

#[test]
fn rejects_manual_probability_with_nondefault_automatic_multiplier() {
    let mut config = config();
    config.set_tail_probability(Some(0.01));
    config.set_tail_probability_multiplier(2.0);

    let error = config
        .validate()
        .expect_err("manual and automatic threshold controls must conflict");

    assert!(error.to_string().contains("cannot be used together"));
}

#[test]
fn rejects_paired_only_filters_when_reads_are_fragments() {
    let mut ignore_gap_config = config();
    ignore_gap_config.unpaired.reads_are_fragments = true;
    ignore_gap_config.set_ignore_gap(true);
    let ignore_gap_error = ignore_gap_config
        .validate()
        .expect_err("unpaired reads do not have an inter-mate gap");

    let mut proper_pair_config = config();
    proper_pair_config.unpaired.reads_are_fragments = true;
    proper_pair_config.set_require_proper_pair(true);
    let proper_pair_error = proper_pair_config
        .validate()
        .expect_err("unpaired reads cannot require a proper pair");

    assert!(ignore_gap_error.to_string().contains("--ignore-gap"));
    assert!(
        proper_pair_error
            .to_string()
            .contains("--require-proper-pair")
    );
}

#[test]
fn rejects_model_core_larger_than_the_minimum_context_span() {
    let mut config = config();
    config.set_stride(500_000);
    config.set_bin_size(250_000);

    let error = config
        .validate()
        .expect_err("a fitting context cannot be smaller than its core");

    assert!(
        error
            .to_string()
            .contains("stride (500000) cannot be greater")
    );
}

#[test]
fn initial_inspection_defaults_and_model_selection_are_serialized() {
    let mut config = config();
    assert_eq!(config.model, OutlierModel::Zinb);
    assert_eq!(config.protection_window_size, 100_000);
    assert_eq!(config.min_protection_positions, 50_000);
    assert_eq!(config.max_outlier_percentage, 5.0);
    assert_eq!(
        config.ambiguous_region_action,
        AmbiguousRegionAction::Preserve
    );
    config.set_model(OutlierModel::Zip);
    let arguments = config.to_cli_args().unwrap();
    assert!(
        arguments
            .windows(2)
            .any(|pair| pair[0] == "--model" && pair[1] == "zip")
    );
    config.set_min_protection_positions(100_001);
    assert!(config.validate().is_err());
    config.set_min_protection_positions(50_000);
    for percentage in [f64::NAN, -0.1, 100.1] {
        config.set_max_outlier_percentage(percentage);
        assert!(config.validate().is_err());
    }
}

#[test]
fn preserved_plot_settings_are_optional_validated_and_serialized() {
    let mut config = config();
    assert!(!config.plot_preserved_regions);
    config.set_plot_preserved_regions(true);
    config.set_preserved_plot_context_multiplier(2.5);
    config.set_preserved_plot_min_span(20_000);
    assert!(config.validate().is_ok());
    let arguments = config.to_cli_args().unwrap();
    assert!(
        arguments
            .iter()
            .any(|argument| argument == "--plot-preserved-regions")
    );
    assert!(
        arguments
            .windows(2)
            .any(|pair| pair[0] == "--preserved-plot-context-multiplier" && pair[1] == "2.5")
    );
    assert!(
        arguments
            .windows(2)
            .any(|pair| pair[0] == "--preserved-plot-min-span" && pair[1] == "20000")
    );
    for multiplier in [f64::NAN, f64::INFINITY, 0.5] {
        config.set_preserved_plot_context_multiplier(multiplier);
        assert!(config.validate().is_err());
    }
    config.set_preserved_plot_context_multiplier(1.0);
    config.set_preserved_plot_min_span(0);
    assert!(config.validate().is_err());
}
