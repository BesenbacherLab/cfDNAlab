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
fn context_defaults_and_weighting_choices_are_serialized() {
    let mut config = config();
    assert_eq!(config.bin_size, 7_000_000);
    assert_eq!(config.context_weighting, ContextWeighting::Triangular);

    for (weighting, expected_argument) in [
        (ContextWeighting::Triangular, "triangular"),
        (ContextWeighting::Uniform, "uniform"),
        (ContextWeighting::TriangularBlend(0.2), "triangular=0.2"),
    ] {
        config.set_context_weighting(weighting);
        let arguments = config.to_cli_args().unwrap();
        assert!(
            arguments
                .windows(2)
                .any(|pair| { pair[0] == "--context-weighting" && pair[1] == expected_argument })
        );
        assert!(
            arguments
                .windows(2)
                .any(|pair| { pair[0] == "--bin-size" && pair[1] == "7000000" })
        );
        config.validate().expect("valid context defaults");
    }
}

#[cfg(feature = "cli")]
#[test]
fn cli_context_defaults_match_programmatic_defaults_and_reject_unknown_weighting() {
    use clap::{Args, FromArgMatches};

    let arguments = ["outliers", "--bam", "input.bam", "--output-dir", "out"];
    let command = OutliersConfig::augment_args(clap::Command::new("outliers"));
    let matches = command.clone().try_get_matches_from(arguments).unwrap();
    let parsed = OutliersConfig::from_arg_matches(&matches).unwrap();
    assert_eq!(parsed.bin_size, 7_000_000);
    assert_eq!(parsed.context_weighting, ContextWeighting::Triangular);

    let mut invalid_arguments = arguments.to_vec();
    invalid_arguments.extend(["--context-weighting", "unknown"]);
    let error = command.try_get_matches_from(invalid_arguments).unwrap_err();
    assert_eq!(error.kind(), clap::error::ErrorKind::ValueValidation);
}

#[test]
fn context_weighting_parses_and_roundtrips_blend_fractions_including_endpoints() {
    for (argument, fraction) in [
        ("triangular=0", 0.0),
        ("triangular=0.2", 0.2),
        ("triangular=1", 1.0),
        ("triangular=1e-3", 0.001),
    ] {
        let weighting: ContextWeighting = argument.parse().unwrap();

        assert_eq!(weighting, ContextWeighting::TriangularBlend(fraction));
        assert_eq!(
            weighting.to_string().parse::<ContextWeighting>().unwrap(),
            weighting
        );
        assert!((weighting.triangular_fraction().unwrap() - fraction).abs() < f64::EPSILON);
    }
}

#[test]
fn context_weighting_rejects_invalid_blends_in_cli_and_programmatic_configs() {
    for argument in [
        "triangular=",
        "triangular=text",
        "triangular=NaN",
        "triangular=inf",
        "triangular=-inf",
        "triangular=-0.01",
        "triangular=1.01",
        "uniform=0.2",
    ] {
        assert!(argument.parse::<ContextWeighting>().is_err(), "{argument}");
    }
    let mut config = config();
    for fraction in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.01, 1.01] {
        config.set_context_weighting(ContextWeighting::TriangularBlend(fraction));

        let error = config.validate().unwrap_err();

        assert!(
            error
                .to_string()
                .contains("context weighting blend fraction")
        );
    }
}

#[cfg(feature = "cli")]
#[test]
fn rendered_blended_context_weighting_parses_through_clap() {
    use clap::{Args, FromArgMatches};

    for fraction in [0.0, 0.2, 1.0] {
        let mut config = config();
        config.set_context_weighting(ContextWeighting::TriangularBlend(fraction));
        let arguments = config.to_cli_args().unwrap();
        let command = OutliersConfig::augment_args(clap::Command::new("outliers"));

        // The standalone argument parser starts at the subcommand name after `cfdna`
        let matches = command
            .try_get_matches_from(arguments.into_iter().skip(1))
            .unwrap();
        let parsed = OutliersConfig::from_arg_matches(&matches).unwrap();

        assert_eq!(
            parsed.context_weighting,
            ContextWeighting::TriangularBlend(fraction)
        );
    }
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
fn automatic_multiplier_must_be_finite_and_positive() {
    let mut config = config();
    for multiplier in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0, -1.0] {
        config.set_tail_probability_multiplier(multiplier);
        let error = config.validate().unwrap_err();
        assert!(error.to_string().contains("tail_probability_multiplier"));
    }
    for multiplier in [0.5, 1.0, 2.5] {
        config.set_tail_probability_multiplier(multiplier);
        config
            .validate()
            .expect("positive finite automatic multiplier");
    }
}

#[test]
fn manual_probability_accepts_its_upper_bound_and_rejects_invalid_values() {
    let mut config = config();
    for probability in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0,
        -0.1,
        0.500_001,
    ] {
        config.set_tail_probability(Some(probability));
        let error = config.validate().unwrap_err();
        assert!(error.to_string().contains("tail_probability"));
    }
    for probability in [1e-12, 0.5] {
        config.set_tail_probability(Some(probability));
        config.validate().expect("manual probability in (0, 0.5]");
    }
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
fn rejects_model_core_larger_than_the_minimum_eligible_context_size() {
    let mut config = config();
    config.set_stride(500_000);
    config.set_bin_size(250_000);

    let error = config
        .validate()
        .expect_err("the eligible context minimum must be at least the core size");

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
    #[cfg(feature = "plotters")]
    config
        .validate()
        .expect("plotting is available in this build");
    #[cfg(not(feature = "plotters"))]
    assert!(
        config
            .validate()
            .unwrap_err()
            .to_string()
            .contains("plotters")
    );
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
    // Context coordinates are included in the inspection table even without PNG generation
    config.set_plot_preserved_regions(false);
    config
        .validate()
        .expect("table-only inspection is available in every build");
    for multiplier in [f64::NAN, f64::INFINITY, 0.5] {
        config.set_preserved_plot_context_multiplier(multiplier);
        assert!(config.validate().is_err());
    }
    config.set_preserved_plot_context_multiplier(1.0);
    config.set_preserved_plot_min_span(0);
    assert!(config.validate().is_err());
}

#[cfg(not(feature = "plotters"))]
#[test]
fn unavailable_plotting_is_rejected_before_opening_the_bam() {
    let mut config = config();
    config.set_plot_preserved_regions(true);

    let error = config
        .validate()
        .expect_err("an explicit PNG request must not be ignored");

    assert!(error.to_string().contains("--plot-preserved-regions"));
    assert!(error.to_string().contains("plotters"));
}

#[cfg(feature = "cli")]
#[test]
fn cli_plot_multiplier_and_percentage_parsers_match_validation_bounds() {
    for value in ["0.5", "NaN", "inf", "-1"] {
        assert!(parse_context_multiplier(value).is_err());
    }
    assert_eq!(parse_context_multiplier("1").unwrap(), 1.0);
    for value in ["-1", "100.1", "NaN", "inf"] {
        assert!(parse_percentage(value).is_err());
    }
    assert_eq!(parse_percentage("0").unwrap(), 0.0);
    assert_eq!(parse_percentage("100").unwrap(), 100.0);
}
