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
    assert!(proper_pair_error.to_string().contains("--require-proper-pair"));
}

#[test]
fn rejects_model_core_larger_than_the_minimum_context_span() {
    let mut config = config();
    config.set_stride(500_000);
    config.set_bin_size(250_000);

    let error = config
        .validate()
        .expect_err("a fitting context cannot be smaller than its core");

    assert!(error.to_string().contains("stride (500000) cannot be greater"));
}
