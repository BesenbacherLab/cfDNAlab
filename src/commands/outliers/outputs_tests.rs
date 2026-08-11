use super::*;
use crate::commands::cli_common::{ChromosomeArgs, IOCArgs};

#[test]
fn global_diagnostics_include_the_maximum_observed_coverage() {
    let observed_counts = [(0, 100), (1, 20), (5, 1)].into_iter().collect();

    let maximum_coverage =
        maximum_observed_coverage(&observed_counts).expect("histogram with an extreme tail");

    assert_eq!(maximum_coverage, 5);
}

#[test]
fn common_metadata_records_threshold_provenance_and_fragment_semantics() {
    let mut config = OutliersConfig::new(
        IOCArgs {
            bam: PathBuf::from("sample.bam"),
            output_dir: PathBuf::from("results"),
            n_threads: 3,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr2".to_string(), "chr1".to_string()]),
            chromosomes_file: None,
        },
    );
    config.set_tail_probability_multiplier(2.5);
    config.set_stride(250_000);
    config.set_bin_size(2_000_000);
    config.set_min_context_fragments(500);
    config.set_blacklist(Some(vec![PathBuf::from("excluded.bed")]));
    let chromosomes = vec!["chr2".to_string(), "chr1".to_string()];
    let metadata = OutputMetadata {
        config: &config,
        chromosomes: &chromosomes,
        tail_probability: 0.000_25,
        eligible_positions: 10_000,
        blacklist_flank: 1_000,
    };
    let mut written = Vec::new();

    write_common_metadata(&mut written, &metadata).expect("valid metadata");
    let text = String::from_utf8(written).expect("UTF-8 metadata");

    assert!(text.contains("# selected_chromosomes=[\"chr2\",\"chr1\"]\n"));
    assert!(text.contains("# tail_probability=0.00025\n"));
    assert!(text.contains("# tail_probability_mode=automatic\n"));
    assert!(text.contains("# tail_probability_multiplier=2.5\n"));
    assert!(text.contains("# automatic_probability_denominator=10000\n"));
    assert!(text.contains("# model_core_size=250000\n"));
    assert!(text.contains("# minimum_context_span=2000000\n"));
    assert!(text.contains("# minimum_context_fragment_support=500\n"));
    assert!(
        text.contains(
            "# fragment_support_definition=accepted_fragments_with_unblacklisted_midpoint\n"
        )
    );
    assert!(text.contains("# raw_coverage_source=uncorrected_mapped_reference_segments\n"));
    assert!(text.contains("# input_blacklists=[\"excluded.bed\"]\n"));
    assert!(text.contains("# blacklist_flank=1000\n"));
}
