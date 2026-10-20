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
fn common_metadata_records_settings_and_analysed_chromosome_lengths() {
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
    let contigs = Contigs {
        contigs: [
            ("chr1".to_string(), (0, 4_000)),
            ("chr2".to_string(), (1, 6_000)),
        ]
        .into_iter()
        .collect(),
    };
    let metadata = OutputMetadata {
        config: &config,
        chromosomes: &chromosomes,
        contigs: &contigs,
        tail_probability: 0.000_25,
        eligible_positions: 10_000,
        blacklist_flank: 1_000,
    };
    let mut written = Vec::new();

    write_common_metadata(&mut written, &metadata).expect("valid metadata");
    let text = String::from_utf8(written).expect("UTF-8 metadata");

    assert!(text.contains("# analysed_chromosomes=[{\"name\":\"chr1\",\"size\":4000},{\"name\":\"chr2\",\"size\":6000}]\n"));
    assert!(!text.contains("selected_chromosomes="));
    assert!(!text.contains("reference_contig_footprint_json="));
    assert!(text.contains("# tail_probability=0.00025\n"));
    assert!(text.contains("# tail_probability_mode=automatic\n"));
    assert!(text.contains("# tail_probability_multiplier=2.5\n"));
    assert!(text.contains("# automatic_probability_denominator=10000\n"));
    assert!(text.contains("# model_core_size=250000\n"));
    assert!(text.contains("# minimum_context_span=2000000\n"));
    assert!(text.contains("# minimum_context_fragment_support=500\n"));
    assert!(text.contains("# raw_coverage_source=uncorrected_mapped_reference_segments\n"));
    assert!(text.contains("# blacklists=[\"excluded.bed\"]\n"));
    assert!(text.contains("# blacklist_flank=1000\n"));
    assert!(text.contains(&format!(
        "# plotting_available={}\n",
        cfg!(feature = "plotters")
    )));
}

#[test]
fn output_manifest_includes_png_only_when_the_build_can_write_it() {
    let paths = OutputPaths::new(Path::new("results"), "sample");
    let mut expected = vec![
        "sample.outliers.keep_weights.tsv",
        "sample.outliers.exact.bed",
        "sample.outliers.flanked.bed",
        "sample.outliers.histograms.tsv.zst",
        "sample.outliers.models.tsv",
        "sample.outliers.global_fit.tsv",
    ];
    if cfg!(feature = "plotters") {
        expected.push("sample.outliers.global_fit.png");
    }
    expected.extend([
        "sample.outliers.audit.tsv",
        "sample.outliers.summary.tsv",
        "sample.outliers.preserved_regions.tsv",
        "sample.outliers.README.txt",
    ]);

    let manifest = paths.all();

    assert_eq!(
        manifest,
        expected
            .iter()
            .map(|name| Path::new("results").join(name))
            .collect::<Vec<_>>()
    );
}
