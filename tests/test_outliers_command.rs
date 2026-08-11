#![cfg(feature = "cmd_outliers")]

use std::{collections::BTreeMap, fs, path::Path};

use anyhow::{Context, Result, ensure};
use cfdnalab::{
    RunOptions,
    run_like_cli::{
        common::{ChromosomeArgs, IOCArgs},
        outliers::{OutlierTarget, OutliersConfig, run_outliers},
    },
    testing::{PairedFragmentSpec, TempBam, TempBamBuilder, read_zst_to_string},
};
use tempfile::TempDir;

fn paired_outlier_fixture() -> Result<TempBam> {
    let mut builder = TempBamBuilder::new()
        .name("paired_outlier")
        .contig("chr1", 1_000)
        .use_record_indexed_read_names();

    // Three non-overlapping baseline fragments create 90 positions at coverage 1.
    for start in [0, 300, 600] {
        builder = builder.paired_fragment(PairedFragmentSpec::new(0, start, 30, 10));
    }
    // Twenty distinct molecules share [100, 130), creating 30 positions at coverage 20.
    for _ in 0..20 {
        builder = builder.paired_fragment(PairedFragmentSpec::new(0, 100, 30, 10));
    }

    builder.build()
}

fn boundary_blacklist_and_fallback_fixture() -> Result<TempBam> {
    let mut builder = TempBamBuilder::new()
        .name("outlier_boundaries")
        .contig("chr1", 1_000_030)
        .contig("chr2", 1_000)
        .use_record_indexed_read_names();

    for start in [100, 300_000, 700_000] {
        builder = builder.paired_fragment(PairedFragmentSpec::new(0, start, 30, 10));
    }
    builder = builder.paired_fragment(PairedFragmentSpec::new(1, 100, 30, 10));
    // The three pileups cross a model-core boundary, a blacklisted core boundary, and a tile
    // boundary, respectively.
    for start in [249_990, 499_990, 999_990] {
        for _ in 0..20 {
            builder = builder.paired_fragment(PairedFragmentSpec::new(0, start, 30, 10));
        }
    }

    builder.build()
}

fn data_lines(text: &str) -> Vec<&str> {
    text.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect()
}

fn base_config(bam: &Path, output_dir: &Path) -> OutliersConfig {
    let mut config = OutliersConfig::new(
        IOCArgs {
            bam: bam.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 2,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    );
    config.set_output_prefix("sample");
    config.set_tile_size(1_000_000);
    config.set_stride(250);
    config.set_bin_size(1_000);
    config.set_min_context_fragments(1);
    config.set_tail_probability(Some(0.05));
    config.set_target(OutlierTarget::Threshold);
    config.set_blacklist_flank(Some(50));
    config.set_min_mapq(0);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 30;
        fragment_lengths.max_fragment_length = 30;
    }
    config
}

#[test]
fn paired_bam_produces_the_hand_derived_models_calls_and_all_outputs() -> Result<()> {
    let bam = paired_outlier_fixture()?;
    let output_dir = TempDir::new()?;
    let config = base_config(bam.bam_path(), output_dir.path());

    let result = run_outliers(&config, RunOptions::new_quiet())?;

    assert_eq!(result.counters.base.counted_fragments, 23);
    assert_eq!(result.output_files.len(), 7);
    let expected_names = [
        "sample.outliers.keep_weights.tsv",
        "sample.outliers.exact.bed",
        "sample.outliers.flanked.bed",
        "sample.outliers.histograms.tsv.zst",
        "sample.outliers.models.tsv",
        "sample.outliers.global_fit.tsv",
        "sample.outliers.global_fit.png",
    ];
    let actual_names = result
        .output_files
        .iter()
        .map(|path| {
            path.file_name()
                .and_then(|name| name.to_str())
                .context("UTF-8 output filename")
        })
        .collect::<Result<Vec<_>>>()?;
    assert_eq!(actual_names, expected_names);
    assert!(result.output_files.iter().all(|path| path.is_file()));

    let histogram_text = read_zst_to_string(&result.output_histograms)?;
    assert!(histogram_text.contains("# selected_chromosomes=[\"chr1\"]\n"));
    assert!(histogram_text.contains(
        "# fragment_support_definition=accepted_fragments_with_unblacklisted_midpoint\n"
    ));
    let histogram_lines = data_lines(&histogram_text);
    assert_eq!(
        histogram_lines[0],
        "chromosome\tstart\tend\tfragment_support\tcoverage\tobserved_positions"
    );
    assert_eq!(
        &histogram_lines[1..],
        &[
            "chr1\t0\t250\t21\t0\t190",
            "chr1\t0\t250\t21\t1\t30",
            "chr1\t0\t250\t21\t20\t30",
            "chr1\t250\t500\t1\t0\t220",
            "chr1\t250\t500\t1\t1\t30",
            "chr1\t500\t750\t1\t0\t220",
            "chr1\t500\t750\t1\t1\t30",
            "chr1\t750\t1000\t0\t0\t250",
        ]
    );

    let model_text = fs::read_to_string(&result.output_models)?;
    let model_lines = data_lines(&model_text);
    let model_header = model_lines[0].split('\t').collect::<Vec<_>>();
    assert_eq!(model_header.len(), 24);
    assert_eq!(model_header[10], "initial_threshold");
    assert_eq!(model_header[15], "final_threshold");
    assert_eq!(model_header[16], "observed_positive_coverage_mean");
    assert_eq!(model_header[17], "observed_positive_coverage_variance");
    assert_eq!(
        model_header[18],
        "underlying_zip_expected_positive_coverage_variance"
    );
    assert_eq!(model_header[19], "positive_coverage_variance_ratio");
    assert_eq!(model_header[20], "observed_initial_tail_positions");
    assert_eq!(model_header[22], "observed_final_tail_positions");
    assert_eq!(model_lines.len(), 5);
    for (core_index, line) in model_lines[1..].iter().enumerate() {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields[0], "chr1");
        assert_eq!(fields[1].parse::<u32>()?, core_index as u32 * 250);
        assert_eq!(fields[2].parse::<u32>()?, (core_index as u32 + 1) * 250);
        assert_eq!(&fields[3..5], &["0", "1000"]);
        assert_eq!(fields[5], "local");
        assert_eq!(fields[6], "1000");
        assert_eq!(fields[7], "23");
        assert_eq!(fields[10], "7");
        assert_eq!(fields[11], "970");
        assert_eq!(fields[15], "2");
        assert_eq!(fields[20], "30");
        assert_eq!(fields[22], "30");
    }

    let global_text = fs::read_to_string(&result.output_global_fit)?;
    assert!(global_text.contains("# positive_coverage_variance_ratio_definition="));
    assert!(global_text.contains("# initial_threshold=7\n"));
    assert!(global_text.contains("# final_threshold=2\n"));
    assert!(global_text.contains("# observed_initial_tail_positions=30\n"));
    assert!(global_text.contains("# observed_final_tail_positions=30\n"));
    let global_lines = data_lines(&global_text);
    assert_eq!(
        global_lines[0],
        "coverage\tobserved_positions\tinitial_fitted_positions\tunderlying_fitted_positions"
    );
    assert_eq!(global_lines.len(), 22);
    for (coverage, line) in global_lines[1..].iter().enumerate() {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields[0].parse::<usize>()?, coverage);
        let observed = fields[1].parse::<u64>()?;
        let expected_observed = match coverage {
            0 => 880,
            1 => 90,
            20 => 30,
            _ => 0,
        };
        assert_eq!(observed, expected_observed);
        ensure!(fields[2].parse::<f64>()?.is_finite());
        ensure!(fields[3].parse::<f64>()?.is_finite());
    }

    let keep_text = fs::read_to_string(&result.output_keep_weights)?;
    assert!(keep_text.contains("# omitted_keep_weight=1.0\n"));
    assert!(keep_text.contains("# fragment_overlap_rule=minimum_keep_weight\n"));
    let keep_lines = data_lines(&keep_text);
    assert_eq!(keep_lines[0], "chromosome\tstart\tend\tkeep_weight");
    let keep_fields = keep_lines[1].split('\t').collect::<Vec<_>>();
    assert_eq!(&keep_fields[..3], &["chr1", "100", "130"]);
    assert!((keep_fields[3].parse::<f64>()? - 0.1).abs() < 1e-12);
    assert_eq!(keep_lines.len(), 2);

    let exact_text = fs::read_to_string(&result.output_exact_blacklist)?;
    assert_eq!(data_lines(&exact_text), ["chr1\t100\t130"]);
    let flanked_text = fs::read_to_string(&result.output_flanked_blacklist)?;
    assert_eq!(data_lines(&flanked_text), ["chr1\t50\t180"]);

    let plot = fs::read(&result.output_plot)?;
    assert!(plot.len() > 8);
    assert_eq!(&plot[..8], b"\x89PNG\r\n\x1a\n");

    Ok(())
}

#[test]
fn boundaries_blacklist_and_short_contig_preserve_calls_and_record_fallback() -> Result<()> {
    let bam = boundary_blacklist_and_fallback_fixture()?;
    let output_dir = TempDir::new()?;
    let blacklist_path = output_dir.path().join("input_blacklist.bed");
    fs::write(&blacklist_path, "chr1\t500000\t500001\nchr2\t0\t1000\n")?;
    let mut config = OutliersConfig::new(
        IOCArgs {
            bam: bam.bam_path().to_path_buf(),
            output_dir: output_dir.path().to_path_buf(),
            n_threads: 2,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
            chromosomes_file: None,
        },
    );
    config.set_output_prefix("boundaries");
    config.set_tile_size(1_000_000);
    config.set_stride(250_000);
    config.set_bin_size(1_000_000);
    config.set_min_context_fragments(1);
    config.set_tail_probability(Some(0.000_001));
    config.set_target(OutlierTarget::Threshold);
    config.set_blacklist_flank(Some(0));
    config.set_blacklist(Some(vec![blacklist_path]));
    config.set_min_mapq(0);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 30;
        fragment_lengths.max_fragment_length = 30;
    }

    let result = run_outliers(&config, RunOptions::new_quiet())?;

    let histogram_text = read_zst_to_string(&result.output_histograms)?;
    let histogram_lines = data_lines(&histogram_text);
    let eligible_positions = histogram_lines[1..]
        .iter()
        .map(|line| {
            line.split('\t')
                .nth(5)
                .context("histogram observed-position column")?
                .parse::<u64>()
                .context("integer observed positions")
        })
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .sum::<u64>();
    // The two selected contigs contain 1,001,030 positions. One chr1 base and all 1,000 chr2
    // bases are blacklisted.
    assert_eq!(eligible_positions, 1_000_029);
    let mut global_counts = BTreeMap::<u32, u64>::new();
    for line in &histogram_lines[1..] {
        let fields = line.split('\t').collect::<Vec<_>>();
        *global_counts.entry(fields[4].parse()?).or_default() += fields[5].parse::<u64>()?;
    }
    assert_eq!(
        global_counts,
        [(0, 999_850), (1, 90), (20, 89)].into_iter().collect()
    );
    assert!(
        histogram_lines
            .iter()
            .any(|line| *line == "chr2\t0\t1000\t0\t0\t0")
    );

    let models_text = fs::read_to_string(&result.output_models)?;
    let model_lines = data_lines(&models_text);
    let chr2_fields = model_lines
        .iter()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|fields| fields[0] == "chr2")
        .context("chr2 model row")?;
    assert_eq!(chr2_fields[5], "global_fallback");
    assert_eq!(&chr2_fields[3..5], &["0", "1000"]);
    assert_eq!(chr2_fields[6], "0");
    assert_eq!(chr2_fields[7], "0");
    assert_eq!(chr2_fields[16], "NaN");
    assert_eq!(chr2_fields[17], "NaN");
    assert_eq!(chr2_fields[19], "NaN");

    let exact_text = fs::read_to_string(&result.output_exact_blacklist)?;
    assert_eq!(
        data_lines(&exact_text),
        [
            "chr1\t249990\t250020",
            "chr1\t499990\t500000",
            "chr1\t500001\t500020",
            "chr1\t999990\t1000020",
        ]
    );
    // The first and last calls prove that core and tile boundaries do not split contiguous calls.
    // The middle two prove that the masked base is neither called nor bridged.
    let keep_text = fs::read_to_string(&result.output_keep_weights)?;
    let keep_lines = data_lines(&keep_text);
    assert_eq!(keep_lines.len(), 5);
    assert!(keep_lines[1..].iter().all(|line| {
        let weight = line
            .split('\t')
            .nth(3)
            .expect("keep-weight column")
            .parse::<f64>()
            .expect("numeric keep weight");
        (weight - 0.1).abs() < 1e-12
    }));

    Ok(())
}

#[test]
fn local_mean_and_zero_targets_write_their_regional_mass_weights() -> Result<()> {
    let bam = paired_outlier_fixture()?;

    let local_output_dir = TempDir::new()?;
    let mut local_config = base_config(bam.bam_path(), local_output_dir.path());
    local_config.set_target(OutlierTarget::LocalMean);
    let local_result = run_outliers(&local_config, RunOptions::new_quiet())?;
    let local_model_text = fs::read_to_string(&local_result.output_models)?;
    let local_model_fields = data_lines(&local_model_text)[1]
        .split('\t')
        .collect::<Vec<_>>();
    let underlying_mean = local_model_fields[14].parse::<f64>()?;
    let local_weight_text = fs::read_to_string(&local_result.output_keep_weights)?;
    let local_weight = data_lines(&local_weight_text)[1]
        .split('\t')
        .nth(3)
        .context("local-mean keep weight")?
        .parse::<f64>()?;
    // Every called position has coverage 20 and the same underlying fitted mean, so the regional
    // mass ratio reduces exactly to underlying_mean / 20.
    assert!((local_weight - underlying_mean / 20.0).abs() < 1e-12);
    assert!(local_weight > 0.004 && local_weight < 0.005);

    let zero_output_dir = TempDir::new()?;
    let mut zero_config = base_config(bam.bam_path(), zero_output_dir.path());
    zero_config.set_target(OutlierTarget::Zero);
    let zero_result = run_outliers(&zero_config, RunOptions::new_quiet())?;
    let zero_weight_text = fs::read_to_string(&zero_result.output_keep_weights)?;
    let zero_weight = data_lines(&zero_weight_text)[1]
        .split('\t')
        .nth(3)
        .context("zero keep weight")?
        .parse::<f64>()?;
    assert_eq!(zero_weight, 0.0);

    Ok(())
}
