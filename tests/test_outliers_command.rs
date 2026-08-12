#![cfg(feature = "cmd_outliers")]

use std::{collections::BTreeMap, fs, path::Path};

use anyhow::{Context, Result, ensure};
use cfdnalab::{
    RunOptions,
    run_like_cli::{
        common::{ChromosomeArgs, IOCArgs, OutlierWeightsArgs},
        fcoverage::{FCoverageConfig, run_fcoverage},
        outliers::{OutlierTarget, OutliersConfig, run_outliers},
    },
    testing::{PairedFragmentSpec, TempBam, TempBamBuilder, read_zst_to_string},
};
use tempfile::TempDir;

const ROUNDTRIP_CHROMOSOME_LENGTH: usize = 10_000;
const ROUNDTRIP_BACKGROUND_START: usize = 500;
const ROUNDTRIP_BACKGROUND_END: usize = 9_500;
const ROUNDTRIP_OUTLIER_START: usize = 4_400;
const ROUNDTRIP_OUTLIER_END: usize = 4_550;
const ROUNDTRIP_FRAGMENT_LENGTH: i64 = 150;
const ROUNDTRIP_BACKGROUND_COVERAGE: f64 = 3.0;
const ROUNDTRIP_RAW_OUTLIER_COVERAGE: f64 = 33.0;

const COMPLEX_CHROMOSOME_LENGTH: usize = 60_000;
const COMPLEX_OUTLIER_START: usize = 30_000;
const COMPLEX_OUTLIER_END: usize = 30_150;
const COMPLEX_MIN_FRAGMENT_LENGTH: u32 = 120;
const COMPLEX_MAX_FRAGMENT_LENGTH: u32 = 200;

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

fn normal_coverage_with_extreme_pileup_fixture() -> Result<TempBam> {
    let mut builder = TempBamBuilder::new()
        .name("outlier_weight_roundtrip")
        .contig("chr1", ROUNDTRIP_CHROMOSOME_LENGTH as u32)
        .use_record_indexed_read_names();

    // Three independently generated fragment layers tile [500, 9500), producing exact background
    // coverage 3. Their different read lengths and sequences exercise real paired-read assembly
    // while preserving the same 150 bp `pos` to `reference_end` fragment spans.
    let background_layers = [(40, b'A', b'T'), (50, b'C', b'G'), (60, b'G', b'C')];
    for (read_length, forward_base, reverse_base) in background_layers {
        for start in (ROUNDTRIP_BACKGROUND_START..ROUNDTRIP_BACKGROUND_END)
            .step_by(ROUNDTRIP_FRAGMENT_LENGTH as usize)
        {
            builder = builder.paired_fragment(
                PairedFragmentSpec::new(0, start as i64, ROUNDTRIP_FRAGMENT_LENGTH, read_length)
                    .bases(forward_base, reverse_base),
            );
        }
    }

    // Thirty additional distinct molecules share one complete background fragment span. Together
    // with the three ordinary fragments, they create coverage 33 only on [4400, 4550).
    for pileup_index in 0..30 {
        let read_length = 45 + (pileup_index % 3) * 5;
        builder = builder.paired_fragment(
            PairedFragmentSpec::new(
                0,
                ROUNDTRIP_OUTLIER_START as i64,
                ROUNDTRIP_FRAGMENT_LENGTH,
                read_length,
            )
            .base_quality(30 + (pileup_index % 10) as u8),
        );
    }

    // These plausible alignments fail the shared MAPQ filter and must not affect either command.
    for start in [650, 2_000, 4_400, 6_800, 8_900] {
        builder = builder.paired_fragment(
            PairedFragmentSpec::new(0, start, ROUNDTRIP_FRAGMENT_LENGTH, 50).mapq(5),
        );
    }

    builder.build()
}

fn complex_poisson_like_coverage_with_extreme_pileup_fixture() -> Result<TempBam> {
    let mut builder = TempBamBuilder::new()
        .name("complex_outlier_weight_roundtrip")
        .contig("chr1", COMPLEX_CHROMOSOME_LENGTH as u32)
        .use_record_indexed_read_names();

    // Independent uniform fragment starts with a mean fragment length near 160 bp produce mean
    // background coverage near 3.2. A fixed local generator keeps the BAM deterministic while
    // avoiding a periodic layer structure. Fragment and read lengths, sequence, and base quality
    // all vary across the background molecules.
    let mut generator_state = 0x4d59_5df4_d0f3_3173_u64;
    for fragment_index in 0..1_200_u64 {
        generator_state = generator_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let fragment_length = COMPLEX_MIN_FRAGMENT_LENGTH as u64
            + ((generator_state >> 32)
                % u64::from(COMPLEX_MAX_FRAGMENT_LENGTH - COMPLEX_MIN_FRAGMENT_LENGTH + 1));

        generator_state = generator_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let maximum_start = COMPLEX_CHROMOSOME_LENGTH as u64 - fragment_length;
        let start = (generator_state >> 32) % (maximum_start + 1);

        generator_state = generator_state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let read_length = 45 + ((generator_state >> 32) % 36);
        let (forward_base, reverse_base) = match fragment_index % 4 {
            0 => (b'A', b'T'),
            1 => (b'C', b'G'),
            2 => (b'G', b'C'),
            _ => (b'T', b'A'),
        };

        builder = builder.paired_fragment(
            PairedFragmentSpec::new(0, start as i64, fragment_length as i64, read_length as i64)
                .base_quality(30 + (fragment_index % 11) as u8)
                .bases(forward_base, reverse_base),
        );
    }

    // These ordinary fragments guarantee that fragment-level weighting produces shoulders outside
    // the called interval. They overlap an outlier boundary but extend into otherwise normal
    // coverage on the left or right.
    builder = builder
        .paired_fragment(PairedFragmentSpec::new(0, 29_900, 150, 55))
        .paired_fragment(PairedFragmentSpec::new(0, 30_100, 150, 65));

    // Forty distinct fragments add an extreme population across [30000, 30150). Background
    // fragments remain present underneath it, so the observed outlier coverage varies rather than
    // forming an isolated constant-coverage block.
    for pileup_index in 0..40 {
        builder = builder.paired_fragment(
            PairedFragmentSpec::new(0, COMPLEX_OUTLIER_START as i64, 150, 45 + pileup_index % 36)
                .base_quality(30 + (pileup_index % 10) as u8),
        );
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

fn roundtrip_fcoverage_config(
    bam: &Path,
    output_dir: &Path,
    output_prefix: &str,
) -> FCoverageConfig {
    let mut config = FCoverageConfig::new(
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
    config.set_output_prefix(output_prefix);
    config.set_tile_size(1_000_000);
    config.set_decimals(6);
    config.set_min_mapq(30);
    config.set_require_proper_pair(false);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = ROUNDTRIP_FRAGMENT_LENGTH as u32;
        fragment_lengths.max_fragment_length = ROUNDTRIP_FRAGMENT_LENGTH as u32;
    }
    config
}

fn complex_roundtrip_fcoverage_config(
    bam: &Path,
    output_dir: &Path,
    output_prefix: &str,
) -> FCoverageConfig {
    let mut config = FCoverageConfig::new(
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
    config.set_output_prefix(output_prefix);
    config.set_tile_size(1_000_000);
    config.set_decimals(6);
    config.set_min_mapq(30);
    config.set_require_proper_pair(false);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = COMPLEX_MIN_FRAGMENT_LENGTH;
        fragment_lengths.max_fragment_length = COMPLEX_MAX_FRAGMENT_LENGTH;
    }
    config
}

fn dense_bedgraph(text: &str) -> Vec<f64> {
    dense_bedgraph_with_length(text, ROUNDTRIP_CHROMOSOME_LENGTH)
}

fn dense_bedgraph_with_length(text: &str, chromosome_length: usize) -> Vec<f64> {
    let mut coverage = vec![0.0; chromosome_length];
    for line in text.lines().filter(|line| !line.is_empty()) {
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields.len(), 4, "unexpected bedGraph row: {line}");
        assert_eq!(fields[0], "chr1");
        let start = fields[1].parse::<usize>().expect("bedGraph start");
        let end = fields[2].parse::<usize>().expect("bedGraph end");
        let value = fields[3].parse::<f64>().expect("bedGraph value");
        coverage[start..end].fill(value);
    }
    coverage
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
fn detected_outlier_weights_flatten_an_extreme_pileup_in_fcoverage() -> Result<()> {
    // Arrange
    // Accepted fragments create the following exact raw coverage:
    //
    // - [0, 500) and [9500, 10000): coverage 0 over 1,000 positions in total
    // - [500, 4400) and [4550, 9500): ordinary coverage 3 over 8,850 positions
    // - [4400, 4550): coverage 33 over 150 positions
    //
    // The five MAPQ-5 decoys are excluded by both commands. The complete accepted-fragment count
    // is therefore 3 layers * 60 fragments + 30 pileup fragments = 210.
    let bam = normal_coverage_with_extreme_pileup_fixture()?;
    let output_dir = TempDir::new()?;

    // Act and assert the unweighted input signal before involving the detector.
    let raw_config = roundtrip_fcoverage_config(bam.bam_path(), output_dir.path(), "roundtrip_raw");
    let raw_result = run_fcoverage(&raw_config, RunOptions::new_quiet())?;
    assert_eq!(raw_result.counters.base.counted_fragments, 210);
    let raw_text = read_zst_to_string(&raw_result.final_out_path)?;
    assert_eq!(
        raw_text,
        "chr1\t500\t4400\t3\n\
         chr1\t4400\t4550\t33\n\
         chr1\t4550\t9500\t3\n"
    );
    let raw_coverage = dense_bedgraph(&raw_text);

    // Detect the pileup programmatically. Every core uses the full 10 kb chromosome as its local
    // context. The automatic tail probability is 1 / 10,000 eligible positions. Local-mean
    // targeting should estimate the ordinary coverage near 2.7 and assign that mass to the call.
    let mut outlier_config = OutliersConfig::new(
        IOCArgs {
            bam: bam.bam_path().to_path_buf(),
            output_dir: output_dir.path().to_path_buf(),
            n_threads: 2,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    );
    outlier_config.set_output_prefix("roundtrip_detected");
    outlier_config.set_tile_size(1_000_000);
    outlier_config.set_stride(1_000);
    outlier_config.set_bin_size(10_000);
    outlier_config.set_min_context_fragments(200);
    outlier_config.set_target(OutlierTarget::LocalMean);
    outlier_config.set_blacklist_flank(Some(0));
    outlier_config.set_min_mapq(30);
    outlier_config.set_require_proper_pair(false);
    {
        let fragment_lengths = outlier_config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = ROUNDTRIP_FRAGMENT_LENGTH as u32;
        fragment_lengths.max_fragment_length = ROUNDTRIP_FRAGMENT_LENGTH as u32;
    }

    let outlier_result = run_outliers(&outlier_config, RunOptions::new_quiet())?;
    assert_eq!(outlier_result.counters.base.counted_fragments, 210);

    // The persisted histograms must reproduce the hand-derived raw signal, not merely contain an
    // unspecified high-coverage bin.
    let histogram_text = read_zst_to_string(&outlier_result.output_histograms)?;
    let mut global_coverage_counts = BTreeMap::<u32, u64>::new();
    for line in data_lines(&histogram_text).into_iter().skip(1) {
        let fields = line.split('\t').collect::<Vec<_>>();
        *global_coverage_counts
            .entry(fields[4].parse::<u32>()?)
            .or_default() += fields[5].parse::<u64>()?;
    }
    assert_eq!(
        global_coverage_counts,
        [(0, 1_000), (3, 8_850), (33, 150)].into_iter().collect()
    );

    // Inspect the actual local model assigned to the core containing the pileup. Its final
    // threshold must separate ordinary coverage 3 from extreme coverage 33. The fitted ZIP mean
    // is close to 2.7 because 90% of retained positions have positive coverage with conditional
    // positive mean 3.
    let model_text = fs::read_to_string(&outlier_result.output_models)?;
    let pileup_model = data_lines(&model_text)
        .into_iter()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|fields| fields[1] == "4000")
        .context("model core containing the roundtrip pileup")?;
    assert_eq!(
        &pileup_model[..8],
        &[
            "chr1", "4000", "5000", "0", "10000", "local", "10000", "210"
        ]
    );
    let underlying_mean = pileup_model[14].parse::<f64>()?;
    let final_threshold = pileup_model[15].parse::<u32>()?;
    assert!(
        (2.69..=2.71).contains(&underlying_mean),
        "unexpected fitted ordinary coverage {underlying_mean}"
    );
    assert!(
        final_threshold > ROUNDTRIP_BACKGROUND_COVERAGE as u32
            && final_threshold < ROUNDTRIP_RAW_OUTLIER_COVERAGE as u32,
        "threshold {final_threshold} should separate coverage 3 from coverage 33"
    );

    let exact_text = fs::read_to_string(&outlier_result.output_exact_blacklist)?;
    assert_eq!(data_lines(&exact_text), ["chr1\t4400\t4550"]);
    let keep_weight_text = fs::read_to_string(&outlier_result.output_keep_weights)?;
    let keep_weight_lines = data_lines(&keep_weight_text);
    assert_eq!(keep_weight_lines.len(), 2);
    let keep_weight_fields = keep_weight_lines[1].split('\t').collect::<Vec<_>>();
    assert_eq!(&keep_weight_fields[..3], &["chr1", "4400", "4550"]);
    let keep_weight = keep_weight_fields[3].parse::<f64>()?;
    let expected_keep_weight = underlying_mean / ROUNDTRIP_RAW_OUTLIER_COVERAGE;
    assert!(
        (keep_weight - expected_keep_weight).abs() <= 1e-12,
        "expected keep weight {expected_keep_weight}, got {keep_weight}"
    );

    // Pass the detector's file directly to fcoverage. All 33 fragments covering the called region
    // receive the same scalar weight, so its raw coverage 33 becomes the fitted local mean. The
    // surrounding ordinary coverage and uncovered chromosome ends must remain unchanged.
    let mut weighted_config =
        roundtrip_fcoverage_config(bam.bam_path(), output_dir.path(), "roundtrip_weighted");
    weighted_config.set_outlier_weights(OutlierWeightsArgs {
        outlier_weights: Some(outlier_result.output_keep_weights.clone()),
    });
    let weighted_result = run_fcoverage(&weighted_config, RunOptions::new_quiet())?;
    assert_eq!(weighted_result.counters.base.counted_fragments, 210);
    let weighted_text = read_zst_to_string(&weighted_result.final_out_path)?;
    let weighted_coverage = dense_bedgraph(&weighted_text);

    let expected_flattened_coverage = ROUNDTRIP_RAW_OUTLIER_COVERAGE * expected_keep_weight;
    for position in 0..ROUNDTRIP_CHROMOSOME_LENGTH {
        let expected = if (ROUNDTRIP_OUTLIER_START..ROUNDTRIP_OUTLIER_END).contains(&position) {
            expected_flattened_coverage
        } else {
            raw_coverage[position]
        };
        assert!(
            (weighted_coverage[position] - expected).abs() <= 2e-6,
            "position {position}: expected {expected}, got {}",
            weighted_coverage[position]
        );
    }
    assert_eq!(
        raw_coverage[ROUNDTRIP_OUTLIER_START],
        ROUNDTRIP_RAW_OUTLIER_COVERAGE
    );
    assert!((weighted_coverage[ROUNDTRIP_OUTLIER_START] - underlying_mean).abs() <= 2e-6);
    assert!(
        weighted_coverage.iter().copied().fold(0.0_f64, f64::max) <= ROUNDTRIP_BACKGROUND_COVERAGE
    );
    let weighted_outlier_mass = weighted_coverage[ROUNDTRIP_OUTLIER_START..ROUNDTRIP_OUTLIER_END]
        .iter()
        .sum::<f64>();
    let expected_target_mass =
        underlying_mean * (ROUNDTRIP_OUTLIER_END - ROUNDTRIP_OUTLIER_START) as f64;
    assert!(
        (weighted_outlier_mass - expected_target_mass).abs() <= 0.000_3,
        "expected regional target mass {expected_target_mass}, got {weighted_outlier_mass}"
    );

    Ok(())
}

#[test]
fn detected_outlier_weights_flatten_an_extreme_with_complex_poisson_like_background() -> Result<()>
{
    // Arrange
    // The background consists of 1,202 irregular paired fragments with independently distributed
    // starts and lengths from 120 through 200 bp. Forty additional fragments create the extreme
    // population. The raw and weighted tracks are generated through the public fcoverage runner,
    // while the weights between them come directly from the public outliers runner.
    let bam = complex_poisson_like_coverage_with_extreme_pileup_fixture()?;
    let output_dir = TempDir::new()?;

    // Act
    let raw_config = complex_roundtrip_fcoverage_config(
        bam.bam_path(),
        output_dir.path(),
        "complex_roundtrip_raw",
    );
    let raw_result = run_fcoverage(&raw_config, RunOptions::new_quiet())?;
    assert_eq!(raw_result.counters.base.counted_fragments, 1_242);
    let raw_text = read_zst_to_string(&raw_result.final_out_path)?;
    let raw_coverage = dense_bedgraph_with_length(&raw_text, COMPLEX_CHROMOSOME_LENGTH);

    // Assert that this is a variable Poisson-like background rather than another constant layer.
    // Excluding the known pileup span, the expected mean is near 3.2 from total fragment span over
    // chromosome length. Independently placed fragments should also have variance near their mean.
    let background_values = raw_coverage
        .iter()
        .enumerate()
        .filter(|(position, _)| !(COMPLEX_OUTLIER_START..COMPLEX_OUTLIER_END).contains(position))
        .map(|(_, &coverage)| coverage)
        .collect::<Vec<_>>();
    let background_mean = background_values.iter().sum::<f64>() / background_values.len() as f64;
    let background_variance = background_values
        .iter()
        .map(|coverage| (coverage - background_mean).powi(2))
        .sum::<f64>()
        / background_values.len() as f64;
    let distinct_background_coverages = background_values
        .iter()
        .map(|coverage| *coverage as u32)
        .collect::<std::collections::BTreeSet<_>>();
    assert!(
        (2.8..=3.6).contains(&background_mean),
        "unexpected background mean {background_mean}"
    );
    assert!(
        (0.65..=1.35).contains(&(background_variance / background_mean)),
        "background variance-to-mean ratio was {}",
        background_variance / background_mean
    );
    assert!(
        distinct_background_coverages.len() >= 8,
        "expected a broad background coverage distribution, got {distinct_background_coverages:?}"
    );

    let mut outlier_config = OutliersConfig::new(
        IOCArgs {
            bam: bam.bam_path().to_path_buf(),
            output_dir: output_dir.path().to_path_buf(),
            n_threads: 2,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    );
    outlier_config.set_output_prefix("complex_roundtrip_detected");
    outlier_config.set_tile_size(1_000_000);
    outlier_config.set_stride(5_000);
    outlier_config.set_bin_size(COMPLEX_CHROMOSOME_LENGTH as u32);
    outlier_config.set_min_context_fragments(1_000);
    outlier_config.set_target(OutlierTarget::LocalMean);
    outlier_config.set_blacklist_flank(Some(0));
    outlier_config.set_min_mapq(30);
    outlier_config.set_require_proper_pair(false);
    {
        let fragment_lengths = outlier_config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = COMPLEX_MIN_FRAGMENT_LENGTH;
        fragment_lengths.max_fragment_length = COMPLEX_MAX_FRAGMENT_LENGTH;
    }

    let outlier_result = run_outliers(&outlier_config, RunOptions::new_quiet())?;
    assert_eq!(outlier_result.counters.base.counted_fragments, 1_242);

    // Every core fits the full chromosome. The threshold used by the pileup core must sit above
    // every ordinary background value and below every value in the added extreme population.
    let model_text = fs::read_to_string(&outlier_result.output_models)?;
    let pileup_model = data_lines(&model_text)
        .into_iter()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .find(|fields| fields[1] == "30000")
        .context("model core containing the complex roundtrip pileup")?;
    assert_eq!(
        &pileup_model[3..8],
        &["0", "60000", "local", "60000", "1242"]
    );
    let underlying_mean = pileup_model[14].parse::<f64>()?;
    let final_threshold = pileup_model[15].parse::<u32>()?;
    let maximum_background_coverage =
        background_values.iter().copied().fold(0.0_f64, f64::max) as u32;
    let minimum_outlier_coverage = raw_coverage[COMPLEX_OUTLIER_START..COMPLEX_OUTLIER_END]
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min) as u32;
    assert!(
        (2.8..=3.6).contains(&underlying_mean),
        "unexpected fitted background mean {underlying_mean}"
    );
    assert!(
        final_threshold > maximum_background_coverage
            && final_threshold <= minimum_outlier_coverage,
        "threshold {final_threshold} did not separate background maximum {maximum_background_coverage} from outlier minimum {minimum_outlier_coverage}"
    );

    let exact_text = fs::read_to_string(&outlier_result.output_exact_blacklist)?;
    assert_eq!(
        data_lines(&exact_text),
        ["chr1\t30000\t30150"],
        "only the added extreme population should be called"
    );
    let keep_weight_text = fs::read_to_string(&outlier_result.output_keep_weights)?;
    let keep_weight_lines = data_lines(&keep_weight_text);
    assert_eq!(keep_weight_lines.len(), 2);
    let keep_weight_fields = keep_weight_lines[1].split('\t').collect::<Vec<_>>();
    assert_eq!(&keep_weight_fields[..3], &["chr1", "30000", "30150"]);
    let keep_weight = keep_weight_fields[3].parse::<f64>()?;
    let raw_outlier_mass = raw_coverage[COMPLEX_OUTLIER_START..COMPLEX_OUTLIER_END]
        .iter()
        .sum::<f64>();
    let expected_target_mass =
        underlying_mean * (COMPLEX_OUTLIER_END - COMPLEX_OUTLIER_START) as f64;
    let expected_keep_weight = expected_target_mass / raw_outlier_mass;
    assert!(
        (keep_weight - expected_keep_weight).abs() <= 1e-12,
        "expected keep weight {expected_keep_weight}, got {keep_weight}"
    );

    let mut weighted_config = complex_roundtrip_fcoverage_config(
        bam.bam_path(),
        output_dir.path(),
        "complex_roundtrip_weighted",
    );
    weighted_config.set_outlier_weights(OutlierWeightsArgs {
        outlier_weights: Some(outlier_result.output_keep_weights.clone()),
    });
    let weighted_result = run_fcoverage(&weighted_config, RunOptions::new_quiet())?;
    assert_eq!(weighted_result.counters.base.counted_fragments, 1_242);
    let weighted_text = read_zst_to_string(&weighted_result.final_out_path)?;
    let weighted_coverage = dense_bedgraph_with_length(&weighted_text, COMPLEX_CHROMOSOME_LENGTH);

    // Every fragment covering a called position necessarily overlaps the call, so the complete
    // called interval is multiplied by the same regional keep weight. Regional mass must therefore
    // reach the local-mean target despite the varying original coverage within the interval.
    for position in COMPLEX_OUTLIER_START..COMPLEX_OUTLIER_END {
        let expected = raw_coverage[position] * keep_weight;
        assert!(
            (weighted_coverage[position] - expected).abs() <= 2e-6,
            "position {position}: expected {expected}, got {}",
            weighted_coverage[position]
        );
    }
    let weighted_outlier_mass = weighted_coverage[COMPLEX_OUTLIER_START..COMPLEX_OUTLIER_END]
        .iter()
        .sum::<f64>();
    assert!(
        (weighted_outlier_mass - expected_target_mass).abs() <= 0.000_3,
        "expected regional target mass {expected_target_mass}, got {weighted_outlier_mass}"
    );

    // Fragment weighting also reduces the portions of overlapping ordinary fragments that extend
    // outside the called interval. Positions farther than the maximum accepted fragment length
    // cannot belong to a fragment touching the call and must remain exactly unchanged.
    assert!(
        (COMPLEX_OUTLIER_START - COMPLEX_MAX_FRAGMENT_LENGTH as usize..COMPLEX_OUTLIER_START)
            .any(|position| weighted_coverage[position] < raw_coverage[position])
    );
    assert!(
        (COMPLEX_OUTLIER_END..COMPLEX_OUTLIER_END + COMPLEX_MAX_FRAGMENT_LENGTH as usize)
            .any(|position| weighted_coverage[position] < raw_coverage[position])
    );
    for position in 0..COMPLEX_CHROMOSOME_LENGTH {
        assert!(
            weighted_coverage[position] <= raw_coverage[position] + 1e-6,
            "weighting increased coverage at position {position}"
        );
        if position < COMPLEX_OUTLIER_START - COMPLEX_MAX_FRAGMENT_LENGTH as usize
            || position >= COMPLEX_OUTLIER_END + COMPLEX_MAX_FRAGMENT_LENGTH as usize
        {
            assert!(
                (weighted_coverage[position] - raw_coverage[position]).abs() <= 1e-6,
                "distant position {position} changed from {} to {}",
                raw_coverage[position],
                weighted_coverage[position]
            );
        }
    }

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
