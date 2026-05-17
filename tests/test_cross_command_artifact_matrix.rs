#![cfg(all(
    feature = "cmd_bam_to_bam",
    feature = "cmd_coverage_weights",
    feature = "cmd_fcoverage",
    feature = "cmd_gc_bias",
    feature = "cmd_lengths",
    feature = "cmd_midpoints",
    feature = "cmd_ref_gc_bias"
))]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::{
    bam_to_bam::{bam_to_bam::run_inner as run_bam_to_bam, config::BamToBamConfig},
    cli_common::{
        ApplyGCArgFileOnly, ApplyGCArgs, AssignToWindowArgs, ChromosomeArgs,
        DistributionWindowsArgs, IOCArgs, ScaleGenomeArgs,
    },
    coverage_weights::{
        config::CoverageWeightsConfig, coverage_weights::run as run_coverage_weights,
    },
    fcoverage::{config::FCoverageConfig, fcoverage::run as run_fcoverage},
    lengths::{config::LengthsConfig, lengths::run as run_lengths},
    midpoints::{
        config::MidpointsConfig, midpoints::run as run_midpoints, smoothing::MidpointSmoothing,
    },
};
use cfdnalab::shared::{indel_mode::IndelMode, io::dot_join};
use fixtures::{
    BamFixture, TwoBitFixture, bam_from_specs, build_real_neutral_gc_package,
    build_real_neutral_gc_package_for_range, paired_fragment, read_length_counts_tsv,
    read_midpoint_zarr_counts, read_zst_to_string, simple_inward_bam, simple_reference_twobit,
    write_bed,
};
use ndarray::{Array2, Array3};
use rust_htslib::bam::{Read, Reader, record::Aux};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use tempfile::TempDir;

const F64_TOL: f64 = 1e-6;
const F32_TOL: f32 = 1e-6;
const EXPECTED_FRAGMENT_AVERAGE: f64 = 2146.0_f64 / 2745.0_f64;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|chr| chr.to_string()).collect()),
        chromosomes_file: None,
    }
}

fn make_real_scaling_config(output_dir: &Path, bam_path: &Path) -> CoverageWeightsConfig {
    let mut cfg = CoverageWeightsConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_output_prefix("coverage".to_string());
    cfg.set_bin_size(40);
    cfg.set_stride(20);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    {
        let fragment_lengths = cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 200;
    }
    cfg
}

#[derive(Debug)]
struct SharedRealArtifacts {
    tempdir: TempDir,
    consumer_bam: PathBuf,
    reference_path: PathBuf,
    scaling_path: PathBuf,
    gc_path: PathBuf,
    bed_path: PathBuf,
}

#[derive(Debug)]
struct SharedRealArtifactCache {
    _producer_bam: BamFixture,
    _consumer_bam_fixture: BamFixture,
    _reference_fixture: TwoBitFixture,
    _fixture_root: TempDir,
    consumer_bam: PathBuf,
    reference_path: PathBuf,
    scaling_path: PathBuf,
    gc_path: PathBuf,
    bed_path: PathBuf,
}

fn shared_real_artifact_cache() -> Result<&'static SharedRealArtifactCache> {
    static CACHE: OnceLock<std::result::Result<SharedRealArtifactCache, String>> = OnceLock::new();
    let cache_result = CACHE.get_or_init(|| {
        (|| -> Result<SharedRealArtifactCache> {
            let scaling_producer_bam = simple_inward_bam()?;
            let consumer_bam = bam_from_specs(
                vec![("chr1".to_string(), 200)],
                vec![paired_fragment(20, 61, 20)],
                Vec::new(),
                "shared_real_artifacts_consumer",
            )?;
            let reference = simple_reference_twobit()?;
            let fixture_root = TempDir::new()?;
            let weights_out_dir = fixture_root.path().join("coverage_weights");
            std::fs::create_dir_all(&weights_out_dir)?;
            let scaling_gc_path = build_real_neutral_gc_package_for_range(
                &scaling_producer_bam.bam,
                &reference.path,
                fixture_root.path(),
                10,
                200,
            )?;

            // Shared fixture reasoning:
            // - The scaling producer is the standard one-fragment fixture [20,80).
            // - Its neutral real GC package is built from the same BAM and repeated ACGT reference with
            //   the same length range that `coverage-weights` is configured to consume, 10..200.
            // - Within that broader package, the only accepted 60 bp fragment still receives GC weight
            //   1.0, so the coverage profile stays the same while the GC path is exercised honestly.
            // - `coverage-weights --gc-file <neutral package>` with bin_size=40 and stride=20 therefore
            //   writes the already hand-derived scaling profile:
            //     [0,20):   37/20
            //     [20,40):  37/45
            //     [40,60):  37/60
            //     [60,80):  37/45
            //     [80,100): 37/15
            // - The consumer is one 61 bp fragment [20,81) on the same repeated ACGT reference.
            // - The real `ref-gc-bias -> gc-bias` chain is neutral for that consumer:
            //   all GC-by-length mass lands in one shared cell, so the final GC weight is exactly 1.0.
            // - Fragment-average consumers therefore see only the scaling average:
            //     (20*(37/45) + 20*(37/60) + 20*(37/45) + 1*(37/15)) / 61 = 2146/2745.
            let mut scaling_cfg =
                make_real_scaling_config(&weights_out_dir, &scaling_producer_bam.bam);
            scaling_cfg.set_gc(ApplyGCArgs {
                gc_file: Some(scaling_gc_path),
                gc_tag: None,
                neutralize_invalid_gc: false,
            });
            scaling_cfg.set_ref_2bit(Some(reference.path.clone()));
            run_coverage_weights(&scaling_cfg)?;

            let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");
            let gc_path = build_real_neutral_gc_package(
                &consumer_bam.bam,
                &reference.path,
                fixture_root.path(),
                61,
            )?;
            let bed_path = fixture_root.path().join("windows.bed");
            write_bed(&bed_path, &[("chr1", 45, 56, "groupA")])?;

            Ok(SharedRealArtifactCache {
                _producer_bam: scaling_producer_bam,
                _consumer_bam_fixture: consumer_bam,
                _reference_fixture: reference,
                _fixture_root: fixture_root,
                consumer_bam: PathBuf::new(),
                reference_path: PathBuf::new(),
                scaling_path,
                gc_path,
                bed_path,
            })
            .map(|mut cache| {
                cache.consumer_bam = cache._consumer_bam_fixture.bam.clone();
                cache.reference_path = cache._reference_fixture.path.clone();
                cache
            })
        })()
        .map_err(|err| format!("{err:#}"))
    });
    match cache_result {
        Ok(cache) => Ok(cache),
        Err(err) => Err(anyhow::anyhow!(
            "failed to build shared real artifact cache: {err}"
        )),
    }
}

fn build_shared_real_artifacts() -> Result<SharedRealArtifacts> {
    let cache = shared_real_artifact_cache()?;
    let tempdir = TempDir::new()?;

    Ok(SharedRealArtifacts {
        tempdir,
        consumer_bam: cache.consumer_bam.clone(),
        reference_path: cache.reference_path.clone(),
        scaling_path: cache.scaling_path.clone(),
        gc_path: cache.gc_path.clone(),
        bed_path: cache.bed_path.clone(),
    })
}

fn read_bam_tags(path: &Path) -> Result<Vec<(Option<f32>, Option<f32>, Option<u32>)>> {
    let mut reader = Reader::from_path(path)?;
    let mut tags = Vec::new();
    for record_result in reader.records() {
        let record = record_result?;
        let gc = match record.aux(b"GC") {
            Ok(Aux::Float(value)) => Some(value),
            _ => None,
        };
        let cov = match record.aux(b"cw") {
            Ok(Aux::Float(value)) => Some(value),
            _ => None,
        };
        let flen = match record.aux(b"fl") {
            Ok(Aux::U32(value)) => Some(value),
            _ => None,
        };
        tags.push((gc, cov, flen));
    }
    Ok(tags)
}

fn assert_close_f64(actual: f64, expected: f64, context: &str) {
    assert!(
        (actual - expected).abs() <= F64_TOL,
        "{context}: expected {expected}, got {actual}"
    );
}

fn assert_close_f32(actual: f32, expected: f32, context: &str) {
    assert!(
        (actual - expected).abs() <= F32_TOL,
        "{context}: expected {expected}, got {actual}"
    );
}

#[test]
fn bam_to_bam_consumes_shared_real_artifacts_with_expected_fragment_tags() -> Result<()> {
    // Arrange:
    // The shared builder gives us:
    // - one real neutral GC package with weight 1.0 for the only accepted fragment
    // - one real scaling TSV whose full-fragment average over [20,81) is 2146/2745
    //
    // `bam-to-bam` writes fragment-level metadata on both mates, so the expected output is:
    // - two records
    // - each tagged with GC=1.0, cw=2146/2745, fl=61
    let artifacts = build_shared_real_artifacts()?;
    let out_bam = artifacts.tempdir.path().join("tagged.bam");

    let mut cfg = BamToBamConfig::new(
        artifacts.consumer_bam.clone(),
        out_bam.clone(),
        base_chromosomes(&["chr1"]),
    );
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_coverage_scaling_factors(Some(artifacts.scaling_path.clone()));
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(artifacts.gc_path.clone()),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(artifacts.reference_path.clone()));
    {
        let fragment_lengths = cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 61;
        fragment_lengths.max_fragment_length = 61;
    }

    // Act
    run_bam_to_bam(&cfg)?;

    // Assert
    let tags = read_bam_tags(&out_bam)?;
    assert_eq!(tags.len(), 2);
    for (mate_index, (gc, cov, flen)) in tags.into_iter().enumerate() {
        assert_eq!(gc, Some(1.0), "mate {mate_index} GC tag");
        assert_close_f32(
            cov.expect("cw tag should be present"),
            EXPECTED_FRAGMENT_AVERAGE as f32,
            &format!("mate {mate_index} cw tag"),
        );
        assert_eq!(flen, Some(61), "mate {mate_index} fl tag");
    }

    Ok(())
}

#[test]
fn lengths_consumes_shared_real_artifacts_with_expected_weighted_count() -> Result<()> {
    // Arrange:
    // The shared real artifacts define exactly one accepted 61 bp fragment.
    // The real GC package is neutral, so `lengths` should only apply the real scaling average:
    //   2146/2745
    // Therefore the 1x1 global output matrix must contain exactly that value.
    let artifacts = build_shared_real_artifacts()?;
    let out_dir = artifacts.tempdir.path().join("lengths");

    let mut cfg = LengthsConfig::new(
        IOCArgs {
            bam: artifacts.consumer_bam.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_indel_mode(IndelMode::Ignore);
    cfg.set_windows(DistributionWindowsArgs::default());
    cfg.set_window_assignment(AssignToWindowArgs::default());
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scaling_factors(Some(artifacts.scaling_path.clone()));
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(artifacts.gc_path.clone()),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(artifacts.reference_path.clone()));
    cfg.set_per_bp_length_bins(61, 61);

    // Act
    run_lengths(&cfg)?;

    // Assert
    let counts_path = out_dir.join(dot_join(&[
        cfg.output_prefix.trim(),
        "length_counts.tsv.zst",
    ]));
    let arr: Array2<f64> = read_length_counts_tsv(&counts_path)?;
    assert_eq!(arr.dim(), (1, 1));
    assert_close_f64(
        arr[(0, 0)],
        EXPECTED_FRAGMENT_AVERAGE,
        "lengths weighted count",
    );

    Ok(())
}

#[test]
fn midpoints_consumes_shared_real_artifacts_with_expected_profile_mass() -> Result<()> {
    // Arrange:
    // The shared consumer fragment is [20,81), so its midpoint is:
    //   20 + floor(61 / 2) = 50
    // One BED window [45,56) therefore receives the mass at profile position:
    //   50 - 45 = 5
    // The real GC package is neutral, so the written midpoint mass is just the real scaling
    // average 2146/2745 at that single position.
    let artifacts = build_shared_real_artifacts()?;
    let out_dir = artifacts.tempdir.path().join("midpoints");

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: artifacts.consumer_bam.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        artifacts.bed_path.clone(),
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![61, 62]);
    cfg.set_smoothing(MidpointSmoothing::None);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    {
        let mut scale_genome = ScaleGenomeArgs::default();
        scale_genome.scaling_factors = Some(artifacts.scaling_path.clone());
        cfg.set_scale_genome(scale_genome);
    }
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(artifacts.gc_path.clone()),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(artifacts.reference_path.clone()));

    // Act
    run_midpoints(&cfg)?;

    // Assert
    let counts_path = out_dir.join("sites.midpoint_profiles.zarr");
    let arr: Array3<f32> = read_midpoint_zarr_counts(&counts_path)?;
    assert_eq!(arr.shape(), &[1, 1, 11]);
    for (position, value) in arr.slice(ndarray::s![0, 0, ..]).iter().enumerate() {
        let expected = if position == 5 {
            EXPECTED_FRAGMENT_AVERAGE as f32
        } else {
            0.0_f32
        };
        assert_close_f32(*value, expected, &format!("midpoint position {position}"));
    }
    assert_close_f32(
        arr.sum(),
        EXPECTED_FRAGMENT_AVERAGE as f32,
        "midpoints total mass",
    );

    Ok(())
}

#[test]
fn fcoverage_consumes_shared_real_artifacts_with_expected_per_base_profile() -> Result<()> {
    // Arrange:
    // `fcoverage` is intentionally different from the fragment-average consumers.
    // It applies the same real scaling artifact per covered base in place, not as one fragment
    // average. The real GC package is neutral, so the bedGraph must contain the raw real scaling
    // runs over the covered fragment [20,81):
    // - [20,40):  37/45
    // - [40,60):  37/60
    // - [60,80):  37/45
    // - [80,81):  37/15
    let artifacts = build_shared_real_artifacts()?;
    let out_dir = artifacts.tempdir.path().join("fcoverage");

    let mut cfg = FCoverageConfig::new(
        IOCArgs {
            bam: artifacts.consumer_bam.clone(),
            output_dir: out_dir.clone(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    cfg.set_output_prefix("testcov");
    cfg.set_tile_size(1_000);
    cfg.set_decimals(6);
    cfg.set_keep_zero_runs(false);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    {
        let mut scale_genome = ScaleGenomeArgs::default();
        scale_genome.scaling_factors = Some(artifacts.scaling_path.clone());
        cfg.set_scale_genome(scale_genome);
    }
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(artifacts.gc_path.clone()),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(artifacts.reference_path.clone()));
    {
        let fragment_lengths = cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 61;
        fragment_lengths.max_fragment_length = 61;
    }

    // Act
    run_fcoverage(&cfg)?;

    // Assert
    let bedgraph_path = out_dir.join("testcov.fcoverage.per_position.bedgraph.zst");
    let text = read_zst_to_string(&bedgraph_path)?;
    let lines: Vec<_> = text.lines().collect();
    let expected = [
        (20_u64, 40_u64, 37.0_f64 / 45.0_f64),
        (40_u64, 60_u64, 37.0_f64 / 60.0_f64),
        (60_u64, 80_u64, 37.0_f64 / 45.0_f64),
        (80_u64, 81_u64, 37.0_f64 / 15.0_f64),
    ];
    assert_eq!(lines.len(), expected.len());
    for (line_index, (line, (expected_start, expected_end, expected_value))) in
        lines.iter().zip(expected).enumerate()
    {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(
            fields.len(),
            4,
            "unexpected bedGraph row {line_index}: {line}"
        );
        assert_eq!(fields[0], "chr1");
        assert_eq!(fields[1].parse::<u64>()?, expected_start);
        assert_eq!(fields[2].parse::<u64>()?, expected_end);
        assert_close_f64(
            fields[3].parse::<f64>()?,
            expected_value,
            &format!("fcoverage row {line_index}"),
        );
    }

    Ok(())
}
