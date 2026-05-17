#![cfg(all(feature = "cmd_midpoints", feature = "cmd_ends"))]

mod fixtures;

use anyhow::{Context, Result};
use cfdnalab::commands::{
    cli_common::{ChromosomeArgs, DistributionWindowsArgs, IOCArgs},
    ends::{
        config::EndsConfig,
        config_structs::{AssignMotifToWindowArgs, ClipStrategy, KmerSource, WindowMotifAssigner},
        ends::run as run_ends,
    },
    midpoints::{
        config::MidpointsConfig, midpoints::run as run_midpoints, smoothing::MidpointSmoothing,
    },
};
use fixtures::{
    bam_from_specs_strict_identity, paired_fragment, read_midpoint_zarr_counts,
    read_midpoint_zarr_i32_1d, read_midpoint_zarr_u32_1d, simple_reference_twobit, write_bed,
};
use ndarray::{Array2, arr3};
use serde_json::Value;
use std::{
    env,
    path::{Path, PathBuf},
    sync::Arc,
};
use zarrs::{array::Array, filesystem::FilesystemStore};

#[test]
// This is ignored because it is not a normal correctness test. It generates a real cfDNAlab
// midpoint Zarr output for the downstream Python/R compatibility workflow, which invokes it
// explicitly from the downstream-zarr GitHub Action.
#[ignore = "generates the Zarr fixture consumed by downstream Python/R compatibility tests"]
fn generate_midpoint_zarr_fixture_with_cfdnalab() -> Result<()> {
    let output_dir = env::var_os("CFDNALAB_DOWNSTREAM_FIXTURE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("downstream_tests/tmp"));
    std::fs::create_dir_all(&output_dir)?;

    let bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 300)],
        // The midpoint set intentionally creates distinct nonzero values across
        // group, length, and position axes. With --bin-size 2 each raw midpoint
        // contributes 0.5 to its final position bin, so repeated fragments create
        // asymmetric values that catch axis swaps and broad downstream assertions.
        vec![
            // LYL1, length bin [30, 50)
            paired_fragment(26, 41, 20),
            paired_fragment(26, 41, 20),
            paired_fragment(28, 41, 20),
            // LYL1, length bin [50, 70)
            paired_fragment(95, 61, 20),
            paired_fragment(95, 61, 20),
            paired_fragment(95, 61, 20),
            paired_fragment(97, 61, 20),
            // LYL1, length bin [70, 100)
            paired_fragment(89, 81, 20),
            paired_fragment(89, 81, 20),
            paired_fragment(89, 81, 20),
            paired_fragment(89, 81, 20),
            // beta-site, length bin [30, 50)
            paired_fragment(141, 41, 20),
            paired_fragment(143, 41, 20),
            paired_fragment(143, 41, 20),
            // beta-site, length bin [50, 70)
            paired_fragment(25, 61, 20),
            paired_fragment(25, 61, 20),
            paired_fragment(25, 61, 20),
            paired_fragment(29, 61, 20),
            // beta-site, length bin [70, 100)
            paired_fragment(123, 81, 20),
            paired_fragment(127, 81, 20),
            paired_fragment(127, 81, 20),
            // gamma_long, length bin [30, 50)
            paired_fragment(55, 41, 20),
            paired_fragment(55, 41, 20),
            paired_fragment(55, 41, 20),
            paired_fragment(55, 41, 20),
            paired_fragment(55, 41, 20),
            // gamma_long, length bin [50, 70)
            paired_fragment(41, 61, 20),
            paired_fragment(47, 61, 20),
            paired_fragment(47, 61, 20),
            paired_fragment(47, 61, 20),
            // gamma_long, length bin [70, 100)
            paired_fragment(35, 81, 20),
            paired_fragment(39, 81, 20),
            paired_fragment(39, 81, 20),
        ],
        Vec::new(),
        "downstream_midpoint_fixture",
    )?;
    let intervals = output_dir.join("tiny.midpoint_intervals.bed");
    write_bed(
        &intervals,
        &[
            ("chr1", 45, 55, "LYL1"),
            ("chr1", 120, 130, "LYL1"),
            ("chr1", 50, 60, "beta-site"),
            ("chr1", 160, 170, "beta-site"),
            ("chr1", 70, 80, "gamma_long"),
            ("chr1", 200, 210, "gamma_long"),
        ],
    )?;

    let mut config = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: output_dir.clone(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        intervals,
    );
    config.set_output_prefix("tiny");
    config.set_length_bins(vec![30, 50, 70, 100]);
    config.set_bin_size(2);
    config.set_smoothing(MidpointSmoothing::None);
    config.set_tile_size(1_000_000);
    config.set_min_mapq(0);
    config.set_require_proper_pair(false);
    config.plot_groups.clear();

    run_midpoints(&config)?;

    let zarr_path = output_dir.join("tiny.midpoint_profiles.zarr");
    let group_index_path = output_dir.join("tiny.group_index.tsv");
    let settings_path = output_dir.join("tiny.midpoint_settings.json");

    assert!(zarr_path.is_dir());
    assert!(group_index_path.is_file());
    assert!(settings_path.is_file());

    let counts = read_midpoint_zarr_counts(&zarr_path)?;
    assert_eq!(counts.shape(), &[3, 3, 5]);
    assert_eq!(
        counts,
        arr3(&[
            [
                [1.0_f32, 0.5, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.5, 0.5, 0.0],
                [0.0, 0.0, 0.0, 0.0, 2.0],
            ],
            [
                [0.5, 1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 1.5, 0.0, 0.5],
                [0.0, 0.5, 0.0, 1.0, 0.0],
            ],
            [
                [0.0, 0.0, 2.5, 0.0, 0.0],
                [0.5, 0.0, 0.0, 1.5, 0.0],
                [0.0, 0.0, 0.5, 0.0, 1.0],
            ],
        ])
    );
    assert_eq!(counts.sum(), 16.5);

    assert_eq!(
        read_midpoint_zarr_i32_1d(&zarr_path, "/group")?,
        vec![0, 1, 2]
    );
    assert_eq!(
        read_midpoint_zarr_i32_1d(&zarr_path, "/length_bin")?,
        vec![0, 1, 2]
    );
    assert_eq!(
        read_midpoint_zarr_i32_1d(&zarr_path, "/position")?,
        vec![0, 1, 2, 3, 4]
    );
    assert_eq!(
        read_midpoint_zarr_i32_1d(&zarr_path, "/position_bin_start_bp")?,
        vec![0, 2, 4, 6, 8]
    );
    assert_eq!(
        read_midpoint_zarr_i32_1d(&zarr_path, "/position_bin_end_bp")?,
        vec![2, 4, 6, 8, 10]
    );
    assert_eq!(
        read_midpoint_zarr_u32_1d(&zarr_path, "/eligible_intervals")?,
        vec![2, 2, 2]
    );
    assert_eq!(
        read_midpoint_zarr_u32_1d(&zarr_path, "/length_start_bp")?,
        vec![30, 50, 70]
    );
    assert_eq!(
        read_midpoint_zarr_u32_1d(&zarr_path, "/length_end_bp")?,
        vec![50, 70, 100]
    );
    let group_metadata: Value =
        serde_json::from_str(&std::fs::read_to_string(zarr_path.join("group/zarr.json"))?)?;
    assert_eq!(group_metadata["attributes"]["label_field"], "group_name");
    assert_eq!(
        group_metadata["attributes"]["labels"],
        serde_json::json!(["LYL1", "beta-site", "gamma_long"])
    );

    let group_index = std::fs::read_to_string(group_index_path)?;
    assert_eq!(
        group_index,
        "group_idx\tgroup_name\teligible_intervals\n0\tLYL1\t2\n1\tbeta-site\t2\n2\tgamma_long\t2\n"
    );

    let settings: Value = serde_json::from_str(&std::fs::read_to_string(settings_path)?)?;
    assert_eq!(
        settings["array_axes"],
        serde_json::json!(["group", "length_bin", "position"])
    );
    assert_eq!(
        settings["length_axis"]["bin_definition"]["edges"],
        serde_json::json!([30, 50, 70, 100])
    );
    assert_eq!(settings["position_axis"]["output_interval_length_bp"], 10);
    assert_eq!(settings["position_axis"]["bin_size_bp"], 2);
    assert_eq!(settings["position_axis"]["n_bins"], 5);
    assert_eq!(settings["smoothing"]["method"], "none");

    Ok(())
}

#[test]
// This is ignored because it is not a normal correctness test. It generates real cfDNAlab
// end-motif Zarr outputs for the downstream Python/R compatibility workflow, which invokes it
// explicitly from the downstream-zarr GitHub Action.
#[ignore = "generates the Zarr fixtures consumed by downstream Python/R compatibility tests"]
fn generate_end_motif_zarr_fixtures_with_cfdnalab() -> Result<()> {
    let output_dir = env::var_os("CFDNALAB_DOWNSTREAM_FIXTURE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("downstream_tests/tmp"));
    std::fs::create_dir_all(&output_dir)?;

    let bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 256)],
        vec![paired_fragment(10, 10, 4)],
        Vec::new(),
        "downstream_end_motif_fixture",
    )?;
    let reference = simple_reference_twobit()?;

    let dense_global_path = run_end_fixture(
        &bam.bam,
        &reference.path,
        &output_dir,
        "tiny_dense_global",
        true,
        None,
    )?;
    let (dense_motifs, dense_counts) = read_dense_end_counts(&dense_global_path)?;
    assert_eq!(dense_motifs, vec!["_A", "_C", "_G", "_T"]);
    assert_eq!(dense_counts.shape(), &[1, 4]);
    assert_eq!(dense_counts[(0, 0)], 1.0);
    assert_eq!(dense_counts[(0, 1)], 0.0);
    assert_eq!(dense_counts[(0, 2)], 1.0);
    assert_eq!(dense_counts[(0, 3)], 0.0);

    let sparse_window_bed = output_dir.join("tiny_ends_windows.bed");
    write_bed(
        &sparse_window_bed,
        &[("chr1", 10, 11, "left"), ("chr1", 19, 20, "right")],
    )?;
    let sparse_window_path = run_end_fixture(
        &bam.bam,
        &reference.path,
        &output_dir,
        "tiny_sparse_windowed",
        false,
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(sparse_window_bed),
            by_grouped_bed: None,
        }),
    )?;
    let (window_motifs, window_counts) = read_sparse_end_counts(&sparse_window_path)?;
    assert_eq!(window_motifs, vec!["_A", "_G"]);
    assert_eq!(window_counts.shape(), &[2, 2]);
    assert_eq!(window_counts[(0, 1)], 1.0);
    assert_eq!(window_counts[(1, 0)], 1.0);
    assert_eq!(window_counts.sum(), 2.0);

    let sparse_grouped_bed = output_dir.join("tiny_ends_grouped.bed");
    write_bed(
        &sparse_grouped_bed,
        &[
            ("chr1", 10, 11, "beta"),
            ("chr1", 19, 20, "alpha"),
            ("chr1", 10, 20, "beta"),
            ("chr1", 30, 31, "gamma"),
        ],
    )?;
    let sparse_grouped_path = run_end_fixture(
        &bam.bam,
        &reference.path,
        &output_dir,
        "tiny_sparse_grouped",
        false,
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(sparse_grouped_bed),
        }),
    )?;
    let (grouped_motifs, grouped_counts) = read_sparse_end_counts(&sparse_grouped_path)?;
    assert_eq!(grouped_motifs, vec!["_A", "_G"]);
    assert_eq!(grouped_counts.shape(), &[3, 2]);
    assert_eq!(grouped_counts[(0, 0)], 1.0);
    assert_eq!(grouped_counts[(0, 1)], 2.0);
    assert_eq!(grouped_counts[(1, 0)], 1.0);
    assert_eq!(grouped_counts.row(2).sum(), 0.0);
    assert_eq!(grouped_counts.sum(), 4.0);
    let group_metadata: Value = serde_json::from_str(&std::fs::read_to_string(
        sparse_grouped_path.join("group/zarr.json"),
    )?)?;
    assert_eq!(
        group_metadata["attributes"]["labels"],
        serde_json::json!(["beta", "alpha", "gamma"])
    );

    Ok(())
}

fn run_end_fixture(
    bam_path: &Path,
    reference_path: &Path,
    output_dir: &Path,
    prefix: &str,
    all_motifs: bool,
    windows: Option<DistributionWindowsArgs>,
) -> Result<PathBuf> {
    let mut config = EndsConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
        1,
        0,
    );
    config.output_prefix = prefix.to_string();
    config.set_ref_2bit(Some(reference_path.to_path_buf()));
    config.source_inside = KmerSource::Reference;
    config.all_motifs = all_motifs;
    config.clip.clip_strategy = ClipStrategy::Aligned;
    config.set_min_mapq(0);
    config.set_tile_size(1_000_000);
    config.set_require_proper_pair(false);
    config.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    if let Some(windows) = windows {
        config.set_windows(windows);
    }
    {
        let lengths = config.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    run_ends(&config)?;
    let zarr_path = output_dir.join(format!("{prefix}.end_motifs.zarr"));
    assert!(zarr_path.is_dir());
    assert!(
        output_dir
            .join(format!("{prefix}.end_settings.json"))
            .is_file()
    );
    Ok(zarr_path)
}

fn read_dense_end_counts(store_path: &Path) -> Result<(Vec<String>, Array2<f64>)> {
    let motifs = read_motif_labels(store_path)?;
    let counts: Vec<f64> = read_zarr_array(store_path, "/counts")?;
    let shape = read_zarr_shape(store_path, "/counts")?;
    let matrix = Array2::from_shape_vec((shape[0], shape[1]), counts)?;
    Ok((motifs, matrix))
}

fn read_sparse_end_counts(store_path: &Path) -> Result<(Vec<String>, Array2<f64>)> {
    let motifs = read_motif_labels(store_path)?;
    let row: Vec<i32> = read_zarr_array(store_path, "/sparse/row")?;
    let motif: Vec<i32> = read_zarr_array(store_path, "/sparse/motif")?;
    let count: Vec<f64> = read_zarr_array(store_path, "/sparse/count")?;
    let shape: Vec<i32> = read_zarr_array(store_path, "/sparse/shape")?;
    let row_count = usize::try_from(shape[0]).context("sparse row count must be non-negative")?;
    let motif_count =
        usize::try_from(shape[1]).context("sparse motif count must be non-negative")?;
    let mut matrix = Array2::<f64>::zeros((row_count, motif_count));
    for ((row, motif), count) in row.into_iter().zip(motif).zip(count) {
        let row = usize::try_from(row).context("sparse row index must be non-negative")?;
        let motif = usize::try_from(motif).context("sparse motif index must be non-negative")?;
        matrix[(row, motif)] = count;
    }
    Ok((motifs, matrix))
}

fn read_motif_labels(store_path: &Path) -> Result<Vec<String>> {
    let motif_index: Vec<i32> = read_zarr_array(store_path, "/motif_index")?;
    let motif_byte: Vec<i32> = read_zarr_array(store_path, "/motif_byte")?;
    let motif_ascii: Vec<u8> = read_zarr_array(store_path, "/motif_ascii")?;
    let motif_width = motif_byte.len();
    assert_eq!(motif_ascii.len(), motif_index.len() * motif_width);
    if motif_width == 0 {
        return Ok(Vec::new());
    }
    motif_ascii
        .chunks_exact(motif_width)
        .map(|bytes| {
            String::from_utf8(bytes.to_vec()).context("motif_ascii row must be valid UTF-8")
        })
        .collect()
}

fn read_zarr_array<T>(store_path: &Path, array_path: &str) -> Result<Vec<T>>
where
    T: zarrs::array::ElementOwned,
{
    let store = Arc::new(FilesystemStore::new(store_path)?);
    let array = Array::open(store, array_path)?;
    Ok(array.retrieve_array_subset(&array.subset_all())?)
}

fn read_zarr_shape(store_path: &Path, array_path: &str) -> Result<Vec<usize>> {
    let store = Arc::new(FilesystemStore::new(store_path)?);
    let array = Array::open(store, array_path)?;
    Ok(array
        .shape()
        .iter()
        .map(|dimension| *dimension as usize)
        .collect())
}
