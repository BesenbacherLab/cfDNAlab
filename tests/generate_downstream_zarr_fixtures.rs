#![cfg(all(
    feature = "cmd_midpoints",
    feature = "cmd_ends",
    feature = "cmd_lengths"
))]

mod fixtures;

use anyhow::{Context, Result};
use cfdnalab::commands::{
    cli_common::{ChromosomeArgs, DistributionWindowsArgs, IOCArgs},
    ends::{
        config::EndsConfig,
        config_structs::{AssignMotifToWindowArgs, ClipStrategy, KmerSource, WindowMotifAssigner},
        ends::run as run_ends,
    },
    lengths::{config::LengthsConfig, lengths::run as run_lengths},
    midpoints::{
        config::MidpointsConfig, midpoints::run as run_midpoints, smoothing::MidpointSmoothing,
    },
};
use fixtures::{
    bam_from_specs_strict_identity, paired_fragment, read_length_counts_text,
    read_midpoint_zarr_counts, read_midpoint_zarr_i32_1d, twobit_from_sequences, write_bed,
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
        read_midpoint_zarr_i32_1d(&zarr_path, "/eligible_intervals")?,
        vec![2, 2, 2]
    );
    assert_eq!(
        read_midpoint_zarr_i32_1d(&zarr_path, "/length_start_bp")?,
        vec![30, 50, 70]
    );
    assert_eq!(
        read_midpoint_zarr_i32_1d(&zarr_path, "/length_end_bp")?,
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

    let mut chr2_fragment = paired_fragment(10, 10, 4);
    chr2_fragment.forward.tid = 1;
    chr2_fragment.forward.mate_tid = Some(1);
    chr2_fragment.reverse.tid = 1;
    chr2_fragment.reverse.mate_tid = Some(1);

    let bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 256), ("chr2".to_string(), 256)],
        vec![paired_fragment(10, 10, 4), chr2_fragment],
        Vec::new(),
        "downstream_end_motif_fixture",
    )?;
    let reference = twobit_from_sequences(
        "downstream_end_reference",
        vec![
            ("chr1".to_string(), "ACGT".repeat(64)),
            ("chr2".to_string(), "ACGT".repeat(64)),
        ],
    )?;

    let dense_global_path = run_end_fixture(
        &bam.bam,
        &reference.path,
        &output_dir,
        "tiny_dense_global",
        &["chr1"],
        1,
        0,
        true,
        None,
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
        &[
            ("chr1", 10, 11, "left"),
            ("chr1", 19, 20, "right"),
            ("chr2", 10, 11, "chr2_left"),
        ],
    )?;
    let sparse_window_path = run_end_fixture(
        &bam.bam,
        &reference.path,
        &output_dir,
        "tiny_sparse_windowed",
        &["chr1", "chr2"],
        1,
        0,
        false,
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(sparse_window_bed),
            by_grouped_bed: None,
        }),
        None,
    )?;
    let (window_motifs, window_counts) = read_sparse_end_counts(&sparse_window_path)?;
    assert_eq!(window_motifs, vec!["_A", "_G"]);
    assert_eq!(window_counts.shape(), &[3, 2]);
    assert_eq!(window_counts[(0, 1)], 1.0);
    assert_eq!(window_counts[(1, 0)], 1.0);
    assert_eq!(window_counts[(2, 1)], 1.0);
    assert_eq!(window_counts.sum(), 3.0);

    let selected_motifs_file = output_dir.join("tiny_ends_selected_motifs.tsv");
    std::fs::write(&selected_motifs_file, "GT_AC\nAC_GT\nTT_TT\n")?;
    let sparse_windowed_selected_motifs_path = run_end_fixture(
        &bam.bam,
        &reference.path,
        &output_dir,
        "tiny_sparse_windowed_selected_motifs",
        &["chr1", "chr2"],
        2,
        2,
        false,
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(output_dir.join("tiny_ends_windows.bed")),
            by_grouped_bed: None,
        }),
        Some(selected_motifs_file),
    )?;
    let (selected_motifs, selected_counts) =
        read_sparse_end_counts(&sparse_windowed_selected_motifs_path)?;
    // Sparse motifs-file output keeps only observed motifs unless --all-motifs is set.
    assert_eq!(selected_motifs, vec!["GT_AC", "AC_GT"]);
    assert_eq!(selected_counts.shape(), &[3, 2]);
    assert_eq!(selected_counts[(0, 1)], 1.0);
    assert_eq!(selected_counts[(1, 0)], 1.0);
    assert_eq!(selected_counts[(2, 1)], 1.0);
    assert_eq!(selected_counts.sum(), 3.0);
    let selected_motif_root_metadata: Value = serde_json::from_str(&std::fs::read_to_string(
        sparse_windowed_selected_motifs_path.join("zarr.json"),
    )?)?;
    assert_eq!(
        selected_motif_root_metadata["attributes"]["motif_axis_kind"],
        serde_json::json!("motif")
    );
    assert!(
        sparse_windowed_selected_motifs_path
            .join("motif_ascii")
            .is_dir()
    );

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
        &["chr1"],
        1,
        0,
        false,
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(sparse_grouped_bed.clone()),
        }),
        None,
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

    let motif_groups_file = output_dir.join("tiny_ends_motif_groups.tsv");
    std::fs::write(
        &motif_groups_file,
        "G\tleft-hit\nA\tright-hit\nC\tunused-hit\n",
    )?;
    let sparse_grouped_motif_groups_path = run_end_fixture(
        &bam.bam,
        &reference.path,
        &output_dir,
        "tiny_sparse_grouped_motif_groups",
        &["chr1"],
        1,
        0,
        false,
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(sparse_grouped_bed),
        }),
        Some(motif_groups_file),
    )?;
    let motif_group_labels = read_motif_group_labels(&sparse_grouped_motif_groups_path)?;
    let motif_group_counts = read_sparse_count_matrix(&sparse_grouped_motif_groups_path)?;
    assert_eq!(motif_group_labels, vec!["left-hit", "right-hit"]);
    assert_eq!(motif_group_counts.shape(), &[3, 2]);
    assert_eq!(motif_group_counts[(0, 0)], 2.0);
    assert_eq!(motif_group_counts[(0, 1)], 1.0);
    assert_eq!(motif_group_counts[(1, 0)], 0.0);
    assert_eq!(motif_group_counts[(1, 1)], 1.0);
    assert_eq!(motif_group_counts.row(2).sum(), 0.0);
    assert_eq!(motif_group_counts.sum(), 4.0);
    let motif_group_root_metadata: Value = serde_json::from_str(&std::fs::read_to_string(
        sparse_grouped_motif_groups_path.join("zarr.json"),
    )?)?;
    assert_eq!(
        motif_group_root_metadata["attributes"]["motif_axis_kind"],
        serde_json::json!("motif_group")
    );

    let wide_motif_groups_file = output_dir.join("tiny_ends_wide_motif_groups.tsv");
    std::fs::write(
        &wide_motif_groups_file,
        "GT_AC\tright-hit-wide\nAC_GT\tleft-hit-wide\nTT_TT\tunused-wide\n",
    )?;
    let sparse_grouped_wide_motif_groups_path = run_end_fixture(
        &bam.bam,
        &reference.path,
        &output_dir,
        "tiny_sparse_grouped_wide_motif_groups",
        &["chr1"],
        2,
        2,
        false,
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(output_dir.join("tiny_ends_grouped.bed")),
        }),
        Some(wide_motif_groups_file),
    )?;
    let wide_motif_group_labels = read_motif_group_labels(&sparse_grouped_wide_motif_groups_path)?;
    let wide_motif_group_counts = read_sparse_count_matrix(&sparse_grouped_wide_motif_groups_path)?;
    assert_eq!(
        wide_motif_group_labels,
        vec!["right-hit-wide", "left-hit-wide"]
    );
    assert_eq!(wide_motif_group_counts.shape(), &[3, 2]);
    assert_eq!(wide_motif_group_counts[(0, 0)], 1.0);
    assert_eq!(wide_motif_group_counts[(0, 1)], 2.0);
    assert_eq!(wide_motif_group_counts[(1, 0)], 1.0);
    assert_eq!(wide_motif_group_counts.row(2).sum(), 0.0);
    assert_eq!(wide_motif_group_counts.sum(), 4.0);
    let wide_motif_group_root_metadata: Value = serde_json::from_str(&std::fs::read_to_string(
        sparse_grouped_wide_motif_groups_path.join("zarr.json"),
    )?)?;
    assert_eq!(
        wide_motif_group_root_metadata["attributes"]["motif_axis_kind"],
        serde_json::json!("motif_group")
    );
    assert!(
        !sparse_grouped_wide_motif_groups_path
            .join("motif_ascii")
            .exists()
    );

    Ok(())
}

#[test]
// This is ignored because it is not a normal correctness test. It generates real cfDNAlab
// length-count TSV outputs for the downstream R compatibility workflow, which invokes it
// explicitly from the downstream fixture GitHub Action.
#[ignore = "generates the length-count fixtures consumed by downstream R compatibility tests"]
fn generate_length_count_tsv_fixtures_with_cfdnalab() -> Result<()> {
    let output_dir = env::var_os("CFDNALAB_DOWNSTREAM_FIXTURE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("downstream_tests/tmp"));
    std::fs::create_dir_all(&output_dir)?;

    let bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 500)],
        vec![
            paired_fragment(10, 35, 20),
            paired_fragment(50, 45, 20),
            paired_fragment(120, 61, 20),
            paired_fragment(130, 65, 20),
            paired_fragment(215, 81, 20),
            paired_fragment(320, 35, 20),
        ],
        Vec::new(),
        "downstream_length_fixture",
    )?;

    let window_bed = output_dir.join("tiny_lengths_windows.bed");
    write_bed(
        &window_bed,
        &[
            ("chr1", 0, 100, "left"),
            ("chr1", 100, 200, "middle"),
            ("chr1", 200, 300, "right"),
            ("chr1", 300, 360, "tail"),
        ],
    )?;
    let grouped_bed = output_dir.join("tiny_lengths_grouped.bed");
    write_bed(
        &grouped_bed,
        &[
            ("chr1", 0, 100, "beta"),
            ("chr1", 100, 200, "alpha"),
            ("chr1", 200, 300, "beta"),
            ("chr1", 300, 360, "gamma"),
            ("chr1", 360, 390, "zero"),
        ],
    )?;
    let blacklist_bed = output_dir.join("tiny_lengths_blacklist.bed");
    write_bed(
        &blacklist_bed,
        &[
            ("chr1", 96, 100, "mask_left"),
            ("chr1", 100, 105, "mask_middle"),
            ("chr1", 200, 210, "mask_right"),
            ("chr1", 300, 315, "mask_tail"),
            ("chr1", 360, 370, "mask_zero"),
        ],
    )?;

    let global_path = run_length_fixture(&bam.bam, &output_dir, "tiny_lengths_global", None, None)?;
    assert_eq!(
        read_length_counts_text(&global_path)?,
        "count_30_50\tcount_50_70\tcount_70_100\n3\t2\t1\n"
    );

    let windowed_path = run_length_fixture(
        &bam.bam,
        &output_dir,
        "tiny_lengths_windowed",
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(window_bed.clone()),
            by_grouped_bed: None,
        }),
        Some(blacklist_bed.clone()),
    )?;
    assert_eq!(
        read_length_counts_text(&windowed_path)?,
        concat!(
            "chrom\tstart\tend\tblacklisted_fraction\tcount_30_50\tcount_50_70\tcount_70_100\n",
            "chr1\t0\t100\t0.04\t2\t0\t0\n",
            "chr1\t100\t200\t0.05\t0\t2\t0\n",
            "chr1\t200\t300\t0.1\t0\t0\t1\n",
            "chr1\t300\t360\t0.25\t1\t0\t0\n",
        )
    );

    let windowed_no_blacklist_path = run_length_fixture(
        &bam.bam,
        &output_dir,
        "tiny_lengths_windowed_no_blacklist",
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: Some(window_bed),
            by_grouped_bed: None,
        }),
        None,
    )?;
    assert_eq!(
        read_length_counts_text(&windowed_no_blacklist_path)?,
        concat!(
            "chrom\tstart\tend\tcount_30_50\tcount_50_70\tcount_70_100\n",
            "chr1\t0\t100\t2\t0\t0\n",
            "chr1\t100\t200\t0\t2\t0\n",
            "chr1\t200\t300\t0\t0\t1\n",
            "chr1\t300\t360\t1\t0\t0\n",
        )
    );

    let grouped_path = run_length_fixture(
        &bam.bam,
        &output_dir,
        "tiny_lengths_grouped",
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed.clone()),
        }),
        Some(blacklist_bed),
    )?;
    assert_eq!(
        read_length_counts_text(&grouped_path)?,
        concat!(
            "group_name\teligible_windows\tblacklisted_fraction\tcount_30_50\tcount_50_70\tcount_70_100\n",
            "beta\t2\t0.07\t2\t0\t1\n",
            "alpha\t1\t0.05\t0\t2\t0\n",
            "gamma\t1\t0.25\t1\t0\t0\n",
            "zero\t1\t0.333\t0\t0\t0\n",
        )
    );

    let grouped_no_blacklist_path = run_length_fixture(
        &bam.bam,
        &output_dir,
        "tiny_lengths_grouped_no_blacklist",
        Some(DistributionWindowsArgs {
            by_size: None,
            by_bed: None,
            by_grouped_bed: Some(grouped_bed),
        }),
        None,
    )?;
    assert_eq!(
        read_length_counts_text(&grouped_no_blacklist_path)?,
        concat!(
            "group_name\teligible_windows\tcount_30_50\tcount_50_70\tcount_70_100\n",
            "beta\t2\t2\t0\t1\n",
            "alpha\t1\t0\t2\t0\n",
            "gamma\t1\t1\t0\t0\n",
            "zero\t1\t0\t0\t0\n",
        )
    );

    Ok(())
}

fn run_length_fixture(
    bam_path: &Path,
    output_dir: &Path,
    prefix: &str,
    windows: Option<DistributionWindowsArgs>,
    blacklist: Option<PathBuf>,
) -> Result<PathBuf> {
    let mut config = LengthsConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    );
    config.output_prefix = prefix.to_string();
    config.set_length_bins(vec![30, 50, 70, 100]);
    config.set_min_mapq(0);
    config.set_require_proper_pair(false);
    config.set_tile_size(1_000_000);
    if let Some(windows) = windows {
        config.set_windows(windows);
    }
    if let Some(blacklist) = blacklist {
        config.blacklist = Some(vec![blacklist]);
    }

    run_lengths(&config)?;
    let counts_path = output_dir.join(format!("{prefix}.length_counts.tsv.zst"));
    assert!(counts_path.is_file());
    assert!(
        output_dir
            .join(format!("{prefix}.length_settings.json"))
            .is_file()
    );
    Ok(counts_path)
}

fn run_end_fixture(
    bam_path: &Path,
    reference_path: &Path,
    output_dir: &Path,
    prefix: &str,
    chromosomes: &[&str],
    k_inside: usize,
    k_outside: usize,
    all_motifs: bool,
    windows: Option<DistributionWindowsArgs>,
    motifs_file: Option<PathBuf>,
) -> Result<PathBuf> {
    let mut config = EndsConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(
                chromosomes
                    .iter()
                    .map(|chromosome| (*chromosome).to_string())
                    .collect(),
            ),
            chromosomes_file: None,
        },
        k_inside,
        k_outside,
    );
    config.output_prefix = prefix.to_string();
    config.set_ref_2bit(Some(reference_path.to_path_buf()));
    config.source_inside = KmerSource::Reference;
    config.all_motifs = all_motifs;
    config.motifs_file = motifs_file;
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
    let matrix = read_sparse_count_matrix(store_path)?;
    Ok((motifs, matrix))
}

fn read_sparse_count_matrix(store_path: &Path) -> Result<Array2<f64>> {
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
    Ok(matrix)
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

fn read_motif_group_labels(store_path: &Path) -> Result<Vec<String>> {
    let motif_metadata: Value = serde_json::from_str(&std::fs::read_to_string(
        store_path.join("motif_index/zarr.json"),
    )?)?;
    assert_eq!(
        motif_metadata["attributes"]["label_field"],
        serde_json::json!("motif_group")
    );
    motif_metadata["attributes"]["labels"]
        .as_array()
        .context("motif_group labels must be a JSON array")?
        .iter()
        .map(|label| {
            label
                .as_str()
                .map(ToString::to_string)
                .context("motif_group label must be a string")
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
