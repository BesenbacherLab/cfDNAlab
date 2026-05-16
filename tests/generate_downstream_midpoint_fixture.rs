#![cfg(feature = "cmd_midpoints")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::{
    cli_common::{ChromosomeArgs, IOCArgs},
    midpoints::{
        config::MidpointsConfig, midpoints::run as run_midpoints, smoothing::MidpointSmoothing,
    },
};
use fixtures::{
    bam_from_specs_strict_identity, paired_fragment, read_midpoint_zarr_counts,
    read_midpoint_zarr_i32_1d, read_midpoint_zarr_u32_1d, write_bed,
};
use ndarray::arr3;
use serde_json::Value;
use std::{env, path::PathBuf};

#[test]
// This is ignored because it is not a normal correctness test. It generates a real cfDNAlab
// midpoint Zarr output for the downstream Python/R compatibility workflow, which invokes it
// explicitly with `cargo test -- --ignored --exact generate_midpoint_zarr_fixture_with_cfdnalab`.
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
            // alpha, length bin [30, 50)
            paired_fragment(26, 41, 20),
            paired_fragment(26, 41, 20),
            paired_fragment(28, 41, 20),
            // alpha, length bin [50, 70)
            paired_fragment(95, 61, 20),
            paired_fragment(95, 61, 20),
            paired_fragment(95, 61, 20),
            paired_fragment(97, 61, 20),
            // alpha, length bin [70, 100)
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
            ("chr1", 45, 55, "alpha"),
            ("chr1", 120, 130, "alpha"),
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
    assert_eq!(
        read_midpoint_zarr_u32_1d(&zarr_path, "/group_name_nbytes")?,
        vec![5, 9, 10]
    );

    let group_index = std::fs::read_to_string(group_index_path)?;
    assert_eq!(
        group_index,
        "group_idx\tgroup_name\teligible_intervals\n0\talpha\t2\n1\tbeta-site\t2\n2\tgamma_long\t2\n"
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
