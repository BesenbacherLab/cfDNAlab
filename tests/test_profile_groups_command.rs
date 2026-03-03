#![cfg(feature = "cmd_midpoints")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs, ScaleGenomeArgs};
use cfdnalab::commands::midpoints::config::MidpointsConfig;
use cfdnalab::commands::midpoints::midpoints::run;
use fixtures::{complex_bam_fixture, write_bed};
use ndarray::Array3;
use ndarray_npy::read_npy;
use std::path::PathBuf;
use tempfile::TempDir;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
        chromosomes_file: None,
    }
}

fn base_midpoints_config_for_length_bins() -> MidpointsConfig {
    MidpointsConfig::new(
        IOCArgs {
            bam: PathBuf::from("dummy.bam"),
            output_dir: PathBuf::from("out"),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        PathBuf::from("intervals.bed"),
    )
}

#[test]
fn length_bin_range_spec_matches_brace_expansion_edges() -> Result<()> {
    // Arrange: Hand-derived expected edges for 100..220 with step 10.
    // The end is an edge (not a counted length), so we expect:
    // 100, 110, 120, ..., 220.
    let expected_edges = vec![
        100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220,
    ];

    let mut edge_list_config = base_midpoints_config_for_length_bins();
    edge_list_config.set_length_bins(expected_edges.clone());

    let mut range_spec_config = base_midpoints_config_for_length_bins();
    range_spec_config.set_length_bins_spec("100:220:10");

    // Act
    let edges_from_edge_list = edge_list_config.resolve_length_bins()?;
    let edges_from_range_spec = range_spec_config.resolve_length_bins()?;

    // Assert
    assert_eq!(edges_from_edge_list, expected_edges);
    assert_eq!(edges_from_range_spec, expected_edges);
    assert_eq!(edges_from_edge_list, edges_from_range_spec);

    Ok(())
}

#[test]
fn length_bin_start_end_list_format_is_rejected() {
    // Arrange: This format was intentionally removed.
    let mut config = base_midpoints_config_for_length_bins();
    config.set_length_bins_spec("30-80,80-150");

    // Act
    let error = config
        .resolve_length_bins()
        .expect_err("start-end list format should fail");

    // Assert
    assert!(
        format!("{error}").contains("explicit start-end lists are not supported"),
        "Unexpected error message: {error}"
    );
}

#[test]
fn midpoint_profiles_written_with_group_index() -> Result<()> {
    let bam = complex_bam_fixture()?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 40, 80, "groupA"),
            ("chr1", 180, 220, "groupA"),
            ("chr2", 20, 60, "groupB"),
            ("chr2", 60, 100, "groupB"),
        ],
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1", "chr2"]),
        bed_path.clone(),
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![20, 60, 120]);
    cfg.set_tile_size(1_000);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());

    run(&cfg)?;

    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    assert!(counts_path.exists());
    let arr: Array3<f32> = read_npy(&counts_path)?;
    assert_eq!(arr.shape(), &[2, 2, 40]); // groups, length bins, window size
    assert!(arr.sum() > 0.0);

    let map_path = temp.path().join("sites.group_index.tsv");
    let map_text = std::fs::read_to_string(&map_path)?;
    assert!(map_text.contains("groupA"));
    assert!(map_text.contains("groupB"));

    Ok(())
}
