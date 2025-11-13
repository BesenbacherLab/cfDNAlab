#![cfg(feature = "cmd_midpoints")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs, ScaleGenomeArgs};
use cfdnalab::commands::profile_groups::config::ProfileGroupsConfig;
use cfdnalab::commands::profile_groups::profile_groups::run;
use fixtures::{complex_bam_fixture, write_bed};
use ndarray::Array3;
use ndarray_npy::read_npy;
use tempfile::TempDir;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|c| c.to_string()).collect()),
        chromosomes_file: None,
    }
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

    let mut cfg = ProfileGroupsConfig::new(
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
