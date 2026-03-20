#![cfg(feature = "cmd_midpoints")]

mod fixtures;

use anyhow::Result;
use cfdnalab::commands::cli_common::{ChromosomeArgs, IOCArgs, ScaleGenomeArgs};
use cfdnalab::commands::midpoints::config::MidpointsConfig;
use cfdnalab::commands::midpoints::midpoints::run;
use fixtures::{FragmentSpec, ReadSpec, bam_from_specs, complex_bam_fixture, write_bed};
use ndarray::Array3;
use ndarray_npy::read_npy;
use std::collections::HashMap;
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

fn paired_fragment_on_tid(
    tid: usize,
    start: i64,
    fragment_len: i64,
    read_len: i64,
) -> FragmentSpec {
    const FLAG_FIRST_MATE: u16 = 0x40;
    const FLAG_SECOND_MATE: u16 = 0x80;
    const FLAG_PROPER_PAIR: u16 = 0x2;
    const FLAG_MATE_REVERSE: u16 = 0x20;

    let reverse_start = start + fragment_len - read_len;
    let insert_size = fragment_len;
    FragmentSpec {
        forward: ReadSpec {
            tid,
            pos: start,
            cigar: vec![('M', read_len as u32)],
            seq: vec![b'A'; read_len as usize],
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: FLAG_FIRST_MATE | FLAG_MATE_REVERSE | FLAG_PROPER_PAIR,
            mate_tid: Some(tid),
            mate_pos: Some(reverse_start),
            insert_size,
        },
        reverse: ReadSpec {
            tid,
            pos: reverse_start,
            cigar: vec![('M', read_len as u32)],
            seq: vec![b'T'; read_len as usize],
            qual: 40,
            is_reverse: true,
            mapq: 60,
            flags: FLAG_SECOND_MATE | FLAG_PROPER_PAIR,
            mate_tid: Some(tid),
            mate_pos: Some(start),
            insert_size: -insert_size,
        },
    }
}

fn read_group_index_map(path: &std::path::Path) -> Result<HashMap<String, usize>> {
    let text = std::fs::read_to_string(path)?;
    let mut out = HashMap::new();
    for line in text.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let idx = fields.next().unwrap().parse::<usize>()?;
        let name = fields.next().unwrap().to_string();
        out.insert(name, idx);
    }
    Ok(out)
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

#[test]
fn midpoint_fetch_narrowing_preserves_tile_halo_near_chromosome_end_on_three_chromosomes()
-> Result<()> {
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 95),
            ("chr2".to_string(), 95),
            ("chr3".to_string(), 95),
        ],
        vec![
            paired_fragment_on_tid(0, 84, 11, 3),
            paired_fragment_on_tid(1, 84, 11, 3),
            paired_fragment_on_tid(2, 84, 11, 3),
        ],
        Vec::new(),
        "midpoints_chrom_end_halo_three_chr",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp.path().join("windows_three_chr_near_end.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 89, 95, "groupA"),
            ("chr2", 89, 95, "groupB"),
            ("chr3", 89, 95, "groupC"),
        ],
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1", "chr2", "chr3"]),
        bed_path.clone(),
    );
    cfg.set_output_prefix("sites");
    cfg.set_length_bins(vec![10, 15]);
    cfg.set_tile_size(40);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());

    // Manual expectations:
    // - Each chromosome ends with a 6 bp site [89,95), which falls in the last tile [80,95).
    // - The only fragment on each chromosome is [84,95), length 11, midpoint 89.
    // - The midpoint lies at window position 89 - 89 = 0, so each group gets one count at
    //   length-bin [10,15) and position 0.
    // - This command-level fixture checks that narrowing to the extreme midpoint sites does not
    //   discard the fetch halo already carried by the last tile near chromosome end.
    // - It does not isolate the separate `halo_bp` argument to the narrowing helper, because the
    //   tile fetch band was already built with the same maximum-fragment-length halo.
    run(&cfg)?;

    let counts_path = temp.path().join("sites.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
    assert_eq!(arr.shape(), &[3, 1, 6]);

    let map_path = temp.path().join("sites.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;

    assert_eq!(group_to_idx.len(), 3);
    assert_eq!(arr.sum(), 3.0);
    for group_name in ["groupA", "groupB", "groupC"] {
        let group_idx = group_to_idx[group_name];
        let row = arr.slice(ndarray::s![group_idx, 0, ..]).to_vec();
        assert_eq!(
            row,
            vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "{group_name} should have exactly one midpoint count at position 0"
        );
    }

    Ok(())
}

#[test]
fn midpoint_fetch_narrowing_reads_all_eligible_fragments_near_chromosome_end_on_three_chromosomes()
-> Result<()> {
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 95),
            ("chr2".to_string(), 95),
            ("chr3".to_string(), 95),
        ],
        vec![
            paired_fragment_on_tid(0, 79, 11, 3),
            paired_fragment_on_tid(0, 80, 11, 3),
            paired_fragment_on_tid(0, 82, 11, 3),
            paired_fragment_on_tid(0, 84, 11, 3),
            paired_fragment_on_tid(1, 79, 11, 3),
            paired_fragment_on_tid(1, 80, 11, 3),
            paired_fragment_on_tid(1, 82, 11, 3),
            paired_fragment_on_tid(1, 84, 11, 3),
            paired_fragment_on_tid(2, 79, 11, 3),
            paired_fragment_on_tid(2, 80, 11, 3),
            paired_fragment_on_tid(2, 82, 11, 3),
            paired_fragment_on_tid(2, 84, 11, 3),
        ],
        Vec::new(),
        "midpoints_chrom_end_fetch_reads_all_eligible",
    )?;
    let temp = TempDir::new()?;
    let bed_path = temp
        .path()
        .join("windows_three_chr_fetch_read_coverage.bed");
    write_bed(
        &bed_path,
        &[
            ("chr1", 85, 95, "groupA"),
            ("chr2", 85, 95, "groupB"),
            ("chr3", 85, 95, "groupC"),
        ],
    )?;

    let mut cfg = MidpointsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: temp.path().to_path_buf(),
            n_threads: 2,
        },
        base_chromosomes(&["chr1", "chr2", "chr3"]),
        bed_path,
    );
    cfg.set_output_prefix("sites_fetch_reads_all");
    cfg.set_length_bins(vec![10, 15]);
    cfg.set_tile_size(40);
    cfg.set_min_mapq(0);
    cfg.set_require_proper_pair(false);
    cfg.set_scale_genome(ScaleGenomeArgs::default());

    // Manual expectations:
    // - Each chromosome has one site [85,95), which lies in the last tile [80,95).
    // - Four fragments are present per chromosome, all length 11:
    //     * [79,90) midpoint 84 -> outside the site, so it must not be counted
    //     * [80,91) midpoint 85 -> counted at site position 0
    //     * [82,93) midpoint 87 -> counted at site position 2
    //     * [84,95) midpoint 89 -> counted at site position 4
    // - The narrowing step therefore has to preserve enough of the tile fetch band to read all
    //   three eligible fragments, not just the one closest to chromosome end.
    // - Each group row must therefore be exactly [1,0,1,0,1,0,0,0,0,0].
    run(&cfg)?;

    let counts_path = temp
        .path()
        .join("sites_fetch_reads_all.midpoint_profiles.npy");
    let arr: Array3<f32> = read_npy(&counts_path)?;
    assert_eq!(arr.shape(), &[3, 1, 10]);

    let map_path = temp.path().join("sites_fetch_reads_all.group_index.tsv");
    let group_to_idx = read_group_index_map(&map_path)?;

    assert_eq!(group_to_idx.len(), 3);
    assert_eq!(arr.sum(), 9.0);
    for group_name in ["groupA", "groupB", "groupC"] {
        let group_idx = group_to_idx[group_name];
        let row = arr.slice(ndarray::s![group_idx, 0, ..]).to_vec();
        assert_eq!(
            row,
            vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "{group_name} should count exactly the three eligible near-end fragments"
        );
    }

    Ok(())
}
