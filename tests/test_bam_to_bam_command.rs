#![cfg(feature = "cmd_bam_to_bam")]

mod fixtures;

use std::{collections::HashMap, fs, path::Path};

use anyhow::Result;
use cfdnalab::commands::{
    bam_to_bam::{bam_to_bam::run_inner, config::BamToBamConfig},
    cli_common::ChromosomeArgs,
};
use fixtures::{FragmentSpec, ReadSpec, bam_from_specs, paired_fragment};
use rust_htslib::bam::{self, Read, record::Aux};
use tempfile::tempdir;

#[test]
fn filters_on_mapping_quality_and_fragment_membership() -> Result<()> {
    let surviving = paired_fragment(50, 160, 40);
    let mut low_mapq = paired_fragment(300, 160, 40);
    low_mapq.forward.mapq = 5;
    low_mapq.reverse.mapq = 30; // One read passes but the entire fragment is removed

    let orphan = orphan_read(500);

    let bam = bam_from_specs(
        vec![("chr1".to_string(), 1000)],
        vec![surviving, low_mapq],
        vec![orphan],
        "bam_to_bam_mapq",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("filtered.bam");

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.min_mapq = 30;

    run_inner(&cfg)?;

    let counts = read_qname_counts(&out_bam)?;
    assert_eq!(
        counts.get("frag0_50"),
        Some(&2),
        "Expected both mates of the surviving fragment"
    );
    assert_eq!(
        counts.len(),
        1,
        "Fragments with low MAPQ or missing mates must be discarded"
    );

    let lengths = read_fragment_lengths(&out_bam)?;
    assert!(
        lengths.iter().all(|&len| len == 160),
        "Surviving reads must carry the FLEN AUX tag"
    );

    Ok(())
}

#[test]
fn blacklisting_removes_fragment_when_single_mate_overlaps() -> Result<()> {
    let safe = paired_fragment(10, 120, 40);
    let with_blacklisted_mate = paired_fragment(200, 120, 40);

    let bam = bam_from_specs(
        vec![("chr1".to_string(), 1000)],
        vec![safe, with_blacklisted_mate],
        Vec::new(),
        "bam_to_bam_blacklist",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("blacklisted.bam");
    let blacklist = work.path().join("blacklist.bed");
    fs::write(&blacklist, "chr1\t280\t290\n")?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.blacklist = Some(vec![blacklist]);

    run_inner(&cfg)?;

    let counts = read_qname_counts(&out_bam)?;
    assert_eq!(
        counts.get("frag0_10"),
        Some(&2),
        "Clean fragments should remain after blacklist filtering"
    );
    assert!(
        !counts.contains_key("frag0_200"),
        "Fragments with a mate in a blacklisted region must be removed entirely"
    );
    assert_eq!(counts.len(), 1);

    Ok(())
}

#[test]
fn blacklisting_removes_fragment_when_gap_is_blacklisted() -> Result<()> {
    let safe = paired_fragment(10, 160, 40);
    let gap_hit = paired_fragment(300, 200, 40);

    let bam = bam_from_specs(
        vec![("chr1".to_string(), 1000)],
        vec![safe, gap_hit],
        Vec::new(),
        "bam_to_bam_blacklist_gap",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("blacklisted_gap.bam");
    let blacklist = work.path().join("blacklist_gap.bed");
    fs::write(&blacklist, "chr1\t380\t390\n")?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.blacklist = Some(vec![blacklist]);

    run_inner(&cfg)?;

    let counts = read_qname_counts(&out_bam)?;
    assert!(
        !counts.contains_key("frag0_300"),
        "Fragments with blacklisted gaps must be removed"
    );
    assert_eq!(
        counts.get("frag0_10"),
        Some(&2),
        "Unrelated fragments should survive"
    );

    Ok(())
}

#[test]
fn respects_chromosome_sort_toggle() -> Result<()> {
    let frag_chr2 = paired_fragment(10, 160, 40);
    let frag_chr10 = fragment_on_tid(paired_fragment(20, 160, 40), 1);
    let chroms = vec![("chr2".to_string(), 500), ("chr10".to_string(), 500)];
    let bam = bam_from_specs(
        chroms,
        vec![frag_chr2, frag_chr10],
        Vec::new(),
        "chrom_sort",
    )?;

    let work = tempdir()?;
    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec!["chr2".to_string(), "chr10".to_string()]),
        chromosomes_file: None,
    };

    let sorted_out = work.path().join("sorted.bam");
    let cfg_sorted = BamToBamConfig::new(bam.bam.clone(), sorted_out.clone(), chrom_args.clone());
    run_inner(&cfg_sorted)?;
    assert_eq!(
        first_record_chrom(&sorted_out)?,
        "chr10",
        "Default sorting should reorder chromosomes lexicographically"
    );

    let unsorted_out = work.path().join("unsorted.bam");
    let mut cfg_unsorted = BamToBamConfig::new(bam.bam.clone(), unsorted_out.clone(), chrom_args);
    cfg_unsorted.skip_chromosome_sort = true;
    run_inner(&cfg_unsorted)?;
    assert_eq!(
        first_record_chrom(&unsorted_out)?,
        "chr2",
        "Skipping chromosome sort should keep the provided order"
    );

    Ok(())
}

#[test]
fn windows_keep_fragment_when_single_read_is_inside() -> Result<()> {
    let target = paired_fragment(20, 200, 40);
    let outside = paired_fragment(400, 200, 40);
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 1000)],
        vec![target, outside],
        Vec::new(),
        "window_single_read",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("window_single.bam");
    let bed = work.path().join("windows_single.bed");
    write_bed(&bed, &[(0, 60)])?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.by_bed = Some(bed);

    run_inner(&cfg)?;
    let counts = read_qname_counts(&out_bam)?;
    assert!(
        counts.contains_key("frag0_20"),
        "Fragments overlapping the window through just one read must survive"
    );
    assert!(
        !counts.contains_key("frag0_400"),
        "Fragments outside the defined windows must be filtered"
    );

    Ok(())
}

#[test]
fn windows_count_overlap_spanning_gap_between_mates() -> Result<()> {
    let target = paired_fragment(0, 200, 40);
    let outside = paired_fragment(400, 200, 40);
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 1000)],
        vec![target, outside],
        Vec::new(),
        "window_gap",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("window_gap.bam");
    let bed = work.path().join("windows_gap.bed");
    write_bed(&bed, &[(80, 90)])?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.by_bed = Some(bed);

    run_inner(&cfg)?;
    let counts = read_qname_counts(&out_bam)?;
    assert!(
        counts.contains_key("frag0_0"),
        "Fragments should be counted when windows intersect the gap between mates"
    );
    assert!(
        !counts.contains_key("frag0_400"),
        "Fragments without any overlap should be removed"
    );

    Ok(())
}

#[test]
fn global_mode_keeps_expected_fragments_across_three_chromosomes() -> Result<()> {
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 1_000),
            ("chr2".to_string(), 1_000),
            ("chr3".to_string(), 1_000),
        ],
        vec![
            paired_fragment(10, 120, 40),
            fragment_on_tid(paired_fragment(30, 120, 40), 1),
            fragment_on_tid(paired_fragment(50, 120, 40), 2),
        ],
        Vec::new(),
        "bam_to_bam_three_chr_global",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("three_chr_global.bam");
    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec![
            "chr1".to_string(),
            "chr2".to_string(),
            "chr3".to_string(),
        ]),
        chromosomes_file: None,
    };
    let mut cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);
    cfg.skip_chromosome_sort = true;

    run_inner(&cfg)?;

    let counts = read_qname_counts(&out_bam)?;
    assert_eq!(
        counts,
        HashMap::from([
            ("frag0_10".to_string(), 2usize),
            ("frag1_30".to_string(), 2usize),
            ("frag2_50".to_string(), 2usize),
        ]),
        "Global mode should keep one complete fragment per chromosome"
    );

    Ok(())
}

#[test]
fn bed_mode_filters_expected_fragments_across_three_chromosomes() -> Result<()> {
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 1_000),
            ("chr2".to_string(), 1_000),
            ("chr3".to_string(), 1_000),
        ],
        vec![
            paired_fragment(10, 120, 40),
            fragment_on_tid(paired_fragment(30, 120, 40), 1),
            fragment_on_tid(paired_fragment(50, 120, 40), 2),
            paired_fragment(400, 120, 40),
            fragment_on_tid(paired_fragment(420, 120, 40), 1),
            fragment_on_tid(paired_fragment(440, 120, 40), 2),
        ],
        Vec::new(),
        "bam_to_bam_three_chr_bed",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("three_chr_bed.bam");
    let bed = work.path().join("three_chr_windows.bed");
    fs::write(&bed, "chr1\t0\t60\nchr2\t0\t80\nchr3\t40\t100\n")?;

    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec![
            "chr1".to_string(),
            "chr2".to_string(),
            "chr3".to_string(),
        ]),
        chromosomes_file: None,
    };
    let mut cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);
    cfg.skip_chromosome_sort = true;
    cfg.by_bed = Some(bed);

    run_inner(&cfg)?;

    let counts = read_qname_counts(&out_bam)?;
    assert_eq!(
        counts,
        HashMap::from([
            ("frag0_10".to_string(), 2usize),
            ("frag1_30".to_string(), 2usize),
            ("frag2_50".to_string(), 2usize),
        ]),
        "BED mode should keep only the fragments overlapping the per-chromosome windows"
    );

    Ok(())
}

#[test]
fn writes_coverage_weight_when_scaling_factors_provided() -> Result<()> {
    let fragment = paired_fragment(0, 100, 40);
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 500)],
        vec![fragment],
        Vec::new(),
        "scaling_weights",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("scaling.bam");
    let scaling = work.path().join("scaling.tsv");
    write_scaling_file(&scaling, "chr1", 500, 2.0)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.scale_genome.scaling_factors = Some(scaling);

    run_inner(&cfg)?;
    let weights = read_tag_values(&out_bam, b"COV")?;
    assert_eq!(
        weights,
        vec![2.0f32, 2.0f32],
        "Both mates should emit the configured scaling factor"
    );

    Ok(())
}

fn base_config(in_bam: &Path, out_bam: &Path) -> BamToBamConfig {
    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec!["chr1".to_string()]),
        chromosomes_file: None,
    };
    let mut cfg = BamToBamConfig::new(in_bam.to_path_buf(), out_bam.to_path_buf(), chrom_args);
    cfg.skip_chromosome_sort = true;
    cfg
}

fn orphan_read(pos: i64) -> ReadSpec {
    ReadSpec {
        tid: 0,
        pos,
        cigar: vec![('M', 40)],
        seq: vec![b'G'; 40],
        qual: 30,
        is_reverse: false,
        mapq: 60,
        flags: 0,
        mate_tid: None,
        mate_pos: None,
        insert_size: 0,
    }
}

fn read_qname_counts(path: &Path) -> Result<HashMap<String, usize>> {
    let mut reader = bam::Reader::from_path(path)?;
    let mut counts: HashMap<String, usize> = HashMap::new();
    for rec in reader.records() {
        let rec = rec?;
        let qname = std::str::from_utf8(rec.qname()).unwrap().to_string();
        *counts.entry(qname).or_default() += 1;
    }
    Ok(counts)
}

fn read_fragment_lengths(path: &Path) -> Result<Vec<u32>> {
    let mut reader = bam::Reader::from_path(path)?;
    let mut values = Vec::new();
    for rec in reader.records() {
        let rec = rec?;
        if let Ok(Aux::U32(value)) = rec.aux(b"FLEN") {
            values.push(value);
        }
    }
    Ok(values)
}

fn read_tag_values(path: &Path, tag: &[u8]) -> Result<Vec<f32>> {
    let mut reader = bam::Reader::from_path(path)?;
    let mut values = Vec::new();
    for rec in reader.records() {
        let rec = rec?;
        if let Ok(Aux::Float(value)) = rec.aux(tag) {
            values.push(value);
        }
    }
    Ok(values)
}

fn fragment_on_tid(mut fragment: FragmentSpec, tid: usize) -> FragmentSpec {
    fragment.forward.tid = tid;
    fragment.reverse.tid = tid;
    fragment.forward.mate_tid = Some(tid);
    fragment.reverse.mate_tid = Some(tid);
    fragment
}

fn first_record_chrom(path: &Path) -> Result<String> {
    let mut reader = bam::Reader::from_path(path)?;
    let header = reader.header().to_owned();
    let rec = reader
        .records()
        .next()
        .expect("BAM should contain records")?;
    let tid = rec.tid() as u32;
    let name = std::str::from_utf8(header.tid2name(tid))
        .unwrap()
        .to_string();
    Ok(name)
}

fn write_bed(path: &Path, windows: &[(u64, u64)]) -> Result<()> {
    let mut contents = String::new();
    for &(start, end) in windows {
        contents.push_str(&format!("chr1\t{}\t{}\n", start, end));
    }
    fs::write(path, contents)?;
    Ok(())
}

fn write_scaling_file(path: &Path, chr: &str, len: u64, factor: f32) -> Result<()> {
    let contents = format!("chromosome\tstart\tend\tscaling_factor\n{chr}\t0\t{len}\t{factor}\n");
    fs::write(path, contents)?;
    Ok(())
}
