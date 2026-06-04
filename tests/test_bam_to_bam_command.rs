#![cfg(feature = "cmd_bam_to_bam")]

// KEEP-IN-TESTS: all active tests in this file cover bam-to-bam command output, errors, or artifacts.

mod fixtures;

use std::{collections::HashMap, fs, path::Path};

use anyhow::{Context, Result};
use cfdnalab::RunOptions;
use cfdnalab::gc_bias::GCCorrectionPackage;
#[cfg(feature = "cmd_coverage_weights")]
use cfdnalab::run_like_cli::coverage_weights::{
    CoverageWeightsConfig, run_coverage_weights as run_coverage_weights_command,
};
#[cfg(feature = "cmd_fragment_count_weights")]
use cfdnalab::run_like_cli::fragment_count_weights::{
    FragmentCountWeightsConfig, run_fragment_count_weights as run_fragment_count_weights_command,
};
use cfdnalab::run_like_cli::{
    bam_to_bam::{BamToBamConfig, BamToBamRunResult, run_bam_to_bam},
    common::{ApplyGCArgFileOnly, ChromosomeArgs},
};
use cfdnalab::testing::{
    Cigar, FragmentSpec, PairedFragmentSpec, ReadSpec, TempBam, bam_from_fragments,
    build_command_produced_gc_correction_package_from_reference_windows,
    single_contig_inward_pair_bam, twobit_from_sequences, twobit_with_single_repeating_contig,
};
use cfdnalab::{
    constants::GC_CORRECTION_SCHEMA_VERSION,
    reference::{ContigFootprintEntry, twobit_contig_footprint},
};
use ndarray::array;
use rust_htslib::bam::{self, Read, record::Aux};
use tempfile::tempdir;

fn run_bam_to_bam_for_test(cfg: &BamToBamConfig) -> Result<BamToBamRunResult> {
    run_bam_to_bam(cfg, RunOptions::new_quiet())
}

#[cfg(feature = "cmd_coverage_weights")]
fn run_coverage_weights(cfg: &CoverageWeightsConfig) -> Result<()> {
    run_coverage_weights_command(cfg, RunOptions::new_quiet()).map(|_| ())
}

#[cfg(feature = "cmd_fragment_count_weights")]
fn run_fragment_count_weights(cfg: &FragmentCountWeightsConfig) -> Result<()> {
    run_fragment_count_weights_command(cfg, RunOptions::new_quiet()).map(|_| ())
}

#[test]
fn filters_on_mapping_quality_and_fragment_membership() -> Result<()> {
    let surviving = PairedFragmentSpec::new(0, 50, 160, 40).build()?;
    let mut low_mapq = PairedFragmentSpec::new(0, 300, 160, 40).build()?;
    low_mapq.forward.mapq = 5;
    low_mapq.reverse.mapq = 30; // One read passes but the entire fragment is removed

    let orphan = orphan_read(500);

    let bam = bam_from_fragments(
        "bam_to_bam_mapq",
        vec![("chr1".to_string(), 1000)],
        vec![surviving, low_mapq],
        vec![orphan],
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("filtered.bam");

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.min_mapq = 30;

    run_bam_to_bam_for_test(&cfg)?;

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
        "Surviving reads must carry the fl AUX tag"
    );

    Ok(())
}

#[test]
fn blacklisting_removes_fragment_when_single_mate_overlaps() -> Result<()> {
    let safe = PairedFragmentSpec::new(0, 10, 120, 40).build()?;
    let with_blacklisted_mate = PairedFragmentSpec::new(0, 200, 120, 40).build()?;

    let bam = bam_from_fragments(
        "bam_to_bam_blacklist",
        vec![("chr1".to_string(), 1000)],
        vec![safe, with_blacklisted_mate],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("blacklisted.bam");
    let blacklist = work.path().join("blacklist.bed");
    fs::write(&blacklist, "chr1\t280\t290\n")?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.blacklist = Some(vec![blacklist]);

    run_bam_to_bam_for_test(&cfg)?;

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
    let safe = PairedFragmentSpec::new(0, 10, 160, 40).build()?;
    let gap_hit = PairedFragmentSpec::new(0, 300, 200, 40).build()?;

    let bam = bam_from_fragments(
        "bam_to_bam_blacklist_gap",
        vec![("chr1".to_string(), 1000)],
        vec![safe, gap_hit],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("blacklisted_gap.bam");
    let blacklist = work.path().join("blacklist_gap.bed");
    fs::write(&blacklist, "chr1\t380\t390\n")?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.blacklist = Some(vec![blacklist]);

    run_bam_to_bam_for_test(&cfg)?;

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
fn by_bed_excludes_chromosomes_without_any_windows() -> Result<()> {
    // Arrange:
    // Build one fragment on each chromosome:
    // - chr1: qname `frag0_20`, span [20, 80)
    // - chr2: qname `frag1_120`, span [120, 180)
    //
    // Then provide a BED file that only contains chr1 [0, 100).
    // In `bam-to-bam`, `--by-bed` is an inclusion filter on fragments overlapping declared BED
    // windows, so chr2 must emit no reads at all.
    let chr1_fragment = PairedFragmentSpec::new(0, 20, 60, 30).build()?;
    let chr2_fragment = fragment_on_tid(PairedFragmentSpec::new(0, 120, 60, 30).build()?, 1);
    let bam = bam_from_fragments(
        "bam_to_bam_missing_chr_bed_windows",
        vec![("chr1".to_string(), 300), ("chr2".to_string(), 300)],
        vec![chr1_fragment, chr2_fragment],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let global_out_bam = work.path().join("global.bam");
    let out_bam = work.path().join("bed_only_chr1.bam");
    let bed = work.path().join("chr1_only.bed");
    fs::write(&bed, "chr1\t0\t100\n")?;

    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
        chromosomes_file: None,
    };
    let mut global_cfg =
        BamToBamConfig::new(bam.bam.clone(), global_out_bam.clone(), chrom_args.clone());
    global_cfg.min_mapq = 0;

    let mut cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);
    cfg.min_mapq = 0;
    cfg.by_bed = Some(bed);

    // Act
    run_bam_to_bam_for_test(&global_cfg)?;
    run_bam_to_bam_for_test(&cfg)?;

    // Assert:
    // First prove the fixture and output path really preserve both chromosomes in global mode.
    // Then show that `--by-bed` removes only the chromosome with no BED rows.
    let global_counts = read_qname_counts(&global_out_bam)?;
    assert_eq!(
        global_counts,
        HashMap::from([
            ("frag0_20".to_string(), 2usize),
            ("frag1_120".to_string(), 2usize),
        ])
    );

    let counts = read_qname_counts(&out_bam)?;
    assert_eq!(counts, HashMap::from([("frag0_20".to_string(), 2usize)]));

    Ok(())
}

#[test]
fn default_min_mapq_matches_explicit_zero_and_differs_from_explicit_thirty() -> Result<()> {
    // Arrange:
    // Build one inward fragment where both mates have MAPQ 20.
    //
    // `bam-to-bam` intentionally defaults to `min_mapq = 0`, so:
    // - default config must keep the fragment
    // - explicit `min_mapq = 0` must do the same
    // - explicit `min_mapq = 30` must remove the fragment entirely
    //
    // Because the command writes both mates of each kept fragment, the expected BAM row counts are:
    // - default / explicit 0: 2 records for qname `frag0_20`
    // - explicit 30:         0 records
    let mut low_mapq = PairedFragmentSpec::new(0, 20, 60, 20).build()?;
    low_mapq.forward.mapq = 20;
    low_mapq.reverse.mapq = 20;
    let bam = bam_from_fragments(
        "bam_to_bam_default_min_mapq",
        vec![("chr1".to_string(), 200)],
        vec![low_mapq],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let default_out = work.path().join("default.bam");
    let explicit_zero_out = work.path().join("explicit_zero.bam");
    let explicit_thirty_out = work.path().join("explicit_thirty.bam");

    let default_cfg = base_config(&bam.bam, &default_out);

    let mut explicit_zero_cfg = base_config(&bam.bam, &explicit_zero_out);
    explicit_zero_cfg.set_min_mapq(0);

    let mut explicit_thirty_cfg = base_config(&bam.bam, &explicit_thirty_out);
    explicit_thirty_cfg.set_min_mapq(30);

    // Act
    run_bam_to_bam_for_test(&default_cfg)?;
    run_bam_to_bam_for_test(&explicit_zero_cfg)?;
    run_bam_to_bam_for_test(&explicit_thirty_cfg)?;

    // Assert
    let default_counts = read_qname_counts(&default_out)?;
    let explicit_zero_counts = read_qname_counts(&explicit_zero_out)?;
    let explicit_thirty_counts = read_qname_counts(&explicit_thirty_out)?;

    assert_eq!(default_counts, explicit_zero_counts);
    assert_eq!(
        default_counts,
        HashMap::from([("frag0_20".to_string(), 2usize)])
    );
    assert!(
        explicit_thirty_counts.is_empty(),
        "raising min_mapq to 30 should remove the MAPQ-20 fragment"
    );

    Ok(())
}

#[test]
fn writes_explicit_chromosomes_in_bam_header_order() -> Result<()> {
    let frag_chr2 = PairedFragmentSpec::new(0, 10, 160, 40).build()?;
    let frag_chr10 = fragment_on_tid(PairedFragmentSpec::new(0, 20, 160, 40).build()?, 1);
    let chroms = vec![("chr2".to_string(), 500), ("chr10".to_string(), 500)];
    let bam = bam_from_fragments(
        "chrom_sort",
        chroms,
        vec![frag_chr2, frag_chr10],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec!["chr10".to_string(), "chr2".to_string()]),
        chromosomes_file: None,
    };

    let out_bam = work.path().join("header_order.bam");
    let cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);
    run_bam_to_bam_for_test(&cfg)?;

    assert_eq!(
        read_record_chromosomes(&out_bam)?,
        vec![
            "chr2".to_string(),
            "chr2".to_string(),
            "chr10".to_string(),
            "chr10".to_string(),
        ],
        "output must follow BAM header target order, not the user-provided chromosome order",
    );
    assert_bam_records_are_coordinate_sorted(&out_bam)?;

    Ok(())
}

#[test]
fn chromosomes_all_follows_bam_header_order() -> Result<()> {
    // Arrange:
    // Build a BAM whose header/contig order is intentionally non-lexicographic:
    //   [chr2, chr10, chr1]
    //
    // With `--chromosomes all`, `bam-to-bam` resolves chromosomes from the BAM header and writes
    // them in that same target-id order so the output remains coordinate-sorted.
    // We place one fragment on each chromosome so the record order reveals that behavior
    // directly. Each kept fragment writes two BAM records, so the expected chromosome sequence is:
    //   [chr2, chr2, chr10, chr10, chr1, chr1]
    let bam = bam_from_fragments(
        "bam_to_bam_all_default_sort",
        vec![
            ("chr2".to_string(), 500),
            ("chr10".to_string(), 500),
            ("chr1".to_string(), 500),
        ],
        vec![
            fragment_on_tid(PairedFragmentSpec::new(0, 10, 120, 40).build()?, 0),
            fragment_on_tid(PairedFragmentSpec::new(0, 20, 120, 40).build()?, 1),
            fragment_on_tid(PairedFragmentSpec::new(0, 30, 120, 40).build()?, 2),
        ],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("all_sorted.bam");
    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec!["all".to_string()]),
        chromosomes_file: None,
    };
    let cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);

    // Act
    run_bam_to_bam_for_test(&cfg)?;

    // Assert
    assert_eq!(
        read_record_chromosomes(&out_bam)?,
        vec![
            "chr2".to_string(),
            "chr2".to_string(),
            "chr10".to_string(),
            "chr10".to_string(),
            "chr1".to_string(),
            "chr1".to_string(),
        ],
        "`--chromosomes all` should follow BAM header target order"
    );
    assert_bam_records_are_coordinate_sorted(&out_bam)?;

    Ok(())
}

#[test]
fn writes_chromosomes_file_selection_in_bam_header_order() -> Result<()> {
    // Arrange:
    // The chromosome file deliberately lists chr1 before chr2, but the BAM header order is chr2,
    // chr10, chr1. The command should keep only the selected subset and write that subset in BAM
    // header order: chr2 first, then chr1. chr10 must not appear because it was not selected.
    let bam = bam_from_fragments(
        "bam_to_bam_chromosome_file_header_order",
        vec![
            ("chr2".to_string(), 500),
            ("chr10".to_string(), 500),
            ("chr1".to_string(), 500),
        ],
        vec![
            fragment_on_tid(PairedFragmentSpec::new(0, 10, 120, 40).build()?, 0),
            fragment_on_tid(PairedFragmentSpec::new(0, 20, 120, 40).build()?, 1),
            fragment_on_tid(PairedFragmentSpec::new(0, 30, 120, 40).build()?, 2),
        ],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("chromosome_file_header_order.bam");
    let chromosomes_file = work.path().join("chromosomes.txt");
    fs::write(&chromosomes_file, "chr1\nchr2\n")?;
    let chrom_args = ChromosomeArgs {
        chromosomes: None,
        chromosomes_file: Some(chromosomes_file),
    };
    let cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);

    // Act
    run_bam_to_bam_for_test(&cfg)?;

    // Assert
    assert_eq!(
        read_record_chromosomes(&out_bam)?,
        vec![
            "chr2".to_string(),
            "chr2".to_string(),
            "chr1".to_string(),
            "chr1".to_string(),
        ],
        "chromosome-file selection should keep the selected subset but write it in BAM header order"
    );
    assert_bam_records_are_coordinate_sorted(&out_bam)?;

    Ok(())
}

#[test]
fn windows_keep_fragment_when_single_read_is_inside() -> Result<()> {
    let target = PairedFragmentSpec::new(0, 20, 200, 40).build()?;
    let outside = PairedFragmentSpec::new(0, 400, 200, 40).build()?;
    let bam = bam_from_fragments(
        "window_single_read",
        vec![("chr1".to_string(), 1000)],
        vec![target, outside],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("window_single.bam");
    let bed = work.path().join("windows_single.bed");
    write_bed4(&bed, &[(0, 60)])?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.by_bed = Some(bed);

    run_bam_to_bam_for_test(&cfg)?;
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
    let target = PairedFragmentSpec::new(0, 0, 200, 40).build()?;
    let outside = PairedFragmentSpec::new(0, 400, 200, 40).build()?;
    let bam = bam_from_fragments(
        "window_gap",
        vec![("chr1".to_string(), 1000)],
        vec![target, outside],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("window_gap.bam");
    let bed = work.path().join("windows_gap.bed");
    write_bed4(&bed, &[(80, 90)])?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.by_bed = Some(bed);

    run_bam_to_bam_for_test(&cfg)?;
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
    let bam = bam_from_fragments(
        "bam_to_bam_three_chr_global",
        vec![
            ("chr1".to_string(), 1_000),
            ("chr2".to_string(), 1_000),
            ("chr3".to_string(), 1_000),
        ],
        vec![
            PairedFragmentSpec::new(0, 10, 120, 40).build()?,
            fragment_on_tid(PairedFragmentSpec::new(0, 30, 120, 40).build()?, 1),
            fragment_on_tid(PairedFragmentSpec::new(0, 50, 120, 40).build()?, 2),
        ],
        Vec::new(),
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
    let cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);

    run_bam_to_bam_for_test(&cfg)?;

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
    let bam = bam_from_fragments(
        "bam_to_bam_three_chr_bed",
        vec![
            ("chr1".to_string(), 1_000),
            ("chr2".to_string(), 1_000),
            ("chr3".to_string(), 1_000),
        ],
        vec![
            PairedFragmentSpec::new(0, 10, 120, 40).build()?,
            fragment_on_tid(PairedFragmentSpec::new(0, 30, 120, 40).build()?, 1),
            fragment_on_tid(PairedFragmentSpec::new(0, 50, 120, 40).build()?, 2),
            PairedFragmentSpec::new(0, 400, 120, 40).build()?,
            fragment_on_tid(PairedFragmentSpec::new(0, 420, 120, 40).build()?, 1),
            fragment_on_tid(PairedFragmentSpec::new(0, 440, 120, 40).build()?, 2),
        ],
        Vec::new(),
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
    cfg.by_bed = Some(bed);

    run_bam_to_bam_for_test(&cfg)?;

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
    let fragment = PairedFragmentSpec::new(0, 0, 100, 40).build()?;
    let bam = bam_from_fragments(
        "scaling_weights",
        vec![("chr1".to_string(), 500)],
        vec![fragment],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("scaling.bam");
    let scaling = work.path().join("scaling.tsv");
    write_scaling_file(&scaling, "chr1", 500, 2.0)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_coverage_scaling_factors(Some(scaling));

    run_bam_to_bam_for_test(&cfg)?;
    let weights = read_tag_values(&out_bam, b"cw")?;
    assert_eq!(
        weights,
        vec![2.0f32, 2.0f32],
        "Both mates should emit the configured scaling factor"
    );

    Ok(())
}

#[test]
fn writes_count_weight_when_count_scaling_factors_provided() -> Result<()> {
    // Arrange:
    // One paired fragment spanning [0, 100) and one chromosome-wide count-scaling factor 0.5.
    //
    // Expected:
    // - both mates receive nw = 0.5
    // - no cw tags are written when only count scaling is configured
    let fragment = PairedFragmentSpec::new(0, 0, 100, 40).build()?;
    let bam = bam_from_fragments(
        "count_scaling_weights",
        vec![("chr1".to_string(), 500)],
        vec![fragment],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("count_scaling.bam");
    let scaling = work.path().join("count_scaling.tsv");
    write_scaling_file(&scaling, "chr1", 500, 0.5)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_count_scaling_factors(Some(scaling));

    // Act
    run_bam_to_bam_for_test(&cfg)?;

    // Assert
    assert_eq!(read_tag_values(&out_bam, b"nw")?, vec![0.5_f32, 0.5_f32]);
    assert!(
        read_tag_values(&out_bam, b"cw")?.is_empty(),
        "cw tags should be absent when only count scaling is configured"
    );

    Ok(())
}

#[test]
fn writes_coverage_and_fragment_count_scaling_to_separate_aux_tags() -> Result<()> {
    // Arrange:
    // Use one simple paired fragment spanning [0, 100) and provide two chromosome-wide scaling
    // files with different constants:
    // - coverage scaling       -> 2.0, expected on tag cw
    // - count-based scaling -> 0.5, expected on tag nw
    //
    // Both mates should carry both tags, so the per-record expectations are:
    // - cw  = 2.0
    // - nw  = 0.5
    // - fl = 100
    let fragment = PairedFragmentSpec::new(0, 0, 100, 40).build()?;
    let bam = bam_from_fragments(
        "dual_scaling_weights",
        vec![("chr1".to_string(), 500)],
        vec![fragment],
        Vec::new(),
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("dual_scaling.bam");
    let coverage_scaling = work.path().join("coverage_scaling.tsv");
    let fragment_count_scaling = work.path().join("fragment_count_scaling.tsv");
    write_scaling_file(&coverage_scaling, "chr1", 500, 2.0)?;
    write_scaling_file(&fragment_count_scaling, "chr1", 500, 0.5)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_coverage_scaling_factors(Some(coverage_scaling));
    cfg.set_count_scaling_factors(Some(fragment_count_scaling));

    // Act
    run_bam_to_bam_for_test(&cfg)?;

    // Assert
    assert_eq!(read_tag_values(&out_bam, b"cw")?, vec![2.0_f32, 2.0_f32]);
    assert_eq!(read_tag_values(&out_bam, b"nw")?, vec![0.5_f32, 0.5_f32]);
    assert_eq!(read_fragment_lengths(&out_bam)?, vec![100_u32, 100_u32]);
    assert_first_record_has_exact_aux_tags(&out_bam, &[b"cw", b"nw", b"fl"])?;
    assert_first_record_lacks_aux_tags(&out_bam, &[b"CO", b"CN", b"FL"])?;

    Ok(())
}

#[test]
fn paired_and_unpaired_fragment_modes_apply_same_full_fragment_scaling_for_same_span() -> Result<()>
{
    // Arrange:
    // We represent the same physical fragment span [20, 80) in two different supported input
    // modes:
    // - paired-end: the `single_contig_inward_pair_bam()` fixture has one fragment with forward [20, 40) and
    //   reverse [60, 80), so fragment span = [20, 80), length = 60
    // - unpaired `reads_are_fragments`: one single read with CIGAR 60M starting at 20, so the
    //   read span is also [20, 80), length = 60
    //
    // `bam-to-bam` computes scaling over the full fragment span in both modes. We therefore use a
    // non-uniform scaling TSV whose hand-derived full-span average is not an integer:
    // - [0, 40)  factor 2.0  -> contributes 20 bp over [20, 40)
    // - [40, 80) factor 1.0  -> contributes 40 bp over [40, 80)
    // - [80,200) factor 1.0  -> not touched
    //
    // Average scaling over [20, 80):
    //   (20 * 2.0 + 40 * 1.0) / 60 = 80 / 60 = 4/3
    //
    // Therefore:
    // - paired mode must emit two records, each tagged with cw = 4/3 and fl = 60
    // - unpaired mode must emit one record tagged with cw = 4/3 and fl = 60
    let paired_bam = single_contig_inward_pair_bam()?;
    let unpaired_bam = bam_from_fragments(
        "bam_to_bam_unpaired_fragment_scaling",
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: 20,
            cigar: vec![Cigar::Match(60)],
            seq: vec![b'A'; 60],
            base_quality: 30,
            is_reverse: false,
            mapq: 40,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
    )?;
    let work = tempdir()?;
    let paired_out = work.path().join("paired_scaled.bam");
    let unpaired_out = work.path().join("unpaired_scaled.bam");
    let scaling = work.path().join("piecewise_scaling.tsv");
    fs::write(
        &scaling,
        "chromosome\tstart\tend\tscaling_factor\n\
chr1\t0\t40\t2.0\n\
chr1\t40\t80\t1.0\n\
chr1\t80\t200\t1.0\n",
    )?;

    let mut paired_cfg = base_config(&paired_bam.bam, &paired_out);
    paired_cfg.set_coverage_scaling_factors(Some(scaling.clone()));
    paired_cfg.set_min_mapq(0);
    {
        let frag = paired_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }

    let mut unpaired_cfg = base_config(&unpaired_bam.bam, &unpaired_out);
    unpaired_cfg.set_coverage_scaling_factors(Some(scaling));
    unpaired_cfg.unpaired.reads_are_fragments = true;
    unpaired_cfg.set_min_mapq(0);
    {
        let frag = unpaired_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }

    // Act
    run_bam_to_bam_for_test(&paired_cfg)?;
    run_bam_to_bam_for_test(&unpaired_cfg)?;

    // Assert
    let paired_cov = read_tag_values(&paired_out, b"cw")?;
    let unpaired_cov = read_tag_values(&unpaired_out, b"cw")?;
    let paired_flen = read_fragment_lengths(&paired_out)?;
    let unpaired_flen = read_fragment_lengths(&unpaired_out)?;
    let expected_cov = 4.0_f32 / 3.0_f32;

    assert_eq!(paired_cov.len(), 2);
    assert_eq!(unpaired_cov.len(), 1);
    for value in paired_cov {
        assert!(
            (value - expected_cov).abs() <= 1e-6,
            "paired output should tag both mates with full-fragment scaling 4/3, got {value}"
        );
    }
    assert!(
        (unpaired_cov[0] - expected_cov).abs() <= 1e-6,
        "unpaired reads_are_fragments output should use the same full-fragment scaling 4/3, got {}",
        unpaired_cov[0]
    );
    assert_eq!(paired_flen, vec![60_u32, 60_u32]);
    assert_eq!(unpaired_flen, vec![60_u32]);

    Ok(())
}

#[test]
fn paired_and_unpaired_fragment_modes_apply_same_full_fragment_count_scaling_for_same_span()
-> Result<()> {
    // Arrange:
    // Reuse the same physical fragment span [20, 80) in paired and unpaired modes and provide
    // a non-uniform count-scaling TSV:
    // - [0, 40)  factor 2.0  -> contributes 20 bp over [20, 40)
    // - [40, 80) factor 1.0  -> contributes 40 bp over [40, 80)
    //
    // Full-fragment average over [20, 80):
    //   (20 * 2.0 + 40 * 1.0) / 60 = 4/3
    //
    // Expected:
    // - paired mode writes nw = 4/3 on both mates
    // - unpaired mode writes nw = 4/3 on the single fragment record
    // - cw is absent in both runs
    let paired_bam = single_contig_inward_pair_bam()?;
    let unpaired_bam = bam_from_fragments(
        "bam_to_bam_unpaired_fragment_count_scaling",
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: 20,
            cigar: vec![Cigar::Match(60)],
            seq: vec![b'A'; 60],
            base_quality: 30,
            is_reverse: false,
            mapq: 40,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
    )?;
    let work = tempdir()?;
    let paired_out = work.path().join("paired_count_scaled.bam");
    let unpaired_out = work.path().join("unpaired_count_scaled.bam");
    let scaling = work.path().join("piecewise_count_scaling.tsv");
    fs::write(
        &scaling,
        "chromosome\tstart\tend\tscaling_factor\n\
chr1\t0\t40\t2.0\n\
chr1\t40\t80\t1.0\n\
chr1\t80\t200\t1.0\n",
    )?;

    let mut paired_cfg = base_config(&paired_bam.bam, &paired_out);
    paired_cfg.set_count_scaling_factors(Some(scaling.clone()));
    paired_cfg.set_min_mapq(0);
    {
        let frag = paired_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }

    let mut unpaired_cfg = base_config(&unpaired_bam.bam, &unpaired_out);
    unpaired_cfg.set_count_scaling_factors(Some(scaling));
    unpaired_cfg.unpaired.reads_are_fragments = true;
    unpaired_cfg.set_min_mapq(0);
    {
        let frag = unpaired_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }

    // Act
    run_bam_to_bam_for_test(&paired_cfg)?;
    run_bam_to_bam_for_test(&unpaired_cfg)?;

    // Assert
    let paired_cnt = read_tag_values(&paired_out, b"nw")?;
    let unpaired_cnt = read_tag_values(&unpaired_out, b"nw")?;
    let paired_flen = read_fragment_lengths(&paired_out)?;
    let unpaired_flen = read_fragment_lengths(&unpaired_out)?;
    // The fragment spans [20,80), so it overlaps:
    // - 20 bp of the [0,40) bin with factor 2.0
    // - 40 bp of the [40,80) bin with factor 1.0
    // Full-fragment average = (20*2 + 40*1) / 60 = 4/3.
    let expected_cnt = 4.0_f32 / 3.0_f32;

    assert_eq!(paired_cnt.len(), 2);
    assert_eq!(unpaired_cnt.len(), 1);
    for value in paired_cnt {
        assert!(
            (value - expected_cnt).abs() <= 1e-6,
            "paired output should tag both mates with full-fragment count scaling 4/3, got {value}"
        );
    }
    assert!(
        (unpaired_cnt[0] - expected_cnt).abs() <= 1e-6,
        "unpaired reads_are_fragments output should use the same full-fragment count scaling 4/3, got {}",
        unpaired_cnt[0]
    );
    assert!(
        read_tag_values(&paired_out, b"cw")?.is_empty(),
        "paired count-scaling-only output should not write cw tags"
    );
    assert!(
        read_tag_values(&unpaired_out, b"cw")?.is_empty(),
        "unpaired count-scaling-only output should not write cw tags"
    );
    assert_eq!(paired_flen, vec![60_u32, 60_u32]);
    assert_eq!(unpaired_flen, vec![60_u32]);

    Ok(())
}

#[test]
fn gc_file_neutralize_invalid_writes_gc_tag_one_on_both_mates() -> Result<()> {
    let bam = single_contig_inward_pair_bam()?;
    let ref_twobit = twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_fallback.bam");
    let gc_path = work.path().join("gc_pkg.zarr");
    build_gc_package(&gc_path, 26, twobit_contig_footprint(&ref_twobit.path)?)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: true,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        // The GC package uses end_offset=26, so the command requires min length > 52.
        frag.min_fragment_length = 53;
        frag.max_fragment_length = 200;
    }

    // Manual expectations:
    // - The fixture contains one paired fragment spanning [20, 80), length 60.
    // - With end_offset=26, the GC-corrected span is only 8 bp long.
    // - The corrector requires at least 10 A/C/G/T bases, so GC lookup fails.
    // - With `neutralize_invalid_gc=true`, `bam-to-bam` keeps the fragment,
    //   `gc_failed_fragments` increments once per fragment,
    //   and both mates receive a `GC` AUX tag of 1.0.
    let counters = run_bam_to_bam_for_test(&cfg)?.counters;

    assert_eq!(counters.base.counted_fragments, 1);
    assert_eq!(counters.gc_failed_fragments, 1);

    let gc_weights = read_tag_values(&out_bam, b"GC")?;
    assert_eq!(
        gc_weights,
        vec![1.0_f32, 1.0_f32],
        "both mates should receive the neutralized GC weight"
    );

    let lengths = read_fragment_lengths(&out_bam)?;
    assert_eq!(
        lengths,
        vec![60_u32, 60_u32],
        "the fragment should still be emitted with its fl tags"
    );

    Ok(())
}

#[test]
fn gc_file_default_behavior_skips_fragment_entirely() -> Result<()> {
    let bam = single_contig_inward_pair_bam()?;
    let ref_twobit = twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_drop_invalid.bam");
    let gc_path = work.path().join("gc_pkg.zarr");
    build_gc_package(&gc_path, 26, twobit_contig_footprint(&ref_twobit.path)?)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        // The only fixture fragment is length 60, so with end_offset=26 the corrected span is
        // only 8 bp. That is below the corrector's minimum 10 A/C/G/T requirement, so the GC
        // lookup fails for the one logical fragment in this BAM.
        frag.min_fragment_length = 53;
        frag.max_fragment_length = 200;
    }

    // Manual expectations:
    // - one fragment is encountered
    // - GC lookup fails once for that fragment
    // - with the default GC behavior, the command must skip the fragment instead of silently
    //   neutralizing it to weight 1.0
    // - therefore the output BAM must contain no records and no AUX tags at all
    let counters = run_bam_to_bam_for_test(&cfg)?.counters;

    assert_eq!(counters.gc_failed_fragments, 1);
    assert_eq!(counters.base.counted_fragments, 0);
    assert!(
        read_qname_counts(&out_bam)?.is_empty(),
        "dropping the only invalid-GC fragment should leave an empty BAM"
    );
    assert!(
        read_tag_values(&out_bam, b"GC")?.is_empty(),
        "no GC tags should be written when every fragment is skipped"
    );
    assert!(
        read_fragment_lengths(&out_bam)?.is_empty(),
        "no fl tags should be written when every fragment is skipped"
    );

    Ok(())
}

#[test]
fn gc_file_and_scaling_factors_write_identical_gc_cov_and_flen_tags_on_both_mates() -> Result<()> {
    let bam = single_contig_inward_pair_bam()?;
    let ref_twobit = twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_and_cov.bam");
    let gc_path = work.path().join("gc_pkg.zarr");
    let scaling_path = work.path().join("scaling.tsv");
    build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;
    write_scaling_file(&scaling_path, "chr1", 200, 4.0_f32 / 3.0_f32)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_coverage_scaling_factors(Some(scaling_path));
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    // Manual expectations:
    // - `single_contig_inward_pair_bam()` contains one fragment spanning [20,80), so fl must be 60
    // - the whole-chrom scaling TSV sets factor 4/3 everywhere, so both mates must receive
    //   identical `cw = 4/3`
    // - the helper GC package has:
    //     length bins [10,60) and [60,200]
    //     GC bins     [0,50) and [50,101]
    //     correction matrix row for length 60 = [2, 10]
    // - the repeated ACGT reference gives GC%=50 for any 60 bp fragment, so the fragment lands in
    //   the second GC bin and both mates must receive `GC = 10`
    // - because `bam-to-bam` tags per fragment, both output records must carry the same
    //   `GC`, `cw`, and `fl` values
    let counters = run_bam_to_bam_for_test(&cfg)?.counters;

    assert_eq!(counters.base.counted_fragments, 1);
    assert_eq!(read_tag_values(&out_bam, b"GC")?, vec![10.0_f32, 10.0_f32]);
    assert_eq!(
        read_tag_values(&out_bam, b"cw")?,
        vec![4.0_f32 / 3.0_f32, 4.0_f32 / 3.0_f32]
    );
    assert_eq!(read_fragment_lengths(&out_bam)?, vec![60_u32, 60_u32]);

    Ok(())
}

#[test]
fn gc_file_and_count_scaling_factors_write_identical_gc_cnt_and_flen_tags_on_both_mates()
-> Result<()> {
    let bam = single_contig_inward_pair_bam()?;
    let ref_twobit = twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_and_cnt.bam");
    let gc_path = work.path().join("gc_pkg.zarr");
    let scaling_path = work.path().join("count_scaling.tsv");
    build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;
    write_scaling_file(&scaling_path, "chr1", 200, 4.0_f32 / 3.0_f32)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_count_scaling_factors(Some(scaling_path));
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    // Manual expectations:
    // - same fragment and GC derivation as the matching coverage-based test
    // - both mates must receive GC = 10 and nw = 4/3
    // - cw stays absent when only count scaling is configured
    let counters = run_bam_to_bam_for_test(&cfg)?.counters;

    assert_eq!(counters.base.counted_fragments, 1);
    assert_eq!(read_tag_values(&out_bam, b"GC")?, vec![10.0_f32, 10.0_f32]);
    assert_eq!(
        read_tag_values(&out_bam, b"nw")?,
        vec![4.0_f32 / 3.0_f32, 4.0_f32 / 3.0_f32]
    );
    assert!(
        read_tag_values(&out_bam, b"cw")?.is_empty(),
        "cw tags should be absent when only count scaling is configured"
    );
    assert_eq!(read_fragment_lengths(&out_bam)?, vec![60_u32, 60_u32]);

    Ok(())
}

#[cfg(feature = "cmd_fragment_count_weights")]
#[test]
fn real_fragment_count_weights_tsv_is_applied_per_fragment_in_bam_to_bam() -> Result<()> {
    // Arrange:
    // Use the same mixed-length fixture as the normalize-genome test:
    // - short fragment [0,20)
    // - long  fragment [20,80)
    //
    // With bin_size == stride == 20, fragment-count-weights yields:
    // - short covered bin scaling factor = 0.5
    // - long covered bins scaling factor = 1.5
    //
    // `bam-to-bam` averages scaling over the full fragment span, so:
    // - short fragment gets nw = 0.5
    // - long fragment gets nw = 1.5
    let bam = mixed_length_fragment_bam("bam_to_bam_real_count_weights")?;
    let work = tempdir()?;
    let weights_out_dir = work.path().join("weights_out");
    fs::create_dir_all(&weights_out_dir)?;

    let mut weights_cfg = FragmentCountWeightsConfig::new(
        cfdnalab::run_like_cli::common::IOCArgs {
            bam: bam.bam.clone(),
            output_dir: weights_out_dir.clone(),
            n_threads: 2,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    );
    weights_cfg.set_bin_size(20);
    weights_cfg.set_stride(20);
    weights_cfg.set_min_mapq(0);
    weights_cfg.set_require_proper_pair(false);
    weights_cfg.set_output_prefix("counts".to_string());
    {
        let frag = weights_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }
    run_fragment_count_weights(&weights_cfg)?;

    let scaling_path = weights_out_dir.join("counts.fragment_counts.scaling_factors.tsv");
    let out_bam = work.path().join("count_scaled_real.bam");
    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_count_scaling_factors(Some(scaling_path));
    cfg.set_min_mapq(0);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }

    // Act
    run_bam_to_bam_for_test(&cfg)?;

    // Assert
    let mut reader = bam::Reader::from_path(&out_bam)?;
    let mut observed = Vec::new();
    for record in reader.records() {
        let record = record?;
        let cnt = match record.aux(b"nw") {
            Ok(Aux::Float(value)) => value,
            other => panic!("expected nw float tag on every read, got {other:?}"),
        };
        let flen = match record.aux(b"fl") {
            Ok(Aux::U32(value)) => value,
            other => panic!("expected fl u32 tag on every read, got {other:?}"),
        };
        observed.push((flen, cnt));
    }
    observed.sort_by(|left, right| left.0.cmp(&right.0));

    assert_eq!(observed.len(), 4);
    for (idx, (flen, cnt)) in observed.iter().enumerate() {
        let (expected_flen, expected_cnt) = if idx < 2 {
            (20_u32, 0.5_f32)
        } else {
            (60_u32, 1.5_f32)
        };
        assert_eq!(*flen, expected_flen);
        assert!(
            (*cnt - expected_cnt).abs() <= 1e-6,
            "unexpected nw for fl {expected_flen}: expected {expected_cnt}, got {cnt}"
        );
    }
    assert!(
        read_tag_values(&out_bam, b"cw")?.is_empty(),
        "cw tags should be absent when consuming fragment-count weights"
    );

    Ok(())
}

#[test]
fn gc_file_rejects_package_when_fragment_length_range_is_outside_supported_range() -> Result<()> {
    // Arrange:
    // The fixture contributes one fragment of length 60. We configure the converter to accept
    // exactly that length, then write the smallest valid GC package that only covers 10..=59.
    //
    // The converter validates the GC package before it starts chromosome iteration, so the
    // expected error is the shared compatibility failure from the GC loader.
    let bam = single_contig_inward_pair_bam()?;
    let ref_twobit = twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_range_error.bam");
    let gc_path = work.path().join("gc_pkg_short.zarr");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![10, 59],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
        correction_matrix: array![[1.0_f64]],
    };
    package.write_zarr(&gc_path)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    // Act
    let err = run_bam_to_bam_for_test(&cfg).expect_err("out-of-range GC package should fail");

    // Assert
    let msg = err.to_string();
    assert!(
        msg.contains("fragment length range [60-60] is outside the range covered by the correction package [10-59]"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn real_ref_gc_bias_then_gc_bias_package_changes_bam_to_bam_in_expected_direction() -> Result<()> {
    // Arrange:
    // Use the same real non-neutral producer workflow as the corresponding `gc-bias` test:
    // - reference: chr1[0,100) all A, chr1[100,200) all C
    // - reference-side windows: [0,91) and [100,191)
    // - with the `ref-gc-bias` fit rule those count starts 0..=81 and 100..=181 only, so the
    //   written reference support is still balanced between pure-A and pure-C starts
    // - sample BAM: one A-only 10 bp fragment and nine C-only 10 bp fragments
    //
    // The resulting real package is hand-derived as:
    // - GC%=0   -> weight 5.0
    // - GC%=100 -> weight 5/9
    //
    // `bam-to-bam` writes one tagged BAM record per input read, so the expected output is:
    // - first mate pair tagged GC=5.0, fl=10
    // - remaining nine mate pairs tagged GC=5/9, fl=10
    let reference = twobit_from_sequences(
        "bam_to_bam_real_non_neutral_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;
    let starts = [10_i64, 110, 120, 130, 140, 150, 160, 170, 180, 190];
    let fragments = starts
        .into_iter()
        .map(|start| PairedFragmentSpec::new(0, start, 10, 5).build())
        .collect::<Result<Vec<_>>>()?;
    let bam = bam_from_fragments(
        "bam_to_bam_real_non_neutral_bam",
        vec![("chr1".to_string(), 200)],
        fragments,
        Vec::new(),
    )?;
    let work = tempdir()?;
    let out_bam = work.path().join("real_non_neutral_gc.bam");
    let gc_path = build_command_produced_gc_correction_package_from_reference_windows(
        &bam.bam,
        &reference.path,
        work.path(),
        10,
        "chr1\t0\t91\nchr1\t100\t191\n",
        // Chromosome length 200 and fragment length 10 give:
        //   200 - 10 + 1 = 191 valid starts.
        191,
    )?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 10;
    }

    // Act
    let counters = run_bam_to_bam_for_test(&cfg)?.counters;

    // Assert
    assert_eq!(counters.base.counted_fragments, 10);
    let gc_tags = read_tag_values(&out_bam, b"GC")?;
    let flens = read_fragment_lengths(&out_bam)?;
    assert_eq!(gc_tags.len(), 20);
    assert_eq!(flens, vec![10_u32; 20]);
    assert_eq!(gc_tags[0], 5.0_f32);
    assert_eq!(gc_tags[1], 5.0_f32);
    for value in gc_tags.iter().skip(2) {
        assert!(
            (*value as f64 - (5.0 / 9.0)).abs() <= 1e-6,
            "expected GC tag 5/9 on C-only fragments, got {value}"
        );
    }

    Ok(())
}

#[test]
fn bed_blacklist_scaling_and_gc_together_keep_only_the_expected_tagged_fragment() -> Result<()> {
    let reference = twobit_from_sequences(
        "bam_to_bam_combined_filters_reference",
        vec![("chr1".to_string(), "ACGT".repeat(75))],
    )?;
    let bam = bam_from_fragments(
        "bam_to_bam_combined_filters_bam",
        vec![("chr1".to_string(), 300)],
        vec![
            PairedFragmentSpec::new(0, 20, 60, 20).build()?,
            PairedFragmentSpec::new(0, 100, 60, 20).build()?,
            PairedFragmentSpec::new(0, 220, 60, 20).build()?,
        ],
        Vec::new(),
    )?;
    let work = tempdir()?;
    let out_bam = work.path().join("combined_filters.bam");
    let gc_path = work.path().join("gc_pkg.zarr");
    let scaling_path = work.path().join("scaling.tsv");
    let bed_path = work.path().join("windows.bed");
    let blacklist_path = work.path().join("blacklist.bed");
    build_gc_package(&gc_path, 0, twobit_contig_footprint(&reference.path)?)?;
    write_scaling_file(&scaling_path, "chr1", 300, 4.0_f32 / 3.0_f32)?;
    write_bed4(&bed_path, &[(0, 180)])?;
    fs::write(&blacklist_path, "chr1\t120\t130\n")?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_coverage_scaling_factors(Some(scaling_path));
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.by_bed = Some(bed_path);
    cfg.blacklist = Some(vec![blacklist_path]);
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    // Manual expectations:
    // - fragment A spans [20,80): inside BED, outside blacklist -> kept
    // - fragment B spans [100,160): inside BED, but overlaps blacklist [120,130) -> removed
    // - fragment C spans [220,280): outside BED -> removed before tagging
    // - the kept fragment still sees the same repeated-ACGT reference semantics as above:
    //     GC%=50 -> `GC = 10`
    //     scaling factor is uniform -> `cw = 4/3`
    //     fragment length -> `fl = 60`
    let counters = run_bam_to_bam_for_test(&cfg)?.counters;

    assert_eq!(counters.base.counted_fragments, 1);
    assert_eq!(counters.blacklisted_fragments, 1);
    assert_eq!(
        read_qname_counts(&out_bam)?,
        HashMap::from([(String::from("frag0_20"), 2_usize)])
    );
    assert_eq!(read_tag_values(&out_bam, b"GC")?, vec![10.0_f32, 10.0_f32]);
    assert_eq!(
        read_tag_values(&out_bam, b"cw")?,
        vec![4.0_f32 / 3.0_f32, 4.0_f32 / 3.0_f32]
    );
    assert_eq!(read_fragment_lengths(&out_bam)?, vec![60_u32, 60_u32]);

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn real_multi_chromosome_coverage_weights_tsv_is_applied_per_chromosome_in_bam_to_bam() -> Result<()>
{
    // Arrange:
    // Reuse the same two-chromosome fixture and hand derivation as the matching `bam-to-frag`
    // workflow and the command-level `coverage-weights` shared-global-mean test.
    //
    // chr1 fragment:
    // - span [20, 80), length 60
    // - avg-overlap profile:
    //   [1/3, 3/4, 1, 3/4, 1/4, 0, ...]
    //
    // chr2 fragment:
    // - span [20, 40), length 20
    // - avg-overlap profile:
    //   [1/3, 1/2, 1/4, 0, ...]
    //
    // Shared non-zero global mean:
    //   chr1 sum = 37/12
    //   chr2 sum = 13/12
    //   total    = 25/6 across 8 non-zero bins
    //   mean     = (25/6) / 8 = 25/48
    //
    // Inverted per-bin scaling factors are therefore:
    //   chr1 [20,40) = (25/48) / (3/4) = 25/36
    //   chr1 [40,60) = (25/48) / 1     = 25/48
    //   chr1 [60,80) = (25/48) / (3/4) = 25/36
    //   chr2 [20,40) = (25/48) / (1/2) = 25/24
    //
    // `bam-to-bam` averages scaling over the full fragment span:
    //   chr1 weight = ((25/36) + (25/48) + (25/36)) / 3 = 275/432
    //   chr2 weight = 25/24
    //
    // It writes one tag set per emitted read, so the final BAM must contain:
    // - two chr1 reads with cw = 275/432 and fl = 60
    // - two chr2 reads with cw = 25/24   and fl = 20
    let mut chr2_fragment = PairedFragmentSpec::new(0, 20, 20, 10).build()?;
    chr2_fragment = fragment_on_tid(chr2_fragment, 1);

    let bam = bam_from_fragments(
        "bam_to_bam_real_multi_chr_scaling",
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![
            PairedFragmentSpec::new(0, 20, 60, 20).build()?,
            chr2_fragment,
        ],
        Vec::new(),
    )?;
    let work = tempdir()?;

    let weights_out_dir = work.path().join("weights_out");
    std::fs::create_dir_all(&weights_out_dir)?;
    let mut weights_cfg = CoverageWeightsConfig::new(
        cfdnalab::run_like_cli::common::IOCArgs {
            bam: bam.bam.clone(),
            output_dir: weights_out_dir.clone(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
            chromosomes_file: None,
        },
    );
    weights_cfg.set_output_prefix("coverage".to_string());
    weights_cfg.set_bin_size(40);
    weights_cfg.set_stride(20);
    weights_cfg.set_min_mapq(0);
    weights_cfg.set_require_proper_pair(false);
    {
        let frag = weights_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }
    run_coverage_weights(&weights_cfg)?;

    let scaling_path = weights_out_dir.join("coverage.coverage.scaling_factors.tsv");
    let out_bam = work.path().join("scaled_multi_chr.bam");
    let mut cfg = BamToBamConfig::new(
        bam.bam.clone(),
        out_bam.clone(),
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
            chromosomes_file: None,
        },
    );
    cfg.set_coverage_scaling_factors(Some(scaling_path));
    cfg.min_mapq = 0;
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    // Act
    let counters = run_bam_to_bam_for_test(&cfg)?.counters;

    // Assert
    assert_eq!(counters.base.counted_fragments, 2);
    let mut reader = bam::Reader::from_path(&out_bam)?;
    let header = reader.header().to_owned();
    let mut observed: Vec<(String, f32, u32)> = Vec::new();
    for record_result in reader.records() {
        let record = record_result?;
        let tid = record.tid() as u32;
        let chromosome = std::str::from_utf8(header.tid2name(tid))?.to_string();
        let scaling = match record.aux(b"cw") {
            Ok(Aux::Float(value)) => value,
            other => panic!("expected cw float tag on every read, got {other:?}"),
        };
        let flen = match record.aux(b"fl") {
            Ok(Aux::U32(value)) => value,
            other => panic!("expected fl u32 tag on every read, got {other:?}"),
        };
        observed.push((chromosome, scaling, flen));
    }
    observed.sort_by(|left, right| left.0.cmp(&right.0));

    let expected_chr1 = 275.0_f32 / 432.0_f32;
    let expected_chr2 = 25.0_f32 / 24.0_f32;

    assert_eq!(observed.len(), 4);
    for (idx, (chromosome, scaling, flen)) in observed.iter().enumerate() {
        let (expected_chrom, expected_scaling, expected_flen) = if idx < 2 {
            ("chr1", expected_chr1, 60_u32)
        } else {
            ("chr2", expected_chr2, 20_u32)
        };
        assert_eq!(chromosome, expected_chrom);
        assert!(
            (*scaling - expected_scaling).abs() <= 1e-6,
            "unexpected scaling for {expected_chrom}: expected {expected_scaling}, got {scaling}"
        );
        assert_eq!(
            *flen, expected_flen,
            "unexpected fl for {expected_chrom}: expected {expected_flen}, got {flen}"
        );
    }

    Ok(())
}

#[test]
fn gc_file_rejects_package_with_schema_version_mismatch() -> Result<()> {
    // Arrange:
    // Build the smallest valid GC correction package shape, but with an intentionally
    // incompatible schema version. `bam-to-bam` should reject it while loading the package,
    // before chromosome processing starts or any output BAM records are written.
    let bam = single_contig_inward_pair_bam()?;
    let ref_twobit = twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_bad_version.bam");
    let gc_path = work.path().join("gc_pkg_bad_version.zarr");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION + 1,
        end_offset: 0,
        length_edges: vec![10, 200],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&ref_twobit.path)?,
        correction_matrix: array![[1.0_f64]],
    };
    package.write_zarr(&gc_path)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

    // Act
    let err = run_bam_to_bam_for_test(&cfg).expect_err("schema version mismatch should fail");

    // Assert
    let msg = err.to_string();
    assert!(
        msg.contains("GC correction package schema version mismatch"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn scaling_tsv_must_cover_requested_chromosome_end_in_bam_to_bam() -> Result<()> {
    // Arrange:
    // `single_contig_inward_pair_bam()` uses chr1 length 200.
    // The shared scaling loader requires bins to cover the whole requested chromosome exactly.
    //
    // This TSV stops at 100:
    //   [0,100) factor 2.0
    // so the converter must reject it before chromosome processing begins.
    let bam = single_contig_inward_pair_bam()?;
    let work = tempdir()?;
    let out_bam = work.path().join("scaling_truncated.bam");
    let scaling_path = work.path().join("truncated_scaling.tsv");
    write_scaling_file(&scaling_path, "chr1", 100, 2.0)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_coverage_scaling_factors(Some(scaling_path));

    // Act
    let err = run_bam_to_bam_for_test(&cfg).expect_err("truncated scaling TSV should fail");

    // Assert:
    // `bam-to-bam` wraps the shared loader with `load scaling factors`, so inspect the full
    // error chain rather than only the top-level context.
    let msg = format!("{err:#}");
    assert!(
        msg.contains("scaling TSV: bins on 'chr1' must end at chrom_len=200 (got end=100)"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn count_scaling_tsv_must_cover_requested_chromosome_end_in_bam_to_bam() -> Result<()> {
    // Arrange:
    // This mirrors the coverage-side truncated TSV regression, but through the separate
    // `--count-scaling-factors` loader path. The artifact contract is the same: bins must cover
    // the full requested chromosome even if the counted fragment lies entirely inside [0,100).
    let bam = single_contig_inward_pair_bam()?;
    let work = tempdir()?;
    let out_bam = work.path().join("count_scaling_truncated.bam");
    let scaling_path = work.path().join("truncated_count_scaling.tsv");
    write_scaling_file(&scaling_path, "chr1", 100, 2.0)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_count_scaling_factors(Some(scaling_path));

    // Act
    let err = run_bam_to_bam_for_test(&cfg).expect_err("truncated count scaling TSV should fail");

    // Assert
    let msg = format!("{err:#}");
    assert!(
        msg.contains("scaling TSV: bins on 'chr1' must end at chrom_len=200 (got end=100)"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

#[test]
fn count_scaling_tsv_with_uncorrected_metadata_rejects_gc_corrected_bam_to_bam_run() -> Result<()> {
    // Arrange:
    // The BAM run uses `--gc-file`, so the current run is explicitly GC-corrected.
    // A count-scaling TSV marked `gc_mode=uncorrected` is therefore a known mismatch and should
    // fail during scaling-map loading before any BAM records are written.
    let bam = single_contig_inward_pair_bam()?;
    let ref_twobit = twobit_with_single_repeating_contig("simple_reference", "chr1", "ACGT", 256)?;
    let work = tempdir()?;
    let out_bam = work.path().join("count_scaling_gc_mismatch.bam");
    let scaling_path = work.path().join("uncorrected_count_scaling.tsv");
    let gc_path = work.path().join("gc_pkg.zarr");
    build_gc_package(&gc_path, 0, twobit_contig_footprint(&ref_twobit.path)?)?;
    fs::write(
        &scaling_path,
        "# gc_mode=uncorrected\nchromosome\tstart\tend\tscaling_factor\nchr1\t0\t200\t1.0\n",
    )?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_count_scaling_factors(Some(scaling_path));
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

    // Act
    let err = run_bam_to_bam_for_test(&cfg)
        .expect_err("uncorrected count scaling must fail in a GC-corrected bam-to-bam run");

    // Assert
    let msg = format!("{err:#}");
    assert!(
        msg.contains("no GC correction"),
        "unexpected error message: {msg}"
    );
    assert!(
        msg.contains("via --gc-file"),
        "unexpected error message: {msg}"
    );

    Ok(())
}

fn base_config(in_bam: &Path, out_bam: &Path) -> BamToBamConfig {
    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec!["chr1".to_string()]),
        chromosomes_file: None,
    };
    BamToBamConfig::new(in_bam.to_path_buf(), out_bam.to_path_buf(), chrom_args)
}

fn orphan_read(pos: i64) -> ReadSpec {
    ReadSpec {
        tid: 0,
        pos,
        cigar: vec![Cigar::Match(40)],
        seq: vec![b'G'; 40],
        base_quality: 30,
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
        if let Ok(Aux::U32(value)) = rec.aux(b"fl") {
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

fn read_record_chromosomes(path: &Path) -> Result<Vec<String>> {
    let mut reader = bam::Reader::from_path(path)?;
    let header = reader.header().to_owned();
    let mut chromosomes = Vec::new();
    for record_result in reader.records() {
        let record = record_result?;
        let tid = record.tid() as u32;
        chromosomes.push(std::str::from_utf8(header.tid2name(tid))?.to_string());
    }
    Ok(chromosomes)
}

fn assert_bam_records_are_coordinate_sorted(path: &Path) -> Result<()> {
    let mut reader = bam::Reader::from_path(path)?;
    let mut previous_coordinate: Option<(i32, i64)> = None;
    for record_result in reader.records() {
        let record = record_result?;
        let current_coordinate = (record.tid(), record.pos());
        if let Some(previous_coordinate) = previous_coordinate {
            assert!(
                previous_coordinate <= current_coordinate,
                "BAM records in {} are not coordinate-sorted: observed {:?} after {:?}",
                path.display(),
                current_coordinate,
                previous_coordinate
            );
        }
        previous_coordinate = Some(current_coordinate);
    }
    Ok(())
}

fn assert_first_record_has_exact_aux_tags(path: &Path, expected_tags: &[&[u8; 2]]) -> Result<()> {
    let aux_tags = read_first_record_aux_tag_names(path)?;
    for expected_tag in expected_tags {
        assert!(
            aux_tags
                .iter()
                .any(|observed_tag| observed_tag.as_slice() == expected_tag.as_slice()),
            "expected first record in {} to contain exact AUX tag {:?}, observed {:?}",
            path.display(),
            std::str::from_utf8(expected_tag.as_slice()).unwrap(),
            aux_tags
        );
    }
    Ok(())
}

fn assert_first_record_lacks_aux_tags(path: &Path, unexpected_tags: &[&[u8; 2]]) -> Result<()> {
    let aux_tags = read_first_record_aux_tag_names(path)?;
    for unexpected_tag in unexpected_tags {
        assert!(
            !aux_tags
                .iter()
                .any(|observed_tag| observed_tag.as_slice() == unexpected_tag.as_slice()),
            "first record in {} should not contain old truncated AUX tag {:?}, observed {:?}",
            path.display(),
            std::str::from_utf8(unexpected_tag.as_slice()).unwrap(),
            aux_tags
        );
    }
    Ok(())
}

fn read_first_record_aux_tag_names(path: &Path) -> Result<Vec<Vec<u8>>> {
    let mut reader = bam::Reader::from_path(path)?;
    let record = reader
        .records()
        .next()
        .context("expected BAM to contain at least one record")??;
    record
        .aux_iter()
        .map(|aux_result| {
            aux_result
                .map(|(tag, _aux_value)| tag.to_vec())
                .map_err(Into::into)
        })
        .collect()
}

fn fragment_on_tid(mut fragment: FragmentSpec, tid: usize) -> FragmentSpec {
    fragment.forward.tid = tid;
    fragment.reverse.tid = tid;
    fragment.forward.mate_tid = Some(tid);
    fragment.reverse.mate_tid = Some(tid);
    fragment
}

#[cfg(feature = "cmd_fragment_count_weights")]
fn mixed_length_fragment_bam(name: &str) -> Result<TempBam> {
    bam_from_fragments(
        name,
        vec![("chr1".to_string(), 100)],
        vec![
            PairedFragmentSpec::new(0, 0, 20, 10).build()?,
            PairedFragmentSpec::new(0, 20, 60, 10).build()?,
        ],
        Vec::new(),
    )
}

fn write_bed4(path: &Path, windows: &[(u64, u64)]) -> Result<()> {
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

fn build_gc_package(
    path: &Path,
    end_offset: u64,
    reference_contig_footprint: Vec<ContigFootprintEntry>,
) -> Result<()> {
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset,
        length_edges: vec![10, 60, 200],
        gc_edges: vec![0, 50, 101],
        length_bin_frequencies: array![1.0_f64, 3.0_f64],
        reference_contig_footprint,
        correction_matrix: array![[1.0_f64, 1.0_f64], [2.0_f64, 10.0_f64]],
    };
    package.write_zarr(path)?;
    Ok(())
}
