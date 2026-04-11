#![cfg(feature = "cmd_bam_to_bam")]

mod fixtures;

use std::{collections::HashMap, fs, path::Path};

use anyhow::Result;
#[cfg(feature = "cmd_coverage_weights")]
use cfdnalab::commands::coverage_weights::{
    config::CoverageWeightsConfig, coverage_weights::run as run_coverage_weights,
};
use cfdnalab::commands::{
    bam_to_bam::{bam_to_bam::run_inner, config::BamToBamConfig},
    cli_common::{ApplyGCArgFileOnly, ChromosomeArgs},
    gc_bias::{GC_CORRECTION_SCHEMA_VERSION, package::GCCorrectionPackage},
};
use fixtures::{
    FragmentSpec, ReadSpec, bam_from_specs, build_real_non_neutral_gc_package, paired_fragment,
    simple_inward_bam, simple_reference_twobit, twobit_from_sequences,
};
use ndarray::array;
use rust_htslib::bam::{self, Read, record::Aux};
use tempfile::tempdir;

#[test]
fn filters_on_mapping_quality_and_fragment_membership() -> Result<()> {
    // Human verification status: unverified
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
    // Human verification status: unverified
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
    // Human verification status: unverified
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
fn by_bed_excludes_chromosomes_without_any_windows() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Build one fragment on each chromosome:
    // - chr1: qname `frag0_20`, span [20, 80)
    // - chr2: qname `frag1_120`, span [120, 180)
    //
    // Then provide a BED file that only contains chr1 [0, 100).
    // In `bam-to-bam`, `--by-bed` is an inclusion filter on fragments overlapping declared BED
    // windows, so chr2 must emit no reads at all.
    let chr1_fragment = paired_fragment(20, 60, 30);
    let chr2_fragment = fragment_on_tid(paired_fragment(120, 60, 30), 1);
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 300), ("chr2".to_string(), 300)],
        vec![chr1_fragment, chr2_fragment],
        Vec::new(),
        "bam_to_bam_missing_chr_bed_windows",
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
    global_cfg.skip_chromosome_sort = true;
    global_cfg.min_mapq = 0;

    let mut cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);
    cfg.skip_chromosome_sort = true;
    cfg.min_mapq = 0;
    cfg.by_bed = Some(bed);

    // Act
    run_inner(&global_cfg)?;
    run_inner(&cfg)?;

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
    // Human verification status: unverified
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
    let mut low_mapq = paired_fragment(20, 60, 20);
    low_mapq.forward.mapq = 20;
    low_mapq.reverse.mapq = 20;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![low_mapq],
        Vec::new(),
        "bam_to_bam_default_min_mapq",
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
    run_inner(&default_cfg)?;
    run_inner(&explicit_zero_cfg)?;
    run_inner(&explicit_thirty_cfg)?;

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
fn respects_chromosome_sort_toggle() -> Result<()> {
    // Human verification status: unverified
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
fn chromosomes_all_default_sorts_lexicographically_instead_of_bam_header_order() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Build a BAM whose header/contig order is intentionally non-lexicographic:
    //   [chr2, chr10, chr1]
    //
    // With `--chromosomes all`, `bam-to-bam` first resolves chromosomes from the BAM header and
    // then, by default, sorts them lexicographically unless `--skip-chromosome-sort` is set.
    //
    // Therefore the output processing order must be:
    //   chr1, chr10, chr2
    //
    // We place one fragment on each chromosome so the emitted record order reveals that behavior
    // directly. Each kept fragment writes two BAM records, so the expected chromosome sequence is:
    //   [chr1, chr1, chr10, chr10, chr2, chr2]
    let bam = bam_from_specs(
        vec![
            ("chr2".to_string(), 500),
            ("chr10".to_string(), 500),
            ("chr1".to_string(), 500),
        ],
        vec![
            fragment_on_tid(paired_fragment(10, 120, 40), 0),
            fragment_on_tid(paired_fragment(20, 120, 40), 1),
            fragment_on_tid(paired_fragment(30, 120, 40), 2),
        ],
        Vec::new(),
        "bam_to_bam_all_default_sort",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("all_sorted.bam");
    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec!["all".to_string()]),
        chromosomes_file: None,
    };
    let cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);

    // Act
    run_inner(&cfg)?;

    // Assert
    let mut reader = bam::Reader::from_path(&out_bam)?;
    let header = reader.header().to_owned();
    let mut record_chromosomes = Vec::new();
    for record_result in reader.records() {
        let record = record_result?;
        let tid = record.tid() as u32;
        record_chromosomes.push(
            std::str::from_utf8(header.tid2name(tid))
                .unwrap()
                .to_string(),
        );
    }
    assert_eq!(
        record_chromosomes,
        vec![
            "chr1".to_string(),
            "chr1".to_string(),
            "chr10".to_string(),
            "chr10".to_string(),
            "chr2".to_string(),
            "chr2".to_string(),
        ],
        "`--chromosomes all` should still follow the command's default lexicographic sorting"
    );

    Ok(())
}

#[test]
fn chromosomes_all_with_skip_sort_follows_bam_header_order() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Reuse the same intentionally non-lexicographic BAM header/contig order:
    //   [chr2, chr10, chr1]
    //
    // With `--chromosomes all`, chromosome resolution starts from that BAM header order.
    // When `--skip-chromosome-sort` is enabled, the command should keep that resolved order
    // instead of applying its default lexicographic reordering.
    //
    // With one fragment on each chromosome, the emitted chromosome sequence must therefore be:
    //   [chr2, chr2, chr10, chr10, chr1, chr1]
    let bam = bam_from_specs(
        vec![
            ("chr2".to_string(), 500),
            ("chr10".to_string(), 500),
            ("chr1".to_string(), 500),
        ],
        vec![
            fragment_on_tid(paired_fragment(10, 120, 40), 0),
            fragment_on_tid(paired_fragment(20, 120, 40), 1),
            fragment_on_tid(paired_fragment(30, 120, 40), 2),
        ],
        Vec::new(),
        "bam_to_bam_all_skip_sort",
    )?;

    let work = tempdir()?;
    let out_bam = work.path().join("all_unsorted.bam");
    let chrom_args = ChromosomeArgs {
        chromosomes: Some(vec!["all".to_string()]),
        chromosomes_file: None,
    };
    let mut cfg = BamToBamConfig::new(bam.bam.clone(), out_bam.clone(), chrom_args);
    cfg.skip_chromosome_sort = true;

    // Act
    run_inner(&cfg)?;

    // Assert
    let mut reader = bam::Reader::from_path(&out_bam)?;
    let header = reader.header().to_owned();
    let mut record_chromosomes = Vec::new();
    for record_result in reader.records() {
        let record = record_result?;
        let tid = record.tid() as u32;
        record_chromosomes.push(
            std::str::from_utf8(header.tid2name(tid))
                .unwrap()
                .to_string(),
        );
    }
    assert_eq!(
        record_chromosomes,
        vec![
            "chr2".to_string(),
            "chr2".to_string(),
            "chr10".to_string(),
            "chr10".to_string(),
            "chr1".to_string(),
            "chr1".to_string(),
        ],
        "`--skip-chromosome-sort` should preserve BAM-header-derived order for `--chromosomes all`"
    );

    Ok(())
}

#[test]
fn windows_keep_fragment_when_single_read_is_inside() -> Result<()> {
    // Human verification status: unverified
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
    // Human verification status: unverified
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
    // Human verification status: unverified
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
    // Human verification status: unverified
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
    // Human verification status: unverified
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

#[test]
fn paired_and_unpaired_fragment_modes_apply_same_full_fragment_scaling_for_same_span() -> Result<()>
{
    // Human verification status: unverified
    // Arrange:
    // We represent the same physical fragment span [20, 80) in two different supported input
    // modes:
    // - paired-end: the `simple_inward_bam()` fixture has one fragment with forward [20, 40) and
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
    // - paired mode must emit two records, each tagged with COV = 4/3 and FLEN = 60
    // - unpaired mode must emit one record tagged with COV = 4/3 and FLEN = 60
    let paired_bam = simple_inward_bam()?;
    let unpaired_bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos: 20,
            cigar: vec![('M', 60)],
            seq: vec![b'A'; 60],
            qual: 30,
            is_reverse: false,
            mapq: 40,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
        "bam_to_bam_unpaired_fragment_scaling",
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
    paired_cfg.scale_genome.scaling_factors = Some(scaling.clone());
    paired_cfg.set_min_mapq(0);
    {
        let frag = paired_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }

    let mut unpaired_cfg = base_config(&unpaired_bam.bam, &unpaired_out);
    unpaired_cfg.scale_genome.scaling_factors = Some(scaling);
    unpaired_cfg.unpaired.reads_are_fragments = true;
    unpaired_cfg.set_min_mapq(0);
    {
        let frag = unpaired_cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 100;
    }

    // Act
    run_inner(&paired_cfg)?;
    run_inner(&unpaired_cfg)?;

    // Assert
    let paired_cov = read_tag_values(&paired_out, b"COV")?;
    let unpaired_cov = read_tag_values(&unpaired_out, b"COV")?;
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
fn gc_file_fallback_writes_gc_tag_one_on_both_mates() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_fallback.bam");
    let gc_path = work.path().join("gc_pkg.npz");
    build_gc_package(&gc_path, 26)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        skip_invalid_gc: false,
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
    // - `bam-to-bam` follows the same fallback semantics as `bam-to-frag` here:
    //   the fragment is kept, `gc_failed_fragments` increments once per fragment,
    //   and both mates receive a `GC` AUX tag of 1.0.
    let counters = run_inner(&cfg)?;

    assert_eq!(counters.base.counted_fragments, 1);
    assert_eq!(counters.gc_failed_fragments, 1);

    let gc_weights = read_tag_values(&out_bam, b"GC")?;
    assert_eq!(
        gc_weights,
        vec![1.0_f32, 1.0_f32],
        "both mates should receive the GC fallback weight"
    );

    let lengths = read_fragment_lengths(&out_bam)?;
    assert_eq!(
        lengths,
        vec![60_u32, 60_u32],
        "the fragment should still be emitted with its FLEN tags"
    );

    Ok(())
}

#[test]
fn gc_file_drop_invalid_skips_fragment_entirely() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_drop_invalid.bam");
    let gc_path = work.path().join("gc_pkg.npz");
    build_gc_package(&gc_path, 26)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        skip_invalid_gc: true,
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
    // - with `skip_invalid_gc=true`, the command must skip the fragment instead of silently
    //   falling back to weight 1.0
    // - therefore the output BAM must contain no records and no AUX tags at all
    let counters = run_inner(&cfg)?;

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
        "no FLEN tags should be written when every fragment is skipped"
    );

    Ok(())
}

#[test]
fn gc_file_and_scaling_factors_write_identical_gc_cov_and_flen_tags_on_both_mates() -> Result<()> {
    // Human verification status: unverified
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_and_cov.bam");
    let gc_path = work.path().join("gc_pkg.npz");
    let scaling_path = work.path().join("scaling.tsv");
    build_gc_package(&gc_path, 0)?;
    write_scaling_file(&scaling_path, "chr1", 200, 4.0_f32 / 3.0_f32)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.scale_genome.scaling_factors = Some(scaling_path);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        skip_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    // Manual expectations:
    // - `simple_inward_bam()` contains one fragment spanning [20,80), so FLEN must be 60
    // - the whole-chrom scaling TSV sets factor 4/3 everywhere, so both mates must receive
    //   identical `COV = 4/3`
    // - the helper GC package has:
    //     length bins [10,60) and [60,200]
    //     GC bins     [0,50) and [50,101]
    //     correction matrix row for length 60 = [2, 10]
    // - the repeated ACGT reference gives GC%=50 for any 60 bp fragment, so the fragment lands in
    //   the second GC bin and both mates must receive `GC = 10`
    // - because `bam-to-bam` tags per fragment, both output records must carry the same
    //   `GC`, `COV`, and `FLEN` values
    let counters = run_inner(&cfg)?;

    assert_eq!(counters.base.counted_fragments, 1);
    assert_eq!(read_tag_values(&out_bam, b"GC")?, vec![10.0_f32, 10.0_f32]);
    assert_eq!(
        read_tag_values(&out_bam, b"COV")?,
        vec![4.0_f32 / 3.0_f32, 4.0_f32 / 3.0_f32]
    );
    assert_eq!(read_fragment_lengths(&out_bam)?, vec![60_u32, 60_u32]);

    Ok(())
}

#[test]
fn gc_file_rejects_package_when_fragment_length_range_is_outside_supported_range() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // The fixture contributes one fragment of length 60. We configure the converter to accept
    // exactly that length, then write the smallest valid GC package that only covers 10..=59.
    //
    // The converter validates the GC package before it starts chromosome iteration, so the
    // expected error is the shared compatibility failure from the GC loader.
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_range_error.bam");
    let gc_path = work.path().join("gc_pkg_short.npz");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![10, 59],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        correction_matrix: array![[1.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        skip_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 60;
        frag.max_fragment_length = 60;
    }

    // Act
    let err = run_inner(&cfg).expect_err("out-of-range GC package should fail");

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
    // Human verification status: unverified
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
    // - first mate pair tagged GC=5.0, FLEN=10
    // - remaining nine mate pairs tagged GC=5/9, FLEN=10
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
        .map(|start| paired_fragment(start, 10, 5))
        .collect();
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        fragments,
        Vec::new(),
        "bam_to_bam_real_non_neutral_bam",
    )?;
    let work = tempdir()?;
    let out_bam = work.path().join("real_non_neutral_gc.bam");
    let gc_path = build_real_non_neutral_gc_package(
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
        skip_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 10;
    }

    // Act
    let counters = run_inner(&cfg)?;

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
    // Human verification status: unverified
    let reference = twobit_from_sequences(
        "bam_to_bam_combined_filters_reference",
        vec![("chr1".to_string(), "ACGT".repeat(75))],
    )?;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 300)],
        vec![
            paired_fragment(20, 60, 20),
            paired_fragment(100, 60, 20),
            paired_fragment(220, 60, 20),
        ],
        Vec::new(),
        "bam_to_bam_combined_filters_bam",
    )?;
    let work = tempdir()?;
    let out_bam = work.path().join("combined_filters.bam");
    let gc_path = work.path().join("gc_pkg.npz");
    let scaling_path = work.path().join("scaling.tsv");
    let bed_path = work.path().join("windows.bed");
    let blacklist_path = work.path().join("blacklist.bed");
    build_gc_package(&gc_path, 0)?;
    write_scaling_file(&scaling_path, "chr1", 300, 4.0_f32 / 3.0_f32)?;
    write_bed(&bed_path, &[(0, 180)])?;
    fs::write(&blacklist_path, "chr1\t120\t130\n")?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.scale_genome.scaling_factors = Some(scaling_path);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        skip_invalid_gc: false,
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
    //     scaling factor is uniform -> `COV = 4/3`
    //     fragment length -> `FLEN = 60`
    let counters = run_inner(&cfg)?;

    assert_eq!(counters.base.counted_fragments, 1);
    assert_eq!(counters.blacklisted_fragments, 1);
    assert_eq!(
        read_qname_counts(&out_bam)?,
        HashMap::from([(String::from("frag0_20"), 2_usize)])
    );
    assert_eq!(read_tag_values(&out_bam, b"GC")?, vec![10.0_f32, 10.0_f32]);
    assert_eq!(
        read_tag_values(&out_bam, b"COV")?,
        vec![4.0_f32 / 3.0_f32, 4.0_f32 / 3.0_f32]
    );
    assert_eq!(read_fragment_lengths(&out_bam)?, vec![60_u32, 60_u32]);

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn real_multi_chromosome_coverage_weights_tsv_is_applied_per_chromosome_in_bam_to_bam() -> Result<()>
{
    // Human verification status: unverified
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
    // - two chr1 reads with COV = 275/432 and FLEN = 60
    // - two chr2 reads with COV = 25/24   and FLEN = 20
    let mut chr2_fragment = paired_fragment(20, 20, 10);
    chr2_fragment = fragment_on_tid(chr2_fragment, 1);

    let bam = bam_from_specs(
        vec![("chr1".to_string(), 200), ("chr2".to_string(), 200)],
        vec![paired_fragment(20, 60, 20), chr2_fragment],
        Vec::new(),
        "bam_to_bam_real_multi_chr_scaling",
    )?;
    let work = tempdir()?;

    let weights_out_dir = work.path().join("weights_out");
    std::fs::create_dir_all(&weights_out_dir)?;
    let mut weights_cfg = CoverageWeightsConfig::new(
        cfdnalab::commands::cli_common::IOCArgs {
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

    let scaling_path = weights_out_dir.join("coverage.scaling_factors.tsv");
    let out_bam = work.path().join("scaled_multi_chr.bam");
    let mut cfg = BamToBamConfig::new(
        bam.bam.clone(),
        out_bam.clone(),
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string(), "chr2".to_string()]),
            chromosomes_file: None,
        },
    );
    cfg.skip_chromosome_sort = true;
    cfg.scale_genome.scaling_factors = Some(scaling_path);
    cfg.min_mapq = 0;
    {
        let frag = cfg.fragment_lengths_mut();
        frag.min_fragment_length = 10;
        frag.max_fragment_length = 200;
    }

    // Act
    let counters = run_inner(&cfg)?;

    // Assert
    assert_eq!(counters.base.counted_fragments, 2);
    let mut reader = bam::Reader::from_path(&out_bam)?;
    let header = reader.header().to_owned();
    let mut observed: Vec<(String, f32, u32)> = Vec::new();
    for record_result in reader.records() {
        let record = record_result?;
        let tid = record.tid() as u32;
        let chromosome = std::str::from_utf8(header.tid2name(tid))?.to_string();
        let scaling = match record.aux(b"COV") {
            Ok(Aux::Float(value)) => value,
            other => panic!("expected COV float tag on every read, got {other:?}"),
        };
        let flen = match record.aux(b"FLEN") {
            Ok(Aux::U32(value)) => value,
            other => panic!("expected FLEN u32 tag on every read, got {other:?}"),
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
            "unexpected FLEN for {expected_chrom}: expected {expected_flen}, got {flen}"
        );
    }

    Ok(())
}

#[test]
fn gc_file_rejects_package_with_schema_version_mismatch() -> Result<()> {
    // Human verification status: unverified
    // Arrange:
    // Build the smallest valid GC correction package shape, but with an intentionally
    // incompatible schema version. `bam-to-bam` should reject it while loading the package,
    // before chromosome processing starts or any output BAM records are written.
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_bad_version.bam");
    let gc_path = work.path().join("gc_pkg_bad_version.npz");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION + 1,
        end_offset: 0,
        length_edges: vec![10, 200],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        correction_matrix: array![[1.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        skip_invalid_gc: false,
    });
    cfg.set_ref_2bit(Some(ref_twobit.path.clone()));

    // Act
    let err = run_inner(&cfg).expect_err("schema version mismatch should fail");

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
    // Human verification status: unverified
    // Arrange:
    // `simple_inward_bam()` uses chr1 length 200.
    // The shared scaling loader requires bins to cover the whole requested chromosome exactly.
    //
    // This TSV stops at 100:
    //   [0,100) factor 2.0
    // so the converter must reject it before chromosome processing begins.
    let bam = simple_inward_bam()?;
    let work = tempdir()?;
    let out_bam = work.path().join("scaling_truncated.bam");
    let scaling_path = work.path().join("truncated_scaling.tsv");
    write_scaling_file(&scaling_path, "chr1", 100, 2.0)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.scale_genome.scaling_factors = Some(scaling_path);

    // Act
    let err = run_inner(&cfg).expect_err("truncated scaling TSV should fail");

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

fn build_gc_package(path: &Path, end_offset: u64) -> Result<()> {
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset,
        length_edges: vec![10, 60, 200],
        gc_edges: vec![0, 50, 101],
        length_bin_frequencies: array![1.0_f64, 3.0_f64],
        correction_matrix: array![[1.0_f64, 1.0_f64], [2.0_f64, 10.0_f64]],
    };
    package.write_npz(path)?;
    Ok(())
}
