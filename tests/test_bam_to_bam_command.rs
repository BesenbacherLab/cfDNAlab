#![cfg(feature = "cmd_bam_to_bam")]

mod fixtures;

use std::{collections::HashMap, fs, path::Path};

use anyhow::Result;
use cfdnalab::commands::{
    bam_to_bam::{bam_to_bam::run_inner, config::BamToBamConfig},
    cli_common::{ApplyGCArgFileOnly, ChromosomeArgs},
    gc_bias::{GC_CORRECTION_SCHEMA_VERSION, package::GCCorrectionPackage},
};
use fixtures::{
    FragmentSpec, ReadSpec, bam_from_specs, paired_fragment, simple_inward_bam,
    simple_reference_twobit,
};
use ndarray::array;
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
    assert_eq!(default_counts, HashMap::from([("frag0_20".to_string(), 2usize)]));
    assert!(
        explicit_thirty_counts.is_empty(),
        "raising min_mapq to 30 should remove the MAPQ-20 fragment"
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
fn chromosomes_all_default_sorts_lexicographically_instead_of_bam_header_order() -> Result<()> {
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

#[test]
fn gc_file_fallback_writes_gc_tag_one_on_both_mates() -> Result<()> {
    let bam = simple_inward_bam()?;
    let ref_twobit = simple_reference_twobit()?;
    let work = tempdir()?;
    let out_bam = work.path().join("gc_fallback.bam");
    let gc_path = work.path().join("gc_pkg.npz");
    build_gc_package(&gc_path, 26)?;

    let mut cfg = base_config(&bam.bam, &out_bam);
    cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        drop_invalid_gc: false,
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
fn gc_file_rejects_package_when_fragment_length_range_is_outside_supported_range() -> Result<()> {
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
        drop_invalid_gc: false,
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
fn gc_file_rejects_package_with_schema_version_mismatch() -> Result<()> {
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
        drop_invalid_gc: false,
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
