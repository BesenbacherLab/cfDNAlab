#![cfg(feature = "cmd_ends")]

mod fixtures;

use anyhow::{Context, Result};
#[cfg(feature = "cmd_gc_bias")]
use cfdnalab::commands::gc_bias::{GC_CORRECTION_SCHEMA_VERSION, package::GCCorrectionPackage};
use cfdnalab::commands::{
    cli_common::{ApplyGCArgs, ChromosomeArgs, IOCArgs, UnpairedArgs, WindowsArgs},
    ends::{
        config::EndsConfig,
        config_structs::{
            AssignMotifToWindowArgs, BaseQualityFilter, ClipStrategy, DEFAULT_MAX_SOFT_CLIPS,
            KmerSource, WindowMotifAssigner,
        },
        ends::run,
    },
};
use cfdnalab::shared::{blacklist::BlacklistStrategy, indel_mode::IndelMotifFilterPolicy};
use fixtures::{
    BamFixture, FragmentSpec, ReadSpec, bam_from_specs, paired_fragment, simple_reference_twobit,
    single_read_bam_with_qualities, twobit_from_sequences, write_bed,
};
#[cfg(feature = "cmd_gc_bias")]
use ndarray::array;
use ndarray::{Array1, Array2};
use ndarray_npy::{NpzReader, read_npy};
use serde_json::{Map, Value, json};
#[cfg(feature = "cli")]
use std::process::Command;
use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::{Path, PathBuf},
};
use tempfile::TempDir;

fn base_chromosomes(chrs: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(chrs.iter().map(|chr| chr.to_string()).collect()),
        chromosomes_file: None,
    }
}

fn base_config(
    bam_path: &Path,
    output_dir: &Path,
    k_inside: usize,
    k_outside: usize,
) -> EndsConfig {
    let mut cfg = EndsConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        k_inside,
        k_outside,
    );
    cfg.output_prefix = "ends".to_string();
    cfg.tile_size = 100;
    cfg.min_mapq = 0;
    cfg.require_proper_pair = false;
    cfg.clip.clip_strategy = ClipStrategy::Aligned;
    cfg
}

fn set_exact_fragment_length(cfg: &mut EndsConfig, len: u32) {
    let lengths = cfg.fragment_lengths_mut();
    lengths.min_fragment_length = len;
    lengths.max_fragment_length = len;
}

#[test]
fn ends_config_new_defaults_to_skip_clip_strategy() -> Result<()> {
    // Arrange
    let out_dir = TempDir::new()?;
    let cfg = EndsConfig::new(
        IOCArgs {
            bam: out_dir.path().join("dummy.bam"),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        1,
        0,
    );

    // Assert
    assert_eq!(cfg.clip.clip_strategy, ClipStrategy::Skip);
    assert_eq!(cfg.clip.max_soft_clips, DEFAULT_MAX_SOFT_CLIPS);
    Ok(())
}

#[test]
fn end_scope_bq_filter_drops_only_the_failing_end() -> Result<()> {
    // Arrange: left end quality is 40 and right end quality is 10.
    //
    // Mental derivation for k_inside=1, k_outside=0, source_inside=read:
    // - left read base `A` contributes `_A`
    // - right read base `G` contributes `_C` after right-end reverse complementation
    // - `min in end >= 30` keeps the left end and drops the right end
    let mut fragment = paired_fragment(10, 10, 4);
    fragment.forward.seq = b"AAAA".to_vec();
    fragment.forward.qual = 40;
    fragment.reverse.seq = b"GGGG".to_vec();
    fragment.reverse.qual = 10;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 256)],
        vec![fragment],
        Vec::new(),
        "ends_bq_end_filter",
    )?;
    let baseline_out_dir = TempDir::new()?;
    let mut baseline_cfg = base_config(&bam.bam, baseline_out_dir.path(), 1, 0);
    baseline_cfg.all_motifs = true;
    set_exact_fragment_length(&mut baseline_cfg, 10);

    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.all_motifs = true;
    set_exact_fragment_length(&mut cfg, 10);
    cfg.bq_filters = vec!["min in end >= 30".parse::<BaseQualityFilter>().unwrap()];

    // Act
    run(&baseline_cfg)?;
    let (baseline_motifs, baseline_matrix) = read_dense_output(baseline_out_dir.path())?;
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: without BQ filters, the same 10 bp fragment is counted at both ends, so this test
    // cannot pass because the default fragment-length filter removed the fragment earlier.
    assert_eq!(baseline_matrix.sum(), 2.0);
    assert_eq!(
        motif_count(&baseline_matrix, &baseline_motifs, 0, "_A"),
        1.0
    );
    assert_eq!(
        motif_count(&baseline_matrix, &baseline_motifs, 0, "_C"),
        1.0
    );
    assert_eq!(matrix.sum(), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_C"), 0.0);
    Ok(())
}

#[test]
fn fragment_scope_bq_filter_drops_the_full_fragment_when_it_fails() -> Result<()> {
    // Arrange: the fragment mean across both kept end qualities is (40 + 10) / 2 = 25.
    let mut fragment = paired_fragment(10, 10, 4);
    fragment.forward.seq = b"AAAA".to_vec();
    fragment.forward.qual = 40;
    fragment.reverse.seq = b"GGGG".to_vec();
    fragment.reverse.qual = 10;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 256)],
        vec![fragment],
        Vec::new(),
        "ends_bq_fragment_filter",
    )?;
    let baseline_out_dir = TempDir::new()?;
    let mut baseline_cfg = base_config(&bam.bam, baseline_out_dir.path(), 1, 0);
    baseline_cfg.all_motifs = true;
    set_exact_fragment_length(&mut baseline_cfg, 10);

    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.all_motifs = true;
    set_exact_fragment_length(&mut cfg, 10);
    cfg.bq_filters = vec![
        "mean in fragment >= 30"
            .parse::<BaseQualityFilter>()
            .unwrap(),
    ];

    // Act
    run(&baseline_cfg)?;
    let (_baseline_motifs, baseline_matrix) = read_dense_output(baseline_out_dir.path())?;
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(baseline_matrix.sum(), 2.0);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn fragment_scope_bq_filters_apply_before_end_scope_filters_remove_one_end() -> Result<()> {
    // Arrange: left end quality 40 passes the end filter, right end quality 10 fails it,
    // but the fragment mean over both ends is still only 25.
    //
    // Mental derivation:
    // - `min in end >= 30` would keep only the left end
    // - `mean in fragment >= 30` fails on the raw candidate fragment
    // - so the fragment must contribute nothing
    let mut fragment = paired_fragment(10, 10, 4);
    fragment.forward.seq = b"AAAA".to_vec();
    fragment.forward.qual = 40;
    fragment.reverse.seq = b"GGGG".to_vec();
    fragment.reverse.qual = 10;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 256)],
        vec![fragment],
        Vec::new(),
        "ends_bq_filter_order",
    )?;
    let baseline_out_dir = TempDir::new()?;
    let mut baseline_cfg = base_config(&bam.bam, baseline_out_dir.path(), 1, 0);
    baseline_cfg.all_motifs = true;
    set_exact_fragment_length(&mut baseline_cfg, 10);

    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.all_motifs = true;
    set_exact_fragment_length(&mut cfg, 10);
    cfg.bq_filters = vec![
        "min in end >= 30".parse::<BaseQualityFilter>().unwrap(),
        "mean in fragment >= 30"
            .parse::<BaseQualityFilter>()
            .unwrap(),
    ];

    // Act
    run(&baseline_cfg)?;
    let (_baseline_motifs, baseline_matrix) = read_dense_output(baseline_out_dir.path())?;
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(baseline_matrix.sum(), 2.0);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn fragment_scope_bq_filter_can_pass_while_end_scope_filter_drops_one_end() -> Result<()> {
    // Arrange: left end quality is 40 and right end quality is 20.
    //
    // Mental derivation for k_inside=1:
    // - fragment mean is (40 + 20) / 2 = 30, so `mean in fragment >= 30` passes
    // - `min in end >= 30` still drops the right end
    // - only the left `_A` motif should remain
    let mut fragment = paired_fragment(10, 10, 4);
    fragment.forward.seq = b"AAAA".to_vec();
    fragment.forward.qual = 40;
    fragment.reverse.seq = b"GGGG".to_vec();
    fragment.reverse.qual = 20;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 256)],
        vec![fragment],
        Vec::new(),
        "ends_bq_fragment_passes_end_drops_one",
    )?;
    let baseline_out_dir = TempDir::new()?;
    let mut baseline_cfg = base_config(&bam.bam, baseline_out_dir.path(), 1, 0);
    baseline_cfg.all_motifs = true;
    set_exact_fragment_length(&mut baseline_cfg, 10);

    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.all_motifs = true;
    set_exact_fragment_length(&mut cfg, 10);
    cfg.bq_filters = vec![
        "mean in fragment >= 30"
            .parse::<BaseQualityFilter>()
            .unwrap(),
        "min in end >= 30".parse::<BaseQualityFilter>().unwrap(),
    ];

    // Act
    run(&baseline_cfg)?;
    let (baseline_motifs, baseline_matrix) = read_dense_output(baseline_out_dir.path())?;
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(baseline_matrix.sum(), 2.0);
    assert_eq!(
        motif_count(&baseline_matrix, &baseline_motifs, 0, "_A"),
        1.0
    );
    assert_eq!(
        motif_count(&baseline_matrix, &baseline_motifs, 0, "_C"),
        1.0
    );
    assert_eq!(matrix.sum(), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_C"), 0.0);
    Ok(())
}

#[test]
fn fragment_scope_bq_filter_can_pass_while_end_scope_filters_drop_both_ends() -> Result<()> {
    // Arrange: both end qualities are 20.
    //
    // Mental derivation for k_inside=1:
    // - fragment mean is (20 + 20) / 2 = 20, so `mean in fragment >= 20` passes
    // - `min in end >= 30` fails on both ends
    // - once both ends are removed, the fragment must contribute nothing
    let mut fragment = paired_fragment(10, 10, 4);
    fragment.forward.seq = b"AAAA".to_vec();
    fragment.forward.qual = 20;
    fragment.reverse.seq = b"GGGG".to_vec();
    fragment.reverse.qual = 20;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 256)],
        vec![fragment],
        Vec::new(),
        "ends_bq_fragment_passes_both_ends_fail",
    )?;
    let baseline_out_dir = TempDir::new()?;
    let mut baseline_cfg = base_config(&bam.bam, baseline_out_dir.path(), 1, 0);
    baseline_cfg.all_motifs = true;
    set_exact_fragment_length(&mut baseline_cfg, 10);

    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.all_motifs = true;
    set_exact_fragment_length(&mut cfg, 10);
    cfg.bq_filters = vec![
        "mean in fragment >= 20"
            .parse::<BaseQualityFilter>()
            .unwrap(),
        "min in end >= 30".parse::<BaseQualityFilter>().unwrap(),
    ];

    // Act
    run(&baseline_cfg)?;
    let (_baseline_motifs, baseline_matrix) = read_dense_output(baseline_out_dir.path())?;
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(baseline_matrix.sum(), 2.0);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn reads_are_fragments_supports_combined_end_and_fragment_bq_filters() -> Result<()> {
    // Arrange: left inside base quality is 40 and right inside base quality is 10.
    //
    // Mental derivation for k_inside=1:
    // - the read spans one fragment, so left and right scores come from the first and last bases
    // - fragment mean is (40 + 10) / 2 = 25, so `mean in fragment >= 25` passes
    // - `min in end >= 30` drops only the right end
    // - the first base is `A`, so only `_A` remains after filtering
    let bam = single_read_bam_with_qualities(
        "ends_unpaired_bq_filters",
        10,
        vec![('M', 10)],
        b"AAAAAAAAAG",
        &[40, 30, 30, 30, 30, 30, 30, 30, 30, 10],
    )?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.all_motifs = true;
    cfg.bq_filters = vec![
        "mean in fragment >= 25"
            .parse::<BaseQualityFilter>()
            .unwrap(),
        "min in end >= 30".parse::<BaseQualityFilter>().unwrap(),
    ];
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.sum(), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_C"), 0.0);
    Ok(())
}

#[test]
fn reads_are_fragments_drop_the_fragment_when_fragment_filter_passes_but_both_ends_fail()
-> Result<()> {
    // Arrange: the first and last base qualities are both 20.
    //
    // Mental derivation for k_inside=1:
    // - fragment mean is (20 + 20) / 2 = 20, so `mean in fragment >= 20` passes
    // - `min in end >= 30` fails for both ends
    // - after both ends are removed, there is nothing left to count
    let bam = single_read_bam_with_qualities(
        "ends_unpaired_bq_fragment_passes_both_ends_fail",
        10,
        vec![('M', 10)],
        b"AAAAAAAAAA",
        &[20, 30, 30, 30, 30, 30, 30, 30, 30, 20],
    )?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.all_motifs = true;
    cfg.bq_filters = vec![
        "mean in fragment >= 20"
            .parse::<BaseQualityFilter>()
            .unwrap(),
        "min in end >= 30".parse::<BaseQualityFilter>().unwrap(),
    ];
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn aligned_bq_filter_uses_aligned_inside_qualities_for_k_inside_gt_one() -> Result<()> {
    // Arrange: 2S10M2S with high-quality clipped bases and low-quality aligned interior.
    //
    // Mental derivation for k_inside=3:
    // - aligned left slice is [30, 10, 10], so `min in end >= 30` fails
    // - aligned right slice is [10, 10, 30], so it also fails
    // - if aligned slicing is implemented correctly, neither end is counted
    let bam = single_read_bam_with_qualities(
        "ends_aligned_bq_k3",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
        &[40, 35, 30, 10, 10, 10, 10, 10, 10, 10, 10, 30, 35, 40],
    )?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 3, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Aligned;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.bq_filters = vec!["min in end >= 30".parse::<BaseQualityFilter>().unwrap()];
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_bq_filter_uses_raw_inside_qualities_for_k_inside_gt_one() -> Result<()> {
    // Arrange: the same 2S10M2S read should now use the raw terminal bases.
    //
    // Mental derivation for k_inside=3:
    // - raw left slice is [40, 35, 30], so `min in end >= 30` passes
    // - raw right slice is [30, 35, 40], so it also passes
    // - if raw-aligned slicing regresses back to aligned slices, the count would drop to zero
    let bam = single_read_bam_with_qualities(
        "ends_raw_aligned_bq_k3",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
        &[40, 35, 30, 10, 10, 10, 10, 10, 10, 10, 10, 30, 35, 40],
    )?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 3, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.bq_filters = vec!["min in end >= 30".parse::<BaseQualityFilter>().unwrap()];
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_bq_filter_uses_raw_inside_qualities_for_k_inside_gt_one() -> Result<()> {
    // Arrange: raw-shifted-boundary uses the same raw quality slices as raw-aligned-boundary but
    // keeps the shifted fragment length of 14 bp.
    let bam = single_read_bam_with_qualities(
        "ends_raw_shifted_bq_k3",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
        &[40, 35, 30, 10, 10, 10, 10, 10, 10, 10, 10, 30, 35, 40],
    )?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 3, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.bq_filters = vec!["min in end >= 30".parse::<BaseQualityFilter>().unwrap()];
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 14;
        lengths.max_fragment_length = 14;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn bq_filter_rejects_reference_backed_inside_bases_with_descriptive_error() -> Result<()> {
    // Arrange
    let bam = simple_paired_fragment_bam("ends_bq_reference_source_error", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.bq_filters = vec!["min in end >= 30".parse::<BaseQualityFilter>().unwrap()];

    // Act
    let err = run(&cfg).expect_err("reference-backed inside bases should reject BQ filters");

    // Assert
    assert!(err.to_string().contains("`--bq-filter`"));
    assert!(err.to_string().contains("`--source-inside reference`"));
    assert!(err.to_string().contains("read base qualities"));
    Ok(())
}

#[test]
fn bq_filter_rejects_zero_inside_bases_with_descriptive_error() -> Result<()> {
    // Arrange
    let bam = simple_paired_fragment_bam("ends_bq_k_inside_zero_error", 10, 10, 4)?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 0, 1);
    cfg.bq_filters = vec!["min in end >= 30".parse::<BaseQualityFilter>().unwrap()];

    // Act
    let err = run(&cfg).expect_err("BQ filters without inside bases should fail");

    // Assert
    assert!(
        err.to_string()
            .contains("`--bq-filter` requires `--k-inside > 0`")
    );
    Ok(())
}

#[test]
fn bq_filter_rejects_missing_base_qualities_with_descriptive_error() -> Result<()> {
    // Arrange: BAM uses 255 placeholders to mean missing qualities.
    let bam = single_read_bam_with_qualities(
        "ends_bq_missing_qualities_error",
        10,
        vec![('M', 10)],
        b"AAAAAAAAAA",
        &[255; 10],
    )?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.all_motifs = true;
    cfg.bq_filters = vec!["min in end >= 30".parse::<BaseQualityFilter>().unwrap()];
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    let err = run(&cfg).expect_err("missing BAM qualities should fail loudly");
    let error_text = format!("{err:#}");

    // Assert
    assert!(error_text.contains("missing base qualities"));
    assert!(error_text.contains("--bq-filter"));
    Ok(())
}

#[test]
fn settings_json_includes_bq_filters_in_command_output() -> Result<()> {
    // Arrange
    let bam = simple_paired_fragment_bam("ends_bq_settings", 10, 10, 4)?;
    let out_dir = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.bq_filters = vec![
        "min in end >= 30".parse::<BaseQualityFilter>().unwrap(),
        "max in fragment < 20".parse::<BaseQualityFilter>().unwrap(),
    ];

    // Act
    run(&cfg)?;
    let settings = read_text_file(&settings_path(out_dir.path()))?;

    // Assert
    assert_eq!(
        parse_json(&settings).get("bq_filters"),
        Some(&json!(["min in end >= 30", "max in fragment < 20"]))
    );
    Ok(())
}

fn simple_paired_fragment_bam(
    name: &str,
    start: i64,
    fragment_len: i64,
    read_len: i64,
) -> Result<BamFixture> {
    bam_from_specs(
        vec![("chr1".to_string(), 256)],
        vec![paired_fragment(start, fragment_len, read_len)],
        Vec::new(),
        name,
    )
}

fn single_read_bam(
    name: &str,
    pos: i64,
    cigar: Vec<(char, u32)>,
    seq: &[u8],
) -> Result<BamFixture> {
    bam_from_specs(
        vec![("chr1".to_string(), 256)],
        Vec::new(),
        vec![ReadSpec {
            tid: 0,
            pos,
            cigar,
            seq: seq.to_vec(),
            qual: 40,
            is_reverse: false,
            mapq: 60,
            flags: 0,
            mate_tid: None,
            mate_pos: None,
            insert_size: 0,
        }],
        name,
    )
}

fn custom_paired_fragment_bam(
    name: &str,
    forward: ReadSpec,
    reverse: ReadSpec,
) -> Result<BamFixture> {
    bam_from_specs(
        vec![("chr1".to_string(), 256)],
        vec![FragmentSpec { forward, reverse }],
        Vec::new(),
        name,
    )
}

fn fragment_on_tid(tid: usize, start: i64, fragment_len: i64, read_len: i64) -> FragmentSpec {
    let mut fragment = paired_fragment(start, fragment_len, read_len);
    fragment.forward.tid = tid;
    fragment.reverse.tid = tid;
    fragment.forward.mate_tid = Some(tid);
    fragment.reverse.mate_tid = Some(tid);
    fragment
}

fn three_chrom_reference_end_fixture(name: &str) -> Result<(BamFixture, fixtures::TwoBitFixture)> {
    let bam = bam_from_specs(
        vec![
            ("chr1".to_string(), 200),
            ("chr2".to_string(), 200),
            ("chr3".to_string(), 200),
        ],
        vec![
            fragment_on_tid(0, 20, 60, 20),
            fragment_on_tid(1, 30, 80, 20),
            fragment_on_tid(2, 40, 100, 20),
        ],
        Vec::new(),
        name,
    )?;

    // Mental derivation for k_inside=1, k_outside=0, source_inside=reference:
    // - chr1 fragment [20,80):  left base A -> "_A", right base A -> revcomp "_T"
    // - chr2 fragment [30,110): left base C -> "_C", right base C -> revcomp "_G"
    // - chr3 fragment [40,140): left base T -> "_T", right base G -> revcomp "_C"
    let reference = twobit_from_sequences(
        &format!("{name}_reference"),
        vec![
            ("chr1".to_string(), "A".repeat(200)),
            ("chr2".to_string(), "C".repeat(200)),
            (
                "chr3".to_string(),
                format!("{}T{}G{}", "A".repeat(40), "A".repeat(98), "A".repeat(60)),
            ),
        ],
    )?;

    Ok((bam, reference))
}

fn dense_output_paths(out_dir: &Path) -> (PathBuf, PathBuf) {
    (
        out_dir.join("ends.end_motifs.npy"),
        out_dir.join("ends.end_motifs.txt"),
    )
}

fn sparse_output_paths(out_dir: &Path) -> (PathBuf, PathBuf) {
    (
        out_dir.join("ends.end_motifs.sparse.npz"),
        out_dir.join("ends.end_motifs.txt"),
    )
}

fn settings_path(out_dir: &Path) -> PathBuf {
    out_dir.join("ends.end_motif_settings.json")
}

fn read_motif_labels(path: &Path) -> Result<Vec<String>> {
    let reader = BufReader::new(File::open(path)?);
    reader
        .lines()
        .collect::<std::io::Result<Vec<_>>>()
        .context("read motif labels")
}

fn read_dense_output(out_dir: &Path) -> Result<(Vec<String>, Array2<f64>)> {
    let (matrix_path, motifs_path) = dense_output_paths(out_dir);
    let motifs = read_motif_labels(&motifs_path)?;
    let matrix: Array2<f64> = read_npy(&matrix_path)?;
    Ok((motifs, matrix))
}

fn read_sparse_output(out_dir: &Path) -> Result<(Vec<String>, Array2<f64>)> {
    let (npz_path, motifs_path) = sparse_output_paths(out_dir);
    let motifs = read_motif_labels(&motifs_path)?;
    let file = File::open(&npz_path)?;
    let mut npz = NpzReader::new(file)?;
    let row: Array1<u64> = npz.by_name("row.npy")?;
    let col: Array1<u64> = npz.by_name("col.npy")?;
    let data: Array1<f64> = npz.by_name("data.npy")?;
    let shape: Array1<i64> = npz.by_name("shape.npy")?;

    let n_rows = shape[0] as usize;
    let n_cols = shape[1] as usize;
    let mut dense = Array2::<f64>::zeros((n_rows, n_cols));
    for ((&r, &c), &value) in row.iter().zip(col.iter()).zip(data.iter()) {
        dense[(r as usize, c as usize)] = value;
    }

    Ok((motifs, dense))
}

fn motif_count(matrix: &Array2<f64>, motifs: &[String], row: usize, motif: &str) -> f64 {
    let column = motifs
        .iter()
        .position(|label| label == motif)
        .unwrap_or_else(|| panic!("missing motif column {motif} in {:?}", motifs));
    matrix[(row, column)]
}

fn expected_combined_1_plus_1_dense_order_without_collapse() -> Vec<String> {
    let bases = ["A", "C", "G", "T"];
    let mut motifs = Vec::new();
    for outside in bases {
        for inside in bases {
            motifs.push(format!("{outside}_{inside}"));
        }
    }
    motifs
}

fn expected_combined_1_plus_1_dense_order_with_collapse() -> Vec<String> {
    let bases = ["A", "C", "G", "T"];
    let mut motifs = Vec::new();
    for outside in ["A", "C"] {
        for inside in bases {
            motifs.push(format!("{outside}_{inside}"));
        }
    }
    motifs
}

fn expected_combined_2_plus_2_dense_order_without_collapse() -> Vec<String> {
    let bases = ["A", "C", "G", "T"];
    let mut motifs = Vec::new();
    for first_outside in bases {
        for second_outside in bases {
            for first_inside in bases {
                for second_inside in bases {
                    motifs.push(format!(
                        "{first_outside}{second_outside}_{first_inside}{second_inside}"
                    ));
                }
            }
        }
    }
    motifs
}

fn expected_combined_2_plus_2_dense_order_with_collapse() -> Vec<String> {
    let bases = ["A", "C", "G", "T"];
    let mut motifs = Vec::new();
    for first_outside in ["A", "C"] {
        for second_outside in bases {
            for first_inside in bases {
                for second_inside in bases {
                    motifs.push(format!(
                        "{first_outside}{second_outside}_{first_inside}{second_inside}"
                    ));
                }
            }
        }
    }
    motifs
}

fn read_text_file(path: &Path) -> Result<String> {
    let mut buf = String::new();
    File::open(path)?.read_to_string(&mut buf)?;
    Ok(buf)
}

fn parse_json(text: &str) -> Value {
    serde_json::from_str(text).expect("settings sidecar should be valid JSON")
}

fn expected_settings_json(
    source_inside: &str,
    clip_strategy: &str,
    window_assignment: &str,
) -> Value {
    let mut expected = Map::new();
    expected.insert("source_inside".to_string(), json!(source_inside));
    expected.insert("clip_strategy".to_string(), json!(clip_strategy));
    expected.insert("window_assignment".to_string(), json!(window_assignment));
    #[cfg(feature = "ends_experimental")]
    expected.insert("collapse_complement".to_string(), json!(false));
    Value::Object(expected)
}

#[cfg(feature = "cli")]
fn cfdna_binary_path() -> Result<String> {
    std::env::var("CARGO_BIN_EXE_cfdna")
        .context("CARGO_BIN_EXE_cfdna is not set for this CLI integration test")
}

#[test]
fn blacklist_any_skips_a_fragment_before_any_end_motifs_are_counted() -> Result<()> {
    // Arrange: fragment [10,20) overlaps the blacklist at [15,16), so blacklist_strategy=Any
    // should exclude the fragment before either end motif is counted.
    let bam = simple_paired_fragment_bam("ends_blacklist_fragment", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let blacklist_bed = out_dir.path().join("blacklist.bed");
    write_bed(&blacklist_bed, &[("chr1", 15, 16, "blk")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.blacklist = Some(vec![blacklist_bed]);
    cfg.blacklist_strategy = BlacklistStrategy::Any;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn blacklist_masking_skips_only_the_reference_backed_end_motif_that_overlaps_a_blacklisted_base()
-> Result<()> {
    // Arrange: fragment [10,20) has left terminal base at 10 and right terminal base at 19.
    // Blacklisting [10,11) masks only the left inside base. Using blacklist_strategy=All keeps
    // the fragment itself, so only the left endpoint motif should disappear.
    let bam = simple_paired_fragment_bam("ends_blacklist_reference_end", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let blacklist_bed = out_dir.path().join("blacklist.bed");
    write_bed(&blacklist_bed, &[("chr1", 10, 11, "blk")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.blacklist = Some(vec![blacklist_bed]);
    cfg.blacklist_strategy = BlacklistStrategy::All;
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert: right terminal base at 19 is T, which orients to "_A". The left "_G" is masked.
    assert_eq!(motifs, vec!["_A"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 1.0);
    Ok(())
}

#[cfg(feature = "cli")]
#[test]
fn cli_statistics_count_one_fragment_and_one_motif_when_only_one_end_survives() -> Result<()> {
    // Arrange: the fragment survives fragment-level filtering, but the blacklist masks only
    // the left terminal base. The right end motif still counts, so the public statistics should
    // report one counted fragment and one counted end motif.
    //
    // Mental derivation:
    // - the fragment itself is kept because `blacklist-strategy=all` and only one base is blacklisted
    // - the left endpoint motif is skipped because its terminal base overlaps the masked base
    // - the right endpoint motif still counts
    // - therefore the public stats must say:
    //   * 1 fragment with one or more counted motifs
    //   * 1 distinct counted end motif across those fragments
    let bam = simple_paired_fragment_bam("ends_cli_stats_one_end", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let blacklist_bed = out_dir.path().join("blacklist.bed");
    write_bed(&blacklist_bed, &[("chr1", 10, 11, "blk")])?;
    let binary = cfdna_binary_path()?;

    // Act
    let output = Command::new(binary)
        .args([
            "ends",
            "--bam",
            bam.bam.to_str().context("bam path is not valid UTF-8")?,
            "--output-dir",
            out_dir
                .path()
                .to_str()
                .context("output dir is not valid UTF-8")?,
            "--chromosomes",
            "chr1",
            "--k-inside",
            "1",
            "--k-outside",
            "0",
            "--ref-2bit",
            reference
                .path
                .to_str()
                .context("reference path is not valid UTF-8")?,
            "--source-inside",
            "reference",
            "--blacklist",
            blacklist_bed
                .to_str()
                .context("blacklist path is not valid UTF-8")?,
            "--blacklist-strategy",
            "all",
            "--assign-by",
            "endpoint",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "10",
            "--min-mapq",
            "0",
            "--tile-size",
            "1000000",
            "--output-prefix",
            "ends",
            "--n-threads",
            "1",
        ])
        .output()
        .context("running cfdna ends CLI")?;

    // Assert
    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).context("stdout is not valid UTF-8")?;
    assert!(stdout.contains("Fragments with one or more counted motifs: 1"));
    assert!(stdout.contains("Distinct counted end motifs across those fragments: 1"));
    Ok(())
}

#[cfg(feature = "cli")]
#[test]
fn cli_statistics_only_count_reads_from_tiles_with_relevant_windows() -> Result<()> {
    // Arrange: two paired fragments on a 2 Mb chromosome with a 1 Mb tile size.
    //
    // Mental derivation:
    // - fragment A starts at 10, so it belongs to tile 0
    // - fragment B starts at 1_500_000, so it belongs to tile 1
    // - the BED file contains only one window around fragment B, so tile 0 has no relevant windows
    //   and is skipped before any BAM reads are scanned there
    // - each paired fragment contributes 2 reads, so:
    //   * whole BAM contains 4 reads
    //   * processed tiles contain only the 2 reads from fragment B
    // - fragment B survives and both of its end motifs count, so public stats should report:
    //   * 2 observed reads in processed tiles
    //   * 1 fragment with counted motifs
    //   * 2 distinct counted end motifs
    let fragment_a = paired_fragment(10, 10, 4);
    let fragment_b = paired_fragment(1_500_000, 10, 4);
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 2_000_100)],
        vec![fragment_a, fragment_b],
        Vec::new(),
        "ends_cli_stats_skipped_tile",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 1_500_000, 1_500_010, "processed_tile_window")],
    )?;
    let binary = cfdna_binary_path()?;

    // Act
    let output = Command::new(binary)
        .args([
            "ends",
            "--bam",
            bam.bam.to_str().context("bam path is not valid UTF-8")?,
            "--output-dir",
            out_dir
                .path()
                .to_str()
                .context("output dir is not valid UTF-8")?,
            "--chromosomes",
            "chr1",
            "--k-inside",
            "1",
            "--k-outside",
            "0",
            "--by-bed",
            windows_bed
                .to_str()
                .context("windows BED path is not valid UTF-8")?,
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "10",
            "--min-mapq",
            "0",
            "--tile-size",
            "1000000",
            "--output-prefix",
            "ends",
            "--n-threads",
            "1",
        ])
        .output()
        .context("running cfdna ends CLI with skipped tile")?;

    // Assert
    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).context("stdout is not valid UTF-8")?;
    assert!(stdout.contains("Note: counts below cover only tiles with relevant output windows"));
    assert!(stdout.contains("Observed reads in processed tiles: 2"));
    assert!(stdout.contains("Initially accepted observed reads: 2"));
    assert!(stdout.contains("Fragments with one or more counted motifs: 1"));
    assert!(stdout.contains("Distinct counted end motifs across those fragments: 2"));
    Ok(())
}

#[cfg(feature = "cli")]
#[test]
fn cli_source_inside_reference_writes_the_expected_dense_reference_backed_counts() -> Result<()> {
    // Arrange: this fixture intentionally makes the read-backed and reference-backed outcomes
    // disagree, so the test fails loudly if `--source-inside reference` is ignored.
    //
    // `paired_fragment(...)` writes forward read bases as `AAAA` and reverse read bases as `TTTT`.
    // For fragment [10,20):
    //
    // - reference-backed inside bases on the ACGT-repeat reference are:
    //   - left inside base at 10 = G -> "_G"
    //   - right inside base at 19 = T -> reverse-complement "_A"
    //   - expected dense counts: `_A = 1`, `_G = 1`
    //
    // - read-backed inside bases would instead be:
    //   - left inside from forward read `AAAA` -> "_A"
    //   - right inside from reverse read `TTTT` -> reverse-complement "_A"
    //   - wrong dense counts if the CLI ignored `--source-inside reference`: `_A = 2`, `_G = 0`
    let bam = simple_paired_fragment_bam("ends_cli_reference_source", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let binary = cfdna_binary_path()?;

    // Act
    let output = Command::new(binary)
        .args([
            "ends",
            "--bam",
            bam.bam.to_str().context("bam path is not valid UTF-8")?,
            "--output-dir",
            out_dir
                .path()
                .to_str()
                .context("output dir is not valid UTF-8")?,
            "--chromosomes",
            "chr1",
            "--k-inside",
            "1",
            "--k-outside",
            "0",
            "--ref-2bit",
            reference
                .path
                .to_str()
                .context("reference path is not valid UTF-8")?,
            "--source-inside",
            "reference",
            "--all-motifs",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "10",
            "--min-mapq",
            "0",
            "--tile-size",
            "1000000",
            "--output-prefix",
            "ends",
            "--n-threads",
            "1",
        ])
        .output()
        .context("running cfdna ends CLI with --source-inside reference")?;

    // Assert
    assert!(
        output.status.success(),
        "stdout:\n{}\n\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let (motifs, matrix) = read_dense_output(out_dir.path())?;
    assert_eq!(motifs, vec!["_A", "_C", "_G", "_T"]);
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn blacklist_masking_still_skips_read_backed_inside_motifs_using_genomic_reference_coordinates()
-> Result<()> {
    // Arrange: unpaired read-fragment [10,20) with read sequence A C G A A C G A A A.
    // - left read-backed motif = "_A"
    // - right read-backed motif = reverse-complement("A") = "_T"
    // Blacklisting [10,11) should drop only the left motif even though inside bases come from
    // the read, because blacklist validation is still genomic.
    let bam = single_read_bam(
        "ends_blacklist_read_end",
        10,
        vec![('M', 10)],
        b"ACGAACGAAA",
    )?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let blacklist_bed = out_dir.path().join("blacklist.bed");
    write_bed(&blacklist_bed, &[("chr1", 10, 11, "blk")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    cfg.blacklist = Some(vec![blacklist_bed]);
    cfg.blacklist_strategy = BlacklistStrategy::All;
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["_T"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 1.0);
    Ok(())
}

#[test]
fn blacklist_masking_skips_outside_bases_independently_of_inside_source() -> Result<()> {
    // Arrange: fragment [11,21) with k_outside=1 and k_inside=1.
    //
    // The left outside base is genomic position 10. Blacklisting [10,11) therefore masks only the
    // left outside base and must skip the full left-end motif regardless of where the inside base
    // comes from.
    //
    // This fixture intentionally makes the surviving right-end label differ between inside-source
    // modes, so we can prove that only the outside masking is shared:
    //
    // - reference-backed inside:
    //   - right inside genomic base at 20 = A
    //   - right outside genomic base at 21 = C
    //   - storage order = "AC", decode RC("AC") = "GT" -> "G_T"
    //
    // - read-backed inside:
    //   - reverse read base is T (from the paired fixture's `TTTT`)
    //   - right outside genomic base at 21 = C
    //   - storage order = "TC", decode RC("TC") = "GA" -> "G_A"
    //
    // Without outside masking, both ends would contribute the same label within each mode, so the
    // surviving label would have count 2. After masking [10,11), the left end must disappear and
    // the surviving right-end label must have count exactly 1 in both modes.
    let bam = simple_paired_fragment_bam("ends_blacklist_outside_source_independent", 11, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir_read = TempDir::new()?;
    let out_dir_reference = TempDir::new()?;
    let blacklist_bed = out_dir_read.path().join("blacklist.bed");
    write_bed(&blacklist_bed, &[("chr1", 10, 11, "mask_left_outside")])?;

    let run_with_source =
        |output_dir: &Path, source_inside: KmerSource| -> Result<(Vec<String>, Array2<f64>)> {
            let mut cfg = base_config(&bam.bam, output_dir, 1, 1);
            cfg.set_ref_2bit(Some(reference.path.clone()));
            cfg.source_inside = source_inside;
            cfg.all_motifs = false;
            cfg.blacklist = Some(vec![blacklist_bed.clone()]);
            cfg.blacklist_strategy = BlacklistStrategy::All;
            {
                let lengths = cfg.fragment_lengths_mut();
                lengths.min_fragment_length = 10;
                lengths.max_fragment_length = 10;
            }

            run(&cfg)?;
            read_sparse_output(output_dir)
        };

    // Act
    let (read_motifs, read_matrix) = run_with_source(out_dir_read.path(), KmerSource::Read)?;
    let (reference_motifs, reference_matrix) =
        run_with_source(out_dir_reference.path(), KmerSource::Reference)?;

    // Assert
    assert_eq!(read_motifs, vec!["G_A"]);
    assert_eq!(read_matrix.shape(), &[1, 1]);
    assert_eq!(read_matrix[(0, 0)], 1.0);

    assert_eq!(reference_motifs, vec!["G_T"]);
    assert_eq!(reference_matrix.shape(), &[1, 1]);
    assert_eq!(reference_matrix[(0, 0)], 1.0);
    Ok(())
}

#[test]
fn scaling_factors_weight_each_counted_end_motif() -> Result<()> {
    // Arrange: one chromosome-wide scaling factor of 2.0 should double both endpoint counts for
    // fragment [10,20), whose reference-backed motifs are "_G" on the left and "_A" on the right.
    let bam = simple_paired_fragment_bam("ends_scaling", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let scaling_path = out_dir.path().join("scaling.tsv");
    std::fs::write(
        &scaling_path,
        "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t256\t2\n",
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_scaling_factors(Some(scaling_path));
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 2.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 2.0);
    assert_eq!(matrix.sum(), 4.0);
    Ok(())
}

#[cfg(feature = "cmd_gc_bias")]
#[test]
fn gc_file_weights_each_counted_end_motif_by_the_fragment_gc_correction() -> Result<()> {
    // Arrange: fragment [10,20) on the ACGT-repeat reference has GC%=50 over 10 bp.
    // The package below assigns weight 3.0 to length bin [10,11) and GC bin [0,51), so both
    // endpoint motifs should each be counted with weight 3.0.
    let bam = simple_paired_fragment_bam("ends_gc_weight", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("gc_package.npz");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![10, 11, 20],
        gc_edges: vec![0, 51, 100],
        length_bin_frequencies: array![1.0_f64, 1.0_f64],
        correction_matrix: array![[3.0_f64, 1.0_f64], [1.0_f64, 1.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 3.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 3.0);
    assert_eq!(matrix.sum(), 6.0);
    Ok(())
}

#[cfg(feature = "cmd_gc_bias")]
#[test]
fn blacklist_gc_and_scaling_weights_combine_to_the_exact_expected_endpoint_counts() -> Result<()> {
    // Arrange: this fixture combines all three weighting/filtering layers in one run.
    //
    // We use four paired fragments of length 10 on one chromosome:
    //
    //   F1: [10,20)   = AAAAAAAAAA     -> GC 0%
    //   F2: [40,50)   = CCCCCCCCCC     -> GC 100%
    //   F3: [70,80)   = AAAAACCCCC     -> GC 50%
    //   F4: [100,110) = GGGGGGGGGG     -> fully blacklisted
    //
    // The scaling bins fully cover the chromosome:
    //
    //   [0,30)   -> 1.0
    //   [30,60)  -> 2.0
    //   [60,75)  -> 1.0
    //   [75,90)  -> 2.0
    //   [90,120) -> 1.0
    //
    // So the fragment-level scaling averages are:
    //
    //   F1: 10 bp inside [0,30)                 -> 1.0
    //   F2: 10 bp inside [30,60)                -> 2.0
    //   F3: 5 bp in [60,75) and 5 bp in [75,90) -> (5*1 + 5*2) / 10 = 1.5
    //   F4: irrelevant because the fragment is dropped by full-fragment blacklisting
    //
    // The blacklist BED has two roles:
    //
    //   [40,41)   masks only the left terminal base of F2
    //   [100,110) fully covers F4, so `blacklist-strategy=all` drops that fragment
    //
    // The GC package is set up so length 10 gets:
    //
    //   GC 0%   -> weight 2.0
    //   GC 50%  -> weight 4.0
    //   GC 100% -> weight 3.0
    //
    // Endpoint windows are one base wide and isolate each expected end:
    //
    //   W0 [10,11)   -> F1 left
    //   W1 [19,20)   -> F1 right
    //   W2 [40,41)   -> F2 left, but masked by blacklist
    //   W3 [49,50)   -> F2 right
    //   W4 [70,71)   -> F3 left
    //   W5 [79,80)   -> F3 right
    //   W6 [100,101) -> F4 left, but whole fragment is blacklisted away
    //
    // Mental derivation of the exact final endpoint weights:
    //
    //   F1 left:  `_A` weight = GC 2.0 * scaling 1.0 = 2.0
    //   F1 right: `_T` weight = GC 2.0 * scaling 1.0 = 2.0
    //             right terminal base is A, which reverse-complements to T
    //
    //   F2 left:  skipped entirely because [40,41) is blacklist-masked
    //   F2 right: `_G` weight = GC 3.0 * scaling 2.0 = 6.0
    //             right terminal base is C, which reverse-complements to G
    //
    //   F3 left:  `_A` weight = GC 4.0 * scaling 1.5 = 6.0
    //   F3 right: `_G` weight = GC 4.0 * scaling 1.5 = 6.0
    //
    //   F4 left:  skipped because the full fragment is blacklisted by `all`
    //
    // Therefore the row sums must be:
    //
    //   [2.0, 2.0, 0.0, 6.0, 6.0, 6.0, 0.0]
    //
    // and the global motif totals must be:
    //
    //   `_A` = 8.0
    //   `_T` = 2.0
    //   `_G` = 12.0
    //   `_C` = 0.0
    //   total = 22.0
    let chromosome_sequence = format!(
        "{}{}{}{}{}{}{}{}{}",
        "T".repeat(10),
        "A".repeat(10),
        "T".repeat(20),
        "C".repeat(10),
        "T".repeat(20),
        "AAAAACCCCC",
        "T".repeat(20),
        "G".repeat(10),
        "T".repeat(10),
    );
    let reference = twobit_from_sequences(
        "ends_blacklist_gc_scaling_reference",
        vec![("chr1".to_string(), chromosome_sequence)],
    )?;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 120)],
        vec![
            paired_fragment(10, 10, 4),
            paired_fragment(40, 10, 4),
            paired_fragment(70, 10, 4),
            paired_fragment(100, 10, 4),
        ],
        Vec::new(),
        "ends_blacklist_gc_scaling",
    )?;
    let out_dir = TempDir::new()?;

    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[
            ("chr1", 10, 11, "f1_left"),
            ("chr1", 19, 20, "f1_right"),
            ("chr1", 40, 41, "f2_left_masked"),
            ("chr1", 49, 50, "f2_right"),
            ("chr1", 70, 71, "f3_left"),
            ("chr1", 79, 80, "f3_right"),
            ("chr1", 100, 101, "f4_left_blacklisted"),
        ],
    )?;

    let blacklist_bed = out_dir.path().join("blacklist.bed");
    write_bed(
        &blacklist_bed,
        &[
            ("chr1", 40, 41, "mask_f2_left"),
            ("chr1", 100, 110, "drop_f4"),
        ],
    )?;

    let scaling_path = out_dir.path().join("scaling.tsv");
    std::fs::write(
        &scaling_path,
        concat!(
            "chromosome\tstart\tend\tscaling_factor\n",
            "chr1\t0\t30\t1\n",
            "chr1\t30\t60\t2\n",
            "chr1\t60\t75\t1\n",
            "chr1\t75\t90\t2\n",
            "chr1\t90\t120\t1\n",
        ),
    )?;

    let gc_path = out_dir.path().join("gc_package.npz");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![10, 11, 20],
        gc_edges: vec![0, 1, 50, 51, 100],
        length_bin_frequencies: array![1.0_f64, 1.0_f64],
        correction_matrix: array![
            [2.0_f64, 1.0_f64, 4.0_f64, 3.0_f64],
            [1.0_f64, 1.0_f64, 1.0_f64, 1.0_f64]
        ],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.blacklist = Some(vec![blacklist_bed]);
    cfg.blacklist_strategy = BlacklistStrategy::All;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    cfg.set_scaling_factors(Some(scaling_path));
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    cfg.tile_size = 1_000_000;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: one row per BED window, one column per single-base motif.
    assert_eq!(matrix.shape(), &[7, 4]);

    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 2.0);
    assert_eq!(matrix.row(0).sum(), 2.0);

    assert_eq!(motif_count(&matrix, &motifs, 1, "_T"), 2.0);
    assert_eq!(matrix.row(1).sum(), 2.0);

    assert_eq!(matrix.row(2).sum(), 0.0);

    assert_eq!(motif_count(&matrix, &motifs, 3, "_G"), 6.0);
    assert_eq!(matrix.row(3).sum(), 6.0);

    assert_eq!(motif_count(&matrix, &motifs, 4, "_A"), 6.0);
    assert_eq!(matrix.row(4).sum(), 6.0);

    assert_eq!(motif_count(&matrix, &motifs, 5, "_G"), 6.0);
    assert_eq!(matrix.row(5).sum(), 6.0);

    assert_eq!(matrix.row(6).sum(), 0.0);

    assert_eq!(motif_count(&matrix, &motifs, 0, "_C"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 0.0);

    let total_a: f64 = matrix
        .column(
            motifs
                .iter()
                .position(|label| label == "_A")
                .context("missing _A motif")?,
        )
        .sum();
    let total_t: f64 = matrix
        .column(
            motifs
                .iter()
                .position(|label| label == "_T")
                .context("missing _T motif")?,
        )
        .sum();
    let total_g: f64 = matrix
        .column(
            motifs
                .iter()
                .position(|label| label == "_G")
                .context("missing _G motif")?,
        )
        .sum();
    let total_c: f64 = matrix
        .column(
            motifs
                .iter()
                .position(|label| label == "_C")
                .context("missing _C motif")?,
        )
        .sum();
    assert_eq!(total_a, 8.0);
    assert_eq!(total_t, 2.0);
    assert_eq!(total_g, 12.0);
    assert_eq!(total_c, 0.0);
    assert_eq!(matrix.sum(), 22.0);
    Ok(())
}

#[cfg(feature = "cmd_gc_bias")]
#[test]
fn default_gc_behavior_skips_fragments_when_gc_correction_cannot_be_computed() -> Result<()> {
    // Arrange: use a reference where the fragment GC window contains only `N`, so GC fraction
    // cannot be computed even though the correction package covers the fragment length. With
    // the default GC behavior, the fragment should be skipped instead of falling back to weight 1.0.
    let bam = simple_paired_fragment_bam("ends_neutralize_invalid_gc", 10, 10, 4)?;
    let reference = twobit_from_sequences(
        "ends_neutralize_invalid_gc_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}{}", "A".repeat(10), "N".repeat(10), "A".repeat(236)),
        )],
    )?;
    let out_dir = TempDir::new()?;
    let gc_path = out_dir.path().join("invalid_gc_package.npz");
    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![10, 11, 20],
        gc_edges: vec![0, 51, 100],
        length_bin_frequencies: array![1.0_f64, 1.0_f64],
        correction_matrix: array![[2.0_f64, 1.0_f64], [1.0_f64, 1.0_f64]],
    };
    package.write_npz(&gc_path)?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn auto_indel_filter_keeps_indel_affected_read_backed_end_motifs() -> Result<()> {
    // Arrange: the forward read has an insertion inside the left aligned 4-base footprint.
    // In auto mode with read-backed inside bases, that end should still be kept.
    let forward = ReadSpec {
        tid: 0,
        pos: 100,
        cigar: vec![('M', 2), ('I', 1), ('M', 5)],
        seq: b"ACGTACGT".to_vec(),
        qual: 40,
        is_reverse: false,
        mapq: 60,
        flags: 0x40 | 0x20 | 0x2,
        mate_tid: Some(0),
        mate_pos: Some(110),
        insert_size: 16,
    };
    let reverse = ReadSpec {
        tid: 0,
        pos: 110,
        cigar: vec![('M', 6)],
        seq: b"AACCGG".to_vec(),
        qual: 40,
        is_reverse: true,
        mapq: 60,
        flags: 0x80 | 0x2,
        mate_tid: Some(0),
        mate_pos: Some(100),
        insert_size: -16,
    };
    let bam = custom_paired_fragment_bam("ends_auto_indel_read", forward, reverse)?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 4, 0);
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    cfg.set_indel_filter(IndelMotifFilterPolicy::Auto);
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 16;
        lengths.max_fragment_length = 16;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert: both ends survive.
    //
    // Mental derivation:
    // - left read-backed inside bases come directly from the forward read prefix: `ACGT`
    // - right read-backed inside bases come from the reverse read suffix `CCGG`
    // - right ends are reverse-complemented on decode, but `CCGG` is palindromic
    // - therefore the two observed labels are exactly `_ACGT` and `_CCGG`
    assert_eq!(matrix.shape(), &[1, 2]);
    assert_eq!(motifs, vec!["_ACGT", "_CCGG"]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_ACGT"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_CCGG"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn auto_indel_filter_skips_indel_affected_reference_backed_end_motifs() -> Result<()> {
    // Arrange: the same pair as above, but now auto mode uses reference-backed inside bases and
    // should therefore skip the indel-affected left end while keeping the right end.
    let forward = ReadSpec {
        tid: 0,
        pos: 100,
        cigar: vec![('M', 2), ('I', 1), ('M', 5)],
        seq: b"ACGTACGT".to_vec(),
        qual: 40,
        is_reverse: false,
        mapq: 60,
        flags: 0x40 | 0x20 | 0x2,
        mate_tid: Some(0),
        mate_pos: Some(110),
        insert_size: 16,
    };
    let reverse = ReadSpec {
        tid: 0,
        pos: 110,
        cigar: vec![('M', 6)],
        seq: b"AACCGG".to_vec(),
        qual: 40,
        is_reverse: true,
        mapq: 60,
        flags: 0x80 | 0x2,
        mate_tid: Some(0),
        mate_pos: Some(100),
        insert_size: -16,
    };
    let bam = custom_paired_fragment_bam("ends_auto_indel_reference", forward, reverse)?;
    let reference = twobit_from_sequences(
        "ends_auto_indel_reference_context",
        vec![(
            "chr1".to_string(),
            format!(
                "{}{}{}",
                "T".repeat(100),
                "TGCAAAAACCCCGGTT",
                "T".repeat(140)
            ),
        )],
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 4, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.set_indel_filter(IndelMotifFilterPolicy::Auto);
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 16;
        lengths.max_fragment_length = 16;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    //
    // Mental derivation on the custom reference:
    // - left reference-backed 4-mer would have been `TGCA`, but that end is indel-affected and skipped
    // - right reference-backed 4-mer is ref[112..116) = `GGTT`
    // - right-end decoding reverse-complements that to `AACC`
    // - therefore the only surviving motif is `_AACC`
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(motifs, vec!["_AACC"]);
    assert_eq!(matrix[(0, 0)], 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn skip_affected_fragment_drops_the_whole_fragment_when_any_end_motif_is_indel_affected()
-> Result<()> {
    // Arrange: the same pair as above has an indel in the left end motif footprint, so
    // skip-affected-fragment must suppress all counting.
    let forward = ReadSpec {
        tid: 0,
        pos: 100,
        cigar: vec![('M', 2), ('I', 1), ('M', 5)],
        seq: b"ACGTACGT".to_vec(),
        qual: 40,
        is_reverse: false,
        mapq: 60,
        flags: 0x40 | 0x20 | 0x2,
        mate_tid: Some(0),
        mate_pos: Some(110),
        insert_size: 16,
    };
    let reverse = ReadSpec {
        tid: 0,
        pos: 110,
        cigar: vec![('M', 6)],
        seq: b"AACCGG".to_vec(),
        qual: 40,
        is_reverse: true,
        mapq: 60,
        flags: 0x80 | 0x2,
        mate_tid: Some(0),
        mate_pos: Some(100),
        insert_size: -16,
    };
    let bam = custom_paired_fragment_bam("ends_skip_affected_fragment", forward, reverse)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 4, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.set_indel_filter(IndelMotifFilterPolicy::SkipAffectedFragment);
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 16;
        lengths.max_fragment_length = 16;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn outside_reference_lookup_uses_preloaded_tile_reference_when_the_motif_extends_left_of_tile_fetch()
-> Result<()> {
    // Arrange: by-BED window [20,21) with max_fragment_length=10 gives tile fetch [10,31).
    // The motif preload now widens from full tile fetch by k_outside on both sides, so with
    // k_outside=11 the loaded reference span starts at 0. Asking for the left outside motif at
    // boundary 20 needs reference bases [9,20), which is left of tile.fetch but still inside the
    // preloaded reference span. On the ACGT-repeat reference, seq[9..20) is
    // C G T A C G T A C G T, so the outside-only label is "CGTACGTACGT_".
    let bam = single_read_bam(
        "ends_exact_reference_fallback",
        20,
        vec![('M', 10)],
        b"AAAAAAAAAA",
    )?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 20, 21, "left")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 0, 11);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["CGTACGTACGT_"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 1.0);
    Ok(())
}

#[test]
fn endpoint_assigns_left_and_right_end_motifs_to_separate_windows() -> Result<()> {
    // Arrange: fragment [10,20) on ACGT-repeat reference.
    // - left terminal base:  seq[10] = G  -> label "_G"
    // - right terminal base: seq[19] = T  -> oriented right-end label "_A"
    let bam = simple_paired_fragment_bam("ends_endpoint_split", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 10, 11, "left"), ("chr1", 19, 20, "right")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["_A", "_C", "_G", "_T"]);
    assert_eq!(matrix.shape(), &[2, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_A"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn midpoint_assigns_both_end_motifs_to_the_midpoint_window() -> Result<()> {
    // Arrange: fragment [10,20) has even midpoint 14 or 15, both inside [14,16).
    // So midpoint assignment should count both end motifs in that one window.
    let bam = simple_paired_fragment_bam("ends_midpoint", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 14, 16, "mid")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Midpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn count_overlap_weights_both_end_motifs_by_the_fragment_overlap_fraction() -> Result<()> {
    // Arrange: fragment [10,20) overlaps window [10,15) by 5 of 10 bp, so each end motif should
    // contribute 0.5 under count-overlap weighting.
    let bam = simple_paired_fragment_bam("ends_count_overlap", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 15, "half")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::CountOverlap,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 0.5);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 0.5);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn cross_tile_fragment_is_counted_once_per_window_when_it_reaches_into_the_next_tile() -> Result<()>
{
    // Arrange: tile size 20 puts [15,35) across two tile cores.
    // The fragment starts in the first tile core, so it should still be counted in both
    // overlapping windows [0,20) and [20,40), but only once overall after tile reduction.
    let bam = simple_paired_fragment_bam("ends_cross_tile_once", 15, 20, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.tile_size = 20;
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: Some(20),
        by_bed: None,
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Any,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 20;
        lengths.max_fragment_length = 20;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: [15,35) overlaps windows 0 and 1 only, and each of those rows should receive the
    // two end motifs exactly once rather than being doubled by neighboring tiles.
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_C"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_T"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_C"), 1.0);
    assert_eq!(matrix.row(0).sum(), 2.0);
    assert_eq!(matrix.row(1).sum(), 2.0);
    Ok(())
}

#[test]
fn all_requires_the_full_fragment_to_overlap_the_window() -> Result<()> {
    // Arrange: fragment [10,20) does not fully overlap [10,19), so "all" should count nothing.
    let bam = simple_paired_fragment_bam("ends_all", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 19, "almost")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::All,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn all_assignment_counts_the_fragment_when_the_window_fully_contains_it() -> Result<()> {
    // Arrange: fragment [10,20) is fully contained in [10,20), so `all` should accept it.
    let bam = simple_paired_fragment_bam("ends_all_accept", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 20, "full")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::All,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: `All` is fragment-centric once the fragment passes the overlap test, so both ends
    // are counted into the same window.
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn proportion_assignment_counts_the_fragment_when_the_requested_fraction_is_met() -> Result<()> {
    // Arrange: fragment [10,20) overlaps [10,15) by 5/10 bp, so proportion=0.5 should accept it.
    let bam = simple_paired_fragment_bam("ends_proportion", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 15, "half")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Proportion(0.5),
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn proportion_assignment_rejects_the_fragment_when_the_requested_fraction_is_not_met() -> Result<()>
{
    // Arrange: fragment [10,20) overlaps [10,14) by 4/10 bp, so proportion=0.5 should reject it.
    let bam = simple_paired_fragment_bam("ends_proportion_reject", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 14, "short")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Proportion(0.5),
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn any_assignment_counts_both_end_motifs_when_any_fragment_base_overlaps() -> Result<()> {
    // Arrange: fragment [10,20) overlaps [19,20) by exactly one base, which should still count
    // both end motifs under "any".
    let bam = simple_paired_fragment_bam("ends_any", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 19, 20, "one_bp")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Any,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn endpoint_counts_both_end_motifs_when_one_window_contains_both_terminal_bases() -> Result<()> {
    // Arrange: one window covering the full fragment contains both endpoint bases, so endpoint
    // assignment should place both motifs in the same row.
    let bam = simple_paired_fragment_bam("ends_endpoint_same_window", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 20, "full")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn collapse_complement_merges_complement_equivalent_single_base_end_motifs() -> Result<()> {
    // Arrange: fragment [11,21) on ACGT-repeat gives:
    // - left terminal base  seq[11] = T -> "_T"
    // - right terminal base seq[20] = A -> oriented right-end also "_T"
    // With complement collapsing enabled, both should map to canonical "_A".
    let bam = simple_paired_fragment_bam("ends_collapse", 11, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.collapse_complement = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: single-base collapse keeps only canonical A/C columns.
    assert_eq!(motifs, vec!["_A", "_C"]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 2.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn dense_all_motifs_output_enumerates_the_full_combined_1_plus_1_universe_without_collapse()
-> Result<()> {
    // Arrange: left endpoint only for fragment [10,20) with one outside and one inside
    // reference base. The observed motif is "C_G", but dense output must still write the full
    // combined 1+1 universe.
    let bam = simple_paired_fragment_bam("ends_dense_combined_no_collapse", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 11, "left")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(
        motifs,
        expected_combined_1_plus_1_dense_order_without_collapse()
    );
    assert_eq!(matrix.shape(), &[1, 16]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "C_G"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "A_A"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "T_T"), 0.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn dense_all_motifs_output_enumerates_the_full_collapsed_combined_1_plus_1_universe() -> Result<()>
{
    // Arrange: two left endpoints produce "G_T" and "C_A", which are same-orientation
    // complements and must both collapse to the canonical label "C_A".
    let reference = twobit_from_sequences(
        "ends_dense_collapsed_1_plus_1_reference",
        vec![(
            "chr1".to_string(),
            format!(
                "{}{}{}{}{}",
                "T".repeat(9),
                "GT",
                "T".repeat(8),
                "CA",
                "T".repeat(235)
            ),
        )],
    )?;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 256)],
        Vec::new(),
        vec![
            ReadSpec {
                tid: 0,
                pos: 10,
                cigar: vec![('M', 10)],
                seq: b"AAAAAAAAAA".to_vec(),
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            },
            ReadSpec {
                tid: 0,
                pos: 20,
                cigar: vec![('M', 10)],
                seq: b"CCCCCCCCCC".to_vec(),
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            },
        ],
        "ends_dense_collapsed_1_plus_1",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 10, 11, "left_a"), ("chr1", 20, 21, "left_b")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.collapse_complement = true;
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(
        motifs,
        expected_combined_1_plus_1_dense_order_with_collapse()
    );
    assert_eq!(matrix.shape(), &[2, 8]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "C_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "C_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "A_A"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "C_T"), 0.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn dense_all_motifs_output_enumerates_the_full_combined_2_plus_2_universe_without_collapse()
-> Result<()> {
    // Arrange: one left endpoint from a 10 bp fragment with a 2+2 reference motif. Dense output
    // must still enumerate the full 4^4 universe, not just the observed "GT_AC" column.
    let reference = twobit_from_sequences(
        "ends_dense_combined_2_plus_2_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}{}", "T".repeat(8), "GTAC", "T".repeat(244)),
        )],
    )?;
    let bam = single_read_bam(
        "ends_dense_combined_2_plus_2",
        10,
        vec![('M', 10)],
        b"AAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 11, "left")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 2, 2);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(
        motifs,
        expected_combined_2_plus_2_dense_order_without_collapse()
    );
    assert_eq!(matrix.shape(), &[1, 256]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "GT_AC"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "CA_TG"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "AA_AA"), 0.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn dense_all_motifs_output_enumerates_the_full_collapsed_combined_2_plus_2_universe() -> Result<()>
{
    // Arrange: two left endpoints produce "GT_AC" and "CA_TG", which are same-orientation
    // complements and must both land in the canonical "CA_TG" column.
    let reference = twobit_from_sequences(
        "ends_dense_collapsed_2_plus_2_reference",
        vec![(
            "chr1".to_string(),
            format!(
                "{}{}{}{}{}",
                "T".repeat(8),
                "GTAC",
                "T".repeat(6),
                "CATG",
                "T".repeat(234)
            ),
        )],
    )?;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 256)],
        Vec::new(),
        vec![
            ReadSpec {
                tid: 0,
                pos: 10,
                cigar: vec![('M', 10)],
                seq: b"AAAAAAAAAA".to_vec(),
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            },
            ReadSpec {
                tid: 0,
                pos: 20,
                cigar: vec![('M', 10)],
                seq: b"CCCCCCCCCC".to_vec(),
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            },
        ],
        "ends_dense_collapsed_2_plus_2",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 10, 11, "left_a"), ("chr1", 20, 21, "left_b")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 2, 2);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.collapse_complement = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(
        motifs,
        expected_combined_2_plus_2_dense_order_with_collapse()
    );
    assert_eq!(matrix.shape(), &[2, 128]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "CA_TG"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "CA_TG"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "AA_AA"), 0.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn collapse_complement_preserves_outside_inside_order_for_odd_length_end_motifs() -> Result<()> {
    // Arrange: two unpaired fragments of length 10 on a custom reference, with k_outside=1 and
    // k_inside=2 so the full motif length is 3.
    //
    // The intended contract is:
    // - decode first into biological 5'->3' `outside || inside` order
    // - then collapse against the same-orientation complement
    //
    // Fragment A spans [10,20):
    // - left full motif uses reference [9,12) = G T A -> "GTA" -> label "G_TA"
    // - right storage uses reference [18,21) = A T G -> revcomp("ATG") = "CAT" -> label "C_AT"
    //
    // Fragment B spans [30,40):
    // - left full motif uses reference [29,32) = C A T -> "CAT" -> label "C_AT"
    // - right storage uses reference [38,41) = A T G -> revcomp("ATG") = "CAT" -> label "C_AT"
    //
    // `GTA` and `CAT` are same-orientation complements:
    // - complement("GTA") = "CAT"
    // - complement("CAT") = "GTA"
    //
    // Lexicographically, "CAT" < "GTA", so all four counted ends must collapse to "C_AT".
    let reference = twobit_from_sequences(
        "ends_odd_length_complement_reference",
        vec![(
            "chr1".to_string(),
            format!(
                "{}{}{}{}{}{}{}{}{}",
                "T".repeat(9),
                "GTA",
                "T".repeat(6),
                "ATG",
                "T".repeat(8),
                "CAT",
                "T".repeat(6),
                "ATG",
                "T".repeat(215)
            ),
        )],
    )?;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 256)],
        Vec::new(),
        vec![
            ReadSpec {
                tid: 0,
                pos: 10,
                cigar: vec![('M', 10)],
                seq: b"AAAAAAAAAA".to_vec(),
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            },
            ReadSpec {
                tid: 0,
                pos: 30,
                cigar: vec![('M', 10)],
                seq: b"CCCCCCCCCC".to_vec(),
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            },
        ],
        "ends_odd_length_complement",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 2, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.collapse_complement = true;
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["C_AT"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 4.0);
    assert_eq!(matrix.sum(), 4.0);
    Ok(())
}

#[test]
fn collapse_complement_preserves_outside_inside_order_for_multi_base_2_plus_2_end_motifs()
-> Result<()> {
    // Arrange: two unpaired fragments of length 10 on a custom reference, with k_outside=2 and
    // k_inside=2 so the full motif length is 4.
    //
    // The intended contract is:
    // - decode first into biological 5'->3' `outside || inside` order
    // - then collapse against the same-orientation complement on the full 4-base motif
    // - only after that split into `<outside>_<inside>`
    //
    // Fragment A spans [10,20):
    // - left full motif uses reference [8,12) = G T A C -> "GTAC" -> label "GT_AC"
    // - right storage uses reference [18,22) = C A T G -> revcomp("CATG") = "CATG"
    //   -> label "CA_TG"
    //
    // Fragment B spans [30,40):
    // - left full motif uses reference [28,32) = T G C A -> "TGCA" -> label "TG_CA"
    // - right storage uses reference [38,42) = A C G T -> revcomp("ACGT") = "ACGT"
    //   -> label "AC_GT"
    //
    // Same-orientation complement pairs:
    // - complement("GTAC") = "CATG", so canonical label is "CA_TG"
    // - complement("TGCA") = "ACGT", so canonical label is "AC_GT"
    //
    // So the four counted ends must collapse into exactly two labels with count 2 each.
    let reference = twobit_from_sequences(
        "ends_2_plus_2_complement_reference",
        vec![(
            "chr1".to_string(),
            format!(
                "{}{}{}{}{}{}{}{}{}",
                "T".repeat(8),
                "GTAC",
                "T".repeat(6),
                "CATG",
                "T".repeat(6),
                "TGCA",
                "T".repeat(6),
                "ACGT",
                "T".repeat(214)
            ),
        )],
    )?;
    let bam = bam_from_specs(
        vec![("chr1".to_string(), 256)],
        Vec::new(),
        vec![
            ReadSpec {
                tid: 0,
                pos: 10,
                cigar: vec![('M', 10)],
                seq: b"AAAAAAAAAA".to_vec(),
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            },
            ReadSpec {
                tid: 0,
                pos: 30,
                cigar: vec![('M', 10)],
                seq: b"CCCCCCCCCC".to_vec(),
                qual: 40,
                is_reverse: false,
                mapq: 60,
                flags: 0,
                mate_tid: None,
                mate_pos: None,
                insert_size: 0,
            },
        ],
        "ends_2_plus_2_complement",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 2, 2);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.collapse_complement = true;
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["AC_GT", "CA_TG"]);
    assert_eq!(matrix.shape(), &[1, 2]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "AC_GT"), 2.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "CA_TG"), 2.0);
    assert_eq!(matrix.sum(), 4.0);
    Ok(())
}

#[test]
fn sparse_output_is_the_default_when_all_motifs_is_disabled() -> Result<()> {
    // Arrange: one fragment with two observed motifs should produce a sparse 1x2 matrix by default.
    let bam = simple_paired_fragment_bam("ends_sparse_default", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;
    let settings = read_text_file(&settings_path(out_dir.path()))?;

    // Assert
    assert!(!dense_output_paths(out_dir.path()).0.exists());
    assert!(sparse_output_paths(out_dir.path()).0.exists());
    assert_eq!(motifs, vec!["_A", "_G"]);
    assert_eq!(matrix.shape(), &[1, 2]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(
        parse_json(&settings),
        expected_settings_json("reference", "aligned", "endpoint")
    );
    Ok(())
}

#[test]
fn dense_all_motifs_output_still_uses_the_same_settings_sidecar() -> Result<()> {
    // Arrange: the sidecar should describe motif semantics, not mirror obvious output format state.
    let bam = simple_paired_fragment_bam("ends_dense_settings", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let settings = read_text_file(&settings_path(out_dir.path()))?;

    // Assert
    assert!(dense_output_paths(out_dir.path()).0.exists());
    assert!(!sparse_output_paths(out_dir.path()).0.exists());
    assert_eq!(
        parse_json(&settings),
        expected_settings_json("reference", "aligned", "endpoint")
    );
    Ok(())
}

#[test]
fn global_mode_accumulates_reference_end_motifs_across_three_chromosomes() -> Result<()> {
    // Arrange:
    // The three-chrom fixture contributes these reference-backed inside-only motifs:
    //
    // - chr1 [20,80):  `_A` and `_T`
    // - chr2 [30,110): `_C` and `_G`
    // - chr3 [40,140): `_T` and `_C`
    //
    // Global mode must collapse all three chromosomes into one output row, so the exact totals are:
    // - `_A` = 1
    // - `_C` = 2
    // - `_G` = 1
    // - `_T` = 2
    let (bam, reference) = three_chrom_reference_end_fixture("ends_three_chr_global")?;
    let out_dir = TempDir::new()?;

    let mut cfg = EndsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1", "chr2", "chr3"]),
        1,
        0,
    );
    cfg.output_prefix = "ends".to_string();
    cfg.tile_size = 50;
    cfg.min_mapq = 0;
    cfg.require_proper_pair = false;
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 60;
        lengths.max_fragment_length = 100;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(motifs, vec!["_A", "_C", "_G", "_T"]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_C"), 2.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 2.0);
    assert_eq!(matrix.sum(), 6.0);
    Ok(())
}

#[test]
fn by_size_windowing_keeps_three_chromosome_rows_in_chromosome_order() -> Result<()> {
    // Arrange:
    // With by_size=200 on three 200 bp chromosomes, the output must have one row per chromosome
    // in the selected chromosome order. The three-chrom fixture yields exact row motifs:
    //
    // - row 0 / chr1: `_A` = 1, `_T` = 1
    // - row 1 / chr2: `_C` = 1, `_G` = 1
    // - row 2 / chr3: `_C` = 1, `_T` = 1
    let (bam, reference) = three_chrom_reference_end_fixture("ends_three_chr_by_size")?;
    let out_dir = TempDir::new()?;

    let mut cfg = EndsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1", "chr2", "chr3"]),
        1,
        0,
    );
    cfg.output_prefix = "ends".to_string();
    cfg.tile_size = 50;
    cfg.min_mapq = 0;
    cfg.require_proper_pair = false;
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: Some(200),
        by_bed: None,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 60;
        lengths.max_fragment_length = 100;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;
    let bins_tsv = read_text_file(&out_dir.path().join("ends.bins.tsv"))?;

    // Assert
    assert_eq!(matrix.shape(), &[3, 4]);
    assert_eq!(motifs, vec!["_A", "_C", "_G", "_T"]);

    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_C"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 1.0);

    assert_eq!(motif_count(&matrix, &motifs, 1, "_A"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_C"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_G"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_T"), 0.0);

    assert_eq!(motif_count(&matrix, &motifs, 2, "_A"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 2, "_C"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 2, "_G"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 2, "_T"), 1.0);
    assert_eq!(matrix.sum(), 6.0);

    let rows: Vec<&str> = bins_tsv.lines().collect();
    assert_eq!(rows.len(), 4);
    assert_eq!(rows[0], "chrom\tstart\tend\tblacklisted_fraction");
    assert_eq!(rows[1], "chr1\t0\t200\t0");
    assert_eq!(rows[2], "chr2\t0\t200\t0");
    assert_eq!(rows[3], "chr3\t0\t200\t0");
    Ok(())
}

#[test]
fn bed_windowing_preserves_bed_row_order_and_skips_selected_chromosomes_without_windows()
-> Result<()> {
    // Arrange:
    // Chromosome selection includes chr1, chr2, chr3, but the BED file only contains chr3 then chr1.
    // chr2 must therefore contribute no output row at all.
    //
    // BED order is intentionally non-chromosomal to test original-index preservation:
    // - row 0 / chr3: `_C` = 1, `_T` = 1
    // - row 1 / chr1: `_A` = 1, `_T` = 1
    let (bam, reference) = three_chrom_reference_end_fixture("ends_three_chr_by_bed_sparse")?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows_three_chr.bed");
    write_bed(
        &windows_bed,
        &[
            ("chr3", 0, 200, "chr3_window"),
            ("chr1", 0, 200, "chr1_window"),
        ],
    )?;

    let mut cfg = EndsConfig::new(
        IOCArgs {
            bam: bam.bam.clone(),
            output_dir: out_dir.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1", "chr2", "chr3"]),
        1,
        0,
    );
    cfg.output_prefix = "ends".to_string();
    cfg.tile_size = 50;
    cfg.min_mapq = 0;
    cfg.require_proper_pair = false;
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 60;
        lengths.max_fragment_length = 100;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;
    let bins_tsv = read_text_file(&out_dir.path().join("ends.bins.tsv"))?;

    // Assert
    assert_eq!(matrix.shape(), &[2, 4]);
    assert_eq!(motifs, vec!["_A", "_C", "_G", "_T"]);

    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_C"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 1.0);

    assert_eq!(motif_count(&matrix, &motifs, 1, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_C"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_G"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_T"), 1.0);
    assert_eq!(matrix.sum(), 4.0);

    let rows: Vec<&str> = bins_tsv.lines().collect();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0], "chrom\tstart\tend\tblacklisted_fraction");
    assert_eq!(rows[1], "chr3\t0\t200\t0");
    assert_eq!(rows[2], "chr1\t0\t200\t0");
    assert!(rows.iter().skip(1).all(|row| !row.starts_with("chr2\t")));
    Ok(())
}

#[test]
fn by_size_and_bed_equivalent_full_chromosome_windows_match_across_three_chromosomes() -> Result<()>
{
    // Arrange:
    // - three chromosomes of length 200
    // - by-size 200 produces one full-chromosome row per chromosome
    // - BED windows [0,200) for each chromosome describe the exact same row partition
    // Under aligned endpoint counting, these two window modes must therefore produce identical
    // motif labels and identical row-wise count matrices.
    let (bam, reference) = three_chrom_reference_end_fixture("ends_three_chr_bed_vs_size")?;
    let by_size_out = TempDir::new()?;
    let bed_out = TempDir::new()?;
    let windows_bed = bed_out.path().join("windows_three_chr_full.bed");
    write_bed(
        &windows_bed,
        &[
            ("chr1", 0, 200, "chr1_window"),
            ("chr2", 0, 200, "chr2_window"),
            ("chr3", 0, 200, "chr3_window"),
        ],
    )?;

    let make_cfg = |out_dir: &Path, windows: WindowsArgs| {
        let mut cfg = EndsConfig::new(
            IOCArgs {
                bam: bam.bam.clone(),
                output_dir: out_dir.to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1", "chr2", "chr3"]),
            1,
            0,
        );
        cfg.output_prefix = "ends".to_string();
        cfg.tile_size = 50;
        cfg.min_mapq = 0;
        cfg.require_proper_pair = false;
        cfg.set_ref_2bit(Some(reference.path.clone()));
        cfg.source_inside = KmerSource::Reference;
        cfg.all_motifs = true;
        cfg.set_windows(windows);
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = 60;
            lengths.max_fragment_length = 100;
        }
        cfg
    };

    let by_size_cfg = make_cfg(
        by_size_out.path(),
        WindowsArgs {
            by_size: Some(200),
            by_bed: None,
        },
    );
    let bed_cfg = make_cfg(
        bed_out.path(),
        WindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed),
        },
    );

    // Act
    run(&by_size_cfg)?;
    run(&bed_cfg)?;

    // Assert
    let (by_size_motifs, by_size_matrix) = read_dense_output(by_size_out.path())?;
    let (bed_motifs, bed_matrix) = read_dense_output(bed_out.path())?;

    assert_eq!(by_size_motifs, vec!["_A", "_C", "_G", "_T"]);
    assert_eq!(by_size_motifs, bed_motifs);
    assert_eq!(by_size_matrix.shape(), &[3, 4]);
    assert_eq!(bed_matrix.shape(), &[3, 4]);
    assert_eq!(by_size_matrix, bed_matrix);
    assert_eq!(by_size_matrix.sum(), 6.0);
    Ok(())
}

#[test]
fn global_mode_counts_both_end_motifs_in_one_output_row() -> Result<()> {
    // Arrange: with no windows configured, the command should produce one global output row.
    let bam = simple_paired_fragment_bam("ends_global", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    assert!(!out_dir.path().join("ends.bins.tsv").exists());
    Ok(())
}

#[test]
fn all_motifs_dense_output_includes_zero_count_columns_for_unobserved_motifs() -> Result<()> {
    // Arrange: the one-fragment case only observes _A and _G, so _C and _T must still be present
    // as explicit zero columns under all-motifs output.
    let bam = simple_paired_fragment_bam("ends_zero_columns", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["_A", "_C", "_G", "_T"]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_C"), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 0.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn all_motifs_dense_output_enumerates_outside_only_labels_when_k_inside_is_zero() -> Result<()> {
    // Arrange: outside-only motifs should still have a fixed dense column universe.
    let bam = simple_paired_fragment_bam("ends_outside_only_dense", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 0, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, _matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["A_", "C_", "G_", "T_"]);
    Ok(())
}

#[test]
fn reference_backed_inside_and_outside_bases_are_combined_into_the_expected_labels() -> Result<()> {
    // Arrange: fragment [10,20) on the ACGT-repeat reference with `k_outside = 1, k_inside = 2`.
    //
    // Mental derivation:
    // - left outside base is seq[9]  = C
    // - left inside bases are seq[10..12) = GT
    //   -> left label = `C_GT`
    //
    // - right inside bases are seq[18..20) = GT
    // - right outside base is seq[20] = A
    // - right-end storage order is `GTA`, which reverse-complements to `TAC`
    //   -> right label = `T_AC`
    let bam = simple_paired_fragment_bam("ends_combined_inside_outside_reference", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 2, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["C_GT", "T_AC"]);
    assert_eq!(matrix.shape(), &[1, 2]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "C_GT"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "T_AC"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn config_docstring_visualization_example_counts_the_expected_motifs_with_and_without_collapse()
-> Result<()> {
    // Arrange: this is the exact example from the `EndsConfig` docstring.
    //
    // Reference indices:
    //   0  1  2  3  4  5  6  7  8  9  10 11 12 13 14
    //   A  T  C  G  T  T  T  T  T  T   T  C  A  T  C
    //
    // One read-as-fragment spans [2,13), so with `k_outside=2` and `k_inside=2`:
    // - left full motif uses reference [0,4)  = A T C G -> `ATCG` -> label `AT_CG`
    // - right storage uses reference [11,15) = C A T C -> `CATC`
    //   and right-end decode reverse-complements that to `GATG` -> label `GA_TG`
    //
    // With `collapse_complement=true`:
    // - complement(`ATCG`) = `TAGC`, so canonical label stays `AT_CG`
    // - complement(`GATG`) = `CTAC`, so canonical label becomes `CT_AC`
    let reference = twobit_from_sequences(
        "ends_config_docstring_visualization_reference",
        vec![("chr1".to_string(), "ATCGTTTTTTTCATC".to_string())],
    )?;
    let bam = single_read_bam(
        "ends_config_docstring_visualization",
        2,
        vec![('M', 11)],
        b"AAAAAAAAAAA",
    )?;

    let out_dir_uncollapsed = TempDir::new()?;
    let mut cfg = base_config(&bam.bam, out_dir_uncollapsed.path(), 2, 2);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 11;
        lengths.max_fragment_length = 11;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir_uncollapsed.path())?;

    // Assert
    assert_eq!(motifs, vec!["AT_CG", "GA_TG"]);
    assert_eq!(matrix.shape(), &[1, 2]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "AT_CG"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "GA_TG"), 1.0);
    assert_eq!(matrix.sum(), 2.0);

    let out_dir_collapsed = TempDir::new()?;
    let mut collapsed_cfg = base_config(&bam.bam, out_dir_collapsed.path(), 2, 2);
    collapsed_cfg.set_ref_2bit(Some(reference.path.clone()));
    collapsed_cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    collapsed_cfg.source_inside = KmerSource::Reference;
    collapsed_cfg.all_motifs = false;
    collapsed_cfg.collapse_complement = true;
    {
        let lengths = collapsed_cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 11;
        lengths.max_fragment_length = 11;
    }

    // Act
    run(&collapsed_cfg)?;
    let (collapsed_motifs, collapsed_matrix) = read_sparse_output(out_dir_collapsed.path())?;

    // Assert
    assert_eq!(collapsed_motifs, vec!["AT_CG", "CT_AC"]);
    assert_eq!(collapsed_matrix.shape(), &[1, 2]);
    assert_eq!(
        motif_count(&collapsed_matrix, &collapsed_motifs, 0, "AT_CG"),
        1.0
    );
    assert_eq!(
        motif_count(&collapsed_matrix, &collapsed_motifs, 0, "CT_AC"),
        1.0
    );
    assert_eq!(collapsed_matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_assignment_uses_the_shifted_assignment_boundaries() -> Result<()> {
    // Arrange: unpaired read-as-fragment with 2S10M2S at pos 10.
    // - aligned interval [10,20)
    // - raw assignment interval [8,22)
    // - endpoint positions 8 and 21
    // The raw terminal bases are T on the left and A on the right, which both orient to "_T".
    let bam = single_read_bam(
        "ends_raw_shifted",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 8, 9, "left_raw"), ("chr1", 21, 22, "right_raw")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 14;
        lengths.max_fragment_length = 14;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[2, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_T"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_assignment_keeps_a_left_window_that_only_raw_reach_can_touch()
-> Result<()> {
    // Arrange: unpaired 2S10M2S at pos 10 with tile size 10.
    //
    // Mental derivation:
    // - aligned interval is [10,20), so the fragment belongs to tile core [10,20)
    // - raw assignment interval is [8,22), so the left endpoint is 8
    // - BED window [8,9) lies left of the tile core and is only relevant because raw clipping
    //   moves the counted endpoint there
    // - the left raw base is T, which decodes to "_T"
    let bam = single_read_bam(
        "ends_raw_left_of_tile_core",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 8, 9, "left_raw_only")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 14;
        lengths.max_fragment_length = 14;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_left_only_window_is_tile_size_invariant() -> Result<()> {
    // Arrange:
    // - unpaired 2S10M2S at pos 10
    // - BED row [8,9) is reachable only through raw left clipping
    // - tile_size=10 forces a cross-tile situation; tile_size=1000 does not
    // The final motif output must be identical across both decompositions.
    let bam = single_read_bam(
        "ends_raw_left_only_tile_invariance",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let windows_bed = out_dir.path().join(format!("windows_{tile_size}.bed"));
        write_bed(&windows_bed, &[("chr1", 8, 9, "left_raw_only")])?;

        let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
        cfg.source_inside = KmerSource::Read;
        cfg.all_motifs = true;
        cfg.tile_size = tile_size;
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed),
        });
        cfg.set_window_assignment(AssignMotifToWindowArgs {
            assign_by: WindowMotifAssigner::Endpoint,
        });
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = 14;
            lengths.max_fragment_length = 14;
        }

        run(&cfg)?;
        outputs.push(read_dense_output(out_dir.path())?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[0].1.shape(), &[1, 4]);
    assert_eq!(motif_count(&outputs[0].1, &outputs[0].0, 0, "_T"), 1.0);
    assert_eq!(outputs[0].1.sum(), 1.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_assignment_does_not_count_a_window_ending_at_the_left_raw_boundary()
-> Result<()> {
    // Arrange: the same 2S10M2S read has raw left boundary 8, so BED window [7,8) touches the
    // boundary but does not contain the left endpoint under half-open semantics.
    let bam = single_read_bam(
        "ends_raw_left_boundary_open",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 7, 8, "touches_left_only")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 14;
        lengths.max_fragment_length = 14;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn aligned_endpoint_assignment_ignores_raw_shifted_boundary_positions() -> Result<()> {
    // Arrange: the same 2S10M2S read uses aligned endpoints [10,19] under aligned clipping, so
    // windows at the raw-shifted positions [8,9) and [21,22) should receive no counts.
    let bam = single_read_bam(
        "ends_aligned_not_raw",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 8, 9, "raw_left"), ("chr1", 21, 22, "raw_right")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Aligned;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[2, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_assignment_uses_aligned_positions_with_raw_bases() -> Result<()> {
    // Arrange: unpaired 2S10M2S at pos 10.
    // - aligned interval [10,20)
    // - raw-aligned-boundary keeps aligned endpoint positions 10 and 19
    // - raw inside bases still use the clipped terminal bases, so both ends orient to "_T"
    let bam = single_read_bam(
        "ends_raw_aligned_endpoint",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[
            ("chr1", 10, 11, "left_aligned"),
            ("chr1", 19, 20, "right_aligned"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[2, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_T"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_assignment_does_not_count_windows_at_shifted_positions()
-> Result<()> {
    // Arrange: the same 2S10M2S read keeps aligned endpoint positions [10,19] in
    // raw-aligned-boundary mode, so windows at the shifted positions [8,9) and [21,22) must stay
    // empty even though the inside bases come from the raw read.
    let bam = single_read_bam(
        "ends_raw_aligned_not_shifted",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[
            ("chr1", 8, 9, "shifted_left"),
            ("chr1", 21, 22, "shifted_right"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[2, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_left_only_window_outside_aligned_reach_is_tile_size_invariant()
-> Result<()> {
    // Arrange:
    // - unpaired 2S10M2S at pos 10 keeps aligned endpoint positions [10,19]
    // - BED row [8,9) is therefore outside aligned reach and must stay empty
    // - tile decomposition must not change that
    let bam = single_read_bam(
        "ends_raw_aligned_left_only_tile_invariance",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let windows_bed = out_dir.path().join(format!("windows_{tile_size}.bed"));
        write_bed(&windows_bed, &[("chr1", 8, 9, "shifted_left_only")])?;

        let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
        cfg.source_inside = KmerSource::Read;
        cfg.all_motifs = true;
        cfg.tile_size = tile_size;
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed),
        });
        cfg.set_window_assignment(AssignMotifToWindowArgs {
            assign_by: WindowMotifAssigner::Endpoint,
        });
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = 10;
            lengths.max_fragment_length = 10;
        }

        run(&cfg)?;
        outputs.push(read_dense_output(out_dir.path())?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[0].1.shape(), &[1, 4]);
    assert_eq!(outputs[0].1.sum(), 0.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_assignment_keeps_an_aligned_right_halo_window() -> Result<()> {
    // Arrange:
    // - unpaired 10M10S at pos 19 has aligned interval [19,29)
    // - raw-aligned-boundary keeps the right endpoint at 28
    // - the right inside base is raw read `A`, which reverse-complements to motif `_T`
    let bam = single_read_bam(
        "ends_raw_aligned_right_halo",
        19,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 28, 29, "right_aligned_only")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_assignment_does_not_count_a_far_right_window_beyond_aligned_reach()
-> Result<()> {
    // Arrange: unpaired 10M10S at pos 10 keeps aligned endpoint positions [10,19] in
    // raw-aligned-boundary mode, so BED row [29,30) must stay empty.
    let bam = single_read_bam(
        "ends_raw_aligned_far_right_window",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 29, 30, "far_right_only")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn aligned_endpoint_assignment_keeps_a_right_halo_only_window_reached_by_an_owned_fragment()
-> Result<()> {
    // Arrange:
    // - unpaired 10M at pos 19 has aligned interval [19,29)
    // - fragment ownership is by aligned start in tile core [10,20)
    // - BED row [28,29) is outside the core but contains the aligned right endpoint
    let bam = single_read_bam(
        "ends_aligned_right_halo",
        19,
        vec![('M', 10)],
        b"AAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 28, 29, "right_halo_only")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Aligned;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.row(0).sum(), 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn aligned_endpoint_assignment_mixed_core_and_right_halo_rows_count_only_the_true_target()
-> Result<()> {
    // Arrange:
    // - unpaired 10M at pos 19 has aligned endpoints at 19 and 28
    // - row [10,11) is a non-target core row
    // - row [28,29) is the true right-end target
    let bam = single_read_bam(
        "ends_aligned_mixed_core_and_halo",
        19,
        vec![('M', 10)],
        b"AAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[
            ("chr1", 10, 11, "core_only"),
            ("chr1", 28, 29, "right_halo_only"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Aligned;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    assert_eq!(matrix.shape(), &[2, 4]);
    assert_eq!(matrix.row(0).sum(), 0.0);
    assert_eq!(matrix.row(1).sum(), 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn skip_endpoint_assignment_keeps_a_right_halo_only_window_when_no_end_is_soft_clipped()
-> Result<()> {
    // Arrange: in skip mode an unclipped fragment should match aligned endpoint assignment.
    let bam = single_read_bam("ends_skip_right_halo", 19, vec![('M', 10)], b"AAAAAAAAAA")?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 28, 29, "right_halo_only")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Skip;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.row(0).sum(), 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn fragment_length_filters_use_the_adjusted_assignment_length_in_raw_shifted_boundary_mode()
-> Result<()> {
    // Arrange: the same 2S10M2S read has aligned length 10 but raw assignment length 14.
    // Keeping only length 10 should therefore exclude the fragment.
    let bam = single_read_bam(
        "ends_adjusted_length_filter",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: the fragment is excluded because its raw assignment length is 14 bp.
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn fragment_length_filters_use_the_aligned_assignment_length_in_raw_aligned_boundary_mode()
-> Result<()> {
    // Arrange: the same 2S10M2S read has aligned assignment length 10 in raw-aligned-boundary
    // mode, so keeping only length 10 should retain the fragment.
    let bam = single_read_bam(
        "ends_raw_aligned_length_filter",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn skip_clipping_skips_only_the_clipped_end_and_keeps_the_unclipped_end() -> Result<()> {
    // Arrange: 2S10M at pos 10 has a clipped left end and an unclipped right end.
    // With skip clipping, only the right end should remain. The aligned terminal base is T,
    // which orients to right-end label "_A".
    let bam = single_read_bam(
        "ends_skip_clipped_end",
        10,
        vec![('S', 2), ('M', 10)],
        b"TTAAAAAAAAAT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Skip;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn skip_clipping_skips_the_fragment_when_both_ends_are_soft_clipped() -> Result<()> {
    // Arrange: 2S10M2S has soft clipping on both ends, so skip clipping should leave no motif.
    let bam = single_read_bam(
        "ends_skip_both_clipped",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Skip;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn max_soft_clips_skips_only_the_over_clipped_end() -> Result<()> {
    // Arrange: 2S10M has two soft-clipped bases on the left end and none on the right end.
    // With max_soft_clips=1, the left end should be skipped while the unclipped right end remains.
    let bam = single_read_bam(
        "ends_max_soft_clips",
        10,
        vec![('S', 2), ('M', 10)],
        b"TTAAAAAAAAAT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.clip.max_soft_clips = 1;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_blacklist_validation_ignores_inside_bases_without_reference_overlap()
-> Result<()> {
    // Arrange: left motif uses two fully clipped bases in raw-aligned-boundary mode.
    // Blacklisting [8,10) covers only those clipped-only positions and must therefore not skip
    // the left endpoint motif, because blacklist validation only applies where the motif overlaps
    // reference.
    let bam = single_read_bam(
        "ends_raw_aligned_blacklist_clipped_prefix",
        10,
        vec![('S', 2), ('M', 10)],
        b"TTAAAAAAAAAT",
    )?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let blacklist_bed = out_dir.path().join("blacklist.bed");
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&blacklist_bed, &[("chr1", 8, 10, "clipped_prefix_only")])?;
    write_bed(&windows_bed, &[("chr1", 10, 11, "left_endpoint")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 2, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    cfg.blacklist = Some(vec![blacklist_bed]);
    cfg.blacklist_strategy = BlacklistStrategy::All;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["_TT"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 1.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_forbids_reference_inside_source() -> Result<()> {
    // Arrange
    let bam = single_read_bam(
        "ends_raw_aligned_reference_forbidden",
        10,
        vec![('S', 2), ('M', 10)],
        b"TTAAAAAAAAAT",
    )?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Reference;

    // Act
    let err = run(&cfg).expect_err("raw-aligned-boundary + reference source should fail");

    // Assert
    assert!(err.to_string().contains(
        "`--clip-strategy raw-aligned-boundary` cannot be combined with `--source-inside reference`"
    ));
    Ok(())
}

#[test]
fn max_soft_clips_keeps_a_raw_shifted_boundary_end_when_the_clip_count_equals_the_threshold()
-> Result<()> {
    // Arrange: with max_soft_clips=2, a 2S10M read should still keep the clipped left end because
    // the documented rule is "higher number of soft-clipped bases than this".
    let bam = single_read_bam(
        "ends_max_soft_clips_equal",
        10,
        vec![('S', 2), ('M', 10)],
        b"TTAAAAAAAAAT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.clip.max_soft_clips = 2;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 12;
        lengths.max_fragment_length = 12;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: both the left raw-clipped end ("_T") and the unclipped right end ("_A") survive.
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_A"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn max_soft_clips_skips_the_fragment_when_both_ends_exceed_the_threshold() -> Result<()> {
    // Arrange: 2S10M2S has two soft-clipped bases on both ends, so max_soft_clips=1 should leave
    // no surviving motifs.
    let bam = single_read_bam(
        "ends_max_soft_clips_both",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.clip.max_soft_clips = 1;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn hard_clipped_fragments_are_discarded_entirely() -> Result<()> {
    // Arrange: hard clipping is documented as an always-on fragment exclusion.
    let bam = single_read_bam(
        "ends_hard_clip",
        10,
        vec![('H', 2), ('M', 10)],
        b"AAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn motif_labels_use_outside_inside_order_when_outside_bases_are_present() -> Result<()> {
    // Arrange: left endpoint only for fragment [10,20) with one outside and one inside reference base.
    // - outside base before left boundary: seq[9]  = C
    // - inside base at left boundary:      seq[10] = G
    // So the final user-facing label should be "C_G".
    let bam = simple_paired_fragment_bam("ends_label_order", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 11, "left")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["C_G"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 1.0);
    Ok(())
}

#[test]
fn right_end_motif_labels_use_outside_inside_order_when_outside_bases_are_present() -> Result<()> {
    // Arrange: right endpoint only for fragment [4,14) with one outside and one inside
    // reference base chosen so the right-end decode is non-palindromic.
    //
    // Reference around the right boundary:
    // - inside base at fragment end-1: seq[13] = T
    // - outside base after fragment:   seq[14] = G
    //
    // Right-end storage order is "TG", and revcomp("TG") = "CA", so the final public label must
    // be "C_A".
    let reference = twobit_from_sequences(
        "ends_right_label_order_reference",
        vec![(
            "chr1".to_string(),
            format!("{}{}{}", "T".repeat(13), "TG", "T".repeat(241)),
        )],
    )?;
    let bam = single_read_bam("ends_right_label_order", 4, vec![('M', 10)], b"AAAAAAAAAA")?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 13, 14, "right")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["C_A"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 1.0);
    Ok(())
}

#[test]
fn outside_only_motifs_are_labeled_with_an_empty_inside_half() -> Result<()> {
    // Arrange: left endpoint only with k_inside=0 and k_outside=1.
    // The base immediately outside the left boundary of [10,20) is seq[9] = C, so the label must
    // be "C_".
    let bam = simple_paired_fragment_bam("ends_outside_only", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 11, "left")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 0, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["C_"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 1.0);
    Ok(())
}

#[test]
fn left_edge_missing_outside_context_drops_the_left_endpoint_motif() -> Result<()> {
    // Arrange: fragment [0,10) has no reference base outside the left boundary, so with
    // k_outside=1 the left endpoint motif should decode to a sentinel and be dropped.
    let bam = simple_paired_fragment_bam("ends_left_edge_sentinel", 0, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 0, 1, "left"), ("chr1", 9, 10, "right")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: the left endpoint row is empty, while the right endpoint still contributes one motif.
    assert_eq!(matrix.row(0).sum(), 0.0);
    assert_eq!(matrix.row(1).sum(), 1.0);
    Ok(())
}

#[test]
fn right_edge_missing_outside_context_drops_the_right_endpoint_motif() -> Result<()> {
    // Arrange: fragment [246,256) ends at the chromosome boundary, so with k_outside=1 the right
    // endpoint has no outside reference context and should be dropped.
    let bam = simple_paired_fragment_bam("ends_right_edge_sentinel", 246, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 246, 247, "left"), ("chr1", 255, 256, "right")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: the right endpoint row is empty, while the left endpoint still contributes one motif.
    assert_eq!(matrix.row(0).sum(), 1.0);
    assert_eq!(matrix.row(1).sum(), 0.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_shifting_still_applies_when_only_outside_bases_are_counted() -> Result<()> {
    // Arrange: with k_inside=0 and raw clipping, the shifted raw boundary should still control
    // endpoint assignment and outside-base extraction.
    //
    // 2S10M2S at pos 10 gives raw assignment interval [8,22), so the left endpoint is 8.
    // The base immediately outside that raw left boundary is seq[7] on the reference, which is T.
    let bam = single_read_bam(
        "ends_raw_outside_only",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 8, 9, "left_raw")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 0, 1);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 14;
        lengths.max_fragment_length = 14;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["T_"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 1.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_outside_only_motifs_use_the_aligned_boundary() -> Result<()> {
    // Arrange: with k_inside=0 and raw-aligned-boundary clipping, outside-base extraction should
    // still use the aligned boundary.
    //
    // 2S10M2S at pos 10 keeps the left endpoint at 10, so the base immediately outside is
    // seq[9] = C in the simple ACGT reference.
    let bam = single_read_bam(
        "ends_raw_aligned_outside_only",
        10,
        vec![('S', 2), ('M', 10), ('S', 2)],
        b"TTAAAAAAAAAAAA",
    )?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 11, "left_aligned")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 0, 1);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["C_"]);
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(matrix[(0, 0)], 1.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_assignment_keeps_a_far_right_window_beyond_aligned_reach()
-> Result<()> {
    // Arrange: unpaired 10M10S at pos 10 with tile size 10.
    //
    // Mental derivation:
    // - aligned interval is [10,20), so the fragment again belongs to tile core [10,20)
    // - raw assignment interval is [10,30), so the right endpoint is 29
    // - BED window [29,30) is outside the tile core and outside the aligned-length-only right
    //   reach from that core, but it is still a valid raw endpoint window
    // - the right raw base is A, which reverse-complements to "_T" on decode
    let bam = single_read_bam(
        "ends_raw_far_right_window",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 29, 30, "right_raw_only")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 20;
        lengths.max_fragment_length = 20;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_far_right_only_window_is_tile_size_invariant() -> Result<()> {
    // Arrange:
    // - unpaired 10M10S at pos 10
    // - BED row [29,30) is reachable only through raw right clipping
    // - tile_size=10 and tile_size=1000 must therefore produce the same final motif row
    let bam = single_read_bam(
        "ends_raw_far_right_tile_invariance",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let windows_bed = out_dir.path().join(format!("windows_{tile_size}.bed"));
        write_bed(&windows_bed, &[("chr1", 29, 30, "right_raw_only")])?;

        let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
        cfg.source_inside = KmerSource::Read;
        cfg.all_motifs = true;
        cfg.tile_size = tile_size;
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed),
        });
        cfg.set_window_assignment(AssignMotifToWindowArgs {
            assign_by: WindowMotifAssigner::Endpoint,
        });
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = 20;
            lengths.max_fragment_length = 20;
        }

        run(&cfg)?;
        outputs.push(read_dense_output(out_dir.path())?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[0].1.shape(), &[1, 4]);
    assert_eq!(motif_count(&outputs[0].1, &outputs[0].0, 0, "_T"), 1.0);
    assert_eq!(outputs[0].1.sum(), 1.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_far_right_only_window_is_tile_size_invariant() -> Result<()> {
    // Arrange:
    // - unpaired 10M10S at pos 10 keeps aligned endpoint positions [10,19]
    // - BED row [29,30) is therefore always outside aligned reach
    let bam = single_read_bam(
        "ends_raw_aligned_far_right_tile_invariance",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let windows_bed = out_dir.path().join(format!("windows_{tile_size}.bed"));
        write_bed(&windows_bed, &[("chr1", 29, 30, "far_right_only")])?;

        let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
        cfg.source_inside = KmerSource::Read;
        cfg.all_motifs = true;
        cfg.tile_size = tile_size;
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed),
        });
        cfg.set_window_assignment(AssignMotifToWindowArgs {
            assign_by: WindowMotifAssigner::Endpoint,
        });
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = 10;
            lengths.max_fragment_length = 10;
        }

        run(&cfg)?;
        outputs.push(read_dense_output(out_dir.path())?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[0].1.shape(), &[1, 4]);
    assert_eq!(outputs[0].1.sum(), 0.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_assignment_does_not_count_a_window_starting_at_the_right_raw_boundary()
-> Result<()> {
    // Arrange: 10M10S at pos 10 has raw assignment interval [10,30), so the counted right-end
    // position is 29. BED window [30,31) starts exactly at the exclusive boundary and must remain
    // empty under half-open semantics.
    let bam = single_read_bam(
        "ends_raw_right_boundary_open",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 30, 31, "touches_right_only")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 20;
        lengths.max_fragment_length = 20;
    }

    // Act
    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[1, 4]);
    assert_eq!(matrix.sum(), 0.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_assignment_must_not_shrink_fetch_to_unrelated_core_windows()
-> Result<()> {
    // Arrange: unpaired 10M10S at pos 19 with tile size 10.
    //
    // Mental derivation:
    // - aligned interval is [19,29), so the fragment is owned by tile core [10,20)
    // - raw assignment interval is [19,39), so the right endpoint is 38
    // - BED window [10,11) overlaps the core but not either endpoint
    // - BED window [38,39) is the true target window for the right raw endpoint
    // The correct output is therefore two rows with counts [0.0, 1.0].
    let bam = single_read_bam(
        "ends_raw_far_right_with_core_window",
        19,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[
            ("chr1", 10, 11, "core_only"),
            ("chr1", 38, 39, "right_raw_only"),
        ],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 20;
        lengths.max_fragment_length = 20;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[2, 4]);
    assert_eq!(matrix.row(0).sum(), 0.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_T"), 1.0);
    assert_eq!(matrix.sum(), 1.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_mixed_core_and_far_right_rows_are_tile_size_invariant()
-> Result<()> {
    // Arrange:
    // - unpaired 10M10S at pos 19
    // - BED row [10,11) is a non-target core row and must stay zero
    // - BED row [38,39) is the true raw right-end target and must carry one motif
    // - with tile_size=10 these rows live in different tiles; with tile_size=1000 they do not
    let bam = single_read_bam(
        "ends_raw_mixed_tile_invariance",
        19,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let windows_bed = out_dir.path().join(format!("windows_{tile_size}.bed"));
        write_bed(
            &windows_bed,
            &[
                ("chr1", 10, 11, "core_only"),
                ("chr1", 38, 39, "right_raw_only"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
        cfg.source_inside = KmerSource::Read;
        cfg.all_motifs = true;
        cfg.tile_size = tile_size;
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed),
        });
        cfg.set_window_assignment(AssignMotifToWindowArgs {
            assign_by: WindowMotifAssigner::Endpoint,
        });
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = 20;
            lengths.max_fragment_length = 20;
        }

        run(&cfg)?;
        outputs.push(read_dense_output(out_dir.path())?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[0].1.shape(), &[2, 4]);
    assert_eq!(outputs[0].1.row(0).sum(), 0.0);
    assert_eq!(motif_count(&outputs[0].1, &outputs[0].0, 1, "_T"), 1.0);
    assert_eq!(outputs[0].1.sum(), 1.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_mixed_core_and_far_right_rows_are_tile_size_invariant()
-> Result<()> {
    // Arrange:
    // - unpaired 10M10S at pos 19 keeps aligned endpoint positions [19,28]
    // - BED rows [10,11) and [38,39) are both non-target rows and must stay zero
    let bam = single_read_bam(
        "ends_raw_aligned_mixed_tile_invariance",
        19,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let windows_bed = out_dir.path().join(format!("windows_{tile_size}.bed"));
        write_bed(
            &windows_bed,
            &[
                ("chr1", 10, 11, "core_only"),
                ("chr1", 38, 39, "far_right_only"),
            ],
        )?;

        let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
        cfg.source_inside = KmerSource::Read;
        cfg.all_motifs = true;
        cfg.tile_size = tile_size;
        cfg.set_windows(WindowsArgs {
            by_size: None,
            by_bed: Some(windows_bed),
        });
        cfg.set_window_assignment(AssignMotifToWindowArgs {
            assign_by: WindowMotifAssigner::Endpoint,
        });
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = 10;
            lengths.max_fragment_length = 10;
        }

        run(&cfg)?;
        outputs.push(read_dense_output(out_dir.path())?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[0].1.shape(), &[2, 4]);
    assert_eq!(outputs[0].1.row(0).sum(), 0.0);
    assert_eq!(outputs[0].1.row(1).sum(), 0.0);
    assert_eq!(outputs[0].1.sum(), 0.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_by_size_counts_the_previous_bin_reached_by_left_raw_clipping()
-> Result<()> {
    // Arrange:
    // - unpaired 2S10M at pos 10 has raw endpoints at 8 and 19
    // - fixed-size windows of 10 bp therefore put the left endpoint in [0,10)
    //   and the right endpoint in [10,20)
    // - ownership still comes from aligned start in tile core [10,20)
    let bam = single_read_bam(
        "ends_raw_by_size_left_previous_bin",
        10,
        vec![('S', 2), ('M', 10)],
        b"TTAAAAAAAAAT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: Some(10),
        by_bed: None,
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 12;
        lengths.max_fragment_length = 12;
    }

    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    assert_eq!(matrix.row(0).sum(), 1.0);
    assert_eq!(matrix.row(1).sum(), 1.0);
    assert_eq!(matrix.row(2).sum(), 0.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_by_size_keeps_both_ends_in_the_aligned_bin_with_left_clipping()
-> Result<()> {
    // Arrange:
    // - unpaired 2S10M at pos 10 keeps aligned endpoint positions 10 and 19
    // - fixed-size windows of 10 bp therefore keep both endpoints in [10,20)
    let bam = single_read_bam(
        "ends_raw_aligned_by_size_left",
        10,
        vec![('S', 2), ('M', 10)],
        b"TTAAAAAAAAAT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: Some(10),
        by_bed: None,
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    assert_eq!(matrix.row(0).sum(), 0.0);
    assert_eq!(matrix.row(1).sum(), 2.0);
    assert_eq!(matrix.row(2).sum(), 0.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_by_size_counts_the_next_bin_reached_by_right_raw_clipping()
-> Result<()> {
    // Arrange:
    // - unpaired 10M10S at pos 10 has raw endpoints at 10 and 29
    // - fixed-size windows of 10 bp therefore put the left endpoint in [10,20)
    //   and the right endpoint in [20,30)
    let bam = single_read_bam(
        "ends_raw_by_size_right_next_bin",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: Some(10),
        by_bed: None,
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 20;
        lengths.max_fragment_length = 20;
    }

    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    assert_eq!(matrix.row(0).sum(), 0.0);
    assert_eq!(matrix.row(1).sum(), 1.0);
    assert_eq!(matrix.row(2).sum(), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_by_size_keeps_both_ends_in_the_aligned_bin_with_right_clipping()
-> Result<()> {
    // Arrange:
    // - unpaired 10M10S at pos 10 keeps aligned endpoint positions 10 and 19
    // - fixed-size windows of 10 bp therefore keep both endpoints in [10,20)
    let bam = single_read_bam(
        "ends_raw_aligned_by_size_right",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.tile_size = 10;
    cfg.set_windows(WindowsArgs {
        by_size: Some(10),
        by_bed: None,
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    run(&cfg)?;
    let (_motifs, matrix) = read_dense_output(out_dir.path())?;

    assert_eq!(matrix.row(0).sum(), 0.0);
    assert_eq!(matrix.row(1).sum(), 2.0);
    assert_eq!(matrix.row(2).sum(), 0.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_by_size_keeps_exact_half_open_boundary_bins_zero() -> Result<()> {
    // Arrange:
    // - 2S10M at pos 10 has raw left endpoint 8, so [7,8) must stay zero while [8,9) counts
    // - 10M10S at pos 10 has raw right endpoint 29, so [30,31) must stay zero while [29,30) counts
    let left_bam = single_read_bam(
        "ends_raw_by_size_left_boundary_open",
        10,
        vec![('S', 2), ('M', 10)],
        b"TTAAAAAAAAAT",
    )?;
    let right_bam = single_read_bam(
        "ends_raw_by_size_right_boundary_open",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let left_out = TempDir::new()?;
    let right_out = TempDir::new()?;

    let make_cfg = |bam_path: &Path, out_dir: &Path, fragment_length: u32| {
        let mut cfg = base_config(bam_path, out_dir, 1, 0);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
        cfg.source_inside = KmerSource::Read;
        cfg.all_motifs = true;
        cfg.tile_size = 10;
        cfg.set_windows(WindowsArgs {
            by_size: Some(1),
            by_bed: None,
        });
        cfg.set_window_assignment(AssignMotifToWindowArgs {
            assign_by: WindowMotifAssigner::Endpoint,
        });
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = fragment_length;
            lengths.max_fragment_length = fragment_length;
        }
        cfg
    };

    run(&make_cfg(&left_bam.bam, left_out.path(), 12))?;
    run(&make_cfg(&right_bam.bam, right_out.path(), 20))?;

    let (_left_motifs, left_matrix) = read_dense_output(left_out.path())?;
    let (_right_motifs, right_matrix) = read_dense_output(right_out.path())?;

    assert_eq!(left_matrix.row(7).sum(), 0.0);
    assert_eq!(left_matrix.row(8).sum(), 1.0);
    assert_eq!(right_matrix.row(29).sum(), 1.0);
    assert_eq!(right_matrix.row(30).sum(), 0.0);
    Ok(())
}

#[test]
fn raw_shifted_boundary_endpoint_by_size_output_is_tile_size_invariant() -> Result<()> {
    // Arrange:
    // - fixed-size windows are 10 bp bins
    // - changing tile_size must not change the final by-size endpoint rows
    let bam = single_read_bam(
        "ends_raw_by_size_tile_invariance",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.clip.clip_strategy = ClipStrategy::RawShiftedBoundary;
        cfg.source_inside = KmerSource::Read;
        cfg.all_motifs = true;
        cfg.tile_size = tile_size;
        cfg.set_windows(WindowsArgs {
            by_size: Some(10),
            by_bed: None,
        });
        cfg.set_window_assignment(AssignMotifToWindowArgs {
            assign_by: WindowMotifAssigner::Endpoint,
        });
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = 20;
            lengths.max_fragment_length = 20;
        }

        run(&cfg)?;
        outputs.push(read_dense_output(out_dir.path())?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[0].1.row(1).sum(), 1.0);
    assert_eq!(outputs[0].1.row(2).sum(), 1.0);
    assert_eq!(outputs[0].1.sum(), 2.0);
    Ok(())
}

#[test]
fn raw_aligned_boundary_endpoint_by_size_output_is_tile_size_invariant() -> Result<()> {
    // Arrange:
    // - unpaired 10M10S at pos 10 keeps aligned endpoint positions [10,19]
    // - with 10 bp bins, both endpoints therefore stay in row 1 regardless of tile size
    let bam = single_read_bam(
        "ends_raw_aligned_by_size_tile_invariance",
        10,
        vec![('M', 10), ('S', 10)],
        b"AAAAAAAAAAAAAAAAAAAA",
    )?;
    let tile_sizes = [10_u32, 1_000_u32];
    let mut outputs = Vec::new();

    for tile_size in tile_sizes {
        let out_dir = TempDir::new()?;
        let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
        cfg.set_unpaired(UnpairedArgs {
            reads_are_fragments: true,
        });
        cfg.clip.clip_strategy = ClipStrategy::RawAlignedBoundary;
        cfg.source_inside = KmerSource::Read;
        cfg.all_motifs = true;
        cfg.tile_size = tile_size;
        cfg.set_windows(WindowsArgs {
            by_size: Some(10),
            by_bed: None,
        });
        cfg.set_window_assignment(AssignMotifToWindowArgs {
            assign_by: WindowMotifAssigner::Endpoint,
        });
        {
            let lengths = cfg.fragment_lengths_mut();
            lengths.min_fragment_length = 10;
            lengths.max_fragment_length = 10;
        }

        run(&cfg)?;
        outputs.push(read_dense_output(out_dir.path())?);
    }

    assert_eq!(outputs[0], outputs[1]);
    assert_eq!(outputs[0].1.row(0).sum(), 0.0);
    assert_eq!(outputs[0].1.row(1).sum(), 2.0);
    assert_eq!(outputs[0].1.row(2).sum(), 0.0);
    assert_eq!(outputs[0].1.sum(), 2.0);
    Ok(())
}

#[test]
fn both_kmer_sizes_zero_is_rejected() -> Result<()> {
    // Arrange: an empty motif is intentionally undefined.
    let bam = simple_paired_fragment_bam("ends_empty_motif", 10, 10, 4)?;
    let out_dir = TempDir::new()?;
    let cfg = base_config(&bam.bam, out_dir.path(), 0, 0);

    // Act
    let err = run(&cfg).expect_err("empty motif definition should be rejected");

    // Assert
    assert!(
        err.to_string()
            .contains("At least one of --k-inside or --k-outside must be > 0")
    );
    Ok(())
}

#[test]
fn settings_json_keeps_the_runtime_fields_needed_to_interpret_output() -> Result<()> {
    // Arrange: this run changes only fields that still belong in the sidecar contract:
    // - source_inside = read
    // - clip_strategy = skip
    // - window_assignment = endpoint
    // - collapse_complement = false
    let bam = single_read_bam(
        "ends_settings_semantics",
        10,
        vec![('S', 2), ('M', 10)],
        b"TTAAAAAAAAAT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Skip;
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let settings = read_text_file(&settings_path(out_dir.path()))?;

    // Assert
    assert_eq!(
        parse_json(&settings),
        expected_settings_json("read", "skip", "endpoint")
    );
    Ok(())
}

#[test]
fn settings_json_formats_proportion_window_assignment_stably() -> Result<()> {
    // Arrange: use a simple exact decimal that should stay readable in the sidecar.
    let bam = single_read_bam(
        "ends_settings_proportion_precision",
        10,
        vec![('M', 10)],
        b"AAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Proportion(0.125),
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let settings = read_text_file(&settings_path(out_dir.path()))?;

    // Assert: exact sidecar contract, not just substring presence.
    assert_eq!(
        parse_json(&settings),
        expected_settings_json("read", "aligned", "proportion=0.125")
    );
    Ok(())
}

#[test]
fn reference_backed_inside_bases_require_ref_2bit() -> Result<()> {
    // Arrange: reference-backed inside extraction is documented to require --ref-2bit.
    let bam = simple_paired_fragment_bam("ends_missing_ref_2bit", 10, 10, 4)?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    let err = run(&cfg).expect_err("reference-backed ends should require --ref-2bit");

    // Assert
    assert!(err.to_string().contains("--ref-2bit"));
    Ok(())
}

#[test]
fn outside_bases_require_ref_2bit_even_when_inside_bases_come_from_reads() -> Result<()> {
    // Arrange: outside motifs always come from the reference, so k_outside>0 still requires --ref-2bit.
    let bam = simple_paired_fragment_bam("ends_missing_ref_for_outside", 10, 10, 4)?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 1);
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    let err = run(&cfg).expect_err("outside bases should require --ref-2bit");

    // Assert
    assert!(err.to_string().contains("--ref-2bit"));
    Ok(())
}

#[test]
fn scaling_factors_must_cover_every_counted_fragment() -> Result<()> {
    // Arrange: fragment [10,20) is counted, but the scaling file only covers [0,10), so there
    // is no overlapping scaling bin at all for that fragment.
    let bam = simple_paired_fragment_bam("ends_scaling_gap", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let scaling_path = out_dir.path().join("scaling.tsv");
    std::fs::write(
        &scaling_path,
        "chromosome\tstart\tend\tscaling_factor\nchr1\t0\t10\t2\n",
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_scaling_factors(Some(scaling_path));
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    let err = run(&cfg).expect_err("incomplete scaling coverage should fail loudly");

    // Assert
    assert!(
        format!("{err:#}")
            .contains("scaling TSV: bins on 'chr1' must end at chrom_len=256 (got end=10)")
    );
    Ok(())
}

#[test]
fn windowed_runs_write_bins_tsv_with_the_selected_windows() -> Result<()> {
    // Arrange: in BED-windowed mode the command should persist the selected windows as TSV.
    let bam = simple_paired_fragment_bam("ends_bins_tsv", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 10, 11, "left"), ("chr1", 19, 20, "right")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let bins_tsv = read_text_file(&out_dir.path().join("ends.bins.tsv"))?;

    // Assert: header plus one row per selected window.
    let rows: Vec<&str> = bins_tsv.lines().collect();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0], "chrom\tstart\tend\tblacklisted_fraction");
    assert!(rows[1].starts_with("chr1\t10\t11\t"));
    assert!(rows[2].starts_with("chr1\t19\t20\t"));
    Ok(())
}

#[test]
fn settings_json_ignores_fragment_length_bounds_but_keeps_motif_definition_fields() -> Result<()> {
    // Arrange:
    // - fragment-length bounds change counting eligibility only, so they should not appear
    // - source_inside stays reference
    // - clip_strategy stays aligned
    // - window_assignment stays endpoint
    // - collapse_complement stays false
    let bam = simple_paired_fragment_bam("ends_settings_lengths", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 11;
    }

    // Act
    run(&cfg)?;
    let settings = read_text_file(&settings_path(out_dir.path()))?;

    // Assert
    assert_eq!(
        parse_json(&settings),
        expected_settings_json("reference", "aligned", "endpoint")
    );
    Ok(())
}

#[test]
fn unpaired_mode_rejects_require_proper_pair() -> Result<()> {
    // Arrange: the command explicitly forbids combining reads-as-fragments with proper-pair filtering.
    let bam = single_read_bam(
        "ends_unpaired_proper_pair",
        10,
        vec![('M', 10)],
        b"AAAAAAAAAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.require_proper_pair = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    let err = run(&cfg).expect_err("reads-as-fragments cannot require proper pairing");

    // Assert
    assert!(
        err.to_string()
            .contains("--require-proper-pair cannot be used with --reads-are-fragments")
    );
    Ok(())
}

#[test]
fn sparse_output_motif_labels_only_include_observed_motifs() -> Result<()> {
    // Arrange: one observed motif pair should not force unobserved motifs into the sparse label file.
    let bam = simple_paired_fragment_bam("ends_sparse_observed_only", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let motifs = read_motif_labels(&sparse_output_paths(out_dir.path()).1)?;

    // Assert
    assert_eq!(motifs, vec!["_A", "_G"]);
    Ok(())
}

#[test]
fn all_motifs_dense_output_enumerates_inside_only_labels_when_k_outside_is_zero() -> Result<()> {
    // Arrange: inside-only motifs should still have a fixed dense universe under all-motifs.
    let bam = simple_paired_fragment_bam("ends_inside_only_dense", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, _matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["_A", "_C", "_G", "_T"]);
    Ok(())
}

#[test]
fn read_backed_inside_only_runs_without_ref_2bit() -> Result<()> {
    // Arrange: when both the inside bases come from reads and k_outside=0, the command should not
    // require a reference genome.
    let bam = single_read_bam("ends_read_only_no_ref", 10, vec![('M', 10)], b"AAAAAAAAAA")?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;

    // Assert
    assert!(sparse_output_paths(out_dir.path()).0.exists());
    assert!(settings_path(out_dir.path()).exists());
    Ok(())
}

#[test]
fn output_prefix_is_applied_to_all_primary_end_outputs() -> Result<()> {
    // Arrange: output prefix should namespace all primary artifacts.
    let bam = simple_paired_fragment_bam("ends_prefixed_outputs", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.output_prefix = "sampleA".to_string();
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;

    // Assert
    assert!(
        out_dir
            .path()
            .join("sampleA.end_motifs.sparse.npz")
            .exists()
    );
    assert!(out_dir.path().join("sampleA.end_motifs.txt").exists());
    assert!(
        out_dir
            .path()
            .join("sampleA.end_motif_settings.json")
            .exists()
    );
    assert!(!out_dir.path().join("end_motifs.sparse.npz").exists());
    assert!(!out_dir.path().join("end_motifs.txt").exists());
    assert!(!out_dir.path().join("end_motif_settings.json").exists());
    Ok(())
}

#[test]
fn default_window_assignment_is_endpoint() -> Result<()> {
    // Arrange: without overriding assign-by, the documented default is endpoint.
    let bam = simple_paired_fragment_bam("ends_default_endpoint", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 10, 11, "left"), ("chr1", 19, 20, "right")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_windows(WindowsArgs {
        by_size: None,
        by_bed: Some(windows_bed),
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert
    assert_eq!(matrix.shape(), &[2, 4]);
    assert_eq!(motif_count(&matrix, &motifs, 0, "_G"), 1.0);
    assert_eq!(motif_count(&matrix, &motifs, 1, "_A"), 1.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn by_size_windowing_writes_bins_tsv() -> Result<()> {
    // Arrange: fixed-size windowing should also persist the resolved window coordinates.
    let bam = simple_paired_fragment_bam("ends_by_size_bins", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.set_windows(WindowsArgs {
        by_size: Some(20),
        by_bed: None,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;
    let bins_tsv = read_text_file(&out_dir.path().join("ends.bins.tsv"))?;

    // Assert
    let rows: Vec<&str> = bins_tsv.lines().collect();
    assert!(!rows.is_empty());
    assert_eq!(rows[0], "chrom\tstart\tend\tblacklisted_fraction");
    assert!(rows.iter().skip(1).all(|row| row.starts_with("chr1\t")));
    Ok(())
}

#[test]
fn output_prefix_is_applied_to_bins_tsv_for_windowed_runs() -> Result<()> {
    // Arrange: prefixed runs should namespace the auxiliary bins TSV too.
    let bam = simple_paired_fragment_bam("ends_prefixed_bins", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.output_prefix = "sampleA".to_string();
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.set_windows(WindowsArgs {
        by_size: Some(20),
        by_bed: None,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;

    // Assert
    assert!(out_dir.path().join("sampleA.bins.tsv").exists());
    assert!(!out_dir.path().join("bins.tsv").exists());
    Ok(())
}

#[test]
fn output_prefix_is_applied_to_dense_all_motifs_outputs() -> Result<()> {
    // Arrange: prefixed all-motifs runs should namespace the dense primary outputs too.
    let bam = simple_paired_fragment_bam("ends_prefixed_dense_outputs", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.output_prefix = "sampleA".to_string();
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;

    // Assert
    assert!(out_dir.path().join("sampleA.end_motifs.npy").exists());
    assert!(out_dir.path().join("sampleA.end_motifs.txt").exists());
    assert!(
        out_dir
            .path()
            .join("sampleA.end_motif_settings.json")
            .exists()
    );
    assert!(!out_dir.path().join("end_motifs.npy").exists());
    assert!(!out_dir.path().join("end_motifs.txt").exists());
    assert!(!out_dir.path().join("end_motif_settings.json").exists());
    Ok(())
}

#[test]
fn empty_output_prefix_writes_unprefixed_primary_outputs() -> Result<()> {
    // Arrange: the documented empty-prefix behavior is to write filenames without a leading prefix.
    let bam = simple_paired_fragment_bam("ends_empty_prefix", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.output_prefix.clear();
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;

    // Assert
    assert!(out_dir.path().join("end_motifs.sparse.npz").exists());
    assert!(out_dir.path().join("end_motifs.txt").exists());
    assert!(out_dir.path().join("end_motif_settings.json").exists());
    Ok(())
}

#[test]
fn empty_output_prefix_writes_unprefixed_bins_tsv_for_windowed_runs() -> Result<()> {
    // Arrange: the empty-prefix contract should also apply to auxiliary window outputs.
    let bam = simple_paired_fragment_bam("ends_empty_prefix_bins", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.output_prefix.clear();
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_inside = KmerSource::Reference;
    cfg.all_motifs = false;
    cfg.set_windows(WindowsArgs {
        by_size: Some(20),
        by_bed: None,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;

    // Assert
    assert!(out_dir.path().join("bins.tsv").exists());
    Ok(())
}

#[test]
fn read_backed_paired_inside_only_runs_without_ref_2bit() -> Result<()> {
    // Arrange: paired runs also should not require a reference when both the inside bases come
    // from reads and k_outside=0.
    let bam = simple_paired_fragment_bam("ends_paired_read_only_no_ref", 10, 10, 4)?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.source_inside = KmerSource::Read;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 10;
        lengths.max_fragment_length = 10;
    }

    // Act
    run(&cfg)?;

    // Assert
    assert!(sparse_output_paths(out_dir.path()).0.exists());
    assert!(settings_path(out_dir.path()).exists());
    Ok(())
}
