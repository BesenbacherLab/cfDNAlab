#![cfg(feature = "cmd_ends")]

mod fixtures;

use anyhow::{Context, Result};
#[cfg(feature = "cmd_gc_bias")]
use cfdnalab::commands::gc_bias::{GC_CORRECTION_SCHEMA_VERSION, package::GCCorrectionPackage};
use cfdnalab::commands::{
    cli_common::{ApplyGCArgs, ChromosomeArgs, IOCArgs, UnpairedArgs, WindowsArgs},
    ends::{
        config::EndsConfig,
        config_structs::{AssignMotifToWindowArgs, ClipStrategy, KmerSource, WindowMotifAssigner},
        ends::run,
    },
};
use cfdnalab::shared::{blacklist::BlacklistStrategy, indel_mode::IndelMotifFilterPolicy};
use fixtures::{
    BamFixture, FragmentSpec, ReadSpec, bam_from_specs, paired_fragment, simple_reference_twobit,
    twobit_from_sequences, write_bed,
};
#[cfg(feature = "cmd_gc_bias")]
use ndarray::array;
use ndarray::{Array1, Array2};
use ndarray_npy::{NpzReader, read_npy};
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
    k_within: usize,
    k_outside: usize,
) -> EndsConfig {
    let mut cfg = EndsConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        k_within,
        k_outside,
    );
    cfg.output_prefix = "ends".to_string();
    cfg.tile_size = 100;
    cfg.min_mapq = 0;
    cfg.require_proper_pair = false;
    cfg.clip.clip_strategy = ClipStrategy::Aligned;
    cfg
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

fn read_text_file(path: &Path) -> Result<String> {
    let mut buf = String::new();
    File::open(path)?.read_to_string(&mut buf)?;
    Ok(buf)
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
    cfg.source_within = KmerSource::Reference;
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
    // Blacklisting [10,11) masks only the left within base. Using blacklist_strategy=All keeps
    // the fragment itself, so only the left endpoint motif should disappear.
    let bam = simple_paired_fragment_bam("ends_blacklist_reference_end", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let blacklist_bed = out_dir.path().join("blacklist.bed");
    write_bed(&blacklist_bed, &[("chr1", 10, 11, "blk")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
            "--k-within",
            "1",
            "--k-outside",
            "0",
            "--ref-2bit",
            reference
                .path
                .to_str()
                .context("reference path is not valid UTF-8")?,
            "--source-within",
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
            out_dir.path().to_str().context("output dir is not valid UTF-8")?,
            "--chromosomes",
            "chr1",
            "--k-within",
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

#[test]
fn blacklist_masking_still_skips_read_backed_within_motifs_using_genomic_reference_coordinates()
-> Result<()> {
    // Arrange: unpaired read-fragment [10,14) with read sequence A C G A.
    // - left read-backed motif = "_A"
    // - right read-backed motif = reverse-complement("A") = "_T"
    // Blacklisting [10,11) should drop only the left motif even though within bases come from
    // the read, because blacklist validation is still genomic.
    let bam = single_read_bam("ends_blacklist_read_end", 10, vec![('M', 4)], b"ACGA")?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let blacklist_bed = out_dir.path().join("blacklist.bed");
    write_bed(&blacklist_bed, &[("chr1", 10, 11, "blk")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Read;
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
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
    cfg.all_motifs = true;
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        drop_invalid_gc: false,
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
fn drop_invalid_gc_skips_fragments_when_gc_correction_cannot_be_computed() -> Result<()> {
    // Arrange: use a reference where the fragment GC window contains only `N`, so GC fraction
    // cannot be computed even though the correction package covers the fragment length. With
    // drop_invalid_gc=true the fragment should be skipped instead of falling back to weight 1.0.
    let bam = simple_paired_fragment_bam("ends_drop_invalid_gc", 10, 10, 4)?;
    let reference = twobit_from_sequences(
        "ends_drop_invalid_gc_reference",
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
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = true;
    cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path),
        gc_tag: None,
        drop_invalid_gc: true,
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
    // In auto mode with read-backed within bases, that end should still be kept.
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
    cfg.source_within = KmerSource::Read;
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

    // Assert: both ends survive, so one global row should contain two observed motifs.
    assert_eq!(matrix.shape(), &[1, 2]);
    assert_eq!(motifs.len(), 2);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn auto_indel_filter_skips_indel_affected_reference_backed_end_motifs() -> Result<()> {
    // Arrange: the same pair as above, but now auto mode uses reference-backed within bases and
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
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 4, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
    assert_eq!(matrix.shape(), &[1, 1]);
    assert_eq!(motifs.len(), 1);
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
    cfg.source_within = KmerSource::Reference;
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
fn outside_reference_lookup_falls_back_to_exact_reference_fetch_when_the_motif_crosses_the_tile_slice()
-> Result<()> {
    // Arrange: by-BED window [10,11) with max_fragment_length=4 gives a tile-local fetch halo of
    // 4 bp, so the loaded slice starts at 6. Asking for k_outside=5 on the left endpoint at 10
    // needs reference bases [5,10), which crosses one base outside the preloaded tile slice and
    // must therefore use the exact fallback path. On the ACGT-repeat reference, seq[5..10) is
    // C G T A C, so the outside-only label is "CGTAC_".
    let bam = single_read_bam("ends_exact_reference_fallback", 10, vec![('M', 4)], b"ACGT")?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 11, "left")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 0, 5);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Read;
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
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_sparse_output(out_dir.path())?;

    // Assert
    assert_eq!(motifs, vec!["CGTAC_"]);
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
fn proportion_assignment_counts_the_fragment_when_the_requested_fraction_is_met() -> Result<()> {
    // Arrange: fragment [10,20) overlaps [10,15) by 5/10 bp, so proportion=0.5 should accept it.
    let bam = simple_paired_fragment_bam("ends_proportion", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 15, "half")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
fn collapse_complement_merges_reverse_complement_equivalent_end_motifs() -> Result<()> {
    // Arrange: fragment [11,21) on ACGT-repeat gives:
    // - left terminal base  seq[11] = T -> "_T"
    // - right terminal base seq[20] = A -> oriented right-end also "_T"
    // With complement collapsing enabled, both should map to canonical "_A".
    let bam = simple_paired_fragment_bam("ends_collapse", 11, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
fn sparse_output_is_the_default_when_all_motifs_is_disabled() -> Result<()> {
    // Arrange: one fragment with two observed motifs should produce a sparse 1x2 matrix by default.
    let bam = simple_paired_fragment_bam("ends_sparse_default", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
    assert!(settings.contains("\"k_within\": 1"));
    assert!(settings.contains("\"k_outside\": 0"));
    assert!(settings.contains("\"source_within\": \"reference\""));
    assert!(!settings.contains("\"output_format\""));
    assert!(!settings.contains("\"all_motifs\""));
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
    cfg.source_within = KmerSource::Reference;
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
    assert!(settings.contains("\"fragment_length_basis\": \"aligned\""));
    assert!(!settings.contains("\"output_format\""));
    assert!(!settings.contains("\"all_motifs\""));
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
    cfg.source_within = KmerSource::Reference;
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
    assert!(!out_dir.path().join("ends.bins.bed").exists());
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
    cfg.source_within = KmerSource::Reference;
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
fn all_motifs_dense_output_enumerates_outside_only_labels_when_k_within_is_zero() -> Result<()> {
    // Arrange: outside-only motifs should still have a fixed dense column universe.
    let bam = simple_paired_fragment_bam("ends_outside_only_dense", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 0, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
fn raw_endpoint_assignment_uses_the_shifted_assignment_boundaries() -> Result<()> {
    // Arrange: unpaired read-as-fragment with 2S4M2S at pos 10.
    // - aligned interval [10,14)
    // - raw assignment interval [8,16)
    // - endpoint positions 8 and 15
    // The raw terminal bases are T on the left and A on the right, which both orient to "_T".
    let bam = single_read_bam(
        "ends_raw_shifted",
        10,
        vec![('S', 2), ('M', 4), ('S', 2)],
        b"TTACGTAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 8, 9, "left_raw"), ("chr1", 15, 16, "right_raw")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Raw;
    cfg.source_within = KmerSource::Read;
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
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
fn aligned_endpoint_assignment_ignores_raw_shifted_boundary_positions() -> Result<()> {
    // Arrange: the same 2S4M2S read uses aligned endpoints [10,13] under aligned clipping, so
    // windows at the raw-shifted positions [8,9) and [15,16) should receive no counts.
    let bam = single_read_bam(
        "ends_aligned_not_raw",
        10,
        vec![('S', 2), ('M', 4), ('S', 2)],
        b"TTACGTAA",
    )?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 8, 9, "raw_left"), ("chr1", 15, 16, "raw_right")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Aligned;
    cfg.source_within = KmerSource::Read;
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
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
fn fragment_length_filters_use_the_aligned_fragment_length_even_in_raw_mode() -> Result<()> {
    // Arrange: the same 2S4M2S read has aligned length 4 but raw assignment length 8.
    // Keeping only length 4 should therefore still retain the fragment.
    let bam = single_read_bam(
        "ends_aligned_length_filter",
        10,
        vec![('S', 2), ('M', 4), ('S', 2)],
        b"TTACGTAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Raw;
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
    }

    // Act
    run(&cfg)?;
    let (motifs, matrix) = read_dense_output(out_dir.path())?;

    // Assert: both raw end motifs survive because the filter is based on the aligned 4 bp span.
    assert_eq!(motif_count(&matrix, &motifs, 0, "_T"), 2.0);
    assert_eq!(matrix.sum(), 2.0);
    Ok(())
}

#[test]
fn drop_clipping_skips_only_the_clipped_end_and_keeps_the_unclipped_end() -> Result<()> {
    // Arrange: 2S4M at pos 10 has a clipped left end and an unclipped right end.
    // With drop-clipping, only the right end should remain. The aligned terminal base is T,
    // which orients to right-end label "_A".
    let bam = single_read_bam(
        "ends_drop_clipped_end",
        10,
        vec![('S', 2), ('M', 4)],
        b"TTACGT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Drop;
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
fn drop_clipping_skips_the_fragment_when_both_ends_are_soft_clipped() -> Result<()> {
    // Arrange: 2S4M2S has soft clipping on both ends, so drop-clipping should leave no motif.
    let bam = single_read_bam(
        "ends_drop_both_clipped",
        10,
        vec![('S', 2), ('M', 4), ('S', 2)],
        b"TTACGTAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Drop;
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
    // Arrange: 2S4M has two soft-clipped bases on the left end and none on the right end.
    // With max_soft_clips=1, the left end should be skipped while the unclipped right end remains.
    let bam = single_read_bam(
        "ends_max_soft_clips",
        10,
        vec![('S', 2), ('M', 4)],
        b"TTACGT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Raw;
    cfg.clip.max_soft_clips = Some(1);
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
fn max_soft_clips_keeps_an_end_when_the_clip_count_equals_the_threshold() -> Result<()> {
    // Arrange: with max_soft_clips=2, a 2S4M read should still keep the clipped left end because
    // the documented rule is "higher number of soft-clipped bases than this".
    let bam = single_read_bam(
        "ends_max_soft_clips_equal",
        10,
        vec![('S', 2), ('M', 4)],
        b"TTACGT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Raw;
    cfg.clip.max_soft_clips = Some(2);
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
    // Arrange: 2S4M2S has two soft-clipped bases on both ends, so max_soft_clips=1 should leave
    // no surviving motifs.
    let bam = single_read_bam(
        "ends_max_soft_clips_both",
        10,
        vec![('S', 2), ('M', 4), ('S', 2)],
        b"TTACGTAA",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Raw;
    cfg.clip.max_soft_clips = Some(1);
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
    let bam = single_read_bam("ends_hard_clip", 10, vec![('H', 2), ('M', 4)], b"ACGT")?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
fn motif_labels_use_outside_within_order_when_outside_bases_are_present() -> Result<()> {
    // Arrange: left endpoint only for fragment [10,20) with one outside and one within reference base.
    // - outside base before left boundary: seq[9]  = C
    // - within base at left boundary:      seq[10] = G
    // So the final user-facing label should be "C_G".
    let bam = simple_paired_fragment_bam("ends_label_order", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 11, "left")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
fn outside_only_motifs_are_labeled_with_an_empty_within_half() -> Result<()> {
    // Arrange: left endpoint only with k_within=0 and k_outside=1.
    // The base immediately outside the left boundary of [10,20) is seq[9] = C, so the label must
    // be "C_".
    let bam = simple_paired_fragment_bam("ends_outside_only", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(&windows_bed, &[("chr1", 10, 11, "left")])?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 0, 1);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
fn raw_boundary_shifting_still_applies_when_only_outside_bases_are_counted() -> Result<()> {
    // Arrange: with k_within=0 and raw clipping, the shifted raw boundary should still control
    // endpoint assignment and outside-base extraction.
    //
    // 2S4M2S at pos 10 gives raw assignment interval [8,16), so the left endpoint is 8.
    // The base immediately outside that raw left boundary is seq[7] on the reference, which is T.
    let bam = single_read_bam(
        "ends_raw_outside_only",
        10,
        vec![('S', 2), ('M', 4), ('S', 2)],
        b"TTACGTAA",
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
    cfg.clip.clip_strategy = ClipStrategy::Raw;
    cfg.source_within = KmerSource::Reference;
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
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
            .contains("At least one of --k-within or --k-outside must be > 0")
    );
    Ok(())
}

#[test]
fn settings_json_records_clip_and_window_assignment_semantics() -> Result<()> {
    // Arrange: the settings sidecar should tell downstream users how the motifs were defined.
    let bam = single_read_bam(
        "ends_settings_semantics",
        10,
        vec![('S', 2), ('M', 4)],
        b"TTACGT",
    )?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.clip.clip_strategy = ClipStrategy::Drop;
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = false;
    cfg.set_window_assignment(AssignMotifToWindowArgs {
        assign_by: WindowMotifAssigner::Endpoint,
    });
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
    }

    // Act
    run(&cfg)?;
    let settings = read_text_file(&settings_path(out_dir.path()))?;

    // Assert
    assert!(settings.contains("\"clip_strategy\": \"drop\""));
    assert!(settings.contains("\"window_assignment\": \"endpoint\""));
    assert!(settings.contains("\"reads_are_fragments\": true"));
    assert!(settings.contains("\"source_within\": \"read\""));
    Ok(())
}

#[test]
fn reference_backed_within_bases_require_ref_2bit() -> Result<()> {
    // Arrange: reference-backed within extraction is documented to require --ref-2bit.
    let bam = simple_paired_fragment_bam("ends_missing_ref_2bit", 10, 10, 4)?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.source_within = KmerSource::Reference;
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
fn outside_bases_require_ref_2bit_even_when_within_bases_come_from_reads() -> Result<()> {
    // Arrange: outside motifs always come from the reference, so k_outside>0 still requires --ref-2bit.
    let bam = simple_paired_fragment_bam("ends_missing_ref_for_outside", 10, 10, 4)?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 1);
    cfg.source_within = KmerSource::Read;
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
fn windowed_runs_write_bins_bed_with_the_selected_windows() -> Result<()> {
    // Arrange: in BED-windowed mode the command should persist the window coordinates it used.
    let bam = simple_paired_fragment_bam("ends_bins_bed", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;
    let windows_bed = out_dir.path().join("windows.bed");
    write_bed(
        &windows_bed,
        &[("chr1", 10, 11, "left"), ("chr1", 19, 20, "right")],
    )?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
    let bins_bed = read_text_file(&out_dir.path().join("ends.bins.bed"))?;

    // Assert: BED rows should preserve the selected windows.
    let rows: Vec<&str> = bins_bed.lines().collect();
    assert_eq!(rows.len(), 2);
    assert!(rows[0].starts_with("chr1\t10\t11\t"));
    assert!(rows[1].starts_with("chr1\t19\t20\t"));
    Ok(())
}

#[test]
fn settings_json_records_fragment_length_bounds() -> Result<()> {
    // Arrange: the sidecar should tell downstream users exactly which aligned fragment lengths were kept.
    let bam = simple_paired_fragment_bam("ends_settings_lengths", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 9;
        lengths.max_fragment_length = 11;
    }

    // Act
    run(&cfg)?;
    let settings = read_text_file(&settings_path(out_dir.path()))?;

    // Assert
    assert!(settings.contains("\"min_fragment_length\": 9"));
    assert!(settings.contains("\"max_fragment_length\": 11"));
    assert!(settings.contains("\"fragment_length_basis\": \"aligned\""));
    Ok(())
}

#[test]
fn unpaired_mode_rejects_require_proper_pair() -> Result<()> {
    // Arrange: the command explicitly forbids combining reads-as-fragments with proper-pair filtering.
    let bam = single_read_bam("ends_unpaired_proper_pair", 10, vec![('M', 4)], b"ACGT")?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.require_proper_pair = true;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
    cfg.source_within = KmerSource::Reference;
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
fn all_motifs_dense_output_enumerates_within_only_labels_when_k_outside_is_zero() -> Result<()> {
    // Arrange: within-only motifs should still have a fixed dense universe under all-motifs.
    let bam = simple_paired_fragment_bam("ends_within_only_dense", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
fn read_backed_within_only_runs_without_ref_2bit() -> Result<()> {
    // Arrange: when both the within bases come from reads and k_outside=0, the command should not
    // require a reference genome.
    let bam = single_read_bam("ends_read_only_no_ref", 10, vec![('M', 4)], b"ACGT")?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_unpaired(UnpairedArgs {
        reads_are_fragments: true,
    });
    cfg.source_within = KmerSource::Read;
    cfg.all_motifs = false;
    {
        let lengths = cfg.fragment_lengths_mut();
        lengths.min_fragment_length = 4;
        lengths.max_fragment_length = 4;
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
fn by_size_windowing_writes_bins_bed() -> Result<()> {
    // Arrange: fixed-size windowing should also persist the resolved window coordinates.
    let bam = simple_paired_fragment_bam("ends_by_size_bins", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
    let bins_bed = read_text_file(&out_dir.path().join("ends.bins.bed"))?;

    // Assert
    assert!(!bins_bed.trim().is_empty());
    assert!(bins_bed.lines().all(|row| row.starts_with("chr1\t")));
    Ok(())
}

#[test]
fn output_prefix_is_applied_to_bins_bed_for_windowed_runs() -> Result<()> {
    // Arrange: prefixed runs should namespace the auxiliary bins BED too.
    let bam = simple_paired_fragment_bam("ends_prefixed_bins", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.output_prefix = "sampleA".to_string();
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
    assert!(out_dir.path().join("sampleA.bins.bed").exists());
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
    cfg.source_within = KmerSource::Reference;
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
    cfg.source_within = KmerSource::Reference;
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
fn empty_output_prefix_writes_unprefixed_bins_bed_for_windowed_runs() -> Result<()> {
    // Arrange: the empty-prefix contract should also apply to auxiliary window outputs.
    let bam = simple_paired_fragment_bam("ends_empty_prefix_bins", 10, 10, 4)?;
    let reference = simple_reference_twobit()?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.output_prefix.clear();
    cfg.set_ref_2bit(Some(reference.path.clone()));
    cfg.source_within = KmerSource::Reference;
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
    assert!(out_dir.path().join("bins.bed").exists());
    Ok(())
}

#[test]
fn read_backed_paired_within_only_runs_without_ref_2bit() -> Result<()> {
    // Arrange: paired runs also should not require a reference when both the within bases come
    // from reads and k_outside=0.
    let bam = simple_paired_fragment_bam("ends_paired_read_only_no_ref", 10, 10, 4)?;
    let out_dir = TempDir::new()?;

    let mut cfg = base_config(&bam.bam, out_dir.path(), 1, 0);
    cfg.source_within = KmerSource::Read;
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
