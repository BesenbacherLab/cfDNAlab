#![cfg(feature = "cmd_frag_to_bam")]

mod fixtures;

use anyhow::{Context, Result};
#[cfg(feature = "cmd_bam_to_frag")]
use cfdnalab::commands::bam_to_frag::{
    bam_to_frag::run_inner as run_bam_to_frag, config::BamToFragConfig,
};
use cfdnalab::commands::cli_common::ChromosomeArgs;
#[cfg(feature = "cmd_bam_to_frag")]
use cfdnalab::commands::cli_common::IOCArgs;
use cfdnalab::commands::frag_to_bam::{
    config::FragToBamConfig, frag_to_bam::run as run_frag_to_bam,
};
use cfdnalab::shared::blacklist::BlacklistStrategy;
#[cfg(feature = "cmd_bam_to_frag")]
use flate2::read::MultiGzDecoder;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};
use std::fs;
#[cfg(feature = "cmd_bam_to_frag")]
use std::io::Read as _;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

#[cfg(feature = "cmd_bam_to_frag")]
use fixtures::{FragmentSpec, ReadSpec, bam_from_specs};

#[derive(Debug, Clone, PartialEq, Eq)]
struct BamRow {
    chromosome: String,
    start: u64,
    end: u64,
    mapq: u8,
    strand: char,
    qname: String,
    flags: u16,
    cigar: String,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
    insert_size: i64,
    is_paired: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct AuxTags {
    gc_weight: Option<f32>,
    scaling_weight: Option<f32>,
    fragment_length: Option<u32>,
}

fn base_chromosomes(chromosome_names: &[&str]) -> ChromosomeArgs {
    ChromosomeArgs {
        chromosomes: Some(
            chromosome_names
                .iter()
                .map(|name| name.to_string())
                .collect(),
        ),
        chromosomes_file: None,
    }
}

fn write_frag_file(path: &Path, lines: &[&str]) -> Result<()> {
    let mut content = String::new();
    for line in lines {
        content.push_str(line);
        content.push('\n');
    }
    fs::write(path, content)?;
    Ok(())
}

fn write_chrom_sizes(path: &Path, rows: &[(&str, u32)]) -> Result<()> {
    let mut content = String::new();
    for (chromosome, length) in rows {
        content.push_str(&format!("{chromosome}\t{length}\n"));
    }
    fs::write(path, content)?;
    Ok(())
}

fn write_blacklist_bed(path: &Path, rows: &[(&str, u64, u64)]) -> Result<()> {
    let mut content = String::new();
    for (chromosome, start, end) in rows {
        content.push_str(&format!("{chromosome}\t{start}\t{end}\n"));
    }
    fs::write(path, content)?;
    Ok(())
}

fn make_config(
    frag_path: PathBuf,
    output_dir: PathBuf,
    chrom_sizes_path: PathBuf,
    chromosomes: ChromosomeArgs,
) -> FragToBamConfig {
    let mut config = FragToBamConfig::new(frag_path, output_dir, chromosomes, chrom_sizes_path);
    config.set_output_prefix("restored");
    // Keep non-length tests independent from the production defaults (30..=1000).
    // Individual tests that validate length behavior override this explicitly.
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 1;
        fragment_lengths.max_fragment_length = 1_000;
    }
    config
}

fn read_bam_rows(path: &Path) -> Result<Vec<BamRow>> {
    let mut reader = bam::Reader::from_path(path)?;
    let header = reader.header().to_owned();
    let mut rows = Vec::new();

    for record_result in reader.records() {
        let record = record_result?;
        let tid = record.tid() as u32;
        let chromosome = std::str::from_utf8(header.tid2name(tid))
            .context("invalid chromosome name in BAM header")?
            .to_string();
        let strand = if record.is_reverse() { '-' } else { '+' };
        rows.push(BamRow {
            chromosome,
            start: record.pos() as u64,
            end: record.reference_end() as u64,
            mapq: record.mapq(),
            strand,
            qname: std::str::from_utf8(record.qname())
                .context("invalid qname")?
                .to_string(),
            flags: record.flags(),
            cigar: format!("{}", record.cigar()),
            sequence: record.seq().as_bytes(),
            qualities: record.qual().to_vec(),
            insert_size: record.insert_size(),
            is_paired: record.is_paired(),
        });
    }

    Ok(rows)
}

fn read_bam_header_chromosomes(path: &Path) -> Result<Vec<String>> {
    let reader = bam::Reader::from_path(path)?;
    let header = reader.header();
    let mut chromosomes = Vec::new();
    for tid in 0..header.target_count() {
        chromosomes.push(
            std::str::from_utf8(header.tid2name(tid))
                .context("invalid chromosome name in BAM header")?
                .to_string(),
        );
    }
    Ok(chromosomes)
}

fn read_aux_tags(path: &Path) -> Result<Vec<AuxTags>> {
    let mut reader = bam::Reader::from_path(path)?;
    let mut tags = Vec::new();
    for record_result in reader.records() {
        let record = record_result?;
        let gc_weight = match record.aux(b"GC") {
            Ok(Aux::Float(value)) => Some(value),
            _ => None,
        };
        let scaling_weight = match record.aux(b"COV") {
            Ok(Aux::Float(value)) => Some(value),
            _ => None,
        };
        let fragment_length = match record.aux(b"FLEN") {
            Ok(Aux::U32(value)) => Some(value),
            _ => None,
        };
        tags.push(AuxTags {
            gc_weight,
            scaling_weight,
            fragment_length,
        });
    }
    Ok(tags)
}

fn assert_optional_f32_eq(actual: Option<f32>, expected: Option<f32>, label: &str) {
    match (actual, expected) {
        (None, None) => {}
        (Some(actual_value), Some(expected_value)) => {
            let diff = (actual_value - expected_value).abs();
            assert!(
                diff <= 1e-6,
                "{} mismatch: expected {}, got {}",
                label,
                expected_value,
                actual_value
            );
        }
        _ => panic!(
            "{} mismatch: expected {:?}, got {:?}",
            label, expected, actual
        ),
    }
}

fn assert_unpaired_full_match_record(
    row: &BamRow,
    expected_chromosome: &str,
    expected_start: u64,
    expected_end: u64,
    expected_mapq: u8,
    expected_strand: char,
    expected_qname: &str,
) {
    let expected_length = expected_end - expected_start;
    assert_eq!(row.chromosome, expected_chromosome);
    assert_eq!(row.start, expected_start);
    assert_eq!(row.end, expected_end);
    assert_eq!(row.mapq, expected_mapq);
    assert_eq!(row.strand, expected_strand);
    assert_eq!(row.qname, expected_qname);
    assert_eq!(row.insert_size, 0);
    assert!(
        !row.is_paired,
        "frag-to-bam output must be unpaired, but record had the paired flag set"
    );

    // We generate one contiguous match over the full fragment span.
    assert_eq!(row.cigar, format!("{expected_length}M"));

    // The converter stores unknown sequence as 'N' and fixed quality 40.
    assert_eq!(row.sequence.len(), expected_length as usize);
    assert!(row.sequence.iter().all(|base| *base == b'N'));
    assert_eq!(row.qualities.len(), expected_length as usize);
    assert!(row.qualities.iter().all(|quality| *quality == 40));

    let expected_flags = if expected_strand == '-' { 0x10 } else { 0 };
    assert_eq!(row.flags, expected_flags);
}

fn run_blacklist_strategy_case(strategy: BlacklistStrategy) -> Result<Vec<u64>> {
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");
    let blacklist_path = input_dir.path().join("blacklist.bed");

    // Four 10bp fragments and one blacklist interval [10,20):
    // - [0,10)   overlap=0
    // - [5,15)   overlap=5/10, midpoint=10
    // - [10,20)  overlap=10/10, midpoint=15
    // - [15,25)  overlap=5/10, midpoint=20 (outside half-open [10,20))
    write_frag_file(
        &frag_path,
        &[
            "chr1\t0\t10\t60\t+",
            "chr1\t5\t15\t60\t+",
            "chr1\t10\t20\t60\t+",
            "chr1\t15\t25\t60\t+",
        ],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;
    write_blacklist_bed(&blacklist_path, &[("chr1", 10, 20)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_blacklist(Some(vec![blacklist_path]));
    config.set_blacklist_strategy(strategy);
    config.set_min_mapq(0);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 1;
        fragment_lengths.max_fragment_length = 100;
    }

    run_frag_to_bam(&config)?;

    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    Ok(rows.iter().map(|row| row.start).collect())
}

#[test]
fn given_valid_frag_when_run_then_writes_expected_unpaired_bam_records() -> Result<()> {
    // Arrange:
    // Two fragments in input order chr1 then chr2.
    // Chrom sizes order is intentionally chr2 then chr1.
    // The writer iterates chrom-sizes order in the second pass, so output row order should be:
    //   1) chr2 fragment [5,9), mapq=30, strand='-'
    //   2) chr1 fragment [10,20), mapq=60, strand='+'
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\t+", "chr2\t5\t9\t30\t-"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr2", 100), ("chr1", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1", "chr2"]),
    );

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr2", 5, 9, 30, '-', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 10, 20, 60, '+', "fragment_2");

    Ok(())
}

#[test]
fn given_filters_and_extra_columns_when_run_then_only_expected_fragments_remain() -> Result<()> {
    // Arrange:
    // - Keep chromosomes: chr1 only.
    // - Keep mapq >= 20.
    // - Keep fragment length in [5, 25].
    // Expected keeps:
    // - chr1 [0,20) mapq 60
    // - chr1 [40,50) mapq 60 with extra trailing columns (ignored by parser)
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chr1\t0\t20\t60\t+",
            "chr1\t10\t25\t10\t+",
            "chr1\t30\t33\t60\t-",
            "chr1\t40\t50\t60\t-\textra_col\tignored_col",
            "chr2\t5\t15\t50\t+",
        ],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100), ("chr2", 100)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_min_mapq(20);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 5;
        fragment_lengths.max_fragment_length = 25;
    }

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr1", 0, 20, 60, '+', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 40, 50, 60, '-', "fragment_2");

    Ok(())
}

#[test]
fn given_length_bounds_when_run_then_only_fragments_within_inclusive_range_are_kept() -> Result<()>
{
    // Arrange:
    // Fragment lengths are:
    // - [0,9)   -> 9   (below min, excluded)
    // - [10,20) -> 10  (equal to min, included)
    // - [30,50) -> 20  (equal to max, included)
    // - [60,81) -> 21  (above max, excluded)
    //
    // With min=10 and max=20 we expect only starts [10,30] to remain.
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chr1\t0\t9\t60\t+",
            "chr1\t10\t20\t60\t+",
            "chr1\t30\t50\t60\t-",
            "chr1\t60\t81\t60\t+",
        ],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 200)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_min_mapq(0);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 20;
    }

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 20, 60, '+', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 30, 50, 60, '-', "fragment_2");

    Ok(())
}

#[test]
fn given_inline_header_with_unsupported_extra_column_when_run_then_returns_clear_error()
-> Result<()> {
    // Arrange:
    // The header contains one extra column `gc`.
    // Only these extra names are supported: gc_weight, scaling_weight, flen.
    // Expected behavior is a validation error before conversion starts.
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tgc",
            "chr1\t10\t20\t60\t+\t0.3",
        ],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    // Act
    let error = run_frag_to_bam(&config).expect_err("expected unsupported extra column error");

    // Assert
    let error_text = format!("{error}");
    assert!(
        error_text.contains("Unsupported frag header column name(s): gc"),
        "unexpected error: {error}"
    );
    assert!(
        error_text.contains("gc_weight, scaling_weight, or flen"),
        "unexpected error: {error}"
    );

    Ok(())
}

#[test]
fn given_unsupported_extra_columns_and_ignore_extras_when_run_then_conversion_succeeds_without_aux_tags()
-> Result<()> {
    // Arrange:
    // Same unsupported header as the previous test, but `--ignore-extras` is enabled.
    // Hand-derived expectation:
    // - The fragment [10,20) passes all filters and becomes one BAM record.
    // - Extra column `gc` is ignored, so GC/COV/FLEN tags are absent.
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tgc",
            "chr1\t10\t20\t60\t+\t0.3",
        ],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_ignore_extras(true);

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 20, 60, '+', "fragment_1");

    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, None);

    Ok(())
}

#[test]
fn given_inline_header_with_unknown_extra_and_allow_unknown_extras_when_run_then_known_extras_are_still_transferred()
-> Result<()> {
    // Arrange:
    // Inline header has unknown `gc` and known `flen`.
    // Hand-derived expectation:
    // - Unknown `gc` is ignored
    // - Known `flen` still maps to FLEN=30
    // - GC and COV tags are absent
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tgc\tflen",
            "chr1\t10\t40\t60\t+\t0.2\t30",
        ],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 200)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_allow_unknown_extras(true);

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 40, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, Some(30));

    Ok(())
}

#[test]
fn given_supported_extra_column_names_when_run_then_gc_cov_and_flen_are_transferred_to_aux_tags()
-> Result<()> {
    // Arrange:
    // Header uses the three supported extra names exactly.
    // Hand-derived tag expectations:
    // Row 1 -> GC=0.25, COV=1.5, FLEN=10
    // Row 2 -> GC=None ("na"), COV=None ("."), FLEN=11
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tgc_weight\tscaling_weight\tflen",
            "chr1\t0\t10\t60\t+\t0.25\t1.5\t10",
            "chr1\t20\t31\t40\t-\tna\t.\t11",
        ],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr1", 0, 10, 60, '+', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 20, 31, 40, '-', "fragment_2");

    assert_eq!(aux_tags.len(), 2);
    assert_optional_f32_eq(aux_tags[0].gc_weight, Some(0.25), "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, Some(1.5), "COV");
    assert_eq!(aux_tags[0].fragment_length, Some(10));
    assert_optional_f32_eq(aux_tags[1].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[1].scaling_weight, None, "COV");
    assert_eq!(aux_tags[1].fragment_length, Some(11));

    Ok(())
}

#[test]
fn given_no_header_and_extra_columns_when_run_then_extra_columns_are_ignored_and_no_aux_tags_are_written()
-> Result<()> {
    // Arrange:
    // No inline header is present, no explicit header is configured, and no companion
    // header file exists for this file name. The parser therefore uses fixed 5-column
    // layout and ignores trailing columns.
    //
    // Hand-derived expectation:
    // - One fragment [10,20) is converted to one BAM record
    // - Trailing values `0.25 1.5 10` are ignored
    // - GC, COV, and FLEN tags are absent
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\t+\t0.25\t1.5\t10"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 20, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, None);

    Ok(())
}

#[test]
fn given_inline_header_with_only_flen_when_run_then_only_flen_aux_tag_is_written() -> Result<()> {
    // Arrange:
    // Inline header defines only one supported extra column (`flen`).
    // Hand-derived expectation:
    // - Row 1 has flen=80 so FLEN=80 is written
    // - Row 2 has flen="." so FLEN is absent
    // - GC and COV remain absent for both rows
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tflen",
            "chr1\t0\t80\t60\t+\t80",
            "chr1\t100\t180\t55\t-\t.",
        ],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 500)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr1", 0, 80, 60, '+', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 100, 180, 55, '-', "fragment_2");

    assert_eq!(aux_tags.len(), 2);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, Some(80));
    assert_optional_f32_eq(aux_tags[1].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[1].scaling_weight, None, "COV");
    assert_eq!(aux_tags[1].fragment_length, None);

    Ok(())
}

#[test]
fn given_explicit_header_with_only_flen_when_run_then_only_flen_aux_tag_is_written() -> Result<()> {
    // Arrange:
    // Frag file has no inline header. We pass an explicit header file that maps column 6 to `flen`.
    // Hand-derived expectation:
    // - One fragment [10,60) with FLEN=50
    // - GC and COV remain absent
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let header_path = input_dir.path().join("explicit_header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t60\t42\t+\t50"])?;
    write_frag_file(
        &header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tflen"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 500)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_frag_header(Some(header_path));

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 60, 42, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, Some(50));

    Ok(())
}

#[test]
fn given_explicit_header_with_unsupported_extra_column_when_run_then_returns_clear_error()
-> Result<()> {
    // Arrange:
    // Explicit header is chosen first and contains unsupported extra name `gc`.
    // Hand-derived expectation is a validation error with the unsupported name.
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let explicit_header_path = input_dir.path().join("explicit_header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\t+\t0.3"])?;
    write_frag_file(
        &explicit_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tgc"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_frag_header(Some(explicit_header_path));

    // Act
    let error = run_frag_to_bam(&config).expect_err("expected unsupported extra column error");

    // Assert
    let error_text = format!("{error}");
    assert!(
        error_text.contains("Unsupported frag header column name(s): gc"),
        "unexpected error: {error}"
    );
    assert!(
        error_text.contains("gc_weight, scaling_weight, or flen"),
        "unexpected error: {error}"
    );

    Ok(())
}

#[test]
fn given_explicit_header_with_unsupported_extra_column_and_ignore_extras_when_run_then_conversion_succeeds_without_aux_tags()
-> Result<()> {
    // Arrange:
    // Explicit header includes unsupported extra name `gc`.
    // With --ignore-extras, the header is accepted and only the core five columns are used.
    // Hand-derived expectation:
    // - One fragment [10,20) is written
    // - GC, COV, and FLEN tags are absent
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let explicit_header_path = input_dir.path().join("explicit_header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\t+\t0.3"])?;
    write_frag_file(
        &explicit_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tgc"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_frag_header(Some(explicit_header_path));
    config.set_ignore_extras(true);

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 20, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, None);

    Ok(())
}

#[test]
fn given_explicit_header_with_unknown_extra_and_allow_unknown_extras_when_run_then_known_extras_are_still_transferred()
-> Result<()> {
    // Arrange:
    // Explicit header has unknown `gc` and known `flen`.
    // Hand-derived expectation:
    // - Unknown `gc` is ignored
    // - Known `flen` maps to FLEN=30
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let explicit_header_path = input_dir.path().join("explicit_header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t5\t35\t60\t-\t0.2\t30"])?;
    write_frag_file(
        &explicit_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tgc\tflen"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 200)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_frag_header(Some(explicit_header_path));
    config.set_allow_unknown_extras(true);

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 5, 35, 60, '-', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, Some(30));

    Ok(())
}

#[test]
fn given_companion_header_with_only_flen_when_run_then_only_flen_aux_tag_is_written() -> Result<()>
{
    // Arrange:
    // No explicit header is configured and frag has no inline header.
    // The companion header path is auto-detected from `<prefix>.frag.tsv`.
    // Hand-derived expectation:
    // - One record [20,70) with FLEN=50 from the 6th column
    // - GC and COV remain absent
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("sample.frag.tsv");
    let companion_header_path = input_dir.path().join("sample.frag.header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t20\t70\t39\t-\t50"])?;
    write_frag_file(
        &companion_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tflen"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 500)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 20, 70, 39, '-', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, Some(50));

    Ok(())
}

#[test]
fn given_companion_header_with_unsupported_extra_column_and_ignore_extras_when_run_then_conversion_succeeds_without_aux_tags()
-> Result<()> {
    // Arrange:
    // Companion header includes unsupported extra name `gc`.
    // With --ignore-extras, conversion should continue using only core columns.
    // Hand-derived expectation:
    // - One fragment [10,20) is written
    // - GC, COV, and FLEN tags are absent
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("sample.frag.tsv");
    let companion_header_path = input_dir.path().join("sample.frag.header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\t+\t0.3"])?;
    write_frag_file(
        &companion_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tgc"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_ignore_extras(true);

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 20, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, None);

    Ok(())
}

#[test]
fn given_companion_header_with_unknown_extra_and_allow_unknown_extras_when_run_then_known_extras_are_still_transferred()
-> Result<()> {
    // Arrange:
    // Companion header has unknown `gc` and known `flen`.
    // Hand-derived expectation:
    // - Unknown `gc` is ignored
    // - Known `flen` maps to FLEN=30
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("sample.frag.tsv");
    let companion_header_path = input_dir.path().join("sample.frag.header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t5\t35\t60\t-\t0.2\t30"])?;
    write_frag_file(
        &companion_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tgc\tflen"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 200)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_allow_unknown_extras(true);

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 5, 35, 60, '-', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");
    assert_eq!(aux_tags[0].fragment_length, Some(30));

    Ok(())
}

#[test]
fn given_companion_header_with_unsupported_extra_column_when_run_then_returns_clear_error()
-> Result<()> {
    // Arrange:
    // No explicit header is set. Companion header is therefore selected.
    // Companion header contains unsupported extra name `gc`, so validation must fail.
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("sample.frag.tsv");
    let companion_header_path = input_dir.path().join("sample.frag.header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\t+\t0.3"])?;
    write_frag_file(
        &companion_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tgc"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    // Act
    let error = run_frag_to_bam(&config).expect_err("expected unsupported extra column error");

    // Assert
    let error_text = format!("{error}");
    assert!(
        error_text.contains("Unsupported frag header column name(s): gc"),
        "unexpected error: {error}"
    );
    assert!(
        error_text.contains("gc_weight, scaling_weight, or flen"),
        "unexpected error: {error}"
    );

    Ok(())
}

#[test]
fn given_companion_and_inline_headers_when_run_then_returns_header_conflict_error() -> Result<()> {
    // Arrange:
    // Both a companion header file and an inline header are present.
    // Hand-derived expectation:
    // - Conversion fails with a clear conflict error
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("sample.frag.tsv");
    let companion_header_path = input_dir.path().join("sample.frag.header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tflen",
            "chr1\t0\t40\t60\t+\t40",
        ],
    )?;
    write_frag_file(
        &companion_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 500)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    // Act
    let error = run_frag_to_bam(&config).expect_err("expected header conflict error");

    // Assert
    let error_text = format!("{error}");
    assert!(
        error_text.contains("Conflicting headers detected"),
        "unexpected error: {error}"
    );
    assert!(
        error_text.contains("companion header file"),
        "unexpected error: {error}"
    );

    Ok(())
}

#[test]
fn given_explicit_and_companion_headers_when_run_then_explicit_header_takes_precedence()
-> Result<()> {
    // Arrange:
    // Companion header is core-only. Explicit header maps the 6th column to flen.
    // Priority order uses explicit `--frag-header` first.
    //
    // Hand-derived expectation:
    // - Fragment converts and FLEN=30 is written
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("sample.frag.tsv");
    let companion_header_path = input_dir.path().join("sample.frag.header.tsv");
    let explicit_header_path = input_dir.path().join("explicit_header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t5\t35\t60\t-\t30"])?;
    write_frag_file(
        &companion_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand"],
    )?;
    write_frag_file(
        &explicit_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tflen"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 500)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_frag_header(Some(explicit_header_path));

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 5, 35, 60, '-', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_eq!(aux_tags[0].fragment_length, Some(30));
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].scaling_weight, None, "COV");

    Ok(())
}

#[test]
fn given_explicit_and_inline_headers_when_run_then_returns_header_conflict_error() -> Result<()> {
    // Arrange:
    // Both --frag-header and inline header are present.
    // Hand-derived expectation:
    // - Conversion fails with a clear conflict error
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("sample.frag.tsv");
    let explicit_header_path = input_dir.path().join("explicit_header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tflen",
            "chr1\t5\t35\t60\t-\t30",
        ],
    )?;
    write_frag_file(
        &explicit_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tflen"],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 500)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_frag_header(Some(explicit_header_path));

    // Act
    let error = run_frag_to_bam(&config).expect_err("expected header conflict error");

    // Assert
    let error_text = format!("{error}");
    assert!(
        error_text.contains("Conflicting headers detected"),
        "unexpected error: {error}"
    );
    assert!(
        error_text.contains("--frag-header"),
        "unexpected error: {error}"
    );

    Ok(())
}

#[test]
fn given_blacklist_any_when_run_then_any_overlap_is_excluded() -> Result<()> {
    // Any overlap excludes [5,15), [10,20), and [15,25), so only [0,10) remains.
    let kept_starts = run_blacklist_strategy_case(BlacklistStrategy::Any)?;
    assert_eq!(kept_starts, vec![0]);
    Ok(())
}

#[test]
fn given_blacklist_all_when_run_then_only_fully_overlapped_fragments_are_excluded() -> Result<()> {
    // Full overlap excludes only [10,20).
    let kept_starts = run_blacklist_strategy_case(BlacklistStrategy::All)?;
    assert_eq!(kept_starts, vec![0, 5, 15]);
    Ok(())
}

#[test]
fn given_blacklist_midpoint_when_run_then_midpoint_overlap_controls_exclusion() -> Result<()> {
    // Midpoints:
    // - [5,15) midpoint=10 -> excluded
    // - [10,20) midpoint=15 -> excluded
    // - [15,25) midpoint=20 -> kept (end-exclusive blacklist interval)
    let kept_starts = run_blacklist_strategy_case(BlacklistStrategy::Midpoint)?;
    assert_eq!(kept_starts, vec![0, 15]);
    Ok(())
}

#[test]
fn given_blacklist_proportion_when_run_then_threshold_controls_exclusion() -> Result<()> {
    // Overlap fractions are [0.0, 0.5, 1.0, 0.5] for the four fragments.
    // With threshold 0.6, only [10,20) is excluded.
    let kept_starts = run_blacklist_strategy_case(BlacklistStrategy::Proportion(0.6))?;
    assert_eq!(kept_starts, vec![0, 5, 15]);
    Ok(())
}

#[test]
fn given_unsorted_fragments_within_chromosome_when_run_then_returns_order_error() -> Result<()> {
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t20\t30\t60\t+", "chr1\t10\t20\t60\t+"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    let error = run_frag_to_bam(&config).expect_err("expected out-of-order error");
    assert!(
        format!("{error}").contains("out of order on chr1"),
        "unexpected error: {error}"
    );
    Ok(())
}

#[test]
fn given_chromosome_reappears_when_run_then_returns_order_error() -> Result<()> {
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chr1\t10\t20\t60\t+",
            "chr2\t5\t15\t60\t+",
            "chr1\t30\t40\t60\t+",
        ],
    )?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100), ("chr2", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1", "chr2"]),
    );

    let error = run_frag_to_bam(&config).expect_err("expected chromosome reappearance error");
    assert!(
        format!("{error}").contains("appears after moving past it"),
        "unexpected error: {error}"
    );
    Ok(())
}

#[test]
fn given_fragment_exceeds_chromosome_when_run_then_returns_bounds_error() -> Result<()> {
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t60\t60\t+"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 50)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    let error = run_frag_to_bam(&config).expect_err("expected bounds error");
    assert!(
        format!("{error}").contains("exceeds chromosome bounds"),
        "unexpected error: {error}"
    );
    Ok(())
}

#[test]
fn given_invalid_strand_when_run_then_returns_parse_error() -> Result<()> {
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\tx"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    let error = run_frag_to_bam(&config).expect_err("expected invalid strand error");
    assert!(
        format!("{error}").contains("Strand must be '+' or '-'"),
        "unexpected error: {error}"
    );
    Ok(())
}

#[test]
fn given_all_fragments_fail_mapq_when_run_then_writes_empty_bam() -> Result<()> {
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t10\t+"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_min_mapq(20);

    // `frag-to-bam` only raises "No fragments passed filters" when no chromosome
    // is observed at all. If chromosomes are observed but all fragments are later
    // filtered (e.g. by mapq), it returns success and writes an empty BAM.
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let rows = read_bam_rows(&output_bam_path)?;
    assert!(rows.is_empty(), "expected no records in output BAM");

    Ok(())
}

#[test]
fn given_all_fragments_fail_chromosome_selection_when_run_then_returns_clear_error() -> Result<()> {
    // The command raises "No fragments passed filters; no BAM to write" when
    // no chromosome survives selection and therefore no chromosome is observed.
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr2\t10\t20\t60\t+"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100), ("chr2", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );

    let error = run_frag_to_bam(&config).expect_err("expected no-fragments error");
    assert!(
        format!("{error}").contains("No fragments passed filters; no BAM to write"),
        "unexpected error: {error}"
    );

    Ok(())
}

#[test]
fn given_chromosomes_all_when_run_then_header_and_output_follow_chrom_sizes_order() -> Result<()> {
    // Arrange:
    // `--chromosomes all` in this command resolves to chrom-sizes order.
    // Chrom sizes file order: chr2, chr1, chr3.
    // Input has fragments for chr1 and chr2 only.
    // Expected:
    // - Header targets: [chr2, chr1, chr3]
    // - Output records ordered by that same chromosome order: chr2 row first, then chr1 row.
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\t+", "chr2\t5\t9\t30\t-"])?;
    write_chrom_sizes(
        &chrom_sizes_path,
        &[("chr2", 100), ("chr1", 100), ("chr3", 100)],
    )?;

    let all_chromosomes = ChromosomeArgs {
        chromosomes: Some(vec!["all".to_string()]),
        chromosomes_file: None,
    };
    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        all_chromosomes,
    );

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_dir.path().join("restored.bam");
    let header_chromosomes = read_bam_header_chromosomes(&output_bam_path)?;
    let rows = read_bam_rows(&output_bam_path)?;

    // Assert
    assert_eq!(
        header_chromosomes,
        vec!["chr2".to_string(), "chr1".to_string(), "chr3".to_string()]
    );
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].chromosome, "chr2");
    assert_eq!(rows[1].chromosome, "chr1");

    Ok(())
}

#[cfg(feature = "cmd_bam_to_frag")]
fn read_gzip_text(path: &Path) -> Result<String> {
    let file = fs::File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut decoder = MultiGzDecoder::new(file);
    let mut text = String::new();
    decoder.read_to_string(&mut text)?;
    Ok(text)
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn given_bam_to_frag_then_frag_to_bam_when_roundtrip_then_restores_all_available_fields()
-> Result<()> {
    // Arrange:
    // We create two inward fragments whose restorable fields are hand-derived:
    //
    // Fragment A (chr1):
    // - forward read is read1 at 100, reverse read is read2 at 160 with 20M each
    // - fragment span: [100, 180), length 80
    // - min mapq: min(55,45)=45
    // - read1 strand: '+'
    //
    // Fragment B (chr2):
    // - forward read is read2 at 200, reverse read is read1 at 260 with 20M each
    // - fragment span: [200, 280), length 80
    // - min mapq: min(35,30)=30
    // - read1 strand: '-' (read1 is reverse)
    //
    // So `bam-to-frag` should emit:
    // - chr1 100 180 45 +
    // - chr2 200 280 30 -
    //
    // Then `frag-to-bam` should restore all available fields:
    // chromosome, start, end (via CIGAR), mapq, strand, and one record per fragment.
    // It cannot restore pair structure, mate fields, original sequence, or original qualities.
    let first_mate_flag: u16 = 0x40;
    let second_mate_flag: u16 = 0x80;
    let proper_pair_flag: u16 = 0x2;
    let mate_reverse_flag: u16 = 0x20;

    let fragment_a = FragmentSpec {
        forward: ReadSpec {
            tid: 0,
            pos: 100,
            cigar: vec![('M', 20)],
            seq: vec![b'A'; 20],
            qual: 30,
            is_reverse: false,
            mapq: 55,
            flags: first_mate_flag | proper_pair_flag | mate_reverse_flag,
            mate_tid: Some(0),
            mate_pos: Some(160),
            insert_size: 80,
        },
        reverse: ReadSpec {
            tid: 0,
            pos: 160,
            cigar: vec![('M', 20)],
            seq: vec![b'T'; 20],
            qual: 30,
            is_reverse: true,
            mapq: 45,
            flags: second_mate_flag | proper_pair_flag,
            mate_tid: Some(0),
            mate_pos: Some(100),
            insert_size: -80,
        },
    };

    let fragment_b = FragmentSpec {
        forward: ReadSpec {
            tid: 1,
            pos: 200,
            cigar: vec![('M', 20)],
            seq: vec![b'C'; 20],
            qual: 30,
            is_reverse: false,
            mapq: 30,
            flags: second_mate_flag | proper_pair_flag | mate_reverse_flag,
            mate_tid: Some(1),
            mate_pos: Some(260),
            insert_size: 80,
        },
        reverse: ReadSpec {
            tid: 1,
            pos: 260,
            cigar: vec![('M', 20)],
            seq: vec![b'G'; 20],
            qual: 30,
            is_reverse: true,
            mapq: 35,
            flags: first_mate_flag | proper_pair_flag,
            mate_tid: Some(1),
            mate_pos: Some(200),
            insert_size: -80,
        },
    };

    let source_bam = bam_from_specs(
        vec![("chr1".to_string(), 500), ("chr2".to_string(), 500)],
        vec![fragment_a, fragment_b],
        Vec::new(),
        "frag_to_bam_roundtrip",
    )?;

    let bam_to_frag_output = TempDir::new()?;
    let frag_to_bam_output = TempDir::new()?;
    let chrom_sizes_path = frag_to_bam_output.path().join("chrom.sizes");
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 500), ("chr2", 500)])?;

    let mut bam_to_frag_config = BamToFragConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: bam_to_frag_output.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1", "chr2"]),
    );
    bam_to_frag_config.set_output_prefix("roundtrip");
    bam_to_frag_config.set_min_mapq(0);
    {
        let fragment_lengths = bam_to_frag_config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 500;
    }

    // Act 1: BAM -> FRAG
    run_bam_to_frag(&bam_to_frag_config)?;
    let frag_path = bam_to_frag_output.path().join("roundtrip.frag.tsv.gz");
    let frag_text = read_gzip_text(&frag_path)?;
    let frag_lines: Vec<&str> = frag_text.lines().collect();

    // Assert 1: hand-derived FRAG rows
    assert_eq!(
        frag_lines,
        vec!["chr1\t100\t180\t45\t+", "chr2\t200\t280\t30\t-"]
    );

    let mut frag_to_bam_config = make_config(
        frag_path,
        frag_to_bam_output.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1", "chr2"]),
    );
    frag_to_bam_config.set_min_mapq(0);
    {
        let fragment_lengths = frag_to_bam_config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 500;
    }

    // Act 2: FRAG -> BAM
    run_frag_to_bam(&frag_to_bam_config)?;
    let restored_bam_path = frag_to_bam_output.path().join("restored.bam");
    let restored_rows = read_bam_rows(&restored_bam_path)?;

    // Assert 2: all available fields are restored exactly
    assert_eq!(restored_rows.len(), 2);
    assert_unpaired_full_match_record(&restored_rows[0], "chr1", 100, 180, 45, '+', "fragment_1");
    assert_unpaired_full_match_record(&restored_rows[1], "chr2", 200, 280, 30, '-', "fragment_2");

    Ok(())
}
