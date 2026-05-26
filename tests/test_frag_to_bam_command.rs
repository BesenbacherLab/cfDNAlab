#![cfg(feature = "cmd_frag_to_bam")]

// KEEP-IN-TESTS: all active tests in this file cover frag-to-bam command output, errors, or end-to-end round trips.

mod fixtures;

use anyhow::{Context, Result};
use cfdnalab::RunOptions;
use cfdnalab::blacklist::BlacklistStrategy;
use cfdnalab::constants::GC_CORRECTION_SCHEMA_VERSION;
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
use cfdnalab::gc_bias::GCCorrectionPackage;
#[cfg(feature = "cmd_bam_to_frag")]
use cfdnalab::reference::twobit_contig_footprint;
#[cfg(all(feature = "cmd_bam_to_bam", feature = "cmd_bam_to_frag"))]
use cfdnalab::run_like_cli::bam_to_bam::{
    BamToBamConfig, run_bam_to_bam as run_bam_to_bam_command,
};
#[cfg(feature = "cmd_bam_to_frag")]
use cfdnalab::run_like_cli::bam_to_frag::{
    BamToFragConfig, run_bam_to_frag as run_bam_to_frag_command,
};
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_lengths"))]
use cfdnalab::run_like_cli::common::AssignToWindowArgs;
use cfdnalab::run_like_cli::common::ChromosomeArgs;
#[cfg(feature = "cmd_bam_to_frag")]
use cfdnalab::run_like_cli::common::IOCArgs;
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
use cfdnalab::run_like_cli::common::{ApplyGCArgFileOnly, ApplyGCArgs};
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_fcoverage"))]
use cfdnalab::run_like_cli::fcoverage::{FCoverageConfig, run_fcoverage as run_fcoverage_command};
use cfdnalab::run_like_cli::frag_to_bam::{
    FragToBamConfig, run_frag_to_bam as run_frag_to_bam_command,
};
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_lengths"))]
use cfdnalab::run_like_cli::lengths::{LengthsConfig, run_lengths as run_lengths_command};
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
use cfdnalab::run_like_cli::midpoints::{
    MidpointSmoothing, MidpointsConfig, run_midpoints as run_midpoints_command,
};
#[cfg(feature = "cmd_bam_to_frag")]
use flate2::read::MultiGzDecoder;
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
use ndarray::Array3;
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
use ndarray::array;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};
use std::fs;
#[cfg(feature = "cmd_bam_to_frag")]
use std::io::Read as _;
use std::path::{Path, PathBuf};
use tempfile::TempDir;

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_lengths"))]
use fixtures::read_length_counts_tsv;
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
use fixtures::read_midpoint_zarr_counts;
#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_fcoverage"))]
use fixtures::read_zst_to_string;
#[cfg(feature = "cmd_bam_to_frag")]
use fixtures::{
    FragmentSpec, ReadSpec, bam_from_specs, bam_from_specs_strict_identity,
    build_real_non_neutral_gc_package, paired_fragment, simple_reference_twobit,
    twobit_from_sequences,
};

fn run_frag_to_bam(config: &FragToBamConfig) -> Result<()> {
    run_frag_to_bam_command(config, RunOptions::new_quiet()).map(|_| ())
}

#[cfg(feature = "cmd_bam_to_frag")]
fn run_bam_to_frag(config: &BamToFragConfig) -> Result<()> {
    run_bam_to_frag_command(config, RunOptions::new_quiet()).map(|_| ())
}

#[cfg(all(feature = "cmd_bam_to_bam", feature = "cmd_bam_to_frag"))]
fn run_bam_to_bam(config: &BamToBamConfig) -> Result<()> {
    run_bam_to_bam_command(config, RunOptions::new_quiet()).map(|_| ())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_fcoverage"))]
fn run_fcoverage(config: &FCoverageConfig) -> Result<()> {
    run_fcoverage_command(config, RunOptions::new_quiet()).map(|_| ())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_lengths"))]
fn run_lengths(config: &LengthsConfig) -> Result<()> {
    run_lengths_command(config, RunOptions::new_quiet()).map(|_| ())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
fn run_midpoints(config: &MidpointsConfig) -> Result<()> {
    run_midpoints_command(config, RunOptions::new_quiet()).map(|_| ())
}

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
    coverage_scaling_weight: Option<f32>,
    count_scaling_weight: Option<f32>,
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

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
fn read_group_index_map(path: &Path) -> Result<std::collections::HashMap<String, usize>> {
    let text = fs::read_to_string(path)?;
    let mut out = std::collections::HashMap::new();
    for line in text.lines().skip(1) {
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let idx = fields
            .next()
            .context("missing group index field")?
            .parse::<usize>()?;
        let name = fields
            .next()
            .context("missing group name field")?
            .to_string();
        out.insert(name, idx);
    }
    Ok(out)
}

fn write_scaling_tsv(path: &Path, rows: &[(&str, u64, u64, f32)]) -> Result<()> {
    let mut content = String::from("chromosome\tstart\tend\tscaling_factor\n");
    for (chromosome, start, end, factor) in rows {
        content.push_str(&format!("{chromosome}\t{start}\t{end}\t{factor}\n"));
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
    // Still respect the shared minimum included fragment length of 10.
    // Individual tests that validate length behavior override this explicitly.
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 1_000;
    }
    config
}

fn output_bam_path(output_dir: &Path, prefix: &str) -> PathBuf {
    output_dir.join(format!("{prefix}.fragments.bam"))
}

fn dot_join(parts: &[&str]) -> String {
    parts.join(".")
}

fn build_bai_for_test_bam(bam_path: &Path) -> Result<PathBuf> {
    let bai_path = bam_path.with_extension("bam.bai");
    bam::index::build(bam_path, None, bam::index::Type::Bai, 1)
        .with_context(|| format!("index bam {}", bam_path.display()))?;
    let target = bam_path.with_extension("bai");
    if bai_path.exists() {
        fs::rename(&bai_path, &target)?;
    }
    Ok(target)
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
        let coverage_scaling_weight = match record.aux(b"cw") {
            Ok(Aux::Float(value)) => Some(value),
            _ => None,
        };
        let count_scaling_weight = match record.aux(b"nw") {
            Ok(Aux::Float(value)) => Some(value),
            _ => None,
        };
        let fragment_length = match record.aux(b"fl") {
            Ok(Aux::U32(value)) => Some(value),
            _ => None,
        };
        tags.push(AuxTags {
            gc_weight,
            coverage_scaling_weight,
            count_scaling_weight,
            fragment_length,
        });
    }
    Ok(tags)
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

    // Five 10bp fragments and one blacklist interval [10,20):
    // - [0,10)   overlap=0
    // - [5,15)   overlap=5/10, central bases 9 and 10
    // - [10,20)  overlap=10/10, central bases 14 and 15
    // - [15,25)  overlap=5/10, central bases 19 and 20
    // - [16,26)  overlap=4/10, central bases 20 and 21
    write_frag_file(
        &frag_path,
        &[
            "chr1\t0\t10\t60\t+",
            "chr1\t5\t15\t60\t+",
            "chr1\t10\t20\t60\t+",
            "chr1\t15\t25\t60\t+",
            "chr1\t16\t26\t60\t+",
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
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    run_frag_to_bam(&config)?;

    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    Ok(rows.iter().map(|row| row.start).collect())
}

#[test]
fn given_valid_frag_when_run_then_writes_expected_unpaired_bam_records() -> Result<()> {
    // Arrange:
    // Two fragments in input order chr1 then chr2.
    // Chrom sizes order is intentionally chr2 then chr1.
    // The writer iterates chrom-sizes order in the second pass, so output row order should be:
    //   1) chr2 fragment [5,15), mapq=30, strand='-'
    //   2) chr1 fragment [10,20), mapq=60, strand='+'
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\t+", "chr2\t5\t15\t30\t-"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr2", 100), ("chr1", 100)])?;

    let config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1", "chr2"]),
    );

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr2", 5, 15, 30, '-', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 10, 20, 60, '+', "fragment_2");

    Ok(())
}

#[test]
fn given_filters_and_extra_columns_when_run_then_only_expected_fragments_remain() -> Result<()> {
    // Arrange:
    // - Keep chromosomes: chr1 only.
    // - Keep mapq >= 20.
    // - Keep fragment length in [10, 25].
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
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 25;
    }

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
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
    // Only these extra names are supported: gc_weight, coverage_scaling_weight,
    // count_scaling_weight, and flen.
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
        error_text.contains("gc_weight, coverage_scaling_weight, count_scaling_weight, or flen"),
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
    // - Extra column `gc` is ignored, so GC/cw/fl tags are absent.
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 20, 60, '+', "fragment_1");

    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
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
    // - Known `flen` still maps to fl=30
    // - GC and cw tags are absent
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 40, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
    assert_eq!(aux_tags[0].fragment_length, Some(30));

    Ok(())
}

#[test]
fn given_supported_extra_column_names_when_run_then_gc_cov_and_flen_are_transferred_to_aux_tags()
-> Result<()> {
    // Arrange:
    // Header uses the three supported extra names exactly.
    // Hand-derived tag expectations:
    // Row 1 -> GC=0.25, cw=1.5, fl=10
    // Row 2 -> GC=None ("na"), cw=None ("."), fl=11
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tgc_weight\tcoverage_scaling_weight\tflen",
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr1", 0, 10, 60, '+', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 20, 31, 40, '-', "fragment_2");

    assert_eq!(aux_tags.len(), 2);
    assert_optional_f32_eq(aux_tags[0].gc_weight, Some(0.25), "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, Some(1.5), "cw");
    assert_eq!(aux_tags[0].fragment_length, Some(10));
    assert_optional_f32_eq(aux_tags[1].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[1].coverage_scaling_weight, None, "cw");
    assert_eq!(aux_tags[1].fragment_length, Some(11));
    assert_first_record_has_exact_aux_tags(&output_bam_path, &[b"GC", b"cw", b"fl"])?;
    assert_first_record_lacks_aux_tags(&output_bam_path, &[b"CO", b"FL"])?;

    Ok(())
}

#[test]
fn given_fragment_count_scaling_column_when_run_then_cnt_aux_tag_is_written() -> Result<()> {
    // Arrange:
    // Header includes the new recognized fragment-count scaling column.
    // Hand-derived tag expectations:
    // Row 1 -> cw=1.5, nw=0.5, fl=10
    // Row 2 -> cw=None, nw=2.0, fl=11
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tcoverage_scaling_weight\tcount_scaling_weight\tflen",
            "chr1\t0\t10\t60\t+\t1.5\t0.5\t10",
            "chr1\t20\t31\t40\t-\t.\t2.0\t11",
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr1", 0, 10, 60, '+', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 20, 31, 40, '-', "fragment_2");

    assert_eq!(aux_tags.len(), 2);
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, Some(1.5), "cw");
    assert_optional_f32_eq(aux_tags[0].count_scaling_weight, Some(0.5), "nw");
    assert_eq!(aux_tags[0].fragment_length, Some(10));
    assert_optional_f32_eq(aux_tags[1].coverage_scaling_weight, None, "cw");
    assert_optional_f32_eq(aux_tags[1].count_scaling_weight, Some(2.0), "nw");
    assert_eq!(aux_tags[1].fragment_length, Some(11));
    assert_first_record_has_exact_aux_tags(&output_bam_path, &[b"cw", b"nw", b"fl"])?;
    assert_first_record_lacks_aux_tags(&output_bam_path, &[b"CO", b"CN", b"FL"])?;

    Ok(())
}

#[test]
fn given_inline_header_with_only_count_scaling_when_run_then_only_cnt_aux_tag_is_written()
-> Result<()> {
    // Arrange:
    // Inline header defines only the count-scaling extra column.
    // Hand-derived expectations:
    // - Row 1 -> nw=0.5
    // - Row 2 -> nw absent because value is "."
    // - GC, cw, and fl remain absent for both rows
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(
        &frag_path,
        &[
            "chromosome\tstart\tend\tmapq\tstrand\tcount_scaling_weight",
            "chr1\t0\t10\t60\t+\t0.5",
            "chr1\t20\t31\t40\t-\t.",
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr1", 0, 10, 60, '+', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 20, 31, 40, '-', "fragment_2");

    assert_eq!(aux_tags.len(), 2);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
    assert_optional_f32_eq(aux_tags[0].count_scaling_weight, Some(0.5), "nw");
    assert_eq!(aux_tags[0].fragment_length, None);
    assert_optional_f32_eq(aux_tags[1].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[1].coverage_scaling_weight, None, "cw");
    assert_optional_f32_eq(aux_tags[1].count_scaling_weight, None, "nw");
    assert_eq!(aux_tags[1].fragment_length, None);

    Ok(())
}

#[test]
fn given_no_header_and_extra_columns_when_run_then_extra_columns_are_ignored_and_no_aux_tags_are_written()
-> Result<()> {
    // Arrange:
    // No inline header is present, no explicit header is configured, and no companion
    // header file exists for this filename. The parser therefore uses fixed 5-column
    // layout and ignores trailing columns.
    //
    // Hand-derived expectation:
    // - One fragment [10,20) is converted to one BAM record
    // - Trailing values `0.25 1.5 10` are ignored
    // - GC, cw, and fl tags are absent
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 20, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
    assert_eq!(aux_tags[0].fragment_length, None);

    Ok(())
}

#[test]
fn given_inline_header_with_only_flen_when_run_then_only_flen_aux_tag_is_written() -> Result<()> {
    // Arrange:
    // Inline header defines only one supported extra column (`flen`).
    // Hand-derived expectation:
    // - Row 1 has flen=80 so fl=80 is written
    // - Row 2 has flen="." so fl is absent
    // - GC and cw remain absent for both rows
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 2);
    assert_unpaired_full_match_record(&rows[0], "chr1", 0, 80, 60, '+', "fragment_1");
    assert_unpaired_full_match_record(&rows[1], "chr1", 100, 180, 55, '-', "fragment_2");

    assert_eq!(aux_tags.len(), 2);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
    assert_eq!(aux_tags[0].fragment_length, Some(80));
    assert_optional_f32_eq(aux_tags[1].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[1].coverage_scaling_weight, None, "cw");
    assert_eq!(aux_tags[1].fragment_length, None);

    Ok(())
}

#[test]
fn given_explicit_header_with_only_flen_when_run_then_only_flen_aux_tag_is_written() -> Result<()> {
    // Arrange:
    // Frag file has no inline header. We pass an explicit header file that maps column 6 to `flen`.
    // Hand-derived expectation:
    // - One fragment [10,60) with fl=50
    // - GC and cw remain absent
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 60, 42, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
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
        error_text.contains("gc_weight, coverage_scaling_weight, count_scaling_weight, or flen"),
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
    // - GC, cw, and fl tags are absent
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 20, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
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
    // - Known `flen` maps to fl=30
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 5, 35, 60, '-', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
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
    // - One record [20,70) with fl=50 from the 6th column
    // - GC and cw remain absent
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 20, 70, 39, '-', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
    assert_eq!(aux_tags[0].fragment_length, Some(50));

    Ok(())
}

#[test]
fn given_companion_header_with_only_count_scaling_when_run_then_only_cnt_aux_tag_is_written()
-> Result<()> {
    // Arrange:
    // No explicit header and no inline header. The companion header therefore defines the
    // count-scaling column.
    // Hand-derived expectations:
    // - one fragment [20,70) gets nw=0.5
    // - GC, cw, and fl remain absent
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("sample.frag.tsv");
    let companion_header_path = input_dir.path().join("sample.frag.header.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t20\t70\t39\t-\t0.5"])?;
    write_frag_file(
        &companion_header_path,
        &["chromosome\tstart\tend\tmapq\tstrand\tcount_scaling_weight"],
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 20, 70, 39, '-', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
    assert_optional_f32_eq(aux_tags[0].count_scaling_weight, Some(0.5), "nw");
    assert_eq!(aux_tags[0].fragment_length, None);

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
    // - GC, cw, and fl tags are absent
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 10, 20, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
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
    // - Known `flen` maps to fl=30
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 5, 35, 60, '-', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
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
        error_text.contains("gc_weight, coverage_scaling_weight, count_scaling_weight, or flen"),
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
    // - Fragment converts and fl=30 is written
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    let aux_tags = read_aux_tags(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 5, 35, 60, '-', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_eq!(aux_tags[0].fragment_length, Some(30));
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");

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
    // Any overlap excludes all fragments except [0,10).
    let kept_starts = run_blacklist_strategy_case(BlacklistStrategy::Any)?;
    assert_eq!(kept_starts, vec![0]);
    Ok(())
}

#[test]
fn given_blacklist_all_when_run_then_only_fully_overlapped_fragments_are_excluded() -> Result<()> {
    // Full overlap excludes only [10,20).
    let kept_starts = run_blacklist_strategy_case(BlacklistStrategy::All)?;
    assert_eq!(kept_starts, vec![0, 5, 15, 16]);
    Ok(())
}

#[test]
fn given_blacklist_midpoint_when_run_then_midpoint_overlap_controls_exclusion() -> Result<()> {
    // Midpoint strategy checks central-base support:
    // - [5,15) central bases 9 and 10 -> excluded
    // - [10,20) central bases 14 and 15 -> excluded
    // - [15,25) central bases 19 and 20 -> excluded because 19 lies in [10,20)
    // - [16,26) central bases 20 and 21 -> kept because [10,20) is end-exclusive
    let kept_starts = run_blacklist_strategy_case(BlacklistStrategy::Midpoint)?;
    assert_eq!(kept_starts, vec![0, 16]);
    Ok(())
}

#[test]
fn given_blacklist_proportion_when_run_then_threshold_controls_exclusion() -> Result<()> {
    // Overlap fractions are [0.0, 0.5, 1.0, 0.5, 0.4] for the five fragments.
    // With threshold 0.6, only [10,20) is excluded.
    let kept_starts = run_blacklist_strategy_case(BlacklistStrategy::Proportion(0.6))?;
    assert_eq!(kept_starts, vec![0, 5, 15, 16]);
    Ok(())
}

#[test]
fn given_touching_short_blacklists_in_separate_files_when_min_size_exceeds_each_then_filter_happens_before_merge()
-> Result<()> {
    // Arrange:
    // One fragment spans [0,10).
    //
    // Two blacklist files contain touching 5 bp intervals:
    //   file A -> [0,5)
    //   file B -> [5,10)
    //
    // We set `blacklist_min_size = 6`.
    //
    // The documented shared-loader contract is:
    // 1. discard intervals shorter than `min_size`
    // 2. then merge the remaining intervals
    //
    // So both 5 bp intervals must be discarded before any touching merge can happen, and the
    // fragment must survive.
    //
    // This specifically distinguishes the intended behavior from the wrong "merge first, then
    // filter" order:
    // - merge-first would create one [0,10) blacklist interval of size 10
    // - that merged interval would then pass the size filter and exclude the fragment
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");
    let blacklist_a = input_dir.path().join("blacklist_a.bed");
    let blacklist_b = input_dir.path().join("blacklist_b.bed");

    write_frag_file(&frag_path, &["chr1\t0\t10\t60\t+"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;
    write_blacklist_bed(&blacklist_a, &[("chr1", 0, 5)])?;
    write_blacklist_bed(&blacklist_b, &[("chr1", 5, 10)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_blacklist(Some(vec![blacklist_a, blacklist_b]));
    config.set_blacklist_min_size(6);
    config.set_blacklist_strategy(BlacklistStrategy::Any);
    config.set_min_mapq(0);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    // Act
    run_frag_to_bam(&config)?;
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;

    // Assert
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 0, 10, 60, '+', "fragment_1");

    Ok(())
}

#[test]
fn given_same_blacklist_premerged_in_one_file_when_min_size_is_met_then_fragment_is_excluded()
-> Result<()> {
    // Arrange:
    // Use the same logical masked region as the previous test, but this time pre-merge it in the
    // BED itself:
    //   [0,10)
    //
    // With `blacklist_min_size = 6`, this single interval is kept because its own length is 10.
    // Under the default `any` strategy it therefore excludes the only fragment [0,10).
    //
    // Paired with the previous test, this locks down the exact shared-loader contract:
    // - split short intervals in separate files are filtered before merge and therefore disappear
    // - an already merged long interval survives the size filter and excludes the fragment
    let input_dir = TempDir::new()?;
    let output_dir = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");
    let blacklist_path = input_dir.path().join("blacklist_merged.bed");

    write_frag_file(&frag_path, &["chr1\t0\t10\t60\t+"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;
    write_blacklist_bed(&blacklist_path, &[("chr1", 0, 10)])?;

    let mut config = make_config(
        frag_path,
        output_dir.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    config.set_blacklist(Some(vec![blacklist_path]));
    config.set_blacklist_min_size(6);
    config.set_blacklist_strategy(BlacklistStrategy::Any);
    config.set_min_mapq(0);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    // Act
    run_frag_to_bam(&config)?;

    // Assert:
    // `frag-to-bam` only raises "No fragments passed filters; no BAM to write" when no
    // chromosome is observed at all. Here chr1 is observed before the blacklist removes the only
    // fragment, so the command succeeds and writes an empty BAM.
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    assert!(
        rows.is_empty(),
        "expected the premerged blacklist to filter out the only fragment, got {rows:?}"
    );

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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
    let rows = read_bam_rows(&output_bam_path)?;
    assert!(rows.is_empty(), "expected no records in output BAM");

    Ok(())
}

#[test]
fn default_min_mapq_matches_explicit_zero_and_differs_from_explicit_twenty() -> Result<()> {
    // Arrange:
    // `frag-to-bam` intentionally defaults to `min_mapq = 0`.
    // Use two frag rows on chr1:
    // - [10,20) with MAPQ 10
    // - [30,45) with MAPQ 30
    //
    // Therefore:
    // - default config must keep both rows
    // - explicit `min_mapq = 0` must match exactly
    // - explicit `min_mapq = 20` must drop only the MAPQ-10 row
    let input_dir = TempDir::new()?;
    let out_default = TempDir::new()?;
    let out_zero = TempDir::new()?;
    let out_twenty = TempDir::new()?;
    let frag_path = input_dir.path().join("input.frag.tsv");
    let chrom_sizes_path = input_dir.path().join("chrom.sizes");

    write_frag_file(&frag_path, &["chr1\t10\t20\t10\t+", "chr1\t30\t45\t30\t-"])?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 100)])?;

    let default_config = make_config(
        frag_path.clone(),
        out_default.path().to_path_buf(),
        chrom_sizes_path.clone(),
        base_chromosomes(&["chr1"]),
    );
    let mut explicit_zero_config = make_config(
        frag_path.clone(),
        out_zero.path().to_path_buf(),
        chrom_sizes_path.clone(),
        base_chromosomes(&["chr1"]),
    );
    explicit_zero_config.set_min_mapq(0);
    let mut explicit_twenty_config = make_config(
        frag_path,
        out_twenty.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    explicit_twenty_config.set_min_mapq(20);

    // Act
    run_frag_to_bam(&default_config)?;
    run_frag_to_bam(&explicit_zero_config)?;
    run_frag_to_bam(&explicit_twenty_config)?;

    // Assert:
    // `frag-to-bam` writes one unpaired BAM record per kept row. The row fields map directly to:
    // - chromosome/start/end
    // - mapq
    // - strand
    // So default and explicit-zero runs must contain both rows, while `min_mapq = 20` keeps only
    // the second row.
    //
    // Qname derivation is separate from input line number: the implementation renumbers kept
    // fragments during the second pass as `fragment_1`, `fragment_2`, ... in write order.
    // Therefore the single surviving row in the `min_mapq = 20` run becomes `fragment_1`.
    let default_rows = read_bam_rows(&output_bam_path(out_default.path(), "restored"))?;
    let explicit_zero_rows = read_bam_rows(&output_bam_path(out_zero.path(), "restored"))?;
    let explicit_twenty_rows = read_bam_rows(&output_bam_path(out_twenty.path(), "restored"))?;

    assert_eq!(default_rows, explicit_zero_rows);
    assert_eq!(default_rows.len(), 2);
    assert_unpaired_full_match_record(&default_rows[0], "chr1", 10, 20, 10, '+', "fragment_1");
    assert_unpaired_full_match_record(&default_rows[1], "chr1", 30, 45, 30, '-', "fragment_2");

    assert_eq!(explicit_twenty_rows.len(), 1);
    assert_unpaired_full_match_record(
        &explicit_twenty_rows[0],
        "chr1",
        30,
        45,
        30,
        '-',
        "fragment_1",
    );

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

    write_frag_file(&frag_path, &["chr1\t10\t20\t60\t+", "chr2\t5\t15\t30\t-"])?;
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
    let output_bam_path = output_bam_path(output_dir.path(), "restored");
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
fn roundtrip_simple_inward_to_unpaired_bam() -> Result<(fixtures::BamFixture, TempDir, PathBuf)> {
    let source_bam = fixtures::simple_inward_bam()?;
    let bam_to_frag_output = TempDir::new()?;
    let frag_to_bam_output = TempDir::new()?;
    let chrom_sizes_path = frag_to_bam_output.path().join("chrom.sizes");
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 200)])?;

    let mut bam_to_frag_config = BamToFragConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: bam_to_frag_output.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    bam_to_frag_config.set_output_prefix("roundtrip");
    bam_to_frag_config.set_min_mapq(0);
    {
        let fragment_lengths = bam_to_frag_config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 200;
    }
    run_bam_to_frag(&bam_to_frag_config)?;

    let frag_path = bam_to_frag_output.path().join("roundtrip.frag.tsv.gz");

    let mut frag_to_bam_config = make_config(
        frag_path,
        frag_to_bam_output.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    frag_to_bam_config.set_min_mapq(0);
    {
        let fragment_lengths = frag_to_bam_config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 200;
    }
    run_frag_to_bam(&frag_to_bam_config)?;

    let restored_bam_path = output_bam_path(frag_to_bam_output.path(), "restored");
    // Downstream counting commands fetch by genomic region and therefore require an index.
    // `frag-to-bam` intentionally does not create one, so this roundtrip helper adds it
    // explicitly to keep the test focused on roundtrip scientific equivalence.
    build_bai_for_test_bam(&restored_bam_path)?;
    Ok((source_bam, frag_to_bam_output, restored_bam_path))
}

#[cfg(feature = "cmd_bam_to_frag")]
fn roundtrip_single_paired_fragment_to_unpaired_bam(
    name: &str,
    fragment_start: i64,
    fragment_len: i64,
    read_len: i64,
) -> Result<(fixtures::BamFixture, TempDir, PathBuf)> {
    let source_bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment(fragment_start, fragment_len, read_len)],
        Vec::new(),
        name,
    )?;
    let bam_to_frag_output = TempDir::new()?;
    let frag_to_bam_output = TempDir::new()?;
    let chrom_sizes_path = frag_to_bam_output.path().join("chrom.sizes");
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 200)])?;

    let mut bam_to_frag_config = BamToFragConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: bam_to_frag_output.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    bam_to_frag_config.set_output_prefix("roundtrip");
    bam_to_frag_config.set_min_mapq(0);
    {
        let fragment_lengths = bam_to_frag_config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 200;
    }
    run_bam_to_frag(&bam_to_frag_config)?;

    let frag_path = bam_to_frag_output.path().join("roundtrip.frag.tsv.gz");

    let mut frag_to_bam_config = make_config(
        frag_path,
        frag_to_bam_output.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    frag_to_bam_config.set_min_mapq(0);
    {
        let fragment_lengths = frag_to_bam_config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 200;
    }
    run_frag_to_bam(&frag_to_bam_config)?;

    let restored_bam_path = output_bam_path(frag_to_bam_output.path(), "restored");
    build_bai_for_test_bam(&restored_bam_path)?;
    Ok((source_bam, frag_to_bam_output, restored_bam_path))
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
    let restored_bam_path = output_bam_path(frag_to_bam_output.path(), "restored");
    let restored_rows = read_bam_rows(&restored_bam_path)?;

    // Assert 2: all available fields are restored exactly
    assert_eq!(restored_rows.len(), 2);
    assert_unpaired_full_match_record(&restored_rows[0], "chr1", 100, 180, 45, '+', "fragment_1");
    assert_unpaired_full_match_record(&restored_rows[1], "chr2", 200, 280, 30, '-', "fragment_2");

    Ok(())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_lengths"))]
#[test]
fn given_bam_to_frag_then_frag_to_bam_when_counting_lengths_then_roundtrip_matches_original()
-> Result<()> {
    // Arrange:
    // `simple_inward_bam` contains exactly one inward fragment on chr1:
    // - forward read covers [20, 40)
    // - reverse read covers [60, 80)
    // - fragment span is therefore [20, 80)
    // - fragment length is 80 - 20 = 60
    //
    // `bam-to-frag` emits the row:
    //   chr1  20  80  60  +
    //
    // `frag-to-bam` restores that as one unpaired BAM record spanning [20, 80).
    //
    // `lengths` defines fragment length as:
    // - paired mode: end(reverse) - start(forward)
    // - reads_are_fragments mode: end(read) - start(read)
    //
    // So both inputs must yield one global-window count in the length-60 bin only.
    let (source_bam, _restored_dir, restored_bam_path) = roundtrip_simple_inward_to_unpaired_bam()?;

    let original_out = TempDir::new()?;
    let restored_out = TempDir::new()?;

    let mut original_cfg = LengthsConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: original_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    original_cfg.set_min_mapq(0);
    original_cfg.set_require_proper_pair(false);
    original_cfg.set_window_assignment(AssignToWindowArgs::default());
    original_cfg.set_per_bp_length_bins(10, 100);

    let mut restored_cfg = LengthsConfig::new(
        IOCArgs {
            bam: restored_bam_path.clone(),
            output_dir: restored_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    restored_cfg.set_min_mapq(0);
    restored_cfg.unpaired.reads_are_fragments = true;
    restored_cfg.set_require_proper_pair(false);
    restored_cfg.set_window_assignment(AssignToWindowArgs::default());
    restored_cfg.set_per_bp_length_bins(10, 100);

    // Act
    run_lengths(&original_cfg)?;
    run_lengths(&restored_cfg)?;

    let original_counts = read_length_counts_tsv(&original_out.path().join(dot_join(&[
        original_cfg.output_prefix.trim(),
        "length_counts.tsv.zst",
    ])))?;
    let restored_counts = read_length_counts_tsv(&restored_out.path().join(dot_join(&[
        restored_cfg.output_prefix.trim(),
        "length_counts.tsv.zst",
    ])))?;

    // Assert:
    // One global window and lengths 10..=100 give shape (1, 91).
    // Length 60 is column 60 - 10 = 50.
    assert_eq!(original_counts.dim(), (1, 91));
    assert_eq!(restored_counts.dim(), (1, 91));
    assert_eq!(original_counts, restored_counts);

    let len60_idx = 60 - 10;
    assert!((original_counts[(0, len60_idx)] - 1.0).abs() < 1e-12);
    for idx in 0..original_counts.ncols() {
        if idx == len60_idx {
            continue;
        }
        assert!(
            original_counts[(0, idx)].abs() < 1e-12,
            "expected only length 60 to be occupied, but column {idx} had {}",
            original_counts[(0, idx)]
        );
    }

    Ok(())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_lengths"))]
#[test]
fn given_bam_to_frag_then_frag_to_bam_when_counting_lengths_with_blacklist_then_roundtrip_matches_original()
-> Result<()> {
    // Arrange:
    // `simple_inward_bam()` contains one fragment spanning [20, 80), length 60.
    // We apply a blacklist interval [25, 35), which overlaps that fragment.
    //
    // `lengths` defaults to blacklist strategy `any`, so a fragment is excluded entirely as soon as
    // any part of its span overlaps the blacklist.
    //
    // Therefore both representations of the same physical fragment:
    // - original paired BAM
    // - roundtripped unpaired BAM in `reads_are_fragments` mode
    // must yield the same global length-count array:
    // - shape (1, 91) for lengths 10..=100
    // - all zeros, because the only fragment is blacklisted away
    let (source_bam, _restored_dir, restored_bam_path) = roundtrip_simple_inward_to_unpaired_bam()?;
    let original_out = TempDir::new()?;
    let restored_out = TempDir::new()?;
    let blacklist_path = original_out.path().join("blacklist.bed");
    write_blacklist_bed(&blacklist_path, &[("chr1", 25, 35)])?;

    let mut original_cfg = LengthsConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: original_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    original_cfg.set_min_mapq(0);
    original_cfg.set_require_proper_pair(false);
    original_cfg.set_window_assignment(AssignToWindowArgs::default());
    original_cfg.blacklist = Some(vec![blacklist_path.clone()]);
    original_cfg.set_per_bp_length_bins(10, 100);

    let mut restored_cfg = LengthsConfig::new(
        IOCArgs {
            bam: restored_bam_path,
            output_dir: restored_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    restored_cfg.set_min_mapq(0);
    restored_cfg.unpaired.reads_are_fragments = true;
    restored_cfg.set_require_proper_pair(false);
    restored_cfg.set_window_assignment(AssignToWindowArgs::default());
    restored_cfg.blacklist = Some(vec![blacklist_path]);
    restored_cfg.set_per_bp_length_bins(10, 100);

    // Act
    run_lengths(&original_cfg)?;
    run_lengths(&restored_cfg)?;

    // Assert
    let original_counts = read_length_counts_tsv(&original_out.path().join(dot_join(&[
        original_cfg.output_prefix.trim(),
        "length_counts.tsv.zst",
    ])))?;
    let restored_counts = read_length_counts_tsv(&restored_out.path().join(dot_join(&[
        restored_cfg.output_prefix.trim(),
        "length_counts.tsv.zst",
    ])))?;

    assert_eq!(original_counts.dim(), (1, 91));
    assert_eq!(restored_counts.dim(), (1, 91));
    assert_eq!(original_counts, restored_counts);
    assert!(
        original_counts.iter().all(|value| value.abs() < 1e-12),
        "the only fragment overlaps the blacklist and should be removed entirely"
    );

    Ok(())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_fcoverage"))]
#[test]
fn given_bam_to_frag_then_frag_to_bam_when_counting_coverage_then_roundtrip_matches_original()
-> Result<()> {
    // Arrange:
    // The same `simple_inward_bam` roundtrip represents exactly one fragment/read spanning [20, 80).
    //
    // `fcoverage` counts:
    // - paired input by the fragment span [20, 80)
    // - unpaired `reads_are_fragments` input by the read span [20, 80)
    //
    // Therefore both inputs must yield the same positional bedGraph:
    //   chr1  20  80  1
    let (source_bam, _restored_dir, restored_bam_path) = roundtrip_simple_inward_to_unpaired_bam()?;

    let original_out = TempDir::new()?;
    let restored_out = TempDir::new()?;

    let mut original_cfg = FCoverageConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: original_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    original_cfg.set_output_prefix("origcov");
    original_cfg.set_tile_size(1_000);
    original_cfg.set_min_mapq(0);
    original_cfg.set_require_proper_pair(false);
    original_cfg.set_ignore_gap(false);
    {
        let fragment_lengths = original_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    let mut restored_cfg = FCoverageConfig::new(
        IOCArgs {
            bam: restored_bam_path.clone(),
            output_dir: restored_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    restored_cfg.set_output_prefix("restoredcov");
    restored_cfg.set_tile_size(1_000);
    restored_cfg.set_min_mapq(0);
    restored_cfg.set_require_proper_pair(false);
    restored_cfg.unpaired.reads_are_fragments = true;
    {
        let fragment_lengths = restored_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    // Act
    run_fcoverage(&original_cfg)?;
    run_fcoverage(&restored_cfg)?;

    let original_text = read_zst_to_string(
        &original_out
            .path()
            .join("origcov.fcoverage.per_position.bedgraph.zst"),
    )?;
    let restored_text = read_zst_to_string(
        &restored_out
            .path()
            .join("restoredcov.fcoverage.per_position.bedgraph.zst"),
    )?;

    // Assert
    assert_eq!(original_text, restored_text);
    assert_eq!(original_text.trim(), "chr1\t20\t80\t1");

    Ok(())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_fcoverage"))]
#[test]
fn given_bam_to_frag_then_frag_to_bam_when_counting_coverage_with_blacklist_then_roundtrip_matches_original()
-> Result<()> {
    // Arrange:
    // `simple_inward_bam()` roundtrips to one unpaired fragment/read span [20, 80).
    // We blacklist [25, 35), which lies inside that covered span.
    //
    // `fcoverage` masks blacklisted positions out of the positional output rather than dropping the
    // whole fragment. So the surviving coverage must split into two runs:
    // - [20, 25) -> 1
    // - [35, 80) -> 1
    //
    // Because both the original paired input and the roundtripped unpaired input describe the same
    // fragment span, they must produce the same masked bedGraph.
    let (source_bam, _restored_dir, restored_bam_path) = roundtrip_simple_inward_to_unpaired_bam()?;
    let original_out = TempDir::new()?;
    let restored_out = TempDir::new()?;
    let blacklist_path = original_out.path().join("blacklist.bed");
    write_blacklist_bed(&blacklist_path, &[("chr1", 25, 35)])?;

    let mut original_cfg = FCoverageConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: original_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    original_cfg.set_output_prefix("origcov");
    original_cfg.set_tile_size(1_000);
    original_cfg.set_min_mapq(0);
    original_cfg.set_require_proper_pair(false);
    original_cfg.set_ignore_gap(false);
    original_cfg.blacklist = Some(vec![blacklist_path.clone()]);
    {
        let fragment_lengths = original_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    let mut restored_cfg = FCoverageConfig::new(
        IOCArgs {
            bam: restored_bam_path,
            output_dir: restored_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    restored_cfg.set_output_prefix("restoredcov");
    restored_cfg.set_tile_size(1_000);
    restored_cfg.set_min_mapq(0);
    restored_cfg.set_require_proper_pair(false);
    restored_cfg.unpaired.reads_are_fragments = true;
    restored_cfg.blacklist = Some(vec![blacklist_path]);
    {
        let fragment_lengths = restored_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    // Act
    run_fcoverage(&original_cfg)?;
    run_fcoverage(&restored_cfg)?;

    let original_text = read_zst_to_string(
        &original_out
            .path()
            .join("origcov.fcoverage.per_position.bedgraph.zst"),
    )?;
    let restored_text = read_zst_to_string(
        &restored_out
            .path()
            .join("restoredcov.fcoverage.per_position.bedgraph.zst"),
    )?;

    // Assert
    assert_eq!(original_text, restored_text);
    assert_eq!(original_text.trim(), "chr1\t20\t25\t1\nchr1\t35\t80\t1");

    Ok(())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_bam_to_bam"))]
#[test]
fn bam_frag_bam_roundtrip_preserves_coverage_tags_for_same_span() -> Result<()> {
    // Arrange:
    // This test compares the same physical fragment represented in two different ways:
    // - the original paired BAM from `simple_inward_bam()`
    // - a restored unpaired BAM after BAM -> FRAG -> BAM roundtrip
    //
    // We test this because FRAG does not preserve pair structure, so the restored BAM will have
    // one unpaired record instead of two mates. Even so, downstream fragment-level tags should
    // stay identical when the underlying span is the same.
    //
    // In both cases the fragment span is [20, 80), length 60.
    //
    // We then run `bam-to-bam` on both representations with a non-uniform scaling TSV:
    // - [0, 40)  factor 2.0  -> contributes 20 bp over [20, 40)
    // - [40, 80) factor 1.0  -> contributes 40 bp over [40, 80)
    // - [80,200) factor 1.0  -> not touched
    //
    // Hand-derived full-fragment average over [20, 80):
    //   (20 * 2.0 + 40 * 1.0) / 60 = 80 / 60 = 4/3
    //
    // `bam-to-bam` always tags each emitted record with fragment-level values, so the scientific
    // contract is:
    // - original paired BAM -> 2 records, each with cw = 4/3 and fl = 60
    // - restored unpaired BAM -> 1 record, with cw = 4/3 and fl = 60
    //
    // Record counts differ because pair structure is not representable through FRAG, but the
    // fragment tags must stay identical for the same physical span.
    let (source_bam, _restored_dir, restored_bam_path) = roundtrip_simple_inward_to_unpaired_bam()?;
    let work = TempDir::new()?;
    let original_out = work.path().join("original_tagged.bam");
    let restored_out = work.path().join("restored_tagged.bam");
    let scaling_path = work.path().join("piecewise_scaling.tsv");
    write_scaling_tsv(
        &scaling_path,
        &[
            ("chr1", 0, 40, 2.0),
            ("chr1", 40, 80, 1.0),
            ("chr1", 80, 200, 1.0),
        ],
    )?;

    let mut original_cfg = BamToBamConfig::new(
        source_bam.bam.clone(),
        original_out.clone(),
        base_chromosomes(&["chr1"]),
    );
    original_cfg.set_coverage_scaling_factors(Some(scaling_path.clone()));
    original_cfg.min_mapq = 0;
    {
        let fragment_lengths = original_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    let mut restored_cfg = BamToBamConfig::new(
        restored_bam_path,
        restored_out.clone(),
        base_chromosomes(&["chr1"]),
    );
    restored_cfg.set_coverage_scaling_factors(Some(scaling_path));
    restored_cfg.min_mapq = 0;
    restored_cfg.unpaired.reads_are_fragments = true;
    {
        let fragment_lengths = restored_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    // Act
    run_bam_to_bam(&original_cfg)?;
    run_bam_to_bam(&restored_cfg)?;

    // Assert
    let original_tags = read_aux_tags(&original_out)?;
    let restored_tags = read_aux_tags(&restored_out)?;
    let expected_cov = Some(4.0_f32 / 3.0_f32);

    assert_eq!(original_tags.len(), 2);
    assert_eq!(restored_tags.len(), 1);
    for tags in &original_tags {
        assert_optional_f32_eq(tags.coverage_scaling_weight, expected_cov, "paired cw");
        assert_optional_f32_eq(tags.gc_weight, None, "paired GC");
        assert_eq!(tags.fragment_length, Some(60));
    }
    assert_optional_f32_eq(
        restored_tags[0].coverage_scaling_weight,
        expected_cov,
        "restored cw",
    );
    assert_optional_f32_eq(restored_tags[0].gc_weight, None, "restored GC");
    assert_eq!(restored_tags[0].fragment_length, Some(60));

    Ok(())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_bam_to_bam"))]
#[test]
fn bam_frag_bam_roundtrip_preserves_count_tags_for_same_span() -> Result<()> {
    // Arrange:
    // This is the count-scaling mirror of the coverage-scaling roundtrip test above.
    // It again compares the same physical fragment represented as:
    // - the original paired BAM
    // - a restored unpaired BAM after BAM -> FRAG -> BAM
    //
    // We test this because FRAG roundtrips intentionally collapse pair structure, but count
    // scaling in downstream `bam-to-bam` should still depend only on the fragment span itself.
    //
    // In both inputs the fragment span is [20, 80), length 60.
    //
    // Use a non-uniform count-scaling TSV:
    // - [0, 40)  factor 2.0  -> contributes 20 bp over [20, 40)
    // - [40, 80) factor 1.0  -> contributes 40 bp over [40, 80)
    //
    // Hand-derived full-fragment average over [20, 80):
    //   (20 * 2.0 + 40 * 1.0) / 60 = 4/3
    //
    // The scientific contract is therefore:
    // - original paired BAM -> 2 records, each with nw = 4/3 and fl = 60
    // - restored unpaired BAM -> 1 record, with nw = 4/3 and fl = 60
    // - cw stays absent in both outputs because only count scaling is configured
    let (source_bam, _restored_dir, restored_bam_path) = roundtrip_simple_inward_to_unpaired_bam()?;
    let work = TempDir::new()?;
    let original_out = work.path().join("original_count_tagged.bam");
    let restored_out = work.path().join("restored_count_tagged.bam");
    let scaling_path = work.path().join("piecewise_count_scaling.tsv");
    write_scaling_tsv(
        &scaling_path,
        &[
            ("chr1", 0, 40, 2.0),
            ("chr1", 40, 80, 1.0),
            ("chr1", 80, 200, 1.0),
        ],
    )?;

    let mut original_cfg = BamToBamConfig::new(
        source_bam.bam.clone(),
        original_out.clone(),
        base_chromosomes(&["chr1"]),
    );
    original_cfg.set_count_scaling_factors(Some(scaling_path.clone()));
    original_cfg.min_mapq = 0;
    {
        let fragment_lengths = original_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    let mut restored_cfg = BamToBamConfig::new(
        restored_bam_path,
        restored_out.clone(),
        base_chromosomes(&["chr1"]),
    );
    restored_cfg.set_count_scaling_factors(Some(scaling_path));
    restored_cfg.min_mapq = 0;
    restored_cfg.unpaired.reads_are_fragments = true;
    {
        let fragment_lengths = restored_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 100;
    }

    // Act
    run_bam_to_bam(&original_cfg)?;
    run_bam_to_bam(&restored_cfg)?;

    // Assert
    let original_tags = read_aux_tags(&original_out)?;
    let restored_tags = read_aux_tags(&restored_out)?;
    let expected_cnt = Some(4.0_f32 / 3.0_f32);

    assert_eq!(original_tags.len(), 2);
    assert_eq!(restored_tags.len(), 1);
    for tags in &original_tags {
        assert_optional_f32_eq(tags.count_scaling_weight, expected_cnt, "paired nw");
        assert_optional_f32_eq(tags.coverage_scaling_weight, None, "paired cw");
        assert_optional_f32_eq(tags.gc_weight, None, "paired GC");
        assert_eq!(tags.fragment_length, Some(60));
    }
    assert_optional_f32_eq(
        restored_tags[0].count_scaling_weight,
        expected_cnt,
        "restored nw",
    );
    assert_optional_f32_eq(
        restored_tags[0].coverage_scaling_weight,
        None,
        "restored cw",
    );
    assert_optional_f32_eq(restored_tags[0].gc_weight, None, "restored GC");
    assert_eq!(restored_tags[0].fragment_length, Some(60));

    Ok(())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
#[test]
fn given_bam_to_frag_then_frag_to_bam_when_counting_midpoints_then_roundtrip_matches_original()
-> Result<()> {
    // Arrange:
    // We use one odd-length fragment so the midpoint is deterministic across both
    // representations:
    // - paired BAM fragment span [20, 81), length 61
    // - roundtripped unpaired BAM read span [20, 81), length 61
    //
    // Hand-derived midpoint:
    //   20 + floor(61 / 2) = 50
    // For one BED window [45, 56), that midpoint lands at profile position:
    //   50 - 45 = 5
    //
    // Therefore both the original paired BAM and the FRAG->BAM roundtrip must yield the same
    // midpoint profile array:
    // - shape [1 group, 1 length-bin, 11 positions]
    // - exactly one count at [0, 0, 5]
    let (source_bam, _restored_dir, restored_bam_path) =
        roundtrip_single_paired_fragment_to_unpaired_bam(
            "frag_to_bam_midpoints_roundtrip",
            20,
            61,
            20,
        )?;
    let original_out = TempDir::new()?;
    let restored_out = TempDir::new()?;
    let intervals = original_out.path().join("sites.bed");
    fs::write(&intervals, "chr1\t45\t56\tgroupA\n")?;

    let make_cfg = |bam_path: &Path, out_dir: &Path, unpaired: bool, prefix: &str| {
        let mut cfg = MidpointsConfig::new(
            IOCArgs {
                bam: bam_path.to_path_buf(),
                output_dir: out_dir.to_path_buf(),
                n_threads: 1,
            },
            base_chromosomes(&["chr1"]),
            intervals.clone(),
        );
        cfg.set_output_prefix(prefix);
        cfg.set_length_bins(vec![61, 62]);
        cfg.set_smoothing(MidpointSmoothing::None);
        cfg.set_tile_size(1_000);
        cfg.set_min_mapq(0);
        cfg.set_require_proper_pair(false);
        cfg.unpaired.reads_are_fragments = unpaired;
        cfg
    };

    let original_cfg = make_cfg(&source_bam.bam, original_out.path(), false, "origmid");
    let restored_cfg = make_cfg(&restored_bam_path, restored_out.path(), true, "restoredmid");

    // Act
    run_midpoints(&original_cfg)?;
    run_midpoints(&restored_cfg)?;

    // Assert
    let original_arr: Array3<f32> =
        read_midpoint_zarr_counts(original_out.path().join("origmid.midpoint_profiles.zarr"))?;
    let restored_arr: Array3<f32> = read_midpoint_zarr_counts(
        restored_out
            .path()
            .join("restoredmid.midpoint_profiles.zarr"),
    )?;
    let original_groups = fs::read_to_string(original_out.path().join("origmid.group_index.tsv"))?;
    let restored_groups =
        fs::read_to_string(restored_out.path().join("restoredmid.group_index.tsv"))?;

    assert_eq!(original_arr, restored_arr);
    assert_eq!(original_arr.shape(), &[1, 1, 11]);
    assert_eq!(original_arr[[0, 0, 5]], 1.0);
    assert_eq!(original_arr.sum(), 1.0);
    assert_eq!(original_groups, restored_groups);
    assert_eq!(
        original_groups.trim(),
        "group_idx\tgroup_name\teligible_intervals\n0\tgroupA\t1"
    );

    Ok(())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
#[test]
fn given_bam_to_frag_gc_weights_then_frag_to_bam_then_midpoints_gc_tag_matches_original_gc_file_weighting()
-> Result<()> {
    // Arrange:
    // Start from one odd-length paired fragment spanning [20, 81), length 61. Its midpoint is
    // deterministic:
    //   20 + floor(61 / 2) = 50
    // and one window [45, 56) therefore receives the fragment at profile position:
    //   50 - 45 = 5
    //
    // We build the smallest GC package that assigns a constant weight 3.0 to every 61 bp
    // fragment, independent of GC percentage:
    // - length_edges = [61, 62]
    // - gc_edges     = [0, 101]
    // - correction_matrix = [[3.0]]
    //
    // Then we compare two logically equivalent released workflows:
    // 1. original paired BAM -> `midpoints --gc-file <pkg>`
    // 2. original paired BAM -> `bam-to-frag --gc-file <pkg>` ->
    //    `frag-to-bam` (auto-detect companion header, restore `GC` tag) ->
    //    `midpoints --gc-tag GC`
    //
    // Because the package weight is constant 3.0 for the only supported fragment length, both
    // workflows must produce the exact same midpoint profile:
    // - shape [1, 1, 11]
    // - exactly 3.0 at position 5
    // - total mass 3.0
    let source_bam = bam_from_specs(
        vec![("chr1".to_string(), 200)],
        vec![paired_fragment(20, 61, 20)],
        Vec::new(),
        "frag_to_bam_gc_roundtrip_source",
    )?;
    let reference = simple_reference_twobit()?;
    let work = TempDir::new()?;
    let bam_to_frag_out = TempDir::new()?;
    let frag_to_bam_out = TempDir::new()?;
    let original_midpoints_out = TempDir::new()?;
    let restored_midpoints_out = TempDir::new()?;
    let gc_path = work.path().join("constant_gc_pkg.zarr");
    let intervals = work.path().join("sites.bed");
    let chrom_sizes_path = frag_to_bam_out.path().join("chrom.sizes");
    fs::write(&intervals, "chr1\t45\t56\tgroupA\n")?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 200)])?;

    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![61, 62],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&reference.path)?,
        correction_matrix: array![[3.0_f64]],
    };
    package.write_zarr(&gc_path)?;

    let mut bam_to_frag_cfg = BamToFragConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: bam_to_frag_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    bam_to_frag_cfg.set_output_prefix("weighted");
    bam_to_frag_cfg.set_min_mapq(0);
    bam_to_frag_cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path.clone()),
        neutralize_invalid_gc: false,
    });
    bam_to_frag_cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let fragment_lengths = bam_to_frag_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 61;
        fragment_lengths.max_fragment_length = 61;
    }

    let mut original_midpoints_cfg = MidpointsConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: original_midpoints_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        intervals.clone(),
    );
    original_midpoints_cfg.set_output_prefix("origsites");
    original_midpoints_cfg.set_length_bins(vec![61, 62]);
    original_midpoints_cfg.set_smoothing(MidpointSmoothing::None);
    original_midpoints_cfg.set_tile_size(1_000);
    original_midpoints_cfg.set_min_mapq(0);
    original_midpoints_cfg.set_require_proper_pair(false);
    original_midpoints_cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path.clone()),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    original_midpoints_cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act 1: create a FRAG file with `gc_weight` and restore it back to BAM with `GC` tags.
    run_bam_to_frag(&bam_to_frag_cfg)?;
    let frag_path = bam_to_frag_out.path().join("weighted.frag.tsv.gz");
    let mut frag_to_bam_cfg = make_config(
        frag_path,
        frag_to_bam_out.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    frag_to_bam_cfg.set_min_mapq(0);
    {
        let fragment_lengths = frag_to_bam_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 61;
        fragment_lengths.max_fragment_length = 61;
    }
    run_frag_to_bam(&frag_to_bam_cfg)?;
    let restored_bam_path = output_bam_path(frag_to_bam_out.path(), "restored");
    build_bai_for_test_bam(&restored_bam_path)?;

    // Act 2: consume the original BAM via `--gc-file` and the restored BAM via `--gc-tag`.
    run_midpoints(&original_midpoints_cfg)?;

    let mut restored_midpoints_cfg = MidpointsConfig::new(
        IOCArgs {
            bam: restored_bam_path,
            output_dir: restored_midpoints_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        intervals,
    );
    restored_midpoints_cfg.set_output_prefix("restoredsites");
    restored_midpoints_cfg.set_length_bins(vec![61, 62]);
    restored_midpoints_cfg.set_smoothing(MidpointSmoothing::None);
    restored_midpoints_cfg.set_tile_size(1_000);
    restored_midpoints_cfg.set_min_mapq(0);
    restored_midpoints_cfg.set_require_proper_pair(false);
    restored_midpoints_cfg.unpaired.reads_are_fragments = true;
    restored_midpoints_cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    });
    run_midpoints(&restored_midpoints_cfg)?;

    // Assert
    let original_arr: Array3<f32> = read_midpoint_zarr_counts(
        original_midpoints_out
            .path()
            .join("origsites.midpoint_profiles.zarr"),
    )?;
    let restored_arr: Array3<f32> = read_midpoint_zarr_counts(
        restored_midpoints_out
            .path()
            .join("restoredsites.midpoint_profiles.zarr"),
    )?;
    let restored_tags = read_aux_tags(&output_bam_path(frag_to_bam_out.path(), "restored"))?;

    assert_eq!(restored_tags.len(), 1);
    assert_optional_f32_eq(restored_tags[0].gc_weight, Some(3.0), "restored GC");
    assert_eq!(original_arr, restored_arr);
    assert_eq!(original_arr.shape(), &[1, 1, 11]);
    assert_eq!(original_arr[[0, 0, 5]], 3.0);
    assert_eq!(original_arr.sum(), 3.0);

    Ok(())
}

#[cfg(all(feature = "cmd_bam_to_frag", feature = "cmd_midpoints"))]
#[test]
fn given_bam_to_frag_real_non_neutral_gc_then_frag_to_bam_then_midpoints_gc_tag_matches_original_gc_file_weighting()
-> Result<()> {
    // Arrange:
    // Use a real non-neutral `ref-gc-bias -> gc-bias` package rather than a handcrafted constant
    // matrix.
    //
    // Reference genome:
    // - chr1[0,100)    = all A
    // - chr1[100,101)  = one-base spacer excluded from the reference BED windows
    // - chr1[101,201)  = all C
    //
    // Real package derivation:
    // - fragment length is fixed at 61
    // - reference windows [0,100) and [101,201) count only starts that both lie in the window
    //   and leave room for the full 61 bp fragment, so the counted support is:
    //     starts 0..=39    -> GC%=0
    //     starts 101..=140 -> GC%=100
    // - sample producer BAM has one A-only fragment and nine C-only fragments
    // - resulting weights are:
    //     GC%=0   -> 5.0
    //     GC%=100 -> 5/9
    //
    // Consumer BAM:
    // - one A-only fragment [10,71), midpoint 40
    // - one C-only fragment [110,171), midpoint 140
    // - windows [35,46) and [135,146) therefore receive both fragments at profile position 5
    //
    // We compare two logically equivalent released workflows:
    // 1. original paired BAM -> `midpoints --gc-file <real pkg>`
    // 2. original paired BAM -> `bam-to-frag --gc-file <real pkg>` ->
    //    `frag-to-bam` -> `midpoints --gc-tag GC`
    //
    // Both workflows must therefore produce the same profile:
    // - groupA: 5.0 at position 5
    // - groupC: 5/9 at position 5
    let reference = twobit_from_sequences(
        "frag_to_bam_real_non_neutral_reference",
        vec![(
            "chr1".to_string(),
            format!("{}N{}", "A".repeat(100), "C".repeat(100)),
        )],
    )?;
    // The producer BAM stacks nine identical C-only fragments at one start. Give each fragment
    // a distinct qname here so the real GC package is built from ten molecules, not one pair
    // that gets aliased by the fixture writer.
    let producer_bam = bam_from_specs_strict_identity(
        vec![("chr1".to_string(), 201)],
        {
            let mut fragments = vec![paired_fragment(10, 61, 20)];
            for _ in 0..9 {
                fragments.push(paired_fragment(110, 61, 20));
            }
            fragments
        },
        Vec::new(),
        "frag_to_bam_real_non_neutral_producer",
    )?;
    let source_bam = bam_from_specs(
        vec![("chr1".to_string(), 201)],
        vec![paired_fragment(10, 61, 20), paired_fragment(110, 61, 20)],
        Vec::new(),
        "frag_to_bam_real_non_neutral_source",
    )?;
    let work = TempDir::new()?;
    let bam_to_frag_out = TempDir::new()?;
    let frag_to_bam_out = TempDir::new()?;
    let original_midpoints_out = TempDir::new()?;
    let restored_midpoints_out = TempDir::new()?;
    let gc_path = build_real_non_neutral_gc_package(
        &producer_bam.bam,
        &reference.path,
        work.path(),
        61,
        "chr1\t0\t100\nchr1\t101\t201\n",
        // Chromosome length 201 and fragment length 61 give 141 valid starts in total. Under the
        // `ref-gc-bias` fit rule these BED rows count exactly 40 pure-A starts and 40 pure-C
        // starts.
        141,
    )?;
    let intervals = work.path().join("sites.bed");
    let chrom_sizes_path = frag_to_bam_out.path().join("chrom.sizes");
    fs::write(&intervals, "chr1\t35\t46\tgroupA\nchr1\t135\t146\tgroupC\n")?;
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 201)])?;

    let mut bam_to_frag_cfg = BamToFragConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: bam_to_frag_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    bam_to_frag_cfg.set_output_prefix("weighted");
    bam_to_frag_cfg.set_min_mapq(0);
    bam_to_frag_cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path.clone()),
        neutralize_invalid_gc: false,
    });
    bam_to_frag_cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let fragment_lengths = bam_to_frag_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 61;
        fragment_lengths.max_fragment_length = 61;
    }

    let mut original_midpoints_cfg = MidpointsConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: original_midpoints_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        intervals.clone(),
    );
    original_midpoints_cfg.set_output_prefix("origsites");
    original_midpoints_cfg.set_length_bins(vec![61, 62]);
    original_midpoints_cfg.set_smoothing(MidpointSmoothing::None);
    original_midpoints_cfg.set_tile_size(1_000);
    original_midpoints_cfg.set_min_mapq(0);
    original_midpoints_cfg.set_require_proper_pair(false);
    original_midpoints_cfg.set_gc(ApplyGCArgs {
        gc_file: Some(gc_path.clone()),
        gc_tag: None,
        neutralize_invalid_gc: false,
    });
    original_midpoints_cfg.set_ref_2bit(Some(reference.path.clone()));

    // Act 1: create a FRAG file with real `gc_weight` values and restore it back to BAM.
    run_bam_to_frag(&bam_to_frag_cfg)?;
    let frag_path = bam_to_frag_out.path().join("weighted.frag.tsv.gz");
    let mut frag_to_bam_cfg = make_config(
        frag_path,
        frag_to_bam_out.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    frag_to_bam_cfg.set_min_mapq(0);
    {
        let fragment_lengths = frag_to_bam_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 61;
        fragment_lengths.max_fragment_length = 61;
    }
    run_frag_to_bam(&frag_to_bam_cfg)?;
    let restored_bam_path = output_bam_path(frag_to_bam_out.path(), "restored");
    build_bai_for_test_bam(&restored_bam_path)?;

    // Act 2: consume the original BAM via `--gc-file` and the restored BAM via `--gc-tag`.
    run_midpoints(&original_midpoints_cfg)?;

    let mut restored_midpoints_cfg = MidpointsConfig::new(
        IOCArgs {
            bam: restored_bam_path,
            output_dir: restored_midpoints_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
        intervals,
    );
    restored_midpoints_cfg.set_output_prefix("restoredsites");
    restored_midpoints_cfg.set_length_bins(vec![61, 62]);
    restored_midpoints_cfg.set_smoothing(MidpointSmoothing::None);
    restored_midpoints_cfg.set_tile_size(1_000);
    restored_midpoints_cfg.set_min_mapq(0);
    restored_midpoints_cfg.set_require_proper_pair(false);
    restored_midpoints_cfg.unpaired.reads_are_fragments = true;
    restored_midpoints_cfg.set_gc(ApplyGCArgs {
        gc_file: None,
        gc_tag: Some("GC".to_string()),
        neutralize_invalid_gc: false,
    });
    run_midpoints(&restored_midpoints_cfg)?;

    // Assert
    let original_arr: Array3<f32> = read_midpoint_zarr_counts(
        original_midpoints_out
            .path()
            .join("origsites.midpoint_profiles.zarr"),
    )?;
    let restored_arr: Array3<f32> = read_midpoint_zarr_counts(
        restored_midpoints_out
            .path()
            .join("restoredsites.midpoint_profiles.zarr"),
    )?;
    let restored_tags = read_aux_tags(&output_bam_path(frag_to_bam_out.path(), "restored"))?;

    assert_eq!(restored_tags.len(), 2);
    assert_optional_f32_eq(restored_tags[0].gc_weight, Some(5.0), "restored GC A-only");
    assert_optional_f32_eq(
        restored_tags[1].gc_weight,
        Some(5.0_f32 / 9.0_f32),
        "restored GC C-only",
    );
    assert_eq!(original_arr, restored_arr);
    assert_eq!(original_arr.shape(), &[2, 1, 11]);
    let group_to_idx = read_group_index_map(
        &original_midpoints_out
            .path()
            .join("origsites.group_index.tsv"),
    )?;
    let group_a_idx = group_to_idx["groupA"];
    let group_c_idx = group_to_idx["groupC"];
    assert!((original_arr[[group_a_idx, 0, 5]] - 5.0).abs() <= 1e-6);
    assert!((original_arr[[group_c_idx, 0, 5]] - (5.0_f32 / 9.0_f32)).abs() <= 1e-6);
    assert!((original_arr.sum() - (5.0_f32 + 5.0_f32 / 9.0_f32)).abs() <= 1e-6);

    Ok(())
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn given_bam_to_frag_with_real_gc_and_scaling_outputs_when_frag_to_bam_runs_then_companion_header_restores_both_aux_tags()
-> Result<()> {
    // Arrange:
    // `simple_inward_bam()` contains one paired fragment spanning [20, 80), length 60.
    //
    // We generate the FRAG input with the real released producer `bam-to-frag` using:
    // - a constant GC package that assigns weight 3.0 to all 60 bp fragments
    // - a uniform scaling TSV with factor 2.0 across chr1
    //
    // Hand-derived expectations:
    // - `bam-to-frag` writes exactly one FRAG row for [20, 80)
    // - the companion header must advertise both extra columns:
    //     gc_weight, coverage_scaling_weight
    // - `frag-to-bam` auto-detects that companion header and restores one unpaired BAM record
    //   with:
    //     GC   = 3.0
    //     cw  = 2.0
    //     fl absent (because `bam-to-frag` does not emit a `flen` column)
    let source_bam = fixtures::simple_inward_bam()?;
    let reference = simple_reference_twobit()?;
    let bam_to_frag_out = TempDir::new()?;
    let frag_to_bam_out = TempDir::new()?;
    let work = TempDir::new()?;
    let gc_path = work.path().join("constant_gc_pkg.zarr");
    let scaling_path = work.path().join("uniform_scaling.tsv");
    let chrom_sizes_path = frag_to_bam_out.path().join("chrom.sizes");
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 200)])?;
    write_scaling_tsv(&scaling_path, &[("chr1", 0, 200, 2.0)])?;

    let package = GCCorrectionPackage {
        version: GC_CORRECTION_SCHEMA_VERSION,
        end_offset: 0,
        length_edges: vec![60, 61],
        gc_edges: vec![0, 101],
        length_bin_frequencies: array![1.0_f64],
        reference_contig_footprint: twobit_contig_footprint(&reference.path)?,
        correction_matrix: array![[3.0_f64]],
    };
    package.write_zarr(&gc_path)?;

    let mut bam_to_frag_cfg = BamToFragConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: bam_to_frag_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    bam_to_frag_cfg.set_output_prefix("weighted");
    bam_to_frag_cfg.set_min_mapq(0);
    bam_to_frag_cfg.set_coverage_scaling_factors(Some(scaling_path));
    bam_to_frag_cfg.set_gc(ApplyGCArgFileOnly {
        gc_file: Some(gc_path),
        neutralize_invalid_gc: false,
    });
    bam_to_frag_cfg.set_ref_2bit(Some(reference.path.clone()));
    {
        let fragment_lengths = bam_to_frag_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 60;
        fragment_lengths.max_fragment_length = 60;
    }

    // Act 1: real producer writes FRAG row plus companion header.
    run_bam_to_frag(&bam_to_frag_cfg)?;
    let frag_path = bam_to_frag_out.path().join("weighted.frag.tsv.gz");
    let header_path = bam_to_frag_out.path().join("weighted.frag.header.tsv");
    let header_text = fs::read_to_string(&header_path)?;

    let mut frag_to_bam_cfg = make_config(
        frag_path,
        frag_to_bam_out.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    frag_to_bam_cfg.set_min_mapq(0);
    {
        let fragment_lengths = frag_to_bam_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 60;
        fragment_lengths.max_fragment_length = 60;
    }

    // Act 2: companion header is auto-detected and restored to BAM AUX tags.
    run_frag_to_bam(&frag_to_bam_cfg)?;
    let output_bam = output_bam_path(frag_to_bam_out.path(), "restored");
    let rows = read_bam_rows(&output_bam)?;
    let aux_tags = read_aux_tags(&output_bam)?;

    // Assert
    assert_eq!(
        header_text.trim(),
        "chromosome\tstart\tend\tmin_mapq\tread1_strand\tgc_weight\tcoverage_scaling_weight"
    );
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 20, 80, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, Some(3.0), "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, Some(2.0), "cw");
    assert_eq!(aux_tags[0].fragment_length, None);

    Ok(())
}

#[cfg(feature = "cmd_bam_to_frag")]
#[test]
fn given_bam_to_frag_with_count_scaling_output_when_frag_to_bam_runs_then_companion_header_restores_cnt_aux_tag()
-> Result<()> {
    // Arrange:
    // `simple_inward_bam()` contains one paired fragment spanning [20, 80), length 60.
    //
    // We generate the FRAG input with the released producer `bam-to-frag` using only a
    // chromosome-wide count-scaling TSV with factor 0.5.
    //
    // Hand-derived expectations:
    // - `bam-to-frag` writes one FRAG row for [20, 80)
    // - the companion header advertises `count_scaling_weight`
    // - `frag-to-bam` auto-detects that companion header and restores one unpaired BAM record
    //   with nw = 0.5
    // - GC, cw, and fl stay absent
    let source_bam = fixtures::simple_inward_bam()?;
    let bam_to_frag_out = TempDir::new()?;
    let frag_to_bam_out = TempDir::new()?;
    let work = TempDir::new()?;
    let scaling_path = work.path().join("uniform_count_scaling.tsv");
    let chrom_sizes_path = frag_to_bam_out.path().join("chrom.sizes");
    write_chrom_sizes(&chrom_sizes_path, &[("chr1", 200)])?;
    write_scaling_tsv(&scaling_path, &[("chr1", 0, 200, 0.5)])?;

    let mut bam_to_frag_cfg = BamToFragConfig::new(
        IOCArgs {
            bam: source_bam.bam.clone(),
            output_dir: bam_to_frag_out.path().to_path_buf(),
            n_threads: 1,
        },
        base_chromosomes(&["chr1"]),
    );
    bam_to_frag_cfg.set_output_prefix("weighted");
    bam_to_frag_cfg.set_min_mapq(0);
    bam_to_frag_cfg.set_count_scaling_factors(Some(scaling_path));
    {
        let fragment_lengths = bam_to_frag_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 60;
        fragment_lengths.max_fragment_length = 60;
    }

    // Act 1: real producer writes FRAG row plus companion header.
    run_bam_to_frag(&bam_to_frag_cfg)?;
    let frag_path = bam_to_frag_out.path().join("weighted.frag.tsv.gz");
    let header_path = bam_to_frag_out.path().join("weighted.frag.header.tsv");
    let header_text = fs::read_to_string(&header_path)?;

    let mut frag_to_bam_cfg = make_config(
        frag_path,
        frag_to_bam_out.path().to_path_buf(),
        chrom_sizes_path,
        base_chromosomes(&["chr1"]),
    );
    frag_to_bam_cfg.set_min_mapq(0);
    {
        let fragment_lengths = frag_to_bam_cfg.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 60;
        fragment_lengths.max_fragment_length = 60;
    }

    // Act 2: companion header is auto-detected and restored to BAM AUX tags.
    run_frag_to_bam(&frag_to_bam_cfg)?;
    let output_bam = output_bam_path(frag_to_bam_out.path(), "restored");
    let rows = read_bam_rows(&output_bam)?;
    let aux_tags = read_aux_tags(&output_bam)?;

    // Assert
    assert_eq!(
        header_text.trim(),
        "chromosome\tstart\tend\tmin_mapq\tread1_strand\tcount_scaling_weight"
    );
    assert_eq!(rows.len(), 1);
    assert_unpaired_full_match_record(&rows[0], "chr1", 20, 80, 60, '+', "fragment_1");
    assert_eq!(aux_tags.len(), 1);
    assert_optional_f32_eq(aux_tags[0].gc_weight, None, "GC");
    assert_optional_f32_eq(aux_tags[0].coverage_scaling_weight, None, "cw");
    assert_optional_f32_eq(aux_tags[0].count_scaling_weight, Some(0.5), "nw");
    assert_eq!(aux_tags[0].fragment_length, None);

    Ok(())
}
