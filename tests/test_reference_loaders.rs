use anyhow::{Result, anyhow};
use cfdnalab::commands::cli_common::{ChromosomeArgs, ContigSource};
use cfdnalab::shared::reference::{
    load_chrom_sizes, read_seq, read_seq_in_range, twobit_contig_lengths, twobit_contig_names,
};
use std::io::{BufWriter, Write};
use tempfile::NamedTempFile;
use twobit::convert::{fasta::FastaReader, to_2bit};

const SAMPLE_FASTA: &str = ">chr1\nACGTACGTNN\n>chr2\nTTAA\n";

fn write_twobit(fasta: &str) -> Result<NamedTempFile> {
    let reader = FastaReader::mem_open(fasta.as_bytes().to_vec())?;
    let file = NamedTempFile::new()?;
    {
        let mut writer = BufWriter::new(file.reopen()?);
        to_2bit(&mut writer, &reader).map_err(|err| anyhow!(err.to_string()))?;
        writer.flush()?;
    }
    Ok(file)
}

#[test]
fn read_seq_reads_full_chromosome() -> Result<()> {
    let twobit = write_twobit(SAMPLE_FASTA)?;
    let seq = read_seq(twobit.path(), "chr1")?;
    assert_eq!(seq, b"ACGTACGTNN");
    Ok(())
}

#[test]
fn read_seq_in_range_reads_slice() -> Result<()> {
    let twobit = write_twobit(SAMPLE_FASTA)?;
    let seq = read_seq_in_range(twobit.path(), "chr1", 2..8)?;
    assert_eq!(seq, b"GTACGT");
    Ok(())
}

#[test]
fn read_seq_roundtrips_full_50bp_sequence_with_terminal_partial_byte() -> Result<()> {
    // Arrange:
    // This sequence is 50 bp long, so the final .2bit packed byte stores only two real bases.
    // twobit previously had a bug that wrote the final base wrongly, so we keep a check of this going forward:
    //   chr1[35..50] = TTTTTCCCCCCCCCC
    // If the writer/reader corrupts the last partial byte, the terminal pure-C 10 bp interval
    // `[40,50)` stops being `CCCCCCCCCC`, and downstream GC tests will observe GC%=90 instead of
    // GC%=100. So this test checks the direct I/O contract with no other command logic involved.
    let expected = format!(
        "{}{}{}{}{}",
        "A".repeat(10),
        "T".repeat(10),
        "C".repeat(5) + &"A".repeat(5),
        "T".repeat(10),
        "C".repeat(10)
    );
    let fasta = format!(">chr1\n{expected}\n");
    let twobit = write_twobit(&fasta)?;

    // Act:
    // Read both the full chromosome and the exact half-open range used by the failing fixture.
    let full = read_seq(twobit.path(), "chr1")?;
    let range = read_seq_in_range(twobit.path(), "chr1", 0..50)?;

    // Assert:
    // Both loader entry points must reproduce the original bytes exactly. We also pin the final
    // 15 bp tail so any future failure immediately shows whether the last packed .2bit byte is
    // being decoded incorrectly.
    assert_eq!(full, expected.as_bytes());
    assert_eq!(range, expected.as_bytes());
    assert_eq!(&full[35..50], b"TTTTTCCCCCCCCCC");
    assert_eq!(&range[35..50], b"TTTTTCCCCCCCCCC");
    Ok(())
}

#[test]
fn read_seq_roundtrips_full_52bp_sequence_when_tail_byte_is_not_partial() -> Result<()> {
    // Arrange:
    // This is the same logical fixture as the 50 bp regression above, but with two unused padding
    // bases appended so the chromosome length becomes divisible by 4. That keeps the important
    // `[40,50)` pure-C interval unchanged while avoiding the known upstream partial-byte bug.
    let expected = format!(
        "{}{}{}{}{}{}",
        "A".repeat(10),
        "T".repeat(10),
        "C".repeat(5) + &"A".repeat(5),
        "T".repeat(10),
        "C".repeat(10),
        "A".repeat(2)
    );
    let fasta = format!(">chr1\n{expected}\n");
    let twobit = write_twobit(&fasta)?;

    // Act
    let full = read_seq(twobit.path(), "chr1")?;
    let range = read_seq_in_range(twobit.path(), "chr1", 0..52)?;

    // Assert
    assert_eq!(full, expected.as_bytes());
    assert_eq!(range, expected.as_bytes());
    assert_eq!(&full[35..50], b"TTTTTCCCCCCCCCC");
    assert_eq!(&range[35..50], b"TTTTTCCCCCCCCCC");
    assert_eq!(&full[50..52], b"AA");
    assert_eq!(&range[50..52], b"AA");
    Ok(())
}

#[test]
fn twobit_contig_lengths_filters_requested_contigs() -> Result<()> {
    let twobit = write_twobit(SAMPLE_FASTA)?;
    let lengths = twobit_contig_lengths(twobit.path(), &["chr1".to_string(), "chr3".to_string()])?;
    assert_eq!(lengths.len(), 1);
    assert_eq!(lengths.get("chr1"), Some(&10));
    assert!(!lengths.contains_key("chr2"));
    assert!(!lengths.contains_key("chr3"));
    Ok(())
}

#[test]
fn twobit_contig_names_preserves_reference_order() -> Result<()> {
    let twobit = write_twobit(">chrB\nACGT\n>chrA\nTTAA\n>chrTiny\nCC\n")?;

    let names = twobit_contig_names(twobit.path())?;

    assert_eq!(names, vec!["chrB", "chrA", "chrTiny"]);
    Ok(())
}

#[test]
fn chromosome_args_all_resolves_from_twobit_source_order() -> Result<()> {
    let twobit = write_twobit(">chrB\nACGT\n>chrA\nTTAA\n")?;
    let args = ChromosomeArgs {
        chromosomes: Some(vec!["all".to_string()]),
        chromosomes_file: None,
    };

    let names = args.resolve_chromosomes(Some(ContigSource::ref_2bit(twobit.path())))?;

    assert_eq!(names, vec!["chrB", "chrA"]);
    Ok(())
}

#[test]
fn load_chrom_sizes_valid_file() -> Result<()> {
    let mut file = NamedTempFile::new()?;
    writeln!(
        file,
        "# comment\n\nchr1\t100\nchr2 200\nchr3\t{}",
        (u32::MAX as u64) + 10
    )?;

    let sizes = load_chrom_sizes(file.path())?;
    assert_eq!(sizes.get("chr1"), Some(&100));
    assert_eq!(sizes.get("chr2"), Some(&200));
    assert_eq!(sizes.get("chr3"), Some(&u32::MAX));

    Ok(())
}

#[test]
fn load_chrom_sizes_invalid_file() -> Result<()> {
    let mut file = NamedTempFile::new()?;
    writeln!(file, "chr1\tnot_a_number")?;

    let err = load_chrom_sizes(file.path()).expect_err("expected parse failure");
    assert!(err.to_string().contains("Invalid size"));
    Ok(())
}
