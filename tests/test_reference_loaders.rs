use anyhow::{Result, anyhow};
use cfdnalab::shared::reference::{
    load_chrom_sizes, read_seq, read_seq_in_range, twobit_contig_lengths,
};
use std::io::{BufWriter, Write};
use tempfile::NamedTempFile;
use twobit::convert::{fasta::FastaReader, to_2bit};

const SAMPLE_FASTA: &str = ">chr1\nACGTACGTNN\n>chr2\nTTAA\n";

fn write_twobit(fasta: &str) -> Result<NamedTempFile> {
    let reader = FastaReader::mem_open(fasta.as_bytes().to_vec())?;
    let mut file = NamedTempFile::new()?;
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
