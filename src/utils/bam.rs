use anyhow::{Context, Result, anyhow};
use rust_htslib::bam::{IndexedReader, Read};
use std::path::Path;

/// Create a BAM file reader for a given chromosome.
///
/// Returns Reader, tid, and chromosome length.
pub fn create_chromosome_reader(bam_path: &Path, chr: &str) -> Result<(IndexedReader, u32, u64)> {
    let reader = IndexedReader::from_path(bam_path).context(format!("opening BAM for {}", chr))?;
    let header = reader.header().to_owned();
    let tid = header
        .tid(chr.as_bytes())
        .ok_or_else(|| anyhow!("{} not in BAM", chr))?;
    let chrom_len = header
        .target_len(tid)
        .ok_or_else(|| anyhow!("No length for {}", chr))? as u64;
    Ok((reader, tid, chrom_len))
}
