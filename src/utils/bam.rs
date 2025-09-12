use anyhow::{Context, Result, anyhow};
use rust_htslib::bam::{IndexedReader, Read, Reader};
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

pub fn bam_header_contigs<P: AsRef<std::path::Path>>(bam_path: P) -> Result<Vec<String>> {
    let reader = Reader::from_path(bam_path)?;
    let header = reader.header();
    let names = header
        .target_names()
        .iter()
        .map(|b| {
            std::str::from_utf8(b)
                .context("non-UTF8 contig name in BAM header")
                .map(|s| s.to_owned())
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(names)
}

// pub fn bam_header_contigs_with_len<P: AsRef<std::path::Path>>(
//     bam_path: P,
// ) -> Result<Vec<(String, u64)>> {
//     let reader = Reader::from_path(bam_path)?;
//     let header = reader.header();

//     let mut result = Vec::new();
//     for (tid, name_bytes) in header.target_names().iter().enumerate() {
//         let name = std::str::from_utf8(name_bytes)
//             .context("non-UTF8 contig name in BAM header")?
//             .to_owned();
//         let len = header
//             .target_len(tid as u32)
//             .context("missing contig length in BAM header")?;
//         result.push((name, len as u64));
//     }
//     Ok(result)
// }
