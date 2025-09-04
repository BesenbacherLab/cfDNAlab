use anyhow::Context;
use std::path::Path;
use twobit::TwoBitFile;

/// Load reference genome sequence for
/// a single chromosome from a 2bit file.
pub fn read_seq(path: &Path, chr: &str) -> anyhow::Result<Vec<u8>> {
    // Open 2bit file
    let mut tb = TwoBitFile::open(path).context("opening 2bit")?;
    // Extract reference sequence
    let seq = tb
        .read_sequence(chr, ..)
        .context(format!("extracting reference seq for {}", chr))?;
    Ok(seq.as_bytes().to_vec())
}
