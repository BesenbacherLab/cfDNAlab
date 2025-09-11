use anyhow::{Context, Result};
use fxhash::{FxHashMap, FxHashSet};
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

/// Return (chrom_name, length) for the requested contigs in a .2bit file
pub fn twobit_contig_lengths<P: AsRef<Path>>(
    path: P,
    chromosomes: &[String],
) -> Result<FxHashMap<String, usize>> {
    let tb = TwoBitFile::open(path)?;
    let mut name_to_size: FxHashMap<String, usize> =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());
    let chromosomes_set: FxHashSet<&str> = chromosomes.iter().map(|s| s.as_str()).collect();
    for (name, size) in tb.chrom_names().iter().zip(tb.chrom_sizes()) {
        if chromosomes_set.contains(name.as_str()) {
            name_to_size.insert(name.to_string(), size);
        }
    }
    Ok(name_to_size)
}
