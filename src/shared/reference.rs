use anyhow::{Context, Result};
use fxhash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::io::BufRead;
use std::{
    fs::{self, File},
    io::BufReader,
    ops::RangeBounds,
    path::{Path, PathBuf},
};
use twobit::{TwoBitFile, TwoBitPhysicalFile};

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct ContigFootprintEntry {
    pub name: String,
    pub size: u64,
}

/// Reusable reader for sequence queries against one 2bit file.
pub(crate) struct ReferenceReader {
    reader: TwoBitPhysicalFile,
}

impl ReferenceReader {
    /// Open a 2bit file and load its sequence index.
    pub(crate) fn open(path: &Path) -> anyhow::Result<Self> {
        let reader = TwoBitFile::open(path)
            .with_context(|| format!("opening 2bit reference {}", path.display()))?;
        Ok(Self { reader })
    }

    /// Load the full sequence for one chromosome.
    pub(crate) fn read_seq(&mut self, chr: &str) -> anyhow::Result<Vec<u8>> {
        let seq = self
            .reader
            .read_sequence(chr, ..)
            .with_context(|| format!("extracting reference seq for {chr}"))?;
        Ok(seq.into_bytes())
    }

    /// Load a sequence range for one chromosome.
    pub(crate) fn read_seq_in_range<R>(&mut self, chr: &str, range: R) -> anyhow::Result<Vec<u8>>
    where
        R: RangeBounds<usize> + Clone,
    {
        let seq = self
            .reader
            .read_sequence(chr, range.clone())
            .with_context(|| {
                format!(
                    "extracting reference seq for {}:{:?}-{:?}",
                    chr,
                    range.start_bound().cloned(),
                    range.end_bound().cloned()
                )
            })?;
        Ok(seq.into_bytes())
    }
}

/// Copy a reference into a command's unique working directory.
#[allow(
    dead_code,
    reason = "single-command feature builds may compile reference helpers without staging a 2bit file"
)]
pub(crate) fn stage_reference_2bit(source: &Path, work_dir: &Path) -> anyhow::Result<PathBuf> {
    let destination = work_dir.join("reference.2bit");
    fs::copy(source, &destination).with_context(|| {
        format!(
            "copying 2bit reference {} to {}",
            source.display(),
            destination.display()
        )
    })?;
    Ok(destination)
}

/// Load reference genome sequence for
/// a single chromosome from a 2bit file.
pub fn read_seq<P: AsRef<Path>>(path: P, chr: &str) -> anyhow::Result<Vec<u8>> {
    let mut reader = ReferenceReader::open(path.as_ref())?;
    reader.read_seq(chr)
}

/// Load reference genome sequence for a range of positions
/// in a single chromosome from a 2bit file.
pub fn read_seq_in_range<R, P: AsRef<Path>>(path: P, chr: &str, range: R) -> anyhow::Result<Vec<u8>>
where
    R: RangeBounds<usize> + Clone,
{
    let mut reader = ReferenceReader::open(path.as_ref())?;
    reader.read_seq_in_range(chr, range)
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

/// Return contig names from a .2bit file in reference order.
pub fn twobit_contig_names<P: AsRef<Path>>(path: P) -> Result<Vec<String>> {
    let path = path.as_ref();
    let tb =
        TwoBitFile::open(path).with_context(|| format!("opening 2bit reference {:?}", path))?;
    Ok(tb.chrom_names())
}

/// Return a stable footprint from 2bit contig names and lengths.
///
/// GC correction packages depend on the reference contig set used when they were built.
/// This records the exact contig names and sizes so downstream commands can warn when a package
/// is applied with a different `--ref-2bit`.
///
/// The footprint intentionally excludes file paths and sequence content. It sorts contigs by
/// `(name, size)` so the value does not depend on the order stored in the 2bit header.
pub fn twobit_contig_footprint<P: AsRef<Path>>(path: P) -> Result<Vec<ContigFootprintEntry>> {
    let tb = TwoBitFile::open(path)?;
    let mut entries: Vec<ContigFootprintEntry> = tb
        .chrom_names()
        .into_iter()
        .zip(tb.chrom_sizes())
        .map(|(name, size)| ContigFootprintEntry {
            name,
            size: size as u64,
        })
        .collect();
    entries.sort_unstable_by(|left, right| {
        left.name.cmp(&right.name).then(left.size.cmp(&right.size))
    });
    Ok(entries)
}

/// Load chromosome sizes from a two-column sizes file or .fai.
///
/// Parameters
/// ----------
/// - path:
///   Path to sizes or FAI.
///
/// Returns
/// -------
/// - sizes:
///   Map of chrom -> size (u32, saturating if > u32::MAX).
pub fn load_chrom_sizes<P: AsRef<Path>>(path: P) -> Result<FxHashMap<String, u32>> {
    let (_, sizes) = load_chrom_sizes_with_order(path)?;
    Ok(sizes)
}

/// Load chromosome sizes *in order* from a two-column sizes file or .fai.
///
/// Parameters
/// ----------
/// - path:
///   Path to sizes or FAI.
///
/// Returns
/// -------
/// - sizes:
///   Map of chrom -> size (u32, saturating if > u32::MAX).
pub fn load_chrom_sizes_with_order<P: AsRef<std::path::Path>>(
    path: P,
) -> Result<(Vec<String>, FxHashMap<String, u32>)> {
    let path = path.as_ref();
    let file = File::open(path).with_context(|| format!("Opening chrom sizes {:?}", path))?;
    let reader = BufReader::with_capacity(1 << 20, file);
    let mut sizes: FxHashMap<String, u32> = FxHashMap::default();
    let mut order: Vec<String> = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // Accept both FAI and two-column TSV
        let parts: Vec<&str> = line.split(['\t', ' ']).collect();
        if parts.len() < 2 {
            continue;
        }
        let name = parts[0].trim();
        if name.is_empty() {
            continue;
        }
        let size: u64 = parts[1]
            .trim()
            .parse()
            .with_context(|| format!("Invalid size for '{}'", name))?;
        if sizes.contains_key(name) {
            anyhow::bail!(
                "Duplicate chromosome '{}' in chrom-sizes file {:?}",
                name,
                path
            );
        }
        order.push(name.to_string());
        sizes.insert(name.to_string(), size.min(u32::MAX as u64) as u32);
    }

    Ok((order, sizes))
}

#[cfg(test)]
mod tests {
    include!("reference_tests.rs");
}
