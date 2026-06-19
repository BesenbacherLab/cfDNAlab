use anyhow::{Context, Result, bail};
use fxhash::{FxHashMap, FxHashSet};
#[cfg(writes_bam_output)]
use rust_htslib::bam;
#[cfg(reads_indexed_bam)]
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::{Read, Reader};
#[cfg(writes_bam_output)]
use std::ffi::OsString;
use std::path::Path;
#[cfg(writes_bam_output)]
use std::path::PathBuf;
use url::Url;

#[cfg(writes_bam_output)]
use crate::shared::thread_pool::default_thread_count;

/// Create a BAM file reader for a given chromosome.
///
/// Returns Reader, tid, and chromosome length.
#[cfg(reads_indexed_bam)]
pub fn create_chromosome_reader(bam_path: &Path, chr: &str) -> Result<(IndexedReader, u32, u64)> {
    let reader = open_indexed_bam_reader(bam_path).context(format!("opening BAM for {}", chr))?;
    let header = reader.header().to_owned();
    let tid = header
        .tid(chr.as_bytes())
        .ok_or_else(|| anyhow::anyhow!("{} not in BAM", chr))?;
    let chrom_len = header
        .target_len(tid)
        .ok_or_else(|| anyhow::anyhow!("No length for {}", chr))? as u64;
    Ok((reader, tid, chrom_len))
}

pub fn open_bam_reader(bam_path: &Path) -> Result<Reader> {
    match bam_input_url(bam_path)? {
        Some(url) => Reader::from_url(&url).with_context(|| format!("opening BAM URL {}", url)),
        None => Reader::from_path(bam_path)
            .with_context(|| format!("opening BAM {}", bam_path.display())),
    }
}

#[cfg(reads_indexed_bam)]
fn open_indexed_bam_reader(bam_path: &Path) -> Result<IndexedReader> {
    match bam_input_url(bam_path)? {
        Some(url) => IndexedReader::from_url(&url)
            .with_context(|| format!("opening indexed BAM URL {}", url)),
        None => IndexedReader::from_path(bam_path)
            .with_context(|| format!("opening indexed BAM {}", bam_path.display())),
    }
}

#[cfg(writes_bam_output)]
pub(crate) fn bam_bai_path(bam_path: &Path) -> Result<PathBuf> {
    let file_name = bam_path
        .file_name()
        .with_context(|| format!("BAM path has no file name: {}", bam_path.display()))?;
    let mut index_file_name = OsString::from(file_name);
    index_file_name.push(".bai");
    Ok(bam_path.with_file_name(index_file_name))
}

#[cfg(writes_bam_output)]
pub(crate) fn build_bam_bai_index(bam_path: &Path) -> Result<PathBuf> {
    let bai_path = bam_bai_path(bam_path)?;
    let indexing_threads = u32::try_from(default_thread_count()).unwrap_or(u32::MAX);
    // `samtools index sample.bam` conventionally creates `sample.bam.bai`. Passing the path
    // explicitly keeps cfDNAlab's generated BAM outputs predictable instead of depending on HTSlib's
    // default index-file naming.
    //
    // These conversion commands do not expose their own thread count. Use the same default policy as
    // the shared CLI thread option: leave one core free when possible, but always use at least one
    // indexing thread.
    bam::index::build(
        bam_path,
        Some(bai_path.as_path()),
        bam::index::Type::Bai,
        indexing_threads,
    )
    .with_context(|| {
        format!(
            "indexing BAM {} to {} with {} thread(s)",
            bam_path.display(),
            bai_path.display(),
            indexing_threads
        )
    })?;
    if !bai_path.exists() {
        bail!(
            "BAM index build completed but {} was not created",
            bai_path.display()
        );
    }
    Ok(bai_path)
}

fn bam_input_url(bam_path: &Path) -> Result<Option<Url>> {
    let Some(raw_path) = bam_path.to_str() else {
        return Ok(None);
    };
    if !raw_path.contains("://") {
        return Ok(None);
    }

    let url = Url::parse(raw_path).with_context(|| format!("parsing BAM URL {}", raw_path))?;
    match url.scheme() {
        "ftp" | "http" | "https" => Ok(Some(url)),
        scheme => bail!("unsupported BAM URL scheme '{}'", scheme),
    }
}

pub fn bam_header_contigs<P: AsRef<std::path::Path>>(bam_path: P) -> Result<Vec<String>> {
    let reader = open_bam_reader(bam_path.as_ref())?;
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

/// (tid, len) for each requested chromosome from the BAM header (no index/reads needed).
pub fn bam_contigs_info<P: AsRef<Path>>(bam_path: P, chromosomes: &[String]) -> Result<Contigs> {
    let rdr = open_bam_reader(bam_path.as_ref())?;
    let hdr = rdr.header().to_owned(); // HeaderView -> clone needed for tid2name() lifetime

    let want: FxHashSet<&str> = chromosomes.iter().map(|s| s.as_str()).collect();
    let n_targets = hdr.target_count() as i32;

    let mut out: FxHashMap<String, (i32, u32)> =
        FxHashMap::with_capacity_and_hasher(want.len(), Default::default());

    for tid in 0..n_targets {
        let name = std::str::from_utf8(hdr.tid2name(tid as u32))
            .context("non-UTF8 contig name in BAM header")?
            .to_string();
        if !want.contains(name.as_str()) {
            continue;
        }
        let len = hdr
            .target_len(tid as u32)
            .ok_or_else(|| anyhow::anyhow!("missing target_len for {}", name))?;
        out.insert(name, (tid, len as u32));
    }

    // Ensure all requested chromosomes were found
    for chr in chromosomes {
        if !out.contains_key(chr) {
            anyhow::bail!("chromosome '{}' not found in BAM header", chr);
        }
    }

    Ok(Contigs { contigs: out })
}

#[derive(Debug, Clone)]
pub struct Contigs {
    /// Chromosome -> (tid, length)
    pub contigs: FxHashMap<String, (i32, u32)>,
}

#[cfg(test)]
mod tests {
    use super::*;

    include!("bam_tests.rs");
}
