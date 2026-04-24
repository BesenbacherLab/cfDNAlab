use anyhow::{Context, Result};
use fxhash::{FxHashMap, FxHashSet};
use std::io::BufRead;
use std::{fs::File, io::BufReader, ops::RangeBounds, path::Path};
use twobit::TwoBitFile;

// FNV-1a constants for a compact, deterministic contig-set fingerprint. This is not a
// cryptographic hash. It is only a low-cost mismatch signal for GC correction packages.
const CONTIG_SIGNATURE_OFFSET_A: u64 = 0xcbf29ce484222325;
const CONTIG_SIGNATURE_OFFSET_B: u64 = 0x84222325cbf29ce4;
const CONTIG_SIGNATURE_PRIME: u64 = 0x100000001b3;

/// Load reference genome sequence for
/// a single chromosome from a 2bit file.
pub fn read_seq<P: AsRef<Path>>(path: P, chr: &str) -> anyhow::Result<Vec<u8>> {
    // Open 2bit file
    let mut tb = TwoBitFile::open(path).context("opening 2bit")?;
    // Extract reference sequence
    let seq = tb
        .read_sequence(chr, ..)
        .context(format!("extracting reference seq for {}", chr))?;
    Ok(seq.as_bytes().to_vec())
}

/// Load reference genome sequence for a range of positions
/// in a single chromosome from a 2bit file.
pub fn read_seq_in_range<R, P: AsRef<Path>>(path: P, chr: &str, range: R) -> anyhow::Result<Vec<u8>>
where
    R: RangeBounds<usize> + Clone,
{
    // Open 2bit file
    let mut tb = TwoBitFile::open(path).context("opening 2bit")?;
    // Extract reference sequence
    let seq = tb.read_sequence(chr, range.clone()).context(format!(
        "extracting reference seq for {}:{:?}-{:?}",
        chr,
        range.start_bound().cloned(),
        range.end_bound().cloned()
    ))?;
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

/// Compute a stable signature from 2bit contig names and lengths.
///
/// GC correction packages depend on the reference contig set used when they were built.
/// This records a compact fingerprint of that set so downstream commands can warn when a
/// package is applied with a different `--ref-2bit`.
///
/// The signature intentionally excludes file paths and sequence content. It sorts contigs by
/// `(name, length)` so the value does not depend on the order stored in the 2bit header, then
/// feeds each name and little-endian length into `update_contig_signature` with delimiter bytes
/// between fields and records. The delimiters keep adjacent fields from being concatenated into
/// the same byte stream.
///
/// This is a mismatch signal, not a cryptographic identity check. Equal signatures mean the
/// contig names and lengths probably match. Different signatures mean the package and reference
/// should be treated as different contig sets.
pub fn twobit_contig_signature<P: AsRef<Path>>(path: P) -> Result<[u64; 2]> {
    let tb = TwoBitFile::open(path)?;
    let mut entries: Vec<(String, usize)> =
        tb.chrom_names().into_iter().zip(tb.chrom_sizes()).collect();
    entries.sort_unstable_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));

    let mut signature = [CONTIG_SIGNATURE_OFFSET_A, CONTIG_SIGNATURE_OFFSET_B];
    for (name, size) in entries {
        update_contig_signature(&mut signature, name.as_bytes());
        update_contig_signature(&mut signature, &[0]);
        update_contig_signature(&mut signature, &(size as u64).to_le_bytes());
        update_contig_signature(&mut signature, &[0xff]);
    }
    Ok(signature)
}

/// Fold bytes into the two-lane contig-set signature.
///
/// The first lane is standard 64-bit FNV-1a, using the Fowler-Noll-Vo offset basis
/// `0xcbf29ce484222325` and prime `0x100000001b3`. For each byte, FNV-1a XORs the byte into the
/// accumulator and then multiplies by the prime with wrapping `u64` arithmetic.
/// The recurrence and constants are from the IETF FNV draft:
/// https://www.ietf.org/archive/id/draft-eastlake-fnv-35.html
///
/// The second lane uses the same recurrence with a different starting value and a one-bit-left
/// rotated byte. That second lane is not a separate published hash. It is just an independent
/// cheap check so a collision in one 64-bit lane is less likely to hide a changed contig set.
///
/// This helper is intentionally streaming. `twobit_contig_signature` can call it for each field
/// and delimiter without allocating one combined byte buffer for all contigs.
fn update_contig_signature(signature: &mut [u64; 2], bytes: &[u8]) {
    for byte in bytes {
        signature[0] ^= u64::from(*byte);
        signature[0] = signature[0].wrapping_mul(CONTIG_SIGNATURE_PRIME);
        signature[1] ^= u64::from(*byte).rotate_left(1);
        signature[1] = signature[1].wrapping_mul(CONTIG_SIGNATURE_PRIME);
    }
}

/// Load chromosome sizes from a two-column sizes file or .fai.
///
/// Parameters
/// ----------
/// - path:
///     Path to sizes or FAI.
///
/// Returns
/// -------
/// - sizes:
///     Map of chrom -> size (u32, saturating if > u32::MAX).
pub fn load_chrom_sizes<P: AsRef<Path>>(path: P) -> Result<FxHashMap<String, u32>> {
    let (_, sizes) = load_chrom_sizes_with_order(path)?;
    Ok(sizes)
}

/// Load chromosome sizes *in order* from a two-column sizes file or .fai.
///
/// Parameters
/// ----------
/// - path:
///     Path to sizes or FAI.
///
/// Returns
/// -------
/// - sizes:
///     Map of chrom -> size (u32, saturating if > u32::MAX).
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
