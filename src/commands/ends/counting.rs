use crate::shared::{
    base::{ZEROISH_F64_TOLERANCE, make_canonical, rev_complement},
    kmers::{kmer_codec::KmerSpec, motifs_file::EncodedMotifKey},
};
use anyhow::{Result, ensure};
use fxhash::FxHashMap;
use serde::{Deserialize, Serialize};

/// Sparse motif counts for one output window.
///
/// This keeps the hot path compact by counting encoded motif keys directly. Decoding to user-facing
/// motif strings happens later during reduction and output writing.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct EndMotifCounts {
    pub(crate) counts: FxHashMap<EncodedMotifKey, f64>,
}

/// Sparse motif counts for all windows touched by one tile.
pub(crate) type EndCountsByWindow = FxHashMap<u64, EndMotifCounts>;

/// Sparse selected-target counts for all windows touched by one tile.
///
/// The outer key is the global output row. The inner key is the motifs-file target index assigned
/// during parsing. Post-processing compacts those target indices when `--all-motifs` is not set.
pub(crate) type SelectedEndCountsByWindow = FxHashMap<u64, FxHashMap<u32, f64>>;

impl EndMotifCounts {
    /// Create an empty sparse count map.
    #[inline]
    #[allow(dead_code)]
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Add a single already-validated weighted motif observation to this window.
    ///
    /// Callers are responsible for checking `EndMotifCounts::should_store_weight` before calling
    /// this helper.
    ///
    /// Parameters
    /// ----------
    /// - `key`:
    ///   Encoded motif identity to increment
    /// - `weight`:
    ///   Weight to add for this observation
    ///
    /// Returns
    /// -------
    /// - `()`:
    ///   The map is updated in place
    #[inline]
    pub(crate) fn incr_weighted(&mut self, key: EncodedMotifKey, weight: f64) {
        *self.counts.entry(key).or_insert(0.0) += weight;
    }

    /// Return whether a weight should create a sparse count entry.
    #[inline]
    pub(crate) fn should_store_weight(weight: f64) -> Result<bool> {
        ensure!(
            weight.is_finite(),
            "sparse end motif weight {weight} is not finite"
        );
        ensure!(
            weight >= -ZEROISH_F64_TOLERANCE,
            "sparse end motif weight {weight} is negative, this is not currently supported"
        );
        Ok(weight > ZEROISH_F64_TOLERANCE)
    }
}

/// Decode one sparse end-motif window into final motif labels.
///
/// Decoding happens on the full joined motif, then optional complement collapsing is applied.
/// Motifs containing `N` are dropped to match the existing k-mer output conventions.
///
/// Parameters
/// ----------
/// - `counts`:
///   Sparse encoded counts for one window
/// - `inside_spec`:
///   Codec spec for the inside half, or `None` when `k_inside = 0`
/// - `outside_spec`:
///   Codec spec for the outside half, or `None` when `k_outside = 0`
/// - `collapse_complement`:
///   Whether reverse-complement-equivalent full motifs should be collapsed
///
/// Returns
/// -------
/// - `FxHashMap<String, f64>`:
///   Final motif labels mapped to their merged counts
pub(crate) fn decode_end_motif_counts(
    counts: &EndMotifCounts,
    inside_spec: Option<&KmerSpec>,
    outside_spec: Option<&KmerSpec>,
    collapse_complement: bool,
) -> FxHashMap<String, f64> {
    let mut decoded = FxHashMap::default();

    for (&key, &value) in &counts.counts {
        // Decode motif such that right-end motifs are reverse-complemented
        // into their 5'->3' orientation
        let full_motif = decode_full_motif(key, inside_spec, outside_spec);
        // Collapse to the same-orientation complement when requested
        // This keeps the outside_inside order contract
        let full_motif = if collapse_complement {
            make_canonical(full_motif, false, false)
        } else {
            full_motif
        };
        if full_motif.contains('N') {
            continue;
        }
        let motif_label = format_end_motif_label(&full_motif, inside_spec, outside_spec);
        *decoded.entry(motif_label).or_insert(0.0) += value;
    }

    decoded
}

/// Format a fully decoded motif as `<outside>_<inside>`.
///
/// The full motif string is expected to be oriented already and ordered as `outside || inside`.
///
/// Parameters
/// ----------
/// - `full_motif`:
///   Fully decoded motif sequence in `outside || inside` order
/// - `inside_spec`:
///   Codec spec for the inside half, or `None` when that half is empty
/// - `outside_spec`:
///   Codec spec for the outside half, or `None` when that half is empty
///
/// Returns
/// -------
/// - `String`:
///   Public motif label in `<outside>_<inside>` form
pub(crate) fn format_end_motif_label(
    full_motif: &str,
    inside_spec: Option<&KmerSpec>,
    outside_spec: Option<&KmerSpec>,
) -> String {
    let inside_len = inside_spec.map_or(0, |spec| spec.k);
    let outside_len = outside_spec.map_or(0, |spec| spec.k);
    assert_eq!(
        full_motif.len(),
        inside_len + outside_len,
        "decoded end motif length did not match configured k-mer sizes"
    );

    let (outside, inside) = full_motif.split_at(outside_len);
    format!("{outside}_{inside}")
}

/// Decode one counted key back into its full motif string.
///
/// The two encoded halves are decoded in storage order first:
/// - left ends: `outside || inside`
/// - right ends: `inside || outside`
///
/// Then the full joined string is reverse-complemented when `reverse_on_decode`
/// is set, so the final motif always runs from the fragment end inward in
/// 5'->3' orientation.
///
/// Parameters
/// ----------
/// - `key`:
///   Encoded motif key to decode
/// - `inside_spec`:
///   Codec spec for the inside half, or `None` when that half is empty
/// - `outside_spec`:
///   Codec spec for the outside half, or `None` when that half is empty
///
/// Returns
/// -------
/// - `String`:
///   Full biological motif sequence before optional complement collapse
pub(crate) fn decode_full_motif(
    key: EncodedMotifKey,
    inside_spec: Option<&KmerSpec>,
    outside_spec: Option<&KmerSpec>,
) -> String {
    let inside = inside_spec.map_or_else(String::new, |spec| spec.decode_kmer(key.inside_code));
    let outside = outside_spec.map_or_else(String::new, |spec| spec.decode_kmer(key.outside_code));

    if key.reverse_on_decode {
        rev_complement(&format!("{inside}{outside}"))
    } else {
        format!("{outside}{inside}")
    }
}

/// Serialized tile entry for one counted motif.
///
/// This is the compact on-disk representation used while merging sparse tile count records.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct TileEndMotifCountEntry {
    pub(crate) inside_code: u64,
    pub(crate) outside_code: u64,
    pub(crate) reverse_on_decode: bool,
    pub(crate) value: f64,
}

impl From<(EncodedMotifKey, f64)> for TileEndMotifCountEntry {
    fn from((key, value): (EncodedMotifKey, f64)) -> Self {
        Self {
            inside_code: key.inside_code,
            outside_code: key.outside_code,
            reverse_on_decode: key.reverse_on_decode,
            value,
        }
    }
}

impl From<&TileEndMotifCountEntry> for EncodedMotifKey {
    fn from(entry: &TileEndMotifCountEntry) -> Self {
        Self {
            inside_code: entry.inside_code,
            outside_code: entry.outside_code,
            reverse_on_decode: entry.reverse_on_decode,
        }
    }
}

/// Serialized sparse counts for one output window in one tile.
///
/// Each tile writes one of these per touched output window so reduction can merge the sparse
/// count records later without reconstructing dense matrices first.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct TileWindowEndCounts {
    pub(crate) original_idx: u64,
    pub(crate) entries: Vec<TileEndMotifCountEntry>,
}

#[cfg(test)]
mod tests {
    include!("counting_tests.rs");
}
