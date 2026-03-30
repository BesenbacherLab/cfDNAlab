use crate::shared::{
    base::{make_canonical, rev_complement},
    kmers::kmer_codec::KmerSpec,
};
use anyhow::Result;
use fxhash::FxHashMap;
use serde::{Deserialize, Serialize};

/// Encoded key for one full end motif before final decoding.
///
/// The two halves use the shared radix-5 codec but are counted together because complement
/// collapsing and final orientation are defined on the full joined motif, not on each half
/// independently.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct EncodedEndMotifKey {
    pub inside_code: u64,
    pub outside_code: u64,
    pub reverse_on_decode: bool,
}

/// Sparse motif counts for one output window.
///
/// This keeps the hot path compact by counting encoded motif keys directly. Decoding to user-facing
/// motif strings happens later during reduction and output writing.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct EndMotifCounts {
    pub counts: FxHashMap<EncodedEndMotifKey, f64>,
}

/// Sparse motif counts for all windows touched by one tile.
pub type EndCountsByWindow = FxHashMap<u64, EndMotifCounts>;

impl EndMotifCounts {
    /// Create an empty sparse end-motif counter.
    ///
    /// This is the standard starting point for one window or one temporary merge target.
    ///
    /// Returns
    /// -------
    /// - `Self`:
    ///   An empty count map
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Create another empty counter with the same logical shape.
    ///
    /// This exists for consistency with other reducers in the codebase where a "zeroed like"
    /// helper is convenient during accumulation.
    ///
    /// Parameters
    /// ----------
    /// - `self`:
    ///   Existing sparse counts used only as a shape hint
    ///
    /// Returns
    /// -------
    /// - `Self`:
    ///   A new empty count map
    #[inline]
    pub fn zeroed_like(&self) -> Self {
        Self::new()
    }

    /// Add one weighted motif observation to this window.
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
    pub fn incr_weighted(&mut self, key: EncodedEndMotifKey, weight: f64) {
        *self.counts.entry(key).or_insert(0.0) += weight;
    }

    /// Merge another sparse window counter into this one.
    ///
    /// Parameters
    /// ----------
    /// - `other`:
    ///   Counts to add into `self`
    ///
    /// Returns
    /// -------
    /// - `Result<()>`:
    ///   `Ok(())` after all weights have been added
    pub fn merge_from(&mut self, other: &Self) -> Result<()> {
        for (key, value) in &other.counts {
            *self.counts.entry(*key).or_insert(0.0) += *value;
        }
        Ok(())
    }

    /// Collapse several sparse window counters into one merged result.
    ///
    /// Parameters
    /// ----------
    /// - `iter`:
    ///   Counts to merge
    ///
    /// Returns
    /// -------
    /// - `Result<Self>`:
    ///   The merged sparse count map
    pub fn collapse<'a, I>(iter: I) -> Result<Self>
    where
        I: IntoIterator<Item = &'a EndMotifCounts>,
    {
        let mut merged = EndMotifCounts::new();
        for counts in iter {
            merged.merge_from(counts)?;
        }
        Ok(merged)
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
pub fn decode_end_motif_counts(
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
pub fn format_end_motif_label(
    full_motif: &str,
    inside_spec: Option<&KmerSpec>,
    outside_spec: Option<&KmerSpec>,
) -> String {
    let inside_len = inside_spec.map_or(0, |spec| spec.k);
    let outside_len = outside_spec.map_or(0, |spec| spec.k);
    debug_assert_eq!(full_motif.len(), inside_len + outside_len);

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
pub fn decode_full_motif(
    key: EncodedEndMotifKey,
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
/// This is the compact on-disk representation used while merging sparse tile payloads.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TileEndMotifCountEntry {
    pub inside_code: u64,
    pub outside_code: u64,
    pub reverse_on_decode: bool,
    pub value: f64,
}

impl From<(EncodedEndMotifKey, f64)> for TileEndMotifCountEntry {
    fn from((key, value): (EncodedEndMotifKey, f64)) -> Self {
        Self {
            inside_code: key.inside_code,
            outside_code: key.outside_code,
            reverse_on_decode: key.reverse_on_decode,
            value,
        }
    }
}

impl From<&TileEndMotifCountEntry> for EncodedEndMotifKey {
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
/// payloads later without reconstructing dense matrices first.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TileWindowEndCounts {
    pub original_idx: u64,
    pub entries: Vec<TileEndMotifCountEntry>,
}

#[cfg(test)]
mod tests {
    include!("counting_tests.rs");
}
