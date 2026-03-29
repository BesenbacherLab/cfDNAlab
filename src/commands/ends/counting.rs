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
    pub within_code: u64,
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
/// - `within_spec`:
///   Codec spec for the within half, or `None` when `k_within = 0`
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
    within_spec: Option<&KmerSpec>,
    outside_spec: Option<&KmerSpec>,
    collapse_complement: bool,
) -> FxHashMap<String, f64> {
    let mut decoded = FxHashMap::default();

    for (&key, &value) in &counts.counts {
        let full_motif = maybe_collapse_full_motif(
            decode_full_motif(key, within_spec, outside_spec),
            collapse_complement,
        );
        if full_motif.contains('N') {
            continue;
        }
        let motif_label = format_end_motif_label(&full_motif, within_spec, outside_spec);
        *decoded.entry(motif_label).or_insert(0.0) += value;
    }

    decoded
}

/// Format a fully decoded motif as `<outside>_<within>`.
///
/// The full motif string is expected to be oriented already and ordered as `outside || within`.
///
/// Parameters
/// ----------
/// - `full_motif`:
///   Fully decoded motif sequence in `outside || within` order
/// - `within_spec`:
///   Codec spec for the within half, or `None` when that half is empty
/// - `outside_spec`:
///   Codec spec for the outside half, or `None` when that half is empty
///
/// Returns
/// -------
/// - `String`:
///   Public motif label in `<outside>_<within>` form
pub fn format_end_motif_label(
    full_motif: &str,
    within_spec: Option<&KmerSpec>,
    outside_spec: Option<&KmerSpec>,
) -> String {
    let within_len = within_spec.map_or(0, |spec| spec.k);
    let outside_len = outside_spec.map_or(0, |spec| spec.k);
    debug_assert_eq!(full_motif.len(), within_len + outside_len);

    let (outside, within) = full_motif.split_at(outside_len);
    format!("{outside}_{within}")
}

/// Decode one counted key back into its full motif string.
///
/// The two encoded halves are decoded in storage order first:
/// - left ends: `outside || within`
/// - right ends: `within || outside`
///
/// Then the full joined string is reverse-complemented when `reverse_on_decode`
/// is set, so the final motif always runs from the fragment end inward in
/// 5'->3' orientation.
///
/// Parameters
/// ----------
/// - `key`:
///   Encoded motif key to decode
/// - `within_spec`:
///   Codec spec for the within half, or `None` when that half is empty
/// - `outside_spec`:
///   Codec spec for the outside half, or `None` when that half is empty
///
/// Returns
/// -------
/// - `String`:
///   Full biological motif sequence before optional complement collapse
pub fn decode_full_motif(
    key: EncodedEndMotifKey,
    within_spec: Option<&KmerSpec>,
    outside_spec: Option<&KmerSpec>,
) -> String {
    let within = within_spec.map_or_else(String::new, |spec| spec.decode_kmer(key.within_code));
    let outside = outside_spec.map_or_else(String::new, |spec| spec.decode_kmer(key.outside_code));

    let storage_order = if key.reverse_on_decode {
        format!("{within}{outside}")
    } else {
        format!("{outside}{within}")
    };

    if key.reverse_on_decode {
        rev_complement(&storage_order)
    } else {
        storage_order
    }
}

/// Apply optional reverse-complement collapsing to a decoded motif string.
///
/// Parameters
/// ----------
/// - `motif`:
///   Fully decoded motif sequence
/// - `collapse_complement`:
///   Whether to canonicalize the motif with its reverse complement
///
/// Returns
/// -------
/// - `String`:
///   The original motif or its canonical representative
#[inline]
pub fn maybe_collapse_full_motif(motif: String, collapse_complement: bool) -> String {
    if collapse_complement {
        make_canonical(motif)
    } else {
        motif
    }
}

/// Serialized tile entry for one counted motif.
///
/// This is the compact on-disk representation used while merging sparse tile payloads.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TileEndMotifCountEntry {
    pub within_code: u64,
    pub outside_code: u64,
    pub reverse_on_decode: bool,
    pub value: f64,
}

impl From<(EncodedEndMotifKey, f64)> for TileEndMotifCountEntry {
    fn from((key, value): (EncodedEndMotifKey, f64)) -> Self {
        Self {
            within_code: key.within_code,
            outside_code: key.outside_code,
            reverse_on_decode: key.reverse_on_decode,
            value,
        }
    }
}

impl From<&TileEndMotifCountEntry> for EncodedEndMotifKey {
    fn from(entry: &TileEndMotifCountEntry) -> Self {
        Self {
            within_code: entry.within_code,
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
    use super::*;
    use crate::shared::kmers::kmer_codec::build_kmer_specs;

    fn spec_for_k(k: u8) -> KmerSpec {
        let specs = build_kmer_specs(&[k]).expect("valid k-mer spec");
        specs[&k].clone()
    }

    #[test]
    fn decode_full_motif_keeps_left_end_in_storage_order() {
        // Arrange: left-end storage order is outside || within.
        let within_spec = spec_for_k(2);
        let outside_spec = spec_for_k(2);
        let key = EncodedEndMotifKey {
            within_code: within_spec.encode_kmer_bytes(b"GT"),
            outside_code: outside_spec.encode_kmer_bytes(b"AC"),
            reverse_on_decode: false,
        };

        // Act
        let motif = decode_full_motif(key, Some(&within_spec), Some(&outside_spec));

        // Assert
        assert_eq!(motif, "ACGT");
    }

    #[test]
    fn decode_full_motif_reverse_complements_right_end() {
        // Arrange: right-end storage order is within || outside, then reverse-complemented.
        // "AAAC" reverse-complements to "GTTT".
        let within_spec = spec_for_k(2);
        let outside_spec = spec_for_k(2);
        let key = EncodedEndMotifKey {
            within_code: within_spec.encode_kmer_bytes(b"AA"),
            outside_code: outside_spec.encode_kmer_bytes(b"AC"),
            reverse_on_decode: true,
        };

        // Act
        let motif = decode_full_motif(key, Some(&within_spec), Some(&outside_spec));

        // Assert
        assert_eq!(motif, "GTTT");
    }

    #[test]
    fn decode_end_motif_counts_collapses_reverse_complements_when_requested() {
        // Arrange: "GT" and "AC" are reverse complements and both canonicalize to "AC".
        let within_spec = spec_for_k(2);
        let mut counts = EndMotifCounts::new();
        counts.incr_weighted(
            EncodedEndMotifKey {
                within_code: within_spec.encode_kmer_bytes(b"GT"),
                outside_code: 0,
                reverse_on_decode: false,
            },
            1.5,
        );
        counts.incr_weighted(
            EncodedEndMotifKey {
                within_code: within_spec.encode_kmer_bytes(b"AC"),
                outside_code: 0,
                reverse_on_decode: false,
            },
            2.0,
        );

        // Act
        let decoded = decode_end_motif_counts(&counts, Some(&within_spec), None, true);

        // Assert
        assert_eq!(decoded.len(), 1);
        assert_eq!(decoded.get("_AC"), Some(&3.5));
    }

    #[test]
    fn decode_end_motif_counts_drops_motifs_with_n() {
        // Arrange: sentinel-N decodes to an N-containing motif and should be dropped.
        let within_spec = spec_for_k(2);
        let mut counts = EndMotifCounts::new();
        counts.incr_weighted(
            EncodedEndMotifKey {
                within_code: within_spec.encode_kmer_bytes(b"AN"),
                outside_code: 0,
                reverse_on_decode: false,
            },
            1.0,
        );

        // Act
        let decoded = decode_end_motif_counts(&counts, Some(&within_spec), None, false);

        // Assert
        assert!(decoded.is_empty());
    }

    #[test]
    fn format_end_motif_label_formats_full_motif_as_outside_within() {
        // Arrange / Act
        let within_spec = spec_for_k(2);
        let outside_spec = spec_for_k(2);
        let label = format_end_motif_label("ACGT", Some(&within_spec), Some(&outside_spec));

        // Assert
        assert_eq!(label, "AC_GT");
    }
}
