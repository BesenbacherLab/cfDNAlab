use crate::shared::{
    base::{ZEROISH_F64_TOLERANCE, make_canonical, rev_complement},
    kmers::kmer_codec::{KmerCodes, KmerSpec, SubspaceKmerSpec, build_left_aligned_codes_for_spec},
};
use anyhow::{Result, ensure};
use fxhash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

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

/// Meaning of the public `motif` axis in end-motif Zarr output.
///
/// The count matrix always uses a numeric column coordinate internally. This enum records whether
/// those numeric columns should be interpreted as concrete motif labels or as motif-file group
/// labels when metadata is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EndMotifColumnKind {
    /// Each column is one concrete `<outside>_<inside>` motif
    Motif,
    /// Each column is one user-defined group from the second motifs-file column
    MotifGroup,
}

/// Precomputed motifs-file lookup used by the selected counting path.
///
/// A motifs file defines two related things:
/// - The output target axis, stored in `labels`
/// - The encoded motif states that should count into each target, stored in `lookup`
///
/// `lookup` is keyed by the same encoded state that tile-local motif counting builds. It includes
/// the `reverse_on_decode` flag, so left-end motifs and reverse-complemented right-end states stay
/// distinct when the motifs file maps them to different targets.
#[derive(Debug, Clone)]
pub(crate) struct SelectedEndMotifLookup {
    /// Output labels in motifs-file target order
    pub(crate) labels: Vec<String>,
    /// Whether `labels` are concrete motifs or user-defined motif groups
    pub(crate) column_kind: EndMotifColumnKind,
    /// Codec spec for the inside half, if `k_inside > 0`
    pub(crate) inside_spec: Option<EndMotifHalfSpec>,
    /// Codec spec for the outside half, if `k_outside > 0`
    pub(crate) outside_spec: Option<EndMotifHalfSpec>,
    /// Encoded end-motif key to original target index in `labels`
    pub(crate) lookup: FxHashMap<EncodedEndMotifKey, u32>,
}

impl SelectedEndMotifLookup {
    /// Return the motifs-file target for an encoded end motif.
    ///
    /// This is intentionally a thin map lookup. Counting without a motifs file never calls it, and
    /// motifs-file counting has already paid the parsing and validation cost before tile processing
    /// starts.
    ///
    /// Parameters
    /// ----------
    /// - `key`:
    ///   Encoded motif state observed for one fragment end
    ///
    /// Returns
    /// -------
    /// - `Option<u32>`:
    ///   Original motifs-file target index when the observed state is selected
    pub(crate) fn target_for(&self, key: EncodedEndMotifKey) -> Option<u32> {
        self.lookup.get(&key).copied()
    }
}

/// Sparse selected-target counts for all windows touched by one tile.
///
/// The outer key is the global output row. The inner key is the motifs-file target index assigned
/// during parsing. Post-processing compacts those target indices when `--all-motifs` is not set.
pub(crate) type SelectedEndCountsByWindow = FxHashMap<u64, FxHashMap<u32, f64>>;

/// Codec used for one inside or outside motif half during tile-local counting.
///
/// Full motif output uses the radix-5 [`KmerSpec`] so sparse encoded keys can be decoded into motif
/// strings during reduction. Motifs-file output also uses radix-5 for halves up to the full
/// radix-5 limit. It switches to [`SubspaceKmerSpec`] only for larger halves, because output labels
/// already come from the motifs file and the full motif universe cannot be represented.
#[derive(Clone, Debug)]
pub(crate) enum EndMotifHalfSpec {
    /// Full radix-5 k-mer space
    Radix5(Arc<KmerSpec>),
    /// Byte-backed selected-k-mer subspace for motifs-file halves above the radix-5 limit
    Subspace(Arc<SubspaceKmerSpec>),
}

impl EndMotifHalfSpec {
    /// Wrap a full-space radix-5 spec in the end-motif codec enum.
    pub(crate) fn from_radix5(spec: KmerSpec) -> Self {
        EndMotifHalfSpec::Radix5(Arc::new(spec))
    }

    /// Wrap a byte-backed selected-subspace spec in the end-motif codec enum.
    pub(crate) fn from_subspace(spec: SubspaceKmerSpec) -> Self {
        EndMotifHalfSpec::Subspace(Arc::new(spec))
    }

    /// Wrap a shared byte-backed selected-subspace spec in the end-motif codec enum.
    pub(crate) fn from_shared_subspace(spec: Arc<SubspaceKmerSpec>) -> Self {
        EndMotifHalfSpec::Subspace(spec)
    }

    /// Return the motif-half length.
    #[inline]
    pub(crate) fn k(&self) -> usize {
        match self {
            EndMotifHalfSpec::Radix5(spec) => spec.k,
            EndMotifHalfSpec::Subspace(spec) => spec.k,
        }
    }

    /// Encode one exact motif-half byte slice.
    #[inline]
    pub(crate) fn encode_kmer_bytes(&self, seq: &[u8]) -> u64 {
        match self {
            EndMotifHalfSpec::Radix5(spec) => spec.encode_kmer_bytes(seq),
            EndMotifHalfSpec::Subspace(spec) => spec.encode_kmer_bytes(seq),
        }
    }

    /// Build per-position codes for a tile-local reference slice.
    #[inline]
    pub(crate) fn build_left_aligned_codes(&self, seq: &[u8]) -> KmerCodes {
        match self {
            EndMotifHalfSpec::Radix5(spec) => build_left_aligned_codes_for_spec(seq, spec),
            EndMotifHalfSpec::Subspace(spec) => spec.build_left_aligned_codes(seq),
        }
    }

    /// Return whether an encoded code is invalid for this motif half.
    #[inline]
    pub(crate) fn code_is_invalid(&self, code: u64) -> bool {
        match self {
            EndMotifHalfSpec::Radix5(spec) => {
                code == spec.sentinel_none() || code == spec.sentinel_n()
            }
            EndMotifHalfSpec::Subspace(spec) => code == spec.sentinel_missing(),
        }
    }

    /// Return the invalid code used when a reference-coordinate lookup has no full motif half.
    #[inline]
    pub(crate) fn missing_reference_code(&self) -> u64 {
        match self {
            EndMotifHalfSpec::Radix5(spec) => spec.sentinel_none(),
            EndMotifHalfSpec::Subspace(spec) => spec.sentinel_missing(),
        }
    }

    /// Return the invalid code used when blacklist masking touches a reference motif half.
    #[inline]
    pub(crate) fn masked_reference_code(&self) -> u64 {
        match self {
            EndMotifHalfSpec::Radix5(spec) => spec.sentinel_n(),
            EndMotifHalfSpec::Subspace(spec) => spec.sentinel_missing(),
        }
    }

    /// Return whether two specs can share one precomputed reference-code vector.
    #[inline]
    pub(crate) fn can_share_reference_codes_with(&self, other: &Self) -> bool {
        matches!(
            (self, other),
            (EndMotifHalfSpec::Radix5(left), EndMotifHalfSpec::Radix5(right)) if left.k == right.k
        ) || matches!(
            (self, other),
            (EndMotifHalfSpec::Subspace(left), EndMotifHalfSpec::Subspace(right)) if Arc::ptr_eq(left, right)
        )
    }
}

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
    pub fn incr_weighted(&mut self, key: EncodedEndMotifKey, weight: f64) {
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
