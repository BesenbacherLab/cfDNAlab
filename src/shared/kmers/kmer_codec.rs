use crate::shared::{
    base::{BASES, encode_base, rev_complement},
    positioning::PositionGroup,
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

/// * `k`    – length
/// * `code` – packed reference code in the narrowest type, promoted to u64
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Kmer {
    pub k: u8,
    pub code: u64,
    pub orientation: KmerOrientation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum KmerOrientation {
    Forward,
    Reverse,
}

impl KmerOrientation {
    /// Get KmerOrientation from a PositionGroup
    ///
    /// Left/Mid => forward, right => reverse
    pub fn from_position_group(group: PositionGroup) -> KmerOrientation {
        match group {
            PositionGroup::Left | PositionGroup::Mid => KmerOrientation::Forward,
            PositionGroup::Right => KmerOrientation::Reverse,
        }
    }
}

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
impl Kmer {
    /// User-readable string representation.
    /// Requires a `KmerSpec` table to know how to decode arbitrary k.
    pub fn to_string(&self, specs: &FxHashMap<u8, KmerSpec>) -> String {
        let motif = specs[&self.k].decode_kmer(self.code);
        match self.orientation {
            KmerOrientation::Forward => motif,
            KmerOrientation::Reverse => rev_complement(&motif),
        }
    }
}

/// The narrowest integer width that can accommodate the code space for a k‑mer
/// length, *plus* the two reserved sentinel values.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum Width {
    U8,
    U16,
    U32,
    U64,
}

/// Per-position code vector stored in the tightest possible type.
#[derive(Debug)]
pub enum KmerCodes {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
    U64(Vec<u64>),
}

impl KmerCodes {
    /// Return the code at position `idx` as `u64`.
    #[inline]
    pub fn get(&self, idx: usize) -> u64 {
        match self {
            KmerCodes::U8(v) => v[idx] as u64,
            KmerCodes::U16(v) => v[idx] as u64,
            KmerCodes::U32(v) => v[idx] as u64,
            KmerCodes::U64(v) => v[idx],
        }
    }
}

/// One fully‑specified encoder/decoder for a particular k.
#[derive(Clone, Debug)]
pub struct KmerSpec {
    /// Window length
    pub k: usize,
    /// Integer width used for storage
    width: Width,
    /// Code used when no full k‑mer is available (chromosome ends)
    sentinel_none: u64,
    /// Code used when the window contains any 'N' base
    sentinel_n: u64,
}

impl KmerSpec {
    /// Build per‑position codes for the provided reference sequence.
    pub fn build_left_aligned_codes(&self, seq: &[u8]) -> Vec<u64> {
        build_left_aligned_codes(seq, self.k, self.sentinel_none, self.sentinel_n)
    }

    /// Encode one exact k-mer-sized slice with the shared radix-5 representation.
    ///
    /// Returns:
    /// - `sentinel_none` when `seq.len() != self.k`
    /// - `sentinel_n` when any base encodes as `N`
    /// - otherwise the radix-5 code for the provided bases
    pub fn encode_kmer_bytes(&self, seq: &[u8]) -> u64 {
        if seq.len() != self.k {
            return self.sentinel_none;
        }

        let mut code: u64 = 0;
        for &base in seq {
            let encoded = encode_base(base) as u64;
            if encoded == 4 {
                return self.sentinel_n;
            }
            code = code * 5 + encoded;
        }
        code
    }

    /// Decode a single code back to its k‑mer string, returning all‑'N' if the
    /// code is one of the sentinels.
    pub fn decode_kmer(&self, code: u64) -> String {
        decode_kmer(code, self.k, self.sentinel_none, self.sentinel_n)
    }

    /// Public accessor for the “no full k‑mer” sentinel.
    pub fn sentinel_none(&self) -> u64 {
        self.sentinel_none
    }

    /// Public accessor for the “contains N” sentinel.
    pub fn sentinel_n(&self) -> u64 {
        self.sentinel_n
    }
}

/// K-mer encoder restricted to a selected k-mer subspace.
///
/// This is useful when callers know the only k-mers they care about before counting starts. Codes
/// are compact indices into first-seen unique normalized k-mers. Unselected, malformed,
/// N-containing, and out-of-bounds k-mers map to `sentinel_missing`.
#[derive(Clone, Debug)]
pub struct SubspaceKmerSpec {
    /// Window length
    pub k: usize,
    /// Integer width used for precomputed selected-code arrays
    width: Width,
    /// Sentinel used when a k-mer is not in the selected subspace
    sentinel_missing: u64,
    /// Encoding strategy chosen from k and the selected k-mer set
    encoding: SubspaceKmerEncoding,
}

/// Encoding strategy for selected k-mer subspaces.
///
/// Small enough k-mers reuse the existing radix-5 encoder and remap full-space codes to compact
/// selected codes. Larger k-mers use byte-slice lookup, avoiding the full-space `u64` limit.
#[derive(Clone, Debug)]
enum SubspaceKmerEncoding {
    Radix5 {
        full_spec: KmerSpec,
        selected_codes_by_radix_code: SelectedCodesByRadixCode,
    },
    Bytes {
        selected_codes_by_kmer_bytes: SelectedCodesByKmerBytes,
    },
}

/// Lookup table for radix-backed subspace encoding.
///
/// The key is the ordinary full-space radix-5 code produced by `KmerSpec`. The value is the compact
/// selected-subspace code assigned from the first-seen unique normalized k-mer order. Missing keys
/// mean the k-mer is valid DNA but not part of the selected subspace.
#[derive(Clone, Debug)]
struct SelectedCodesByRadixCode(FxHashMap<u64, u64>);

/// Lookup table for byte-backed subspace encoding.
///
/// The key is the normalized uppercase ACGT k-mer byte string. The value is the compact
/// selected-subspace code assigned from the first-seen unique normalized k-mer order. This variant
/// is used when k is too large for the full-space radix-5 `u64` representation.
#[derive(Clone, Debug)]
struct SelectedCodesByKmerBytes(FxHashMap<Box<[u8]>, u64>);

impl SubspaceKmerSpec {
    /// Return the selected-code sentinel.
    #[inline]
    pub fn sentinel_missing(&self) -> u64 {
        self.sentinel_missing
    }

    /// Return the selected code for one exact k-mer.
    ///
    /// The returned code is either a compact index into the selected k-mer subspace or
    /// `sentinel_missing` when the input is not exactly length `k`, contains a non-ACGT base, or is
    /// not selected.
    #[inline]
    pub fn encode_kmer_bytes(&self, seq: &[u8]) -> u64 {
        match &self.encoding {
            SubspaceKmerEncoding::Radix5 {
                full_spec,
                selected_codes_by_radix_code,
            } => {
                let full_code = full_spec.encode_kmer_bytes(seq);
                selected_codes_by_radix_code
                    .0
                    .get(&full_code)
                    .copied()
                    .unwrap_or(self.sentinel_missing)
            }
            SubspaceKmerEncoding::Bytes {
                selected_codes_by_kmer_bytes,
            } => {
                let Some(normalized) = normalize_acgt_kmer(seq, self.k) else {
                    return self.sentinel_missing;
                };
                selected_codes_by_kmer_bytes
                    .0
                    .get(normalized.as_ref())
                    .copied()
                    .unwrap_or(self.sentinel_missing)
            }
        }
    }

    /// Build left-aligned selected codes for every reference position.
    ///
    /// The result has the same length as `seq`. Positions without a complete selected k-mer contain
    /// `sentinel_missing`. The packed dtype is chosen from the subspace size, not from the full
    /// k-mer universe.
    pub fn build_left_aligned_codes(&self, seq: &[u8]) -> KmerCodes {
        let raw = match &self.encoding {
            SubspaceKmerEncoding::Radix5 {
                full_spec,
                selected_codes_by_radix_code,
            } => full_spec
                .build_left_aligned_codes(seq)
                .into_iter()
                .map(|full_code| {
                    selected_codes_by_radix_code
                        .0
                        .get(&full_code)
                        .copied()
                        .unwrap_or(self.sentinel_missing)
                })
                .collect(),
            SubspaceKmerEncoding::Bytes { .. } => build_left_aligned_subspace_codes(self, seq),
        };

        pack_codes(raw, self.width)
    }
}

/// Construct a `KmerSpec` for each k.
///
/// * Duplicate sizes result in an error.
pub fn build_kmer_specs(kmer_sizes: &[u8]) -> Result<FxHashMap<u8, KmerSpec>> {
    let mut seen = FxHashSet::default();
    let mut specs = FxHashMap::default();

    for &k in kmer_sizes {
        if k < 1 {
            bail!("Illegal k-mer size {k}. Must be positive.");
        }
        // TODO: Calculate actual limit possible!
        if k > 27 {
            bail!("k-mer size {k} is too large. Highest allowed k is 27");
        }
        if !seen.insert(k) {
            bail!("Duplicate k-mer size {k}");
        }
        let (width, sentinel_none, sentinel_n) =
            choose_width(k as usize).context(format!("calculating dtype for k={:?}", k))?;
        specs.insert(
            k,
            KmerSpec {
                k: k as usize,
                width,
                sentinel_none,
                sentinel_n,
            },
        );
    }
    Ok(specs)
}

/// Build a selected-subspace k-mer spec.
///
/// Selected codes are assigned in first-seen order after normalizing to uppercase ACGT and dropping
/// duplicate k-mers. For k values that fit the existing radix-5 representation, this reuses
/// `KmerSpec` internally. Larger k values fall back to byte-slice lookup while keeping the same
/// public selected-code behavior.
///
/// Parameters
/// ----------
/// - `k`:
///   K-mer length for every selected entry
/// - `selected_kmers`:
///   Selected ACGT k-mers in desired code order. Duplicates after normalization are ignored.
///
/// Returns
/// -------
/// - `Result<SubspaceKmerSpec>`:
///   A selected-subspace encoder with compact codes and an adaptive precompute dtype
pub fn build_subspace_kmer_spec<T>(k: usize, selected_kmers: &[T]) -> Result<SubspaceKmerSpec>
where
    T: AsRef<[u8]>,
{
    ensure!(k > 0, "subspace k-mer size must be positive");
    ensure!(
        !selected_kmers.is_empty(),
        "subspace k-mer spec requires at least one selected k-mer"
    );

    let mut normalized_selected_kmers = Vec::with_capacity(selected_kmers.len());
    let mut seen_kmers = FxHashSet::default();
    for (input_index, raw_kmer) in selected_kmers.iter().enumerate() {
        let normalized = normalize_acgt_kmer(raw_kmer.as_ref(), k)
            .with_context(|| format!("invalid selected k-mer at index {input_index}"))?;
        if seen_kmers.insert(normalized.clone()) {
            normalized_selected_kmers.push(normalized);
        }
    }

    let (width, sentinel_missing) = choose_subspace_width(normalized_selected_kmers.len())?;
    let full_spec = build_radix5_spec_if_supported(k)?;
    let encoding = match full_spec {
        Some(full_spec) => {
            let mut selected_codes_by_radix_code = FxHashMap::default();
            for (code_index, normalized) in normalized_selected_kmers.iter().enumerate() {
                let full_code = full_spec.encode_kmer_bytes(normalized.as_ref());
                selected_codes_by_radix_code.insert(full_code, code_index as u64);
            }
            SubspaceKmerEncoding::Radix5 {
                full_spec,
                selected_codes_by_radix_code: SelectedCodesByRadixCode(
                    selected_codes_by_radix_code,
                ),
            }
        }
        None => {
            let mut selected_codes_by_kmer_bytes = FxHashMap::default();
            for (code_index, normalized) in normalized_selected_kmers.into_iter().enumerate() {
                selected_codes_by_kmer_bytes.insert(normalized, code_index as u64);
            }
            SubspaceKmerEncoding::Bytes {
                selected_codes_by_kmer_bytes: SelectedCodesByKmerBytes(
                    selected_codes_by_kmer_bytes,
                ),
            }
        }
    };

    Ok(SubspaceKmerSpec {
        k,
        width,
        sentinel_missing,
        encoding,
    })
}

/// Build one kmer code vector for every `KmerSpec` and store it in a map keyed by `k`.
///
/// The vector is kept in the narrowest width dictated by `spec.width`.
/// This preserves the RAM benefit of the width-selection logic.
///
/// The hash map key is always the `k` value of the corresponding spec.
///
/// Example:
/// ```rust,ignore
/// use pairbase::pairbase::kmer_codec::build_codes_per_k;
/// let codes_by_k = build_codes_per_k(&seq_bytes, kmer_specs);
/// let trinuc_codes = &codes_by_k[&3];
/// let dinuc_codes  = &codes_by_k[&2];
/// ```
pub fn build_left_aligned_codes_per_k(
    seq: &[u8],
    specs: &FxHashMap<u8, KmerSpec>,
) -> FxHashMap<u8, KmerCodes> {
    let mut map = FxHashMap::default();

    for (k, spec) in specs {
        map.insert(
            *k,
            pack_codes(spec.build_left_aligned_codes(seq), spec.width),
        );
    }

    map
}

/* ------------------------------------------------------------------------- */
/*  Internal helpers                                                         */
/* ------------------------------------------------------------------------- */

/// Decide which integer width is sufficient for the code space of this k.
/// The top two codes of the chosen width are reserved as sentinels.
pub fn choose_width(k: usize) -> Result<(Width, u64, u64)> {
    // `u128` is used so that 5^k never overflows during width selection.
    // Even for k = 27 we have 5^k ≈ 7.4e18 < 2^128, so the calculation is safe.
    // The value is then compared to the MAX of each smaller integer type.
    let max_real_code = 5u128.pow(k as u32) - 1; // Highest real code (no sentinels)

    macro_rules! fits_in {
        ($ty:ty) => {
            max_real_code <= (<$ty>::MAX as u128 - 2)
        };
    }

    if fits_in!(u8) {
        Ok((Width::U8, u8::MAX as u64, (u8::MAX - 1) as u64))
    } else if fits_in!(u16) {
        Ok((Width::U16, u16::MAX as u64, (u16::MAX - 1) as u64))
    } else if fits_in!(u32) {
        Ok((Width::U32, u32::MAX as u64, (u32::MAX - 1) as u64))
    } else if fits_in!(u64) {
        Ok((Width::U64, u64::MAX, u64::MAX - 1))
    } else {
        bail!("k is too large to fit in u64 while keeping sentinel space")
    }
}

/// Choose storage for compact subspace codes plus one missing sentinel.
fn choose_subspace_width(n_selected: usize) -> Result<(Width, u64)> {
    ensure!(
        n_selected > 0,
        "subspace k-mer spec requires at least one selected k-mer"
    );

    if n_selected < u8::MAX as usize {
        Ok((Width::U8, u8::MAX as u64))
    } else if n_selected < u16::MAX as usize {
        Ok((Width::U16, u16::MAX as u64))
    } else if n_selected < u32::MAX as usize {
        Ok((Width::U32, u32::MAX as u64))
    } else {
        Ok((Width::U64, u64::MAX))
    }
}

/// Build a radix-5 spec when k fits the existing packed `u64` representation.
fn build_radix5_spec_if_supported(k: usize) -> Result<Option<KmerSpec>> {
    if !radix5_fits_u64_with_sentinels(k) {
        return Ok(None);
    }
    let (width, sentinel_none, sentinel_n) =
        choose_width(k).with_context(|| format!("calculating dtype for k={k}"))?;
    Ok(Some(KmerSpec {
        k,
        width,
        sentinel_none,
        sentinel_n,
    }))
}

/// Return whether radix-5 can store every real k-mer plus two sentinels in `u64`.
fn radix5_fits_u64_with_sentinels(k: usize) -> bool {
    let mut code_space = 1u128;
    for _ in 0..k {
        let Some(next_code_space) = code_space.checked_mul(5) else {
            return false;
        };
        if next_code_space > u64::MAX as u128 - 1 {
            return false;
        }
        code_space = next_code_space;
    }
    true
}

/// Normalize one exact selected k-mer to uppercase ACGT bytes.
fn normalize_acgt_kmer(seq: &[u8], k: usize) -> Option<Box<[u8]>> {
    if seq.len() != k {
        return None;
    }

    let mut normalized = Vec::with_capacity(seq.len());
    for &base in seq {
        let encoded = encode_base(base);
        if encoded == 4 {
            return None;
        }
        normalized.push(BASES[encoded as usize] as u8);
    }
    Some(normalized.into_boxed_slice())
}

/// Build selected codes for byte-backed subspaces.
/// TODO: Make a rolling-encoder version of this to reduce hashing work
fn build_left_aligned_subspace_codes(spec: &SubspaceKmerSpec, seq: &[u8]) -> Vec<u64> {
    let chrom_len = seq.len();
    if spec.k > chrom_len {
        return vec![spec.sentinel_missing; chrom_len];
    }

    let mut out = Vec::with_capacity(chrom_len);
    for start in 0..=chrom_len - spec.k {
        out.push(spec.encode_kmer_bytes(&seq[start..start + spec.k]));
    }
    out.extend(std::iter::repeat_n(spec.sentinel_missing, spec.k - 1));
    out
}

/// Pack promoted `u64` codes into the selected storage width.
fn pack_codes(raw: Vec<u64>, width: Width) -> KmerCodes {
    match width {
        Width::U8 => KmerCodes::U8(raw.into_iter().map(|code| code as u8).collect()),
        Width::U16 => KmerCodes::U16(raw.into_iter().map(|code| code as u16).collect()),
        Width::U32 => KmerCodes::U32(raw.into_iter().map(|code| code as u32).collect()),
        Width::U64 => KmerCodes::U64(raw),
    }
}

/// Build radix-5 codes for every **left-aligned** k-mer in `seq`.
/// * `sentinel_none` – code for positions where **no** complete k-mer exists
/// * `sentinel_n`   – code for any window that contains an 'N'
///
/// The result length always equals `seq.len()`.
fn build_left_aligned_codes(seq: &[u8], k: usize, sentinel_none: u64, sentinel_n: u64) -> Vec<u64> {
    let chrom_len = seq.len();

    // No complete window fits at all
    if k > chrom_len {
        return vec![sentinel_none; chrom_len];
    }

    // Output will always be exactly chrom_len long
    let mut out = Vec::with_capacity(chrom_len);

    // Rolling-hash helpers
    let highest_place = 5u64.pow((k - 1) as u32); // weight of the left-most digit
    let mut code: u64 = 0; // radix-5 value of current window
    let mut n_in_window: u32 = 0; // 'N' counter in current window

    // First full k-mer window
    for i in 0..k {
        let val = encode_base(seq[i]) as u64;
        if val == 4 {
            n_in_window += 1;
        }
        code = code * 5 + val;
    }
    out.push(if n_in_window > 0 { sentinel_n } else { code });

    // Slide the window through the chromosome
    for i in k..chrom_len {
        // outgoing (left-most) base
        let val_left = encode_base(seq[i - k]) as u64;
        if val_left == 4 {
            n_in_window -= 1;
        }
        code -= val_left * highest_place;

        // shift the remaining k-1 digits one place left (×5)
        code *= 5;

        // incoming (right-most) base
        let val_right = encode_base(seq[i]) as u64;
        if val_right == 4 {
            n_in_window += 1;
        }
        code += val_right;

        out.push(if n_in_window > 0 { sentinel_n } else { code });
    }

    // Pad the tail where no full window fits
    // (exactly k-1 positions)
    out.extend(std::iter::repeat_n(sentinel_none, k - 1));

    debug_assert_eq!(out.len(), chrom_len);
    out
}

/// Decode a code to its k‑mer string, returning 'N'×k for sentinels.
fn decode_kmer(code: u64, k: usize, sentinel_none: u64, sentinel_n: u64) -> String {
    if code == sentinel_none || code == sentinel_n {
        return "N".repeat(k);
    }
    let mut tmp = code;
    let mut buf = vec!['N'; k];
    for pos in (0..k).rev() {
        buf[pos] = BASES[(tmp % 5) as usize];
        tmp /= 5;
    }
    buf.into_iter().collect()
}

#[cfg(test)]
mod tests {
    include!("kmer_codec_tests.rs");
}
