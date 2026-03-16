use crate::{
    commands::fragment_kmers::positions::PositionGroup,
    shared::base::{BASES, encode_base, rev_complement},
};
use anyhow::{Context, Result, bail};
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
    /// Human-readable string representation.
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
    /// Code used when the window contains any ‘N’ base
    sentinel_n: u64,
}

impl KmerSpec {
    /// Build per‑position codes for the provided reference sequence.
    pub fn build_left_aligned_codes(&self, seq: &[u8]) -> Vec<u64> {
        build_left_aligned_codes(seq, self.k, self.sentinel_none, self.sentinel_n)
    }

    /// Decode a single code back to its k‑mer string, returning all‑‘N’ if the
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
        // Generic builder returns Vec<u64>
        let raw: Vec<u64> = spec.build_left_aligned_codes(seq);

        // Down-cast into the tightest variant
        let packed = match spec.width {
            Width::U8 => KmerCodes::U8(raw.into_iter().map(|c| c as u8).collect()),
            Width::U16 => KmerCodes::U16(raw.into_iter().map(|c| c as u16).collect()),
            Width::U32 => KmerCodes::U32(raw.into_iter().map(|c| c as u32).collect()),
            Width::U64 => KmerCodes::U64(raw),
        };

        map.insert(*k, packed);
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

/// Build radix-5 codes for every **left-aligned** k-mer in `seq`.
/// * `sentinel_none` – code for positions where **no** complete k-mer exists
/// * `sentinel_n`   – code for any window that contains an ‘N’
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
    let mut n_in_window: u32 = 0; // ‘N’ counter in current window

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

/// Decode a code to its k‑mer string, returning ‘N’×k for sentinels.
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
