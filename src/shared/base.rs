/// The 5 bases: `A, C, G, T, N`
pub(crate) const BASES: [char; 5] = ['A', 'C', 'G', 'T', 'N'];

/// Shared zeroish tolerance for values that originate from `f32` arithmetic.
pub(crate) const ZEROISH_F32_TOLERANCE: f32 = 2.0 * f32::EPSILON;

/// Shared zeroish tolerance for values that originate from `f64` arithmetic.
pub(crate) const ZEROISH_F64_TOLERANCE: f64 = 2.0 * f64::EPSILON;

/// Clamp tiny finite `f32` values to exact zero using the shared `f32` tolerance.
///
/// This keeps reported derived statistics stable when arithmetic leaves behind
/// very small positive or negative roundoff residues.
#[inline]
#[allow(dead_code)]
pub(crate) fn clamp_close_to_zero_f32(value: f32) -> f32 {
    if value.is_finite() && value.abs() <= ZEROISH_F32_TOLERANCE {
        0.0
    } else {
        value
    }
}

/// Clamp tiny finite `f64` values to exact zero using the shared `f64` tolerance.
///
/// This is the right helper when the value itself originates from `f64` arithmetic.
#[inline]
#[allow(dead_code)]
pub(crate) fn clamp_close_to_zero_f64(value: f64) -> f64 {
    if value.is_finite() && value.abs() <= ZEROISH_F64_TOLERANCE {
        0.0
    } else {
        value
    }
}

/// Clamp tiny finite `f64` values to exact zero using the shared `f32` tolerance.
///
/// Use this when the reported `f64` value is derived from `f32`-originating coverage and should
/// therefore keep the same zeroish threshold as the underlying coverage representation.
#[inline]
#[allow(dead_code)]
pub(crate) fn clamp_close_to_zero_f64_with_f32_threshold(value: f64) -> f64 {
    if value.is_finite() && value.abs() <= ZEROISH_F32_TOLERANCE as f64 {
        0.0
    } else {
        value
    }
}

/// Encode a single nucleotide into its base‑5 digit.
///
/// - A or a -> 0  
/// - C or c -> 1  
/// - G or g -> 2  
/// - T or t -> 3  
/// - anything else -> 4
#[inline(always)]
pub(crate) fn encode_base(b: u8) -> u8 {
    LUT[b as usize]
}

/// Static ASCII->radix-5 lookup table.
/// 0 = A, 1 = C, 2 = G, 3 = T, 4 = N/other
static LUT: [u8; 256] = {
    const N: u8 = 4;
    let mut t = [N; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

// Non-ACGT bytes also map to N, including ambiguity codes, lowercase masked bases, and any other unexpected byte

/// Get the complement of a single nucleotide base.
///
/// - A or a -> T  
/// - C or c -> G  
/// - G or g -> C  
/// - T or t -> A
/// - N or n -> N  
/// - anything else -> identity (return `b`)
#[inline]
pub(crate) fn complement(b: char) -> char {
    match b {
        'A' | 'a' => 'T',
        'T' | 't' => 'A',
        'C' | 'c' => 'G',
        'G' | 'g' => 'C',
        'N' | 'n' => 'N',
        _ => b,
    }
}

/// Reverse-complement of a plain sequence, e.g. "AC" -> "GT"
pub(crate) fn rev_complement(seq: &str) -> String {
    seq.chars().rev().map(complement).collect()
}

/// Complement of a plain sequence, e.g. "AC" -> "TG"
pub(crate) fn complement_seq(seq: &str) -> String {
    seq.chars().map(complement).collect()
}

/// Return the canonical form of `kmer`.
///
/// When `reverse=true`, k-mers are compared against their reverse complement.
/// Otherwise, just their direct complement.
///
/// When `odd_by_center=true`, **odd-length** k-mers are compared by their the middle base,
/// keeping the k-mer as-is if it is `A`, `C`, or `N`, and otherwise returning the
/// complement.
///
/// Otherwise, k-mers are compared against their complement,
/// returning the lexicographically smaller of the two.
#[inline]
pub(crate) fn make_canonical(kmer: String, reverse: bool, odd_by_center: bool) -> String {
    let len = kmer.len();

    if odd_by_center && len % 2 == 1 {
        let mid = kmer.as_bytes()[len / 2].to_ascii_uppercase();
        if mid == b'G' || mid == b'T' {
            if reverse {
                return rev_complement(&kmer);
            } else {
                return complement_seq(&kmer);
            }
        }
        return kmer;
    }

    let compl = if reverse {
        rev_complement(&kmer)
    } else {
        complement_seq(&kmer)
    };
    if kmer <= compl { kmer } else { compl }
}

#[cfg(test)]
mod tests {
    include!("base_tests.rs");
}
