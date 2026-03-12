/// The 5 bases: `A, C, G, T, N`
pub const BASES: [char; 5] = ['A', 'C', 'G', 'T', 'N'];

/// Encode a single nucleotide into its base‑5 digit.
///
/// - A or a -> 0  
/// - C or c -> 1  
/// - G or g -> 2  
/// - T or t -> 3  
/// - anything else -> 4
#[inline(always)]
pub fn encode_base(b: u8) -> u8 {
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
pub fn complement(b: char) -> char {
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
pub fn rev_complement(seq: &str) -> String {
    seq.chars().rev().map(complement).collect()
}

/// Return the canonical form of `kmer`.
///
/// Even-length k-mers are compared against their reverse complement,
/// returning the lexicographically smaller of the two. For odd-length
/// k-mers we only inspect the middle base, keeping the k-mer as-is if it
/// is `A`, `C`, or `N`, and otherwise returning the reverse complement.
#[inline]
pub fn make_canonical(kmer: String) -> String {
    let len = kmer.len();

    if len % 2 == 1 {
        let mid = kmer.as_bytes()[len / 2].to_ascii_uppercase();
        if mid == b'G' || mid == b'T' {
            return rev_complement(&kmer);
        }
        return kmer;
    }

    let rc = rev_complement(&kmer);
    if kmer <= rc { kmer } else { rc }
}
