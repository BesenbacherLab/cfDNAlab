// -- Ref sequence position blacklisting --

/// Byte used for blacklisted bases in the reference sequence
pub const BLACKLIST_BYTE: u8 = b'X';

/// Mask every base that falls inside a blacklist interval with `BLACKLIST_BYTE`.
///
/// * `seq`: mutable byte slice of the reference chromosome  
/// * `intervals`: merged, **sorted**, non-overlapping `[start, end)` pairs  
///
/// Runs in **O(total interval length)** – no per-base scanning.
pub fn apply_blacklist_mask_to_seq(seq: &mut [u8], intervals: &[(u64, u64)]) {
    for &(start, end) in intervals {
        let s = start as usize;
        let e = end as usize;
        // Silent bounds-check: some BEDs can extend past chromosome end
        if s >= seq.len() {
            break;
        }
        let e = e.min(seq.len());
        seq[s..e].fill(BLACKLIST_BYTE);
    }
}
