// -- Ref sequence position blacklisting --

/// Byte used for blacklisted bases in the reference sequence
pub const BLACKLIST_BYTE: u8 = b'X';

/// Mask every base that falls inside a blacklist interval with `BLACKLIST_BYTE`.
///
/// * `seq`: mutable byte slice of the reference chromosome.
/// * `intervals`: merged, **sorted**, non-overlapping `[start, end)` pairs.
/// * `start_from`: The 0th element of `seq` represents this position in the chromosome.
///
/// Runs in **O(total interval length)** – no per-base scanning.
pub fn apply_blacklist_mask_to_seq(seq: &mut [u8], intervals: &[(u64, u64)], start_from: u64) {
    if seq.is_empty() || intervals.is_empty() {
        return;
    }

    let seq_len = seq.len();
    let seq_end = start_from.saturating_add(seq_len as u64);

    // Skip intervals that end before this slice starts
    let mut idx = 0;
    while idx < intervals.len() && intervals[idx].1 <= start_from {
        idx += 1;
    }

    for &(start, end) in &intervals[idx..] {
        if start >= seq_end {
            break;
        }

        let mask_start = start.saturating_sub(start_from) as usize;

        let capped_end = end.min(seq_end);
        if capped_end <= start_from {
            continue;
        }

        let mask_end = capped_end.saturating_sub(start_from) as usize;
        if mask_end <= mask_start {
            continue;
        }

        seq[mask_start..mask_end].fill(BLACKLIST_BYTE);
    }
}
