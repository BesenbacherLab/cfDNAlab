use crate::commands::prepare_windows::config::{OobPolicy, PrepareConfig};

/// Apply resize or flank transform to a window.
///
/// The function returns final coordinates and ensures they are valid with respect
/// to chromosome bounds according to the out-of-bounds policy when a size
/// transform is applied.
///
/// Resizing centers the new window on the midpoint of the original interval.
/// When the input length and target size have different parity (odd/even),
/// there are two equally centered placements. In that case, the code chooses
/// left or right with a deterministic hash based on the midpoint, input length,
/// target size, and optional seed.
///
/// Size Examples
/// -------------
///
/// Interval size 6, resize 4: [011110] -> unique placement
///
/// Interval size 6, resize 3: [001110] or [011100] -> left or right choice
///
/// Interval size 5, resize 4: [11110] or [01111] -> left or right choice
///
/// Interval size 5, resize 3: [01110] -> unique placement
///
/// Parameters
/// ----------
/// - start:
///     Original start.
/// - end:
///     Original end (exclusive).
/// - chrom_size_bp:
///     Chromosome size in base pairs, if trimming or dropping is requested.
/// - cfg:
///     Configuration with resize/flank/oob policy.
///
/// Returns
/// -------
/// - out_start:
///     Transformed start.
/// - out_end:
///     Transformed end (exclusive).
pub fn apply_size_transform(
    start: u32,
    end: u32,
    chrom_size_bp: Option<u32>,
    cfg: &PrepareConfig,
) -> Option<(u32, u32)> {
    let (mut out_start, mut out_end) = (start, end);
    let mut transformed = false;

    if let Some(size) = cfg.resize {
        transformed = true;
        // Center on midpoint, then resolve left or right when parity (odd/even) makes centering ambiguous
        let length = end.saturating_sub(start);
        let midpoint = start + (length / 2);
        let half = size / 2;
        let parity_matches = (length % 2) == (size % 2);
        if parity_matches {
            out_start = midpoint.saturating_sub(half);
            out_end = out_start.saturating_add(size);
        } else {
            // Resolve left or right deterministically to avoid bias
            let decision_seed = cfg.seed.unwrap_or(0);
            let hash =
                fxhash::hash64(&(midpoint as u64, length as u64, size as u64, decision_seed));
            // Low bit acts as a deterministic coin flip for left or right choice
            let choose_left = (hash & 1) == 0;
            let midpoint_i = midpoint as i64;
            let half_i = half as i64;
            let offset = if length % 2 == 0 {
                if choose_left { -1 } else { 0 }
            } else if choose_left {
                0
            } else {
                1
            };
            let start_i = midpoint_i - half_i + offset;
            out_start = if start_i < 0 { 0 } else { start_i as u32 };
            out_end = out_start.saturating_add(size);
        }
    } else if let Some(flanks) = cfg.flank.as_ref() {
        transformed = true;
        let left = flanks[0];
        let right = flanks[1];
        // Allow zero and directionality; negative values are clipped later by policy
        let new_start = (start as i64) - (left as i64);
        let new_end = (end as i64) + (right as i64);
        out_start = if new_start < 0 { 0 } else { new_start as u32 };
        out_end = if new_end < 0 { 0 } else { new_end as u32 };
        if out_end < out_start {
            return None;
        }
    }

    if !transformed {
        return Some((out_start, out_end));
    }

    // Out-of-bounds handling
    match cfg.oob {
        OobPolicy::Allow => Some((out_start, out_end)),
        OobPolicy::Trim => {
            let size = chrom_size_bp.expect("chromosome sizes required for trim/drop policies");
            let s = out_start.min(size);
            let e = out_end.min(size);
            if e <= s { None } else { Some((s, e)) }
        }
        OobPolicy::Drop => {
            let size = chrom_size_bp.expect("chromosome sizes required for trim/drop policies");
            if out_start >= size || out_end > size {
                None
            } else if out_end <= out_start {
                None
            } else {
                Some((out_start, out_end))
            }
        }
    }
}
