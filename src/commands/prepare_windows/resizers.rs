use crate::commands::prepare_windows::config::{OobPolicy, PrepareConfig};

/// Apply resize or flank transform to a window.
///
/// The function returns final coordinates and ensures they are valid with respect
/// to chromosome bounds according to the out-of-bounds policy.
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

    if let Some(size) = cfg.resize {
        // Center on midpoint. Odd sizes center exactly; even sizes choose a side randomly
        // to avoid bias if you configured a seed; otherwise, round left deterministically
        // by defaulting to left if no seed was given.
        let midpoint = start + ((end - start) / 2);
        let half = size / 2;
        if size % 2 == 1 {
            out_start = midpoint.saturating_sub(half);
            out_end = midpoint.saturating_add(half) + 1;
        } else {
            // Even size: choose left or right deterministically using a hash that
            // incorporates the midpoint, target size, and optional seed.
            let decision_seed = cfg.seed.unwrap_or(0);
            let hash = fxhash::hash64(&(midpoint as u64, size as u64, decision_seed));
            if (hash & 1) == 0 {
                out_start = midpoint.saturating_sub(half);
                out_end = midpoint.saturating_add(half);
            } else {
                out_start = midpoint.saturating_sub(half.saturating_sub(1));
                out_end = midpoint.saturating_add(half + 1);
            }
        }
    } else if let Some(flanks) = cfg.flank.as_ref() {
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
