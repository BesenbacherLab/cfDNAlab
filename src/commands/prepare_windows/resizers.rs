use anyhow::{Result, bail};

use crate::commands::prepare_windows::config::{OobPolicy, PrepareConfig};
use crate::shared::interval::Interval;

/// Apply resize or flank transform to a window.
///
/// The function returns final coordinates and ensures they are valid with respect
/// to chromosome bounds according to the out-of-bounds policy when `chrom_size_bp`
/// is specified (underflow always checked).
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
/// - out:
///     Transformed checked interval, or `None` when the chosen policy drops it.
pub(crate) fn apply_size_transform(
    start: u32,
    end: u32,
    chrom_size_bp: Option<u32>,
    cfg: &PrepareConfig,
) -> Result<Option<Interval<u32>>> {
    let mut transformed = false;
    let mut intended_start_i: i64 = start as i64;
    let mut intended_end_i: i64 = end as i64;

    if let Some(size) = cfg.resize {
        transformed = true;
        // Center on midpoint, then resolve left or right when parity (odd/even) makes centering ambiguous
        let length = end.saturating_sub(start);
        let midpoint = start + (length / 2);
        let half = size / 2;
        let parity_matches = (length % 2) == (size % 2);
        if parity_matches {
            intended_start_i = midpoint.saturating_sub(half) as i64;
            intended_end_i = intended_start_i + size as i64;
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
            intended_start_i = midpoint_i - half_i + offset;
            intended_end_i = intended_start_i + size as i64;
        }
    } else if let Some(flanks) = cfg.flank.as_ref() {
        transformed = true;
        let left = flanks[0];
        let right = flanks[1];
        // Allow zero and directionality. Negative values are clipped later by policy
        let new_start = (start as i64) - (left as i64);
        let new_end = (end as i64) + (right as i64);
        intended_start_i = new_start;
        intended_end_i = new_end;
        if intended_end_i < intended_start_i {
            return Ok(None);
        }
    }

    let underflow = intended_start_i < 0;

    // Handle underflow (drop or trim)
    if underflow {
        match cfg.oob {
            OobPolicy::Allow => {
                eprintln!("Warning: window underflowed chromosome start. Dropping.");
                return Ok(None);
            }
            OobPolicy::Drop => {
                return Ok(None);
            }
            OobPolicy::Trim => {
                intended_start_i = 0;
            }
        }
    }

    if matches!(cfg.oob, OobPolicy::Allow) {
        // No check is requested
        return Ok(Some(Interval::new(
            intended_start_i as u32,
            intended_end_i as u32,
        )?));
    }

    // Resizing/flanking without chromosome sizes is not allowed: we cannot
    // enforce OOB policies or preserve centering correctly
    if transformed && chrom_size_bp.is_none() {
        bail!(
            "resize/flank requested without chromosome sizes. Provide chrom sizes to enforce OOB policy"
        );
    }

    let overflow = chrom_size_bp
        .map(|size| intended_end_i > size as i64)
        .unwrap_or(false);

    // When bounds are unknown, don't check overflow
    if chrom_size_bp.is_none() || !overflow {
        return Ok(Some(Interval::new(
            intended_start_i as u32,
            intended_end_i as u32,
        )?));
    }

    let size = chrom_size_bp.expect("chromosome sizes required for trim/drop policies");
    match cfg.oob {
        OobPolicy::Trim => {
            let trimmed_start = intended_start_i as u32;
            let trimmed_end = intended_end_i.min(size as i64) as u32;
            if trimmed_end <= trimmed_start {
                Ok(None)
            } else {
                Ok(Some(Interval::new(trimmed_start, trimmed_end)?))
            }
        }
        OobPolicy::Drop => Ok(None),
        _ => unreachable!("OobPolicy::Allow already handled"),
    }
}

#[cfg(test)]
mod tests {
    include!("resizers_tests.rs");
}
