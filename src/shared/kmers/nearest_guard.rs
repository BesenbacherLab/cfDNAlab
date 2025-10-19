#[derive(Debug, Clone, Copy)]
pub struct NearestGuardBounds {
    /// Inclusive 0-based start index for forward (left) k-mers.
    pub max_forward_start: u64,
    /// Inclusive 0-based anchor index (last base) for reverse (right) k-mers.
    pub min_reverse_anchor: u64,
}

impl NearestGuardBounds {
    /// Inclusive 0-based start index for reverse k-mers.
    pub fn min_reverse_start(&self, k_span: u32) -> u64 {
        let span = k_span as u64;
        if span == 0 {
            0
        } else {
            self.min_reverse_anchor
                .saturating_sub(span.saturating_sub(1))
        }
    }
}

/// Calculate midpoint guard bounds for the `nearest` frame.
///
/// Returns `None` when the k-mer span exceeds the fragment length or when either
/// operand is zero.
pub fn nearest_guard_bounds(length: u32, k_span: u32) -> Option<NearestGuardBounds> {
    if length == 0 || k_span == 0 || k_span > length {
        return None;
    }

    let len = length as u64;
    let span = k_span as u64;
    let half = len / 2;

    let mut max_forward_start = half.saturating_sub(span);
    let mut min_reverse_anchor = if len % 2 == 1 {
        half.saturating_add(span)
    } else {
        half.saturating_add(span.saturating_sub(1))
    };

    // Clamp to valid domain (full length)
    let latest_start = len.saturating_sub(span);
    max_forward_start = max_forward_start.min(latest_start);
    min_reverse_anchor = min_reverse_anchor.min(len.saturating_sub(1));

    Some(NearestGuardBounds {
        max_forward_start,
        min_reverse_anchor,
    })
}
