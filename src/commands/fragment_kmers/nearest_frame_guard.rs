use crate::shared::kmers::nearest_guard::nearest_guard_bounds;

/// Midpoint guard used by the `nearest` frame to prevent k-mers from straddling the fold.
#[derive(Debug, Clone, Copy)]
pub struct NearestFrameGuard {
    forward_max_start: u64,
    reverse_min_anchor: u64,
    reverse_min_start: u64,
}

impl NearestFrameGuard {
    pub fn by_flag(apply_nearest_guard: bool, fragment_length: u32, k_span: u32) -> Option<Self> {
        if apply_nearest_guard && let Some(bounds) = nearest_guard_bounds(fragment_length, k_span) {
            return Some(Self {
                forward_max_start: bounds.max_forward_start,
                reverse_min_anchor: bounds.min_reverse_anchor,
                reverse_min_start: bounds.min_reverse_start(k_span),
            });
        }
        None
    }

    #[inline]
    pub fn allows_forward(&self, start_offset_0: u64) -> bool {
        start_offset_0 <= self.forward_max_start
    }

    #[inline]
    pub fn allows_reverse_anchor(&self, anchor_offset_0: u64) -> bool {
        anchor_offset_0 >= self.reverse_min_anchor
    }

    #[inline]
    pub fn allows_reverse_start(&self, start_offset_0: u64) -> bool {
        start_offset_0 >= self.reverse_min_start
    }
}
