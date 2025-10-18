use crate::commands::visualize_positions::ReferenceFrame;
use crate::shared::kmers::nearest_guard::nearest_guard_bounds;

#[derive(Debug, Clone, Copy)]
pub struct KmerFrameGuard {
    forward_max_start: u64,
    reverse_min_anchor: u64,
    reverse_min_start: u64,
}

impl KmerFrameGuard {
    pub fn new(frame: ReferenceFrame, fragment_length: u32, k_span: u32) -> Self {
        if matches!(frame, ReferenceFrame::Nearest) {
            if let Some(bounds) = nearest_guard_bounds(fragment_length, k_span) {
                return Self {
                    forward_max_start: bounds.max_forward_start,
                    reverse_min_anchor: bounds.min_reverse_anchor,
                    reverse_min_start: bounds.min_reverse_start(k_span),
                };
            }
        }
        Self {
            forward_max_start: u64::MAX,
            reverse_min_anchor: 0,
            reverse_min_start: 0,
        }
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
