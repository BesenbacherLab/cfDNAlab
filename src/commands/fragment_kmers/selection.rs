use crate::commands::fragment_kmers::nearest_frame_guard::NearestFrameGuard;
use crate::commands::fragment_kmers::positions::{
    AllowedWindows, PositionOrientation, PositionSelection, in_any_run,
};

#[derive(Debug, Clone)]
pub enum SelectionDecision {
    IncludeForward {
        start_offset_0: u64,
    },
    IncludeReverse {
        start_offset_0: u64,
        anchor_offset_0: u64,
    },
    SkipAdvance,
}

pub fn evaluate_selection(
    selection: PositionSelection,
    windows: &AllowedWindows,
    guard: Option<&NearestFrameGuard>,
    k_span: u64,
    offset: u64,
    forward_range: Option<(u64, u64)>,
    reverse_range: Option<(u64, u64)>,
) -> SelectionDecision {
    match selection.orientation() {
        PositionOrientation::Forward => {
            if !in_any_run(offset, &windows.forward_starts) {
                return SelectionDecision::SkipAdvance;
            }
            let Some((forward_min, forward_max)) = forward_range else {
                return SelectionDecision::SkipAdvance;
            };
            if offset < forward_min || offset > forward_max {
                return SelectionDecision::SkipAdvance;
            }
            if let Some(guard) = guard {
                if !guard.allows_forward(offset) {
                    return SelectionDecision::SkipAdvance;
                }
            }
            SelectionDecision::IncludeForward {
                start_offset_0: offset,
            }
        }
        PositionOrientation::Reverse => {
            if !in_any_run(offset, &windows.reverse_anchors) {
                return SelectionDecision::SkipAdvance;
            }
            let Some((reverse_min, reverse_max)) = reverse_range else {
                return SelectionDecision::SkipAdvance;
            };
            if offset < reverse_min || offset > reverse_max {
                return SelectionDecision::SkipAdvance;
            }
            if let Some(guard) = guard {
                if !guard.allows_reverse_anchor(offset) {
                    return SelectionDecision::SkipAdvance;
                }
            }
            let start_offset_0 = offset.saturating_sub(k_span.saturating_sub(1));
            if let Some(guard) = guard {
                if !guard.allows_reverse_start(start_offset_0) {
                    return SelectionDecision::SkipAdvance;
                }
            }
            SelectionDecision::IncludeReverse {
                start_offset_0,
                anchor_offset_0: offset,
            }
        }
    }
}
