use crate::utils::blacklist::strategy::BlacklistStrategy;

pub fn is_blacklisted(
    blacklist_intervals: &[(u64, u64)],
    blacklist_strategy: BlacklistStrategy,
    start: u64,
    end: u64,
    look_back: u64,
    ptr: &mut usize,
) -> bool {
    // Determine blacklist status
    let in_blacklist = match blacklist_strategy {
        BlacklistStrategy::All => is_all(&blacklist_intervals, start, end, look_back, ptr),
        BlacklistStrategy::Any => is_any(&blacklist_intervals, start, end, look_back, ptr),
        BlacklistStrategy::Midpoint => {
            is_midpoint(&blacklist_intervals, start, end, look_back, ptr)
        }
        BlacklistStrategy::Proportion(th) => {
            is_proportion(&blacklist_intervals, start, end, look_back, ptr, th)
        }
    };
    in_blacklist
}

/// Advance `ptr` to skip any intervals ending before `start`, then
/// sum up how many bases of [start,end) overlap the intervals.
/// `ptr` is left of the first interval that might overlap the next bin.
///
/// Intervals must be sorted by start position and be non‐overlapping.
pub fn compute_blacklist_overlap(
    intervals: &[(u64, u64)],
    start: u64,
    end: u64,
    look_back: u64,
    ptr: &mut usize,
) -> f64 {
    // Skip intervals that end at or before the bin start minus `look_back`
    while *ptr < intervals.len() && intervals[*ptr].1 <= start.saturating_sub(look_back) {
        *ptr += 1;
    }
    // Sum all overlap lengths for this bin
    let mut covered = 0;
    let mut i = *ptr;
    while i < intervals.len() && intervals[i].0 < end {
        let (s, e) = intervals[i];
        covered += e.min(end).saturating_sub(s.max(start));
        i += 1;
    }
    covered as f64 / (end - start) as f64
}

/// Check if the full fragment lies within an interval
pub fn is_all(
    intervals: &[(u64, u64)],
    start: u64,
    end: u64,
    look_back: u64,
    ptr: &mut usize,
) -> bool {
    // Skip any intervals that end entirely before our fragment start minus `look_back`
    while *ptr < intervals.len() && intervals[*ptr].1 <= start.saturating_sub(look_back) {
        *ptr += 1;
    }
    // If there's an interval here, check full overlap of [start,end)
    if let Some(&(s, e)) = intervals.get(*ptr) {
        s <= start && e >= end
    } else {
        false
    }
}

/// Check if fragment midpoint lies within an interval
pub fn is_midpoint(
    intervals: &[(u64, u64)],
    start: u64,
    end: u64,
    look_back: u64,
    ptr: &mut usize,
) -> bool {
    // Compute fragment midpoint
    let mid = start + (end - start) / 2;
    // Skip any intervals that end entirely before the fragment start minus `look_back`
    while *ptr < intervals.len() && intervals[*ptr].1 <= start.saturating_sub(look_back) {
        *ptr += 1;
    }
    // Iterate over every interval that could cover the midpoint
    let mut i = *ptr;
    while i < intervals.len() && intervals[i].0 <= mid {
        let (s, e) = intervals[i];
        if s <= mid && mid < e {
            return true;
        }
        i += 1;
    }
    false
}

/// Returns true if at least `thr` proportion of [start,end) is covered
/// by `intervals`.  `ptr` only skips intervals whose end ≤ start,
/// so we never “lose” an interval that might still overlap the next fragment.
pub fn is_proportion(
    intervals: &[(u64, u64)],
    start: u64,
    end: u64,
    look_back: u64,
    ptr: &mut usize,
    thr: f64,
) -> bool {
    // Skip any intervals that end entirely before our fragment start minus `look_back`
    while *ptr < intervals.len() && intervals[*ptr].1 <= start.saturating_sub(look_back) {
        *ptr += 1;
    }
    // Now *ptr is the first interval that might overlap [start,end)
    let mut covered = 0f64;
    let mut i = *ptr;
    // Sum overlaps only through intervals whose .0 < end
    while i < intervals.len() && intervals[i].0 < end {
        let (s, e) = intervals[i];
        covered += (e.min(end).saturating_sub(s.max(start)) as f64) / ((end - start) as f64);
        if covered >= thr {
            // Early stopping if we reached the threshold
            break;
        }
        i += 1;
    }
    // ptr remains at the first possible overlapping interval
    covered >= thr
}

/// Returns true if at least 1 base is covered by `intervals`.
/// `ptr` only skips intervals whose end ≤ start, so we never
/// “lose” an interval that might still overlap the next fragment.
pub fn is_any(
    intervals: &[(u64, u64)],
    start: u64,
    end: u64,
    look_back: u64,
    ptr: &mut usize,
) -> bool {
    // Skip any intervals that end entirely before our fragment start minus `look_back`
    while *ptr < intervals.len() && intervals[*ptr].1 <= start.saturating_sub(look_back) {
        *ptr += 1;
    }
    // Now *ptr is the first interval that might overlap [start,end)
    let mut covered_bases = 0u16;
    let mut i = *ptr;
    // Sum overlaps only through intervals whose .0 < end
    while i < intervals.len() && intervals[i].0 < end {
        let (s, e) = intervals[i];
        covered_bases += e.min(end).saturating_sub(s.max(start)) as u16;
        if covered_bases > 0 {
            // Early stopping if we reached the threshold
            break;
        }
        i += 1;
    }
    // ptr remains at the first possible overlapping interval
    covered_bases > 0
}
