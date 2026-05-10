use crate::shared::blacklist::strategy::BlacklistStrategy;
use crate::shared::interval::Interval;

pub fn is_blacklisted(
    blacklist_intervals: &[Interval<u64>],
    blacklist_strategy: BlacklistStrategy,
    interval: Interval<u64>,
    look_back: u64,
    ptr: &mut usize,
) -> bool {
    // Determine blacklist status

    match blacklist_strategy {
        BlacklistStrategy::All => is_all(blacklist_intervals, interval, look_back, ptr),
        BlacklistStrategy::Any => is_any(blacklist_intervals, interval, look_back, ptr),
        BlacklistStrategy::Midpoint => is_midpoint(blacklist_intervals, interval, look_back, ptr),
        BlacklistStrategy::Proportion(th) => {
            is_proportion(blacklist_intervals, interval, look_back, ptr, th)
        }
    }
}

/// Advance `ptr` to skip any intervals ending before `start`, then
/// sum up how many bases of [start,end) overlap the intervals.
/// `ptr` is left of the first interval that might overlap the next bin.
///
/// Intervals must be sorted by start position and be non‐overlapping.
pub fn compute_blacklist_overlap(
    intervals: &[Interval<u64>],
    interval: Interval<u64>,
    look_back: u64,
    ptr: &mut usize,
) -> f64 {
    let (start, end) = interval.as_tuple();
    // Skip intervals that end at or before the bin start minus `look_back`
    while *ptr < intervals.len() && intervals[*ptr].end() <= start.saturating_sub(look_back) {
        *ptr += 1;
    }
    // Sum all overlap lengths for this bin
    let mut covered = 0;
    let mut i = *ptr;
    while i < intervals.len() && intervals[i].start() < end {
        let blacklist_interval = intervals[i];
        let (blacklist_start, blacklist_end) = blacklist_interval.as_tuple();
        covered += blacklist_end
            .min(end)
            .saturating_sub(blacklist_start.max(start));
        i += 1;
    }
    covered as f64 / (interval.len() as f64)
}

/// Check if the full fragment lies within an interval
pub fn is_all(
    intervals: &[Interval<u64>],
    interval: Interval<u64>,
    look_back: u64,
    ptr: &mut usize,
) -> bool {
    // Skip any intervals that end entirely before our fragment start minus `look_back`
    while *ptr < intervals.len()
        && intervals[*ptr].end() <= interval.start().saturating_sub(look_back)
    {
        *ptr += 1;
    }
    // If there's an interval here, check full overlap of [start,end)
    if let Some(blacklist_interval) = intervals.get(*ptr) {
        blacklist_interval.contains_interval(interval)
    } else {
        false
    }
}

/// Check if the fragment midpoint support overlaps a blacklist interval.
///
/// Odd-length fragments have one central base. Even-length fragments have two
/// central bases in discrete base coordinates, so either central base is enough
/// to blacklist the fragment.
pub fn is_midpoint(
    intervals: &[Interval<u64>],
    interval: Interval<u64>,
    look_back: u64,
    ptr: &mut usize,
) -> bool {
    let (start, end) = interval.as_tuple();
    let length = end - start;
    let right_center = start + length / 2;
    let left_center = if length % 2 == 0 {
        right_center - 1
    } else {
        right_center
    };

    // Skip any intervals that end entirely before the fragment start minus `look_back`
    while *ptr < intervals.len() && intervals[*ptr].end() <= start.saturating_sub(look_back) {
        *ptr += 1;
    }
    // Iterate over every interval that could cover one of the central bases
    let mut i = *ptr;
    while i < intervals.len() && intervals[i].start() <= right_center {
        if intervals[i].contains_point(left_center) || intervals[i].contains_point(right_center) {
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
    intervals: &[Interval<u64>],
    interval: Interval<u64>,
    look_back: u64,
    ptr: &mut usize,
    thr: f64,
) -> bool {
    let (start, end) = interval.as_tuple();
    // Skip any intervals that end entirely before our fragment start minus `look_back`
    while *ptr < intervals.len() && intervals[*ptr].end() <= start.saturating_sub(look_back) {
        *ptr += 1;
    }
    // Now *ptr is the first interval that might overlap [start,end)
    let mut covered = 0f64;
    let mut i = *ptr;
    // Sum overlaps only through intervals whose .0 < end
    while i < intervals.len() && intervals[i].start() < end {
        let blacklist_interval = intervals[i];
        covered += blacklist_interval
            .intersection(interval)
            .map(|shared_interval| shared_interval.len() as f64 / interval.len() as f64)
            .unwrap_or(0.0);
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
    intervals: &[Interval<u64>],
    interval: Interval<u64>,
    look_back: u64,
    ptr: &mut usize,
) -> bool {
    let (start, end) = interval.as_tuple();
    // Skip any intervals that end entirely before our fragment start minus `look_back`
    while *ptr < intervals.len() && intervals[*ptr].end() <= start.saturating_sub(look_back) {
        *ptr += 1;
    }
    // Now *ptr is the first interval that might overlap [start,end)
    let mut covered_bases = 0u16;
    let mut i = *ptr;
    // Sum overlaps only through intervals whose .0 < end
    while i < intervals.len() && intervals[i].start() < end {
        let blacklist_interval = intervals[i];
        covered_bases += blacklist_interval
            .intersection(interval)
            .map(|shared_interval| shared_interval.len() as u16)
            .unwrap_or(0);
        if covered_bases > 0 {
            // Early stopping if we reached the threshold
            break;
        }
        i += 1;
    }
    // ptr remains at the first possible overlapping interval
    covered_bases > 0
}

#[cfg(test)]
mod tests {
    include!("overlaps_tests.rs");
}
