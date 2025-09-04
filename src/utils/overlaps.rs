/// Window index and bounds for a single overlap.
#[derive(Debug)]
pub struct OverlappedWindow {
    /// Window index.
    pub idx: usize,
    /// Window start (inclusive).
    pub win_start: u64,
    /// Window end (exclusive).
    pub win_end: u64,
}

/// Collection of windows that overlap a given interval.
#[derive(Debug)]
pub struct OverlappingWindows {
    /// Each window touched by the interval.
    pub windows: Vec<OverlappedWindow>,
    /// Interval start (inclusive).
    pub interval_start: u64,
    /// Interval end (exclusive).
    pub interval_end: u64,
}

/// Return an iterator of 0-based indices `k` for fixed-size bins
/// `[k*bin_size, (k+1)*bin_size)` touched by `[interval_start, interval_end)`.
///
/// Parameters
/// ----------
/// - `interval_start`: Interval start (inclusive).
/// - `interval_end`:   Interval end (exclusive); must be `> interval_start`.
/// - `bin_size`:       Bin size in bases (half-open bins).
///
/// Returns
/// -------
/// Iterator over touched bin indices (0-based).
///
/// Example
/// -------
/// `interval(95, 105, 100) => 0, 1`
#[inline]
pub fn create_overlapping_bins_by_size(
    interval_start: u64,
    interval_end: u64,
    bin_size: u64,
) -> impl Iterator<Item = u64> {
    debug_assert!(interval_end > interval_start, "empty interval");
    let first = interval_start / bin_size;
    let last_excl = (interval_end.saturating_sub(1)) / bin_size + 1;
    first..last_excl
}

/// Half-open interval overlap test. Returns `true` when the intervals share
/// at least one position.
///
/// Parameters
/// ----------
/// - `a_start`, `a_end`: Interval A (start inclusive, end exclusive).
/// - `b_start`, `b_end`: Interval B (start inclusive, end exclusive).
///
/// Returns
/// ---------
/// `true` if the intervals intersect.
#[inline]
pub fn half_open_intervals_overlap(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> bool {
    debug_assert!(a_end >= a_start && b_end >= b_start, "misordered interval");
    b_end > a_start && b_start < a_end
}

/// Enumerate all windows touched by a half-open interval.
///
/// Notes
/// -----
/// - Coordinates are half-open: `[start, end)`.
/// - If `windows` (BED-mode) is provided, it must be **sorted by start** (non-decreasing).
/// - The moving pointer `wd_ptr` is advanced past windows whose end is
///   left of `interval_start - look_back`, enabling streaming even when later
///   intervals may start left of the current one.
/// - In fixed-size mode (`by_size = Some(bin_size)`), bins are not clamped to
///   `chrom_len` (clamp downstream if required).
///
/// Parameters
/// ----------
/// - `chrom_len`: Chromosome length (used only for the global window case).
/// - `wd_ptr`:    Moving pointer into `windows` (BED-mode) for streaming scans.
/// - `windows`:   Optional BED-like windows as `(start, end, original_idx)`.
///                Returned `OverlappedWindow.idx` is the **scan index (`bin_idx`)**, not `original_idx`.
///                Technically, the `original_idx` is ignored and can be any u64.
/// - `by_size`:   If `Some(bin_size)`, use fixed-size bins; otherwise use `windows` if provided.
/// - `interval_start`, `interval_end`: Interval coordinates (start inclusive, end exclusive).
/// - `look_back`: Max look-back distance for advancing `wd_ptr` (e.g., max fragment length).
///
/// Returns
/// -------
/// `Some(OverlappingWindows)` when at least one window is hit; otherwise `None`.
#[inline]
pub fn find_overlapping_windows(
    chrom_len: u64,
    wd_ptr: &mut usize,
    windows: Option<&[(u64, u64, u64)]>, // (start, end, original_idx)
    by_size: Option<u64>,                // bin size for size‑mode
    interval_start: u64,
    interval_end: u64,
    look_back: u64,
) -> Option<OverlappingWindows> {
    // Build window list according to mode
    let mut overlaps = OverlappingWindows {
        windows: Vec::new(),
        interval_start: interval_start,
        interval_end: interval_end,
    };

    // Size‑mode bins
    if let Some(bin_size) = by_size {
        for bin_idx in create_overlapping_bins_by_size(interval_start, interval_end, bin_size) {
            overlaps.windows.push(OverlappedWindow {
                idx: bin_idx as usize,
                win_start: bin_idx * bin_size,
                win_end: bin_idx * bin_size + bin_size,
            });
        }

    // BED‑mode windows
    } else if let Some(window_list) = windows {
        // Skip any intervals that end entirely before the interval start (minus `look_back`)
        // Note that `interval_start` may not be the most left interval position in the outer stash
        while *wd_ptr < window_list.len()
            && window_list[*wd_ptr].1 <= interval_start.saturating_sub(look_back)
        {
            *wd_ptr += 1;
        }
        let mut bin_idx = *wd_ptr;
        while bin_idx < window_list.len() && window_list[bin_idx].0 < interval_end {
            let (win_start, win_end, _) = window_list[bin_idx];
            if half_open_intervals_overlap(interval_start, interval_end, win_start, win_end) {
                overlaps.windows.push(OverlappedWindow {
                    idx: bin_idx,
                    win_start,
                    win_end,
                });
            }
            bin_idx += 1;
        }

    // Global chromosome‑wide window
    } else {
        overlaps.windows.push(OverlappedWindow {
            idx: 0,
            win_start: 0,
            win_end: chrom_len,
        });
    }

    if overlaps.windows.is_empty() {
        return None;
    }

    Some(overlaps)
}
