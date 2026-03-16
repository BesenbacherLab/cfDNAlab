use crate::shared::interval::Interval;
use crate::{Error, Result};

/// A single window hit for one queried interval.
///
/// Use this when you need both the touched window span and the fraction of the
/// queried interval that fell inside that window. The `idx` field is the scan
/// index used by the current window source, so it is useful for looking up
/// per-window state but should not be treated as a stable genomic identifier.
#[derive(Debug)]
pub struct OverlappingWindow {
    /// Window index.
    pub idx: usize,
    /// Window interval as a checked half-open span.
    pub interval: Interval<u64>,
    /// Overlap fraction (overlap_bp / fragment_length_bp)
    pub overlap_fraction: f32,
}

impl OverlappingWindow {
    /// Create one overlap record for a window hit.
    ///
    /// This validates that the overlap fraction is within `[0.0, 1.0]`. The
    /// interval itself must already be a checked non-empty half-open span.
    ///
    /// Parameters
    /// ----------
    /// - `idx`:
    ///   Index of the touched window in the current scan order.
    /// - `interval`:
    ///   Checked non-empty window interval.
    /// - `overlap_fraction`:
    ///   Fraction of the queried interval covered by this window.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   A validated overlap record.
    pub fn new(idx: usize, interval: Interval<u64>, overlap_fraction: f32) -> Result<Self> {
        if !(0.0..=1.0).contains(&overlap_fraction) {
            return Err(Error::OverlapFractionOutOfBounds { overlap_fraction });
        }
        Ok(Self {
            idx,
            interval,
            overlap_fraction,
        })
    }

    /// Return the inclusive window start coordinate.
    ///
    /// This is a convenience accessor for code that still works with separate
    /// start and end coordinates.
    #[inline]
    pub fn start(&self) -> u64 {
        self.interval.start()
    }

    /// Return the exclusive window end coordinate.
    ///
    /// This is a convenience accessor for code that still works with separate
    /// start and end coordinates.
    #[inline]
    pub fn end(&self) -> u64 {
        self.interval.end()
    }

    /// Replace the stored overlap fraction.
    ///
    /// Use this when overlap fractions are computed in more than one pass and
    /// the window record should be reused. The new fraction must stay within
    /// `[0.0, 1.0]`.
    ///
    /// Parameters
    /// ----------
    /// - `new_fraction`:
    ///   New overlap fraction for this window.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   `Ok(())` when the fraction is valid.
    pub fn set_overlap_fraction(&mut self, new_fraction: f32) -> Result<()> {
        if !(0.0..=1.0).contains(&new_fraction) {
            return Err(Error::OverlapFractionOutOfBounds {
                overlap_fraction: new_fraction,
            });
        }
        self.overlap_fraction = new_fraction;
        Ok(())
    }
}

/// All windows touched by one queried interval.
///
/// This groups the original queried interval together with the windows it hit.
/// It is the natural return type for window assignment logic where one query
/// can overlap zero, one, or many windows.
#[derive(Debug)]
pub struct OverlappingWindows {
    /// Each window touched by the interval.
    pub windows: Vec<OverlappingWindow>,
    /// Queried interval as a checked half-open span.
    pub interval: Interval<u64>,
}

impl OverlappingWindows {
    /// Create an empty overlap collection for one queried interval.
    ///
    /// Parameters
    /// ----------
    /// - `interval`:
    ///   Checked non-empty interval that will be compared to windows.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Empty collection ready to receive window hits.
    pub fn new(interval: Interval<u64>) -> Self {
        Self {
            windows: Vec::new(),
            interval,
        }
    }

    /// Return the inclusive queried interval start coordinate.
    ///
    /// This is a convenience accessor for code that still works with separate
    /// start and end coordinates.
    #[inline]
    pub fn query_start(&self) -> u64 {
        self.interval.start()
    }

    /// Return the exclusive queried interval end coordinate.
    ///
    /// This is a convenience accessor for code that still works with separate
    /// start and end coordinates.
    #[inline]
    pub fn query_end(&self) -> u64 {
        self.interval.end()
    }
}

/// Return the fixed-size bins touched by a queried interval.
///
/// Use this when windows are implied by a constant bin size rather than an
/// explicit BED-like list. The returned range contains 0-based bin indices for
/// every bin that shares at least one base with the input interval.
///
/// The input interval must be non-empty. Bin size must be greater than zero.
///
/// Parameters
/// ----------
/// - `interval_start`:
///   Interval start coordinate, inclusive.
/// - `interval_end`:
///   Interval end coordinate, exclusive. Must be greater than `interval_start`.
/// - `bin_size`:
///   Bin size in bases. Must be greater than zero.
///
/// Returns
/// -------
/// - `out`:
///   Half-open range of touched 0-based bin indices.
#[inline]
pub fn create_overlapping_bins_by_size(
    interval_start: u64,
    interval_end: u64,
    bin_size: u64,
) -> Result<std::ops::Range<u64>> {
    if bin_size == 0 {
        return Err(Error::InvalidBinSize { bin_size });
    }
    let interval = Interval::new(interval_start, interval_end)?;
    let first = interval.start() / bin_size;
    let last_excl = (interval.end().saturating_sub(1)) / bin_size + 1;
    Ok(first..last_excl)
}

/// Check whether two checked half-open intervals overlap.
///
/// This is the basic geometric predicate used by the window assignment code.
/// Both intervals are already validated, so the helper only answers whether
/// they share at least one base.
///
/// Parameters
/// ----------
/// - `interval_a`:
///   First checked non-empty interval.
/// - `interval_b`:
///   Second checked non-empty interval.
///
/// Returns
/// -------
/// - `out`:
///   `true` when the intervals intersect.
#[inline]
pub fn half_open_intervals_overlap(interval_a: Interval<u64>, interval_b: Interval<u64>) -> bool {
    interval_b.end() > interval_a.start() && interval_b.start() < interval_a.end()
}

/// Find the windows hit by one queried interval.
///
/// This helper supports the three windowing modes used by the crate:
/// fixed-size bins, explicit BED-like windows, and a single chromosome-wide
/// window. It returns the windows that overlap the queried interval together
/// with overlap fractions measured relative to the queried interval.
///
/// In BED mode the `wd_ptr` pointer is updated in place so repeated calls can
/// stream through sorted windows without rescanning from the beginning. That
/// means the `windows` input must be sorted by start coordinate.
///
/// Parameters
/// ----------
/// - `chrom_len`:
///   Chromosome length. Window ends are clamped to this coordinate.
/// - `wd_ptr`:
///   Moving pointer into `windows` for streaming BED-mode scans.
/// - `windows`:
///   Optional BED-like windows as `(start, end, original_idx)`. The returned
///   `OverlappingWindow.idx` is the scan index, not `original_idx`.
/// - `by_size`:
///   Fixed bin size. When present, fixed-size mode is used instead of
///   `windows`.
/// - `interval_start`:
///   Queried interval start coordinate, inclusive.
/// - `interval_end`:
///   Queried interval end coordinate, exclusive.
/// - `min_overlap_fraction`:
///   Minimum fraction of the queried interval that must overlap a window for
///   that window to be retained.
/// - `look_back`:
///   Maximum distance used when advancing `wd_ptr` in BED mode.
///
/// Returns
/// -------
/// - `out`:
///   `Some(OverlappingWindows)` when at least one window is hit, otherwise
///   `None`.
#[inline]
pub fn find_overlapping_windows(
    chrom_len: u64,
    wd_ptr: &mut usize,
    windows: Option<&[(u64, u64, u64)]>, // (start, end, original_idx)
    by_size: Option<u64>,                // bin size for size‑mode
    interval_start: u64,
    interval_end: u64,
    min_overlap_fraction: f64,
    look_back: u64,
) -> Result<Option<OverlappingWindows>> {
    let query_interval = Interval::new(interval_start, interval_end)?;

    // Build window list according to mode
    let mut overlaps = OverlappingWindows::new(query_interval);

    // Size‑mode bins
    if let Some(bin_size) = by_size {
        for bin_idx in create_overlapping_bins_by_size(interval_start, interval_end, bin_size)? {
            let window_start = bin_idx * bin_size;
            let window_end = (bin_idx * bin_size + bin_size).min(chrom_len);
            if window_end <= window_start {
                continue;
            }
            let window_interval = Interval::new(window_start, window_end)?;
            let overlap_proportion = fraction_overlap_of_a(query_interval, window_interval);
            if (overlap_proportion as f64) < min_overlap_fraction {
                continue;
            }
            overlaps.windows.push(OverlappingWindow::new(
                bin_idx as usize,
                window_interval,
                overlap_proportion,
            )?);
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
            let (win_start, mut win_end, _) = window_list[bin_idx];
            win_end = win_end.min(chrom_len);
            if win_end <= win_start {
                bin_idx += 1;
                continue;
            }
            let window_interval = Interval::new(win_start, win_end)?;
            if half_open_intervals_overlap(query_interval, window_interval) {
                let overlap_proportion = fraction_overlap_of_a(query_interval, window_interval);
                if (overlap_proportion as f64) >= min_overlap_fraction {
                    overlaps.windows.push(OverlappingWindow::new(
                        bin_idx,
                        window_interval,
                        overlap_proportion,
                    )?);
                }
            }
            bin_idx += 1;
        }

    // Global chromosome‑wide window
    } else {
        overlaps.windows.push(OverlappingWindow::new(
            0,
            Interval::new(0, chrom_len)?,
            1.0,
        )?);
    }

    if overlaps.windows.is_empty() {
        return Ok(None);
    }

    Ok(Some(overlaps))
}

/// Compute how much of `interval_a` is covered by `interval_b`.
///
/// The returned value is measured relative to `interval_a`, not symmetrically
/// across both intervals. This is useful when the queried interval is the
/// denominator and windows are being filtered by minimum overlap fraction.
///
/// Parameters
/// ----------
/// - `interval_a`:
///   Checked non-empty interval used as the denominator.
/// - `interval_b`:
///   Checked non-empty interval whose overlap with `interval_a` is measured.
///
/// Returns
/// -------
/// - `out`:
///   Fraction of `interval_a` covered by `interval_b`, in `[0.0, 1.0]`.
#[inline]
pub fn fraction_overlap_of_a(interval_a: Interval<u64>, interval_b: Interval<u64>) -> f32 {
    let overlap_bp = overlap_len(interval_a, interval_b) as f64;
    let interval_a_len = (interval_a.end() - interval_a.start()) as f64;
    (overlap_bp / interval_a_len) as f32
}

/// Compute the number of overlapping bases shared by two intervals.
///
/// This is a low-level helper for overlap-based calculations. It returns zero
/// when the intervals do not intersect.
///
/// Parameters
/// ----------
/// - `interval_a`:
///   First checked non-empty interval.
/// - `interval_b`:
///   Second checked non-empty interval.
///
/// Returns
/// -------
/// - `out`:
///   Number of shared bases.
#[inline]
pub fn overlap_len(interval_a: Interval<u64>, interval_b: Interval<u64>) -> u64 {
    let start = interval_a.start().max(interval_b.start());
    let end = interval_a.end().min(interval_b.end());
    end.saturating_sub(start)
}
