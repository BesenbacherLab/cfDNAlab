use crate::shared::interval::IndexedInterval;
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
    pub overlap_fraction: f64,
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
    pub fn new(idx: usize, interval: Interval<u64>, overlap_fraction: f64) -> Result<Self> {
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
    pub fn set_overlap_fraction(&mut self, new_fraction: f64) -> Result<()> {
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
/// - `interval`:
///   Checked queried interval.
/// - `bin_size`:
///   Bin size in bases. Must be greater than zero.
///
/// Returns
/// -------
/// - `out`:
///   Half-open range of touched 0-based bin indices.
#[inline]
pub fn create_overlapping_bins_by_size(
    interval: Interval<u64>,
    bin_size: u64,
) -> Result<std::ops::Range<u64>> {
    if bin_size == 0 {
        return Err(Error::InvalidBinSize { bin_size });
    }
    let first = interval.start() / bin_size;
    let last_excl = (interval.end().saturating_sub(1)) / bin_size + 1;
    Ok(first..last_excl)
}

/// Check whether two checked half-open intervals share at least one base.
///
/// Half-open intervals include their start coordinate and exclude their end coordinate, so
/// `[0, 10)` and `[10, 20)` touch but do not overlap. Both inputs are already validated
/// non-empty intervals.
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
    interval_a.intersects(interval_b)
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
/// - `query_interval`:
///   Queried checked non-empty interval.
/// - `min_overlap_fraction`:
///   Minimum one-way overlap fraction for retaining a window, computed as
///   `overlap_bases / query_interval.len()`. This is not the fraction of the
///   window covered by the query. Must be in `[0.0, 1.0]`.
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
    windows: Option<&[IndexedInterval<u64>]>,
    by_size: Option<u64>, // bin size for size‑mode
    query_interval: Interval<u64>,
    min_overlap_fraction: f64,
    look_back: u64,
) -> Result<Option<OverlappingWindows>> {
    if !(0.0..=1.0).contains(&min_overlap_fraction) {
        return Err(Error::OverlapFractionOutOfBounds {
            overlap_fraction: min_overlap_fraction,
        });
    }

    // Build window list according to mode
    let mut overlaps = OverlappingWindows::new(query_interval);

    // Size‑mode bins
    if let Some(bin_size) = by_size {
        for bin_idx in create_overlapping_bins_by_size(query_interval, bin_size)? {
            let window_start = bin_idx * bin_size;
            let window_end = (bin_idx * bin_size + bin_size).min(chrom_len);
            if window_end <= window_start {
                continue;
            }
            let window_interval = Interval::new(window_start, window_end)?;
            let overlap_proportion = fraction_overlap_of_a(query_interval, window_interval);
            if overlap_proportion < min_overlap_fraction {
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
        // Note that `query_interval.start()` may not be the most left interval position in the
        // outer stash
        while *wd_ptr < window_list.len()
            && window_list[*wd_ptr].end() <= query_interval.start().saturating_sub(look_back)
        {
            *wd_ptr += 1;
        }
        let mut bin_idx = *wd_ptr;
        while bin_idx < window_list.len() && window_list[bin_idx].start() < query_interval.end() {
            let window = window_list[bin_idx];
            let win_start = window.start();
            let win_end = window.end().min(chrom_len);
            if win_end <= win_start {
                bin_idx += 1;
                continue;
            }
            let window_interval = Interval::new(win_start, win_end)?;
            if half_open_intervals_overlap(query_interval, window_interval) {
                let overlap_proportion = fraction_overlap_of_a(query_interval, window_interval);
                if overlap_proportion >= min_overlap_fraction {
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

/// Find window overlaps for many same-width intervals on one chromosome.
///
/// Each lookup asks for the windows touched by `[query_start, query_start + query_width)`, clipped
/// to `chrom_len`. This is useful when scanning k-mer starts or other fixed-width intervals in
/// nondecreasing coordinate order.
///
/// Global and fixed-size windows are handled directly. BED mode gives the same answers as
/// [`find_overlapping_windows`] with `look_back = 0`, but keeps a small cache for adjacent query
/// starts. The cache is valid until a BED window enters or leaves the candidate set. Each cached
/// candidate is still rechecked for the current query start, because its overlap fraction can change
/// while the candidate set stays the same. BED-mode callers must provide query starts in
/// nondecreasing order.
///
/// For full-width cached BED queries, a window passes when its overlapping bases meet the requested
/// fraction of `query_width`. A threshold of `0.0` means any actual overlap. Zero-overlap windows
/// are never returned. Queries clipped at the chromosome end use the generic overlap finder so
/// their fractions are measured against the clipped query length.
#[cfg(any(feature = "cmd_ref_kmers", test))]
#[derive(Debug)]
pub(crate) struct FixedWidthOverlapCursor<'a> {
    /// Chromosome length used to clip query intervals and BED window ends.
    chrom_len: u64,
    /// Optional sorted BED-like windows for BED mode.
    ///
    /// The slice must use chromosome coordinates. In BED mode, returned `OverlappingWindow.idx`
    /// values are scan positions in this slice, matching [`find_overlapping_windows`].
    windows: Option<&'a [IndexedInterval<u64>]>,
    /// Optional fixed window size for size mode.
    ///
    /// When this is present, the cursor does not use the BED cache. Fixed-size windows can be found
    /// directly from the query interval.
    by_size: Option<u64>,
    /// Width of the unclipped query interval.
    ///
    /// Ref-kmers uses this as the k-mer size. Caching is only used when the clipped query still has
    /// this full width.
    query_width: u64,
    /// Minimum overlap fraction used to decide whether a window is returned.
    ///
    /// This is the one-way fraction of the query covered by a window, not the fraction of the
    /// window covered by the query. Cached full-width BED queries compare
    /// `overlap_bases / query_width`. Clipped queries and fixed-size windows use the current query
    /// interval length, matching [`find_overlapping_windows`].
    min_overlap_fraction: f64,
    /// Forward scan pointer into `windows` for BED mode.
    ///
    /// This is the first BED window not known to end before the current query start. It only moves
    /// forward, so BED-mode query starts must be requested in nondecreasing order.
    wd_ptr: usize,
    /// Whether `cached_windows` and `next_candidate_change_query_start` describe the current BED
    /// window slice.
    cache_ready: bool,
    /// Candidate BED windows for the current query-start range.
    ///
    /// Each cached window stores the query-start ranges where it passes the overlap threshold and
    /// where its overlap length is constant.
    cached_windows: Vec<CachedFixedWidthWindow>,
    /// First query start where the BED candidate slice may change.
    ///
    /// This is a query-start coordinate, not a genomic window coordinate. The cache is valid while
    /// `query_start < next_candidate_change_query_start`. At this boundary, a BED window may enter
    /// or leave the candidate slice, so the cursor must rescan with `wd_ptr`.
    next_candidate_change_query_start: u64,
    /// Number of BED cache refreshes, used by tests to verify cache reuse.
    refresh_count: usize,
}

#[cfg(any(feature = "cmd_ref_kmers", test))]
impl<'a> FixedWidthOverlapCursor<'a> {
    /// Create a cursor for fixed-width overlap queries.
    ///
    /// `starting_wd_ptr` has the same meaning as the `wd_ptr` argument to
    /// [`find_overlapping_windows`]. BED-mode query starts must be requested in nondecreasing
    /// coordinate order because the cursor uses a forward-only window pointer.
    ///
    /// Parameters
    /// ----------
    /// - `chrom_len`:
    ///   Chromosome length used to clip query intervals and BED window ends.
    /// - `windows`:
    ///   Optional sorted BED-like windows in chromosome coordinates. Use `None` for global or
    ///   fixed-size mode.
    /// - `by_size`:
    ///   Fixed window size for size mode. Use `None` for BED or global mode.
    /// - `query_width`:
    ///   Width of each unclipped query interval. Ref-kmers passes the k-mer size here.
    /// - `min_overlap_fraction`:
    ///   Minimum one-way overlap fraction for keeping a window. Full-width cached BED queries
    ///   measure this as `overlap_bases / query_width`. Clipped queries measure it against the
    ///   clipped query length. This is not the fraction of the window covered by the query. Must be
    ///   in `[0.0, 1.0]`.
    /// - `starting_wd_ptr`:
    ///   Initial BED scan position, usually the first window precomputed for the current tile.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Cursor initialized for repeated fixed-width overlap lookups.
    pub(crate) fn new(
        chrom_len: u64,
        windows: Option<&'a [IndexedInterval<u64>]>,
        by_size: Option<u64>,
        query_width: u64,
        min_overlap_fraction: f64,
        starting_wd_ptr: usize,
    ) -> Result<Self> {
        Interval::new(0, query_width)?;
        if !(0.0..=1.0).contains(&min_overlap_fraction) {
            return Err(Error::OverlapFractionOutOfBounds {
                overlap_fraction: min_overlap_fraction,
            });
        }
        if let Some(bin_size) = by_size {
            if bin_size == 0 {
                return Err(Error::InvalidBinSize { bin_size });
            }
        }

        Ok(Self {
            chrom_len,
            windows,
            by_size,
            query_width,
            min_overlap_fraction,
            wd_ptr: starting_wd_ptr,
            cache_ready: false,
            cached_windows: Vec::new(),
            next_candidate_change_query_start: u64::MAX,
            refresh_count: 0,
        })
    }

    /// Return overlaps for the fixed-width query interval beginning at `query_start`.
    ///
    /// Parameters
    /// ----------
    /// - `query_start`:
    ///   Inclusive chromosome coordinate for the query start. The query end is computed as
    ///   `query_start + query_width` and clipped at `chrom_len`.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Overlapping windows with fractions relative to the queried interval, or `None` when no
    ///   window passes the configured threshold.
    pub(crate) fn find_overlaps(&mut self, query_start: u64) -> Result<Option<OverlappingWindows>> {
        let Some(query_interval) = self.query_interval(query_start)? else {
            return Ok(None);
        };

        if let Some(bin_size) = self.by_size {
            return self.find_by_size(query_interval, bin_size);
        }

        if let Some(windows) = self.windows {
            return self.find_bed(query_interval, windows);
        }

        self.find_global(query_interval)
    }

    /// Return how many times the BED cache has been rebuilt.
    ///
    /// This is test-only instrumentation for checking that adjacent fixed-width lookups reuse the
    /// cache.
    #[cfg(test)]
    pub(crate) fn refresh_count(&self) -> usize {
        self.refresh_count
    }

    /// Build the checked query interval for a requested start.
    ///
    /// The interval is clipped to the chromosome length. Starts at or beyond the chromosome end
    /// have no interval and return `None`. Ref-kmers filters full k-mers before calling this cursor,
    /// so clipped intervals are mainly for generic fixed-width callers.
    fn query_interval(&self, query_start: u64) -> Result<Option<Interval<u64>>> {
        let query_end = query_start
            .saturating_add(self.query_width)
            .min(self.chrom_len);
        if query_end <= query_start {
            return Ok(None);
        }

        Ok(Some(Interval::new(query_start, query_end)?))
    }

    /// Find overlaps when windows are fixed-size bins.
    ///
    /// Fixed-size bins are computed directly from the query interval, so this path does not use the
    /// BED cache. Overlap fractions are measured as `overlap_bases / query_interval.len()`.
    fn find_by_size(
        &self,
        query_interval: Interval<u64>,
        bin_size: u64,
    ) -> Result<Option<OverlappingWindows>> {
        let mut overlaps = OverlappingWindows::new(query_interval);

        for bin_idx in create_overlapping_bins_by_size(query_interval, bin_size)? {
            let window_start = bin_idx * bin_size;
            let window_end = (bin_idx * bin_size + bin_size).min(self.chrom_len);
            if window_end <= window_start {
                continue;
            }
            let window_interval = Interval::new(window_start, window_end)?;
            let overlap_fraction = fraction_overlap_of_a(query_interval, window_interval);
            if overlap_fraction < self.min_overlap_fraction {
                continue;
            }
            overlaps.windows.push(OverlappingWindow::new(
                bin_idx as usize,
                window_interval,
                overlap_fraction,
            )?);
        }

        if overlaps.windows.is_empty() {
            Ok(None)
        } else {
            Ok(Some(overlaps))
        }
    }

    /// Find overlaps against explicit BED-like windows.
    ///
    /// Full-width queries use cached candidate windows. Clipped queries delegate to
    /// [`find_overlapping_windows`] so their overlap fraction uses the clipped query length as the
    /// denominator.
    fn find_bed(
        &mut self,
        query_interval: Interval<u64>,
        windows: &[IndexedInterval<u64>],
    ) -> Result<Option<OverlappingWindows>> {
        if query_interval.len() != self.query_width {
            let mut wd_ptr = self.wd_ptr;
            return find_overlapping_windows(
                self.chrom_len,
                &mut wd_ptr,
                Some(windows),
                None,
                query_interval,
                self.min_overlap_fraction,
                0,
            );
        }

        if !self.cache_valid_for_query_start(query_interval.start()) {
            self.refresh_bed_cache(query_interval, windows)?;
        }

        // Build this query's result from the reusable candidate windows. The cache tells us which
        // BED windows can matter in this query-start range. It does not decide whether each window
        // passes the threshold at this exact start.
        let mut overlaps = OverlappingWindows::new(query_interval);
        for cached_window in &self.cached_windows {
            // Cached windows are candidates for this query-start range, not guaranteed hits for
            // every start in the range. Re-check the accepted range for the current start
            let Some(overlap_bases) = cached_window.overlap_bases_at(query_interval.start()) else {
                continue;
            };

            // Full-width cached BED queries use `query_width` as the denominator. Clipped queries
            // use the generic finder above, so they do not reach this branch
            let overlap_fraction = overlap_bases as f64 / self.query_width as f64;
            overlaps.windows.push(OverlappingWindow::new(
                cached_window.idx,
                cached_window.interval,
                overlap_fraction,
            )?);
        }

        if overlaps.windows.is_empty() {
            Ok(None)
        } else {
            Ok(Some(overlaps))
        }
    }

    /// Return the chromosome-wide window hit for global mode.
    ///
    /// Global mode has a single row for the whole chromosome, so every non-empty query interval
    /// contributes with overlap fraction `1.0`.
    fn find_global(&self, query_interval: Interval<u64>) -> Result<Option<OverlappingWindows>> {
        let mut overlaps = OverlappingWindows::new(query_interval);
        overlaps.windows.push(OverlappingWindow::new(
            0,
            Interval::new(0, self.chrom_len)?,
            1.0,
        )?);
        Ok(Some(overlaps))
    }

    /// Return whether the current BED cache can answer a query beginning at `query_start`.
    ///
    /// `next_candidate_change_query_start` uses the same coordinate system as `query_start`. The
    /// query end does not take part in this check because the cache boundary has already been
    /// converted to the first query start where the candidate slice can change.
    fn cache_valid_for_query_start(&self, query_start: u64) -> bool {
        self.cache_ready && query_start < self.next_candidate_change_query_start
    }

    /// Rebuild cached BED candidates for the current query-start range.
    ///
    /// This advances the forward BED pointer, collects windows that overlap the current query, and
    /// precomputes the query-start ranges where each candidate passes the overlap threshold.
    fn refresh_bed_cache(
        &mut self,
        query_interval: Interval<u64>,
        windows: &[IndexedInterval<u64>],
    ) -> Result<()> {
        // Drop windows that cannot overlap this query or any later query start
        while self.wd_ptr < windows.len() && windows[self.wd_ptr].end() <= query_interval.start() {
            self.wd_ptr += 1;
        }

        // Keep all windows whose start is before the current query end. Together with `wd_ptr`,
        // this makes the half-open candidate slice `[wd_ptr, candidate_end_idx)`
        let mut candidate_end_idx = self.wd_ptr;
        while candidate_end_idx < windows.len()
            && windows[candidate_end_idx].start() < query_interval.end()
        {
            candidate_end_idx += 1;
        }

        self.cached_windows.clear();
        let required_overlap = required_overlap_bases(self.query_width, self.min_overlap_fraction);

        // Convert the next unseen BED start into the first query start that can reach it. If there
        // is no unseen window, the sentinel keeps this out of the minimum
        let mut next_candidate_change_query_start = first_query_start_reaching_window(
            windows
                .get(candidate_end_idx)
                .map(IndexedInterval::start)
                .unwrap_or(u64::MAX),
            self.query_width,
        );

        for window_idx in self.wd_ptr..candidate_end_idx {
            let window = windows[window_idx];

            // Existing candidates leave the slice when the query start reaches their clipped end
            // The cache is valid only until the earliest candidate entry or exit query start
            next_candidate_change_query_start =
                next_candidate_change_query_start.min(later_query_start_or_never(
                    window.end().min(self.chrom_len),
                    query_interval.start(),
                ));

            // Convert this BED window into query-start ranges that can be reused while the
            // candidate slice is unchanged. Windows too short for the threshold are excluded here
            if let Some(cached_window) = CachedFixedWidthWindow::new(
                window_idx,
                window,
                self.chrom_len,
                self.query_width,
                required_overlap,
            )? {
                self.cached_windows.push(cached_window);
            }
        }

        self.next_candidate_change_query_start = next_candidate_change_query_start;
        self.cache_ready = true;
        self.refresh_count += 1;
        Ok(())
    }
}

/// Cached query-start ranges for one candidate BED window and one fixed query width.
///
/// The accepted range controls whether this window can be returned for a query start. The
/// constant-overlap range records the starts where the overlap length is unchanged. Outside that
/// range, moving the query by one base changes the overlap by one base.
#[cfg(any(feature = "cmd_ref_kmers", test))]
#[derive(Debug, Clone, Copy)]
struct CachedFixedWidthWindow {
    idx: usize,
    interval: Interval<u64>,
    query_width: u64,
    /// First query start whose overlap reaches the configured minimum.
    accepted_start: u64,
    /// First query start after the window falls below the configured minimum.
    accepted_end_exclusive: u64,
    /// First query start where overlap length is constant.
    constant_overlap_start: u64,
    /// First query start after overlap length stops being constant.
    constant_overlap_end_exclusive: u64,
    /// Overlap length inside the constant-overlap range.
    ///
    /// This is `query_width` for windows at least as wide as the query and `window_len` for shorter
    /// windows.
    constant_overlap_bases: u64,
}

#[cfg(any(feature = "cmd_ref_kmers", test))]
impl CachedFixedWidthWindow {
    /// Precompute query-start ranges for a candidate BED window.
    ///
    /// The cached window is returned only when it can reach `required_overlap` bases with a
    /// full-width query. Empty or chromosome-clipped-away windows return `None`.
    fn new(
        idx: usize,
        window: IndexedInterval<u64>,
        chrom_len: u64,
        query_width: u64,
        required_overlap: u64,
    ) -> Result<Option<Self>> {
        let window_start = window.start();
        let window_end = window.end().min(chrom_len);
        if window_end <= window_start {
            return Ok(None);
        }

        let window_len = window_end - window_start;
        let max_overlap = query_width.min(window_len);
        if required_overlap == 0 || required_overlap > max_overlap {
            return Ok(None);
        }

        let accepted_start = window_start.saturating_sub(query_width - required_overlap);
        let accepted_end_exclusive = window_end - required_overlap + 1;
        if accepted_end_exclusive <= accepted_start {
            return Ok(None);
        }

        // The constant-overlap range is full query overlap for long windows and full window
        // coverage for short windows
        let (constant_overlap_start, constant_overlap_end_exclusive) = if window_len >= query_width
        {
            (window_start, window_end - query_width + 1)
        } else {
            (window_end.saturating_sub(query_width), window_start + 1)
        };

        Ok(Some(Self {
            idx,
            interval: Interval::new(window_start, window_end)?,
            query_width,
            accepted_start,
            accepted_end_exclusive,
            constant_overlap_start,
            constant_overlap_end_exclusive,
            constant_overlap_bases: max_overlap,
        }))
    }

    /// Return the overlap length at `query_start` when the threshold is met.
    ///
    /// Starts outside the accepted range return `None`. Inside the accepted range, the returned
    /// length is capped at the maximum possible overlap for this query and window.
    fn overlap_bases_at(&self, query_start: u64) -> Option<u64> {
        if query_start < self.accepted_start || query_start >= self.accepted_end_exclusive {
            return None;
        }

        // Outside the constant-overlap range, moving the query by one base changes the overlap by
        // one base. Inside that range, use the precomputed constant length
        let overlap = if query_start < self.constant_overlap_start {
            query_start
                .saturating_add(self.query_width)
                .saturating_sub(self.interval.start())
        } else if query_start < self.constant_overlap_end_exclusive {
            self.constant_overlap_bases
        } else {
            self.interval.end().saturating_sub(query_start)
        };

        Some(overlap.min(self.constant_overlap_bases))
    }
}

/// Convert an overlap-fraction threshold into a required number of bases.
///
/// A threshold of `0.0` means any actual overlap. Other values are expected to have already been
/// validated as finite values in `(0.0, 1.0]`.
#[cfg(any(feature = "cmd_ref_kmers", test))]
fn required_overlap_bases(query_width: u64, min_overlap_fraction: f64) -> u64 {
    if min_overlap_fraction == 0.0 {
        return 1;
    }

    (min_overlap_fraction * query_width as f64).ceil() as u64
}

/// Return the first query start where a BED window can enter the candidate set.
///
/// For a full-width query `[query_start, query_start + query_width)`, a window beginning at
/// `window_start` first overlaps when `query_start + query_width > window_start`. Solving that
/// inequality gives `window_start - query_width + 1`, saturated at zero. `u64::MAX` is used as the
/// no-future-window sentinel.
#[cfg(any(feature = "cmd_ref_kmers", test))]
fn first_query_start_reaching_window(window_start: u64, query_width: u64) -> u64 {
    if window_start == u64::MAX {
        u64::MAX
    } else {
        window_start.saturating_sub(query_width.saturating_sub(1))
    }
}

/// Return `candidate_change_query_start` only when it lies after the current query start.
///
/// Candidate changes at or before the current start have already happened for this lookup and are
/// ignored by returning `u64::MAX`.
#[cfg(any(feature = "cmd_ref_kmers", test))]
fn later_query_start_or_never(candidate_change_query_start: u64, current_query_start: u64) -> u64 {
    if candidate_change_query_start > current_query_start {
        candidate_change_query_start
    } else {
        u64::MAX
    }
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
pub fn fraction_overlap_of_a(interval_a: Interval<u64>, interval_b: Interval<u64>) -> f64 {
    let overlap_bp = overlap_len(interval_a, interval_b) as f64;
    let interval_a_len = interval_a.len() as f64;
    overlap_bp / interval_a_len
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
    interval_a
        .intersection(interval_b)
        .map(|interval| interval.len())
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    include!("overlaps_tests.rs");
}
