use crate::shared::interval::IndexedInterval;
use crate::shared::interval::Interval;
#[cfg(uses_bed_window_tier_helpers)]
use crate::shared::tiled_run::{Tile, TileWindowSpan, precompute_tile_window_spans};
use crate::{Error, Result};
#[cfg(uses_bed_window_tier_helpers)]
use fxhash::FxHashMap;
#[cfg(uses_bed_window_tier_helpers)]
use rayon::prelude::*;

/// Default minimum length for treating a BED-like window as broad.
///
/// Broad windows get their own tile-local pointer so they cannot keep the narrow-window pointer
/// pinned behind nested short intervals.
#[cfg(uses_bed_window_tier_helpers)]
pub(crate) const DEFAULT_BROAD_WINDOW_MIN_BP: u64 = 100_000;

/// A single window hit for one queried interval.
///
/// Use this when you need both the touched window span and the fraction of the
/// queried interval that fell inside that window. The `idx` field is the scan
/// index used by the current window source, so it is useful for looking up
/// per-window state but should not be treated as a stable genomic identifier.
/// BED-like overlap finders also carry `output_idx`, which is the row identity
/// used for output aggregation.
#[derive(Debug)]
pub struct OverlappingWindow {
    /// Source-window position or bin index used by the current window source.
    pub idx: usize,
    /// Optional stable output identity for BED-like windows.
    pub output_idx: Option<u64>,
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
        Self::new_with_output_idx(idx, interval, overlap_fraction, None)
    }

    /// Create one overlap record for a BED-like window hit with an output identity.
    ///
    /// `idx` stays the `all_windows` position used for local count arrays. `output_idx` is the
    /// original BED row id for ordinary BED windows and the group index for grouped BED windows.
    pub(crate) fn new_with_output_idx(
        idx: usize,
        interval: Interval<u64>,
        overlap_fraction: f64,
        output_idx: Option<u64>,
    ) -> Result<Self> {
        if !(0.0..=1.0).contains(&overlap_fraction) {
            return Err(Error::OverlapFractionOutOfBounds { overlap_fraction });
        }
        Ok(Self {
            idx,
            output_idx,
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

/// One BED-like window with both local `all_windows` identity and output identity.
///
/// `all_windows_idx` is the position in the chromosome-local `all_windows` list. It is used for
/// tile-local count arrays. `output_idx` is the identity carried by the original `IndexedInterval`:
/// the BED row id for ordinary BED input and the group index for grouped BED.
#[cfg(uses_bed_window_tier_helpers)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BedWindowEntry {
    pub(crate) interval: Interval<u64>,
    pub(crate) all_windows_idx: usize,
    pub(crate) output_idx: u64,
}

#[cfg(uses_bed_window_tier_helpers)]
impl BedWindowEntry {
    /// Build a BED entry from a window and its chromosome-local `all_windows` position.
    #[inline]
    pub(crate) fn from_indexed_interval(
        all_windows_idx: usize,
        window: IndexedInterval<u64>,
    ) -> Self {
        Self {
            interval: window.interval,
            all_windows_idx,
            output_idx: window.idx(),
        }
    }

    /// Return the inclusive start coordinate.
    #[inline]
    pub(crate) fn start(&self) -> u64 {
        self.interval.start()
    }

    /// Return the exclusive end coordinate.
    #[inline]
    pub(crate) fn end(&self) -> u64 {
        self.interval.end()
    }

    /// Return the window length in bases.
    #[inline]
    pub(crate) fn len(&self) -> u64 {
        self.interval.len()
    }

    #[inline]
    fn with_interval(self, interval: Interval<u64>) -> Self {
        Self { interval, ..self }
    }
}

#[cfg(uses_bed_window_tier_helpers)]
impl AsRef<Interval<u64>> for BedWindowEntry {
    #[inline]
    fn as_ref(&self) -> &Interval<u64> {
        &self.interval
    }
}

/// Size class for a BED window tier.
#[cfg(uses_bed_window_tier_helpers)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BedWindowTierKind {
    /// Windows with length greater than or equal to `broad_window_min_bp`.
    Broad,
    /// Windows shorter than `broad_window_min_bp`.
    Narrow,
}

/// One start-sorted BED window tier used by the tile-local overlap finder.
#[cfg(uses_bed_window_tier_helpers)]
#[derive(Debug, Clone)]
pub(crate) struct BedWindowTier {
    /// Size class used for tier-specific fast paths.
    pub(crate) kind: BedWindowTierKind,
    pub(crate) windows: Vec<BedWindowEntry>,
}

/// Chromosome-local BED windows used by `lengths` and `ends`.
///
/// `all_windows` stores the full chromosome-local, start-sorted BED window list used for count
/// arrays and fetch narrowing. `tiers` stores split views of the same entries in independent
/// start-sorted scan lists: currently broad windows and narrow windows. Each tier carries its size
/// class, so fast paths can depend on `BedWindowTierKind` rather than a vector position. Entries
/// inside the tiers keep their `all_windows_idx`, so overlap results can be mapped back to the full
/// list even though each tier advances through its own pointer.
#[cfg(uses_bed_window_tier_helpers)]
#[derive(Debug, Clone, Default)]
pub(crate) struct ChromosomeBedWindows {
    pub(crate) all_windows: Vec<BedWindowEntry>,
    pub(crate) tiers: Vec<BedWindowTier>,
}

#[cfg(uses_bed_window_tier_helpers)]
impl ChromosomeBedWindows {
    /// Build the current broad/narrow tier views from start-sorted BED-like windows.
    ///
    /// `broad_window_min_bp` is the length threshold for the size tiers. Windows with length
    /// greater than or equal to this threshold go into the broad tier. Shorter windows go into the
    /// narrow tier. Later tile-local splitting reuses these tiers and does not check the threshold
    /// again.
    pub(crate) fn from_indexed_windows(
        windows: &[IndexedInterval<u64>],
        broad_window_min_bp: u64,
    ) -> Self {
        let all_windows: Vec<BedWindowEntry> = windows
            .iter()
            .copied()
            .enumerate()
            .map(|(all_windows_idx, window)| {
                BedWindowEntry::from_indexed_interval(all_windows_idx, window)
            })
            .collect();

        let mut broad_windows = Vec::new();
        let mut narrow_windows = Vec::new();
        for window in all_windows.iter().copied() {
            // This is the broad/narrow threshold check for BED tiering
            if window.len() >= broad_window_min_bp {
                broad_windows.push(window);
            } else {
                narrow_windows.push(window);
            }
        }

        Self {
            all_windows,
            tiers: vec![
                BedWindowTier {
                    kind: BedWindowTierKind::Broad,
                    windows: broad_windows,
                },
                BedWindowTier {
                    kind: BedWindowTierKind::Narrow,
                    windows: narrow_windows,
                },
            ],
        }
    }
}

/// Build chromosome-local BED window collections in parallel.
///
/// Each input value provides a sorted interval slice. Construction is independent by chromosome,
/// so callers should initialize the shared Rayon pool with their requested thread count before
/// calling this function.
#[cfg(uses_bed_window_tier_helpers)]
pub(crate) fn build_bed_windows_by_chr<T>(
    windows_by_chr: &FxHashMap<String, T>,
    broad_window_min_bp: u64,
) -> FxHashMap<String, ChromosomeBedWindows>
where
    T: AsRef<[IndexedInterval<u64>]> + Sync,
{
    windows_by_chr
        .par_iter()
        .map(|(chr, windows)| {
            (
                chr.clone(),
                ChromosomeBedWindows::from_indexed_windows(windows.as_ref(), broad_window_min_bp),
            )
        })
        .collect()
}

/// Cached tile spans for `all_windows` and each tier.
#[cfg(uses_bed_window_tier_helpers)]
#[derive(Clone, Debug)]
pub(crate) struct TileBedWindowSpans {
    #[cfg_attr(not(uses_tile_bed_overlap_context), allow(dead_code))]
    pub(crate) all_windows_span: Option<TileWindowSpan>,
    pub(crate) tier_spans: Vec<Option<TileWindowSpan>>,
}

/// Tile-local view of chromosome BED windows.
///
/// The windows stay in chromosome-level storage. The spans select the `all_windows` and tier
/// ranges that can matter for one tile.
#[cfg(uses_bed_window_tier_helpers)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct TileBedWindowView<'a> {
    pub(crate) chromosome_windows: &'a ChromosomeBedWindows,
    pub(crate) spans: &'a TileBedWindowSpans,
}

/// Precompute tile-local `all_windows` and tier spans for BED-like windows.
#[cfg(uses_bed_window_tier_helpers)]
pub(crate) fn precompute_tile_bed_window_spans(
    tiles: &[Tile],
    bed_windows_by_chr: &FxHashMap<String, ChromosomeBedWindows>,
    left_halo: u64,
    right_halo: u64,
) -> Vec<TileBedWindowSpans> {
    let all_windows_spans = precompute_tile_window_spans(
        tiles,
        |chr| {
            bed_windows_by_chr
                .get(chr)
                .map(|windows| windows.all_windows.as_slice())
                .unwrap_or(&[])
        },
        left_halo,
        right_halo,
    );

    let tier_count = bed_windows_by_chr
        .values()
        .map(|windows| windows.tiers.len())
        .max()
        .unwrap_or(0);

    let mut tile_spans: Vec<TileBedWindowSpans> = all_windows_spans
        .into_iter()
        .map(|all_windows_span| TileBedWindowSpans {
            all_windows_span,
            tier_spans: vec![None; tier_count],
        })
        .collect();

    for tier_idx in 0..tier_count {
        let tier_spans = precompute_tile_window_spans(
            tiles,
            |chr| {
                bed_windows_by_chr
                    .get(chr)
                    .and_then(|windows| windows.tiers.get(tier_idx))
                    .map(|tier| tier.windows.as_slice())
                    .unwrap_or(&[])
            },
            left_halo,
            right_halo,
        );

        for (tile_span, tier_span) in tile_spans.iter_mut().zip(tier_spans) {
            tile_span.tier_spans[tier_idx] = tier_span;
        }
    }

    tile_spans
}

#[cfg(uses_tile_bed_overlap_context)]
struct TileBedTierCursor {
    windows: Vec<BedWindowEntry>,
    wd_ptr: usize,
}

/// Tile-local BED candidates split into always-hit windows and scanned windows.
///
/// `always_hit_windows` contains broad BED rows whose clipped interval covers the full tile
/// assignment range. Every query assigned to the tile lies inside that range, so these broad
/// windows do not need per-query scanning. `scanned_windows_by_size_tier` contains all other
/// candidates, still grouped by the broad/narrow size tier chosen by `broad_window_min_bp`.
#[cfg(uses_bed_window_tier_helpers)]
#[derive(Debug)]
struct TileBedWindowSplit {
    always_hit_windows: Vec<BedWindowEntry>,
    scanned_windows_by_size_tier: Vec<Vec<BedWindowEntry>>,
}

/// Split already-tiered tile candidates into always-hit windows and windows to scan.
///
/// Only broad windows are checked for the always-hit path. In production tile runs, narrow windows
/// are not expected to cover the full tile assignment range, so testing every narrow candidate for
/// that condition is dead work. Narrow tiers are copied directly into
/// `scanned_windows_by_size_tier`.
///
/// Empty size tiers are valid and remain empty in `scanned_windows_by_size_tier`. A broad row is
/// moved to `always_hit_windows` only after its end has been clipped to `chrom_len`, because
/// overlap rows must never expose bases beyond the chromosome.
#[cfg(uses_bed_window_tier_helpers)]
fn split_tile_bed_candidates_into_always_hit_and_scanned<'a, I>(
    chrom_len: u64,
    windows_by_size_tier: I,
    tile_assignment_envelope: Interval<u64>,
) -> Result<TileBedWindowSplit>
where
    I: IntoIterator<Item = (BedWindowTierKind, &'a [BedWindowEntry])>,
{
    let mut always_hit_windows = Vec::new();
    let mut scanned_windows_by_size_tier = Vec::new();

    for (tier_kind, windows_in_size_tier) in windows_by_size_tier {
        if tier_kind != BedWindowTierKind::Broad {
            scanned_windows_by_size_tier.push(windows_in_size_tier.to_vec());
            continue;
        }

        let mut scanned_windows = Vec::with_capacity(windows_in_size_tier.len());
        for window in windows_in_size_tier {
            let Some(clipped_window_interval) = clamp_bed_window_to_chrom(*window, chrom_len)?
            else {
                continue;
            };
            if clipped_window_interval.contains_interval(tile_assignment_envelope) {
                // This window contains every query assigned to the tile
                always_hit_windows.push(window.with_interval(clipped_window_interval));
            } else {
                // Keep original bounds so pointer retirement matches the generic BED finder
                scanned_windows.push(*window);
            }
        }

        scanned_windows_by_size_tier.push(scanned_windows);
    }

    Ok(TileBedWindowSplit {
        always_hit_windows,
        scanned_windows_by_size_tier,
    })
}

/// Split tile-local BED candidates once so overlap finders can share the same tier selection.
///
/// `spans.tier_spans` must have an entry for every chromosome BED tier. An entry may be `None`,
/// which means that tier has no candidate windows for this tile. A missing entry would silently
/// drop a whole tier, so it is treated as an internal span-cache error.
#[cfg(uses_bed_window_tier_helpers)]
fn split_tile_bed_windows(
    chrom_len: u64,
    chromosome_windows: &ChromosomeBedWindows,
    spans: &TileBedWindowSpans,
    tile_assignment_envelope: Interval<u64>,
) -> Result<TileBedWindowSplit> {
    if spans.tier_spans.len() < chromosome_windows.tiers.len() {
        return Err(Error::InvalidBedWindowTierSpanCount {
            tier_count: chromosome_windows.tiers.len(),
            span_count: spans.tier_spans.len(),
        });
    }

    let mut windows_by_size_tier = Vec::with_capacity(chromosome_windows.tiers.len());
    for (tier_idx, tier) in chromosome_windows.tiers.iter().enumerate() {
        let candidate_span = spans.tier_spans.get(tier_idx).copied().flatten();
        windows_by_size_tier.push((
            tier.kind,
            span_slice(tier.windows.as_slice(), candidate_span, "BED size tier")?,
        ));
    }

    split_tile_bed_candidates_into_always_hit_and_scanned(
        chrom_len,
        windows_by_size_tier,
        tile_assignment_envelope,
    )
}

/// Tile-local BED overlap finder with independent tier pointers.
///
/// Use this only for BED-like count windows. It returns the same candidate windows as the generic
/// BED finder, but each tier is scanned through a separate pointer so long nested intervals do not
/// pin unrelated windows. Returned `OverlappingWindow.idx` values are chromosome-local positions
/// in `all_windows`. `OverlappingWindow.output_idx` carries the BED row id or grouped BED group
/// index.
///
/// The returned window order is an implementation detail and may differ from
/// [`find_overlapping_windows`]. Callers must treat overlap rows as a set keyed by `idx` plus
/// interval, not as an ordered stream.
#[cfg(uses_tile_bed_overlap_context)]
pub(crate) struct TileBedOverlapContext {
    chrom_len: u64,
    always_hit_windows: Vec<BedWindowEntry>,
    tier_cursors: Vec<TileBedTierCursor>,
}

#[cfg(uses_tile_bed_overlap_context)]
impl TileBedOverlapContext {
    /// Build the tile-local overlap context from tiered BED-like windows.
    ///
    /// `tile_assignment_envelope` must contain every query interval that can be submitted to this
    /// context for the tile. Broad candidate windows that contain the full range are stored as
    /// always-hit windows and are not scanned per fragment. Narrow candidates are scanned directly,
    /// because they are not expected to cover full tile assignment ranges in production tile runs.
    pub(crate) fn new(
        chrom_len: u64,
        chromosome_windows: &ChromosomeBedWindows,
        spans: &TileBedWindowSpans,
        tile_assignment_envelope: Interval<u64>,
    ) -> Result<Self> {
        let split = split_tile_bed_windows(
            chrom_len,
            chromosome_windows,
            spans,
            tile_assignment_envelope,
        )?;
        let tier_cursors = split
            .scanned_windows_by_size_tier
            .into_iter()
            .map(|windows| TileBedTierCursor { windows, wd_ptr: 0 })
            .collect();

        Ok(Self {
            chrom_len,
            always_hit_windows: split.always_hit_windows,
            tier_cursors,
        })
    }

    /// Find BED-like windows hit by one query interval.
    pub(crate) fn find_overlapping_windows(
        &mut self,
        query_interval: Interval<u64>,
        min_overlap_fraction: f64,
        look_back: u64,
    ) -> Result<Option<OverlappingWindows>> {
        if !(0.0..=1.0).contains(&min_overlap_fraction) {
            return Err(Error::OverlapFractionOutOfBounds {
                overlap_fraction: min_overlap_fraction,
            });
        }

        let mut overlaps = OverlappingWindows::new(query_interval);
        for window in &self.always_hit_windows {
            overlaps
                .windows
                .push(OverlappingWindow::new_with_output_idx(
                    window.all_windows_idx,
                    window.interval,
                    1.0,
                    Some(window.output_idx),
                )?);
        }

        for tier_cursor in &mut self.tier_cursors {
            append_bed_window_tier_overlaps(
                self.chrom_len,
                &mut tier_cursor.wd_ptr,
                tier_cursor.windows.as_slice(),
                query_interval,
                min_overlap_fraction,
                look_back,
                &mut overlaps,
            )?;
        }

        if overlaps.windows.is_empty() {
            Ok(None)
        } else {
            Ok(Some(overlaps))
        }
    }
}

#[cfg(uses_bed_window_tier_helpers)]
fn span_slice<'a>(
    items: &'a [BedWindowEntry],
    span: Option<TileWindowSpan>,
    split_name: &'static str,
) -> Result<&'a [BedWindowEntry]> {
    let Some(span) = span else {
        return Ok(&[]);
    };
    if span.first_idx > span.last_idx_exclusive || span.last_idx_exclusive > items.len() {
        return Err(Error::InvalidBedWindowSplitSpan {
            split_name,
            start: span.first_idx,
            end: span.last_idx_exclusive,
            len: items.len(),
        });
    }
    Ok(&items[span.first_idx..span.last_idx_exclusive])
}

#[cfg(uses_bed_window_tier_helpers)]
fn clamp_bed_window_to_chrom(
    window: BedWindowEntry,
    chrom_len: u64,
) -> Result<Option<Interval<u64>>> {
    let window_start = window.start();
    let window_end = window.end().min(chrom_len);
    if window_end <= window_start {
        return Ok(None);
    }
    Ok(Some(Interval::new(window_start, window_end)?))
}

#[cfg(uses_bed_window_tier_helpers)]
fn append_bed_window_tier_overlaps(
    chrom_len: u64,
    wd_ptr: &mut usize,
    windows: &[BedWindowEntry],
    query_interval: Interval<u64>,
    min_overlap_fraction: f64,
    look_back: u64,
    overlaps: &mut OverlappingWindows,
) -> Result<()> {
    while *wd_ptr < windows.len()
        && windows[*wd_ptr].end() <= query_interval.start().saturating_sub(look_back)
    {
        *wd_ptr += 1;
    }

    let mut scan_window_idx = *wd_ptr;
    while scan_window_idx < windows.len() && windows[scan_window_idx].start() < query_interval.end()
    {
        let window = windows[scan_window_idx];
        let Some(window_interval) = clamp_bed_window_to_chrom(window, chrom_len)? else {
            scan_window_idx += 1;
            continue;
        };
        if query_interval.intersects(window_interval) {
            let overlap_fraction = fraction_overlap_of_a(query_interval, window_interval);
            if overlap_fraction >= min_overlap_fraction {
                overlaps
                    .windows
                    .push(OverlappingWindow::new_with_output_idx(
                        window.all_windows_idx,
                        window_interval,
                        overlap_fraction,
                        Some(window.output_idx),
                    )?);
            }
        }
        scan_window_idx += 1;
    }

    Ok(())
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

/// Source of windows for repeated fixed-width overlap lookups.
///
/// The cursor owns this source so BED mode can keep cache state next to the tier it describes.
/// Global and fixed-size sources do not need cache state. BED sources are built from tile-local
/// BED window views, which keeps the ref-kmers hot loop from rescanning chromosome-wide background
/// windows for every k-mer start.
#[cfg(feature = "cmd_ref_kmers")]
#[derive(Debug)]
pub(crate) enum FixedWidthWindowSource {
    /// One chromosome-wide output row.
    ///
    /// Every non-empty query interval maps to row 0 with overlap fraction `1.0`.
    Global,
    /// Fixed-size chromosome-local bins.
    ///
    /// Bin indices are chromosome-local. Callers add the chromosome row offset when converting
    /// overlaps into output rows.
    FixedSize(u64),
    /// BED-like windows split into always-hit rows and independently cached size tiers.
    ///
    /// `always_hit_windows` contains broad tile-local BED windows that fully contain the tile
    /// assignment envelope, so every query submitted to this cursor intersects them. `tiers`
    /// contains the remaining tile-local BED windows split by size, with one forward pointer and
    /// cache per tier.
    Bed {
        /// Broad BED windows that contain every query interval assigned to this tile.
        always_hit_windows: Vec<BedWindowEntry>,
        /// Start-sorted scanned tiers, each with its own cache.
        tiers: Vec<FixedWidthBedTier>,
    },
}

#[cfg(feature = "cmd_ref_kmers")]
impl FixedWidthWindowSource {
    /// Build a BED source from the tile-local tier spans precomputed by the command runner.
    ///
    /// `tile_assignment_envelope` must contain every query interval that can be submitted to the
    /// returned source. Broad BED windows containing that whole range are moved to the always-hit
    /// list. Narrow windows and non-covering broad windows stay in their size tier and are scanned
    /// through that tier's independent cache. This avoids checking narrow windows for a full-range
    /// condition they are not expected to satisfy in production tile runs.
    #[cfg(uses_bed_window_tier_helpers)]
    pub(crate) fn bed_from_tile_view(
        chrom_len: u64,
        bed_window_view: TileBedWindowView<'_>,
        tile_assignment_envelope: Interval<u64>,
    ) -> Result<Self> {
        let split = split_tile_bed_windows(
            chrom_len,
            bed_window_view.chromosome_windows,
            bed_window_view.spans,
            tile_assignment_envelope,
        )?;
        Ok(Self::from_bed_split(split))
    }

    fn from_bed_split(split: TileBedWindowSplit) -> Self {
        Self::Bed {
            always_hit_windows: split.always_hit_windows,
            tiers: split
                .scanned_windows_by_size_tier
                .into_iter()
                .filter(|windows| !windows.is_empty())
                .map(FixedWidthBedTier::new)
                .collect(),
        }
    }
}

#[cfg(feature = "cmd_ref_kmers")]
#[derive(Debug)]
pub(crate) struct FixedWidthBedTier {
    /// Tile-local BED windows for this size tier, sorted by start coordinate.
    windows: Vec<BedWindowEntry>,
    /// Forward pointer and cached candidate windows for this tier.
    cache: FixedWidthBedCache,
}

#[cfg(feature = "cmd_ref_kmers")]
impl FixedWidthBedTier {
    fn new(windows: Vec<BedWindowEntry>) -> Self {
        Self {
            windows,
            cache: FixedWidthBedCache::new(),
        }
    }
}

#[cfg(feature = "cmd_ref_kmers")]
#[derive(Debug, Default)]
struct FixedWidthBedCache {
    /// Forward scan pointer into the tier's start-sorted BED windows.
    ///
    /// This is the first window not known to end before the current query start. It only moves
    /// forward, so BED-mode query starts must be requested in nondecreasing order.
    wd_ptr: usize,
    /// Whether `cached_windows` and `next_candidate_change_query_start` describe this tier.
    cache_ready: bool,
    /// Candidate BED windows for the current query-start range.
    ///
    /// Each cached window stores the query-start ranges where it passes the overlap threshold and
    /// where its overlap length is constant.
    cached_windows: Vec<CachedFixedWidthWindow>,
    /// First query start where the tier candidate slice may change.
    ///
    /// This is a query-start coordinate, not a genomic window coordinate. The cache is valid while
    /// `query_start < next_candidate_change_query_start`. At this boundary, a BED window may enter
    /// or leave the candidate slice, so this tier must rescan from `wd_ptr`.
    next_candidate_change_query_start: u64,
    /// Number of cache refreshes, used by tests to verify cache reuse.
    refresh_count: usize,
}

#[cfg(feature = "cmd_ref_kmers")]
impl FixedWidthBedCache {
    fn new() -> Self {
        Self {
            next_candidate_change_query_start: u64::MAX,
            ..Self::default()
        }
    }

    /// Mark the candidate cache stale without moving the forward pointer.
    ///
    /// Clipped queries use a direct tier scan so their fractions use the clipped query length as
    /// the denominator. The next full-width query rebuilds the cache from the current pointer.
    fn invalidate(&mut self) {
        self.cache_ready = false;
        self.cached_windows.clear();
        self.next_candidate_change_query_start = u64::MAX;
    }

    /// Return whether the current tier cache can answer a query beginning at `query_start`.
    ///
    /// `next_candidate_change_query_start` uses the same coordinate system as `query_start`. The
    /// query end does not take part in this check because the cache boundary has already been
    /// converted to the first query start where the candidate slice can change.
    fn cache_valid_for_query_start(&self, query_start: u64) -> bool {
        self.cache_ready && query_start < self.next_candidate_change_query_start
    }

    /// Rebuild cached BED candidates for one tier and the current query-start range.
    ///
    /// This advances the tier's forward pointer, collects windows that overlap the current
    /// full-width query, and precomputes the query-start ranges where each candidate passes the
    /// overlap threshold.
    fn refresh(
        &mut self,
        chrom_len: u64,
        windows: &[BedWindowEntry],
        query_interval: Interval<u64>,
        query_width: u64,
        min_overlap_fraction: f64,
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
        let required_overlap = required_overlap_bases(query_width, min_overlap_fraction);
        // Convert the next unseen BED start into the first query start that can reach it. If there
        // is no unseen window, the sentinel keeps this out of the minimum
        let mut next_candidate_change_query_start = first_query_start_reaching_window(
            windows
                .get(candidate_end_idx)
                .map(BedWindowEntry::start)
                .unwrap_or(u64::MAX),
            query_width,
        );

        for window in &windows[self.wd_ptr..candidate_end_idx] {
            // Existing candidates leave the slice when the query start reaches their clipped end.
            // The cache is valid only until the earliest candidate entry or exit query start
            next_candidate_change_query_start = next_candidate_change_query_start.min(
                later_query_start_or_never(window.end().min(chrom_len), query_interval.start()),
            );

            // Convert this BED window into query-start ranges that can be reused while the
            // candidate slice is unchanged. Windows too short for the threshold are excluded here
            if let Some(cached_window) = CachedFixedWidthWindow::new(
                window.all_windows_idx,
                Some(window.output_idx),
                window.start(),
                window.end(),
                chrom_len,
                query_width,
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

/// Find window overlaps for many same-width intervals on one chromosome.
///
/// Each lookup asks for the windows touched by `[query_start, query_start + query_width)`, clipped
/// to `chrom_len`. This is useful when scanning k-mer starts or other fixed-width intervals in
/// nondecreasing coordinate order.
///
/// Global and fixed-size sources are handled directly. BED mode gives the same overlap set as
/// [`find_overlapping_windows`] with `look_back = 0`, but the source is split before scanning:
/// broad windows covering the full tile assignment envelope are emitted from an always-hit list,
/// and the remaining windows keep independent caches for each size tier. That prevents broad
/// background windows from pinning the pointer for shorter nested windows. BED-mode callers must
/// provide query starts in nondecreasing order.
///
/// For full-width cached BED queries, a window passes when its overlapping bases meet the requested
/// fraction of `query_width`. A threshold of `0.0` means any actual overlap. Zero-overlap windows
/// are never returned. Queries clipped at the chromosome end scan the same tiers without using the
/// cache so their fractions are measured against the clipped query length. BED overlaps carry
/// `output_idx`, which is the original BED row id for ordinary BED input and the group index for
/// grouped BED input.
#[cfg(feature = "cmd_ref_kmers")]
#[derive(Debug)]
pub(crate) struct FixedWidthOverlapCursor {
    /// Chromosome length used to clip query intervals and BED window ends.
    chrom_len: u64,
    /// Window source for this cursor.
    ///
    /// The source owns any BED tier caches needed by the cursor. Global and fixed-size modes keep
    /// no mutable source state.
    source: FixedWidthWindowSource,
    /// Width of the unclipped query interval.
    ///
    /// Ref-kmers uses this as the k-mer assignment width. Caching is only used when the clipped
    /// query still has this full width.
    query_width: u64,
    /// Minimum one-way overlap fraction used to decide whether a window is returned.
    ///
    /// This is the fraction of the query covered by a window, not the fraction of the window
    /// covered by the query. Cached full-width BED queries compare `overlap_bases / query_width`.
    /// Clipped queries and fixed-size windows use the current query interval length, matching
    /// [`find_overlapping_windows`].
    min_overlap_fraction: f64,
}

#[cfg(feature = "cmd_ref_kmers")]
impl FixedWidthOverlapCursor {
    /// Create a cursor for fixed-width overlap queries.
    ///
    /// Parameters
    /// ----------
    /// - `chrom_len`:
    ///   Chromosome length used to clip query intervals and BED window ends.
    /// - `source`:
    ///   Global, fixed-size, or prepared tile-local BED window source.
    /// - `query_width`:
    ///   Width of each unclipped query interval. Ref-kmers passes the k-mer assignment width here.
    /// - `min_overlap_fraction`:
    ///   Minimum one-way overlap fraction for keeping a window. Full-width cached BED queries
    ///   measure this as `overlap_bases / query_width`. Clipped queries measure it against the
    ///   clipped query length. This is not the fraction of the window covered by the query. Must be
    ///   in `[0.0, 1.0]`.
    ///
    /// Returns
    /// -------
    /// - `out`:
    ///   Cursor initialized for repeated fixed-width overlap lookups.
    pub(crate) fn new(
        chrom_len: u64,
        source: FixedWidthWindowSource,
        query_width: u64,
        min_overlap_fraction: f64,
    ) -> Result<Self> {
        Interval::new(0, query_width)?;
        if !(0.0..=1.0).contains(&min_overlap_fraction) {
            return Err(Error::OverlapFractionOutOfBounds {
                overlap_fraction: min_overlap_fraction,
            });
        }
        if let FixedWidthWindowSource::FixedSize(bin_size) = &source {
            if *bin_size == 0 {
                return Err(Error::InvalidBinSize {
                    bin_size: *bin_size,
                });
            }
        }

        Ok(Self {
            chrom_len,
            source,
            query_width,
            min_overlap_fraction,
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

        match &mut self.source {
            FixedWidthWindowSource::Global => Self::find_global(self.chrom_len, query_interval),
            FixedWidthWindowSource::FixedSize(bin_size) => Self::find_by_size(
                self.chrom_len,
                self.min_overlap_fraction,
                query_interval,
                *bin_size,
            ),
            FixedWidthWindowSource::Bed {
                always_hit_windows,
                tiers,
            } => Self::find_bed(
                self.chrom_len,
                self.query_width,
                self.min_overlap_fraction,
                always_hit_windows,
                tiers,
                query_interval,
            ),
        }
    }

    /// Return how many times the BED cache has been rebuilt.
    ///
    /// This is test-only instrumentation for checking that adjacent fixed-width lookups reuse the
    /// cache.
    #[cfg(test)]
    pub(crate) fn refresh_count(&self) -> usize {
        match &self.source {
            FixedWidthWindowSource::Bed { tiers, .. } => {
                tiers.iter().map(|tier| tier.cache.refresh_count).sum()
            }
            FixedWidthWindowSource::Global | FixedWidthWindowSource::FixedSize(_) => 0,
        }
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
        chrom_len: u64,
        min_overlap_fraction: f64,
        query_interval: Interval<u64>,
        bin_size: u64,
    ) -> Result<Option<OverlappingWindows>> {
        let mut overlaps = OverlappingWindows::new(query_interval);

        for bin_idx in create_overlapping_bins_by_size(query_interval, bin_size)? {
            let window_start = bin_idx * bin_size;
            let window_end = (bin_idx * bin_size + bin_size).min(chrom_len);
            if window_end <= window_start {
                continue;
            }
            let window_interval = Interval::new(window_start, window_end)?;
            let overlap_fraction = fraction_overlap_of_a(query_interval, window_interval);
            if overlap_fraction < min_overlap_fraction {
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
    /// Full-width queries use cached candidate windows. Clipped queries scan the tier windows
    /// directly so their overlap fraction uses the clipped query length as the denominator.
    fn find_bed(
        chrom_len: u64,
        query_width: u64,
        min_overlap_fraction: f64,
        always_hit_windows: &[BedWindowEntry],
        tiers: &mut [FixedWidthBedTier],
        query_interval: Interval<u64>,
    ) -> Result<Option<OverlappingWindows>> {
        let mut overlaps = OverlappingWindows::new(query_interval);
        Self::append_always_hit_windows(
            always_hit_windows,
            query_interval,
            min_overlap_fraction,
            &mut overlaps,
        )?;

        if query_interval.len() != query_width {
            for tier in tiers {
                tier.cache.invalidate();
                append_bed_window_tier_overlaps(
                    chrom_len,
                    &mut tier.cache.wd_ptr,
                    tier.windows.as_slice(),
                    query_interval,
                    min_overlap_fraction,
                    0,
                    &mut overlaps,
                )?;
            }
            return if overlaps.windows.is_empty() {
                Ok(None)
            } else {
                Ok(Some(overlaps))
            };
        }

        for tier in tiers {
            if !tier
                .cache
                .cache_valid_for_query_start(query_interval.start())
            {
                tier.cache.refresh(
                    chrom_len,
                    tier.windows.as_slice(),
                    query_interval,
                    query_width,
                    min_overlap_fraction,
                )?;
            }

            for cached_window in &tier.cache.cached_windows {
                let Some(overlap_bases) = cached_window.overlap_bases_at(query_interval.start())
                else {
                    continue;
                };

                let overlap_fraction = overlap_bases as f64 / query_width as f64;
                overlaps
                    .windows
                    .push(OverlappingWindow::new_with_output_idx(
                        cached_window.idx,
                        cached_window.interval,
                        overlap_fraction,
                        cached_window.output_idx,
                    )?);
            }
        }

        if overlaps.windows.is_empty() {
            Ok(None)
        } else {
            Ok(Some(overlaps))
        }
    }

    fn append_always_hit_windows(
        always_hit_windows: &[BedWindowEntry],
        query_interval: Interval<u64>,
        min_overlap_fraction: f64,
        overlaps: &mut OverlappingWindows,
    ) -> Result<()> {
        for window in always_hit_windows {
            if !query_interval.intersects(window.interval) {
                continue;
            }
            let overlap_fraction = fraction_overlap_of_a(query_interval, window.interval);
            if overlap_fraction < min_overlap_fraction {
                continue;
            }
            overlaps
                .windows
                .push(OverlappingWindow::new_with_output_idx(
                    window.all_windows_idx,
                    window.interval,
                    overlap_fraction,
                    Some(window.output_idx),
                )?);
        }
        Ok(())
    }

    /// Return the chromosome-wide window hit for global mode.
    ///
    /// Global mode has a single row for the whole chromosome, so every non-empty query interval
    /// contributes with overlap fraction `1.0`.
    fn find_global(
        chrom_len: u64,
        query_interval: Interval<u64>,
    ) -> Result<Option<OverlappingWindows>> {
        let mut overlaps = OverlappingWindows::new(query_interval);
        overlaps.windows.push(OverlappingWindow::new(
            0,
            Interval::new(0, chrom_len)?,
            1.0,
        )?);
        Ok(Some(overlaps))
    }
}

/// Cached query-start ranges for one candidate BED window and one fixed query width.
///
/// The accepted range controls whether this window can be returned for a query start. The
/// constant-overlap range records the starts where the overlap length is unchanged. Outside that
/// range, moving the query by one base changes the overlap by one base.
#[cfg(feature = "cmd_ref_kmers")]
#[derive(Debug, Clone, Copy)]
struct CachedFixedWidthWindow {
    idx: usize,
    output_idx: Option<u64>,
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

#[cfg(feature = "cmd_ref_kmers")]
impl CachedFixedWidthWindow {
    /// Precompute query-start ranges for a candidate BED window.
    ///
    /// The cached window is returned only when it can reach `required_overlap` bases with a
    /// full-width query. Empty or chromosome-clipped-away windows return `None`.
    fn new(
        idx: usize,
        output_idx: Option<u64>,
        window_start: u64,
        window_end: u64,
        chrom_len: u64,
        query_width: u64,
        required_overlap: u64,
    ) -> Result<Option<Self>> {
        let window_end = window_end.min(chrom_len);
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
            output_idx,
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
#[cfg(feature = "cmd_ref_kmers")]
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
#[cfg(feature = "cmd_ref_kmers")]
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
#[cfg(feature = "cmd_ref_kmers")]
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
