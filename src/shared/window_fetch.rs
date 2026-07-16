use anyhow::Result;

#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_lengths"
))]
use crate::commands::cli_common::WindowSpec;
#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_lengths"
))]
use crate::shared::tiled_run::clamp_fetch_to_window_span;
#[cfg(any(
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_midpoints"
))]
use crate::shared::tiled_run::overlapping_windows_for_tile;
use crate::shared::{
    interval::{IndexedInterval, Interval},
    tiled_run::{Tile, TileWindowSpan},
};

/// Helpers for turning already-selected windows into an aligned BAM fetch interval.
///
/// This module separates three different questions that used to get mixed together:
///
/// 1. Which BED windows are relevant for a tile?
/// 2. Can those relevant BED windows be converted into an aligned-coordinate window extent?
/// 3. After that extent is known, how should it be clamped onto `tile.fetch`?
///
/// Here "window extent" means the smallest aligned interval running from the minimum relevant
/// window start to the maximum relevant window end before any extra halo is applied.
///
/// The second question matters because BAM fetching happens in aligned reference coordinates,
/// while command-level window relevance may be decided in a different coordinate model. Some
/// commands can safely say "these BED starts and ends bound every aligned read I need, plus halo".
/// Other commands select windows from shifted, clipped, midpoint, or fragment-reach logic where
/// BED coordinates alone are not enough to prove a smaller aligned fetch is safe.
///
/// The supported BED fetch policies are:
///
/// - `CoreOverlap`:
///   BED relevance is defined only by overlap between the BED interval and the tile core. This is
///   the right model when the command is asking "which windows physically intersect this tile's
///   owned reference-coordinate region?", without using fragment start, endpoint, midpoint, or
///   clipped-boundary reach to make windows relevant.
/// - `CandidateWindowExtent`:
///   BED relevance has already been decided elsewhere, and BAM fetch narrowing should use the
///   min/max window extent of that candidate window set. The caller must supply a halo that is large
///   enough to preserve every aligned read that could make those windows relevant.
/// - `KeepTileFetch`:
///   BED windows may still be relevant for counting, but the BED coordinates are not safe inputs
///   to aligned BAM fetch narrowing. In that case the caller keeps the full aligned `tile.fetch`
///   band.
#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_lengths"
))]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum BedFetchPolicy {
    #[cfg(any(feature = "cmd_fcoverage", feature = "cmd_fragment_kmers"))]
    CoreOverlap,
    #[cfg(any(feature = "cmd_ends", feature = "cmd_lengths"))]
    CandidateWindowExtent,
    #[allow(dead_code)]
    KeepTileFetch,
}

/// Determine the fetch span for a tile based on the active window strategy.
///
/// Global mode uses the full tile fetch range. Fixed-size window mode narrows
/// the fetch span to the first and last windows touching the tile core. BED
/// mode uses the precomputed window bounds for the tile and returns `None`
/// when the tile does not intersect any BED windows at all.
///
/// Parameters
/// ----------
/// - `tile`:
///   Tile describing the chromosome, core span, and fetch span
/// - `tile_window_span`:
///   Cached min and max window bounds for the tile in BED mode
/// - `windows_chr`:
///   Chromosome BED windows as `(start, end, idx)` tuples in BED mode
/// - `window_opt`:
///   Window specification selecting global, fixed-size, or BED mode
/// - `chrom_len`:
///   Chromosome length used to clamp fetch coordinates
/// - `halo_bp`:
///   Extra bases to keep on both sides of the active window extent before
///   clamping back onto the tile fetch interval
/// - `bed_fetch_policy`:
///   Explicit BED-mode policy describing whether BED windows narrow fetch by tile/core overlap,
///   by aligned fragment reach, or not at all
///
/// Returns
/// -------
/// - `Result<Option<Interval<u64>>>`:
///   Checked absolute fetch interval, or `None` when no windows apply
#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_lengths"
))]
pub(crate) fn fetch_span_for_tile(
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    windows_chr: Option<&[IndexedInterval<u64>]>,
    window_opt: &WindowSpec,
    chrom_len: u64,
    halo_bp: u64,
    bed_fetch_policy: BedFetchPolicy,
) -> Result<Option<Interval<u64>>> {
    match window_opt {
        WindowSpec::Global => {
            let fetch_start = tile.fetch_start() as u64;
            let fetch_end = (tile.fetch_end().min(chrom_len as u32)) as u64;
            if fetch_start >= fetch_end {
                return Ok(None);
            }
            Ok(Some(Interval::new(fetch_start, fetch_end)?))
        }
        WindowSpec::Size(window_bp) => {
            let core_start = tile.core_start() as u64;
            let core_end = (tile.core_end() as u64).min(chrom_len);
            if core_start >= chrom_len || core_end == 0 {
                return Ok(None);
            }
            let window_idx_start = core_start / window_bp;
            let window_idx_end = (core_end.saturating_sub(1)) / window_bp;
            let window_start = window_idx_start * window_bp;
            let window_end = ((window_idx_end + 1) * window_bp).min(chrom_len);
            let window_span = Interval::new(window_start, window_end)?;
            Ok(clamp_fetch_to_window_span(
                tile,
                chrom_len,
                window_span,
                halo_bp,
            )?)
        }
        WindowSpec::Bed(_) => {
            let Some(windows_chr) = windows_chr else {
                return Ok(None);
            };
            match bed_fetch_policy {
                #[cfg(any(feature = "cmd_fcoverage", feature = "cmd_fragment_kmers"))]
                BedFetchPolicy::CoreOverlap => {
                    let Some(window_span) = window_derived_fetch_extent_for_core_overlap(
                        windows_chr,
                        tile,
                        tile_window_span,
                    )?
                    else {
                        return Ok(None);
                    };
                    Ok(clamp_fetch_to_window_span(
                        tile,
                        chrom_len,
                        window_span,
                        halo_bp,
                    )?)
                }
                #[cfg(any(feature = "cmd_ends", feature = "cmd_lengths"))]
                BedFetchPolicy::CandidateWindowExtent => {
                    let Some(window_span) =
                        window_derived_fetch_extent_for_candidates(windows_chr, tile_window_span)?
                    else {
                        return Ok(None);
                    };
                    Ok(clamp_fetch_to_window_span(
                        tile,
                        chrom_len,
                        window_span,
                        halo_bp,
                    )?)
                }
                BedFetchPolicy::KeepTileFetch => full_tile_fetch_span(tile, chrom_len),
            }
        }
    }
}

/// Determine the fetch span for an already-selected BED candidate set.
///
/// This is the BED-only equivalent of `fetch_span_for_tile(..., CandidateWindowExtent)`, but it
/// accepts any window record that exposes a checked `Interval<u64>`.
#[cfg(any(feature = "cmd_ends", feature = "cmd_lengths"))]
pub(crate) fn fetch_span_for_bed_candidates<T>(
    tile: &Tile,
    candidate_span: Option<&TileWindowSpan>,
    windows_chr: &[T],
    chrom_len: u64,
    halo_bp: u64,
) -> Result<Option<Interval<u64>>>
where
    T: AsRef<Interval<u64>>,
{
    let Some(window_span) =
        window_derived_fetch_extent_for_candidates(windows_chr, candidate_span)?
    else {
        return Ok(None);
    };
    Ok(clamp_fetch_to_window_span(
        tile,
        chrom_len,
        window_span,
        halo_bp,
    )?)
}

/// Derive the aligned min/max window extent from BED windows that truly overlap the tile core.
///
/// The helper iterates over the core-overlapping windows, using the cached candidate span when
/// provided, and tracks the minimum start and maximum end across those true overlaps. When no BED
/// window intersects the tile core, it returns `None`.
///
/// Coordinate space:
/// - consumes BED window coordinates
/// - returns an aligned-coordinate window extent
///
/// Fragment ownership rule:
/// - none; this helper is for commands whose BED relevance is defined by core overlap
///
/// Counting or assignment interval assumption:
/// - none beyond core-overlap BED ownership
///
/// Aligned fetch narrowing:
/// - allowed
///
/// This helper must derive the leftmost and rightmost bounds from windows that actually overlap
/// the tile core, even when a wider cached candidate span is available.
#[cfg(any(
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_midpoints"
))]
pub(crate) fn window_derived_fetch_extent_for_core_overlap(
    windows_chr: &[IndexedInterval<u64>],
    tile: &Tile,
    candidate_span: Option<&TileWindowSpan>,
) -> Result<Option<Interval<u64>>> {
    let mut iter = overlapping_windows_for_tile(windows_chr, tile, candidate_span);
    let Some(first) = iter.next() else {
        return Ok(None);
    };

    let mut min_start = first.start();
    let mut max_end = first.end();
    for window in iter {
        min_start = min_start.min(window.start());
        max_end = max_end.max(window.end());
    }

    Ok(Some(Interval::new(min_start, max_end)?))
}

/// Derive the aligned min/max window extent from an already-selected BED candidate set.
///
/// Coordinate space:
/// - consumes BED window coordinates
/// - returns an aligned-coordinate window extent
///
/// Fragment ownership rule:
/// - not interpreted here; the caller has already selected the candidate windows
///
/// Counting or assignment interval assumption:
/// - not interpreted here; the caller must prove that this BED extent plus the supplied halo is a
///   safe aligned fetch interval for the active counting coordinate choice
///
/// Aligned fetch narrowing:
/// - allowed only because the caller has already proven that the candidate-window extent plus the
///   supplied halo preserves every aligned read needed for counting
///
/// This helper must derive the min/max directly from the candidate span and must not reapply a
/// core-overlap filter.
#[cfg(any(feature = "cmd_ends", feature = "cmd_lengths"))]
pub(crate) fn window_derived_fetch_extent_for_candidates<T>(
    windows_chr: &[T],
    candidate_span: Option<&TileWindowSpan>,
) -> Result<Option<Interval<u64>>>
where
    T: AsRef<Interval<u64>>,
{
    let Some(candidate_span) = candidate_span else {
        return Ok(None);
    };
    if candidate_span.is_empty() {
        return Ok(None);
    }

    let start_idx = candidate_span.first_idx.min(windows_chr.len());
    let end_idx = candidate_span.last_idx_exclusive.min(windows_chr.len());
    if start_idx >= end_idx {
        return Ok(None);
    }

    let first_window = windows_chr[start_idx].as_ref();
    let mut min_start = first_window.start();
    let mut max_end = first_window.end();
    for window in &windows_chr[start_idx + 1..end_idx] {
        let window_interval = window.as_ref();
        min_start = min_start.min(window_interval.start());
        max_end = max_end.max(window_interval.end());
    }

    Ok(Some(Interval::new(min_start, max_end)?))
}

/// Return the full aligned BAM fetch span already stored on the tile.
///
/// Coordinate space:
/// - consumes the tile's aligned fetch band
/// - returns a final aligned BAM fetch span
///
/// Fragment ownership rule:
/// - not interpreted here; the caller has already decided that BED windows must not narrow fetch
///
/// Counting or assignment interval assumption:
/// - BED relevance may use different coordinates from aligned BAM fetch
///
/// Aligned fetch narrowing:
/// - not performed here
///
/// This is the generic "do not narrow fetch from BED windows" policy.
#[cfg(any(
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_lengths"
))]
pub(crate) fn full_tile_fetch_span(tile: &Tile, chrom_len: u64) -> Result<Option<Interval<u64>>> {
    let fetch_start = tile.fetch_start() as u64;
    let fetch_end = (tile.fetch_end().min(chrom_len as u32)) as u64;
    if fetch_start >= fetch_end {
        return Ok(None);
    }
    Ok(Some(Interval::new(fetch_start, fetch_end)?))
}

#[cfg(test)]
mod tests {
    include!("window_fetch_tests.rs");
}
