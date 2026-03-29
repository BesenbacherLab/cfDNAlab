use anyhow::Result;

use crate::{
    commands::cli_common::WindowSpec,
    shared::{
        interval::{IndexedInterval, Interval},
        tiled_run::{
            Tile, TileWindowSpan, clamp_fetch_to_window_span, tile_window_min_max,
        },
    },
};

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
///   Extra bases to keep on both sides of the active window span before
///   clamping back onto the tile fetch interval
///
/// Returns
/// -------
/// - `Result<Option<Interval<u64>>>`:
///   Checked absolute fetch interval, or `None` when no windows apply
pub fn fetch_span_for_tile(
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    windows_chr: Option<&[IndexedInterval<u64>]>,
    window_opt: &WindowSpec,
    chrom_len: u64,
    halo_bp: u64,
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
            let Some(window_span) = tile_window_min_max(windows_chr, tile, tile_window_span)? else {
                return Ok(None);
            };
            Ok(clamp_fetch_to_window_span(
                tile,
                chrom_len,
                window_span,
                halo_bp,
            )?)
        }
    }
}
