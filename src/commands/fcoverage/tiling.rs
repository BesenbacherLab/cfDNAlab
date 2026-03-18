use std::io::{BufWriter, Write};

use anyhow::{Context, Result};

use crate::{
    commands::fcoverage::window_results::CoverageWindowAction,
    shared::interval::Interval,
    shared::tiled_run::{
        Tile, TileMode, TileWindowSpan, clamp_fetch_to_window_span, parse_tile_index,
        tile_window_min_max,
    },
};

/// Concatenates per-tile positional files into a single merged output.
///
/// The routine scans each chromosome, orders the matching tile files by index, and stream-copies
/// their contents into the destination writer without re-encoding, allowing pre-compressed chunks
/// to remain untouched.
///
/// # Parameters
/// - `temp_dir`: Directory containing the per-tile files.
/// - `out_dir`: Directory where the merged file should be written.
/// - `chromosomes`: Chromosome names that determine merge order.
/// - `per_tile_prefix`: Prefix used in the per-tile file names.
/// - `final_name`: File name for the merged output.
///
/// # Returns
/// Path to the merged file on success.
pub fn merge_positional_tiles(
    temp_dir: &std::path::Path,
    out_dir: &std::path::Path,
    chromosomes: &[String],
    per_tile_prefix: &str, // e.g. "coverage.pos" (whole-genome) or "coverage.pos.win" (windowed)
    final_name: &str,      // e.g. "coverage.per_position.tsv"
) -> Result<std::path::PathBuf> {
    let final_path = out_dir.join(final_name);
    let mut out = BufWriter::new(
        std::fs::File::create(&final_path)
            .with_context(|| format!("Creating merged output: {}", final_path.display()))?,
    );

    for chr in chromosomes {
        // Collect tile files for this chromosome from temp_dir
        let mut chr_files: Vec<(u32, std::path::PathBuf)> = Vec::new();
        for entry in std::fs::read_dir(temp_dir)
            .with_context(|| format!("Listing temp_dir: {}", temp_dir.display()))?
        {
            let path = entry?.path();
            if !path.is_file() {
                continue;
            }
            let fname = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
            // Expect "{per_tile_prefix}.{chr}.{index}.tsv"
            if fname.starts_with(per_tile_prefix)
                && fname.contains(&format!(".{chr}."))
                && let Some(idx) = parse_tile_index(fname)
            {
                chr_files.push((idx, path));
            }
        }

        // Sort by tile index to preserve genomic order within chr
        chr_files.sort_by_key(|(i, _)| *i);

        // Stream copy each tile into the final file
        for (_idx, path) in chr_files {
            let mut f = std::fs::File::open(&path)
                .with_context(|| format!("Opening tile file: {}", path.display()))?;
            std::io::copy(&mut f, &mut out).with_context(|| {
                format!(
                    "Copying from {} into {}",
                    path.display(),
                    final_path.display()
                )
            })?;
        }
    }

    out.flush().context("Flushing merged output")?;
    Ok(final_path)
}

/// Joins already-compressed per-tile final outputs while preserving frame boundaries.
///
/// A compressed header frame is written first, followed by each tile frame in genomic order, so
/// the resulting file stays a valid zstd concatenation stream suitable for downstream tools.
///
/// # Parameters
/// - `temp_dir`: Directory containing the compressed per-tile final files.
/// - `out_dir`: Directory where the merged file will be placed.
/// - `chromosomes`: Chromosome names that dictate processing order.
/// - `per_tile_prefix`: Prefix shared by the per-tile files.
/// - `final_name`: File name of the merged artifact.
/// - `header_line`: Plain-text header to encode as its own compressed frame.
///
/// # Returns
/// Path to the merged file on success.
pub fn concat_aligned_size_tile_finals(
    temp_dir: &std::path::Path,
    out_dir: &std::path::Path,
    chromosomes: &[String],
    per_tile_prefix: &str, // e.g., "<prefix>.fin"
    final_name: &str,      // e.g., "<prefix>.avg.tsv.zst"
    header_line: &str,     // single header line without trailing newline
) -> Result<std::path::PathBuf> {
    let final_path = out_dir.join(final_name);
    let mut out = BufWriter::new(
        std::fs::File::create(&final_path)
            .with_context(|| format!("Creating {}", final_path.display()))?,
    );

    // Write a compressed header frame first (so we never touch tile frames).
    let mut header_bytes = header_line.as_bytes().to_vec();
    header_bytes.push(b'\n');
    let header_frame =
        zstd::encode_all(&header_bytes[..], 3).context("Compressing header frame")?;
    out.write_all(&header_frame)?;

    // Then append each tile's compressed frame in genomic order
    for chr in chromosomes {
        // Collect tile files for this chromosome
        let mut chr_files: Vec<(u32, std::path::PathBuf)> = Vec::new();
        for entry in std::fs::read_dir(temp_dir)
            .with_context(|| format!("Listing {}", temp_dir.display()))?
        {
            let path = entry?.path();
            if !path.is_file() {
                continue;
            }
            let fname = path.file_name().and_then(|s| s.to_str()).unwrap_or("");
            if fname.starts_with(per_tile_prefix)
                && fname.contains(&format!(".{chr}."))
                && let Some(idx) = parse_tile_index(fname)
            {
                chr_files.push((idx, path));
            }
        }
        chr_files.sort_by_key(|(i, _)| *i);

        // Copy bytes verbatim (frame concatenation)
        for (_i, p) in chr_files {
            let mut f =
                std::fs::File::open(&p).with_context(|| format!("Opening {}", p.display()))?;
            std::io::copy(&mut f, &mut out).with_context(|| {
                format!("Copying {} into {}", p.display(), final_path.display())
            })?;
        }
    }

    out.flush()?;
    Ok(final_path)
}

/// Shrinks a tile's fetch region to the range implied by the overlapping windows.
///
/// Depending on the tile mode, the function either keeps the original fetch span or intersects it
/// with the min/max window bounds (expanded by halos) so that downstream fetches only read the
/// necessary bases.
///
/// # Parameters
/// - `tile`: Tile whose fetch interval may be reduced.
/// - `tile_span`: Optional cached window span for the tile.
/// - `mode`: Output mode describing whether windows are used.
/// - `chrom_len`: Length of the chromosome in bases.
/// - `halo_bp`: Extra bases to keep on both sides of the overlapping window span so
///   fragments that extend outside the window itself can still be reconstructed.
///
/// # Returns
/// Checked absolute fetch interval, or `None` when no fetch is needed.
pub fn adapt_fetch_to_extreme_windows(
    tile: &Tile,
    tile_span: Option<&TileWindowSpan>,
    mode: &TileMode<'_>,
    chrom_len: u32,
    halo_bp: u64,
) -> Result<Option<Interval<u64>>> {
    let chrom_len_u64 = chrom_len as u64;

    // Decide the fetch interval based on mode/windows.
    // For whole-genome positional: use the full tile fetch band.
    // For windowed runs: restrict to the overlapping window span widened by a fragment-sized halo,
    // then intersect it with the tile’s existing fetch band.
    match mode {
        TileMode::Positional { windows: None, .. } => Ok(Some(Interval::new(
            tile.fetch_start() as u64,
            tile.fetch_end() as u64,
        )?)),
        TileMode::Positional {
            windows: Some(wchr),
            ..
        } => {
            let Some(window_span) = tile_window_min_max(wchr, tile, tile_span)? else {
                return Ok(None);
            };
            Ok(clamp_fetch_to_window_span(
                tile,
                chrom_len_u64,
                window_span,
                halo_bp,
            )?)
        }
        TileMode::AggregatesByBed { windows: wchr, .. } => {
            let Some(window_span) = tile_window_min_max(wchr, tile, tile_span)? else {
                return Ok(None);
            };
            Ok(clamp_fetch_to_window_span(
                tile,
                chrom_len_u64,
                window_span,
                halo_bp,
            )?)
        }
        TileMode::AggregatesBySize { window_bp, .. } => {
            let core_start = tile.core_start() as u64;
            let core_end = tile.core_end() as u64;
            if core_start >= chrom_len_u64 {
                return Ok(None);
            }
            let window_size_bp = *window_bp;
            let first_window_idx = core_start / window_size_bp;
            let last_window_idx = (core_end.saturating_sub(1)) / window_size_bp;
            let min_window_start = first_window_idx * window_size_bp;
            let max_window_end = ((last_window_idx + 1) * window_size_bp).min(chrom_len_u64);
            let window_span = Interval::new(min_window_start, max_window_end)?;
            Ok(clamp_fetch_to_window_span(
                tile,
                chrom_len_u64,
                window_span,
                halo_bp,
            )?)
        }
    }
}

/// Computes the aggregate coverage statistics over a tile-local range.
///
/// The function pulls values from prefix sums whenever available and only scans the mask slice when
/// necessary, yielding both the summed coverage and the count of allowed versus blacklisted bases.
///
/// # Parameters
/// - `local_start_idx`: Inclusive start index inside the tile-local arrays.
/// - `local_end_idx`: Exclusive end index inside the tile-local arrays.
/// - `masked`: Whether masked mode is enabled.
/// - `ps_all`: Prefix sums over all bases.
/// - `ps_allow`: Optional prefix sums over allowed bases.
/// - `cnt_allow`: Optional prefix sums over the count of allowed bases.
/// - `mask`: Optional mask where `1` marks blacklisted bases (used when `cnt_allow` is absent).
///
/// # Returns
/// A triple `(sum, allowed_bases, blacklisted_bases)` for the requested span.
///
/// # Panics
/// The caller must ensure the indices are within bounds of the prefix sum arrays.
#[inline]
pub fn coverage_sum_and_counts(
    local_start_idx: usize,
    local_end_idx: usize,
    masked: bool,
    ps_all: &[f64],
    ps_allow: Option<&[f64]>,
    cnt_allow: Option<&[u32]>,
    mask: Option<&[u8]>,
) -> (f64, u64, u64) {
    let sum = if masked {
        if let Some(pa) = ps_allow {
            pa[local_end_idx] - pa[local_start_idx]
        } else {
            ps_all[local_end_idx] - ps_all[local_start_idx]
        }
    } else {
        ps_all[local_end_idx] - ps_all[local_start_idx]
    };

    let span = (local_end_idx - local_start_idx) as u64;

    let allowed = if masked {
        if let Some(cnt) = cnt_allow {
            (cnt[local_end_idx] - cnt[local_start_idx]) as u64
        } else if let Some(m) = mask {
            let mut ok = 0u64;
            for i in local_start_idx..local_end_idx {
                if m[i] == 0 {
                    ok += 1;
                }
            }
            ok
        } else {
            span
        }
    } else {
        span
    };

    let blacklisted = span - allowed;
    (sum, allowed, blacklisted)
}

/// Converts accumulated coverage statistics into a final window value.
///
/// Depending on the requested action the result is either an average (over allowed or full span)
/// or the raw total; zero denominators yield zero to avoid NaNs.
///
/// # Parameters
/// - `sum`: Accumulated coverage sum.
/// - `allowed_positions`: Number of unmasked bases in the window.
/// - `unmasked_span_bp`: Full span length in bases for unmasked mode.
/// - `masked`: Whether the masked mode is active.
/// - `mode`: Window action describing how to interpret the aggregates.
///
/// # Returns
/// The final window value after applying the requested action.
#[inline]
pub fn finalize_value(
    sum: f64,
    allowed_positions: u64,
    unmasked_span_bp: u64, // end-start when unmasked mode; ignored when masked mode
    masked: bool,
    mode: &CoverageWindowAction,
) -> f64 {
    match mode {
        CoverageWindowAction::Average => {
            if masked {
                if allowed_positions == 0 {
                    0.0
                } else {
                    sum / allowed_positions as f64
                }
            } else {
                if unmasked_span_bp == 0 {
                    0.0
                } else {
                    sum / unmasked_span_bp as f64
                }
            }
        }
        CoverageWindowAction::Total => sum,
        _ => unreachable!(),
    }
}

/// Absolute/local overlap between a requested interval and a tile core.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoreClip {
    /// Inclusive start index in tile-core local coordinates.
    pub local_start_idx: usize,
    /// Exclusive end index in tile-core local coordinates.
    pub local_end_idx: usize,
    /// Absolute interval after clipping the requested interval to the tile core.
    pub clipped_abs_interval: Interval<u64>,
}

/// Clips an absolute interval to the tile core and converts the overlap to core-local indices.
///
/// The helper returns both tile-core local indices and the clipped absolute interval so callers
/// can reuse whichever representation is needed.
///
/// # Parameters
/// - `absolute_interval`: Absolute interval to clip against the tile core.
/// - `core_interval`: Absolute tile-core interval used as the clipping target.
///
/// # Returns
/// Checked clipped overlap when the intervals intersect, otherwise `None`.
#[inline]
pub fn clip_interval_to_core_and_localize(
    absolute_interval: Interval<u64>,
    core_interval: Interval<u64>,
) -> Result<Option<CoreClip>> {
    let Some(clipped_abs_interval) = absolute_interval.intersection(core_interval) else {
        return Ok(None);
    };
    let local_start_idx = (clipped_abs_interval.start() - core_interval.start()) as usize;
    let local_end_idx = (clipped_abs_interval.end() - core_interval.start()) as usize;
    Ok(Some(CoreClip {
        local_start_idx,
        local_end_idx,
        clipped_abs_interval,
    }))
}
