use crate::utils::coverage::coverage_prefix::Coverage;
use crate::utils::fragment::minimal_fragment::Fragment;
use crate::utils::fragment::segment_fragment::FragmentWithSegments;
use crate::utils::{bam::Contigs, coverage::window_results::CoverageWindowAction};
use anyhow::{Context, Result};
use rand::{Rng, distr::Alphanumeric};
use std::fmt;
use std::io::{BufWriter, Write};

/// A processing tile for one chromosome
#[derive(Debug, Clone)]
pub struct Tile {
    pub chr: String,
    pub tid: i32,
    pub index: u32,       // 0-based index within chromosome
    pub core_start: u32,  // Inclusive
    pub core_end: u32,    // Exclusive
    pub fetch_start: u32, // Inclusive (core expanded by halo)
    pub fetch_end: u32,   // Exclusive
}

/// Half-open window indices covering a tile core.
///
/// The range `[first_idx, last_idx_exclusive)` selects the portion of the chromosome-specific
/// window slice whose starts fall before the tile core end. Windows that end before the tile core
/// start are excluded when the span is constructed, so streaming from `first_idx` is safe.
#[derive(Clone, Copy, Debug, Default)]
pub struct TileWindowSpan {
    pub first_idx: usize,
    pub last_idx_exclusive: usize,
}

impl TileWindowSpan {
    /// Reports whether the cached window span is empty.
    ///
    /// The span contains no usable windows when the start and end indices are identical,
    /// meaning every candidate was pruned beforehand.
    ///
    /// # Returns
    /// `true` when `first_idx == last_idx_exclusive`, otherwise `false`.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.first_idx == self.last_idx_exclusive
    }
}

/// Prepares the per-tile window spans that downstream code can reuse directly.
///
/// The routine advances twin pointers across the start-sorted window list of each chromosome as it
/// iterates through tiles, trimming windows that end too early and extending the range until starts
/// escape the tile core. This preserves streaming behaviour and avoids re-scanning the same
/// windows for neighbouring tiles.
///
/// # Parameters
/// - `tiles`: All tiles sorted by chromosome and core position.
/// - `windows_for_chr`: Closure that retrieves the start-sorted windows `(start, end, idx)` for the
///   requested chromosome.
///
/// # Returns
/// A vector the same length as `tiles` containing the optional `[first, last)` window index span
/// for each tile.
pub fn precompute_tile_window_spans<'a, F>(
    tiles: &[Tile],
    mut windows_for_chr: F,
) -> Vec<Option<TileWindowSpan>>
where
    F: FnMut(&str) -> &'a [(u64, u64, u64)],
{
    let mut spans: Vec<Option<TileWindowSpan>> = vec![None; tiles.len()];
    let mut tile_idx = 0usize;

    while tile_idx < tiles.len() {
        let chr = tiles[tile_idx].chr.as_str();
        let windows = windows_for_chr(chr);

        // Capture the range of tiles that share this chromosome so we can reuse the same
        // streaming window pointers across them
        let chr_tile_start = tile_idx;
        while tile_idx < tiles.len() && tiles[tile_idx].chr == chr {
            tile_idx += 1;
        }
        let chr_tile_end = tile_idx;

        if windows.is_empty() {
            continue;
        }

        let windows_len = windows.len();
        let mut w_left = 0usize;
        let mut w_right = 0usize;

        for idx in chr_tile_start..chr_tile_end {
            let tile = &tiles[idx];
            let core_start = tile.core_start as u64;
            let core_end = tile.core_end as u64;

            // Discard windows that end before the tile core
            while w_left < windows_len && windows[w_left].1 <= core_start {
                w_left += 1;
            }

            if w_right < w_left {
                w_right = w_left;
            }

            // Extend to cover every window whose start lies inside the tile core span
            while w_right < windows_len && windows[w_right].0 < core_end {
                w_right += 1;
            }

            if w_left == windows_len {
                // No windows remain for this chromosome, so the rest of the tiles cannot overlap
                spans[idx] = None;
                for span in spans[idx + 1..chr_tile_end].iter_mut() {
                    *span = None;
                }
                break;
            }

            // `w_left` now points at the first surviving window candidate for this tile
            let first_candidate = &windows[w_left];

            // Check if the earliest remaining window begins at/after the core end, so later ones do too
            if first_candidate.0 >= core_end {
                spans[idx] = None;
                continue;
            }

            // The iterator helpers perform the precise overlap check (end > core_start &&
            // start < core_end) so the recorded range only needs to bound the candidate windows
            spans[idx] = Some(TileWindowSpan {
                first_idx: w_left,
                last_idx_exclusive: w_right,
            });
        }
    }

    spans
}

/// Iterator over the windows that overlap a tile core.
///
/// The underlying slice is filtered on-the-fly to skip windows whose span does not intersect the
/// tile, so callers do not need to duplicate the overlap predicates.
pub struct TileWindowsIter<'a> {
    windows: &'a [(u64, u64, u64)],
    next_idx: usize,
    end_idx: usize,
    core_start: u64,
    core_end: u64,
}

impl<'a> Iterator for TileWindowsIter<'a> {
    type Item = &'a (u64, u64, u64);

    /// Produces the next window that intersects the stored tile core.
    ///
    /// Cached index bounds may include neighbouring candidates, so this method checks overlap on
    /// the fly and discards false positives before yielding.
    ///
    /// # Returns
    /// `Some(window)` when another overlapping window is available, otherwise `None` at the end of
    /// the span.
    fn next(&mut self) -> Option<Self::Item> {
        while self.next_idx < self.end_idx {
            let window = &self.windows[self.next_idx];
            self.next_idx += 1;
            // Only return windows that truly intersect the tile core, even if the cached span
            // contains neighbouring candidates with the same chromosome ordering
            if window.1 > self.core_start && window.0 < self.core_end {
                return Some(window);
            }
        }
        None
    }
}

/// Locates the window slice that could overlap the provided core without using cached spans.
///
/// The helper performs the same two-pointer scan as the cached precomputation but confined to the
/// given window list, making it useful when a tile span was not precomputed.
///
/// # Parameters
/// - `windows`: Start-sorted window triples `(start, end, idx)` for the chromosome.
/// - `core_start`: Inclusive core start position in absolute coordinates.
/// - `core_end`: Exclusive core end position in absolute coordinates.
///
/// # Returns
/// A pair `(left, right)` giving the half-open window index range whose members may overlap.
fn span_bounds_without_cache(
    windows: &[(u64, u64, u64)],
    core_start: u64,
    core_end: u64,
) -> (usize, usize) {
    let mut left = 0usize;
    while left < windows.len() && windows[left].1 <= core_start {
        left += 1;
    }

    let mut right = left;
    while right < windows.len() && windows[right].0 < core_end {
        right += 1;
    }

    (left, right)
}

/// Provides an iterator over the windows that overlap a given tile core.
///
/// It reuses a cached span when available, falling back to a local scan otherwise, and clamps the
/// resulting indices to the source slice before constructing the iterator.
///
/// # Parameters
/// - `windows`: Start-sorted window triples `(start, end, idx)` for the chromosome.
/// - `tile`: Tile whose core boundaries determine the overlap test.
/// - `span`: Optional cached span previously produced by `precompute_tile_window_spans`.
///
/// # Returns
/// A `TileWindowsIter` positioned to stream the overlapping windows.
pub fn overlapping_windows_for_tile<'a>(
    windows: &'a [(u64, u64, u64)],
    tile: &Tile,
    span: Option<&TileWindowSpan>,
) -> TileWindowsIter<'a> {
    let core_start = tile.core_start as u64;
    let core_end = tile.core_end as u64;

    let (start_idx, end_idx) = match span {
        Some(span) if !span.is_empty() => (span.first_idx, span.last_idx_exclusive),
        _ => span_bounds_without_cache(windows, core_start, core_end),
    };

    TileWindowsIter {
        windows,
        next_idx: start_idx.min(windows.len()),
        end_idx: end_idx.min(windows.len()),
        core_start,
        core_end,
    }
}

/// Finds the extreme start and end among windows that overlap the tile core.
///
/// The function iterates over the overlapping windows (using the cached span when provided) and
/// tracks the minimum start and maximum end, yielding `None` when no windows intersect the core.
///
/// # Parameters
/// - `windows`: Start-sorted window triples `(start, end, idx)` for the chromosome.
/// - `tile`: Tile whose core is used for overlap checks.
/// - `span`: Optional cached span for fast window access.
///
/// # Returns
/// `Some((min_start, max_end))` when at least one window overlaps; otherwise `None`.
pub fn tile_window_min_max(
    windows: &[(u64, u64, u64)],
    tile: &Tile,
    span: Option<&TileWindowSpan>,
) -> Option<(u64, u64)> {
    let mut iter = overlapping_windows_for_tile(windows, tile, span);
    let first = iter.next()?;
    let mut min_start = first.0;
    let mut max_end = first.1;

    for &(start, end, _) in iter {
        if start < min_start {
            min_start = start;
        }
        if end > max_end {
            max_end = end;
        }
    }

    Some((min_start, max_end))
}

/// Tightens a tile's fetch bounds to the observed window span while respecting halos.
///
/// The narrowed span subtracts the left/right halo from the minimum/maximum window edges and then
/// clamps the result back onto the tile fetch interval and chromosome length, discarding empty
/// ranges.
///
/// # Parameters
/// - `tile`: Tile providing the original fetch and core coordinates.
/// - `chrom_len`: Total length of the chromosome in bases.
/// - `min_ws`: Minimum window start observed among overlaps.
/// - `max_we`: Maximum window end observed among overlaps.
///
/// # Returns
/// `Some((start, end))` as absolute fetch limits when a non-empty span remains; otherwise `None`.
#[inline]
pub fn clamp_fetch_to_window_span(
    tile: &Tile,
    chrom_len: u64,
    min_ws: u64,
    max_we: u64,
) -> Option<(i64, i64)> {
    if min_ws >= max_we {
        return None;
    }

    let left_halo = (tile.core_start as u64).saturating_sub(tile.fetch_start as u64);
    let right_halo = (tile.fetch_end as u64).saturating_sub(tile.core_end as u64);

    let narrowed_start = min_ws.saturating_sub(left_halo);
    let narrowed_end = max_we.saturating_add(right_halo);

    let start_u64 = narrowed_start.max(tile.fetch_start as u64);
    let end_u64 = narrowed_end.min(tile.fetch_end as u64).min(chrom_len);

    (start_u64 < end_u64).then(|| (start_u64 as i64, end_u64 as i64))
}

/// Builds non-overlapping core tiles and their fetch halos for the requested contigs.
///
/// Tiles are generated sequentially across each chromosome. When an alignment size is supplied and
/// a large enough number of bins fits into the tile length, the core width is rounded down to the
/// nearest multiple to keep window edges aligned; otherwise the original tile length is used.
///
/// # Parameters
/// - `chromosomes`: Chromosome names in the order tiles should be produced.
/// - `contigs`: Mapping from chromosome name to contig metadata (tid and length).
/// - `tile_bp`: Desired tile core length in bases.
/// - `halo_bp`: Halo length applied to both sides of each tile.
/// - `align_bp`: Optional bin size that cores should align to when feasible.
///
/// # Returns
/// A tuple `(tiles, guaranteed_aligned)` where `tiles` contains every generated `Tile` and
/// `guaranteed_aligned` flags whether cores were aligned to `align_bp`.
pub fn build_tiles(
    chromosomes: &[String],
    contigs: &Contigs,
    tile_bp: u32,
    halo_bp: u32,
    align_bp: Option<u64>,
) -> anyhow::Result<(Vec<Tile>, bool)> {
    let mut tiles = Vec::new();

    // Decide the effective core size in bases (possibly aligned)
    let (effective_tile_bp, guaranteed_aligned) = match align_bp {
        Some(bin_size) if bin_size > 0 && bin_size <= tile_bp as u64 => {
            if (tile_bp as u64) % bin_size == 0 {
                (tile_bp, true)
            } else if (tile_bp as u64) / bin_size >= 10 {
                let k = (tile_bp as u64) / bin_size;
                ((k * bin_size) as u32, true) // Never drop below one bin
            } else {
                (tile_bp as u32, false)
            }
        }
        _ => (tile_bp as u32, false),
    };

    for chr in chromosomes {
        let &(tid, chrom_len_u32) = contigs
            .contigs
            .get(chr)
            .ok_or_else(|| anyhow::anyhow!("missing contig for '{}'", chr))?;
        let chrom_len = chrom_len_u32 as u32;

        let mut start = 0u32;
        let mut idx = 0u32;
        while start < chrom_len {
            let core_end = (start.saturating_add(effective_tile_bp)).min(chrom_len);

            // Halos do not need to be aligned; they are just fetch guards.
            let fetch_start = start.saturating_sub(halo_bp);
            let fetch_end = (core_end.saturating_add(halo_bp)).min(chrom_len);

            tiles.push(Tile {
                chr: chr.clone(),
                tid: tid as i32,
                index: idx,
                core_start: start,
                core_end,
                fetch_start,
                fetch_end,
            });

            idx += 1;
            start = core_end;
        }
    }

    // Just in case we decide to move start of cores in the future
    // We'll have an extensive debug test
    #[cfg(debug_assertions)]
    if let Some(bs) = align_bp {
        if guaranteed_aligned {
            for t in &tiles {
                // Starts/ends of *cores* (not final chromosome end) line up on the grid
                debug_assert_eq!((t.core_start as u64) % bs, 0);
            }
        }
    }

    Ok((tiles, guaranteed_aligned))
}

/// Adds a fragment's coverage contribution into the tile-local accumulator.
///
/// Segmented fragments are processed segment by segment, while simple fragments are clipped once;
/// in both cases the coordinates are translated into the tile's local frame before they are added.
/// The caller must provide a coverage array sized to the tile core.
///
/// # Parameters
/// - `cp`: Tile-local coverage structure to update.
/// - `fragment`: Fragment carrying absolute coordinates and optional segments.
/// - `weight`: Weight applied when inserting the fragment.
/// - `core_start`: Inclusive start of the tile core in absolute coordinates.
/// - `core_end`: Exclusive end of the tile core in absolute coordinates.
///
/// # Returns
/// `Ok(())` when the contribution is applied successfully, or an error bubbling from the coverage
/// accumulator.
#[inline]
pub fn add_fragment_clipped_to_core(
    cp: &mut Coverage,
    fragment: &FragmentWithSegments,
    weight: f32,
    core_start: u32,
    core_end: u32,
) -> Result<()> {
    // Use explicit segments if present
    if let Some(segments) = &fragment.segments {
        for &(seg_start_abs, seg_end_abs) in segments {
            let s = seg_start_abs.max(core_start);
            let e = seg_end_abs.min(core_end);
            if s < e {
                // Skips fragments completely outside tile
                // Shift to tile-local coordinates
                let local = Fragment {
                    tid: fragment.tid,
                    start: s - core_start,
                    end: e - core_start,
                };
                cp.add_fragment_weighted(local, weight)?;
            }
        }
    } else {
        // No explicit segments -> treat as one span (this already encodes your include_inter_mate_gap policy)
        let s = fragment.start.max(core_start);
        let e = fragment.end.min(core_end);
        if s < e {
            // Skips fragments completely outside tile
            // Shift to tile-local coordinates
            let local = Fragment {
                tid: fragment.tid,
                start: s - core_start,
                end: e - core_start,
            };
            cp.add_fragment_weighted(local, weight)?;
        }
    }
    Ok(())
}

/// What the tile should write
pub enum TileMode<'w> {
    /// Whole positional coverage for the core,
    /// or windowed positional coverage without index (unique positions)
    Positional {
        windows: Option<&'w [(u64, u64, u64)]>, // Per-chr windows if provided
        out_path: std::path::PathBuf,           // Per-tile file path
        indexed: bool,                          // Whether to save index
    },
    AggregatesByBed {
        windows: &'w [(u64, u64, u64)],    // Per-chr windows
        masked: bool,                      // Use masked counts/sums
        partials_out: std::path::PathBuf, // Cross-boundary windows (idx, sum, allowed, blacklisted)
        cross_idx_out: std::path::PathBuf, // Sidecar listing crossers
    },
    AggregatesBySize {
        window_bp: u64,                    // Fixed window size in bases
        masked: bool,                      // Use masked counts/sums
        finals_out: std::path::PathBuf,    // Final windows that do not need reducing
        partials_out: std::path::PathBuf, // Cross-boundary windows (idx, sum, allowed, blacklisted)
        cross_idx_out: std::path::PathBuf, // Sidecar listing crossers
        guaranteed_aligned: bool, // Tiles and window_bp align, so write per-tile FINALS (no reducer needed)
    },
}

/// Filters a chromosome's windows down to those that touch the tile core.
///
/// The iterator simply checks for half-open interval overlap between each window and the core
/// bounds expressed as absolute coordinates.
///
/// # Parameters
/// - `windows_chr`: Chromosome-specific windows `(start, end, idx)` in start order.
/// - `core_start`: Inclusive tile core start in absolute bases.
/// - `core_end`: Exclusive tile core end in absolute bases.
///
/// # Returns
/// An iterator yielding references to the overlapping windows.
#[inline]
pub fn windows_overlapping_core<'a>(
    windows_chr: &'a [(u64, u64, u64)],
    core_start: u32,
    core_end: u32,
) -> impl Iterator<Item = &'a (u64, u64, u64)> {
    let cs = core_start as u64;
    let ce = core_end as u64;
    windows_chr
        .iter()
        .filter(move |&&(ws, we, _idx)| we > cs && ws < ce)
}

/// Extracts the tile index suffix from a coverage file name.
///
/// The search proceeds right-to-left and returns the first segment that contains only ASCII digits,
/// making it tolerant to multi-part extensions such as `.tsv.zst`.
///
/// # Parameters
/// - `file_name`: File name to inspect.
///
/// # Returns
/// `Some(index)` when a numeric segment is found; otherwise `None`.
pub fn parse_tile_index(file_name: &str) -> Option<u32> {
    for seg in file_name.rsplit('.') {
        if !seg.is_empty() && seg.chars().all(|c| c.is_ascii_digit()) {
            return seg.parse().ok();
        }
    }
    None
}

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
            if fname.starts_with(per_tile_prefix) && fname.contains(&format!(".{chr}.")) {
                if let Some(idx) = parse_tile_index(fname) {
                    chr_files.push((idx, path));
                }
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
/// A compressed header frame is emitted before concatenating the tile frames in genomic order, so
/// the resulting file is still a valid zstd concatenation stream suitable for downstream tools.
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
            if fname.starts_with(per_tile_prefix) && fname.contains(&format!(".{chr}.")) {
                if let Some(idx) = parse_tile_index(fname) {
                    chr_files.push((idx, path));
                }
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
///
/// # Returns
/// `Some((start, end))` giving the adjusted fetch coordinates, or `None` when no fetch is needed.
pub fn adapt_fetch_to_extreme_windows(
    tile: &Tile,
    tile_span: Option<&TileWindowSpan>,
    mode: &TileMode<'_>,
    chrom_len: u32,
) -> Option<(i64, i64)> {
    let chrom_len_u64 = chrom_len as u64;

    // Decide the fetch interval based on mode/windows.
    // For whole-genome positional: use the full tile fetch band.
    // For windowed runs: restrict to [min_overlapping_window, max_overlapping_window] ± halo,
    // intersected with the tile’s existing fetch band.
    match mode {
        TileMode::Positional { windows: None, .. } => {
            Some((tile.fetch_start as i64, tile.fetch_end as i64))
        }
        TileMode::Positional {
            windows: Some(wchr),
            ..
        } => {
            let (min_ws, max_we) = tile_window_min_max(wchr, tile, tile_span)?;
            clamp_fetch_to_window_span(tile, chrom_len_u64, min_ws, max_we)
        }
        TileMode::AggregatesByBed { windows: wchr, .. } => {
            let (min_ws, max_we) = tile_window_min_max(wchr, tile, tile_span)?;
            clamp_fetch_to_window_span(tile, chrom_len_u64, min_ws, max_we)
        }
        TileMode::AggregatesBySize { window_bp, .. } => {
            let cs = tile.core_start as u64;
            let ce = tile.core_end as u64;
            if cs >= chrom_len_u64 {
                return None;
            }
            let owned_window_bp = *window_bp;
            let k_lo = cs / owned_window_bp;
            let k_hi = (ce.saturating_sub(1)) / owned_window_bp;
            let min_ws = k_lo * owned_window_bp;
            let max_we = ((k_hi + 1) * owned_window_bp).min(chrom_len_u64);
            clamp_fetch_to_window_span(tile, chrom_len_u64, min_ws, max_we)
        }
    }
}

/// Generates a random alphanumeric suffix of the requested length.
///
/// This helper is used for temporary directory creation where collisions must be unlikely.
///
/// # Parameters
/// - `n`: Number of characters to sample.
///
/// # Returns
/// A string consisting of `n` random ASCII letters or digits.
fn random_suffix(n: usize) -> String {
    rand::rng()
        .sample_iter(&Alphanumeric)
        .take(n)
        .map(char::from)
        .collect()
}

/// Creates a unique temporary directory within the output tree.
///
/// The function attempts a handful of random suffixes before falling back to a timestamp-based
/// name, ensuring directories can be created even under heavy parallelism.
///
/// # Parameters
/// - `base_out`: Root directory that should contain the temporary directory.
/// - `prefix`: Human-readable prefix used when building the directory name.
///
/// # Returns
/// Path to the created temporary directory.
pub fn make_temp_dir(
    base_out: &std::path::Path,
    prefix: &str,
) -> anyhow::Result<std::path::PathBuf> {
    // Try a few times just in case
    for _ in 0..8 {
        let suffix = random_suffix(10);
        let p = base_out.join(format!("tmp.{prefix}.{suffix}"));
        if !p.exists() {
            std::fs::create_dir_all(&p)?;
            return Ok(p);
        }
    }
    // Fallback: timestamped
    let ts = chrono::Utc::now().timestamp_millis();
    let p = base_out.join(format!("tmp.{prefix}.{ts}"));
    std::fs::create_dir_all(&p)?;
    Ok(p)
}

/// Rounds a floating-point value to a fixed number of decimal places.
///
/// When the requested precision is non-positive the value is rounded to the nearest integer.
///
/// # Parameters
/// - `x`: Value to round.
/// - `decimals`: Number of decimal places to preserve.
///
/// # Returns
/// The rounded value.
pub fn round_to(x: f64, decimals: i32) -> f64 {
    if decimals <= 0 {
        return x.round();
    }
    let f = 10f64.powi(decimals);
    (x * f).round() / f
}

/// Rounds using a precomputed scaling factor for repeated operations.
///
/// This variant avoids recomputing the power-of-ten factor when many values share the same
/// precision requirement.
///
/// # Parameters
/// - `x`: Value to round.
/// - `factor`: Precomputed `10f64.powi(decimals)` for the desired precision.
///
/// # Returns
/// The rounded value.
pub fn round_to_with_precomputed_factor(x: f64, factor: f64) -> f64 {
    if factor == 1.0 {
        return x.round();
    }
    (x * factor).round() / factor
}

/// Lightweight adapter for printing rounded numbers without heap allocations.
///
/// The stored value and decimal precision are used to format coverage values compactly in hot
/// loops.
pub struct CompactNumber {
    pub v: f64,
    pub decimals: i32,
}

impl fmt::Display for CompactNumber {
    /// Formats the number using the stored decimal precision without allocating.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Normalize negative zero (can apparently happen after rounding tiny negatives)
        // After this, values will be non-zero
        if self.v == 0.0 {
            return f.write_str("0");
        }
        if self.decimals <= 0 {
            // Integer path: Avoid float fmt; Round then print int
            let r = self.v.round();
            // After integer rounding, also normalize -0 just in case
            if r == 0.0 {
                return f.write_str("0");
            }
            // Write directly; no heap
            return write!(f, "{:.0}", r);
        }
        // Stack buffer; Write fixed decimals
        let mut buf = arrayvec::ArrayString::<64>::new();
        // Write with N decimals
        let _ = fmt::write(
            &mut buf,
            format_args!("{:.*}", self.decimals as usize, self.v),
        );
        // Trim zeros and trailing dot
        while buf.as_bytes().last() == Some(&b'0') {
            buf.pop();
        }
        if buf.as_bytes().last() == Some(&b'.') {
            buf.pop();
        }
        f.write_str(buf.as_str())
    }
}

/// Writes BedGraph segments for a window of coverage values.
///
/// Consecutive bases with the same rounded value are merged into runs, any masked positions are
/// omitted entirely, and absolute coordinates are reconstructed from the tile origin.
///
/// # Parameters
/// - `chr`: Chromosome name to print.
/// - `cov`: Tile-local coverage values.
/// - `mask`: Optional mask slice where `1` marks blacklisted bases.
/// - `local_start_idx`: Inclusive start index inside `cov`.
/// - `local_end_idx`: Exclusive end index inside `cov`.
/// - `tile_abs_start`: Absolute coordinate of `cov[0]`.
/// - `decimals`: Number of decimals to keep when grouping runs.
/// - `keep_zero_runs`: Whether runs with value zero should be emitted.
/// - `out`: Writer receiving the BedGraph lines.
///
/// # Returns
/// `Ok(())` on success, or the underlying I/O error.
pub fn emit_bedgraph_runs<W: Write>(
    chr: &str,
    cov: &[f32],
    mask: Option<&[u8]>,    // 1 = blacklisted(masked), 0 = allowed
    local_start_idx: usize, // Local start (inclusive)
    local_end_idx: usize,   // Local end (exclusive)
    tile_abs_start: u64,    // Absolute position of index 0 in `cov` (tile.core_start)
    decimals: i32,          // Decimals to round coverage
    keep_zero_runs: bool,   // Whether to write zero-runs
    out: &mut W,
) -> Result<()> {
    if local_start_idx >= local_end_idx {
        return Ok(());
    }

    visit_runs_in_window(
        cov,
        mask,
        local_start_idx,
        local_end_idx,
        decimals,
        keep_zero_runs,
        |run_lo, run_hi, value| {
            let s_abs = tile_abs_start + run_lo as u64;
            let e_abs = tile_abs_start + run_hi as u64;
            // Ignore write errors here; bubbled up by caller on flush
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}",
                chr,
                s_abs,
                e_abs,
                CompactNumber { v: value, decimals },
            );
        },
    );

    Ok(())
}

/// Writes run-length encoded coverage for a single window in TSV form.
///
/// The helper mirrors `emit_bedgraph_runs` but optionally appends the window's original index to
/// each line when provided, which is needed for downstream grouping workflows.
///
/// # Parameters
/// - `chr`: Chromosome name to print.
/// - `cov`: Tile-local coverage values.
/// - `mask`: Optional mask slice where `1` marks blacklisted bases.
/// - `local_start_idx`: Inclusive start index inside `cov`.
/// - `local_end_idx`: Exclusive end index inside `cov`.
/// - `tile_abs_start`: Absolute coordinate of `cov[0]`.
/// - `orig_idx`: Optional original window index to append.
/// - `decimals`: Number of decimals to keep when grouping runs.
/// - `keep_zero_runs`: Whether runs with value zero should be emitted.
/// - `out`: Writer receiving the TSV lines.
///
/// # Returns
/// `Ok(())` on success, or the underlying I/O error.
pub fn emit_windowed_runs<W: Write>(
    chr: &str,
    cov: &[f32],
    mask: Option<&[u8]>,    // 1 = blacklisted(masked), 0 = allowed
    local_start_idx: usize, // Local start (inclusive)
    local_end_idx: usize,   // Local end (exclusive)
    tile_abs_start: u64,    // Absolute position of index 0 in `cov` (tile.core_start)
    orig_idx: Option<u64>,  // Window’s original index
    decimals: i32,          // Decimals to round coverage
    keep_zero_runs: bool,   // Whether to write zero-runs
    out: &mut W,
) -> Result<()> {
    if local_start_idx >= local_end_idx {
        return Ok(());
    }
    visit_runs_in_window(
        cov,
        mask,
        local_start_idx,
        local_end_idx,
        decimals,
        keep_zero_runs,
        |run_lo, run_hi, value| {
            let s_abs = tile_abs_start + run_lo as u64;
            let e_abs = tile_abs_start + run_hi as u64;
            let _ = if let Some(idx) = orig_idx {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    chr,
                    s_abs,
                    e_abs,
                    CompactNumber { v: value, decimals },
                    idx
                )
            } else {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}",
                    chr,
                    s_abs,
                    e_abs,
                    CompactNumber { v: value, decimals },
                )
            };
        },
    );

    Ok(())
}

/// Iterates over contiguous runs of equal rounded coverage within a slice.
///
/// Masked indices are skipped so that the visitor sees only unmasked stretches. Rounding is
/// applied before comparing values, ensuring that small floating-point perturbations do not split
/// runs unnecessarily.
///
/// # Parameters
/// - `cov`: Tile-local coverage values.
/// - `mask`: Optional mask where `1` marks blacklisted bases.
/// - `local_start_idx`: Inclusive start index inside `cov`.
/// - `local_end_idx`: Exclusive end index inside `cov`.
/// - `decimals`: Number of decimals to keep when grouping runs.
/// - `keep_zero_runs`: Whether runs with value zero should be visited.
/// - `on_run`: Visitor called with `(run_start, run_end, rounded_value)`.
#[inline]
fn visit_runs_in_window(
    cov: &[f32],
    mask: Option<&[u8]>,
    local_start_idx: usize,
    local_end_idx: usize,
    decimals: i32,
    keep_zero_runs: bool,
    on_run: impl FnMut(usize, usize, f64),
) {
    let m = mask.unwrap_or(&[]);
    let m_has_elements = !m.is_empty();
    if m_has_elements {
        visit_runs_masked(
            cov,
            m,
            local_start_idx,
            local_end_idx,
            decimals,
            keep_zero_runs,
            on_run,
        )
    } else {
        visit_runs_unmasked(
            cov,
            local_start_idx,
            local_end_idx,
            decimals,
            keep_zero_runs,
            on_run,
        )
    }
}

/// Visits runs when no masking is applied.
///
/// Values are rounded using the provided precision and adjacent equal values are merged before the
/// visitor callback is invoked.
///
/// # Parameters
/// - `cov`: Tile-local coverage values.
/// - `local_start_idx`: Inclusive start index inside `cov`.
/// - `local_end_idx`: Exclusive end index inside `cov`.
/// - `decimals`: Number of decimals to keep when grouping runs.
/// - `keep_zero_runs`: Whether runs with value zero should be visited.
/// - `on_run`: Visitor called with `(run_start, run_end, rounded_value)`.
#[inline]
fn visit_runs_unmasked(
    cov: &[f32],
    local_start_idx: usize,
    local_end_idx: usize,
    decimals: i32,
    keep_zero_runs: bool,
    mut on_run: impl FnMut(usize, usize, f64),
) {
    let mut i = local_start_idx;
    let rounding_factor = if decimals <= 0 {
        1.0
    } else {
        10f64.powi(decimals)
    };

    while i < local_end_idx {
        // Start run
        let run_start_idx = i;

        let value0 = round_to_with_precomputed_factor(cov[i] as f64, rounding_factor);

        // Extend run
        let mut j = i + 1;
        while j < local_end_idx {
            let vj = round_to_with_precomputed_factor(cov[j] as f64, rounding_factor);
            if vj != value0 {
                break;
            }
            j += 1;
        }

        // Optionally drop zero runs
        if value0 == 0.0 && !keep_zero_runs {
            i = j;
            continue;
        }

        on_run(run_start_idx, j, value0);
        i = j;
    }
}

/// Visits runs while respecting a binary mask that excludes certain bases.
///
/// The visitor is skipped whenever the mask marks the base as blacklisted, effectively splitting
/// runs around masked positions.
///
/// # Parameters
/// - `cov`: Tile-local coverage values.
/// - `m`: Mask where `1` marks blacklisted bases.
/// - `local_start_idx`: Inclusive start index inside `cov`.
/// - `local_end_idx`: Exclusive end index inside `cov`.
/// - `decimals`: Number of decimals to keep when grouping runs.
/// - `keep_zero_runs`: Whether runs with value zero should be visited.
/// - `on_run`: Visitor called with `(run_start, run_end, rounded_value)`.
#[inline]
fn visit_runs_masked(
    cov: &[f32],
    m: &[u8],
    local_start_idx: usize,
    local_end_idx: usize,
    decimals: i32,
    keep_zero_runs: bool,
    mut on_run: impl FnMut(usize, usize, f64),
) {
    let mut i = local_start_idx;
    let rounding_factor = if decimals <= 0 {
        1.0
    } else {
        10f64.powi(decimals)
    };

    while i < local_end_idx {
        // Skip masked base
        if m[i] == 1 {
            i += 1;
            continue;
        }

        // Start run
        let run_start_idx = i;

        let value0 = round_to_with_precomputed_factor(cov[i] as f64, rounding_factor);

        // Extend run
        let mut j = i + 1;
        while j < local_end_idx {
            if m[j] == 1 {
                break;
            }
            let vj = round_to_with_precomputed_factor(cov[j] as f64, rounding_factor);
            if vj != value0 {
                break;
            }
            j += 1;
        }

        // Optionally drop zero runs
        if value0 == 0.0 && !keep_zero_runs {
            i = j;
            continue;
        }

        on_run(run_start_idx, j, value0);
        i = j;
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

/// Clips an absolute interval to the tile core and converts it to local indices.
///
/// The helper returns both the local indices and the clipped absolute bounds so callers can reuse
/// whichever representation is needed.
///
/// # Parameters
/// - `abs_start`: Inclusive absolute start of the interval.
/// - `abs_end`: Exclusive absolute end of the interval.
/// - `core_start`: Inclusive absolute start of the tile core.
/// - `core_end`: Exclusive absolute end of the tile core.
///
/// # Returns
/// `Some((local_start_idx, local_end_idx, clipped_start, clipped_end))` when the intervals
/// intersect; otherwise `None`.
#[inline]
pub fn intersect_abs_with_core_to_local(
    abs_start: u64,
    abs_end: u64,
    core_start: u32,
    core_end: u32,
) -> Option<(usize, usize, u32, u32)> {
    let s_abs = (abs_start as u32).max(core_start);
    let e_abs = (abs_end as u32).min(core_end);
    if e_abs <= s_abs {
        return None;
    }
    let local_start_idx = (s_abs - core_start) as usize;
    let local_end_idx = (e_abs - core_start) as usize;
    Some((local_start_idx, local_end_idx, s_abs, e_abs))
}
