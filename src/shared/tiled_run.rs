use crate::shared::bam::Contigs;
use rand::{Rng, distr::Alphanumeric};

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
/// escape the tile core expanded by optional halos. This preserves streaming behaviour and avoids
/// re-scanning the same windows for neighbouring tiles.
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
    left_halo: u64,
    right_halo: u64,
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
            let left_bound = core_start.saturating_sub(left_halo);
            let right_bound = core_end.saturating_add(right_halo);

            // Discard windows that end before the left bound (core minus halo)
            while w_left < windows_len && windows[w_left].1 <= left_bound {
                w_left += 1;
            }

            if w_right < w_left {
                w_right = w_left;
            }

            // Extend to cover every window whose start lies inside the right bound (core plus halo)
            while w_right < windows_len && windows[w_right].0 < right_bound {
                w_right += 1;
            }

            if w_left == windows_len {
                // No windows remain for this chromosome, so the rest of the tiles cannot overlap
                // TODO: Already initialized with None?
                spans[idx] = None;
                for span in spans[idx + 1..chr_tile_end].iter_mut() {
                    *span = None;
                }
                break;
            }

            // `w_left` now points at the first surviving window candidate for this tile
            let first_candidate = &windows[w_left];

            // Check if the earliest remaining window begins at/after the right bound, so later ones do too
            if first_candidate.0 >= right_bound {
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

            // Halos do not need to be aligned, they are just fetch guards
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
