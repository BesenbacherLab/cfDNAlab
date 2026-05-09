use crate::shared::bam::Contigs;
use crate::shared::interval::{IndexedInterval, Interval};
use crate::shared::io::dot_join;
use anyhow::{Context, ensure};
use rand::{Rng, distr::Alphanumeric};
use std::io::ErrorKind;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};
use tracing::warn;

const TEMP_DIR_CLEANUP_TARGET: &str = "temp-dir-cleanup";
static TEMP_DIR_CTRL_C_CLEANUP_STARTED: AtomicBool = AtomicBool::new(false);
static TEMP_DIR_CTRL_C_HANDLER_INSTALLED: OnceLock<()> = OnceLock::new();
static TEMP_DIR_CTRL_C_REGISTRY: OnceLock<Mutex<Vec<PathBuf>>> = OnceLock::new();

/// A processing tile for one chromosome
#[derive(Debug, Clone)]
pub struct Tile {
    pub chr: String,
    pub tid: i32,
    pub index: u32, // 0-based index within chromosome
    pub core: Interval<u32>,
    pub fetch: Interval<u32>,
}

impl Tile {
    /// Create a tile from checked half-open core and fetch intervals.
    ///
    /// The fetch interval must fully cover the tile core.
    pub fn new(
        chr: String,
        tid: i32,
        index: u32,
        core: Interval<u32>,
        fetch: Interval<u32>,
    ) -> crate::Result<Self> {
        if !fetch.contains_interval(core) {
            return Err(crate::Error::TileFetchDoesNotCoverCore);
        }
        Ok(Self {
            chr,
            tid,
            index,
            core,
            fetch,
        })
    }

    /// Create a tile from raw half-open core and fetch bounds.
    ///
    /// This is a convenience constructor for call sites that still work with
    /// coordinates. It validates the bounds as intervals and then delegates to
    /// the typed constructor.
    pub fn from_coords(
        chr: String,
        tid: i32,
        index: u32,
        core_start: u32,
        core_end: u32,
        fetch_start: u32,
        fetch_end: u32,
    ) -> crate::Result<Self> {
        let core = Interval::new(core_start, core_end)?;
        let fetch = Interval::new(fetch_start, fetch_end)?;
        Self::new(chr, tid, index, core, fetch)
    }

    #[inline]
    pub fn core_start(&self) -> u32 {
        self.core.start()
    }

    #[inline]
    pub fn core_end(&self) -> u32 {
        self.core.end()
    }

    #[inline]
    pub fn fetch_start(&self) -> u32 {
        self.fetch.start()
    }

    #[inline]
    pub fn fetch_end(&self) -> u32 {
        self.fetch.end()
    }

    /// Ensure a BAM reader opened for this tile resolved the same chromosome tid.
    ///
    /// `Tile::tid` is stored as `i32` because it follows the tiling inputs, while
    /// rust-htslib returns BAM tids as `u32`. Keep the conversion and comparison
    /// together so callers cannot accidentally wrap negative tids with `as u32`.
    pub fn ensure_matches_bam_tid(&self, bam_tid: u32) -> anyhow::Result<()> {
        let tile_tid =
            u32::try_from(self.tid).context("tile tid is negative for BAM tid comparison")?;
        ensure!(
            bam_tid == tile_tid,
            "BAM tid mismatch for chromosome {}: tile tid is {}, BAM tid is {}",
            self.chr,
            tile_tid,
            bam_tid
        );
        Ok(())
    }
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
    F: FnMut(&str) -> &'a [IndexedInterval<u64>],
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
            let core_start = tile.core_start() as u64;
            let core_end = tile.core_end() as u64;
            let left_bound = core_start.saturating_sub(left_halo);
            let right_bound = core_end.saturating_add(right_halo);

            (w_left, w_right) =
                advance_window_span_bounds(windows, w_left, w_right, left_bound, right_bound);

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
            if first_candidate.start() >= right_bound {
                spans[idx] = None;
                continue;
            }

            // Some commands consume this cached span directly as their production candidate set
            // (`lengths`, `ends`, `gc_bias`), while core-overlap commands tighten it later with
            // iterator-level overlap checks. So the half-open boundary logic here must itself stay
            // correct; later filtering is model-dependent and not universal.
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
    windows: &'a [IndexedInterval<u64>],
    next_idx: usize,
    end_idx: usize,
    core_start: u64,
    core_end: u64,
}

impl<'a> Iterator for TileWindowsIter<'a> {
    type Item = &'a IndexedInterval<u64>;

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
            if window.end() > self.core_start && window.start() < self.core_end {
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
    windows: &[IndexedInterval<u64>],
    core_start: u64,
    core_end: u64,
) -> (usize, usize) {
    advance_window_span_bounds(windows, 0, 0, core_start, core_end)
}

/// Advance cached half-open window-span bounds through a start-sorted window list.
///
/// The caller supplies previously valid `left` and `right` scan positions from an earlier tile on
/// the same chromosome. Because tiles are processed left-to-right and the requested bounds move
/// monotonically forward, those positions are safe lower bounds for the next scan. The function
/// advances them just far enough to satisfy the new half-open bounds:
///
/// - advance `left` while windows end at or before `left_bound`
/// - advance `right` while windows start before `right_bound`
///
/// Passing `left = 0` and `right = 0` reproduces the uncached single-tile scan used by
/// `span_bounds_without_cache(...)`.
///
/// Parameters
/// ----------
/// - `windows`:
///   Start-sorted window triples `(start, end, idx)` for one chromosome
/// - `left`:
///   Previous lower-bound scan position for the first surviving candidate window
/// - `right`:
///   Previous lower-bound scan position for the exclusive right edge of the candidate span
/// - `left_bound`:
///   New inclusive left pruning bound. Windows with `end <= left_bound` are discarded
/// - `right_bound`:
///   New exclusive right inclusion bound. Windows with `start < right_bound` stay in the span
///
/// Returns
/// -------
/// - `(usize, usize)`:
///   Updated half-open candidate-window index span `(left, right)`
fn advance_window_span_bounds(
    windows: &[IndexedInterval<u64>],
    mut left: usize,
    mut right: usize,
    left_bound: u64,
    right_bound: u64,
) -> (usize, usize) {
    while left < windows.len() && windows[left].end() <= left_bound {
        left += 1;
    }

    if right < left {
        right = left;
    }
    while right < windows.len() && windows[right].start() < right_bound {
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
    windows: &'a [IndexedInterval<u64>],
    tile: &Tile,
    span: Option<&TileWindowSpan>,
) -> TileWindowsIter<'a> {
    let core_start = tile.core_start() as u64;
    let core_end = tile.core_end() as u64;

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

/// Return the cached candidate-window span for a core-overlap tile/window model.
///
/// Coordinate space:
/// - consumes BED window coordinates
/// - returns a BED-window index span
///
/// Fragment ownership rule:
/// - none; this helper does not reason about fragment ownership
///
/// Counting or assignment interval assumption:
/// - none; relevance is defined only by BED/core overlap
///
/// Aligned fetch narrowing:
/// - not performed here
///
/// This helper answers only:
/// - which BED windows overlap the tile core?
///
/// It is not valid for fragment-reach commands such as `lengths`, `ends`, or `gc_bias`.
pub fn candidate_window_span_for_tile_core_overlap(
    windows: &[IndexedInterval<u64>],
    tile: &Tile,
) -> Option<TileWindowSpan> {
    let core_start = tile.core_start() as u64;
    let core_end = tile.core_end() as u64;
    let (first_idx, last_idx_exclusive) = span_bounds_without_cache(windows, core_start, core_end);
    (first_idx < last_idx_exclusive).then_some(TileWindowSpan {
        first_idx,
        last_idx_exclusive,
    })
}

/// Return the cached candidate-window span for a fragment-reach tile/window model.
///
/// Coordinate space:
/// - consumes BED window coordinates
/// - returns a BED-window index span
///
/// Fragment ownership rule:
/// - fragment is owned iff its aligned start lies in the tile core
///
/// Counting or assignment interval assumption:
/// - the caller must choose the left and right reach values that correspond to the command's
///   actual counting or assignment interval
///
/// Aligned fetch narrowing:
/// - not performed here
///
/// This helper answers only:
/// - which BED windows could receive counts from tile-owned fragments under the supplied reach?
///
/// It is not valid for future commands that use a different ownership rule, such as "fragment end
/// lies in the tile core", unless they define a separate helper or prove the same reach model.
pub fn candidate_window_span_for_tile_fragment_reach(
    windows: &[IndexedInterval<u64>],
    tile: &Tile,
    left_reach_bp: u64,
    right_reach_bp: u64,
) -> Option<TileWindowSpan> {
    let left_bound = (tile.core_start() as u64).saturating_sub(left_reach_bp);
    let right_bound = (tile.core_end() as u64).saturating_add(right_reach_bp);
    let (first_idx, last_idx_exclusive) =
        span_bounds_without_cache(windows, left_bound, right_bound);
    (first_idx < last_idx_exclusive).then_some(TileWindowSpan {
        first_idx,
        last_idx_exclusive,
    })
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
/// - `window_span`: Minimum-to-maximum interval covering the overlapping windows.
/// - `halo_bp`: Extra bases to keep on both sides of the observed window span before clamping
///   back onto the tile fetch interval.
///
/// # Returns
/// `Some(interval)` as absolute fetch limits when a non-empty span remains, otherwise `None`.
#[inline]
pub fn clamp_fetch_to_window_span(
    tile: &Tile,
    chrom_len: u64,
    window_span: Interval<u64>,
    halo_bp: u64,
) -> crate::Result<Option<Interval<u64>>> {
    let min_ws = window_span.start();
    let max_we = window_span.end();

    let tile_left_halo = (tile.core_start() as u64).saturating_sub(tile.fetch_start() as u64);
    let tile_right_halo = (tile.fetch_end() as u64).saturating_sub(tile.core_end() as u64);
    let effective_left_halo = tile_left_halo.max(halo_bp);
    let effective_right_halo = tile_right_halo.max(halo_bp);

    let narrowed_start = min_ws.saturating_sub(effective_left_halo);
    let narrowed_end = max_we.saturating_add(effective_right_halo);

    let start_u64 = narrowed_start.max(tile.fetch_start() as u64);
    let end_u64 = narrowed_end.min(tile.fetch_end() as u64).min(chrom_len);

    if start_u64 >= end_u64 {
        return Ok(None);
    }

    Ok(Some(Interval::new(start_u64, end_u64)?))
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
            if (tile_bp as u64).is_multiple_of(bin_size) {
                (tile_bp, true)
            } else if (tile_bp as u64) / bin_size >= 10 {
                let k = (tile_bp as u64) / bin_size;
                ((k * bin_size) as u32, true) // Never drop below one bin
            } else {
                (tile_bp, false)
            }
        }
        _ => (tile_bp, false),
    };

    for chr in chromosomes {
        let &(tid, chrom_len_u32) = contigs
            .contigs
            .get(chr)
            .ok_or_else(|| anyhow::anyhow!("missing contig for '{}'", chr))?;
        let chrom_len = chrom_len_u32;

        let mut start = 0u32;
        let mut idx = 0u32;
        while start < chrom_len {
            let core_end = (start.saturating_add(effective_tile_bp)).min(chrom_len);

            // Halos do not need to be aligned, they are just fetch guards
            let fetch_start = start.saturating_sub(halo_bp);
            let fetch_end = (core_end.saturating_add(halo_bp)).min(chrom_len);

            tiles.push(Tile::from_coords(
                chr.clone(),
                tid,
                idx,
                start,
                core_end,
                fetch_start,
                fetch_end,
            )?);

            idx += 1;
            start = core_end;
        }
    }

    // Just in case we decide to move start of cores in the future
    // We'll have an extensive debug test
    #[cfg(debug_assertions)]
    if let Some(alignment_bp) = align_bp
        && guaranteed_aligned
    {
        for t in &tiles {
            // Starts/ends of *cores* (not final chromosome end) line up on the grid
            debug_assert_eq!((t.core_start() as u64) % alignment_bp, 0);
        }
    }

    Ok((tiles, guaranteed_aligned))
}

/// What the tile should write
pub enum TileMode<'w> {
    /// Whole positional coverage for the core,
    /// or windowed positional coverage without index (unique positions)
    Positional {
        windows: Option<&'w [IndexedInterval<u64>]>, // Per-chr windows if provided
        out_path: std::path::PathBuf,                // Per-tile file path
        indexed: bool,                               // Whether to save index
    },
    AggregatesByBed {
        windows: &'w [IndexedInterval<u64>], // Per-chr windows
        masked: bool,                        // Use masked counts/sums
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
pub fn windows_overlapping_core(
    windows_chr: &[IndexedInterval<u64>],
    core_start: u32,
    core_end: u32,
) -> impl Iterator<Item = &IndexedInterval<u64>> {
    let core_start_abs = core_start as u64;
    let core_end_abs = core_end as u64;
    windows_chr
        .iter()
        .filter(move |window| window.end() > core_start_abs && window.start() < core_end_abs)
}

/// Extracts the tile index suffix from a coverage filename.
///
/// The search proceeds right-to-left and returns the first segment that contains only ASCII digits,
/// making it tolerant to multi-part extensions such as `.tsv.zst`.
///
/// # Parameters
/// - `file_name`: Filename to inspect.
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
/// - `prefix`: User-readable prefix used when building the directory name.
///
/// # Returns
/// Path to the created temporary directory.
pub fn make_temp_dir(base_out: &Path, prefix: &str) -> anyhow::Result<PathBuf> {
    // Try a few times just in case
    for _ in 0..8 {
        let suffix = random_suffix(10);
        let p = base_out.join(dot_join(&["tmp", prefix, &suffix]));
        if !p.exists() {
            std::fs::create_dir_all(&p)?;
            return Ok(p);
        }
    }
    // Fallback: timestamped
    let ts = chrono::Utc::now().timestamp_millis();
    let p = base_out.join(dot_join(&["tmp", prefix, &ts.to_string()]));
    std::fs::create_dir_all(&p)?;
    Ok(p)
}

/// Guard for per-run temporary directories.
///
/// Commands keep this value in scope for as long as tile files may be needed. The directory is
/// removed when the guard is dropped, so early returns clean up the same way as successful runs.
/// Call `remove()` when cleanup failure should be reported on the normal success path.
pub struct TempDirGuard {
    path: PathBuf,
    removed: bool,
}

impl TempDirGuard {
    /// Creates and guards a unique temporary directory inside `base_out`.
    pub fn new(base_out: &Path, prefix: &str) -> anyhow::Result<Self> {
        install_temp_dir_ctrl_c_cleanup_handler();
        let path = make_temp_dir(base_out, prefix)?;
        Ok(Self::from_existing_path(path))
    }

    /// Guards an existing temporary directory path.
    pub fn from_existing_path(path: PathBuf) -> Self {
        install_temp_dir_ctrl_c_cleanup_handler();
        register_temp_dir_for_ctrl_c_cleanup(&path);
        Self {
            path,
            removed: false,
        }
    }

    /// Returns the guarded temporary directory path.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Removes the guarded directory and disables drop-time cleanup after success.
    pub fn remove(&mut self) -> anyhow::Result<()> {
        if self.removed {
            return Ok(());
        }

        match std::fs::remove_dir_all(&self.path) {
            Ok(()) => {
                self.removed = true;
                unregister_temp_dir_for_ctrl_c_cleanup(&self.path);
                Ok(())
            }
            Err(err) if err.kind() == ErrorKind::NotFound => {
                self.removed = true;
                unregister_temp_dir_for_ctrl_c_cleanup(&self.path);
                Ok(())
            }
            Err(err) => Err(err.into()),
        }
    }
}

impl Drop for TempDirGuard {
    fn drop(&mut self) {
        if self.removed {
            return;
        }

        match std::fs::remove_dir_all(&self.path) {
            Ok(()) => {
                self.removed = true;
                unregister_temp_dir_for_ctrl_c_cleanup(&self.path);
            }
            Err(err) if err.kind() == ErrorKind::NotFound => {
                self.removed = true;
                unregister_temp_dir_for_ctrl_c_cleanup(&self.path);
            }
            Err(err) => {
                warn!(
                    target: TEMP_DIR_CLEANUP_TARGET,
                    "warning: failed to remove temp dir {}: {}",
                    self.path.display(),
                    err
                );
            }
        }
    }
}

fn temp_dir_ctrl_c_registry() -> &'static Mutex<Vec<PathBuf>> {
    TEMP_DIR_CTRL_C_REGISTRY.get_or_init(|| Mutex::new(Vec::new()))
}

fn install_temp_dir_ctrl_c_cleanup_handler() {
    TEMP_DIR_CTRL_C_HANDLER_INSTALLED.get_or_init(|| {
        if let Err(err) = ctrlc::set_handler(|| {
            if TEMP_DIR_CTRL_C_CLEANUP_STARTED.swap(true, Ordering::SeqCst) {
                return;
            }

            let paths = match temp_dir_ctrl_c_registry().lock() {
                Ok(paths) => paths.clone(),
                Err(poisoned) => poisoned.into_inner().clone(),
            };
            for path in paths {
                match std::fs::remove_dir_all(&path) {
                    Ok(()) => {}
                    Err(err) if err.kind() == ErrorKind::NotFound => {}
                    Err(err) => {
                        // Write directly to stderr because the handler exits immediately after cleanup
                        eprintln!(
                            "Warning: failed to remove temporary directory {}: {}",
                            path.display(),
                            err
                        );
                    }
                }
            }

            process::exit(130);
        }) {
            warn!(
                target: TEMP_DIR_CLEANUP_TARGET,
                "warning: failed to install Ctrl+C temp-dir cleanup handler: {}",
                err
            );
        }
    });
}

fn register_temp_dir_for_ctrl_c_cleanup(path: &Path) {
    let mut paths = match temp_dir_ctrl_c_registry().lock() {
        Ok(paths) => paths,
        Err(poisoned) => poisoned.into_inner(),
    };
    paths.push(path.to_path_buf());
}

fn unregister_temp_dir_for_ctrl_c_cleanup(path: &Path) {
    let mut paths = match temp_dir_ctrl_c_registry().lock() {
        Ok(paths) => paths,
        Err(poisoned) => poisoned.into_inner(),
    };
    paths.retain(|registered_path| registered_path != path);
}

#[cfg(test)]
mod tests {
    include!("tiled_run_tests.rs");
}
