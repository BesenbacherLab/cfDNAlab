use crate::commands::cli_common::WindowSpec;
use crate::commands::gc_bias::counting::{GCCounts, GCPrefixes};
use crate::shared::{
    bam::Contigs,
    bed::Windows,
    tiled_run::{Tile, TileWindowSpan, overlapping_windows_for_tile},
};
use anyhow::{Result, anyhow, bail, ensure};
use fxhash::FxHashMap;
use std::path::PathBuf;

/// Summary of windowing configuration.
pub struct WindowStats {
    pub avg_span: f64,
    pub total_windows: u64,
}

/// Compute window statistics without extra passes over the window set.
///
/// Returns the mean span per window and the total number of windows for the
/// configured mode (BED, fixed-size, or global).
pub fn compute_window_stats(
    window_opt: &WindowSpec,
    windows_map: Option<&FxHashMap<String, Windows>>,
    contigs: &Contigs,
    chromosomes: &[String],
) -> Result<WindowStats> {
    match window_opt {
        WindowSpec::Bed(_) => {
            let map = windows_map
                .ok_or_else(|| anyhow!("BED window spec requires loaded windows"))?;
            let mut total_len: u64 = 0;
            let mut count: u64 = 0;
            for chr in chromosomes {
                if let Some(ws) = map.get(chr) {
                    for (s, e, _) in ws.as_slice() {
                        total_len = total_len.saturating_add(e.saturating_sub(*s));
                        count += 1;
                    }
                }
            }
            ensure!(count > 0, "No windows to compute average span from");
            Ok(WindowStats {
                avg_span: total_len as f64 / count as f64,
                total_windows: count,
            })
        }
        WindowSpec::Global => {
            let mut total_len: u64 = 0;
            for chr in chromosomes {
                if let Some((_, len)) = contigs.contigs.get(chr) {
                    total_len = total_len.saturating_add(*len as u64);
                }
            }
            ensure!(
                total_len > 0,
                "Chromosome lengths unavailable for global window span"
            );
            Ok(WindowStats {
                avg_span: total_len as f64,
                total_windows: 1,
            })
        }
        WindowSpec::Size(win_size) => {
            ensure!(*win_size > 0, "Window size must be positive");
            let mut total_windows: u64 = 0;
            let mut total_span: u64 = 0;
            for chr in chromosomes {
                let Some((_, len_u32)) = contigs.contigs.get(chr) else {
                    bail!("Missing contig length for {}", chr);
                };
                let len = *len_u32 as u64;
                if len == 0 {
                    continue;
                }
                let full = len / *win_size;
                let rem = len % *win_size;
                let windows_for_chr = full + u64::from(rem > 0);
                total_windows = total_windows.saturating_add(windows_for_chr);
                total_span = total_span.saturating_add(len);
            }
            ensure!(total_windows > 0, "No windows computed for --by-size");
            Ok(WindowStats {
                avg_span: total_span as f64 / total_windows as f64,
                total_windows,
            })
        }
    }
}

#[derive(Clone, Debug)]
pub struct WindowState {
    pub idx: u64,
    pub start: u64,
    pub end: u64,
    pub contained: bool,
    pub counts: GCCounts,
    pub has_counts: bool,
    pub weight: usize,
    pub crossing_file: Option<PathBuf>,
}

impl WindowState {
    pub fn new(
        idx: u64,
        start: u64,
        end: u64,
        contained: bool,
        template: &GCCounts,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            idx,
            start,
            end,
            contained,
            counts: template.zeroed_like()?,
            has_counts: false,
            weight: 0,
            crossing_file: None,
        })
    }

    pub fn reset(
        &mut self,
        idx: u64,
        start: u64,
        end: u64,
        contained: bool,
        template: &GCCounts,
    ) -> anyhow::Result<()> {
        self.idx = idx;
        self.start = start;
        self.end = end;
        self.contained = contained;
        if self.counts.buffer_len() != template.buffer_len() {
            self.counts = template.zeroed_like()?;
        } else {
            self.counts.clear();
        }
        self.has_counts = false;
        self.weight = 0;
        self.crossing_file = None;
        Ok(())
    }
}

pub fn fixed_size_window_bounds(idx: u64, window_bp: u64, chrom_len: u64) -> (u64, u64) {
    let start = idx.saturating_mul(window_bp);
    let end = start.saturating_add(window_bp).min(chrom_len);
    (start, end)
}

#[inline]
pub fn window_state_from_idx(
    idx: u64,
    window_bp: u64,
    chrom_len: u64,
    core_start: u64,
    core_end: u64,
    template: &GCCounts,
) -> Result<WindowState> {
    let (start, end) = fixed_size_window_bounds(idx, window_bp, chrom_len);
    let contained = start >= core_start && end <= core_end;
    WindowState::new(idx, start, end, contained, template)
}

/// Compute the number of ACGT bases in the current window using prefix sums.
///
/// Uses the overlap between the window and the provided sequence to index the prefix
/// array at the window bounds, then subtracts those prefix values to recover the count
/// in O(1) without per-base scanning.
///
/// If a window is completely outside the available sequence (metadata shorter than the
/// tiling or truncated input), this throws an error instead of silently treating the
/// window as empty.
///
/// Parameters
/// ----------
/// - `buf`:
///     Window state updated in place with ACGT counts.
///
/// - `gc_prefixes`:
///     Prefix sums where each entry stores total ACGT up to that index.
///
/// - `seq_start`:
///     Genomic start coordinate associated with the prefix sums.
///
/// - `seq_end`:
///     Genomic end coordinate associated with the prefix sums.
///
/// Returns
/// -------
/// - Updates `buf.counts.num_acgt_out_of` with `(acgt_in_window, length_used)`.
/// - Returns an error when the window does not overlap the available sequence.
pub fn compute_window_acgt(
    buf: &mut WindowState,
    gc_prefixes: &GCPrefixes,
    seq_start: u64,
    seq_end: u64,
) -> Result<()> {
    let observed_start = buf.start.max(seq_start).min(seq_end);
    let observed_end = buf.end.min(seq_end);
    ensure!(
        observed_end > observed_start,
        "Window [{}, {}) does not overlap sequence [{}, {})",
        buf.start,
        buf.end,
        seq_start,
        seq_end
    );
    let start_local = (observed_start - seq_start) as usize;
    let end_local = (observed_end - seq_start) as usize;
    // Prefix arrays are sized to the loaded sequence range. A failure here means the
    // window's overlap extends beyond the available prefix data, which is a data error.
    ensure!(
        end_local < gc_prefixes.acgt.len(),
        "Window end index {} exceeds prefix length {}",
        end_local,
        gc_prefixes.acgt.len().saturating_sub(1)
    );
    let acgt_count = gc_prefixes.acgt[end_local] - gc_prefixes.acgt[start_local];
    let observed_len = (end_local - start_local) as u64;
    buf.counts.num_acgt_out_of = (acgt_count as u64, observed_len);
    Ok(())
}

pub fn overlap_length(a_start: u64, a_end: u64, b_start: u64, b_end: u64) -> u64 {
    let start = a_start.max(b_start);
    let end = a_end.min(b_end);
    end.saturating_sub(start)
}

#[derive(Debug)]
pub struct PreparedTileWindows {
    pub windows: Vec<WindowState>,
    pub streaming_buffers: Option<(u64, WindowState, WindowState)>,
    pub skip_tile: bool,
}

/// Prepares window buffers for a tile based on the configured window specification.
///
/// Builds either a vector of per-window states for BED or global windows, or two streaming buffers
/// for fixed-size windows so callers avoid reallocating buffers while scanning fragments. Returns
/// a skip flag when the chromosome carries no windows, allowing the caller to exit early without
/// extra work.
///
/// Parameters
/// ----------
/// - `window_opt`:
///     Window configuration indicating whether to use BED, fixed-size, or global windows.
/// - `windows_opt`:
///     Optional slice of per-chromosome windows `(start, end, idx)` when using BED windows.
/// - `tile`:
///     Tile whose core and fetch bounds determine which windows are relevant.
/// - `tile_window_span`:
///     Optional cached span that bounds candidate windows for the tile.
/// - `chrom_len`:
///     Total chromosome length used to cap fixed-size window ends.
/// - `template`:
///     GC count template used to initialize each window buffer.
///
/// Returns
/// -------
/// - `PreparedTileWindows`:
///     Contains the ready-to-use window buffers, any streaming pair, and a skip flag for empty BED
///     chromosomes.
pub fn prepare_tile_windows(
    window_opt: &WindowSpec,
    windows_opt: Option<&[(u64, u64, u64)]>,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    chrom_len: u64,
    template: &GCCounts,
) -> Result<PreparedTileWindows> {
    let mut windows: Vec<WindowState> = Vec::new();
    let mut streaming_buffers: Option<(u64, WindowState, WindowState)> = None;
    let mut skip_tile = false;

    let core_start = tile.core_start as u64;
    let core_end = tile.core_end as u64;

    match window_opt {
        WindowSpec::Bed(_) => {
            let win_slice = windows_opt.ok_or_else(|| {
                anyhow!(
                    "Window specification is windowed, but no windows provided for chromosome {}",
                    &tile.chr
                )
            })?;
            if win_slice.is_empty() {
                skip_tile = true;
            } else {
                let capacity = tile_window_span
                    .map(|span| span.last_idx_exclusive.saturating_sub(span.first_idx))
                    .unwrap_or(win_slice.len());
                windows.reserve(capacity);
                for &(window_start, window_end, window_idx) in
                    overlapping_windows_for_tile(win_slice, tile, tile_window_span)
                {
                    let contained = window_start >= core_start && window_end <= core_end;
                    windows.push(WindowState::new(
                        window_idx,
                        window_start,
                        window_end,
                        contained,
                        template,
                    )?);
                }
            }
        }
        WindowSpec::Size(window_bp) => {
            let current_idx = core_start / *window_bp;
            let current = window_state_from_idx(
                current_idx,
                *window_bp,
                chrom_len,
                core_start,
                core_end,
                template,
            )?;

            let next_idx = current_idx + 1;
            let next = window_state_from_idx(
                next_idx,
                *window_bp,
                chrom_len,
                core_start,
                core_end,
                template,
            )?;

            streaming_buffers = Some((*window_bp, current, next));
        }
        WindowSpec::Global => {
            windows.push(WindowState::new(0, core_start, core_end, true, template)?);
        }
    }

    Ok(PreparedTileWindows {
        windows,
        streaming_buffers,
        skip_tile,
    })
}
