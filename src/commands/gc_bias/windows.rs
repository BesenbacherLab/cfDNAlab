use crate::commands::cli_common::WindowSpec;
use crate::commands::gc_bias::counting::{GCCounts, GCPrefixes};
use crate::shared::{
    bam::Contigs,
    bed::Windows,
    interval::{IndexedInterval, Interval},
    tiled_run::{Tile, TileWindowSpan, overlapping_windows_for_tile},
};
use anyhow::{Context, Result, anyhow, bail, ensure};
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
            let map =
                windows_map.ok_or_else(|| anyhow!("BED window spec requires loaded windows"))?;
            let mut total_len: u64 = 0;
            let mut count: u64 = 0;
            for chr in chromosomes {
                if let Some(ws) = map.get(chr) {
                    for window in ws.as_slice() {
                        total_len =
                            total_len.saturating_add(window.end().saturating_sub(window.start()));
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

/// Mutable state for one GC-bias counting interval.
///
/// The struct holds the interval being accumulated, whether that interval is fully contained in
/// the tile core, the current counts, and bookkeeping used while windows are streamed, finalized,
/// and optionally spilled for later reduction.
#[derive(Clone, Debug)]
pub struct WindowState {
    // Stable window index in the current mode
    pub idx: u64,
    // Checked genomic interval represented by this state
    pub interval: Interval<u64>,
    // Whether the interval lies fully inside the tile core
    pub contained: bool,
    // Counts accumulated for this interval
    pub counts: GCCounts,
    // Whether any fragment has contributed counts to this state
    pub has_counts: bool,
    // Number of finalized contributions merged into this state
    pub weight: usize,
    // Optional path to spill data written for later reduction
    pub crossing_file: Option<PathBuf>,
}

impl WindowState {
    pub fn new(
        idx: u64,
        interval: Interval<u64>,
        contained: bool,
        template: &GCCounts,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            idx,
            interval,
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
        interval: Interval<u64>,
        contained: bool,
        template: &GCCounts,
    ) -> anyhow::Result<()> {
        self.idx = idx;
        self.interval = interval;
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

    #[inline]
    pub fn start(&self) -> u64 {
        self.interval.start()
    }

    #[inline]
    pub fn end(&self) -> u64 {
        self.interval.end()
    }
}

pub fn fixed_size_window_interval(
    idx: u64,
    window_bp: u64,
    chrom_len: u64,
) -> crate::Result<Interval<u64>> {
    let start = idx.saturating_mul(window_bp);
    if start >= chrom_len {
        return Err(crate::Error::InvalidFixedWindowIndex {
            idx,
            start,
            chrom_len,
        });
    }
    let end = start.saturating_add(window_bp).min(chrom_len);
    Interval::new(start, end)
}

#[inline]
pub fn window_state_from_idx(
    idx: u64,
    window_bp: u64,
    chrom_len: u64,
    core_interval: Interval<u64>,
    template: &GCCounts,
) -> Result<WindowState> {
    let interval = fixed_size_window_interval(idx, window_bp, chrom_len)?;
    let contained = core_interval.contains_interval(interval);
    WindowState::new(idx, interval, contained, template)
}

/// Build the rolling fixed-size window buffers for one tile.
///
/// Fixed-size GC bias counting keeps two buffers alive at a time: the current
/// window and the next window. This helper creates those checked buffers from
/// the checked tile core interval so callers do not repeat the index math.
///
/// The two-buffer setup works for most tiles because a fragment can only overlap
/// the current fixed window and the one immediately after it. Near chromosome
/// end, the "next" window may be the last partial window, and there may be no
/// valid non-empty window after that. In those cases, the missing next window
/// is represented as `None`.
fn prepare_fixed_size_streaming_buffers(
    window_bp: u64,
    chrom_len: u64,
    core_interval: Interval<u64>,
    template: &GCCounts,
) -> Result<(WindowState, Option<WindowState>)> {
    let current_idx = core_interval.start() / window_bp;
    let current =
        window_state_from_idx(current_idx, window_bp, chrom_len, core_interval, template)?;

    let next_idx = current_idx + 1;
    let next = if next_idx.saturating_mul(window_bp) < chrom_len {
        Some(window_state_from_idx(
            next_idx,
            window_bp,
            chrom_len,
            core_interval,
            template,
        )?)
    } else {
        None
    };

    Ok((current, next))
}

/// Advance the rolling fixed-size window buffers by one window.
///
/// After the old current window has been finalized, the old next window becomes
/// current and the recycled allocation is reset to the following window when it
/// exists.
///
/// This is the only place where the streaming path asks for the synthetic
/// window after the current `next` buffer. That means the chromosome-end edge
/// case also lives here: if `next` already points at the last partial window,
/// the following window should become `None` instead of constructing a fake
/// end-of-chromosome interval.
pub(crate) fn advance_fixed_size_streaming_buffers(
    current: WindowState,
    next: WindowState,
    window_bp: u64,
    chrom_len: u64,
    core_interval: Interval<u64>,
    template: &GCCounts,
) -> Result<(WindowState, Option<WindowState>)> {
    let mut recycled = current;
    let current = next;
    let next_idx = current.idx + 1;
    if next_idx.saturating_mul(window_bp) >= chrom_len {
        return Ok((current, None));
    }

    let next_interval = fixed_size_window_interval(next_idx, window_bp, chrom_len)?;
    let next_contained = core_interval.contains_interval(next_interval);
    recycled.reset(next_idx, next_interval, next_contained, template)?;
    Ok((current, Some(recycled)))
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
/// - `observed_interval`:
///     Genomic interval whose support should be measured for this window.
///
/// - `sequence_interval`:
///     Genomic interval associated with the prefix sums.
///
/// Returns
/// -------
/// - Updates `buf.counts.num_acgt_out_of` with `(acgt_in_window, length_used)`.
/// - Returns an error when the observed interval does not overlap the available sequence.
pub fn set_window_acgt_in_observed_interval(
    buf: &mut WindowState,
    gc_prefixes: &GCPrefixes,
    observed_interval: Interval<u64>,
    sequence_interval: Interval<u64>,
) -> Result<()> {
    let (seq_start, seq_end) = sequence_interval.as_tuple();
    let observed_interval = observed_interval
        .clip_to(sequence_interval)
        .ok_or_else(|| {
            anyhow!(
                "Observed interval [{}, {}) for window [{}, {}) does not overlap sequence [{}, {})",
                observed_interval.start(),
                observed_interval.end(),
                buf.start(),
                buf.end(),
                seq_start,
                seq_end
            )
        })?;
    let observed_local = observed_interval.shift_left(seq_start)?.try_to_usize()?;
    let acgt_count = gc_prefixes.acgt_count(observed_local).with_context(|| {
        format!(
            "counting ACGT support for observed interval [{}, {}) within window [{}, {}): \
                 sequence interval [{}, {}), local interval [{}, {}), prefix length {}",
            observed_interval.start(),
            observed_interval.end(),
            buf.start(),
            buf.end(),
            seq_start,
            seq_end,
            observed_local.start(),
            observed_local.end(),
            gc_prefixes.acgt.len().saturating_sub(1)
        )
    })?;
    let observed_len = observed_local.len() as u64;
    buf.counts.num_acgt_out_of = (acgt_count as u64, observed_len);
    Ok(())
}

pub fn overlap_length(a: Interval<u64>, b: Interval<u64>) -> u64 {
    a.intersection(b).map_or(0, |shared| shared.len())
}

#[derive(Debug)]
pub struct PreparedTileWindows {
    pub windows: Vec<WindowState>,
    pub streaming_buffers: Option<(u64, WindowState, Option<WindowState>)>,
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
    windows_opt: Option<&[IndexedInterval<u64>]>,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    chrom_len: u64,
    template: &GCCounts,
) -> Result<PreparedTileWindows> {
    let mut windows: Vec<WindowState> = Vec::new();
    let mut streaming_buffers: Option<(u64, WindowState, Option<WindowState>)> = None;
    let mut skip_tile = false;

    let core_interval = tile.core.try_to_u64()?;
    let core_start = tile.core_start() as u64;
    let core_end = tile.core_end() as u64;

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
                for window in overlapping_windows_for_tile(win_slice, tile, tile_window_span) {
                    let window_start = window.start();
                    let window_end = window.end();
                    let window_idx = window.idx();
                    let contained = window_start >= core_start && window_end <= core_end;
                    windows.push(WindowState::new(
                        window_idx,
                        window.interval,
                        contained,
                        template,
                    )?);
                }
            }
        }
        WindowSpec::Size(window_bp) => {
            let (current, next) = prepare_fixed_size_streaming_buffers(
                *window_bp,
                chrom_len,
                core_interval,
                template,
            )?;
            streaming_buffers = Some((*window_bp, current, next));
        }
        WindowSpec::Global => {
            windows.push(WindowState::new(
                0,
                Interval::new(core_start, core_end)?,
                true,
                template,
            )?);
        }
    }

    Ok(PreparedTileWindows {
        windows,
        streaming_buffers,
        skip_tile,
    })
}

#[cfg(test)]
mod tests {
    include!("windows_tests.rs");
}
