use std::str::FromStr;

use anyhow::Result;

use crate::shared::coverage::Coverage;

/// What to do per window
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PeaksWindowAction {
    #[default]
    Stats,
    OnlyIncludeThesePositionsUnique,
    OnlyIncludeThesePositionsIndexed,
}

// For the CLI
impl FromStr for PeaksWindowAction {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "stats" {
            Ok(PeaksWindowAction::Stats)
        } else if s == "unique-positions" {
            Ok(PeaksWindowAction::OnlyIncludeThesePositionsUnique)
        } else if s == "indexed-positions" {
            Ok(PeaksWindowAction::OnlyIncludeThesePositionsIndexed)
        } else {
            Err("Use 'stats', 'indexed-positions', or 'unique-positions'".into())
        }
    }
}

/// Collection of statistics about WPS peaks.
#[derive(Debug, Clone)]
pub struct PeakStats {
    pub count: u32,
    pub avg_distance: f32,
    pub median_distance: u32,
}

/// Per-window payload
#[derive(Debug, Clone)]
pub enum WindowPeaksValue {
    /// Average coverage in the window
    Stats(PeakStats),
    /// Positional coverage for every base in the window, left->right
    Positions(Vec<f32>),
}

/// One window's result (keeps original ordering info)
#[derive(Debug, Clone)]
pub struct WindowPeaksResult {
    pub start: u64,
    pub end: u64,
    pub original_idx: u64,
    pub value: WindowPeaksValue,
    pub num_blacklisted_pos: Option<u32>,
}

/// Top-level result for a run with or without windows
#[derive(Debug, Clone)]
pub enum PeaksOutput {
    /// Results for each input window
    PerWindow {
        action: PeaksWindowAction,
        results: Vec<WindowPeaksResult>,
    },
    /// No windows given -> return positional peaks for the whole sequence
    WholePositional {
        /// Start offset, typically 0
        start: u64,
        /// End offset, typically `length`
        end: u64,
        /// Per-base coverage, left->right
        values: Vec<f32>,
    },
}

// TODO: This is not implemented
/// Compute peak outputs for windows or whole-chromosome positions
///
/// Parameters
/// ----------
/// - cp: Coverage with coverage finalized and indexes buildable
/// - windows: Optional triplets of `(start, end, original_idx)`
/// - action: What to return per window
/// - nan_blacklisted: Set blacklisted positions to `f32::NAN` and exclude when computing sums/averages
///
/// Returns
/// -------
/// - out: `WindowPeaksResult` with either per-window results or whole positional coverage
pub fn compute_window_peaks_outputs(
    cp: &mut Coverage, // TODO: Assume not here?
    windows: Option<&[(u64, u64, u64)]>,
    action: PeaksWindowAction,
    nan_blacklisted: bool,
) -> Result<PeaksOutput> {
    // Require finalized coverage (do not finalize here to keep behavior explicit)
    if cp.coverage().is_none() {
        anyhow::bail!(
            "coverage not finalized; call finalize_coverage() before compute_window_outputs"
        )
    }

    // No windows (None or empty) -> positional coverage for entire sequence
    if windows.map_or(true, |w| w.is_empty()) {
        let cov = cp.coverage_in_window(0, cp.length(), nan_blacklisted)?;
        return Ok(PeaksOutput::WholePositional {
            start: 0,
            end: cp.length() as u64,
            values: cov.to_vec(),
        });
    }

    // Safe unwrap after empty-check
    let windows = windows.unwrap();

    // Bounds check once up front
    let len_u64 = cp.length() as u64;
    for &(s, e, _) in windows {
        if s > e || e > len_u64 {
            anyhow::bail!("window [{s}..{e}) out of bounds for length {len_u64}");
        }
    }

    match action {
        PeaksWindowAction::Average => {
            // Build (or reuse) indexes explicitly for clarity
            cp.build_indexes(true)?;

            let spans: Vec<(u32, u32)> = windows
                .iter()
                .map(|&(s, e, _)| (s as u32, e as u32))
                .collect();
            let avgs = cp.bulk_avg_coverage(&spans, nan_blacklisted, false)?;

            let mut results = Vec::with_capacity(windows.len());
            for (&(s, e, idx), &avg) in windows.iter().zip(avgs.iter()) {
                let bl = cp.blacklist_mask().map(|mask| {
                    let a = s as usize;
                    let b = e as usize;
                    mask[a..b].iter().map(|&m| (m == 1) as u32).sum()
                });
                results.push(WindowPeaksResult {
                    start: s,
                    end: e,
                    original_idx: idx,
                    value: WindowPeaksValue::Average(avg),
                    num_blacklisted_pos: bl,
                });
            }

            Ok(PeaksOutput::PerWindow { action, results })
        }
        PeaksWindowAction::Total => {
            cp.build_indexes(true)?;
            let spans: Vec<(u32, u32)> = windows
                .iter()
                .map(|&(s, e, _)| (s as u32, e as u32))
                .collect();
            let sums = cp.bulk_sum_coverage(&spans, nan_blacklisted, false)?;

            let mut results = Vec::with_capacity(windows.len());
            for (&(s, e, idx), &sum) in windows.iter().zip(sums.iter()) {
                let bl = cp.blacklist_mask().map(|mask| {
                    let a = s as usize;
                    let b = e as usize;
                    mask[a..b].iter().map(|&m| (m == 1) as u32).sum()
                });
                results.push(WindowPeaksResult {
                    start: s,
                    end: e,
                    original_idx: idx,
                    value: WindowPeaksValue::Total(sum),
                    num_blacklisted_pos: bl,
                });
            }

            Ok(PeaksOutput::PerWindow { action, results })
        }
        PeaksWindowAction::OnlyIncludeThesePositionsUnique
        | PeaksWindowAction::OnlyIncludeThesePositionsIndexed => {
            // Positional coverage per window, with optional NaNs for blacklisted sites
            let cov = cp.coverage_in_window(0, cp.length(), nan_blacklisted)?;

            let mut results = Vec::with_capacity(windows.len());
            for &(s, e, idx) in windows {
                let a = s as usize;
                let b = e as usize;

                // If you need to exclude blacklisted positions here, map through `bl_mask`
                // For now we return all positions in-window
                let vals = cov[a..b].to_vec();

                results.push(WindowPeaksResult {
                    start: s,
                    end: e,
                    original_idx: idx,
                    value: WindowPeaksValue::Positions(vals),
                    num_blacklisted_pos: None,
                });
            }

            Ok(PeaksOutput::PerWindow { action, results })
        }
    }
}
