use std::str::FromStr;

use anyhow::Result;

use crate::shared::coverage::Coverage;

/// What to do per window
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CoverageWindowAction {
    #[default]
    Average,
    Total,
    OnlyIncludeThesePositionsUnique,
    OnlyIncludeThesePositionsIndexed,
}

// For the CLI
impl FromStr for CoverageWindowAction {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "average" {
            Ok(CoverageWindowAction::Average)
        } else if s == "total" {
            Ok(CoverageWindowAction::Total)
        } else if s == "unique-positions" {
            Ok(CoverageWindowAction::OnlyIncludeThesePositionsUnique)
        } else if s == "indexed-positions" {
            Ok(CoverageWindowAction::OnlyIncludeThesePositionsIndexed)
        } else {
            Err("Use 'average', 'total', 'indexed-positions', or 'unique-positions'".into())
        }
    }
}

/// Per-window payload
#[derive(Debug, Clone)]
pub enum WindowValue {
    /// Average coverage in the window
    Average(f32),
    /// Total coverage in the window
    Total(f64),
    /// Positional coverage for every base in the window, left->right
    Positions(Vec<f32>),
}

/// One window's result (keeps original ordering info)
#[derive(Debug, Clone)]
pub struct WindowResult {
    pub start: u64,
    pub end: u64,
    pub original_idx: u64,
    pub value: WindowValue,
    pub num_blacklisted_pos: Option<u32>,
}

/// Top-level result for a run with or without windows
#[derive(Debug, Clone)]
pub enum CoverageOutput {
    /// Results for each input window
    PerWindow {
        action: CoverageWindowAction,
        results: Vec<WindowResult>,
    },
    /// No windows given -> return positional coverage for the whole sequence
    WholePositional {
        /// Start offset, typically 0
        start: u64,
        /// End offset, typically `length`
        end: u64,
        /// Per-base coverage, left->right
        values: Vec<f32>,
    },
}

/// Compute outputs for windows or whole-chromosome positions
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
/// - out: `CoverageOutput` with either per-window results or whole positional coverage
pub fn compute_window_outputs(
    cp: &mut Coverage,
    windows: Option<&[(u64, u64, u64)]>,
    action: CoverageWindowAction,
    nan_blacklisted: bool,
) -> Result<CoverageOutput> {
    // Require finalized coverage (do not finalize here to keep behavior explicit)
    if cp.coverage().is_none() {
        anyhow::bail!(
            "coverage not finalized; call finalize_coverage() before compute_window_outputs"
        )
    }

    // No windows (None or empty) -> positional coverage for entire sequence
    if windows.is_none_or(|w| w.is_empty()) {
        let cov = cp.coverage_in_window(0, cp.length(), nan_blacklisted)?;
        return Ok(CoverageOutput::WholePositional {
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
        CoverageWindowAction::Average => {
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
                results.push(WindowResult {
                    start: s,
                    end: e,
                    original_idx: idx,
                    value: WindowValue::Average(avg),
                    num_blacklisted_pos: bl,
                });
            }

            Ok(CoverageOutput::PerWindow { action, results })
        }
        CoverageWindowAction::Total => {
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
                results.push(WindowResult {
                    start: s,
                    end: e,
                    original_idx: idx,
                    value: WindowValue::Total(sum),
                    num_blacklisted_pos: bl,
                });
            }

            Ok(CoverageOutput::PerWindow { action, results })
        }
        CoverageWindowAction::OnlyIncludeThesePositionsUnique
        | CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
            // Positional coverage per window, with optional NaNs for blacklisted sites
            let cov = cp.coverage_in_window(0, cp.length(), nan_blacklisted)?;

            let mut results = Vec::with_capacity(windows.len());
            for &(s, e, idx) in windows {
                let a = s as usize;
                let b = e as usize;

                // If you need to exclude blacklisted positions here, map through `bl_mask`
                // For now we return all positions in-window
                let vals = cov[a..b].to_vec();

                results.push(WindowResult {
                    start: s,
                    end: e,
                    original_idx: idx,
                    value: WindowValue::Positions(vals),
                    num_blacklisted_pos: None,
                });
            }

            Ok(CoverageOutput::PerWindow { action, results })
        }
    }
}
