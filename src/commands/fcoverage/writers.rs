use crate::shared::formatters::{CompactNumber, round_to_with_precomputed_factor};
use crate::shared::interval::Interval;
use anyhow::Result;
use std::io::Write;

/// Write a final aggregate row: `chromosome  start  end  value  blacklisted_positions`
#[inline]
pub fn write_final_row<W: Write>(
    w: &mut W,
    chr: &str,
    interval: Interval<u64>,
    value: f64,
    blacklisted_positions: u64,
    decimals: i32,
) -> anyhow::Result<()> {
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}",
        chr,
        interval.start(),
        interval.end(),
        CompactNumber { v: value, decimals },
        blacklisted_positions
    )?;
    Ok(())
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
/// - `keep_zero_runs`: Whether runs with values that rounds to zero should still be written.
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
            let run_start_abs = tile_abs_start + run_lo as u64;
            let run_end_abs = tile_abs_start + run_hi as u64;
            // Ignore write errors here; bubbled up by caller on flush
            let _ = writeln!(
                out,
                "{}\t{}\t{}\t{}",
                chr,
                run_start_abs,
                run_end_abs,
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
/// - `keep_zero_runs`: Whether runs with values that rounds to zero should still be written.
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
    orig_idx: Option<u64>,  // Window's original index
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
            let run_start_abs = tile_abs_start + run_lo as u64;
            let run_end_abs = tile_abs_start + run_hi as u64;
            let _ = if let Some(idx) = orig_idx {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    chr,
                    run_start_abs,
                    run_end_abs,
                    CompactNumber { v: value, decimals },
                    idx
                )
            } else {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}",
                    chr,
                    run_start_abs,
                    run_end_abs,
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
