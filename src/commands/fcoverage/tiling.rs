use std::io::{BufRead, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::{
    commands::cli_common::WindowSpec,
    commands::fcoverage::window_results::CoverageWindowAction,
    shared::coverage::Coverage,
    shared::formatters::{CompactNumber, round_to},
    shared::interval::Interval,
    shared::io::open_text_reader,
    shared::tiled_run::{Tile, TileMode, TileWindowSpan, clamp_fetch_to_window_span},
    shared::window_fetch::{
        BedFetchPolicy, fetch_span_for_tile, full_tile_fetch_span,
        window_derived_fetch_extent_for_core_overlap,
    },
    shared::writers::open_zstd_auto_writer,
};

/// Kind of simple tile temp file returned by tile processing.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub(crate) enum TileTempFileKind {
    Positional,
    SizeFinal,
}

/// Returned path for one simple tile output.
///
/// The merge code uses this explicit path instead of inferring files from temporary directory
/// names. `kind` keeps positional outputs and already-finalized fixed-size outputs distinct while
/// sharing the same carrier fields. Final order is still requested chromosome order, then tile
/// index within chromosome.
#[derive(Debug, Clone)]
pub(crate) struct TileTempFile {
    pub kind: TileTempFileKind,
    pub chromosome: String,
    pub tile_index: u32,
    pub path: PathBuf,
}

fn sorted_tile_outputs_for_chromosome<'a>(
    tile_outputs: &'a [TileTempFile],
    kind: TileTempFileKind,
    chromosome: &str,
) -> Vec<&'a TileTempFile> {
    let mut chromosome_outputs = tile_outputs
        .iter()
        .filter(|output| output.kind == kind && output.chromosome == chromosome)
        .collect::<Vec<_>>();
    chromosome_outputs.sort_by_key(|output| output.tile_index);
    chromosome_outputs
}

/// Concatenates returned positional tile outputs into one final positional file.
///
/// The function streams tile frames in requested chromosome order and tile-index order. It never
/// scans the temp directory, so decoys or stale files with matching names cannot affect the final
/// output. Without restore-mean scaling, compressed tile frames are copied byte-for-byte.
///
/// Parameters
/// ----------
/// - `out_dir`:
///     Directory where the merged file is written.
/// - `chromosomes`:
///     Requested chromosome order for the final output.
/// - `tile_outputs`:
///     Returned positional tile paths from tile processing.
/// - `final_name`:
///     Filename for the merged output.
///
/// Returns
/// -------
/// - `PathBuf`:
///     Path to the merged final output.
fn merge_positional_tile_outputs(
    out_dir: &std::path::Path,
    chromosomes: &[String],
    tile_outputs: &[TileTempFile],
    final_name: &str,
) -> Result<std::path::PathBuf> {
    let final_path = out_dir.join(final_name);
    let mut out = BufWriter::new(
        std::fs::File::create(&final_path)
            .with_context(|| format!("Creating merged output: {}", final_path.display()))?,
    );

    for chromosome in chromosomes {
        for output in sorted_tile_outputs_for_chromosome(
            tile_outputs,
            TileTempFileKind::Positional,
            chromosome,
        ) {
            let mut tile_file = std::fs::File::open(&output.path)
                .with_context(|| format!("Opening tile file: {}", output.path.display()))?;
            std::io::copy(&mut tile_file, &mut out).with_context(|| {
                format!(
                    "Copying from {} into {}",
                    output.path.display(),
                    final_path.display()
                )
            })?;
        }
    }

    out.flush().context("Flushing merged output")?;
    Ok(final_path)
}

fn merge_scaled_positional_tile_outputs(
    out_dir: &std::path::Path,
    chromosomes: &[String],
    tile_outputs: &[TileTempFile],
    final_name: &str,
    multiplier: f64,
    indexed: bool,
    decimals: i32,
    n_threads: usize,
) -> Result<std::path::PathBuf> {
    let final_path = out_dir.join(final_name);
    let mut out = open_zstd_auto_writer(&final_path, 3, Some(n_threads as u32))?;
    let mut line = String::new();

    for chromosome in chromosomes {
        for output in sorted_tile_outputs_for_chromosome(
            tile_outputs,
            TileTempFileKind::Positional,
            chromosome,
        ) {
            let mut reader = open_text_reader(&output.path)?;
            loop {
                line.clear();
                if reader.read_line(&mut line)? == 0 {
                    break;
                }

                let raw = line.trim_end_matches('\n').trim_end_matches('\r');
                if raw.is_empty() {
                    continue;
                }

                let mut cols = raw.split('\t');
                let chr_col = cols.next().ok_or_else(|| {
                    anyhow::anyhow!("Missing chromosome column in {}", output.path.display())
                })?;
                let start_col = cols.next().ok_or_else(|| {
                    anyhow::anyhow!("Missing start column in {}", output.path.display())
                })?;
                let end_col = cols.next().ok_or_else(|| {
                    anyhow::anyhow!("Missing end column in {}", output.path.display())
                })?;
                let value_col = cols.next().ok_or_else(|| {
                    anyhow::anyhow!("Missing value column in {}", output.path.display())
                })?;
                let value = value_col.parse::<f64>().with_context(|| {
                    format!(
                        "Parsing positional value '{}' in {}",
                        value_col,
                        output.path.display()
                    )
                })?;
                let scaled_value = round_to(value * multiplier, decimals);

                if indexed {
                    let idx_col = cols.next().ok_or_else(|| {
                        anyhow::anyhow!("Missing window index column in {}", output.path.display())
                    })?;
                    anyhow::ensure!(
                        cols.next().is_none(),
                        "Unexpected extra columns in indexed positional tile {}",
                        output.path.display()
                    );
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{}",
                        chr_col,
                        start_col,
                        end_col,
                        CompactNumber {
                            v: scaled_value,
                            decimals
                        },
                        idx_col
                    )?;
                } else {
                    anyhow::ensure!(
                        cols.next().is_none(),
                        "Unexpected extra columns in positional tile {}",
                        output.path.display()
                    );
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}",
                        chr_col,
                        start_col,
                        end_col,
                        CompactNumber {
                            v: scaled_value,
                            decimals
                        }
                    )?;
                }
            }
        }
    }

    out.flush()
        .context("Flushing scaled merged positional output")?;
    Ok(final_path)
}

/// Merge positional tile outputs, optionally applying the late restore-mean multiplier.
///
/// Restore-mean is a final merge concern because the multiplier is only known after all tiles have
/// been counted. Row identity, chromosome order, tile order, and indexed positional columns stay
/// the same as the unscaled positional merge.
pub(crate) fn merge_positional_tile_outputs_with_optional_scaling(
    out_dir: &std::path::Path,
    chromosomes: &[String],
    tile_outputs: &[TileTempFile],
    final_name: &str,
    restore_mean_multiplier: Option<f64>,
    indexed: bool,
    decimals: i32,
    n_threads: usize,
) -> Result<std::path::PathBuf> {
    if let Some(multiplier) = restore_mean_multiplier {
        merge_scaled_positional_tile_outputs(
            out_dir,
            chromosomes,
            tile_outputs,
            final_name,
            multiplier,
            indexed,
            decimals,
            n_threads,
        )
    } else {
        merge_positional_tile_outputs(out_dir, chromosomes, tile_outputs, final_name)
    }
}

/// Joins returned aligned fixed-size final outputs while preserving compressed frame boundaries.
///
/// A compressed header frame is written first, followed by each returned tile frame in genomic
/// order. The tile payloads are copied verbatim so the final file stays a valid zstd frame
/// concatenation and does not re-derive values that were already finalized during tile processing.
///
/// Parameters
/// ----------
/// - `out_dir`:
///     Directory where the merged file is written.
/// - `chromosomes`:
///     Requested chromosome order for the final output.
/// - `tile_outputs`:
///     Returned final tile paths from aligned fixed-size tile processing.
/// - `final_name`:
///     Filename for the merged artifact.
/// - `header_line`:
///     Plain-text header to encode as its own compressed frame.
///
/// Returns
/// -------
/// - `PathBuf`:
///     Path to the merged final output.
pub(crate) fn concat_aligned_size_tile_final_outputs(
    out_dir: &std::path::Path,
    chromosomes: &[String],
    tile_outputs: &[TileTempFile],
    final_name: &str,
    header_line: &str,
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

    for chromosome in chromosomes {
        for output in sorted_tile_outputs_for_chromosome(
            tile_outputs,
            TileTempFileKind::SizeFinal,
            chromosome,
        ) {
            let mut tile_file = std::fs::File::open(&output.path)
                .with_context(|| format!("Opening {}", output.path.display()))?;
            std::io::copy(&mut tile_file, &mut out).with_context(|| {
                format!(
                    "Copying {} into {}",
                    output.path.display(),
                    final_path.display()
                )
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
    // then intersect it with the tile's existing fetch band.
    match mode {
        TileMode::Positional { windows: None, .. } => full_tile_fetch_span(tile, chrom_len_u64),
        TileMode::Positional {
            windows: Some(wchr),
            ..
        } => {
            let Some(window_span) =
                window_derived_fetch_extent_for_core_overlap(wchr, tile, tile_span)?
            else {
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
            let Some(window_span) =
                window_derived_fetch_extent_for_core_overlap(wchr, tile, tile_span)?
            else {
                return Ok(None);
            };
            Ok(clamp_fetch_to_window_span(
                tile,
                chrom_len_u64,
                window_span,
                halo_bp,
            )?)
        }
        TileMode::AggregatesBySize { window_bp, .. } => fetch_span_for_tile(
            tile,
            None,
            None,
            &WindowSpec::Size(*window_bp),
            chrom_len_u64,
            halo_bp,
            BedFetchPolicy::CoreOverlap,
        ),
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

/// Prefix arrays used to derive summary statistics over coverage slices.
#[derive(Debug, Clone)]
pub struct CoverageSummaryPrefixes {
    pub sum_of_squares_all: Vec<f64>,
    pub sum_of_squares_unmasked: Option<Vec<f64>>,
    pub nonzero_all: Vec<u64>,
    pub nonzero_unmasked: Option<Vec<u64>>,
}

pub fn build_summary_prefixes(cp: &Coverage) -> Result<CoverageSummaryPrefixes> {
    let coverage = cp.coverage().ok_or_else(|| {
        anyhow::anyhow!("coverage must remain available while building summary-stat prefixes")
    })?;
    let mask = cp.blacklist_mask();

    let mut sum_of_squares_all = Vec::with_capacity(coverage.len() + 1);
    let mut nonzero_all = Vec::with_capacity(coverage.len() + 1);
    let mut sum_of_squares_unmasked = None;
    let mut nonzero_unmasked = None;

    sum_of_squares_all.push(0.0);
    nonzero_all.push(0);

    let mut running_sum_of_squares_all = 0.0;
    let mut running_nonzero_all = 0u64;

    if let Some(mask) = mask {
        let mut sum_of_squares_unmasked_prefix = Vec::with_capacity(coverage.len() + 1);
        let mut nonzero_unmasked_prefix = Vec::with_capacity(coverage.len() + 1);
        sum_of_squares_unmasked_prefix.push(0.0);
        nonzero_unmasked_prefix.push(0);

        let mut running_sum_of_squares_unmasked = 0.0;
        let mut running_nonzero_unmasked = 0u64;

        for (index, &value_f32) in coverage.iter().enumerate() {
            let value = value_f32 as f64;
            let squared_value = value * value;
            let is_nonzero = value > 0.0;

            running_sum_of_squares_all += squared_value;
            if is_nonzero {
                running_nonzero_all += 1;
            }

            if mask[index] == 0 {
                running_sum_of_squares_unmasked += squared_value;
                if is_nonzero {
                    running_nonzero_unmasked += 1;
                }
            }

            sum_of_squares_all.push(running_sum_of_squares_all);
            nonzero_all.push(running_nonzero_all);
            sum_of_squares_unmasked_prefix.push(running_sum_of_squares_unmasked);
            nonzero_unmasked_prefix.push(running_nonzero_unmasked);
        }

        sum_of_squares_unmasked = Some(sum_of_squares_unmasked_prefix);
        nonzero_unmasked = Some(nonzero_unmasked_prefix);
    } else {
        for &value_f32 in coverage {
            let value = value_f32 as f64;
            let squared_value = value * value;

            running_sum_of_squares_all += squared_value;
            if value > 0.0 {
                running_nonzero_all += 1;
            }

            sum_of_squares_all.push(running_sum_of_squares_all);
            nonzero_all.push(running_nonzero_all);
        }
    }

    Ok(CoverageSummaryPrefixes {
        sum_of_squares_all,
        sum_of_squares_unmasked,
        nonzero_all,
        nonzero_unmasked,
    })
}

/// Raw summary statistics over one tile-local slice.
#[derive(Debug, Clone, Copy)]
pub struct CoverageSummarySliceStats {
    pub coverage_sum: f64,
    pub eligible_positions: u64,
    pub blacklisted_positions: u64,
    pub nonzero_positions: u64,
    pub coverage_sum_of_squares: f64,
}

/// Compute raw summary statistics over a tile-local range.
///
/// This extends `coverage_sum_and_counts` with the second raw moment and the count of eligible
/// nonzero bases so later reducers can derive variance, SD, coverage fraction, and grouped
/// correlations without revisiting per-base coverage.
#[inline]
pub fn coverage_summary_and_counts(
    local_start_idx: usize,
    local_end_idx: usize,
    masked: bool,
    ps_all: &[f64],
    ps_allow: Option<&[f64]>,
    cnt_allow: Option<&[u32]>,
    mask: Option<&[u8]>,
    summary_prefixes: &CoverageSummaryPrefixes,
) -> CoverageSummarySliceStats {
    let (coverage_sum, eligible_positions, blacklisted_positions) = coverage_sum_and_counts(
        local_start_idx,
        local_end_idx,
        masked,
        ps_all,
        ps_allow,
        cnt_allow,
        mask,
    );

    let coverage_sum_of_squares = if masked {
        if let Some(prefix) = summary_prefixes.sum_of_squares_unmasked.as_ref() {
            prefix[local_end_idx] - prefix[local_start_idx]
        } else {
            summary_prefixes.sum_of_squares_all[local_end_idx]
                - summary_prefixes.sum_of_squares_all[local_start_idx]
        }
    } else {
        summary_prefixes.sum_of_squares_all[local_end_idx]
            - summary_prefixes.sum_of_squares_all[local_start_idx]
    };

    let nonzero_positions = if masked {
        if let Some(prefix) = summary_prefixes.nonzero_unmasked.as_ref() {
            prefix[local_end_idx] - prefix[local_start_idx]
        } else {
            summary_prefixes.nonzero_all[local_end_idx]
                - summary_prefixes.nonzero_all[local_start_idx]
        }
    } else {
        summary_prefixes.nonzero_all[local_end_idx] - summary_prefixes.nonzero_all[local_start_idx]
    };

    CoverageSummarySliceStats {
        coverage_sum,
        eligible_positions,
        blacklisted_positions,
        nonzero_positions,
        coverage_sum_of_squares,
    }
}

/// Converts accumulated coverage statistics into a final window value.
///
/// Depending on the requested action the result is either an average over eligible positions or
/// the raw total. Averages with no eligible positions are undefined and return `NaN`.
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
        CoverageWindowAction::Average | CoverageWindowAction::AverageOnUniqueBases => {
            if masked {
                if allowed_positions == 0 {
                    f64::NAN
                } else {
                    sum / allowed_positions as f64
                }
            } else {
                if unmasked_span_bp == 0 {
                    f64::NAN
                } else {
                    sum / unmasked_span_bp as f64
                }
            }
        }
        CoverageWindowAction::Total | CoverageWindowAction::TotalOnUniqueBases => sum,
        CoverageWindowAction::SummaryStats | CoverageWindowAction::SummaryStatsOnUniqueBases => {
            unreachable!("summary-stats uses raw aggregate rows instead of finalize_value()")
        }
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

#[cfg(test)]
mod tests {
    include!("tiling_tests.rs");
}
