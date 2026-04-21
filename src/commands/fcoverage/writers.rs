use crate::shared::formatters::{CompactNumber, round_to_with_precomputed_factor};
use crate::shared::interval::Interval;
use crate::shared::writers::open_zstd_auto_writer;
use anyhow::{Result, bail};
use fxhash::FxHashMap;
use std::fmt;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::commands::fcoverage::reducer::{
    reduce_aggregates_by_size_with_cross_index_for_chr,
    reduce_bed_basic_with_cross_index_for_chr_rows,
    reduce_aggregates_by_size_with_cross_index_for_chr_rows, reduce_bed_with_cross_index_for_chr,
    reduce_bed_with_cross_index_for_chr_rows,
};
use crate::commands::fcoverage::tiling::concat_aligned_size_tile_finals;
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use crate::shared::base::{ZEROISH_F32_TOLERANCE, ZEROISH_F64_TOLERANCE};
use crate::shared::bed::{GroupedCoverageLayout, Windows};

/// Treat an `f64` mean this close to zero as numerically zero for CV reporting.
///
/// This is intentionally a numeric cleanup threshold, not a domain/reporting threshold.
/// If we ever want to suppress CV at larger low-coverage means, that should be a separate,
/// explicitly named policy.
const ZEROISH_COVERAGE_MEAN_F64: f64 = ZEROISH_F64_TOLERANCE;

/// Largest finite CV that we still print as a raw number in summary-stat outputs.
///
/// A CV above this bound is still mathematically finite, but the exact value is usually not
/// informative in a compact TSV. We therefore render it as `>1e6` to make it obvious that the
/// relative dispersion was extreme without pretending that the exact magnitude is useful.
const MAX_DISPLAYABLE_COVERAGE_CV: f64 = 1.0e6;

#[derive(Debug, Clone, Copy, Default)]
struct GroupedAggregateAccum {
    span_positions: u64,
    blacklisted_positions: u64,
    eligible_positions: u64,
    nonzero_positions: u64,
    coverage_sum: f64,
    coverage_sum_of_squares: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct SummaryStatsRow {
    pub span_positions: u64,
    pub blacklisted_positions: u64,
    pub eligible_positions: u64,
    pub nonzero_positions: u64,
    pub coverage_sum: f64,
    pub coverage_sum_of_squares: f64,
    pub mean_coverage: f64,
    pub total_coverage: f64,
    pub variance_coverage: f64,
    pub sd_coverage: f64,
    pub coefficient_of_variation_coverage: f64,
    pub covered_fraction: f64,
}

/// Formatter for the CV column only.
///
/// Most summary-stat fields should stay fully numeric, but CV has a special readability issue:
/// when the mean is extremely small, the exact finite ratio can become so large that it no
/// longer helps humans compare windows. This formatter keeps ordinary values numeric and only
/// switches to the sentinel string `>1e6` for truly extreme finite CV values.
struct CoverageCoefficientOfVariation {
    value: f64,
    decimals: i32,
}

impl fmt::Display for CoverageCoefficientOfVariation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.value.is_finite() && self.value > MAX_DISPLAYABLE_COVERAGE_CV {
            return f.write_str(&coverage_cv_overflow_label());
        }
        CompactNumber {
            v: self.value,
            decimals: self.decimals,
        }
        .fmt(f)
    }
}

/// Build the sentinel label used when the CV exceeds the display threshold.
///
/// The label is derived from the same constant that drives the numeric comparison so the
/// printed text cannot silently drift away from the active threshold. If the threshold
/// constant is ever changed to an invalid value, we panic rather than write a misleading label.
fn coverage_cv_overflow_label() -> String {
    if !MAX_DISPLAYABLE_COVERAGE_CV.is_finite() || MAX_DISPLAYABLE_COVERAGE_CV <= 0.0 {
        panic!(
            "MAX_DISPLAYABLE_COVERAGE_CV must be finite and > 0, got {}",
            MAX_DISPLAYABLE_COVERAGE_CV
        );
    }

    let scientific = format!("{:e}", MAX_DISPLAYABLE_COVERAGE_CV);
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("scientific notation formatter must include an exponent");
    let exponent_value: i32 = exponent
        .parse()
        .expect("scientific notation exponent must parse as i32");

    let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
    format!(">{}e{}", mantissa, exponent_value)
}

/// Build derived summary statistics from exact raw aggregate values.
pub fn derive_summary_stats(
    span_positions: u64,
    blacklisted_positions: u64,
    eligible_positions: u64,
    nonzero_positions: u64,
    coverage_sum: f64,
    coverage_sum_of_squares: f64,
) -> Result<SummaryStatsRow> {
    let total_coverage = coverage_sum;

    if eligible_positions == 0 {
        return Ok(SummaryStatsRow {
            span_positions,
            blacklisted_positions,
            eligible_positions,
            nonzero_positions,
            coverage_sum: total_coverage,
            coverage_sum_of_squares,
            mean_coverage: f64::NAN,
            total_coverage,
            variance_coverage: f64::NAN,
            sd_coverage: f64::NAN,
            coefficient_of_variation_coverage: f64::NAN,
            covered_fraction: f64::NAN,
        });
    }

    let eligible_positions_f64 = eligible_positions as f64;
    let mean_coverage = coverage_sum / eligible_positions_f64;
    let variance_coverage = derive_nonnegative_variance_coverage(
        eligible_positions,
        coverage_sum_of_squares,
        mean_coverage,
    )?;
    let sd_coverage = variance_coverage.sqrt();
    let coefficient_of_variation_coverage =
        derive_coefficient_of_variation_coverage(mean_coverage, sd_coverage);
    let covered_fraction = nonzero_positions as f64 / eligible_positions_f64;

    Ok(SummaryStatsRow {
        span_positions,
        blacklisted_positions,
        eligible_positions,
        nonzero_positions,
        coverage_sum,
        coverage_sum_of_squares,
        mean_coverage,
        total_coverage,
        variance_coverage,
        sd_coverage,
        coefficient_of_variation_coverage,
        covered_fraction,
    })
}

/// Derive a CV value that is safe and readable to report.
///
/// The CV itself is `sd / mean`, so it is mathematically undefined when the mean is zero.
/// We also treat means within `2 * f64::EPSILON` as zero because values that small are not a
/// meaningful coverage signal and should not create enormous ratios from denominator noise.
///
/// This helper does not cap large finite CV values. It preserves the real numeric result and
/// leaves the `>1e6` presentation policy to the dedicated display wrapper used only when
/// writing the TSV rows.
fn derive_coefficient_of_variation_coverage(mean_coverage: f64, sd_coverage: f64) -> f64 {
    if !mean_coverage.is_finite() || mean_coverage.abs() <= ZEROISH_COVERAGE_MEAN_F64 {
        f64::NAN
    } else {
        sd_coverage / mean_coverage
    }
}

/// Derive a variance value that is safe to send into `sqrt`.
///
/// Coverage values are stored as `f32` and later promoted to `f64` for summary statistics.
/// The variance formula used here is `E[x^2] - E[x]^2`, which subtracts two nearby floating-
/// point values and can therefore leave behind a tiny negative residue even when the exact
/// mathematical variance is zero.
///
/// We only repair the narrow case that is impossible in exact arithmetic but expected from
/// floating-point cancellation:
/// - Tiny finite negative values within the shared `f32` tolerance are snapped to exact zero
/// - Larger finite negative values are treated as invariant violations and returned as errors
/// - Positive values, even tiny ones, are kept because a very small positive variance is valid
///
/// The last point is important: `sqrt` has no special problem with a tiny positive variance.
/// For example, `sqrt(1e-9)` is about `3.16e-5`, which is still an ordinary finite result.
/// What we need to prevent is taking `sqrt` of a materially negative value and silently
/// producing `NaN` in the output.
fn derive_nonnegative_variance_coverage(
    eligible_positions: u64,
    coverage_sum_of_squares: f64,
    mean_coverage: f64,
) -> Result<f64> {
    let eligible_positions_f64 = eligible_positions as f64;
    let raw_variance = coverage_sum_of_squares / eligible_positions_f64 - mean_coverage.powi(2);

    if raw_variance.is_finite() && raw_variance < 0.0 {
        // Tiny negative variance is a known cancellation artifact from `E[x^2] - E[x]^2`
        // after `f32`-originating coverage was promoted to `f64`. Keep the tolerance tied
        // to the source precision so we do not silently widen what counts as recoverable noise
        if raw_variance.abs() <= ZEROISH_F32_TOLERANCE as f64 {
            return Ok(0.0);
        }

        bail!(
            "derived a negative variance_coverage {} below the allowed cancellation tolerance {}. eligible_positions={}, coverage_sum_of_squares={}, mean_coverage={}",
            raw_variance,
            ZEROISH_F32_TOLERANCE,
            eligible_positions,
            coverage_sum_of_squares,
            mean_coverage,
        );
    }

    Ok(raw_variance)
}

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

/// Write a headerless summary-stats row for BED or fixed-size outputs.
pub fn write_summary_stats_row<W: Write>(
    w: &mut W,
    chr: &str,
    interval: Interval<u64>,
    stats: SummaryStatsRow,
    decimals: i32,
) -> Result<()> {
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        chr,
        interval.start(),
        interval.end(),
        stats.span_positions,
        stats.blacklisted_positions,
        stats.eligible_positions,
        stats.nonzero_positions,
        CompactNumber {
            v: stats.coverage_sum,
            decimals
        },
        CompactNumber {
            v: stats.coverage_sum_of_squares,
            decimals
        },
        CompactNumber {
            v: stats.mean_coverage,
            decimals
        },
        CompactNumber {
            v: stats.total_coverage,
            decimals
        },
        CompactNumber {
            v: stats.variance_coverage,
            decimals
        },
        CompactNumber {
            v: stats.sd_coverage,
            decimals
        },
        CoverageCoefficientOfVariation {
            value: stats.coefficient_of_variation_coverage,
            decimals,
        },
        CompactNumber {
            v: stats.covered_fraction,
            decimals
        },
    )?;
    Ok(())
}

/// Write a headerless grouped aggregate row.
pub fn write_grouped_value_row<W: Write>(
    w: &mut W,
    group_idx: u64,
    span_positions: u64,
    blacklisted_positions: u64,
    eligible_positions: u64,
    value: f64,
    decimals: i32,
) -> Result<()> {
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}",
        group_idx,
        span_positions,
        blacklisted_positions,
        eligible_positions,
        CompactNumber { v: value, decimals }
    )?;
    Ok(())
}

/// Write a headerless grouped summary-stats row.
pub fn write_grouped_summary_stats_row<W: Write>(
    w: &mut W,
    group_idx: u64,
    stats: SummaryStatsRow,
    decimals: i32,
) -> Result<()> {
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        group_idx,
        stats.span_positions,
        stats.blacklisted_positions,
        stats.eligible_positions,
        stats.nonzero_positions,
        CompactNumber {
            v: stats.coverage_sum,
            decimals
        },
        CompactNumber {
            v: stats.coverage_sum_of_squares,
            decimals
        },
        CompactNumber {
            v: stats.mean_coverage,
            decimals
        },
        CompactNumber {
            v: stats.total_coverage,
            decimals
        },
        CompactNumber {
            v: stats.variance_coverage,
            decimals
        },
        CompactNumber {
            v: stats.sd_coverage,
            decimals
        },
        CoverageCoefficientOfVariation {
            value: stats.coefficient_of_variation_coverage,
            decimals,
        },
        CompactNumber {
            v: stats.covered_fraction,
            decimals
        },
    )?;
    Ok(())
}

fn aggregate_value_header(action: CoverageWindowAction) -> &'static str {
    match action {
        CoverageWindowAction::Average | CoverageWindowAction::AverageOnUniqueBases => {
            "mean_coverage"
        }
        CoverageWindowAction::Total | CoverageWindowAction::TotalOnUniqueBases => "total_coverage",
        CoverageWindowAction::SummaryStats | CoverageWindowAction::SummaryStatsOnUniqueBases => {
            unreachable!("summary-stats uses a dedicated header")
        }
        CoverageWindowAction::OnlyIncludeThesePositionsUnique
        | CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
            unreachable!("positional actions do not use aggregate headers")
        }
    }
}

pub(crate) fn summary_stats_header() -> &'static str {
    "chromosome\tstart\tend\tspan_positions\tblacklisted_positions\teligible_positions\tnonzero_positions\tcoverage_sum\tcoverage_sum_of_squares\tmean_coverage\ttotal_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\tcovered_fraction"
}

fn grouped_summary_stats_header() -> &'static str {
    "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\tnonzero_positions\tcoverage_sum\tcoverage_sum_of_squares\tmean_coverage\ttotal_coverage\tvariance_coverage\tsd_coverage\tcoefficient_of_variation_coverage\tcovered_fraction"
}

/// Write one final BED aggregate output file across all requested chromosomes.
///
/// This is the top-level writer for ordinary `--by-bed` aggregate outputs after tile counting has
/// finished. It opens the final compressed output file, writes the correct header for the selected
/// action, and then streams one chromosome at a time through the matching reducer path.
///
/// Summary-stats uses the raw-row reducer and derives the final statistics here, because the
/// reducer itself should stay responsible only for exact additive raw moments. `average` and
/// `total` use the lighter reducer that can write final numeric values directly.
///
/// Parameters
/// ----------
/// - `final_path`:
///     Final compressed output path to create.
/// - `temp_dir`:
///     Directory that holds the per-tile partial files.
/// - `partials_prefix`:
///     Shared filename prefix used to discover partial rows for this output mode.
/// - `windows_by_chr`:
///     BED windows keyed by chromosome. These windows must match the identities written during
///     tile counting.
/// - `chromosomes`:
///     Chromosome order to reduce and write.
/// - `masked`:
///     Whether blacklist masking was active. Needed for the non-summary reducer path.
/// - `action`:
///     Requested aggregate mode for these BED windows.
/// - `decimals`:
///     Decimal precision used in the final output.
/// - `n_threads`:
///     Compression worker count for the output writer.
///
/// Returns
/// -------
/// - `()`:
///     Completes once every requested chromosome has been reduced and written to `final_path`.
pub(crate) fn write_bed_aggregate_output(
    final_path: &Path,
    temp_dir: &Path,
    partials_prefix: &str,
    windows_by_chr: &FxHashMap<String, Windows>,
    chromosomes: &[String],
    masked: bool,
    action: CoverageWindowAction,
    decimals: i32,
    n_threads: usize,
) -> Result<()> {
    let mut writer = open_zstd_auto_writer(final_path, 3, Some(n_threads as u32))?;
    if action.is_summary_stats() {
        // Summary-stats keeps the reducer focused on exact additive raw rows.
        // We derive the human-facing metrics here, right before writing the final TSV row.
        writeln!(writer, "{}", summary_stats_header())?;
        for chromosome in chromosomes {
            if let Some(windows_for_chr) = windows_by_chr.get(chromosome) {
                reduce_bed_with_cross_index_for_chr_rows(
                    chromosome,
                    temp_dir,
                    partials_prefix,
                    windows_for_chr.as_slice(),
                    |row| {
                        let stats = derive_summary_stats(
                            row.interval.len(),
                            row.blacklisted_positions,
                            row.eligible_positions,
                            row.nonzero_positions,
                            row.coverage_sum,
                            row.coverage_sum_of_squares,
                        )?;
                        write_summary_stats_row(
                            &mut writer,
                            chromosome,
                            row.interval,
                            stats,
                            decimals,
                        )
                    },
                )?;
            }
        }
    } else {
        // `average` and `total` can use the lighter reducer path because they do not need raw
        // moments such as `nonzero_positions` or `coverage_sum_of_squares`.
        writeln!(
            writer,
            "chromosome\tstart\tend\t{}\tblacklisted_positions",
            aggregate_value_header(action)
        )?;
        for chromosome in chromosomes {
            if let Some(windows_for_chr) = windows_by_chr.get(chromosome) {
                reduce_bed_with_cross_index_for_chr(
                    chromosome,
                    temp_dir,
                    partials_prefix,
                    windows_for_chr.as_slice(),
                    masked,
                    action,
                    decimals,
                    &mut writer,
                )?;
            }
        }
    }
    writer.flush()?;
    Ok(())
}

/// Write one final grouped BED aggregate output file across all requested chromosomes.
///
/// Grouped BED outputs are reduced in two stages:
/// 1. reduce each chromosome's segment-level partial rows back into exact per-segment totals
/// 2. fold those segment totals into per-group accumulators using `segment_idx -> group_idx`
///
/// This writer owns the second stage. It preloads one accumulator per group with that group's
/// total span, walks the segment reducers chromosome by chromosome, merges each reduced segment
/// into its group, and finally writes one row per `group_idx` in deterministic order.
///
/// Summary-stats keeps the full raw moments through the group fold and derives the final
/// statistics only once the complete grouped row is known. Non-summary grouped actions keep only
/// the fields needed for final `average` or `total`.
///
/// Parameters
/// ----------
/// - `final_path`:
///     Final compressed output path to create.
/// - `temp_dir`:
///     Directory that holds the per-tile partial files.
/// - `partials_prefix`:
///     Shared filename prefix used to discover partial rows for this output mode.
/// - `grouped_layout`:
///     Precomputed grouped segment layout, including `segment_idx -> group_idx` and stable group
///     metadata.
/// - `chromosomes`:
///     Chromosome order to reduce and write.
/// - `action`:
///     Requested grouped aggregate mode.
/// - `decimals`:
///     Decimal precision used in the final output.
/// - `n_threads`:
///     Compression worker count for the output writer.
///
/// Returns
/// -------
/// - `()`:
///     Completes once every reduced segment has been folded into its group and all grouped rows
///     have been written to `final_path`.
pub(crate) fn write_grouped_bed_aggregate_output(
    final_path: &Path,
    temp_dir: &Path,
    partials_prefix: &str,
    grouped_layout: &GroupedCoverageLayout,
    chromosomes: &[String],
    action: CoverageWindowAction,
    decimals: i32,
    n_threads: usize,
) -> Result<()> {
    // Start every group accumulator with its known span. Segment reduction will fill in the
    // coverage-derived fields later, but span positions come from the grouped layout itself.
    let mut grouped_accums: FxHashMap<u64, GroupedAggregateAccum> = FxHashMap::default();
    for (group_idx, span_positions) in &grouped_layout.group_span_positions {
        grouped_accums.insert(
            *group_idx,
            GroupedAggregateAccum {
                span_positions: *span_positions,
                ..Default::default()
            },
        );
    }

    if action.is_summary_stats() {
        // First reduce segment-level summary rows, then fold those raw moments into the owning
        // group. This keeps the segment reducer and grouped writer responsibilities cleanly split.
        for chromosome in chromosomes {
            if let Some(windows_for_chr) = grouped_layout.segments_by_chr.get(chromosome) {
                reduce_bed_with_cross_index_for_chr_rows(
                    chromosome,
                    temp_dir,
                    partials_prefix,
                    windows_for_chr.as_slice(),
                    |row| {
                        let group_idx = *grouped_layout
                            .segment_idx_to_group_idx
                            .get(&row.idx)
                            .ok_or_else(|| {
                                anyhow::anyhow!(
                                    "grouped reducer is missing a group mapping for segment_idx {}",
                                    row.idx
                                )
                            })?;
                        let group_accum = grouped_accums.get_mut(&group_idx).ok_or_else(|| {
                            anyhow::anyhow!("missing grouped accumulator for group_idx {}", group_idx)
                        })?;
                        // Summary-stats fields are additive across the group's reduced segments.
                        group_accum.blacklisted_positions += row.blacklisted_positions;
                        group_accum.eligible_positions += row.eligible_positions;
                        group_accum.nonzero_positions += row.nonzero_positions;
                        group_accum.coverage_sum += row.coverage_sum;
                        group_accum.coverage_sum_of_squares += row.coverage_sum_of_squares;
                        Ok(())
                    },
                )?;
            }
        }
    } else {
        // The lighter grouped path folds only the fields needed for `average` and `total`.
        for chromosome in chromosomes {
            if let Some(windows_for_chr) = grouped_layout.segments_by_chr.get(chromosome) {
                reduce_bed_basic_with_cross_index_for_chr_rows(
                    chromosome,
                    temp_dir,
                    partials_prefix,
                    windows_for_chr.as_slice(),
                    |row| {
                        let group_idx = *grouped_layout
                            .segment_idx_to_group_idx
                            .get(&row.idx)
                            .ok_or_else(|| {
                                anyhow::anyhow!(
                                    "grouped reducer is missing a group mapping for segment_idx {}",
                                    row.idx
                                )
                            })?;
                        let group_accum = grouped_accums.get_mut(&group_idx).ok_or_else(|| {
                            anyhow::anyhow!("missing grouped accumulator for group_idx {}", group_idx)
                        })?;
                        // Non-summary grouped outputs never need the summary-only raw moments.
                        group_accum.blacklisted_positions += row.blacklisted_positions;
                        group_accum.eligible_positions += row.eligible_positions;
                        group_accum.coverage_sum += row.coverage_sum;
                        Ok(())
                    },
                )?;
            }
        }
    }

    // Always write groups in stable ascending `group_idx` order, regardless of chromosome order
    // or the order in which segment reductions completed.
    let mut sorted_group_indices: Vec<u64> =
        grouped_layout.group_idx_to_name.keys().copied().collect();
    sorted_group_indices.sort_unstable();

    let mut writer = open_zstd_auto_writer(final_path, 3, Some(n_threads as u32))?;
    if action.is_summary_stats() {
        writeln!(writer, "{}", grouped_summary_stats_header())?;
        for group_idx in sorted_group_indices {
            // Missing groups default to all-zero accumulators. That only matters for groups that
            // exist in the layout but received no reduced segment contributions in this run.
            let accum = grouped_accums.get(&group_idx).copied().unwrap_or_default();
            let stats = derive_summary_stats(
                accum.span_positions,
                accum.blacklisted_positions,
                accum.eligible_positions,
                accum.nonzero_positions,
                accum.coverage_sum,
                accum.coverage_sum_of_squares,
            )?;
            write_grouped_summary_stats_row(&mut writer, group_idx, stats, decimals)?;
        }
    } else {
        writeln!(
            writer,
            "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\t{}",
            aggregate_value_header(action)
        )?;
        for group_idx in sorted_group_indices {
            let accum = grouped_accums.get(&group_idx).copied().unwrap_or_default();
            // Group-level `average` and `total` are derived only after all reduced segments have
            // been folded into the group accumulator.
            let value = match action {
                CoverageWindowAction::Average | CoverageWindowAction::AverageOnUniqueBases => {
                    if accum.eligible_positions == 0 {
                        0.0
                    } else {
                        accum.coverage_sum / accum.eligible_positions as f64
                    }
                }
                CoverageWindowAction::Total | CoverageWindowAction::TotalOnUniqueBases => {
                    accum.coverage_sum
                }
                _ => unreachable!("summary-stats handled above"),
            };
            write_grouped_value_row(
                &mut writer,
                group_idx,
                accum.span_positions,
                accum.blacklisted_positions,
                accum.eligible_positions,
                value,
                decimals,
            )?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Write one final fixed-size aggregate output file across all requested chromosomes.
///
/// This is the top-level writer for `--by-size` aggregate outputs. It has two distinct paths:
/// - when tile and bin boundaries are guaranteed to align, it concatenates pre-finalized per-tile
///   files directly for speed
/// - otherwise, it opens the final output and runs the matching cross-tile reducer per chromosome
///
/// Summary-stats again derives the final human-facing metrics only after exact raw rows have been
/// reduced. Non-summary outputs use the lighter reducer path.
///
/// Parameters
/// ----------
/// - `final_path`:
///     Final compressed output path to create.
/// - `temp_dir`:
///     Directory that holds per-tile finals or partial files.
/// - `partials_prefix`:
///     Shared filename prefix used to discover cross-tile partial rows when reduction is needed.
/// - `finals_prefix`:
///     Shared filename prefix used to discover already-finalized aligned tile outputs.
/// - `chromosomes`:
///     Chromosome order to reduce or concatenate.
/// - `contigs`:
///     Contig metadata used to recover chromosome lengths for the reducer path.
/// - `masked`:
///     Whether blacklist masking was active. Needed for the non-summary reducer path.
/// - `action`:
///     Requested aggregate mode for the fixed-size bins.
/// - `decimals`:
///     Decimal precision used in the final output.
/// - `n_threads`:
///     Compression worker count for the output writer.
/// - `tile_and_window_boundaries_align`:
///     Whether every tile final already matches complete fixed-size bins and can therefore be
///     concatenated directly without any cross-tile reduction.
///
/// Returns
/// -------
/// - `()`:
///     Completes once all requested chromosomes have been concatenated or reduced into
///     `final_path`.
pub(crate) fn write_size_aggregate_output(
    final_path: &Path,
    temp_dir: &Path,
    partials_prefix: &str,
    finals_prefix: &str,
    chromosomes: &[String],
    contigs: &crate::shared::bam::Contigs,
    masked: bool,
    action: CoverageWindowAction,
    decimals: i32,
    n_threads: usize,
    tile_and_window_boundaries_align: bool,
) -> Result<()> {
    if tile_and_window_boundaries_align {
        // Fast path: every tile already wrote complete fixed-size bins, so we can concatenate the
        // ready-made outputs instead of reopening the cross-tile reducers.
        let header = if action.is_summary_stats() {
            summary_stats_header().to_string()
        } else {
            format!(
                "chromosome\tstart\tend\t{}\tblacklisted_positions",
                aggregate_value_header(action)
            )
        };
        let _ = concat_aligned_size_tile_finals(
            temp_dir,
            &final_path.parent().map(PathBuf::from).ok_or_else(|| {
                anyhow::anyhow!("missing parent directory for {}", final_path.display())
            })?,
            chromosomes,
            finals_prefix,
            final_path
                .file_name()
                .and_then(|name| name.to_str())
                .ok_or_else(|| anyhow::anyhow!("invalid final output filename"))?,
            header.as_str(),
        )?;
        return Ok(());
    }

    // General path: some bins crossed tile boundaries, so we have to reduce the partial rows back
    // into complete fixed-size bins chromosome by chromosome.
    let mut writer = open_zstd_auto_writer(final_path, 3, Some(n_threads as u32))?;
    if action.is_summary_stats() {
        writeln!(writer, "{}", summary_stats_header())?;
        for chromosome in chromosomes {
            let chrom_len = contigs
                .contigs
                .get(chromosome)
                .map(|&(_, len)| len as u64)
                .ok_or_else(|| {
                    anyhow::anyhow!("Chromosome '{}' not found in contig map", chromosome)
                })?;
            // Summary-stats keeps the reducer focused on exact raw moments, then derives the final
            // metrics here once the full bin row is known.
            reduce_aggregates_by_size_with_cross_index_for_chr_rows(
                chromosome,
                temp_dir,
                partials_prefix,
                chrom_len,
                |row| {
                    let stats = derive_summary_stats(
                        row.interval.len(),
                        row.blacklisted_positions,
                        row.eligible_positions,
                        row.nonzero_positions,
                        row.coverage_sum,
                        row.coverage_sum_of_squares,
                    )?;
                    write_summary_stats_row(&mut writer, chromosome, row.interval, stats, decimals)
                },
            )?;
        }
    } else {
        writeln!(
            writer,
            "chromosome\tstart\tend\t{}\tblacklisted_positions",
            aggregate_value_header(action)
        )?;
        for chromosome in chromosomes {
            let chrom_len = contigs
                .contigs
                .get(chromosome)
                .map(|&(_, len)| len as u64)
                .ok_or_else(|| {
                    anyhow::anyhow!("Chromosome '{}' not found in contig map", chromosome)
                })?;
            // `average` and `total` can use the lighter reducer because they do not need the
            // summary-only raw moments.
            reduce_aggregates_by_size_with_cross_index_for_chr(
                chromosome,
                temp_dir,
                partials_prefix,
                masked,
                action,
                chrom_len,
                decimals,
                &mut writer,
            )?;
        }
    }
    writer.flush()?;
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
pub fn write_bedgraph_runs<W: Write>(
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
/// The helper mirrors `write_bedgraph_runs` but optionally appends the window's original index to
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
pub fn write_windowed_runs<W: Write>(
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

#[cfg(test)]
mod tests {
    include!("writers_tests.rs");
}
