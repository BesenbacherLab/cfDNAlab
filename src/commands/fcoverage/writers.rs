use crate::shared::formatters::{CompactNumber, round_to, round_to_with_precomputed_factor};
use crate::shared::interval::Interval;
use crate::shared::writers::open_zstd_auto_writer;
use anyhow::{Result, bail};
use fxhash::FxHashMap;
use std::fmt;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::commands::fcoverage::reducer::{
    ReducedAggregateRow, TileAggregateTempFiles, reduce_bed_rows, reduce_size_rows,
};
use crate::commands::fcoverage::tiling::{
    TileTempFile, concat_aligned_size_tile_final_outputs, finalize_value,
};
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

/// Fold one reduced segment row into its owning grouped BED accumulator.
///
/// This is intentionally writer-side logic rather than reducer logic. By the time this helper runs,
/// the reducer has already rebuilt exact per-segment raw rows. Grouping is a second stage that:
/// - looks up `segment_idx -> group_idx`
/// - adds the reduced segment's raw fields into the group's totals
/// - keeps `span_positions` owned by the grouped layout rather than by the reducer
fn fold_reduced_segment_into_group(
    grouped_layout: &GroupedCoverageLayout,
    grouped_accums: &mut FxHashMap<u64, GroupedAggregateAccum>,
    row: ReducedAggregateRow,
) -> Result<()> {
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

    // Reduced segment rows already have stable identity and exact additive raw fields.
    // Group folding should therefore stay as simple field-wise addition, regardless of whether
    // the source reducer came from basic or summary temp files.
    group_accum.blacklisted_positions += row.blacklisted_positions;
    group_accum.eligible_positions += row.eligible_positions;
    group_accum.nonzero_positions += row.nonzero_positions;
    group_accum.coverage_sum += row.coverage_sum;
    group_accum.coverage_sum_of_squares += row.coverage_sum_of_squares;
    Ok(())
}

/// Summary statistics for one final output row.
///
/// Reducers and grouped folding keep internal raw fields named for their additive role, such as
/// `coverage_sum`. This final row uses the public TSV names instead, so the writer cannot
/// accidentally expose both an internal raw name and a user-facing alias for the same value.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SummaryStatsRow {
    pub(crate) span_positions: u64,
    pub(crate) blacklisted_positions: u64,
    pub(crate) eligible_positions: u64,
    pub(crate) nonzero_positions: u64,
    pub(crate) covered_fraction: f64,
    pub(crate) total_coverage: f64,
    pub(crate) total_squared_coverage: f64,
    pub(crate) average_coverage: f64,
    pub(crate) variance_coverage: f64,
    pub(crate) sd_coverage: f64,
    pub(crate) coefficient_of_variation_coverage: f64,
}

/// Formatter for the CV column only.
///
/// Most summary-stat fields should stay fully numeric, but CV has a special readability issue:
/// when the mean is extremely small, the exact finite ratio can become so large that it no
/// longer helps users compare windows. This formatter keeps ordinary values numeric and only
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
pub(crate) fn derive_summary_stats(
    span_positions: u64,
    blacklisted_positions: u64,
    eligible_positions: u64,
    nonzero_positions: u64,
    total_coverage: f64,
    total_squared_coverage: f64,
) -> Result<SummaryStatsRow> {
    if eligible_positions == 0 {
        return Ok(SummaryStatsRow {
            span_positions,
            blacklisted_positions,
            eligible_positions,
            nonzero_positions,
            covered_fraction: f64::NAN,
            total_coverage,
            total_squared_coverage,
            average_coverage: f64::NAN,
            variance_coverage: f64::NAN,
            sd_coverage: f64::NAN,
            coefficient_of_variation_coverage: f64::NAN,
        });
    }

    let eligible_positions_f64 = eligible_positions as f64;
    let average_coverage = total_coverage / eligible_positions_f64;
    let variance_coverage = derive_nonnegative_variance_coverage(
        eligible_positions,
        total_squared_coverage,
        average_coverage,
    )?;
    let sd_coverage = variance_coverage.sqrt();
    let coefficient_of_variation_coverage =
        derive_coefficient_of_variation_coverage(average_coverage, sd_coverage);
    let covered_fraction = nonzero_positions as f64 / eligible_positions_f64;

    Ok(SummaryStatsRow {
        span_positions,
        blacklisted_positions,
        eligible_positions,
        nonzero_positions,
        covered_fraction,
        total_coverage,
        total_squared_coverage,
        average_coverage,
        variance_coverage,
        sd_coverage,
        coefficient_of_variation_coverage,
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
fn derive_coefficient_of_variation_coverage(average_coverage: f64, sd_coverage: f64) -> f64 {
    if !average_coverage.is_finite() || average_coverage.abs() <= ZEROISH_COVERAGE_MEAN_F64 {
        f64::NAN
    } else {
        sd_coverage / average_coverage
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
    total_squared_coverage: f64,
    average_coverage: f64,
) -> Result<f64> {
    let eligible_positions_f64 = eligible_positions as f64;
    let raw_variance = total_squared_coverage / eligible_positions_f64 - average_coverage.powi(2);

    if raw_variance.is_finite() && raw_variance < 0.0 {
        // Tiny negative variance is a known cancellation artifact from `E[x^2] - E[x]^2`
        // after `f32`-originating coverage was promoted to `f64`. Keep the tolerance tied
        // to the source precision so we do not silently widen what counts as recoverable noise
        if raw_variance.abs() <= ZEROISH_F32_TOLERANCE as f64 {
            return Ok(0.0);
        }

        bail!(
            "derived a negative variance_coverage {} below the allowed cancellation tolerance {}. eligible_positions={}, total_squared_coverage={}, average_coverage={}",
            raw_variance,
            ZEROISH_F32_TOLERANCE,
            eligible_positions,
            total_squared_coverage,
            average_coverage,
        );
    }

    Ok(raw_variance)
}

/// Write the shared summary-stat value columns that follow the row identity columns.
///
/// BED/fixed-size rows and grouped rows differ only in their leading identity fields. The
/// scientific value columns after that are identical, so this helper keeps the formatting logic in one
/// place without hiding which writer owns which leading columns.
fn write_summary_stats_fields<W: Write>(
    w: &mut W,
    stats: SummaryStatsRow,
    decimals: i32,
) -> Result<()> {
    write!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        stats.span_positions,
        stats.blacklisted_positions,
        stats.eligible_positions,
        stats.nonzero_positions,
        CompactNumber {
            v: stats.covered_fraction,
            decimals
        },
        CompactNumber {
            v: stats.total_coverage,
            decimals
        },
        CompactNumber {
            v: stats.total_squared_coverage,
            decimals
        },
        CompactNumber {
            v: stats.average_coverage,
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
    )?;
    Ok(())
}

/// Derive summary statistics from one reduced raw row and write the final TSV row.
///
/// Reducers intentionally stop at exact additive raw values. This helper keeps the writer-side
/// summary path compact while preserving that separation: derive summary statistics only once the
/// full reduced row is known, then write them in the public output schema.
fn write_reduced_summary_stats_row<W: Write>(
    w: &mut W,
    chr: &str,
    row: ReducedAggregateRow,
    decimals: i32,
) -> Result<()> {
    let stats = derive_summary_stats(
        row.interval.len(),
        row.blacklisted_positions,
        row.eligible_positions,
        row.nonzero_positions,
        row.coverage_sum,
        row.coverage_sum_of_squares,
    )?;
    write_summary_stats_row(w, chr, row.interval, stats, decimals)
}

fn scale_reduced_row(
    row: ReducedAggregateRow,
    restore_mean_multiplier: Option<f64>,
) -> ReducedAggregateRow {
    let Some(multiplier) = restore_mean_multiplier else {
        return row;
    };

    ReducedAggregateRow {
        coverage_sum: row.coverage_sum * multiplier,
        coverage_sum_of_squares: row.coverage_sum_of_squares * multiplier.powi(2),
        ..row
    }
}

fn write_non_summary_reduced_row<W: Write>(
    writer: &mut W,
    chromosome: &str,
    row: ReducedAggregateRow,
    masked: bool,
    action: CoverageWindowAction,
    decimals: i32,
    restore_mean_multiplier: Option<f64>,
) -> Result<()> {
    let row = scale_reduced_row(row, restore_mean_multiplier);
    let value = finalize_value(
        row.coverage_sum,
        row.eligible_positions,
        row.interval.len(),
        masked,
        &action,
    );
    write_final_row(
        writer,
        chromosome,
        row.interval,
        round_to(value, decimals),
        row.blacklisted_positions,
        decimals,
    )
}

/// Write a final aggregate row: `chromosome  start  end  value  blacklisted_positions`
#[inline]
pub(crate) fn write_final_row<W: Write>(
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
pub(crate) fn write_summary_stats_row<W: Write>(
    w: &mut W,
    chr: &str,
    interval: Interval<u64>,
    stats: SummaryStatsRow,
    decimals: i32,
) -> Result<()> {
    write!(w, "{}\t{}\t{}\t", chr, interval.start(), interval.end(),)?;
    write_summary_stats_fields(w, stats, decimals)?;
    writeln!(w)?;
    Ok(())
}

/// Write a headerless grouped aggregate row.
pub(crate) fn write_grouped_value_row<W: Write>(
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
pub(crate) fn write_grouped_summary_stats_row<W: Write>(
    w: &mut W,
    group_idx: u64,
    stats: SummaryStatsRow,
    decimals: i32,
) -> Result<()> {
    write!(w, "{}\t", group_idx)?;
    write_summary_stats_fields(w, stats, decimals)?;
    writeln!(w)?;
    Ok(())
}

fn aggregate_value_header(action: CoverageWindowAction, signal_label: &str) -> String {
    match action {
        CoverageWindowAction::Average | CoverageWindowAction::AverageOnUniqueBases => {
            format!("average_{signal_label}")
        }
        CoverageWindowAction::Total | CoverageWindowAction::TotalOnUniqueBases => {
            format!("total_{signal_label}")
        }
        CoverageWindowAction::SummaryStats | CoverageWindowAction::SummaryStatsOnUniqueBases => {
            unreachable!("summary-stats uses a dedicated header")
        }
        CoverageWindowAction::OnlyIncludeThesePositionsUnique
        | CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
            unreachable!("positional actions do not use aggregate headers")
        }
    }
}

pub(crate) fn summary_stats_header(signal_label: &str) -> String {
    format!(
        "chromosome\tstart\tend\tspan_positions\tblacklisted_positions\teligible_positions\tnonzero_positions\tcovered_fraction\ttotal_{signal_label}\ttotal_squared_{signal_label}\taverage_{signal_label}\tvariance_{signal_label}\tsd_{signal_label}\tcoefficient_of_variation_{signal_label}"
    )
}

fn grouped_summary_stats_header(signal_label: &str) -> String {
    format!(
        "group_idx\tspan_positions\tblacklisted_positions\teligible_positions\tnonzero_positions\tcovered_fraction\ttotal_{signal_label}\ttotal_squared_{signal_label}\taverage_{signal_label}\tvariance_{signal_label}\tsd_{signal_label}\tcoefficient_of_variation_{signal_label}"
    )
}

/// Convert one grouped raw accumulator into the final non-summary reported value.
///
/// Group folding keeps exact additive raw totals until the very end. This helper makes the final
/// `average` vs `total` policy explicit in one place and keeps the grouped write loop focused on
/// output order and formatting.
fn finalize_grouped_value(action: CoverageWindowAction, accum: GroupedAggregateAccum) -> f64 {
    match action {
        CoverageWindowAction::Average | CoverageWindowAction::AverageOnUniqueBases => {
            if accum.eligible_positions == 0 {
                f64::NAN
            } else {
                accum.coverage_sum / accum.eligible_positions as f64
            }
        }
        CoverageWindowAction::Total | CoverageWindowAction::TotalOnUniqueBases => {
            accum.coverage_sum
        }
        _ => unreachable!("summary-stats handled separately"),
    }
}

/// Look up one chromosome length from the contig map and convert it to `u64`.
///
/// Fixed-size reduction needs the true chromosome end so the final bin can be clipped after
/// cross-tile reduction. Keeping the lookup here avoids repeating the same error handling in each
/// writer branch.
fn chromosome_length_u64(contigs: &crate::shared::bam::Contigs, chromosome: &str) -> Result<u64> {
    contigs
        .contigs
        .get(chromosome)
        .map(|&(_, len)| len as u64)
        .ok_or_else(|| anyhow::anyhow!("Chromosome '{}' not found in contig map", chromosome))
}

fn aggregate_tile_outputs_for_chromosome<'a>(
    tile_outputs_by_chr: &'a FxHashMap<String, Vec<TileAggregateTempFiles>>,
    chromosome: &str,
    output_kind: &str,
) -> Result<&'a [TileAggregateTempFiles]> {
    let tile_outputs = tile_outputs_by_chr.get(chromosome).ok_or_else(|| {
        anyhow::anyhow!(
            "No returned {} aggregate tile outputs for chromosome '{}'",
            output_kind,
            chromosome
        )
    })?;
    anyhow::ensure!(
        !tile_outputs.is_empty(),
        "Returned {} aggregate tile output list is empty for chromosome '{}'",
        output_kind,
        chromosome
    );
    Ok(tile_outputs.as_slice())
}

/// Write one final BED aggregate output across all requested chromosomes.
///
/// This is the top-level writer for ordinary `--by-bed` aggregate outputs after tile counting has
/// returned explicit temp paths. It opens the final compressed output file, writes the selected
/// header, and reduces each chromosome through the BED reducer using only the returned tile paths.
///
/// Summary-stats keeps the reducer focused on exact additive raw rows. The writer derives mean,
/// nonzero, and dispersion columns immediately before writing each final TSV row. `average` and
/// `total` use the same raw reducer but finalize only the requested scalar value.
///
/// Parameters
/// ----------
/// - `final_path`:
///   Final compressed output path to create.
/// - `tile_outputs_by_chr`:
///   Returned aggregate tile paths grouped by chromosome.
/// - `windows_by_chr`:
///   BED windows keyed by chromosome. These identities must match the partial rows.
/// - `chromosomes`:
///   Chromosome order to reduce and write.
/// - `masked`:
///   Whether blacklist masking was active. Needed for scalar finalization.
/// - `action`:
///   Requested aggregate mode for these BED windows.
/// - `decimals`:
///   Decimal precision used in the final output.
/// - `n_threads`:
///   Compression worker count for the output writer.
/// - `signal_label`:
///   Public signal name used inside aggregate value column names.
/// - `restore_mean_multiplier`:
///   Optional late multiplier applied to raw coverage sums before finalization.
pub(crate) fn write_bed_aggregate_output(
    final_path: &Path,
    tile_outputs_by_chr: &FxHashMap<String, Vec<TileAggregateTempFiles>>,
    windows_by_chr: &FxHashMap<String, Windows>,
    chromosomes: &[String],
    masked: bool,
    action: CoverageWindowAction,
    decimals: i32,
    n_threads: usize,
    signal_label: &str,
    restore_mean_multiplier: Option<f64>,
) -> Result<()> {
    let mut writer = open_zstd_auto_writer(final_path, 3, Some(n_threads as u32))?;
    if action.is_summary_stats() {
        // Summary-stats derives final columns from exact reduced sums and moments.
        writeln!(writer, "{}", summary_stats_header(signal_label))?;
        for chromosome in chromosomes {
            let Some(windows_for_chr) = windows_by_chr.get(chromosome) else {
                continue;
            };
            if windows_for_chr.is_empty() {
                continue;
            }
            let tile_outputs =
                aggregate_tile_outputs_for_chromosome(tile_outputs_by_chr, chromosome, "BED")?;
            reduce_bed_rows(
                chromosome,
                tile_outputs,
                windows_for_chr.as_slice(),
                true,
                |row| {
                    write_reduced_summary_stats_row(
                        &mut writer,
                        chromosome,
                        scale_reduced_row(row, restore_mean_multiplier),
                        decimals,
                    )
                },
            )?;
        }
    } else {
        // Average and total use the same raw reducer, then write only the requested scalar.
        writeln!(
            writer,
            "chromosome\tstart\tend\t{}\tblacklisted_positions",
            aggregate_value_header(action, signal_label)
        )?;
        for chromosome in chromosomes {
            let Some(windows_for_chr) = windows_by_chr.get(chromosome) else {
                continue;
            };
            if windows_for_chr.is_empty() {
                continue;
            }
            let tile_outputs =
                aggregate_tile_outputs_for_chromosome(tile_outputs_by_chr, chromosome, "BED")?;
            reduce_bed_rows(
                chromosome,
                tile_outputs,
                windows_for_chr.as_slice(),
                false,
                |row| {
                    write_non_summary_reduced_row(
                        &mut writer,
                        chromosome,
                        row,
                        masked,
                        action,
                        decimals,
                        restore_mean_multiplier,
                    )
                },
            )?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Write one final grouped BED aggregate output across all requested chromosomes.
///
/// Grouped BED output is a second pass over ordinary BED reduction. The reducer first reconstructs
/// exact per-segment rows from the returned tile paths. Those rows are then added to their owning
/// groups using the grouped layout's `segment_idx -> group_idx` map.
///
/// Group span is not recomputed from observed segment rows. It is initialized from
/// `group_span_positions` so groups with no reduced coverage rows still keep their represented
/// span and write zero coverage-derived fields. Final rows are written in ascending `group_idx`.
///
/// Basic and summary reducers share the same fold because basic partial rows are expanded into the
/// same in-memory accumulator shape with summary-only fields set to zero. Summary-stats derives
/// final statistics after all segments in a group have been folded.
///
/// `signal_label` only changes final public column names. The grouped reducer still folds the same
/// internal additive values.
pub(crate) fn write_grouped_bed_aggregate_output(
    final_path: &Path,
    tile_outputs_by_chr: &FxHashMap<String, Vec<TileAggregateTempFiles>>,
    grouped_layout: &GroupedCoverageLayout,
    chromosomes: &[String],
    action: CoverageWindowAction,
    decimals: i32,
    n_threads: usize,
    signal_label: &str,
    restore_mean_multiplier: Option<f64>,
) -> Result<()> {
    // Start each group with the span declared by the grouped layout. Segment rows fill in only the
    // coverage-derived fields.
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

    // Reduce chromosome-level segment rows back to exact raw values, then fold those reduced
    // segments into their owning groups. Basic and summary outputs use the same accumulation loop
    // because basic reducers fill summary-only fields with zeroes in memory.
    for chromosome in chromosomes {
        let Some(windows_for_chr) = grouped_layout.segments_by_chr.get(chromosome) else {
            continue;
        };
        if windows_for_chr.is_empty() {
            continue;
        }
        // A chromosome with grouped segments must have returned tile outputs. Missing outputs
        // would turn an upstream mismatch into zero coverage.
        let tile_outputs =
            aggregate_tile_outputs_for_chromosome(tile_outputs_by_chr, chromosome, "grouped BED")?;

        if action.is_summary_stats() {
            reduce_bed_rows(
                chromosome,
                tile_outputs,
                windows_for_chr.as_slice(),
                true,
                |row| {
                    fold_reduced_segment_into_group(
                        grouped_layout,
                        &mut grouped_accums,
                        scale_reduced_row(row, restore_mean_multiplier),
                    )
                },
            )?;
        } else {
            reduce_bed_rows(
                chromosome,
                tile_outputs,
                windows_for_chr.as_slice(),
                false,
                |row| {
                    fold_reduced_segment_into_group(
                        grouped_layout,
                        &mut grouped_accums,
                        scale_reduced_row(row, restore_mean_multiplier),
                    )
                },
            )?;
        }
    }

    // Output order is the stable group index order from the grouped layout, not chromosome order
    // or hash-map iteration order.
    let mut sorted_group_indices: Vec<u64> =
        grouped_layout.group_idx_to_name.keys().copied().collect();
    sorted_group_indices.sort_unstable();

    let mut writer = open_zstd_auto_writer(final_path, 3, Some(n_threads as u32))?;
    if action.is_summary_stats() {
        writeln!(writer, "{}", grouped_summary_stats_header(signal_label))?;
        for group_idx in sorted_group_indices {
            // Groups without segment contributions still write one row. Their span was seeded
            // from the layout above, while coverage-derived fields remain zero.
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
            aggregate_value_header(action, signal_label)
        )?;
        for group_idx in sorted_group_indices {
            let accum = grouped_accums.get(&group_idx).copied().unwrap_or_default();
            // Finalize average or total only after every reduced segment has been folded into the
            // group accumulator.
            let value = finalize_grouped_value(action, accum);
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

/// Write one final fixed-size aggregate output across all requested chromosomes.
///
/// This writer has two paths. When tile and bin boundaries are guaranteed to align and no
/// restore-mean multiplier is needed, every tile already wrote complete final rows, so the writer
/// concatenates those returned final paths directly. Otherwise it opens the final output and runs
/// the fixed-size reducer per chromosome using the returned partial paths.
///
/// Restore-mean intentionally bypasses the aligned fast path. In that mode aligned tiles write
/// exact raw partial rows, and the reducer path applies the late multiplier with the same final
/// semantics as non-aligned runs. Missing cross-index files still mean one contribution per bin.
///
/// Summary-stats derives mean, nonzero, and dispersion columns only after exact raw rows have been
/// reduced. Non-summary outputs finalize `average` or `total` from the reduced raw sums.
///
/// `signal_label` only changes final public column names. The fixed-size reducers still use the
/// same internal additive values for coverage and length-normalized fragment mass.
pub(crate) fn write_size_aggregate_output(
    final_path: &Path,
    tile_outputs_by_chr: &FxHashMap<String, Vec<TileAggregateTempFiles>>,
    final_tile_outputs: &[TileTempFile],
    chromosomes: &[String],
    contigs: &crate::shared::bam::Contigs,
    masked: bool,
    action: CoverageWindowAction,
    decimals: i32,
    n_threads: usize,
    tile_and_window_boundaries_align: bool,
    signal_label: &str,
    restore_mean_multiplier: Option<f64>,
) -> Result<()> {
    if tile_and_window_boundaries_align && restore_mean_multiplier.is_none() {
        // Fast path: every tile already wrote complete fixed-size bins, so we can concatenate the
        // ready-made outputs instead of reopening the cross-tile reducers.
        //
        // `restore-mean` intentionally bypasses this fast path and falls through to the shared
        // reducer/finalizer below. In that mode, aligned tiles already wrote exact raw rows to
        // `partials_out`, and the ordinary reducer path can finalize them directly with the
        // existing "missing cross-index means one contribution" rule.

        let header = if action.is_summary_stats() {
            summary_stats_header(signal_label)
        } else {
            format!(
                "chromosome\tstart\tend\t{}\tblacklisted_positions",
                aggregate_value_header(action, signal_label)
            )
        };
        let _ = concat_aligned_size_tile_final_outputs(
            &final_path.parent().map(PathBuf::from).ok_or_else(|| {
                anyhow::anyhow!("missing parent directory for {}", final_path.display())
            })?,
            chromosomes,
            final_tile_outputs,
            final_path
                .file_name()
                .and_then(|name| name.to_str())
                .ok_or_else(|| anyhow::anyhow!("invalid final output filename"))?,
            header.as_str(),
        )?;
        return Ok(());
    }

    // General path: reduce raw partial rows back into complete fixed-size bins chromosome by
    // chromosome. This is required when bins cross tile boundaries.
    //
    // This same path also handles aligned `restore-mean` runs. In that case there is still no
    // cross-tile reduction to do, but reusing the reducer/finalizer avoids a second raw-row
    // writer path with the same final semantics.
    let mut writer = open_zstd_auto_writer(final_path, 3, Some(n_threads as u32))?;
    if action.is_summary_stats() {
        // Summary-stats derives final columns from exact reduced sums and moments.
        writeln!(writer, "{}", summary_stats_header(signal_label))?;
        for chromosome in chromosomes {
            let chrom_len = chromosome_length_u64(contigs, chromosome)?;
            if chrom_len == 0 {
                continue;
            }
            let tile_outputs =
                aggregate_tile_outputs_for_chromosome(tile_outputs_by_chr, chromosome, "size")?;
            reduce_size_rows(chromosome, tile_outputs, chrom_len, true, |row| {
                write_reduced_summary_stats_row(
                    &mut writer,
                    chromosome,
                    scale_reduced_row(row, restore_mean_multiplier),
                    decimals,
                )
            })?;
        }
    } else {
        // Average and total use the same raw reducer, then write only the requested scalar.
        writeln!(
            writer,
            "chromosome\tstart\tend\t{}\tblacklisted_positions",
            aggregate_value_header(action, signal_label)
        )?;
        for chromosome in chromosomes {
            let chrom_len = chromosome_length_u64(contigs, chromosome)?;
            if chrom_len == 0 {
                continue;
            }
            let tile_outputs =
                aggregate_tile_outputs_for_chromosome(tile_outputs_by_chr, chromosome, "size")?;
            reduce_size_rows(chromosome, tile_outputs, chrom_len, false, |row| {
                write_non_summary_reduced_row(
                    &mut writer,
                    chromosome,
                    row,
                    masked,
                    action,
                    decimals,
                    restore_mean_multiplier,
                )
            })?;
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
pub(crate) fn write_bedgraph_runs<W: Write>(
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
            writeln!(
                out,
                "{}\t{}\t{}\t{}",
                chr,
                run_start_abs,
                run_end_abs,
                CompactNumber { v: value, decimals },
            )?;
            Ok(())
        },
    )
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
pub(crate) fn write_windowed_runs<W: Write>(
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
            if let Some(idx) = orig_idx {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    chr,
                    run_start_abs,
                    run_end_abs,
                    CompactNumber { v: value, decimals },
                    idx
                )?;
            } else {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}",
                    chr,
                    run_start_abs,
                    run_end_abs,
                    CompactNumber { v: value, decimals },
                )?;
            };
            Ok(())
        },
    )
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
///   Any error stops iteration and is returned to the caller.
#[inline]
fn visit_runs_in_window(
    cov: &[f32],
    mask: Option<&[u8]>,
    local_start_idx: usize,
    local_end_idx: usize,
    decimals: i32,
    keep_zero_runs: bool,
    on_run: impl FnMut(usize, usize, f64) -> Result<()>,
) -> Result<()> {
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
///   Any error stops iteration and is returned to the caller.
#[inline]
fn visit_runs_unmasked(
    cov: &[f32],
    local_start_idx: usize,
    local_end_idx: usize,
    decimals: i32,
    keep_zero_runs: bool,
    mut on_run: impl FnMut(usize, usize, f64) -> Result<()>,
) -> Result<()> {
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

        on_run(run_start_idx, j, value0)?;
        i = j;
    }
    Ok(())
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
///   Any error stops iteration and is returned to the caller.
#[inline]
fn visit_runs_masked(
    cov: &[f32],
    m: &[u8],
    local_start_idx: usize,
    local_end_idx: usize,
    decimals: i32,
    keep_zero_runs: bool,
    mut on_run: impl FnMut(usize, usize, f64) -> Result<()>,
) -> Result<()> {
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

        on_run(run_start_idx, j, value0)?;
        i = j;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("writers_tests.rs");
}
