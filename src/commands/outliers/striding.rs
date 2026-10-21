use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use rayon::prelude::*;

use crate::{
    commands::outliers::{
        config::{ContextWeighting, OutlierModel},
        distribution::{FitDiagnostics, TwoStageCoverageModel, diagnose_fit, fit_two_stage},
        failure::AnalysisFailures,
        model::CoverageCounts,
    },
    shared::{interval::Interval, triangular::TriangularWeightCache},
};

const MAXIMUM_EXACT_RAW_COVERAGE: f32 = 16_777_215.0;
/// A blended fitting frequency has one million units per eligible base before rounding.
const BLENDED_FREQUENCY_SCALE: f64 = 1_000_000.0;

/// Raw eligible-position coverage histogram for a fixed-size model core.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CoreHistogram {
    /// Fixed-size model core covered by this histogram.
    pub(crate) interval: Interval<u32>,
    /// Eligible-position counts keyed by raw positional coverage.
    pub(crate) counts: CoverageCounts,
    /// Accepted fragments assigned once to this core by a midpoint outside the blacklists.
    pub(crate) fragment_support: u64,
}

impl CoreHistogram {
    pub(crate) fn new(interval: Interval<u32>) -> Self {
        Self {
            interval,
            counts: CoverageCounts::new(),
            fragment_support: 0,
        }
    }

    pub(crate) fn add_coverage(
        &mut self,
        coverage: &[f32],
        blacklist_mask: Option<&[u8]>,
    ) -> Result<()> {
        if let Some(mask) = blacklist_mask {
            ensure!(
                mask.len() == coverage.len(),
                "coverage length {} does not match blacklist-mask length {}",
                coverage.len(),
                mask.len()
            );
        }

        for (position, &coverage_value) in coverage.iter().enumerate() {
            if blacklist_mask.is_some_and(|mask| mask[position] != 0) {
                continue;
            }
            ensure!(
                coverage_value.is_finite() && coverage_value >= 0.0,
                "raw coverage must be finite and nonnegative, got {}",
                coverage_value
            );
            let rounded_coverage = coverage_value.round();
            ensure!(
                (coverage_value - rounded_coverage).abs() <= 1e-4,
                "raw outlier coverage must be integer-valued, got {}",
                coverage_value
            );
            ensure!(
                rounded_coverage <= MAXIMUM_EXACT_RAW_COVERAGE,
                "raw coverage {} reaches the f32 exact-integer limit of {}",
                rounded_coverage,
                MAXIMUM_EXACT_RAW_COVERAGE as u32 + 1
            );
            let coverage = rounded_coverage as u32;
            let positions = self.counts.entry(coverage).or_default();
            *positions = positions
                .checked_add(1)
                .context("coverage histogram count overflow")?;
        }
        Ok(())
    }

    pub(crate) fn add_histogram(&mut self, other: &Self) -> Result<()> {
        for (&coverage, &positions_to_add) in &other.counts {
            let positions = self.counts.entry(coverage).or_default();
            *positions = positions
                .checked_add(positions_to_add)
                .context("coverage histogram count overflow")?;
        }
        self.fragment_support = self
            .fragment_support
            .checked_add(other.fragment_support)
            .context("fragment support overflow")?;
        Ok(())
    }

    pub(crate) fn eligible_positions(&self) -> Result<u64> {
        self.counts.values().try_fold(0_u64, |total, &positions| {
            total
                .checked_add(positions)
                .context("eligible-position count overflow")
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ModelSource {
    Local,
    GlobalFallback,
}

impl ModelSource {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Local => "local",
            Self::GlobalFallback => "global_fallback",
        }
    }
}

/// Model and fit diagnostics assigned to a fixed-size model core.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CoreModel {
    /// Fixed-size chromosome interval over which this model supplies T2 and correction targets.
    pub(crate) core: Interval<u32>,
    /// Broader chromosome-local interval assembled to fit the model or assess fallback support.
    pub(crate) context: Interval<u32>,
    /// Eligible positions in `context`, even when the applied model was fitted globally.
    pub(crate) eligible_positions: u64,
    /// Accepted fragments assigned to `context` by a midpoint outside the blacklists.
    pub(crate) fragment_support: u64,
    /// Whether parameters were fitted from `context` or reused from the global histogram.
    pub(crate) source: ModelSource,
    pub(crate) model: TwoStageCoverageModel,
    /// Comparison of this context's original histogram with the applied model.
    pub(crate) diagnostics: FitDiagnostics,
}

/// Local models and the pooled fit, present only when needed for chromosome fallback.
#[derive(Debug)]
pub(crate) struct FittedModels {
    pub(crate) global: Option<TwoStageCoverageModel>,
    /// Original core-index slots. Fully blacklisted cores have no fitted model.
    pub(crate) by_chromosome: FxHashMap<String, Vec<Option<CoreModel>>>,
}

/// Sum all model-core histograms without another BAM scan.
pub(crate) fn sum_global_histogram(
    histograms_by_chromosome: &FxHashMap<String, Vec<CoreHistogram>>,
) -> Result<CoreHistogram> {
    let mut global = CoreHistogram::new(Interval::new(0_u32, 1_u32)?);
    for chromosome_histograms in histograms_by_chromosome.values() {
        for histogram in chromosome_histograms {
            global.add_histogram(histogram)?;
        }
    }
    ensure!(
        global.eligible_positions()? > 0,
        "no eligible positions remained after chromosome selection and blacklisting"
    );
    Ok(global)
}

/// Fit cores with eligible bases from adaptively expanded, chromosome-local histogram contexts.
///
/// Fully blacklisted cores retain an empty slot and require neither local nor global fitting.
/// The global model is used only when the complete chromosome cannot meet both the required
/// eligible-position count and fragment support. Eligible zeros count toward context size.
/// A complete local context whose two-stage fit fails returns that error rather than silently
/// changing its model source. Context weighting changes fitting frequencies in both stages,
/// while support counts and comparisons with observed coverage remain unweighted.
pub(crate) fn fit_core_models(
    chromosomes: &[String],
    histograms_by_chromosome: &FxHashMap<String, Vec<CoreHistogram>>,
    global_histogram: &CoreHistogram,
    family: OutlierModel,
    tail_probability: f64,
    minimum_eligible_positions: u32,
    minimum_fragment_support: u64,
    context_weighting: ContextWeighting,
) -> Result<FittedModels> {
    let triangular_fraction = context_weighting.triangular_fraction()?;
    let mut needs_global = false;
    let mut eligible_positions_by_chromosome = FxHashMap::default();
    for chromosome in chromosomes {
        let cores = histograms_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing core histograms for {chromosome}"))?;
        // Histograms stay immutable during fitting. Count eligibility once per core and reuse
        // it for fallback, empty-slot checks, context growth, and recorded context support
        let eligible_positions = cores
            .iter()
            .map(CoreHistogram::eligible_positions)
            .collect::<Result<Vec<_>>>()?;
        needs_global |= needs_global_fallback(
            cores,
            &eligible_positions,
            minimum_eligible_positions,
            minimum_fragment_support,
        )?;
        eligible_positions_by_chromosome.insert(chromosome.clone(), eligible_positions);
    }
    let mut messages = Vec::new();
    let mut failed_histograms = String::from(
        "# GLOBAL pools all selected chromosomes and has no genomic coordinates (NA)\n# coverage is raw positional coverage, observed_positions counts unweighted eligible bases\n# fitting_frequency includes local distance weighting and blend discretization, whose common normalization cancels in the fit\nchromosome\tcore_start\tcore_end\tcontext_start\tcontext_end\tfragments_by_midpoint\tcoverage\tobserved_positions\tfitting_frequency\n",
    );
    let global_model = if needs_global {
        match fit_two_stage(&global_histogram.counts, tail_probability, family) {
            Ok(model) => Some(model),
            Err(error) => {
                messages.push(format!("required global {} fit, alpha={tail_probability}, eligible positions {}, fragment support {}: {error:#}", family.as_str(), global_histogram.eligible_positions()?, global_histogram.fragment_support));
                append_failed_histogram(
                    &mut failed_histograms,
                    "GLOBAL",
                    global_histogram.interval,
                    global_histogram.interval,
                    global_histogram,
                    None,
                );
                None
            }
        }
    } else {
        None
    };
    let mut models_by_chromosome = FxHashMap::default();
    let triangular_weight_cache = TriangularWeightCache::default();
    for chromosome in chromosomes {
        let cores = &histograms_by_chromosome[chromosome];
        let eligible_positions_by_core = &eligible_positions_by_chromosome[chromosome];
        // Collect every independent local result before returning failures. Histogram processing
        // is cheap compared with another BAM pass and supplies all failing contexts at once
        let results = (0..cores.len())
            .into_par_iter()
            .map(|core_index| {
                let core = cores[core_index].interval;
                if eligible_positions_by_core[core_index] == 0 {
                    // Keep this grid slot without constructing or fitting an unused context
                    return Ok((core, core, 0, CoreHistogram::new(core), None, None));
                }
                let (mut histogram, context, eligible_positions, complete) = build_context(
                    cores,
                    eligible_positions_by_core,
                    core_index,
                    minimum_eligible_positions,
                    minimum_fragment_support,
                )?;
                let mut weighted_counts = if complete && triangular_fraction > 0.0 {
                    let triangular_counts = triangular_context_counts(
                        cores,
                        core_index,
                        context,
                        &triangular_weight_cache,
                    )
                    .with_context(|| {
                        format!(
                            "core {:?}, context {:?}: constructing triangular fitting frequencies",
                            core.as_tuple(),
                            context.as_tuple()
                        )
                    })?;
                    Some(blend_context_counts(
                        &histogram.counts,
                        triangular_counts,
                        eligible_positions,
                        triangular_fraction,
                    )?)
                } else {
                    None
                };
                let fit = if complete {
                    Some(fit_two_stage(
                        weighted_counts.as_ref().unwrap_or(&histogram.counts),
                        tail_probability,
                        family,
                    ))
                } else {
                    global_model.map(Ok)
                };
                let result = fit.map(|fit| -> Result<CoreModel> {
                    let mut model = fit?;
                    if complete {
                        // Integer kernel multipliers are fitting frequencies, not extra bases
                        // Report the actual positions below this weighted fit's T1
                        model.second_fit_retained_positions = histogram
                            .counts
                            .range(..model.initial_threshold)
                            .try_fold(0_u64, |total, (_, &positions)| {
                                total
                                    .checked_add(positions)
                                    .context("retained-position count overflow")
                            })?;
                    }
                    Ok(CoreModel {
                        core,
                        context,
                        eligible_positions,
                        fragment_support: histogram.fragment_support,
                        source: if complete {
                            ModelSource::Local
                        } else {
                            ModelSource::GlobalFallback
                        },
                        diagnostics: diagnose_fit(&histogram.counts, model)?,
                        model,
                    })
                });
                // Retain coverage counts only for failed contexts while collecting all results
                // Successful models already carry their support and numerical diagnostics
                if !matches!(&result, Some(Err(_))) {
                    histogram.counts.clear();
                    weighted_counts = None;
                }
                Ok((
                    core,
                    context,
                    eligible_positions,
                    histogram,
                    weighted_counts,
                    result,
                ))
            })
            .collect::<Vec<Result<_>>>();
        let mut chromosome_models = Vec::with_capacity(cores.len());
        for result in results {
            match result {
                Ok((
                    core,
                    context,
                    eligible_positions,
                    histogram,
                    weighted_counts,
                    Some(Err(error)),
                )) => {
                    messages.push(format!("{chromosome} core {:?}, context {:?}, {} fit, alpha={tail_probability}, eligible positions {}, fragment support {}: {error:#}", core.as_tuple(), context.as_tuple(), family.as_str(), eligible_positions, histogram.fragment_support));
                    append_failed_histogram(
                        &mut failed_histograms,
                        chromosome,
                        core,
                        context,
                        &histogram,
                        weighted_counts.as_ref(),
                    );
                    chromosome_models.push(None);
                }
                Ok((_, _, _, _, _, Some(Ok(model)))) => chromosome_models.push(Some(model)),
                // A fully blacklisted core or an already recorded required global failure
                Ok((_, _, _, _, _, None)) => chromosome_models.push(None),
                Err(error) => {
                    messages.push(format!("{chromosome} context construction: {error:#}"));
                    chromosome_models.push(None);
                }
            }
        }
        models_by_chromosome.insert(chromosome.clone(), chromosome_models);
    }
    if !messages.is_empty() {
        return Err(AnalysisFailures {
            stage: "coverage model fitting after the first BAM pass",
            messages,
            histogram: Some(failed_histograms),
            audit: None,
        }
        .into());
    }
    Ok(FittedModels {
        global: global_model,
        by_chromosome: models_by_chromosome,
    })
}

/// A full chromosome supplies the maximum possible local eligible-position and midpoint support.
///
/// A wholly blacklisted chromosome has no positions needing a model, so it cannot require fallback.
fn needs_global_fallback(
    cores: &[CoreHistogram],
    eligible_positions_by_core: &[u64],
    minimum_eligible_positions: u32,
    minimum_fragment_support: u64,
) -> Result<bool> {
    ensure!(!cores.is_empty(), "chromosome has no model-core histograms");
    let eligible_positions =
        eligible_positions_by_core
            .iter()
            .try_fold(0_u64, |total, &positions| {
                total
                    .checked_add(positions)
                    .context("chromosome eligible-position count overflow")
            })?;
    if eligible_positions == 0 {
        return Ok(false);
    }
    let support = cores.iter().try_fold(0_u64, |total, core| {
        total
            .checked_add(core.fragment_support)
            .context("chromosome fragment support overflow")
    })?;
    Ok(eligible_positions < u64::from(minimum_eligible_positions)
        || support < minimum_fragment_support)
}

/// Preserve sparse failing-context counts without requiring decompression or a whole-run scan.
fn append_failed_histogram(
    output: &mut String,
    chromosome: &str,
    core: Interval<u32>,
    context: Interval<u32>,
    histogram: &CoreHistogram,
    fitting_counts: Option<&CoverageCounts>,
) {
    use std::fmt::Write;
    for (coverage, positions) in &histogram.counts {
        let coordinates = if chromosome == "GLOBAL" {
            "NA\tNA\tNA\tNA".to_string()
        } else {
            format!(
                "{}\t{}\t{}\t{}",
                core.start(),
                core.end(),
                context.start(),
                context.end()
            )
        };
        writeln!(
            output,
            "{chromosome}\t{coordinates}\t{}\t{coverage}\t{positions}\t{}",
            histogram.fragment_support,
            fitting_counts.unwrap_or(&histogram.counts)[coverage]
        )
        .expect("writing to a string cannot fail");
    }
}

/// Pool core histograms using a triangle centered on the core whose model is being estimated.
///
/// If the farthest context core is `radius` grid steps away, its multiplier is 1 and the
/// central multiplier is `radius + 1`. At chromosome ends retain that same distance rule
/// rather than moving the peak to the context midpoint. Masked cores retain their grid
/// distance but contribute no observations. Partial cores contribute only their actual bases.
///
/// Integer multipliers preserve every observed coverage level without rounding frequencies.
/// Dividing all frequencies by their sum gives the normalized triangle. That common factor
/// cancels from ZIP estimates and from the average likelihood optimized for ZINB.
fn triangular_context_counts(
    cores: &[CoreHistogram],
    core_index: usize,
    context: Interval<u32>,
    weight_cache: &TriangularWeightCache,
) -> Result<CoverageCounts> {
    let first_index = cores.partition_point(|core| core.interval.end() <= context.start());
    let end_index = cores.partition_point(|core| core.interval.start() < context.end());
    ensure!(
        first_index <= core_index && core_index < end_index,
        "triangular context does not contain its center core"
    );
    let radius = (core_index - first_index).max(end_index - 1 - core_index);
    // Skip the unavailable left part of the shared profile without moving its peak
    let weights_start = radius - (core_index - first_index);
    let weights = weight_cache.get(radius)?;
    let mut counts = CoverageCounts::new();
    for (core, &weight) in cores[first_index..end_index]
        .iter()
        .zip(&weights[weights_start..])
    {
        let multiplier = u64::try_from(weight)?;
        for (&coverage, &positions) in &core.counts {
            let fitting_frequency = positions
                .checked_mul(multiplier)
                .context("triangular fitting-frequency multiplication overflow")?;
            let combined_frequency = counts.entry(coverage).or_default();
            *combined_frequency = combined_frequency
                .checked_add(fitting_frequency)
                .context("triangular fitting-frequency sum overflow")?;
        }
    }
    Ok(counts)
}

/// Mix uniform and triangular histograms after normalizing their total fitting frequencies.
///
/// For fraction `fraction`, each effective frequency is `(1-fraction)*raw +
/// fraction*triangular*eligible_positions/triangular_total`. This keeps the requested mixture
/// independent of masking or the absolute triangle multipliers. Both stages then fit this
/// same histogram, with the second retaining only coverage below the context's T1.
///
/// Intermediate blends use one million integer units per effective base. Round frequencies
/// only after pooling, preserving each positive coverage frequency with at least one unit.
/// The rounding error is at most one millionth of a base per coverage level. Coverage levels
/// themselves never change. Exact endpoints reuse the original integer histograms.
fn blend_context_counts(
    uniform_counts: &CoverageCounts,
    mut triangular_counts: CoverageCounts,
    eligible_positions: u64,
    fraction: f64,
) -> Result<CoverageCounts> {
    if fraction == 0.0 {
        return Ok(uniform_counts.clone());
    }
    if fraction == 1.0 {
        return Ok(triangular_counts);
    }
    let triangular_total = triangular_counts
        .values()
        .try_fold(0_u64, |total, &frequency| {
            total
                .checked_add(frequency)
                .context("triangular fitting-frequency total overflow")
        })?;
    ensure!(
        eligible_positions > 0 && triangular_total > 0,
        "cannot blend empty context histograms"
    );
    let triangular_normalization = eligible_positions as f64 / triangular_total as f64;
    let mut discretized_total = 0_u64;
    for (&coverage, frequency) in &mut triangular_counts {
        let observed_positions = *uniform_counts
            .get(&coverage)
            .context("triangular coverage level missing from raw context histogram")?;
        let blended_frequency = ((1.0 - fraction) * observed_positions as f64
            + fraction * triangular_normalization * *frequency as f64)
            * BLENDED_FREQUENCY_SCALE;
        ensure!(
            blended_frequency.is_finite()
                && blended_frequency >= 0.0
                && blended_frequency < u64::MAX as f64,
            "invalid blended fitting frequency {blended_frequency} at coverage {coverage}"
        );
        // Never erase a rare observed level merely because its distance weight is small
        *frequency = (blended_frequency.round() as u64).max(u64::from(observed_positions > 0));
        discretized_total = discretized_total
            .checked_add(*frequency)
            .context("blended fitting-frequency total overflow")?;
    }
    Ok(triangular_counts)
}

/// Grow a contiguous chromosome context until both eligible bases and midpoint fragments suffice.
///
/// Expand by whole cores on each available side. Masked positions retain their genomic span but
/// contribute no observations. Eligible positions with zero coverage still count toward the minimum.
/// The immutable per-core counts are summed incrementally, avoiding repeated histogram scans.
fn build_context(
    cores: &[CoreHistogram],
    eligible_positions_by_core: &[u64],
    core_index: usize,
    minimum_eligible_positions: u32,
    minimum_fragment_support: u64,
) -> Result<(CoreHistogram, Interval<u32>, u64, bool)> {
    let mut left_index = core_index;
    let mut right_index = core_index;
    let mut combined = CoreHistogram::new(cores[core_index].interval);
    combined.add_histogram(&cores[core_index])?;
    let mut eligible_positions = eligible_positions_by_core[core_index];

    loop {
        let context = Interval::new(
            cores[left_index].interval.start(),
            cores[right_index].interval.end(),
        )?;
        let has_required_positions = eligible_positions >= u64::from(minimum_eligible_positions);
        let has_required_fragments = combined.fragment_support >= minimum_fragment_support;
        if has_required_positions && has_required_fragments {
            return Ok((combined, context, eligible_positions, true));
        }

        let mut expanded = false;
        if left_index > 0 {
            left_index -= 1;
            combined.add_histogram(&cores[left_index])?;
            eligible_positions = eligible_positions
                .checked_add(eligible_positions_by_core[left_index])
                .context("context eligible-position count overflow")?;
            expanded = true;
        }
        if right_index + 1 < cores.len() {
            right_index += 1;
            combined.add_histogram(&cores[right_index])?;
            eligible_positions = eligible_positions
                .checked_add(eligible_positions_by_core[right_index])
                .context("context eligible-position count overflow")?;
            expanded = true;
        }
        if !expanded {
            return Ok((combined, context, eligible_positions, false));
        }
    }
}

#[cfg(test)]
mod tests {
    include!("striding_tests.rs");
}
