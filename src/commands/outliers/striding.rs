use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use rayon::prelude::*;

use crate::{
    commands::outliers::{
        config::OutlierModel,
        distribution::{FitDiagnostics, TwoStageCoverageModel, diagnose_fit, fit_two_stage},
        failure::AnalysisFailures,
        model::CoverageCounts,
    },
    shared::interval::Interval,
};

const MAXIMUM_EXACT_RAW_COVERAGE: f32 = 16_777_215.0;

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
/// The global model is used only when the complete chromosome cannot meet the required span and
/// fragment support. A complete local context whose two-stage fit fails returns that error rather
/// than silently changing its model source.
pub(crate) fn fit_core_models(
    chromosomes: &[String],
    histograms_by_chromosome: &FxHashMap<String, Vec<CoreHistogram>>,
    global_histogram: &CoreHistogram,
    family: OutlierModel,
    tail_probability: f64,
    minimum_context_span: u32,
    minimum_fragment_support: u64,
) -> Result<FittedModels> {
    let mut needs_global = false;
    for chromosome in chromosomes {
        let cores = histograms_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing core histograms for {chromosome}"))?;
        needs_global |=
            needs_global_fallback(cores, minimum_context_span, minimum_fragment_support)?;
    }
    let mut messages = Vec::new();
    let mut failed_histograms = String::from(
        "# GLOBAL pools all selected chromosomes and has no genomic coordinates (NA)\n# coverage is raw positional coverage, observed_positions counts eligible bases\nchromosome\tcore_start\tcore_end\tcontext_start\tcontext_end\tfragments_by_midpoint\tcoverage\tobserved_positions\n",
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
                );
                None
            }
        }
    } else {
        None
    };
    let mut models_by_chromosome = FxHashMap::default();
    for chromosome in chromosomes {
        let cores = &histograms_by_chromosome[chromosome];
        // Collect every independent local result before returning failures. Histogram processing
        // is cheap compared with another BAM pass and supplies all failing contexts at once
        let results = (0..cores.len())
            .into_par_iter()
            .map(|core_index| {
                let core = cores[core_index].interval;
                if cores[core_index].eligible_positions()? == 0 {
                    // Keep this grid slot without constructing or fitting an unused context
                    return Ok((core, core, CoreHistogram::new(core), None));
                }
                let (mut histogram, context, complete) = build_context(
                    cores,
                    core_index,
                    minimum_context_span,
                    minimum_fragment_support,
                )?;
                let fit = if complete {
                    Some(fit_two_stage(&histogram.counts, tail_probability, family))
                } else {
                    global_model.map(Ok)
                };
                let result = fit.map(|fit| -> Result<CoreModel> {
                    let model = fit?;
                    Ok(CoreModel {
                        core,
                        context,
                        eligible_positions: histogram.eligible_positions()?,
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
                }
                Ok((core, context, histogram, result))
            })
            .collect::<Vec<Result<_>>>();
        let mut chromosome_models = Vec::with_capacity(cores.len());
        for result in results {
            match result {
                Ok((core, context, histogram, Some(Err(error)))) => {
                    messages.push(format!("{chromosome} core {:?}, context {:?}, {} fit, alpha={tail_probability}, eligible positions {}, fragment support {}: {error:#}", core.as_tuple(), context.as_tuple(), family.as_str(), histogram.eligible_positions()?, histogram.fragment_support));
                    append_failed_histogram(
                        &mut failed_histograms,
                        chromosome,
                        core,
                        context,
                        &histogram,
                    );
                    chromosome_models.push(None);
                }
                Ok((_, _, _, Some(Ok(model)))) => chromosome_models.push(Some(model)),
                // A fully blacklisted core or an already recorded required global failure
                Ok((_, _, _, None)) => chromosome_models.push(None),
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

/// A full chromosome supplies the maximum possible local span and midpoint support.
///
/// A wholly blacklisted chromosome has no positions needing a model, so it cannot require fallback.
fn needs_global_fallback(
    cores: &[CoreHistogram],
    minimum_span: u32,
    minimum_support: u64,
) -> Result<bool> {
    let first = cores
        .first()
        .context("chromosome has no model-core histograms")?;
    let last = cores
        .last()
        .context("chromosome has no model-core histograms")?;
    if !cores
        .iter()
        .any(|core| core.counts.values().any(|&positions| positions > 0))
    {
        return Ok(false);
    }
    let span = first.interval.expand_to_include(last.interval).len();
    let support = cores.iter().try_fold(0_u64, |total, core| {
        total
            .checked_add(core.fragment_support)
            .context("chromosome fragment support overflow")
    })?;
    Ok(span < minimum_span || support < minimum_support)
}

/// Preserve sparse failing-context counts without requiring decompression or a whole-run scan.
fn append_failed_histogram(
    output: &mut String,
    chromosome: &str,
    core: Interval<u32>,
    context: Interval<u32>,
    histogram: &CoreHistogram,
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
            "{chromosome}\t{coordinates}\t{}\t{coverage}\t{positions}",
            histogram.fragment_support
        )
        .expect("writing to a string cannot fail");
    }
}

fn build_context(
    cores: &[CoreHistogram],
    core_index: usize,
    minimum_context_span: u32,
    minimum_fragment_support: u64,
) -> Result<(CoreHistogram, Interval<u32>, bool)> {
    let mut left_index = core_index;
    let mut right_index = core_index;
    let mut combined = CoreHistogram::new(cores[core_index].interval);
    combined.add_histogram(&cores[core_index])?;

    loop {
        let context = Interval::new(
            cores[left_index].interval.start(),
            cores[right_index].interval.end(),
        )?;
        let has_required_span = context.len() >= minimum_context_span;
        let has_required_fragments = combined.fragment_support >= minimum_fragment_support;
        if has_required_span && has_required_fragments {
            return Ok((combined, context, true));
        }

        let mut expanded = false;
        if left_index > 0 {
            left_index -= 1;
            combined.add_histogram(&cores[left_index])?;
            expanded = true;
        }
        if right_index + 1 < cores.len() {
            right_index += 1;
            combined.add_histogram(&cores[right_index])?;
            expanded = true;
        }
        if !expanded {
            return Ok((combined, context, false));
        }
    }
}

#[cfg(test)]
mod tests {
    include!("striding_tests.rs");
}
