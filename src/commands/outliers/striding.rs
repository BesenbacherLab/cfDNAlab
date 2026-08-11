use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use rayon::prelude::*;

use crate::{
    commands::outliers::model::{
        CoverageCounts, TwoStageZipModel, ZipFitDiagnostics, diagnose_zip_fit, fit_two_stage_zip,
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
    /// Accepted fragments with an unblacklisted midpoint, assigned once to this core.
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
    pub(crate) core: Interval<u32>,
    pub(crate) context: Interval<u32>,
    pub(crate) eligible_positions: u64,
    pub(crate) fragment_support: u64,
    pub(crate) source: ModelSource,
    pub(crate) zip: TwoStageZipModel,
    pub(crate) diagnostics: ZipFitDiagnostics,
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

/// Fit every core from an adaptively expanded, chromosome-local histogram context.
///
/// The global model is used only when the complete chromosome cannot meet the required span and
/// fragment support. A complete local context whose two-stage fit fails returns that error rather
/// than silently changing its model source.
pub(crate) fn fit_core_models(
    chromosomes: &[String],
    histograms_by_chromosome: &FxHashMap<String, Vec<CoreHistogram>>,
    global_model: TwoStageZipModel,
    tail_probability: f64,
    minimum_context_span: u32,
    minimum_fragment_support: u64,
) -> Result<FxHashMap<String, Vec<CoreModel>>> {
    let mut models_by_chromosome =
        FxHashMap::with_capacity_and_hasher(histograms_by_chromosome.len(), Default::default());

    for chromosome in chromosomes {
        let cores = histograms_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing core histograms for chromosome '{chromosome}'"))?;
        ensure!(
            !cores.is_empty(),
            "chromosome '{}' has no model-core histograms",
            chromosome
        );
        let chromosome_models = (0..cores.len())
            .into_par_iter()
            .map(|core_index| {
                let (context_histogram, context, context_complete) = build_context(
                    cores,
                    core_index,
                    minimum_context_span,
                    minimum_fragment_support,
                )?;
                let local_fit = context_complete
                    .then(|| fit_two_stage_zip(&context_histogram.counts, tail_probability))
                    .transpose()?;
                let (source, zip) = match local_fit {
                    Some(model) => (ModelSource::Local, model),
                    None => (ModelSource::GlobalFallback, global_model),
                };
                let diagnostics = diagnose_zip_fit(&context_histogram.counts, zip)?;
                Ok(CoreModel {
                    core: cores[core_index].interval,
                    context,
                    eligible_positions: context_histogram.eligible_positions()?,
                    fragment_support: context_histogram.fragment_support,
                    source,
                    zip,
                    diagnostics,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        models_by_chromosome.insert(chromosome.clone(), chromosome_models);
    }

    Ok(models_by_chromosome)
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
