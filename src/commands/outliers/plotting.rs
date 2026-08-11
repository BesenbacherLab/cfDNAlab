use std::path::Path;

use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use plotters::prelude::*;

use crate::commands::outliers::{
    model::{CoverageCounts, TwoStageZipModel},
    striding::CoreModel,
};

#[derive(Debug, Clone, Copy, PartialEq)]
struct ThresholdSummary {
    minimum: u32,
    mean_across_cores: f64,
    maximum: u32,
}

#[derive(Debug)]
struct LocalModelDiagnostics {
    initial_thresholds: Vec<f64>,
    final_thresholds: Vec<f64>,
    underlying_means: Vec<f64>,
    chromosome_boundaries: Vec<usize>,
    initial_summary: ThresholdSummary,
    final_summary: ThresholdSummary,
}

/// Write global fit diagnostics together with the local models used during calling.
///
/// The upper panels show the complete original global histogram and both global ZIP fits. The
/// global thresholds are sample-level references only. Local threshold summaries and the lower
/// spatial panel describe the per-core models that actually control outlier calling.
pub(crate) fn write_global_fit_plot(
    output_path: &Path,
    observed_counts: &CoverageCounts,
    global_model: TwoStageZipModel,
    chromosomes: &[String],
    models_by_chromosome: &FxHashMap<String, Vec<CoreModel>>,
) -> Result<()> {
    let maximum_observed_coverage = observed_counts
        .last_key_value()
        .map(|(&coverage, _)| coverage)
        .context("global ZIP plot has no observed coverage values")?;
    ensure!(
        maximum_observed_coverage > 0,
        "global ZIP plot has no positive coverage values"
    );
    let local_diagnostics = collect_local_model_diagnostics(chromosomes, models_by_chromosome)?;
    let maximum_threshold = global_model
        .initial_threshold
        .max(global_model.final_threshold)
        .max(local_diagnostics.initial_summary.maximum)
        .max(local_diagnostics.final_summary.maximum);
    let maximum_diagnostic_coverage = maximum_observed_coverage.max(maximum_threshold);
    let detailed_coverage_limit = maximum_threshold
        .saturating_add((maximum_threshold / 3).max(5))
        .min(maximum_diagnostic_coverage)
        .max(1);
    let total_positions = observed_counts.values().copied().sum::<u64>() as f64;
    let maximum_log_count = (0..=maximum_diagnostic_coverage)
        .flat_map(|coverage| {
            let observed = observed_counts.get(&coverage).copied().unwrap_or(0) as f64;
            [
                (observed + 1.0).log10(),
                (total_positions * global_model.initial.probability_mass(coverage) + 1.0).log10(),
                (total_positions * global_model.underlying_fit.probability_mass(coverage) + 1.0)
                    .log10(),
            ]
        })
        .fold(0.0_f64, f64::max)
        .max(1.0);

    let root = BitMapBackend::new(output_path, (1600, 1500)).into_drawing_area();
    root.fill(&WHITE)?;
    let panels = root.split_evenly((3, 1));

    draw_detailed_fit_panel(
        &panels[0],
        observed_counts,
        total_positions,
        global_model,
        &local_diagnostics,
        detailed_coverage_limit,
        maximum_log_count,
    )?;
    draw_complete_histogram_panel(
        &panels[1],
        observed_counts,
        total_positions,
        global_model,
        &local_diagnostics,
        maximum_diagnostic_coverage,
        maximum_log_count,
    )?;
    draw_local_model_panel(&panels[2], &local_diagnostics)?;

    root.present().with_context(|| {
        format!(
            "writing outlier diagnostic plot to {}",
            output_path.display()
        )
    })?;
    Ok(())
}

fn collect_local_model_diagnostics(
    chromosomes: &[String],
    models_by_chromosome: &FxHashMap<String, Vec<CoreModel>>,
) -> Result<LocalModelDiagnostics> {
    let mut initial_thresholds = Vec::new();
    let mut final_thresholds = Vec::new();
    let mut underlying_means = Vec::new();
    let mut chromosome_boundaries = Vec::new();

    for chromosome in chromosomes {
        let models = models_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing models for chromosome '{chromosome}'"))?;
        if !initial_thresholds.is_empty() && !models.is_empty() {
            chromosome_boundaries.push(initial_thresholds.len());
        }
        for model in models {
            initial_thresholds.push(model.zip.initial_threshold as f64);
            final_thresholds.push(model.zip.final_threshold as f64);
            underlying_means.push(model.zip.underlying_fit.mean());
        }
    }

    ensure!(
        !initial_thresholds.is_empty(),
        "cannot plot local ZIP diagnostics without core models"
    );
    let initial_summary = summarize_thresholds(&initial_thresholds);
    let final_summary = summarize_thresholds(&final_thresholds);
    Ok(LocalModelDiagnostics {
        initial_thresholds,
        final_thresholds,
        underlying_means,
        chromosome_boundaries,
        initial_summary,
        final_summary,
    })
}

fn summarize_thresholds(thresholds: &[f64]) -> ThresholdSummary {
    let minimum = thresholds.iter().copied().fold(f64::INFINITY, f64::min);
    let maximum = thresholds.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    ThresholdSummary {
        minimum: minimum as u32,
        mean_across_cores: thresholds.iter().sum::<f64>() / thresholds.len() as f64,
        maximum: maximum as u32,
    }
}

fn draw_detailed_fit_panel(
    area: &DrawingArea<BitMapBackend<'_>, plotters::coord::Shift>,
    observed_counts: &CoverageCounts,
    total_positions: f64,
    global_model: TwoStageZipModel,
    local_diagnostics: &LocalModelDiagnostics,
    maximum_coverage: u32,
    maximum_log_count: f64,
) -> Result<()> {
    let x_axis_end = maximum_coverage as f64 + 1.0;
    let y_axis_end = maximum_log_count * 1.05;
    let mut chart = ChartBuilder::on(area)
        .caption(
            "Global raw-coverage ZIP fit and threshold ranges",
            ("sans-serif", 30),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(0.0_f64..x_axis_end, 0.0_f64..y_axis_end)?;
    chart
        .configure_mesh()
        .x_desc("Raw positional coverage")
        .y_desc("log10(position count + 1)")
        .draw()?;

    let observed_series = (0..=maximum_coverage).map(|coverage| {
        let count = observed_counts.get(&coverage).copied().unwrap_or(0) as f64;
        (coverage as f64, (count + 1.0).log10())
    });
    chart
        .draw_series(LineSeries::new(observed_series, BLACK.stroke_width(3)))?
        .label("Observed")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], BLACK.stroke_width(3)));

    let initial_series = (0..=maximum_coverage).map(|coverage| {
        (
            coverage as f64,
            (total_positions * global_model.initial.probability_mass(coverage) + 1.0).log10(),
        )
    });
    chart
        .draw_series(LineSeries::new(initial_series, BLUE.stroke_width(2)))?
        .label("Global initial ZIP reference")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], BLUE.stroke_width(2)));

    let final_series = (0..=maximum_coverage).map(|coverage| {
        (
            coverage as f64,
            (total_positions * global_model.underlying_fit.probability_mass(coverage) + 1.0)
                .log10(),
        )
    });
    chart
        .draw_series(LineSeries::new(final_series, RED.stroke_width(2)))?
        .label("Global tail-excluded ZIP reference")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], RED.stroke_width(2)));

    chart
        .draw_series(std::iter::once(PathElement::new(
            [
                (global_model.initial_threshold as f64, 0.0),
                (global_model.initial_threshold as f64, y_axis_end),
            ],
            BLUE.mix(0.45).stroke_width(1),
        )))?
        .label(format!(
            "Global reference T1 = {}",
            global_model.initial_threshold
        ))
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], BLUE.mix(0.45)));
    chart
        .draw_series(std::iter::once(PathElement::new(
            [
                (global_model.final_threshold as f64, 0.0),
                (global_model.final_threshold as f64, y_axis_end),
            ],
            RED.mix(0.45).stroke_width(1),
        )))?
        .label(format!(
            "Global reference T2 = {}",
            global_model.final_threshold
        ))
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], RED.mix(0.45)));

    for threshold in [
        local_diagnostics.initial_summary.minimum as f64,
        local_diagnostics.initial_summary.maximum as f64,
    ] {
        chart.draw_series(std::iter::once(PathElement::new(
            [(threshold, 0.0), (threshold, y_axis_end)],
            CYAN.mix(0.45).stroke_width(1),
        )))?;
    }
    chart
        .draw_series(std::iter::once(PathElement::new(
            [
                (local_diagnostics.initial_summary.mean_across_cores, 0.0),
                (
                    local_diagnostics.initial_summary.mean_across_cores,
                    y_axis_end,
                ),
            ],
            CYAN.stroke_width(3),
        )))?
        .label(format!(
            "Used T1 across cores min/mean/max = {}/{:.1}/{}",
            local_diagnostics.initial_summary.minimum,
            local_diagnostics.initial_summary.mean_across_cores,
            local_diagnostics.initial_summary.maximum
        ))
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], CYAN.stroke_width(3)));

    for threshold in [
        local_diagnostics.final_summary.minimum as f64,
        local_diagnostics.final_summary.maximum as f64,
    ] {
        chart.draw_series(std::iter::once(PathElement::new(
            [(threshold, 0.0), (threshold, y_axis_end)],
            MAGENTA.mix(0.45).stroke_width(1),
        )))?;
    }
    chart
        .draw_series(std::iter::once(PathElement::new(
            [
                (local_diagnostics.final_summary.mean_across_cores, 0.0),
                (
                    local_diagnostics.final_summary.mean_across_cores,
                    y_axis_end,
                ),
            ],
            MAGENTA.stroke_width(3),
        )))?
        .label(format!(
            "Used T2 across cores min/mean/max = {}/{:.1}/{}",
            local_diagnostics.final_summary.minimum,
            local_diagnostics.final_summary.mean_across_cores,
            local_diagnostics.final_summary.maximum
        ))
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], MAGENTA.stroke_width(3)));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .draw()?;
    Ok(())
}

fn draw_complete_histogram_panel(
    area: &DrawingArea<BitMapBackend<'_>, plotters::coord::Shift>,
    observed_counts: &CoverageCounts,
    total_positions: f64,
    global_model: TwoStageZipModel,
    local_diagnostics: &LocalModelDiagnostics,
    maximum_coverage: u32,
    maximum_log_count: f64,
) -> Result<()> {
    let transformed_maximum = ((maximum_coverage as f64) + 1.0).log10();
    let y_axis_end = maximum_log_count * 1.05;
    let mut chart = ChartBuilder::on(area)
        .caption("Complete original histogram", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(0.0_f64..transformed_maximum, 0.0_f64..y_axis_end)?;
    chart
        .configure_mesh()
        .x_desc("log10(raw positional coverage + 1)")
        .y_desc("log10(position count + 1)")
        .draw()?;

    let transform_coverage = |coverage: u32| ((coverage as f64) + 1.0).log10();
    let observed_series = (0..=maximum_coverage).filter_map(|coverage| {
        let count = observed_counts.get(&coverage).copied().unwrap_or(0);
        (count > 0).then(|| (transform_coverage(coverage), (count as f64 + 1.0).log10()))
    });
    chart.draw_series(
        observed_series.map(|point| Circle::new(point, 2, BLACK.mix(0.75).filled())),
    )?;

    let initial_series = (0..=maximum_coverage).map(|coverage| {
        (
            transform_coverage(coverage),
            (total_positions * global_model.initial.probability_mass(coverage) + 1.0).log10(),
        )
    });
    chart.draw_series(LineSeries::new(
        initial_series,
        BLUE.mix(0.7).stroke_width(2),
    ))?;
    let final_series = (0..=maximum_coverage).map(|coverage| {
        (
            transform_coverage(coverage),
            (total_positions * global_model.underlying_fit.probability_mass(coverage) + 1.0)
                .log10(),
        )
    });
    chart.draw_series(LineSeries::new(final_series, RED.mix(0.7).stroke_width(2)))?;

    for (threshold, style) in [
        (
            global_model.initial_threshold as f64,
            BLUE.mix(0.45).stroke_width(1),
        ),
        (
            global_model.final_threshold as f64,
            RED.mix(0.45).stroke_width(1),
        ),
        (
            local_diagnostics.initial_summary.minimum as f64,
            CYAN.mix(0.45).stroke_width(1),
        ),
        (
            local_diagnostics.initial_summary.mean_across_cores,
            CYAN.stroke_width(3),
        ),
        (
            local_diagnostics.initial_summary.maximum as f64,
            CYAN.mix(0.45).stroke_width(1),
        ),
        (
            local_diagnostics.final_summary.minimum as f64,
            MAGENTA.mix(0.45).stroke_width(1),
        ),
        (
            local_diagnostics.final_summary.mean_across_cores,
            MAGENTA.stroke_width(3),
        ),
        (
            local_diagnostics.final_summary.maximum as f64,
            MAGENTA.mix(0.45).stroke_width(1),
        ),
    ] {
        let transformed_threshold = (threshold + 1.0).log10();
        chart.draw_series(std::iter::once(PathElement::new(
            [
                (transformed_threshold, 0.0),
                (transformed_threshold, y_axis_end),
            ],
            style,
        )))?;
    }
    Ok(())
}

fn draw_local_model_panel(
    area: &DrawingArea<BitMapBackend<'_>, plotters::coord::Shift>,
    diagnostics: &LocalModelDiagnostics,
) -> Result<()> {
    let number_of_cores = diagnostics.initial_thresholds.len();
    let maximum_value = diagnostics
        .initial_thresholds
        .iter()
        .chain(&diagnostics.final_thresholds)
        .chain(&diagnostics.underlying_means)
        .copied()
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let mut chart = ChartBuilder::on(area)
        .caption("Local models used during calling", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(70)
        .build_cartesian_2d(
            0.0_f64..number_of_cores as f64,
            0.0_f64..maximum_value * 1.05,
        )?;
    chart
        .configure_mesh()
        .x_desc("Fixed-size model cores in chromosome order")
        .y_desc("Raw positional coverage")
        .draw()?;

    for &boundary in &diagnostics.chromosome_boundaries {
        chart.draw_series(std::iter::once(PathElement::new(
            [
                (boundary as f64, 0.0),
                (boundary as f64, maximum_value * 1.05),
            ],
            BLACK.mix(0.15).stroke_width(1),
        )))?;
    }

    let initial_points = diagnostics
        .initial_thresholds
        .iter()
        .enumerate()
        .map(|(index, &threshold)| (index as f64 + 0.5, threshold));
    chart
        .draw_series(initial_points.map(|point| Circle::new(point, 2, CYAN.filled())))?
        .label("Used T1")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], CYAN.stroke_width(2)));
    let final_points = diagnostics
        .final_thresholds
        .iter()
        .enumerate()
        .map(|(index, &threshold)| (index as f64 + 0.5, threshold));
    chart
        .draw_series(final_points.map(|point| Circle::new(point, 2, MAGENTA.filled())))?
        .label("Used T2")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], MAGENTA.stroke_width(2)));
    let mean_points = diagnostics
        .underlying_means
        .iter()
        .enumerate()
        .map(|(index, &mean)| (index as f64 + 0.5, mean));
    chart
        .draw_series(mean_points.map(|point| Circle::new(point, 2, GREEN.filled())))?
        .label("Used underlying ZIP mean")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 30, y)], GREEN.stroke_width(2)));
    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .draw()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("plotting_tests.rs");
}
