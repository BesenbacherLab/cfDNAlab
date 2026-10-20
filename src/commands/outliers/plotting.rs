use std::path::Path;

use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};

use crate::commands::outliers::{
    distribution::TwoStageCoverageModel, model::CoverageCounts, striding::CoreModel,
};

const INITIAL_THRESHOLD_COLOR: RGBColor = RGBColor(0, 130, 160);
const FINAL_THRESHOLD_COLOR: RGBColor = RGBColor(180, 0, 180);

#[derive(Debug, Clone, Copy, PartialEq)]
struct ThresholdSummary {
    minimum: u32,
    mean_across_cores: f64,
    maximum: u32,
}

#[derive(Debug)]
struct LocalModelDiagnostics {
    initial_thresholds: Vec<Option<f64>>,
    final_thresholds: Vec<Option<f64>>,
    underlying_means: Vec<Option<f64>>,
    active_component_means: Vec<Option<f64>>,
    chromosome_boundaries: Vec<usize>,
    chromosome_labels: Vec<(usize, String)>,
    initial_summary: ThresholdSummary,
    final_summary: ThresholdSummary,
}

/// Write global fit diagnostics together with the local models used during calling.
///
/// Both histogram panels show the same raw eligible-base counts, with linear and logarithmic
/// coverage axes. Fitted curves appear only when a chromosome needs global fallback. Their
/// thresholds are sample-level references. The lower panels show the cutoffs and means used for
/// calling on separate scales so extreme cutoffs cannot hide variation in the fitted means.
/// Fully blacklisted cores remain gaps at their original positions and do not enter summaries.
pub(crate) fn write_global_fit_plot(
    output_path: &Path,
    observed_counts: &CoverageCounts,
    global_model: Option<TwoStageCoverageModel>,
    chromosomes: &[String],
    models_by_chromosome: &FxHashMap<String, Vec<Option<CoreModel>>>,
) -> Result<()> {
    let maximum_observed_coverage = observed_counts
        .last_key_value()
        .map(|(&coverage, _)| coverage)
        .context("global coverage plot has no observed coverage values")?;
    ensure!(
        maximum_observed_coverage > 0,
        "global coverage plot has no positive coverage values"
    );
    let local_diagnostics = collect_local_model_diagnostics(chromosomes, models_by_chromosome)?;
    let maximum_threshold = global_model
        .map(|model| model.initial_threshold.max(model.final_threshold))
        .unwrap_or(0)
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
                global_model
                    .map(|model| {
                        (total_positions * model.initial.probability_mass(coverage) + 1.0).log10()
                    })
                    .unwrap_or(0.0),
                global_model
                    .map(|model| {
                        (total_positions * model.underlying_fit.probability_mass(coverage) + 1.0)
                            .log10()
                    })
                    .unwrap_or(0.0),
            ]
        })
        .fold(0.0_f64, f64::max)
        .max(1.0);

    let root = BitMapBackend::new(output_path, (1600, 2000)).into_drawing_area();
    root.fill(&WHITE)?;
    let panels = root.split_evenly((4, 1));

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
    draw_local_model_panel(
        &panels[2],
        &local_diagnostics,
        "Local calling cutoffs",
        [
            (
                "T1: initial exclusion cutoff",
                &local_diagnostics.initial_thresholds,
                INITIAL_THRESHOLD_COLOR,
            ),
            (
                "T2: final calling cutoff",
                &local_diagnostics.final_thresholds,
                FINAL_THRESHOLD_COLOR,
            ),
        ],
    )?;
    draw_local_model_panel(
        &panels[3],
        &local_diagnostics,
        "Local fitted coverage means",
        [
            (
                "Mean including structural zeros",
                &local_diagnostics.underlying_means,
                RGBColor(0, 130, 0),
            ),
            (
                "Active-component mean",
                &local_diagnostics.active_component_means,
                BLUE,
            ),
        ],
    )?;

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
    models_by_chromosome: &FxHashMap<String, Vec<Option<CoreModel>>>,
) -> Result<LocalModelDiagnostics> {
    let mut initial_thresholds = Vec::new();
    let mut final_thresholds = Vec::new();
    let mut underlying_means = Vec::new();
    let mut active_component_means = Vec::new();
    let mut chromosome_boundaries = Vec::new();
    let mut chromosome_labels = Vec::new();

    for chromosome in chromosomes {
        let models = models_by_chromosome
            .get(chromosome)
            .with_context(|| format!("missing models for chromosome '{chromosome}'"))?;
        if !initial_thresholds.is_empty() && !models.is_empty() {
            chromosome_boundaries.push(initial_thresholds.len());
        }
        chromosome_labels.push((initial_thresholds.len(), chromosome.clone()));
        for model in models {
            // Keep one point slot per grid core, with no value for a fully blacklisted core
            let model = model.as_ref();
            initial_thresholds.push(model.map(|model| model.model.initial_threshold as f64));
            final_thresholds.push(model.map(|model| model.model.final_threshold as f64));
            underlying_means.push(model.map(|model| model.model.underlying_fit.mean()));
            active_component_means
                .push(model.map(|model| model.model.underlying_fit.component_mean));
        }
    }

    let initial_summary = summarize_thresholds(&initial_thresholds)?;
    let final_summary = summarize_thresholds(&final_thresholds)?;
    Ok(LocalModelDiagnostics {
        initial_thresholds,
        final_thresholds,
        underlying_means,
        active_component_means,
        chromosome_boundaries,
        chromosome_labels,
        initial_summary,
        final_summary,
    })
}

/// Summarize actual fitted cutoffs without treating fully blacklisted cores as observations.
fn summarize_thresholds(thresholds: &[Option<f64>]) -> Result<ThresholdSummary> {
    let values = thresholds.iter().flatten().copied();
    let number_of_models = values.clone().count();
    ensure!(
        number_of_models > 0,
        "cannot plot local coverage diagnostics without fitted core models"
    );
    let minimum = values.clone().fold(f64::INFINITY, f64::min);
    let maximum = values.clone().fold(f64::NEG_INFINITY, f64::max);
    Ok(ThresholdSummary {
        minimum: minimum as u32,
        mean_across_cores: values.sum::<f64>() / number_of_models as f64,
        maximum: maximum as u32,
    })
}

fn draw_detailed_fit_panel(
    area: &DrawingArea<BitMapBackend<'_>, plotters::coord::Shift>,
    observed_counts: &CoverageCounts,
    total_positions: f64,
    global_model: Option<TwoStageCoverageModel>,
    local_diagnostics: &LocalModelDiagnostics,
    maximum_coverage: u32,
    maximum_log_count: f64,
) -> Result<()> {
    let x_axis_end = maximum_coverage as f64 + 1.0;
    let y_axis_end = maximum_log_count * 1.05;
    let mut chart = ChartBuilder::on(area)
        .caption(
            "Raw coverage histogram (linear coverage axis)",
            ("sans-serif", 32),
        )
        .margin(20)
        .x_label_area_size(75)
        .y_label_area_size(100)
        .build_cartesian_2d(0.0_f64..x_axis_end, 0.0_f64..y_axis_end)?;
    chart
        .configure_mesh()
        .label_style(("sans-serif", 22))
        .axis_desc_style(("sans-serif", 24))
        .max_light_lines(0)
        .bold_line_style(BLACK.mix(0.12))
        .x_desc("Raw positional fragment coverage")
        .y_desc("log10(eligible bases + 1)")
        .draw()?;

    let observed_series = (0..=maximum_coverage).map(|coverage| {
        let count = observed_counts.get(&coverage).copied().unwrap_or(0) as f64;
        (coverage as f64, (count + 1.0).log10())
    });
    chart
        .draw_series(LineSeries::new(observed_series, BLACK.stroke_width(3)))?
        .label("Raw eligible-base counts")
        .legend(|(left, top)| {
            PathElement::new([(left, top), (left + 30, top)], BLACK.stroke_width(3))
        });

    if let Some(global_model) = global_model {
        let initial_series = (0..=maximum_coverage).map(|coverage| {
            (
                coverage as f64,
                (total_positions * global_model.initial.probability_mass(coverage) + 1.0).log10(),
            )
        });
        chart
            .draw_series(LineSeries::new(initial_series, BLUE.stroke_width(2)))?
            .label(format!(
                "Global initial {} (dispersion {:.4})",
                global_model.family.as_str().to_uppercase(),
                global_model.initial.dispersion
            ))
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
            .label(format!(
                "Global tail-excluded {} (dispersion {:.4})",
                global_model.family.as_str().to_uppercase(),
                global_model.underlying_fit.dispersion
            ))
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
    }

    for threshold in [
        local_diagnostics.initial_summary.minimum as f64,
        local_diagnostics.initial_summary.maximum as f64,
    ] {
        chart.draw_series(std::iter::once(PathElement::new(
            [(threshold, 0.0), (threshold, y_axis_end)],
            INITIAL_THRESHOLD_COLOR.mix(0.45).stroke_width(1),
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
            INITIAL_THRESHOLD_COLOR.stroke_width(3),
        )))?
        .label(format!(
            "T1 cutoff: mean {:.1} (thick), min/max {}/{} (thin)",
            local_diagnostics.initial_summary.mean_across_cores,
            local_diagnostics.initial_summary.minimum,
            local_diagnostics.initial_summary.maximum
        ))
        .legend(|(left, top)| {
            PathElement::new(
                [(left, top), (left + 30, top)],
                INITIAL_THRESHOLD_COLOR.stroke_width(3),
            )
        });

    for threshold in [
        local_diagnostics.final_summary.minimum as f64,
        local_diagnostics.final_summary.maximum as f64,
    ] {
        chart.draw_series(std::iter::once(PathElement::new(
            [(threshold, 0.0), (threshold, y_axis_end)],
            FINAL_THRESHOLD_COLOR.mix(0.45).stroke_width(1),
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
            FINAL_THRESHOLD_COLOR.stroke_width(3),
        )))?
        .label(format!(
            "T2 cutoff: mean {:.1} (thick), min/max {}/{} (thin)",
            local_diagnostics.final_summary.mean_across_cores,
            local_diagnostics.final_summary.minimum,
            local_diagnostics.final_summary.maximum
        ))
        .legend(|(left, top)| {
            PathElement::new(
                [(left, top), (left + 30, top)],
                FINAL_THRESHOLD_COLOR.stroke_width(3),
            )
        });

    chart
        .configure_series_labels()
        .label_font(("sans-serif", 22))
        .position(SeriesLabelPosition::UpperRight)
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .draw()?;
    Ok(())
}

fn draw_complete_histogram_panel(
    area: &DrawingArea<BitMapBackend<'_>, plotters::coord::Shift>,
    observed_counts: &CoverageCounts,
    total_positions: f64,
    global_model: Option<TwoStageCoverageModel>,
    local_diagnostics: &LocalModelDiagnostics,
    maximum_coverage: u32,
    maximum_log_count: f64,
) -> Result<()> {
    let transformed_maximum = ((maximum_coverage as f64) + 1.0).log10();
    let y_axis_end = maximum_log_count * 1.05;
    let mut chart = ChartBuilder::on(area)
        .caption(
            "Raw coverage histogram (log coverage axis)",
            ("sans-serif", 32),
        )
        .margin(20)
        .x_label_area_size(75)
        .y_label_area_size(100)
        .build_cartesian_2d(0.0_f64..transformed_maximum, 0.0_f64..y_axis_end)?;
    chart
        .configure_mesh()
        .label_style(("sans-serif", 22))
        .axis_desc_style(("sans-serif", 24))
        .max_light_lines(0)
        .bold_line_style(BLACK.mix(0.12))
        .x_desc("log10(raw positional fragment coverage + 1)")
        .y_desc("log10(eligible bases + 1)")
        .draw()?;

    let transform_coverage = |coverage: u32| ((coverage as f64) + 1.0).log10();
    let observed_series = (0..=maximum_coverage).filter_map(|coverage| {
        let count = observed_counts.get(&coverage).copied().unwrap_or(0);
        (count > 0).then(|| (transform_coverage(coverage), (count as f64 + 1.0).log10()))
    });
    chart
        .draw_series(observed_series.map(|point| Circle::new(point, 2, BLACK.mix(0.75).filled())))?
        .label("Raw eligible-base counts")
        .legend(|(left, top)| Circle::new((left + 15, top), 3, BLACK.filled()));

    if let Some(global_model) = global_model {
        let initial_series = (0..=maximum_coverage).map(|coverage| {
            (
                transform_coverage(coverage),
                (total_positions * global_model.initial.probability_mass(coverage) + 1.0).log10(),
            )
        });
        chart
            .draw_series(LineSeries::new(
                initial_series,
                BLUE.mix(0.7).stroke_width(2),
            ))?
            .label("Global initial fit / T1")
            .legend(|(left, top)| PathElement::new([(left, top), (left + 30, top)], BLUE));
        let final_series = (0..=maximum_coverage).map(|coverage| {
            (
                transform_coverage(coverage),
                (total_positions * global_model.underlying_fit.probability_mass(coverage) + 1.0)
                    .log10(),
            )
        });
        chart
            .draw_series(LineSeries::new(final_series, RED.mix(0.7).stroke_width(2)))?
            .label("Global underlying fit / T2")
            .legend(|(left, top)| PathElement::new([(left, top), (left + 30, top)], RED));

        for (threshold, color) in [
            (global_model.initial_threshold, BLUE),
            (global_model.final_threshold, RED),
        ] {
            let position = transform_coverage(threshold);
            chart.draw_series(std::iter::once(PathElement::new(
                [(position, 0.0), (position, y_axis_end)],
                color.mix(0.45),
            )))?;
        }
    }
    for (threshold, style) in [
        (
            local_diagnostics.initial_summary.minimum as f64,
            INITIAL_THRESHOLD_COLOR.mix(0.45).stroke_width(1),
        ),
        (
            local_diagnostics.initial_summary.mean_across_cores,
            INITIAL_THRESHOLD_COLOR.stroke_width(3),
        ),
        (
            local_diagnostics.initial_summary.maximum as f64,
            INITIAL_THRESHOLD_COLOR.mix(0.45).stroke_width(1),
        ),
        (
            local_diagnostics.final_summary.minimum as f64,
            FINAL_THRESHOLD_COLOR.mix(0.45).stroke_width(1),
        ),
        (
            local_diagnostics.final_summary.mean_across_cores,
            FINAL_THRESHOLD_COLOR.stroke_width(3),
        ),
        (
            local_diagnostics.final_summary.maximum as f64,
            FINAL_THRESHOLD_COLOR.mix(0.45).stroke_width(1),
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
    for (label, color) in [
        (
            format!(
                "T1 cutoff: mean {:.1} (thick), min/max {}/{} (thin)",
                local_diagnostics.initial_summary.mean_across_cores,
                local_diagnostics.initial_summary.minimum,
                local_diagnostics.initial_summary.maximum
            ),
            INITIAL_THRESHOLD_COLOR,
        ),
        (
            format!(
                "T2 cutoff: mean {:.1} (thick), min/max {}/{} (thin)",
                local_diagnostics.final_summary.mean_across_cores,
                local_diagnostics.final_summary.minimum,
                local_diagnostics.final_summary.maximum
            ),
            FINAL_THRESHOLD_COLOR,
        ),
    ] {
        chart
            .draw_series(std::iter::empty::<PathElement<(f64, f64)>>())?
            .label(label)
            .legend(move |(left, top)| {
                PathElement::new([(left, top), (left + 30, top)], color.stroke_width(3))
            });
    }
    chart
        .configure_series_labels()
        .label_font(("sans-serif", 22))
        .position(SeriesLabelPosition::UpperRight)
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .draw()?;
    Ok(())
}

/// Draw two related model quantities on their own coverage scale.
///
/// Chromosome names are centered over their cores and alternate between two rows to leave
/// room for neighboring short chromosomes. Both tracks retain the same core order and boundaries.
/// Legends sit below the axes so they cannot obscure data or chromosome labels.
fn draw_local_model_panel(
    area: &DrawingArea<BitMapBackend<'_>, plotters::coord::Shift>,
    diagnostics: &LocalModelDiagnostics,
    title: &str,
    series: [(&str, &[Option<f64>], RGBColor); 2],
) -> Result<()> {
    let (chart_area, legend_area) = area.split_vertically(area.dim_in_pixel().1.saturating_sub(45));
    let (legend_width, legend_height) = legend_area.dim_in_pixel();
    let legend_column_width = legend_width / series.len() as u32;
    let legend_text_style = ("sans-serif", 22)
        .into_font()
        .color(&BLACK)
        .pos(Pos::new(HPos::Left, VPos::Center));
    let number_of_cores = diagnostics.initial_thresholds.len();
    let maximum_value = series
        .iter()
        .flat_map(|(_, values, _)| values.iter().flatten().copied())
        .fold(0.0_f64, f64::max)
        .max(f64::EPSILON);
    // Leave two rows above the data for chromosome labels, independent of the coverage scale
    let y_axis_end = maximum_value * 1.35;
    let mut chart = ChartBuilder::on(&chart_area)
        .caption(title, ("sans-serif", 32))
        .margin(20)
        .x_label_area_size(75)
        .y_label_area_size(100)
        .build_cartesian_2d(0.0_f64..number_of_cores as f64, 0.0_f64..y_axis_end)?;
    chart
        .configure_mesh()
        .label_style(("sans-serif", 22))
        .axis_desc_style(("sans-serif", 24))
        .max_light_lines(0)
        .disable_x_mesh()
        .bold_line_style(BLACK.mix(0.12))
        .x_desc("Fixed-size model cores in chromosome order")
        .y_desc("Raw positional fragment coverage")
        .draw()?;

    for &boundary in &diagnostics.chromosome_boundaries {
        chart.draw_series(std::iter::once(PathElement::new(
            [(boundary as f64, 0.0), (boundary as f64, y_axis_end)],
            BLACK.mix(0.15).stroke_width(1),
        )))?;
    }

    for (label_index, (start, chromosome)) in diagnostics.chromosome_labels.iter().enumerate() {
        let end = diagnostics
            .chromosome_labels
            .get(label_index + 1)
            .map(|(next_start, _)| *next_start)
            .unwrap_or(number_of_cores);
        if *start == end {
            continue;
        }
        let center = (*start as f64 + end as f64) / 2.0;
        let label_height = if label_index % 2 == 0 { 1.26 } else { 1.12 };
        chart.draw_series(std::iter::once(Text::new(
            chromosome.clone(),
            (center, maximum_value * label_height),
            ("sans-serif", 22)
                .into_font()
                .color(&BLACK)
                .pos(Pos::new(HPos::Center, VPos::Center)),
        )))?;
    }
    for (series_index, (label, values, color)) in series.into_iter().enumerate() {
        // A larger first marker remains visible when both fitted quantities coincide
        let radius = if series_index == 0 { 3 } else { 2 };
        chart.draw_series(values.iter().enumerate().filter_map(|(index, &value)| {
            value.map(|value| Circle::new((index as f64 + 0.5, value), radius, color.filled()))
        }))?;

        let label_width = legend_area.estimate_text_size(label, &legend_text_style)?.0 as i32;
        let column_start = series_index as i32 * legend_column_width as i32;
        let legend_left = column_start + (legend_column_width as i32 - label_width - 40) / 2;
        let legend_top = legend_height as i32 / 2;
        legend_area.draw(&PathElement::new(
            [(legend_left, legend_top), (legend_left + 30, legend_top)],
            color.stroke_width(2),
        ))?;
        legend_area.draw(&Text::new(
            label,
            (legend_left + 40, legend_top),
            legend_text_style.clone(),
        ))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("plotting_tests.rs");
}
