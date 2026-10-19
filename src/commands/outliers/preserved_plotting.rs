//! Coverage summaries and optional PNG views around preserved reporting groups.

use super::{outliers::CalledRegion, preserved_regions::PreservedRegion, striding::CoreModel};
use crate::shared::interval::Interval;
use anyhow::{Context, Result, ensure};
use plotters::{coord::Shift, prelude::*};
use std::path::Path;

const MAXIMUM_PLOT_BINS: u32 = 2_000;

/// Coverage over one reference bin. Blacklisted positions do not contribute to any metric.
#[derive(Debug)]
struct CoverageBin {
    interval: Interval<u32>,
    coverage_sum: f64,
    eligible_positions: u64,
    maximum: f64,
    /// Position of the first observed peak, so mask shading cannot hide a relocated bin peak.
    maximum_position: u32,
}

impl CoverageBin {
    /// An empty bin is unavailable, whereas an observed all-zero bin has mean zero.
    fn mean(&self) -> Option<f64> {
        (self.eligible_positions > 0).then(|| self.coverage_sum / self.eligible_positions as f64)
    }
}

/// Fixed genomic bins spanning a plot context, accumulated across disjoint processing chunks.
pub(super) struct CoverageBins {
    interval: Interval<u32>,
    bin_size: u32,
    bins: Vec<CoverageBin>,
}

impl CoverageBins {
    /// Choose enough reference bases per bin to bound figure complexity without discarding peaks.
    pub(super) fn for_plot(interval: Interval<u32>) -> Result<Self> {
        Self::new(interval, interval.len().div_ceil(MAXIMUM_PLOT_BINS))
    }

    fn new(interval: Interval<u32>, bin_size: u32) -> Result<Self> {
        ensure!(bin_size > 0, "plot bin size must be greater than zero");
        let mut bins = Vec::new();
        let mut start = interval.start();
        while start < interval.end() {
            let end = start.saturating_add(bin_size).min(interval.end());
            bins.push(CoverageBin {
                interval: Interval::new(start, end)?,
                coverage_sum: 0.0,
                eligible_positions: 0,
                maximum: 0.0,
                maximum_position: start,
            });
            start = end;
        }
        Ok(Self {
            interval,
            bin_size,
            bins,
        })
    }

    /// Add each unblacklisted base once, including genuine zeros, independent of chunk boundaries.
    pub(super) fn add_coverage(
        &mut self,
        core: Interval<u32>,
        coverage: &[f32],
        mask: Option<&[u8]>,
    ) -> Result<()> {
        ensure!(
            self.interval.contains_interval(core) && core.len() as usize == coverage.len(),
            "plot coverage does not match its requested genomic interval"
        );
        ensure!(
            mask.is_none_or(|values| values.len() == coverage.len()),
            "plot coverage and blacklist mask have different lengths"
        );
        for (offset, &value) in coverage.iter().enumerate() {
            if mask.is_some_and(|values| values[offset] != 0) {
                continue;
            }
            ensure!(
                value.is_finite() && value >= 0.0,
                "invalid raw coverage {value} while plotting preserved region"
            );
            let position = core.start() + offset as u32;
            let bin_index = (position - self.interval.start()) / self.bin_size;
            let bin = &mut self.bins[bin_index as usize];
            bin.coverage_sum += f64::from(value);
            if bin.eligible_positions == 0 || f64::from(value) > bin.maximum {
                bin.maximum = f64::from(value);
                bin.maximum_position = position;
            }
            bin.eligible_positions += 1;
        }
        Ok(())
    }
}

/// Shared evidence for the context and close-up panels of a preserved-region figure.
pub(super) struct PreservedPlotData<'a> {
    pub(super) region: &'a PreservedRegion,
    pub(super) coverage: &'a CoverageBins,
    pub(super) models: &'a [CoreModel],
    pub(super) blacklist: &'a [Interval<u64>],
    pub(super) calls: &'a [CalledRegion],
}

/// Save an overview and a closer view with observed coverage, fitted baseline, cutoff, and masks.
///
/// The blue track is the mean among unblacklisted bases in each reference bin. Grey peak marks
/// retain narrow elevations lost by averaging. Masked spans are drawn as gaps, never as zeros.
pub(super) fn write_preserved_region_plot(path: &Path, data: PreservedPlotData<'_>) -> Result<()> {
    let root = BitMapBackend::new(path, (1600, 900)).into_drawing_area();
    root.fill(&WHITE)?;
    let titled = root.titled(
        &format!(
            "{} | {}:{}-{} | {} calls, {} called bases",
            data.region.region_id,
            data.region.chromosome,
            data.region.interval.start(),
            data.region.interval.end(),
            data.region.call_indices.len(),
            data.region.called_bases
        ),
        ("sans-serif", 21),
    )?;
    let panels = titled.split_evenly((2, 1));
    draw_coverage_panel(
        &panels[0],
        data.region.plot_interval,
        "Surrounding coverage",
        &data,
    )?;
    let padding = (data.region.interval.len() / 10).max(100);
    let zoom = Interval::new(
        data.region.interval.start().saturating_sub(padding),
        data.region.interval.end().saturating_add(padding),
    )?
    .intersection(data.region.plot_interval)
    .context("preserved plot close-up is outside its context")?;
    draw_coverage_panel(&panels[1], zoom, "Closer view of the grouped calls", &data)?;
    root.present()
        .with_context(|| format!("writing preserved-region plot {}", path.display()))?;
    Ok(())
}

/// Draw separate bin segments so unavailable coverage never connects across a masked gap.
fn draw_coverage_panel(
    area: &DrawingArea<BitMapBackend<'_>, Shift>,
    interval: Interval<u32>,
    title: &str,
    data: &PreservedPlotData<'_>,
) -> Result<()> {
    let bins = data
        .coverage
        .bins
        .iter()
        .filter(|bin| bin.interval.intersects(interval))
        .collect::<Vec<_>>();
    let models = data
        .models
        .iter()
        .filter(|model| model.core.intersects(interval))
        .collect::<Vec<_>>();
    let maximum = bins
        .iter()
        .filter(|bin| bin.eligible_positions > 0)
        .map(|bin| {
            let mean = bin.mean().unwrap_or(0.0);
            if interval.contains_point(bin.maximum_position) {
                bin.maximum
            } else {
                mean
            }
        })
        .chain(models.iter().flat_map(|model| {
            [
                f64::from(model.model.final_threshold),
                model.model.underlying_fit.component_mean,
            ]
        }))
        .fold(1.0, f64::max)
        * 1.1;
    let coordinate = |position: u32| f64::from(position) / 1000.0;
    let mut chart = ChartBuilder::on(area)
        .caption(
            format!("{title} ({} bp coverage bins)", data.coverage.bin_size),
            ("sans-serif", 17),
        )
        .margin(12)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(
            coordinate(interval.start())..coordinate(interval.end()),
            0.0..maximum,
        )?;
    let precision = if interval.len() < 10_000 {
        3
    } else if interval.len() < 100_000 {
        2
    } else {
        1
    };
    chart
        .configure_mesh()
        .x_labels(8)
        .x_label_formatter(&|position| format!("{position:.precision$}"))
        .x_desc(format!("Position on {} (kb)", data.region.chromosome))
        .y_desc("Fragment coverage")
        .draw()?;

    chart
        .draw_series(data.region.call_indices.iter().filter_map(|&index| {
            data.calls[index]
                .interval
                .intersection(interval)
                .map(|overlap| {
                    Rectangle::new(
                        [
                            (coordinate(overlap.start()), 0.0),
                            (coordinate(overlap.end()), maximum),
                        ],
                        RGBColor(240, 190, 60).mix(0.23).filled(),
                    )
                })
        }))?
        .label("Preserved calls")
        .legend(|(legend_x, legend_y)| {
            Rectangle::new(
                [(legend_x, legend_y - 4), (legend_x + 18, legend_y + 4)],
                RGBColor(240, 190, 60).mix(0.5).filled(),
            )
        });
    chart
        .draw_series(bins.iter().filter_map(|bin| {
            let mean = bin.mean()?;
            if !interval.contains_point(bin.maximum_position) {
                return None;
            }
            let position = coordinate(bin.maximum_position);
            Some(PathElement::new(
                [(position, mean), (position, bin.maximum)],
                RGBColor(160, 160, 160),
            ))
        }))?
        .label("Bin peak")
        .legend(|(legend_x, legend_y)| {
            PathElement::new(
                [(legend_x, legend_y), (legend_x + 18, legend_y)],
                RGBColor(160, 160, 160),
            )
        });
    chart
        .draw_series(bins.iter().filter_map(|bin| {
            let mean = bin.mean()?;
            let overlap = bin.interval.intersection(interval)?;
            Some(PathElement::new(
                [
                    (coordinate(overlap.start()), mean),
                    (coordinate(overlap.end()), mean),
                ],
                BLUE.stroke_width(2),
            ))
        }))?
        .label("Mean coverage")
        .legend(|(legend_x, legend_y)| {
            PathElement::new(
                [(legend_x, legend_y), (legend_x + 18, legend_y)],
                BLUE.stroke_width(2),
            )
        });
    chart
        .draw_series(models.iter().filter_map(|model| {
            let overlap = model.core.intersection(interval)?;
            let mean = model.model.underlying_fit.component_mean;
            Some(PathElement::new(
                [
                    (coordinate(overlap.start()), mean),
                    (coordinate(overlap.end()), mean),
                ],
                GREEN.stroke_width(2),
            ))
        }))?
        .label("Background mean")
        .legend(|(legend_x, legend_y)| {
            PathElement::new(
                [(legend_x, legend_y), (legend_x + 18, legend_y)],
                GREEN.stroke_width(2),
            )
        });
    chart
        .draw_series(models.iter().filter_map(|model| {
            let overlap = model.core.intersection(interval)?;
            let threshold = f64::from(model.model.final_threshold);
            Some(PathElement::new(
                [
                    (coordinate(overlap.start()), threshold),
                    (coordinate(overlap.end()), threshold),
                ],
                RED.stroke_width(2),
            ))
        }))?
        .label("Outlier cutoff")
        .legend(|(legend_x, legend_y)| {
            PathElement::new(
                [(legend_x, legend_y), (legend_x + 18, legend_y)],
                RED.stroke_width(2),
            )
        });
    // Paint masks last to hide both observed tracks and fitted lines at excluded positions
    let reference_interval = interval.try_to_u64()?;
    let masks = data
        .blacklist
        .iter()
        .filter_map(|mask| {
            let overlap = mask.intersection(reference_interval)?;
            Some(Rectangle::new(
                [
                    (overlap.start() as f64 / 1000.0, 0.0),
                    (overlap.end() as f64 / 1000.0, maximum),
                ],
                RGBColor(225, 225, 225).filled(),
            ))
        })
        .collect::<Vec<_>>();
    if !masks.is_empty() {
        chart
            .draw_series(masks)?
            .label("Blacklisted")
            .legend(|(legend_x, legend_y)| {
                Rectangle::new(
                    [(legend_x, legend_y - 4), (legend_x + 18, legend_y + 4)],
                    RGBColor(225, 225, 225).filled(),
                )
            });
    }
    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.9))
        .border_style(BLACK)
        .position(SeriesLabelPosition::UpperRight)
        .label_font(("sans-serif", 12))
        .draw()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    include!("preserved_plotting_tests.rs");
}
