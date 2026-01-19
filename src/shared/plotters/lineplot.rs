use anyhow::Result;
use plotters::{coord::Shift, prelude::*};
use std::path::Path;

/// Render a quick-look line plot to a PNG file for fast QC.
///
/// Draws labeled axes with a small padding around both ranges so the
/// series is not clipped, and lets callers control the canvas size for
/// embedding in reports. Fails fast when inputs are empty or have
/// mismatched lengths to keep caller errors clear.
///
/// Parameters
/// ----------
/// - `out_path`:
///     Destination path for the image file.
/// - `title`:
///     Plot title shown above the chart.
/// - `x_label`:
///     Label for the x axis.
/// - `y_label`:
///     Label for the y axis.
/// - `x_values`:
///     X coordinates for the series.
/// - `y_values`:
///     Y coordinates for the series. Must match `x_values` length.
/// - `width`:
///     Canvas width in pixels.
/// - `height`:
///     Canvas height in pixels.
///
/// Returns
/// -------
/// - `Result<()>`:
///     Ok when the plot is written.
pub fn write_line_plot_png<P: AsRef<Path>>(
    out_path: P,
    title: &str,
    x_label: &str,
    y_label: &str,
    x_values: &[f64],
    y_values: &[f64],
    width: u32,
    height: u32,
) -> Result<()> {
    let drawing_area = BitMapBackend::new(out_path.as_ref(), (width, height)).into_drawing_area();
    draw_line_plot(&drawing_area, title, x_label, y_label, x_values, y_values)
}

/// Render a quick-look line plot to an SVG file for fast QC.
///
/// Matches the PNG plot defaults while keeping the vector format for
/// downstream editing or embedding, and lets callers control the canvas
/// size.
///
/// Parameters
/// ----------
/// - `out_path`:
///     Destination path for the image file.
/// - `title`:
///     Plot title shown above the chart.
/// - `x_label`:
///     Label for the x axis.
/// - `y_label`:
///     Label for the y axis.
/// - `x_values`:
///     X coordinates for the series.
/// - `y_values`:
///     Y coordinates for the series. Must match `x_values` length.
/// - `width`:
///     Canvas width in pixels.
/// - `height`:
///     Canvas height in pixels.
///
/// Returns
/// -------
/// - `Result<()>`:
///     Ok when the plot is written.
pub fn write_line_plot_svg<P: AsRef<Path>>(
    out_path: P,
    title: &str,
    x_label: &str,
    y_label: &str,
    x_values: &[f64],
    y_values: &[f64],
    width: u32,
    height: u32,
) -> Result<()> {
    let drawing_area = SVGBackend::new(out_path.as_ref(), (width, height)).into_drawing_area();
    draw_line_plot(&drawing_area, title, x_label, y_label, x_values, y_values)
}

fn draw_line_plot<DB: DrawingBackend>(
    drawing_area: &DrawingArea<DB, Shift>,
    title: &str,
    x_label: &str,
    y_label: &str,
    x_values: &[f64],
    y_values: &[f64],
) -> Result<()>
where
    DB::ErrorType: 'static + std::error::Error + Send + Sync,
{
    // Fail fast on malformed data to give clear feedback to callers
    anyhow::ensure!(!x_values.is_empty(), "x axis values are empty");
    anyhow::ensure!(x_values.len() == y_values.len(), "x/y length mismatch");

    // Set range of `x` values

    let x_min = x_values.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = x_values.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    // Pad x range when all values are identical so Plotters has a non-zero domain
    let x_span = (x_max - x_min).abs();
    let x_pad = if x_span > 0.0 && x_span.is_finite() {
        x_span * 0.05
    } else {
        1.0
    };
    let x_range = (x_min - x_pad)..(x_max + x_pad);

    // Set range of `y` values

    let y_min = y_values.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = y_values.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let mut y_low = if y_min.is_finite() { y_min } else { 0.0 };
    let mut y_high = if y_max.is_finite() { y_max } else { 1.0 };

    // Ensure a non-zero span even if all values are identical
    if (y_high - y_low).abs() < f64::EPSILON {
        let pad = (y_low.abs() * 0.05).max(1.0);
        y_low -= pad;
        y_high += pad;
    }

    let y_span = (y_high - y_low).abs();
    let y_pad = if y_span > 0.0 && y_span.is_finite() {
        y_span * 0.05
    } else {
        1.0
    };
    let y_range = (y_low - y_pad)..(y_high + y_pad);

    // Keep the background clear and axis labels readable for quick QC
    drawing_area.fill(&WHITE)?;

    // Reserve a modest fixed strip for axis labels
    let x_label_area = 52;
    let y_label_area = 52;

    let mut chart = ChartBuilder::on(drawing_area)
        .caption(title, ("sans-serif", 22))
        .margin(20)
        .x_label_area_size(x_label_area)
        .y_label_area_size(y_label_area)
        .build_cartesian_2d(x_range, y_range)?;

    chart
        .configure_mesh()
        .axis_desc_style(("sans-serif", 22))
        .x_desc(x_label)
        .y_desc(y_label)
        .draw()?;

    chart.draw_series(LineSeries::new(
        x_values.iter().copied().zip(y_values.iter().copied()),
        &BLUE,
    ))?;

    drawing_area.present()?;
    Ok(())
}
