use anyhow::{Result, ensure};
use ndarray::Array2;
use plotters::{
    coord::Shift,
    prelude::*,
    style::text_anchor::{HPos, Pos, VPos},
};
use std::borrow::Cow;
use std::path::Path;

use super::histogram::HistogramSpec;

/// Output formats supported by the heatmap writer.
pub enum HeatmapFormat {
    Png,
    Svg,
}

/// Upsampling algorithm for heatmap rendering.
pub enum HeatmapUpsample {
    Nearest,
    Bilinear,
}

// Height reserved below the heatmap for the legend, lowering gives the heatmap more vertical room
const LEGEND_HEIGHT: u32 = 40;
// Top margin around the heatmap plot, lowering pulls the plot closer to the top histogram
const HEATMAP_MARGIN_TOP: u32 = 4;
// Right margin around the heatmap plot, lowering pulls the plot closer to the right histogram
const HEATMAP_MARGIN_RIGHT: u32 = 4;
// Left margin around the heatmap plot, keep unchanged to preserve axis label spacing
const HEATMAP_MARGIN_LEFT: u32 = 10;
// Bottom margin around the heatmap plot, keep unchanged to preserve x-axis spacing
const HEATMAP_MARGIN_BOTTOM: u32 = 10;
// Vertical space for x-axis labels and ticks on the heatmap, reducing pulls the plot toward the legend
const HEATMAP_X_LABEL_AREA: u32 = 52;
// Horizontal space for y-axis labels and ticks on the heatmap, reducing pulls the plot toward the left edge
const HEATMAP_Y_LABEL_AREA: u32 = 62;
// Height reserved for the title when a top histogram exists, lowering moves the histogram upward
const TITLE_HEIGHT_WITH_TOP_HIST: u32 = 50;
// Vertical draw height for the top histogram bars, lowering makes the bars shorter without affecting the gap
const TOP_HIST_HEIGHT: u32 = 70;
// Padding between the bottom of the top histogram and the start of the heatmap, set to zero for no extra gap
const TOP_HIST_GAP_BELOW: u32 = 0;
// Minimum height guaranteed for the heatmap after carving the top histogram and its gap, raising forces the top panel to shrink first
const MIN_HEATMAP_HEIGHT_AFTER_TOP: u32 = 140;
// Desired width for the right histogram panel, lowering gives more width to the heatmap
const RIGHT_PANEL_TARGET_WIDTH: u32 = 70;
// Minimum width the heatmap keeps when allocating the right panel, raising shrinks the right panel when space is tight
const MIN_HEATMAP_WIDTH_AFTER_RIGHT: u32 = 200;
// Padding between the heatmap and the right histogram panel, set to zero for no gap
const RIGHT_PANEL_GAP: u32 = 0;
// Margin above histogram bars inside their panels, set to zero to let bars touch the panel top
const HIST_MARGIN_TOP: u32 = 0;
// Margin below histogram bars inside their panels, set to zero to let bars touch the panel bottom
const HIST_MARGIN_BOTTOM: u32 = 0;
// Margin to the left of histogram bars inside their panels, set to zero to let bars touch the panel left edge
const HIST_MARGIN_LEFT: u32 = 0;
// Margin to the right of histogram bars inside their panels, set to zero to let bars touch the panel right edge
const HIST_MARGIN_RIGHT: u32 = 0;
// Space for x-axis labels on histogram charts, set to zero to remove histogram x-axis
const HIST_X_LABEL_AREA: u32 = 0;
// Space for y-axis labels on histogram charts, set to zero to remove histogram y-axis
const HIST_Y_LABEL_AREA: u32 = 0;
// Enable debug backgrounds for histogram plot areas (top: green, right: red)
const DEBUG_HIST_BACKGROUNDS: bool = false;
// Fill color for histogram bars
const HIST_BAR_COLOR: RGBColor = RGBColor(161, 174, 177);

/// Render a heatmap from a matrix to an image.
///
/// Draws each cell as a filled rectangle spanning the provided axis edges,
/// normalizes colors between the finite minimum and maximum values (or the
/// supplied range), and skips non-finite entries so callers can mask
/// unsupported regions. Supports both PNG and SVG backends so you can choose
/// between quick raster snapshots and vector output for reports.
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
/// - `values`:
///     Matrix to render where rows align with `y_edges` and columns with `x_edges`.
/// - `x_edges`:
///     Optional x-axis boundaries for each column.
///     Used to map matrix column indices onto real-world units (e.g., GC bins -> GC percent edges).
///     Length must be `values.ncols() + 1` when provided.
/// - `y_edges`:
///     Optional y-axis boundaries for each row.
///     Used to map matrix row indices onto real-world units (e.g., fragment length bins -> fragment length edges).
///     Length must be `values.nrows() + 1` when provided.
/// - `val_min`:
///     Optional lower bound for color scaling. Defaults to the minimum finite value.
/// - `val_max`:
///     Optional upper bound for color scaling. Defaults to the maximum finite value.
/// - `val_center`:
///     Optional center value for a diverging scale. Values below use a cool gradient toward
///     the center and values above use a warm gradient.
/// - `min_color`:
///     Optional color for the minimum value. Defaults depend on the palette: pink when
///     using a diverging center, black otherwise.
/// - `max_color`:
///     Optional color for the maximum value. Defaults to green.
/// - `symmetric_diverging`:
///     When true, uses the maximum absolute distance from the center to scale both
///     sides of a diverging palette so gradients are symmetric.
/// - `upsample_factor`:
///     Bilinear upsampling factor applied to the matrix before plotting to reduce visible blockiness. Use 1 to disable.
/// - `upsample_method`:
///     Algorithm used when upsampling: nearest (crisp blocks) or bilinear (smooth).
/// - `width`:
///     Canvas width in pixels.
/// - `height`:
///     Canvas height in pixels.
/// - `format`:
///     Output format, choose PNG for quick looks or SVG for vector editing.
///
/// Returns
/// -------
/// - `Result<()>`:
///     Ok when the plot is written.
pub fn write_heatmap<P: AsRef<Path>>(
    out_path: P,
    title: &str,
    x_label: &str,
    y_label: &str,
    values: &Array2<f64>,
    x_edges: Option<&[f64]>,
    y_edges: Option<&[f64]>,
    val_min: Option<f64>,
    val_max: Option<f64>,
    val_center: Option<f64>,
    min_color: Option<RGBColor>,
    max_color: Option<RGBColor>,
    symmetric_diverging: bool,
    upsample_factor: usize,
    upsample_method: HeatmapUpsample,
    width: u32,
    height: u32,
    format: HeatmapFormat,
) -> Result<()> {
    let (x_edges, y_edges, values_to_plot, mean_val, min_val, max_val) = prepare_heatmap_inputs(
        values,
        x_edges,
        y_edges,
        val_min,
        val_max,
        val_center,
        upsample_factor,
        upsample_method,
    )?;

    let out_path = out_path.as_ref();
    match format {
        HeatmapFormat::Png => {
            let drawing_area = BitMapBackend::new(out_path, (width, height)).into_drawing_area();
            draw_heatmap(
                &drawing_area,
                Some(title),
                x_label,
                y_label,
                &x_edges,
                &y_edges,
                &values_to_plot,
                min_val,
                max_val,
                mean_val,
                val_center,
                min_color,
                max_color,
                symmetric_diverging,
            )
        }
        HeatmapFormat::Svg => {
            let drawing_area = SVGBackend::new(out_path, (width, height)).into_drawing_area();
            draw_heatmap(
                &drawing_area,
                Some(title),
                x_label,
                y_label,
                &x_edges,
                &y_edges,
                &values_to_plot,
                min_val,
                max_val,
                mean_val,
                val_center,
                min_color,
                max_color,
                symmetric_diverging,
            )
        }
    }
}

/// Render a heatmap with optional top and right histograms for marginal mass.
///
/// Splits the canvas into panels for the histograms and heatmap, keeping the
/// legend attached to the heatmap area. When a histogram is omitted the space
/// is returned to the heatmap so defaults remain unchanged.
///
/// Parameters
/// ----------
/// - `x_hist`:
///     Optional histogram to render above the heatmap using the x-axis scale.
/// - `y_hist`:
///     Optional histogram to render to the right of the heatmap using the y-axis scale.
pub fn write_heatmap_with_histograms<P: AsRef<Path>>(
    out_path: P,
    title: &str,
    x_label: &str,
    y_label: &str,
    values: &Array2<f64>,
    x_edges: Option<&[f64]>,
    y_edges: Option<&[f64]>,
    x_hist: Option<&HistogramSpec>,
    y_hist: Option<&HistogramSpec>,
    val_min: Option<f64>,
    val_max: Option<f64>,
    val_center: Option<f64>,
    min_color: Option<RGBColor>,
    max_color: Option<RGBColor>,
    symmetric_diverging: bool,
    upsample_factor: usize,
    upsample_method: HeatmapUpsample,
    width: u32,
    height: u32,
    format: HeatmapFormat,
) -> Result<()> {
    let (x_edges, y_edges, values_to_plot, mean_val, min_val, max_val) = prepare_heatmap_inputs(
        values,
        x_edges,
        y_edges,
        val_min,
        val_max,
        val_center,
        upsample_factor,
        upsample_method,
    )?;

    let out_path = out_path.as_ref();
    match format {
        HeatmapFormat::Png => draw_heatmap_with_layout(
            BitMapBackend::new(out_path, (width, height)).into_drawing_area(),
            title,
            x_label,
            y_label,
            &x_edges,
            &y_edges,
            &values_to_plot,
            min_val,
            max_val,
            mean_val,
            val_center,
            min_color,
            max_color,
            symmetric_diverging,
            x_hist,
            y_hist,
        ),
        HeatmapFormat::Svg => draw_heatmap_with_layout(
            SVGBackend::new(out_path, (width, height)).into_drawing_area(),
            title,
            x_label,
            y_label,
            &x_edges,
            &y_edges,
            &values_to_plot,
            min_val,
            max_val,
            mean_val,
            val_center,
            min_color,
            max_color,
            symmetric_diverging,
            x_hist,
            y_hist,
        ),
    }
}

/// Draw the heatmap and optional legend on the provided drawing area.
///
/// Builds axes, fills each cell using the chosen palette, and places the legend
/// beneath the plot when there is vertical room.
///
/// Parameters
/// ----------
/// - `drawing_area`:
///     Target drawing area from the backend.
/// - `title`:
///     Plot title.
/// - `x_label`:
///     Label for the x axis.
/// - `y_label`:
///     Label for the y axis.
/// - `x_edges`:
///     X boundaries for each column.
/// - `y_edges`:
///     Y boundaries for each row.
/// - `values`:
///     Matrix of values to render.
/// - `min_val`:
///     Lower bound for color scaling.
/// - `max_val`:
///     Upper bound for color scaling.
/// - `mean_val`:
///     Optional mean value shown in the legend.
/// - `center_val`:
///     Optional diverging center.
/// - `min_color`:
///     Optional color for the minimum value. Defaults depend on the palette: pink when
///     using a diverging center, pink otherwise.
/// - `max_color`:
///     Optional color for the maximum value. Defaults to green.
/// - `symmetric_diverging`:
///     When true, uses the maximum absolute distance from the center to scale both sides so the diverging gradients share a common curve.
///
/// Returns
/// -------
/// - `Result<()>`:
///     Ok when drawing finishes.
fn draw_heatmap<DB: DrawingBackend>(
    drawing_area: &DrawingArea<DB, Shift>,
    title: Option<&str>,
    x_label: &str,
    y_label: &str,
    x_edges: &[f64],
    y_edges: &[f64],
    values: &Array2<f64>,
    min_val: f64,
    max_val: f64,
    mean_val: Option<f64>,
    center_val: Option<f64>,
    min_color: Option<RGBColor>,
    max_color: Option<RGBColor>,
    symmetric_diverging: bool,
) -> Result<()>
where
    DB::ErrorType: 'static + std::error::Error + Send + Sync,
{
    drawing_area.fill(&WHITE)?;

    let (_, area_h) = drawing_area.dim_in_pixel();
    let (plot_area, legend_area) = if area_h > LEGEND_HEIGHT {
        let (upper, lower) = drawing_area.split_vertically(area_h - LEGEND_HEIGHT);
        (upper, Some(lower))
    } else {
        (drawing_area.clone(), None)
    };

    let x_range = *x_edges.first().unwrap()..*x_edges.last().unwrap();
    let y_range = *y_edges.first().unwrap()..*y_edges.last().unwrap();

    let mut base_builder = ChartBuilder::on(&plot_area);
    let mut builder = base_builder
        .margin_top(HEATMAP_MARGIN_TOP)
        .margin_right(HEATMAP_MARGIN_RIGHT)
        .margin_bottom(HEATMAP_MARGIN_BOTTOM)
        .margin_left(HEATMAP_MARGIN_LEFT)
        .x_label_area_size(HEATMAP_X_LABEL_AREA)
        .y_label_area_size(HEATMAP_Y_LABEL_AREA);
    if let Some(t) = title {
        builder = builder.caption(t, ("sans-serif", 22));
    }
    let mut chart = builder.build_cartesian_2d(x_range, y_range)?;

    chart
        .configure_mesh()
        .axis_desc_style(("sans-serif", 22))
        .x_desc(x_label)
        .y_desc(y_label)
        .disable_mesh()
        .draw()?;

    for (row_idx, row) in values.outer_iter().enumerate() {
        let y0 = y_edges[row_idx];
        let y1 = y_edges[row_idx + 1];
        for (col_idx, &value) in row.iter().enumerate() {
            if !value.is_finite() {
                continue;
            }
            let x0 = x_edges[col_idx];
            let x1 = x_edges[col_idx + 1];
            let color = color_for_value(
                value,
                min_val,
                max_val,
                center_val,
                min_color,
                max_color,
                symmetric_diverging,
            );
            chart.draw_series(std::iter::once(Rectangle::new(
                [(x0, y0), (x1, y1)],
                color.filled(),
            )))?;
        }
    }

    if let Some(legend_area) = legend_area {
        draw_color_legend(
            &legend_area,
            min_val,
            max_val,
            mean_val,
            center_val,
            min_color,
            max_color,
            symmetric_diverging,
        )?;
    }

    plot_area.present()?;
    Ok(())
}

fn draw_heatmap_with_layout<DB: DrawingBackend>(
    root_area: DrawingArea<DB, Shift>,
    title: &str,
    x_label: &str,
    y_label: &str,
    x_edges: &[f64],
    y_edges: &[f64],
    values_to_plot: &Array2<f64>,
    min_val: f64,
    max_val: f64,
    mean_val: Option<f64>,
    val_center: Option<f64>,
    min_color: Option<RGBColor>,
    max_color: Option<RGBColor>,
    symmetric_diverging: bool,
    x_hist: Option<&HistogramSpec>,
    y_hist: Option<&HistogramSpec>,
) -> Result<()>
where
    DB::ErrorType: 'static + std::error::Error + Send + Sync,
{
    root_area.fill(&WHITE)?;

    // First, reserve space for a title (when a top histogram is present)
    let mut work_area = root_area.clone();
    let mut title_area = None;
    if x_hist.is_some() {
        let title_height = TITLE_HEIGHT_WITH_TOP_HIST;
        let (_, h) = work_area.dim_in_pixel();
        if h > title_height {
            let split = work_area.split_vertically(title_height);
            title_area = Some(split.0);
            work_area = split.1;
        }
    }

    // Reserve space for the right histogram so the top panel width matches the heatmap
    let mut main_area = work_area.clone();
    let mut right_area = None;
    if y_hist.is_some() {
        let available_w = main_area.dim_in_pixel().0;
        let desired = RIGHT_PANEL_TARGET_WIDTH;
        let total_right = desired + RIGHT_PANEL_GAP;
        let right_width =
            total_right.min(available_w.saturating_sub(MIN_HEATMAP_WIDTH_AFTER_RIGHT));
        if right_width > 0 {
            // First split off the combined gap+panel, then drop the gap so only the panel remains
            let split = main_area.split_horizontally(available_w.saturating_sub(right_width));
            main_area = split.0;
            let mut panel_block = split.1;
            if RIGHT_PANEL_GAP > 0 {
                let (_, panel_only) = panel_block.split_horizontally(RIGHT_PANEL_GAP);
                panel_block = panel_only;
            }
            right_area = Some(panel_block);
        }
    }

    // Then carve off the top histogram from the main area so widths stay aligned
    let mut heatmap_area = main_area.clone();
    let mut top_area = None;
    let mut top_height = 0;
    if x_hist.is_some() {
        let (_, main_h) = main_area.dim_in_pixel();
        let total_needed = TOP_HIST_HEIGHT + TOP_HIST_GAP_BELOW;
        top_height = total_needed.min(main_h.saturating_sub(MIN_HEATMAP_HEIGHT_AFTER_TOP));
        if top_height > 0 {
            let split = heatmap_area.split_vertically(top_height);
            let panel = split.0;
            heatmap_area = split.1;

            // Split the panel into the histogram draw area and a gap below it
            let panel_h = panel.dim_in_pixel().1;
            let hist_h = TOP_HIST_HEIGHT.min(panel_h);
            let (hist_area, _gap) = panel.split_vertically(hist_h);
            top_area = Some(hist_area);
        }
    }

    // Align the right histogram vertical span with the heatmap by removing the top histogram height
    if let Some(area) = right_area.take() {
        if top_height > 0 {
            let (_, area_h) = area.dim_in_pixel();
            let effective_top = top_height.min(area_h);
            let (_, lower) = area.split_vertically(effective_top);
            right_area = Some(lower);
        } else {
            right_area = Some(area);
        }
    }

    // Draw title above the top histogram when present
    if let Some(area) = title_area {
        let (w, h) = area.dim_in_pixel();
        let text_style = ("sans-serif", 22)
            .into_text_style(&area)
            .pos(Pos::new(HPos::Center, VPos::Center))
            .color(&BLACK);
        area.draw(&Text::new(
            title.to_string(),
            (w as i32 / 2, h as i32 / 2),
            text_style,
        ))?;
    }

    draw_heatmap(
        &heatmap_area,
        if x_hist.is_some() { None } else { Some(title) },
        x_label,
        y_label,
        x_edges,
        y_edges,
        values_to_plot,
        min_val,
        max_val,
        mean_val,
        val_center,
        min_color,
        max_color,
        symmetric_diverging,
    )?;

    let x_limits = (*x_edges.first().unwrap(), *x_edges.last().unwrap());
    let y_limits = (*y_edges.first().unwrap(), *y_edges.last().unwrap());

    if let Some(area) = top_area {
        if let Some(hist) = x_hist {
            draw_histogram_top(&area, hist, x_limits)?;
        }
    }
    if let Some(area) = right_area {
        if let Some(hist) = y_hist {
            // Match right histogram height to the heatmap plot box (exclude legend and x-axis area) and tint for debugging
            let heatmap_plot_h = heatmap_area.dim_in_pixel().1.saturating_sub(
                LEGEND_HEIGHT + HEATMAP_MARGIN_TOP + HEATMAP_MARGIN_BOTTOM + HEATMAP_X_LABEL_AREA,
            );
            let (_, area_h) = area.dim_in_pixel();
            let top_pad = HEATMAP_MARGIN_TOP;
            let target_h = heatmap_plot_h.min(area_h.saturating_sub(top_pad));
            // Carve off a top padding band (kept white) and draw the sidebar just below it without shortening height
            let (pad_area, rest) = area.split_vertically(top_pad);
            pad_area.fill(&WHITE)?;
            let (aligned_area, _) = rest.split_vertically(target_h);
            if DEBUG_HIST_BACKGROUNDS {
                aligned_area.fill(&RED)?;
            }
            draw_histogram_right(&aligned_area, hist, y_limits)?;
        }
    }

    root_area.present()?;
    Ok(())
}

fn prepare_heatmap_inputs<'a>(
    values: &'a Array2<f64>,
    x_edges: Option<&[f64]>,
    y_edges: Option<&[f64]>,
    val_min: Option<f64>,
    val_max: Option<f64>,
    val_center: Option<f64>,
    upsample_factor: usize,
    upsample_method: HeatmapUpsample,
) -> Result<(
    Vec<f64>,
    Vec<f64>,
    Cow<'a, Array2<f64>>,
    Option<f64>,
    f64,
    f64,
)> {
    let upsample_factor = upsample_factor.max(1);
    let mut x_edges = resolve_edges("x", x_edges, values.ncols())?;
    let mut y_edges = resolve_edges("y", y_edges, values.nrows())?;
    let mut values_to_plot: Cow<'a, Array2<f64>> = Cow::Borrowed(values);
    if upsample_factor > 1 {
        values_to_plot = Cow::Owned(match upsample_method {
            HeatmapUpsample::Nearest => upsample_nearest(values, upsample_factor),
            HeatmapUpsample::Bilinear => upsample_bilinear(values, upsample_factor),
        });
        x_edges = subdivide_edges(&x_edges, upsample_factor)?;
        y_edges = subdivide_edges(&y_edges, upsample_factor)?;
    }
    let mean_val = find_finite_mean(&values_to_plot);

    let (data_min, data_max) = find_finite_min_max(&values_to_plot);
    ensure!(
        data_min.is_finite() && data_max.is_finite(),
        "heatmap values are empty"
    );

    let min_val = val_min.unwrap_or(data_min);
    let max_val = val_max.unwrap_or(data_max);
    ensure!(
        min_val.is_finite() && max_val.is_finite(),
        "heatmap limits must be finite"
    );
    ensure!(max_val > min_val, "heatmap max must be greater than min");
    if let Some(center) = val_center {
        ensure!(
            center > min_val && center < max_val,
            "heatmap center must be within (min, max)"
        );
    }

    Ok((x_edges, y_edges, values_to_plot, mean_val, min_val, max_val))
}

fn draw_histogram_top<DB: DrawingBackend>(
    area: &DrawingArea<DB, Shift>,
    hist: &HistogramSpec,
    x_range: (f64, f64),
) -> Result<()>
where
    DB::ErrorType: 'static + std::error::Error + Send + Sync,
{
    // Keep panel white and only tint the plotting region for debugging alignment
    area.fill(&WHITE)?;
    let (x_min, x_max) = x_range;
    let max_y = hist.max().max(1.0);
    let mut chart = ChartBuilder::on(area)
        // Align horizontally with the heatmap by matching its left/right insets
        .margin_top(HIST_MARGIN_TOP)
        .margin_bottom(HIST_MARGIN_BOTTOM)
        .margin_left(HEATMAP_MARGIN_LEFT + HEATMAP_Y_LABEL_AREA)
        .margin_right(HEATMAP_MARGIN_RIGHT)
        // No axes for the histogram, so keep label areas at zero
        .x_label_area_size(0)
        .y_label_area_size(0)
        .build_cartesian_2d(x_min..x_max, 0.0..max_y)?;

    if DEBUG_HIST_BACKGROUNDS {
        chart.plotting_area().fill(&GREEN)?;
    }

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_labels(0)
        .y_labels(0)
        .axis_style(&WHITE)
        .draw()?;

    let bar_style = ShapeStyle {
        color: HIST_BAR_COLOR.to_rgba(),
        filled: true,
        stroke_width: 0,
    };
    for (idx, &count) in hist.counts.iter().enumerate() {
        let left = hist.edges[idx];
        let right = hist.edges[idx + 1];
        // Clamp bars so the histogram stops exactly at the heatmap bounds
        let clamped_left = left.max(x_min);
        let clamped_right = right.min(x_max);
        if clamped_right <= clamped_left {
            continue;
        }
        chart.draw_series(std::iter::once(Rectangle::new(
            [(clamped_left, 0.0), (clamped_right, count)],
            bar_style,
        )))?;
    }

    Ok(())
}

fn draw_histogram_right<DB: DrawingBackend>(
    area: &DrawingArea<DB, Shift>,
    hist: &HistogramSpec,
    y_range: (f64, f64),
) -> Result<()>
where
    DB::ErrorType: 'static + std::error::Error + Send + Sync,
{
    let (y_min, y_max) = y_range;
    let max_x = hist.max().max(1.0);
    let mut chart = ChartBuilder::on(area)
        .margin(0)
        .margin_top(HIST_MARGIN_TOP)
        .margin_bottom(HIST_MARGIN_BOTTOM)
        .margin_left(HIST_MARGIN_LEFT)
        .margin_right(HIST_MARGIN_RIGHT)
        .x_label_area_size(HIST_X_LABEL_AREA)
        .y_label_area_size(HIST_Y_LABEL_AREA)
        .build_cartesian_2d(0.0..max_x, y_min..y_max)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_labels(0)
        .y_labels(0)
        .axis_style(&WHITE)
        .draw()?;

    let bar_style = ShapeStyle {
        color: HIST_BAR_COLOR.to_rgba(),
        filled: true,
        stroke_width: 0,
    };
    for (idx, &count) in hist.counts.iter().enumerate() {
        let bottom = hist.edges[idx];
        let top = hist.edges[idx + 1];
        // Clamp bars so the histogram ends flush with the heatmap range
        let clamped_bottom = bottom.max(y_min);
        let clamped_top = top.min(y_max);
        if clamped_top <= clamped_bottom {
            continue;
        }
        chart.draw_series(std::iter::once(Rectangle::new(
            [(0.0, clamped_bottom), (count, clamped_top)],
            bar_style,
        )))?;
    }

    Ok(())
}

/// Find finite min and max in a matrix.
///
/// Skips non-finite entries so masked regions do not affect limits.
///
/// Parameters
/// ----------
/// - `values`:
///     Matrix to scan.
///
/// Returns
/// -------
/// - `(f64, f64)`:
///     Finite minimum and maximum.
fn find_finite_min_max(values: &Array2<f64>) -> (f64, f64) {
    let mut min_val = f64::INFINITY;
    let mut max_val = f64::NEG_INFINITY;
    for &v in values.iter() {
        if v.is_finite() {
            if v < min_val {
                min_val = v;
            }
            if v > max_val {
                max_val = v;
            }
        }
    }
    (min_val, max_val)
}

/// Compute the mean of finite entries in a matrix.
///
/// Ignores non-finite values and returns `None` when no finite entries exist.
///
/// Parameters
/// ----------
/// - `values`:
///     Matrix to average.
///
/// Returns
/// -------
/// - `Option<f64>`:
///     Mean of finite entries, or None when empty.
fn find_finite_mean(values: &Array2<f64>) -> Option<f64> {
    let mut sum = 0.0;
    let mut count = 0usize;
    for &v in values.iter() {
        if v.is_finite() {
            sum += v;
            count += 1;
        }
    }
    if count > 0 {
        Some(sum / count as f64)
    } else {
        None
    }
}

/// Resolve axis edges, defaulting to contiguous indices when missing.
///
/// Validates caller-supplied edges to ensure they match the expected length.
///
/// Parameters
/// ----------
/// - `name`:
///     Axis label used in error messages.
/// - `edges`:
///     Optional user-provided edge vector.
/// - `len`:
///     Number of bins along the axis.
///
/// Returns
/// -------
/// - `Vec<f64>`:
///     Validated edge vector.
fn resolve_edges(name: &str, edges: Option<&[f64]>, len: usize) -> Result<Vec<f64>> {
    if let Some(edges) = edges {
        ensure!(
            edges.len() == len + 1,
            "{}_edges length must be len + 1 (len={})",
            name,
            len
        );
        return Ok(edges.to_vec());
    }
    Ok((0..=len).map(|i| i as f64).collect())
}

/// Upsample by nearest-neighbor (pixel replication).
///
/// Simply repeats each cell `factor` times along both axes to preserve crisp
/// boundaries.
fn upsample_nearest(values: &Array2<f64>, factor: usize) -> Array2<f64> {
    let factor = factor.max(1);
    if factor == 1 {
        return values.clone();
    }
    let (rows, cols) = values.dim();
    let mut out = Array2::<f64>::zeros((rows * factor, cols * factor));
    for r in 0..rows {
        for c in 0..cols {
            let v = values[(r, c)];
            let row_start = r * factor;
            let col_start = c * factor;
            for rr in row_start..row_start + factor {
                for cc in col_start..col_start + factor {
                    out[(rr, cc)] = v;
                }
            }
        }
    }
    out
}

/// Bilinearly upsample a matrix by an integer factor.
///
/// Expands each cell smoothly so higher output resolutions do not appear blocky.
/// Uses standard bilinear interpolation: blend along x within the nearest two
/// source rows, then blend those row results along y.
///
/// Parameters
/// ----------
/// - `values`:
///     Input matrix.
/// - `factor`:
///     Integer upsampling factor. Values below 1 are treated as 1.
///
/// Returns
/// -------
/// - `Array2<f64>`:
///     Upsampled matrix.
fn upsample_bilinear(values: &Array2<f64>, factor: usize) -> Array2<f64> {
    let factor = factor.max(1);
    if factor == 1 {
        return values.clone();
    }
    let (rows, cols) = values.dim();
    if rows == 0 || cols == 0 {
        return values.clone();
    }
    let new_rows = rows * factor;
    let new_cols = cols * factor;
    let mut out = Array2::<f64>::zeros((new_rows, new_cols));

    for r in 0..new_rows {
        // Map the output row back to the source grid in floating point
        let src_y = (r as f64) / factor as f64;
        // Nearest source rows above and below
        let y0 = src_y.floor().max(0.0) as usize;
        let y1 = (y0 + 1).min(rows - 1);
        // Fractional position between y0 and y1
        let ty = (src_y - y0 as f64).clamp(0.0, 1.0);

        for c in 0..new_cols {
            // Map the output column back to the source grid in floating point
            let src_x = (c as f64) / factor as f64;
            // Nearest source columns left and right
            let x0 = src_x.floor().max(0.0) as usize;
            let x1 = (x0 + 1).min(cols - 1);
            // Fractional position between x0 and x1
            let tx = (src_x - x0 as f64).clamp(0.0, 1.0);

            // Source cell values at the four surrounding corners
            let v00 = values[(y0, x0)];
            let v01 = values[(y0, x1)];
            let v10 = values[(y1, x0)];
            let v11 = values[(y1, x1)];

            // Interpolate horizontally on the top and bottom edges
            let v0 = v00 * (1.0 - tx) + v01 * tx;
            let v1 = v10 * (1.0 - tx) + v11 * tx;

            // Interpolate vertically between the two edges
            out[(r, c)] = v0 * (1.0 - ty) + v1 * ty;
        }
    }
    out
}

/// Subdivide axis edges to align with an upsampled matrix.
///
/// Inserts evenly spaced intermediate edges within each original interval.
///
/// Parameters
/// ----------
/// - `edges`:
///     Original axis edge vector.
/// - `factor`:
///     Upsampling factor. Values below 1 are treated as 1.
///
/// Returns
/// -------
/// - `Vec<f64>`:
///     Refined edge vector.
fn subdivide_edges(edges: &[f64], factor: usize) -> Result<Vec<f64>> {
    ensure!(edges.len() >= 2, "edges must contain at least two points");
    let factor = factor.max(1);
    if factor == 1 {
        return Ok(edges.to_vec());
    }
    let mut out = Vec::with_capacity((edges.len() - 1) * factor + 1);
    for window in edges.windows(2) {
        let start = window[0];
        let end = window[1];
        let step = (end - start) / factor as f64;
        for i in 0..factor {
            out.push(start + step * i as f64);
        }
    }
    out.push(*edges.last().unwrap());
    Ok(out)
}

/// Map a value to a color with optional diverging center.
///
/// Uses a Blue-Yellow-Red diverging palette when a center is set, or a single Yellow-Red
/// gradient otherwise.
///
/// Parameters
/// ----------
/// - `value`:
///     Value to map.
/// - `min_val`:
///     Lower bound for scaling.
/// - `max_val`:
///     Upper bound for scaling.
/// - `center_val`:
///     Optional diverging center.
/// - `symmetric_diverging`:
///     When true, scales both sides using the maximum absolute distance from the center.
///
/// Returns
/// -------
/// - `RGBColor`:
///     Color for the value.
fn color_for_value(
    value: f64,
    min_val: f64,
    max_val: f64,
    center_val: Option<f64>,
    min_color: Option<RGBColor>,
    max_color: Option<RGBColor>,
    symmetric_diverging: bool,
) -> RGBColor {
    // Color palettes:
    // Diverging: pink: ff00f6 (255,0,246), black: 000000 (0,0,0), green: 0cff00 (12,255,0)
    // Single: black: 000000 (0,0,0), green: 0cff00 (12,255,0)

    let (default_min, default_max) = if center_val.is_some() {
        (RGBColor(255, 0, 246), RGBColor(12, 255, 0))
    } else {
        (RGBColor(0, 0, 0), RGBColor(12, 255, 0))
    };
    let center_color = RGBColor(0, 0, 0);
    let min_color = min_color.unwrap_or(default_min);
    let max_color = max_color.unwrap_or(default_max);

    if let Some(center) = center_val {
        if symmetric_diverging {
            let span = (max_val - center)
                .abs()
                .max((center - min_val).abs())
                .max(f64::EPSILON);
            let norm = ((value - center) / span).clamp(-1.0, 1.0);
            return if norm < 0.0 {
                interpolate_rgb(center_color, min_color, -norm)
            } else {
                interpolate_rgb(center_color, max_color, norm)
            };
        }

        if value <= center {
            let norm = ((value - min_val) / (center - min_val).max(f64::EPSILON)).clamp(0.0, 1.0);
            return interpolate_rgb(min_color, center_color, norm);
        } else {
            let norm = ((value - center) / (max_val - center).max(f64::EPSILON)).clamp(0.0, 1.0);
            return interpolate_rgb(center_color, max_color, norm);
        }
    }
    let norm = ((value - min_val) / (max_val - min_val).max(f64::EPSILON)).clamp(0.0, 1.0);
    interpolate_rgb(min_color, max_color, norm)
}

/// Linearly interpolate between two RGB colors.
///
/// Parameters
/// ----------
/// - `start`:
///     Start color.
/// - `end`:
///     End color.
/// - `t`:
///     Position in [0, 1].
///
/// Returns
/// -------
/// - `RGBColor`:
///     Interpolated color.
fn interpolate_rgb(start: RGBColor, end: RGBColor, t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let r = start.0 as f64 + (end.0 as f64 - start.0 as f64) * t;
    let g = start.1 as f64 + (end.1 as f64 - start.1 as f64) * t;
    let b = start.2 as f64 + (end.2 as f64 - start.2 as f64) * t;
    RGBColor(r as u8, g as u8, b as u8)
}

/// Draw a horizontal color legend with bordered swatches and labels.
///
/// Lays out min, max, center, and mean entries when provided, using the same
/// palette as the plot.
///
/// Parameters
/// ----------
/// - `legend_area`:
///     Drawing area reserved for the legend.
/// - `min_val`:
///     Minimum value label.
/// - `max_val`:
///     Maximum value label.
/// - `mean_val`:
///     Optional mean label.
/// - `center_val`:
///     Optional diverging center label.
/// - `min_color`:
///     Optional color for the minimum value. Defaults follow the heatmap palette.
/// - `max_color`:
///     Optional color for the maximum value. Defaults follow the heatmap palette.
/// - `symmetric_diverging`:
///     When true, legend swatches use the symmetric diverging scaling.
///
/// Returns
/// -------
/// - `Result<()>`:
///     Ok when rendering succeeds.
fn draw_color_legend<DB: DrawingBackend>(
    legend_area: &DrawingArea<DB, Shift>,
    min_val: f64,
    max_val: f64,
    mean_val: Option<f64>,
    center_val: Option<f64>,
    min_color: Option<RGBColor>,
    max_color: Option<RGBColor>,
    symmetric_diverging: bool,
) -> Result<()>
where
    DB::ErrorType: 'static + std::error::Error + Send + Sync,
{
    let (_, area_h) = legend_area.dim_in_pixel();
    legend_area.fill(&WHITE)?;

    let swatch_w: i32 = 32;
    let swatch_h: i32 = 18;
    let bottom_pad: i32 = 20; // Space reserved below the legend content
    let h_pad: i32 = 80; // Aligns with heatmap
    let x0: i32 = h_pad;
    let y0: i32 = area_h as i32 - swatch_h - bottom_pad;

    let mut items = vec![("min", min_val), ("max", max_val)];
    if let Some(center) = center_val {
        items.insert(1, ("center", center));
    }
    if let Some(mean) = mean_val {
        items.push(("mean", mean));
    }

    let mut x_cursor = x0;
    for (label, value) in items.iter() {
        let color = color_for_value(
            *value,
            min_val,
            max_val,
            center_val,
            min_color,
            max_color,
            symmetric_diverging,
        );
        let fill_style = ShapeStyle {
            color: color.to_rgba(),
            filled: true,
            stroke_width: 0,
        };
        legend_area.draw(&Rectangle::new(
            [(x_cursor, y0), (x_cursor + swatch_w, y0 + swatch_h)],
            fill_style,
        ))?;
        let border_style = ShapeStyle {
            color: BLACK.to_rgba(),
            filled: false,
            stroke_width: 1,
        };
        legend_area.draw(&Rectangle::new(
            [(x_cursor, y0), (x_cursor + swatch_w, y0 + swatch_h)],
            border_style,
        ))?;

        let text = format!("{}: {:.2}", label, value);
        let text_x = x_cursor + swatch_w + 8;
        let text_y = y0 + swatch_h / 2;
        let text_style = ("sans-serif", 16)
            .into_text_style(legend_area)
            .pos(Pos::new(HPos::Left, VPos::Center))
            .color(&BLACK);
        legend_area.draw(&Text::new(text, (text_x, text_y), text_style))?;

        let step = swatch_w + h_pad + 90;
        x_cursor += step;
    }

    Ok(())
}
