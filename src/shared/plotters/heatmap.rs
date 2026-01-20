use anyhow::{Result, ensure};
use ndarray::Array2;
use plotters::{coord::Shift, prelude::*};
use std::path::Path;

/// Output formats supported by the heatmap writer.
pub enum HeatmapFormat {
    Png,
    Svg,
}

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
    width: u32,
    height: u32,
    format: HeatmapFormat,
) -> Result<()> {
    let x_edges = resolve_edges("x", x_edges, values.ncols())?;
    let y_edges = resolve_edges("y", y_edges, values.nrows())?;
    let mean_val = find_finite_mean(values);

    let (data_min, data_max) = find_finite_min_max(values);
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

    let out_path = out_path.as_ref();
    match format {
        HeatmapFormat::Png => {
            let drawing_area = BitMapBackend::new(out_path, (width, height)).into_drawing_area();
            draw_heatmap(
                &drawing_area,
                title,
                x_label,
                y_label,
                &x_edges,
                &y_edges,
                values,
                min_val,
                max_val,
                mean_val,
                val_center,
            )
        }
        HeatmapFormat::Svg => {
            let drawing_area = SVGBackend::new(out_path, (width, height)).into_drawing_area();
            draw_heatmap(
                &drawing_area,
                title,
                x_label,
                y_label,
                &x_edges,
                &y_edges,
                values,
                min_val,
                max_val,
                mean_val,
                val_center,
            )
        }
    }
}

fn draw_heatmap<DB: DrawingBackend>(
    drawing_area: &DrawingArea<DB, Shift>,
    title: &str,
    x_label: &str,
    y_label: &str,
    x_edges: &[f64],
    y_edges: &[f64],
    values: &Array2<f64>,
    min_val: f64,
    max_val: f64,
    mean_val: Option<f64>,
    center_val: Option<f64>,
) -> Result<()>
where
    DB::ErrorType: 'static + std::error::Error + Send + Sync,
{
    drawing_area.fill(&WHITE)?;

    let legend_height: u32 = 90;
    let (area_w, area_h) = drawing_area.dim_in_pixel();
    let (plot_area, legend_area) = if area_h > legend_height {
        let (upper, lower) = drawing_area.split_vertically(area_h - legend_height);
        (upper, Some(lower))
    } else {
        (drawing_area.clone(), None)
    };

    let x_range = *x_edges.first().unwrap()..*x_edges.last().unwrap();
    let y_range = *y_edges.first().unwrap()..*y_edges.last().unwrap();

    let x_label_area = 52;
    let y_label_area = 62;

    let mut chart = ChartBuilder::on(&plot_area)
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
            let color = color_for_value(value, min_val, max_val, center_val);
            chart.draw_series(std::iter::once(Rectangle::new(
                [(x0, y0), (x1, y1)],
                color.filled(),
            )))?;
        }
    }

    if let Some(legend_area) = legend_area {
        draw_color_legend(&legend_area, min_val, max_val, mean_val, center_val)?;
    }

    plot_area.present()?;
    Ok(())
}

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

fn color_for_value(value: f64, min_val: f64, max_val: f64, center_val: Option<f64>) -> RGBColor {
    if let Some(center) = center_val {
        if value <= center {
            let norm = ((value - min_val) / (center - min_val).max(f64::EPSILON)).clamp(0.0, 1.0);
            return interpolate_rgb(RGBColor(49, 130, 189), RGBColor(255, 255, 255), norm);
        } else {
            let norm = ((value - center) / (max_val - center).max(f64::EPSILON)).clamp(0.0, 1.0);
            return interpolate_rgb(RGBColor(255, 255, 255), RGBColor(203, 24, 29), norm);
        }
    }
    let norm = ((value - min_val) / (max_val - min_val).max(f64::EPSILON)).clamp(0.0, 1.0);
    interpolate_plasma(norm)
}

fn interpolate_rgb(start: RGBColor, end: RGBColor, t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let r = start.0 as f64 + (end.0 as f64 - start.0 as f64) * t;
    let g = start.1 as f64 + (end.1 as f64 - start.1 as f64) * t;
    let b = start.2 as f64 + (end.2 as f64 - start.2 as f64) * t;
    RGBColor(r as u8, g as u8, b as u8)
}

fn interpolate_plasma(t: f64) -> RGBColor {
    // Lightweight approximation of matplotlib plasma for quick contrast
    let t = t.clamp(0.0, 1.0);
    let r = (241.0 * t + 12.0 * (1.0 - t)) as u8;
    let g = (103.0 * t + 7.0 * (1.0 - t)) as u8;
    let b = (33.0 * t + 134.0 * (1.0 - t)) as u8;
    RGBColor(r, g, b)
}

fn draw_color_legend<DB: DrawingBackend>(
    legend_area: &DrawingArea<DB, Shift>,
    min_val: f64,
    max_val: f64,
    mean_val: Option<f64>,
    center_val: Option<f64>,
) -> Result<()>
where
    DB::ErrorType: 'static + std::error::Error + Send + Sync,
{
    let (area_w, area_h) = legend_area.dim_in_pixel();
    legend_area.fill(&WHITE)?;

    let swatch_w: i32 = 32;
    let swatch_h: i32 = 18;
    let v_pad: i32 = 8;
    let h_pad: i32 = 12;
    let x0: i32 = h_pad;
    let y0: i32 = v_pad;

    let mut items = vec![("min", min_val), ("max", max_val)];
    if let Some(center) = center_val {
        items.insert(1, ("center", center));
    }
    if let Some(mean) = mean_val {
        items.push(("mean", mean));
    }

    for (idx, (label, value)) in items.iter().enumerate() {
        let y = y0 + idx as i32 * (swatch_h + v_pad);
        let color = color_for_value(*value, min_val, max_val, center_val);
        let swatch_style = ShapeStyle {
            color: color.to_rgba(),
            filled: true,
            stroke_width: 1,
        };
        legend_area.draw(&Rectangle::new(
            [(x0, y), (x0 + swatch_w, y + swatch_h)],
            swatch_style,
        ))?;

        let text = format!("{}: {:.4}", label, value);
        let text_x = x0 + swatch_w + h_pad;
        let text_y = y + swatch_h - 2;
        legend_area.draw(&Text::new(
            text,
            (text_x, text_y),
            ("sans-serif", 16).into_font().color(&BLACK),
        ))?;
    }

    // Border for the legend region to distinguish it
    let content_height = y0 + (items.len() as i32) * (swatch_h + v_pad) + v_pad / 2;
    let border = Rectangle::new(
        [
            (0, 0),
            (area_w as i32 - 1, content_height.min(area_h as i32 - 1)),
        ],
        ShapeStyle::from(&BLACK).stroke_width(1),
    );
    legend_area.draw(&border)?;

    Ok(())
}
