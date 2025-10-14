use std::fmt::Write;

use super::model::{Anchor, LengthVisualization, Track, VizConfig};

/// Render the visualization as an SVG string.
pub fn render_svg(results: &[LengthVisualization], config: &VizConfig) -> String {
    let width = config.width as f64;
    let mut height_estimate = 20.0;
    for viz in results {
        height_estimate += 18.0; // header line
        height_estimate += viz.tracks.len() as f64 * 24.0;
        if viz.all_tracks_empty() {
            height_estimate += 16.0;
        }
        height_estimate += 16.0; // block spacing
    }
    height_estimate += 20.0;
    let height = height_estimate.max(config.height as f64);

    let mut svg = String::new();
    writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{:.0}" height="{:.0}" font-family="monospace" font-size="12">"##,
        width, height
    )
    .ok();

    let mut y_cursor = 30.0;
    for viz in results {
        let header = format!(
            "L={} | anchor={} | positions={} | step={} | bases={}",
            viz.fragment_length,
            config.anchor.as_str(),
            config.positions_input,
            config.step.get(),
            config.bases.as_str()
        );
        writeln!(
            svg,
            r##"<text x="12" y="{:.1}" fill="#111">{}</text>"##,
            y_cursor, header
        )
        .ok();
        if let Some(label) = &config.label {
            writeln!(
                svg,
                r##"<text x="12" y="{:.1}" fill="#555">label: {}</text>"##,
                y_cursor + 14.0,
                label
            )
            .ok();
            y_cursor += 14.0;
        }
        y_cursor += 18.0;

        for track in &viz.tracks {
            draw_track_svg(
                &mut svg,
                track,
                viz.fragment_length,
                config,
                width,
                y_cursor,
            );
            y_cursor += 24.0;
        }

        if viz.all_tracks_empty() {
            writeln!(
                svg,
                r##"<text x="12" y="{:.1}" fill="#b91c1c">no positions selected for L={}</text>"##,
                y_cursor, viz.fragment_length
            )
            .ok();
            y_cursor += 20.0;
        }

        y_cursor += 16.0;
    }

    svg.push_str("</svg>\n");
    svg
}

fn draw_track_svg(
    svg: &mut String,
    track: &Track,
    fragment_length: u32,
    config: &VizConfig,
    full_width: f64,
    baseline_y: f64,
) {
    let margin_left = (full_width * 0.18).max(70.0).min(full_width * 0.4);
    let margin_right = 16.0;
    let bar_left = margin_left;
    let bar_width = (full_width - margin_left - margin_right).max(1.0);
    let bar_height = 10.0;
    let bar_top = baseline_y;
    let text_y = baseline_y + bar_height + 12.0;

    writeln!(
        svg,
        r##"<text x="12" y="{:.1}" fill="#111">{}</text>"##,
        baseline_y + bar_height - 2.0,
        track.name
    )
    .ok();

    writeln!(
        svg,
        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="#e5e7eb" stroke="#94a3b8" stroke-width="0.5"/>"##,
        bar_left,
        bar_top,
        bar_width,
        bar_height
    )
    .ok();

    for (start, end) in contiguous_segments(&track.selected_indices) {
        let x0 = value_to_px(start as f64, track, bar_left, bar_width);
        let x1 = value_to_px(end as f64, track, bar_left, bar_width);
        let width = (x1 - x0).abs().max(1.0);
        writeln!(
            svg,
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="#2563eb"/>"##,
            x0.min(x1),
            bar_top,
            width,
            bar_height
        )
        .ok();
    }

    if config.show_half {
        if config.anchor == Anchor::Nearest {
            let half = fragment_length / 2;
            if half > 0 {
                draw_marker(
                    svg,
                    '^',
                    half as f64,
                    track,
                    bar_left,
                    bar_top,
                    bar_width,
                    bar_height,
                );
            }
        } else if config.anchor == Anchor::Span {
            let half = fragment_length / 2;
            if half > 0 {
                draw_marker(
                    svg,
                    '^',
                    half as f64,
                    track,
                    bar_left,
                    bar_top,
                    bar_width,
                    bar_height,
                );
            }
        }
    }

    if config.show_mid {
        match config.anchor {
            Anchor::Mid => draw_marker(
                svg, '*', 0.0, track, bar_left, bar_top, bar_width, bar_height,
            ),
            Anchor::Span => {
                let mid = (track.axis.start as f64 + track.axis.end as f64) / 2.0;
                draw_marker(
                    svg, '*', mid, track, bar_left, bar_top, bar_width, bar_height,
                );
            }
            _ => {}
        }
    }

    writeln!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" fill="#475569">{}</text>"##,
        bar_left,
        text_y,
        format!("{}..{}", track.axis.start, track.axis.end)
    )
    .ok();
}

fn contiguous_segments(indices: &[i32]) -> Vec<(i32, i32)> {
    if indices.is_empty() {
        return Vec::new();
    }
    let mut segments = Vec::new();
    let mut start = indices[0];
    let mut prev = indices[0];
    for &value in &indices[1..] {
        if value == prev + 1 {
            prev = value;
        } else {
            segments.push((start, prev));
            start = value;
            prev = value;
        }
    }
    segments.push((start, prev));
    segments
}

fn value_to_px(value: f64, track: &Track, left: f64, width: f64) -> f64 {
    let axis_start = track.axis.start as f64;
    let axis_end = track.axis.end as f64;
    if axis_end <= axis_start {
        return left;
    }
    let ratio = ((value - axis_start) / (axis_end - axis_start)).clamp(0.0, 1.0);
    left + ratio * width
}

fn draw_marker(
    svg: &mut String,
    symbol: char,
    value: f64,
    track: &Track,
    bar_left: f64,
    bar_top: f64,
    bar_width: f64,
    bar_height: f64,
) {
    let x = value_to_px(value, track, bar_left, bar_width);
    writeln!(
        svg,
        r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#1f2937" stroke-width="1"/>"##,
        x,
        bar_top,
        x,
        bar_top + bar_height,
    )
    .ok();
    writeln!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" fill="#1f2937" text-anchor="middle">{}</text>"##,
        x,
        bar_top - 2.0,
        symbol
    )
    .ok();
}
