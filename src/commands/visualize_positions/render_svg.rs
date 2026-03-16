use std::cmp::Ordering;
use std::fmt::Write;

use crate::commands::fragment_kmers::positions::ReferenceFrame;

use super::model::{LengthVisualization, Track, VizConfig};

const CHAR_WIDTH: f64 = 7.0;
const MARKER_BAND: f64 = 12.0;
const BAR_HEIGHT: f64 = 10.0;
const INDEX_BAND: f64 = 28.0;
const INDEX_LABEL_PAD: f64 = 4.0;
const LABEL_BAND: f64 = 12.0;
const LABEL_COLUMN_PADDING: f64 = 60.0;
const FRAGMENT_PADDING: f64 = 30.0;

fn svg_track_label(track: &Track, config: &VizConfig) -> String {
    if config.position_specs[0].frame == ReferenceFrame::Nearest && track.name == "nearest" {
        let max_val = track.axis.end.max(track.axis.start);
        format!("{} (max distance {})", track.name, max_val)
    } else {
        track.name.clone()
    }
}

/// Render the visualization as an SVG string.
pub fn render_svg(results: &[LengthVisualization], config: &VizConfig) -> String {
    let width = config.width as f64;
    let mut height_estimate = 20.0;
    let per_track_height = track_block_height(config);
    for viz in results {
        height_estimate += 18.0; // header line
        height_estimate += viz.tracks.len() as f64 * per_track_height;
        if viz.all_tracks_empty() {
            height_estimate += 16.0;
        }
        height_estimate += FRAGMENT_PADDING; // block spacing
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
        let mut header = format!(
            "L={} | frame={} | positions={} | step={} | bases={} | mismatches={}",
            viz.fragment_length,
            config
                .position_specs
                .iter()
                .map(|ps| ps.frame.as_str().to_string())
                .collect::<Vec<String>>()
                .join(","),
            config
                .position_specs
                .iter()
                .map(|ps| ps.positions_string.clone())
                .collect::<Vec<String>>()
                .join(","),
            config
                .position_specs
                .iter()
                .map(|ps| ps.step.get().to_string())
                .collect::<Vec<String>>()
                .join(","),
            config.bases.as_str(),
            config.mismatch_bases_from.as_str()
        );
        if let Some(kmer_sizes) = &config.kmer_sizes
            && !kmer_sizes.is_empty()
        {
            let list = kmer_sizes
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<_>>()
                .join(",");
            write!(header, " | k-mer-sizes={}", list).ok();
        }
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

        let max_label_chars = viz
            .tracks
            .iter()
            .map(|track| svg_track_label(track, config).chars().count())
            .max()
            .unwrap_or(0);

        for track in &viz.tracks {
            let advance = draw_track_svg(
                &mut svg,
                track,
                viz.fragment_length,
                config,
                width,
                y_cursor,
                max_label_chars,
            );
            y_cursor += advance;
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

        y_cursor += FRAGMENT_PADDING;
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
    label_char_width: usize,
) -> f64 {
    let track_label = svg_track_label(track, config);
    let effective_chars = label_char_width.max(track_label.chars().count());
    let label_space = effective_chars as f64 * CHAR_WIDTH + LABEL_COLUMN_PADDING;
    let max_margin = (full_width * 0.4).max(LABEL_COLUMN_PADDING);
    let margin_right = 16.0;
    let mut margin_left = label_space.clamp(LABEL_COLUMN_PADDING, max_margin);
    if margin_left > full_width - margin_right - 10.0 {
        margin_left = (full_width - margin_right - 10.0).max(12.0);
    }
    let bar_left = margin_left;
    let bar_width = (full_width - margin_left - margin_right).max(1.0);
    let bar_top = baseline_y + MARKER_BAND;
    let bar_height = BAR_HEIGHT;

    writeln!(
        svg,
        r##"<text x="12" y="{:.1}" fill="#111">{}</text>"##,
        bar_top + bar_height - 2.0,
        track_label
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

    let markers = axis_markers(track, fragment_length, config);
    for (value, symbol) in markers {
        draw_axis_marker(svg, symbol, value, track, bar_left, bar_top, bar_width);
    }

    let axis_bottom = bar_top + bar_height;
    if config.show_index {
        draw_tick_marks(svg, track, bar_left, axis_bottom, bar_width);
    }

    let text_y = if config.show_index {
        axis_bottom + INDEX_BAND + INDEX_LABEL_PAD
    } else {
        axis_bottom + LABEL_BAND
    };

    writeln!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" fill="#475569">{}</text>"##,
        bar_left,
        text_y,
        format!("{}..{}", track.axis.start, track.axis.end)
    )
    .ok();

    track_block_height(config)
}

fn track_block_height(config: &VizConfig) -> f64 {
    MARKER_BAND
        + BAR_HEIGHT
        + if config.show_index {
            INDEX_BAND + INDEX_LABEL_PAD
        } else {
            LABEL_BAND
        }
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

fn draw_axis_marker(
    svg: &mut String,
    symbol: char,
    value: f64,
    track: &Track,
    bar_left: f64,
    bar_top: f64,
    bar_width: f64,
) {
    let x = value_to_px(value, track, bar_left, bar_width);
    let top = (bar_top - MARKER_BAND + 6.0).max(4.0);
    let bottom = bar_top - 2.0;
    writeln!(
        svg,
        r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#1f2937" stroke-width="1"/>"##,
        x, top, x, bottom
    )
    .ok();
    writeln!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" fill="#1f2937" text-anchor="middle">{}</text>"##,
        x,
        top - 2.0,
        symbol
    )
    .ok();
}

fn draw_tick_marks(
    svg: &mut String,
    track: &Track,
    bar_left: f64,
    axis_bottom: f64,
    bar_width: f64,
) {
    let start = track.axis.start;
    let end = track.axis.end;
    if start > end {
        return;
    }

    let mut candidates = Vec::new();
    for value in start..=end {
        if should_mark_tick(value, start, end) {
            candidates.push((value, tick_priority(value, start, end)));
        }
    }
    if candidates.is_empty() {
        return;
    }

    for (value, _) in &candidates {
        let x = value_to_px(*value as f64, track, bar_left, bar_width);
        writeln!(
            svg,
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#1f2937" stroke-width="1"/>"##,
            x,
            axis_bottom,
            x,
            axis_bottom + 6.0
        )
        .ok();
    }

    let mut placed_labels = Vec::new();
    let mut occupied: Vec<(f64, f64)> = Vec::new();
    candidates.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    for (value, _) in candidates {
        let x = value_to_px(value as f64, track, bar_left, bar_width);
        let label = value.to_string();
        let text_width = label.len() as f64 * CHAR_WIDTH;
        let half = text_width / 2.0;
        let left = x - half;
        let right = x + half;
        if occupied
            .iter()
            .any(|&(occupied_left, occupied_right)| left < occupied_right && right > occupied_left)
        {
            continue;
        }
        occupied.push((left, right));
        placed_labels.push((x, label));
    }
    placed_labels.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    for (x, label) in placed_labels {
        writeln!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" fill="#334155" font-size="11" text-anchor="middle">{}</text>"##,
            x,
            axis_bottom + 20.0,
            label
        )
        .ok();
    }
}

fn should_mark_tick(value: i32, start: i32, end: i32) -> bool {
    value == start || value == end || value % 10 == 0
}

fn tick_priority(value: i32, start: i32, end: i32) -> u8 {
    if value == start || value == end { 2 } else { 1 }
}

fn axis_markers(track: &Track, fragment_length: u32, config: &VizConfig) -> Vec<(f64, char)> {
    let mut markers = Vec::new();
    if config.show_half {
        match config.position_specs[0].frame {
            ReferenceFrame::Nearest => {
                let half = (fragment_length / 2) as f64;
                if half > 0.0 {
                    push_marker(&mut markers, half.max(track.axis.start as f64), '^');
                }
            }
            ReferenceFrame::Left | ReferenceFrame::Right | ReferenceFrame::PerEnd => {
                let half = (fragment_length / 2) as f64;
                if half > 0.0 {
                    push_marker(&mut markers, half.max(track.axis.start as f64), '^');
                }
            }
            _ => {}
        }
    }
    if config.show_mid {
        let mid = match config.position_specs[0].frame {
            ReferenceFrame::Mid => Some(0.0),
            ReferenceFrame::Nearest => Some((fragment_length / 2) as f64),
            ReferenceFrame::Left | ReferenceFrame::Right | ReferenceFrame::PerEnd => {
                Some((track.axis.start as f64 + track.axis.end as f64) / 2.0)
            }
        };
        if let Some(value) = mid {
            push_marker(&mut markers, value, '*');
        }
    }
    markers
}

fn push_marker(markers: &mut Vec<(f64, char)>, value: f64, symbol: char) {
    let priority = marker_priority(symbol);
    for (existing_value, existing_symbol) in markers.iter_mut() {
        if (*existing_value - value).abs() <= 1.0 {
            if marker_priority(*existing_symbol) < priority {
                *existing_symbol = symbol;
                *existing_value = value;
            }
            return;
        }
    }
    markers.push((value, symbol));
}

fn marker_priority(symbol: char) -> u8 {
    match symbol {
        '*' => 2,
        '^' => 1,
        _ => 0,
    }
}
