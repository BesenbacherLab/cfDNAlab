use std::fmt::Write;

use super::model::{LengthVisualization, ReferenceFrame, Track, VizConfig};

/// Render the visualization as ASCII art.
pub fn render_ascii(results: &[LengthVisualization], config: &VizConfig) -> String {
    let mut output = String::new();
    for (idx, viz) in results.iter().enumerate() {
        if idx > 0 {
            output.push('\n');
        }
        write_header(&mut output, viz, config);
        if let Some(first_track) = viz.tracks.first() {
            let markers = axis_markers(first_track, viz.fragment_length, config);
            let marker_columns = marker_columns(first_track, config.width, &markers);

            let mut axis_chars: Vec<char> = build_ruler(config.width).chars().collect();
            for &(column, symbol) in &marker_columns {
                if column < axis_chars.len() {
                    axis_chars[column] = symbol;
                }
            }
            output.push_str("axis  : ");
            output.push_str(&axis_chars.into_iter().collect::<String>());
            output.push('\n');

            if config.show_index {
                let (ticks, labels) = build_tick_lines(first_track, config.width);
                let mut ticks_chars: Vec<char> = ticks.chars().collect();
                for &(column, _) in &marker_columns {
                    if column < ticks_chars.len() {
                        ticks_chars[column] = '|';
                    }
                }
                output.push_str("ticks : ");
                output.push_str(&ticks_chars.into_iter().collect::<String>());
                output.push('\n');
                output.push_str("index : ");
                output.push_str(&labels);
                output.push('\n');
            }
        }

        for track in &viz.tracks {
            let bar = build_track_bar(track, config);
            output.push_str(&format!("{:>6}: ", track.name));
            output.push_str(&bar);
            output.push('\n');
        }

        if viz.all_tracks_empty() {
            output.push_str("note  : no positions selected for L=");
            output.push_str(&viz.fragment_length.to_string());
            output.push('\n');
        }
    }
    output
}

fn write_header(buffer: &mut String, viz: &LengthVisualization, config: &VizConfig) {
    let mut line = format!(
        "L={}  | frame={}  positions={}  step={}  bases={}  overlap={}",
        viz.fragment_length,
        config.frame.as_str(),
        config.positions_input,
        config.step.get(),
        config.bases.as_str(),
        config.mismatch_bases_from.as_str()
    );
    if let Some(label) = &config.label {
        write!(line, "  label={}", label).ok();
    }
    buffer.push_str(&line);
    buffer.push('\n');
}

fn build_ruler(width: usize) -> String {
    if width == 0 {
        return String::new();
    }
    let mut chars = vec!['-'; width];
    if let Some(first) = chars.first_mut() {
        *first = '|';
    }
    if width > 1 {
        chars[width - 1] = '|';
    }
    chars.into_iter().collect()
}

pub fn build_tick_lines(track: &Track, width: usize) -> (String, String) {
    if width == 0 {
        return (String::new(), String::new());
    }
    let mut ticks = vec![' '; width];
    let mut labels = vec![' '; width];
    let mut chosen: Vec<Option<(i32, u8)>> = vec![None; width];

    let start = track.axis.start;
    let end = track.axis.end;
    if start <= end {
        for value in start..=end {
            if should_mark_tick(value, start, end) {
                let column = value_to_column(value as f64, start as f64, end as f64, width);
                if column < width {
                    let priority = tick_priority(value, start, end);
                    let slot = &mut chosen[column];
                    let should_replace = match slot {
                        None => true,
                        Some((existing_value, existing_priority)) => {
                            priority > *existing_priority
                                || (priority == *existing_priority
                                    && value == end
                                    && *existing_value != end)
                        }
                    };
                    if should_replace {
                        *slot = Some((value, priority));
                    }
                }
            }
        }
    }

    for (column, slot) in chosen.into_iter().enumerate() {
        if let Some((value, _)) = slot {
            ticks[column] = '|';
            place_label(&mut labels, column, value);
        }
    }

    (ticks.into_iter().collect(), labels.into_iter().collect())
}

fn should_mark_tick(value: i32, start: i32, end: i32) -> bool {
    value == start || value == end || value % 10 == 0
}

fn tick_priority(value: i32, start: i32, end: i32) -> u8 {
    if value == start || value == end { 2 } else { 1 }
}

fn place_label(line: &mut [char], column: usize, value: i32) {
    let label = value.to_string();
    let len = label.len();
    if len == 0 || len > line.len() {
        return;
    }
    let start = column.saturating_sub(len.saturating_sub(1));
    for (idx, ch) in label.chars().enumerate() {
        let pos = start + idx;
        if pos < line.len() {
            line[pos] = ch;
        }
    }
}

fn build_track_bar(track: &Track, config: &VizConfig) -> String {
    if config.width == 0 {
        return String::new();
    }
    let mut cells = vec!['.'; config.width];
    let axis_start = track.axis.start as f64;
    let axis_end = track.axis.end as f64;
    for &index in &track.selected_indices {
        let column = value_to_column(index as f64, axis_start, axis_end, config.width);
        if column < cells.len() {
            cells[column] = '#';
        }
    }

    cells.into_iter().collect()
}

pub fn value_to_column(value: f64, axis_start: f64, axis_end: f64, width: usize) -> usize {
    if width == 0 {
        return 0;
    }
    if axis_end <= axis_start {
        return 0;
    }
    let span = axis_end - axis_start;
    let mut ratio = (value - axis_start) / span;
    if ratio.is_nan() {
        ratio = 0.0;
    }
    ratio = ratio.clamp(0.0, 1.0);
    let max_index = (width - 1) as f64;
    let scaled = ratio * max_index;
    scaled.round().clamp(0.0, max_index) as usize
}

fn axis_markers(track: &Track, fragment_length: u32, config: &VizConfig) -> Vec<(f64, char)> {
    let mut markers = Vec::new();
    if config.show_half {
        match config.frame {
            ReferenceFrame::Nearest => {
                let half = (fragment_length / 2) as f64;
                if half > 0.0 {
                    markers.push((half.max(track.axis.start as f64), '^'));
                }
            }
            ReferenceFrame::Left | ReferenceFrame::Right | ReferenceFrame::PerEnd => {
                let half = (fragment_length / 2) as f64;
                if half > 0.0 {
                    markers.push((half.max(track.axis.start as f64), '^'));
                }
            }
            _ => {}
        }
    }
    if config.show_mid {
        let mid = match config.frame {
            ReferenceFrame::Mid => Some(0.0),
            ReferenceFrame::Nearest => Some((fragment_length / 2) as f64),
            ReferenceFrame::Left | ReferenceFrame::Right | ReferenceFrame::PerEnd => {
                Some((track.axis.start as f64 + track.axis.end as f64) / 2.0)
            }
        };
        if let Some(value) = mid {
            markers.push((value, '*'));
        }
    }
    markers
}

fn marker_columns(track: &Track, width: usize, markers: &[(f64, char)]) -> Vec<(usize, char)> {
    if width == 0 {
        return Vec::new();
    }
    let axis_start = track.axis.start as f64;
    let axis_end = track.axis.end as f64;
    if axis_end <= axis_start {
        return Vec::new();
    }
    markers
        .iter()
        .filter_map(|(value, symbol)| {
            if *value < axis_start || *value > axis_end {
                None
            } else {
                let column = value_to_column(*value, axis_start, axis_end, width);
                Some((column, *symbol))
            }
        })
        .collect()
}
