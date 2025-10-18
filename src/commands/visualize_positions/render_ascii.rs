use std::fmt::Write;

use super::model::{LengthVisualization, ReferenceFrame, Track, VizConfig};

/// Render the visualization as ASCII art.
pub fn render_ascii(results: &[LengthVisualization], config: &VizConfig) -> String {
    let mut output = String::new();
    for (idx, viz) in results.iter().enumerate() {
        if idx > 0 {
            output.push('\n');
        }
        let mut label_width = viz
            .tracks
            .iter()
            .map(|track| track.name.len())
            .max()
            .unwrap_or(0);
        for static_label in ["axis", "ticks", "index", "note"] {
            label_width = label_width.max(static_label.len());
        }
        let reference_axis = viz.tracks.first().map(|track| track.axis.clone());
        if let Some(reference_axis) = &reference_axis {
            for track in &viz.tracks {
                if track.axis.start != reference_axis.start || track.axis.end != reference_axis.end
                {
                    label_width = label_width.max(axis_label_for_track(track, config).len());
                    if config.show_index {
                        label_width = label_width.max(ticks_label_for_track(track).len());
                        label_width = label_width.max(index_label_for_track(track).len());
                    }
                }
            }
        }
        write_header(&mut output, viz, config, label_width);

        if let Some(first_track) = viz.tracks.first() {
            let markers = axis_markers(first_track, viz.fragment_length, config);
            let marker_columns = marker_columns(first_track, config.width, &markers);

            let mut axis_chars: Vec<char> = build_ruler(config.width).chars().collect();
            for &(column, symbol) in &marker_columns {
                if column < axis_chars.len() {
                    axis_chars[column] = symbol;
                }
            }
            write!(output, "{:>width$}: ", "axis", width = label_width).ok();
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
                write!(output, "{:>width$}: ", "ticks", width = label_width).ok();
                output.push_str(&ticks_chars.into_iter().collect::<String>());
                output.push('\n');
                write!(output, "{:>width$}: ", "index", width = label_width).ok();
                output.push_str(&labels);
                output.push('\n');
            }
        }

        for track in &viz.tracks {
            if let Some(reference_axis) = &reference_axis {
                if track.axis.start != reference_axis.start || track.axis.end != reference_axis.end
                {
                    write_track_axis(&mut output, track, viz.fragment_length, config, label_width);
                }
            }
            let bar = build_track_bar(track, config);
            write!(output, "{:>width$}: ", track.name, width = label_width).ok();
            output.push_str(&bar);
            if config.frame == ReferenceFrame::Nearest && track.name == "nearest" {
                write!(
                    output,
                    " | max distance {}",
                    track.axis.end.max(track.axis.start)
                )
                .ok();
            }
            output.push('\n');
        }

        if viz.all_tracks_empty() {
            write!(output, "{:>width$}: ", "note", width = label_width).ok();
            output.push_str("no positions selected for L=");
            output.push_str(&viz.fragment_length.to_string());
            output.push('\n');
        }
    }
    output
}

fn write_header(
    buffer: &mut String,
    viz: &LengthVisualization,
    config: &VizConfig,
    label_width: usize,
) {
    let mut line = format!("L={}", viz.fragment_length);
    let target_width = label_width + 2;
    if line.len() < target_width {
        line.push_str(&" ".repeat(target_width - line.len()));
    } else {
        line.push(' ');
    }
    write!(
        line,
        "| frame={}  positions={}  step={}  bases={}  mismatches={}",
        config.frame.as_str(),
        config.positions_input,
        config.step.get(),
        config.bases.as_str(),
        config.mismatch_bases_from.as_str()
    )
    .ok();
    if let Some(orders) = &config.orders {
        if !orders.is_empty() {
            let list = orders
                .iter()
                .map(|order| order.to_string())
                .collect::<Vec<_>>()
                .join(",");
            write!(line, "  orders={}", list).ok();
        }
    }
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

    let mut placements = Vec::new();
    for (column, slot) in chosen.into_iter().enumerate() {
        if let Some((value, priority)) = slot {
            ticks[column] = '|';
            placements.push((column, value, priority));
        }
    }

    placements.sort_by(|a, b| b.2.cmp(&a.2).then(a.0.cmp(&b.0)));
    for (column, value, _) in placements {
        try_place_label(&mut labels, column, value);
    }

    (ticks.into_iter().collect(), labels.into_iter().collect())
}

fn should_mark_tick(value: i32, start: i32, end: i32) -> bool {
    value == start || value == end || value % 10 == 0
}

fn tick_priority(value: i32, start: i32, end: i32) -> u8 {
    if value == start || value == end { 2 } else { 1 }
}

fn try_place_label(line: &mut [char], column: usize, value: i32) -> bool {
    let label = value.to_string();
    let len = label.len();
    if len == 0 || len > line.len() {
        return false;
    }
    let start = column.saturating_sub(len.saturating_sub(1));
    if start >= line.len() {
        return false;
    }
    let end = start + len;
    if end > line.len() {
        return false;
    }
    if (start..end).any(|pos| line[pos] != ' ') {
        return false;
    }
    for (idx, ch) in label.chars().enumerate() {
        let pos = start + idx;
        if pos < line.len() {
            line[pos] = ch;
        }
    }
    true
}

fn build_track_bar(track: &Track, config: &VizConfig) -> String {
    if config.width == 0 {
        return String::new();
    }
    let mut cells = vec!['.'; config.width];
    if track.selected_indices.is_empty() {
        return cells.into_iter().collect();
    }

    if config.width == 1 {
        cells[0] = '#';
        return cells.into_iter().collect();
    }

    let axis_start = track.axis.start as f64;
    let axis_end = track.axis.end as f64;
    if axis_end <= axis_start {
        for &index in &track.selected_indices {
            let column = value_to_column(index as f64, axis_start, axis_end, config.width);
            if column < cells.len() {
                cells[column] = '#';
            }
        }
        return cells.into_iter().collect();
    }

    let mut run_start: Option<i32> = None;
    let mut previous_value: Option<i32> = None;
    for &value in &track.selected_indices {
        if run_start.is_none() {
            run_start = Some(value);
            previous_value = Some(value);
            continue;
        }
        if let Some(prev) = previous_value {
            if value <= prev {
                fill_run_columns(
                    &mut cells,
                    run_start.unwrap(),
                    prev,
                    axis_start,
                    axis_end,
                    config.width,
                );
                run_start = Some(value);
            } else if value == prev + 1 {
                // Continue the current contiguous run
            } else {
                fill_run_columns(
                    &mut cells,
                    run_start.unwrap(),
                    prev,
                    axis_start,
                    axis_end,
                    config.width,
                );
                run_start = Some(value);
            }
        }
        previous_value = Some(value);
    }

    if let (Some(start), Some(end)) = (run_start, previous_value) {
        fill_run_columns(&mut cells, start, end, axis_start, axis_end, config.width);
    }

    cells.into_iter().collect()
}

fn write_track_axis(
    output: &mut String,
    track: &Track,
    fragment_length: u32,
    config: &VizConfig,
    label_width: usize,
) {
    let label = axis_label_for_track(track, config);
    let markers = axis_markers(track, fragment_length, config);
    let marker_columns = marker_columns(track, config.width, &markers);

    let mut axis_chars: Vec<char> = build_ruler(config.width).chars().collect();
    for &(column, symbol) in &marker_columns {
        if column < axis_chars.len() {
            axis_chars[column] = symbol;
        }
    }
    write!(output, "{:>width$}: ", label, width = label_width).ok();
    output.push_str(&axis_chars.into_iter().collect::<String>());
    output.push('\n');

    if config.show_index {
        let (ticks, labels) = build_tick_lines(track, config.width);
        let mut ticks_chars: Vec<char> = ticks.chars().collect();
        for &(column, _) in &marker_columns {
            if column < ticks_chars.len() {
                ticks_chars[column] = '|';
            }
        }
        let ticks_label = ticks_label_for_track(track);
        write!(output, "{:>width$}: ", ticks_label, width = label_width).ok();
        output.push_str(&ticks_chars.into_iter().collect::<String>());
        output.push('\n');
        let index_label = index_label_for_track(track);
        write!(output, "{:>width$}: ", index_label, width = label_width).ok();
        output.push_str(&labels);
        output.push('\n');
    }
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

fn fill_run_columns(
    cells: &mut [char],
    run_start: i32,
    run_end: i32,
    axis_start: f64,
    axis_end: f64,
    width: usize,
) {
    if width == 0 {
        return;
    }
    let start = run_start.min(run_end);
    let end = run_start.max(run_end);
    if width == 1 {
        cells[0] = '#';
        return;
    }
    if axis_end <= axis_start {
        let column = value_to_column(start as f64, axis_start, axis_end, width);
        if column < cells.len() {
            cells[column] = '#';
        }
        return;
    }

    let start_col = value_to_column(start as f64, axis_start, axis_end, width);
    let end_col = value_to_column(end as f64, axis_start, axis_end, width);
    let (lower, upper) = if start_col <= end_col {
        (start_col, end_col)
    } else {
        (end_col, start_col)
    };
    for column in lower..=upper {
        if column < cells.len() {
            cells[column] = '#';
        }
    }
}

fn axis_label_for_track(track: &Track, config: &VizConfig) -> String {
    if config.frame == ReferenceFrame::Nearest && track.name == "nearest" {
        let max_val = track.axis.end.max(track.axis.start);
        format!("axis({} max={})", track.name, max_val)
    } else {
        format!("axis({})", track.name)
    }
}

fn ticks_label_for_track(track: &Track) -> String {
    format!("ticks({})", track.name)
}

fn index_label_for_track(track: &Track) -> String {
    format!("index({})", track.name)
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
