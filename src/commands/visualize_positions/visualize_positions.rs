use crate::commands::visualize_positions::{
    BasesFrom, LengthVisualization, ReadClamp, ReferenceFrame, Style,
    build_nearest_guard_overlays, build_tracks_for_length, render_ascii, render_svg,
};
use crate::commands::visualize_positions::config::VisualizePositionsConfig;
use anyhow::Result;
use std::fs;
use std::io::{self, Write};

/// Execute the visualize-selected-region command.
pub fn run(cfg: &VisualizePositionsConfig) -> Result<()> {
    let viz_cfg = cfg.build()?;

    let mut results: Vec<LengthVisualization> = Vec::new();
    let clamp_mode = match viz_cfg.bases {
        BasesFrom::NearestRead => ReadClamp::Nearest,
        BasesFrom::Reads => ReadClamp::Both,
        _ => ReadClamp::None,
    };

    for &length in &viz_cfg.fragment_lengths {
        let mut viz = build_tracks_for_length(
            length,
            viz_cfg.frame,
            &viz_cfg.positions,
            viz_cfg.step,
            clamp_mode,
        );

        if viz_cfg.frame == ReferenceFrame::Nearest {
            if let Some(orders) = &viz_cfg.orders {
                if !orders.is_empty() {
                    let fragment_track = viz
                        .tracks
                        .iter()
                        .find(|track| track.name == "fragment")
                        .cloned();
                    let nearest_track = viz
                        .tracks
                        .iter()
                        .find(|track| track.name == "nearest")
                        .cloned();
                    if let (Some(fragment_track), Some(nearest_track)) =
                        (fragment_track, nearest_track)
                    {
                        let overlays = build_nearest_guard_overlays(
                            length,
                            &fragment_track,
                            &nearest_track,
                            orders,
                        );
                        viz.tracks.extend(overlays);
                    }
                }
            }
        }

        results.push(viz);
    }

    let rendered = match viz_cfg.style {
        Style::Ascii => render_ascii(&results, &viz_cfg),
        Style::Svg => render_svg(&results, &viz_cfg),
    };

    if let Some(path) = &viz_cfg.output {
        fs::write(path, rendered)?;
    } else {
        let mut stdout = io::stdout().lock();
        stdout.write_all(rendered.as_bytes())?;
    }

    Ok(())
}
