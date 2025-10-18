use crate::commands::visualize_positions::config::VisualizePositionsConfig;
use crate::commands::visualize_positions::{
    BasesFrom, LengthVisualization, ReadClamp, Style, build_kmer_start_overlays,
    build_tracks_for_length, render_ascii, render_svg,
};
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

        if let Some(kmer_sizes) = &viz_cfg.kmer_sizes {
            if !kmer_sizes.is_empty() {
                let base_tracks = viz.tracks.clone();
                let overlays = build_kmer_start_overlays(
                    viz_cfg.frame,
                    length,
                    &viz_cfg.positions,
                    viz_cfg.step,
                    &base_tracks,
                    kmer_sizes,
                );
                viz.tracks.extend(overlays);
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
