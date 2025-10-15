use std::fs;
use std::io::{self, Write};

use anyhow::Result;

use crate::pos_kmer_viz::{
    LengthVisualization, Style, build_tracks_for_length, render_ascii, render_svg,
};

use super::config::VisualizeSelectedRegionConfig;

/// Execute the visualize-selected-region command.
pub fn run(cfg: &VisualizeSelectedRegionConfig) -> Result<()> {
    let viz_cfg = cfg.build()?;

    let mut results: Vec<LengthVisualization> = Vec::new();
    for &length in &viz_cfg.fragment_lengths {
        let viz = build_tracks_for_length(length, viz_cfg.frame, &viz_cfg.positions, viz_cfg.step);
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
