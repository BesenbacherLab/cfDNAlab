use std::num::NonZeroUsize;
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
#[cfg(feature = "cli")]
use clap::Parser;

use crate::commands::{
    cli_common::FragmentPositionSelectionArgs,
    visualize_positions::{
        BasesFrom, ReferenceFrame, Style, VizConfig, parse_lengths, parse_positions,
    },
};

/// `fragment-kmers` helper: Draw which fragment bases will be counted for a given frame and range setup.
///
/// Use this helper to prototype the “where to count” arguments (`--frame`, `--positions`, `--step`, `--bases-from`, `--mismatch-bases-from`),
/// before you run `cfdna fragment-kmers` on a BAM file. For every fragment length you request, the selected
/// bases are rendered as ASCII or SVG, so you can check the correct positions are counted at. The
/// command is geometry-only: no BAM or reference reads are touched while you iterate.
///
/// Describe your selections with the **1-based inclusive** grammar (`A..B`, `A..-B`, `..half`, `5..half-3`,
/// `-60..60` (`mid`-frame-only), and friends) and the diagram will show the regions counted by
/// `cfdna fragment-kmers`, assuming the same arguments are passed.
#[cfg_attr(feature = "cli", derive(Parser, Clone))]
pub struct VisualizePositionsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    pub position_selection: FragmentPositionSelectionArgs,

    /// Explicit fragment lengths to sketch (comma-separated) `[integers]`.
    ///
    /// Handy for inspecting a bespoke menu of lengths (e.g., `90,123,200`). Conflicts with `--length-range`.
    #[cfg_attr(
        feature = "cli",
        arg(
            long,
            value_delimiter = ',',
            num_args = 1..,
            value_parser = clap::value_parser!(u32).range(1..),
            conflicts_with = "length_range",
            help_heading = "Visualization"
        )
    )]
    pub lengths: Option<Vec<u32>>,

    /// Generate a ladder of fragment lengths to sketch `[MIN:MAX[:STEP]]`.
    ///
    /// The default step is 20 when omitted (e.g., `80:200:20`). Conflicts with `--lengths`.
    /// Omit both `--lengths` and `--length-range` to fall back to `100:220:20`.
    #[cfg_attr(
        feature = "cli",
        arg(
            long,
            conflicts_with = "lengths",
            help_heading = "Visualization",
            value_parser = clap::builder::NonEmptyStringValueParser::new()
        )
    )]
    pub length_range: Option<String>,

    /// Rendering backend for the diagram `[ascii|svg]`.
    ///
    /// ASCII is compact and stdout-friendly. SVG produces a figure for slides or docs.
    #[cfg_attr(
        feature = "cli",
        arg(long, value_enum, default_value_t = Style::Ascii, help_heading = "Visualization")
    )]
    pub style: Style,

    /// Width of the plotted track in characters (ASCII) or pixels (SVG) `[integer > 0]`.
    ///
    /// Wider widths give finer horizontal resolution.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub width: Option<usize>,

    /// Height of the SVG canvas in pixels (ignored for ASCII output) `[integer > 0]`.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub height: Option<u32>,

    /// Optional file to write the visualization to `[path]`.
    ///
    /// Without `--output`, both ASCII and SVG are printed to stdout so you can pipe or preview inline.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub output: Option<PathBuf>,

    /// Free-form label appended to the header to annotate the sketch `[string]`.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub label: Option<String>,

    /// Show numeric tick marks alongside the ASCII axis `[flag]`.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub show_index: bool,

    /// Mark the halfway distance (`floor(length/2)` from the frame origin) with `^` on the axis `[flag]`.
    ///
    /// We deliberately use `floor(length/2)` rather than the mathematical midpoint so that “half” always refers to the
    /// largest part entirely contained within the first half of the fragment. This keeps ranges such as `..half` and
    /// `half+1..` disjoint, leaving the exact midpoint base (for odd lengths) to the separate `*` marker. For
    /// `nearest`, the preview line shows the full fragment and places `^` that far along the folded track. For the
    /// linear frames (`left`, `right`, `per-end`) the marker appears exactly `floor(length/2)` bases away from the
    /// selected end, reinforcing the “halfway distance” interpretation even when the true center is the next base.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub show_half: bool,

    /// Hide the fragment midpoint marker (`*`) on the axis `[flag]`.
    ///
    /// The midpoint marker denotes the exact center base that sits between both halves. It is drawn by default whenever
    /// the frame exposes this conceptual center (`mid` at 0, `left`/`right`/`per-end` at the central coordinate, `nearest`
    /// where the fold reaches its apex). Use this flag to suppress it.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub hide_mid: bool,
}

impl VisualizePositionsConfig {
    pub fn build(&self) -> Result<VizConfig> {
        let step = NonZeroUsize::new(self.position_selection.step)
            .ok_or_else(|| anyhow!("--step must be at least 1 (example: --step 3)"))?;

        let fragment_lengths = if let Some(list) = &self.lengths {
            list.clone()
        } else {
            let ladder_spec = self.length_range.as_deref().unwrap_or("100:220:20");
            parse_lengths(None, Some(ladder_spec)).context("failed to parse fragment lengths")?
        };

        if fragment_lengths.is_empty() {
            return Err(anyhow!(
                "no fragment lengths provided; use --lengths or --length-range"
            ));
        }

        let positions = parse_positions(
            self.position_selection.frame,
            &self.position_selection.positions,
        )
        .with_context(|| {
            format!(
                "invalid --positions \"{}\" for frame {}",
                self.position_selection.positions,
                self.position_selection.frame.as_str()
            )
        })?;

        if self.position_selection.frame == ReferenceFrame::Mid
            && self.position_selection.bases_from == BasesFrom::NearestRead
        {
            return Err(anyhow!(
                "`--bases-from nearest-read` is incompatible with the `mid` frame. Choose a different bases-from mode."
            ));
        }

        let width = self.width.unwrap_or(100);
        if width == 0 {
            return Err(anyhow!("--width must be positive (example: --width 120)"));
        }

        let height = self.height.unwrap_or(120);
        if height == 0 {
            return Err(anyhow!("--height must be positive (example: --height 160)"));
        }

        Ok(VizConfig {
            frame: self.position_selection.frame,
            positions,
            positions_input: self.position_selection.positions.clone(),
            step,
            bases: self.position_selection.bases_from,
            mismatch_bases_from: self.position_selection.mismatch_bases_from,
            fragment_lengths,
            style: self.style,
            width,
            height,
            output: self.output.clone(),
            label: self.label.clone(),
            show_index: self.show_index,
            show_half: self.show_half,
            show_mid: !self.hide_mid,
        })
    }
}
