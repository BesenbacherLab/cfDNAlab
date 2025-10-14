use std::num::NonZeroUsize;
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
#[cfg(feature = "cli")]
use clap::Parser;

use crate::pos_kmer_viz::{Anchor, Bases, Style, VizConfig, parse_lengths, parse_positions};

/// Draw which fragment bases will be counted for a given anchor and range setup.
///
/// Use this helper to prototype the “where to count” arguments (`--anchor`, `--positions`, `--step`, `--bases`),
/// before you run the counter on a BAM file. For every fragment length you request it renders the selected
/// bases as ASCII or SVG so you can sanity-check trims, interior windows, and where coverage concentrates. The
/// command is geometry-only: no BAM or reference reads are touched while you iterate.
///
/// Describe your selections with the shared 1-based inclusive grammar (`A-B`, `A:-B`, `:half`, `5..half-3`,
/// `-60..+60`, and friends) and the diagram will match the counting engine once it consumes the same arguments.
#[cfg_attr(feature = "cli", derive(Parser, Clone))]
pub struct VisualizeSelectedRegionConfig {
    /// Choose the reference frame that interprets every other flag `[left|right|per-end|nearest|mid|span]`.
    ///
    /// `left` counts bases from the forward 5′ end, `right` from the reverse 5′ end, `per-end` renders both
    /// of those tracks side-by-side, `nearest` folds the fragment so distances grow away from the closer end
    /// (`half` refers to floor(length/2)), `mid` centres the axis at 0 to emphasise symmetry, and `span`
    /// walks linearly from the left 5′ end to the right 5′ end. In every case you describe bases using
    /// 1-based inclusive indices relative to the chosen anchor.
    #[cfg_attr(
        feature = "cli",
        arg(long, value_enum, help_heading = "Region Selection")
    )]
    pub anchor: Anchor,

    /// Describe which bases remain after anchoring `[string]`.
    ///
    /// Ranges are written with 1-based inclusive bounds inside that frame: `A-B`, `A:`, `:B`, and `A:-B`
    /// for end-anchored systems (`left`, `right`, `per-end`, `span`). For the `nearest` anchor the grammar
    /// folds the fragment, so `:half` spans from distance 1 up to the fold point and `A..half-K` reads as
    /// “start at distance `A` and stop `K` bases before the fold” (e.g., `5..half-3` covers distances 5 through
    /// floor(length/2) − 3). Symmetric anchors like `mid` use forms such as `-60..+60`. Put differently:
    /// `1-10` keeps the first ten bases, `10:-10` trims both ends, `:half` reaches the fold point in `nearest`,
    /// `5..half-3` keeps only the interior band, and `-60..+60` sketches a ±60 bp window around the midpoint.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Region Selection"))]
    pub positions: String,

    /// Downsample after selection by keeping every Nth index `[integer ≥ 1]`.
    ///
    /// Applied independently to each track in anchor order (e.g., per-end left and right both stride through
    /// their own selections). Leave at 1 to keep every base.
    #[cfg_attr(
        feature = "cli",
        arg(long, default_value_t = 1, help_heading = "Region Selection")
    )]
    pub step: usize,

    /// Label the axis using read or reference coordinates `[read|reference]`.
    ///
    /// The visualization always uses geometry from the reference. When you pick `read`, the legend reminds
    /// you that downstream reports will show read-space positions; otherwise it notes that reference bases
    /// frame the interpretation.
    #[cfg_attr(
        feature = "cli",
        arg(long, value_enum, help_heading = "Region Selection")
    )]
    pub bases_from: Bases,

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
            help_heading = "Region Selection"
        )
    )]
    pub lengths: Option<Vec<u32>>,

    /// Generate a ladder of fragment lengths to sketch `[MIN:MAX[:STEP]]`.
    ///
    /// The default step is 10 when omitted (e.g., `80:200:20`). Conflicts with `--lengths`.
    /// Omit both `--lengths` and `--length-range` to fall back to `100:220:20`.
    #[cfg_attr(
        feature = "cli",
        arg(
            long,
            conflicts_with = "lengths",
            help_heading = "Region Selection",
            value_parser = clap::builder::NonEmptyStringValueParser::new()
        )
    )]
    pub length_range: Option<String>,

    /// Rendering backend for the diagram `[ascii|svg]`.
    ///
    /// ASCII is compact and stdout-friendly; SVG produces a figure for slides or docs.
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

    /// Mark the halfway point with `^` (ASCII) or a vertical line (SVG) `[flag]`.
    ///
    /// For `nearest`, the marker lands on `floor(length/2)` - the furthest folded distance before the ends meet.
    /// For `span`, it marks the halfway distance from the left 5′ end (length/2). Other anchors do not display it.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub show_half: bool,

    /// Mark the exact midpoint with `*` when the anchor exposes it (`mid` or `span`) `[flag]`.
    ///
    /// On the `mid` anchor this labels the `q=0` origin; on `span` it sits at the geometrical centre between
    /// the two ends. Combine with `--show-half` on `span` to see both the halfway distance and the true midpoint.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub show_mid: bool,
}

impl VisualizeSelectedRegionConfig {
    pub fn build(&self) -> Result<VizConfig> {
        let step = NonZeroUsize::new(self.step)
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

        let positions = parse_positions(self.anchor, &self.positions).with_context(|| {
            format!(
                "invalid --positions \"{}\" for anchor {}",
                self.positions,
                self.anchor.as_str()
            )
        })?;

        let width = self.width.unwrap_or(100);
        if width == 0 {
            return Err(anyhow!("--width must be positive (example: --width 120)"));
        }

        let height = self.height.unwrap_or(120);
        if height == 0 {
            return Err(anyhow!("--height must be positive (example: --height 160)"));
        }

        Ok(VizConfig {
            anchor: self.anchor,
            positions,
            positions_input: self.positions.clone(),
            step,
            bases: self.bases_from,
            fragment_lengths,
            style: self.style,
            width,
            height,
            output: self.output.clone(),
            label: self.label.clone(),
            show_index: self.show_index,
            show_half: self.show_half,
            show_mid: self.show_mid,
        })
    }
}
