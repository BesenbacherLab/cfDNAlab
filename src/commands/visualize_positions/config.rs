use std::num::NonZeroUsize;
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
#[cfg(feature = "cli")]
use clap::Parser;

use crate::pos_kmer_viz::{
    BasesFrom, OverlapResolution, ReferenceFrame, Style, VizConfig, parse_lengths, parse_positions,
};

/// `fragment-kmers` helper: Draw which fragment bases will be counted for a given frame and range setup.
///
/// Use this helper to prototype the “where to count” arguments (`--frame`, `--positions`, `--step`, `--bases-from`, `--overlap-resolution`),
/// before you run `cfdna fragment-kmers` on a BAM file. For every fragment length you request, the selected
/// bases are rendered as ASCII or SVG, so you can check the correct positions are counted at. The
/// command is geometry-only: no BAM or reference reads are touched while you iterate.
///
/// Describe your selections with the **1-based inclusive** grammar (`A..B`, `A..-B`, `..half`, `5..half-3`,
/// `-60..60` (`mid`-frame-only), and friends) and the diagram will show the regions counted by
/// `cfdna fragment-kmers`, assuming the same arguments are passed.
#[cfg_attr(feature = "cli", derive(Parser, Clone))]
pub struct VisualizeSelectedRegionConfig {
    /// Choose the reference frame that interprets every other region selection argument `[left|right|per-end|nearest|mid]`.
    ///
    /// Note: `--positions` describe positions to count at relative to the chosen frame.
    /// Some frames are only relevant when `fragment-kmers` return positionally indexed counts.
    ///
    /// - **`left`** counts bases from the forward 5' end. Indices increase along the fragment and
    ///   k-mers are counted in the forward-orientation.
    ///
    /// - **`right`** counts bases from the reverse 5' end. Indices decrease along the fragment and
    ///   k-mers are counted in the reverse-orientation with **complemented** bases.
    ///
    /// - **`per-end`** counts both `left` and `right` simultaneously, producing two sets of k-mer counts.
    ///   The `step` start can differ per side.
    ///
    /// - **`nearest`** folds the fragment around the midpoint so distances grow away from the nearest end.
    ///   The positional keyword `half` represents the midpoint (and maximum position).
    ///   Bases contributed by the reverse 5' side are complemented.
    ///
    /// - **`mid`** centres the axis on the midpoint, allowing selections around zero with negative/positive offsets.
    ///   K-mers are counted in the forward-orientation.
    #[cfg_attr(
        feature = "cli",
        arg(long, value_enum, default_value_t = ReferenceFrame::Left, help_heading = "Region Selection")
    )]
    pub frame: ReferenceFrame,

    /// Describe which positions to count at relative to the selected frame `[string]`.
    ///
    /// Indices are **1-based inclusive**, why e.g. `1..10` would start at the first position and end at the tenth position (included).
    ///
    /// The allowed shapes depend on `--frame`:
    ///
    /// - **`left`**, **`right`**, **`per-end`**: use `A..B`, `A..`, `..B`, `A..-B`, `..half`, or `A..half-K`.
    ///   For example, `1..10` keeps the first ten bases, `10..-10` trims both ends, and `..half-5`
    ///   includes bases from the start up to five before the fragment midpoint. Open intervals like `A..`
    ///   include every coordinate from `A` to the end of the frame.
    ///
    /// - **`nearest`** (folded 1..length/2): use `A..B`, `A..`, `..B`, `..half`, or `A..half-K`. Here, `half` expands to the
    ///   largest folded distance (ties are randomly assigned for even-length fragments), ensuring the centre base is
    ///   maximally counted once. Forms like `10..-10` are rejected for this frame.
    ///
    /// - **`mid`** (centered at 0): use `-M..N`, `-M..`, or `..N`. E.g. `-10..10` for the 20 bases around the midpoint.
    #[cfg_attr(
        feature = "cli",
        arg(long, help_heading = "Region Selection", allow_hyphen_values = true)
    )]
    pub positions: String,

    /// Downsample after selection by keeping every Nth index `[integer >= 1]`.
    ///
    /// Applied independently to each track in frame order (e.g., per-end left and right both stride through
    /// their own selections). Leave at 1 to keep every base.
    ///
    /// For the `mid` frame, zero is treated as the origin of the stride: when the chosen range includes the
    /// midpoint, it is always retained and every `step`th offset is kept symmetrically
    /// (`-2*step`, `-step`, `0`, `step`, `2*step`, ...). Ranges that exclude the origin fall back to the default stride.
    #[cfg_attr(
        feature = "cli",
        arg(long, default_value_t = 1, help_heading = "Region Selection")
    )]
    pub step: usize,

    /// Choose which coordinate source defines the counted positions `[prefer-read|read|reference|nearest-read]`.
    ///
    /// - `prefer-reads`: Use read-space coordinates whenever an observed base covers the requested position
    ///   and fall back to the reference span when reads don't cover the positions.
    ///
    /// - `reads`: Only count positions the reads cover.
    ///
    /// - `reference`: Always use the reference span, even when reads do not cover those bases.
    ///
    /// - `nearest-read`: Clamp the selection to the read that corresponds to the frame origin (e.g., the
    ///   left read for the `left` frame).
    #[cfg_attr(
        feature = "cli",
        arg(
            long,
            value_enum,
            default_value_t = BasesFrom::PreferReads,
            help_heading = "Region Selection"
        )
    )]
    pub bases_from: BasesFrom,

    /// Resolve overlapping read mismatches when preferring read bases `[nearest-read|base-quality|reference]`.
    ///
    /// - `nearest-read`: Take the base from whichever read is closest to the frame origin. **NOTE**: Incompatible with `--frame mid`.
    ///
    /// - `base-quality`: Take the base with the highest quality score.
    ///
    /// - `reference`: Ignore the reads and fall back to the reference base for that coordinate.
    #[cfg_attr(
        feature = "cli",
        arg(
            long = "overlap-resolution",
            value_enum,
            default_value_t = OverlapResolution::NearestRead,
            help_heading = "Region Selection"
        )
    )]
    pub overlap_resolution: OverlapResolution,

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
    /// The default step is 10 when omitted (e.g., `80:200:20`). Conflicts with `--lengths`.
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

    /// Mark the halfway distance (length/2 from the frame origin) with `^` on the axis `[flag]`.
    ///
    /// For `nearest`, the preview line highlights the full fragment and marks `length/2` before the folded track.
    /// For linear frames (`left`, `right`, `per-end`), the mark lands at `length/2` from their respective
    /// origin. This differs from the fragment midpoint (`*`), which is the conceptual centre point of the fragment.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub show_half: bool,

    /// Hide the fragment midpoint marker (`*`) on the axis `[flag]`.
    ///
    /// The midpoint marker is drawn by default whenever the frame exposes the conceptual centre (`mid` at 0,
    /// `left`/`right`/`per-end` at the halfway coordinate, `nearest` at the folded maximum). Use this flag to suppress it.
    #[cfg_attr(feature = "cli", arg(long, help_heading = "Visualization"))]
    pub hide_mid: bool,
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

        let positions = parse_positions(self.frame, &self.positions).with_context(|| {
            format!(
                "invalid --positions \"{}\" for frame {}",
                self.positions,
                self.frame.as_str()
            )
        })?;

        if self.frame == ReferenceFrame::Mid && self.bases_from == BasesFrom::NearestRead {
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
            frame: self.frame,
            positions,
            positions_input: self.positions.clone(),
            step,
            bases: self.bases_from,
            overlap_resolution: self.overlap_resolution,
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
