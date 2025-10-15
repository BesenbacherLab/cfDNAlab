use std::num::NonZeroUsize;
use std::path::PathBuf;

#[cfg(feature = "cli")]
use clap::ValueEnum;

/// Enumeration of the available coordinate frames used to interpret positional selections.
///
/// The variants mirror the CLI keyword semantics documented in AGENTS.md.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceFrame {
    Span,
    Left,
    Right,
    #[cfg_attr(feature = "cli", value(alias = "per-end"))]
    PerEnd,
    Nearest,
    Mid,
}

impl ReferenceFrame {
    pub fn as_str(self) -> &'static str {
        match self {
            ReferenceFrame::Left => "left",
            ReferenceFrame::Right => "right",
            ReferenceFrame::PerEnd => "per-end",
            ReferenceFrame::Nearest => "nearest",
            ReferenceFrame::Mid => "mid",
            ReferenceFrame::Span => "span",
        }
    }
}

/// Whether the user wants to reason about read or reference coordinates.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BasesFrom {
    /// Prefer observed read coordinates, but fall back to the reference span when a read is missing.
    PreferRead,
    Read,
    Reference,
}

impl BasesFrom {
    pub fn as_str(self) -> &'static str {
        match self {
            BasesFrom::PreferRead => "prefer-read",
            BasesFrom::Read => "read",
            BasesFrom::Reference => "reference",
        }
    }
}

impl Default for BasesFrom {
    fn default() -> Self {
        BasesFrom::PreferRead
    }
}

/// Available rendering backends for the CLI.
#[cfg_attr(feature = "cli", derive(ValueEnum))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Style {
    Ascii,
    Svg,
}

impl Style {
    pub fn as_str(self) -> &'static str {
        match self {
            Style::Ascii => "ascii",
            Style::Svg => "svg",
        }
    }
}

/// Axis bounds used for rendering a track.
///
/// The axis is inclusive on both ends because the selections operate over
/// discrete base indices.
#[derive(Debug, Clone)]
pub struct AxisBounds {
    pub start: i32,
    pub end: i32,
}

impl AxisBounds {
    pub fn new(start: i32, end: i32) -> Self {
        Self { start, end }
    }

    pub fn length(&self) -> i32 {
        self.end - self.start
    }
}

/// A single visualization track (one logical coordinate system).
#[derive(Debug, Clone)]
pub struct Track {
    pub name: String,
    pub axis: AxisBounds,
    pub selected_indices: Vec<i32>,
}

impl Track {
    pub fn is_empty(&self) -> bool {
        self.selected_indices.is_empty()
    }
}

/// Per-fragment visualization data.
#[derive(Debug, Clone)]
pub struct LengthVisualization {
    pub fragment_length: u32,
    pub tracks: Vec<Track>,
}

impl LengthVisualization {
    pub fn all_tracks_empty(&self) -> bool {
        self.tracks.iter().all(Track::is_empty)
    }
}

/// Parsed representation of the CLI configuration.
#[derive(Debug, Clone)]
pub struct VizConfig {
    pub frame: ReferenceFrame,
    pub positions: PositionsSpec,
    pub positions_input: String,
    pub step: NonZeroUsize,
    pub bases: BasesFrom,
    pub fragment_lengths: Vec<u32>,
    pub style: Style,
    pub width: usize,
    pub height: u32,
    pub output: Option<PathBuf>,
    pub label: Option<String>,
    pub show_index: bool,
    pub show_half: bool,
    pub show_mid: bool,
}

/// Range grammar for frames that index strictly from one end.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LinearRange {
    /// Closed inclusive range `A..B`.
    Closed { start: u32, end: u32 },
    /// Open-right range `A:`.
    From { start: u32 },
    /// Open-left range `:B`.
    To { end: u32 },
    /// Opposite-end trimmed range `A..-B`.
    TrimOtherEnd { start: u32, other_end_trim: u32 },
}

/// Range grammar used with the `nearest` frame.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NearestRange {
    Closed { start: u32, end: u32 },
    From { start: u32 },
    ToHalf { minus: u32 },
    FromToHalf { start: u32, minus: u32 },
}

/// Range grammar used with the `mid` frame.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MidRange {
    Closed { neg: u32, pos: u32 },
    LeftOpen { neg: u32 },
    RightOpen { pos: u32 },
}

/// The position specification tagged to allow frame-specific dispatch.
#[derive(Debug, Clone)]
pub enum PositionsSpec {
    Linear(LinearRange),
    Nearest(NearestRange),
    Mid(MidRange),
}
