use std::path::PathBuf;

#[cfg(feature = "cli")]
use clap::ValueEnum;

use crate::commands::fragment_kmers::parse::PositionalSelectionSpec;
pub use crate::commands::fragment_kmers::positions::{
    BasesFrom, MismatchBasesFrom, ReferenceFrame,
};

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
    pub position_specs: Vec<PositionalSelectionSpec>,
    pub bases: BasesFrom,
    pub mismatch_bases_from: MismatchBasesFrom,
    pub kmer_sizes: Option<Vec<u8>>,
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
