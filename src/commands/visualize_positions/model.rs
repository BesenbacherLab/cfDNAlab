use std::path::PathBuf;

#[cfg(feature = "cli")]
use clap::ValueEnum;

use crate::commands::fragment_kmers::parse::PositionalSelectionSpec;
use crate::shared::positioning::{BasesFrom, MismatchBasesFrom};
use crate::shared::visualization::Track;

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
