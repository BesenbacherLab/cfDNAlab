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
pub(crate) struct LengthVisualization {
    pub(crate) fragment_length: u32,
    pub(crate) tracks: Vec<Track>,
}

impl LengthVisualization {
    pub(crate) fn all_tracks_empty(&self) -> bool {
        self.tracks.iter().all(Track::is_empty)
    }
}

/// Parsed representation of the CLI configuration.
#[derive(Debug, Clone)]
pub(crate) struct VizConfig {
    pub(crate) position_specs: Vec<PositionalSelectionSpec>,
    pub(crate) bases: BasesFrom,
    pub(crate) mismatch_bases_from: MismatchBasesFrom,
    pub(crate) kmer_sizes: Option<Vec<u8>>,
    pub(crate) fragment_lengths: Vec<u32>,
    pub(crate) style: Style,
    pub(crate) width: usize,
    pub(crate) height: u32,
    pub(crate) output: Option<PathBuf>,
    pub(crate) label: Option<String>,
    pub(crate) show_index: bool,
    pub(crate) show_half: bool,
    pub(crate) show_mid: bool,
}
