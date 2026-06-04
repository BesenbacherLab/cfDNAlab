use std::str::FromStr;

use crate::shared::interval::{IndexedInterval, Interval};

/// What to do per window
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PeaksWindowAction {
    #[default]
    Stats,
    OnlyIncludeThesePositionsUnique,
    OnlyIncludeThesePositionsIndexed,
}

// For the CLI
impl FromStr for PeaksWindowAction {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        if s == "stats" {
            Ok(PeaksWindowAction::Stats)
        } else if s == "unique-positions" {
            Ok(PeaksWindowAction::OnlyIncludeThesePositionsUnique)
        } else if s == "indexed-positions" {
            Ok(PeaksWindowAction::OnlyIncludeThesePositionsIndexed)
        } else {
            Err("Use 'stats', 'indexed-positions', or 'unique-positions'".into())
        }
    }
}

/// Collection of statistics about WPS peaks.
#[derive(Debug, Clone)]
pub struct PeakStats {
    pub count: u32,
    pub avg_distance: f32,
    pub median_distance: u32,
}

/// Per-window WPS peak value
#[derive(Debug, Clone)]
pub enum WindowPeaksValue {
    /// Average coverage in the window
    Stats(PeakStats),
    /// Positional coverage for every base in the window, left->right
    Positions(Vec<f32>),
}

/// One window's result (keeps original ordering info)
#[derive(Debug, Clone)]
pub struct WindowPeaksResult {
    pub window: IndexedInterval<u64>,
    pub value: WindowPeaksValue,
    pub num_blacklisted_pos: Option<u32>,
}

impl WindowPeaksResult {
    #[inline]
    pub fn start(&self) -> u64 {
        self.window.start()
    }

    #[inline]
    pub fn end(&self) -> u64 {
        self.window.end()
    }

    #[inline]
    pub fn original_idx(&self) -> u64 {
        self.window.idx()
    }
}

/// Top-level result for a run with or without windows
#[derive(Debug, Clone)]
pub enum PeaksOutput {
    /// Results for each input window
    PerWindow {
        action: PeaksWindowAction,
        results: Vec<WindowPeaksResult>,
    },
    /// No windows given -> return positional peaks for the whole sequence
    WholePositional {
        /// Covered span, typically `[0, length)`
        interval: Interval<u64>,
        /// Per-base coverage, left->right
        values: Vec<f32>,
    },
}
