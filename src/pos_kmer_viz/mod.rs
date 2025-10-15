//! Support code for the positional k-mer visualisation command.
//!
//! The public API is consumed by the `visualize-selected-region`
//! subcommand so tests and other tooling can exercise the parsing,
//! selection, and rendering primitives without going through Clap.
pub mod model;
pub mod parse;
pub mod render_ascii;
pub mod render_svg;
pub mod select;

pub use model::{
    BasesFrom, LengthVisualization, LinearRange, MidRange, NearestRange, MismatchBasesFrom,
    PositionsSpec, ReferenceFrame, Style, Track, VizConfig,
};
pub use parse::{RangeParseError, parse_lengths, parse_positions};
pub use render_ascii::render_ascii;
pub use render_svg::render_svg;
pub use select::{ReadClamp, build_tracks_for_length};
