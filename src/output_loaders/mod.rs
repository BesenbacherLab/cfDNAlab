//! Public loaders for files produced by cfDNAlab commands.
//!
//! These APIs read command outputs from disk and return loader structs with
//! parsed metadata and count values. They are separate from `run_like_cli`:
//! command runners produce files, while `output_loaders` opens files that have
//! already been written.
//!
//! Loader entry points:
//!
//! ```text
//! load_lengths_output(path)   -> LengthsOutput
//! load_fcoverage_output(path) -> FCoverageOutput
//! load_ends_output(path)      -> EndsOutput
//! load_midpoints_output(path) -> MidpointsOutput
//! ```

mod common;

#[cfg(feature = "cmd_ends")]
mod ends;

#[cfg(feature = "cmd_fcoverage")]
mod fcoverage;

#[cfg(feature = "cmd_lengths")]
mod lengths;

#[cfg(feature = "cmd_midpoints")]
mod midpoints;

pub use common::{DenseMatrix, LengthBin, WindowRow};

#[cfg(feature = "cmd_midpoints")]
pub use common::DenseArray3;

#[cfg(feature = "cmd_ends")]
pub use ends::{
    EndMotifAxisKind, EndMotifCountSelection, EndMotifCountsData, EndMotifGroupRow,
    EndMotifRowMetadata, EndMotifRowMode, EndMotifSparseCountLookup, EndMotifSparseCounts,
    EndMotifSparseEntry, EndMotifStorageMode, EndMotifWindowMode, EndsOutput, EndsSelector,
    load_ends_output,
};

#[cfg(feature = "cmd_fcoverage")]
pub use fcoverage::{
    FCoverageCoefficientOfVariation, FCoverageData, FCoverageGroupRow, FCoverageOutput,
    FCoverageRowMetadata, FCoverageRowMode, FCoverageSelection, FCoverageSelector, FCoverageSignal,
    FCoverageSummaryStats, FCoverageSummaryStatsSelection, FCoverageValueMode,
    FCoverageValueSelection, FCoverageWindowRow, load_fcoverage_output,
    load_fcoverage_output_with_group_index,
};

#[cfg(feature = "cmd_lengths")]
pub use lengths::{
    LengthCountSelection, LengthGroupRow, LengthOutputMode, LengthRowMetadata, LengthsOutput,
    LengthsSelector, load_lengths_output,
};

#[cfg(feature = "cmd_midpoints")]
pub use midpoints::{
    MidpointCountSelection, MidpointGroupRow, MidpointPositionBin, MidpointsOutput,
    MidpointsSelector, load_midpoints_output,
};
