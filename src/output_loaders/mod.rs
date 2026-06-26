//! Public loaders for files produced by cfDNAlab commands.
//!
//! These APIs read command outputs from disk and return loader structs with
//! parsed metadata and count values. They are separate from `run_like_cli`:
//! command runners produce files, while `output_loaders` opens files that have
//! already been written.
//!
//! Loader entry points are compiled with the matching command feature:
//!
//! ```text
//! cmd_lengths   -> load_lengths_output(path)
//! cmd_fcoverage -> load_fcoverage_output(path)
//! cmd_ends      -> load_ends_output(path)
//! cmd_midpoints -> load_midpoints_output(path)
//! cmd_ref_kmers -> load_ref_kmers_output(path)
//! ```
//!
//! Public loader methods return `OutputLoaderResult<T>`. Error messages include
//! the path, line, array, or selector context needed to fix malformed cfDNAlab
//! outputs.

mod common;
mod error;

#[cfg(feature = "cmd_ends")]
mod ends;

#[cfg(feature = "cmd_fcoverage")]
mod fcoverage;

#[cfg(feature = "cmd_lengths")]
mod lengths;

#[cfg(feature = "cmd_midpoints")]
mod midpoints;

#[cfg(feature = "cmd_ref_kmers")]
mod ref_kmers;

pub use common::{DenseMatrix, LengthBin, WindowRow};
pub use error::{OutputLoaderError, OutputLoaderResult};

#[cfg(feature = "cmd_midpoints")]
pub use common::DenseArray3;

#[cfg(feature = "cmd_ends")]
pub use ends::{
    EndMotifAxisKind, EndMotifCountSelection, EndMotifCountsData, EndMotifGroupRow,
    EndMotifOutputMetadata, EndMotifRowMetadata, EndMotifRowMode, EndMotifSparseCountLookup,
    EndMotifSparseCounts, EndMotifSparseEntry, EndMotifStorageMode, EndMotifWindowMode, EndsOutput,
    EndsSelector, load_ends_output,
};

#[cfg(feature = "cmd_fcoverage")]
pub use fcoverage::{
    FCoverageAggregationBasis, FCoverageCoefficientOfVariation, FCoverageData, FCoverageDataMode,
    FCoverageFilenameMetadata, FCoverageGroupRow, FCoverageLengthNormalization, FCoverageOutput,
    FCoverageOutputMetadata, FCoverageRowMetadata, FCoverageRowMode, FCoverageSelection,
    FCoverageSelector, FCoverageSignal, FCoverageSummaryStats, FCoverageSummaryStatsSelection,
    FCoverageValueMode, FCoverageValueSelection, FCoverageWindowRow, load_fcoverage_output,
    load_fcoverage_output_with_group_index,
};

#[cfg(feature = "cmd_lengths")]
pub use lengths::{
    LengthCountSelection, LengthGroupRow, LengthOutputMetadata, LengthOutputMode,
    LengthRowMetadata, LengthsOutput, LengthsSelector, load_lengths_output,
};

#[cfg(feature = "cmd_midpoints")]
pub use midpoints::{
    MidpointCountSelection, MidpointGroupRow, MidpointPositionBin, MidpointsOutput,
    MidpointsOutputMetadata, MidpointsSelector, load_midpoints_output,
};

#[cfg(feature = "cmd_ref_kmers")]
pub use ref_kmers::{
    RefKmerFrequencyData, RefKmerFrequencySelection, RefKmerGroupRow, RefKmerMotifAxisKind,
    RefKmerOutputMetadata, RefKmerRowMetadata, RefKmerRowMode, RefKmerSparseCountEntry,
    RefKmerSparseFrequencies, RefKmerSparseFrequencyEntry, RefKmerStorageMode, RefKmerWindowMode,
    RefKmersOutput, RefKmersSelector, load_ref_kmers_output,
};
