//! LIONHEART-derived average overlapping fragment length normalization.
//!
//! The command performs one tiled genomic sweep, stores additive sufficient statistics, runs two
//! in-memory mixture fits, and writes a Zarr lookup package. fcoverage applies that package with a
//! bounded rolling adaptor that preserves the normal fragment iterator's exact order.

pub(crate) mod config;
pub(crate) mod inference;
pub(crate) mod model;
pub(crate) mod overlapping_lengths_correction;
pub(crate) mod package;
pub(crate) mod reducer;
mod scipy_bfgs;
pub(crate) mod tiling;
