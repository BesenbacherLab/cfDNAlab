//! Writers for small genomic scaling-factor inputs.
//!
//! cfDNAlab scaling-factor files are TSV files with `chromosome`, `start`,
//! `end`, and `scaling_factor` columns. These helpers write tiny inputs for
//! tests that exercise scaling-aware commands. Coordinates are zero-based and
//! half-open.

use anyhow::Result;
use std::{fs::File, io::Write, path::Path};

/// One row in a cfDNAlab scaling-factor TSV input.
///
/// Rows are written unchanged by `write_scaling_factors_tsv`. This lets tests
/// construct both valid and intentionally invalid files at the call site. Use
/// command-level validation when the test is about rejecting bad scaling input.
///
/// The row represents a zero-based half-open interval on `chromosome` with a
/// multiplicative scaling factor applied by commands that support regional
/// scaling. The helper does not require the interval to be ordered or the
/// scaling factor to be finite.
#[derive(Clone, Debug, PartialEq)]
pub struct ScalingFactorRow {
    /// Chromosome or contig name.
    pub chromosome: String,
    /// Zero-based half-open start coordinate.
    pub start: u64,
    /// Zero-based half-open end coordinate.
    pub end: u64,
    /// Multiplicative scaling factor.
    pub scaling_factor: f32,
}

impl ScalingFactorRow {
    /// Create a scaling-factor row.
    ///
    /// `start` and `end` are zero-based half-open coordinates. The scaling
    /// factor is stored unchanged and later written with Rust's standard display
    /// formatting for `f32`.
    pub fn new(chromosome: impl Into<String>, start: u64, end: u64, scaling_factor: f32) -> Self {
        Self {
            chromosome: chromosome.into(),
            start,
            end,
            scaling_factor,
        }
    }
}

/// Write a cfDNAlab scaling-factor TSV input.
///
/// The file starts with the header expected by cfDNAlab commands:
/// `chromosome`, `start`, `end`, and `scaling_factor`. Rows are written in the
/// order supplied by the caller.
///
/// Empty `rows` creates a header-only file. The helper does not sort, merge,
/// deduplicate, validate interval order, or reject unusual numeric scaling
/// factors.
pub fn write_scaling_factors_tsv<P: AsRef<Path>>(path: P, rows: &[ScalingFactorRow]) -> Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "chromosome\tstart\tend\tscaling_factor")?;
    for row in rows {
        writeln!(
            file,
            "{}\t{}\t{}\t{}",
            row.chromosome, row.start, row.end, row.scaling_factor
        )?;
    }
    Ok(())
}
