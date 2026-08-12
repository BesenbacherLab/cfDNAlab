//! Writers for small fragment-level outlier keep-weight inputs.
//!
//! cfDNAlab outlier-weight files are sparse TSV files with `chromosome`, `start`, `end`, and
//! `keep_weight` columns. Coordinates are zero-based and half-open. Positions omitted from the
//! file have implicit weight `1.0`.

use anyhow::Result;
use std::{fs::File, io::Write, path::Path};

/// A row in a cfDNAlab fragment-level outlier keep-weight TSV.
///
/// Rows are written unchanged so tests can construct valid and intentionally invalid inputs.
/// Command-level validation remains responsible for interval ordering, overlap, contig bounds,
/// and the allowed weight range.
#[derive(Clone, Debug, PartialEq)]
pub struct OutlierWeightRow {
    /// Chromosome or contig name.
    pub chromosome: String,
    /// Zero-based half-open start coordinate.
    pub start: u64,
    /// Zero-based half-open end coordinate.
    pub end: u64,
    /// Fragment keep weight assigned by any positive overlap.
    pub keep_weight: f64,
}

impl OutlierWeightRow {
    /// Create an outlier keep-weight row without validating it.
    pub fn new(chromosome: impl Into<String>, start: u64, end: u64, keep_weight: f64) -> Self {
        Self {
            chromosome: chromosome.into(),
            start,
            end,
            keep_weight,
        }
    }
}

/// Write a sparse fragment-level outlier keep-weight TSV.
///
/// Rows are written in caller order after metadata describing the implicit identity weight and
/// minimum-over-overlaps rule. Empty `rows` creates a valid header-only sparse file.
pub fn write_outlier_weights_tsv<P: AsRef<Path>>(path: P, rows: &[OutlierWeightRow]) -> Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "# omitted_keep_weight=1.0")?;
    writeln!(file, "# fragment_overlap_rule=minimum_keep_weight")?;
    writeln!(file, "chromosome\tstart\tend\tkeep_weight")?;
    for row in rows {
        writeln!(
            file,
            "{}\t{}\t{}\t{}",
            row.chromosome, row.start, row.end, row.keep_weight
        )?;
    }
    Ok(())
}
