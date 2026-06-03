//! Writers for small BED inputs.
//!
//! These helpers are intentionally small and explicit. They write tab-separated
//! BED rows without adding sorting, merging, or coordinate conversion. Callers
//! are expected to pass zero-based half-open coordinates, matching ordinary BED
//! conventions and cfDNAlab command inputs.

use anyhow::Result;
use std::{fs::File, io::Write, path::Path};

/// One BED4 row.
///
/// `chrom`, `start`, and `end` are written unchanged. `name` is the fourth BED
/// column and is useful for grouped-window command tests where expected counts
/// are grouped by window name.
///
/// This type does not validate interval ordering, empty names, or chromosome
/// naming. Some tests intentionally create invalid files to check command
/// errors, and this helper should not hide those inputs before the command sees
/// them.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bed4Row {
    /// Chromosome or contig name.
    pub chrom: String,
    /// Zero-based half-open start coordinate.
    pub start: u64,
    /// Zero-based half-open end coordinate.
    pub end: u64,
    /// BED name field.
    pub name: String,
}

impl Bed4Row {
    /// Create a BED4 row.
    ///
    /// Coordinates are zero-based and half-open. Values are stored as provided
    /// and written unchanged by `write_bed4`. In particular, `start > end` is
    /// allowed here so command-level validation can be tested explicitly.
    pub fn new(chrom: impl Into<String>, start: u64, end: u64, name: impl Into<String>) -> Self {
        Self {
            chrom: chrom.into(),
            start,
            end,
            name: name.into(),
        }
    }
}

/// Write BED4 rows to a tab-separated file.
///
/// The output has no header. Each row is written as `chrom`, `start`, `end`,
/// and `name`. Use this for small command inputs where the expected windows are
/// easier to derive directly from the call site than from a committed fixture
/// file.
///
/// Rows are written in caller-supplied order. Empty `rows` creates an empty
/// file. The helper does not sort, merge, deduplicate, clamp, or validate
/// coordinates.
pub fn write_bed4<P: AsRef<Path>>(path: P, rows: &[Bed4Row]) -> Result<()> {
    let mut file = File::create(path)?;
    for row in rows {
        writeln!(
            file,
            "{}\t{}\t{}\t{}",
            row.chrom, row.start, row.end, row.name
        )?;
    }
    Ok(())
}
