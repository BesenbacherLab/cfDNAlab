use anyhow::Result;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::utils::coverage::tiled_run::format_number_simplify;

/// Open a zstd-compressed buffered writer once per file
///
/// Use when you plan to write many lines to the same file
/// Keeps the concrete encoder type hidden behind `Box<dyn Write>` which is fine for IO-bound paths
#[inline]
pub fn open_zstd_auto_writer<P: AsRef<Path>>(
    path: P,
    level: i32,
    n_threads: Option<u32>,
) -> Result<BufWriter<Box<dyn Write>>> {
    let file = std::fs::File::create(path)?;
    let mut enc = zstd::Encoder::new(file, level)?; // Create once
    if let Some(num_threads) = n_threads {
        enc.multithread(num_threads).ok();
    }
    let sink: Box<dyn Write> = Box::new(enc.auto_finish()); // Auto-finishing wrapper
    Ok(BufWriter::new(sink))
}

/// Write a final aggregate row: `chromosome  start  end  value  blacklisted_positions`
#[inline]
pub fn write_final_row<W: Write>(
    w: &mut W,
    chr: &str,
    start: u64,
    end: u64,
    value: f64,
    blacklisted_positions: u64,
    decimals: i32,
) -> anyhow::Result<()> {
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}",
        chr,
        start,
        end,
        format_number_simplify(value, decimals),
        blacklisted_positions
    )?;
    Ok(())
}
