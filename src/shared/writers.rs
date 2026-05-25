use anyhow::{Context, Result};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Open a zstd-compressed buffered writer once per file
///
/// Use when you plan to write many lines to the same file.
///
/// Keeps the concrete encoder type hidden behind `Box<dyn Write>` which is fine for IO-bound paths
pub(crate) fn open_zstd_auto_writer<P: AsRef<Path>>(
    path: P,
    level: i32,
    n_threads: Option<u32>,
) -> Result<BufWriter<Box<dyn Write>>> {
    // Borrow once, reuse the &Path
    let path_ref = path.as_ref();

    let file = std::fs::File::create(path_ref)
        .with_context(|| format!("creating {}", path_ref.display()))?;

    let mut enc = zstd::Encoder::new(file, level)?; // Create once
    if let Some(num_threads) = n_threads.filter(|&num_threads| num_threads > 1) {
        enc.multithread(num_threads).ok();
    }
    let sink: Box<dyn Write> = Box::new(enc.auto_finish()); // Auto-finishing wrapper
    Ok(BufWriter::with_capacity(512 * 1024, sink))
}
