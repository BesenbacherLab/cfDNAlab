use crate::shared::temp_chrom_names::temp_chrom_token;
use anyhow::{Context, Result};
use fxhash::FxHashMap;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

/// Buffered writer paired with a chromosome-scoped temporary BED file.
///
/// Streaming each chromosome through its own writer keeps memory usage low and
/// lets the pipeline defer global policies until the final concatenation pass.
///
/// The writer wraps a 1 MiB `BufWriter<File>` pointing at an opaque-token path under the run's
/// temp directory. Callers borrow the writer when writing rows and later reuse the stored path
/// during concatenation.
#[derive(Debug)]
pub(crate) struct ChromTempWriter {
    path: PathBuf,
    writer: BufWriter<File>,
}

impl ChromTempWriter {
    #[inline]
    pub(crate) fn writer(&mut self) -> &mut BufWriter<File> {
        &mut self.writer
    }
}

/// Ensure a chromosome-specific temp writer exists and return it.
///
/// The helper creates the writer lazily the first time a chromosome appears and
/// reuses it for subsequent chunks.
///
/// It assigns an opaque chromosome token, opens the file inside `temp_dir`, and wraps it in a
/// 1 MiB buffer to keep I/O efficient.
///
/// # Parameters
/// - `chrom`: chromosome identifier (e.g., `chr7`).
/// - `temp_dir`: per-run temporary directory.
/// - `temp_writers`: map of already-open writers.
///
/// # Returns
/// Mutable reference to the `ChromTempWriter` for `chrom`.
pub(crate) fn ensure_temp_writer_for_chrom<'a>(
    chrom: &str,
    temp_dir: &Path,
    temp_writers: &'a mut FxHashMap<String, ChromTempWriter>,
) -> Result<&'a mut ChromTempWriter> {
    if !temp_writers.contains_key(chrom) {
        let token = temp_chrom_token(temp_writers.len());
        let path = temp_dir.join(format!("chrom.{token}.bed.tmp"));
        let file = File::create(&path)
            .with_context(|| format!("creating temp file for chromosome {}", chrom))?;
        let writer = BufWriter::with_capacity(1 << 20, file);
        temp_writers.insert(chrom.to_string(), ChromTempWriter { path, writer });
    }
    Ok(temp_writers.get_mut(chrom).unwrap())
}

/// Flush all chromosome writers and surface their file paths.
///
/// After streaming finishes we flush each buffer to disk and collect the temp
/// file paths so the concatenation pass can replay them.
///
/// Each writer is flushed in place, the path is cloned, and the map entry is
/// removed; the returned vector preserves the map's insertion order.
///
/// # Parameters
/// - `temp_writers`: map of chromosome -> writer.
///
/// # Returns
/// Vector of `(chromosome, path)` pairs ready for concatenation.
pub(crate) fn finalize_temp_writers(
    temp_writers: &mut FxHashMap<String, ChromTempWriter>,
) -> Result<Vec<(String, PathBuf)>> {
    let mut entries: Vec<(String, PathBuf)> = Vec::with_capacity(temp_writers.len());
    for (chrom, chrom_writer) in temp_writers.iter_mut() {
        chrom_writer.writer.flush()?;
        entries.push((chrom.clone(), chrom_writer.path.clone()));
    }
    temp_writers.clear();
    Ok(entries)
}
