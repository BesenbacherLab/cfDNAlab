use crate::commands::prepare_windows::{
    config::PrepareConfig,
    labels::{LabelKey, LabelSchema, build_tuple_compositions, render_label_for_key},
    prepare_windows::Window,
};
use crate::shared::io::{TextWriter, create_text_writer, stdout_text_writer};
use anyhow::{Context, Result};
use fxhash::FxHashMap;
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

/// Buffered writer paired with a chromosome-scoped temporary BED file.
///
/// Streaming each chromosome through its own writer keeps memory usage low and
/// lets the pipeline defer global policies until the final concatenation pass.
///
/// The writer wraps a 1 MiB `BufWriter<File>` pointing at a sanitized path under
/// the run's temp directory. Callers borrow the writer when writing rows and
/// later reuse the stored path during concatenation.
#[derive(Debug)]
pub struct ChromTempWriter {
    path: PathBuf,
    writer: BufWriter<File>,
}

impl ChromTempWriter {
    #[inline]
    pub fn writer(&mut self) -> &mut BufWriter<File> {
        &mut self.writer
    }
}

/// Serialize finalized windows as minimal BED-like rows.
///
/// Windows without an assigned group become compact three-column BED rows, while
/// grouped windows receive a fourth column so the metadata survives the write.
///
/// The function iterates the slice in order, writing either
/// `chrom<sep>start<sep>end` or `chrom<sep>start<sep>end<sep>group`. The caller
/// provides the separator (typically `\t`).
///
/// # Parameters
/// - `writer`: destination implementing [`Write`].
/// - `windows`: window slice to serialize.
/// - `separator`: delimiter to place between columns.
/// - `out_labels`: label keys to write after coordinates.
/// - `label_schema`: resolved label compositions.
///
/// # Returns
/// `Ok(())` on success or an error if writing fails.
pub fn write_windows<W: Write>(
    writer: &mut W,
    windows: &[Window],
    separator: char,
    out_labels: &[LabelKey],
    label_schema: &LabelSchema,
) -> Result<()> {
    for w in windows {
        let tuple_compositions = if label_schema.compositions().is_empty() {
            Vec::new()
        } else {
            build_tuple_compositions(&w.label_tuples, label_schema)
        };
        write!(
            writer,
            "{}{sep}{}{sep}{}",
            w.chrom.as_ref(),
            w.resized_start,
            w.resized_end,
            sep = separator
        )?;
        for key in out_labels {
            let label =
                render_label_for_key(&w.label_tuples, &tuple_compositions, key, label_schema);
            write!(writer, "{sep}{}", label, sep = separator)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Ensure a chromosome-specific temp writer exists and return it.
///
/// The helper creates the writer lazily the first time a chromosome appears and
/// reuses it for subsequent chunks.
///
/// It sanitizes the chromosome name, opens the file inside `temp_dir`, and wraps
/// it in a 1 MiB buffer to keep I/O efficient.
///
/// # Parameters
/// - `chrom`: chromosome identifier (e.g., `chr7`).
/// - `temp_dir`: per-run temporary directory.
/// - `temp_writers`: map of already-open writers.
///
/// # Returns
/// Mutable reference to the `ChromTempWriter` for `chrom`.
pub fn ensure_temp_writer_for_chrom<'a>(
    chrom: &str,
    temp_dir: &Path,
    temp_writers: &'a mut FxHashMap<String, ChromTempWriter>,
) -> Result<&'a mut ChromTempWriter> {
    if !temp_writers.contains_key(chrom) {
        let sanitized = chrom.replace('/', "_");
        // Temp file name is "chrom.<sanitized>.bed.tmp", for example "chrom.chr1.bed.tmp" or "chrom.chr1_KI270706v1_random.bed.tmp"
        let path = temp_dir.join(format!("chrom.{sanitized}.bed.tmp"));
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
pub fn finalize_temp_writers(
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

/// Concatenate temp outputs into the final writer.
///
/// Replays each chromosome's temp file back to back.
///
/// Output is buffered (stdout or file). Each temp file is streamed line by line,
/// preserving columns exactly as they were written, and temp entries are processed
/// in lexicographic chromosome order.
///
/// # Parameters
/// - `cfg`: resolved configuration.
/// - `temp_entries`: `(chromosome, temp_path)` pairs returned by
///   [`finalize_temp_writers`].
///
/// # Returns
/// `Ok(())` on success or an error if reading/writing fails.
pub fn concatenate_temps(cfg: &PrepareConfig, temp_entries: &[(String, PathBuf)]) -> Result<()> {
    let mut out: TextWriter = if cfg.output.as_os_str() == "-" {
        stdout_text_writer()
    } else {
        create_text_writer(&cfg.output)?
    };

    // Concatenate in lexicographic chrom order by default
    let mut entries: Vec<&(String, PathBuf)> = temp_entries.iter().collect();
    entries.sort_by(|a, b| a.0.cmp(&b.0));

    for (_chrom, path) in entries {
        let file = File::open(path)?;
        let reader = BufReader::with_capacity(1 << 20, file);

        for line_res in reader.lines() {
            let line = line_res?;
            if line.is_empty() {
                continue;
            }
            writeln!(out, "{}", line)?;
        }
    }

    out.finish()?;
    Ok(())
}
