use crate::commands::prepare_windows::{config::PrepareConfig, prepare_windows::FinalWindow};
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
/// lets the pipeline defer global policies (such as `min_per_group`) until the
/// final concatenation pass.
///
/// The writer wraps a 1 MiB `BufWriter<File>` pointing at a sanitized path under
/// the run's temp directory. Callers borrow the writer when emitting rows and
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
/// The function iterates the slice in order, emitting either
/// `chrom<sep>start<sep>end` or `chrom<sep>start<sep>end<sep>group`. The caller
/// provides the separator (typically `\t`).
///
/// # Parameters
/// - `writer`: destination implementing [`Write`].
/// - `windows`: window slice to serialize.
/// - `separator`: delimiter to place between columns.
///
/// # Returns
/// `Ok(())` on success or an error if writing fails.
pub fn write_windows<W: Write>(
    writer: &mut W,
    windows: &[FinalWindow],
    separator: char,
) -> Result<()> {
    for w in windows {
        if w.group.is_empty() {
            writeln!(
                writer,
                "{}{sep}{}{sep}{}",
                w.chrom,
                w.start,
                w.end,
                sep = separator
            )?;
        } else {
            writeln!(
                writer,
                "{}{sep}{}{sep}{}{sep}{}",
                w.chrom,
                w.start,
                w.end,
                w.group,
                sep = separator
            )?;
        }
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

/// Concatenate temp outputs into the final writer, honoring `min_per_group`.
///
/// Replays each chromosome's temp file back to back, filtering out group labels
/// that fail the minimum-count threshold when requested.
///
/// Output is buffered (stdout or file). Each temp file is streamed line by line;
/// the optional group column determines whether to emit a three- or four-column
/// row, and temp entries are processed in lexicographic chromosome order.
///
/// # Parameters
/// - `cfg`: resolved configuration.
/// - `temp_entries`: `(chromosome, temp_path)` pairs returned by
///   [`finalize_temp_writers`].
/// - `global_group_counts`: counts per group after spacing/merging.
///
/// # Returns
/// `Ok(())` on success or an error if reading/writing fails.
pub fn concatenate_temps_enforcing_min_per_group(
    cfg: &PrepareConfig,
    temp_entries: &[(String, PathBuf)],
    global_group_counts: &FxHashMap<String, u64>,
) -> Result<()> {
    // Open output
    let mut out: Box<dyn Write> = if let Some(path) = &cfg.output {
        if path.as_os_str() == "-" {
            Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout()))
        } else {
            Box::new(BufWriter::with_capacity(
                1 << 20,
                File::create(path).with_context(|| format!("Creating output {:?}", path))?,
            ))
        }
    } else {
        Box::new(BufWriter::with_capacity(1 << 20, std::io::stdout()))
    };

    // Helper to decide if a group passes the threshold
    let min_required = cfg.min_per_group;

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
            if cfg.group_cols.is_empty() && cfg.distance_bins.is_none() {
                // No group anywhere, write as-is (3 columns)
                writeln!(out, "{}", line)?;
                continue;
            }

            // If group column exists, it is the 4th field
            // We do not add a header; this is a minimal BED-like.
            let mut parts = line.split(cfg.separator);
            let chrom_name = parts.next().unwrap_or_default();
            let start = parts.next().unwrap_or_default();
            let end = parts.next().unwrap_or_default();
            let group = parts.next().unwrap_or_default();

            if group.is_empty() {
                // If there is no group after all processing, write the 3 cols only
                writeln!(
                    out,
                    "{}{}{}{}{}",
                    chrom_name, cfg.separator, start, cfg.separator, end
                )?;
            } else if let Some(min_n) = min_required {
                let count = *global_group_counts.get(group).unwrap_or(&0);
                if count >= min_n as u64 {
                    writeln!(
                        out,
                        "{}{sep}{}{sep}{}{sep}{}",
                        chrom_name,
                        start,
                        end,
                        group,
                        sep = cfg.separator
                    )?;
                }
            } else {
                writeln!(
                    out,
                    "{}{sep}{}{sep}{}{sep}{}",
                    chrom_name,
                    start,
                    end,
                    group,
                    sep = cfg.separator
                )?;
            }
        }
    }

    Ok(())
}
