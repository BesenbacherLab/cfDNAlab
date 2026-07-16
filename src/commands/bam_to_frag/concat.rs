use anyhow::{Context, Result};
use flate2::{Compression, write::GzEncoder};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use zstd::stream::read::Decoder;

/// Concatenate *.frag.tsv.zst files into a single gzip-compressed *.frag.tsv.gz, streaming.
///
/// Parameters
/// ----------
/// - input_paths:
///   Paths to input .frag.tsv.zst files, in the desired order.
/// - output_path:
///   Path to the resulting .frag.tsv.gz file to write.
/// - has_header:
///   If true, write the first file entirely and drop the first line (header) from subsequent files.
///
/// Returns
/// -------
/// - result:
///   Ok on success, or an error explaining what failed.
pub(crate) fn concat_frag_zst_to_gzip(
    input_paths: &[PathBuf],
    output_path: &Path,
    has_header: bool,
) -> Result<()> {
    ensure_nonempty(input_paths)?;

    // Open output encoder with buffering
    let out_file = File::create(output_path)
        .with_context(|| format!("Cannot create output: {}", output_path.display()))?;
    let out_buf = BufWriter::new(out_file);
    let mut gz = GzEncoder::new(out_buf, Compression::default());

    for (idx, in_path) in input_paths.iter().enumerate() {
        // Open and stream-decode zstd
        let in_file = File::open(in_path)
            .with_context(|| format!("Cannot open input: {}", in_path.display()))?;
        let in_buf = BufReader::new(in_file);
        let zstd_dec = Decoder::new(in_buf)
            .with_context(|| format!("Cannot zstd-decode: {}", in_path.display()))?;

        if has_header && idx > 0 {
            // Skip first line only for subsequent files
            let mut lines_reader = BufReader::new(zstd_dec);
            skip_first_line(&mut lines_reader)
                .with_context(|| format!("Failed skipping header in: {}", in_path.display()))?;
            io::copy(&mut lines_reader, &mut gz).with_context(|| {
                format!(
                    "Stream copy failed after header skip: {}",
                    in_path.display()
                )
            })?;
        } else {
            // Copy raw decoded bytes
            let mut decoded = BufReader::new(zstd_dec);
            io::copy(&mut decoded, &mut gz)
                .with_context(|| format!("Stream copy failed for: {}", in_path.display()))?;
        }
    }

    // Finish gzip stream
    gz.try_finish().context("Failed to finish gzip stream")?;
    Ok(())
}

fn ensure_nonempty(input_paths: &[PathBuf]) -> Result<()> {
    if input_paths.is_empty() {
        anyhow::bail!("No input files provided")
    }
    Ok(())
}

fn skip_first_line<R: BufRead>(r: &mut R) -> io::Result<()> {
    let mut _discard = String::new();
    r.read_line(&mut _discard).map(|_| ())
}
