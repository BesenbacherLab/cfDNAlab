use anyhow::{Context, Result};
use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::Path,
};
use zstd::Decoder as ZstdDecoder;
use zstd::Encoder as ZstdEncoder;

const BUF_CAP: usize = 1 << 20;
const DEFAULT_ZSTD_LEVEL: i32 = 3;

/// Join dot-separated name segments while skipping empty parts.
///
/// This keeps output naming consistent across commands when the optional output
/// prefix is omitted. For example, `["sample", "length_counts.npy"]` becomes
/// `sample.length_counts.npy`, while `["", "length_counts.npy"]` becomes
/// `length_counts.npy`.
pub fn dot_join(parts: &[&str]) -> String {
    parts
        .iter()
        .map(|part| part.trim())
        .filter(|part| !part.is_empty())
        .collect::<Vec<_>>()
        .join(".")
}

/// Open a text reader that transparently handles `.gz`, `.bgz`, `.zst`, or plain files.
///
/// The caller is responsible for handling "-" or stdin separately.
pub fn open_text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase());

    match ext.as_deref() {
        Some("gz") | Some("bgz") => {
            let file = File::open(path).with_context(|| format!("Opening {}", path.display()))?;
            let decoder = MultiGzDecoder::new(file);
            Ok(Box::new(BufReader::with_capacity(BUF_CAP, decoder)))
        }
        Some("zst") | Some("zstd") => {
            let file = File::open(path).with_context(|| format!("Opening {}", path.display()))?;
            let decoder = ZstdDecoder::new(file)
                .with_context(|| format!("Opening zstd decoder for {}", path.display()))?;
            Ok(Box::new(BufReader::with_capacity(BUF_CAP, decoder)))
        }
        _ => {
            let file = File::open(path).with_context(|| format!("Opening {}", path.display()))?;
            Ok(Box::new(BufReader::with_capacity(BUF_CAP, file)))
        }
    }
}

enum WriterInner {
    Stdout(BufWriter<io::Stdout>),
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
    Zstd(BufWriter<Box<dyn Write>>),
}

/// Writer that finishes compression streams when dropped via [`finish`](TextWriter::finish).
pub struct TextWriter {
    inner: WriterInner,
}

impl TextWriter {
    fn new(inner: WriterInner) -> Self {
        Self { inner }
    }

    /// Finalize the underlying stream and flush any buffered bytes.
    pub fn finish(self) -> Result<()> {
        match self.inner {
            WriterInner::Stdout(mut w) => {
                w.flush()?;
                Ok(())
            }
            WriterInner::Plain(mut w) => {
                w.flush()?;
                Ok(())
            }
            WriterInner::Gzip(mut enc) => {
                enc.flush()?;
                enc.try_finish()?;
                Ok(())
            }
            WriterInner::Zstd(mut w) => {
                w.flush()?;
                Ok(())
            }
        }
    }
}

impl Write for TextWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.inner {
            WriterInner::Stdout(w) => w.write(buf),
            WriterInner::Plain(w) => w.write(buf),
            WriterInner::Gzip(w) => w.write(buf),
            WriterInner::Zstd(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.inner {
            WriterInner::Stdout(w) => w.flush(),
            WriterInner::Plain(w) => w.flush(),
            WriterInner::Gzip(w) => w.flush(),
            WriterInner::Zstd(w) => w.flush(),
        }
    }
}

/// Construct a writer suitable for stdout.
pub fn stdout_text_writer() -> TextWriter {
    TextWriter::new(WriterInner::Stdout(BufWriter::with_capacity(
        BUF_CAP,
        io::stdout(),
    )))
}

/// Create a writer that compresses according to the file extension (`.gz`, `.zst`).
pub fn create_text_writer(path: &Path) -> Result<TextWriter> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase());

    match ext.as_deref() {
        Some("gz") | Some("bgz") => {
            let file =
                File::create(path).with_context(|| format!("Creating {}", path.display()))?;
            let buf = BufWriter::with_capacity(BUF_CAP, file);
            Ok(TextWriter::new(WriterInner::Gzip(GzEncoder::new(
                buf,
                Compression::default(),
            ))))
        }
        Some("zst") | Some("zstd") => {
            let file =
                File::create(path).with_context(|| format!("Creating {}", path.display()))?;
            let encoder = ZstdEncoder::new(file, DEFAULT_ZSTD_LEVEL)?;
            let sink: Box<dyn Write> = Box::new(encoder.auto_finish());
            Ok(TextWriter::new(WriterInner::Zstd(
                BufWriter::with_capacity(BUF_CAP, sink),
            )))
        }
        _ => {
            let file =
                File::create(path).with_context(|| format!("Creating {}", path.display()))?;
            Ok(TextWriter::new(WriterInner::Plain(
                BufWriter::with_capacity(BUF_CAP, file),
            )))
        }
    }
}
