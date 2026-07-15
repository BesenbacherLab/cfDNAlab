use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
#[cfg(writes_text_outputs)]
use flate2::{Compression, write::GzEncoder};
use fxhash::FxHashMap;
#[cfg(writes_text_outputs)]
use std::io::{BufWriter, Write};
use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read},
    path::{Path, PathBuf},
    sync::mpsc::{Receiver, SyncSender, TryRecvError, sync_channel},
    thread::{self, JoinHandle},
};
use zstd::Decoder as ZstdDecoder;
#[cfg(writes_text_outputs)]
use zstd::Encoder as ZstdEncoder;

const BUF_CAP: usize = 1 << 20;
const BACKGROUND_READ_CHUNK_SIZE: usize = 4 << 20;
const BACKGROUND_READ_QUEUE_SIZE: usize = 2;
#[cfg(writes_text_outputs)]
const DEFAULT_ZSTD_LEVEL: i32 = 3;
const REPLACEABLE_DIRECTORY_EXTENSIONS: &[&str] = &["zarr"];

/// Join dot-separated name segments while skipping empty parts.
///
/// This keeps output naming consistent across commands when the optional output
/// prefix is omitted. For example, `["sample", "length_counts.npy"]` becomes
/// `sample.length_counts.npy`, while `["", "length_counts.npy"]` becomes
/// `length_counts.npy`.
pub(crate) fn dot_join(parts: &[&str]) -> String {
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
pub(crate) fn open_text_reader(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase());
    let file = File::open(path).with_context(|| format!("Opening {}", path.display()))?;

    #[cfg(target_os = "linux")]
    advise_sequential_access(&file, path);

    match ext.as_deref() {
        Some("gz") | Some("bgz") => {
            let decoder = MultiGzDecoder::new(file);
            Ok(Box::new(BufReader::with_capacity(BUF_CAP, decoder)))
        }
        Some("zst") | Some("zstd") => {
            let decoder = ZstdDecoder::new(file)
                .with_context(|| format!("Opening zstd decoder for {}", path.display()))?;
            Ok(Box::new(BufReader::with_capacity(BUF_CAP, decoder)))
        }
        _ => Ok(Box::new(BufReader::with_capacity(BUF_CAP, file))),
    }
}

/// Open a text reader that loads and decompresses bytes on a background thread.
///
/// The returned reader presents the same sequential byte stream as [`open_text_reader`]. A small
/// bounded queue lets file reading and decompression overlap with parsing without allowing the
/// reader to buffer the full input in memory.
pub(crate) fn open_text_reader_in_background(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let source = open_text_reader(path)?;
    Ok(Box::new(BackgroundTextReader::new(source)?))
}

enum BackgroundReadMessage {
    Data(Vec<u8>),
    Error(io::Error),
    End,
}

/// Sequential reader backed by a bounded queue filled from another thread.
struct BackgroundTextReader {
    messages: Option<Receiver<BackgroundReadMessage>>,
    recycled_buffers: SyncSender<Vec<u8>>,
    current_buffer: Vec<u8>,
    current_position: usize,
    reached_end: bool,
    worker: Option<JoinHandle<()>>,
}

impl BackgroundTextReader {
    fn new(mut source: Box<dyn BufRead + Send>) -> Result<Self> {
        let (message_sender, message_receiver) = sync_channel(BACKGROUND_READ_QUEUE_SIZE);
        let (recycled_buffer_sender, recycled_buffer_receiver) =
            sync_channel(BACKGROUND_READ_QUEUE_SIZE);
        let worker = thread::Builder::new()
            .name("cfdnalab-text-reader".to_string())
            .spawn(move || {
                loop {
                    let mut buffer = match recycled_buffer_receiver.try_recv() {
                        Ok(buffer) => buffer,
                        Err(TryRecvError::Empty) => vec![0; BACKGROUND_READ_CHUNK_SIZE],
                        Err(TryRecvError::Disconnected) => return,
                    };
                    buffer.resize(BACKGROUND_READ_CHUNK_SIZE, 0);

                    match source.read(&mut buffer) {
                        Ok(0) => {
                            let _ = message_sender.send(BackgroundReadMessage::End);
                            return;
                        }
                        Ok(bytes_read) => {
                            buffer.truncate(bytes_read);
                            if message_sender
                                .send(BackgroundReadMessage::Data(buffer))
                                .is_err()
                            {
                                return;
                            }
                        }
                        Err(error) => {
                            let _ = message_sender.send(BackgroundReadMessage::Error(error));
                            return;
                        }
                    }
                }
            })
            .context("starting background text reader")?;

        Ok(Self {
            messages: Some(message_receiver),
            recycled_buffers: recycled_buffer_sender,
            current_buffer: Vec::new(),
            current_position: 0,
            reached_end: false,
            worker: Some(worker),
        })
    }

    fn receive_next_buffer(&mut self) -> io::Result<bool> {
        if !self.current_buffer.is_empty() {
            let previous_buffer = std::mem::take(&mut self.current_buffer);
            let _ = self.recycled_buffers.try_send(previous_buffer);
            self.current_position = 0;
        }

        let message = self
            .messages
            .as_ref()
            .expect("background reader message receiver missing")
            .recv()
            .map_err(|_| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "background text reader stopped before reporting end of input",
                )
            })?;
        match message {
            BackgroundReadMessage::Data(buffer) => {
                self.current_buffer = buffer;
                Ok(true)
            }
            BackgroundReadMessage::Error(error) => {
                self.reached_end = true;
                Err(error)
            }
            BackgroundReadMessage::End => {
                self.reached_end = true;
                Ok(false)
            }
        }
    }
}

impl Read for BackgroundTextReader {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        if output.is_empty() {
            return Ok(0);
        }

        let bytes_to_copy = {
            let available = self.fill_buf()?;
            let bytes_to_copy = available.len().min(output.len());
            output[..bytes_to_copy].copy_from_slice(&available[..bytes_to_copy]);
            bytes_to_copy
        };
        self.consume(bytes_to_copy);
        Ok(bytes_to_copy)
    }
}

impl BufRead for BackgroundTextReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if !self.reached_end
            && self.current_position == self.current_buffer.len()
            && !self.receive_next_buffer()?
        {
            return Ok(&[]);
        }
        Ok(&self.current_buffer[self.current_position..])
    }

    fn consume(&mut self, amount: usize) {
        self.current_position = self
            .current_position
            .saturating_add(amount)
            .min(self.current_buffer.len());
    }
}

impl Drop for BackgroundTextReader {
    fn drop(&mut self) {
        // Disconnect the queue before joining so a producer blocked on a full queue can stop
        self.messages.take();
        if let Some(worker) = self.worker.take()
            && worker.join().is_err()
            && !thread::panicking()
        {
            tracing::warn!("Background text reader thread panicked");
        }
    }
}

/// Ask Linux to use a larger readahead window for a sequential text stream.
///
/// This is only a performance hint. A failure should not prevent the caller from reading an
/// otherwise valid file, but it is reported so an unsupported or unexpected platform setup is
/// visible to the user.
#[cfg(target_os = "linux")]
fn advise_sequential_access(file: &File, path: &Path) {
    use std::os::fd::AsRawFd;

    // posix_fadvise returns the error number directly instead of setting errno
    let error_code =
        unsafe { libc::posix_fadvise(file.as_raw_fd(), 0, 0, libc::POSIX_FADV_SEQUENTIAL) };
    if error_code != 0 {
        let error = std::io::Error::from_raw_os_error(error_code);
        tracing::warn!(
            "Could not enable sequential file readahead for {}: {}",
            path.display(),
            error
        );
    }
}

#[cfg(writes_text_outputs)]
enum WriterInner {
    #[cfg(feature = "cmd_prepare_windows")]
    Stdout(BufWriter<io::Stdout>),
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
    Zstd(BufWriter<Box<dyn Write>>),
}

/// Writer that finishes compression streams when dropped via [`finish`](TextWriter::finish).
#[cfg(writes_text_outputs)]
pub(crate) struct TextWriter {
    inner: WriterInner,
}

#[cfg(writes_text_outputs)]
impl TextWriter {
    fn new(inner: WriterInner) -> Self {
        Self { inner }
    }

    /// Finalize the underlying stream and flush any buffered bytes.
    pub(crate) fn finish(self) -> Result<()> {
        match self.inner {
            #[cfg(feature = "cmd_prepare_windows")]
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

#[cfg(writes_text_outputs)]
impl Write for TextWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.inner {
            #[cfg(feature = "cmd_prepare_windows")]
            WriterInner::Stdout(w) => w.write(buf),
            WriterInner::Plain(w) => w.write(buf),
            WriterInner::Gzip(w) => w.write(buf),
            WriterInner::Zstd(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.inner {
            #[cfg(feature = "cmd_prepare_windows")]
            WriterInner::Stdout(w) => w.flush(),
            WriterInner::Plain(w) => w.flush(),
            WriterInner::Gzip(w) => w.flush(),
            WriterInner::Zstd(w) => w.flush(),
        }
    }
}

/// Construct a writer suitable for stdout.
#[cfg(feature = "cmd_prepare_windows")]
pub(crate) fn stdout_text_writer() -> TextWriter {
    TextWriter::new(WriterInner::Stdout(BufWriter::with_capacity(
        BUF_CAP,
        io::stdout(),
    )))
}

/// Create a writer that compresses according to the file extension (`.gz`, `.zst`).
#[cfg(writes_text_outputs)]
pub(crate) fn create_text_writer(path: &Path) -> Result<TextWriter> {
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

/* Helpers for writing final outputs to temp files before moving them into place.
 */

/// Tracks final output files written through a command temp directory.
pub(crate) struct FinalOutputFiles {
    temp_dir: PathBuf,
    paths: FxHashMap<PathBuf, PathBuf>,
}

impl FinalOutputFiles {
    pub(crate) fn new(temp_dir: &Path) -> Result<Self> {
        Ok(Self {
            temp_dir: create_final_output_temp_dir(temp_dir)?,
            paths: FxHashMap::default(),
        })
    }

    pub(crate) fn temp_path_for(&self, final_path: &Path) -> Result<PathBuf> {
        output_path_in_dir(&self.temp_dir, final_path)
    }

    #[cfg(any(
        feature = "cmd_ends",
        feature = "cmd_fcoverage",
        feature = "cmd_fragment_kmers",
        feature = "cmd_gc_bias",
        feature = "cmd_wps"
    ))]
    pub(crate) fn temp_dir(&self) -> &Path {
        &self.temp_dir
    }

    pub(crate) fn record(&mut self, temp_path: PathBuf, final_path: PathBuf) -> Result<()> {
        if self.paths.contains_key(&temp_path) {
            anyhow::bail!("duplicate temporary output path: {}", temp_path.display());
        }
        if self.paths.values().any(|path| path == &final_path) {
            anyhow::bail!("duplicate final output path: {}", final_path.display());
        }

        self.paths.insert(temp_path, final_path);
        Ok(())
    }

    /// Record files written in the final-output temp directory under their final filenames.
    ///
    /// Each temp path is mapped to `output_dir/<file name>`. Use this after passing
    /// `FinalOutputFiles::temp_dir` to a writer that returns the files it created.
    #[cfg(any(
        all(feature = "cmd_gc_bias", feature = "plotters"),
        feature = "cmd_fragment_kmers"
    ))]
    pub(crate) fn record_temp_files_with_same_names_in(
        &mut self,
        temp_paths: impl IntoIterator<Item = PathBuf>,
        output_dir: &Path,
    ) -> Result<()> {
        for temp_path in temp_paths {
            let file_name = temp_path.file_name().with_context(|| {
                format!(
                    "temporary output path has no filename: {}",
                    temp_path.display()
                )
            })?;
            let final_path = output_dir.join(file_name);
            self.record(temp_path, final_path)?;
        }
        Ok(())
    }

    /// Move recorded temp files to their final paths one by one.
    ///
    /// Each recorded artifact has already been fully written in the temp directory, so callers do
    /// not expose half-written individual files. Directory outputs replace an existing directory
    /// at the final path before the rename only for explicitly supported directory-backed formats,
    /// such as Zarr. This helper is still best-effort across multiple artifacts: if a later move
    /// fails, earlier artifacts may already be visible at their final paths and the error is
    /// returned to the caller.
    pub(crate) fn move_into_place(self) -> Result<()> {
        for (temp_path, final_path) in self.paths {
            move_output_file_into_place(&temp_path, &final_path)?;
        }
        Ok(())
    }
}

/// Create the temp subdirectory used for final output files.
fn create_final_output_temp_dir(temp_dir: &Path) -> Result<PathBuf> {
    let final_output_temp_dir = temp_dir.join("final_outputs");
    std::fs::create_dir_all(&final_output_temp_dir).with_context(|| {
        format!(
            "creating final output temp directory {}",
            final_output_temp_dir.display()
        )
    })?;
    Ok(final_output_temp_dir)
}

/// Return the final output filename inside a directory.
///
/// The caller should pass a directory under the output directory, so the final rename does not copy
/// across filesystems.
fn output_path_in_dir(directory: &Path, final_path: &Path) -> Result<PathBuf> {
    let file_name = final_path
        .file_name()
        .with_context(|| format!("output path has no filename: {}", final_path.display()))?
        .to_string_lossy();
    Ok(directory.join(file_name.as_ref()))
}

/// Move a fully written output artifact into place.
///
/// File outputs rely on platform rename behavior. Directory outputs cannot be renamed over an
/// existing directory, so an existing final directory is removed first after validating that the
/// path has an explicitly supported directory-backed output extension. Replacing a file with a
/// directory is treated as an error because that usually means the output contract changed or the
/// destination path is wrong.
fn move_output_file_into_place(temp_path: &Path, final_path: &Path) -> Result<()> {
    if temp_path.is_dir() && final_path.exists() {
        if final_path.is_dir() {
            ensure_safe_directory_replacement_path(final_path)?;
            std::fs::remove_dir_all(final_path).with_context(|| {
                format!(
                    "removing previous output directory {}",
                    final_path.display()
                )
            })?;
        } else {
            anyhow::bail!(
                "cannot replace output file {} with output directory {}",
                final_path.display(),
                temp_path.display()
            );
        }
    }
    std::fs::rename(temp_path, final_path).with_context(|| {
        format!(
            "moving output file {} to {}",
            temp_path.display(),
            final_path.display()
        )
    })
}

fn ensure_safe_directory_replacement_path(final_path: &Path) -> Result<()> {
    anyhow::ensure!(
        final_path.file_name().is_some(),
        "refusing to replace output directory {} because it has no final path component",
        final_path.display()
    );
    anyhow::ensure!(
        has_replaceable_directory_extension(final_path),
        "refusing to replace output directory {} because its extension is not one of: {}",
        final_path.display(),
        REPLACEABLE_DIRECTORY_EXTENSIONS.join(", ")
    );
    let canonical_path = final_path
        .canonicalize()
        .with_context(|| format!("canonicalizing output directory {}", final_path.display()))?;
    anyhow::ensure!(
        canonical_path.file_name().is_some(),
        "refusing to replace output directory {} because it resolves to {}",
        final_path.display(),
        canonical_path.display()
    );
    Ok(())
}

fn has_replaceable_directory_extension(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .map(|extension| {
            REPLACEABLE_DIRECTORY_EXTENSIONS
                .iter()
                .any(|allowed| extension.eq_ignore_ascii_case(allowed))
        })
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    include!("io_tests.rs");
}
