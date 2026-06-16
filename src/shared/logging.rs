use anyhow::{Context, Result};
use chrono::Local;
use rand::{Rng, distr::Alphanumeric};
use std::{
    fs::{self, File, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
    sync::{Arc, Mutex, OnceLock},
};
use tracing::Level;
use tracing_subscriber::{Layer, filter::filter_fn, fmt, layer::SubscriberExt};

/// Shared logging argument used by commands that opt into tracing-based CLI output.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LoggingArgs {
    /// Logging destination `[stdout|quiet|file|file=<path>]`
    ///
    /// `stdout` keeps the normal run narrative on standard output.
    ///
    /// `quiet` suppresses the normal run narrative and progress bars, while warnings
    /// and errors still go to `stderr`.
    ///
    /// `file` writes the normal run narrative to an auto-generated log file under
    /// the command output directory.
    ///
    /// `file=<path>` writes the normal run narrative to the exact path you provide.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "stdout",
            value_parser = parse_log_spec,
            help_heading = "Logging"
        )
    )]
    pub log: LogSpec,
}

/// Parsed logging mode for a top-level command.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum LogSpec {
    #[default]
    Stdout,
    Quiet,
    File(Option<PathBuf>),
}

/// Parse the compact `--log` grammar shared by tracing-enabled commands.
pub fn parse_log_spec(value: &str) -> Result<LogSpec, String> {
    match value {
        "stdout" => Ok(LogSpec::Stdout),
        "quiet" => Ok(LogSpec::Quiet),
        "file" => Ok(LogSpec::File(None)),
        _ => {
            if let Some(path) = value.strip_prefix("file=") {
                if path.is_empty() {
                    return Err(
                        "invalid --log value 'file=': expected a non-empty path after 'file='"
                            .to_string(),
                    );
                }
                Ok(LogSpec::File(Some(PathBuf::from(path))))
            } else {
                Err(format!(
                    "invalid --log value '{value}'. Expected one of: stdout, quiet, file, file=<path>"
                ))
            }
        }
    }
}

#[derive(Clone)]
enum PrimaryOutput {
    Stdout,
    Quiet,
    File(Arc<Mutex<File>>),
}

static PRIMARY_OUTPUT: OnceLock<PrimaryOutput> = OnceLock::new();

/// Initialize tracing and the shared primary output sink for one CLI invocation.
///
/// The primary sink carries the normal run narrative and explicit summary blocks.
/// Warnings and errors always stay on `stderr`.
pub fn init_cli_logging(
    command_name: &str,
    log_spec: &LogSpec,
    default_output_dir: Option<&Path>,
) -> Result<()> {
    let resolved_log_path = resolve_log_path(command_name, log_spec, default_output_dir)?;
    let primary_output = match resolved_log_path {
        Some(path) => {
            if let Some(parent) = path.parent() {
                fs::create_dir_all(parent)
                    .with_context(|| format!("creating log directory {}", parent.display()))?;
            }
            let file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&path)
                .with_context(|| format!("opening log file {}", path.display()))?;
            PrimaryOutput::File(Arc::new(Mutex::new(file)))
        }
        None => match log_spec {
            LogSpec::Stdout => PrimaryOutput::Stdout,
            LogSpec::Quiet => PrimaryOutput::Quiet,
            LogSpec::File(_) => unreachable!("file mode must resolve to a concrete path"),
        },
    };

    install_tracing(primary_output.clone()).context("initializing tracing subscriber")?;
    PRIMARY_OUTPUT
        .set(primary_output)
        .map_err(|_| anyhow::anyhow!("primary CLI output sink was already initialized"))?;
    Ok(())
}

/// Write a preformatted block to the primary sink without appending a newline.
pub fn write_primary(text: &str) {
    match PRIMARY_OUTPUT.get() {
        Some(PrimaryOutput::Stdout) | None => {
            let mut stdout = io::stdout().lock();
            stdout
                .write_all(text.as_bytes())
                .expect("writing primary CLI output");
            stdout.flush().expect("flushing primary CLI output");
        }
        Some(PrimaryOutput::Quiet) => {}
        Some(PrimaryOutput::File(file)) => {
            let mut file = file.lock().expect("locking primary log file");
            file.write_all(text.as_bytes())
                .expect("writing primary CLI output");
            file.flush().expect("flushing primary CLI output");
        }
    }
}

/// Write one logical line to the primary sink.
pub fn write_primary_line(line: &str) {
    write_primary(&format!("{line}\n"));
}

/// Duplicate a top-level stderr line into the log file when file logging is active.
pub fn duplicate_stderr_line_to_file(line: &str) {
    if let Some(PrimaryOutput::File(file)) = PRIMARY_OUTPUT.get() {
        let mut file = file.lock().expect("locking primary log file");
        writeln!(file, "{line}").expect("writing mirrored stderr line to log file");
        file.flush()
            .expect("flushing mirrored stderr line to log file");
    }
}

/// Return whether the current primary sink should use terminal-oriented formatting.
pub fn primary_uses_terminal_formatting() -> bool {
    matches!(PRIMARY_OUTPUT.get(), Some(PrimaryOutput::Stdout) | None)
}

/// Return whether progress reporting should remain enabled.
pub fn progress_enabled() -> bool {
    !matches!(PRIMARY_OUTPUT.get(), Some(PrimaryOutput::Quiet))
}

fn resolve_log_path(
    command_name: &str,
    log_spec: &LogSpec,
    default_output_dir: Option<&Path>,
) -> Result<Option<PathBuf>> {
    match log_spec {
        LogSpec::Stdout | LogSpec::Quiet => Ok(None),
        LogSpec::File(Some(path)) => Ok(Some(path.clone())),
        LogSpec::File(None) => {
            let base_dir = match default_output_dir {
                Some(path) => path.to_path_buf(),
                None => std::env::current_dir().context("reading current working directory")?,
            };
            let logs_dir = base_dir.join("logs");
            let timestamp = Local::now().format("%Y-%m-%d_%H-%M-%S");
            let suffix = random_suffix(8);
            Ok(Some(
                logs_dir.join(format!("{command_name}_{timestamp}_{suffix}.log")),
            ))
        }
    }
}

fn random_suffix(length: usize) -> String {
    rand::rng()
        .sample_iter(&Alphanumeric)
        .take(length)
        .map(char::from)
        .collect()
}

fn install_tracing(primary_output: PrimaryOutput) -> Result<()> {
    let stderr_layer = fmt::layer()
        .with_ansi(false)
        .without_time()
        .with_target(true)
        .with_level(true)
        .with_writer(io::stderr)
        .with_filter(filter_fn(|metadata| {
            matches!(*metadata.level(), Level::WARN | Level::ERROR)
        }));

    match primary_output {
        PrimaryOutput::Quiet => {
            let subscriber = tracing_subscriber::registry().with(stderr_layer);
            tracing::subscriber::set_global_default(subscriber)
                .context("setting global tracing subscriber")?;
        }
        PrimaryOutput::Stdout => {
            let primary_layer = fmt::layer()
                .with_ansi(false)
                .without_time()
                .with_target(true)
                .with_level(false)
                .with_writer(PrimaryMakeWriter::new(PrimaryOutput::Stdout))
                .with_filter(filter_fn(|metadata| {
                    matches!(*metadata.level(), Level::INFO)
                }));
            let subscriber = tracing_subscriber::registry()
                .with(stderr_layer)
                .with(primary_layer);
            tracing::subscriber::set_global_default(subscriber)
                .context("setting global tracing subscriber")?;
        }
        PrimaryOutput::File(file) => {
            let primary_layer = fmt::layer()
                .with_ansi(false)
                .without_time()
                .with_target(true)
                .with_level(true)
                .with_writer(PrimaryMakeWriter::new(PrimaryOutput::File(file)))
                .with_filter(filter_fn(|metadata| {
                    matches!(*metadata.level(), Level::INFO | Level::WARN | Level::ERROR)
                }));
            let subscriber = tracing_subscriber::registry()
                .with(stderr_layer)
                .with(primary_layer);
            tracing::subscriber::set_global_default(subscriber)
                .context("setting global tracing subscriber")?;
        }
    }

    Ok(())
}

#[derive(Clone)]
struct PrimaryMakeWriter {
    output: PrimaryOutput,
}

impl PrimaryMakeWriter {
    fn new(output: PrimaryOutput) -> Self {
        Self { output }
    }
}

impl<'a> tracing_subscriber::fmt::writer::MakeWriter<'a> for PrimaryMakeWriter {
    type Writer = PrimaryWriter;

    fn make_writer(&'a self) -> Self::Writer {
        match &self.output {
            PrimaryOutput::Stdout => PrimaryWriter::Stdout(io::stdout()),
            PrimaryOutput::Quiet => PrimaryWriter::Sink(io::sink()),
            PrimaryOutput::File(file) => PrimaryWriter::File(file.clone()),
        }
    }
}

enum PrimaryWriter {
    Stdout(io::Stdout),
    Sink(io::Sink),
    File(Arc<Mutex<File>>),
}

impl Write for PrimaryWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Stdout(stdout) => stdout.write(buf),
            Self::Sink(sink) => sink.write(buf),
            Self::File(file) => file
                .lock()
                .map_err(|_| io::Error::other("failed to lock log file"))?
                .write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Stdout(stdout) => stdout.flush(),
            Self::Sink(sink) => sink.flush(),
            Self::File(file) => file
                .lock()
                .map_err(|_| io::Error::other("failed to lock log file"))?
                .flush(),
        }
    }
}
