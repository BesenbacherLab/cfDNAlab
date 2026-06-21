use std::path::{Path, PathBuf};

/// Controls reporting side effects for a programmatic command run.
///
/// These options do not change the scientific computation. They only decide whether the command
/// prints statistics, shows progress bars, and writes status messages while it runs.
///
/// Programmatic runners do not initialize cfDNAlab's CLI logging runtime. If an application wants
/// to collect cfDNAlab status messages, it should install its own `tracing` subscriber once for
/// the process and enable the relevant fields here. Use `RunOptions::new_quiet()` for library-style
/// calls that should not write progress bars, status messages, equivalent CLI commands, or final
/// statistics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RunOptions {
    /// Print the final command statistics after a successful run.
    pub report_statistics: bool,
    /// Show progress bars for long-running tiled work.
    pub show_progress: bool,
    /// Write status messages such as input loading and output writing progress.
    pub log_statuses: bool,
    /// Write the equivalent full CLI command before the command starts heavy work.
    pub log_equivalent_cli: bool,
}

impl RunOptions {
    /// Build a reporting profile with all reporting side effects enabled.
    ///
    /// This enables statistics, progress bars, status messages, and equivalent CLI logging.
    /// In direct Rust calls, status messages go through the caller's `tracing` subscriber if one
    /// is installed, while statistics are written to the primary output sink.
    pub fn new_cli() -> Self {
        Self {
            report_statistics: true,
            show_progress: true,
            log_statuses: true,
            log_equivalent_cli: true,
        }
    }

    /// Build a reporting profile with all reporting side effects disabled.
    ///
    /// This preserves the same command computation and output files while suppressing cfDNAlab's
    /// own progress bars, status messages, equivalent CLI command, and final statistics.
    pub fn new_quiet() -> Self {
        Self {
            report_statistics: false,
            show_progress: false,
            log_statuses: false,
            log_equivalent_cli: false,
        }
    }
}

impl Default for RunOptions {
    fn default() -> Self {
        Self::new_cli()
    }
}

/// Writes an info-level status event when command status logging is enabled.
///
/// Command runners use this for user-facing progress milestones such as loading inputs,
/// reducing temporary files, and writing outputs. This keeps the `RunOptions` check near the
/// status event while preserving the normal `tracing::info!` syntax, including `target:` and
/// formatting arguments.
///
/// Formatting arguments stay inside the guard, so they are only evaluated when
/// `options.log_statuses` is true.
///
/// Parameters
/// ----------
/// - `$options`:
///   A `RunOptions` value or expression with a `log_statuses` field.
/// - `$arg`:
///   The arguments forwarded to `tracing::info!` when status logging is enabled.
#[cfg(any(
    feature = "cmd_bam_to_bam",
    feature = "cmd_bam_to_frag",
    feature = "cmd_coverage_weights",
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_gc_bias",
    feature = "cmd_lengths",
    feature = "cmd_midpoints",
    feature = "cmd_wps",
    feature = "cmd_wps_peaks",
))]
macro_rules! status_info {
    ($options:expr, $($arg:tt)+) => {{
        if ($options).log_statuses {
            tracing::info!($($arg)+);
        }
    }};
}

#[cfg(any(
    feature = "cmd_bam_to_bam",
    feature = "cmd_bam_to_frag",
    feature = "cmd_coverage_weights",
    feature = "cmd_ends",
    feature = "cmd_fcoverage",
    feature = "cmd_fragment_kmers",
    feature = "cmd_gc_bias",
    feature = "cmd_lengths",
    feature = "cmd_midpoints",
    feature = "cmd_wps",
    feature = "cmd_wps_peaks",
))]
pub(crate) use status_info;

/// Format the equivalent CLI command as a status log message.
///
/// The tracing formatter already prefixes the event with the command target,
/// so keep the command on the same line as the label for predictable assertions
/// and simple line-oriented logs.
pub(crate) fn equivalent_cli_log_message(command: &str) -> String {
    format!("Equivalent CLI: {command}")
}

/// Shared interface for command-specific run results.
///
/// Each command returns its own result type because output files and counters differ by command.
/// This trait exposes the common parts so callers can write generic code for result inspection
/// without losing access to command-specific fields.
pub trait CommandRunResult {
    /// Counter type reported by the command.
    type Counters;

    /// Return the counters collected during the run.
    fn counters(&self) -> &Self::Counters;
    /// Return all final output files that the command produced.
    fn output_files(&self) -> &[PathBuf];
    /// Return the main output file when the command has a single primary artifact.
    fn primary_output(&self) -> Option<&Path>;
}
