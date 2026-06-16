use std::path::{Path, PathBuf};

/// Controls reporting side effects for a programmatic command run.
///
/// These options do not change the scientific computation. They only decide whether the command
/// prints statistics, shows progress bars, and writes status messages while it runs.
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
    /// This enables statistics, progress bars, and status messages.
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
    /// This disables printing and progress side effects while preserving the same command
    /// computation and output files.
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
