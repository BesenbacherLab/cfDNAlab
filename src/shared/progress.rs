use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use std::io::{self, IsTerminal};
#[allow(dead_code)]
use std::time::Duration;

const DEFAULT_BAR_TEMPLATE: &str = "       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}";
#[allow(dead_code)]
const DEFAULT_SPINNER_TEMPLATE: &str = "{spinner} {msg} [{elapsed_precise}]";
#[allow(dead_code)]
const DEFAULT_SPINNER_TICK_INTERVAL: Duration = Duration::from_millis(100);

/// Creates terminal-aware progress bars and spinners for command runners.
///
/// This keeps progress reporting consistent across commands while avoiding
/// progress redraw noise in redirected logs and pipelines. The factory checks
/// whether `stderr` is a terminal because `indicatif` draws to `stderr` by
/// default.
///
/// Use [`ProgressFactory::with_enabled`] with the command's `RunOptions.show_progress`
/// value.
/// When drawing is disabled, the factory still returns normal `ProgressBar`
/// values, they simply use a hidden draw target.
///
/// Examples
/// --------
/// ```ignore
/// let progress = ProgressFactory::new();
/// let bar = progress.default_bar(total_tiles as u64);
/// bar.set_message("Counting tiles");
/// ```
///
/// ```ignore
/// let progress = ProgressFactory::with_enabled(options.show_progress);
/// let spinner = progress.default_spinner();
/// spinner.set_message("Reading streamed input");
/// ```
#[derive(Clone, Copy, Debug)]
pub(crate) struct ProgressFactory {
    enabled: bool,
}

impl Default for ProgressFactory {
    fn default() -> Self {
        Self::new()
    }
}

impl ProgressFactory {
    /// Create a progress factory that draws only when `stderr` is a terminal.
    pub(crate) fn new() -> Self {
        Self::with_enabled(true)
    }

    /// Create a progress factory that never draws progress output.
    #[allow(dead_code)]
    pub(crate) fn hidden() -> Self {
        Self { enabled: false }
    }

    /// Create a progress factory that honors the explicit progress setting and TTY state.
    pub(crate) fn with_enabled(enabled: bool) -> Self {
        Self {
            enabled: enabled && io::stderr().is_terminal(),
        }
    }

    /// Create a bar with the repo-wide default template.
    pub(crate) fn default_bar(&self, len: u64) -> ProgressBar {
        self.bar_with_style(len, Self::default_bar_style())
    }

    /// Create a spinner with the repo-wide default template and tick interval.
    #[allow(dead_code)]
    pub(crate) fn default_spinner(&self) -> ProgressBar {
        self.spinner_with_style(Self::default_spinner_style())
    }

    /// Create a bar with a caller-provided style.
    pub(crate) fn bar_with_style(&self, len: u64, style: ProgressStyle) -> ProgressBar {
        let bar = ProgressBar::with_draw_target(Some(len), self.draw_target());
        bar.set_style(style);
        bar
    }

    /// Create a spinner with a caller-provided style.
    ///
    /// Visible spinners automatically start steady ticking. Hidden spinners
    /// stay inert so quiet and non-interactive runs do not spend time
    /// redrawing.
    #[allow(dead_code)]
    pub(crate) fn spinner_with_style(&self, style: ProgressStyle) -> ProgressBar {
        let spinner = ProgressBar::with_draw_target(None, self.draw_target());
        spinner.set_style(style);
        if self.enabled {
            spinner.enable_steady_tick(DEFAULT_SPINNER_TICK_INTERVAL);
        }
        spinner
    }

    /// Build the default bar style used across the CLI.
    pub(crate) fn default_bar_style() -> ProgressStyle {
        ProgressStyle::default_bar()
            .template(DEFAULT_BAR_TEMPLATE)
            .expect("hardcoded progress template")
    }

    /// Build the default spinner style used across the CLI.
    #[allow(dead_code)]
    pub(crate) fn default_spinner_style() -> ProgressStyle {
        ProgressStyle::default_spinner()
            .template(DEFAULT_SPINNER_TEMPLATE)
            .expect("hardcoded progress template")
    }

    fn draw_target(&self) -> ProgressDrawTarget {
        if self.enabled {
            ProgressDrawTarget::stderr()
        } else {
            ProgressDrawTarget::hidden()
        }
    }
}
