use crate::shared::logging;

#[cfg(feature = "cli")]
use crate::cli_app::{CLI_SEPARATOR_WIDTH, plain_terminal_signature, terminal_signature};

#[cfg(not(feature = "cli"))]
const CLI_SEPARATOR_WIDTH: usize = 48;

/// Print the top-level command banner to the configured primary sink.
pub(crate) fn print_command_banner(command_name: &str) {
    let signature = command_banner_signature();
    logging::write_primary(&signature);
    logging::write_primary_line(&format!("Command: cfdna {command_name}"));
    logging::write_primary_line(&"─".repeat(CLI_SEPARATOR_WIDTH));
}

/// Print the closing command separator to the configured primary sink.
pub(crate) fn print_command_footer() {
    logging::write_primary_line(&"─".repeat(CLI_SEPARATOR_WIDTH));
}

/// Write one logical line to the configured primary sink.
pub(crate) fn write_primary_line(line: &str) {
    logging::write_primary_line(line);
}

#[cfg(feature = "cli")]
fn command_banner_signature() -> String {
    if logging::primary_uses_terminal_formatting() {
        terminal_signature()
    } else {
        plain_terminal_signature()
    }
}

#[cfg(not(feature = "cli"))]
fn command_banner_signature() -> String {
    format!(
        "\n{}\n\n  cfDNAlab\n\n{}\n",
        "_".repeat(CLI_SEPARATOR_WIDTH),
        "─".repeat(CLI_SEPARATOR_WIDTH)
    )
}
