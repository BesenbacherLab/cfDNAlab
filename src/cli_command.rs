use crate::Result;
use std::ffi::OsString;

pub(crate) mod helpers;

use helpers::shell_quote;

/// Render a filled command config as an equivalent `cfdna` CLI invocation.
///
/// `to_cli_args()` returns argument tokens. Use this when passing the command
/// to another process API. `to_cli_string()` is for display in logs, reports,
/// and provenance text.
pub trait ToCliCommand {
    /// Return the full argument vector, including `cfdna` and the subcommand.
    fn to_cli_args(&self) -> Result<Vec<OsString>>;

    /// Return a shell-quoted command string for display.
    fn to_cli_string(&self) -> Result<String> {
        Ok(self
            .to_cli_args()?
            .iter()
            .map(shell_quote)
            .collect::<Vec<_>>()
            .join(" "))
    }
}
