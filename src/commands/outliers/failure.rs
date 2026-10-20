//! Human-readable diagnostics under filenames that never indicate a successful run.

use std::{
    fmt,
    fs::File,
    io::{BufWriter, Write},
    time::{SystemTime, UNIX_EPOCH},
};

use anyhow::{Context, Result};

use super::config::OutliersConfig;
use crate::{
    ToCliCommand,
    shared::io::{FinalOutputFiles, dot_join},
};

/// Collected failures and the already available evidence needed to diagnose them.
#[derive(Debug)]
pub(super) struct AnalysisFailures {
    pub(super) stage: &'static str,
    pub(super) messages: Vec<String>,
    /// Plain sparse TSV, with a header explaining each retained context's coverage counts.
    pub(super) histogram: Option<String>,
    /// Failed merged regions and their contributing evidence, never a normal call audit.
    pub(super) audit: Option<String>,
}

impl fmt::Display for AnalysisFailures {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{}: {} failure(s)\n{}",
            self.stage,
            self.messages.len(),
            self.messages.join("\n")
        )
    }
}

impl std::error::Error for AnalysisFailures {}

/// Save the current attempt's error artifacts only after each has been completely written.
///
/// The text report identifies its companion files, so artifacts from earlier attempts cannot
/// be mistaken for evidence from this attempt. Normal results from previous runs stay intact.
pub(super) fn write_failure_report(
    config: &OutliersConfig,
    error: &anyhow::Error,
) -> Result<std::path::PathBuf> {
    let attempt = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let directory = tempfile::tempdir_in(&config.ioc.output_dir)?;
    let mut outputs = FinalOutputFiles::new(directory.path())?;
    let failures = error.downcast_ref::<AnalysisFailures>();
    let mut companions = Vec::new();
    for (suffix, contents) in [
        (
            "error_histogram.tsv",
            failures.and_then(|failure| failure.histogram.as_deref()),
        ),
        (
            "error_audit.tsv",
            failures.and_then(|failure| failure.audit.as_deref()),
        ),
    ] {
        if let Some(contents) = contents {
            let path = config.ioc.output_dir.join(dot_join(&[
                config.output_prefix.trim(),
                "outliers",
                suffix,
            ]));
            let temporary = outputs.temp_path_for(&path)?;
            let mut writer = BufWriter::new(File::create(&temporary)?);
            writeln!(
                writer,
                "# failed_attempt_unix_nanoseconds={attempt}\n# normal_results_withheld=true"
            )?;
            writer.write_all(contents.as_bytes())?;
            writer.flush()?;
            drop(writer);
            outputs.record(temporary, path.clone())?;
            companions.push(path);
        }
    }
    let path = config.ioc.output_dir.join(dot_join(&[
        config.output_prefix.trim(),
        "outliers.error.txt",
    ]));
    let temporary = outputs.temp_path_for(&path)?;
    let mut writer = BufWriter::new(File::create(&temporary)?);
    writeln!(
        writer,
        "OUTLIERS RUN FAILED\nAttempt (Unix nanoseconds): {attempt}\nVersion: {}\n\nAnalysis failures publish no new normal results. A failure while publishing completed files can leave some new results in place. Previous results may also remain after a failed attempt.\n\n{error:#}\n",
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(
        writer,
        "Settings:\n{}\n",
        config
            .to_cli_string()
            .context("formatting failed run settings")?
    )?;
    if companions.is_empty() {
        writeln!(
            writer,
            "No context histogram or ambiguity audit was available at this stage."
        )?;
    } else {
        writeln!(writer, "Evidence from this attempt:")?;
        for companion in companions {
            writeln!(writer, "{}", companion.display())?;
        }
    }
    writeln!(
        writer,
        "\nRetries replace error files with matching names. Use the attempt identifier to match companions. Error files are not success outputs and may remain after a later successful run."
    )?;
    writer.flush()?;
    drop(writer);
    outputs.record(temporary, path.clone())?;
    outputs.move_into_place()?;
    Ok(path)
}

#[cfg(test)]
mod tests {
    include!("failure_tests.rs");
}
