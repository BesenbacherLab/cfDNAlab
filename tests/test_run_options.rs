#![cfg(feature = "cmd_fcoverage")]

use anyhow::Result;
use cfdnalab::RunOptions;
use cfdnalab::run_like_cli::common::{ChromosomeArgs, IOCArgs};
use cfdnalab::run_like_cli::fcoverage::{FCoverageConfig, run_fcoverage};
use cfdnalab::testing::single_contig_inward_pair_bam;
use std::io::{self, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};
use tempfile::TempDir;
use tracing_subscriber::fmt::writer::MakeWriter;

/// Shared byte buffer that captures formatted tracing events.
///
/// `tracing_subscriber` asks its `MakeWriter` for a fresh `Write` handle when it formats an event.
/// This type plays both roles. `make_writer()` returns a cheap clone, and every clone points at the
/// same `Arc<Mutex<Vec<u8>>>`, so all captured events end up in one buffer that the test can read
/// after the public runner returns.
#[derive(Clone, Default)]
struct CapturedTraceBuffer {
    bytes: Arc<Mutex<Vec<u8>>>,
}

impl CapturedTraceBuffer {
    /// Return all captured tracing output as UTF-8 text.
    ///
    /// The `tracing_subscriber` formatter writes valid UTF-8 for the messages used here. Invalid
    /// bytes would indicate a broken test writer or formatter configuration, not a command result.
    fn text(&self) -> String {
        let bytes = self.bytes.lock().expect("captured tracing output lock");
        String::from_utf8(bytes.clone()).expect("captured tracing output should be UTF-8")
    }
}

impl Write for CapturedTraceBuffer {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        self.bytes
            .lock()
            .expect("captured tracing output lock")
            .extend_from_slice(buffer);
        Ok(buffer.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl<'writer> MakeWriter<'writer> for CapturedTraceBuffer {
    type Writer = CapturedTraceBuffer;

    fn make_writer(&'writer self) -> Self::Writer {
        self.clone()
    }
}

/// Build the smallest real `fcoverage` command config needed to reach status logging.
///
/// The fixture has one inward paired-end fragment on `chr1`. The fragment length lower bound stays
/// at the project minimum of 10 bp so the test exercises the normal public runner path.
fn fcoverage_config(bam_path: &Path, output_dir: &Path, output_prefix: &str) -> FCoverageConfig {
    let mut config = FCoverageConfig::new(
        IOCArgs {
            bam: bam_path.to_path_buf(),
            output_dir: output_dir.to_path_buf(),
            n_threads: 1,
        },
        ChromosomeArgs {
            chromosomes: Some(vec!["chr1".to_string()]),
            chromosomes_file: None,
        },
    );
    config.set_output_prefix(output_prefix);
    config.set_min_mapq(0);
    {
        let fragment_lengths = config.fragment_lengths_mut();
        fragment_lengths.min_fragment_length = 10;
        fragment_lengths.max_fragment_length = 200;
    }
    config
}

/// Run `fcoverage` through the public API while collecting application-owned tracing output.
///
/// This is the downstream Rust-side behavior `RunOptions` is meant to support. The test owns the
/// subscriber, while cfDNAlab only decides whether to write status and equivalent CLI events.
fn capture_fcoverage_logs(config: &FCoverageConfig, options: RunOptions) -> Result<String> {
    let captured_trace_buffer = CapturedTraceBuffer::default();
    let subscriber = tracing_subscriber::fmt()
        .with_ansi(false)
        .without_time()
        .with_target(true)
        .with_level(false)
        .with_max_level(tracing::Level::INFO)
        .with_writer(captured_trace_buffer.clone())
        .finish();

    tracing::subscriber::with_default(subscriber, || run_fcoverage(config, options).map(|_| ()))?;
    Ok(captured_trace_buffer.text())
}

/// Run the public `fcoverage` API with a fresh output directory for one reporting profile.
///
/// Fresh directories keep repeated downstream calls independent and avoid output-file collisions.
fn capture_fcoverage_logs_for_options(
    bam_path: &Path,
    output_prefix: &str,
    options: RunOptions,
) -> Result<String> {
    let output_dir = TempDir::new()?;
    let config = fcoverage_config(bam_path, output_dir.path(), output_prefix);
    capture_fcoverage_logs(&config, options)
}

#[test]
fn downstream_run_options_control_status_and_equivalent_cli_tracing() -> Result<()> {
    // Arrange: downstream code uses the public fixture and public config types, then owns the
    // tracing subscriber around each public runner call.
    let bam_fixture = single_contig_inward_pair_bam()?;

    let quiet_logs =
        capture_fcoverage_logs_for_options(&bam_fixture.bam, "quiet", RunOptions::new_quiet())?;
    assert!(
        quiet_logs.trim().is_empty(),
        "quiet RunOptions should not write tracing status or equivalent CLI output.\nlogs:\n{quiet_logs}"
    );

    let status_only_logs = capture_fcoverage_logs_for_options(
        &bam_fixture.bam,
        "status_only",
        RunOptions {
            log_statuses: true,
            ..RunOptions::new_quiet()
        },
    )?;
    assert!(
        status_only_logs.contains("fcoverage: Counting per tile"),
        "status logging should include fcoverage milestones.\nlogs:\n{status_only_logs}"
    );
    assert!(
        status_only_logs.contains("fcoverage: Merging temporary tile files to final output"),
        "status logging should include the reduction milestone.\nlogs:\n{status_only_logs}"
    );
    assert!(
        !status_only_logs.contains("Equivalent CLI:"),
        "log_equivalent_cli=false should suppress equivalent CLI output.\nlogs:\n{status_only_logs}"
    );

    let equivalent_cli_only_logs = capture_fcoverage_logs_for_options(
        &bam_fixture.bam,
        "equivalent_cli_only",
        RunOptions {
            log_equivalent_cli: true,
            ..RunOptions::new_quiet()
        },
    )?;
    assert!(
        equivalent_cli_only_logs.contains("fcoverage: Equivalent CLI: cfdna fcoverage"),
        "log_equivalent_cli=true should write the rendered command.\nlogs:\n{equivalent_cli_only_logs}"
    );
    assert!(
        !equivalent_cli_only_logs.contains("fcoverage: Counting per tile"),
        "log_statuses=false should suppress status milestones.\nlogs:\n{equivalent_cli_only_logs}"
    );

    Ok(())
}
