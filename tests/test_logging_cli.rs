#![cfg(feature = "cli")]

// KEEP-IN-TESTS: CLI logging tests exercise binary behavior and user-visible output.

mod fixtures;

use anyhow::{Context, Result, bail};
use cfdnalab::testing::single_contig_inward_pair_bam;
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::TempDir;

fn binary_name(base_name: &str) -> String {
    if cfg!(windows) {
        format!("{base_name}.exe")
    } else {
        base_name.to_string()
    }
}

fn cfdna_bin_path() -> Result<PathBuf> {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_cfdna") {
        return Ok(PathBuf::from(path));
    }

    let current_exe = std::env::current_exe().context("failed to read current test binary path")?;
    let deps_dir = current_exe
        .parent()
        .context("failed to derive deps directory from current test binary path")?;
    let target_dir = deps_dir
        .parent()
        .context("failed to derive target directory from deps path")?;

    let direct_path = target_dir.join(binary_name("cfdna"));
    if direct_path.is_file() {
        return Ok(direct_path);
    }

    let mut hashed_candidates = Vec::new();
    for entry in std::fs::read_dir(deps_dir).with_context(|| {
        format!(
            "failed to list candidate binaries in {}",
            deps_dir.display()
        )
    })? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let file_name = match path.file_name().and_then(OsStr::to_str) {
            Some(name) => name,
            None => continue,
        };
        let extension = path.extension().and_then(OsStr::to_str);
        let looks_like_hashed_binary = file_name.starts_with("cfdna-");
        let is_makefile_dep = extension == Some("d");
        if looks_like_hashed_binary && !is_makefile_dep {
            hashed_candidates.push(path);
        }
    }
    hashed_candidates.sort_by_key(|path| {
        std::fs::metadata(path)
            .and_then(|meta| meta.modified())
            .ok()
            .unwrap_or(std::time::SystemTime::UNIX_EPOCH)
    });
    if let Some(path) = hashed_candidates.into_iter().last() {
        return Ok(path);
    }

    bail!(
        "Could not locate cfdna binary. Tried CARGO_BIN_EXE_cfdna, {}, and hashed binaries under {}",
        direct_path.display(),
        deps_dir.display()
    );
}

fn command_output(command_name: &str, args: &[&str]) -> Result<std::process::Output> {
    Command::new(cfdna_bin_path()?)
        .arg(command_name)
        .args(args)
        .output()
        .with_context(|| format!("failed running cfdna {command_name} {}", args.join(" ")))
}

fn assert_success(output: &std::process::Output, command_desc: &str) {
    let stdout_text = String::from_utf8_lossy(&output.stdout);
    let stderr_text = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "Expected {command_desc} to succeed.\nstdout:\n{stdout_text}\nstderr:\n{stderr_text}"
    );
}

fn path_text(path: &Path) -> String {
    path.to_string_lossy().to_string()
}

#[cfg(feature = "cmd_fcoverage")]
#[test]
fn fcoverage_log_stdout_routes_normal_messages_to_stdout() -> Result<()> {
    // Arrange: use the minimal one-fragment fixture so the command completes quickly while still
    // reaching the banner, lifecycle logging, and final statistics block.
    let bam_fixture = single_contig_inward_pair_bam()?;
    let out_dir = TempDir::new()?;
    let bam_path = path_text(&bam_fixture.bam);
    let out_path = path_text(out_dir.path());

    // Act
    let output = command_output(
        "fcoverage",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
            "--output-prefix",
            "log_stdout",
            "--log",
            "stdout",
        ],
    )?;
    assert_success(&output, "cfdna fcoverage --log stdout");

    // Assert
    let stdout = String::from_utf8(output.stdout).context("stdout is not valid UTF-8")?;
    let stderr = String::from_utf8(output.stderr).context("stderr is not valid UTF-8")?;
    assert!(stdout.contains("Command: cfdna fcoverage"));
    assert!(stdout.contains("fcoverage: Counting per tile"));
    assert!(stdout.contains("Statistics"));
    assert!(
        stderr.trim().is_empty(),
        "Expected no stderr output without warnings.\nstderr:\n{stderr}"
    );

    Ok(())
}

#[cfg(feature = "cmd_fcoverage")]
#[test]
fn fcoverage_log_quiet_suppresses_normal_cli_output() -> Result<()> {
    // Arrange: the same tiny fixture should produce no warnings, so quiet mode can be asserted as
    // fully silent in this non-TTY integration test.
    let bam_fixture = single_contig_inward_pair_bam()?;
    let out_dir = TempDir::new()?;
    let bam_path = path_text(&bam_fixture.bam);
    let out_path = path_text(out_dir.path());

    // Act
    let output = command_output(
        "fcoverage",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
            "--output-prefix",
            "log_quiet",
            "--log",
            "quiet",
        ],
    )?;
    assert_success(&output, "cfdna fcoverage --log quiet");

    // Assert
    let stdout = String::from_utf8(output.stdout).context("stdout is not valid UTF-8")?;
    let stderr = String::from_utf8(output.stderr).context("stderr is not valid UTF-8")?;
    assert!(
        stdout.trim().is_empty(),
        "Expected quiet mode to suppress stdout.\nstdout:\n{stdout}"
    );
    assert!(
        stderr.trim().is_empty(),
        "Expected quiet mode to produce no stderr output for a warning-free run.\nstderr:\n{stderr}"
    );

    Ok(())
}

#[cfg(feature = "cmd_fcoverage")]
#[test]
fn fcoverage_log_file_writes_to_auto_generated_log_under_output_dir() -> Result<()> {
    // Arrange: plain `--log file` should pick `<output_dir>/logs/` and use a generated file name.
    let bam_fixture = single_contig_inward_pair_bam()?;
    let out_dir = TempDir::new()?;
    let bam_path = path_text(&bam_fixture.bam);
    let out_path = path_text(out_dir.path());

    // Act
    let output = command_output(
        "fcoverage",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
            "--output-prefix",
            "log_file",
            "--log",
            "file",
        ],
    )?;
    assert_success(&output, "cfdna fcoverage --log file");

    // Assert
    let stdout = String::from_utf8(output.stdout).context("stdout is not valid UTF-8")?;
    let stderr = String::from_utf8(output.stderr).context("stderr is not valid UTF-8")?;
    assert!(
        stdout.trim().is_empty(),
        "Expected file mode to avoid stdout.\nstdout:\n{stdout}"
    );
    assert!(
        stderr.trim().is_empty(),
        "Expected file mode to keep stderr empty for a warning-free run.\nstderr:\n{stderr}"
    );

    let logs_dir = out_dir.path().join("logs");
    let mut entries: Vec<PathBuf> = std::fs::read_dir(&logs_dir)
        .with_context(|| format!("reading {}", logs_dir.display()))?
        .map(|entry| entry.map(|value| value.path()))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    entries.sort();
    assert_eq!(
        entries.len(),
        1,
        "Expected exactly one auto-generated log file under {}",
        logs_dir.display()
    );

    let log_path = &entries[0];
    let log_name = log_path
        .file_name()
        .and_then(OsStr::to_str)
        .unwrap_or_default();
    assert!(
        log_name.starts_with("fcoverage_") && log_name.ends_with(".log"),
        "Expected an auto-generated fcoverage log name, got {}",
        log_name
    );

    let log_text = std::fs::read_to_string(log_path)
        .with_context(|| format!("reading {}", log_path.display()))?;
    assert!(log_text.contains("Command: cfdna fcoverage"));
    assert!(log_text.contains("fcoverage: Counting per tile"));
    assert!(log_text.contains("Statistics"));

    Ok(())
}

#[cfg(feature = "cmd_coverage_weights")]
#[test]
fn coverage_weights_log_stdout_keeps_one_top_level_banner_and_shows_nested_fcoverage_logs()
-> Result<()> {
    // Arrange: this command reuses internal fcoverage. The log output should show both command
    // phases, but only one top-level banner because only the binary owns that presentation.
    let bam_fixture = single_contig_inward_pair_bam()?;
    let out_dir = TempDir::new()?;
    let bam_path = path_text(&bam_fixture.bam);
    let out_path = path_text(out_dir.path());

    // Act
    let output = command_output(
        "coverage-weights",
        &[
            "--bam",
            bam_path.as_str(),
            "--output-dir",
            out_path.as_str(),
            "--chromosomes",
            "chr1",
            "--n-threads",
            "1",
            "--bin-size",
            "40",
            "--stride",
            "20",
            "--min-mapq",
            "0",
            "--min-fragment-length",
            "10",
            "--max-fragment-length",
            "200",
            "--output-prefix",
            "nested_logs",
            "--log",
            "stdout",
        ],
    )?;
    assert_success(&output, "cfdna coverage-weights --log stdout");

    // Assert
    let stdout = String::from_utf8(output.stdout).context("stdout is not valid UTF-8")?;
    let stderr = String::from_utf8(output.stderr).context("stderr is not valid UTF-8")?;
    assert!(stdout.contains("Command: cfdna coverage-weights"));
    assert!(stdout.contains("coverage-weights: Calling internal fcoverage"));
    assert!(stdout.contains("fcoverage: Counting per tile"));
    assert!(!stdout.contains("Command: cfdna fcoverage"));
    assert_eq!(
        stdout.matches("Command: cfdna ").count(),
        1,
        "Expected exactly one top-level command banner.\nstdout:\n{stdout}"
    );
    assert!(
        stderr.trim().is_empty(),
        "Expected no stderr output without warnings.\nstderr:\n{stderr}"
    );

    let scaling_path = out_dir
        .path()
        .join("nested_logs.coverage.scaling_factors.tsv");
    assert!(scaling_path.exists(), "Expected {}", scaling_path.display());

    Ok(())
}
