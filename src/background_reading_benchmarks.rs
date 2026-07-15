//! Temporary benchmarks for foreground versus background text reading.
//!
//! These are ignored tests so they can call crate-internal loaders without exposing benchmark-only
//! APIs publicly. Run them in release mode on representative files:
//!
//! ```text
//! CFDNALAB_BENCH_BED=/path/windows.bed \
//! CFDNALAB_BENCH_FRAG=/path/fragments.frag.tsv.zst \
//! CFDNALAB_BENCH_LENGTHS_TSV=/path/length_counts.tsv.zst \
//! CFDNALAB_BENCH_FCOVERAGE_TSV=/path/fcoverage.average.tsv.zst \
//! cargo test --release --features testing background_reading \
//!     -- --ignored --nocapture --test-threads=1
//! ```
//!
//! `CFDNALAB_BENCH_RUNS` controls repetitions and defaults to 6. Set
//! `CFDNALAB_BENCH_BLACKLIST_BED` to benchmark a different BED as a blacklist. Set
//! `CFDNALAB_BENCH_FCOVERAGE_GROUP_INDEX` when the fcoverage input requires its group-index TSV.

use anyhow::{Context, Result, ensure};
use std::{
    env,
    hint::black_box,
    path::PathBuf,
    time::{Duration, Instant},
};

const DEFAULT_BENCHMARK_RUNS: usize = 6;

#[derive(Debug)]
struct DurationSummary {
    minimum_seconds: f64,
    maximum_seconds: f64,
    mean_seconds: f64,
    median_seconds: f64,
}

fn benchmark_runs() -> Result<usize> {
    let runs = env::var("CFDNALAB_BENCH_RUNS")
        .map(|value| {
            value
                .parse::<usize>()
                .context("CFDNALAB_BENCH_RUNS must be a positive integer")
        })
        .unwrap_or(Ok(DEFAULT_BENCHMARK_RUNS))?;
    ensure!(
        runs >= 2,
        "CFDNALAB_BENCH_RUNS must be at least 2 so each mode is measured multiple times"
    );
    Ok(runs)
}

pub(crate) fn required_path(variable_name: &str) -> Result<PathBuf> {
    let value = env::var(variable_name)
        .with_context(|| format!("set {variable_name} to the benchmark input path"))?;
    let path = PathBuf::from(value);
    ensure!(
        path.is_file(),
        "{variable_name} does not point to a file: {}",
        path.display()
    );
    Ok(path)
}

pub(crate) fn optional_path(variable_name: &str) -> Result<Option<PathBuf>> {
    match env::var(variable_name) {
        Ok(value) => {
            let path = PathBuf::from(value);
            ensure!(
                path.is_file(),
                "{variable_name} does not point to a file: {}",
                path.display()
            );
            Ok(Some(path))
        }
        Err(env::VarError::NotPresent) => Ok(None),
        Err(error) => Err(error).with_context(|| format!("reading {variable_name}")),
    }
}

fn summarize(durations: &[Duration]) -> DurationSummary {
    let mut seconds: Vec<f64> = durations.iter().map(Duration::as_secs_f64).collect();
    seconds.sort_by(f64::total_cmp);
    let middle = seconds.len() / 2;
    let median_seconds = if seconds.len().is_multiple_of(2) {
        (seconds[middle - 1] + seconds[middle]) / 2.0
    } else {
        seconds[middle]
    };

    DurationSummary {
        minimum_seconds: seconds[0],
        maximum_seconds: seconds[seconds.len() - 1],
        mean_seconds: seconds.iter().sum::<f64>() / seconds.len() as f64,
        median_seconds,
    }
}

fn print_summary(label: &str, mode: &str, durations: &[Duration]) {
    let summary = summarize(durations);
    eprintln!(
        "{label} [{mode}] runs={} max={:.3}s min={:.3}s mean={:.3}s median={:.3}s",
        durations.len(),
        summary.maximum_seconds,
        summary.minimum_seconds,
        summary.mean_seconds,
        summary.median_seconds,
    );
}

pub(crate) fn compare_read_modes<Output>(
    label: &str,
    mut operation: impl FnMut(bool) -> Result<Output>,
) -> Result<()> {
    let multiple_threads_are_available =
        std::thread::available_parallelism().is_ok_and(|thread_count| thread_count.get() > 1);
    ensure!(
        multiple_threads_are_available,
        "background-reading benchmarks require more than one available execution thread"
    );
    let runs = benchmark_runs()?;
    let mut foreground_durations = Vec::with_capacity(runs);
    let mut background_durations = Vec::with_capacity(runs);

    for run_index in 0..runs {
        let modes = if run_index.is_multiple_of(2) {
            [false, true]
        } else {
            [true, false]
        };
        for read_in_background in modes {
            let started = Instant::now();
            let output = operation(read_in_background)?;
            let elapsed = started.elapsed();
            black_box(&output);
            drop(output);

            if read_in_background {
                background_durations.push(elapsed);
            } else {
                foreground_durations.push(elapsed);
            }
        }
    }

    print_summary(label, "foreground", &foreground_durations);
    print_summary(label, "background", &background_durations);
    Ok(())
}

#[cfg(loads_grouped_bed)]
#[test]
#[ignore = "requires CFDNALAB_BENCH_BED and measures wall-clock performance"]
fn benchmark_grouped_bed_loader() -> Result<()> {
    let path = required_path("CFDNALAB_BENCH_BED")?;
    compare_read_modes("grouped BED loader", |read_in_background| {
        crate::shared::bed::load_grouped_windows_from_bed(
            &path,
            None,
            false,
            None,
            None,
            read_in_background,
        )
    })
}

#[test]
#[ignore = "requires CFDNALAB_BENCH_BED and measures wall-clock performance"]
fn benchmark_plain_bed_loader() -> Result<()> {
    let path = required_path("CFDNALAB_BENCH_BED")?;
    compare_read_modes("plain BED loader", |read_in_background| {
        crate::shared::bed::load_windows_from_bed(&path, None, None, None, read_in_background)
    })
}

#[test]
#[ignore = "requires a benchmark BED and measures wall-clock performance"]
fn benchmark_blacklist_loader() -> Result<()> {
    let path = match optional_path("CFDNALAB_BENCH_BLACKLIST_BED")? {
        Some(path) => path,
        None => required_path("CFDNALAB_BENCH_BED")?,
    };
    compare_read_modes("blacklist loader", |read_in_background| {
        crate::shared::blacklist::load_blacklists(&[path.as_path()], 1, 0, None, read_in_background)
    })
}

#[cfg(feature = "cmd_frag_to_bam")]
#[test]
#[ignore = "requires CFDNALAB_BENCH_FRAG and measures wall-clock performance"]
fn benchmark_frag_to_bam_input_loader() -> Result<()> {
    let path = required_path("CFDNALAB_BENCH_FRAG")?;
    compare_read_modes("frag-to-bam input loader", |read_in_background| {
        crate::commands::frag_to_bam::frag_to_bam::benchmark_frag_input_loading(
            &path,
            read_in_background,
        )
    })
}
