use crate::{
    command_run::{CommandRunResult, RunOptions},
    commands::{
        cli_common::{ensure_output_dir, resolve_chromosomes_and_contigs, validate_output_prefix},
        counters::FCoverageCounters,
        coverage_weights::scaling_weights_config::ScalingWeightsArgs,
        coverage_weights::striding::{
            StrideBin, fill_triangular_overlap, normalize_weighted_average_overlap_by_global_mean,
        },
        fcoverage::{
            config::{FCoverageConfig, LengthNormalizationMode},
            fcoverage::{FCoverageRunResult, run_fcoverage},
            window_results::CoverageWindowAction,
        },
        run_statistics::{
            DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions, GCStatisticsSummary,
            print_fragment_run_statistics,
        },
    },
    shared::{
        interval::Interval,
        io::{FinalOutputFiles, dot_join, open_text_reader},
        tiled_run::TempDirGuard,
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use std::{
    fs::File,
    io::{BufRead, BufWriter, Write},
    path::Path,
    time::Instant,
};
use tracing::info;

const FCOVERAGE_INTERMEDIATE_DECIMALS: u8 = 12;

/// Result from a scaling-weights command.
///
/// This result is shared by `coverage-weights` and `fragment-count-weights`. Both commands run
/// `fcoverage` internally, smooth stride-level signal, and write a scaling-factor TSV.
#[derive(Debug)]
pub struct ScalingWeightsRunResult {
    /// Counters from the internal `fcoverage` run.
    pub counters: FCoverageCounters,
    /// Final scaling-factor TSV written by the command.
    pub output_scaling_factors: std::path::PathBuf,
    /// Result from the internal `fcoverage` producer run.
    pub source_fcoverage: FCoverageRunResult,
    /// Final output files produced by the command.
    pub output_files: Vec<std::path::PathBuf>,
}

impl CommandRunResult for ScalingWeightsRunResult {
    type Counters = FCoverageCounters;

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[std::path::PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&std::path::Path> {
        Some(self.output_scaling_factors.as_path())
    }
}

/// Result from `coverage-weights`.
///
/// This is the scaling-weights result specialized to average coverage.
pub type CoverageWeightsRunResult = ScalingWeightsRunResult;

#[derive(Clone, Copy)]
pub(crate) enum ScalingWeightsCommand {
    Coverage,
    FragmentCount,
}

impl ScalingWeightsCommand {
    fn fcoverage_window_action(self) -> CoverageWindowAction {
        match self {
            Self::Coverage => CoverageWindowAction::Average,
            Self::FragmentCount => CoverageWindowAction::Total,
        }
    }

    fn fcoverage_value_header(self) -> &'static str {
        match self {
            Self::Coverage => "average_coverage",
            Self::FragmentCount => "total_fragment_mass",
        }
    }

    fn output_file_name(self) -> &'static str {
        match self {
            Self::Coverage => "coverage.scaling_factors.tsv",
            Self::FragmentCount => "fragment_counts.scaling_factors.tsv",
        }
    }

    fn output_value_headers(self) -> (&'static str, &'static str) {
        match self {
            Self::Coverage => ("stride_average_coverage", "smoothed_coverage"),
            Self::FragmentCount => ("stride_fragment_mass", "smoothed_fragment_mass"),
        }
    }

    fn normalization_message(self, global_mean: f32, stride: u32) -> String {
        match self {
            Self::Coverage => format!(
                "Normalized smoothed coverage to global mean: {global_mean} (average coverage per eligible base)"
            ),
            Self::FragmentCount => format!(
                "Normalized smoothed fragment counts to global mean: {global_mean} (unit fragments per {stride} bp stride)"
            ),
        }
    }

    fn info(self, message: &str) {
        match self {
            Self::Coverage => info!(target: "coverage-weights", "{message}"),
            Self::FragmentCount => info!(target: "fragment-count-weights", "{message}"),
        }
    }
}

/// Run the `coverage-weights` command.
///
/// This command estimates broad coverage structure and writes genomic scaling factors that can be
/// reused by other commands. It runs `fcoverage` internally, reads stride-level coverage values,
/// smooths them, and normalizes the smoothed signal to a global mean.
///
/// Reporting is controlled by `options`. `report_statistics` prints the final summary and
/// `log_statuses` controls status messages. This command does not use progress bars.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `coverage-weights` command.
/// - `options`:
///     Reporting controls for statistics and status logs.
///
/// Returns
/// -------
/// - `Ok(CoverageWeightsRunResult)`:
///     The scaling-factor path, internal `fcoverage` result, and counters.
///
/// Errors
/// ------
/// Returns an error if internal `fcoverage` fails, the intermediate TSV is malformed, or the final
/// scaling output cannot be written.
pub fn run_coverage_weights(
    opt: &crate::commands::coverage_weights::config::CoverageWeightsConfig,
    options: RunOptions,
) -> Result<CoverageWeightsRunResult> {
    run_with_fcoverage(
        &opt.shared,
        false,
        ScalingWeightsCommand::Coverage,
        Some(opt.ignore_gap),
        options,
    )
}

/// Shared implementation for coverage-based and count-normalized scaling weights.
///
/// `fragment-count-weights` calls this with `normalize_by_length = true`
/// and reads total unit fragment mass from the internal `fcoverage` output.
///
/// Parameters
/// ----------
/// - `opt`:
///     Shared scaling-weight configuration.
/// - `normalize_by_length`:
///     Whether each fragment contributes unit mass across its span.
/// - `command`:
///     Command flavor used to select output naming, headers, and logging target.
/// - `source_ignore_gap`:
///     Optional `fcoverage` gap-handling override for the internal producer run.
/// - `options`:
///     Reporting controls for statistics and status logs.
///
/// Returns
/// -------
/// - `Ok(ScalingWeightsRunResult)`:
///     The scaling-factor path, internal `fcoverage` result, and counters.
///
/// Errors
/// ------
/// Returns an error if the shared configuration is invalid, internal `fcoverage` fails, the
/// intermediate TSV is malformed, or the final scaling output cannot be written.
pub(crate) fn run_with_fcoverage(
    opt: &ScalingWeightsArgs,
    normalize_by_length: bool,
    command: ScalingWeightsCommand,
    source_ignore_gap: Option<bool>,
    options: RunOptions,
) -> Result<ScalingWeightsRunResult> {
    let start_time = Instant::now();
    let (chromosomes, _contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    opt.check_bin_sizes()?;
    opt.fragment_lengths.validate()?;
    opt.gc.validate(opt.ref_2bit.as_deref())?;
    if source_ignore_gap.unwrap_or(false) && opt.unpaired.reads_are_fragments {
        bail!("--ignore-gap cannot be used with --reads-are-fragments");
    }
    validate_output_prefix(opt.output_prefix.trim())?;

    // Keep all intermediate files under the user-chosen output directory so disk usage stays
    // within the filesystem location the user already selected for results.
    ensure_output_dir(&opt.ioc.output_dir)?;
    let final_temp_dir_guard = TempDirGuard::new(
        &opt.ioc.output_dir,
        &dot_join(&[opt.output_prefix.as_str(), "scaling_weights_final"]),
    )
    .context("creating final output temp directory")?;
    let mut final_outputs = FinalOutputFiles::new(final_temp_dir_guard.path())?;

    let fcoverage_output_dir_guard = TempDirGuard::new(
        &opt.ioc.output_dir,
        &dot_join(&[opt.output_prefix.as_str(), "coverage_weights_source"]),
    )
    .context("creating internal fcoverage output directory")?;
    let fcoverage_output_dir = fcoverage_output_dir_guard.path().to_path_buf();

    let fcoverage_cfg = build_fcoverage_stride_config(
        opt,
        &fcoverage_output_dir,
        normalize_by_length,
        command,
        source_ignore_gap.unwrap_or(false),
    );

    if options.log_statuses {
        command.info("Calling internal fcoverage");
    }
    let fcoverage_result = run_fcoverage(
        &fcoverage_cfg,
        RunOptions {
            report_statistics: false,
            show_progress: false,
            log_statuses: options.log_statuses,
        },
    )
    .context("running internal fcoverage")?;
    if options.log_statuses {
        command.info("Reading internal fcoverage output");
    }

    let mut bins_by_chr = load_stride_bins_from_fcoverage_tsv(
        &fcoverage_result,
        chromosomes.as_slice(),
        command.fcoverage_value_header(),
        opt.stride,
    )?;

    for chromosome in &chromosomes {
        let bins = bins_by_chr
            .get_mut(chromosome)
            .with_context(|| format!("missing stride bins for chromosome '{}'", chromosome))?;
        fill_triangular_overlap(bins, opt.bin_size, opt.stride);
    }

    let mean_weighted_average_overlap =
        normalize_weighted_average_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

    if options.log_statuses {
        command.info(&command.normalization_message(mean_weighted_average_overlap, opt.stride));
    }

    if options.log_statuses {
        command.info("Writing stride-bin coordinates and scaling factors to disk");
    }
    let file_name = dot_join(&[opt.output_prefix.as_str(), command.output_file_name()]);
    let final_output_path = opt.ioc.output_dir.join(&file_name);
    let temp_output_path = final_outputs.temp_path_for(&final_output_path)?;
    let mut tsv_writer =
        BufWriter::new(File::create(&temp_output_path).context("creating scaling-factors TSV")?);
    writeln!(
        tsv_writer,
        "# gc_mode={}",
        crate::shared::scale_genome::scaling_gc_mode_for_run(
            opt.gc.gc_file.is_some(),
            opt.gc.gc_tag.is_some(),
        )
        .as_metadata_value()
    )
    .context("writing TSV metadata")?;
    if let Some(ignore_gap) = source_ignore_gap {
        writeln!(tsv_writer, "# ignore_gap={ignore_gap}").context("writing TSV metadata")?;
    }
    let (stride_value_header, smoothed_value_header) = command.output_value_headers();
    writeln!(
        tsv_writer,
        "chromosome\tstart\tend\t{stride_value_header}\t{smoothed_value_header}\tscaling_factor"
    )
    .context("writing TSV header")?;

    for chromosome in chromosomes {
        let bins = bins_by_chr
            .get(&chromosome)
            .with_context(|| format!("missing bins for chromosome: {}", chromosome))?;

        for bin in bins {
            writeln!(
                tsv_writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                chromosome,
                bin.start(),
                bin.end(),
                bin.stride_value,
                bin.smoothed_value,
                bin.scaling_factor
            )
            .context("writing TSV row")?;
        }
    }

    tsv_writer.flush().context("flushing scaling-factors TSV")?;
    drop(tsv_writer);
    final_outputs.record(temp_output_path, final_output_path.clone())?;
    final_outputs.move_into_place()?;
    if options.log_statuses {
        command.info(&format!("Saved output to: {}", final_output_path.display()));
    }

    let global_counter = fcoverage_result.counters;
    let elapsed = start_time.elapsed();
    if options.report_statistics {
        print_fragment_run_statistics(
            &global_counter.base,
            elapsed,
            FragmentRunStatisticsOptions {
                include_section_header: false,
                notes: &[],
                labels: DEFAULT_FRAGMENT_STATISTICS_LABELS,
                blacklist_excluded_fragments: None,
                gc: (opt.gc.gc_file.is_some() || opt.gc.gc_tag.is_some()).then_some(
                    GCStatisticsSummary {
                        neutralize_invalid_gc: opt.gc.neutralize_invalid_gc,
                        failed_fragments: global_counter.gc_failed_fragments,
                        missing_tags: opt
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(global_counter.gc_missing_tags),
                        out_of_range_tags: opt
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(global_counter.gc_out_of_range_tags),
                    },
                ),
            },
            std::iter::empty::<&str>(),
        );
    }
    Ok(ScalingWeightsRunResult {
        counters: global_counter,
        output_scaling_factors: final_output_path.clone(),
        source_fcoverage: fcoverage_result,
        output_files: vec![final_output_path],
    })
}

fn build_fcoverage_stride_config(
    opt: &ScalingWeightsArgs,
    output_dir: &Path,
    normalize_by_length: bool,
    command: ScalingWeightsCommand,
    ignore_gap: bool,
) -> FCoverageConfig {
    let mut cfg = FCoverageConfig::new(
        crate::commands::cli_common::IOCArgs {
            bam: opt.ioc.bam.clone(),
            output_dir: output_dir.to_path_buf(),
            n_threads: opt.ioc.n_threads,
        },
        opt.chromosomes.clone(),
    );

    cfg.set_unpaired(opt.unpaired.clone());
    cfg.set_normalize_by_length_mode(if normalize_by_length {
        LengthNormalizationMode::UnitMass
    } else {
        LengthNormalizationMode::Off
    });
    cfg.set_output_prefix(opt.output_prefix.clone());
    cfg.set_decimals(FCOVERAGE_INTERMEDIATE_DECIMALS);
    cfg.set_per_window(command.fcoverage_window_action());
    cfg.set_ignore_gap(ignore_gap);
    cfg.set_tile_size(opt.tile_size);
    cfg.set_windows(crate::commands::cli_common::DistributionWindowsArgs {
        by_size: Some(opt.stride as u64),
        by_bed: None,
        by_grouped_bed: None,
    });
    cfg.set_fragment_lengths(opt.fragment_lengths.clone());
    cfg.set_min_mapq(opt.min_mapq);
    cfg.set_require_proper_pair(opt.require_proper_pair);
    cfg.set_blacklist(opt.blacklist.clone());
    cfg.set_gc(opt.gc.clone());
    cfg.set_ref_2bit(opt.ref_2bit.clone());
    cfg
}

fn load_stride_bins_from_fcoverage_tsv(
    fcoverage_result: &FCoverageRunResult,
    chromosomes: &[String],
    value_header: &str,
    stride: u32,
) -> Result<FxHashMap<String, Vec<StrideBin>>> {
    ensure!(
        stride > 0,
        "stride must be greater than zero when loading fcoverage stride bins"
    );
    let path = &fcoverage_result.final_out_path;
    let mut reader = open_text_reader(path)?;
    let mut line = String::new();

    line.clear();
    if reader.read_line(&mut line)? == 0 {
        bail!("{}: empty file; header required", path.display());
    }
    let header = line.trim_end();
    let expected_header = format!("chromosome\tstart\tend\t{value_header}\tblacklisted_positions");
    ensure!(
        header == expected_header,
        "{}: unexpected fcoverage header: '{}'",
        path.display(),
        header
    );

    let mut bins_by_chr =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            continue;
        }

        let cols: Vec<&str> = trimmed.split('\t').collect();
        ensure!(
            cols.len() == 5,
            "{}: expected 5 tab-separated columns, got {} in line '{}'",
            path.display(),
            cols.len(),
            trimmed
        );

        let chromosome = cols[0].to_string();
        let start: u32 = cols[1]
            .parse()
            .with_context(|| format!("{}: invalid start '{}'", path.display(), cols[1]))?;
        let end: u32 = cols[2]
            .parse()
            .with_context(|| format!("{}: invalid end '{}'", path.display(), cols[2]))?;
        let mut stride_value: f32 = cols[3].parse().with_context(|| {
            format!("{}: invalid {} '{}'", path.display(), value_header, cols[3])
        })?;
        let blacklisted_positions: u64 = cols[4].parse().with_context(|| {
            format!(
                "{}: invalid blacklisted_positions '{}'",
                path.display(),
                cols[4]
            )
        })?;
        let interval = Interval::new(start, end).with_context(|| {
            format!(
                "{}: invalid stride-bin interval {}..{} for chromosome '{}'",
                path.display(),
                start,
                end,
                chromosome
            )
        })?;
        let span_positions = u64::from(interval.len());
        ensure!(
            blacklisted_positions <= span_positions,
            "{}: blacklisted_positions {} exceeds row span {}..{}",
            path.display(),
            blacklisted_positions,
            start,
            end
        );
        let eligible_positions = (span_positions - blacklisted_positions) as u32;
        if eligible_positions == 0 {
            stride_value = f32::NAN;
        }
        let support_ratio = (eligible_positions as f64) / (stride as f64);

        bins_by_chr
            .entry(chromosome)
            .or_insert_with(Vec::new)
            .push(StrideBin {
                interval,
                eligible_positions,
                support_ratio,
                stride_value,
                smoothed_value: 0.0,
                scaling_factor: 0.0,
            });
    }

    for chromosome in chromosomes {
        let bins = bins_by_chr
            .get(chromosome)
            .with_context(|| format!("{}: missing chromosome '{}'", path.display(), chromosome))?;
        ensure!(
            !bins.is_empty(),
            "{}: chromosome '{}' had no stride bins",
            path.display(),
            chromosome
        );
        ensure!(
            bins[0].start() == 0,
            "{}: chromosome '{}' did not start at 0",
            path.display(),
            chromosome
        );
        for pair in bins.windows(2) {
            ensure!(
                pair[0].end() == pair[1].start(),
                "{}: chromosome '{}' had non-contiguous stride bins at {}..{} and {}..{}",
                path.display(),
                chromosome,
                pair[0].start(),
                pair[0].end(),
                pair[1].start(),
                pair[1].end()
            );
        }
    }

    Ok(bins_by_chr)
}

#[cfg(test)]
mod tests {
    include!("coverage_weights_tests.rs");
}
