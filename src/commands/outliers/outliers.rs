use std::{
    path::{Path, PathBuf},
    sync::Arc,
    time::Instant,
};

use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use tracing::info;

use crate::{
    ToCliCommand,
    command_run::{CommandRunResult, RunOptions, status_info},
    commands::{
        cli_common::{
            ensure_output_dir, load_blacklist_map, resolve_chromosomes_and_contigs,
            validate_output_prefix,
        },
        counters::FCoverageCounters,
        fcoverage::fcoverage::add_fragment_clipped_to_core,
        outliers::{
            config::{OutlierTarget, OutliersConfig},
            model::fit_two_stage_zip,
            outputs::{
                OutputMetadata, OutputPaths, write_call_outputs, write_diagnostic_outputs,
                write_histogram_output,
            },
            striding::{CoreHistogram, CoreModel, fit_core_models, sum_global_histogram},
        },
        run_statistics::{
            DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions,
            TILE_DOUBLE_COUNT_NOTE, print_fragment_run_statistics,
        },
    },
    shared::{
        bam::{Contigs, create_chromosome_reader},
        coverage::Coverage,
        fragment::segment_fragment::FragmentWithSegments,
        fragment_iterators::fragments_with_segments_from_bam,
        interval::Interval,
        io::{FinalOutputFiles, dot_join},
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        thread_pool::init_global_pool,
        tiled_run::{RunTempDirs, Tile, build_tiles},
    },
};

const COMMAND_TARGET: &str = "outliers";

/// Result from `cfdna outliers`.
#[derive(Debug)]
pub struct OutliersRunResult {
    /// Fragment counters from the final tiled BAM scan.
    pub counters: FCoverageCounters,
    /// Sparse fractional fragment weights intended for the planned downstream weight loader.
    pub output_keep_weights: PathBuf,
    /// Exact called regions for blacklist workflows that should not add a flank.
    pub output_exact_blacklist: PathBuf,
    /// Conservative blacklist with nearby positions added around each called region.
    pub output_flanked_blacklist: PathBuf,
    /// Per-core raw data for auditing coverage and fragment support.
    pub output_histograms: PathBuf,
    /// Local fits, actual thresholds, fallback use, and fit diagnostics for every model core.
    pub output_models: PathBuf,
    /// Machine-readable global observed and fitted coverage distribution.
    pub output_global_fit: PathBuf,
    /// Visual assessment of the complete global histogram and per-core models used during calling.
    pub output_plot: PathBuf,
    /// All final output files produced by the command.
    pub output_files: Vec<PathBuf>,
}

impl CommandRunResult for OutliersRunResult {
    type Counters = FCoverageCounters;

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&Path> {
        Some(self.output_keep_weights.as_path())
    }
}

#[derive(Debug)]
struct TileCoverage {
    coverage: Coverage,
    counters: FCoverageCounters,
    first_model_core_index: usize,
    fragment_support: Vec<u64>,
}

/// Histograms and counters produced for a processing tile during the first BAM pass.
#[derive(Debug)]
struct HistogramPassTileResult {
    chromosome: String,
    core_histograms: Vec<(usize, CoreHistogram)>,
    counters: FCoverageCounters,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct CalledRegion {
    pub(super) interval: Interval<u32>,
    pub(super) observed_mass: f64,
    pub(super) target_mass: f64,
}

impl CalledRegion {
    pub(super) fn keep_weight(&self) -> f64 {
        if self.target_mass == 0.0 {
            0.0
        } else {
            (self.target_mass / self.observed_mass).min(1.0)
        }
    }
}

/// Called regions and counters produced for a processing tile during the second BAM pass.
#[derive(Debug)]
struct CallingPassTileResult {
    chromosome: String,
    regions: Vec<CalledRegion>,
    counters: FCoverageCounters,
}

/// Detect high raw-coverage outliers and calculate regional fragment keep weights.
///
/// Coverage is scanned in parallel tiles twice. The first pass builds histograms for fixed-size
/// model cores and fits local two-stage ZIP models. The second pass applies the finalized core
/// thresholds to original raw coverage and joins qualifying positions across core and tile
/// boundaries.
///
/// For a called region, `observed_mass` is the sum of its original positional coverage and
/// `target_mass` is the sum of its selected target coverage. The written regional weight is
/// `min(1, target_mass / observed_mass)`. Positions omitted from the sparse TSV have implicit
/// weight `1.0`. A future downstream loader will assign a fragment the minimum weight among regions
/// overlapping its complete `pos` to `reference_end` span. No downstream command applies these
/// weights yet.
///
/// The final ZIP parameters describe an underlying untruncated distribution even though they are
/// estimated with a likelihood conditioned on coverage below `T1`. Calls use the inclusive tail
/// `P(X >= T2)`. The model does not account for positive-count overdispersion or spatial dependence,
/// so the diagnostic outputs are required model checks rather than optional presentation files.
///
/// Reporting is controlled by `options`. Counters describe the second pass because each pass uses
/// the same fragment filters and tiled fetch pattern.
///
/// Parameters
/// ----------
/// - `config`:
///   Fully resolved command configuration.
/// - `options`:
///   Reporting controls for statistics, progress bars, status logs, and equivalent CLI output.
///
/// Returns
/// -------
/// - `Ok(OutliersRunResult)`:
///   Final output paths and counters from the calling pass.
///
/// Errors
/// ------
/// Returns an error when inputs are invalid, a ZIP model cannot be fitted, raw coverage is not
/// integer-valued, or an output cannot be completed.
pub fn run_outliers(config: &OutliersConfig, options: RunOptions) -> Result<OutliersRunResult> {
    let start_time = Instant::now();
    config.validate()?;
    validate_output_prefix(config.output_prefix.trim())?;

    if options.log_equivalent_cli {
        let command_text = config.to_cli_string()?;
        let message = crate::command_run::equivalent_cli_log_message(&command_text);
        info!(target: COMMAND_TARGET, "{message}");
    }

    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&config.chromosomes, config.ioc.bam.as_path())?;
    ensure_output_dir(&config.ioc.output_dir)?;
    init_global_pool(config.ioc.n_threads)?;

    if config.blacklist.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Loading blacklists");
    }
    let blacklist_map = Arc::new(load_blacklist_map(
        config.blacklist.as_ref(),
        1,
        0,
        &chromosomes,
        config.ioc.n_threads > 1,
    )?);

    let halo = config.fragment_lengths.max_fragment_length;
    let (tiles, _) = build_tiles(
        &chromosomes,
        &contigs,
        config.tile_size,
        halo,
        Some(config.stride as u64),
    )?;
    ensure!(
        !tiles.is_empty(),
        "no tiles were created for the selected chromosomes"
    );

    let work_root = config
        .temp
        .temp_dir
        .as_deref()
        .unwrap_or(&config.ioc.output_dir);
    let run_temp_dirs = RunTempDirs::new(
        work_root,
        &config.ioc.output_dir,
        &dot_join(&[COMMAND_TARGET, config.output_prefix.trim()]),
    )
    .context("creating outliers temporary directories")?;
    let mut final_outputs = FinalOutputFiles::new(run_temp_dirs.final_output_dir())?;
    let output_paths = OutputPaths::new(&config.ioc.output_dir, config.output_prefix.trim());

    status_info!(options, target: COMMAND_TARGET, "First BAM pass: building raw-coverage histograms");
    let first_progress = Arc::new(
        ProgressFactory::with_enabled(options.show_progress).default_bar(tiles.len() as u64),
    );
    let histogram_pass_tile_results = tiles
        .par_iter()
        .map(|tile| {
            let blacklist = blacklist_map
                .get(&tile.chr)
                .map(Vec::as_slice)
                .unwrap_or(&[]);
            let result = process_histogram_pass_tile(config, tile, blacklist);
            first_progress.inc(1);
            result
        })
        .collect::<Result<Vec<_>>>()?;
    first_progress.finish_and_clear();

    let (histograms_by_chromosome, _first_pass_counters) = merge_histogram_pass_results(
        &chromosomes,
        &contigs,
        config.stride,
        histogram_pass_tile_results,
    )?;
    let global_histogram = sum_global_histogram(&histograms_by_chromosome)?;
    let eligible_positions = global_histogram.eligible_positions()?;
    let tail_probability = config
        .tail_probability
        .unwrap_or_else(|| config.tail_probability_multiplier / eligible_positions as f64);
    ensure!(
        tail_probability.is_finite() && tail_probability > 0.0 && tail_probability <= 0.5,
        "effective tail probability must be in (0, 0.5], got {}",
        tail_probability
    );
    let blacklist_flank = config
        .blacklist_flank
        .unwrap_or(config.fragment_lengths.max_fragment_length);
    let output_metadata = OutputMetadata {
        config,
        chromosomes: &chromosomes,
        tail_probability,
        eligible_positions,
        blacklist_flank,
    };

    status_info!(options, target: COMMAND_TARGET, "Writing per-core coverage histograms");
    write_histogram_output(
        &output_metadata,
        &histograms_by_chromosome,
        &output_paths.histograms,
        &mut final_outputs,
    )?;

    status_info!(options, target: COMMAND_TARGET, "Fitting global and local ZIP models");
    let global_model = fit_two_stage_zip(&global_histogram.counts, tail_probability)
        .context("fitting the global two-stage ZIP model")?;
    let models_by_chromosome = fit_core_models(
        &chromosomes,
        &histograms_by_chromosome,
        global_model,
        tail_probability,
        config.bin_size,
        config.min_context_fragments,
    )?;

    status_info!(options, target: COMMAND_TARGET, "Writing ZIP diagnostics");
    write_diagnostic_outputs(
        &output_metadata,
        &global_histogram,
        global_model,
        &models_by_chromosome,
        &output_paths,
        &mut final_outputs,
    )?;

    status_info!(options, target: COMMAND_TARGET, "Second BAM pass: calling outlier regions");
    let second_progress = Arc::new(
        ProgressFactory::with_enabled(options.show_progress).default_bar(tiles.len() as u64),
    );
    let calling_pass_tile_results = tiles
        .par_iter()
        .map(|tile| {
            let blacklist = blacklist_map
                .get(&tile.chr)
                .map(Vec::as_slice)
                .unwrap_or(&[]);
            let chromosome_models = models_by_chromosome
                .get(&tile.chr)
                .with_context(|| format!("missing fitted models for chromosome '{}'", tile.chr))?;
            let result = process_calling_pass_tile(config, tile, blacklist, chromosome_models);
            second_progress.inc(1);
            result
        })
        .collect::<Result<Vec<_>>>()?;
    second_progress.finish_and_clear();

    let (regions_by_chromosome, counters) =
        merge_calling_pass_results(&chromosomes, calling_pass_tile_results)?;
    write_call_outputs(
        &output_metadata,
        &contigs,
        &regions_by_chromosome,
        &output_paths,
        &mut final_outputs,
    )?;

    final_outputs.move_into_place()?;
    let output_files = output_paths.all();
    status_info!(options, target: COMMAND_TARGET, "Saved outlier outputs");

    if options.report_statistics {
        let called_regions = regions_by_chromosome.values().map(Vec::len).sum::<usize>();
        let extra_statistics = [
            format!("Eligible positions: {eligible_positions}"),
            format!("Tail probability: {tail_probability}"),
            format!("Called outlier regions: {called_regions}"),
        ];
        print_fragment_run_statistics(
            &counters.base,
            start_time.elapsed(),
            FragmentRunStatisticsOptions {
                include_section_header: false,
                notes: &[TILE_DOUBLE_COUNT_NOTE],
                labels: DEFAULT_FRAGMENT_STATISTICS_LABELS,
                blacklist_excluded_fragments: None,
                gc: None,
            },
            extra_statistics.iter().map(String::as_str),
        );
    }

    Ok(OutliersRunResult {
        counters,
        output_keep_weights: output_paths.keep_weights,
        output_exact_blacklist: output_paths.exact_blacklist,
        output_flanked_blacklist: output_paths.flanked_blacklist,
        output_histograms: output_paths.histograms,
        output_models: output_paths.models,
        output_global_fit: output_paths.global_fit,
        output_plot: output_paths.plot,
        output_files,
    })
}

/// Process one parallelization tile during the histogram-building BAM pass.
fn process_histogram_pass_tile(
    config: &OutliersConfig,
    tile: &Tile,
    blacklist: &[Interval<u64>],
) -> Result<HistogramPassTileResult> {
    let tile_coverage = build_tile_coverage(config, tile, blacklist, true)?;
    let coverage = tile_coverage
        .coverage
        .coverage()
        .context("tile coverage missing after finalization")?;
    let mask = tile_coverage.coverage.blacklist_mask();
    let mut histograms = Vec::new();
    let mut absolute_start = tile.core_start();

    while absolute_start < tile.core_end() {
        let model_core_index = (absolute_start / config.stride) as usize;
        let overlap_end = ((model_core_index as u64 + 1) * config.stride as u64)
            .min(tile.core_end() as u64) as u32;
        let local_start = (absolute_start - tile.core_start()) as usize;
        let local_end = (overlap_end - tile.core_start()) as usize;
        let mut histogram = CoreHistogram::new(Interval::new(absolute_start, overlap_end)?);
        histogram.add_coverage(
            &coverage[local_start..local_end],
            mask.map(|values| &values[local_start..local_end]),
        )?;
        let support_offset = model_core_index
            .checked_sub(tile_coverage.first_model_core_index)
            .context("model-core support index precedes tile support range")?;
        histogram.fragment_support = tile_coverage
            .fragment_support
            .get(support_offset)
            .copied()
            .unwrap_or(0);
        histograms.push((model_core_index, histogram));
        absolute_start = overlap_end;
    }

    Ok(HistogramPassTileResult {
        chromosome: tile.chr.clone(),
        core_histograms: histograms,
        counters: tile_coverage.counters,
    })
}

/// Build raw positional coverage for one processing tile in either BAM pass.
fn build_tile_coverage(
    config: &OutliersConfig,
    tile: &Tile,
    blacklist: &[Interval<u64>],
    collect_fragment_support: bool,
) -> Result<TileCoverage> {
    let (mut reader, bam_tid, chromosome_length) =
        create_chromosome_reader(&config.ioc.bam, &tile.chr)?;
    tile.ensure_matches_bam_tid(bam_tid)?;
    let (fetch_start, fetch_end) = tile.fetch.try_to_i64()?.as_tuple();
    reader
        .fetch((tile.tid, fetch_start, fetch_end))
        .with_context(|| format!("fetching {} {}-{}", tile.chr, fetch_start, fetch_end))?;

    let mut coverage = Coverage::new(tile.core.len());
    let mut counters = FCoverageCounters::default();
    let fragment_lengths = config.fragment_lengths.clone();
    let fragment_filter =
        move |fragment: &FragmentWithSegments| fragment_lengths.contains(fragment.len());
    let include_read: Box<dyn Fn(&Record) -> bool + Send + Sync> = if config
        .unpaired
        .reads_are_fragments
    {
        let minimum_mapping_quality = config.min_mapq;
        Box::new(move |record| default_include_read_unpaired(record, minimum_mapping_quality))
    } else {
        let minimum_mapping_quality = config.min_mapq;
        let require_proper_pair = config.require_proper_pair;
        Box::new(move |record| {
            default_include_read_paired_end(record, require_proper_pair, minimum_mapping_quality)
        })
    };
    let mut fragments = fragments_with_segments_from_bam(
        reader
            .records()
            .map(|record| record.map_err(anyhow::Error::from)),
        move |record| include_read(record),
        1,
        !config.ignore_gap,
        None,
        fragment_filter,
        config.unpaired.reads_are_fragments,
    )
    .with_local_counters();

    let first_model_core_index = (tile.core_start() / config.stride) as usize;
    let last_model_core_index = ((tile.core_end() - 1) / config.stride) as usize;
    let mut fragment_support = if collect_fragment_support {
        vec![0_u64; last_model_core_index - first_model_core_index + 1]
    } else {
        Vec::new()
    };

    for fragment_result in fragments.by_ref() {
        let fragment = fragment_result.context("reading fragment")?;
        let counted = add_fragment_clipped_to_core(&mut coverage, &fragment, 1.0, tile.core)?;
        if !counted {
            continue;
        }
        counters.base.counted_fragments += 1;

        if collect_fragment_support {
            let midpoint = fragment.start() + fragment.len() / 2;
            if !tile.core.contains_point(midpoint) || position_is_blacklisted(blacklist, midpoint) {
                continue;
            }
            let model_core_index = (midpoint / config.stride) as usize;
            let support_index = model_core_index - first_model_core_index;
            fragment_support[support_index] = fragment_support[support_index]
                .checked_add(1)
                .context("fragment support overflow")?;
        }
    }
    counters.add_from_snapshot(fragments.counters_snapshot());
    coverage.finalize_coverage(true);
    apply_blacklist_mask(&mut coverage, tile, chromosome_length, blacklist)?;

    Ok(TileCoverage {
        coverage,
        counters,
        first_model_core_index,
        fragment_support,
    })
}

fn position_is_blacklisted(blacklist: &[Interval<u64>], position: u32) -> bool {
    let position = position as u64;
    let candidate = blacklist.partition_point(|interval| interval.end() <= position);
    blacklist
        .get(candidate)
        .is_some_and(|interval| interval.contains_point(position))
}

fn apply_blacklist_mask(
    coverage: &mut Coverage,
    tile: &Tile,
    chromosome_length: u64,
    blacklist: &[Interval<u64>],
) -> Result<()> {
    ensure!(
        tile.core_end() as u64 <= chromosome_length,
        "tile core exceeds chromosome length"
    );
    if blacklist.is_empty() {
        return Ok(());
    }

    let core = tile.core.try_to_u64()?;
    let first_overlap = blacklist.partition_point(|interval| interval.end() <= core.start());
    let mut local_intervals = Vec::new();
    for interval in &blacklist[first_overlap..] {
        if interval.start() >= core.end() {
            break;
        }
        if let Some(overlap) = interval.intersection(core) {
            local_intervals.push(overlap.shift_left(core.start())?);
        }
    }
    coverage.set_blacklist_mask(&local_intervals)
}

fn initialize_histograms(
    chromosomes: &[String],
    contigs: &Contigs,
    stride: u32,
) -> Result<FxHashMap<String, Vec<CoreHistogram>>> {
    let mut histograms_by_chromosome =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());
    for chromosome in chromosomes {
        let &(_, chromosome_length) = contigs
            .contigs
            .get(chromosome)
            .with_context(|| format!("missing contig metadata for '{chromosome}'"))?;
        let mut cores = Vec::new();
        let mut start = 0_u32;
        while start < chromosome_length {
            let end = start.saturating_add(stride).min(chromosome_length);
            cores.push(CoreHistogram::new(Interval::new(start, end)?));
            start = end;
        }
        histograms_by_chromosome.insert(chromosome.clone(), cores);
    }
    Ok(histograms_by_chromosome)
}

fn merge_histogram_pass_results(
    chromosomes: &[String],
    contigs: &Contigs,
    stride: u32,
    histogram_pass_tile_results: Vec<HistogramPassTileResult>,
) -> Result<(FxHashMap<String, Vec<CoreHistogram>>, FCoverageCounters)> {
    let mut histograms_by_chromosome = initialize_histograms(chromosomes, contigs, stride)?;
    let mut counters = FCoverageCounters::default();
    for tile_result in histogram_pass_tile_results {
        counters += tile_result.counters;
        let chromosome_histograms = histograms_by_chromosome
            .get_mut(&tile_result.chromosome)
            .with_context(|| {
                format!(
                    "missing histogram destination for chromosome '{}'",
                    tile_result.chromosome
                )
            })?;
        for (model_core_index, histogram) in tile_result.core_histograms {
            chromosome_histograms
                .get_mut(model_core_index)
                .with_context(|| {
                    format!(
                        "model-core histogram index {} is out of range for chromosome '{}'",
                        model_core_index, tile_result.chromosome
                    )
                })?
                .add_histogram(&histogram)?;
        }
    }
    Ok((histograms_by_chromosome, counters))
}

/// Process one parallelization tile during the outlier-calling BAM pass.
fn process_calling_pass_tile(
    config: &OutliersConfig,
    tile: &Tile,
    blacklist: &[Interval<u64>],
    chromosome_models: &[CoreModel],
) -> Result<CallingPassTileResult> {
    let tile_coverage = build_tile_coverage(config, tile, blacklist, false)?;
    let coverage = tile_coverage
        .coverage
        .coverage()
        .context("tile coverage missing during outlier calling")?;
    let mask = tile_coverage.coverage.blacklist_mask();
    let mut regions = Vec::new();
    let mut current_start: Option<u32> = None;
    let mut observed_mass = 0.0_f64;
    let mut target_mass = 0.0_f64;

    for (local_position, &coverage_value) in coverage.iter().enumerate() {
        let absolute_position = tile.core_start() + local_position as u32;
        let core_index = (absolute_position / config.stride) as usize;
        let model = chromosome_models.get(core_index).with_context(|| {
            format!(
                "missing core model {} for chromosome '{}' position {}",
                core_index, tile.chr, absolute_position
            )
        })?;
        let eligible = mask.is_none_or(|values| values[local_position] == 0);
        let rounded_coverage = coverage_value.round();
        ensure!(
            coverage_value.is_finite()
                && coverage_value >= 0.0
                && (coverage_value - rounded_coverage).abs() <= 1e-4,
            "raw coverage must be a finite nonnegative integer, got {}",
            coverage_value
        );
        let is_outlier = eligible && rounded_coverage >= model.zip.final_threshold as f32;

        if is_outlier {
            current_start.get_or_insert(absolute_position);
            observed_mass += rounded_coverage as f64;
            target_mass += match config.target {
                OutlierTarget::LocalMean => model.zip.underlying_fit.mean(),
                OutlierTarget::Threshold => model.zip.final_threshold as f64,
                OutlierTarget::Zero => 0.0,
            };
        } else if let Some(start) = current_start.take() {
            regions.push(CalledRegion {
                interval: Interval::new(start, absolute_position)?,
                observed_mass,
                target_mass,
            });
            observed_mass = 0.0;
            target_mass = 0.0;
        }
    }
    if let Some(start) = current_start {
        regions.push(CalledRegion {
            interval: Interval::new(start, tile.core_end())?,
            observed_mass,
            target_mass,
        });
    }

    Ok(CallingPassTileResult {
        chromosome: tile.chr.clone(),
        regions,
        counters: tile_coverage.counters,
    })
}

fn merge_calling_pass_results(
    chromosomes: &[String],
    calling_pass_tile_results: Vec<CallingPassTileResult>,
) -> Result<(FxHashMap<String, Vec<CalledRegion>>, FCoverageCounters)> {
    let mut regions_by_chromosome: FxHashMap<String, Vec<CalledRegion>> =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());
    let mut counters = FCoverageCounters::default();
    for chromosome in chromosomes {
        regions_by_chromosome.insert(chromosome.clone(), Vec::new());
    }

    for tile_result in calling_pass_tile_results {
        counters += tile_result.counters;
        let destination = regions_by_chromosome
            .get_mut(&tile_result.chromosome)
            .with_context(|| {
                format!(
                    "missing outlier-region destination for chromosome '{}'",
                    tile_result.chromosome
                )
            })?;
        for region in tile_result.regions {
            if let Some(previous) = destination.last_mut() {
                ensure!(
                    previous.interval.end() <= region.interval.start(),
                    "outlier tile regions overlap unexpectedly on chromosome '{}'",
                    tile_result.chromosome
                );
                if previous.interval.end() == region.interval.start() {
                    previous.interval = previous.interval.expand_to_include(region.interval);
                    previous.observed_mass += region.observed_mass;
                    previous.target_mass += region.target_mass;
                    continue;
                }
            }
            destination.push(region);
        }
    }
    Ok((regions_by_chromosome, counters))
}

#[cfg(test)]
mod tests {
    include!("outliers_tests.rs");
}
