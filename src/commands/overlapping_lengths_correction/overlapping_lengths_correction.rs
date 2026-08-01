//! Command orchestration for tiled statistics collection, model fitting, and staged output.

use std::{path::PathBuf, sync::Arc, time::Instant};

use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use rayon::prelude::*;
use tracing::{info, warn};

use crate::{
    command_run::{CommandRunResult, RunOptions, status_info},
    commands::{
        cli_common::{
            ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs, validate_output_prefix,
        },
        counters::FCoverageCounters,
        gc_bias::correct::load_gc_corrector,
        run_statistics::{
            DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions, GCStatisticsSummary,
            TILE_DOUBLE_COUNT_NOTE, print_fragment_run_statistics,
        },
    },
    shared::{
        interval::Interval,
        io::{FinalOutputFiles, dot_join},
        progress::ProgressFactory,
        reference::{ReferenceReader, stage_reference_2bit},
        scale_genome::{ScalingBin, scaling_gc_mode_for_run},
        thread_pool::init_global_pool,
        tiled_run::{RunTempDirs, build_tiles},
    },
};

use super::{
    config::{
        DEFAULT_MAX_FRAGMENT_LENGTH, DEFAULT_MIN_FRAGMENT_LENGTH,
        OverlappingLengthsCorrectionConfig,
    },
    model::fit_overlapping_length_model,
    package::OverlappingLengthsCorrectionPackage,
    reducer::{OverlappingLengthStatistics, build_length_bin_edges},
    tiling::{OverlappingLengthTileResult, process_tile},
};

const COMMAND_TARGET: &str = "overlap-length-model";

/// Result of fitting an average overlapping fragment length normalization model.
#[derive(Debug)]
pub struct OverlappingLengthsCorrectionRunResult {
    /// Fragment and GC-filtering counters accumulated across all tiles.
    pub counters: FCoverageCounters,
    /// Final path of the written Zarr model package.
    pub package_path: PathBuf,
    /// Number of covered, non-blacklisted genomic bases used for fitting.
    pub eligible_covered_bases: u64,
    /// Complete output list used by the shared command-run interface.
    output_files: Vec<PathBuf>,
}

impl CommandRunResult for OverlappingLengthsCorrectionRunResult {
    type Counters = FCoverageCounters;

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&std::path::Path> {
        Some(&self.package_path)
    }
}

/// Fit and write the sample's average overlapping fragment length normalization model.
///
/// This is the programmatic command entry point. It performs a single tiled BAM sweep, merges the
/// additive per-tile statistics, performs both LIONHEART model fits in memory, and writes a
/// self-contained Zarr package. The second fit does not reread the BAM.
///
/// Parameters
/// ----------
/// - `config`:
///   Input, filtering, normalization model, tiling, and output settings.
/// - `options`:
///   Shared controls for progress reporting, equivalent-command logging, and statistics output.
///
/// Returns
/// -------
/// - `OverlappingLengthsCorrectionRunResult`:
///   Counters, the final package path, and the number of bases used for fitting.
pub fn run_overlapping_lengths_correction(
    config: &OverlappingLengthsCorrectionConfig,
    options: RunOptions,
) -> Result<OverlappingLengthsCorrectionRunResult> {
    let start_time = Instant::now();
    let result = execute(config, options)?;
    if options.report_statistics {
        print_fragment_run_statistics(
            &result.counters.base,
            start_time.elapsed(),
            FragmentRunStatisticsOptions {
                include_section_header: true,
                notes: &[TILE_DOUBLE_COUNT_NOTE],
                labels: DEFAULT_FRAGMENT_STATISTICS_LABELS,
                blacklist_excluded_fragments: None,
                gc: (config.gc.gc_file.is_some() || config.gc.gc_tag.is_some()).then_some(
                    GCStatisticsSummary {
                        failed_fragments: result.counters.gc_failed_fragments,
                        neutralize_invalid_gc: config.gc.neutralize_invalid_gc,
                        missing_tags: config
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(result.counters.gc_missing_tags),
                        out_of_range_tags: config
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(result.counters.gc_out_of_range_tags),
                    },
                ),
            },
            [format!(
                "Eligible covered bases used for fitting: {}",
                result.eligible_covered_bases
            )],
        );
    }
    Ok(result)
}

/// Validate inputs, run parallel tile collection, fit the model, and publish the package.
///
/// Keeping orchestration here leaves the public entry point responsible only for timing and
/// optional statistics reporting.
fn execute(
    config: &OverlappingLengthsCorrectionConfig,
    options: RunOptions,
) -> Result<OverlappingLengthsCorrectionRunResult> {
    // Validate scientific combinations before creating output or temporary directories
    let fragment_lengths = config.fragment_lengths();
    fragment_lengths.validate()?;
    config.gc.validate(config.ref_2bit.as_deref())?;
    ensure!(
        config.length_bin_size > 0,
        "--length-bin-size must be positive"
    );
    if config.unpaired.reads_are_fragments && config.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    if config.unpaired.reads_are_fragments && config.ignore_gap {
        bail!("--ignore-gap cannot be used with --reads-are-fragments");
    }
    validate_output_prefix(config.output_prefix.trim())?;
    if config.min_fragment_length != DEFAULT_MIN_FRAGMENT_LENGTH
        || config.max_fragment_length != DEFAULT_MAX_FRAGMENT_LENGTH
    {
        warn!(
            target: COMMAND_TARGET,
            "The overlapping fragment length normalization model was developed primarily for 100-220 bp fragments; the configured {}-{} bp range is experimental",
            config.min_fragment_length,
            config.max_fragment_length
        );
    }
    if options.log_equivalent_cli {
        let command = crate::ToCliCommand::to_cli_string(config)?;
        info!(target: COMMAND_TARGET, "{}", crate::command_run::equivalent_cli_log_message(&command));
    }

    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&config.chromosomes, config.ioc.bam.as_path())?;

    // Stage outputs so an interrupted fit cannot leave a package that appears complete
    ensure_output_dir(&config.ioc.output_dir)?;
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
    .context("create overlapping fragment length model temporary directories")?;
    let mut final_outputs = FinalOutputFiles::new(run_temp_dirs.final_output_dir())?;

    // Each Rayon worker opens its own reader against this local staged reference
    let staged_reference = if config.gc.gc_file.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Copying reference 2bit to the temporary directory");
        Some(stage_reference_2bit(
            config
                .ref_2bit
                .as_deref()
                .context("--ref-2bit is required with --gc-file")?,
            run_temp_dirs.work_dir(),
        )?)
    } else {
        None
    };
    if config.blacklist.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Loading blacklists");
    }
    // Blacklists and scaling tracks are chromosome-indexed once and borrowed by tile workers
    let blacklist_map = load_blacklist_map(
        config.blacklist.as_ref(),
        1,
        0,
        &chromosomes,
        config.ioc.n_threads > 1,
    )?;
    if config.scale_genome.scaling_factors.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<ScalingBin>> = load_scaling_map(
        &config.scale_genome,
        &chromosomes,
        &contigs,
        scaling_gc_mode_for_run(config.gc.gc_file.is_some(), config.gc.gc_tag.is_some()),
        Some(config.ignore_gap),
    )?;
    if config.gc.gc_file.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        config.gc.gc_file.as_ref(),
        staged_reference.as_ref(),
        config.min_fragment_length,
        config.max_fragment_length,
    )?;
    // Bin edges are constructed once so every tile produces merge-compatible vectors
    let length_bin_edges = build_length_bin_edges(
        config.min_fragment_length,
        config.max_fragment_length,
        config.length_bin_size,
    )?;

    init_global_pool(config.ioc.n_threads)?;
    let (tiles, _) = build_tiles(
        &chromosomes,
        &contigs,
        config.tile_size,
        config.max_fragment_length,
        None,
    )?;
    let progress = ProgressFactory::with_enabled(options.show_progress);
    let progress_bar = Arc::new(progress.default_bar(tiles.len() as u64));
    status_info!(options, target: COMMAND_TARGET, "Collecting overlapping fragment length statistics per tile");
    let gc_tag = config.gc.gc_tag.as_deref();
    // A maximum-fragment-length halo supplies complete fragments at every tile-core boundary
    let tile_results = tiles
        .par_iter()
        .map_init(
            || staged_reference.as_deref().map(ReferenceReader::open),
            |reference_reader_result, tile| -> Result<OverlappingLengthTileResult> {
                let reference_reader = match reference_reader_result {
                    Some(reader_result) => Some(
                        reader_result
                            .as_mut()
                            .map_err(|error| anyhow::anyhow!("{error:#}"))?,
                    ),
                    None => None,
                };
                let blacklist = blacklist_map
                    .get(&tile.chr)
                    .map(Vec::as_slice)
                    .unwrap_or(&[] as &[Interval<u64>]);
                let scaling = scaling_map
                    .get(&tile.chr)
                    .map(Vec::as_slice)
                    .unwrap_or(&[] as &[ScalingBin]);
                let result = process_tile(
                    config,
                    tile,
                    length_bin_edges.clone(),
                    blacklist,
                    scaling,
                    gc_corrector.clone(),
                    gc_tag,
                    reference_reader,
                )?;
                progress_bar.inc(1);
                Ok(result)
            },
        )
        .collect::<Result<Vec<_>>>()?;
    progress_bar.finish_with_message("| Finished collecting statistics");

    // Only additive statistics cross thread boundaries, never per-position arrays
    let mut statistics = OverlappingLengthStatistics::new(length_bin_edges)?;
    let mut counters = FCoverageCounters::default();
    for tile_result in tile_results {
        statistics.merge(tile_result.statistics)?;
        counters += tile_result.counters;
    }
    status_info!(options, target: COMMAND_TARGET, "Fitting the two overlapping fragment length mixture models");
    // Both optimization passes operate on the merged sufficient statistics in memory
    let model = fit_overlapping_length_model(&statistics)?;
    let package = OverlappingLengthsCorrectionPackage::from_model(config, &statistics, model);

    let output_name = dot_join(&[config.output_prefix.trim(), "overlap_length_model.zarr"]);
    let final_path = config.ioc.output_dir.join(&output_name);
    let temporary_path = final_outputs.temp_path_for(&final_path)?;
    // Publish the complete directory only after every Zarr array and attribute has been written
    package.write_zarr(&temporary_path)?;
    final_outputs.record(temporary_path, final_path.clone())?;
    final_outputs.move_into_place()?;

    Ok(OverlappingLengthsCorrectionRunResult {
        counters,
        package_path: final_path.clone(),
        eligible_covered_bases: statistics.eligible_covered_bases,
        output_files: vec![final_path],
    })
}
