use crate::command_run::{CommandRunResult, RunOptions};
use crate::commands::fcoverage::config::COVERAGE_SIGNAL_LABEL;
use crate::commands::fcoverage::reducer::TileAggregateTempFiles;
use crate::commands::fcoverage::tiling::{
    TileTempFile, TileTempFileKind, adapt_fetch_to_extreme_windows, finalize_value,
    merge_positional_tile_outputs_with_optional_scaling,
};
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use crate::commands::fcoverage::writers::{
    write_bed_aggregate_output, write_bedgraph_runs, write_final_row, write_size_aggregate_output,
    write_windowed_runs,
};
use crate::commands::gc_bias::correct::{GCCorrector, load_gc_corrector};
use crate::commands::gc_bias::counting::build_gc_prefixes;
use crate::commands::wps::config::{WPSConfig, WPSSharedConfig};
use crate::shared::formatters::round_to;
use crate::shared::fragment::minimal_fragment::Fragment;
use crate::shared::fragment_iterators::fragments_from_bam;
use crate::shared::gc_tag::ClassifiedGCTagWeight;
use crate::shared::interval::{IndexedInterval, Interval};
use crate::shared::io::{FinalOutputFiles, dot_join};
use crate::shared::progress::ProgressFactory;
use crate::shared::read::{default_include_read_paired_end, default_include_read_unpaired};
use crate::shared::reference::read_seq_in_range;
use crate::shared::scale_genome::{ScalingBin, apply_scaling_to_coverage_in_place};
use crate::shared::temp_chrom_names::TempChromNameMap;
use crate::shared::tiled_run::{
    TempDirGuard, Tile, TileMode, TileWindowSpan, build_tiles, overlapping_windows_for_tile,
    precompute_tile_window_spans,
};
use crate::shared::writers::open_zstd_auto_writer;
use crate::{
    commands::cli_common::{
        WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
        resolve_chromosomes_and_contigs, validate_output_prefix,
    },
    commands::counters::WPSCounters,
    commands::run_statistics::{
        DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions, GCStatisticsSummary,
        TILE_DOUBLE_COUNT_NOTE, print_fragment_run_statistics,
    },
    shared::{
        bam::create_chromosome_reader, bed::load_windows_from_bed, thread_pool::init_global_pool,
        windowing::ensure_plain_bed_windows_not_empty,
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::io::Write;
use std::{path::PathBuf, sync::Arc, time::Instant};
use tracing::info;

const COMMAND_TARGET: &str = "wps";

/// Result from `wps`.
///
/// The command writes positional or windowed protection scores. The result records the
/// primary output path, all final output files, and the fragment counters from WPS calculation.
#[derive(Debug)]
pub struct WPSRunResult {
    /// Fragment and filtering counters collected during the run.
    pub counters: WPSCounters,
    /// Main WPS output path.
    pub output_path: PathBuf,
    /// Final output files produced by the command.
    pub output_files: Vec<PathBuf>,
}

impl CommandRunResult for WPSRunResult {
    type Counters = WPSCounters;

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&std::path::Path> {
        Some(self.output_path.as_path())
    }
}

#[derive(Debug, Clone)]
struct WpsTileResult {
    counters: WPSCounters,
    temp_output: Option<WpsTileTempOutput>,
}

#[derive(Debug, Clone)]
pub(crate) enum WpsTileTempOutput {
    Positional {
        chromosome: String,
        tile_index: u32,
        path: PathBuf,
    },
    AggregatesByBed {
        chromosome: String,
        tile_index: u32,
        partials_path: PathBuf,
        cross_index_path: Option<PathBuf>,
    },
    AggregatesBySize {
        chromosome: String,
        tile_index: u32,
        partials_path: PathBuf,
        cross_index_path: Option<PathBuf>,
    },
    SizeFinal {
        chromosome: String,
        tile_index: u32,
        path: PathBuf,
    },
}

impl WpsTileTempOutput {
    fn is_bed_aggregate(&self) -> bool {
        matches!(self, Self::AggregatesByBed { .. })
    }

    fn is_size_aggregate(&self) -> bool {
        matches!(self, Self::AggregatesBySize { .. })
    }
}

/// Run the `wps` command.
///
/// This command computes windowed protection scores from fragment coverage. It resolves
/// chromosomes, prepares optional windows, blacklists, and scaling data, processes tiles in
/// parallel, applies smoothing, and writes positional or aggregated WPS outputs.
///
/// Reporting is controlled by `options`. `report_statistics` prints the final summary and
/// `show_progress` controls progress bars. This command does not use status logs.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `wps` command.
/// - `options`:
///     Reporting controls for statistics and progress bars.
///
/// Returns
/// -------
/// - `Ok(WPSRunResult)`:
///     Counters and output paths for the completed run.
///
/// Errors
/// ------
/// Returns an error if the BAM cannot be read, auxiliary files are invalid, or writing outputs
/// fails at any stage.
pub fn run_wps(opt: &WPSConfig, options: RunOptions) -> Result<WPSRunResult> {
    let start_time = Instant::now();
    if opt.shared_args.unpaired.reads_are_fragments && opt.shared_args.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let (chromosomes, contigs) = resolve_chromosomes_and_contigs(
        &opt.shared_args.chromosomes,
        &opt.shared_args.ioc.bam.as_path(),
    )?;
    let prefix = opt.shared_args.output_prefix.trim();
    validate_output_prefix(prefix)?;
    let window_opt = opt.shared_args.windows.resolve_windows();
    let windowed = matches!(window_opt, WindowSpec::Bed(_) | WindowSpec::Size(_));

    if windowed {
        ensure!(
            opt.per_window.is_some(),
            "when using --by-bed/--by-size, please also specify --per-window"
        );
    }

    ensure!(
        opt.shared_args.min_fragment_length >= opt.shared_args.window_size,
        "min-fragment-length ({}) must be >= window-size ({})",
        opt.shared_args.min_fragment_length,
        opt.shared_args.window_size
    );

    ensure!(
        opt.shared_args.window_size <= opt.shared_args.max_fragment_length,
        "window-size ({}) must be <= max-fragment-length ({})",
        opt.shared_args.window_size,
        opt.shared_args.max_fragment_length
    );
    opt.shared_args
        .gc
        .validate(opt.shared_args.ref_2bit.as_deref())?;

    let per_window_wps_action = opt.per_window;

    if let Some(action) = per_window_wps_action {
        ensure!(
            matches!(
                action,
                CoverageWindowAction::Average
                    | CoverageWindowAction::Total
                    | CoverageWindowAction::OnlyIncludeThesePositionsUnique
                    | CoverageWindowAction::OnlyIncludeThesePositionsIndexed
            ),
            "for WPS, --per-window can only be 'average', 'total', 'unique-positions', or 'indexed-positions'"
        );
    }

    // Create output directory
    ensure_output_dir(&opt.shared_args.ioc.output_dir)?;

    if opt.shared_args.blacklist.is_some() {
        info!(target: COMMAND_TARGET, "Loading blacklists");
    }
    // We don't want WPS scores that were biased from neighbouring blacklisted regions
    // So we don't use positions where any fragments could also touch a blacklisted region
    let blacklist_halo =
        (opt.shared_args.max_fragment_length + (opt.shared_args.window_size + 1) / 2) as u64;
    let blacklist_map = load_blacklist_map(
        opt.shared_args.blacklist.as_ref(),
        1,
        blacklist_halo,
        &chromosomes,
    )?;

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            info!(target: COMMAND_TARGET, "Loading window coordinates");
            let wds = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            ensure_plain_bed_windows_not_empty(&wds)?;
            if matches!(
                per_window_wps_action,
                Some(CoverageWindowAction::OnlyIncludeThesePositionsUnique)
            ) {
                // Merge in-place to avoid double memory-usage
                info!(
                    target: COMMAND_TARGET,
                    "Merging overlapping/touching windows"
                );
                // Take ownership so we can remove entries by chromosome
                let mut wds_owned: FxHashMap<String, crate::shared::bed::Windows> = wds;
                let mut out: FxHashMap<String, crate::shared::bed::Windows> =
                    FxHashMap::with_capacity_and_hasher(wds_owned.len(), Default::default());
                let mut next_idx: u64 = 0;

                // Use the user-provided `chromosomes` order to assign indices deterministically
                for chr in &chromosomes {
                    if let Some(ws) = wds_owned.remove(chr) {
                        // Flatten in-place
                        let (flat, next) = ws.into_flattened_reindexed(next_idx);
                        next_idx = next;
                        out.insert(chr.clone(), flat);
                    }
                }
                Some(out)
            } else {
                Some(wds)
            }
        }
        _ => None,
    };

    // Load genomic scaling factors
    if opt.shared_args.scale_genome.scaling_factors.is_some() {
        info!(target: COMMAND_TARGET, "Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<ScalingBin>> = load_scaling_map(
        &opt.shared_args.scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(
            opt.shared_args.gc.gc_file.is_some(),
            opt.shared_args.gc.gc_tag.is_some(),
        ),
        None,
    )?;

    // Load GC correction package if specified
    if opt.shared_args.gc.gc_file.is_some() {
        info!(target: COMMAND_TARGET, "Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.shared_args.gc.gc_file.as_ref(),
        opt.shared_args.ref_2bit.as_ref(),
        opt.shared_args.min_fragment_length,
        opt.shared_args.max_fragment_length,
    )?;

    let has_scaling_or_correction = opt.shared_args.scale_genome.scaling_factors.is_some()
        || opt.shared_args.gc.gc_file.is_some();

    // Some actions cannot be used with `--by-size`
    if matches!(window_opt, WindowSpec::Size(_))
        && matches!(
            per_window_wps_action,
            Some(CoverageWindowAction::OnlyIncludeThesePositionsUnique)
                | Some(CoverageWindowAction::OnlyIncludeThesePositionsIndexed)
        )
    {
        anyhow::bail!("in --by-size mode, --per-window can only be 'average' or 'total'");
    }

    // Build temporary directory
    let temp_dir_guard = TempDirGuard::new(&opt.shared_args.ioc.output_dir, prefix)
        .context("create per-run temp dir")?;
    let temp_dir = temp_dir_guard.path();
    let mut final_outputs = FinalOutputFiles::new(temp_dir)?;

    // Window size when --by-size (otherwise None)
    let by_size_bp: Option<u64> = match &window_opt {
        WindowSpec::Size(bp) => Some(*bp as u64),
        _ => None,
    };

    // Build tiles
    let halo_bp: u32 = opt
        .shared_args
        .max_fragment_length
        .saturating_add(opt.shared_args.window_size); // Extend fetch so windows see complete fragments
    let (tiles, tile_and_window_boundaries_align) = build_tiles(
        &chromosomes,
        &contigs,
        opt.shared_args.tile_size,
        halo_bp,
        by_size_bp,
    )?;
    let temp_chrom_name_map = TempChromNameMap::from_contigs(&chromosomes)?;

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(
        &tiles,
        |chr| {
            windows_lookup
                .and_then(|m| m.get(chr).map(|w| w.as_slice()))
                .unwrap_or(&[])
        },
        0,
        0,
    ));

    // Where per-tile files go
    let positional_prefix = dot_join(&[prefix, "pos"]);
    let partials_prefix = dot_join(&[prefix, "part"]);
    let finals_prefix = dot_join(&[prefix, "fin"]);

    // Faster to convert to &str once
    let positional_prefix = positional_prefix.as_str();
    let partials_prefix = partials_prefix.as_str();
    let finals_prefix = finals_prefix.as_str();

    // Create filenames of final outputs
    let final_bedgraph_pos_name = dot_join(&[prefix, "wps.per_position.bedgraph.zst"]);
    let final_tsv_pos_name = dot_join(&[prefix, "wps.per_position_per_window.tsv.zst"]);
    let final_average_name = dot_join(&[prefix, "wps.average.tsv.zst"]);
    let final_total_name = dot_join(&[prefix, "wps.total.tsv.zst"]);

    let per_window_action = if windowed {
        Some(per_window_wps_action.context("windowed WPS runs require a per-window action")?)
    } else {
        None
    };

    // Get decimals to use
    let decimals_to_use: i32 = match per_window_action {
        Some(action) => match action {
            CoverageWindowAction::Average | CoverageWindowAction::Total => {
                opt.shared_args.decimals as i32
            }
            CoverageWindowAction::OnlyIncludeThesePositionsUnique
            | CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                if has_scaling_or_correction {
                    opt.shared_args.decimals as i32
                } else {
                    0
                }
            }
            _ => {
                unreachable!("unsupported WPS per-window action must be rejected during validation")
            }
        },
        None => {
            if has_scaling_or_correction {
                opt.shared_args.decimals as i32
            } else {
                0
            }
        }
    };

    let total_tiles = tiles.len();

    // Create progress bar
    let progress = ProgressFactory::with_enabled(options.show_progress);
    let pb = Arc::new(progress.default_bar(total_tiles as u64));

    // Configure global thread‐pool size
    init_global_pool(opt.shared_args.ioc.n_threads as usize)?;

    let mut global_counter = WPSCounters::default();

    info!(target: COMMAND_TARGET, "Calculating WPS per tile");
    pb.set_position(0);

    let tile_window_spans_for_threads = tile_window_spans.clone();

    let tile_results: Vec<WpsTileResult> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<WpsTileResult> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            // Per-chrom projections
            let windows_chr: Option<&[IndexedInterval<u64>]> = windows_map
                .as_ref()
                .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice()));
            let blacklist_chr: &[Interval<u64>] = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr: &[ScalingBin] = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            let tile_result = if matches!(window_opt, WindowSpec::Bed(_))
                && windows_chr.map_or(true, |windows| windows.is_empty())
            {
                WpsTileResult {
                    counters: WPSCounters::default(),
                    temp_output: None,
                }
            } else {
                // Decide tile mode and filename
                let (action_prefix, extensions) = if windowed {
                    let action = per_window_action
                        .context("windowed WPS runs require a per-window action")?;
                    match action {
                        // We need
                        CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                            (positional_prefix, "tsv.zst")
                        }
                        CoverageWindowAction::OnlyIncludeThesePositionsUnique => {
                            (positional_prefix, "bedgraph.zst")
                        }
                        CoverageWindowAction::Average | CoverageWindowAction::Total => {
                            (partials_prefix, "tsv.zst")
                        }
                        _ => unreachable!(
                            "unsupported WPS per-window action must be rejected during validation"
                        ),
                    }
                } else {
                    // Whole positional coverage
                    (positional_prefix, "bedgraph.zst")
                };

                let chr_token = temp_chrom_name_map.token_for(tile.chr.as_str())?;
                let out_path = temp_dir.join(format!(
                    "{prefix}.{chr}.{idx}.{extensions}",
                    prefix = action_prefix,
                    chr = chr_token,
                    idx = tile.index
                ));

                // Windowed tmp outputs for faster reducer later on
                let partials_out = temp_dir.join(format!(
                    "{prefix}.{chr}.{idx}.{extensions}",
                    prefix = partials_prefix,
                    chr = chr_token,
                    idx = tile.index
                ));
                let cross_idx_out = temp_dir.join(format!(
                    "{prefix}.cross.{chr}.{idx}.cross.zst", // Needs this extension!
                    prefix = partials_prefix,
                    chr = chr_token,
                    idx = tile.index
                ));
                let finals_out = temp_dir.join(format!(
                    "{prefix}.{chr}.{idx}.{extensions}",
                    prefix = finals_prefix,
                    chr = chr_token,
                    idx = tile.index
                ));

                let mode = if !windowed {
                    TileMode::Positional {
                        windows: None,
                        out_path,
                        indexed: false,
                    }
                } else {
                    match (
                        &window_opt,
                        per_window_action
                            .context("windowed WPS runs require a per-window action")?,
                    ) {
                        (
                            WindowSpec::Bed(_),
                            CoverageWindowAction::OnlyIncludeThesePositionsUnique,
                        ) => TileMode::Positional {
                            windows: windows_chr,
                            out_path,
                            indexed: false,
                        },
                        (
                            WindowSpec::Bed(_),
                            CoverageWindowAction::OnlyIncludeThesePositionsIndexed,
                        ) => TileMode::Positional {
                            windows: windows_chr,
                            out_path,
                            indexed: true,
                        },
                        (
                            WindowSpec::Bed(_),
                            CoverageWindowAction::Average | CoverageWindowAction::Total,
                        ) => TileMode::AggregatesByBed {
                            windows: windows_chr.context(
                                "BED aggregate tile reached processing without any windows",
                            )?,
                            masked: true,
                            partials_out,
                            cross_idx_out,
                        },
                        (
                            WindowSpec::Size(size),
                            CoverageWindowAction::Average | CoverageWindowAction::Total,
                        ) => TileMode::AggregatesBySize {
                            window_bp: *size,
                            masked: true,
                            finals_out,
                            partials_out,
                            cross_idx_out,
                            guaranteed_aligned: tile_and_window_boundaries_align,
                        },
                        _ => {
                            anyhow::bail!(
                                "Got illegal combination of --by-size/--by-bed and --per-window."
                            )
                        }
                    }
                };

                let (counters, temp_output, _, _) = wps_for_tile(
                    &opt.shared_args,
                    &opt.per_window,
                    opt.keep_zero_runs,
                    tile,
                    tile_span.as_ref(),
                    blacklist_chr,
                    scaling_chr,
                    gc_corrector.clone(),
                    mode,
                    decimals_to_use,
                    0,
                    false,
                )?;
                WpsTileResult {
                    counters,
                    temp_output,
                }
            };
            pb.inc(1);
            Ok(tile_result)
        })
        .collect::<anyhow::Result<_>>()?;

    if options.show_progress {
        pb.finish_with_message("| Finished counting");
    } else {
        pb.finish_and_clear();
    }

    // Collect counters and the temp paths returned by tile processing.
    let mut tile_temp_outputs = Vec::new();
    for tile_result in tile_results {
        global_counter += tile_result.counters;
        if let Some(temp_output) = tile_result.temp_output {
            tile_temp_outputs.push(temp_output);
        }
    }

    info!(
        target: COMMAND_TARGET,
        "Merging temporary tile files to final output"
    );

    // Merge tile temp files into a completed output under the final-output temp directory
    // The merged file moves into output_dir only after the full writer succeeds

    let final_out_path = if let Some(action) = per_window_action {
        match action {
            CoverageWindowAction::OnlyIncludeThesePositionsUnique => {
                // Windowed positional (unique and non-indexed)
                let final_path = opt
                    .shared_args
                    .ioc
                    .output_dir
                    .join(&final_bedgraph_pos_name);
                let temp_path = merge_positional_tile_outputs_with_optional_scaling(
                    final_outputs.temp_dir(),
                    &chromosomes,
                    &wps_positional_tile_outputs(&tile_temp_outputs),
                    final_bedgraph_pos_name.as_str(),
                    None,
                    false,
                    decimals_to_use,
                    opt.shared_args.ioc.n_threads as usize,
                )?;
                final_outputs.record(temp_path, final_path.clone())?;
                final_path
            }
            CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                // Windowed positional with orig_idx column
                let final_path = opt.shared_args.ioc.output_dir.join(&final_tsv_pos_name);
                let temp_path = merge_positional_tile_outputs_with_optional_scaling(
                    final_outputs.temp_dir(),
                    &chromosomes,
                    &wps_positional_tile_outputs(&tile_temp_outputs),
                    final_tsv_pos_name.as_str(),
                    None,
                    true,
                    decimals_to_use,
                    opt.shared_args.ioc.n_threads as usize,
                )?;
                final_outputs.record(temp_path, final_path.clone())?;
                final_path
            }
            CoverageWindowAction::Average | CoverageWindowAction::Total => {
                // Per-chrom reduce of partials into final aggregates
                let final_path = opt.shared_args.ioc.output_dir.join(match action {
                    CoverageWindowAction::Average => final_average_name.as_str(),
                    CoverageWindowAction::Total => final_total_name.as_str(),
                    _ => unreachable!(),
                });
                let temp_final_path = final_outputs.temp_path_for(&final_path)?;

                // Reduce by window source
                match &window_opt {
                    WindowSpec::Bed(_) => {
                        let win_map = windows_map
                            .as_ref()
                            .context("BED WPS reduction requires loaded windows")?;
                        write_bed_aggregate_output(
                            &temp_final_path,
                            &collect_wps_aggregate_tile_outputs_by_chromosome(
                                &tile_temp_outputs,
                                WpsTileTempOutput::is_bed_aggregate,
                            )?,
                            win_map,
                            &chromosomes,
                            true,
                            action,
                            decimals_to_use,
                            opt.shared_args.ioc.n_threads as usize,
                            COVERAGE_SIGNAL_LABEL,
                            None,
                        )?;
                    }
                    WindowSpec::Size(_) => {
                        write_size_aggregate_output(
                            &temp_final_path,
                            &collect_wps_aggregate_tile_outputs_by_chromosome(
                                &tile_temp_outputs,
                                WpsTileTempOutput::is_size_aggregate,
                            )?,
                            &wps_size_final_tile_outputs(&tile_temp_outputs),
                            &chromosomes,
                            &contigs,
                            true,
                            action,
                            decimals_to_use,
                            opt.shared_args.ioc.n_threads as usize,
                            tile_and_window_boundaries_align,
                            COVERAGE_SIGNAL_LABEL,
                            None,
                        )?;
                    }
                    _ => unreachable!(),
                }

                final_outputs.record(temp_final_path, final_path.clone())?;
                final_path
            }
            _ => {
                unreachable!("unsupported WPS per-window action must be rejected during validation")
            }
        }
    } else {
        // Whole-genome positional coverage
        let final_path = opt
            .shared_args
            .ioc
            .output_dir
            .join(&final_bedgraph_pos_name);
        let temp_path = merge_positional_tile_outputs_with_optional_scaling(
            final_outputs.temp_dir(),
            &chromosomes,
            &wps_positional_tile_outputs(&tile_temp_outputs),
            final_bedgraph_pos_name.as_str(),
            None,
            false,
            decimals_to_use,
            opt.shared_args.ioc.n_threads as usize,
        )?;
        final_outputs.record(temp_path, final_path.clone())?;
        final_path
    };

    final_outputs.move_into_place()?;

    info!(
        target: COMMAND_TARGET,
        "Saved output to: {}",
        final_out_path.display()
    );

    let elapsed = start_time.elapsed();
    if options.report_statistics {
        print_fragment_run_statistics(
            &global_counter.base,
            elapsed,
            FragmentRunStatisticsOptions {
                include_section_header: true,
                notes: &[TILE_DOUBLE_COUNT_NOTE],
                labels: DEFAULT_FRAGMENT_STATISTICS_LABELS,
                blacklist_excluded_fragments: None,
                gc: (opt.shared_args.gc.gc_file.is_some() || opt.shared_args.gc.gc_tag.is_some())
                    .then_some(GCStatisticsSummary {
                        neutralize_invalid_gc: opt.shared_args.gc.neutralize_invalid_gc,
                        failed_fragments: global_counter.gc_failed_fragments,
                        missing_tags: opt
                            .shared_args
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(global_counter.gc_missing_tags),
                        out_of_range_tags: opt
                            .shared_args
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(global_counter.gc_out_of_range_tags),
                    }),
            },
            std::iter::empty::<&str>(),
        );
    }
    Ok(WPSRunResult {
        counters: global_counter,
        output_path: final_out_path.clone(),
        output_files: vec![final_out_path],
    })
}

fn wps_positional_tile_outputs(tile_outputs: &[WpsTileTempOutput]) -> Vec<TileTempFile> {
    tile_outputs
        .iter()
        .filter_map(|tile_output| match tile_output {
            WpsTileTempOutput::Positional {
                chromosome,
                tile_index,
                path,
            } => Some(TileTempFile {
                kind: TileTempFileKind::Positional,
                chromosome: chromosome.clone(),
                tile_index: *tile_index,
                path: path.clone(),
            }),
            _ => None,
        })
        .collect()
}

fn wps_size_final_tile_outputs(tile_outputs: &[WpsTileTempOutput]) -> Vec<TileTempFile> {
    tile_outputs
        .iter()
        .filter_map(|tile_output| match tile_output {
            WpsTileTempOutput::SizeFinal {
                chromosome,
                tile_index,
                path,
            } => Some(TileTempFile {
                kind: TileTempFileKind::SizeFinal,
                chromosome: chromosome.clone(),
                tile_index: *tile_index,
                path: path.clone(),
            }),
            _ => None,
        })
        .collect()
}

fn collect_wps_aggregate_tile_outputs_by_chromosome(
    tile_outputs: &[WpsTileTempOutput],
    include_output: fn(&WpsTileTempOutput) -> bool,
) -> Result<FxHashMap<String, Vec<TileAggregateTempFiles>>> {
    let mut by_chromosome: FxHashMap<String, Vec<TileAggregateTempFiles>> = FxHashMap::default();

    for tile_output in tile_outputs.iter().filter(|output| include_output(output)) {
        match tile_output {
            WpsTileTempOutput::AggregatesByBed {
                chromosome,
                tile_index,
                partials_path,
                cross_index_path,
            }
            | WpsTileTempOutput::AggregatesBySize {
                chromosome,
                tile_index,
                partials_path,
                cross_index_path,
            } => {
                by_chromosome
                    .entry(chromosome.clone())
                    .or_default()
                    .push(TileAggregateTempFiles {
                        tile_index: *tile_index,
                        partials_path: partials_path.clone(),
                        cross_index_path: cross_index_path.clone(),
                    })
            }
            WpsTileTempOutput::Positional { .. } | WpsTileTempOutput::SizeFinal { .. } => {}
        }
    }

    for (chromosome, outputs) in &by_chromosome {
        let mut tile_indices = fxhash::FxHashSet::default();
        for output in outputs {
            anyhow::ensure!(
                tile_indices.insert(output.tile_index),
                "duplicate WPS aggregate tile index {} for chromosome '{}'",
                output.tile_index,
                chromosome
            );
        }
    }

    Ok(by_chromosome)
}

/// Process one tile: pair reads, build coverage, and write/return outputs for this tile
pub(crate) fn wps_for_tile(
    opt: &WPSSharedConfig,
    per_window_wps_action: &Option<CoverageWindowAction>,
    keep_zero_runs: bool,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_chr: &[Interval<u64>],
    scaling_chr: &[ScalingBin],
    gc_corrector_opt: Option<GCCorrector>,
    mode: TileMode,
    decimals: i32,
    extra_halo_bp: u32,
    return_wps_instead: bool, // Don't save and aggregate, just return the WPS values
) -> Result<(
    WPSCounters,
    Option<WpsTileTempOutput>,
    Option<Vec<f32>>,
    Option<Vec<u8>>,
)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    tile.ensure_matches_bam_tid(tid_check)?;

    let mut counter = WPSCounters::default();

    let gc_prefixes_opt = if gc_corrector_opt.is_some() {
        let ref_2bit = match opt.ref_2bit.as_ref() {
            Some(r) => r,
            None => bail!("When GC correction is specified, --ref-2bit must also be specified"),
        };
        let seq_bytes = read_seq_in_range(
            &ref_2bit,
            &tile.chr,
            // NOTE: Need for full fetch span to get GC of overlapping fragments!
            (tile.fetch_start() as usize)..(tile.fetch_end() as usize),
        )?;
        Some(build_gc_prefixes(&seq_bytes))
    } else {
        None
    };

    // Adapt the fetch coordinates to the present windows (*in genomic-windowed mode!*)
    // When no windows are present, skip this tile
    let Some(fetch_span) = adapt_fetch_to_extreme_windows(
        tile,
        tile_window_span,
        &mode,
        chrom_len as u32,
        opt.max_fragment_length as u64,
    )?
    else {
        return Ok((counter, None, None, None));
    };
    let (fetch_from, fetch_to) = fetch_span.try_to_i64()?.as_tuple();

    reader
        .fetch((tile.tid as i32, fetch_from, fetch_to))
        .context(format!("fetch {} {}-{}", &tile.chr, fetch_from, fetch_to))?;

    let window_size = opt.window_size;
    let window_left = window_size / 2;
    let window_right = window_size - window_left;
    let context_left = window_left.saturating_add(extra_halo_bp);
    let context_right = window_right.saturating_add(extra_halo_bp);
    let window_left_i64 = window_left as i64;
    let window_right_i64 = window_right as i64;

    // Dilate the tile core so edge positions have full contexts
    // Outputs are trimmed back to the core
    // The difference buffers live on a dilated span that guarantees each core position
    // sees a complete window. We later trim the outputs back to the original core
    let dilated_start_abs = tile.core_start().saturating_sub(context_left);
    let dilated_end_abs = ((tile.core_end() as u64) + context_right as u64).min(chrom_len) as u32;
    if dilated_start_abs >= dilated_end_abs {
        return Ok((counter, None, None, None));
    }

    let dilated_span_len = (dilated_end_abs - dilated_start_abs) as usize; // Length of the dilated buffer (exclusive end)
    if dilated_span_len == 0 {
        return Ok((counter, None, None, None));
    }

    // Offsets of the original core within the dilated span. These values are measured relative
    // to `dilated_start_abs`, so they represent indices into the dilated buffers rather than
    // absolute genomic coordinates.
    let core_start_offset = (tile.core_start() - dilated_start_abs) as usize;
    let core_end_offset_exclusive = (tile.core_end() - dilated_start_abs) as usize;
    let dilated_start_i64 = dilated_start_abs as i64;
    let dilated_end_i64 = dilated_end_abs as i64;
    let dilated_start_abs_u64 = dilated_start_abs as u64; // Absolute coordinate of dilated buffer origin
    let core_start_abs = tile.core_start() as u64;
    let core_end_abs = tile.core_end() as u64;

    // Difference buffers sized to the dilated span plus sentinel
    // We keep them as f32 because weights are f32 and accumulation is stable over these spans
    let mut span_diff = vec![0f32; dilated_span_len + 1];
    let mut end_diff = vec![0f32; dilated_span_len + 1];

    let min_len = opt.min_fragment_length;
    let max_len = opt.max_fragment_length;

    let fragment_filter = move |frag: &Fragment| {
        let len = frag.len();
        len >= min_len && len <= max_len
    };

    let gc_tag_bytes = opt.gc.gc_tag.as_deref().map(|t| t.as_bytes().to_vec());
    let mut iter = if opt.unpaired.reads_are_fragments {
        let min_mapq = opt.min_mapq;
        let include_read_fn = move |r: &Record| default_include_read_unpaired(r, min_mapq);
        fragments_from_bam(
            reader.records().map(|r| r.map_err(anyhow::Error::from)),
            include_read_fn,
            gc_tag_bytes.as_deref(),
            fragment_filter,
            true,
        )
        .with_local_counters()
    } else {
        let min_mapq = opt.min_mapq;
        let require_proper_pair = opt.require_proper_pair;
        let include_read_fn =
            move |r: &Record| default_include_read_paired_end(r, require_proper_pair, min_mapq);
        fragments_from_bam(
            reader.records().map(|r| r.map_err(anyhow::Error::from)),
            include_read_fn,
            gc_tag_bytes.as_deref(),
            fragment_filter,
            false,
        )
        .with_local_counters()
    };

    let get_gc_weight = {
        let gc_corrector = gc_corrector_opt.as_ref();
        let gc_prefixes = gc_prefixes_opt.as_ref();
        move |fragment: &Fragment, fetch_start: u32| -> Result<Option<f64>> {
            match (gc_corrector, gc_prefixes) {
                (Some(corrector), Some(prefixes)) => {
                    let fetch_relative_fragment = fragment
                        .interval
                        .try_to_u64()?
                        .shift_left(fetch_start as u64)?;
                    corrector.correct_fragment(fetch_relative_fragment, prefixes)
                }
                _ => Ok(None),
            }
        }
    };

    let correct_gc = opt.gc.gc_file.is_some();
    let fetch_start = tile.fetch_start();
    let fetch_end = tile.fetch_end();

    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;

        if fragment.start() < fetch_start || fragment.end() > fetch_end {
            // Fragment won't overlap the counting region (assuming correct max_fragment_length+window halo!)
            continue;
        }

        // Get GC correction weight
        let gc_weight = if opt.gc.gc_tag.is_some() {
            match fragment.gc_tag.classify()? {
                ClassifiedGCTagWeight::Usable(weight) => weight as f64,
                ClassifiedGCTagWeight::Missing => {
                    counter.gc_failed_fragments += 1;
                    counter.gc_missing_tags += 1;
                    if opt.gc.neutralize_invalid_gc {
                        1.0
                    } else {
                        continue;
                    }
                }
                ClassifiedGCTagWeight::Invalid { out_of_range } => {
                    counter.gc_failed_fragments += 1;
                    if out_of_range {
                        counter.gc_out_of_range_tags += 1;
                    }
                    if opt.gc.neutralize_invalid_gc {
                        1.0
                    } else {
                        continue;
                    }
                }
            }
        } else {
            let gc_weight_opt = get_gc_weight(&fragment, fetch_start)?;
            match (gc_weight_opt, correct_gc) {
                (Some(w), true) => w,
                (None, true) => {
                    // Tried but failed to make a GC correction weight
                    counter.gc_failed_fragments += 1;
                    if opt.gc.neutralize_invalid_gc {
                        1.0
                    } else {
                        continue;
                    }
                }
                (None, false) => 1.0, // No correction
                (Some(_), false) => unreachable!(),
            }
        };

        let fragment_start = fragment.start() as i64;
        let fragment_end = fragment.end() as i64;

        // Full-span contribution: window stays entirely inside the fragment (edges may align)
        let was_pushed_1 = push_range(
            &mut span_diff,
            dilated_start_i64,
            dilated_end_i64,
            fragment_start + window_left_i64,
            fragment_end - window_right_i64 + 1,
            gc_weight as f32,
        );
        // Left endpoint contribution: window must still contain the fragment start
        let was_pushed_2 = push_range(
            &mut end_diff,
            dilated_start_i64,
            dilated_end_i64,
            fragment_start - window_right_i64 + 1,
            fragment_start + window_left_i64,
            gc_weight as f32,
        );
        // Right endpoint contribution: window must still contain the fragment end (exclusive)
        let was_pushed_3 = push_range(
            &mut end_diff,
            dilated_start_i64,
            dilated_end_i64,
            fragment_end - window_right_i64 + 1,
            fragment_end + window_left_i64 - 1,
            gc_weight as f32,
        );

        if was_pushed_1 || was_pushed_2 || was_pushed_3 {
            counter.base.counted_fragments += 1;
        }
    }

    let overlap_vals = finalize_diff(&mut span_diff);
    let end_vals = finalize_diff(&mut end_diff);
    let mut wps_values = overlap_vals;
    for (value, end_value) in wps_values.iter_mut().zip(end_vals.into_iter()) {
        *value -= end_value;
    }

    if !scaling_chr.is_empty() {
        // Scaling operates on the combined WPS signal (it may apply GC/other corrections)
        apply_scaling_to_coverage_in_place(&mut wps_values, dilated_start_abs, scaling_chr);
    }

    // Build a positional mask over the dilated span so blacklist and chromosome-edge positions
    // never contribute to protected/end counts
    let dilated_interval = Interval::new(dilated_start_abs as u64, dilated_end_abs as u64)?;
    let mask = build_mask_for_core(
        dilated_interval,
        blacklist_chr,
        chrom_len,
        window_left,
        window_right,
    );
    let mask_slice = mask.as_deref();

    counter.add_from_snapshot(iter.counters_snapshot());

    if return_wps_instead {
        return Ok((counter, None, Some(wps_values), mask));
    }

    let mut prefix_all_cache: Option<Vec<f64>> = None;
    let mut prefix_allowed_cache: Option<Vec<f64>> = None;
    let mut prefix_allowed_cnt_cache: Option<Vec<u32>> = None;

    let temp_output = match mode {
        TileMode::Positional {
            windows,
            out_path,
            indexed,
        } => {
            let mut positional_writer = open_zstd_auto_writer(&out_path, 3, None)?;

            match windows {
                None => {
                    write_bedgraph_runs(
                        &tile.chr,
                        &wps_values,
                        mask_slice,
                        core_start_offset,
                        core_end_offset_exclusive,
                        dilated_start_abs_u64,
                        decimals,
                        keep_zero_runs,
                        &mut positional_writer,
                    )?;
                }
                Some(windows_for_chr) => {
                    for window in
                        overlapping_windows_for_tile(windows_for_chr, tile, tile_window_span)
                    {
                        let window_start_abs = window.start();
                        let window_end_abs = window.end();
                        let original_idx = window.idx();
                        let local_overlap = match clip_window_to_core_and_localize(
                            window_start_abs,
                            window_end_abs,
                            tile.core_start(),
                            tile.core_end(),
                            dilated_start_abs,
                        )? {
                            Some(v) => v,
                            None => continue,
                        };

                        if indexed {
                            write_windowed_runs(
                                &tile.chr,
                                &wps_values,
                                mask_slice,
                                local_overlap.local_start_idx,
                                local_overlap.local_end_idx,
                                dilated_start_abs_u64,
                                Some(original_idx),
                                decimals,
                                keep_zero_runs,
                                &mut positional_writer,
                            )?;
                        } else {
                            write_windowed_runs(
                                &tile.chr,
                                &wps_values,
                                mask_slice,
                                local_overlap.local_start_idx,
                                local_overlap.local_end_idx,
                                dilated_start_abs_u64,
                                None,
                                decimals,
                                keep_zero_runs,
                                &mut positional_writer,
                            )?;
                        }
                    }
                }
            }

            positional_writer.flush()?;
            WpsTileTempOutput::Positional {
                chromosome: tile.chr.clone(),
                tile_index: tile.index,
                path: out_path,
            }
        }

        TileMode::AggregatesByBed {
            windows,
            masked: _,
            partials_out,
            cross_idx_out,
        } => {
            let masked_mode = mask_slice.is_some();
            let ps_all_slice = {
                if prefix_all_cache.is_none() {
                    prefix_all_cache = Some(build_prefix(&wps_values));
                }
                prefix_all_cache
                    .as_ref()
                    .context("WPS prefix cache missing after initialization")?
                    .as_slice()
            };
            let (ps_allowed_slice, ps_allowed_cnt_slice) = if masked_mode {
                let mask_slice_ref = mask_slice
                    .context("masked WPS aggregate reduction requires a blacklist mask slice")?;
                if prefix_allowed_cache.is_none() {
                    let (allowed_prefix, allowed_count_prefix) =
                        build_allowed_prefix(&wps_values, mask_slice_ref);
                    prefix_allowed_cache = Some(allowed_prefix);
                    prefix_allowed_cnt_cache = Some(allowed_count_prefix);
                }
                (
                    prefix_allowed_cache.as_ref().map(|v| v.as_slice()),
                    prefix_allowed_cnt_cache.as_ref().map(|v| v.as_slice()),
                )
            } else {
                (None, None)
            };

            let mut partials_writer = open_zstd_auto_writer(&partials_out, 3, None)?;
            let mut cross_index_writer = open_zstd_auto_writer(&cross_idx_out, 3, None)?;

            for window in overlapping_windows_for_tile(windows, tile, tile_window_span) {
                let window_start_abs = window.start();
                let window_end_abs = window.end();
                let original_idx = window.idx();
                let local_overlap = clip_window_to_core_and_localize(
                    window_start_abs,
                    window_end_abs,
                    tile.core_start(),
                    tile.core_end(),
                    dilated_start_abs,
                )?;
                let Some(local_overlap) = local_overlap else {
                    continue;
                };

                let crosses_boundary =
                    !(window_start_abs >= core_start_abs && window_end_abs <= core_end_abs);

                let (sum, allowed, blacklisted) = wps_sum_and_counts(
                    local_overlap.local_start_idx,
                    local_overlap.local_end_idx,
                    masked_mode,
                    ps_all_slice,
                    ps_allowed_slice,
                    ps_allowed_cnt_slice,
                    mask_slice,
                );

                writeln!(
                    partials_writer,
                    "{}\t{}\t{}\t{}",
                    original_idx, sum, allowed, blacklisted
                )?;
                if crosses_boundary {
                    writeln!(cross_index_writer, "{}", original_idx)?;
                }
            }

            partials_writer.flush()?;
            cross_index_writer.flush()?;
            WpsTileTempOutput::AggregatesByBed {
                chromosome: tile.chr.clone(),
                tile_index: tile.index,
                partials_path: partials_out,
                cross_index_path: Some(cross_idx_out),
            }
        }

        TileMode::AggregatesBySize {
            window_bp,
            masked: _,
            finals_out,
            partials_out,
            cross_idx_out,
            guaranteed_aligned,
        } => {
            let action = per_window_wps_action
                .context("aggregate WPS tile mode requires a per-window action")?;
            let masked_mode = mask_slice.is_some();
            let ps_all_slice = {
                if prefix_all_cache.is_none() {
                    prefix_all_cache = Some(build_prefix(&wps_values));
                }
                prefix_all_cache
                    .as_ref()
                    .context("WPS prefix cache missing after initialization")?
                    .as_slice()
            };
            let (ps_allowed_slice, ps_allowed_cnt_slice) = if masked_mode {
                let mask_slice_ref = mask_slice
                    .context("masked WPS aggregate reduction requires a blacklist mask slice")?;
                if prefix_allowed_cache.is_none() {
                    let (pa, cnt) = build_allowed_prefix(&wps_values, mask_slice_ref);
                    prefix_allowed_cache = Some(pa);
                    prefix_allowed_cnt_cache = Some(cnt);
                }
                (
                    prefix_allowed_cache.as_ref().map(|v| v.as_slice()),
                    prefix_allowed_cnt_cache.as_ref().map(|v| v.as_slice()),
                )
            } else {
                (None, None)
            };

            let tile_core_start_abs = core_start_abs;
            let tile_core_end_abs = core_end_abs;
            let first_bin_index = tile_core_start_abs / window_bp;
            let last_bin_index = (tile_core_end_abs.saturating_sub(1)) / window_bp;

            if guaranteed_aligned {
                let mut finals_writer = open_zstd_auto_writer(&finals_out, 3, None)?;

                for bin_index in first_bin_index..=last_bin_index {
                    let bin_start = bin_index * window_bp;
                    let bin_end = (bin_index + 1) * window_bp;

                    let local_overlap = match clip_window_to_core_and_localize(
                        bin_start,
                        bin_end,
                        tile.core_start(),
                        tile.core_end(),
                        dilated_start_abs,
                    )? {
                        Some(v) => v,
                        None => continue,
                    };

                    let (sum, allowed, blacklisted) = wps_sum_and_counts(
                        local_overlap.local_start_idx,
                        local_overlap.local_end_idx,
                        masked_mode,
                        ps_all_slice,
                        ps_allowed_slice,
                        ps_allowed_cnt_slice,
                        mask_slice,
                    );

                    let unmasked_span_bp = local_overlap.clipped_abs_interval.len();
                    let value =
                        finalize_value(sum, allowed, unmasked_span_bp, masked_mode, &action);
                    let value = round_to(value, decimals);

                    write_final_row(
                        &mut finals_writer,
                        &tile.chr,
                        local_overlap.clipped_abs_interval,
                        value,
                        blacklisted,
                        decimals,
                    )?;
                }

                finals_writer.flush()?;
                WpsTileTempOutput::SizeFinal {
                    chromosome: tile.chr.clone(),
                    tile_index: tile.index,
                    path: finals_out,
                }
            } else {
                let mut partials_writer = open_zstd_auto_writer(&partials_out, 3, None)?;
                let mut cross_index_writer = open_zstd_auto_writer(&cross_idx_out, 3, None)?;

                for bin_index in first_bin_index..=last_bin_index {
                    let bin_start = bin_index * window_bp;
                    let bin_end = (bin_index + 1) * window_bp;

                    let local_overlap = match clip_window_to_core_and_localize(
                        bin_start,
                        bin_end,
                        tile.core_start(),
                        tile.core_end(),
                        dilated_start_abs,
                    )? {
                        Some(v) => v,
                        None => continue,
                    };

                    let (sum, allowed, blacklisted) = wps_sum_and_counts(
                        local_overlap.local_start_idx,
                        local_overlap.local_end_idx,
                        masked_mode,
                        ps_all_slice,
                        ps_allowed_slice,
                        ps_allowed_cnt_slice,
                        mask_slice,
                    );

                    writeln!(
                        partials_writer,
                        "{}\t{}\t{}\t{}\t{}",
                        bin_start, bin_end, sum, allowed, blacklisted
                    )?;

                    let fully_inside =
                        (bin_start >= tile_core_start_abs) && (bin_end <= tile_core_end_abs);
                    if !fully_inside {
                        writeln!(cross_index_writer, "{}", bin_start)?;
                    }
                }

                partials_writer.flush()?;
                cross_index_writer.flush()?;
                WpsTileTempOutput::AggregatesBySize {
                    chromosome: tile.chr.clone(),
                    tile_index: tile.index,
                    partials_path: partials_out,
                    cross_index_path: Some(cross_idx_out),
                }
            }
        }
    };

    Ok((counter, Some(temp_output), None, None))
}

// TODO: Add title to docstring
/// Returns whether the range was pushed.
fn push_range(
    diff: &mut [f32],
    core_start: i64,
    core_end: i64,
    raw_start: i64,
    raw_end: i64,
    weight: f32,
) -> bool {
    let start = raw_start.max(core_start);
    let end = raw_end.min(core_end);
    if start >= end {
        return false;
    }
    let from = (start - core_start) as usize;
    let to = (end - core_start) as usize;
    diff[from] += weight;
    diff[to] -= weight;
    true
}

fn finalize_diff(diff: &mut [f32]) -> Vec<f32> {
    let mut acc = 0.0f32;
    let len = diff.len().saturating_sub(1);
    let mut out = Vec::with_capacity(len);
    for i in 0..len {
        acc += diff[i];
        out.push(acc);
    }
    out
}

/// Build a mask over the dilated tile span marking blacklisted bases and centers
/// whose WPS window would exceed chromosome bounds.
fn build_mask_for_core(
    dilated_interval: Interval<u64>,
    blacklist_intervals: &[Interval<u64>],
    chromosome_length: u64,
    left_span: u32,
    right_span: u32,
) -> Option<Vec<u8>> {
    let dilated_span_len = dilated_interval.len() as usize;
    let dilated_start_abs = dilated_interval.start();
    let mut mask = vec![0u8; dilated_span_len];
    let mut has_masked_positions = false;

    for interval in blacklist_intervals {
        let Some(clipped_interval) = interval.clip_to(dilated_interval) else {
            continue;
        };

        // Treat blacklist intervals as half-open [start, end) so the end offset stays exclusive
        let start_offset = (clipped_interval.start() - dilated_start_abs) as usize;
        let end_offset_exclusive = (clipped_interval.end() - dilated_start_abs) as usize;
        if end_offset_exclusive > start_offset {
            mask[start_offset..end_offset_exclusive].fill(1);
            has_masked_positions = true;
        }
    }

    let left_span_u64 = left_span as u64;
    let right_span_u64 = right_span as u64;

    for offset in 0..dilated_span_len {
        let center_abs = dilated_start_abs + offset as u64;
        if center_abs < left_span_u64 || center_abs + right_span_u64 > chromosome_length {
            if mask[offset] == 0 {
                mask[offset] = 1;
            }
            has_masked_positions = true;
        }
    }

    if has_masked_positions {
        Some(mask)
    } else {
        None
    }
}

fn build_prefix(values: &[f32]) -> Vec<f64> {
    let mut prefix_sums = Vec::with_capacity(values.len() + 1);
    let mut running_total = 0.0f64;
    prefix_sums.push(0.0);
    for &value in values {
        running_total += value as f64;
        prefix_sums.push(running_total);
    }
    prefix_sums
}

fn build_allowed_prefix(values: &[f32], mask: &[u8]) -> (Vec<f64>, Vec<u32>) {
    let mut prefix_sum_allowed = Vec::with_capacity(values.len() + 1);
    let mut prefix_count_allowed = Vec::with_capacity(values.len() + 1);
    let mut running_sum_allowed = 0.0f64;
    let mut running_count_allowed = 0u32;
    prefix_sum_allowed.push(0.0);
    prefix_count_allowed.push(0);
    for (position_index, &value) in values.iter().enumerate() {
        if mask.get(position_index).copied().unwrap_or(0) == 0 {
            running_sum_allowed += value as f64;
            running_count_allowed += 1;
        }
        prefix_sum_allowed.push(running_sum_allowed);
        prefix_count_allowed.push(running_count_allowed);
    }
    (prefix_sum_allowed, prefix_count_allowed)
}

fn wps_sum_and_counts(
    local_start_idx: usize,
    local_end_idx: usize,
    masked: bool,
    prefix_all: &[f64],
    prefix_allowed: Option<&[f64]>,
    prefix_allowed_count: Option<&[u32]>,
    mask: Option<&[u8]>,
) -> (f64, u64, u64) {
    let window_span_len = (local_end_idx - local_start_idx) as u64;
    let total_sum = prefix_all[local_end_idx] - prefix_all[local_start_idx];

    if !masked {
        return (total_sum, window_span_len, 0);
    }

    let allowed_sum = if let Some(pa) = prefix_allowed {
        pa[local_end_idx] - pa[local_start_idx]
    } else {
        total_sum
    };

    let allowed_count = if let Some(cnt) = prefix_allowed_count {
        (cnt[local_end_idx] - cnt[local_start_idx]) as u64
    } else if let Some(m) = mask {
        let mut allowed_positions = 0u64;
        for position_index in local_start_idx..local_end_idx {
            if m[position_index] == 0 {
                allowed_positions += 1;
            }
        }
        allowed_positions
    } else {
        window_span_len
    };

    let blacklisted_positions = window_span_len.saturating_sub(allowed_count);
    (allowed_sum, allowed_count, blacklisted_positions)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CoreClip {
    // Inclusive start index in dilated tile-local coordinates
    local_start_idx: usize,
    // Exclusive end index in dilated tile-local coordinates
    local_end_idx: usize,
    // Absolute interval after clipping the requested window to the tile core
    clipped_abs_interval: Interval<u64>,
}

#[inline]
fn clip_window_to_core_and_localize(
    abs_start: u64,
    abs_end: u64,
    core_start: u32,
    core_end: u32,
    dilated_start: u32,
) -> Result<Option<CoreClip>> {
    let core_start_u64 = core_start as u64;
    let core_end_u64 = core_end as u64;
    let start = abs_start.max(core_start_u64);
    let end = abs_end.min(core_end_u64);
    if start >= end {
        return Ok(None);
    }
    let local_start = (start as u32 - dilated_start) as usize;
    let local_end = (end as u32 - dilated_start) as usize;
    Ok(Some(CoreClip {
        local_start_idx: local_start,
        local_end_idx: local_end,
        clipped_abs_interval: Interval::new(start, end)?,
    }))
}

#[cfg(test)]
mod tests {
    include!("wps_tests.rs");
}
