use crate::commands::fcoverage::config::FCoverageConfig;
use crate::commands::fcoverage::tiling::{
    adapt_fetch_to_extreme_windows, build_summary_prefixes, clip_interval_to_core_and_localize,
    coverage_sum_and_counts, coverage_summary_and_counts, finalize_value,
    merge_positional_tiles_with_optional_scaling,
};
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use crate::commands::fcoverage::writers::{
    derive_summary_stats, write_bed_aggregate_output, write_bedgraph_runs, write_final_row,
    write_grouped_bed_aggregate_output, write_size_aggregate_output, write_summary_stats_row,
    write_windowed_runs,
};
use crate::commands::gc_bias::correct::{GCCorrector, load_gc_corrector};
use crate::commands::gc_bias::counting::build_gc_prefixes;
use crate::shared::coverage::{Coverage, clamp_finite_coverage_below_to_zero};
use crate::shared::formatters::round_to;
use crate::shared::fragment::minimal_fragment::Fragment;
use crate::shared::fragment::segment_fragment::FragmentWithSegments;
use crate::shared::fragment_iterators::fragments_with_segments_from_bam;
use crate::shared::gc_tag::{ClassifiedGCTagWeight, MIN_REASONABLE_GC_WEIGHT};
use crate::shared::interval::{IndexedInterval, Interval};
use crate::shared::io::dot_join;
use crate::shared::progress::ProgressFactory;
use crate::shared::read::{default_include_read_paired_end, default_include_read_unpaired};
use crate::shared::reference::read_seq_in_range;
use crate::shared::scale_genome::apply_scaling_to_coverage_in_place;
use crate::shared::tiled_run::{
    Tile, TileMode, TileWindowSpan, build_tiles, make_temp_dir, overlapping_windows_for_tile,
    precompute_tile_window_spans,
};
use crate::{
    commands::cli_common::{
        DistributionWindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
        resolve_chromosomes_and_contigs,
    },
    commands::counters::FCoverageCounters,
    commands::run_statistics::{
        DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions, GCStatisticsSummary,
        TILE_DOUBLE_COUNT_NOTE, print_fragment_run_statistics,
    },
    shared::{
        bam::create_chromosome_reader,
        bed::{
            GroupedCoverageLayout, build_grouped_coverage_layout, load_grouped_windows_from_bed,
            load_windows_from_bed, write_group_idx_to_name_tsv,
        },
        thread_pool::init_global_pool,
        writers::open_zstd_auto_writer,
    },
};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::io::Write;
use std::path::PathBuf;
use std::{sync::Arc, time::Instant};
use tracing::{info, warn};

const COMMAND_TARGET: &str = "fcoverage";

/// Result of an internal `fcoverage` run.
///
/// This is used by other commands that reuse the tiled counting and final by-size
/// reduction without wanting `fcoverage`'s outer statistics wrapper.
pub struct FCoverageRunResult {
    pub counters: FCoverageCounters,
    pub mean_normalization_length: Option<f64>,
    pub final_out_path: PathBuf,
}

/// Execute the fragment coverage pipeline end-to-end.
///
/// Implementation details:
/// - Resolves chromosomes, prepares IO state, then iterates tiles in parallel using Rayon.
/// - Collects per-tile coverage into temporary artefacts before merging them into the final
///   positional or aggregated outputs.
/// - Applies fragment length, blacklist, and optional scaling filters during iteration.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `fcoverage` command.
///
/// Returns:
/// - `Ok(())` when positional and/or windowed outputs are written successfully.
///
/// Errors:
/// - Returns an error if the BAM cannot be read, auxiliary files are invalid, or writing outputs
///   fails at any stage.
pub fn run(opt: &FCoverageConfig) -> Result<()> {
    let start_time = Instant::now();

    let run_result = run_inner(opt)?;
    let global_counter = run_result.counters;

    let elapsed = start_time.elapsed();
    let mut extra_statistics = Vec::new();
    if let Some(mean_normalization_length) = run_result.mean_normalization_length {
        extra_statistics.push(format!(
            "Mean normalization length: {}",
            mean_normalization_length
        ));
    }

    print_fragment_run_statistics(
        &global_counter.base,
        elapsed,
        FragmentRunStatisticsOptions {
            include_section_header: true,
            notes: &[TILE_DOUBLE_COUNT_NOTE],
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
        extra_statistics.iter().map(String::as_str),
    );
    Ok(())
}

pub fn run_inner(opt: &FCoverageConfig) -> Result<FCoverageRunResult> {
    opt.fragment_lengths.validate()?;
    opt.gc.validate(opt.ref_2bit.as_deref())?;
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    if opt.unpaired.reads_are_fragments && opt.ignore_gap {
        bail!("--ignore-gap cannot be used with --reads-are-fragments");
    }
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    let window_opt = opt.windows.resolve_windows();
    let prefix = opt.output_prefix.trim();

    // Create output directory
    ensure_output_dir(&opt.ioc.output_dir)?;

    if opt.blacklist.is_some() {
        info!(target: COMMAND_TARGET, "Loading blacklists");
    }
    let blacklist_map = load_blacklist_map(opt.blacklist.as_ref(), 1, 0, &chromosomes)?;

    let grouped_windowed = matches!(window_opt, DistributionWindowSpec::GroupedBed(_));
    let unique_base_grouped = opt.per_window.is_unique_base_grouped_action();
    match (&window_opt, opt.per_window) {
        (
            DistributionWindowSpec::Size(_),
            CoverageWindowAction::OnlyIncludeThesePositionsUnique
            | CoverageWindowAction::OnlyIncludeThesePositionsIndexed
            | CoverageWindowAction::AverageOnUniqueBases
            | CoverageWindowAction::TotalOnUniqueBases
            | CoverageWindowAction::SummaryStatsOnUniqueBases,
        ) => {
            bail!(
                "in --by-size mode, --per-window can only be 'average', 'total', or 'summary-stats'"
            );
        }
        (
            DistributionWindowSpec::Bed(_),
            CoverageWindowAction::AverageOnUniqueBases
            | CoverageWindowAction::TotalOnUniqueBases
            | CoverageWindowAction::SummaryStatsOnUniqueBases,
        ) => {
            bail!(
                "'*-on-unique-bases' actions require --by-grouped-bed because they change grouped row semantics"
            );
        }
        (
            DistributionWindowSpec::GroupedBed(_),
            CoverageWindowAction::OnlyIncludeThesePositionsUnique
            | CoverageWindowAction::OnlyIncludeThesePositionsIndexed,
        ) => {
            bail!("grouped BED supports only aggregate outputs in fcoverage");
        }
        (DistributionWindowSpec::Global, CoverageWindowAction::Total) => {
            bail!(
                "without windowing, --per-window total is not supported because fcoverage writes positional coverage; use --by-size, --by-bed, or --by-grouped-bed"
            );
        }
        (
            DistributionWindowSpec::Global,
            CoverageWindowAction::SummaryStats
            | CoverageWindowAction::AverageOnUniqueBases
            | CoverageWindowAction::TotalOnUniqueBases
            | CoverageWindowAction::SummaryStatsOnUniqueBases,
        ) => {
            bail!("the requested --per-window mode requires windowed or grouped inputs");
        }
        _ => {}
    }
    // Load the selected windows
    let mut windows_map: Option<FxHashMap<String, crate::shared::bed::Windows>> = None;
    let mut grouped_layout: Option<GroupedCoverageLayout> = None;

    match &window_opt {
        DistributionWindowSpec::Bed(bed) => {
            info!(target: COMMAND_TARGET, "Loading window coordinates");
            let wds = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            if matches!(
                opt.per_window,
                CoverageWindowAction::OnlyIncludeThesePositionsUnique
            ) {
                // Merge in-place to avoid double memory usage
                info!(target: COMMAND_TARGET, "Merging overlapping/touching windows");
                // Take ownership so we can remove entries by chromosome
                let mut wds_owned = wds;
                let mut flattened_windows =
                    FxHashMap::with_capacity_and_hasher(wds_owned.len(), Default::default());

                let mut next_idx: u64 = 0;
                // Use the user-provided chromosome order to assign indices deterministically
                for chromosome in &chromosomes {
                    if let Some(windows_for_chr) = wds_owned.remove(chromosome) {
                        // Flatten in-place
                        let (flattened, next) = windows_for_chr.into_flattened_reindexed(next_idx);
                        next_idx = next;
                        flattened_windows.insert(chromosome.clone(), flattened);
                    }
                }
                windows_map = Some(flattened_windows);
            } else {
                windows_map = Some(wds);
            }
        }
        DistributionWindowSpec::GroupedBed(bed) => {
            info!(target: COMMAND_TARGET, "Loading grouped window coordinates");
            let (grouped_windows_by_chr, group_idx_to_name) =
                load_grouped_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            grouped_layout = Some(build_grouped_coverage_layout(
                &grouped_windows_by_chr,
                &group_idx_to_name,
                &chromosomes,
                unique_base_grouped,
            )?);
            windows_map = Some(
                grouped_layout
                    .as_ref()
                    .expect("grouped coverage layout must exist")
                    .segments_by_chr
                    .clone(),
            );
        }
        DistributionWindowSpec::Size(_) | DistributionWindowSpec::Global => {}
    }

    // Load genomic scaling factors
    if opt.scale_genome.scaling_factors.is_some() {
        info!(target: COMMAND_TARGET, "Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> = load_scaling_map(
        &opt.scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(
            opt.gc.gc_file.is_some(),
            opt.gc.gc_tag.is_some(),
        ),
        Some(opt.ignore_gap),
    )?;

    // Load GC correction package if specified
    if opt.gc.gc_file.is_some() {
        info!(target: COMMAND_TARGET, "Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.gc.gc_file.as_ref(),
        opt.fragment_lengths.min_fragment_length,
        opt.fragment_lengths.max_fragment_length,
    )?;

    // Decide mode once
    let windowed = matches!(
        window_opt,
        DistributionWindowSpec::Bed(_)
            | DistributionWindowSpec::Size(_)
            | DistributionWindowSpec::GroupedBed(_)
    );
    let masked = opt.blacklist.is_some();
    let has_scaling_or_correction = opt.scale_genome.scaling_factors.is_some()
        || opt.gc.gc_file.is_some()
        || opt.gc.gc_tag.is_some()
        || opt.uses_length_normalization();

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;

    // Window size when --by-size (otherwise None)
    let by_size_bp: Option<u64> = match &window_opt {
        DistributionWindowSpec::Size(bp) => Some(*bp),
        _ => None,
    };

    // Build tiles
    let halo_bp: u32 = opt.fragment_lengths.max_fragment_length; // Safe halo for pairing/segments
    let (tiles, tile_and_window_boundaries_align) =
        build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, by_size_bp)?;

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

    let length_norm_prefix = match opt.normalize_by_length_mode {
        crate::commands::fcoverage::config::LengthNormalizationMode::Off => "",
        crate::commands::fcoverage::config::LengthNormalizationMode::UnitMass => {
            "length_normalized"
        }
        crate::commands::fcoverage::config::LengthNormalizationMode::RestoreMean => {
            "length_normalized.restored_mean"
        }
    };

    // Create filenames of final outputs
    let final_bedgraph_pos_name = dot_join(&[
        prefix,
        length_norm_prefix,
        "fcoverage.per_position.bedgraph.zst",
    ]);
    let final_tsv_pos_name = dot_join(&[
        prefix,
        length_norm_prefix,
        "fcoverage.per_position_per_window.tsv.zst",
    ]);
    let final_aggregate_name = dot_join(&[
        prefix,
        length_norm_prefix,
        &format!("fcoverage.{}.tsv.zst", opt.per_window.action_file_stem()),
    ]);

    // Get decimals to use
    let decimals_to_use: i32 = if windowed {
        match opt.per_window {
            CoverageWindowAction::Average
            | CoverageWindowAction::Total
            | CoverageWindowAction::SummaryStats
            | CoverageWindowAction::AverageOnUniqueBases
            | CoverageWindowAction::TotalOnUniqueBases
            | CoverageWindowAction::SummaryStatsOnUniqueBases => opt.decimals as i32,
            CoverageWindowAction::OnlyIncludeThesePositionsUnique
            | CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                if has_scaling_or_correction {
                    opt.decimals as i32
                } else {
                    0
                }
            }
        }
    } else {
        if has_scaling_or_correction {
            opt.decimals as i32
        } else {
            0
        }
    };

    let total_tiles = tiles.len();

    // Create progress bar
    let progress = ProgressFactory::new();
    let pb = Arc::new(progress.default_bar(total_tiles as u64));

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    let mut global_counter = FCoverageCounters::default();

    info!(target: COMMAND_TARGET, "Counting per tile");

    pb.set_position(0);

    let tile_window_spans_for_threads = tile_window_spans.clone();
    let gc_tag = opt.gc.gc_tag.as_deref();

    let tile_results: Vec<FCoverageCounters> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<FCoverageCounters> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let windows_chr: Option<&[IndexedInterval<u64>]> = windows_map
                .as_ref()
                .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice()));
            let blacklist_chr: &[Interval<u64>] = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr: &[(u64, u64, f32)] = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            let counter = if matches!(
                window_opt,
                DistributionWindowSpec::Bed(_) | DistributionWindowSpec::GroupedBed(_)
            ) && windows_chr.map_or(true, |windows| windows.is_empty())
            {
                // Skip this no-window tile and just return an empty counter
                FCoverageCounters::default()
            } else {
                // Decide tile mode and filename
                let (action_prefix, extensions) = if windowed {
                    match opt.per_window {
                        CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                            (positional_prefix, "tsv.zst")
                        }
                        CoverageWindowAction::OnlyIncludeThesePositionsUnique => {
                            (positional_prefix, "bedgraph.zst")
                        }
                        CoverageWindowAction::Average
                        | CoverageWindowAction::Total
                        | CoverageWindowAction::SummaryStats
                        | CoverageWindowAction::AverageOnUniqueBases
                        | CoverageWindowAction::TotalOnUniqueBases
                        | CoverageWindowAction::SummaryStatsOnUniqueBases => {
                            (partials_prefix, "tsv.zst")
                        }
                    }
                } else {
                    // Whole positional coverage
                    (positional_prefix, "bedgraph.zst")
                };

                let out_path = temp_dir.join(format!(
                    "{prefix}.{chr}.{idx}.{extensions}",
                    prefix = action_prefix,
                    chr = tile.chr,
                    idx = tile.index
                ));

                // Windowed tmp outputs for faster reducer later on
                let partials_out = temp_dir.join(format!(
                    "{prefix}.{chr}.{idx}.{extensions}",
                    prefix = partials_prefix,
                    chr = tile.chr,
                    idx = tile.index
                ));
                let cross_idx_out = temp_dir.join(format!(
                    "{prefix}.cross.{chr}.{idx}.cross.zst", // Needs this extension!
                    prefix = partials_prefix,
                    chr = tile.chr,
                    idx = tile.index
                ));
                let finals_out = temp_dir.join(format!(
                    "{prefix}.{chr}.{idx}.{extensions}",
                    prefix = finals_prefix,
                    chr = tile.chr,
                    idx = tile.index
                ));

                let mode = if !windowed {
                    TileMode::Positional {
                        windows: None,
                        out_path,
                        indexed: false,
                    }
                } else {
                    match (&window_opt, opt.per_window) {
                        (
                            DistributionWindowSpec::Bed(_),
                            CoverageWindowAction::OnlyIncludeThesePositionsUnique,
                        ) => TileMode::Positional {
                            windows: windows_chr,
                            out_path,
                            indexed: false,
                        },
                        (
                            DistributionWindowSpec::Bed(_),
                            CoverageWindowAction::OnlyIncludeThesePositionsIndexed,
                        ) => TileMode::Positional {
                            windows: windows_chr,
                            out_path,
                            indexed: true,
                        },
                        (
                            DistributionWindowSpec::Bed(_) | DistributionWindowSpec::GroupedBed(_),
                            CoverageWindowAction::Average
                            | CoverageWindowAction::Total
                            | CoverageWindowAction::SummaryStats
                            | CoverageWindowAction::AverageOnUniqueBases
                            | CoverageWindowAction::TotalOnUniqueBases
                            | CoverageWindowAction::SummaryStatsOnUniqueBases,
                        ) => TileMode::AggregatesByBed {
                            windows: windows_chr.context(
                                "BED aggregate tile reached processing without any windows",
                            )?,
                            masked,
                            partials_out,
                            cross_idx_out,
                        },
                        (
                            DistributionWindowSpec::Size(size),
                            CoverageWindowAction::Average
                            | CoverageWindowAction::Total
                            | CoverageWindowAction::SummaryStats,
                        ) => TileMode::AggregatesBySize {
                            window_bp: *size,
                            masked,
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

                let tile_output_decimals = match &mode {
                    TileMode::Positional { .. }
                        if opt.restores_mean_after_length_normalization() =>
                    {
                        // Positional tile files are re-read and scaled during the late
                        // restore-mean merge, so keep extra precision here only
                        decimals_to_use.max(12)
                    }
                    _ => decimals_to_use,
                };

                process_tile(
                    opt,
                    tile,
                    tile_span.as_ref(),
                    blacklist_chr,
                    scaling_chr,
                    gc_corrector.clone(), // Quite small memory footprint
                    gc_tag,
                    mode,
                    tile_output_decimals,
                )?
            };
            pb.inc(1);
            Ok(counter)
        })
        .collect::<anyhow::Result<_>>()?; // Short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    // Release per-tile inputs before merging outputs
    drop(tile_window_spans_for_threads);
    drop(tile_window_spans);
    drop(tiles);
    drop(blacklist_map);
    drop(scaling_map);
    drop(gc_corrector);

    // Collect counters
    for counter in tile_results {
        global_counter += counter;
    }

    if opt.uses_length_normalization() {
        debug_assert!(
            global_counter.tile_owned_normalization_fragments
                <= global_counter.base.counted_fragments,
            "tile-owned normalization fragments must not exceed counted_fragments"
        );
        debug_assert_eq!(
            global_counter.base.counted_fragments == 0,
            global_counter.tile_owned_normalization_fragments == 0,
            "counted_fragments and tile-owned normalization fragments should agree on zero/non-zero"
        );
    }

    let mean_normalization_length = if opt.uses_length_normalization()
        && global_counter.tile_owned_normalization_fragments > 0
    {
        Some(
            global_counter.tile_owned_normalization_length_sum as f64
                / global_counter.tile_owned_normalization_fragments as f64,
        )
    } else {
        None
    };
    let restore_mean_multiplier = if opt.restores_mean_after_length_normalization() {
        mean_normalization_length
    } else {
        None
    };
    if opt.restores_mean_after_length_normalization() && global_counter.base.counted_fragments == 0
    {
        warn!(
            target: COMMAND_TARGET,
            "restore-mean requested, but mean_normalization_length was undefined because no fragments were counted; leaving the output in unit-mass scale"
        );
    }

    info!(target: COMMAND_TARGET, "Merging temporary tile files to final output");

    // Merge temporary output files and
    // reduce windows present in multiple tiles

    let final_out_path = if !windowed {
        // Whole-genome positional coverage
        merge_positional_tiles_with_optional_scaling(
            &temp_dir,
            &opt.ioc.output_dir,
            &chromosomes,
            positional_prefix,
            final_bedgraph_pos_name.as_str(),
            restore_mean_multiplier,
            false,
            decimals_to_use,
            opt.ioc.n_threads,
        )?
    } else {
        match opt.per_window {
            CoverageWindowAction::OnlyIncludeThesePositionsUnique => {
                // Windowed positional (unique and non-indexed)
                merge_positional_tiles_with_optional_scaling(
                    &temp_dir,
                    &opt.ioc.output_dir,
                    &chromosomes,
                    positional_prefix,
                    final_bedgraph_pos_name.as_str(),
                    restore_mean_multiplier,
                    false,
                    decimals_to_use,
                    opt.ioc.n_threads,
                )?
            }
            CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                // Windowed positional with orig_idx column
                merge_positional_tiles_with_optional_scaling(
                    &temp_dir,
                    &opt.ioc.output_dir,
                    &chromosomes,
                    positional_prefix,
                    final_tsv_pos_name.as_str(),
                    restore_mean_multiplier,
                    true,
                    decimals_to_use,
                    opt.ioc.n_threads,
                )?
            }
            CoverageWindowAction::Average
            | CoverageWindowAction::Total
            | CoverageWindowAction::SummaryStats
            | CoverageWindowAction::AverageOnUniqueBases
            | CoverageWindowAction::TotalOnUniqueBases
            | CoverageWindowAction::SummaryStatsOnUniqueBases => {
                let final_path = opt.ioc.output_dir.join(final_aggregate_name.as_str());

                match &window_opt {
                    DistributionWindowSpec::Bed(_) => write_bed_aggregate_output(
                        &final_path,
                        &temp_dir,
                        partials_prefix,
                        windows_map
                            .as_ref()
                            .context("BED aggregate reduction requires loaded windows")?,
                        &chromosomes,
                        masked,
                        opt.per_window,
                        decimals_to_use,
                        opt.ioc.n_threads,
                        restore_mean_multiplier,
                    )?,
                    DistributionWindowSpec::GroupedBed(_) => write_grouped_bed_aggregate_output(
                        &final_path,
                        &temp_dir,
                        partials_prefix,
                        grouped_layout
                            .as_ref()
                            .context("grouped aggregate reduction requires grouped layout")?,
                        &chromosomes,
                        opt.per_window,
                        decimals_to_use,
                        opt.ioc.n_threads,
                        restore_mean_multiplier,
                    )?,
                    DistributionWindowSpec::Size(_) => write_size_aggregate_output(
                        &final_path,
                        &temp_dir,
                        partials_prefix,
                        finals_prefix,
                        &chromosomes,
                        &contigs,
                        masked,
                        opt.per_window,
                        decimals_to_use,
                        opt.ioc.n_threads,
                        tile_and_window_boundaries_align,
                        restore_mean_multiplier,
                    )?,
                    DistributionWindowSpec::Global => {
                        bail!("unexpected global aggregate path in windowed fcoverage output")
                    }
                }

                if grouped_windowed {
                    let group_index_path = opt
                        .ioc
                        .output_dir
                        .join(dot_join(&[prefix, "group_index.tsv"]));
                    write_group_idx_to_name_tsv(
                        group_index_path,
                        &grouped_layout
                            .as_ref()
                            .context("grouped outputs require group index metadata")?
                            .group_idx_to_name,
                    )?;
                }

                final_path
            }
        }
    };

    info!(
        target: COMMAND_TARGET,
        "Saved output to: {}",
        final_out_path.display()
    );

    if let Err(e) = std::fs::remove_dir_all(&temp_dir) {
        warn!(
            target: COMMAND_TARGET,
            "failed to remove temp dir {}: {}",
            temp_dir.display(),
            e
        );
    }

    Ok(FCoverageRunResult {
        counters: global_counter,
        mean_normalization_length,
        final_out_path,
    })
}

/// Process one tile: pair reads, build coverage, and write outputs for this tile
fn process_tile(
    opt: &FCoverageConfig,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_chr: &[Interval<u64>],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<GCCorrector>,
    gc_tag: Option<&str>,
    mode: TileMode,
    decimals: i32,
) -> Result<FCoverageCounters> {
    // Open a fresh BAM reader for this thread
    let (mut reader, _tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    debug_assert!(_tid_check == tile.tid as u32);

    // Counters
    let mut counter = FCoverageCounters::default();

    let gc_prefixes_opt = if gc_corrector_opt.is_some() {
        let ref_2bit = match opt.ref_2bit.as_ref() {
            Some(r) => r,
            None => bail!("When GC correction is specified, --ref-2bit must also be specified"),
        };
        let seq_bytes = read_seq_in_range(
            ref_2bit,
            &tile.chr,
            // NOTE: Need for full fetch span to get GC of overlapping fragments!
            (tile.fetch_start() as usize)..(tile.fetch_end() as usize),
        )?;
        Some(build_gc_prefixes(&seq_bytes))
    } else {
        None
    };

    // Adapt the fetch coordinates to the present windows (*in windowed mode!*)
    // When no windows are present, skip this tile
    let Some(fetch_span) = adapt_fetch_to_extreme_windows(
        tile,
        tile_window_span,
        &mode,
        chrom_len as u32,
        opt.fragment_lengths.max_fragment_length as u64,
    )?
    else {
        return Ok(counter);
    };
    let (fetch_from, fetch_to) = fetch_span.try_to_i64()?.as_tuple();

    reader
        .fetch((tile.tid, fetch_from, fetch_to))
        .context(format!("fetch {} {}-{}", &tile.chr, fetch_from, fetch_to))?;

    // Prepare CP for tile core length
    let core_len = tile.core_end() - tile.core_start();
    let mut cp = Coverage::new(core_len);

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let lengths = opt.fragment_lengths.clone();
        move |f: &FragmentWithSegments| lengths.contains(f.len())
    };

    let gc_tag_bytes = gc_tag.map(|t| t.as_bytes().to_vec());

    // Create fragment iterator
    let unpaired = opt.unpaired.reads_are_fragments;
    let include_read_fn: Box<dyn Fn(&Record) -> bool + Send + Sync> = if unpaired {
        let min_mapq = opt.min_mapq;
        Box::new(move |r: &Record| default_include_read_unpaired(r, min_mapq))
    } else {
        let min_mapq = opt.min_mapq;
        let require_proper_pair = opt.require_proper_pair;
        Box::new(move |r: &Record| {
            default_include_read_paired_end(r, require_proper_pair, min_mapq)
        })
    };

    let mut iter = fragments_with_segments_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        move |rec| include_read_fn(rec),
        1,
        !opt.ignore_gap,
        gc_tag_bytes.as_deref(),
        fragment_filter,
        unpaired,
    )
    .with_local_counters();

    // Iterate fragments and add coverage
    // Separate branches for with/without GC correction
    if let Some(gc_corrector) = gc_corrector_opt {
        let gc_prefixes = gc_prefixes_opt
            .as_ref()
            .context("GC prefix sums missing despite GC correction being enabled")?;
        let fetch_start = tile.fetch_start();
        let fetch_end = tile.fetch_end();
        for fragment_res in iter.by_ref() {
            let fragment = fragment_res.context("reading fragment")?;
            let normalization_length = normalization_length_for_fragment(&fragment, opt)?;
            let base_weight = calculate_base_weight(normalization_length);

            if fragment.start() < fetch_start || fragment.end() > fetch_end {
                // Fragment won't overlap the tile core (assuming correct max_fragment_length halo!)
                // Note that more fragments (smaller than max_fragment_length) could be outside the tiles
                continue;
            }

            let fetch_relative_fragment = fragment
                .interval
                .try_to_u64()?
                .shift_left(fetch_start as u64)?;

            let gc_weight =
                match gc_corrector.correct_fragment(fetch_relative_fragment, &gc_prefixes)? {
                    Some(weight) => weight,
                    None => {
                        counter.gc_failed_fragments += 1;
                        if opt.gc.neutralize_invalid_gc {
                            1.0
                        } else {
                            continue;
                        }
                    }
                };

            // Clip and add to tile core coverage (segments respected)
            let was_counted = add_fragment_clipped_to_core(
                &mut cp,
                &fragment,
                base_weight * gc_weight,
                tile.core,
            )?;

            if was_counted {
                counter.base.counted_fragments += 1;
                if let Some(normalization_length) = normalization_length
                    && fragment_is_owned_by_tile_for_normalization_stats(&fragment, tile)
                {
                    counter.tile_owned_normalization_fragments += 1;
                    counter.tile_owned_normalization_length_sum += normalization_length as u64;
                }
            }
        }
    } else if gc_tag.is_some() {
        for fragment_res in iter.by_ref() {
            let fragment = fragment_res.context("reading fragment")?;
            let normalization_length = normalization_length_for_fragment(&fragment, opt)?;
            let base_weight = calculate_base_weight(normalization_length);

            let gc_weight = match fragment.gc_tag.classify()? {
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
            };

            let was_counted = add_fragment_clipped_to_core(
                &mut cp,
                &fragment,
                base_weight * gc_weight,
                tile.core,
            )?;

            if was_counted {
                counter.base.counted_fragments += 1;
                if let Some(normalization_length) = normalization_length
                    && fragment_is_owned_by_tile_for_normalization_stats(&fragment, tile)
                {
                    counter.tile_owned_normalization_fragments += 1;
                    counter.tile_owned_normalization_length_sum += normalization_length as u64;
                }
            }
        }
    } else {
        for fragment_res in iter.by_ref() {
            let fragment = fragment_res.context("reading fragment")?;
            let normalization_length = normalization_length_for_fragment(&fragment, opt)?;
            let base_weight = calculate_base_weight(normalization_length);

            // Clip and add to tile core coverage (segments respected)
            let was_counted =
                add_fragment_clipped_to_core(&mut cp, &fragment, base_weight, tile.core)?;

            if was_counted {
                counter.base.counted_fragments += 1;
                if let Some(normalization_length) = normalization_length
                    && fragment_is_owned_by_tile_for_normalization_stats(&fragment, tile)
                {
                    counter.tile_owned_normalization_fragments += 1;
                    counter.tile_owned_normalization_length_sum += normalization_length as u64;
                }
            }
        }
    }

    // Clear up memory before finalizing coverage
    drop(gc_prefixes_opt);

    // Finalize coverage
    cp.finalize_coverage(true);

    // Clamp almost-zero coverages to zero to avoid any f32 roundoff error
    // NOTE: Must come before scaling!
    // NOTE: If we add other normalizations, we must consider its effect here!
    if let Some(cov_mut) = cp.coverage_mut() {
        clamp_finite_coverage_below_to_zero(cov_mut, internal_residual_coverage_floor(opt));
    }

    // Apply per-bin scaling (in-place)
    if !scaling_chr.is_empty()
        && let Some(cov_mut) = cp.coverage_mut()
    {
        apply_scaling_to_coverage_in_place(cov_mut, tile.core_start(), scaling_chr);
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    match mode {
        TileMode::Positional {
            windows,
            out_path,
            indexed,
        } => {
            // Add blacklist late and clipped to the tile core to minimize memory
            // Use binary search to jump to the first overlapping interval
            add_clipped_blacklist_to_cp(&mut cp, tile, !blacklist_chr.is_empty(), blacklist_chr)?;

            // Prepare compressed writer (zstd) for this tile
            let mut w = open_zstd_auto_writer(&out_path, 3, None)?;

            let cov = cp
                .coverage()
                .context("tile coverage missing after finalization")?;
            let mask = cp.blacklist_mask();

            // Write tile data to disk

            match windows {
                None => {
                    // Whole positional coverage for the tile core
                    write_bedgraph_runs(
                        &tile.chr,
                        cov,
                        mask,
                        0,
                        cov.len(),
                        tile.core_start() as u64,
                        decimals,
                        opt.keep_zero_runs,
                        &mut w,
                    )?;
                }
                Some(win_chr) => {
                    for window in overlapping_windows_for_tile(win_chr, tile, tile_window_span) {
                        let window_start = window.start();
                        let window_end = window.end();
                        // Keep the original window index from the BED input so the
                        // indexed positional output and downstream reducers stay aligned
                        let original_idx = window.idx();
                        let core_interval =
                            Interval::new(tile.core_start() as u64, tile.core_end() as u64)?;
                        let local_overlap = if let Some(local_overlap) =
                            clip_interval_to_core_and_localize(
                                Interval::new(window_start, window_end)?,
                                core_interval,
                            )? {
                            local_overlap
                        } else {
                            continue;
                        };

                        if indexed {
                            write_windowed_runs(
                                &tile.chr,
                                cov,
                                mask,
                                local_overlap.local_start_idx,
                                local_overlap.local_end_idx,
                                tile.core_start() as u64,
                                Some(original_idx),
                                decimals,
                                opt.keep_zero_runs,
                                &mut w,
                            )?;
                        } else {
                            write_windowed_runs(
                                &tile.chr,
                                cov,
                                mask,
                                local_overlap.local_start_idx,
                                local_overlap.local_end_idx,
                                tile.core_start() as u64,
                                None,
                                decimals,
                                opt.keep_zero_runs,
                                &mut w,
                            )?;
                        }
                    }
                }
            }

            w.flush()?;
        }

        TileMode::AggregatesByBed {
            windows,
            masked,
            partials_out,
            cross_idx_out,
        } => {
            // Add blacklist (clipped) and build indexes once
            add_clipped_blacklist_to_cp(&mut cp, tile, masked, blacklist_chr)?;
            cp.build_indexes(false)?;

            // Borrow indexes and mask once
            let psum_all = cp
                .psum_all_ref()
                .ok_or_else(|| anyhow::anyhow!("psum_all missing"))?;
            let psum_unmasked = cp.psum_unmasked_ref();
            let psum_cnt_unmasked = cp.psum_unmasked_count_ref();
            let mask: Option<&[u8]> = cp.blacklist_mask();
            let wants_summary_stats = opt.per_window.is_summary_stats();
            let summary_prefixes = if wants_summary_stats {
                Some(build_summary_prefixes(&cp)?)
            } else {
                None
            };

            // Writers: compressed partials and cross sidecar
            let mut w_part = open_zstd_auto_writer(&partials_out, 3, None)?;
            let mut w_cross = open_zstd_auto_writer(&cross_idx_out, 3, None)?;

            let core_start_abs = tile.core_start() as u64;
            let core_end_abs = tile.core_end() as u64;

            // Walk only windows overlapping the core (already start-sorted)
            for window in overlapping_windows_for_tile(windows, tile, tile_window_span) {
                let window_start_abs = window.start();
                let window_end_abs = window.end();
                // This is the original window index, not just the current loop position
                let idx = window.idx();
                let core_interval =
                    Interval::new(tile.core_start() as u64, tile.core_end() as u64)?;
                let local_overlap = if let Some(local_overlap) = clip_interval_to_core_and_localize(
                    Interval::new(window_start_abs, window_end_abs)?,
                    core_interval,
                )? {
                    local_overlap
                } else {
                    continue;
                };

                // Classify as internal (fully inside core) vs boundary (crosses tile core boundary)
                let crosses_boundary =
                    !(window_start_abs >= core_start_abs && window_end_abs <= core_end_abs);

                // Always write a partial row. The reducer merges contributions by `orig_idx`, but
                // ordinary BED windows preserve their original file indices, so callers should not
                // assume the final row order is increasing `orig_idx` unless the windows were
                // explicitly reindexed into coordinate order upstream
                // Internal windows won't appear in the cross-index -> reducer expects 1 contribution
                // Boundary windows will appear in each crossed tile's cross-index -> reducer expects N
                if wants_summary_stats {
                    let raw_stats = coverage_summary_and_counts(
                        local_overlap.local_start_idx,
                        local_overlap.local_end_idx,
                        masked,
                        psum_all,
                        psum_unmasked,
                        psum_cnt_unmasked,
                        mask,
                        summary_prefixes
                            .as_ref()
                            .expect("summary prefixes must exist for summary-stats"),
                    );
                    writeln!(
                        w_part,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        idx,
                        raw_stats.coverage_sum,
                        raw_stats.eligible_positions,
                        raw_stats.blacklisted_positions,
                        raw_stats.nonzero_positions,
                        raw_stats.coverage_sum_of_squares
                    )?;
                } else {
                    // Average/total reducers only need the first raw moment plus position counts,
                    // so keep the per-tile partials narrow in the common non-summary path.
                    let (coverage_sum, eligible_positions, blacklisted_positions) =
                        coverage_sum_and_counts(
                            local_overlap.local_start_idx,
                            local_overlap.local_end_idx,
                            masked,
                            psum_all,
                            psum_unmasked,
                            psum_cnt_unmasked,
                            mask,
                        );
                    writeln!(
                        w_part,
                        "{}\t{}\t{}\t{}",
                        idx, coverage_sum, eligible_positions, blacklisted_positions
                    )?;
                }
                if crosses_boundary {
                    // Cross-index lists the window's orig_idx for the reducer
                    writeln!(w_cross, "{}", idx)?;
                }
            }

            w_part.flush()?;
            w_cross.flush()?;
        }

        TileMode::AggregatesBySize {
            window_bp,
            masked,
            finals_out,
            partials_out,
            cross_idx_out,
            guaranteed_aligned,
        } => {
            // Add blacklist late and clipped to the tile core to minimize memory
            // Use binary search to jump to the first overlapping interval
            add_clipped_blacklist_to_cp(&mut cp, tile, masked, blacklist_chr)?;

            // Build prefix-sum indexes for fast per-window queries
            cp.build_indexes(false)?;

            // Own copies of the prefix arrays and optional mask to avoid long-lived borrows
            let psum_all = cp
                .psum_all_ref()
                .ok_or_else(|| anyhow::anyhow!("psum_all missing"))?;
            let psum_unmasked = cp.psum_unmasked_ref();
            let psum_cnt_unmasked = cp.psum_unmasked_count_ref();
            let mask: Option<&[u8]> = cp.blacklist_mask();
            let wants_summary_stats = opt.per_window.is_summary_stats();
            let summary_prefixes = if wants_summary_stats {
                Some(build_summary_prefixes(&cp)?)
            } else {
                None
            };

            // Determine the fixed-size windows that overlap the tile core
            let core_start_abs = tile.core_start() as u64;
            let core_end_abs = tile.core_end() as u64;
            let first_bin_idx = core_start_abs / window_bp;
            let last_bin_idx = (core_end_abs.saturating_sub(1)) / window_bp;

            if guaranteed_aligned && !opt.restores_mean_after_length_normalization() {
                // FAST PATH: Every bin that touches the core is fully contained in this core
                // so there is no cross-tile reduction later. Scalar actions can write their
                // final value immediately, while summary-stats can write their final derived row
                // immediately from the exact raw aggregates for that bin

                let mut w_fin = open_zstd_auto_writer(&finals_out, 3, None)?;

                for bin_idx in first_bin_idx..=last_bin_idx {
                    let bin_start = bin_idx * window_bp;
                    let bin_end = (bin_idx + 1) * window_bp;

                    // Intersect with core (alignment ensures this equals the bin for non-terminal tiles).
                    let core_interval =
                        Interval::new(tile.core_start() as u64, tile.core_end() as u64)?;
                    let local_overlap = if let Some(local_overlap) =
                        clip_interval_to_core_and_localize(
                            Interval::new(bin_start, bin_end)?,
                            core_interval,
                        )? {
                        local_overlap
                    } else {
                        continue;
                    };

                    // Write the full bin coordinates, not just the overlap inside this tile core.
                    // Aligned tiles guarantee one final row per bin, but the last bin on a
                    // chromosome may still need clipping at the chromosome end
                    let bin_interval = Interval::new(bin_start, bin_end.min(core_end_abs))?;

                    if wants_summary_stats {
                        let raw_stats = coverage_summary_and_counts(
                            local_overlap.local_start_idx,
                            local_overlap.local_end_idx,
                            masked,
                            psum_all,
                            psum_unmasked,
                            psum_cnt_unmasked,
                            mask,
                            summary_prefixes
                                .as_ref()
                                .expect("summary prefixes must exist for summary-stats"),
                        );
                        let stats = derive_summary_stats(
                            bin_interval.len(),
                            raw_stats.blacklisted_positions,
                            raw_stats.eligible_positions,
                            raw_stats.nonzero_positions,
                            raw_stats.coverage_sum,
                            raw_stats.coverage_sum_of_squares,
                        )?;
                        write_summary_stats_row(
                            &mut w_fin,
                            &tile.chr,
                            bin_interval,
                            stats,
                            decimals,
                        )?;
                    } else {
                        let (coverage_sum, eligible_positions, blacklisted_positions) =
                            coverage_sum_and_counts(
                                local_overlap.local_start_idx,
                                local_overlap.local_end_idx,
                                masked,
                                psum_all,
                                psum_unmasked,
                                psum_cnt_unmasked,
                                mask,
                            );
                        let unmasked_span_bp = local_overlap.clipped_abs_interval.len();
                        let value = finalize_value(
                            coverage_sum,
                            eligible_positions,
                            unmasked_span_bp,
                            masked,
                            &opt.per_window,
                        );
                        let value = round_to(value, decimals);

                        write_final_row(
                            &mut w_fin,
                            &tile.chr,
                            bin_interval,
                            value,
                            blacklisted_positions,
                            decimals,
                        )?;
                    }
                }

                w_fin.flush()?;
            } else {
                let mut w_part = open_zstd_auto_writer(&partials_out, 3, None)?;
                let mut w_cross = if guaranteed_aligned {
                    None
                } else {
                    Some(open_zstd_auto_writer(&cross_idx_out, 3, None)?)
                };

                for bin_idx in first_bin_idx..=last_bin_idx {
                    let bin_start = bin_idx * window_bp;
                    let bin_end = (bin_idx + 1) * window_bp;

                    // Intersect with core (alignment ensures this equals the bin for non-terminal tiles).
                    let core_interval =
                        Interval::new(tile.core_start() as u64, tile.core_end() as u64)?;
                    let local_overlap = if let Some(local_overlap) =
                        clip_interval_to_core_and_localize(
                            Interval::new(bin_start, bin_end)?,
                            core_interval,
                        )? {
                        local_overlap
                    } else {
                        continue;
                    };

                    // Write the full bin start/end, not just this tile's overlap.
                    // The reducer merges partial rows by `bin_start` from the cross-index sidecar,
                    // so this row must keep the original bin coordinates even when this tile only
                    // contributes part of that bin.
                    if wants_summary_stats {
                        let raw_stats = coverage_summary_and_counts(
                            local_overlap.local_start_idx,
                            local_overlap.local_end_idx,
                            masked,
                            psum_all,
                            psum_unmasked,
                            psum_cnt_unmasked,
                            mask,
                            summary_prefixes
                                .as_ref()
                                .expect("summary prefixes must exist for summary-stats"),
                        );
                        writeln!(
                            w_part,
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            bin_start,
                            bin_end,
                            raw_stats.coverage_sum,
                            raw_stats.eligible_positions,
                            raw_stats.blacklisted_positions,
                            raw_stats.nonzero_positions,
                            raw_stats.coverage_sum_of_squares
                        )?;
                    } else {
                        // Keep non-summary by-size partials to the fields the scalar reducer
                        // actually consumes so large runs do not write dead summary-stat columns.
                        let (coverage_sum, eligible_positions, blacklisted_positions) =
                            coverage_sum_and_counts(
                                local_overlap.local_start_idx,
                                local_overlap.local_end_idx,
                                masked,
                                psum_all,
                                psum_unmasked,
                                psum_cnt_unmasked,
                                mask,
                            );
                        writeln!(
                            w_part,
                            "{}\t{}\t{}\t{}\t{}",
                            bin_start,
                            bin_end,
                            coverage_sum,
                            eligible_positions,
                            blacklisted_positions
                        )?;
                    }

                    // Mark cross-boundary bins (not fully inside the core) so reducer expects >1 contributions
                    let fully_inside = (bin_start >= core_start_abs) && (bin_end <= core_end_abs);
                    if !fully_inside {
                        let w_cross = w_cross.as_mut().context(
                            "guaranteed_aligned by-size partial writing encountered a cross-boundary bin without a cross-index writer",
                        )?;
                        writeln!(w_cross, "{}", bin_start)?;
                    }
                }
                w_part.flush()?;
                if let Some(mut w_cross) = w_cross {
                    w_cross.flush()?;
                }
            }
        }
    }

    Ok(counter)
}

/// Return whether this fragment is owned by the current tile for normalization statistics.
///
/// Coverage still uses every overlapping tile core contribution. This helper is only for the
/// run-level normalization-length statistics, where each fragment must be counted once even when
/// fetch halos make it visible in neighboring tiles. We assign ownership by fragment start because
/// the current tile builder partitions each chromosome into contiguous, non-overlapping cores, so a
/// counted fragment start belongs to exactly one tile core.
///
/// It must not be reused as a generic fragment-counting rule without re-checking that tiling
/// invariant and the scientific model.
#[inline]
fn fragment_is_owned_by_tile_for_normalization_stats(
    fragment: &FragmentWithSegments,
    tile: &Tile,
) -> bool {
    let fragment_start = fragment.start();
    fragment_start >= tile.core_start() && fragment_start < tile.core_end()
}

/// Add blacklist clipped to tile coordinates to the coverage prefix object
fn add_clipped_blacklist_to_cp(
    cp: &mut Coverage,
    tile: &Tile,
    masked: bool,
    blacklist_chr: &[Interval<u64>],
) -> Result<()> {
    // Add blacklist late and clipped to the tile core to minimize memory
    // Use binary search to jump to the first overlapping interval
    if masked && !blacklist_chr.is_empty() {
        let core_start_abs = tile.core_start() as u64;
        let core_end_abs = tile.core_end() as u64;

        // Find first interval with end > core_start
        let mut i = blacklist_chr.partition_point(|interval| interval.end() <= core_start_abs);

        if i < blacklist_chr.len() {
            let mut clipped: Vec<Interval<u64>> = Vec::new();

            // Walk only the intervals that can overlap the tile core
            while i < blacklist_chr.len() {
                let blacklist_start = blacklist_chr[i].start();
                let blacklist_end = blacklist_chr[i].end();
                if blacklist_start >= core_end_abs {
                    break; // Remaining intervals start after the core band
                }

                let overlap_start_abs = blacklist_start.max(core_start_abs);
                let overlap_end_abs = blacklist_end.min(core_end_abs);
                if overlap_start_abs < overlap_end_abs {
                    // Convert to tile‐local coordinates
                    let local_start = (overlap_start_abs as u32) - tile.core_start();
                    let local_end = (overlap_end_abs as u32) - tile.core_start();
                    clipped.push(Interval::new(local_start as u64, local_end as u64)?);
                }
                i += 1;
            }

            if !clipped.is_empty() {
                cp.set_blacklist_mask(&clipped)?;
            }
        }
    }

    Ok(())
}

/// Adds a fragment's coverage contribution into the tile-local accumulator.
///
/// Segmented fragments are processed segment by segment, while simple fragments are clipped once;
/// in both cases the coordinates are translated into the tile's local frame before they are added.
/// The caller must provide a coverage array sized to the tile core.
///
/// # Parameters
/// - `cp`: Tile-local coverage structure to update.
/// - `fragment`: Fragment carrying absolute coordinates and optional segments.
/// - `weight`: Weight applied when inserting the fragment.
/// - `core_start`: Inclusive start of the tile core in absolute coordinates.
/// - `core_end`: Exclusive end of the tile core in absolute coordinates.
///
/// # Returns
/// - `Ok(true)` if the fragment contributes at least one base to the tile core.
/// - `Ok(false)` if every segment falls outside the core after clipping.
/// - An error when the coverage accumulator rejects the update.
#[inline]
pub fn add_fragment_clipped_to_core(
    cp: &mut Coverage,
    fragment: &FragmentWithSegments,
    weight: f64,
    core_interval: Interval<u32>,
) -> Result<bool> {
    // Use explicit segments if present
    let mut counted = false;
    let core_start = core_interval.start();
    let to_core_local = |interval: Interval<u32>| -> Result<Interval<u32>> {
        interval.shift_left(core_start).map_err(anyhow::Error::from)
    };
    if let Some(segments) = &fragment.segments {
        for segment in segments {
            let Some(clipped_interval) = segment.clip_to(core_interval) else {
                // Skips fragments completely outside tile
                continue;
            };
            // Shift to tile-local coordinates
            let local_interval = to_core_local(clipped_interval)?;
            let local = Fragment {
                tid: fragment.tid,
                interval: local_interval,
                gc_tag: Default::default(),
            };
            cp.add_fragment_weighted(local, weight)?;
            counted = true;
        }
    } else {
        // No explicit segments -> treat as one span (this already encodes your include_inter_mate_gap policy)
        // Skips fragments completely outside tile
        if let Some(clipped_interval) = fragment.interval.clip_to(core_interval) {
            // Shift to tile-local coordinates
            let local_interval = to_core_local(clipped_interval)?;
            let local = Fragment {
                tid: fragment.tid,
                interval: local_interval,
                gc_tag: Default::default(),
            };

            cp.add_fragment_weighted(local, weight)?;
            counted = true;
        }
    }
    Ok(counted)
}

/// Compute the intrinsic normalization length for one fragment when length normalization is active.
///
/// Technical details:
/// - Plain fragments use the full fragment span length `[start, end)`.
/// - Segment-aware fragments use the summed length of their counted reference segments, so deleted
///   or skipped reference bases do not dilute the fragment's total mass below `1.0`.
/// - GC correction and genomic scaling are applied later as multiplicative factors on top of the
///   resulting base weight.
///
/// Parameters
/// ----------
/// - `fragment`:
///     Fragment span and optional counted reference segments for this molecule.
/// - `opt`:
///     Command configuration that controls whether length normalization is enabled.
///
/// Returns
/// -------
/// - `normalization_length`:
///     `Some(denominator)` when length normalization is active, otherwise `None`.
///
/// Errors
/// ------
/// Returns an error if `--normalize-by-length` is enabled but the fragment has zero countable
/// length. That would indicate an internal inconsistency in the counted spans.
#[inline]
fn normalization_length_for_fragment(
    fragment: &FragmentWithSegments,
    opt: &FCoverageConfig,
) -> Result<Option<u32>> {
    if !opt.uses_length_normalization() {
        return Ok(None);
    }

    let normalization_length = if let Some(segments) = &fragment.segments {
        segments.iter().map(|segment| segment.len()).sum::<u32>()
    } else {
        fragment.len()
    };

    if normalization_length == 0 {
        bail!("normalize-by-length encountered a fragment with zero countable length");
    }
    Ok(Some(normalization_length))
}

/// Compute the intrinsic per-base fragment weight before GC correction or genomic scaling.
#[inline]
fn calculate_base_weight(normalization_length: Option<u32>) -> f64 {
    normalization_length.map_or(1.0, |length| 1.0 / length as f64)
}

/// Smallest positive intrinsic per-base weight a counted fragment base can have in this run.
///
/// This is evaluated before GC correction or genomic scaling.
#[inline]
fn minimum_positive_base_weight(opt: &FCoverageConfig) -> f64 {
    if opt.uses_length_normalization() {
        1.0 / opt.fragment_lengths.max_fragment_length as f64
    } else {
        1.0
    }
}

/// Smallest positive GC multiplier that can be considered usable in this run.
#[inline]
fn minimum_positive_gc_weight(opt: &FCoverageConfig) -> f64 {
    if opt.gc.gc_file.is_some() || opt.gc.gc_tag.is_some() {
        MIN_REASONABLE_GC_WEIGHT as f64
    } else {
        1.0
    }
}

/// Smallest real positive pre-scaling support that one counted position can receive.
///
/// This must be updated whenever a new pre-scaling weighting or normalization can lower the
/// per-position support below the current bound. The tests intentionally exercise the current
/// GC and length-normalization combinations so future changes have to revisit this derivation.
#[inline]
fn minimum_positive_pre_scaling_support(opt: &FCoverageConfig) -> f64 {
    minimum_positive_base_weight(opt) * minimum_positive_gc_weight(opt)
}

/// Internal cleanup floor for fake support created by floating-point add/subtract residue.
///
/// The floor stays strictly below the smallest real positive pre-scaling support for the active
/// argument combination, so it can only remove values that should be impossible in exact
/// arithmetic for this run.
#[inline]
fn internal_residual_coverage_floor(opt: &FCoverageConfig) -> f32 {
    (minimum_positive_pre_scaling_support(opt) / 2.0) as f32
}

#[cfg(test)]
mod tests {
    include!("fcoverage_tests.rs");
}
