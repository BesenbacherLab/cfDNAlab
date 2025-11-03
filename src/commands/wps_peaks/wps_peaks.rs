use crate::commands::fcoverage::reducer::{
    reduce_aggregates_by_size_with_cross_index_for_chr, reduce_bed_with_cross_index_for_chr,
};
use crate::commands::fcoverage::tiling::{concat_aligned_size_tile_finals, merge_positional_tiles};
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use crate::commands::wps::wps::wps_for_tile;
use crate::commands::wps_peaks::config::WPSPeaksConfig;
use crate::shared::tiled_run::{
    Tile, TileMode, TileWindowSpan, build_tiles, make_temp_dir, precompute_tile_window_spans,
};
use crate::shared::writers::open_zstd_auto_writer;
use crate::{
    commands::cli_common::{
        WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
        resolve_chromosomes_and_contigs,
    },
    commands::counters::FCoverageCounters,
    shared::{bed::load_windows_from_bed, thread_pool::init_global_pool},
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::io::Write;
use std::{sync::Arc, time::Instant};

/// Execute the windowed protection scores pipeline end-to-end.
///
/// Implementation details:
/// - Resolves chromosomes, prepares IO state, then iterates tiles in parallel.
/// - Collects per-tile scores into temporary artefacts before merging them into the final
///   positional or aggregated outputs.
/// - Applies fragment length, blacklist, and optional scaling filters during iteration.
/// - Applies smoothing to calculated WPS values.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `wps` command.
///
/// Returns:
/// - `Ok(())` when positional and/or windowed outputs are written successfully.
///
/// Errors:
/// - Returns an error if the BAM cannot be read, auxiliary files are invalid, or writing outputs
///   fails at any stage.
pub fn run(opt: &WPSPeaksConfig) -> Result<()> {
    let start_time = Instant::now();
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.shared_args.chromosomes, &opt.shared_args.ioc)?;
    let window_opt = opt.shared_args.windows.resolve_windows();
    let prefix = opt.shared_args.output_prefix.trim();

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

    let windowed = matches!(window_opt, WindowSpec::Bed(_) | WindowSpec::Size(_));
    ensure!(
        !windowed || opt.per_window.is_some(),
        "when using --by-bed/--by-size, please also specify --per-window"
    );
    let per_window_wps_action = opt.per_window;

    // Create output directory
    ensure_output_dir(&opt.shared_args.ioc.output_dir)?;

    if opt.shared_args.blacklist.is_some() {
        println!("Start: Loading blacklists");
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
            println!("Start: Loading window coordinates");
            let wds = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None)?;
            if matches!(
                per_window_wps_action,
                Some(CoverageWindowAction::OnlyIncludeThesePositionsUnique)
            ) {
                // Merge in-place to avoid double memory-usage
                println!("Start: Merging overlapping/touching windows");
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
        println!("Start: Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> =
        load_scaling_map(&opt.shared_args.scale_genome, &chromosomes, &contigs)?;

    let has_scaling = opt.shared_args.scale_genome.scaling_factors.is_some();

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
    let temp_dir = make_temp_dir(&opt.shared_args.ioc.output_dir, prefix)
        .context("create per-run temp dir")?;

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

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(&tiles, |chr| {
        windows_lookup
            .and_then(|m| m.get(chr).map(|w| w.as_slice()))
            .unwrap_or(&[])
    }));

    // Where per-tile files go
    let positional_prefix = format!("{prefix}.pos");
    let partials_prefix = format!("{prefix}.part");
    let finals_prefix = format!("{prefix}.fin");

    // Faster to convert to &str once
    let positional_prefix = positional_prefix.as_str();
    let partials_prefix = partials_prefix.as_str();
    let finals_prefix = finals_prefix.as_str();

    // Create filenames of final outputs
    let final_bedgraph_pos_name = format!("{prefix}.wps.per_position.bedgraph.zst");
    let final_tsv_pos_name = format!("{prefix}.wps.per_position_per_window.tsv.zst");
    let final_avg_name = format!("{prefix}.wps.avg.tsv.zst");
    let final_total_name = format!("{prefix}.wps.total.tsv.zst");

    // Get decimals to use
    let decimals_to_use: i32 = if windowed {
        match per_window_wps_action.expect("per-window action required when windowed") {
            CoverageWindowAction::Average | CoverageWindowAction::Total => {
                opt.shared_args.decimals as i32
            }
            CoverageWindowAction::OnlyIncludeThesePositionsUnique
            | CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                if has_scaling {
                    opt.shared_args.decimals as i32
                } else {
                    0
                }
            }
        }
    } else {
        if has_scaling {
            opt.shared_args.decimals as i32
        } else {
            0
        }
    };

    let total_tiles = tiles.len();

    // Create progress bar
    let pb = Arc::new(ProgressBar::new(total_tiles as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    // Configure global thread‐pool size
    init_global_pool(opt.shared_args.ioc.n_threads as usize)?;

    let mut global_counter = FCoverageCounters::default();

    println!("Start: Counting per tile");

    pb.set_position(0);

    let tile_window_spans_for_threads = tile_window_spans.clone();

    let tile_results: Vec<FCoverageCounters> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<FCoverageCounters> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            // Per-chrom projections
            let windows_chr: Option<&[(u64, u64, u64)]> = windows_map
                .as_ref()
                .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice()));
            let blacklist_chr: &[(u64, u64)] = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr: &[(u64, u64, f32)] = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            // Decide tile mode and file name
            let (action_prefix, extensions) = if windowed {
                let action =
                    per_window_wps_action.expect("per-window action required when windowed");
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
                match (
                    &window_opt,
                    per_window_wps_action.expect("per-window action required when windowed"),
                ) {
                    (WindowSpec::Bed(_), CoverageWindowAction::OnlyIncludeThesePositionsUnique) => {
                        TileMode::Positional {
                            windows: windows_chr,
                            out_path,
                            indexed: false,
                        }
                    }
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
                    ) => {
                        let wchr = windows_chr.expect("windows required for aggregates");
                        TileMode::AggregatesByBed {
                            windows: wchr,
                            masked: true,
                            partials_out,
                            cross_idx_out,
                        }
                    }
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

            let counter = peaks_per_tile(
                &opt,
                tile,
                tile_span.as_ref(),
                blacklist_chr,
                scaling_chr,
                mode,
                decimals_to_use,
            )?;
            pb.inc(1);
            Ok(counter)
        })
        .collect::<anyhow::Result<_>>()?;

    pb.finish_with_message("| Finished counting");

    // Collect counters
    for counter in tile_results {
        global_counter += counter;
    }

    println!("Start: Merging temporary tile files to final output");

    // Merge temporary output files and
    // reduce windows present in multiple tiles

    let final_out_path = if !windowed {
        // Whole-genome positional coverage
        merge_positional_tiles(
            &temp_dir,
            &opt.shared_args.ioc.output_dir,
            &chromosomes,
            positional_prefix,
            final_bedgraph_pos_name.as_str(),
        )?
    } else {
        let action = per_window_wps_action.expect("per-window action required when windowed");
        match action {
            CoverageWindowAction::OnlyIncludeThesePositionsUnique => {
                // Windowed positional (unique and non-indexed)
                merge_positional_tiles(
                    &temp_dir,
                    &opt.shared_args.ioc.output_dir,
                    &chromosomes,
                    positional_prefix,
                    final_bedgraph_pos_name.as_str(),
                )?
            }
            CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                // Windowed positional with orig_idx column
                merge_positional_tiles(
                    &temp_dir,
                    &opt.shared_args.ioc.output_dir,
                    &chromosomes,
                    positional_prefix,
                    final_tsv_pos_name.as_str(),
                )?
            }
            CoverageWindowAction::Average | CoverageWindowAction::Total => {
                // Per-chrom reduce of partials into final aggregates
                let final_path = opt.shared_args.ioc.output_dir.join(match action {
                    CoverageWindowAction::Average => final_avg_name.as_str(),
                    CoverageWindowAction::Total => final_total_name.as_str(),
                    _ => unreachable!(),
                });

                // Header value-column name
                let value_col = match action {
                    CoverageWindowAction::Average => "avg_coverage",
                    CoverageWindowAction::Total => "total_coverage",
                    _ => unreachable!(),
                };

                let header = format!(
                    "chromosome\tstart\tend\t{}\tblacklisted_positions",
                    value_col
                );

                // Reduce by window source
                match &window_opt {
                    WindowSpec::Bed(_) => {
                        let mut positional_writer = open_zstd_auto_writer(
                            &final_path,
                            3,
                            Some(opt.shared_args.ioc.n_threads as u32),
                        )?;

                        // Write header
                        writeln!(positional_writer, "{}", header)?;

                        let win_map = windows_map
                            .as_ref()
                            .expect("windows_map present for --by-bed");
                        for chr in &chromosomes {
                            if let Some(wchr) = win_map.get(chr) {
                                reduce_bed_with_cross_index_for_chr(
                                    chr,
                                    &temp_dir,
                                    partials_prefix,
                                    wchr.as_slice(),
                                    true,
                                    action,
                                    decimals_to_use,
                                    &mut positional_writer,
                                )?;
                            }
                        }
                        positional_writer.flush()?;
                    }
                    WindowSpec::Size(_) => {
                        if tile_and_window_boundaries_align {
                            let _ = concat_aligned_size_tile_finals(
                                &temp_dir,
                                &opt.shared_args.ioc.output_dir,
                                &chromosomes,
                                finals_prefix,
                                match action {
                                    CoverageWindowAction::Average => final_avg_name.as_str(),
                                    CoverageWindowAction::Total => final_total_name.as_str(),
                                    _ => unreachable!(),
                                },
                                &header,
                            )?;
                        } else {
                            let mut positional_writer = open_zstd_auto_writer(
                                &final_path,
                                3,
                                Some(opt.shared_args.ioc.n_threads as u32),
                            )?;

                            // Write header
                            writeln!(positional_writer, "{}", header)?;

                            for chr in &chromosomes {
                                let chrom_len = contigs
                                    .contigs
                                    .get(chr)
                                    .map(|&(_, len)| len as u64)
                                    .expect("missing contig length");
                                reduce_aggregates_by_size_with_cross_index_for_chr(
                                    chr,
                                    &temp_dir,
                                    partials_prefix,
                                    true,
                                    action,
                                    chrom_len,
                                    decimals_to_use,
                                    &mut positional_writer,
                                )?;
                            }
                            positional_writer.flush()?;
                        }
                    }
                    _ => unreachable!(),
                }

                final_path
            }
        }
    };
    println!("Saved output to: {:?}", final_out_path);

    let keep_temp = false; // TODO: Make cli arg behind a feature for dev purposes?
    if !keep_temp {
        if let Err(e) = std::fs::remove_dir_all(&temp_dir) {
            eprintln!(
                "warning: failed to remove temp dir {}: {}",
                temp_dir.display(),
                e
            );
        }
    } else {
        eprintln!("kept temp tiles in {}", temp_dir.display());
    }
    println!("");
    println!("Statistics");
    println!("----------");

    // Print summary statistics and execution time
    let elapsed = start_time.elapsed();
    println!("  Total reads: {}", global_counter.base.total_reads);
    println!(
        "  Initially accepted reads: {} ({:.2}%, forward: {}, reverse: {})",
        global_counter.base.accepted_forward + global_counter.base.accepted_reverse,
        (global_counter.base.accepted_forward + global_counter.base.accepted_reverse) as f64
            / global_counter.base.total_reads as f64
            * 100.0,
        global_counter.base.accepted_forward,
        global_counter.base.accepted_reverse
    );
    // if opt.shared_args.gc.bin_by_gc {
    //     println!("GC-excluded reads: {}", global_counter.base.gc_excl);
    // }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.base.counted_fragments
    );
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

pub fn peaks_per_tile(
    opt: &WPSPeaksConfig,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_chr: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    mode: TileMode,
    decimals: i32,
) -> Result<FCoverageCounters> {
    let (counter, wps_values, mask) = wps_for_tile(
        &opt.shared_args,
        &None,
        false,
        tile,
        tile_window_span,
        blacklist_chr,
        scaling_chr,
        mode,
        decimals,
        true,
    )?;

    /*
    TODO:
    1) Normalize medians in 1kb windows, 2) smoothe with SavGol filter, 3) call peaks, 4) save peaks and/or calculate stats (interpeak-distance, counts, etc) (e.g. per window)
    */

    Ok(counter)
}
