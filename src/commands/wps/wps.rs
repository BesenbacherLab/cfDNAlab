use crate::commands::fcoverage::reducer::{
    reduce_aggregates_by_size_with_cross_index_for_chr, reduce_bed_with_cross_index_for_chr,
};
use crate::commands::fcoverage::tiling::{
    adapt_fetch_to_extreme_windows, concat_aligned_size_tile_finals, finalize_value,
    merge_positional_tiles,
};
use crate::commands::fcoverage::window_results::CoverageWindowAction;
use crate::commands::fcoverage::writers::{
    emit_bedgraph_runs, emit_windowed_runs, write_final_row,
};
use crate::commands::gc_bias::correct::{GCCorrector, load_gc_corrector};
use crate::commands::gc_bias::counting::build_gc_prefixes;
use crate::commands::wps::config::{WPSConfig, WPSSharedConfig};
use crate::shared::formatters::round_to;
use crate::shared::fragment::minimal_fragment::Fragment;
use crate::shared::fragment_iterator::fragments_from_bam;
use crate::shared::read::{default_include_read_paired_end, default_include_read_unpaired};
use crate::shared::reference::read_seq_in_range;
use crate::shared::scale_genome::apply_scaling_to_coverage_in_place;
use crate::shared::tiled_run::{
    Tile, TileMode, TileWindowSpan, build_tiles, make_temp_dir, overlapping_windows_for_tile,
    precompute_tile_window_spans,
};
use crate::shared::writers::open_zstd_auto_writer;
use crate::{
    commands::cli_common::{
        WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
        resolve_chromosomes_and_contigs,
    },
    commands::counters::FCoverageCounters,
    shared::{
        bam::create_chromosome_reader, bed::load_windows_from_bed, thread_pool::init_global_pool,
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
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
pub fn run(opt: &WPSConfig) -> Result<()> {
    let start_time = Instant::now();
    if opt.shared_args.unpaired.reads_are_fragments && opt.shared_args.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let (chromosomes, contigs) = resolve_chromosomes_and_contigs(
        &opt.shared_args.chromosomes,
        &opt.shared_args.ioc.bam.as_path(),
    )?;
    let prefix = opt.shared_args.output_prefix.trim();
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
            let wds = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
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

    // Load GC correction package if specified
    if opt.shared_args.gc.gc_file.is_some() {
        println!("Start: Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.shared_args.gc.gc_file.as_ref(),
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
                if has_scaling_or_correction {
                    opt.shared_args.decimals as i32
                } else {
                    0
                }
            }
        }
    } else {
        if has_scaling_or_correction {
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

    println!("Start: Calculating WPS per tile");
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

            let (counter, _, _) = wps_for_tile(
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
        "  Note: A few reads/fragments may be counted twice in the statistics (only) around the parallelization tiles."
    );
    println!(
        "  Initially accepted reads: {} ({:.2}%, forward: {}, reverse: {})",
        global_counter.base.accepted_forward + global_counter.base.accepted_reverse,
        (global_counter.base.accepted_forward + global_counter.base.accepted_reverse) as f64
            / global_counter.base.total_reads as f64
            * 100.0,
        global_counter.base.accepted_forward,
        global_counter.base.accepted_reverse
    );
    if opt.shared_args.gc.gc_file.is_some() || opt.shared_args.gc.gc_tag.is_some() {
        let gc_fail_action = if opt.shared_args.gc.drop_invalid_gc {
            "fragment skipped"
        } else {
            "fragment counted with weight 1.0"
        };
        println!(
            "  GC correction failures ({}): {}",
            gc_fail_action, global_counter.gc_failed_fragments
        );
    }
    if opt.shared_args.gc.gc_tag.is_some() && global_counter.gc_out_of_range_tags > 0 {
        println!(
            "  GC tag values outside [0, {:.0}] treated as invalid: {}",
            crate::shared::gc_tag::MAX_REASONABLE_GC_WEIGHT,
            global_counter.gc_out_of_range_tags
        );
    }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.base.counted_fragments
    );
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

/// Process one tile: pair reads, build coverage, and write/return outputs for this tile
pub fn wps_for_tile(
    opt: &WPSSharedConfig,
    per_window_wps_action: &Option<CoverageWindowAction>,
    keep_zero_runs: bool,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_chr: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<GCCorrector>,
    mode: TileMode,
    decimals: i32,
    extra_halo_bp: u32,
    return_wps_instead: bool, // Don't save and aggregate, just return the WPS values
) -> Result<(FCoverageCounters, Option<Vec<f32>>, Option<Vec<u8>>)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    debug_assert!(tid_check == tile.tid as u32);

    let mut counter = FCoverageCounters::default();

    let gc_prefixes_opt = if gc_corrector_opt.is_some() {
        let ref_2bit = match opt.ref_2bit.as_ref() {
            Some(r) => r,
            None => bail!("When GC correction is specified, --ref-2bit must also be specified"),
        };
        let seq_bytes = read_seq_in_range(
            &ref_2bit,
            &tile.chr,
            // NOTE: Need for full fetch span to get GC of overlapping fragments!
            (tile.fetch_start as usize)..(tile.fetch_end as usize),
        )?;
        Some(build_gc_prefixes(&seq_bytes))
    } else {
        None
    };

    // Adapt the fetch coordinates to the present windows (*in genomic-windowed mode!*)
    // When no windows are present, skip this tile
    let Some((fetch_from, fetch_to)) =
        adapt_fetch_to_extreme_windows(tile, tile_window_span, &mode, chrom_len as u32)
    else {
        return Ok((counter, None, None));
    };

    reader
        .fetch((tile.tid as i32, fetch_from as i64, fetch_to as i64))
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
    let dilated_start_abs = tile.core_start.saturating_sub(context_left);
    let dilated_end_abs = ((tile.core_end as u64) + context_right as u64).min(chrom_len) as u32;
    if dilated_start_abs >= dilated_end_abs {
        return Ok((counter, None, None));
    }

    let dilated_span_len = (dilated_end_abs - dilated_start_abs) as usize; // Length of the dilated buffer (exclusive end)
    if dilated_span_len == 0 {
        return Ok((counter, None, None));
    }

    // Offsets of the original core within the dilated span. These values are measured relative
    // to `dilated_start_abs`, so they represent indices into the dilated buffers rather than
    // absolute genomic coordinates.
    let core_start_offset = (tile.core_start - dilated_start_abs) as usize;
    let core_end_offset_exclusive = (tile.core_end - dilated_start_abs) as usize;
    let dilated_start_i64 = dilated_start_abs as i64;
    let dilated_end_i64 = dilated_end_abs as i64;
    let dilated_start_abs_u64 = dilated_start_abs as u64; // Absolute coordinate of dilated buffer origin
    let core_start_abs = tile.core_start as u64;
    let core_end_abs = tile.core_end as u64;

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
                    let rel_start = (fragment.start - fetch_start) as u64;
                    let rel_end = (fragment.end - fetch_start) as u64;
                    corrector.correct_fragment(rel_start, rel_end, prefixes)
                }
                _ => Ok(None),
            }
        }
    };

    let correct_gc = opt.gc.gc_file.is_some();
    let fetch_start = tile.fetch_start;
    let fetch_end = tile.fetch_end;

    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;

        if fragment.start < fetch_start || fragment.end > fetch_end {
            // Fragment won't overlap the counting region (assuming correct max_fragment_length+window halo!)
            continue;
        }

        // Get GC correction weight
        let gc_weight = if opt.gc.gc_tag.is_some() {
            if fragment.gc_tag.had_invalid {
                counter.gc_failed_fragments += 1;
                if fragment.gc_tag.was_out_of_range {
                    counter.gc_out_of_range_tags += 1;
                }
                if opt.gc.drop_invalid_gc {
                    continue;
                } else {
                    1.0
                }
            } else if let Some(tag_w) = fragment.gc_tag.weight {
                tag_w as f64
            } else {
                counter.gc_failed_fragments += 1;
                if opt.gc.drop_invalid_gc {
                    continue;
                } else {
                    1.0
                }
            }
        } else {
            let gc_weight_opt = get_gc_weight(&fragment, fetch_start)?;
            match (gc_weight_opt, correct_gc) {
                (Some(w), true) => w,
                (None, true) => {
                    // Tried but failed to make a GC correction weight
                    counter.gc_failed_fragments += 1;
                    if opt.gc.drop_invalid_gc {
                        continue;
                    } else {
                        1.0
                    }
                }
                (None, false) => 1.0, // No correction
                (Some(_), false) => unreachable!(),
            }
        };

        let fragment_start = fragment.start as i64;
        let fragment_end = fragment.end as i64;

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
    let mask = build_mask_for_core(
        dilated_start_abs,
        dilated_end_abs,
        blacklist_chr,
        chrom_len,
        window_left,
        window_right,
    );
    let mask_slice = mask.as_deref();

    counter.add_from_snapshot(iter.counters_snapshot());

    if return_wps_instead {
        return Ok((counter, Some(wps_values), mask));
    }

    let mut prefix_all_cache: Option<Vec<f64>> = None;
    let mut prefix_allowed_cache: Option<Vec<f64>> = None;
    let mut prefix_allowed_cnt_cache: Option<Vec<u32>> = None;

    match mode {
        TileMode::Positional {
            windows,
            out_path,
            indexed,
        } => {
            let mut positional_writer = open_zstd_auto_writer(&out_path, 3, None)?;

            match windows {
                None => {
                    emit_bedgraph_runs(
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
                    for &(window_start_abs, window_end_abs, original_idx) in
                        overlapping_windows_for_tile(windows_for_chr, tile, tile_window_span)
                    {
                        let (local_start_idx, local_end_idx, _, _) = match clip_window_to_core(
                            window_start_abs,
                            window_end_abs,
                            tile.core_start,
                            tile.core_end,
                            dilated_start_abs,
                        ) {
                            Some(v) => v,
                            None => continue,
                        };

                        if indexed {
                            emit_windowed_runs(
                                &tile.chr,
                                &wps_values,
                                mask_slice,
                                local_start_idx,
                                local_end_idx,
                                dilated_start_abs_u64,
                                Some(original_idx),
                                decimals,
                                keep_zero_runs,
                                &mut positional_writer,
                            )?;
                        } else {
                            emit_windowed_runs(
                                &tile.chr,
                                &wps_values,
                                mask_slice,
                                local_start_idx,
                                local_end_idx,
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
                prefix_all_cache.as_ref().unwrap().as_slice()
            };
            let (ps_allowed_slice, ps_allowed_cnt_slice) = if masked_mode {
                let mask_slice_ref = mask_slice.expect("mask slice present");
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

            for &(window_start_abs, window_end_abs, original_idx) in
                overlapping_windows_for_tile(windows, tile, tile_window_span)
            {
                let clipped = clip_window_to_core(
                    window_start_abs,
                    window_end_abs,
                    tile.core_start,
                    tile.core_end,
                    dilated_start_abs,
                );
                let Some((local_start_idx, local_end_idx, _, _)) = clipped else {
                    continue;
                };

                let crosses_boundary =
                    !(window_start_abs >= core_start_abs && window_end_abs <= core_end_abs);

                let (sum, allowed, blacklisted) = wps_sum_and_counts(
                    local_start_idx,
                    local_end_idx,
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
        }

        TileMode::AggregatesBySize {
            window_bp,
            masked: _,
            finals_out,
            partials_out,
            cross_idx_out,
            guaranteed_aligned,
        } => {
            let action =
                per_window_wps_action.expect("per-window action required when using aggregates");
            let masked_mode = mask_slice.is_some();
            let ps_all_slice = {
                if prefix_all_cache.is_none() {
                    prefix_all_cache = Some(build_prefix(&wps_values));
                }
                prefix_all_cache.as_ref().unwrap().as_slice()
            };
            let (ps_allowed_slice, ps_allowed_cnt_slice) = if masked_mode {
                let mask_slice_ref = mask_slice.expect("mask slice present");
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

                    let (local_start_idx, local_end_idx, clipped_start, clipped_end) =
                        match clip_window_to_core(
                            bin_start,
                            bin_end,
                            tile.core_start,
                            tile.core_end,
                            dilated_start_abs,
                        ) {
                            Some(v) => v,
                            None => continue,
                        };

                    let (sum, allowed, blacklisted) = wps_sum_and_counts(
                        local_start_idx,
                        local_end_idx,
                        masked_mode,
                        ps_all_slice,
                        ps_allowed_slice,
                        ps_allowed_cnt_slice,
                        mask_slice,
                    );

                    let unmasked_span_bp = (clipped_end - clipped_start) as u64;
                    let value =
                        finalize_value(sum, allowed, unmasked_span_bp, masked_mode, &action);
                    let value = round_to(value, decimals);

                    write_final_row(
                        &mut finals_writer,
                        &tile.chr,
                        clipped_start,
                        clipped_end,
                        value,
                        blacklisted,
                        decimals,
                    )?;
                }

                finals_writer.flush()?;
            } else {
                let mut partials_writer = open_zstd_auto_writer(&partials_out, 3, None)?;
                let mut cross_index_writer = open_zstd_auto_writer(&cross_idx_out, 3, None)?;

                for bin_index in first_bin_index..=last_bin_index {
                    let bin_start = bin_index * window_bp;
                    let bin_end = (bin_index + 1) * window_bp;

                    let (local_start_idx, local_end_idx, _clipped_start, _clipped_end) =
                        match clip_window_to_core(
                            bin_start,
                            bin_end,
                            tile.core_start,
                            tile.core_end,
                            dilated_start_abs,
                        ) {
                            Some(v) => v,
                            None => continue,
                        };

                    let (sum, allowed, blacklisted) = wps_sum_and_counts(
                        local_start_idx,
                        local_end_idx,
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
            }
        }
    }

    return Ok((counter, None, None));
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
    dilated_start: u32,
    dilated_end: u32,
    blacklist_intervals: &[(u64, u64)],
    chromosome_length: u64,
    left_span: u32,
    right_span: u32,
) -> Option<Vec<u8>> {
    let dilated_span_len = (dilated_end - dilated_start) as usize;
    if dilated_span_len == 0 {
        return None;
    }

    let dilated_start_abs = dilated_start as u64;
    let dilated_end_abs = dilated_end as u64;
    let mut mask = vec![0u8; dilated_span_len];
    let mut has_masked_positions = false;

    for &(interval_start, interval_end) in blacklist_intervals {
        if interval_end <= dilated_start_abs || interval_start >= dilated_end_abs {
            continue;
        }
        let clipped_start = interval_start.max(dilated_start_abs);
        let clipped_end = interval_end.min(dilated_end_abs);
        if clipped_start >= clipped_end {
            continue;
        }

        // Treat blacklist intervals as half-open [start, end) so the end offset stays exclusive
        let start_offset = (clipped_start - dilated_start_abs) as usize;
        let end_offset_exclusive = (clipped_end - dilated_start_abs) as usize;
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

#[inline]
fn clip_window_to_core(
    abs_start: u64,
    abs_end: u64,
    core_start: u32,
    core_end: u32,
    dilated_start: u32,
) -> Option<(usize, usize, u64, u64)> {
    let core_start_u64 = core_start as u64;
    let core_end_u64 = core_end as u64;
    let start = abs_start.max(core_start_u64);
    let end = abs_end.min(core_end_u64);
    if start >= end {
        return None;
    }
    let local_start = (start as u32 - dilated_start) as usize;
    let local_end = (end as u32 - dilated_start) as usize;
    Some((local_start, local_end, start, end))
}
