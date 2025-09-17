use anyhow::{Context, Result};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::io::Write;
use std::{fs::create_dir_all, path::PathBuf, sync::Arc, time::Instant};

use crate::cli_common::ScaleGenomeArgs;
use crate::utils::bam::bam_contigs_info;
use crate::utils::coverage::reducer::{
    reduce_aggregates_by_size_with_cross_index_for_chr, reduce_bed_with_cross_index_for_chr,
};
use crate::utils::coverage::scale_genome::{apply_scaling_in_place, load_scaling_factors_tsv};
use crate::utils::coverage::tiled_run::{
    adapt_fetch_to_extreme_windows, build_tiles, concat_aligned_size_tile_finals,
    coverage_sum_and_counts, emit_bedgraph_runs, emit_windowed_runs, finalize_value,
    intersect_abs_with_core_to_local, make_temp_dir, merge_positional_tiles, round_to,
};
use crate::utils::coverage::writer::{open_zstd_auto_writer, write_final_row};
use crate::utils::fragment::segment_fragment::FragmentWithSegments;
use crate::utils::fragment_iterator::fragments_with_segments_from_bam;
use crate::{
    cli_common::{ChromosomeArgs, FragmentLengthArgs, IOCArgs, WindowSpec, WindowsArgs},
    counters::FCoverageCounters,
    utils::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::load_blacklists,
        coverage::{
            coverage_prefix::Coverage,
            tiled_run::{Tile, TileMode, add_fragment_clipped_to_core, windows_overlapping_core},
            window_results::CoverageWindowAction,
        },
    },
};

// Support:
//  - bigwig

/// Count positional **fragment** coverage across the genome.
///
/// Only paired-end fragments with both reads present are counted. By default,
/// the entire fragment span `[start(forward), end(reverse))` is counted, except for
/// deletions and skipped regions that are not covered by the other read.
///
/// ## Windowing
///
/// When specifying windows (`--by-bed` or `--by-size`), one of the following outputs
/// is possible:
///
///  - Get the average coverage per window (default).
///
///  - Get the total coverage per window.
///
///  - Get the positional coverage for the included windows only (`--by-bed` *only*).
///    Excludes all positions that do not overlap a window from the output.
///    Choose between:
///     1) Indexed: Adds the original window index as an output column and keeps duplicate positions.
///     2) Unique: Overlapping windows are merged to reduce duplicate positions.
///         
///
/// Without windowing, positional coverage for the selected chromosomes are outputted.
///
/// ## Blacklisting
///
/// Positions in blacklisted regions are set to `f32::NaN` (and thus not included in sums or averages).
/// Set `--nan-policy` to change how these positions are handled in the output (positional coverage outputs only).
///
/// ## Temporary files
///
/// We write temporary files to a `<output-dir>/tmp.<output-prefix>.<random>` directory to reduce memory.
/// This directory is deleted at the end of the run. If the software is disrupted, the directory
/// may be left behind.
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct FCoverageConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    ioc: IOCArgs,

    /// Prefix for output files (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   <prefix>.per_position.bedgraph.zst
    ///   <prefix>.per_position_per_window.tsv.zst
    ///   <prefix>.avg.tsv.zst
    ///   <prefix>.total.tsv.zst
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value = "coverage", help_heading = "Core")
    )]
    pub output_prefix: String,

    /// Decimals to round coverage to when writing `[integer]`
    ///
    /// **NOTE**: When floating point precision is not needed,
    /// all coverages are integers, we remove all decimal points!
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "2", value_parser = clap::value_parser!(u8).range(0..), help_heading="Core"))]
    pub decimals: u8,

    /// Output zero-coverage runs in positional coverage outputs `[flag]`
    ///
    /// By default, only covered positions are written to the output.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub keep_zero_runs: bool,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    /// What to return per window `[string]`
    ///
    /// Possible values:
    ///
    ///     - "average": Get the average coverage per window (default).
    ///
    ///     - "total": Get the total coverage per window.
    ///
    ///     - "unique-positions": Get the positional coverage for the included windows only (`--by-bed` *only*).
    ///         Overlapping windows are merged to reduce duplicate positions.
    ///         Excludes all positions that do not overlap a window from the output.
    ///
    ///     - "indexed-positions": Get the positional coverage for the included windows only (`--by-bed` *only*).
    ///         Adds the original window index as an output column and keeps duplicate positions.
    ///         Excludes all positions that do not overlap a window from the output.
    ///
    /// **NOTE**: Ignored when no windows are specified.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "average",
            value_parser,
            ignore_case = true,
            help_heading = "Core"
        )
    )]
    pub per_window: CoverageWindowAction,

    /// Ignore inter-mate gap `[flag]`
    ///
    /// Disable counting of the gap between reads (i.e., `[forward.end, reverse.start)`)
    /// when the two reads do not overlap.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub ignore_gap: bool,

    #[cfg_attr(feature = "cli", clap(flatten))]
    windows: WindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    scale_genome: ScaleGenomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    fragment_lengths: FragmentLengthArgs,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads `[flag]`
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    // TODO: Consider whether blacklist is "filtering" in tools like this?
    /// Optional BED file(s) with blacklisted regions `[path]`
    #[cfg_attr(
        feature = "cli", clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,
    // #[cfg_attr(feature = "cli", clap(flatten))]
    // gc: GCArgs,

    // #[cfg_attr(feature = "cli", clap(flatten))]
    // two_bit: TwoBitArgs,
}

/// Whether to include the read or continue
fn include_read(rec: &Record, opt: &FCoverageConfig) -> bool {
    !(rec.is_unmapped()
        || rec.is_mate_unmapped()
        || rec.tid() != rec.mtid()
        || rec.is_secondary()
        || rec.is_supplementary()
        || rec.is_duplicate()
        || rec.is_quality_check_failed()
        || (opt.require_proper_pair && !rec.is_proper_pair())
        || rec.mapq() < opt.min_mapq) as bool
}

pub fn run(opt: FCoverageConfig) -> Result<()> {
    let start_time = Instant::now();
    let chromosomes = opt
        .chromosomes
        .resolve_chromosomes(Some(&opt.ioc.bam.as_path()))?;
    let window_opt = opt.windows.resolve_windows();
    let prefix = opt.output_prefix.trim();
    let contigs = bam_contigs_info(&opt.ioc.bam, &chromosomes)?;

    // Create output directory
    create_dir_all(&opt.ioc.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = if let Some(beds) = &opt.blacklist {
        println!("Start: Loading blacklists");
        load_blacklists(beds, 1, &chromosomes)?
    } else {
        FxHashMap::default()
    };

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            println!("Start: Loading window coordinates");
            let wds = load_windows_from_bed(bed, &chromosomes, None)?;
            if matches!(
                opt.per_window,
                CoverageWindowAction::OnlyIncludeThesePositionsUnique
            ) {
                // Merge in-place to avoid double memory-usage
                println!("Start: Merging overlapping/touching windows");
                // Take ownership so we can remove entries by chromosome
                let mut wds_owned: FxHashMap<String, crate::utils::bed::Windows> = wds;
                let mut out: FxHashMap<String, crate::utils::bed::Windows> =
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
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> =
        if let Some(path) = &opt.scale_genome.scaling_factors {
            println!("Start: Loading scaling factors");
            load_scaling_factors_tsv(path, &chromosomes, &contigs)?
        } else {
            FxHashMap::with_hasher(Default::default())
        };

    // Decide mode once
    let windowed = matches!(window_opt, WindowSpec::Bed(_) | WindowSpec::Size(_));
    let masked = opt.blacklist.is_some();
    let has_scaling = opt.scale_genome.scaling_factors.is_some();

    // Some actions cannot be used with `--by-size`
    if matches!(window_opt, WindowSpec::Size(_))
        && matches!(
            opt.per_window,
            CoverageWindowAction::OnlyIncludeThesePositionsUnique
                | CoverageWindowAction::OnlyIncludeThesePositionsIndexed
        )
    {
        anyhow::bail!("in --by-size mode, --per-window can only be 'average' or 'total'");
    }

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;

    // Window size when --by-size (otherwise None)
    let by_size_bp: Option<u64> = match &window_opt {
        WindowSpec::Size(bp) => Some(*bp as u64),
        _ => None,
    };

    // Build tiles
    let halo_bp: u32 = opt.fragment_lengths.max_fragment_length; // safe halo for pairing/segments
    let (tiles, tile_and_window_boundaries_align) =
        build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, by_size_bp)?;

    // Where per-tile files go
    let positional_prefix = format!("{prefix}.pos");
    let partials_prefix = format!("{prefix}.part");
    let finals_prefix = format!("{prefix}.fin");

    // Faster to convert to &str once
    let positional_prefix = positional_prefix.as_str();
    let partials_prefix = partials_prefix.as_str();
    let finals_prefix = finals_prefix.as_str();

    // Create filenames of final outputs
    let final_bedgraph_pos_name = format!("{prefix}.per_position.bedgraph.zst");
    let final_tsv_pos_name = format!("{prefix}.per_position_per_window.tsv.zst");
    let final_avg_name = format!("{prefix}.avg.tsv.zst");
    let final_total_name = format!("{prefix}.total.tsv.zst");

    // Get decimals to use
    let decimals_to_use: i32 = if windowed {
        match opt.per_window {
            CoverageWindowAction::Average | CoverageWindowAction::Total => opt.decimals as i32,
            CoverageWindowAction::OnlyIncludeThesePositionsUnique
            | CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                if has_scaling {
                    opt.decimals as i32
                } else {
                    0
                }
            }
        }
    } else {
        if has_scaling { opt.decimals as i32 } else { 0 }
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
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.ioc.n_threads as usize)
        .build_global()
        .context("building Rayon thread pool")?;

    let mut global_counter = FCoverageCounters::default();

    println!("Start: Counting per tile");

    pb.set_position(0);

    let tile_results: Vec<FCoverageCounters> = tiles
        .par_iter()
        .map(|tile| -> Result<FCoverageCounters> {
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
                match opt.per_window {
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
                match (&window_opt, opt.per_window) {
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
                            masked,
                            partials_out,
                            cross_idx_out,
                        }
                    }
                    (
                        WindowSpec::Size(size),
                        CoverageWindowAction::Average | CoverageWindowAction::Total,
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

            let ctr = process_tile(
                &opt,
                tile,
                blacklist_chr,
                scaling_chr,
                mode,
                decimals_to_use,
            )?;
            pb.inc(1);
            Ok(ctr)
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
            &opt.ioc.output_dir,
            &chromosomes,
            positional_prefix,
            final_bedgraph_pos_name.as_str(),
        )?
    } else {
        match opt.per_window {
            CoverageWindowAction::OnlyIncludeThesePositionsUnique => {
                // Windowed positional (unique and non-indexed)
                merge_positional_tiles(
                    &temp_dir,
                    &opt.ioc.output_dir,
                    &chromosomes,
                    positional_prefix,
                    final_bedgraph_pos_name.as_str(),
                )?
            }
            CoverageWindowAction::OnlyIncludeThesePositionsIndexed => {
                // Windowed positional with orig_idx column
                merge_positional_tiles(
                    &temp_dir,
                    &opt.ioc.output_dir,
                    &chromosomes,
                    positional_prefix,
                    final_tsv_pos_name.as_str(),
                )?
            }
            CoverageWindowAction::Average | CoverageWindowAction::Total => {
                // Per-chrom reduce of partials into final aggregates
                let final_path = opt.ioc.output_dir.join(match opt.per_window {
                    CoverageWindowAction::Average => final_avg_name.as_str(),
                    CoverageWindowAction::Total => final_total_name.as_str(),
                    _ => unreachable!(),
                });

                // Header value-column name
                let value_col = match opt.per_window {
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
                        let mut w =
                            open_zstd_auto_writer(&final_path, 3, Some(opt.ioc.n_threads as u32))?;

                        // Write header
                        writeln!(w, "{}", header)?;

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
                                    masked,
                                    opt.per_window,
                                    decimals_to_use,
                                    &mut w,
                                )?;
                            }
                        }
                        w.flush()?;
                    }
                    WindowSpec::Size(_) => {
                        if tile_and_window_boundaries_align {
                            let _ = concat_aligned_size_tile_finals(
                                &temp_dir,
                                &opt.ioc.output_dir,
                                &chromosomes,
                                finals_prefix,
                                match opt.per_window {
                                    CoverageWindowAction::Average => final_avg_name.as_str(),
                                    CoverageWindowAction::Total => final_total_name.as_str(),
                                    _ => unreachable!(),
                                },
                                &header,
                            )?;
                        } else {
                            let mut w = open_zstd_auto_writer(
                                &final_path,
                                3,
                                Some(opt.ioc.n_threads as u32),
                            )?;

                            // Write header
                            writeln!(w, "{}", header)?;

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
                                    masked,
                                    opt.per_window,
                                    chrom_len,
                                    decimals_to_use,
                                    &mut w,
                                )?;
                            }
                            w.flush()?;
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
    println!("  Total reads: {}", global_counter.total_reads);
    println!(
        "  Initially accepted reads: {} ({:.2}%, forward: {}, reverse: {})",
        global_counter.accepted_forward + global_counter.accepted_reverse,
        (global_counter.accepted_forward + global_counter.accepted_reverse) as f64
            / global_counter.total_reads as f64
            * 100.0,
        global_counter.accepted_forward,
        global_counter.accepted_reverse
    );
    // if opt.gc.bin_by_gc {
    //     println!("GC-excluded reads: {}", global_counter.gc_excl);
    // }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.counted_fragments
    );
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

/// Process one tile: pair reads, build coverage, and write outputs for this tile
fn process_tile(
    opt: &FCoverageConfig,
    tile: &Tile,
    blacklist_chr: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    mode: TileMode,
    decimals: i32,
) -> Result<FCoverageCounters> {
    // Open a fresh BAM reader for this thread
    let (mut reader, _tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    debug_assert!(_tid_check == tile.tid as u32);

    // Counters
    let mut counter = FCoverageCounters::default();

    // Adapt the fetch coordinates to the present windows (*in windowed mode!*)
    // When no windows are present, skip this tile
    let Some((fetch_from, fetch_to)) =
        adapt_fetch_to_extreme_windows(&tile, &mode, chrom_len as u32)
    else {
        return Ok(counter);
    };

    reader
        .fetch((tile.tid as i32, fetch_from as i64, fetch_to as i64))
        .context(format!("fetch {} {}-{}", &tile.chr, fetch_from, fetch_to))?;

    // Prepare CP for tile core length
    let core_len = tile.core_end - tile.core_start;
    let mut cp = Coverage::new(core_len);

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let min_len = opt.fragment_lengths.min_fragment_length;
        let max_len = opt.fragment_lengths.max_fragment_length;
        move |f: &FragmentWithSegments| {
            let len = f.len();
            len >= min_len && len <= max_len
        }
    };

    // Wrap to use opt
    let include_read_fn = {
        let opt = (*opt).clone();
        move |r: &Record| include_read(r, &opt)
    };

    let mut iter = fragments_with_segments_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        include_read_fn,
        1,
        !opt.ignore_gap,
        fragment_filter,
    )
    .with_local_counters();

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;
        counter.counted_fragments += 1;

        // Clip and add to tile core coverage (segments respected)
        add_fragment_clipped_to_core(&mut cp, &fragment, 1.0, tile.core_start, tile.core_end)?;
    }

    // Finalize coverage
    cp.finalize_coverage(true);

    // Apply per-bin scaling (in-place)
    if !scaling_chr.is_empty() {
        if let Some(cov_mut) = cp.coverage_mut() {
            apply_scaling_in_place(cov_mut, tile.core_start, scaling_chr);
        }
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
            add_clipped_blacklist_to_cp(&mut cp, &tile, !blacklist_chr.is_empty(), blacklist_chr)?;

            // Prepare compressed writer (zstd) for this tile
            let mut w = open_zstd_auto_writer(&out_path, 3, None)?;

            let cov = cp.coverage().expect("coverage present");
            let mask = cp.blacklist_mask();

            // Write tile data to disk

            match windows {
                None => {
                    // Whole positional coverage for the tile core
                    emit_bedgraph_runs(
                        &tile.chr,
                        cov,
                        mask,
                        0,
                        cov.len(),
                        tile.core_start as u64,
                        decimals,
                        opt.keep_zero_runs,
                        &mut w,
                    )?;
                }
                Some(win_chr) => {
                    // Only include windows that overlap the tile core
                    for &(window_start, window_end, original_idx) in
                        windows_overlapping_core(win_chr, tile.core_start, tile.core_end)
                    {
                        let (local_start_idx, local_end_idx) =
                            if let Some((local_start_idx, local_end_idx, _, _)) =
                                intersect_abs_with_core_to_local(
                                    window_start,
                                    window_end,
                                    tile.core_start,
                                    tile.core_end,
                                )
                            {
                                (local_start_idx, local_end_idx)
                            } else {
                                continue;
                            };

                        if indexed {
                            emit_windowed_runs(
                                &tile.chr,
                                cov,
                                mask,
                                local_start_idx,
                                local_end_idx,
                                tile.core_start as u64,
                                Some(original_idx),
                                decimals,
                                opt.keep_zero_runs,
                                &mut w,
                            )?;
                        } else {
                            emit_windowed_runs(
                                &tile.chr,
                                cov,
                                mask,
                                local_start_idx,
                                local_end_idx,
                                tile.core_start as u64,
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
            add_clipped_blacklist_to_cp(&mut cp, &tile, masked, blacklist_chr)?;
            cp.build_indexes(true)?;

            // Borrow indexes and mask once
            let psum_all = cp
                .psum_all_ref()
                .ok_or_else(|| anyhow::anyhow!("psum_all missing"))?;
            let psum_allowed = cp.psum_allowed_ref();
            let psum_cnt_allowed = cp.psum_allowed_count_ref();
            let mask: Option<&[u8]> = cp.blacklist_mask();

            // Writers: compressed partials and cross sidecar
            let mut w_part = open_zstd_auto_writer(&partials_out, 3, None)?;
            let mut w_cross = open_zstd_auto_writer(&cross_idx_out, 3, None)?;

            let cs = tile.core_start as u64;
            let ce = tile.core_end as u64;

            // Walk only windows overlapping the core (already start-sorted)
            for &(ws, we, idx) in windows_overlapping_core(windows, tile.core_start, tile.core_end)
            {
                let (local_start_idx, local_end_idx) =
                    if let Some((local_start_idx, local_end_idx, _, _)) =
                        intersect_abs_with_core_to_local(ws, we, tile.core_start, tile.core_end)
                    {
                        (local_start_idx, local_end_idx)
                    } else {
                        continue;
                    };

                // Classify as internal (fully inside core) vs boundary (crosses tile core boundary)
                let crosses_boundary = !(ws >= cs && we <= ce);

                // Sum coverage, respecting masked mode
                let (sum, allowed, blacklisted) = coverage_sum_and_counts(
                    local_start_idx,
                    local_end_idx,
                    masked,
                    psum_all,
                    psum_allowed,
                    psum_cnt_allowed,
                    mask,
                );

                // Always write a partial row; reducer will emit in orig_idx order
                // Internal windows won’t appear in the cross-index -> reducer expects 1 contribution
                // Boundary windows will appear in each crossed tile’s cross-index -> reducer expects N
                writeln!(w_part, "{}\t{}\t{}\t{}", idx, sum, allowed, blacklisted)?;
                if crosses_boundary {
                    // Cross-index lists the window’s orig_idx for the reducer
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
            add_clipped_blacklist_to_cp(&mut cp, &tile, masked, blacklist_chr)?;

            // Build prefix-sum indexes for fast per-window queries
            cp.build_indexes(true)?;

            // Own copies of the prefix arrays and optional mask to avoid long-lived borrows
            let psum_all = cp
                .psum_all_ref()
                .ok_or_else(|| anyhow::anyhow!("psum_all missing"))?;
            let psum_allowed = cp.psum_allowed_ref();
            let psum_cnt_allowed = cp.psum_allowed_count_ref();
            let mask: Option<&[u8]> = cp.blacklist_mask();

            // Determine the fixed-size windows that overlap the tile core
            let cs = tile.core_start as u64;
            let ce = tile.core_end as u64;
            let k_lo = cs / window_bp;
            let k_hi = (ce.saturating_sub(1)) / window_bp;

            if guaranteed_aligned {
                // FAST PATH: Every bin that touches the core is fully contained in this core
                // We compute the FINAL value here and write it once. No reducer later

                let mut w_fin = open_zstd_auto_writer(&finals_out, 3, None)?;

                for k in k_lo..=k_hi {
                    let bin_start = k * window_bp;
                    let bin_end = (k + 1) * window_bp;

                    // Intersect with core (alignment ensures this equals the bin for non-terminal tiles).
                    let (local_start_idx, local_end_idx, s_abs, e_abs) =
                        if let Some((local_start_idx, local_end_idx, s_abs, e_abs)) =
                            intersect_abs_with_core_to_local(
                                bin_start,
                                bin_end,
                                tile.core_start,
                                tile.core_end,
                            )
                        {
                            (local_start_idx, local_end_idx, s_abs, e_abs)
                        } else {
                            continue;
                        };

                    // Sum coverage, respecting masked mode
                    let (sum, allowed, blacklisted) = coverage_sum_and_counts(
                        local_start_idx,
                        local_end_idx,
                        masked,
                        psum_all,
                        psum_allowed,
                        psum_cnt_allowed,
                        mask,
                    );

                    // Compute final value now
                    let unmasked_span_bp = (e_abs - s_abs) as u64;
                    let value =
                        finalize_value(sum, allowed, unmasked_span_bp, masked, &opt.per_window);
                    let value = round_to(value, decimals);

                    // Emit FINAL row: chromosome  start  end  value  blacklisted_positions
                    write_final_row(
                        &mut w_fin,
                        &tile.chr,
                        bin_start,
                        bin_end.min(ce),
                        value,
                        blacklisted,
                        decimals,
                    )?;
                }

                w_fin.flush()?;
            } else {
                let mut w_part = open_zstd_auto_writer(&partials_out, 3, None)?;
                let mut w_cross = open_zstd_auto_writer(&cross_idx_out, 3, None)?;

                for k in k_lo..=k_hi {
                    let bin_start = k * window_bp;
                    let bin_end = (k + 1) * window_bp;

                    // Intersect with core (alignment ensures this equals the bin for non-terminal tiles).
                    let (local_start_idx, local_end_idx) =
                        if let Some((local_start_idx, local_end_idx, _, _)) =
                            intersect_abs_with_core_to_local(
                                bin_start,
                                bin_end,
                                tile.core_start,
                                tile.core_end,
                            )
                        {
                            (local_start_idx, local_end_idx)
                        } else {
                            continue;
                        };

                    // Sum coverage, respecting masked mode
                    let (sum, allowed, blacklisted) = coverage_sum_and_counts(
                        local_start_idx,
                        local_end_idx,
                        masked,
                        psum_all,
                        psum_allowed,
                        psum_cnt_allowed,
                        mask,
                    );

                    // PARTIAL row: start  end  sum  allowed  blacklisted
                    writeln!(
                        w_part,
                        "{}\t{}\t{}\t{}\t{}",
                        bin_start, bin_end, sum, allowed, blacklisted
                    )?;

                    // Mark cross-boundary bins (not fully inside the core) so reducer expects >1 contributions
                    let fully_inside = (bin_start >= cs) && (bin_end <= ce);
                    if !fully_inside {
                        writeln!(w_cross, "{}", bin_start)?;
                    }
                }
                w_part.flush()?;
                w_cross.flush()?;
            }
        }
    }

    Ok(counter)
}

/// Add blacklist clipped to tile coordinates to the coverage prefix object
fn add_clipped_blacklist_to_cp(
    cp: &mut Coverage,
    tile: &Tile,
    masked: bool,
    blacklist_chr: &[(u64, u64)],
) -> Result<()> {
    // Add blacklist late and clipped to the tile core to minimize memory
    // Use binary search to jump to the first overlapping interval
    if masked && !blacklist_chr.is_empty() {
        let cs = tile.core_start as u64;
        let ce = tile.core_end as u64;

        // Find first interval with end > core_start
        let mut i = blacklist_chr.partition_point(|&(_bs, be)| be <= cs);

        if i < blacklist_chr.len() {
            let mut clipped: Vec<(u64, u64)> = Vec::new();

            // Walk only the intervals that can overlap [cs, ce)
            while i < blacklist_chr.len() {
                let (bs, be) = blacklist_chr[i];
                if bs >= ce {
                    break; // Remaining intervals start after the core band
                }

                let s_abs = bs.max(cs);
                let e_abs = be.min(ce);
                if s_abs < e_abs {
                    // Convert to tile‐local coordinates
                    let s = (s_abs as u32) - tile.core_start;
                    let e = (e_abs as u32) - tile.core_start;
                    clipped.push((s as u64, e as u64));
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
