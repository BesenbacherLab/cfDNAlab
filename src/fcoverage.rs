use anyhow::{Context, Result};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::io::Write;
use std::{collections::HashMap, fs::create_dir_all, path::PathBuf, sync::Arc, time::Instant};

use crate::utils::coverage::tiled_run::{
    build_tiles, merge_positional_tiles, reduce_aggregates_for_chr,
};
use crate::{
    cli_common::{ChromosomeArgs, FragmentLengthArgs, IOCArgs, WindowSpec, WindowsArgs},
    counters::FCoverageCounters,
    utils::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::load_blacklists,
        coverage::{
            coverage_prefix::CoveragePrefix,
            tiled_run::{Tile, TileMode, add_fragment_clipped_to_core, windows_overlapping_core},
            window_results::CoverageWindowAction,
            writers::NanPolicy,
        },
        fragment::segment_fragment::{SegmentedReadInfo, collect_fragment_with_segments},
    },
};

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
///  - Get the positional coverage for the included windows only. I.e.,
///    exclude all positions that do not overlap a window from the output.
///
/// Without windowing, positional coverage for the selected chromosomes are outputted.
///
/// ## Blacklisting
///
/// Positions in blacklisted regions are set to `f32::NaN` (and thus not included in sums or averages).
/// Set `--nan-policy` to change how these positions are handled in the output (positional coverage outputs only).
#[cfg_attr(feature = "cli", derive(clap::Args))]
pub struct FCoverageConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    ioc: IOCArgs,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory.
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
    ///     - "positions": Get the positional coverage for the included windows only. I.e.,
    ///                    exclude all positions that do not overlap a window from the output.
    ///
    /// NOTE: Ignored when no windows are specified.
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

    /// How to write coverage in blacklisted positions in position-coverage outputs `[string]`
    ///
    /// Possible values:
    ///
    ///     - "drop": Drop the row from the output (default).
    ///
    ///     - "nan": Write an literal NaN string.
    ///
    ///     - "empty": Leave the cell empty.
    ///
    /// NOTE: Ignored when no blacklist(s) are specified or the output is window-aggregates.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "drop",
            value_parser,
            ignore_case = true,
            help_heading = "Core"
        )
    )]
    pub nan_policy: NanPolicy,

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

    // Create output directory
    create_dir_all(&opt.ioc.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = if let Some(beds) = &opt.blacklist {
        println!("Start: Loading blacklists");
        load_blacklists(beds, 1, &chromosomes)?
    } else {
        HashMap::new()
    };

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            println!("Start: Loading window coordinates");
            Some(load_windows_from_bed(bed, &chromosomes, None)?)
        }
        _ => None,
    };

    let halo_bp: u32 = opt.fragment_lengths.max_fragment_length; // safe halo for pairing/segments

    let tiles = build_tiles(&opt.ioc.bam, &chromosomes, opt.tile_size, halo_bp)?;

    // Where per-tile files go
    let positional_prefix = "coverage.pos";
    let partials_prefix = "coverage.part";

    // Decide mode once
    let windowed = matches!(window_opt, WindowSpec::Bed(_));
    let masked = opt.blacklist.is_some();

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

    println!("Start: Counting per chromosome");

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

            // Decide tile mode and file name
            let action_prefix = if windowed {
                match opt.per_window {
                    CoverageWindowAction::OnlyIncludeThesePositions => positional_prefix,
                    CoverageWindowAction::Average | CoverageWindowAction::Total => partials_prefix,
                }
            } else {
                // Whole positional coverage
                positional_prefix
            };

            let out_path = opt.ioc.output_dir.join(format!(
                "{prefix}.{chr}.{idx}.tsv",
                prefix = action_prefix,
                chr = tile.chr,
                idx = tile.index
            ));

            let mode = if !windowed {
                TileMode::Positional {
                    windows: None,
                    out_path,
                }
            } else {
                match opt.per_window {
                    CoverageWindowAction::OnlyIncludeThesePositions => TileMode::Positional {
                        windows: windows_chr,
                        out_path,
                    },
                    CoverageWindowAction::Average | CoverageWindowAction::Total => {
                        let wchr = windows_chr.expect("windows required for aggregates");
                        TileMode::Aggregates {
                            windows: wchr,
                            masked,
                            out_path,
                        }
                    }
                }
            };

            let ctr = process_tile(&opt, tile, blacklist_chr, mode)?;
            pb.inc(1);
            Ok(ctr)
        })
        .collect::<anyhow::Result<_>>()?;

    pb.finish_with_message("| Finished counting");

    // Collect counters
    for counter in tile_results {
        global_counter += counter;
    }

    // Merge temporary output files and
    // reduce windows present in multiple tiles

    let final_out_path = if !windowed {
        // Whole-genome positional coverage
        merge_positional_tiles(
            &opt.ioc.output_dir,
            &chromosomes,
            positional_prefix,
            "coverage.per_positions.tsv",
        )?
    } else {
        match opt.per_window {
            CoverageWindowAction::OnlyIncludeThesePositions => {
                // Windowed positional with orig_idx column
                merge_positional_tiles(
                    &opt.ioc.output_dir,
                    &chromosomes,
                    positional_prefix,
                    "coverage.per_positions.tsv",
                )?
            }
            CoverageWindowAction::Average | CoverageWindowAction::Total => {
                // Per-chrom reduce of partials into final aggregates
                let final_path = opt.ioc.output_dir.join(match opt.per_window {
                    CoverageWindowAction::Average => "coverage.avg.tsv",
                    CoverageWindowAction::Total => "coverage.total.tsv",
                    _ => unreachable!(),
                });
                let mut w = std::io::BufWriter::new(std::fs::File::create(&final_path)?);
                // Header

                let value_col = match opt.per_window {
                    CoverageWindowAction::Average => "avg_coverage",
                    CoverageWindowAction::Total => "total_coverage",
                    _ => unreachable!(),
                };
                writeln!(
                    w,
                    "chromosome\tstart\tend\t{}\tblacklisted_positions",
                    value_col
                )?;

                if let Some(win_map) = &windows_map {
                    for chr in &chromosomes {
                        if let Some(wchr) = win_map.get(chr) {
                            reduce_aggregates_for_chr(
                                chr,
                                &opt.ioc.output_dir,
                                partials_prefix,
                                wchr.as_slice(),
                                masked,
                                &mut w,
                            )?;
                        }
                    }
                } else {
                    anyhow::bail!("Windows required for aggregates")
                }
                w.flush()?;
                final_path
            }
        }
    };

    println!("Saved output to: {:?}", final_out_path);

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
    println!(
        "Out-of-length-range-excluded fragments: {}",
        global_counter.illegal_length_fragments
    );
    // if opt.gc.bin_by_gc {
    //     println!("GC-excluded reads: {}", global_counter.gc_excl);
    // }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.counted_fragments
    );
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

/// Process one tile: pair reads, build coverage, and write outputs for this tile
fn process_tile(
    opt: &FCoverageConfig,
    tile: &Tile,
    blacklist_chr: &[(u64, u64)],
    mode: TileMode,
) -> Result<FCoverageCounters> {
    // Open a fresh BAM reader for this thread
    let (mut reader, _tid_check, _len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    debug_assert!(_tid_check == tile.tid as u32);

    reader
        .fetch((
            tile.tid as i32,
            tile.fetch_start as i64,
            tile.fetch_end as i64,
        ))
        .context(format!(
            "fetch {} {}-{}",
            &tile.chr, tile.fetch_start, tile.fetch_end
        ))?;

    // Prepare CP for tile core length
    let core_len = tile.core_end - tile.core_start;
    let mut cp = CoveragePrefix::initialize_coverage_prefix(core_len);

    // Mate-pair stash keyed by qname
    let mut stash: FxHashMap<Vec<u8>, SegmentedReadInfo> = FxHashMap::default();

    // Counters
    let mut counter = FCoverageCounters::default();

    // Iterate BAM records limited to fetch
    for res in reader.records() {
        let rec = res.context("reading bam record")?;
        counter.total_reads += 1;

        // Filter read using your existing policy
        if !include_read(&rec, opt) {
            continue;
        }

        match rec.is_reverse() {
            true => counter.accepted_reverse += 1,
            false => counter.accepted_forward += 1,
        }

        // Pair with mate if present in the stash
        if let Some(mate) = stash.remove(rec.qname()) {
            // Build a fragment with or without inter-mate gap
            let maybe_frag = collect_fragment_with_segments(
                &SegmentedReadInfo::from(&rec),
                &mate,
                1, // trigger_min_gap_bp (you can surface this in config)
                !opt.ignore_gap,
            );
            let Some(fragment) = maybe_frag else {
                continue;
            };

            counter.collected_fragments += 1;

            // Length filter on the final fragment span
            let fragment_length = fragment.len();
            if fragment_length < opt.fragment_lengths.min_fragment_length
                || fragment_length > opt.fragment_lengths.max_fragment_length
            {
                counter.illegal_length_fragments += 1;
                continue;
            }

            counter.counted_fragments += 1;

            // Clip and add to tile core coverage (segments respected)
            add_fragment_clipped_to_core(&mut cp, &fragment, 1.0, tile.core_start, tile.core_end)?;
        } else {
            // Stash for later pairing
            stash.insert(rec.qname().to_vec(), SegmentedReadInfo::from(&rec));
        }
    }

    // Finalize coverage
    cp.finalize_coverage();

    match mode {
        TileMode::Positional { windows, out_path } => {
            // We need a mask whenever there are blacklist intervals for this chromosome
            let need_mask = !blacklist_chr.is_empty();
            if need_mask && !blacklist_chr.is_empty() {
                // Clip and add blacklists late to minimize memory
                let mut clipped: Vec<(u64, u64)> = Vec::new();
                for &(bs, be) in blacklist_chr {
                    if be <= tile.core_start as u64 || bs >= tile.core_end as u64 {
                        continue;
                    }
                    let s = (bs as u32).max(tile.core_start) - tile.core_start;
                    let e = (be as u32).min(tile.core_end) - tile.core_start;
                    if s < e {
                        clipped.push((s as u64, e as u64));
                    }
                }
                if !clipped.is_empty() {
                    cp.initialize_blacklist_prefix();
                    cp.add_blacklist_many_to_prefix(&clipped)?;
                    cp.finalize_blacklist_prefix();
                }
            }

            let mut w = std::io::BufWriter::new(std::fs::File::create(out_path)?);
            let cov = cp.coverage().expect("coverage present");
            let mask = cp.blacklist_mask();

            // Write tile data to disk

            match windows {
                None => {
                    // Whole positional coverage for the tile core

                    for i in 0..(core_len as usize) {
                        let pos = tile.core_start as u64 + i as u64;

                        if let Some(m) = mask {
                            if m[i] == 1 {
                                // Blacklisted site -> act per policy
                                match opt.nan_policy {
                                    NanPolicy::DropRow => continue,
                                    NanPolicy::WriteLiteralNaN => {
                                        writeln!(w, "{}\t{}\t{}\tNaN", tile.chr, pos, pos + 1)?;
                                        continue;
                                    }
                                    NanPolicy::WriteEmptyCell => {
                                        writeln!(w, "{}\t{}\t{}\t", tile.chr, pos, pos + 1)?;
                                        continue;
                                    }
                                }
                            }
                        }

                        // Not masked
                        writeln!(w, "{}\t{}\t{}\t{}", tile.chr, pos, pos + 1, cov[i])?;
                    }
                }
                Some(win_chr) => {
                    // Only include windows that overlap the tile core
                    for &(window_start, window_end, original_idx) in
                        windows_overlapping_core(win_chr, tile.core_start, tile.core_end)
                    {
                        let s = (window_start as u32).max(tile.core_start);
                        let e = (window_end as u32).min(tile.core_end);
                        let a = (s - tile.core_start) as usize;
                        let b = (e - tile.core_start) as usize;

                        for i in a..b {
                            let pos = tile.core_start as u64 + i as u64;

                            if let Some(m) = mask {
                                if m[i] == 1 {
                                    match opt.nan_policy {
                                        NanPolicy::DropRow => continue,
                                        NanPolicy::WriteLiteralNaN => {
                                            writeln!(
                                                w,
                                                "{}\t{}\t{}\tNaN\t{}",
                                                tile.chr,
                                                pos,
                                                pos + 1,
                                                original_idx
                                            )?;
                                            continue;
                                        }
                                        NanPolicy::WriteEmptyCell => {
                                            writeln!(
                                                w,
                                                "{}\t{}\t{}\t\t{}",
                                                tile.chr,
                                                pos,
                                                pos + 1,
                                                original_idx
                                            )?;
                                            continue;
                                        }
                                    }
                                }
                            }

                            // If no NaNs, write coverage
                            writeln!(
                                w,
                                "{}\t{}\t{}\t{}\t{}",
                                tile.chr,
                                pos,
                                pos + 1,
                                cov[i],
                                original_idx
                            )?;
                        }
                    }
                }
            }

            w.flush()?;
        }

        TileMode::Aggregates {
            windows,
            masked,
            out_path,
        } => {
            // Build indexes once for this tile
            if masked && !blacklist_chr.is_empty() {
                // Same late-mask approach
                let mut clipped: Vec<(u64, u64)> = Vec::new();
                for &(bs, be) in blacklist_chr {
                    if be <= tile.core_start as u64 || bs >= tile.core_end as u64 {
                        continue;
                    }
                    let s = (bs as u32).max(tile.core_start) - tile.core_start;
                    let e = (be as u32).min(tile.core_end) - tile.core_start;
                    if s < e {
                        clipped.push((s as u64, e as u64));
                    }
                }
                if !clipped.is_empty() {
                    cp.initialize_blacklist_prefix();
                    cp.add_blacklist_many_to_prefix(&clipped)?;
                    cp.finalize_blacklist_prefix();
                }
            }
            cp.build_query_index()?;

            // Own everything we need so we don't hold borrows on `cp`
            let psum_all_owned = cp.get_psum_all().ok_or_else(|| {
                anyhow::anyhow!("psum_all missing; build_query_index() should have populated it")
            })?;
            let psum_allowed_owned = cp.get_psum_allowed();
            let cnt_allowed_ps_owned = cp.get_psum_allowed_count();
            let mask_owned: Option<Vec<u8>> = cp.blacklist_mask().map(|m| m.to_vec());

            // Use slices for indexing
            let psum_all: &[f64] = &psum_all_owned;
            let psum_allowed: Option<&[f64]> = psum_allowed_owned.as_deref();
            let cnt_allowed_ps: Option<&[u32]> = cnt_allowed_ps_owned.as_deref();
            let mask: Option<&[u8]> = mask_owned.as_deref();

            // Write per-tile partials: idx, sum, allowed_count, blacklisted_count
            let mut w = std::io::BufWriter::new(std::fs::File::create(out_path)?);

            for &(window_start, window_end, original_idx) in
                windows_overlapping_core(windows, tile.core_start, tile.core_end)
            {
                let s = (window_start as u32).max(tile.core_start);
                let e = (window_end as u32).min(tile.core_end);
                let a_us = (s - tile.core_start) as usize;
                let b_us = (e - tile.core_start) as usize;

                // Sum coverage via prefix sums (avoid calling cp.sum_coverage here)
                let sum = if masked {
                    if let Some(pa) = psum_allowed {
                        pa[b_us] - pa[a_us]
                    } else {
                        // No blacklist present -> allowed == all
                        psum_all[b_us] - psum_all[a_us]
                    }
                } else {
                    psum_all[b_us] - psum_all[a_us]
                };

                // Allowed positions count
                let allowed: u64 = if masked {
                    if let Some(cnt) = cnt_allowed_ps {
                        (cnt[b_us] - cnt[a_us]) as u64
                    } else if let Some(m) = mask {
                        let mut ok = 0u64;
                        for i in a_us..b_us {
                            if m[i] == 0 {
                                ok += 1;
                            }
                        }
                        ok
                    } else {
                        (b_us - a_us) as u64
                    }
                } else {
                    (b_us - a_us) as u64
                };

                let blacklisted = (b_us - a_us) as u64 - allowed;

                writeln!(w, "{}\t{}\t{}\t{}", original_idx, sum, allowed, blacklisted)?;
            }
            w.flush()?;
        }
    }

    Ok(counter)
}
