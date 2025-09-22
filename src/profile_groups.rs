use anyhow::{Context, Result};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{fs::create_dir_all, path::PathBuf, sync::Arc, time::Instant};

use crate::{
    cli_common::{ChromosomeArgs, IOCArgs, ScaleGenomeArgs},
    counters::ProfileGroupsCounters,
    utils::{
        bam::{bam_contigs_info, create_chromosome_reader},
        blacklist::{BlacklistStrategy, is_blacklisted, load_blacklists},
        coverage::{
            scale_genome::{compute_window_scaling_over_fragment, load_scaling_factors_tsv},
            tiled_run::{Tile, build_tiles, make_temp_dir},
        },
        fragment::minimal_fragment::Fragment,
        fragment_iterator::fragments_from_bam,
        grouped_bed::{
            ensure_uniform_window_len, load_grouped_windows_from_bed, write_group_idx_to_name_tsv,
        },
        overlaps::find_overlapping_windows,
        profiling::{
            counting_by_group::ProfileGroupsCounts, midpoint::midpoint_random_even_with_thread_rng,
        },
        read::default_include_read,
    },
};

// Handle deletions?

/// Count fragment lengths in a BAM-file.
///
/// The fragment span is defined as `[end(reverse), start(forward)]`. // TODO: exclusive??
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct ProfileGroupsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    ioc: IOCArgs,

    /// Prefix for output files (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.midpoint_profiles.npy`.
    #[cfg_attr(
        feature = "cli",
        clap(long, short = 'x', default_value = "sites", help_heading = "Core")
    )]
    pub output_prefix: String,

    /// The grouped fixed-size intervals to count within `[path]`
    ///
    /// A BED file of genomic intervals and their respective group names.
    ///
    /// Must be sorted by the `chromosome` and `start` coordinates, and
    /// all intervals must have the same length.
    ///
    /// Sites with the same group name are collapsed to a single profile.
    ///
    /// Columns: `chromosome, start, end, group_name`.
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'w',
            long,
            value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub intervals: PathBuf,

    /// Edges of fragment length bins to count in `[path]`
    ///
    /// The last edge is *exclusive*.
    ///
    /// Example: `--length-bins 20 80 150 220 500 1000`.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            value_parser = clap::value_parser!(u32).range(1..),
            num_args = 2.., // at least two edges per occurrence
            default_values_t = [20_u32, 1001_u32],
            help_heading = "Core"
        )
    )]
    pub length_bins: Vec<u32>,

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

    #[cfg_attr(feature = "cli", clap(flatten))]
    chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    scale_genome: ScaleGenomeArgs,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads `[flag]`
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// **NOTE**: It may be an advantage to instead remove intervals that overlap
    /// blacklisted regions from the BED file.
    #[cfg_attr(
        feature = "cli", clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum size of blacklist intervals to load (bp) `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-min-size",
            default_value = "1",
            help_heading = "Filtering"
        )
    )]
    pub blacklist_min_size: u64,

    /// The fragment positions that should overlap blacklisted regions for it to be excluded `[string]`
    ///
    /// Possible values:
    ///     "any", "all", "midpoint", or "proportion=<threshold>"
    ///
    /// Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-strategy",
            default_value = "any",
            ignore_case = true,
            help_heading = "Filtering"
        )
    )]
    pub blacklist_strategy: BlacklistStrategy,
    // #[cfg_attr(feature = "cli", clap(flatten))]
    // gc: GCArgs,

    // #[cfg_attr(feature = "cli", clap(flatten))]
    // two_bit: TwoBitArgs,
}

pub fn run(opt: ProfileGroupsConfig) -> Result<()> {
    let start_time = Instant::now();
    let chromosomes = opt
        .chromosomes
        .resolve_chromosomes(Some(&opt.ioc.bam.as_path()))?;
    let prefix = opt.output_prefix.trim();
    let contigs = bam_contigs_info(&opt.ioc.bam, &chromosomes)?;

    // Create output directory
    create_dir_all(&opt.ioc.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = if let Some(beds) = &opt.blacklist {
        println!("Start: Loading blacklists");
        load_blacklists(beds, opt.blacklist_min_size, &chromosomes)?
    } else {
        FxHashMap::default()
    };

    // Load sites from BED file
    println!("Start: Loading fixed-size intervals");
    let (windows_map, group_idx_to_name) =
        load_grouped_windows_from_bed(opt.intervals.clone(), &chromosomes, None)?;
    let num_groups = group_idx_to_name.len();
    let total_windows: usize = windows_map.values().map(|gw| gw.len()).sum();
    println!(
        "       Num. chromosomes: {:?} | Num. windows: {:?} | Num. groups: {:?}",
        windows_map.keys().len(),
        total_windows,
        num_groups,
    );

    // Ensure all windows have the same length
    let window_size = ensure_uniform_window_len(&windows_map)?;

    // Load genomic scaling factors
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> =
        if let Some(path) = &opt.scale_genome.scaling_factors {
            println!("Start: Loading scaling factors");
            load_scaling_factors_tsv(path, &chromosomes, &contigs)?
        } else {
            FxHashMap::with_hasher(Default::default())
        };

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;

    // Prepare length bins
    let mut length_bins = opt.length_bins.clone();
    length_bins.sort_unstable();
    let num_length_bins = length_bins.len();
    let max_fragment_length = length_bins[num_length_bins - 1];

    // Build tiles
    let halo_bp: u32 = max_fragment_length; // Safe halo for pairing
    let (tiles, _) = build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, None)?;
    let total_tiles = tiles.len();

    // Where per-tile files go
    let tmp_prefix = format!("{prefix}.midpoint_profiles.tile");
    let tmp_prefix = tmp_prefix.as_str();

    // Create filenames for final output
    let final_counts_path = opt
        .ioc
        .output_dir
        .join(format!("{prefix}.midpoint_profiles.npy"));
    let map_path = opt
        .ioc
        .output_dir
        .join(format!("{}.group_index.tsv", prefix));

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

    // Prepare per-bin counts and metadata
    let mut global_counter = ProfileGroupsCounters::default();

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let tile_results: Vec<(ProfileGroupsCounters, Option<PathBuf>)> = tiles
        .par_iter()
        .map(|tile| -> Result<(_, _)> {
            // Per-chrom projections
            let windows_chr: &[(u64, u64, u64)] = windows_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let blacklist_chr: &[(u64, u64)] = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr: &[(u64, u64, f32)] = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            // Windowed tmp outputs for faster reducer later on
            let tile_counts_out = temp_dir.join(format!(
                "{prefix}.{chr}.{idx}.npy",
                prefix = tmp_prefix,
                chr = tile.chr,
                idx = tile.index
            ));

            let out = process_tile(
                &opt,
                tile,
                tile_counts_out,
                window_size,
                num_groups,
                &length_bins,
                windows_chr,
                blacklist_chr,
                scaling_chr,
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    let mut all_tmp_count_paths: Vec<PathBuf> = Vec::with_capacity(total_tiles);
    // Collect results (in chromosome order) back into the global vectors
    for (counter, tmp_counts_path) in tile_results {
        if let Some(tmp_path) = tmp_counts_path {
            all_tmp_count_paths.push(tmp_path);
        }
        global_counter += counter;
    }

    println!("Start: Merging temporary tile files to final output");

    // Initialize count array and load+fill with tmp counts
    let mut all_counts = ProfileGroupsCounts::new(window_size, num_groups, length_bins.to_vec());
    all_counts.add_from_npy_1d_files_parallel(all_tmp_count_paths)?;
    let all_counts_3d_arr = all_counts.view_ndarray3_group_len_pos();

    println!("Start: Writing final counts to: {:?}", &final_counts_path);
    // Write final counts to output_dir
    write_npy(&final_counts_path, &all_counts_3d_arr).context("Write final fail")?;

    println!("Start: Writing group index to: {:?}", &map_path);
    write_group_idx_to_name_tsv(map_path, &group_idx_to_name)?;

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
    println!(
        "  Blacklist-excluded fragments: {}",
        global_counter.blacklisted_fragments
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

fn process_tile(
    opt: &ProfileGroupsConfig,
    tile: &Tile,
    tile_counts_out: PathBuf,
    window_size: usize,
    num_groups: usize,
    length_bins: &[u32],
    windows: &[(u64, u64, u64)],
    blacklist_intervals: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
) -> anyhow::Result<(ProfileGroupsCounters, Option<PathBuf>)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, _tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    debug_assert!(_tid_check == tile.tid as u32);

    // Initialize counters (default -> 0s)
    let mut counter = ProfileGroupsCounters::default();

    // Replace scaling factor with unused index for overlap finder
    let scaling_with_bin_idx: Vec<(u64, u64, u64)> =
        scaling_chr.iter().map(|(s, e, _)| (*s, *e, 0u64)).collect();

    // Adapt the fetch coordinates to the present windows
    // When no windows are present, skip this tile
    let Some((core_overlapping_windows, fetch_from, fetch_to)) =
        get_overlapping_sites_and_adapt_fetch_to_extremes(windows, &tile, chrom_len as u32)
    else {
        return Ok((counter, None));
    };

    reader
        .fetch((tile.tid as i32, fetch_from as i64, fetch_to as i64))
        .context(format!("fetch {} {}-{}", &tile.chr, fetch_from, fetch_to))?;

    // Extract min/max fragment lengths
    let min_fragment_length = length_bins[0];
    let max_fragment_length = length_bins[length_bins.len() - 1] - 1;

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let min_len = min_fragment_length;
        let max_len = max_fragment_length;
        move |f: &Fragment| {
            let len = f.len();
            len >= min_len && len <= max_len
        }
    };

    // Wrap to use opt
    let include_read_fn = {
        let opt = (*opt).clone();
        move |r: &Record| default_include_read(r, opt.require_proper_pair, opt.min_mapq)
    };

    // Initialize count array
    let mut counts = ProfileGroupsCounts::new(window_size, num_groups, length_bins.to_vec());

    // Streaming pointers and single fetch for this chr
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = 0; // Genomic window
    let mut sf_ptr = 0; // Scaling factor bin

    // Create fragment iterator
    let mut iter = fragments_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        include_read_fn,
        fragment_filter,
    )
    .with_local_counters();

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;
        let fragment_length = fragment.len();

        // Determine blacklist status
        let in_blacklist = is_blacklisted(
            blacklist_intervals,
            opt.blacklist_strategy.clone(),
            fragment.start.into(),
            fragment.end.into(),
            max_fragment_length as u64,
            &mut bl_ptr,
        );
        if in_blacklist {
            counter.blacklisted_fragments += 1;
            continue;
        }

        // Determine fragment midpoint
        // Uses random rounding for even-sized fragments to avoid bias
        let midpoint = midpoint_random_even_with_thread_rng(fragment.start, fragment_length);

        // Only keep fragments with midpoints within the tile
        if midpoint < tile.core_start || midpoint >= tile.core_end {
            continue;
        }

        // Find all overlapping count-windows
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            Some(&core_overlapping_windows),
            None,
            midpoint.into(),
            (midpoint + 1).into(),
            0.99, // "Full" 1bp overlap but avoid rounding error
            max_fragment_length.into(),
        )?;
        let overlapping_windows = if let Some(overlaps) = overlapping_windows {
            overlaps
        } else {
            continue;
        };

        counter.counted_fragments += 1;

        // Find all overlapping scaling-factor bins
        // And count up the weight
        if !scaling_chr.is_empty() {
            // Find overlapping scaling-bins
            let overlapping_scaling_bins = find_overlapping_windows(
                chrom_len,
                &mut sf_ptr,
                Some(&scaling_with_bin_idx),
                None,
                fragment.start.into(), // Full fragment
                fragment.end.into(),
                1. / (max_fragment_length as f64 + 1.0), // Any overlap
                max_fragment_length.into(),
            )?
            .context("unwrapping overlapping scaling bins")?; // Should always find >= 1 bin

            // Extract the indices of the overlapping bins
            let overlapping_scaling_bin_indices: Vec<usize> = overlapping_scaling_bins
                .windows
                .iter()
                .map(|w| w.idx)
                .collect();

            // Calculate the weight per overlapping count-window
            let overlap_weights = compute_window_scaling_over_fragment(
                &overlapping_windows,
                &overlapping_scaling_bin_indices,
                scaling_chr,
            )?;

            // Count up the weight per overlapping count-window
            for (overlapped_window_idx, scaling_weight, _) in overlap_weights {
                let (window_start, _, group_idx) = core_overlapping_windows[overlapped_window_idx];
                let window_position = midpoint - window_start as u32;
                debug_assert!(
                    (window_start as u32) <= midpoint
                        && midpoint < (core_overlapping_windows[overlapped_window_idx].1 as u32),
                    "midpoint not inside window: midpoint={} window=({},{})",
                    midpoint,
                    window_start,
                    core_overlapping_windows[overlapped_window_idx].1
                );
                counts.incr_weighted(
                    window_position as usize,
                    group_idx as usize,
                    fragment_length as usize,
                    scaling_weight,
                )?;
            }
        } else {
            // When no scaling, increment counter by the overlap fraction for each window / bin
            for overlapped_window in overlapping_windows.windows {
                let overlapped_window_idx = overlapped_window.idx;
                let (window_start, _, group_idx) = core_overlapping_windows[overlapped_window_idx];
                let window_position = midpoint - window_start as u32;
                debug_assert!(
                    (window_start as u32) <= midpoint
                        && midpoint < (core_overlapping_windows[overlapped_window_idx].1 as u32),
                    "midpoint not inside window: midpoint={} window=({},{})",
                    midpoint,
                    window_start,
                    core_overlapping_windows[overlapped_window_idx].1
                );
                counts.incr_weighted(
                    window_position as usize,
                    group_idx as usize,
                    fragment_length as usize,
                    1.0,
                )?;
            }
        }
    }

    // Write tile counts to temp dir
    let arr1 = counts.as_ndarray1();
    write_npy(&tile_counts_out, &arr1).context("Write final fail")?;

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    Ok((counter, Some(tile_counts_out)))
}

pub fn get_overlapping_sites_and_adapt_fetch_to_extremes<'a>(
    windows: &'a [(u64, u64, u64)],
    tile: &Tile,
    chrom_len: u32,
) -> Option<(Vec<(u64, u64, u64)>, i64, i64)> {
    let cs = tile.core_start as u64;
    let ce = tile.core_end as u64;

    // Collect indices of overlaps, track extremes
    let mut idxs: Vec<usize> = Vec::new();
    let mut found = false;
    let mut min_ws: u64 = u64::MAX;
    let mut max_we: u64 = 0;

    // For fixed length, start-sorted windows, we can break early once ws >= ce
    for (i, &(ws, we, _)) in windows.iter().enumerate() {
        if ws >= ce {
            break; // Nothing further can overlap [cs, ce)
        }
        if we > cs && ws < ce {
            found = true;
            if ws < min_ws {
                min_ws = ws;
            }
            if we > max_we {
                max_we = we;
            }
            idxs.push(i);
        }
    }

    if !found {
        return None;
    }

    // Use the tile's actual halo (already clamped to chromosome edges in your tile)
    let left_halo = tile.core_start.saturating_sub(tile.fetch_start);
    let right_halo = tile.fetch_end.saturating_sub(tile.core_end);

    // Proposed narrower fetch band from window span ± halo
    let narrowed_start = (min_ws as u32).saturating_sub(left_halo);
    let narrowed_end = (max_we as u32).saturating_add(right_halo);

    // Intersect with tile's original fetch band and clamp to chrom length
    let start_u32 = narrowed_start.max(tile.fetch_start);
    let end_u32 = narrowed_end.min(tile.fetch_end).min(chrom_len);
    if start_u32 >= end_u32 {
        return None;
    }

    // Build vector of the overlappign sites from the indices
    let overlapping_sites = idxs.iter().map(move |i| windows[*i]).collect();

    Some((overlapping_sites, start_u32 as i64, end_u32 as i64))
}
