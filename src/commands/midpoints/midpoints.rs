use crate::{
    commands::{
        cli_common::{
            ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        counters::ProfileGroupsCounters,
        gc_bias::{
            correct::{GCCorrector, load_gc_corrector},
            counting::build_gc_prefixes,
        },
        midpoints::{
            config::MidpointsConfig, counting_by_group::ProfileGroupsCounts,
            windows::ensure_uniform_window_len,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::{load_grouped_windows_from_bed, write_group_idx_to_name_tsv},
        blacklist::is_blacklisted,
        fragment::minimal_fragment::Fragment,
        fragment_iterator::fragments_from_bam,
        interval::IndexedInterval,
        midpoint::midpoint_random_even_with_thread_rng,
        overlaps::find_overlapping_windows,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq_in_range,
        scale_genome::compute_window_scaling_over_fragment,
        thread_pool::init_global_pool,
        tiled_run::{
            Tile, TileWindowSpan, build_tiles, clamp_fetch_to_window_span, make_temp_dir,
            overlapping_windows_for_tile, precompute_tile_window_spans,
        },
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{path::PathBuf, sync::Arc, time::Instant};

// Handle deletions?

/// Execute the grouped midpoint profiling pipeline end-to-end.
///
/// Implementation details:
/// - Resolves chromosomes, loads grouped BED windows, and prepares optional blacklist and scaling
///   data before spawning parallel tiles.
/// - Streams fragments through per-tile accumulators, writing temporary `.npy` slices that are
///   merged into a final 3D array and companion group index.
/// - Applies fragment length, blacklist, and scaling filters during aggregation so downstream tools
///   can consume ready-to-use profiles.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `profile-groups` command.
///
/// Returns:
/// - `Ok(())` when the output `npy` and group-index files are written successfully.
///
/// Errors:
/// - Returns an error if any input cannot be read, the grouped BED is invalid, or writing the
///   outputs fails.
pub fn run(opt: &MidpointsConfig) -> Result<()> {
    let start_time = Instant::now();
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    let prefix = opt.output_prefix.trim();

    // Create output directory
    ensure_output_dir(&opt.ioc.output_dir)?;

    // Load blacklist intervals if provided
    if opt.blacklist.is_some() {
        println!("Start: Loading blacklists");
    }
    let blacklist_map = load_blacklist_map(
        opt.blacklist.as_ref(),
        opt.blacklist_min_size,
        0,
        &chromosomes,
    )?;

    // Load sites from BED file
    println!("Start: Loading fixed-size intervals");
    let (windows_map, group_idx_to_name) = load_grouped_windows_from_bed(
        opt.intervals.clone(),
        Some(chromosomes.as_slice()),
        None,
        None,
    )?;
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
    let mut indexed_windows_map: FxHashMap<String, Vec<IndexedInterval<u64>>> = FxHashMap::default();
    for (chromosome, grouped_windows) in &windows_map {
        let indexed_windows = grouped_windows
            .as_slice()
            .iter()
            .map(|&(start, end, group_idx)| IndexedInterval::new(start, end, group_idx))
            .collect::<crate::Result<Vec<_>>>()?;
        indexed_windows_map.insert(chromosome.clone(), indexed_windows);
    }

    // Parse and validate fragment-length bins once so all tiles use the same edges
    let length_bins = opt.resolve_length_bins()?;
    let num_length_bins = length_bins.len();
    let min_fragment_length = length_bins[0];
    let max_fragment_length = length_bins[num_length_bins - 1] - 1;

    // Load genomic scaling factors
    if opt.scale_genome.scaling_factors.is_some() {
        println!("Start: Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> =
        load_scaling_map(&opt.scale_genome, &chromosomes, &contigs)?;

    // Load GC correction package if specified
    if opt.gc.gc_file.is_some() {
        println!("Start: Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.gc.gc_file.as_ref(),
        min_fragment_length,
        max_fragment_length,
    )?;

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;

    // Build tiles
    let halo_bp: u32 = max_fragment_length; // Safe halo for pairing
    let (tiles, _) = build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, None)?;
    let total_tiles = tiles.len();

    let windows_lookup = &indexed_windows_map;
    let tile_window_spans = Arc::new(precompute_tile_window_spans(
        &tiles,
        |chr| windows_lookup.get(chr).map(|w| w.as_slice()).unwrap_or(&[]),
        0,
        0,
    ));

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
            .expect("hardcoded progress template"),
    );

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    // Prepare per-bin counts and metadata
    let mut global_counter = ProfileGroupsCounters::default();

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let tile_window_spans_for_threads = tile_window_spans.clone();
    let gc_tag = opt.gc.gc_tag.as_deref();

    let tile_results: Vec<(ProfileGroupsCounters, Option<PathBuf>)> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<(_, _)> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            // Per-chrom projections
            let windows_chr: &[IndexedInterval<u64>] = indexed_windows_map
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
                opt,
                tile,
                tile_counts_out,
                window_size,
                num_groups,
                &length_bins,
                windows_chr,
                tile_span.as_ref(),
                blacklist_chr,
                scaling_chr,
                gc_corrector.clone(),
                gc_tag,
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

    #[cfg(feature = "plotters")]
    {
        use crate::commands::midpoints::plotting::plot_midpoint_profiles;

        println!("Start: Plotting selected groups' midpoint profiles");

        plot_midpoint_profiles(
            prefix,
            &opt.ioc.output_dir,
            &opt.plot_groups,
            &length_bins,
            &group_idx_to_name,
            all_counts_3d_arr,
        )?;
    }

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

    println!();
    println!("Statistics");
    println!("----------");
    println!(
        "  Note: A few reads/fragments may be counted twice in the statistics (only) around the parallelization tiles."
    );

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
    println!(
        "  Blacklist-excluded fragments: {}",
        global_counter.blacklisted_fragments
    );
    if opt.gc.gc_file.is_some() || opt.gc.gc_tag.is_some() {
        let gc_fail_action = if opt.gc.drop_invalid_gc {
            "fragment skipped"
        } else {
            "fragment counted with weight 1.0"
        };
        println!(
            "  GC correction failures ({}): {}",
            gc_fail_action, global_counter.gc_failed_fragments
        );
        if opt.gc.gc_tag.is_some() && global_counter.gc_out_of_range_tags > 0 {
            println!(
                "  GC tag values outside [0, {:.0}] treated as invalid: {}",
                crate::shared::gc_tag::MAX_REASONABLE_GC_WEIGHT,
                global_counter.gc_out_of_range_tags
            );
        }
    }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.base.counted_fragments
    );
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn process_tile(
    opt: &MidpointsConfig,
    tile: &Tile,
    tile_counts_out: PathBuf,
    window_size: usize,
    num_groups: usize,
    length_bins: &[u32],
    windows: &[IndexedInterval<u64>],
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_intervals: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<GCCorrector>,
    gc_tag: Option<&str>,
) -> anyhow::Result<(ProfileGroupsCounters, Option<PathBuf>)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, _tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    debug_assert!(_tid_check == tile.tid as u32);

    // Initialize counters (default -> 0s)
    let mut counter = ProfileGroupsCounters::default();

    let gc_prefixes_opt = if gc_corrector_opt.is_some() {
        let ref_2bit = match opt.ref_2bit.as_ref() {
            Some(r) => r,
            None => bail!("When GC correction is specified, --ref-2bit must also be specified"),
        };
        let seq_bytes = read_seq_in_range(
            ref_2bit,
            &tile.chr,
            // NOTE: Need for full fetch span to get GC of overlapping fragments!
            (tile.fetch_start as usize)..(tile.fetch_end as usize),
        )?;
        Some(build_gc_prefixes(&seq_bytes))
    } else {
        None
    };

    // Replace scaling factor with unused index for overlap finder
    let scaling_with_bin_idx: Vec<IndexedInterval<u64>> = scaling_chr
        .iter()
        .map(|(start, end, _)| IndexedInterval::new(*start, *end, 0_u64))
        .collect::<crate::Result<_>>()?;

    // Adapt the fetch coordinates to the present windows
    // When no windows are present, skip this tile
    let Some((core_overlapping_windows, fetch_from, fetch_to)) =
        get_overlapping_sites_and_adapt_fetch_to_extremes(
            windows,
            tile_window_span,
            tile,
            chrom_len as u32,
        )
    else {
        return Ok((counter, None));
    };

    reader
        .fetch((tile.tid, fetch_from, fetch_to))
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

    // Initialize count array
    let mut counts = ProfileGroupsCounts::new(window_size, num_groups, length_bins.to_vec());

    // Streaming pointers
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = tile_window_span
        .and_then(|span| (!span.is_empty()).then_some(span.first_idx))
        .unwrap_or(0);
    let mut sf_ptr = 0; // Scaling factor bin

    // Create fragment iterator
    let gc_tag_bytes = gc_tag.map(|t| t.as_bytes().to_vec());
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

    let correct_gc_from_file = opt.gc.gc_file.is_some();
    let fetch_start = tile.fetch_start;

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;
        let fragment_length = fragment.len();

        // Determine blacklist status
        let in_blacklist = is_blacklisted(
            blacklist_intervals,
            opt.blacklist_strategy,
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

        // Get GC correction weight
        // NOTE: Must come after filtering for midpoints lying within the core!
        let gc_weight = if gc_tag.is_some() {
            // Tag mode: trust tag if valid, otherwise treat as failure
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
            // File-based correction path
            let gc_weight_opt = get_gc_weight(&fragment, fetch_start)?;
            match (gc_weight_opt, correct_gc_from_file) {
                (Some(w), true) => w,
                (None, true) => {
                    counter.gc_failed_fragments += 1;
                    if opt.gc.drop_invalid_gc {
                        continue;
                    } else {
                        1.0
                    }
                }
                (None, false) => 1.0, // No correction
                (Some(_), false) => bail!("unexpected GC weight when GC correction is disabled"),
            }
        };

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

        counter.base.counted_fragments += 1;

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
                let window = core_overlapping_windows[overlapped_window_idx];
                let window_start = window.start();
                let window_end = window.end();
                let group_idx = window.idx();
                let midpoint_u64 = midpoint as u64;
                ensure!(
                    window_start <= midpoint_u64 && midpoint_u64 < window_end,
                    "midpoint not inside window: midpoint={} window=({},{})",
                    midpoint,
                    window_start,
                    window_end
                );
                let window_position = usize::try_from(midpoint_u64 - window_start)
                    .context("window position does not fit in usize")?;
                counts.incr_weighted(
                    window_position,
                    group_idx as usize,
                    fragment_length as usize,
                    scaling_weight * gc_weight,
                )?;
            }
        } else {
            // When no scaling, increment counter by the overlap fraction for each window / bin
            for overlapped_window in overlapping_windows.windows {
                let overlapped_window_idx = overlapped_window.idx;
                let window = core_overlapping_windows[overlapped_window_idx];
                let window_start = window.start();
                let window_end = window.end();
                let group_idx = window.idx();
                let midpoint_u64 = midpoint as u64;
                ensure!(
                    window_start <= midpoint_u64 && midpoint_u64 < window_end,
                    "midpoint not inside window: midpoint={} window=({},{})",
                    midpoint,
                    window_start,
                    window_end
                );
                let window_position = usize::try_from(midpoint_u64 - window_start)
                    .context("window position does not fit in usize")?;
                counts.incr_weighted(
                    window_position,
                    group_idx as usize,
                    fragment_length as usize,
                    gc_weight,
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

// Windows: start, end, **`group_idx`** (not original_idx as in other commands).
pub fn get_overlapping_sites_and_adapt_fetch_to_extremes(
    windows: &[IndexedInterval<u64>],
    tile_span: Option<&TileWindowSpan>,
    tile: &Tile,
    chrom_len: u32,
) -> Option<(Vec<IndexedInterval<u64>>, i64, i64)> {
    let overlapping_sites: Vec<IndexedInterval<u64>> =
        overlapping_windows_for_tile(windows, tile, tile_span)
            .copied()
            .collect();

    if overlapping_sites.is_empty() {
        return None;
    }

    let min_ws = overlapping_sites.iter().map(|window| window.start()).min().unwrap();
    let max_we = overlapping_sites.iter().map(|window| window.end()).max().unwrap();

    let (fetch_from, fetch_to) =
        clamp_fetch_to_window_span(tile, chrom_len as u64, min_ws, max_we)?;

    Some((overlapping_sites, fetch_from, fetch_to))
}
