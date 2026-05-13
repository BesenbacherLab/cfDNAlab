use crate::shared::gc_tag::ClassifiedGCTagWeight;
use crate::shared::io::FinalOutputFiles;
use crate::{
    commands::{
        cli_common::{
            ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs, validate_output_prefix,
        },
        counters::ProfileGroupsCounters,
        gc_bias::{
            correct::{GCCorrector, load_gc_corrector},
            counting::build_gc_prefixes,
        },
        midpoints::{
            config::MidpointsConfig,
            counting_by_group::{ProfileGroupsCounts, SparseProfileGroupsCounts},
            group_index::{eligible_interval_counts_by_group, write_midpoint_group_index_tsv},
            postprocess::{ProfileLayout, postprocess_profile},
            settings::write_midpoint_profile_settings_json,
            windows::{ensure_uniform_window_len, prepare_count_windows},
        },
        run_statistics::{
            DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions, GCStatisticsSummary,
            TILE_DOUBLE_COUNT_NOTE, print_fragment_run_statistics,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_grouped_windows_from_bed,
        blacklist::is_blacklisted,
        fragment::minimal_fragment::Fragment,
        fragment_iterators::fragments_from_bam,
        interval::{IndexedInterval, Interval},
        io::dot_join,
        length_axis::LengthAxis,
        midpoint::midpoint_random_even_for_fragment,
        overlaps::find_overlapping_windows,
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq_in_range,
        scale_genome::{ScalingBin, compute_per_window_scaling_over_fragment},
        temp_chrom_names::TempChromNameMap,
        thread_pool::init_global_pool,
        tiled_run::{
            TempDirGuard, Tile, TileWindowSpan, build_tiles, clamp_fetch_to_window_span,
            overlapping_windows_for_tile, precompute_tile_window_spans,
        },
        window_fetch::window_derived_fetch_extent_for_core_overlap,
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{path::PathBuf, sync::Arc, time::Instant};
use tracing::{info, warn};

const COMMAND_TARGET: &str = "midpoints";
const DENSE_PROFILE_SIZE_WARNING_BYTES: usize = 50 * 1_000_000_000;

/// Execute the grouped midpoint profiling pipeline end-to-end.
///
/// The command produces dense midpoint profiles for grouped BED intervals. Internally, tile
/// workers count into sparse accumulators and write sparse `.npz` temporary files. After all tiles
/// finish, those sparse partial files are merged into one dense `ProfileGroupsCounts` and written
/// as the public `.midpoint_profiles.npy` output with axes `(group, length_bin, position)`.
///
/// Implementation details:
///
/// - Resolves chromosomes, loads grouped BED windows, and prepares optional blacklist and scaling
///   data before spawning parallel tiles.
/// - Streams fragments through per-tile accumulators, writing sparse temporary partial files that
///   are merged into a final 3D array and companion group index.
/// - Applies fragment length, blacklist, and scaling filters during aggregation so downstream tools
///   can consume ready-to-use profiles.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `midpoints` command.
///
/// Returns
/// -------
/// - `Ok(())`:
///     The output `npy` and group-index files were written successfully.
///
/// Errors
/// ------
/// - Returns an error if any input cannot be read, the grouped BED is invalid, or writing the
///   outputs fails.
pub fn run(opt: &MidpointsConfig) -> Result<()> {
    let start_time = Instant::now();
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    opt.gc.validate(opt.ref_2bit.as_deref())?;
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    let prefix = opt.output_prefix.trim();
    validate_output_prefix(prefix)?;

    // Create output directory
    ensure_output_dir(&opt.ioc.output_dir)?;

    // Load blacklist intervals if provided
    if opt.blacklist.is_some() {
        info!(target: COMMAND_TARGET, "Loading blacklists");
    }
    let blacklist_map = load_blacklist_map(
        opt.blacklist.as_ref(),
        opt.blacklist_min_size,
        0,
        &chromosomes,
    )?;

    // Load grouped fixed-size windows from BED
    info!(target: COMMAND_TARGET, "Loading fixed-size intervals");
    let (windows_map, group_idx_to_name) = load_grouped_windows_from_bed(
        opt.intervals.clone(),
        Some(chromosomes.as_slice()),
        None,
        None,
    )?;
    let num_groups = group_idx_to_name.len();
    #[cfg(feature = "plotters")]
    ensure!(
        opt.plot_groups.iter().all(|&idx| idx < num_groups),
        "--plot-groups: group index is out of bounds. There are {} groups (0-based).",
        num_groups
    );
    let total_windows: usize = windows_map.values().map(|gw| gw.len()).sum();
    info!(
        target: COMMAND_TARGET,
        "  Num. chromosomes: {:?} | Num. windows: {:?} | Num. groups: {:?}",
        windows_map.keys().len(),
        total_windows,
        num_groups,
    );

    // Ensure all windows have the same length
    let output_window_size = ensure_uniform_window_len(&windows_map)?;

    // Parse and validate fragment length bins once so all tiles share one lookup table
    let length_axis = Arc::new(LengthAxis::new(opt.resolve_length_bins()?)?);
    let min_fragment_length = length_axis.min_fragment_length();
    let max_fragment_length = length_axis.max_fragment_length();

    // Resolve counted and final profile dimensions
    let profile_layout = ProfileLayout::resolve(output_window_size, opt.bin_size, opt.smooth)?;

    // The grouped BED loader preserves group ids in IndexedInterval.idx. Moving the inner vectors
    // avoids cloning millions of intervals before tiling
    let indexed_windows_map: FxHashMap<String, Vec<IndexedInterval<u64>>> = windows_map
        .into_iter()
        .map(|(chromosome, grouped_windows)| (chromosome, grouped_windows.into_inner()))
        .collect();
    let interval_blacklist_margin = ((u64::from(max_fragment_length) + 1) / 2)
        .checked_add(u64::from(profile_layout.smoothing_flank))
        .context("interval blacklist margin overflow")?;
    let use_blacklist_prefilter = opt.blacklist.is_some() && !opt.keep_blacklisted_intervals;
    let (indexed_windows_map, interval_stats) = prepare_count_windows(
        indexed_windows_map,
        &contigs,
        &blacklist_map,
        profile_layout.smoothing_flank,
        interval_blacklist_margin,
        use_blacklist_prefilter,
    )?;
    ensure!(
        interval_stats.retained_for_counting > 0,
        "No midpoint intervals remain after filtering. Blacklist prefiltering dropped {} interval(s).",
        interval_stats.dropped_by_blacklist_prefilter
    );
    let group_eligible_interval_counts =
        eligible_interval_counts_by_group(&indexed_windows_map, &group_idx_to_name);

    info!(
        target: COMMAND_TARGET,
        "  Intervals after chromosome filtering: {} | blacklist-prefiltered: {} | retained: {}",
        interval_stats.loaded_after_chromosome_filtering,
        interval_stats.dropped_by_blacklist_prefilter,
        interval_stats.retained_for_counting,
    );

    let dense_profile_entries = num_groups
        .checked_mul(length_axis.num_bins())
        .and_then(|size| size.checked_mul(profile_layout.flanked_length))
        .context("dense midpoint profile shape overflow")?;
    let dense_profile_bytes = dense_profile_entries
        .checked_mul(std::mem::size_of::<f32>())
        .context("dense midpoint profile byte size overflow")?;
    let dense_profile_gb = dense_profile_bytes as f64 / 1_000_000_000.0;
    let final_profile_entries = num_groups
        .checked_mul(length_axis.num_bins())
        .and_then(|size| size.checked_mul(profile_layout.output_positions))
        .context("final midpoint profile shape overflow")?;
    let final_profile_bytes = final_profile_entries
        .checked_mul(std::mem::size_of::<f32>())
        .context("final midpoint profile byte size overflow")?;
    let final_profile_gb = final_profile_bytes as f64 / 1_000_000_000.0;
    info!(
        target: COMMAND_TARGET,
        "Dense midpoint counting shape: groups={} length_bins={} positions={} | approx {:.2} GB",
        num_groups,
        length_axis.num_bins(),
        profile_layout.flanked_length,
        dense_profile_gb,
    );
    if profile_layout.output_positions != profile_layout.flanked_length {
        info!(
            target: COMMAND_TARGET,
            "Final midpoint output shape: groups={} length_bins={} positions={} | approx {:.2} GB",
            num_groups,
            length_axis.num_bins(),
            profile_layout.output_positions,
            final_profile_gb,
        );
    }
    if final_profile_bytes >= DENSE_PROFILE_SIZE_WARNING_BYTES {
        warn!(
            target: COMMAND_TARGET,
            "Dense midpoint output is large: approx {:.2} GB for {} entries",
            final_profile_gb,
            final_profile_entries,
        );
    }

    // Load genomic scaling factors
    if opt.scale_genome.scaling_factors.is_some() {
        info!(target: COMMAND_TARGET, "Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<ScalingBin>> = load_scaling_map(
        &opt.scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(
            opt.gc.gc_file.is_some(),
            opt.gc.gc_tag.is_some(),
        ),
        None,
    )?;

    // Load GC correction package if specified
    if opt.gc.gc_file.is_some() {
        info!(target: COMMAND_TARGET, "Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.gc.gc_file.as_ref(),
        opt.ref_2bit.as_ref(),
        min_fragment_length,
        max_fragment_length,
    )?;

    // Build temporary directory
    let temp_dir_guard =
        TempDirGuard::new(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;
    let temp_dir = temp_dir_guard.path();

    // Build tiles with a pairing halo wide enough for any accepted fragment
    let halo_bp: u32 = max_fragment_length; // Safe halo for pairing
    let (tiles, _) = build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, None)?;
    let temp_chrom_name_map = TempChromNameMap::from_contigs(&chromosomes)?;
    let total_tiles = tiles.len();

    let windows_lookup = &indexed_windows_map;
    let tile_window_spans = Arc::new(precompute_tile_window_spans(
        &tiles,
        |chr| windows_lookup.get(chr).map(|w| w.as_slice()).unwrap_or(&[]),
        0,
        0,
    ));

    // Per-tile sparse partial files live in the temp directory
    let tmp_prefix = dot_join(&[prefix, "midpoint_profiles.tile"]);
    let tmp_prefix = tmp_prefix.as_str();

    // Create filenames for final output
    let final_counts_path = opt
        .ioc
        .output_dir
        .join(dot_join(&[prefix, "midpoint_profiles.npy"]));
    let map_path = opt
        .ioc
        .output_dir
        .join(dot_join(&[prefix, "group_index.tsv"]));
    let settings_path = opt
        .ioc
        .output_dir
        .join(dot_join(&[prefix, "midpoint_profile_settings.json"]));

    let progress = ProgressFactory::new();
    let pb = Arc::new(progress.default_bar(total_tiles as u64));

    // Configure global thread-pool size
    init_global_pool(opt.ioc.n_threads)?;

    info!(target: COMMAND_TARGET, "Counting per chromosome");

    pb.set_position(0);

    let tile_window_spans_for_threads = tile_window_spans.clone();
    let gc_tag = opt.gc.gc_tag.as_deref();

    let tile_results: Vec<(ProfileGroupsCounters, Option<PathBuf>)> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<(_, _)> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            // Borrow chromosome-local data for this tile worker
            let windows_chr: &[IndexedInterval<u64>] = indexed_windows_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let blacklist_chr: &[Interval<u64>] = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr: &[ScalingBin] = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            // Sparse tile partial file path. Empty tiles skip writing this file
            let chr_token = temp_chrom_name_map.token_for(tile.chr.as_str())?;
            let tile_counts_out = temp_dir.join(format!(
                "{prefix}.{chr}.{idx}.npz",
                prefix = tmp_prefix,
                chr = chr_token,
                idx = tile.index
            ));

            let out = process_tile(
                opt,
                tile,
                tile_counts_out,
                profile_layout.flanked_length,
                num_groups,
                Arc::clone(&length_axis),
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

    // Initialize global counter for accumulation across tiles
    let mut global_counter = ProfileGroupsCounters::default();

    // Collect temp paths and counters after Rayon returns to the main thread
    let mut all_tmp_count_paths: Vec<PathBuf> = Vec::with_capacity(total_tiles);
    for (counter, tmp_counts_path) in tile_results {
        if let Some(tmp_path) = tmp_counts_path {
            all_tmp_count_paths.push(tmp_path);
        }
        global_counter += counter;
    }

    info!(
        target: COMMAND_TARGET,
        "Merging temporary tile files to final output"
    );

    // Allocate the single final dense profile and merge sparse tile partial files into it
    let mut all_counts = ProfileGroupsCounts::new(
        profile_layout.flanked_length,
        num_groups,
        Arc::clone(&length_axis),
    );
    all_counts.add_from_sparse_npz_files_parallel(all_tmp_count_paths)?;
    let merged_counts_view = all_counts.view_ndarray3_group_len_pos();
    let transformed_counts = postprocess_profile(merged_counts_view, profile_layout)?;
    let final_counts_view = transformed_counts
        .as_ref()
        .map_or(merged_counts_view, |counts| counts.view());

    // Write every final output to the temp directory before moving any of them into place
    // This keeps failed writes from leaving a mix of old and new final files
    let mut final_outputs = FinalOutputFiles::new(temp_dir)?;

    let temp_counts_path = final_outputs.temp_path_for(&final_counts_path)?;
    info!(
        target: COMMAND_TARGET,
        "Writing final counts to temp file {}",
        temp_counts_path.display()
    );
    write_npy(&temp_counts_path, &final_counts_view)
        .with_context(|| format!("writing final counts to {}", temp_counts_path.display()))?;
    final_outputs.record(temp_counts_path, final_counts_path)?;

    let temp_map_path = final_outputs.temp_path_for(&map_path)?;
    info!(
        target: COMMAND_TARGET,
        "Writing group index to temp file {}",
        temp_map_path.display()
    );
    write_midpoint_group_index_tsv(
        &temp_map_path,
        &group_idx_to_name,
        &group_eligible_interval_counts,
    )?;
    final_outputs.record(temp_map_path, map_path)?;

    let temp_settings_path = final_outputs.temp_path_for(&settings_path)?;
    info!(
        target: COMMAND_TARGET,
        "Writing midpoint profile settings to temp file {}",
        temp_settings_path.display()
    );
    write_midpoint_profile_settings_json(
        &temp_settings_path,
        opt,
        &length_axis,
        profile_layout,
        interval_blacklist_margin,
        use_blacklist_prefilter,
    )?;
    final_outputs.record(temp_settings_path, settings_path)?;

    final_outputs.move_into_place()?;

    #[cfg(feature = "plotters")]
    {
        use crate::commands::midpoints::plotting::plot_midpoint_profiles;

        info!(
            target: COMMAND_TARGET,
            "Plotting selected groups' midpoint profiles"
        );

        plot_midpoint_profiles(
            prefix,
            &opt.ioc.output_dir,
            &opt.plot_groups,
            length_axis.edges(),
            &group_idx_to_name,
            final_counts_view,
        )?;
    }

    let elapsed = start_time.elapsed();
    let mut extra_statistics = vec![
        format!(
            "Intervals after chromosome filtering: {}",
            interval_stats.loaded_after_chromosome_filtering
        ),
        format!(
            "Blacklist-prefiltered intervals: {}",
            interval_stats.dropped_by_blacklist_prefilter
        ),
        format!(
            "Intervals retained for counting: {}",
            interval_stats.retained_for_counting
        ),
    ];
    if use_blacklist_prefilter {
        extra_statistics.push(format!(
            "Interval blacklist prefilter margin: {} bp",
            interval_blacklist_margin
        ));
    }

    print_fragment_run_statistics(
        &global_counter.base,
        elapsed,
        FragmentRunStatisticsOptions {
            include_section_header: true,
            notes: &[TILE_DOUBLE_COUNT_NOTE],
            labels: DEFAULT_FRAGMENT_STATISTICS_LABELS,
            blacklist_excluded_fragments: Some(global_counter.blacklisted_fragments),
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

/// Count midpoint contributions for one genomic tile.
///
/// This function is the per-tile worker body used by the parallel `run` loop. It opens a
/// chromosome-scoped BAM reader, narrows the fetch span to windows that can contribute to this
/// tile, streams fragments, and writes one sparse midpoint partial file when the tile has nonzero
/// counts.
///
/// The main counting flow is:
///
/// 1. Build per-tile helper state, including optional GC prefixes and scaling intervals.
/// 2. Narrow the BAM fetch to the extrema of windows overlapping the tile core.
/// 3. Stream paired or unpaired fragments from the BAM reader.
/// 4. Apply blacklist, midpoint-core, GC, and optional scaling filters.
/// 5. Convert each midpoint to a window-relative position and group index.
/// 6. Accumulate weighted counts in `SparseProfileGroupsCounts`.
/// 7. Write a sorted sparse `.npz` partial file, unless the tile had no counts.
///
/// Parameters
/// ----------
/// - `opt`:
///     Command configuration. The worker reads filters, IO paths, and correction options from it.
/// - `tile`:
///     Genomic tile whose core owns midpoint counts and whose fetch band includes pairing halo.
/// - `tile_counts_out`:
///     Destination path for the sparse tile partial file if the tile produces any counts.
/// - `window_size`:
///     Number of midpoint positions per grouped BED window.
/// - `num_groups`:
///     Number of group ids represented in the grouped BED input.
/// - `length_axis`:
///     Shared length-bin lookup used by the sparse accumulator.
/// - `windows`:
///     Chromosome-local grouped windows, sorted by coordinate.
/// - `tile_window_span`:
///     Optional precomputed slice range for windows that can overlap this tile.
/// - `blacklist_intervals`:
///     Chromosome-local blacklist intervals used before midpoint counting.
/// - `scaling_chr`:
///     Chromosome-local genomic scaling bins. Empty means no scaling is applied.
/// - `gc_corrector_opt`:
///     Optional file-based GC corrector for fragment-level weights.
/// - `gc_tag`:
///     Optional BAM tag name for tag-based GC weights.
///
/// Returns
/// -------
/// - `out`:
///     Tile-local run counters and an optional sparse temp-file path. The path is `None` when the
///     tile had no overlapping windows or no counted midpoint cells.
fn process_tile(
    opt: &MidpointsConfig,
    tile: &Tile,
    tile_counts_out: PathBuf,
    window_size: usize,
    num_groups: usize,
    length_axis: Arc<LengthAxis>,
    windows: &[IndexedInterval<u64>],
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_intervals: &[Interval<u64>],
    scaling_chr: &[ScalingBin],
    gc_corrector_opt: Option<GCCorrector>,
    gc_tag: Option<&str>,
) -> anyhow::Result<(ProfileGroupsCounters, Option<PathBuf>)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, _tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    tile.ensure_matches_bam_tid(_tid_check)?;

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
            // NOTE: Need the full fetch span to get GC of overlapping fragments!
            (tile.fetch_start() as usize)..(tile.fetch_end() as usize),
        )?;
        Some(build_gc_prefixes(&seq_bytes))
    } else {
        None
    };

    // The overlap finder only needs checked BED-like intervals here
    //
    // In BED mode, `find_overlapping_windows(...)` stores the scan position in the supplied slice
    // as `OverlappingWindow.idx`. It does not read `IndexedInterval.idx`. Because this temporary
    // list is built in the same order as `scaling_chr`, those scan positions already line up with
    // the chromosome-local indices used later to index back into `scaling_chr`
    //
    // So the carried `IndexedInterval.idx` value is intentionally a placeholder
    let scaling_with_bin_idx: Vec<IndexedInterval<u64>> = scaling_chr
        .iter()
        .map(|b| IndexedInterval::from_interval(b.interval, 0_u64))
        .collect();

    // Extract min/max fragment lengths once so fetch shrinking and fragment filtering share the
    // same bounds
    let min_fragment_length = length_axis.min_fragment_length();
    let max_fragment_length = length_axis.max_fragment_length();

    // Narrow the BAM fetch to windows that overlap the tile core. Tiles without contributing
    // windows return immediately and do not produce a sparse temp file
    let Some((core_overlapping_windows, fetch_span)) =
        get_overlapping_sites_and_adapt_fetch_to_extremes(
            windows,
            tile_window_span,
            tile,
            chrom_len as u32,
            max_fragment_length as u64,
        )?
    else {
        return Ok((counter, None));
    };
    let (fetch_from, fetch_to) = fetch_span.try_to_i64()?.as_tuple();

    reader
        .fetch((tile.tid, fetch_from, fetch_to))
        .context(format!("fetch {} {}-{}", &tile.chr, fetch_from, fetch_to))?;

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

    // Initialize sparse count map
    let mut counts = SparseProfileGroupsCounts::new(window_size, num_groups, length_axis);

    // Streaming pointers let sorted interval scans resume near the previous fragment
    let mut bl_ptr = 0; // Blacklist interval
    // `core_overlapping_windows` is compacted per tile, so this pointer is tile-local
    let mut wd_ptr = 0;
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
        // File-based GC correction shifts fragment coordinates into the fetched reference slice
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

    let correct_gc_from_file = opt.gc.gc_file.is_some();
    let fetch_start = tile.fetch_start();

    // Iterate fragments and count up midpoints
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;
        let fragment_length = fragment.len();

        // Determine fragment midpoint. Even-sized fragments use deterministic coordinate-derived
        // random rounding so tie positions are not always rounded in the same direction
        let midpoint =
            midpoint_random_even_for_fragment(&tile.chr, fragment.start(), fragment_length);

        // Keep only the fragments with midpoints within the tile
        if midpoint < tile.core_start() || midpoint >= tile.core_end() {
            continue;
        }

        // Determine blacklist status
        let in_blacklist = is_blacklisted(
            blacklist_intervals,
            opt.blacklist_strategy,
            fragment.interval.try_to_u64()?,
            max_fragment_length as u64,
            &mut bl_ptr,
        );
        if in_blacklist {
            counter.blacklisted_fragments += 1;
            continue;
        }

        // Get GC correction weight
        // NOTE: Must come after filtering for midpoints lying within the core!
        let gc_weight = if gc_tag.is_some() {
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
            // File-based correction path
            let gc_weight_opt = get_gc_weight(&fragment, fetch_start)?;
            match (gc_weight_opt, correct_gc_from_file) {
                (Some(w), true) => w,
                (None, true) => {
                    counter.gc_failed_fragments += 1;
                    if opt.gc.neutralize_invalid_gc {
                        1.0
                    } else {
                        continue;
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
            Interval::new(midpoint.into(), (midpoint + 1).into())?,
            0.99, // "Full" 1bp overlap but avoid roundoff error
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
                fragment.interval.try_to_u64()?, // Full fragment
                1. / (max_fragment_length as f64 + 1.0), // Any overlap
                max_fragment_length.into(),
            )
            .with_context(|| format!("finding overlapping scaling bins on chr {}", tile.chr))?
            .with_context(|| {
                format!(
                    "no overlapping scaling bins found for fragment {}:{}-{}. Scaling factors must cover every counted base on every counted chromosome",
                    tile.chr,
                    fragment.start(),
                    fragment.end()
                )
            })?;

            // Extract the indices of the overlapping bins
            let overlapping_scaling_bin_indices: Vec<usize> = overlapping_scaling_bins
                .windows
                .iter()
                .map(|w| w.idx)
                .collect();

            // Calculate the weight per overlapping count-window
            let overlap_weights = compute_per_window_scaling_over_fragment(
                fragment.interval.try_to_u64()?,
                &overlapping_windows,
                &overlapping_scaling_bin_indices,
                scaling_chr,
            )?;

            // Count up the weight per overlapping count-window
            for window_scaling in overlap_weights {
                let overlapped_window_idx = window_scaling.window_idx;
                let scaling_weight = window_scaling.scaling_weight;
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
            // When no scaling, increment counter by the GC weight for each window / bin
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

    // Add read and pairing counters captured by the fragment iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    // Empty sparse accumulators do not produce temp files. The final merge only receives paths for
    // tiles with at least one observed midpoint cell
    if counts.is_empty() {
        return Ok((counter, None));
    }

    // Write a sorted sparse partial file to the temporary directory
    counts
        .write_npz(&tile_counts_out)
        .context("Write midpoint tile counts fail")?;

    Ok((counter, Some(tile_counts_out)))
}

/// Collect midpoint sites that overlap a tile core and narrow the fetch span to their extremes.
///
/// Midpoints groups windows by `group_idx`, so the returned `IndexedInterval` values carry group
/// identifiers rather than the original BED row order used in some other commands. The helper
/// keeps only the sites overlapping the tile, computes the minimum start and maximum end across
/// those sites, and clamps the resulting fetch interval back onto the tile fetch band. When no
/// site overlaps the tile, the caller can skip the tile entirely.
///
/// Parameters
/// ----------
/// - `windows`:
///     Start-sorted midpoint windows for the current chromosome. Their `idx` field stores the
///     midpoint group index.
/// - `tile_span`:
///     Optional cached index range for windows that can overlap the tile.
/// - `tile`:
///     Tile whose core determines which sites are kept and whose fetch band bounds the result.
/// - `chrom_len`:
///     Chromosome length used to clamp the final fetch interval.
/// - `halo_bp`:
///     Extra bases to keep on both sides of the extreme overlapping sites before clamping back to
///     the tile fetch band. Callers pass the maximum fragment length here so fragment-overlapping
///     reads can still be reconstructed near tile and chromosome edges.
///
/// Returns
/// -------
/// - `out`:
///     `Some((sites, fetch_span))` when at least one midpoint site overlaps the tile and a
///     non-empty fetch interval remains after clamping. `None` when the tile has no overlapping
///     sites or the clamped fetch interval is empty.
pub fn get_overlapping_sites_and_adapt_fetch_to_extremes(
    windows: &[IndexedInterval<u64>],
    tile_span: Option<&TileWindowSpan>,
    tile: &Tile,
    chrom_len: u32,
    halo_bp: u64,
) -> Result<Option<(Vec<IndexedInterval<u64>>, Interval<u64>)>> {
    let reserve_hint = tile_span
        .map(|span| span.last_idx_exclusive.saturating_sub(span.first_idx))
        .unwrap_or(0)
        .min(windows.len());
    let mut overlapping_sites = Vec::with_capacity(reserve_hint);
    for site in overlapping_windows_for_tile(windows, tile, tile_span) {
        overlapping_sites.push(*site);
    }
    if overlapping_sites.is_empty() {
        return Ok(None);
    }

    let window_span = window_derived_fetch_extent_for_core_overlap(windows, tile, tile_span)?
        .context("midpoint helper found overlapping sites but no core-overlap window extent")?;

    let Some(fetch_span) =
        clamp_fetch_to_window_span(tile, chrom_len as u64, window_span, halo_bp)?
    else {
        return Ok(None);
    };

    Ok(Some((overlapping_sites, fetch_span)))
}
