use crate::{
    commands::{
        cli_common::{
            WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        counters::FragmentKmersCounters,
        fragment_kmers::{
            config::*,
            nearest_frame_guard::NearestFrameGuard,
            parse::*,
            positional_output::*,
            positions::*,
            selection::{SelectionDecision, evaluate_selection},
            tiling::*,
            windows::*,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::{apply_blacklist_mask_to_seq, apply_mask::BLACKLIST_BYTE, is_blacklisted},
        fragment::segment_kmer_fragment::FragmentWithKmerSegments,
        fragment_iterator::fragments_with_kmer_segments_from_bam,
        io::create_text_writer,
        kmers::{
            kmer_codec::{KmerCodes, KmerSpec, build_kmer_specs, build_left_aligned_codes_per_k},
            process_counts::{DecodedCounts, prepare_decoded_counts, split_and_decode_counts},
            write::write_decoded_counts_matrix,
        },
        overlaps::find_overlapping_windows,
        read::default_include_read,
        reference::read_seq_in_range,
        scale_genome::apply_scaling_to_coverage_in_place,
        thread_pool::init_global_pool,
        tiled_run::{
            Tile, TileWindowSpan, build_tiles, make_temp_dir, precompute_tile_window_spans,
        },
    },
};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{convert::TryInto, io::Write, path::Path, sync::Arc, time::Instant};

/// Execute the fragment kmers counting pipeline end-to-end.
///
/// Implementation details:
/// - Resolves chromosomes, prepares optional windows/blacklists/scaling data, and then processes
///   each chromosome in parallel tiles using Rayon.
/// - Streams fragments through per-window accumulators, enumerating the requested k-mers inside
///   every counted window and writing dense (or optional sparse) count matrices plus motif lists.
/// - Applies fragment-length, blacklist, indel, scaling, and strand handling policies consistently
///   across threads.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `fragment-kmers` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &FragmentKmersConfig) -> Result<()> {
    let start_time = Instant::now();
    let global_counter = run_inner(opt)?;

    if !opt.shared_args.quiet {
        println!();
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
        println!(
            "  Blacklist-excluded fragments: {}",
            global_counter.blacklisted_fragments
        );
        // if opt.gc.bin_by_gc {
        //     println!("GC-excluded reads: {}", global_counter.base.gc_excl);
        // }
        println!(
            "  Fragments counted one or more times: {}",
            global_counter.base.counted_fragments
        );
        println!("----------");
        println!("Elapsed time: {:.2?}", elapsed);
    }
    Ok(())
}

pub fn run_inner(opt: &FragmentKmersConfig) -> Result<FragmentKmersCounters> {
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.shared_args.chromosomes, &opt.shared_args.ioc)?;
    let window_opt = opt.shared_args.windows.resolve_windows();
    let position_specs = opt
        .shared_args
        .position_selection
        .clone()
        .into_positional_specs()?;
    let prefix = opt.shared_args.output_prefix.trim();
    let quiet = opt.shared_args.quiet;

    // Create output directory
    ensure_output_dir(&opt.shared_args.ioc.output_dir)?;

    // Load blacklist intervals if provided
    if opt.shared_args.blacklist.is_some() && !quiet {
        println!("Start: Loading blacklists");
    }
    let blacklist_map = load_blacklist_map(
        opt.shared_args.blacklist.as_ref(),
        opt.shared_args.blacklist_min_size,
        0,
        &chromosomes,
    )?;

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            if !quiet {
                println!("Start: Loading window coordinates");
            }
            Some(load_windows_from_bed(
                bed,
                Some(chromosomes.as_slice()),
                None,
            )?)
        }
        _ => None,
    };

    let kmer_specs: FxHashMap<u8, KmerSpec> = build_kmer_specs(&opt.kmer_sizes)?;

    let positional_cache = {
        if opt.shared_args.base_selection.bases_from != BasesFrom::Reference {
            bail!("position selection currently supports bases-from=reference only");
        }
        // Parse each positions specification
        let position_specs = position_specs
            .iter()
            .map(|ps| parse_positions(ps))
            .collect::<Result<Vec<_>, _>>()?;

        let kmer_sizes: Vec<u8> = kmer_specs.keys().cloned().collect();

        Arc::new(PositionSelectionCache::new(
            position_specs,
            &kmer_sizes,
            opt.shared_args.fragment_lengths.min_fragment_length,
            opt.shared_args.fragment_lengths.max_fragment_length,
        )?)
    };

    // Load genomic scaling factors
    if opt.shared_args.scale_genome.scaling_factors.is_some() && !quiet {
        println!("Start: Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> =
        load_scaling_map(&opt.shared_args.scale_genome, &chromosomes, &contigs)?;

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.shared_args.ioc.output_dir, prefix)
        .context("create per-run temp dir")?;

    // Window size when --by-size (otherwise None)
    let by_size_bp: Option<u64> = match &window_opt {
        WindowSpec::Size(bp) => Some(*bp as u64),
        _ => None,
    };

    // Build tiles
    let halo_bp: u32 = opt.shared_args.fragment_lengths.max_fragment_length; // Safe halo for pairing/segments
    let (tiles, _tile_and_window_boundaries_align) = build_tiles(
        &chromosomes,
        &contigs,
        opt.shared_args.tile_size,
        halo_bp,
        by_size_bp,
    )?;

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(&tiles, |chr| {
        windows_lookup
            .and_then(|m| m.get(chr))
            .map(|w| w.as_slice())
            .unwrap_or(&[])
    }));

    // Compute per-chromosome window offsets and overall window count. In BED mode these offsets are
    // zero because windows already carry their global `original_idx` values.
    let (total_windows, chr_offsets_map) =
        compute_window_offsets(&window_opt, &chromosomes, &contigs, windows_map.as_ref())?;
    let chr_offsets = Arc::new(chr_offsets_map);

    let total_tiles = tiles.len();
    let temp_dir = Arc::new(temp_dir);

    // Create progress bar
    let pb = Arc::new(ProgressBar::new(total_tiles as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    if quiet {
        pb.set_draw_target(ProgressDrawTarget::hidden());
    }

    // Configure global thread‐pool size
    init_global_pool(opt.shared_args.ioc.n_threads as usize)?;

    if !quiet {
        println!("Start: Counting per chromosome");
    }

    pb.set_position(0);

    let tile_window_spans_for_threads = tile_window_spans.clone();
    let positional_cache_for_threads = positional_cache.clone();

    let tile_results: Vec<TileResult> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<TileResult> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let counts_path = temp_dir.join(format!(
                "{prefix}.{chr}.{idx}.counts.bin",
                prefix = prefix,
                chr = tile.chr.as_str(),
                idx = tile.index
            ));

            let window_ctx = WindowContext {
                spec: &window_opt,
                windows: windows_map
                    .as_ref()
                    .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice())),
                chr_idx_offset: *chr_offsets.get(&tile.chr).unwrap_or(&0),
            };

            let blacklist_chr = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            let position_cache = positional_cache_for_threads.clone();
            let out = process_tile(
                opt,
                tile,
                &kmer_specs,
                position_cache,
                &window_ctx,
                tile_span.as_ref(),
                blacklist_chr,
                scaling_chr,
                counts_path.as_path(),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    if !quiet {
        pb.finish_with_message("| Finished counting");
    } else {
        pb.finish_and_clear();
    }

    if !quiet {
        println!("Start: Reducing per-tile counts");
    }

    let mut global_counter = FragmentKmersCounters::default();
    let mut tile_results_by_chr: FxHashMap<String, Vec<TileResult>> = FxHashMap::default();

    for tile_result in tile_results {
        global_counter += tile_result.counter;
        tile_results_by_chr
            .entry(tile_result.chr.clone())
            .or_default()
            .push(tile_result);
    }

    let mut payloads: Vec<Vec<TileWindowCounts>> = Vec::with_capacity(tile_results_by_chr.len());
    for chr in &chromosomes {
        if let Some(chr_tile_results) = tile_results_by_chr.remove(chr) {
            payloads.push(reduce_chromosome_tile_results(chr_tile_results)?);
        }
    }
    if !tile_results_by_chr.is_empty() {
        let unexpected_chr = tile_results_by_chr.keys().next().unwrap();
        bail!(
            "tile results produced for unexpected chromosome '{}'",
            unexpected_chr
        );
    }

    let total_windows_usize: usize = total_windows
        .try_into()
        .context("number of windows exceeds addressable size")?;

    if opt.positional_counts {
        let positional_bins = merge_tile_counts_positional(payloads, total_windows_usize)?;

        let mut positional_decoded: Vec<FxHashMap<PositionDescriptor, DecodedCounts>> =
            Vec::with_capacity(positional_bins.len());
        let mut flattened: Vec<DecodedCounts> = Vec::new();

        for window_counts in positional_bins {
            let mut decoded_map: FxHashMap<PositionDescriptor, DecodedCounts> =
                FxHashMap::default();
            decoded_map.reserve(window_counts.len());
            for (descriptor, counts) in window_counts {
                let decoded = split_and_decode_counts(&counts, &kmer_specs);
                flattened.push(decoded.clone());
                decoded_map.insert(descriptor, decoded);
            }
            positional_decoded.push(decoded_map);
        }

        if flattened.is_empty() {
            flattened.push(DecodedCounts {
                counts: FxHashMap::default(),
            });
        }

        let (_, motifs_by_k) = prepare_decoded_counts(&flattened, opt.canonical, &kmer_specs);

        if !quiet {
            println!("Start: Writing positional counts to disk");
        }
        write_positional_output(
            &positional_decoded,
            &motifs_by_k,
            &kmer_specs,
            &opt.shared_args.ioc.output_dir,
            &opt.shared_args.output_prefix,
            opt.save_sparse,
        )?;
    } else {
        let all_bins = merge_tile_counts(payloads, total_windows_usize, &kmer_specs)?;

        // Prepare counts to get correct motifs (collapsed, N-filtered, etc.)
        let (prepared_counts, motifs_by_k) =
            prepare_decoded_counts(&all_bins, opt.canonical, &kmer_specs);

        // Write final counts to output_dir
        if !quiet {
            println!("Start: Writing counts to disk");
        }
        write_decoded_counts_matrix(
            &prepared_counts,
            &kmer_specs,
            &motifs_by_k,
            &opt.shared_args.ioc.output_dir,
            &opt.shared_args.output_prefix,
            opt.save_sparse,
        )?;
    }

    // Build bin metadata when windowed
    let bin_info = build_bin_info(
        &window_opt,
        &chromosomes,
        &contigs,
        windows_map.as_ref(),
        &blacklist_map,
        chr_offsets.as_ref(),
    )?;

    // Write window coordinates as BED file to output_dir
    // Write bins BED file
    if !matches!(window_opt, WindowSpec::Global) {
        if !quiet {
            println!("Start: Writing window coordinates to disk");
        }
        let bins_path = opt.shared_args.ioc.output_dir.join("bins.bed");
        let mut bed_writer = create_text_writer(&bins_path).context("Create bed fail")?;
        for (chr, start, end, _, overlap_perc) in &bin_info {
            writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
                .context("Write bed line fail")?;
        }
        bed_writer.finish().context("Finalize bins.bed writer")?;
    }

    Ok(global_counter)
}

/// Process a single tile: stream fragments, accumulate per-window counts, and persist results.
fn process_tile(
    opt: &FragmentKmersConfig,
    tile: &Tile,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    position_cache: Arc<PositionSelectionCache>,
    window_ctx: &WindowContext,
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_intervals: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    counts_path: &Path,
) -> anyhow::Result<TileResult> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) =
        create_chromosome_reader(&opt.shared_args.ioc.bam, &tile.chr)?;

    let fetch_span = determine_fetch_span(tile, window_ctx, tile_window_span, chrom_len);
    let Some((fetch_from, fetch_to)) = fetch_span else {
        return Ok(TileResult {
            chr: tile.chr.clone(),
            counts_path: None,
            counter: FragmentKmersCounters::default(),
        });
    };

    // Extend the reference slice to include k-mers at the right tile edge
    let max_k: u32 = kmer_specs.keys().copied().max().unwrap_or(1) as u32;
    let seq_end_abs = (tile.core_end as u64)
        .saturating_add((max_k as u64).saturating_sub(1))
        .min(chrom_len) as usize;

    let mut seq_bytes = read_seq_in_range(
        &opt.shared_args.ref_genome.ref_2bit,
        &tile.chr,
        (tile.core_start as usize)..(seq_end_abs),
    )?;

    apply_blacklist_mask_to_seq(&mut seq_bytes, &blacklist_intervals, tile.core_start as u64);

    // Scaled weights to count up
    let positional_scaling_weights = if !scaling_chr.is_empty() {
        let mut scaling_weights = vec![1.0; seq_bytes.len()];
        apply_scaling_to_coverage_in_place(
            &mut scaling_weights,
            tile.core_start as u32,
            scaling_chr,
        );
        // "Blacklist" positions with scaling factors of 0, so they don't get counted
        for (base, weight) in seq_bytes.iter_mut().zip(&scaling_weights) {
            if *weight == 0.0 {
                *base = BLACKLIST_BYTE;
            }
        }
        Some(scaling_weights)
    } else {
        None
    };

    // Prepare left-aligned kmer-codes for each kmer-size
    let positional_codes_by_k: FxHashMap<u8, KmerCodes> =
        build_left_aligned_codes_per_k(&seq_bytes, kmer_specs);

    // Sparse map keyed by original window index -> kmer counts
    let mut counts_by_window: FxHashMap<u64, FxHashMap<CountKey, f64>> = FxHashMap::default();

    // Streaming pointers and single fetch for this chr
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = tile_window_span
        .and_then(|span| (!span.is_empty()).then_some(span.first_idx))
        .unwrap_or(0);

    reader
        .fetch((tid, fetch_from, fetch_to))
        .context(format!("fetch {}", &tile.chr))?;

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let lengths = opt.shared_args.fragment_lengths.clone();
        move |f: &FragmentWithKmerSegments| lengths.contains(f.len())
    };

    // Wrap to use opt
    let include_read_fn = {
        let opt = (*opt).clone();
        move |r: &Record| {
            default_include_read(
                r,
                opt.shared_args.require_proper_pair,
                opt.shared_args.min_mapq,
            )
        }
    };

    // Create fragment iterator
    let mut iter = fragments_with_kmer_segments_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        include_read_fn,
        opt.shared_args.indel_mode,
        !opt.shared_args.ignore_gap,
        0,
        fragment_filter,
    )
    .with_local_counters();

    // Initialize counters (default -> 0s)
    let mut counter = FragmentKmersCounters::default();

    let store_positions = opt.positional_counts;
    let has_nearest_frame = position_cache
        .present_frames
        .contains(&ReferenceFrame::Nearest);

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;

        // Determine blacklist status
        let in_blacklist = is_blacklisted(
            blacklist_intervals,
            opt.shared_args.blacklist_strategy.clone(),
            fragment.start.into(),
            fragment.end.into(),
            opt.shared_args.fragment_lengths.max_fragment_length as u64,
            &mut bl_ptr,
        );
        if in_blacklist {
            counter.blacklisted_fragments += 1;
            continue;
        }

        let cache = position_cache.as_ref();
        let (first, last) = cache
            // Use smallest possible k to include all positions in interval for overlap!
            .bounds(fragment.len(), 1u8)
            .expect("non-empty offsets must have bounds");
        let interval_start = fragment.start as u64 + first as u64;
        let interval_end = fragment.start as u64 + last as u64 + 1;
        if interval_start >= interval_end {
            continue;
        }

        // Find all overlapping count-windows
        debug_assert!(interval_start >= fragment.start as u64);
        let lookback_distance = opt.shared_args.fragment_lengths.max_fragment_length as u64
            + (interval_start - fragment.start as u64);
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            window_ctx.windows_slice(),
            opt.shared_args.windows.by_size,
            interval_start,
            interval_end,
            1. / (opt.shared_args.fragment_lengths.max_fragment_length as f64 + 1.0), // Any overlap
            lookback_distance,
        )?;
        let overlapping_windows = if let Some(overlaps) = overlapping_windows {
            overlaps
        } else {
            continue;
        };

        counter.base.counted_fragments += 1;

        for overlapped_window in overlapping_windows.windows {
            let original_idx = window_ctx.original_idx(overlapped_window.idx);
            let counts = counts_by_window
                .entry(original_idx)
                .or_insert_with(FxHashMap::default);
            count_kmers_at_positions(
                &fragment,
                cache,
                store_positions,
                &positional_codes_by_k,
                kmer_specs,
                counts,
                positional_scaling_weights.as_deref(),
                tile.core_start,
                tile.core_end,
                has_nearest_frame,
            );
        }
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    let mut payload: Vec<TileWindowCounts> = counts_by_window
        .into_iter()
        .filter_map(|(original_idx, hm)| {
            if hm.is_empty() {
                return None;
            }
            let mut entries: Vec<TileKmerCountEntry> = Vec::with_capacity(hm.len());
            for (key, value) in hm {
                entries.push(TileKmerCountEntry::from((key, value)));
            }
            Some(TileWindowCounts {
                original_idx,
                entries,
            })
        })
        .collect();
    payload.sort_unstable_by_key(|w| w.original_idx);

    serialize_tile_counts(counts_path, &payload)?;

    Ok(TileResult {
        chr: tile.chr.clone(),
        counts_path: Some(counts_path.to_path_buf()),
        counter,
    })
}

pub fn count_kmers_at_positions(
    fragment: &FragmentWithKmerSegments,
    cache: &PositionSelectionCache,
    store_positions: bool,
    positional_codes_by_k: &FxHashMap<u8, KmerCodes>,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    counts: &mut FxHashMap<CountKey, f64>,
    weights: Option<&[f32]>,
    tile_core_start: u32,
    tile_core_end: u32,
    apply_nearest_guard: bool,
) {
    // We perform comparisons in absolute genome coordinates first, then translate
    // back to fragment-relative offsets only after clipping
    let fragment_start = fragment.start as u64;
    let tile_start = tile_core_start as u64;
    let tile_end = tile_core_end as u64;

    // We walk the requested k values independently. Each k has its own positional
    // encoding table, so processing them in isolation keeps the hot loop simple
    for (&k, _) in kmer_specs {
        let codes = positional_codes_by_k
            .get(&k)
            .expect("missing positional codes for requested k");
        let k_span = k as u64;
        let selections = match cache.offsets(fragment.len(), k) {
            Some(slice) if !slice.is_empty() => slice,
            _ => {
                continue;
            }
        };
        if selections.is_empty() {
            // Some frames filter out every position for a fragment of a given length
            return;
        }

        // Selections are sorted by offset. We stream through them once per k
        // using a single cursor so the overall complexity stays linear in the number
        // of usable offsets
        let mut offset_cursor = 0usize;

        // Precompute midpoint guards for Nearest: forbid crossing the true midpoint.
        // Rule:
        // - If there IS a physical midpoint (odd length), exclude that base entirely:
        //     forward:  start + (k-1) <= mid-1  -> start <= mid - k
        //     reverse:  start >= mid+1          -> anchor(offset) >= mid + (k-1) + 1 = mid + k
        // - If there is NO physical midpoint (even length), pick the base nearest each side's start:
        //     left  boundary = L/2 - 1 (0-based), right boundary = L/2
        //     forward:  start + (k-1) <= left_boundary   -> start <= (L/2) - k
        //     reverse:  start >= right_boundary          -> anchor(offset) >= (L/2) + (k-1)
        let nearest_guard =
            NearestFrameGuard::by_flag(apply_nearest_guard, fragment.len() as u32, k as u32);

        // Fragments may be gapped by indels, so we examine each contiguous segment
        // and clip it to the tile coordinates before accepting offsets
        'segments: for &(seg_start_raw, seg_end_raw) in &fragment.segments {
            let seg_start = seg_start_raw as u64;
            let seg_end = seg_end_raw as u64;
            if seg_start >= seg_end {
                continue;
            }

            if seg_end <= tile_start {
                // Segment lies completely before the tile
                continue;
            }

            if offset_cursor >= selections.len() {
                // Offsets arrive sorted and we only advance the cursor, therefore exhausting the
                // selections list guarantees no later segment has a usable position
                break;
            }

            if seg_start >= tile_end {
                // Segments are emitted in genomic order, so hitting the tile boundary means the rest
                // start beyond the core and cannot contribute any k-mers
                break;
            }

            let forward_min_abs = seg_start.max(tile_start);
            let effective_seg_end = seg_end.min(tile_end.saturating_add(k_span.saturating_sub(1)));
            let Some(last_start_abs) = effective_seg_end.checked_sub(k_span) else {
                continue;
            };

            // Forward oriented kmers start at the position we count
            // Only offsets whose entire span stays inside both the segment and tile are valid
            let forward_range = if last_start_abs >= forward_min_abs {
                Some((
                    forward_min_abs.saturating_sub(fragment_start),
                    last_start_abs.saturating_sub(fragment_start),
                ))
            } else {
                None
            };

            let reverse_anchor_min_abs = seg_start
                .saturating_add(k_span.saturating_sub(1))
                .max(tile_start.saturating_add(k_span.saturating_sub(1)));
            let reverse_anchor_max_abs = seg_end.min(tile_end);
            let reverse_range = if reverse_anchor_max_abs == 0 {
                None
            } else {
                let max_inclusive = reverse_anchor_max_abs.saturating_sub(1);
                if max_inclusive >= reverse_anchor_min_abs {
                    // Reverse oriented kmers are indexed by their last base
                    // After clipping to the tile and segment we backtrack by k-1
                    // bases to locate the start
                    Some((
                        reverse_anchor_min_abs.saturating_sub(fragment_start),
                        max_inclusive.saturating_sub(fragment_start),
                    ))
                } else {
                    None
                }
            };

            if forward_range.is_none() && reverse_range.is_none() {
                // Clipping removed every valid orientation for this segment
                continue;
            }

            // Build an inclusive span that covers whichever orientations survived
            let segment_range_start = match (forward_range, reverse_range) {
                (Some((fwd_min, _)), Some((rev_min, _))) => fwd_min.min(rev_min),
                (Some((fwd_min, _)), None) => fwd_min,
                (None, Some((rev_min, _))) => rev_min,
                (None, None) => unreachable!(),
            };

            let segment_range_end = match (forward_range, reverse_range) {
                (Some((_, fwd_max)), Some((_, rev_max))) => fwd_max.max(rev_max),
                (Some((_, fwd_max)), None) => fwd_max,
                (None, Some((_, rev_max))) => rev_max,
                (None, None) => unreachable!(),
            };

            while offset_cursor < selections.len()
                && (selections[offset_cursor].offset() as u64) < segment_range_start
            {
                // The cursor still points before the current segment window, so fast-forward it
                offset_cursor += 1;
            }

            let mut idx = offset_cursor;
            while idx < selections.len() {
                let selection = selections[idx];
                let offset = selection.offset() as u64;
                let offset_i32 = selection.offset() as i32;
                if offset > segment_range_end {
                    // Remaining selections start after this segment range, so move to next segment
                    break;
                }

                let idx_abs = fragment_start + offset;
                if idx_abs >= tile_end {
                    // Offsets are ordered, so reaching the right edge of the tile means we are done for this k
                    break 'segments;
                }
                if idx_abs < tile_start {
                    // Offset still lies before the tile core; advance to the next candidate
                    idx += 1;
                    continue;
                }

                let decision = evaluate_selection(
                    selection,
                    nearest_guard.as_ref(),
                    k_span,
                    offset,
                    forward_range,
                    reverse_range,
                );

                match decision {
                    SelectionDecision::SkipAdvance => {
                        idx += 1;
                        continue;
                    }
                    SelectionDecision::IncludeForward { .. } => {
                        let start_local = match idx_abs.checked_sub(tile_start) {
                            Some(val) => val as usize,
                            None => {
                                idx += 1;
                                continue;
                            }
                        };

                        // Ensure the forward k-mer stays within this contiguous segment
                        // idx_abs is the start. Require idx_abs + (k-1) < seg_end
                        if idx_abs.saturating_add(k_span.saturating_sub(1)) >= seg_end {
                            idx += 1;
                            continue;
                        }

                        // We look up weights using the same tile-relative index as the positional codes
                        let weight = match weights {
                            Some(w) => unsafe { *w.get_unchecked(start_local) as f64 },
                            None => 1.0,
                        };

                        // Record the forward kmer code emitted at this start position
                        let key = CountKey {
                            k,
                            code: codes.get(start_local),
                            position: store_positions.then_some(offset_i32),
                            group: selection.group(),
                        };
                        *counts.entry(key).or_insert(0.0) += weight;
                    }
                    SelectionDecision::IncludeReverse {
                        start_offset_0,
                        anchor_offset_0,
                    } => {
                        let kmer_start_abs = fragment_start + start_offset_0;
                        if kmer_start_abs < tile_start || kmer_start_abs < seg_start {
                            idx += 1;
                            continue;
                        }

                        let start_local = match kmer_start_abs.checked_sub(tile_start) {
                            Some(val) => val as usize,
                            None => {
                                idx += 1;
                                continue;
                            }
                        };

                        let end_local =
                            match (fragment_start + anchor_offset_0).checked_sub(tile_start) {
                                Some(val) => val as usize,
                                None => {
                                    idx += 1;
                                    continue;
                                }
                            };

                        let weight = match weights {
                            Some(w) => unsafe { *w.get_unchecked(end_local) as f64 },
                            None => 1.0,
                        };

                        let key = CountKey {
                            k,
                            code: codes.get(start_local),
                            position: store_positions.then_some(offset_i32),
                            group: selection.group(),
                        };
                        *counts.entry(key).or_insert(0.0) += weight;
                    }
                }

                idx += 1;
            }

            // Carry the cursor forward so the next segment starts scanning from the last visited offset
            offset_cursor = idx;
        }
    }
}
