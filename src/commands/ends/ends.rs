use crate::{
    commands::{
        cli_common::{
            WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        counters::EndsCounters,
        ends::{
            config::EndsConfig,
            config_structs::WindowMotifAssigner,
            counting::{EndCountsByWindow, decode_end_motif_counts},
            motifs::{
                build_optional_kmer_spec, build_tile_motif_context, count_fragment_in_window,
            },
            output::{
                build_all_end_motif_order, collect_end_motif_order,
                ensure_all_motifs_enumeration_size, ensure_dense_end_motif_output_size,
            },
            tiling::{
                TileResult, build_tile_payload, deserialize_tile_counts, merge_tile_payload,
                serialize_tile_counts,
            },
            write::{write_end_motif_outputs, write_end_settings_json},
        },
        gc_bias::{
            correct::{GCCorrector, load_gc_corrector},
            counting::build_gc_prefixes,
        },
        lengths::tiling::fetch_span_for_tile,
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::is_blacklisted,
        fragment::ends_fragment::FragmentWithEnds,
        fragment_iterators::fragments_with_ends_from_bam,
        interval::{IndexedInterval, Interval},
        io::{create_text_writer, dot_join},
        midpoint::midpoint_random_even_with_thread_rng,
        overlaps::find_overlapping_windows,
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq_in_range,
        scale_genome::{compute_window_scaling_over_fragment, compute_window_scaling_over_overlap},
        thread_pool::init_global_pool,
        tiled_run::{
            Tile, TileWindowSpan, build_tiles, make_temp_dir, precompute_tile_window_spans,
        },
        windowing::{WindowContext, build_bin_info, compute_window_offsets},
    },
};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{convert::TryInto, io::Write, path::Path, sync::Arc, time::Instant};

/// Execute the end-motif counting pipeline end-to-end.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `ends` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &EndsConfig) -> Result<()> {
    let start_time = Instant::now();
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    if opt.k_within == 0 && opt.k_outside == 0 {
        bail!("At least one of --k-within or --k-outside must be > 0");
    }
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    let window_opt = opt.windows.resolve_windows();
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

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            println!("Start: Loading window coordinates");
            Some(load_windows_from_bed(
                bed,
                Some(chromosomes.as_slice()),
                None,
                None,
            )?)
        }
        _ => None,
    };

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
        opt.fragment_lengths.min_fragment_length,
        opt.fragment_lengths.max_fragment_length,
    )?;

    let halo_bp = opt.fragment_lengths.max_fragment_length;
    let align_bp = match &window_opt {
        WindowSpec::Size(bp) => Some(*bp),
        _ => None,
    };

    // Build tiles (core plus halo)
    let (tiles, _) = build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, align_bp)?;

    let progress = ProgressFactory::new();
    let pb = Arc::new(progress.default_bar(tiles.len() as u64));

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(
        &tiles,
        |chr| {
            windows_lookup
                .and_then(|m| m.get(chr).map(|w| w.as_slice()))
                .unwrap_or(&[])
        },
        0,
        // We use fragments starting in a tile, so we need fragment-overlapping windows starting after the tile
        opt.fragment_lengths.max_fragment_length as u64,
    ));
    let tile_window_spans_for_threads = tile_window_spans.clone();
    // Window rows are global across chromosomes. For fixed-size windows we therefore need a
    // per-chromosome row offset to turn chromosome-local overlap indices into global output rows.
    // BED windows already carry their own original indices, so their offsets stay at zero.
    let (total_windows, chr_offsets_map) =
        compute_window_offsets(&window_opt, &chromosomes, &contigs, windows_map.as_ref())?;
    let chr_offsets = Arc::new(chr_offsets_map);
    let chr_offsets_for_threads = chr_offsets.clone();

    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;
    let counts_prefix = &dot_join(&[prefix, "counts"]);

    println!("Start: Counting per tile");

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    pb.set_position(0);

    let tile_results: Vec<TileResult> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<Option<TileResult>> {
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
            // `find_overlapping_windows` reports chromosome-local indices. `process_tile` needs the
            // chromosome-specific offset so it can emit globally stable window ids into tile payloads.
            let chr_window_idx_offset = *chr_offsets_for_threads.get(&tile.chr).unwrap_or(&0);

            let tile_result = process_tile(
                opt,
                tile,
                tile_span.as_ref(),
                windows_chr,
                chr_window_idx_offset,
                &window_opt,
                blacklist_chr,
                scaling_chr,
                gc_corrector.clone(),
                &temp_dir,
                counts_prefix,
            )?;
            pb.inc(1);
            Ok(tile_result)
        })
        .collect::<Result<Vec<_>>>()? // Short-circuits on the first Err
        .into_iter()
        .flatten()
        .collect();

    pb.finish_with_message("| Finished counting");

    // Release per-tile inputs before merging outputs
    drop(chr_offsets_for_threads);
    drop(tile_window_spans_for_threads);
    drop(tile_window_spans);
    drop(tiles);
    drop(scaling_map);
    drop(gc_corrector);

    // Collect counters
    let mut global_counter = EndsCounters::default();
    for tile_out in &tile_results {
        global_counter += tile_out.counter;
    }

    println!("Start: Reducing temporary tile files");

    let within_spec = build_optional_kmer_spec(opt.k_within, "within")?;
    let outside_spec = build_optional_kmer_spec(opt.k_outside, "outside")?;
    // Start from an all-empty output matrix shape, then fill only the windows that were actually
    // observed in the reduced sparse payloads.
    let mut all_bins = vec![FxHashMap::default(); total_windows as usize];
    let mut reduced_counts: EndCountsByWindow = FxHashMap::default();
    for tile_result in &tile_results {
        merge_tile_payload(
            &mut reduced_counts,
            deserialize_tile_counts(&tile_result.counts_path)?,
        )?;
    }

    // Decode each populated window independently. This is an easy parallel boundary because the
    // windows no longer interact after tile reduction; we only need a final serial pass to place
    // each decoded map into its global output row.
    let decoded_bins: Vec<(usize, FxHashMap<String, f64>)> = reduced_counts
        .into_par_iter()
        .map(
            |(original_idx, counts)| -> Result<(usize, FxHashMap<String, f64>)> {
                let decoded = decode_end_motif_counts(
                    &counts,
                    within_spec.as_ref(),
                    outside_spec.as_ref(),
                    opt.collapse_complement,
                );
                let idx: usize = original_idx
                    .try_into()
                    .context("window index does not fit in usize")?;
                Ok((idx, decoded))
            },
        )
        .collect::<Result<_>>()?;

    for (idx, decoded) in decoded_bins {
        if idx >= all_bins.len() {
            bail!(
                "reduced window index {} is out of bounds for {} output windows",
                idx,
                all_bins.len()
            );
        }
        all_bins[idx] = decoded;
    }

    let bin_info = build_bin_info(
        &window_opt,
        &chromosomes,
        &contigs,
        windows_map.as_ref(),
        &blacklist_map,
        chr_offsets.as_ref(),
    )?;
    // `all_motifs` switches the final output from "observed motifs only" to a dense fixed motif
    // universe. The dense size checks happen before we allocate or enumerate that full universe.
    if opt.all_motifs {
        ensure_all_motifs_enumeration_size(opt.k_within, opt.k_outside, all_bins.len())?;
    }
    let motif_order = if opt.all_motifs {
        build_all_end_motif_order(
            within_spec.as_ref(),
            outside_spec.as_ref(),
            opt.collapse_complement,
        )?
    } else {
        collect_end_motif_order(&all_bins)
    };
    let write_dense_output = opt.all_motifs;
    if write_dense_output {
        ensure_dense_end_motif_output_size(all_bins.len(), motif_order.len())?;
    }
    write_end_motif_outputs(
        &opt.ioc.output_dir,
        prefix,
        &all_bins,
        &motif_order,
        write_dense_output,
    )?;

    drop(blacklist_map);

    let keep_temp = false;
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

    write_end_settings_json(&opt.ioc.output_dir, prefix, opt)?;

    // Write window coordinates as BED file to output_dir
    // Write bins BED file
    if !matches!(window_opt, WindowSpec::Global) {
        println!("Start: Writing window coordinates to disk");
        let bins_path = opt.ioc.output_dir.join(dot_join(&[prefix, "bins.bed"]));
        let mut bed_writer = create_text_writer(&bins_path).context("Create bed fail")?;
        for (chr, start, end, _, overlap_perc) in &bin_info {
            writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
                .context("Write bed line fail")?;
        }
        bed_writer.finish().context("Finalize bins.bed writer")?;
    }

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
    if opt.gc.gc_file.is_some() {
        let gc_fail_action = if opt.gc.drop_invalid_gc {
            "fragment skipped"
        } else {
            "fragment counted with weight 1.0"
        };
        println!(
            "  GC correction failures ({}): {}",
            gc_fail_action, global_counter.gc_failed_fragments
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

/// Count all end motifs owned by one tile and write its sparse payload to disk.
///
/// This function does the tile-local heavy lifting for `ends`: it streams
/// fragments from BAM, applies fragment-level filters and weights, assigns each
/// owned fragment to candidate windows, and accumulates sparse motif counts for
/// those windows before serializing the tile result.
///
/// Parameters
/// ----------
/// - `opt`:
///   Full command configuration
/// - `tile`:
///   Tile currently being processed
/// - `tile_window_span`:
///   Precomputed window span relevant to this tile, if any
/// - `windows_chr`:
///   Window intervals for this chromosome, when not in global mode
/// - `chr_window_idx_offset`:
///   Per-chromosome row offset for fixed-size windows
/// - `window_opt`:
///   Resolved windowing mode
/// - `blacklist_intervals`:
///   Merged blacklist intervals for this chromosome
/// - `scaling_chr`:
///   Genomic scaling bins for this chromosome
/// - `gc_corrector_opt`:
///   Optional GC corrector shared from the outer runner
/// - `temp_dir`:
///   Temporary directory for tile payloads
/// - `counts_prefix`:
///   Prefix used when naming serialized tile payloads
///
/// Returns
/// -------
/// - `Result<Option<TileResult>>`:
///   `None` when the tile has no relevant windows, otherwise the serialized sparse tile result
fn process_tile(
    opt: &EndsConfig,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    windows_chr: Option<&[IndexedInterval<u64>]>,
    chr_window_idx_offset: u64,
    window_opt: &WindowSpec,
    blacklist_intervals: &[Interval<u64>],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<GCCorrector>,
    temp_dir: &Path,
    counts_prefix: &str,
) -> Result<Option<TileResult>> {
    // One BAM reader per tile
    let (mut reader, _tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    debug_assert_eq!(_tid_check, tile.tid as u32);

    // Counters
    let mut counter = EndsCounters::default();
    let counts_path = temp_dir.join(format!(
        "{prefix}.{chr}.{idx}.counts.bin",
        prefix = counts_prefix,
        chr = tile.chr.as_str(),
        idx = tile.index
    ));

    // Build GC prefixes for the full tile fetch span
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

    // Narrow the BAM fetch to the part of the tile that can still contribute to the current
    // windows. In global/by-size modes this usually stays close to the tile fetch span; in BED
    // mode it can shrink substantially.
    let Some(fetch_span) = fetch_span_for_tile(
        tile,
        tile_window_span,
        windows_chr,
        window_opt,
        chrom_len,
        opt.fragment_lengths.max_fragment_length as u64,
    )?
    else {
        // Skip tiles with no relevant windows
        return Ok(None);
    };
    let (fetch_from, fetch_to) = fetch_span.try_to_i64()?.as_tuple();
    let motif_context =
        build_tile_motif_context(opt, tile, fetch_span, chrom_len, blacklist_intervals)?;
    let window_context = WindowContext {
        spec: window_opt,
        windows: windows_chr,
        chr_idx_offset: chr_window_idx_offset,
    };

    reader
        .fetch((tile.tid, fetch_from, fetch_to))
        .context(format!("fetch {} {}-{}", &tile.chr, fetch_from, fetch_to))?;

    let mut counts_by_window: EndCountsByWindow = FxHashMap::default();

    // Fraction of a fragment that must overlap with a window to consider that window as a
    // candidate. Endpoint mode still uses the fragment assignment interval here; the actual
    // left/right terminal-base checks happen later in `count_fragment_in_window(...)`.
    let min_overlap_fraction: f64 = match opt.window_assignment.assign_by {
        WindowMotifAssigner::Any
        | WindowMotifAssigner::Endpoint
        | WindowMotifAssigner::CountOverlap => {
            1. / (2. * opt.fragment_lengths.max_fragment_length as f64 + 1.0)
        } // 2x to allow for raw-clipping-mode expansion and +1 to avoid rounding error issues
        WindowMotifAssigner::All | WindowMotifAssigner::Midpoint => {
            1.0 - (1. / (2. * opt.fragment_lengths.max_fragment_length as f64 + 1.0))
        } // 1.0 but just below to avoid rounding errors
        WindowMotifAssigner::Proportion(p) => p,
    };

    // The overlap finder only needs checked BED-like intervals here.
    //
    // In BED mode, `find_overlapping_windows(...)` returns scan positions in the supplied slice as
    // `OverlappingWindow.idx`; it does not use `IndexedInterval.idx`. This temporary list is built
    // in the same order as `scaling_chr`, so those scan positions are already the correct
    // chromosome-local indices for indexing back into `scaling_chr` later.
    //
    // So the carried `IndexedInterval.idx` value is intentionally a placeholder.
    let scaling_with_bin_idx: Vec<IndexedInterval<u64>> = scaling_chr
        .iter()
        .map(|(start, end, _)| IndexedInterval::new(*start, *end, 0_u64))
        .collect::<crate::Result<_>>()?;

    // Function for filtering fragments after pairing
    let fragment_filter = {
        let lengths = opt.fragment_lengths.clone();
        move |f: &FragmentWithEnds| lengths.contains(f.len())
    };

    // Create fragment iterator with per-tile filtering and optional GC tag handling
    let unpaired = opt.unpaired.reads_are_fragments;
    let max_soft_clips = opt
        .clip
        .max_soft_clips
        .map(|value| {
            value
                .try_into()
                .context("max_soft_clips does not fit in u32 for fragment iteration")
        })
        .transpose()?;
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
    let mut iter = fragments_with_ends_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        move |rec| include_read_fn(rec),
        opt.clip.clip_strategy,
        opt.source_within,
        opt.indel_filter,
        opt.k_within,
        max_soft_clips,
        opt.gc.gc_tag.as_deref().map(str::as_bytes),
        fragment_filter,
        unpaired,
    )
    .with_local_counters();

    let get_gc_weight = {
        let gc_corrector = gc_corrector_opt.as_ref();
        let gc_prefixes = gc_prefixes_opt.as_ref();
        let fetch_start = tile.fetch_start();
        move |fragment: &FragmentWithEnds| -> Result<Option<f64>> {
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

    // Streaming pointers
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = tile_window_span
        .and_then(|span| (!span.is_empty()).then_some(span.first_idx))
        .unwrap_or(0);
    let mut sf_ptr = 0; // Scaling factor bin

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;
        let fragment_assignment_length = fragment.assignment_len();

        // Only count fragments whose start is inside the core to prevent double counting across tiles
        if fragment.start() < tile.core_start() || fragment.start() >= tile.core_end() {
            continue;
        }

        // Determine blacklist status
        let in_blacklist = is_blacklisted(
            blacklist_intervals,
            opt.blacklist_strategy,
            fragment.interval.try_to_u64()?,
            opt.fragment_lengths.max_fragment_length as u64,
            &mut bl_ptr,
        );
        if in_blacklist {
            counter.blacklisted_fragments += 1;
            continue;
        }

        // First find candidate windows from fragment geometry alone. We intentionally do this
        // before motif extraction so all later work only happens for windows that can actually
        // receive counts.
        let query_interval = match opt.window_assignment.assign_by {
            WindowMotifAssigner::Midpoint => {
                let midpoint = midpoint_random_even_with_thread_rng(
                    fragment.assignment_start(),
                    fragment_assignment_length,
                );
                Interval::new(midpoint.into(), (midpoint + 1).into())?
            }
            WindowMotifAssigner::Any
            | WindowMotifAssigner::All
            | WindowMotifAssigner::Proportion(_)
            | WindowMotifAssigner::CountOverlap
            | WindowMotifAssigner::Endpoint => fragment.assignment_interval.try_to_u64()?,
        };
        let by_size = match window_opt {
            WindowSpec::Size(bp) => Some(*bp),
            _ => None,
        };
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            windows_chr,
            by_size,
            query_interval,
            min_overlap_fraction,
            opt.fragment_lengths.max_fragment_length.into(),
        )?;
        let overlapping_windows = if let Some(overlaps) = overlapping_windows {
            overlaps
        } else {
            continue;
        };

        // GC correction is fragment-level, so the same GC weight is reused for every window and
        // every end motif produced from this fragment.
        let gc_weight_opt = get_gc_weight(&fragment)?;
        let gc_weight = match (gc_weight_opt, correct_gc) {
            (Some(w), true) => w,
            (None, true) => {
                // Tried but failed to make a GC correction weight for the current fragment
                // Fall back to no correction or skip
                counter.gc_failed_fragments += 1;
                if opt.gc.drop_invalid_gc {
                    continue;
                } else {
                    1.0
                }
            }
            (None, false) => 1.0, // No correction
            (Some(_), false) => bail!("unexpected GC weight when GC correction is disabled"),
        };

        counter.base.counted_fragments += 1;

        if !scaling_chr.is_empty() {
            // Find overlapping scaling-bins
            let overlapping_scaling_bins = find_overlapping_windows(
                chrom_len,
                &mut sf_ptr,
                Some(&scaling_with_bin_idx),
                None,
                fragment.interval.try_to_u64()?, // Full fragment
                1. / (2. * opt.fragment_lengths.max_fragment_length as f64 + 1.0), // Any overlap without rounding error issues
                opt.fragment_lengths.max_fragment_length.into(),
            )
            .with_context(|| format!("finding overlapping scaling bins on chr {}", tile.chr))?
            .context("no overlapping scaling bins found")?; // Should always find >= 1 bin

            // Extract the indices of the overlapping bins
            let overlapping_scaling_bin_indices: Vec<usize> = overlapping_scaling_bins
                .windows
                .iter()
                .map(|w| w.idx)
                .collect();

            // Calculate the weight per overlapping count-window
            // NOTE: `compute_window_scaling_over_fragment` always returns
            // an overlap fraction of 1.0 (count full fragment)!
            let overlap_weights = match opt.window_assignment.assign_by {
                WindowMotifAssigner::CountOverlap => compute_window_scaling_over_overlap(
                    &overlapping_windows,
                    &overlapping_scaling_bin_indices,
                    scaling_chr,
                )?,
                _ => compute_window_scaling_over_fragment(
                    fragment.interval.try_to_u64()?,
                    &overlapping_windows,
                    &overlapping_scaling_bin_indices,
                    scaling_chr,
                )?,
            };

            // Count up the weight per overlapping count-window. `count_fragment_in_window(...)`
            // still decides whether each end is actually counted in the current window.
            let overlapping_window_intervals: FxHashMap<usize, Interval<u64>> = overlapping_windows
                .windows
                .iter()
                .map(|window| (window.idx, window.interval))
                .collect();
            for (overlapped_window_idx, scaling_weight, overlap_fraction_to_count) in
                overlap_weights
            {
                let original_idx = window_context.original_idx(overlapped_window_idx);
                let window_interval = *overlapping_window_intervals
                    .get(&overlapped_window_idx)
                    .expect("missing overlap interval for scaled count window");
                let count_weight = overlap_fraction_to_count * scaling_weight * gc_weight;
                count_fragment_in_window(
                    &mut counts_by_window,
                    original_idx,
                    window_interval,
                    &fragment,
                    count_weight,
                    &motif_context,
                    opt.source_within,
                    opt.window_assignment.assign_by,
                )?;
            }
        } else {
            // Without genomic scaling, each candidate window gets either weight 1.0 or the raw
            // overlap fraction, depending on the assignment mode.
            for overlapped_window in overlapping_windows.windows {
                let original_idx = window_context.original_idx(overlapped_window.idx);
                let count_weight = match opt.window_assignment.assign_by {
                    WindowMotifAssigner::CountOverlap => overlapped_window.overlap_fraction as f64,
                    _ => 1.0f64,
                } * gc_weight;
                count_fragment_in_window(
                    &mut counts_by_window,
                    original_idx,
                    overlapped_window.interval,
                    &fragment,
                    count_weight,
                    &motif_context,
                    opt.source_within,
                    opt.window_assignment.assign_by,
                )?;
            }
        }
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    let payload = build_tile_payload(counts_by_window);
    serialize_tile_counts(&counts_path, &payload)?;

    Ok(Some(TileResult {
        chr: tile.chr.clone(),
        counts_path,
        counter,
    }))
}
