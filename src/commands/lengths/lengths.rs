use crate::{
    commands::{
        cli_common::{
            WindowAssigner, WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        counters::LengthsCounters,
        gc_bias::{
            correct::{LengthAgnosticGCCorrector, load_length_agnostic_gc_corrector},
            counting::build_gc_prefixes,
        },
        lengths::{
            config::LengthsConfig,
            counting::{LengthCounts, stack_length_counts},
            tiling::{reduce_partials_for_chr, write_cross_npy, write_partials_npz},
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::is_blacklisted,
        fragment::indel_counting_fragment::FragmentWithIndelCounts,
        fragment_iterators::fragments_with_indel_counts_from_bam,
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
        window_fetch::{BedFetchPolicy, fetch_span_for_tile},
        windowing::{build_bin_info, compute_window_offsets},
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{io::Write, path::Path, sync::Arc, time::Instant};

// Map orig_idx -> counts plus containment flag for this tile
#[derive(Clone)]
struct TileCounts {
    counts: LengthCounts,
    contained: bool,
}

#[derive(Clone)]
struct TileOutputs {
    counters: LengthsCounters,
    global_counts: Option<(String, LengthCounts)>,
}

/// Execute the fragment-length counting pipeline end-to-end.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `lengths` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &LengthsConfig) -> Result<()> {
    let start_time = Instant::now();
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
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
    let gc_corrector = load_length_agnostic_gc_corrector(
        opt.gc.gc_file.as_ref(),
        &opt.gc_length_weighting,
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

    // Reusable length-bin template so every tile/window counter shares identical bounds and avoids repeated allocations
    // Cloned with `zeroed_like` when building per-window `TileCounts`, which guarantees merge compatibility during reduction
    let template_counts = LengthCounts::new(
        opt.fragment_lengths.min_fragment_length as usize,
        opt.fragment_lengths.max_fragment_length as usize,
    );

    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;
    let partials_prefix = &dot_join(&[prefix, "part"]);
    let cross_prefix = &dot_join(&[prefix, "cross"]);

    println!("Start: Counting per tile");

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    pb.set_position(0);

    let tile_results: Vec<TileOutputs> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<TileOutputs> {
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

            let counter = process_tile(
                opt,
                tile,
                tile_span.as_ref(),
                windows_chr,
                &window_opt,
                blacklist_chr,
                scaling_chr,
                gc_corrector.clone(),
                &template_counts,
                &temp_dir,
                partials_prefix,
                cross_prefix,
            )?;
            pb.inc(1);
            Ok(counter)
        })
        .collect::<Result<_>>()?; // Short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    // Release per-tile inputs before merging outputs
    drop(tile_window_spans_for_threads);
    drop(tile_window_spans);
    drop(tiles);
    drop(scaling_map);
    drop(gc_corrector);

    // Collect counters
    let mut global_counter = LengthsCounters::default();
    for tile_out in &tile_results {
        global_counter += tile_out.counters;
    }

    println!("Start: Reducing temporary tile files");

    let mut all_bins: Vec<LengthCounts> = Vec::new();

    match &window_opt {
        WindowSpec::Global => {
            let mut counts_by_chr: FxHashMap<String, LengthCounts> = FxHashMap::default();
            for tile_out in &tile_results {
                if let Some((chr, counts)) = &tile_out.global_counts {
                    let entry = counts_by_chr
                        .entry(chr.clone())
                        .or_insert_with(|| template_counts.zeroed_like());
                    entry.merge_from(counts)?;
                }
            }
            for chr in &chromosomes {
                let counts = counts_by_chr
                    .remove(chr)
                    .context("Global mode missing counts for chromosome")?;
                all_bins.push(counts);
            }
            all_bins = vec![LengthCounts::collapse(&all_bins)?];
        }
        WindowSpec::Size(window_bp) => {
            for chr in &chromosomes {
                let chrom_len = contigs
                    .contigs
                    .get(chr)
                    .map(|&(_, len)| len as u64)
                    .context("missing contig length")?;
                let n_windows = chrom_len.div_ceil(*window_bp) as usize;
                let counts = reduce_partials_for_chr(
                    chr,
                    &temp_dir,
                    partials_prefix,
                    cross_prefix,
                    n_windows,
                    &template_counts,
                )?;
                ensure!(
                    counts.len() == n_windows,
                    "Expected {} windows for {} but got {}",
                    n_windows,
                    chr,
                    counts.len()
                );
                all_bins.extend(counts);
            }
        }
        WindowSpec::Bed(_) => {
            let win_map = windows_map
                .as_ref()
                .context("windows_map missing for BED mode")?;
            for chr in &chromosomes {
                let Some(wchr) = win_map.get(chr) else {
                    continue;
                };
                let wchr_slice = wchr.as_slice();
                if wchr_slice.is_empty() {
                    continue;
                }
                let counts = reduce_partials_for_chr(
                    chr,
                    &temp_dir,
                    partials_prefix,
                    cross_prefix,
                    wchr_slice.len(),
                    &template_counts,
                )?;
                ensure!(
                    counts.len() == wchr_slice.len(),
                    "Expected {} windows for {} but got {}",
                    wchr_slice.len(),
                    chr,
                    counts.len()
                );
                all_bins.extend(counts);
            }
        }
    }

    let mut bin_info: Vec<(String, u64, u64, u64, f64)> =
        if matches!(window_opt, WindowSpec::Global) {
            Vec::new()
        } else {
            let (_total_windows, chr_offsets) =
                compute_window_offsets(&window_opt, &chromosomes, &contigs, windows_map.as_ref())?;
            build_bin_info(
                &window_opt,
                &chromosomes,
                &contigs,
                windows_map.as_ref(),
                &blacklist_map,
                &chr_offsets,
            )?
        };

    drop(blacklist_map);

    // Sort by original index (when given a bed file)
    if matches!(window_opt, WindowSpec::Bed(_)) {
        println!("Start: Reordering counts by original window index in BED file");

        // Zip into a single Vec to allow sorting together
        let mut paired: Vec<_> = bin_info.into_iter().zip(all_bins).collect(); // (BinInfo, DecodedCounts)

        // Sort primarily by original window index
        paired.sort_unstable_by_key(|(info, _)| info.3);

        // Unzip back out if you need separate Vecs again
        (bin_info, all_bins) = paired.into_iter().unzip();
    }

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

    // Write final counts to output_dir
    write_npy(
        opt.ioc
            .output_dir
            .join(dot_join(&[prefix, "length_counts.npy"])),
        &stack_length_counts(&all_bins),
    )
    .context("Write final fail")?;

    // Write the min+max fragment length settings
    let settings_path = opt
        .ioc
        .output_dir
        .join(dot_join(&[prefix, "fragment_length_settings.json"]));
    let mut settings_writer =
        create_text_writer(&settings_path).context("Create fragment length settings file")?;
    writeln!(
        settings_writer,
        "{{\"min_fragment_length\":{},\"max_fragment_length\":{}}}",
        opt.fragment_lengths.min_fragment_length, opt.fragment_lengths.max_fragment_length
    )
    .context("Write fragment length settings")?;
    settings_writer
        .finish()
        .context("Finalize fragment length settings writer")?;

    // Plot the global fragment length distribution as a line plot for quick QC
    #[cfg(feature = "plotters")]
    {
        println!("Start: Plotting overall length distribution");

        use crate::shared::plotters::lineplot::write_line_plot_png;

        if all_bins.is_empty() {
            println!("Skipping overall length plot because no bins were produced");
        } else {
            let mut global_counts = vec![0f64; all_bins[0].counts.len()];
            for length_counts in &all_bins {
                for (total, count) in global_counts.iter_mut().zip(length_counts.counts.iter()) {
                    *total += *count;
                }
            }

            let total_counts: f64 = global_counts.iter().sum();
            if total_counts > 0.0 {
                for value in &mut global_counts {
                    *value /= total_counts;
                }
            }

            let x_values: Vec<f64> = (all_bins[0].length_min..=all_bins[0].length_max)
                .map(|len| len as f64)
                .collect();

            let plot_path = opt
                .ioc
                .output_dir
                .join(dot_join(&[prefix, "fragment_lengths_overall.png"]));

            write_line_plot_png(
                &plot_path,
                "Fragment length distribution (summed/global)",
                "Fragment length (bp)",
                "Density",
                &x_values,
                &global_counts,
                1600,
                1000,
            )
            .with_context(|| format!("writing fragment length plot to {}", plot_path.display()))?;
        }
    }

    // Write window coordinates plus overlap metadata as TSV to output_dir
    if !matches!(window_opt, WindowSpec::Global) {
        println!("Start: Writing window coordinates to disk");
        let bins_path = opt.ioc.output_dir.join(dot_join(&[prefix, "bins.tsv"]));
        let mut tsv_writer = create_text_writer(&bins_path).context("Create bins TSV fail")?;
        writeln!(tsv_writer, "chrom\tstart\tend\tblacklisted_fraction")
            .context("Write bins TSV header fail")?;
        for (chr, start, end, _, blacklist_overlap_fraction) in &bin_info {
            writeln!(
                tsv_writer,
                "{}\t{}\t{}\t{}",
                chr, start, end, blacklist_overlap_fraction
            )
            .context("Write bins TSV row fail")?;
        }
        tsv_writer.finish().context("Finalize bins.tsv writer")?;
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

fn process_tile(
    opt: &LengthsConfig,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    windows_chr: Option<&[IndexedInterval<u64>]>,
    window_opt: &WindowSpec,
    blacklist_intervals: &[Interval<u64>],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<LengthAgnosticGCCorrector>,
    template: &LengthCounts,
    temp_dir: &Path,
    partials_prefix: &str,
    cross_prefix: &str,
) -> Result<TileOutputs> {
    // One BAM reader per tile
    let (mut reader, _tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    debug_assert_eq!(_tid_check, tile.tid as u32);

    // Counters
    let mut counter = LengthsCounters::default();

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

    // Adapt the fetch coordinates to the present windows (*in windowed mode!*)
    let Some(fetch_span) = fetch_span_for_tile(
        tile,
        tile_window_span,
        windows_chr,
        window_opt,
        chrom_len,
        opt.fragment_lengths.max_fragment_length as u64,
        BedFetchPolicy::CandidateWindowExtent,
    )?
    else {
        // Skip tiles with no relevant windows
        return Ok(TileOutputs {
            counters: counter,
            global_counts: None,
        });
    };
    let (fetch_from, fetch_to) = fetch_span.try_to_i64()?.as_tuple();

    reader
        .fetch((tile.tid, fetch_from, fetch_to))
        .context(format!("fetch {} {}-{}", &tile.chr, fetch_from, fetch_to))?;

    // Preallocate per-tile window counters
    // Keep indices aligned with global scan order so downstream merging works without remapping
    // Use Option to skip BED windows that cannot be hit by any fragment starting in this tile
    // Track the first and last index covered to translate back to global coordinates
    let (counts_start_idx, counts_end_idx_exclusive, mut counts_by_idx): (
        usize,
        usize,
        Vec<Option<TileCounts>>,
    ) = match window_opt {
        // Global mode has exactly one window covering the chromosome
        WindowSpec::Global => (
            0,
            1,
            vec![Some(TileCounts {
                counts: template.zeroed_like(),
                contained: false,
            })],
        ),

        // Fixed-size mode: allocate only the bins that a fragment from this tile can reach
        WindowSpec::Size(window_bp) => {
            // Total bins on the chromosome
            let chrom_bin_count = chrom_len.div_ceil(*window_bp) as usize;
            // Leftmost bin whose start is at or before the core start
            // (may begin before the core when cores are not aligned)
            let min_bin_idx = (tile.core_start() as u64 / *window_bp) as usize;
            // Furthest coordinate a fragment starting in this tile can reach
            let max_reachable_end = (tile.core_end() as u64)
                .saturating_add(opt.fragment_lengths.max_fragment_length as u64)
                .min(chrom_len);
            // One past the last bin that could overlap that reach
            let max_bin_idx_exclusive = if max_reachable_end == 0 {
                min_bin_idx
            } else {
                (((max_reachable_end - 1) / *window_bp) + 1) as usize
            }
            .min(chrom_bin_count);

            let span_len = max_bin_idx_exclusive.saturating_sub(min_bin_idx);
            let mut counts = Vec::with_capacity(span_len);
            for idx in min_bin_idx..max_bin_idx_exclusive {
                let start = idx as u64 * *window_bp;
                let end = (start + *window_bp).min(chrom_len);
                // Contained means the bin sits fully inside the tile core
                let contained = start >= tile.core_start() as u64 && end <= tile.core_end() as u64;
                counts.push(Some(TileCounts {
                    counts: template.zeroed_like(),
                    contained,
                }));
            }
            (min_bin_idx, max_bin_idx_exclusive, counts)
        }

        // BED mode: reuse the precomputed span but "drop" windows ending before the core start
        WindowSpec::Bed(_) => {
            let span = tile_window_span.context(
                "BED length counting requires a cached tile window span after fetch-span selection",
            )?;
            let wchr = windows_chr.context("BED length counting requires loaded windows")?;
            let span_len = span.last_idx_exclusive.saturating_sub(span.first_idx);
            let mut counts = Vec::with_capacity(span_len);
            for idx in span.first_idx..span.last_idx_exclusive {
                let window = wchr[idx];
                let win_start = window.start();
                let win_end = window.end();
                // Windows fully to the left of the core cannot be hit because every counted fragment
                // starts inside the core. We store None to preserve the global index while skipping
                // both counting work and output rows for those windows
                if win_end <= tile.core_start() as u64 {
                    counts.push(None);
                    continue;
                }
                // Contained flags windows fully inside the core
                let contained =
                    win_start >= tile.core_start() as u64 && win_end <= tile.core_end() as u64;
                counts.push(Some(TileCounts {
                    counts: template.zeroed_like(),
                    contained,
                }));
            }
            (span.first_idx, span.last_idx_exclusive, counts)
        }
    };

    // Fraction of a fragment that must overlap with a window to assign to that window
    let min_overlap_fraction: f64 = match opt.window_assignment.assign_by {
        WindowAssigner::Any | WindowAssigner::CountOverlap => {
            1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0)
        } // +1 to avoid rounding error issues
        WindowAssigner::All | WindowAssigner::Midpoint => {
            1.0 - (1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0))
        } // 1.0 but just below to avoid rounding errors
        WindowAssigner::Proportion(p) => p,
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
        move |f: &FragmentWithIndelCounts| lengths.contains(f.len_indel_adjusted())
    };

    // Create fragment iterator with per-tile filtering and optional GC tag handling
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
    let mut iter = fragments_with_indel_counts_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        move |rec| include_read_fn(rec),
        opt.indel_mode,
        fragment_filter,
        unpaired,
    )
    .with_local_counters();

    let get_gc_weight = {
        let gc_corrector = gc_corrector_opt.as_ref();
        let gc_prefixes = gc_prefixes_opt.as_ref();
        let fetch_start = tile.fetch_start();
        move |fragment: &FragmentWithIndelCounts| -> Result<Option<f64>> {
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

        // Calculate fragment length
        // Note: Only "adjusted" when `--indel-mode` asks for it
        let fragment_length = fragment.len_indel_adjusted();

        // Find all overlapping count-windows

        // Calculate what part needs to overlap to some degree
        let query_interval = match opt.window_assignment.assign_by {
            WindowAssigner::Midpoint => {
                let midpoint =
                    midpoint_random_even_with_thread_rng(fragment.start(), fragment_length);
                Interval::new(midpoint.into(), (midpoint + 1).into())?
            }
            WindowAssigner::Any
            | WindowAssigner::All
            | WindowAssigner::Proportion(_)
            | WindowAssigner::CountOverlap => fragment.interval.try_to_u64()?,
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

        // Get GC correction weight
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
                1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0), // Any overlap
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
                WindowAssigner::CountOverlap => compute_window_scaling_over_overlap(
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

            // Count up the weight per overlapping count-window
            for (overlapped_window_idx, scaling_weight, overlap_fraction_to_count) in
                overlap_weights
            {
                let vec_idx = overlapped_window_idx - counts_start_idx;
                if vec_idx >= counts_by_idx.len() {
                    bail!(
                        "Overlapping window idx {} outside [{}..{}) on {}",
                        overlapped_window_idx,
                        counts_start_idx,
                        counts_end_idx_exclusive,
                        tile.chr
                    );
                }
                if let Some(entry) = counts_by_idx[vec_idx].as_mut() {
                    entry.counts.incr_weighted(
                        fragment_length as usize,
                        overlap_fraction_to_count * scaling_weight * gc_weight,
                    );
                }
            }
        } else {
            // When no scaling, increment counter by 1.0 or by the overlap fraction
            for overlapped_window in overlapping_windows.windows {
                let vec_idx = overlapped_window.idx - counts_start_idx;
                if vec_idx >= counts_by_idx.len() {
                    bail!(
                        "Overlapping window idx {} outside [{}..{}) on {}",
                        overlapped_window.idx,
                        counts_start_idx,
                        counts_end_idx_exclusive,
                        tile.chr
                    );
                }
                if let Some(entry) = counts_by_idx[vec_idx].as_mut() {
                    let count_weight = match opt.window_assignment.assign_by {
                        WindowAssigner::CountOverlap => overlapped_window.overlap_fraction as f64,
                        _ => 1.0f64,
                    };
                    entry
                        .counts
                        .incr_weighted(fragment_length as usize, count_weight * gc_weight);
                }
            }
        }
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    // Prepare outputs
    let mut window_idxs_chr: Vec<u64> = Vec::with_capacity(counts_by_idx.len());
    let mut counts: Vec<LengthCounts> = Vec::with_capacity(counts_by_idx.len());
    let mut contained_flags: Vec<bool> = Vec::with_capacity(counts_by_idx.len());
    let mut crossing_window_idxs_chr: Vec<u64> = Vec::new();
    for (offset, tile_counts_opt) in counts_by_idx.into_iter().enumerate() {
        if let Some(tile_counts) = tile_counts_opt {
            let idx = (counts_start_idx + offset) as u64;
            window_idxs_chr.push(idx);
            counts.push(tile_counts.counts);
            contained_flags.push(tile_counts.contained);
            // Aligned tile/window boundaries do not make non-contained bins tile-exclusive.
            // A fragment can start near the right edge of this core and still contribute to a
            // downstream bin that is fully contained in the next tile. Every non-contained row
            // must therefore stay in the cross-index so the reducer knows to merge it.
            if !tile_counts.contained {
                crossing_window_idxs_chr.push(idx);
            }
        }
    }

    if matches!(window_opt, WindowSpec::Global) {
        debug_assert_eq!(counts.len(), 1);
        let chr_counts = counts
            .into_iter()
            .next()
            .unwrap_or_else(|| template.zeroed_like());
        return Ok(TileOutputs {
            counters: counter,
            global_counts: Some((tile.chr.clone(), chr_counts)),
        });
    }

    let _ = write_partials_npz(
        temp_dir,
        partials_prefix,
        &tile.chr,
        tile.index,
        &window_idxs_chr,
        &contained_flags,
        &counts,
    )?;
    let _ = write_cross_npy(
        temp_dir,
        cross_prefix,
        &tile.chr,
        tile.index,
        &crossing_window_idxs_chr,
    )?;

    Ok(TileOutputs {
        counters: counter,
        global_counts: None,
    })
}
