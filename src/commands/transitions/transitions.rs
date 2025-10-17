use crate::{
    commands::{
        cli_common::{
            WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        counters::FragmentKmersCounters,
        fragment_kmers::{config::*, tiling::*, windows::*},
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::{apply_blacklist_mask_to_seq, apply_mask::BLACKLIST_BYTE, is_blacklisted},
        fragment::segment_kmer_fragment::FragmentWithKmerSegments,
        fragment_iterator::fragments_with_kmer_segments_from_bam,
        io::create_text_writer,
        kmers::{
            kmer_codec::{
                Kmer, KmerCodes, KmerSpec, build_kmer_specs, build_left_aligned_codes_per_k,
            },
            process_counts::prepare_decoded_counts,
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
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{convert::TryInto, io::Write, path::Path, sync::Arc, time::Instant};

/// Execute the base transition probability counting pipeline end-to-end.
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
/// - `opt`: Fully resolved configuration for the `transitions` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &FragmentKmersConfig) -> Result<()> {
    let start_time = Instant::now();
    let (chromosomes, contigs) = resolve_chromosomes_and_contigs(&opt.chromosomes, &opt.ioc)?;
    let window_opt = opt.windows.resolve_windows();
    let prefix = opt.output_prefix.trim();

    // Create output directory
    ensure_output_dir(&opt.ioc.output_dir)?;

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;

    // TODO: Run fragment-kmers and get counts - pass all relevant args
    // TODO: output should be a temporary directory that we create here first!
    FragmentKmersConfig::new(ioc, ref_genome, chromosomes);
    // TODO: Refactor fragment-kmers so we can reuse without printing stats twice - e.g. make a version that returns the counter and paths to outputs! (that fragment_kmers::run wraps)

    println!("Start: Reducing per-tile counts");

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

    let all_bins = merge_tile_counts(payloads, total_windows_usize, &kmer_specs)?;

    // Prepare counts to get correct motifs (collapsed, N-filtered, etc.)
    let (prepared_counts, motifs_by_k) =
        prepare_decoded_counts(&all_bins, opt.canonical, &kmer_specs);

    // Build bin metadata when windowed
    let bin_info = build_bin_info(
        &window_opt,
        &chromosomes,
        &contigs,
        windows_map.as_ref(),
        &blacklist_map,
        chr_offsets.as_ref(),
    )?;

    // Write final counts to output_dir
    println!("Start: Writing counts to disk");
    write_decoded_counts_matrix(
        &prepared_counts,
        &kmer_specs,
        &motifs_by_k,
        &opt.ioc.output_dir,
        &opt.output_prefix,
        opt.save_sparse,
    )?;

    // Write window coordinates as BED file to output_dir
    // Write bins BED file
    if !matches!(window_opt, WindowSpec::Global) {
        println!("Start: Writing window coordinates to disk");
        let bins_path = opt.ioc.output_dir.join("bins.bed");
        let mut bed_writer = create_text_writer(&bins_path).context("Create bed fail")?;
        for (chr, start, end, _, overlap_perc) in &bin_info {
            writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
                .context("Write bed line fail")?;
        }
        bed_writer.finish().context("Finalize bins.bed writer")?;
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
    Ok(())
}

/// Process a single tile: stream fragments, accumulate per-window counts, and persist results.
fn process_tile(
    opt: &FragmentKmersConfig,
    tile: &Tile,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    window_ctx: &WindowContext,
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_intervals: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    counts_path: &Path,
) -> anyhow::Result<TileResult> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;

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
        &opt.ref_genome.ref_2bit,
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
    let mut counts_by_window: FxHashMap<u64, FxHashMap<Kmer, f64>> = FxHashMap::default();

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
        let lengths = opt.fragment_lengths.clone();
        move |f: &FragmentWithKmerSegments| lengths.contains(f.len())
    };

    // Wrap to use opt
    let include_read_fn = {
        let opt = (*opt).clone();
        move |r: &Record| default_include_read(r, opt.require_proper_pair, opt.min_mapq)
    };

    // Create fragment iterator
    let mut iter = fragments_with_kmer_segments_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        include_read_fn,
        opt.indel_mode,
        !opt.ignore_gap,
        opt.end_offset,
        fragment_filter,
    )
    .with_local_counters();

    // Initialize counters (default -> 0s)
    let mut counter = FragmentKmersCounters::default();

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;

        // Determine blacklist status
        let in_blacklist = is_blacklisted(
            blacklist_intervals,
            opt.blacklist_strategy.clone(),
            fragment.start.into(),
            fragment.end.into(),
            opt.fragment_lengths.max_fragment_length as u64,
            &mut bl_ptr,
        );
        if in_blacklist {
            counter.blacklisted_fragments += 1;
            continue;
        }

        // Find all overlapping count-windows
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            window_ctx.windows_slice(),
            opt.windows.by_size,
            (fragment.start + opt.end_offset).into(), // Should only get fragments where this is okay
            (fragment.end - opt.end_offset).into(),
            1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0), // Any overlap
            (opt.fragment_lengths.max_fragment_length + opt.end_offset).into(),
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
            count_kmers_in_segments_clipped(
                &fragment,
                &positional_codes_by_k,
                kmer_specs,
                counts,
                positional_scaling_weights.as_deref(),
                tile.core_start,
                tile.core_end,
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
            for (kmer, value) in hm {
                entries.push(TileKmerCountEntry {
                    k: kmer.k,
                    code: kmer.code,
                    value,
                });
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

/// Count kmers within the fragment’s usable segments, respecting tile core boundaries.
fn count_kmers_in_segments_clipped(
    fragment: &FragmentWithKmerSegments,
    positional_codes_by_k: &FxHashMap<u8, KmerCodes>,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    counts: &mut FxHashMap<Kmer, f64>,
    weights: Option<&[f32]>,
    tile_core_start: u32,
    tile_core_end: u32,
) {
    for (&k, _) in kmer_specs {
        let codes = positional_codes_by_k
            .get(&k)
            .expect("missing positional codes for requested k");
        let k_span = k as u32;

        for &(seg_start, seg_end) in &fragment.segments {
            let seg_start = seg_start.max(tile_core_start);

            // Allow k-mers that START within the core to use up to (k-1) bases past core_end
            let effective_seg_end =
                seg_end.min(tile_core_end.saturating_add(k_span).saturating_sub(1));

            if seg_start >= seg_end.min(tile_core_end) {
                continue;
            }

            let Some(last_start) = effective_seg_end.checked_sub(k_span) else {
                continue;
            };
            if last_start < seg_start {
                continue;
            }

            for idx_abs in seg_start..=last_start {
                // Count only starts inside the core to avoid double counting across tiles
                if idx_abs >= tile_core_end {
                    break;
                }
                let idx_local = (idx_abs - tile_core_start) as usize;
                let w = weights.map_or(1.0, |weights| unsafe { *weights.get_unchecked(idx_local) });

                *counts
                    .entry(Kmer {
                        k,
                        code: codes.get(idx_local),
                    })
                    .or_insert(0.) += w as f64;
            }
        }
    }
}
