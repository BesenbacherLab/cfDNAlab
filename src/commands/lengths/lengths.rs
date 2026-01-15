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
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::{compute_blacklist_overlap, is_blacklisted},
        fragment::indel_counting_fragment::FragmentWithIndelCounts,
        fragment_iterator::{
            fragments_with_indel_counts_from_bam, fragments_with_indel_counts_from_single_end_bam,
        },
        io::create_text_writer,
        midpoint::midpoint_random_even_with_thread_rng,
        overlaps::find_overlapping_windows,
        read::{default_include_read_paired_end, default_include_read_single_end},
        reference::read_seq,
        scale_genome::{compute_window_scaling_over_fragment, compute_window_scaling_over_overlap},
        thread_pool::init_global_pool,
    },
};
use anyhow::{Context, Result, anyhow, bail};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{io::Write, sync::Arc, time::Instant};

/// Execute the fragment-length counting pipeline end-to-end.
///
/// Implementation details:
/// - Resolves chromosomes, prepares optional windows/blacklists/scaling data, and then processes
///   each chromosome in parallel tiles using Rayon.
/// - Streams fragments through per-window accumulators, writing `npy` arrays (and optional BED
///   metadata) summarising the length distribution per window.
/// - Applies fragment-length, blacklist, and assignment policies consistently across threads.
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
    if opt.single_end && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --single-end");
    }
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, &opt.ioc.bam.as_path())?;
    let window_opt = opt.windows.resolve_windows();
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

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

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads as usize)?;

    // Prepare per-bin counts and metadata
    let mut all_bins = Vec::new();
    let mut bin_info = Vec::new();
    let mut global_counter = LengthsCounters::default();

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let results: Vec<(
        Vec<LengthCounts>,
        Option<Vec<(String, u64, u64, u64, f64)>>,
        LengthsCounters,
    )> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _, _)> {
            let out = process_chrom(
                &chr,
                opt,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                &window_opt,
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                scaling_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                gc_corrector.clone(),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    // Collect results (in chromosome order) back into the global vectors
    for (counts_by_bin, bin_vec, counter) in results {
        all_bins.extend(counts_by_bin);
        if !matches!(window_opt, WindowSpec::Global) {
            bin_info.extend(bin_vec.unwrap());
        }
        global_counter += counter;
    }

    // Convert to single `LengthCounts` for global
    // Keep wrapped in vector to simplify writer
    let mut all_bins = if matches!(window_opt, WindowSpec::Global) {
        vec![LengthCounts::collapse(&all_bins)?]
    } else {
        all_bins
    };

    // Sort by original index (when given a bed file)
    if matches!(window_opt, WindowSpec::Bed(_)) {
        println!("Start: Reordering counts by original window index in BED file");

        // Zip into a single Vec to allow sorting together
        let mut paired: Vec<_> = bin_info.into_iter().zip(all_bins.into_iter()).collect(); // (BinInfo, DecodedCounts)

        // Sort primarily by original window index
        paired.sort_unstable_by_key(|(info, _)| info.3);

        // Unzip back out if you need separate Vecs again
        (bin_info, all_bins) = paired.into_iter().unzip();
    }

    // Write final counts to output_dir
    write_npy(
        &opt.ioc.output_dir.join("all_length_counts.npy"),
        &stack_length_counts(&all_bins),
    )
    .context("Write final fail")?;

    // Write the min+max fragment length settings
    let settings_path = opt.ioc.output_dir.join("fragment_length_settings.json");
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
    if opt.gc.gc_file.is_some() {
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
    }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.base.counted_fragments
    );
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn process_chrom(
    chr: &str,
    opt: &LengthsConfig,
    windows: Option<&[(u64, u64, u64)]>,
    window_opt: &WindowSpec,
    blacklist_intervals: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<LengthAgnosticGCCorrector>,
) -> anyhow::Result<(
    Vec<LengthCounts>,
    Option<Vec<(String, u64, u64, u64, f64)>>,
    LengthsCounters,
)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = LengthsCounters::default();

    let gc_prefixes_opt = if gc_corrector_opt.is_some() {
        let ref_2bit = opt.ref_2bit.as_ref().ok_or_else(|| {
            anyhow!("When GC correction is specified, --ref-2bit must also be specified")
        })?;
        let seq_bytes = read_seq(&ref_2bit, chr)?;
        Some(build_gc_prefixes(&seq_bytes))
    } else {
        None
    };

    let num_bins = match window_opt {
        WindowSpec::Bed(_) => windows.unwrap().len(),
        WindowSpec::Size(s) => ((chrom_len + s - 1) / s) as usize,
        WindowSpec::Global => 1,
    };

    // Initialize count arrays
    let mut counts_by_bin = vec![
        LengthCounts::new(
            opt.fragment_lengths.min_fragment_length as usize,
            opt.fragment_lengths.max_fragment_length as usize
        );
        num_bins
    ];

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

    // Replace scaling factor with unused index (for compatibility with overlap finder)
    let scaling_with_bin_idx: Vec<(u64, u64, u64)> =
        scaling_chr.iter().map(|(s, e, _)| (*s, *e, 0u64)).collect();

    // Get coordinates to fetch reads from and to
    let (fetch_from, fetch_to) = match window_opt {
        WindowSpec::Bed(_) => {
            let wn = windows.unwrap();
            let fetch_start = wn[0].0 as i64;
            let fetch_end = wn.iter().map(|w| w.1).max().unwrap() as i64;
            (
                (fetch_start - opt.fragment_lengths.max_fragment_length as i64).max(0i64),
                (fetch_end + opt.fragment_lengths.max_fragment_length as i64).min(chrom_len as i64),
            )
        }
        _ => (0i64, chrom_len as i64),
    };

    reader
        .fetch((tid, fetch_from, fetch_to))
        .context(format!("fetch {}", chr))?;

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let lengths = opt.fragment_lengths.clone();
        move |f: &FragmentWithIndelCounts| lengths.contains(f.len_indel_adjusted())
    };

    // Create fragment iterator
    let mut iter = if opt.single_end {
        let min_mapq = opt.min_mapq;
        let include_read_fn = move |r: &Record| default_include_read_single_end(r, min_mapq);
        fragments_with_indel_counts_from_single_end_bam(
            reader.records().map(|r| r.map_err(anyhow::Error::from)),
            include_read_fn,
            opt.indel_mode,
            fragment_filter,
        )
        .with_local_counters()
    } else {
        let require_proper_pair = opt.require_proper_pair;
        let min_mapq = opt.min_mapq;
        let include_read_fn =
            move |r: &Record| default_include_read_paired_end(r, require_proper_pair, min_mapq);
        fragments_with_indel_counts_from_bam(
            reader.records().map(|r| r.map_err(anyhow::Error::from)),
            include_read_fn,
            opt.indel_mode,
            fragment_filter,
        )
        .with_local_counters()
    };

    let get_gc_weight = {
        let gc_corrector = gc_corrector_opt.as_ref();
        let gc_prefixes = gc_prefixes_opt.as_ref();
        move |fragment: &FragmentWithIndelCounts| -> Result<Option<f64>> {
            match (gc_corrector, gc_prefixes) {
                (Some(corrector), Some(prefixes)) => {
                    corrector.correct_fragment(fragment.start as u64, fragment.end as u64, prefixes)
                }
                _ => Ok(None),
            }
        }
    };

    let correct_gc = opt.gc.gc_file.is_some();

    // Streaming pointers
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = 0; // Genomic window
    let mut sf_ptr = 0; // Scaling factor bin

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;
        let fragment_length = fragment.len_indel_adjusted(); // Only adjusted when --indel-mode asks for it

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

        // Find all overlapping windows
        let (interval_start, interval_end) = match opt.window_assignment.assign_by {
            WindowAssigner::Midpoint => {
                let midpoint =
                    midpoint_random_even_with_thread_rng(fragment.start, fragment_length);
                (midpoint, midpoint + 1)
            }
            WindowAssigner::Any
            | WindowAssigner::All
            | WindowAssigner::Proportion(_)
            | WindowAssigner::CountOverlap => (fragment.start, fragment.end),
        };

        // Find all overlapping count-windows
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            windows,
            opt.windows.by_size,
            interval_start.into(),
            interval_end.into(),
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
            (Some(_), false) => unreachable!(),
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
                1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0), // Any overlap
                opt.fragment_lengths.max_fragment_length.into(),
            )
            .with_context(|| format!("finding overlapping scaling bins on chr {chr}"))?
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
                    &overlapping_windows,
                    &overlapping_scaling_bin_indices,
                    scaling_chr,
                )?,
            };

            // Count up the weight per overlapping count-window
            for (overlapped_window_idx, scaling_weight, overlap_fraction_to_count) in
                overlap_weights
            {
                counts_by_bin[overlapped_window_idx].incr_weighted(
                    fragment_length as usize,
                    overlap_fraction_to_count * scaling_weight * gc_weight,
                );
            }
        } else {
            // When no scaling, increment counter by 1.0 or by the overlap fraction
            for overlapped_window in overlapping_windows.windows {
                let count_weight = match opt.window_assignment.assign_by {
                    WindowAssigner::CountOverlap => overlapped_window.overlap_fraction as f64,
                    _ => 1.0f64,
                };
                counts_by_bin[overlapped_window.idx]
                    .incr_weighted(fragment_length as usize, count_weight * gc_weight);
            }
        }
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    let bin_info = if let Some(size) = opt.windows.by_size {
        // Build bin information for chromosome
        // chrom,start,end,blacklist_overlap
        let mut bl_ptr = 0;
        let mut bin_info = Vec::with_capacity(num_bins);
        for b in 0..num_bins {
            let start = b as u64 * size;
            let end = (start + size).min(chrom_len);
            let overlap_perc =
                compute_blacklist_overlap(blacklist_intervals, start, end, 0u64, &mut bl_ptr);
            // Note: b (index) is a placeholder that is removed later
            bin_info.push((chr.to_string(), start, end, b as u64, overlap_perc));
        }
        Some(bin_info)
    } else if opt.windows.by_bed.is_some() {
        // build bin_info from the exact BED windows
        let mut bl_ptr = 0;
        let windows = windows.unwrap();
        let mut bin_info = Vec::with_capacity(num_bins);
        for (_b, (wstart, wend, original_win_idx)) in windows.iter().cloned().enumerate() {
            let overlap_perc =
                compute_blacklist_overlap(blacklist_intervals, wstart, wend, 0u64, &mut bl_ptr);
            bin_info.push((
                chr.to_string(),
                wstart,
                wend,
                original_win_idx as u64,
                overlap_perc,
            ));
        }
        Some(bin_info)
    } else {
        None
    };

    Ok((counts_by_bin, bin_info, counter))
}
