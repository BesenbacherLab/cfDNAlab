use anyhow::{Context, Result};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{
    collections::HashMap,
    fs::{File, create_dir_all},
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Arc,
    time::Instant,
};

use crate::{
    cli_common::{
        AssignToWindowArgs, ChromosomeArgs, IOCArgs, WindowAssigner, WindowSpec, WindowsArgs,
    },
    counters::LengthsCounters,
    utils::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::{BlackStrategy, compute_blacklist_overlap, is_blacklisted, load_blacklists},
        fragment::{MinimalReadInfo, collect_fragment},
        lengths::counting::{LengthCounts, stack_length_counts},
        overlaps::find_overlapping_windows,
    },
};

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[cfg_attr(
    feature = "cli",
    clap(
        group = clap::ArgGroup::new("windows")
            .required(true)
            .args(&["by_size", "by_bed"])
    )
)]
struct LengthsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    windows: WindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    window_assignment: AssignToWindowArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    chromosomes: ChromosomeArgs,

    /// Minimum fragment length to include (default: 20) [integer]
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20", value_parser = value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub min_fragment_length: u32,

    /// Maximum fragment length to include (default: 600) [integer]
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "600", value_parser = value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub max_fragment_length: u32,

    /// Minimum mapping quality to include (default: 30) [integer]
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads [flag]
    #[cfg_attr(feature = "cli", clap(long))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions [path]
    #[cfg_attr(
        feature = "cli",clap(short = 'b', long, value_parser, num_args = 1.., action = ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum size of blacklist intervals to load (bp) [integer]
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

    /// Blacklist strategy: "any", "full", "midpoint", or "proportion=<threshold>" [string]
    ///
    /// Example of proportion: `--blacklist_strategy proportion=0.2` (no space around `=`)
    #[cfg_attr(
        feature = "cli",clap(
        long,
        alias = "bl-strategy",
        value_parser = value_parser!(BlackStrategy),
        default_value = "any", help_heading = "Filtering"
    ))]
    pub blacklist_strategy: BlackStrategy,
    // #[cfg_attr(feature = "cli", clap(flatten))]
    // gc: GCArgs,

    // #[cfg_attr(feature = "cli", clap(flatten))]
    // two_bit: TwoBitArgs,
}

/// Whether to include the read or continue
fn include_read(rec: &Record, opt: &LengthsConfig) -> bool {
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

pub fn run(opt: LengthsConfig) -> Result<()> {
    let start_time = Instant::now();
    let chromosomes = opt.chromosomes.resolve_chromosomes()?;
    let window_opt = opt.windows.resolve_windows();
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    // Create output directory
    create_dir_all(&opt.ioc.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = if let Some(beds) = &opt.blacklist {
        println!("Start: Loading blacklists");
        load_blacklists(beds, opt.blacklist_min_size, &chromosomes)?
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

    // Cconfigure global thread‐pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.ioc.n_threads as usize)
        .build_global()
        .context("building Rayon thread pool")?;

    // Prepare per-bin counts and metadata
    let mut all_bins = Vec::new();
    let mut bin_info = Vec::new();
    let mut global_counter = LengthsCounters::default();

    // Main loop: process each autosome
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
                &opt,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                &window_opt,
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
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
        &opt.ioc.output_dir.join("all_counts.npy"),
        &stack_length_counts(&all_bins),
    )
    .context("Write final fail")?;

    // Write window coordinates as BED file to output_dir
    // Write bins BED file
    if !matches!(window_opt, WindowSpec::Global) {
        println!("Start: Writing window coordinates to disk");
        let mut bed_writer = BufWriter::new(
            File::create(&opt.ioc.output_dir.join("bins.bed")).context("Create bed fail")?,
        );
        for (chr, start, end, _, overlap_perc) in &bin_info {
            writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
                .context("Write bed line fail")?;
        }
    }

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
        "Blacklist-excluded fragments: {}",
        global_counter.blacklisted_fragments
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

fn process_chrom(
    chr: &str,
    opt: &LengthsConfig,
    windows: Option<&[(u64, u64, u64)]>,
    window_opt: &WindowSpec,
    blacklist_intervals: &[(u64, u64)],
) -> anyhow::Result<(
    Vec<LengthCounts>,
    Option<Vec<(String, u64, u64, u64, f64)>>,
    LengthsCounters,
)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = LengthsCounters::default();

    let num_bins = match window_opt {
        WindowSpec::Bed(_) => windows.unwrap().len(),
        WindowSpec::Size(s) => ((chrom_len + s - 1) / s) as usize,
        WindowSpec::Global => 1,
    };

    // Initialize count arrays
    let mut counts_by_bin = vec![
        LengthCounts::new(
            opt.min_fragment_length as usize,
            opt.max_fragment_length as usize
        );
        num_bins
    ];

    // Streaming pointers and single fetch for this chr
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = 0; // Genomic window

    // Stash for keeping reads until mates arrive
    let mut stash = FxHashMap::<Vec<u8>, MinimalReadInfo>::default();

    // Get coordinates to fetch reads from and to
    let (fetch_from, fetch_to) = match window_opt {
        WindowSpec::Bed(_) => {
            let wn = windows.unwrap();
            let fetch_start = wn[0].0 as i64;
            let fetch_end = wn.iter().map(|w| w.1).max().unwrap() as i64;
            (fetch_start.max(0i64), fetch_end.min(chrom_len as i64))
        }
        _ => (0i64, chrom_len as i64),
    };

    reader
        .fetch((tid, fetch_from, fetch_to))
        .context(format!("fetch {}", chr))?;

    // Loop over records and count
    for res in reader.records() {
        let rec = res.context("reading bam record")?;
        counter.total_reads += 1;

        if rec.tid() != tid as i32 || !include_read(&rec, opt) {
            continue;
        }

        match rec.is_reverse() {
            true => counter.accepted_reverse += 1,
            false => counter.accepted_forward += 1,
        }

        if let Some(mate) = stash.remove(rec.qname()) {
            // Combine reads to a fragment
            let fragment = if let Some(f) = collect_fragment(&MinimalReadInfo::from(&rec), &mate) {
                f
            } else {
                continue;
            };

            counter.collected_fragments += 1;

            // Determine blacklist status
            let in_blacklist = is_blacklisted(
                blacklist_intervals,
                opt.blacklist_strategy.clone(),
                fragment.start.into(),
                fragment.end.into(),
                opt.max_fragment_length as u64,
                &mut bl_ptr,
            );
            if in_blacklist {
                counter.blacklisted_fragments += 1;
            }

            // Check length is within allowed range
            let fragment_length = fragment.len();
            if fragment_length < opt.min_fragment_length
                || fragment_length > opt.max_fragment_length
            {
                counter.illegal_length_fragments += 1;
                continue;
            }

            // Find all overlapping windows
            let (interval_start, interval_end) = match opt.window_assignment.assign_by {
                WindowAssigner::Midpoint => {
                    let mid = fragment.start + (fragment_length / 2);
                    (mid, mid + 1)
                }
                WindowAssigner::Overlap => (fragment.start, fragment.end),
            };
            let overlapping_windows = find_overlapping_windows(
                chrom_len,
                &mut wd_ptr,
                windows,
                opt.windows.by_size,
                interval_start.into(),
                interval_end.into(),
                opt.max_fragment_length.into(),
            );
            let overlapping_windows = if let Some(overlaps) = overlapping_windows {
                overlaps
            } else {
                continue;
            };

            // Increment counter for each window / bin
            for overlapped_window in overlapping_windows.windows.iter() {
                counts_by_bin[overlapped_window.idx].incr(fragment_length as usize);
            }
        } else {
            // Stash read if new qname
            stash.insert(rec.qname().to_vec(), MinimalReadInfo::from(&rec));
        }
    }

    let bin_info = if let Some(size) = opt.windows.by_size {
        // Build bin information for chromosome
        // chrom,start,end,total_count,blacklist_overlap
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
