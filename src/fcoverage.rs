use anyhow::{Context, Result};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{collections::HashMap, fs::create_dir_all, path::PathBuf, sync::Arc, time::Instant};

use crate::{
    cli_common::{ChromosomeArgs, FragmentLengthArgs, IOCArgs, WindowSpec, WindowsArgs},
    counters::FCoverageCounters,
    utils::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::load_blacklists,
        coverage::{
            coverage_prefix::CoveragePrefix,
            window_results::{CoverageOutput, CoverageWindowAction, compute_window_outputs},
            writers::{NanPolicy, write_outputs_auto_with_prefix},
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
/// Without windowing, positional coverage for the entire genome (selected chromosomes) are outputted.
///
/// ## Blacklisting
///
/// Positions in blacklisted regions are set to `f32::NaN` (and thus not included in sums or averages).
/// Set `--nan-policy` to change how these positions are handled in the output (positional coverage outputs only).
#[cfg_attr(feature = "cli", derive(clap::Args))]
pub struct FCoverageConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    windows: WindowsArgs,

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

    // Configure global thread‐pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.ioc.n_threads as usize)
        .build_global()
        .context("building Rayon thread pool")?;

    // Prepare output containers
    let mut out_by_chrom =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());

    let mut global_counter = FCoverageCounters::default();

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let results: Vec<(String, CoverageOutput, FCoverageCounters)> = chromosomes
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
    for (chr, coverage_output, counter) in results {
        out_by_chrom.insert(chr, coverage_output);
        global_counter += counter;
    }

    // Write outputs
    // decide output filename however you like
    let out_path = opt.ioc.output_dir.join("coverage");

    // Default to dropping masked rows for positional writers
    write_outputs_auto_with_prefix(&out_path, &out_by_chrom, opt.nan_policy)?;

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

fn process_chrom(
    chr: &str,
    opt: &FCoverageConfig,
    windows: Option<&[(u64, u64, u64)]>,
    window_opt: &WindowSpec,
    blacklist_intervals: &[(u64, u64)],
) -> anyhow::Result<(String, CoverageOutput, FCoverageCounters)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = FCoverageCounters::default();

    // Stash for keeping reads until mates arrive
    let mut stash = FxHashMap::<Vec<u8>, SegmentedReadInfo>::default();

    // Get coordinates to fetch reads from and to
    let (fetch_from, fetch_to) = match (window_opt, opt.per_window) {
        (WindowSpec::Bed(_), CoverageWindowAction::OnlyIncludeThesePositions) => {
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

    // Initialize coverage counter
    let mut cp = CoveragePrefix::initialize_coverage_prefix(chrom_len as u32);

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
            let fragment = if let Some(f) = collect_fragment_with_segments(
                &SegmentedReadInfo::from(&rec),
                &mate,
                1,
                !opt.ignore_gap,
            ) {
                f
            } else {
                continue;
            };

            counter.collected_fragments += 1;

            // Check length is within allowed range
            let fragment_length = fragment.len();
            if fragment_length < opt.fragment_lengths.min_fragment_length
                || fragment_length > opt.fragment_lengths.max_fragment_length
            {
                counter.illegal_length_fragments += 1;
                continue;
            }

            counter.counted_fragments += 1;

            // Add to coverage prefix
            // TODO: Update weight with GC etc. later!
            cp.add_fragment_with_segments(fragment, 1.0)?;
        } else {
            // Stash read if new qname
            stash.insert(rec.qname().to_vec(), SegmentedReadInfo::from(&rec));
        }
    }

    // Add blacklist
    if !blacklist_intervals.is_empty() {
        cp.initialize_blacklist_prefix();
        cp.add_blacklist_many_to_prefix(blacklist_intervals)?;
        cp.finalize_blacklist_prefix();
    }

    // Get ready to extract average coverage per stride-bin
    cp.finalize_coverage();
    cp.build_query_index()?;

    // Extract outputs
    let window_outputs =
        compute_window_outputs(&mut cp, windows, opt.per_window, opt.blacklist.is_some())?;

    Ok((chr.to_string(), window_outputs, counter))
}
