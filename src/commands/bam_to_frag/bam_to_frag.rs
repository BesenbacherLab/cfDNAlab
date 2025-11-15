use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{fs, io::Write, path::PathBuf, sync::Arc, time::Instant};

use crate::{
    commands::{
        bam_to_frag::{
            concat::concat_frag_zst_to_gzip, config::BamToFragConfig,
            sorted_writer::Entry as WindowEntry, sorted_writer::WindowSorter,
        },
        cli_common::{ensure_output_dir, load_blacklist_map, resolve_chromosomes_and_contigs},
        counters::BamToFragCounters,
    },
    shared::{
        bam::create_chromosome_reader, blacklist::is_blacklisted,
        fragment::frag_file_fragment::FragFileFragment,
        fragment_iterator::fragments_with_frag_file_info_from_bam, read::default_include_read,
        thread_pool::init_global_pool, tiled_run::make_temp_dir, writers::open_zstd_auto_writer,
    },
};

/// Execute the bam-to-frag conversion.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `bam-to-frag` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &BamToFragConfig) -> Result<()> {
    let start_time = Instant::now();
    let global_counter = run_inner(opt)?;

    let quiet = false;

    if !quiet {
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
        println!(
            "  Fragments included: {}",
            global_counter.base.counted_fragments
        );
        println!("----------");
        println!("Elapsed time: {:.2?}", elapsed);
    }
    Ok(())
}

pub fn run_inner(opt: &BamToFragConfig) -> Result<BamToFragCounters> {
    let (chromosomes, _) = resolve_chromosomes_and_contigs(&opt.chromosomes, &opt.ioc.bam.as_path())?;
    let prefix = opt.output_prefix.trim();

    let quiet = false;

    // Create output directory
    ensure_output_dir(&opt.ioc.output_dir)?;

    // Load blacklist intervals if provided
    if opt.blacklist.is_some() && !quiet {
        println!("Start: Loading blacklists");
    }
    let blacklist_map = load_blacklist_map(
        opt.blacklist.as_ref(),
        opt.blacklist_min_size,
        0,
        &chromosomes,
    )?;

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;
    let output_file: PathBuf = opt.ioc.output_dir.join(format!("{prefix}.frag.tsv.gz"));

    // Create progress bar
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    if quiet {
        pb.set_draw_target(ProgressDrawTarget::hidden());
    }

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads as usize)?;

    if !quiet {
        println!("Start: Converting per chromosome");
    }

    pb.set_position(0);

    let results: Vec<(PathBuf, BamToFragCounters)> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _)> {
            let out = process_chrom(
                &chr,
                opt,
                &temp_dir,
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished conversion");

    let mut global_counter = BamToFragCounters::default();
    let mut chromosome_paths: Vec<PathBuf> = Vec::with_capacity(chromosomes.len());

    // Collect results (in chromosome order) back into the global vectors
    for (path, counter) in results {
        global_counter += counter;
        chromosome_paths.push(path);
    }

    // Concatenate chromosome-wise temp files
    concat_frag_zst_to_gzip(&chromosome_paths, &output_file, false)?;

    // Remove temporary directory once final outputs are written
    fs::remove_dir_all(&temp_dir).context("remove temp directory")?;

    Ok(global_counter)
}

fn process_chrom(
    chr: &str,
    opt: &BamToFragConfig,
    temp_dir: &PathBuf,
    blacklist_intervals: &[(u64, u64)],
) -> anyhow::Result<(PathBuf, BamToFragCounters)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, _) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = BamToFragCounters::default();

    reader.fetch(tid).context(format!("fetch {}", chr))?;

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let lengths = opt.fragment_lengths.clone();
        move |f: &FragFileFragment| lengths.contains(f.len())
    };

    // Wrap to use opt
    let include_read_fn = {
        let opt = (*opt).clone();
        move |r: &Record| default_include_read(r, opt.require_proper_pair, opt.min_mapq)
    };

    // Create fragment iterator
    let mut iter = fragments_with_frag_file_info_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        include_read_fn,
        fragment_filter,
    )
    .with_local_counters();

    let mut bl_ptr = 0; // Blacklist interval

    let out_path = temp_dir.join(format!("{chr}.frag.tsv.zst"));

    let mut writer = open_zstd_auto_writer(&out_path, 3, Some(1))?;

    // Write using a bounded window sorter to ensure (start,end)-sorted output
    let mut sorter = WindowSorter::new(opt.fragment_lengths.max_fragment_length);

    // Iterate fragments
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

        counter.base.counted_fragments += 1;

        // Push into windowed sorter
        // That flushes the previous (sorted) entries on the fly
        let line = format!(
            "{}\t{}\t{}\t{}\t{}\n",
            chr, fragment.start, fragment.end, fragment.min_mapq, fragment.read1_strand
        );
        sorter.push(
            WindowEntry {
                start: fragment.start,
                end: fragment.end,
                line,
            },
            &mut writer,
        )?;
    }

    // Flush any fragments still buffered in the sorter tail
    sorter.flush_all(&mut writer)?;
    writer.flush()?;

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    Ok((out_path, counter))
}
