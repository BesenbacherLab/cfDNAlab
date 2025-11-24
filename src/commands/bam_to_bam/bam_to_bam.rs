use anyhow::{Context, Result};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rust_htslib::bam::{self, Format, Header, Read, Record, ext::BamRecordExtensions};
use std::{sync::Arc, time::Instant};

use crate::{
    commands::{
        bam_to_bam::config::BamToBamConfig,
        bam_to_bam::sorted_writer::{RecordEntry, RecordTags, WindowSorter},
        cli_common::{
            WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        counters::BamToFragCounters,
    },
    shared::{
        bam::create_chromosome_reader, bed::load_windows_from_bed, blacklist::is_blacklisted,
        fragment::with_records_fragment::WithRecordsFragment,
        fragment_iterator::fragments_with_records_from_bam, overlaps::find_overlapping_windows,
        read::default_include_read, scale_genome::compute_window_scaling_over_fragment,
    },
};

/// Execute the bam-to-bam filtering.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `bam-to-bam` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &BamToBamConfig) -> Result<()> {
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

pub fn run_inner(opt: &BamToBamConfig) -> Result<BamToFragCounters> {
    let (mut chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, &opt.in_bam.as_path())?;
    if !opt.skip_chromosome_sort {
        chromosomes.sort();
    }
    let window_opt = opt.resolve_windows();

    let quiet = false;

    // Create output directory
    let output_dir = opt
        .out_bam
        .parent()
        .expect("`--out-bam` did not contain a parent directory.");
    ensure_output_dir(&output_dir)?;

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

    if !quiet {
        println!("Start: Converting per chromosome");
    }

    pb.set_position(0);

    let header = {
        let reader = bam::Reader::from_path(&opt.in_bam).context("opening BAM to read header")?;
        Header::from_template(reader.header())
    };
    let mut writer = bam::Writer::from_path(&opt.out_bam, &header, Format::Bam)
        .context("creating BAM writer")?;

    let results: Vec<BamToFragCounters> = chromosomes
        .iter()
        .map(|chr| -> Result<_> {
            let out = process_chrom(
                &chr,
                opt,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                scaling_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                &mut writer,
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished conversion");

    let mut global_counter = BamToFragCounters::default();

    // Collect counters
    for counter in results {
        global_counter += counter;
    }

    Ok(global_counter)
}

fn process_chrom(
    chr: &str,
    opt: &BamToBamConfig,
    windows: Option<&[(u64, u64, u64)]>,
    blacklist_intervals: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    writer: &mut bam::Writer,
) -> anyhow::Result<BamToFragCounters> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.in_bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = BamToFragCounters::default();

    // Replace scaling factor with unused index
    let scaling_with_bin_idx: Vec<(u64, u64, u64)> =
        scaling_chr.iter().map(|(s, e, _)| (*s, *e, 0u64)).collect();

    // Get coordinates to fetch reads from and to
    let (fetch_from, fetch_to) = if windows.is_some() {
        let wn = windows.unwrap();
        let fetch_start = wn[0].0 as i64;
        let fetch_end = wn.iter().map(|w| w.1).max().unwrap() as i64;
        (
            (fetch_start - opt.fragment_lengths.max_fragment_length as i64).max(0i64),
            (fetch_end + opt.fragment_lengths.max_fragment_length as i64).min(chrom_len as i64),
        )
    } else {
        (0i64, chrom_len as i64)
    };

    reader
        .fetch((tid, fetch_from, fetch_to))
        .context(format!("fetch {}", chr))?;

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let lengths = opt.fragment_lengths.clone();
        move |f: &WithRecordsFragment| lengths.contains(f.len())
    };

    // Wrap to use opt
    let include_read_fn = {
        let opt = (*opt).clone();
        move |r: &Record| default_include_read(r, opt.require_proper_pair, opt.min_mapq)
    };

    // Create fragment iterator
    let mut iter = fragments_with_records_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        include_read_fn,
        fragment_filter,
    )
    .with_local_counters();

    // Streaming pointers
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = 0; // Genomic window
    let mut sf_ptr = 0; // Scaling factor bin

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

        // Find all overlapping count-windows
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            windows,
            None,
            fragment.start.into(),
            fragment.end.into(),
            1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0), // Any overlap
            opt.fragment_lengths.max_fragment_length.into(),
        )?;
        let overlapping_windows = if let Some(overlaps) = overlapping_windows {
            overlaps
        } else {
            continue;
        };

        // Find all overlapping scaling-factor bins
        // And count up the weight
        let scaling_weight = if !scaling_chr.is_empty() {
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
            let scaling_weight = compute_window_scaling_over_fragment(
                &overlapping_windows,
                &overlapping_scaling_bin_indices,
                scaling_chr,
            )?
            .pop()
            .map(|(_, w, _)| w)
            .expect("no overlapping scaling bins found");

            Some(scaling_weight)
        } else {
            None
        };

        let fragment_length = fragment.len();
        let tags = Arc::new(RecordTags {
            fragment_length,
            coverage_weight: scaling_weight.map(|w| w as f32),
            gc_weight: None,
        });

        counter.base.counted_fragments += 1;

        let WithRecordsFragment {
            forward_record,
            reverse_record,
            ..
        } = fragment;

        sorter.push(
            RecordEntry {
                start: forward_record.pos() as u32,
                end: forward_record.reference_end() as u32,
                record: forward_record,
                tags: tags.clone(),
            },
            writer,
        )?;

        // Push reverse read
        sorter.push(
            RecordEntry {
                start: reverse_record.pos() as u32,
                end: reverse_record.reference_end() as u32,
                record: reverse_record,
                tags,
            },
            writer,
        )?;
    }

    // Flush any fragments still buffered in the sorter tail
    sorter.flush_all(writer)?;

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    Ok(counter)
}
