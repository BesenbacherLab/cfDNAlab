use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use rust_htslib::bam::{self, Format, Header, Read, Record, ext::BamRecordExtensions};
use std::{sync::Arc, time::Instant};

use crate::{
    commands::{
        bam_to_bam::{
            config::BamToBamConfig,
            sorted_writer::{RecordEntry, RecordTags, WindowSorter},
        },
        cli_common::{
            WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        counters::BamToFragCounters,
        gc_bias::{
            correct::{GCCorrector, load_gc_corrector},
            counting::build_gc_prefixes,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::is_blacklisted,
        fragment::with_records_fragment::WithRecordsFragment,
        fragment_iterator::fragments_with_records_from_bam,
        interval::{IndexedInterval, Interval},
        overlaps::find_overlapping_windows,
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq,
        scale_genome::compute_window_scaling_over_fragment,
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
            "  Fragments included: {}",
            global_counter.base.counted_fragments
        );
        println!("----------");
        println!("Elapsed time: {:.2?}", elapsed);
    }
    Ok(())
}

pub fn run_inner(opt: &BamToBamConfig) -> Result<BamToFragCounters> {
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let (mut chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.in_bam.as_path())?;
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
    ensure_output_dir(output_dir)?;

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

    // Load GC correction package if specified
    if opt.gc.gc_file.is_some() {
        println!("Start: Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.gc.gc_file.as_ref(),
        opt.fragment_lengths.min_fragment_length,
        opt.fragment_lengths.max_fragment_length,
    )?;

    // Create progress bar
    let progress = ProgressFactory::with_enabled(!quiet);
    let pb = Arc::new(progress.default_bar(chromosomes.len() as u64));

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
                chr,
                opt,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                scaling_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                gc_corrector.clone(),
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
    windows: Option<&[IndexedInterval<u64>]>,
    blacklist_intervals: &[Interval<u64>],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<GCCorrector>,
    writer: &mut bam::Writer,
) -> anyhow::Result<BamToFragCounters> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.in_bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = BamToFragCounters::default();

    let gc_prefixes_opt = if gc_corrector_opt.is_some() {
        let ref_2bit = match opt.ref_2bit.as_ref() {
            Some(r) => r,
            None => bail!("When GC correction is specified, --ref-2bit must also be specified"),
        };
        let seq_bytes = read_seq(ref_2bit, chr)?;
        Some(build_gc_prefixes(&seq_bytes))
    } else {
        None
    };

    // Replace scaling factor with unused index
    let scaling_with_bin_idx: Vec<IndexedInterval<u64>> = scaling_chr
        .iter()
        .map(|(start, end, _)| IndexedInterval::new(*start, *end, 0_u64))
        .collect::<crate::Result<_>>()?;

    // Get coordinates to fetch reads from and to
    let (fetch_from, fetch_to) = if windows.is_some() {
        let wn = windows.unwrap();
        let fetch_start = wn[0].start() as i64;
        let fetch_end = wn.iter().map(|window| window.end()).max().unwrap() as i64;
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

    // Create fragment iterator
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

    let mut iter = fragments_with_records_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        move |rec| include_read_fn(rec),
        fragment_filter,
        unpaired,
    )
    .with_local_counters();

    let get_gc_weight = {
        let gc_corrector = gc_corrector_opt.as_ref();
        let gc_prefixes = gc_prefixes_opt.as_ref();
        move |fragment: &WithRecordsFragment| -> Result<Option<f64>> {
            match (gc_corrector, gc_prefixes) {
                (Some(corrector), Some(prefixes)) => {
                    corrector.correct_fragment(fragment.interval.try_to_u64()?, prefixes)
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

    // Write using a bounded window sorter to ensure (start,end)-sorted output
    let mut sorter = WindowSorter::new(opt.fragment_lengths.max_fragment_length);

    // Iterate fragments
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;

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

        // Find all overlapping count-windows
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            windows,
            None,
            fragment.interval.try_to_u64()?,
            1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0), // Any overlap
            opt.fragment_lengths.max_fragment_length.into(),
        )?;
        let overlapping_windows = if let Some(overlaps) = overlapping_windows {
            overlaps
        } else {
            continue;
        };

        let gc_weight = match (get_gc_weight(&fragment)?, correct_gc) {
            (Some(w), true) => Some(w),
            (None, true) => {
                counter.gc_failed_fragments += 1;
                Some(1.0)
            }
            (None, false) => None,
            (Some(_), false) => bail!("unexpected GC weight when GC correction is disabled"),
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
                fragment.interval.try_to_u64()?, // Full fragment
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
            gc_weight: gc_weight.map(|w| w as f32),
        });

        counter.base.counted_fragments += 1;

        let WithRecordsFragment {
            single_record,
            forward_record,
            reverse_record,
            ..
        } = fragment;

        if opt.unpaired.reads_are_fragments {
            let single_record = single_record
                .expect("Single record must exist in unpaired (--reads-are-fragments) mode");
            sorter.push(
                RecordEntry {
                    interval: Interval::new(
                        single_record.pos() as u32,
                        single_record.reference_end() as u32,
                    )?,
                    record: single_record,
                    tags: tags.clone(),
                },
                writer,
            )?;
        } else {
            let forward_record =
                forward_record.expect("Forward record must exist in paired-end mode");
            let reverse_record =
                reverse_record.expect("Reverse record must exist in paired-end mode");

            sorter.push(
                RecordEntry {
                    interval: Interval::new(
                        forward_record.pos() as u32,
                        forward_record.reference_end() as u32,
                    )?,
                    record: forward_record,
                    tags: tags.clone(),
                },
                writer,
            )?;

            // Push reverse read
            sorter.push(
                RecordEntry {
                    interval: Interval::new(
                        reverse_record.pos() as u32,
                        reverse_record.reference_end() as u32,
                    )?,
                    record: reverse_record,
                    tags,
                },
                writer,
            )?;
        }
    }

    // Flush any fragments still buffered in the sorter tail
    sorter.flush_all(writer)?;

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    Ok(counter)
}
