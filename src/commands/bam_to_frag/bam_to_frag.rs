use crate::{
    commands::{
        bam_to_frag::{
            concat::concat_frag_zst_to_gzip,
            config::BamToFragConfig,
            sorted_writer::{Entry as WindowEntry, WindowSorter},
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
        fragment::frag_file_fragment::FragFileFragment,
        fragment_iterator::fragments_with_frag_file_info_from_bam,
        overlaps::find_overlapping_windows,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq,
        scale_genome::compute_window_scaling_over_fragment,
        thread_pool::init_global_pool,
        tiled_run::make_temp_dir,
        writers::open_zstd_auto_writer,
    },
};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{fs, io::Write, path::PathBuf, sync::Arc, time::Instant};
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

pub fn run_inner(opt: &BamToFragConfig) -> Result<BamToFragCounters> {
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    let prefix = opt.output_prefix.trim();
    let window_opt = opt.resolve_windows();

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

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;
    let output_file: PathBuf = opt.ioc.output_dir.join(format!("{prefix}.frag.tsv.gz"));
    let output_header_file: PathBuf = opt.ioc.output_dir.join(format!("{prefix}.frag.header.tsv"));

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
    init_global_pool(opt.ioc.n_threads)?;

    if !quiet {
        println!("Start: Converting per chromosome");
    }

    pb.set_position(0);

    let results: Vec<(PathBuf, BamToFragCounters)> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _)> {
            let out = process_chrom(
                chr,
                opt,
                &temp_dir,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                scaling_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                gc_corrector.clone(),
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
    if !quiet {
        println!("Start: Concatenating chromosome-wise frag files");
    }
    concat_frag_zst_to_gzip(&chromosome_paths, &output_file, false)?;

    // Remove temporary directory once final outputs are written
    fs::remove_dir_all(&temp_dir).context("remove temp directory")?;

    // Create text line
    if !quiet {
        println!("Start: Writing a header file");
    }
    let header = match (
        opt.gc.gc_file.is_some(),
        opt.scale_genome.scaling_factors.is_some(),
    ) {
        (true, true) => {
            "chromosome\tstart\tend\tmin_mapq\tread1_strand\tgc_weight\tscaling_weight\n"
        }
        (true, false) => "chromosome\tstart\tend\tmin_mapq\tread1_strand\tgc_weight\n",
        (false, true) => "chromosome\tstart\tend\tmin_mapq\tread1_strand\tscaling_weight\n",
        (false, false) => "chromosome\tstart\tend\tmin_mapq\tread1_strand\n",
    };

    fs::write(&output_header_file, header).with_context(|| {
        format!(
            "Failed writing fragment header to {}",
            output_header_file.display()
        )
    })?;

    Ok(global_counter)
}

fn process_chrom(
    chr: &str,
    opt: &BamToFragConfig,
    temp_dir: &PathBuf,
    windows: Option<&[(u64, u64, u64)]>,
    blacklist_intervals: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<GCCorrector>,
) -> anyhow::Result<(PathBuf, BamToFragCounters)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = BamToFragCounters::default();

    // TODO: Consider tiling the function to decrease memory from the prefixes
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
        move |f: &FragFileFragment| lengths.contains(f.len())
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

    let mut iter = fragments_with_frag_file_info_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        move |rec| include_read_fn(rec),
        fragment_filter,
        unpaired,
    )
    .with_local_counters();

    let get_gc_weight = {
        let gc_corrector = gc_corrector_opt.as_ref();
        let gc_prefixes = gc_prefixes_opt.as_ref();
        move |fragment: &FragFileFragment| -> Result<Option<f64>> {
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
            opt.blacklist_strategy,
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

        let gc_weight = match (get_gc_weight(&fragment)?, correct_gc) {
            (Some(w), true) => Some(w),
            (None, true) => {
                counter.gc_failed_fragments += 1;
                Some(1.0)
            }
            (None, false) => None,
            (Some(_), false) => unreachable!(),
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

        counter.base.counted_fragments += 1;

        // Create text line
        let line = match (gc_weight, scaling_weight) {
            (Some(gc_w), Some(sf_w)) => format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                chr,
                fragment.start,
                fragment.end,
                fragment.min_mapq,
                fragment.read1_strand,
                gc_w,
                sf_w
            ),
            (Some(gc_w), None) => format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                chr, fragment.start, fragment.end, fragment.min_mapq, fragment.read1_strand, gc_w,
            ),
            (None, Some(sf_w)) => format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                chr, fragment.start, fragment.end, fragment.min_mapq, fragment.read1_strand, sf_w
            ),
            (None, None) => format!(
                "{}\t{}\t{}\t{}\t{}\n",
                chr, fragment.start, fragment.end, fragment.min_mapq, fragment.read1_strand,
            ),
        };

        // Push into windowed sorter
        // That flushes the previous (sorted) entries on the fly
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
