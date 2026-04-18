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
        run_statistics::{
            DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions,
            FragmentStatisticsLabels, GCStatisticsSummary, print_fragment_run_statistics,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::is_blacklisted,
        fragment::frag_file_fragment::FragFileFragment,
        fragment_iterators::fragments_with_frag_file_info_from_bam,
        interval::{IndexedInterval, Interval},
        io::dot_join,
        overlaps::find_overlapping_windows,
        progress::ProgressFactory,
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
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{fs, io::Write, path::PathBuf, sync::Arc, time::Instant};
use tracing::info;

const COMMAND_TARGET: &str = "bam-to-frag";
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
    let elapsed = start_time.elapsed();
    print_fragment_run_statistics(
        &global_counter.base,
        elapsed,
        FragmentRunStatisticsOptions {
            include_section_header: true,
            notes: &[],
            labels: FragmentStatisticsLabels {
                counted_fragments: "Fragments included",
                ..DEFAULT_FRAGMENT_STATISTICS_LABELS
            },
            blacklist_excluded_fragments: Some(global_counter.blacklisted_fragments),
            gc: opt.gc.gc_file.is_some().then_some(GCStatisticsSummary {
                neutralize_invalid_gc: opt.gc.neutralize_invalid_gc,
                failed_fragments: global_counter.gc_failed_fragments,
                missing_tags: None,
                out_of_range_tags: None,
            }),
        },
        std::iter::empty::<&str>(),
    );
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

    // Create output directory
    ensure_output_dir(&opt.ioc.output_dir)?;

    // Load blacklist intervals if provided
    if opt.blacklist.is_some() {
        info!(target: COMMAND_TARGET, "Loading blacklists");
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
            info!(target: COMMAND_TARGET, "Loading window coordinates");
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
    let coverage_scale_genome = opt.coverage_scale_genome_args();
    if coverage_scale_genome.scaling_factors.is_some() {
        info!(target: COMMAND_TARGET, "Loading coverage scaling factors");
    }
    let coverage_scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> = load_scaling_map(
        &coverage_scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(opt.gc.gc_file.is_some(), false),
    )?;
    let count_scale_genome = opt.count_scale_genome_args();
    if count_scale_genome.scaling_factors.is_some() {
        info!(target: COMMAND_TARGET, "Loading count-based scaling factors");
    }
    let count_scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> = load_scaling_map(
        &count_scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(opt.gc.gc_file.is_some(), false),
    )?;

    // Load GC correction package if specified
    if opt.gc.gc_file.is_some() {
        info!(target: COMMAND_TARGET, "Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.gc.gc_file.as_ref(),
        opt.fragment_lengths.min_fragment_length,
        opt.fragment_lengths.max_fragment_length,
    )?;

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;
    let output_file: PathBuf = opt.ioc.output_dir.join(dot_join(&[prefix, "frag.tsv.gz"]));
    let output_header_file: PathBuf = opt
        .ioc
        .output_dir
        .join(dot_join(&[prefix, "frag.header.tsv"]));

    // Create progress bar
    let progress = ProgressFactory::new();
    let pb = Arc::new(progress.default_bar(chromosomes.len() as u64));

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    info!(target: COMMAND_TARGET, "Converting per chromosome");

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
                coverage_scaling_map
                    .get(chr)
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]),
                count_scaling_map
                    .get(chr)
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]),
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
    info!(
        target: COMMAND_TARGET,
        "Concatenating chromosome-wise frag files"
    );
    concat_frag_zst_to_gzip(&chromosome_paths, &output_file, false)?;

    // Remove temporary directory once final outputs are written
    fs::remove_dir_all(&temp_dir).context("remove temp directory")?;

    // Create text line
    info!(target: COMMAND_TARGET, "Writing a header file");
    let mut header = String::from("chromosome\tstart\tend\tmin_mapq\tread1_strand");
    for extra_column in [
        opt.gc.gc_file.is_some().then_some("gc_weight"),
        opt.coverage_scaling_factors
            .is_some()
            .then_some("coverage_scaling_weight"),
        opt.count_scaling_factors
            .is_some()
            .then_some("count_scaling_weight"),
    ]
    .into_iter()
    .flatten()
    {
        header.push('\t');
        header.push_str(extra_column);
    }
    header.push('\n');

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
    windows: Option<&[IndexedInterval<u64>]>,
    blacklist_intervals: &[Interval<u64>],
    coverage_scaling_chr: &[(u64, u64, f32)],
    count_scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<GCCorrector>,
) -> anyhow::Result<(PathBuf, BamToFragCounters)> {
    let out_path = temp_dir.join(format!("{chr}.frag.tsv.zst"));

    if matches!(opt.resolve_windows(), WindowSpec::Bed(_))
        && windows.is_none_or(|window_slice| window_slice.is_empty())
    {
        let mut writer = open_zstd_auto_writer(&out_path, 3, Some(1))?;
        writer.flush()?;
        return Ok((out_path, BamToFragCounters::default()));
    }

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

    // The overlap finder only needs checked BED-like intervals here.
    //
    // In BED mode, `find_overlapping_windows(...)` uses the scan position in this ordered slice as
    // `OverlappingWindow.idx`; it does not use `IndexedInterval.idx`. Because this temporary list
    // is built in the same order as the scaling bins, those scan positions already match the
    // chromosome-local indices needed for indexing back into those bins.
    //
    // So the carried `IndexedInterval.idx` value is intentionally a placeholder.
    let coverage_scaling_with_bin_idx: Vec<IndexedInterval<u64>> = coverage_scaling_chr
        .iter()
        .map(|(start, end, _)| IndexedInterval::new(*start, *end, 0_u64))
        .collect::<crate::Result<_>>()?;
    let count_scaling_with_bin_idx: Vec<IndexedInterval<u64>> = count_scaling_chr
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
    let mut coverage_sf_ptr = 0; // Coverage scaling factor bin
    let mut count_sf_ptr = 0; // Count-based scaling factor bin

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
                if opt.gc.neutralize_invalid_gc {
                    Some(1.0)
                } else {
                    continue;
                }
            }
            (None, false) => None,
            (Some(_), false) => bail!("unexpected GC weight when GC correction is disabled"),
        };

        // Find all overlapping scaling-factor bins
        // And count up the weight
        let coverage_weight = if !coverage_scaling_chr.is_empty() {
            // Find overlapping scaling-bins
            let overlapping_scaling_bins = find_overlapping_windows(
                chrom_len,
                &mut coverage_sf_ptr,
                Some(&coverage_scaling_with_bin_idx),
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
                fragment.interval.try_to_u64()?,
                &overlapping_windows,
                &overlapping_scaling_bin_indices,
                coverage_scaling_chr,
            )?
            .pop()
            .map(|(_, w, _)| w)
            .expect("no overlapping scaling bins found");

            Some(scaling_weight)
        } else {
            None
        };
        let count_weight = if !count_scaling_chr.is_empty() {
            let overlapping_scaling_bins = find_overlapping_windows(
                chrom_len,
                &mut count_sf_ptr,
                Some(&count_scaling_with_bin_idx),
                None,
                fragment.interval.try_to_u64()?,
                1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0),
                opt.fragment_lengths.max_fragment_length.into(),
            )
            .with_context(|| format!("finding overlapping count-based scaling bins on chr {chr}"))?
            .context("no overlapping count-based scaling bins found")?;

            let overlapping_scaling_bin_indices: Vec<usize> = overlapping_scaling_bins
                .windows
                .iter()
                .map(|w| w.idx)
                .collect();

            let scaling_weight = compute_window_scaling_over_fragment(
                fragment.interval.try_to_u64()?,
                &overlapping_windows,
                &overlapping_scaling_bin_indices,
                count_scaling_chr,
            )?
            .pop()
            .map(|(_, w, _)| w)
            .expect("no overlapping count-based scaling bins found");

            Some(scaling_weight)
        } else {
            None
        };

        counter.base.counted_fragments += 1;

        // Create text line
        let mut line = format!(
            "{}\t{}\t{}\t{}\t{}",
            chr,
            fragment.start(),
            fragment.end(),
            fragment.min_mapq,
            fragment.read1_strand,
        );
        for extra_value in [gc_weight, coverage_weight, count_weight]
            .into_iter()
            .flatten()
        {
            line.push('\t');
            line.push_str(&extra_value.to_string());
        }
        line.push('\n');

        // Push into windowed sorter
        // That flushes the previous (sorted) entries on the fly
        sorter.push(
            WindowEntry {
                interval: fragment.interval,
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
