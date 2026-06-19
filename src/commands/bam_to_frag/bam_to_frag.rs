use crate::{
    command_run::{CommandRunResult, RunOptions, status_info},
    commands::{
        bam_to_frag::{
            concat::concat_frag_zst_to_gzip,
            config::BamToFragConfig,
            sorted_writer::{Entry as WindowEntry, WindowSorter},
        },
        cli_common::{
            WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs, validate_output_prefix,
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
        io::{FinalOutputFiles, dot_join},
        overlaps::find_overlapping_windows,
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq,
        scale_genome::{ScalingBin, compute_per_window_scaling_over_fragment},
        temp_chrom_names::TempChromNameMap,
        thread_pool::init_global_pool,
        tiled_run::TempDirGuard,
        windowing::ensure_plain_bed_windows_not_empty,
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

/// Result from `bam-to-frag`.
///
/// The command writes a fragment table and a matching header file. The result also exposes the
/// counters used for fragment filtering and output statistics.
#[derive(Debug)]
pub struct BamToFragRunResult {
    /// Fragment and filtering counters collected during the run.
    pub counters: BamToFragCounters,
    /// Final fragment table path written by the command.
    pub output_frag: PathBuf,
    /// Header path that describes the fragment table columns.
    pub output_header: PathBuf,
    /// Final output files produced by the command.
    pub output_files: Vec<PathBuf>,
}

impl CommandRunResult for BamToFragRunResult {
    type Counters = BamToFragCounters;

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&std::path::Path> {
        Some(&self.output_frag)
    }
}

/// Run the `bam-to-frag` command.
///
/// This is the programmatic entry point for converting BAM records into the fragment table format
/// used by downstream cfDNAlab commands. It applies the configured fragment filters, optional GC
/// correction, optional genomic scaling, and writes the table plus header artifacts.
///
/// Reporting is controlled by `options`. `report_statistics` prints the final summary,
/// `show_progress` controls progress bars, and `log_statuses` controls status messages.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `bam-to-frag` command.
/// - `options`:
///     Reporting controls for statistics, progress bars, and status logs.
///
/// Returns
/// -------
/// - `Ok(BamToFragRunResult)`:
///     Counters and output paths for the completed run.
///
/// Errors
/// ------
/// Returns an error when the configuration is invalid, an input cannot be read, or any output file
/// cannot be written.
pub fn run_bam_to_frag(opt: &BamToFragConfig, options: RunOptions) -> Result<BamToFragRunResult> {
    let start_time = Instant::now();
    let run_result = execute_bam_to_frag(opt, options)?;
    let elapsed = start_time.elapsed();
    if options.report_statistics {
        print_fragment_run_statistics(
            &run_result.counters.base,
            elapsed,
            FragmentRunStatisticsOptions {
                include_section_header: true,
                notes: &[],
                labels: FragmentStatisticsLabels {
                    counted_fragments: "Fragments included",
                    ..DEFAULT_FRAGMENT_STATISTICS_LABELS
                },
                blacklist_excluded_fragments: Some(run_result.counters.blacklisted_fragments),
                gc: opt.gc.gc_file.is_some().then_some(GCStatisticsSummary {
                    neutralize_invalid_gc: opt.gc.neutralize_invalid_gc,
                    failed_fragments: run_result.counters.gc_failed_fragments,
                    missing_tags: None,
                    out_of_range_tags: None,
                }),
            },
            std::iter::empty::<&str>(),
        );
    }
    Ok(run_result)
}

fn execute_bam_to_frag(opt: &BamToFragConfig, options: RunOptions) -> Result<BamToFragRunResult> {
    opt.fragment_lengths.validate()?;
    opt.gc.validate(opt.ref_2bit.as_deref())?;
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let prefix = opt.output_prefix.trim();
    validate_output_prefix(prefix)?;
    let window_opt = opt.resolve_windows();
    if options.log_equivalent_cli {
        let command = crate::ToCliCommand::to_cli_string(opt)?;
        let message = crate::command_run::equivalent_cli_log_message(&command);
        info!(target: COMMAND_TARGET, "{message}");
    }
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    let temp_chrom_name_map = TempChromNameMap::from_contigs(&chromosomes)?;

    // Create output directory
    ensure_output_dir(&opt.ioc.output_dir)?;

    // Load blacklist intervals if provided
    if opt.blacklist.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Loading blacklists");
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
            status_info!(options, target: COMMAND_TARGET, "Loading window coordinates");
            let windows = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            ensure_plain_bed_windows_not_empty(&windows)?;
            Some(windows)
        }
        _ => None,
    };

    // Load genomic scaling factors
    let coverage_scale_genome = opt.coverage_scale_genome_args();
    if coverage_scale_genome.scaling_factors.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Loading coverage scaling factors");
    }
    let coverage_scaling_map: FxHashMap<String, Vec<ScalingBin>> = load_scaling_map(
        &coverage_scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(opt.gc.gc_file.is_some(), false),
        None,
    )?;
    let count_scale_genome = opt.count_scale_genome_args();
    if count_scale_genome.scaling_factors.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Loading count-based scaling factors");
    }
    let count_scaling_map: FxHashMap<String, Vec<ScalingBin>> = load_scaling_map(
        &count_scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(opt.gc.gc_file.is_some(), false),
        None,
    )?;

    // Load GC correction package if specified
    if opt.gc.gc_file.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.gc.gc_file.as_ref(),
        opt.ref_2bit.as_ref(),
        opt.fragment_lengths.min_fragment_length,
        opt.fragment_lengths.max_fragment_length,
    )?;

    // Build temporary directory
    let temp_dir_guard =
        TempDirGuard::new(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;
    let temp_dir = temp_dir_guard.path().to_path_buf();
    let mut final_outputs = FinalOutputFiles::new(temp_dir_guard.path())?;
    let output_file: PathBuf = opt.ioc.output_dir.join(dot_join(&[prefix, "frag.tsv.gz"]));
    let output_header_file: PathBuf = opt
        .ioc
        .output_dir
        .join(dot_join(&[prefix, "frag.header.tsv"]));

    // Create progress bar
    let progress = ProgressFactory::with_enabled(options.show_progress);
    let pb = Arc::new(progress.default_bar(chromosomes.len() as u64));

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    status_info!(options, target: COMMAND_TARGET, "Converting per chromosome");

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
                &temp_chrom_name_map,
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
    status_info!(
        options,
        target: COMMAND_TARGET,
        "Concatenating chromosome-wise frag files"
    );
    let temp_output_file = final_outputs.temp_path_for(&output_file)?;
    concat_frag_zst_to_gzip(&chromosome_paths, &temp_output_file, false)?;
    final_outputs.record(temp_output_file, output_file.clone())?;

    // Create text line
    status_info!(options, target: COMMAND_TARGET, "Writing a header file");
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

    let temp_output_header_file = final_outputs.temp_path_for(&output_header_file)?;
    fs::write(&temp_output_header_file, header).with_context(|| {
        format!(
            "Failed writing fragment header to {}",
            temp_output_header_file.display()
        )
    })?;

    final_outputs.record(temp_output_header_file, output_header_file.clone())?;

    // Keep the final fragment file and header hidden while either write can still fail
    // Move both completed files into output_dir together after the data and header are ready
    final_outputs.move_into_place()?;

    Ok(BamToFragRunResult {
        counters: global_counter,
        output_frag: output_file.clone(),
        output_header: output_header_file.clone(),
        output_files: vec![output_file, output_header_file],
    })
}

fn process_chrom(
    chr: &str,
    opt: &BamToFragConfig,
    temp_dir: &PathBuf,
    windows: Option<&[IndexedInterval<u64>]>,
    blacklist_intervals: &[Interval<u64>],
    coverage_scaling_chr: &[ScalingBin],
    count_scaling_chr: &[ScalingBin],
    gc_corrector_opt: Option<GCCorrector>,
    temp_chrom_name_map: &TempChromNameMap,
) -> anyhow::Result<(PathBuf, BamToFragCounters)> {
    let out_path = temp_chrom_name_map.path_with_suffix(temp_dir, chr, "frag.tsv.zst")?;

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
        .map(|b| IndexedInterval::from_interval(b.interval, 0_u64))
        .collect();
    let count_scaling_with_bin_idx: Vec<IndexedInterval<u64>> = count_scaling_chr
        .iter()
        .map(|b| IndexedInterval::from_interval(b.interval, 0_u64))
        .collect();

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
            // NOTE: `compute_per_window_scaling_over_fragment` always returns
            // an overlap fraction of 1.0 (count full fragment)!
            let scaling_weight = compute_per_window_scaling_over_fragment(
                fragment.interval.try_to_u64()?,
                &overlapping_windows,
                &overlapping_scaling_bin_indices,
                coverage_scaling_chr,
            )?
            .pop()
            .map(|window_scaling| window_scaling.scaling_weight)
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

            let scaling_weight = compute_per_window_scaling_over_fragment(
                fragment.interval.try_to_u64()?,
                &overlapping_windows,
                &overlapping_scaling_bin_indices,
                count_scaling_chr,
            )?
            .pop()
            .map(|window_scaling| window_scaling.scaling_weight)
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
