use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use rust_htslib::bam::{self, Format, Header, Read, Record, ext::BamRecordExtensions};
use std::{sync::Arc, time::Instant};
use tracing::info;

use crate::{
    command_run::{CommandRunResult, RunOptions, status_info},
    commands::{
        bam_to_bam::{
            config::BamToBamConfig,
            sorted_writer::{RecordEntry, RecordTags, WindowSorter},
        },
        cli_common::{
            WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        counters::BamToBamCounters,
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
        bam::{bam_bai_path, build_bam_bai_index, create_chromosome_reader, open_bam_reader},
        bed::load_windows_from_bed,
        blacklist::is_blacklisted,
        fragment::with_records_fragment::WithRecordsFragment,
        fragment_iterators::fragments_with_records_from_bam,
        interval::{IndexedInterval, Interval},
        io::FinalOutputFiles,
        overlaps::find_overlapping_windows,
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq,
        scale_genome::{ScalingBin, compute_per_window_scaling_over_fragment},
        tiled_run::TempDirGuard,
        windowing::ensure_plain_bed_windows_not_empty,
    },
};

const COMMAND_TARGET: &str = "bam-to-bam";

/// Result from `bam-to-bam`.
///
/// The command writes one filtered or annotated BAM file plus its `.bam.bai` index. The result keeps
/// the fragment counters and final output paths together so library callers do not need to
/// reconstruct file names from the configuration.
#[derive(Debug)]
pub struct BamToBamRunResult {
    /// Fragment and filtering counters collected during the run.
    pub counters: BamToBamCounters,
    /// Final BAM path written by the command.
    pub output_bam: std::path::PathBuf,
    /// Final output files produced by the command: the BAM followed by its `.bam.bai` index.
    pub output_files: Vec<std::path::PathBuf>,
}

impl CommandRunResult for BamToBamRunResult {
    type Counters = BamToBamCounters;

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[std::path::PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&std::path::Path> {
        Some(&self.output_bam)
    }
}

/// Run the `bam-to-bam` command.
///
/// This is the programmatic entry point for the same filtering and annotation behavior exposed by
/// the CLI. It reads a BAM, applies the configured fragment filters, optional GC correction,
/// optional genomic scaling, and writes a new BAM.
///
/// Reporting is controlled by `options`. `report_statistics` prints the final summary,
/// `show_progress` controls progress bars, and `log_statuses` controls status messages.
///
/// Parameters
/// ----------
/// - `opt`:
///   Fully resolved configuration for the `bam-to-bam` command.
/// - `options`:
///   Reporting controls for statistics, progress bars, and status logs.
///
/// Returns
/// -------
/// - `Ok(BamToBamRunResult)`:
///   Counters and output paths for the completed run.
///
/// Errors
/// ------
/// Returns an error when the configuration is invalid, an input cannot be read, or the output BAM
/// cannot be written.
pub fn run_bam_to_bam(opt: &BamToBamConfig, options: RunOptions) -> Result<BamToBamRunResult> {
    let start_time = Instant::now();
    let global_counter = execute_bam_to_bam(opt, options)?;
    let elapsed = start_time.elapsed();
    if options.report_statistics {
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
    }
    Ok(BamToBamRunResult {
        counters: global_counter,
        output_bam: opt.out_bam.clone(),
        output_files: vec![opt.out_bam.clone(), bam_bai_path(&opt.out_bam)?],
    })
}

fn execute_bam_to_bam(opt: &BamToBamConfig, options: RunOptions) -> Result<BamToBamCounters> {
    opt.fragment_lengths.validate()?;
    opt.gc.validate(opt.ref_2bit.as_deref())?;
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    if options.log_equivalent_cli {
        let command = crate::ToCliCommand::to_cli_string(opt)?;
        let message = crate::command_run::equivalent_cli_log_message(&command);
        info!(target: COMMAND_TARGET, "{message}");
    }
    let (mut chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.in_bam.as_path())?;
    // Preserve the selected subset, but write it in the input BAM header order.
    // BAM coordinate sorting follows header order, not chromosome-name string order.
    sort_chromosomes_by_bam_header_order(&mut chromosomes, &contigs)?;
    let window_opt = opt.resolve_windows();

    // Create output directory
    let output_dir = opt
        .out_bam
        .parent()
        .expect("`--out-bam` did not contain a parent directory.");
    ensure_output_dir(output_dir)?;
    let temp_dir_guard =
        TempDirGuard::new(output_dir, COMMAND_TARGET).context("create per-run temp dir")?;
    let mut final_outputs = FinalOutputFiles::new(temp_dir_guard.path())?;

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

    // Create progress bar
    let progress = ProgressFactory::with_enabled(options.show_progress);
    let pb = Arc::new(progress.default_bar(chromosomes.len() as u64));

    status_info!(options, target: COMMAND_TARGET, "Converting per chromosome");

    pb.set_position(0);

    let header = {
        let reader = open_bam_reader(&opt.in_bam).context("opening BAM to read header")?;
        Header::from_template(reader.header())
    };
    let temp_out_bam = final_outputs.temp_path_for(&opt.out_bam)?;

    // Write the output BAM under the run temp directory while chromosome writes are ongoing
    // Move it to the requested BAM path only after the writer has closed successfully
    let mut writer = bam::Writer::from_path(&temp_out_bam, &header, Format::Bam)
        .context("creating BAM writer")?;

    let results: Vec<BamToBamCounters> = chromosomes
        .iter()
        .map(|chr| -> Result<_> {
            let out = process_chrom(
                chr,
                opt,
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
                &mut writer,
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished conversion");

    let mut global_counter = BamToBamCounters::default();

    // Collect counters
    for counter in results {
        global_counter += counter;
    }
    drop(writer);

    // HTSlib can only index a complete BAM. Build the BAI while both artifacts are still in the
    // command temp directory, then move the BAM and BAI into place together through `FinalOutputFiles`.
    let temp_out_bai = build_bam_bai_index(&temp_out_bam)?;
    let out_bai = bam_bai_path(&opt.out_bam)?;
    final_outputs.record(temp_out_bam, opt.out_bam.clone())?;
    final_outputs.record(temp_out_bai, out_bai)?;
    final_outputs.move_into_place()?;

    Ok(global_counter)
}

fn sort_chromosomes_by_bam_header_order(
    chromosomes: &mut Vec<String>,
    contigs: &crate::shared::bam::Contigs,
) -> Result<()> {
    // Chromosome selection can come from defaults, explicit CLI values, or a file. Those sources do
    // not necessarily match the BAM header order, so sort the already-selected subset by header tid.
    for chromosome in chromosomes.iter() {
        contigs.contigs.get(chromosome).with_context(|| {
            format!("missing BAM contig metadata for selected chromosome '{chromosome}'")
        })?;
    }

    chromosomes.sort_by_key(|chromosome| {
        contigs
            .contigs
            .get(chromosome)
            .map(|(target_id, _chromosome_length)| *target_id)
            .expect("target IDs were loaded for every selected chromosome")
    });
    Ok(())
}

fn process_chrom(
    chr: &str,
    opt: &BamToBamConfig,
    windows: Option<&[IndexedInterval<u64>]>,
    blacklist_intervals: &[Interval<u64>],
    coverage_scaling_chr: &[ScalingBin],
    count_scaling_chr: &[ScalingBin],
    gc_corrector_opt: Option<GCCorrector>,
    writer: &mut bam::Writer,
) -> anyhow::Result<BamToBamCounters> {
    if matches!(opt.resolve_windows(), WindowSpec::Bed(_))
        && windows.is_none_or(|window_slice| window_slice.is_empty())
    {
        return Ok(BamToBamCounters::default());
    }

    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.in_bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = BamToBamCounters::default();

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
    let (fetch_from, fetch_to) = if let Some(wn) = windows {
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
    let mut coverage_sf_ptr = 0; // Coverage scaling factor bin
    let mut count_sf_ptr = 0; // Count-based scaling factor bin

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

        let fragment_length = fragment.len();
        let tags = Arc::new(RecordTags {
            fragment_length,
            coverage_weight: coverage_weight.map(|w| w as f32),
            fragment_count_weight: count_weight.map(|w| w as f32),
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
