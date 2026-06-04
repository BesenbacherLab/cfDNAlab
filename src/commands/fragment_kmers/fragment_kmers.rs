//! Runner for fragment-kmers extraction from a BAM file.
//!
//! The intended positional selection logic is specified in the
//! `positional_selection_logic.md` document.

use crate::shared::gc_tag::ClassifiedGCTagWeight;
use crate::{
    command_run::{CommandRunResult, RunOptions},
    commands::{
        cli_common::{
            WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs, validate_output_prefix,
        },
        counters::FragmentKmersCounters,
        fragment_kmers::{
            config::*,
            nearest_frame_guard::NearestFrameGuard,
            parse::*,
            positional_output::*,
            positions::*,
            selection::{SelectionDecision, evaluate_selection},
            tiling::*,
        },
        gc_bias::{
            correct::{GCCorrector, load_gc_corrector},
            counting::build_gc_prefixes,
        },
        run_statistics::{
            DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions, GCStatisticsSummary,
            TILE_DOUBLE_COUNT_NOTE, print_fragment_run_statistics,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::{apply_blacklist_mask_to_seq, apply_mask::BLACKLIST_BYTE, is_blacklisted},
        fragment::segment_kmer_fragment::FragmentWithKmerSegments,
        fragment_iterators::fragments_with_kmer_segments_from_bam,
        interval::Interval,
        io::{FinalOutputFiles, dot_join},
        kmers::{
            kmer_codec::{KmerCodes, KmerSpec, build_kmer_specs, build_left_aligned_codes_per_k},
            process_counts::{DecodedCounts, prepare_decoded_counts, split_and_decode_counts},
            write::write_decoded_counts_matrix,
        },
        overlaps::find_overlapping_windows,
        positioning::{BasesFrom, ReferenceFrame},
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq_in_range,
        scale_genome::{ScalingBin, apply_scaling_to_coverage_in_place},
        temp_chrom_names::TempChromNameMap,
        thread_pool::init_global_pool,
        tiled_run::{
            TempDirGuard, Tile, TileWindowSpan, build_tiles, precompute_tile_window_spans,
        },
        window_fetch::{BedFetchPolicy, fetch_span_for_tile},
        windowing::{
            WindowContext, build_bin_info, compute_window_offsets,
            ensure_plain_bed_windows_not_empty, write_bin_info_tsv,
        },
    },
};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{convert::TryInto, path::Path, sync::Arc, time::Instant};
use tracing::info;

const COMMAND_TARGET: &str = "fragment-kmers";

/// Result from `fragment-kmers`.
///
/// The command writes k-mer count matrices and motif metadata for the selected fragment positions.
/// The result records all final artifacts and the counters from fragment processing.
#[derive(Debug)]
pub struct FragmentKmersRunResult {
    /// Fragment and filtering counters collected during the run.
    pub counters: FragmentKmersCounters,
    /// Final output files produced by the command.
    pub output_files: Vec<std::path::PathBuf>,
}

impl CommandRunResult for FragmentKmersRunResult {
    type Counters = FragmentKmersCounters;

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[std::path::PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&std::path::Path> {
        self.output_files.first().map(std::path::PathBuf::as_path)
    }
}

/// Run the `fragment-kmers` command.
///
/// This command counts k-mers at configured fragment-relative positions. It resolves chromosomes,
/// prepares optional windows, blacklists, and scaling data, then processes tiles in parallel and
/// writes count matrices plus motif metadata.
///
/// Reporting is controlled by `options`. `report_statistics` prints the final summary,
/// `show_progress` controls progress bars, and `log_statuses` controls status messages.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `fragment-kmers` command.
/// - `options`:
///     Reporting controls for statistics, progress bars, and status logs.
///
/// Returns
/// -------
/// - `Ok(FragmentKmersRunResult)`:
///     Counters and output paths for the completed run.
///
/// Errors
/// ------
/// Returns an error when the configuration is invalid, an input cannot be read, or any output file
/// cannot be written.
pub fn run_fragment_kmers(
    opt: &FragmentKmersConfig,
    options: RunOptions,
) -> Result<FragmentKmersRunResult> {
    let start_time = Instant::now();
    opt.shared_args.fragment_lengths.validate()?;
    opt.shared_args
        .gc
        .validate(Some(opt.shared_args.ref_genome.ref_2bit.as_path()))?;
    validate_output_prefix(opt.shared_args.output_prefix.trim())?;
    let run_result = execute_fragment_kmers(opt, options)?;
    let global_counter = run_result.counters;
    let elapsed = start_time.elapsed();
    if options.report_statistics {
        print_fragment_run_statistics(
            &global_counter.base,
            elapsed,
            FragmentRunStatisticsOptions {
                include_section_header: true,
                notes: &[TILE_DOUBLE_COUNT_NOTE],
                labels: DEFAULT_FRAGMENT_STATISTICS_LABELS,
                blacklist_excluded_fragments: Some(global_counter.blacklisted_fragments),
                gc: (opt.shared_args.gc.gc_file.is_some() || opt.shared_args.gc.gc_tag.is_some())
                    .then_some(GCStatisticsSummary {
                        neutralize_invalid_gc: opt.shared_args.gc.neutralize_invalid_gc,
                        failed_fragments: global_counter.gc_failed_fragments,
                        missing_tags: opt
                            .shared_args
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(global_counter.gc_missing_tags),
                        out_of_range_tags: opt
                            .shared_args
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(global_counter.gc_out_of_range_tags),
                    }),
            },
            std::iter::empty::<&str>(),
        );
    }
    Ok(run_result)
}

fn execute_fragment_kmers(
    opt: &FragmentKmersConfig,
    options: RunOptions,
) -> Result<FragmentKmersRunResult> {
    opt.shared_args.fragment_lengths.validate()?;
    opt.shared_args
        .gc
        .validate(Some(opt.shared_args.ref_genome.ref_2bit.as_path()))?;
    if opt.shared_args.unpaired.reads_are_fragments && opt.shared_args.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let (chromosomes, contigs) = resolve_chromosomes_and_contigs(
        &opt.shared_args.chromosomes,
        opt.shared_args.ioc.bam.as_path(),
    )?;
    let window_opt = opt.shared_args.windows.resolve_windows();
    let position_specs = opt
        .shared_args
        .position_selection
        .clone()
        .into_positional_specs()?;
    let prefix = opt.shared_args.output_prefix.trim();
    validate_output_prefix(prefix)?;

    // Create output directory
    ensure_output_dir(&opt.shared_args.ioc.output_dir)?;

    // Load blacklist intervals if provided
    if opt.shared_args.blacklist.is_some() && options.log_statuses {
        info!(target: COMMAND_TARGET, "Loading blacklists");
    }
    let blacklist_map = load_blacklist_map(
        opt.shared_args.blacklist.as_ref(),
        opt.shared_args.blacklist_min_size,
        0,
        &chromosomes,
    )?;

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            if options.log_statuses {
                info!(target: COMMAND_TARGET, "Loading window coordinates");
            }
            let windows = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            ensure_plain_bed_windows_not_empty(&windows)?;
            Some(windows)
        }
        _ => None,
    };

    let kmer_specs: FxHashMap<u8, KmerSpec> = build_kmer_specs(&opt.kmer_sizes)?;

    let positional_cache = {
        if opt.shared_args.base_selection.bases_from != BasesFrom::Reference {
            bail!("position selection currently supports bases-from=reference only");
        }
        // Parse each positions specification
        let position_specs = position_specs
            .iter()
            .map(parse_positions)
            .collect::<Result<Vec<_>, _>>()?;

        let kmer_sizes: Vec<u8> = kmer_specs.keys().cloned().collect();

        Arc::new(PositionSelectionCache::new(
            position_specs,
            &kmer_sizes,
            opt.shared_args.fragment_lengths.min_fragment_length,
            opt.shared_args.fragment_lengths.max_fragment_length,
        )?)
    };

    // Load genomic scaling factors
    if opt.shared_args.scale_genome.scaling_factors.is_some() && options.log_statuses {
        info!(target: COMMAND_TARGET, "Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<ScalingBin>> = load_scaling_map(
        &opt.shared_args.scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(
            opt.shared_args.gc.gc_file.is_some(),
            opt.shared_args.gc.gc_tag.is_some(),
        ),
        Some(opt.shared_args.ignore_gap),
    )?;

    // Load GC correction package if specified
    if opt.shared_args.gc.gc_file.is_some() && options.log_statuses {
        info!(target: COMMAND_TARGET, "Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.shared_args.gc.gc_file.as_ref(),
        Some(&opt.shared_args.ref_genome.ref_2bit),
        opt.shared_args.fragment_lengths.min_fragment_length,
        opt.shared_args.fragment_lengths.max_fragment_length,
    )?;

    // Build temporary directory
    let temp_dir_guard = TempDirGuard::new(&opt.shared_args.ioc.output_dir, prefix)
        .context("create per-run temp dir")?;
    let mut final_outputs = FinalOutputFiles::new(temp_dir_guard.path())?;

    // Window size when --by-size (otherwise None)
    let by_size_bp: Option<u64> = match &window_opt {
        WindowSpec::Size(bp) => Some(*bp),
        _ => None,
    };

    // Build tiles
    let halo_bp: u32 = opt.shared_args.fragment_lengths.max_fragment_length; // Safe halo for pairing/segments
    let (tiles, _tile_and_window_boundaries_align) = build_tiles(
        &chromosomes,
        &contigs,
        opt.shared_args.tile_size,
        halo_bp,
        by_size_bp,
    )?;
    let temp_chrom_name_map = TempChromNameMap::from_contigs(&chromosomes)?;

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(
        &tiles,
        |chr| {
            windows_lookup
                .and_then(|m| m.get(chr))
                .map(|w| w.as_slice())
                .unwrap_or(&[])
        },
        0,
        0,
    ));

    // Compute per-chromosome window offsets and overall window count. In BED mode these offsets are
    // zero because windows already carry their global `original_idx` values.
    let (total_windows, chr_offsets_map) =
        compute_window_offsets(&window_opt, &chromosomes, &contigs, windows_map.as_ref())?;
    let chr_offsets = Arc::new(chr_offsets_map);

    let total_tiles = tiles.len();
    let temp_dir = Arc::new(temp_dir_guard.path().to_path_buf());
    let mut output_files = Vec::new();

    // Create progress bar
    let progress = ProgressFactory::with_enabled(options.show_progress);
    let pb = Arc::new(progress.default_bar(total_tiles as u64));

    // Configure global thread‐pool size
    init_global_pool(opt.shared_args.ioc.n_threads)?;

    if options.log_statuses {
        info!(target: COMMAND_TARGET, "Counting per chromosome");
    }

    pb.set_position(0);

    let tile_window_spans_for_threads = tile_window_spans.clone();
    let positional_cache_for_threads = positional_cache.clone();

    let tile_results: Vec<TileResult> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<TileResult> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let chr_token = temp_chrom_name_map.token_for(tile.chr.as_str())?;
            let counts_path = temp_dir.join(format!(
                "{prefix}.{chr}.{idx}.counts.bin",
                prefix = prefix,
                chr = chr_token,
                idx = tile.index
            ));

            let window_ctx = WindowContext {
                spec: &window_opt,
                windows: windows_map
                    .as_ref()
                    .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice())),
                chr_idx_offset: *chr_offsets.get(&tile.chr).unwrap_or(&0),
            };

            let blacklist_chr = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            let position_cache = positional_cache_for_threads.clone();
            let out = process_tile(
                opt,
                tile,
                &kmer_specs,
                position_cache,
                &window_ctx,
                tile_span.as_ref(),
                blacklist_chr,
                scaling_chr,
                gc_corrector.clone(),
                counts_path.as_path(),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    if options.show_progress {
        pb.finish_with_message("| Finished counting");
    } else {
        pb.finish_and_clear();
    }

    if options.log_statuses {
        info!(target: COMMAND_TARGET, "Reducing per-tile counts");
    }

    let mut global_counter = FragmentKmersCounters::default();
    let mut tile_results_by_chr: FxHashMap<String, Vec<TileResult>> = FxHashMap::default();

    for tile_result in tile_results {
        global_counter += tile_result.counter;
        tile_results_by_chr
            .entry(tile_result.chr.clone())
            .or_default()
            .push(tile_result);
    }

    let mut tile_count_batches: Vec<Vec<TileWindowCounts>> =
        Vec::with_capacity(tile_results_by_chr.len());
    for chr in &chromosomes {
        if let Some(chr_tile_results) = tile_results_by_chr.remove(chr) {
            tile_count_batches.push(reduce_chromosome_tile_results(chr_tile_results)?);
        }
    }
    if !tile_results_by_chr.is_empty() {
        let unexpected_chr = tile_results_by_chr.keys().next().unwrap();
        bail!(
            "tile results produced for unexpected chromosome '{}'",
            unexpected_chr
        );
    }

    let total_windows_usize: usize = total_windows
        .try_into()
        .context("number of windows exceeds addressable size")?;

    if opt.positional_counts {
        let positional_bins =
            merge_tile_counts_positional(tile_count_batches, total_windows_usize)?;

        let mut positional_decoded: Vec<FxHashMap<PositionDescriptor, DecodedCounts>> =
            Vec::with_capacity(positional_bins.len());
        let mut flattened: Vec<DecodedCounts> = Vec::new();

        for window_counts in positional_bins {
            let mut decoded_map: FxHashMap<PositionDescriptor, DecodedCounts> =
                FxHashMap::default();
            decoded_map.reserve(window_counts.len());
            for (descriptor, counts) in window_counts {
                let decoded = split_and_decode_counts(&counts, &kmer_specs);
                flattened.push(decoded.clone());
                decoded_map.insert(descriptor, decoded);
            }
            positional_decoded.push(decoded_map);
        }

        if flattened.is_empty() {
            flattened.push(DecodedCounts {
                counts: FxHashMap::default(),
            });
        }

        let (_, motifs_by_k) = prepare_decoded_counts(&flattened, opt.canonical, &kmer_specs);

        if options.log_statuses {
            info!(target: COMMAND_TARGET, "Writing positional counts to disk");
        }
        let temp_output_paths = write_positional_output(
            &positional_decoded,
            &motifs_by_k,
            &kmer_specs,
            final_outputs.temp_dir(),
            &opt.shared_args.output_prefix,
            opt.save_sparse,
        )?;
        output_files.extend(final_paths_for_same_named_temp_files(
            &temp_output_paths,
            &opt.shared_args.ioc.output_dir,
        )?);
        // These files were written to final_outputs.temp_dir() with their final filenames
        // Record each one as output_dir/<file name>, then move all outputs at the end
        final_outputs.record_temp_files_with_same_names_in(
            temp_output_paths,
            &opt.shared_args.ioc.output_dir,
        )?;
    } else {
        let all_bins = merge_tile_counts(tile_count_batches, total_windows_usize, &kmer_specs)?;

        // Prepare counts to get correct motifs (collapsed, N-filtered, etc.)
        let (prepared_counts, motifs_by_k) =
            prepare_decoded_counts(&all_bins, opt.canonical, &kmer_specs);

        // Write counts to the temp folder first
        // They move into output_dir after all requested output files have been written
        if options.log_statuses {
            info!(target: COMMAND_TARGET, "Writing counts to disk");
        }
        let temp_output_paths = write_decoded_counts_matrix(
            &prepared_counts,
            &kmer_specs,
            &motifs_by_k,
            final_outputs.temp_dir(),
            &opt.shared_args.output_prefix,
            opt.save_sparse,
        )?;
        output_files.extend(final_paths_for_same_named_temp_files(
            &temp_output_paths,
            &opt.shared_args.ioc.output_dir,
        )?);
        // These files were written to final_outputs.temp_dir() with their final filenames
        // Record each one as output_dir/<file name>, then move all outputs at the end
        final_outputs.record_temp_files_with_same_names_in(
            temp_output_paths,
            &opt.shared_args.ioc.output_dir,
        )?;
    }

    // Build bin metadata when windowed
    let bin_info = build_bin_info(
        &window_opt,
        &chromosomes,
        &contigs,
        windows_map.as_ref(),
        &blacklist_map,
        chr_offsets.as_ref(),
    )?;

    // Write window coordinates plus overlap metadata to the same temp folder as the count outputs
    if !matches!(window_opt, WindowSpec::Global) {
        if options.log_statuses {
            info!(target: COMMAND_TARGET, "Writing window coordinates to disk");
        }
        let bins_path = opt
            .shared_args
            .ioc
            .output_dir
            .join(dot_join(&[prefix, "bins.tsv"]));
        let temp_bins_path = final_outputs.temp_path_for(&bins_path)?;
        write_bin_info_tsv(&temp_bins_path, &bin_info)?;
        final_outputs.record(temp_bins_path, bins_path.clone())?;
        output_files.push(bins_path);
    }

    final_outputs.move_into_place()?;

    Ok(FragmentKmersRunResult {
        counters: global_counter,
        output_files,
    })
}

fn final_paths_for_same_named_temp_files(
    temp_paths: &[std::path::PathBuf],
    output_dir: &Path,
) -> Result<Vec<std::path::PathBuf>> {
    temp_paths
        .iter()
        .map(|temp_path| {
            let file_name = temp_path.file_name().with_context(|| {
                format!(
                    "temporary output path has no filename: {}",
                    temp_path.display()
                )
            })?;
            Ok(output_dir.join(file_name))
        })
        .collect()
}

/// Process a single tile: stream fragments, accumulate per-window counts, and persist results.
fn process_tile(
    opt: &FragmentKmersConfig,
    tile: &Tile,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    position_cache: Arc<PositionSelectionCache>,
    window_ctx: &WindowContext,
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_intervals: &[Interval<u64>],
    scaling_chr: &[ScalingBin],
    gc_corrector_opt: Option<GCCorrector>,
    counts_path: &Path,
) -> anyhow::Result<TileResult> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) =
        create_chromosome_reader(&opt.shared_args.ioc.bam, &tile.chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = FragmentKmersCounters::default();

    let gc_prefixes_opt = if gc_corrector_opt.is_some() {
        // NOTE: This sequence uses fetch coordinates for enabling GC correction
        // for all fetched fragments. It is only used to build the GC and ACGT prefixes.
        // The later sequence used for kmer extraction uses tile core coordinates.
        let seq_bytes = read_seq_in_range(
            &opt.shared_args.ref_genome.ref_2bit,
            &tile.chr,
            // NOTE: Need for full fetch span to get GC of overlapping fragments!
            (tile.fetch_start() as usize)..(tile.fetch_end() as usize),
        )?;
        Some(build_gc_prefixes(&seq_bytes))
    } else {
        None
    };

    let Some(fetch_span) = fetch_span_for_tile(
        tile,
        tile_window_span,
        window_ctx.windows_slice(),
        window_ctx.spec,
        chrom_len,
        opt.shared_args.fragment_lengths.max_fragment_length as u64,
        BedFetchPolicy::CoreOverlap,
    )?
    else {
        return Ok(TileResult {
            chr: tile.chr.clone(),
            counts_path: None,
            counter: FragmentKmersCounters::default(),
        });
    };
    let (fetch_from, fetch_to) = fetch_span.try_to_i64()?.as_tuple();

    // Extend the reference slice to include k-mers at the right tile edge
    let max_k: u32 = kmer_specs.keys().copied().max().unwrap_or(1) as u32;
    let seq_end_abs = (tile.core_end() as u64)
        .saturating_add((max_k as u64).saturating_sub(1))
        .min(chrom_len) as usize;

    let mut seq_bytes = read_seq_in_range(
        &opt.shared_args.ref_genome.ref_2bit,
        &tile.chr,
        (tile.core_start() as usize)..(seq_end_abs),
    )?;

    apply_blacklist_mask_to_seq(
        &mut seq_bytes,
        blacklist_intervals,
        tile.core_start() as u64,
    );

    // Scaled weights to count up
    let positional_scaling_weights = if !scaling_chr.is_empty() {
        let mut scaling_weights = vec![1.0; seq_bytes.len()];
        apply_scaling_to_coverage_in_place(&mut scaling_weights, tile.core_start(), scaling_chr);
        // "Blacklist" positions with scaling factors of 0, so they don't get counted
        for (base, weight) in seq_bytes.iter_mut().zip(&scaling_weights) {
            if *weight == 0.0 {
                *base = BLACKLIST_BYTE;
            }
        }
        Some(scaling_weights)
    } else {
        None
    };

    // Prepare left-aligned kmer-codes for each kmer-size
    let positional_codes_by_k: FxHashMap<u8, KmerCodes> =
        build_left_aligned_codes_per_k(&seq_bytes, kmer_specs);

    // Sparse map keyed by original window index -> kmer counts
    let mut counts_by_window: FxHashMap<u64, FxHashMap<CountKey, f64>> = FxHashMap::default();

    // Streaming pointers
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = tile_window_span
        .and_then(|span| (!span.is_empty()).then_some(span.first_idx))
        .unwrap_or(0);

    reader
        .fetch((tid, fetch_from, fetch_to))
        .context(format!("fetch {}", &tile.chr))?;

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let lengths = opt.shared_args.fragment_lengths.clone();
        move |f: &FragmentWithKmerSegments| lengths.contains(f.len())
    };

    // Create fragment iterator
    let gc_tag_bytes = opt
        .shared_args
        .gc
        .gc_tag
        .as_deref()
        .map(|t| t.as_bytes().to_vec());
    let unpaired = opt.shared_args.unpaired.reads_are_fragments;
    let include_read_fn: Box<dyn Fn(&Record) -> bool + Send + Sync> = if unpaired {
        let min_mapq = opt.shared_args.min_mapq;
        Box::new(move |r: &Record| default_include_read_unpaired(r, min_mapq))
    } else {
        let min_mapq = opt.shared_args.min_mapq;
        let require_proper_pair = opt.shared_args.require_proper_pair;
        Box::new(move |r: &Record| {
            default_include_read_paired_end(r, require_proper_pair, min_mapq)
        })
    };
    let mut iter = fragments_with_kmer_segments_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        move |rec| include_read_fn(rec),
        opt.shared_args.indel_mode,
        !opt.shared_args.ignore_gap,
        0,
        gc_tag_bytes.as_deref(),
        fragment_filter,
        unpaired,
    )
    .with_local_counters();

    let get_gc_weight = {
        let gc_corrector = gc_corrector_opt.as_ref();
        let gc_prefixes = gc_prefixes_opt.as_ref();
        move |fragment: &FragmentWithKmerSegments, fetch_start: u32| -> Result<Option<f64>> {
            match (gc_corrector, gc_prefixes) {
                (Some(corrector), Some(prefixes)) => {
                    let fetch_relative_fragment = fragment
                        .interval
                        .try_to_u64()?
                        .shift_left(fetch_start as u64)?;
                    corrector.correct_fragment(fetch_relative_fragment, prefixes)
                }
                _ => Ok(None),
            }
        }
    };

    let correct_gc = opt.shared_args.gc.gc_file.is_some();
    let fetch_start = tile.fetch_start();

    let store_positions = opt.positional_counts;
    let has_nearest_frame = position_cache
        .present_frames
        .contains(&ReferenceFrame::Nearest);

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;

        // Determine blacklist status
        let in_blacklist = is_blacklisted(
            blacklist_intervals,
            opt.shared_args.blacklist_strategy,
            fragment.interval.try_to_u64()?,
            opt.shared_args.fragment_lengths.max_fragment_length as u64,
            &mut bl_ptr,
        );
        if in_blacklist {
            counter.blacklisted_fragments += 1;
            continue;
        }

        // Get GC correction weight
        let gc_weight = if opt.shared_args.gc.gc_tag.is_some() {
            match fragment.gc_tag.classify()? {
                ClassifiedGCTagWeight::Usable(weight) => weight as f64,
                ClassifiedGCTagWeight::Missing => {
                    counter.gc_failed_fragments += 1;
                    counter.gc_missing_tags += 1;
                    if opt.shared_args.gc.neutralize_invalid_gc {
                        1.0
                    } else {
                        continue;
                    }
                }
                ClassifiedGCTagWeight::Invalid { out_of_range } => {
                    counter.gc_failed_fragments += 1;
                    if out_of_range {
                        counter.gc_out_of_range_tags += 1;
                    }
                    if opt.shared_args.gc.neutralize_invalid_gc {
                        1.0
                    } else {
                        continue;
                    }
                }
            }
        } else {
            let gc_weight_opt = get_gc_weight(&fragment, fetch_start)?;
            match (gc_weight_opt, correct_gc) {
                (Some(w), true) => w,
                (None, true) => {
                    // Tried but failed to make a GC correction weight
                    counter.gc_failed_fragments += 1;
                    if opt.shared_args.gc.neutralize_invalid_gc {
                        1.0
                    } else {
                        continue;
                    }
                }
                (None, false) => 1.0, // No correction
                (Some(_), false) => unreachable!(),
            }
        };

        // TODO: Does first, last need to be recalculated every iteration?
        let cache = position_cache.as_ref();
        let (first, last) = cache
            // Use smallest possible k to include all positions in interval for overlap
            .bounds(fragment.len(), cache.offsets.keys().copied().min().unwrap())
            .expect("non-empty offsets must have bounds");
        let interval_start = fragment.start() as u64 + first as u64;
        let interval_end = fragment.start() as u64 + last as u64 + 1;
        if interval_start >= interval_end {
            continue;
        }

        // Find all overlapping count-windows
        debug_assert!(interval_start >= fragment.start() as u64);
        let lookback_distance = opt.shared_args.fragment_lengths.max_fragment_length as u64
            + (interval_start - fragment.start() as u64);
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            window_ctx.windows_slice(),
            opt.shared_args.windows.by_size,
            Interval::new(interval_start, interval_end)?,
            1. / (opt.shared_args.fragment_lengths.max_fragment_length as f64 + 1.0), // Any overlap
            lookback_distance,
        )?;
        let overlapping_windows = if let Some(overlaps) = overlapping_windows {
            overlaps
        } else {
            continue;
        };

        counter.base.counted_fragments += 1;

        for overlapped_window in overlapping_windows.windows {
            let original_idx = window_ctx.original_idx(overlapped_window.idx);
            let counts = counts_by_window.entry(original_idx).or_default();
            count_kmers_at_positions(
                &fragment,
                cache,
                store_positions,
                &positional_codes_by_k,
                kmer_specs,
                counts,
                positional_scaling_weights.as_deref(),
                gc_weight,
                tile.core_start(),
                tile.core_end(),
                has_nearest_frame,
            );
        }
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    let mut count_records: Vec<TileWindowCounts> = counts_by_window
        .into_iter()
        .filter_map(|(original_idx, hm)| {
            if hm.is_empty() {
                return None;
            }
            let mut entries: Vec<TileKmerCountEntry> = Vec::with_capacity(hm.len());
            for (key, value) in hm {
                entries.push(TileKmerCountEntry::from((key, value)));
            }
            Some(TileWindowCounts {
                original_idx,
                entries,
            })
        })
        .collect();
    count_records.sort_unstable_by_key(|w| w.original_idx);

    serialize_tile_counts(counts_path, &count_records)?;

    Ok(TileResult {
        chr: tile.chr.clone(),
        counts_path: Some(counts_path.to_path_buf()),
        counter,
    })
}

pub(crate) fn count_kmers_at_positions(
    fragment: &FragmentWithKmerSegments,
    cache: &PositionSelectionCache,
    store_positions: bool,
    positional_codes_by_k: &FxHashMap<u8, KmerCodes>,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    counts: &mut FxHashMap<CountKey, f64>,
    scaling_weights: Option<&[f32]>,
    gc_weight: f64,
    tile_core_start: u32,
    tile_core_end: u32,
    apply_nearest_guard: bool,
) {
    // We perform comparisons in absolute genome coordinates first, then translate
    // back to fragment-relative offsets only after clipping
    let fragment_start = fragment.start() as u64;
    let tile_start = tile_core_start as u64;
    let tile_end = tile_core_end as u64;

    // We walk the requested k values independently. Each k has its own positional
    // encoding table, so processing them in isolation keeps the hot loop simple
    for &k in kmer_specs.keys() {
        let codes = positional_codes_by_k
            .get(&k)
            .expect("missing positional codes for requested k");
        let k_span = k as u64;
        let selections = match cache.offsets(fragment.len(), k) {
            Some(slice) if !slice.is_empty() => slice,
            _ => {
                continue;
            }
        };
        if selections.is_empty() {
            // Some frames filter out every position for a fragment of a given length
            continue;
        }

        // In count_kmers_at_positions, fetch windows once per k/fragment length
        let windows = match cache.windows(fragment.len(), k) {
            Some(w) => w,
            None => continue,
        };

        // Selections are sorted by offset. We stream through them once per k
        // using a single cursor so the overall complexity stays linear in the number
        // of usable offsets
        let mut offset_cursor = 0usize;

        // Precompute midpoint guards for Nearest: forbid crossing the true midpoint.
        // Rule:
        // - If there IS a physical midpoint (odd length), exclude that base entirely:
        //     forward:  start + (k-1) <= mid-1  -> start <= mid - k
        //     reverse:  start >= mid+1          -> anchor(offset) >= mid + (k-1) + 1 = mid + k
        // - If there is NO physical midpoint (even length), pick the base nearest each side's start:
        //     left  boundary = L/2 - 1 (0-based), right boundary = L/2
        //     forward:  start + (k-1) <= left_boundary   -> start <= (L/2) - k
        //     reverse:  start >= right_boundary          -> anchor(offset) >= (L/2) + (k-1)
        let nearest_guard =
            NearestFrameGuard::by_flag(apply_nearest_guard, fragment.len(), k as u32);

        // Fragments may be gapped by indels, so we examine each contiguous segment
        // and clip it to the tile coordinates before accepting offsets
        'segments: for segment in &fragment.segments {
            let seg_start = segment.start() as u64;
            let seg_end = segment.end() as u64;

            if seg_end <= tile_start {
                // Segment lies completely before the tile
                continue;
            }

            if offset_cursor >= selections.len() {
                // Offsets arrive sorted and we only advance the cursor, therefore exhausting the
                // selections list guarantees no later segment has a usable position
                break;
            }

            if seg_start >= tile_end {
                // Segments are emitted in genomic order, so hitting the tile boundary means the rest
                // start beyond the core and cannot contribute any k-mers
                break;
            }

            let forward_min_abs = seg_start.max(tile_start);
            let effective_seg_end = seg_end.min(tile_end.saturating_add(k_span.saturating_sub(1)));
            let Some(last_start_abs) = effective_seg_end.checked_sub(k_span) else {
                continue;
            };

            // Forward oriented kmers start at the position we count
            // Only offsets whose entire span stays inside both the segment and tile are valid
            let forward_range = if last_start_abs >= forward_min_abs {
                Some((
                    forward_min_abs.saturating_sub(fragment_start),
                    last_start_abs.saturating_sub(fragment_start),
                ))
            } else {
                None
            };

            let reverse_anchor_min_abs = seg_start
                .saturating_add(k_span.saturating_sub(1))
                .max(tile_start.saturating_add(k_span.saturating_sub(1)));
            let reverse_anchor_max_abs = seg_end.min(tile_end);
            let reverse_range = if reverse_anchor_max_abs == 0 {
                None
            } else {
                let max_inclusive = reverse_anchor_max_abs.saturating_sub(1);
                if max_inclusive >= reverse_anchor_min_abs {
                    // Reverse oriented kmers are indexed by their last base
                    // After clipping to the tile and segment we backtrack by k-1
                    // bases to locate the start
                    Some((
                        reverse_anchor_min_abs.saturating_sub(fragment_start),
                        max_inclusive.saturating_sub(fragment_start),
                    ))
                } else {
                    None
                }
            };

            if forward_range.is_none() && reverse_range.is_none() {
                // Clipping removed every valid orientation for this segment
                continue;
            }

            // Build an inclusive span that covers whichever orientations survived
            let segment_range_start = match (forward_range, reverse_range) {
                (Some((fwd_min, _)), Some((rev_min, _))) => fwd_min.min(rev_min),
                (Some((fwd_min, _)), None) => fwd_min,
                (None, Some((rev_min, _))) => rev_min,
                (None, None) => unreachable!(),
            };

            let segment_range_end = match (forward_range, reverse_range) {
                (Some((_, fwd_max)), Some((_, rev_max))) => fwd_max.max(rev_max),
                (Some((_, fwd_max)), None) => fwd_max,
                (None, Some((_, rev_max))) => rev_max,
                (None, None) => unreachable!(),
            };

            while offset_cursor < selections.len()
                && (selections[offset_cursor].offset() as u64) < segment_range_start
            {
                // The cursor still points before the current segment window, so fast-forward it
                offset_cursor += 1;
            }

            let mut idx = offset_cursor;
            while idx < selections.len() {
                let selection = selections[idx];
                let offset = selection.offset() as u64;
                let offset_i32 = selection.offset() as i32;
                if offset > segment_range_end {
                    // Remaining selections start after this segment range, so move to next segment
                    break;
                }

                let idx_abs = fragment_start + offset;
                if idx_abs >= tile_end {
                    // Offsets are ordered, so reaching the right edge of the tile means we are done for this k
                    break 'segments;
                }
                if idx_abs < tile_start {
                    // Offset still lies before the tile core; advance to the next candidate
                    idx += 1;
                    continue;
                }

                let decision = evaluate_selection(
                    selection,
                    windows,
                    nearest_guard.as_ref(),
                    k_span,
                    offset,
                    forward_range,
                    reverse_range,
                );

                match decision {
                    SelectionDecision::SkipAdvance => {
                        idx += 1;
                        continue;
                    }
                    SelectionDecision::IncludeForward { .. } => {
                        let start_local = match idx_abs.checked_sub(tile_start) {
                            Some(val) => val as usize,
                            None => {
                                idx += 1;
                                continue;
                            }
                        };

                        // Ensure the forward k-mer stays within this contiguous segment
                        // idx_abs is the start. Require idx_abs + (k-1) < seg_end
                        if idx_abs.saturating_add(k_span.saturating_sub(1)) >= seg_end {
                            idx += 1;
                            continue;
                        }

                        // We look up weights using the same tile-relative index as the positional codes
                        let weight = match scaling_weights {
                            Some(w) => unsafe { *w.get_unchecked(start_local) as f64 },
                            None => 1.0,
                        } * gc_weight;

                        // Record the forward kmer code emitted at this start position
                        let key = CountKey {
                            k,
                            code: codes.get(start_local),
                            position: store_positions.then_some(offset_i32),
                            group: selection.group(),
                        };
                        *counts.entry(key).or_insert(0.0) += weight;
                    }
                    SelectionDecision::IncludeReverse {
                        start_offset_0,
                        anchor_offset_0,
                    } => {
                        let kmer_start_abs = fragment_start + start_offset_0;
                        if kmer_start_abs < tile_start || kmer_start_abs < seg_start {
                            idx += 1;
                            continue;
                        }

                        let start_local = match kmer_start_abs.checked_sub(tile_start) {
                            Some(val) => val as usize,
                            None => {
                                idx += 1;
                                continue;
                            }
                        };

                        let end_local =
                            match (fragment_start + anchor_offset_0).checked_sub(tile_start) {
                                Some(val) => val as usize,
                                None => {
                                    idx += 1;
                                    continue;
                                }
                            };

                        let weight = match scaling_weights {
                            Some(w) => unsafe { *w.get_unchecked(end_local) as f64 },
                            None => 1.0,
                        } * gc_weight;

                        let key = CountKey {
                            k,
                            code: codes.get(start_local),
                            position: store_positions.then_some(offset_i32),
                            group: selection.group(),
                        };
                        *counts.entry(key).or_insert(0.0) += weight;
                    }
                }

                idx += 1;
            }

            // Carry the cursor forward so the next segment starts scanning from the last visited offset
            offset_cursor = idx;
        }
    }
}
