use crate::shared::gc_tag::ClassifiedGCTagWeight;
use crate::{
    commands::{
        cli_common::{
            DistributionWindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        counters::EndsCounters,
        ends::{
            config::EndsConfig,
            config_structs::{ClipStrategy, KmerSource, WindowMotifAssigner},
            counting::{EndCountsByWindow, decode_end_motif_counts},
            motifs::{
                CountedEndFlags, build_optional_kmer_spec, build_tile_motif_context,
                count_fragment_in_window, motif_extraction_ref_2bit_requirement_message,
                motif_extraction_requires_reference, motif_reference_span_for_tile,
            },
            output::{
                build_all_end_motif_order, collect_end_motif_order,
                ensure_all_motifs_enumeration_size, ensure_dense_end_motif_output_size,
            },
            tiling::{
                TileResult, build_tile_payload, deserialize_tile_counts, merge_tile_payload,
                serialize_tile_counts,
            },
            write::{write_end_motif_outputs, write_end_settings_json},
        },
        gc_bias::{
            correct::{GCCorrector, load_gc_corrector},
            counting::build_gc_prefixes,
        },
        run_statistics::{
            FragmentRunStatisticsOptions, FragmentStatisticsLabels, GCStatisticsSummary,
            print_fragment_run_statistics,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::{load_grouped_windows_from_bed, load_windows_from_bed},
        blacklist::is_blacklisted,
        fragment::ends_fragment::FragmentWithEnds,
        fragment_iterators::fragments_with_ends_from_bam,
        interval::{IndexedInterval, Interval},
        io::dot_join,
        midpoint::midpoint_random_even_with_thread_rng,
        overlaps::find_overlapping_windows,
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq_in_range,
        scale_genome::{compute_window_scaling_over_fragment, compute_window_scaling_over_overlap},
        thread_pool::init_global_pool,
        tiled_run::{
            TempDirGuard, Tile, TileWindowSpan, build_tiles, precompute_tile_window_spans,
        },
        window_fetch::{BedFetchPolicy, fetch_span_for_tile},
        windowing::{
            WindowContext, build_bin_info, compute_window_offsets,
            ensure_plain_bed_windows_not_empty, write_bin_info_tsv,
            write_group_index_with_blacklist_tsv,
        },
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{convert::TryInto, path::Path, sync::Arc, time::Instant};
use tracing::{info, warn};

const COMMAND_TARGET: &str = "ends";

fn clip_strategy_name(strategy: ClipStrategy) -> &'static str {
    match strategy {
        ClipStrategy::Aligned => "aligned",
        ClipStrategy::RawAlignedBoundary => "raw-aligned-boundary",
        ClipStrategy::RawShiftedBoundary => "raw-shifted-boundary",
        ClipStrategy::Skip => "skip",
    }
}

fn outside_kmer_clip_strategy_warning(
    k_outside: usize,
    clip_strategy: ClipStrategy,
) -> Option<String> {
    if k_outside == 0 || matches!(clip_strategy, ClipStrategy::Skip) {
        return None;
    }

    Some(format!(
        "`--k-outside > 0` with `--clip-strategy {}` will likely add more noise than signal when soft clipping is present, as it is hard to determine where the outside motif actually lies on the reference. Prefer `--clip-strategy skip` for outside-base counting",
        clip_strategy_name(clip_strategy)
    ))
}

/// Execute the end-motif counting pipeline end-to-end.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `ends` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &EndsConfig) -> Result<()> {
    let start_time = Instant::now();
    opt.fragment_lengths.validate()?;
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    if opt.k_inside == 0 && opt.k_outside == 0 {
        bail!("At least one of --k-inside or --k-outside must be > 0");
    }
    if !opt.bq_filter.is_empty() {
        if opt.k_inside == 0 {
            bail!(
                "`--bq-filter` requires `--k-inside > 0` because it scores the inside read bases"
            );
        }
        if matches!(opt.source_inside, KmerSource::Reference) {
            bail!(
                "`--bq-filter` cannot be combined with `--source-inside reference` because reference-backed inside bases do not have read base qualities"
            );
        }
    }
    if matches!(opt.clip.clip_strategy, ClipStrategy::RawAlignedBoundary)
        && matches!(opt.source_inside, KmerSource::Reference)
    {
        bail!(
            "`--clip-strategy raw-aligned-boundary` cannot be combined with `--source-inside reference`"
        );
    }
    if opt.ref_2bit.is_none() && motif_extraction_requires_reference(opt, opt.blacklist.is_some()) {
        bail!(motif_extraction_ref_2bit_requirement_message());
    }
    if let Some(warning_message) =
        outside_kmer_clip_strategy_warning(opt.k_outside, opt.clip.clip_strategy)
    {
        warn!(target: COMMAND_TARGET, "{warning_message}");
    }
    opt.gc.validate(opt.ref_2bit.as_deref())?;
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    let window_opt = opt.windows.resolve_windows();
    let fetch_window_opt = window_opt.as_fetch_window_spec();
    let prefix = opt.output_prefix.trim();

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
        DistributionWindowSpec::Bed(bed) => {
            info!(target: COMMAND_TARGET, "Loading window coordinates");
            let windows = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            ensure_plain_bed_windows_not_empty(&windows)?;
            Some(windows)
        }
        _ => None,
    };
    let (grouped_windows_map, group_idx_to_name) = match &window_opt {
        DistributionWindowSpec::GroupedBed(bed) => {
            info!(target: COMMAND_TARGET, "Loading grouped window coordinates");
            let (windows_map, group_idx_to_name) =
                load_grouped_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            ensure!(
                !group_idx_to_name.is_empty(),
                "grouped BED file did not contain any valid windows on the selected chromosomes"
            );
            (Some(windows_map), Some(group_idx_to_name))
        }
        _ => (None, None),
    };
    let indexed_windows_map: Option<FxHashMap<String, Vec<IndexedInterval<u64>>>> =
        if let Some(windows_map) = windows_map.as_ref() {
            Some(
                windows_map
                    .iter()
                    .map(|(chr, windows)| (chr.clone(), windows.as_slice().to_vec()))
                    .collect(),
            )
        } else {
            grouped_windows_map.as_ref().map(|windows_map| {
                windows_map
                    .iter()
                    .map(|(chr, windows)| (chr.clone(), windows.as_slice().to_vec()))
                    .collect()
            })
        };

    // Load genomic scaling factors
    if opt.scale_genome.scaling_factors.is_some() {
        info!(target: COMMAND_TARGET, "Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> = load_scaling_map(
        &opt.scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(
            opt.gc.gc_file.is_some(),
            opt.gc.gc_tag.is_some(),
        ),
        None,
    )?;

    // Load GC correction package if specified
    if opt.gc.gc_file.is_some() {
        info!(target: COMMAND_TARGET, "Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.gc.gc_file.as_ref(),
        opt.ref_2bit.as_ref(),
        opt.fragment_lengths.min_fragment_length,
        opt.fragment_lengths.max_fragment_length,
    )?;

    let halo_bp = opt.fragment_lengths.max_fragment_length;
    let align_bp = match &window_opt {
        DistributionWindowSpec::Size(bp) => Some(*bp),
        _ => None,
    };

    // Build tiles (core plus halo)
    let (tiles, _) = build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, align_bp)?;

    let progress = ProgressFactory::new();
    let pb = Arc::new(progress.default_bar(tiles.len() as u64));

    let tile_span_left_halo = if opt.clip.clip_strategy.uses_shifted_boundary() {
        opt.clip.max_soft_clips as u64
    } else {
        0
    };
    let tile_span_right_halo = opt.fragment_lengths.max_fragment_length as u64;

    let windows_lookup = indexed_windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(
        &tiles,
        |chr| {
            windows_lookup
                .and_then(|m| m.get(chr).map(|w| w.as_slice()))
                .unwrap_or(&[])
        },
        tile_span_left_halo,
        // We use fragments starting in a tile, so we need fragment-overlapping windows starting after the tile
        tile_span_right_halo,
    ));
    let tile_window_spans_for_threads = tile_window_spans.clone();
    // Window rows are global across chromosomes. For fixed-size windows we therefore need a
    // per-chromosome row offset to turn chromosome-local overlap indices into global output rows.
    // BED windows already carry their own original indices, so their offsets stay at zero.
    let (total_windows, chr_offsets_map): (u64, FxHashMap<String, u64>) = match &window_opt {
        DistributionWindowSpec::GroupedBed(_) => (
            group_idx_to_name
                .as_ref()
                .context("group_idx_to_name missing for grouped BED mode")?
                .len() as u64,
            chromosomes.iter().map(|chr| (chr.clone(), 0_u64)).collect(),
        ),
        _ => compute_window_offsets(
            &fetch_window_opt,
            &chromosomes,
            &contigs,
            windows_map.as_ref(),
        )?,
    };
    let chr_offsets = Arc::new(chr_offsets_map);
    let chr_offsets_for_threads = chr_offsets.clone();

    let temp_dir_guard =
        TempDirGuard::new(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;
    let temp_dir = temp_dir_guard.path();
    let counts_prefix = &dot_join(&[prefix, "counts"]);
    let inside_spec = build_optional_kmer_spec(opt.k_inside, "inside")?;
    let outside_spec = build_optional_kmer_spec(opt.k_outside, "outside")?;

    info!(target: COMMAND_TARGET, "Counting per tile");

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    pb.set_position(0);

    let tile_results: Vec<TileResult> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<Option<TileResult>> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let windows_chr: Option<&[IndexedInterval<u64>]> = indexed_windows_map
                .as_ref()
                .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice()));
            let blacklist_chr: &[Interval<u64>] = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr: &[(u64, u64, f32)] = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            // `find_overlapping_windows` reports chromosome-local indices. `process_tile` needs the
            // chromosome-specific offset so it can emit globally stable window ids into tile payloads.
            let chr_window_idx_offset = *chr_offsets_for_threads.get(&tile.chr).unwrap_or(&0);

            let tile_result = process_tile(
                opt,
                tile,
                tile_span.as_ref(),
                windows_chr,
                chr_window_idx_offset,
                &window_opt,
                blacklist_chr,
                scaling_chr,
                gc_corrector.clone(),
                temp_dir,
                counts_prefix,
                inside_spec.as_ref(),
                outside_spec.as_ref(),
            )?;
            pb.inc(1);
            Ok(tile_result)
        })
        .collect::<Result<Vec<_>>>()? // Short-circuits on the first Err
        .into_iter()
        .flatten()
        .collect();

    pb.finish_with_message("| Finished counting");

    // Release per-tile inputs before merging outputs
    drop(chr_offsets_for_threads);
    drop(tile_window_spans_for_threads);
    drop(tile_window_spans);
    drop(tiles);
    drop(scaling_map);
    drop(gc_corrector);

    // Collect counters
    let mut global_counter = EndsCounters::default();
    for tile_out in &tile_results {
        global_counter += tile_out.counter;
    }

    info!(target: COMMAND_TARGET, "Reducing temporary tile files");

    // Start from an all-empty output matrix shape, then fill only the windows that were actually
    // observed in the reduced sparse payloads.
    let mut all_bins = vec![FxHashMap::default(); total_windows as usize];
    let mut reduced_counts: EndCountsByWindow = FxHashMap::default();
    for tile_result in &tile_results {
        merge_tile_payload(
            &mut reduced_counts,
            deserialize_tile_counts(&tile_result.counts_path)?,
        )?;
    }

    // Decode each populated window independently. This is an easy parallel boundary because the
    // windows no longer interact after tile reduction; we only need a final serial pass to place
    // each decoded map into its global output row.
    let decoded_bins: Vec<(usize, FxHashMap<String, f64>)> = reduced_counts
        .into_par_iter()
        .map(
            |(original_idx, counts)| -> Result<(usize, FxHashMap<String, f64>)> {
                let decoded = decode_end_motif_counts(
                    &counts,
                    inside_spec.as_ref(),
                    outside_spec.as_ref(),
                    opt.collapse_complement,
                );
                let idx: usize = original_idx
                    .try_into()
                    .context("window index does not fit in usize")?;
                Ok((idx, decoded))
            },
        )
        .collect::<Result<_>>()?;

    for (idx, decoded) in decoded_bins {
        if idx >= all_bins.len() {
            bail!(
                "reduced window index {} is out of bounds for {} output windows",
                idx,
                all_bins.len()
            );
        }
        all_bins[idx] = decoded;
    }

    let bin_info = if matches!(&window_opt, DistributionWindowSpec::GroupedBed(_)) {
        Vec::new()
    } else {
        build_bin_info(
            &fetch_window_opt,
            &chromosomes,
            &contigs,
            windows_map.as_ref(),
            &blacklist_map,
            chr_offsets.as_ref(),
        )?
    };
    // `all_motifs` switches the final output from "observed motifs only" to a dense fixed motif
    // universe. The dense size checks happen before we allocate or enumerate that full universe.
    if opt.all_motifs {
        ensure_all_motifs_enumeration_size(opt.k_inside, opt.k_outside, all_bins.len())?;
    }
    let motif_order = if opt.all_motifs {
        build_all_end_motif_order(
            inside_spec.as_ref(),
            outside_spec.as_ref(),
            opt.collapse_complement,
        )?
    } else {
        collect_end_motif_order(&all_bins)
    };
    let write_dense_output = opt.all_motifs;
    if write_dense_output {
        ensure_dense_end_motif_output_size(all_bins.len(), motif_order.len())?;
    }
    write_end_motif_outputs(
        &opt.ioc.output_dir,
        prefix,
        &all_bins,
        &motif_order,
        write_dense_output,
    )?;

    write_end_settings_json(&opt.ioc.output_dir, prefix, opt)?;

    // Write window coordinates plus overlap metadata as TSV to output_dir
    match &window_opt {
        DistributionWindowSpec::GroupedBed(_) => {
            info!(target: COMMAND_TARGET, "Writing group metadata to disk");
            let group_idx_to_name = group_idx_to_name
                .as_ref()
                .context("group_idx_to_name missing when writing grouped outputs")?;
            let grouped_windows_map = grouped_windows_map
                .as_ref()
                .context("grouped windows missing when writing grouped outputs")?;
            let group_index_path = opt
                .ioc
                .output_dir
                .join(dot_join(&[prefix, "group_index.tsv"]));
            write_group_index_with_blacklist_tsv(
                group_index_path,
                group_idx_to_name,
                &chromosomes,
                grouped_windows_map,
                &blacklist_map,
                opt.blacklist.is_some(),
            )?;
        }
        DistributionWindowSpec::Global => {}
        _ => {
            info!(target: COMMAND_TARGET, "Writing window coordinates to disk");
            let bins_path = opt.ioc.output_dir.join(dot_join(&[prefix, "bins.tsv"]));
            write_bin_info_tsv(bins_path, &bin_info)?;
        }
    }

    drop(blacklist_map);

    let elapsed = start_time.elapsed();
    print_fragment_run_statistics(
        &global_counter.base,
        elapsed,
        FragmentRunStatisticsOptions {
            include_section_header: true,
            notes: &["Note: counts below cover only tiles with relevant output windows"],
            labels: FragmentStatisticsLabels {
                total_reads: "Observed reads in processed tiles",
                accepted_reads: "Initially accepted reads",
                counted_fragments: "Fragments with one or more counted motifs",
            },
            blacklist_excluded_fragments: Some(global_counter.blacklisted_fragments),
            gc: (opt.gc.gc_file.is_some() || opt.gc.gc_tag.is_some()).then_some(
                GCStatisticsSummary {
                    neutralize_invalid_gc: opt.gc.neutralize_invalid_gc,
                    failed_fragments: global_counter.gc_failed_fragments,
                    missing_tags: opt
                        .gc
                        .gc_tag
                        .is_some()
                        .then_some(global_counter.gc_missing_tags),
                    out_of_range_tags: opt
                        .gc
                        .gc_tag
                        .is_some()
                        .then_some(global_counter.gc_out_of_range_tags),
                },
            ),
        },
        [format!(
            "Distinct counted end motifs across those fragments: {}",
            global_counter.counted_motifs
        )],
    );
    Ok(())
}

/// Update the `ends` statistics for one fully processed fragment.
///
/// Fragment-level statistics should reflect emitted motif counts, not just
/// fragments that survived earlier filters. This helper applies the final
/// per-fragment flags after all candidate windows have been processed.
///
/// Parameters
/// ----------
/// - `counter`:
///   Tile-local counters updated in place
/// - `counted_end_flags`:
///   Which distinct end motifs were actually counted for the fragment
fn record_counted_fragment_stats(counter: &mut EndsCounters, counted_end_flags: CountedEndFlags) {
    if counted_end_flags.any_counted() {
        counter.base.counted_fragments += 1;
        counter.counted_motifs += counted_end_flags.counted_motif_total();
    }
}

/// Count all end motifs owned by one tile and write its sparse payload to disk.
///
/// This function does the tile-local heavy lifting for `ends`: it streams
/// fragments from BAM, applies fragment-level filters and weights, assigns each
/// owned fragment to candidate windows, and accumulates sparse motif counts for
/// those windows before serializing the tile result.
///
/// Parameters
/// ----------
/// - `opt`:
///   Full command configuration
/// - `tile`:
///   Tile currently being processed
/// - `tile_window_span`:
///   Precomputed window span relevant to this tile, if any
/// - `windows_chr`:
///   Window intervals for this chromosome, when not in global mode
/// - `chr_window_idx_offset`:
///   Per-chromosome row offset for fixed-size windows
/// - `window_opt`:
///   Resolved windowing mode
/// - `blacklist_intervals`:
///   Merged blacklist intervals for this chromosome
/// - `scaling_chr`:
///   Genomic scaling bins for this chromosome
/// - `gc_corrector_opt`:
///   Optional GC corrector shared from the outer runner
/// - `temp_dir`:
///   Temporary directory for tile payloads
/// - `counts_prefix`:
///   Prefix used when naming serialized tile payloads
/// - `inside_spec`:
///   Shared codec spec for the inside half, or `None` when `k_inside = 0`
/// - `outside_spec`:
///   Shared codec spec for the outside half, or `None` when `k_outside = 0`
///
/// Returns
/// -------
/// - `Result<Option<TileResult>>`:
///   `None` when the tile has no relevant windows, otherwise the serialized sparse tile result
fn process_tile(
    opt: &EndsConfig,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    windows_chr: Option<&[IndexedInterval<u64>]>,
    chr_window_idx_offset: u64,
    window_opt: &DistributionWindowSpec,
    blacklist_intervals: &[Interval<u64>],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<GCCorrector>,
    temp_dir: &Path,
    counts_prefix: &str,
    inside_spec: Option<&crate::shared::kmers::kmer_codec::KmerSpec>,
    outside_spec: Option<&crate::shared::kmers::kmer_codec::KmerSpec>,
) -> Result<Option<TileResult>> {
    let fetch_window_opt = window_opt.as_fetch_window_spec();
    // One BAM reader per tile
    let (mut reader, _tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    debug_assert_eq!(_tid_check, tile.tid as u32);

    let max_fragment_length = opt.fragment_lengths.max_fragment_length;

    // Counters
    let mut counter = EndsCounters::default();
    let counts_path = temp_dir.join(format!(
        "{prefix}.{chr}.{idx}.counts.bin",
        prefix = counts_prefix,
        chr = tile.chr.as_str(),
        idx = tile.index
    ));

    // Build GC prefixes for the full tile fetch span
    let gc_prefixes_opt = if gc_corrector_opt.is_some() {
        let ref_2bit = match opt.ref_2bit.as_ref() {
            Some(r) => r,
            None => bail!("When GC correction is specified, --ref-2bit must also be specified"),
        };
        let seq_bytes = read_seq_in_range(
            ref_2bit,
            &tile.chr,
            // NOTE: Need for full fetch span to get GC of overlapping fragments!
            (tile.fetch_start() as usize)..(tile.fetch_end() as usize),
        )?;
        Some(build_gc_prefixes(&seq_bytes))
    } else {
        None
    };

    // Narrow the BAM fetch to the part of the tile that can still contribute to the current
    // windows. In global/by-size modes this usually stays close to the tile fetch span; in BED
    // mode it can shrink substantially.
    let bed_fetch_halo_bp = opt.fragment_lengths.max_fragment_length as u64;
    let Some(fetch_span) = fetch_span_for_tile(
        tile,
        tile_window_span,
        windows_chr,
        &fetch_window_opt,
        chrom_len,
        bed_fetch_halo_bp,
        BedFetchPolicy::CandidateWindowExtent,
    )?
    else {
        // Skip tiles with no relevant windows
        return Ok(None);
    };
    let (fetch_from, fetch_to) = fetch_span.try_to_i64()?.as_tuple();
    let reference_span = motif_reference_span_for_tile(
        tile,
        chrom_len,
        opt.clip.clip_strategy,
        opt.clip.max_soft_clips,
        opt.k_outside,
    )?;
    let motif_context = build_tile_motif_context(
        opt,
        tile,
        reference_span,
        chrom_len,
        blacklist_intervals,
        inside_spec,
        outside_spec,
    )?;
    let window_context = WindowContext {
        spec: &fetch_window_opt,
        windows: windows_chr,
        chr_idx_offset: chr_window_idx_offset,
    };

    reader
        .fetch((tile.tid, fetch_from, fetch_to))
        .context(format!("fetch {} {}-{}", &tile.chr, fetch_from, fetch_to))?;

    let mut counts_by_window: EndCountsByWindow = FxHashMap::default();

    // Fraction of a fragment that must overlap with a window to consider that window as a
    // candidate. Endpoint mode still uses the fragment assignment interval here; the actual
    // left/right terminal-base checks happen later in `count_fragment_in_window(...)`.
    let min_overlap_fraction: f64 = match opt.window_assignment.assign_by {
        WindowMotifAssigner::Any
        | WindowMotifAssigner::Endpoint
        | WindowMotifAssigner::CountOverlap => 1. / (max_fragment_length as f64 + 1.0), // +1 to avoid rounding error issues
        WindowMotifAssigner::All | WindowMotifAssigner::Midpoint => {
            1.0 - (1. / (max_fragment_length as f64 + 1.0))
        } // 1.0 but just below to avoid rounding errors
        WindowMotifAssigner::Proportion(p) => p,
    };

    // The overlap finder only needs checked BED-like intervals here.
    //
    // In BED mode, `find_overlapping_windows(...)` returns scan positions in the supplied slice as
    // `OverlappingWindow.idx`; it does not use `IndexedInterval.idx`. This temporary list is built
    // in the same order as `scaling_chr`, so those scan positions are already the correct
    // chromosome-local indices for indexing back into `scaling_chr` later.
    //
    // So the carried `IndexedInterval.idx` value is intentionally a placeholder.
    let scaling_with_bin_idx: Vec<IndexedInterval<u64>> = scaling_chr
        .iter()
        .map(|(start, end, _)| IndexedInterval::new(*start, *end, 0_u64))
        .collect::<crate::Result<_>>()?;

    // Function for filtering fragments after pairing
    let fragment_filter = {
        let lengths = opt.fragment_lengths.clone();
        move |fragment: &FragmentWithEnds| lengths.contains(fragment.assignment_len())
    };

    // Create fragment iterator with per-tile filtering and optional GC tag handling
    let unpaired = opt.unpaired.reads_are_fragments;
    let max_soft_clips: u32 = opt
        .clip
        .max_soft_clips
        .try_into()
        .context("max_soft_clips does not fit in u32 for fragment iteration")?;
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
    let mut iter = fragments_with_ends_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        move |rec| include_read_fn(rec),
        opt.clip.clip_strategy,
        opt.source_inside,
        opt.indel_filter,
        opt.k_inside,
        max_soft_clips,
        &opt.bq_filter,
        opt.gc.gc_tag.as_deref().map(str::as_bytes),
        fragment_filter,
        unpaired,
    )
    .with_local_counters();

    let get_gc_weight = {
        let gc_corrector = gc_corrector_opt.as_ref();
        let gc_prefixes = gc_prefixes_opt.as_ref();
        let fetch_start = tile.fetch_start();
        move |fragment: &FragmentWithEnds| -> Result<Option<f64>> {
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

    let correct_gc = opt.gc.gc_file.is_some();

    // Streaming pointers
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = tile_window_span
        .and_then(|span| (!span.is_empty()).then_some(span.first_idx))
        .unwrap_or(0);
    let mut sf_ptr = 0; // Scaling factor bin

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;

        // Only count fragments whose start is inside the core to prevent double counting across tiles
        if fragment.start() < tile.core_start() || fragment.start() >= tile.core_end() {
            continue;
        }

        // Determine blacklist status (based on aligned fragment coordinates)
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

        // First find candidate windows from fragment geometry alone. We intentionally do this
        // before motif extraction so all later work only happens for windows that can actually
        // receive counts.
        let query_interval = match opt.window_assignment.assign_by {
            WindowMotifAssigner::Midpoint => {
                let fragment_assignment_length = fragment.assignment_len();
                let midpoint = midpoint_random_even_with_thread_rng(
                    fragment.assignment_start(),
                    fragment_assignment_length,
                );
                Interval::new(midpoint.into(), (midpoint + 1).into())?
            }
            WindowMotifAssigner::Any
            | WindowMotifAssigner::All
            | WindowMotifAssigner::Proportion(_)
            | WindowMotifAssigner::CountOverlap
            | WindowMotifAssigner::Endpoint => fragment.assignment_interval.try_to_u64()?,
        };
        let by_size = match window_opt {
            DistributionWindowSpec::Size(bp) => Some(*bp),
            _ => None,
        };
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            windows_chr,
            by_size,
            query_interval,
            min_overlap_fraction,
            max_fragment_length.into(),
        )?;
        let overlapping_windows = if let Some(overlaps) = overlapping_windows {
            overlaps
        } else {
            continue;
        };

        // GC correction is fragment-level, so the same GC weight is reused for every window and
        // every end motif produced from this fragment.
        let gc_weight = if opt.gc.gc_tag.is_some() {
            match fragment.gc_tag.classify()? {
                ClassifiedGCTagWeight::Usable(weight) => weight as f64,
                ClassifiedGCTagWeight::Missing => {
                    counter.gc_failed_fragments += 1;
                    counter.gc_missing_tags += 1;
                    if opt.gc.neutralize_invalid_gc {
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
                    if opt.gc.neutralize_invalid_gc {
                        1.0
                    } else {
                        continue;
                    }
                }
            }
        } else {
            let gc_weight_opt = get_gc_weight(&fragment)?;
            match (gc_weight_opt, correct_gc) {
                (Some(w), true) => w,
                (None, true) => {
                    // Tried but failed to make a GC correction weight for the current fragment
                    if opt.gc.neutralize_invalid_gc {
                        counter.gc_failed_fragments += 1;
                        1.0
                    } else {
                        counter.gc_failed_fragments += 1;
                        continue;
                    }
                }
                (None, false) => 1.0, // No correction
                (Some(_), false) => bail!("unexpected GC weight when GC correction is disabled"),
            }
        };
        let mut counted_end_flags = CountedEndFlags::default();

        if !scaling_chr.is_empty() {
            // Find overlapping scaling-bins
            let overlapping_scaling_bins = find_overlapping_windows(
                chrom_len,
                &mut sf_ptr,
                Some(&scaling_with_bin_idx),
                None,
                fragment.interval.try_to_u64()?, // Full fragment
                1. / (max_fragment_length as f64 + 1.0), // Any overlap without rounding error issues
                max_fragment_length.into(),
            )
            .with_context(|| format!("finding overlapping scaling bins on chr {}", tile.chr))?
            .with_context(|| {
                format!(
                    "no overlapping scaling bins found for fragment {}:{}-{}. Scaling factors must cover every counted base on every counted chromosome",
                    tile.chr,
                    fragment.start(),
                    fragment.end()
                )
            })?;

            // Extract the indices of the overlapping bins
            let overlapping_scaling_bin_indices: Vec<usize> = overlapping_scaling_bins
                .windows
                .iter()
                .map(|w| w.idx)
                .collect();

            // Calculate the weight per overlapping count-window
            // NOTE: `compute_window_scaling_over_fragment` always returns
            // an overlap fraction of 1.0 (count full fragment)!
            let overlap_weights = match opt.window_assignment.assign_by {
                WindowMotifAssigner::CountOverlap => compute_window_scaling_over_overlap(
                    &overlapping_windows,
                    &overlapping_scaling_bin_indices,
                    scaling_chr,
                )?,
                _ => compute_window_scaling_over_fragment(
                    fragment.interval.try_to_u64()?,
                    &overlapping_windows,
                    &overlapping_scaling_bin_indices,
                    scaling_chr,
                )?,
            };

            // Count up the weight per overlapping count-window. `count_fragment_in_window(...)`
            // still decides whether each end is actually counted in the current window.
            let overlapping_window_intervals: FxHashMap<usize, Interval<u64>> = overlapping_windows
                .windows
                .iter()
                .map(|window| (window.idx, window.interval))
                .collect();
            for (overlapped_window_idx, scaling_weight, overlap_fraction_to_count) in
                overlap_weights
            {
                let original_idx = window_context.original_idx(overlapped_window_idx);
                let window_interval = *overlapping_window_intervals
                    .get(&overlapped_window_idx)
                    .expect("missing overlap interval for scaled count window");
                let count_weight = overlap_fraction_to_count * scaling_weight * gc_weight;
                counted_end_flags.merge(count_fragment_in_window(
                    &mut counts_by_window,
                    original_idx,
                    window_interval,
                    &fragment,
                    count_weight,
                    &motif_context,
                    opt.source_inside,
                    opt.window_assignment.assign_by,
                )?);
            }
        } else {
            // Without genomic scaling, each candidate window gets either weight 1.0 or the raw
            // overlap fraction, depending on the assignment mode.
            for overlapped_window in overlapping_windows.windows {
                let original_idx = window_context.original_idx(overlapped_window.idx);
                let count_weight = match opt.window_assignment.assign_by {
                    WindowMotifAssigner::CountOverlap => overlapped_window.overlap_fraction,
                    _ => 1.0f64,
                } * gc_weight;
                counted_end_flags.merge(count_fragment_in_window(
                    &mut counts_by_window,
                    original_idx,
                    overlapped_window.interval,
                    &fragment,
                    count_weight,
                    &motif_context,
                    opt.source_inside,
                    opt.window_assignment.assign_by,
                )?);
            }
        }

        record_counted_fragment_stats(&mut counter, counted_end_flags);
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    let payload = build_tile_payload(counts_by_window);
    serialize_tile_counts(&counts_path, &payload)?;

    Ok(Some(TileResult {
        chr: tile.chr.clone(),
        counts_path,
        counter,
    }))
}

#[cfg(test)]
mod tests {
    include!("ends_tests.rs");
}
