use crate::{
    commands::{
        cli_common::{
            DistributionWindowSpec, WindowAssigner, ensure_output_dir, load_blacklist_map,
            load_scaling_map, resolve_chromosomes_and_contigs, validate_output_prefix,
        },
        counters::LengthsCounters,
        gc_bias::{
            correct::{LengthAgnosticGCCorrector, load_length_agnostic_gc_corrector},
            counting::build_gc_prefixes,
        },
        lengths::{
            config::{LengthsConfig, validate_gc_length_trim_rare},
            counting::{LengthAxis, LengthCounts},
            tiling::{reduce_partials_for_chr, write_cross_npy, write_partials_npz},
            writer::{LengthCountRowMetadata, write_length_counts_tsv, write_length_settings_json},
        },
        run_statistics::{
            DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions, GCStatisticsSummary,
            print_fragment_run_statistics,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        base::ZEROISH_F64_TOLERANCE,
        bed::{load_grouped_windows_from_bed, load_windows_from_bed},
        blacklist::is_blacklisted,
        clip_mode::ClipMode,
        fragment::indel_counting_fragment::FragmentWithIndelCounts,
        fragment_iterators::fragments_with_indel_counts_from_bam,
        indel_mode::IndelMode,
        interval::{IndexedInterval, Interval},
        io::{FinalOutputFiles, dot_join},
        midpoint::midpoint_random_even_for_fragment,
        overlaps::find_overlapping_windows,
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq_in_range,
        scale_genome::{
            ScalingBin, build_reference_based_scaling_overlaps_for_assignment_overlaps,
            compute_per_window_scaling_over_fragment_for_selected_windows,
            compute_per_window_scaling_over_overlap,
        },
        temp_chrom_names::TempChromNameMap,
        thread_pool::init_global_pool,
        tiled_run::{
            TempDirGuard, Tile, TileWindowSpan, build_tiles, precompute_tile_window_spans,
        },
        window_fetch::{BedFetchPolicy, fetch_span_for_tile},
        windowing::{
            WindowBinInfo, build_bin_info, compute_window_offsets,
            ensure_plain_bed_windows_not_empty,
        },
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{
    path::{Path, PathBuf},
    sync::Arc,
    time::Instant,
};
use tracing::info;

const COMMAND_TARGET: &str = "lengths";

// Map original window index to counts plus containment flag for this tile
#[derive(Clone)]
struct TileCounts {
    counts: LengthCounts,
    contained: bool,
}

#[derive(Clone)]
struct TileOutputs {
    counters: LengthsCounters,
    chr: String,
    global_counts: Option<(String, LengthCounts)>,
    grouped_counts: Option<FxHashMap<u64, LengthCounts>>,
    partial_path: Option<PathBuf>,
    cross_path: Option<PathBuf>,
}

#[inline]
fn should_store_length_weight(weight: f64) -> Result<bool> {
    ensure!(
        weight.is_finite(),
        "length count weight {weight} is not finite"
    );
    ensure!(
        weight >= -ZEROISH_F64_TOLERANCE,
        "length count weight {weight} is negative, this is not currently supported"
    );
    Ok(weight > ZEROISH_F64_TOLERANCE)
}

fn reduce_tile_partials_for_chr(
    chr: &str,
    tile_results: &[TileOutputs],
    n_windows: usize,
    template_counts: &LengthCounts,
) -> Result<Vec<LengthCounts>> {
    let mut partial_paths: Vec<PathBuf> = Vec::new();
    let mut cross_paths: Vec<PathBuf> = Vec::new();

    for tile_out in tile_results.iter().filter(|tile_out| tile_out.chr == chr) {
        if let Some(partial_path) = &tile_out.partial_path {
            partial_paths.push(partial_path.clone());
        }
        if let Some(cross_path) = &tile_out.cross_path {
            cross_paths.push(cross_path.clone());
        }
    }

    reduce_partials_for_chr(
        chr,
        partial_paths.as_slice(),
        cross_paths.as_slice(),
        n_windows,
        template_counts,
    )
}

/// Maximum reference-coordinate reach a counted fragment can require.
///
/// Length bins bound the adjusted output length. When indel adjustment is enabled, deleted
/// reference bases can make the aligned span longer than that adjusted length. When clip
/// adjustment is enabled, the soft-clip cap bounds how far assignment can move outside the
/// aligned reference span.
fn configured_max_fragment_reach_bp(opt: &LengthsConfig, length_axis: &LengthAxis) -> u32 {
    let max_soft_clip_bases = if matches!(opt.clip_mode, ClipMode::Adjust) {
        opt.max_soft_clips as u32
    } else {
        0
    };
    let max_deletion_bases = if matches!(opt.indel_mode, IndelMode::Adjust) {
        opt.max_deletion_bases as u32
    } else {
        0
    };

    length_axis
        .max_fragment_length()
        .saturating_add(max_soft_clip_bases.max(max_deletion_bases))
}

fn reorder_bed_outputs_by_original_index(
    bin_info: &mut Vec<WindowBinInfo>,
    all_bins: &mut Vec<LengthCounts>,
) -> Result<()> {
    ensure!(
        bin_info.len() == all_bins.len(),
        "BED metadata entries ({}) did not match length count vectors ({})",
        bin_info.len(),
        all_bins.len()
    );

    let mut paired: Vec<_> = std::mem::take(bin_info)
        .into_iter()
        .zip(std::mem::take(all_bins))
        .collect();
    paired.sort_unstable_by_key(|(info, _)| info.output_index);

    let (sorted_bin_info, sorted_bins) = paired.into_iter().unzip();
    *bin_info = sorted_bin_info;
    *all_bins = sorted_bins;
    Ok(())
}

/// Execute the fragment length counting pipeline end-to-end.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `lengths` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: &LengthsConfig) -> Result<()> {
    let start_time = Instant::now();
    let length_axis = Arc::new(LengthAxis::new(opt.resolve_length_bins()?)?);
    validate_gc_length_trim_rare(opt.gc_length_trim_rare)?;
    opt.gc.validate(opt.ref_2bit.as_deref())?;
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    let window_opt = opt.windows.resolve_windows();
    let fetch_window_opt = window_opt.as_fetch_window_spec();
    let prefix = opt.output_prefix.trim();
    validate_output_prefix(prefix)?;

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
            let (windows_map, group_idx_to_name, _strand_detection) =
                load_grouped_windows_from_bed(
                    bed,
                    Some(chromosomes.as_slice()),
                    false,
                    None,
                    None,
                )?;
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
                    .map(|(chr, windows)| (chr.clone(), windows.windows_as_slice().to_vec()))
                    .collect()
            })
        };

    // Load genomic scaling factors
    if opt.scale_genome.scaling_factors.is_some() {
        info!(target: COMMAND_TARGET, "Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<ScalingBin>> = load_scaling_map(
        &opt.scale_genome,
        &chromosomes,
        &contigs,
        crate::shared::scale_genome::scaling_gc_mode_for_run(opt.gc.gc_file.is_some(), false),
        None,
    )?;

    // Load GC correction package if specified
    if opt.gc.gc_file.is_some() {
        info!(target: COMMAND_TARGET, "Loading GC correction matrix");
    }
    let gc_corrector = load_length_agnostic_gc_corrector(
        opt.gc.gc_file.as_ref(),
        opt.ref_2bit.as_ref(),
        &opt.gc_length_weighting,
        opt.gc_length_range,
        opt.gc_length_trim_rare,
        length_axis.min_fragment_length(),
        length_axis.max_fragment_length(),
    )?;

    // BAM fetches are in aligned reference coordinates. For deletion-adjusted lengths, the
    // *aligned* reference span can be larger than the maximum fragment length and still be valid,
    // so fetches use the bounded fragment reach rather than only the maximum adjusted output length.
    let max_fragment_reach_bp = configured_max_fragment_reach_bp(opt, length_axis.as_ref());
    // Fragment ownership is based on the aligned fragment start. The left reach only needs to
    // cover assignment coordinates before that owned start, which can only come from left soft
    // clipping. Deletion-adjusted fragments still use the same aligned start, while their allowed
    // reference span to the right is covered by `max_fragment_reach_bp`.
    let left_assignment_reach_bp = if matches!(opt.clip_mode, ClipMode::Adjust) {
        opt.max_soft_clips as u64
    } else {
        0
    };
    let align_bp = match &window_opt {
        DistributionWindowSpec::Size(bp) => Some(*bp),
        _ => None,
    };

    // Build tiles (core plus halo)
    let (tiles, _) = build_tiles(
        &chromosomes,
        &contigs,
        opt.tile_size,
        max_fragment_reach_bp,
        align_bp,
    )?;
    let temp_chrom_name_map = TempChromNameMap::from_contigs(&chromosomes)?;

    let progress = ProgressFactory::new();
    let pb = Arc::new(progress.default_bar(tiles.len() as u64));

    let windows_lookup = indexed_windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(
        &tiles,
        |chr| {
            windows_lookup
                .and_then(|m| m.get(chr).map(|w| w.as_slice()))
                .unwrap_or(&[])
        },
        left_assignment_reach_bp,
        // We use fragments starting in a tile, so we need windows reachable by that assignment
        // interval to the right of the tile as well
        max_fragment_reach_bp as u64,
    ));
    let tile_window_spans_for_threads = tile_window_spans.clone();

    // Reusable length-bin template so every tile/window counter shares identical bounds and avoids repeated allocations
    // Cloned with `zeroed_like` when building per-window `TileCounts`, which guarantees merge compatibility during reduction
    let template_counts = LengthCounts::new(Arc::clone(&length_axis));

    let temp_dir_guard =
        TempDirGuard::new(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;
    let temp_dir = temp_dir_guard.path();
    let mut final_outputs = FinalOutputFiles::new(temp_dir)?;
    let partials_prefix = &dot_join(&[prefix, "part"]);
    let cross_prefix = &dot_join(&[prefix, "cross"]);

    info!(target: COMMAND_TARGET, "Counting per tile");

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    pb.set_position(0);

    let tile_results: Vec<TileOutputs> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<TileOutputs> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let windows_chr: Option<&[IndexedInterval<u64>]> = indexed_windows_map
                .as_ref()
                .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice()));
            let blacklist_chr: &[Interval<u64>] = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr: &[ScalingBin] = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            let counter = process_tile(
                opt,
                tile,
                tile_span.as_ref(),
                windows_chr,
                &window_opt,
                blacklist_chr,
                scaling_chr,
                gc_corrector.clone(),
                &length_axis,
                &template_counts,
                temp_dir,
                partials_prefix,
                cross_prefix,
                &temp_chrom_name_map,
            )?;
            pb.inc(1);
            Ok(counter)
        })
        .collect::<Result<_>>()?; // Short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    // Release per-tile inputs before merging outputs
    drop(tile_window_spans_for_threads);
    drop(tile_window_spans);
    drop(tiles);
    drop(scaling_map);
    drop(gc_corrector);

    // Collect counters
    let mut global_counter = LengthsCounters::default();
    for tile_out in &tile_results {
        global_counter += tile_out.counters;
    }

    match &window_opt {
        DistributionWindowSpec::GroupedBed(_) => {
            info!(target: COMMAND_TARGET, "Merging grouped counts across tiles");
        }
        _ => {
            info!(target: COMMAND_TARGET, "Reducing temporary tile files");
        }
    }

    let mut all_bins: Vec<LengthCounts> = Vec::new();

    match &window_opt {
        DistributionWindowSpec::Global => {
            let mut counts_by_chr: FxHashMap<String, LengthCounts> = FxHashMap::default();
            for tile_out in &tile_results {
                if let Some((chr, counts)) = &tile_out.global_counts {
                    let entry = counts_by_chr
                        .entry(chr.clone())
                        .or_insert_with(|| template_counts.zeroed_like());
                    entry.merge_from(counts)?;
                }
            }
            for chr in &chromosomes {
                let counts = counts_by_chr
                    .remove(chr)
                    .context("Global mode missing counts for chromosome")?;
                all_bins.push(counts);
            }
            all_bins = vec![LengthCounts::collapse(&all_bins)?];
        }
        DistributionWindowSpec::Size(window_bp) => {
            for chr in &chromosomes {
                let chrom_len = contigs
                    .contigs
                    .get(chr)
                    .map(|&(_, len)| len as u64)
                    .context("missing contig length")?;
                let n_windows = chrom_len.div_ceil(*window_bp) as usize;
                let counts = reduce_tile_partials_for_chr(
                    chr,
                    tile_results.as_slice(),
                    n_windows,
                    &template_counts,
                )?;
                ensure!(
                    counts.len() == n_windows,
                    "Expected {} windows for {} but got {}",
                    n_windows,
                    chr,
                    counts.len()
                );
                all_bins.extend(counts);
            }
        }
        DistributionWindowSpec::Bed(_) => {
            let win_map = windows_map
                .as_ref()
                .context("windows_map missing for BED mode")?;
            for chr in &chromosomes {
                let Some(wchr) = win_map.get(chr) else {
                    continue;
                };
                let wchr_slice = wchr.as_slice();
                if wchr_slice.is_empty() {
                    continue;
                }
                let counts = reduce_tile_partials_for_chr(
                    chr,
                    tile_results.as_slice(),
                    wchr_slice.len(),
                    &template_counts,
                )?;
                ensure!(
                    counts.len() == wchr_slice.len(),
                    "Expected {} windows for {} but got {}",
                    wchr_slice.len(),
                    chr,
                    counts.len()
                );
                all_bins.extend(counts);
            }
        }
        DistributionWindowSpec::GroupedBed(_) => {
            let num_groups = group_idx_to_name
                .as_ref()
                .context("group_idx_to_name missing for grouped BED mode")?
                .len();
            all_bins = (0..num_groups)
                .map(|_| template_counts.zeroed_like())
                .collect();
            for tile_out in &tile_results {
                let Some(grouped_counts) = &tile_out.grouped_counts else {
                    continue;
                };
                for (group_idx, counts) in grouped_counts {
                    let entry = all_bins.get_mut(*group_idx as usize).with_context(|| {
                        format!("group index {group_idx} outside allocated grouped output")
                    })?;
                    entry.merge_from(counts)?;
                }
            }
        }
    }

    let mut bin_info: Vec<WindowBinInfo> = match &window_opt {
        DistributionWindowSpec::Global | DistributionWindowSpec::GroupedBed(_) => Vec::new(),
        _ => {
            let (_total_windows, chr_offsets) = compute_window_offsets(
                &fetch_window_opt,
                &chromosomes,
                &contigs,
                windows_map.as_ref(),
            )?;
            build_bin_info(
                &fetch_window_opt,
                &chromosomes,
                &contigs,
                windows_map.as_ref(),
                &blacklist_map,
                &chr_offsets,
            )?
        }
    };

    // Sort by original index (when given a bed file)
    if matches!(&window_opt, DistributionWindowSpec::Bed(_)) {
        info!(
            target: COMMAND_TARGET,
            "Reordering counts by original window index in BED file"
        );

        reorder_bed_outputs_by_original_index(&mut bin_info, &mut all_bins)?;
    }

    info!(target: COMMAND_TARGET, "Writing final files");

    // Write every final output to the temp directory before moving any of them into place
    // This keeps failed writes from leaving a mix of old and new final files
    let final_counts_path = opt
        .ioc
        .output_dir
        .join(dot_join(&[prefix, "length_counts.tsv.zst"]));
    let temp_counts_path = final_outputs.temp_path_for(&final_counts_path)?;
    let row_metadata = match &window_opt {
        DistributionWindowSpec::Global => LengthCountRowMetadata::Global,
        DistributionWindowSpec::Size(_) | DistributionWindowSpec::Bed(_) => {
            LengthCountRowMetadata::Windows {
                windows: &bin_info,
                include_blacklisted_fraction: opt.blacklist.is_some(),
            }
        }
        DistributionWindowSpec::GroupedBed(_) => {
            let group_idx_to_name = group_idx_to_name
                .as_ref()
                .context("group_idx_to_name missing when writing grouped outputs")?;
            let grouped_windows_map = grouped_windows_map
                .as_ref()
                .context("grouped windows missing when writing grouped outputs")?;
            LengthCountRowMetadata::Groups {
                group_idx_to_name,
                chromosomes: &chromosomes,
                grouped_windows_map,
                blacklist_map: &blacklist_map,
                include_blacklisted_fraction: opt.blacklist.is_some(),
            }
        }
    };
    write_length_counts_tsv(
        &temp_counts_path,
        &all_bins,
        length_axis.as_ref(),
        opt.decimals,
        row_metadata,
    )
    .with_context(|| {
        format!(
            "writing final length counts to temp file {}",
            temp_counts_path.display()
        )
    })?;
    final_outputs.record(temp_counts_path, final_counts_path)?;

    let settings_path = opt
        .ioc
        .output_dir
        .join(dot_join(&[prefix, "length_settings.json"]));
    let temp_settings_path = final_outputs.temp_path_for(&settings_path)?;
    write_length_settings_json(&temp_settings_path, opt, &window_opt, &length_axis)?;
    final_outputs.record(temp_settings_path, settings_path)?;

    final_outputs.move_into_place()?;

    // Plot the global fragment length distribution after the machine-readable outputs are complete
    #[cfg(feature = "plotters")]
    {
        info!(target: COMMAND_TARGET, "Plotting overall length distribution");

        use crate::shared::plotters::lineplot::write_line_plot_png;

        if all_bins.is_empty() {
            info!(
                target: COMMAND_TARGET,
                "Skipping overall length plot because no bins were produced"
            );
        } else {
            let mut global_counts = vec![0f64; all_bins[0].counts.len()];
            for length_counts in &all_bins {
                for (total, count) in global_counts.iter_mut().zip(length_counts.counts.iter()) {
                    *total += *count;
                }
            }

            let total_counts: f64 = global_counts.iter().sum();
            if total_counts > 0.0 {
                for value in &mut global_counts {
                    *value /= total_counts;
                }
            }

            let x_values: Vec<f64> = length_axis
                .edges()
                .windows(2)
                .map(|edge_pair| (edge_pair[0] as f64 + edge_pair[1] as f64 - 1.0) / 2.0)
                .collect();
            let x_label = if length_axis.is_single_bp_bins() {
                "Fragment length (bp)"
            } else {
                "Fragment length bin midpoint (bp)"
            };
            let y_label = if length_axis.is_single_bp_bins() {
                "Density"
            } else {
                "Fraction of counted mass per bin"
            };

            let plot_path = opt
                .ioc
                .output_dir
                .join(dot_join(&[prefix, "fragment_lengths_overall.png"]));

            write_line_plot_png(
                &plot_path,
                "Fragment length distribution (summed/global)",
                x_label,
                y_label,
                &x_values,
                &global_counts,
                1600,
                1000,
            )
            .with_context(|| format!("writing fragment length plot to {}", plot_path.display()))?;
        }
    }

    drop(blacklist_map);

    let elapsed = start_time.elapsed();
    print_fragment_run_statistics(
        &global_counter.base,
        elapsed,
        FragmentRunStatisticsOptions {
            include_section_header: true,
            notes: &[],
            labels: DEFAULT_FRAGMENT_STATISTICS_LABELS,
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

fn process_tile(
    opt: &LengthsConfig,
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    windows_chr: Option<&[IndexedInterval<u64>]>,
    window_opt: &DistributionWindowSpec,
    blacklist_intervals: &[Interval<u64>],
    scaling_chr: &[ScalingBin],
    gc_corrector_opt: Option<LengthAgnosticGCCorrector>,
    length_axis: &Arc<LengthAxis>,
    template: &LengthCounts,
    temp_dir: &Path,
    partials_prefix: &str,
    cross_prefix: &str,
    temp_chrom_name_map: &TempChromNameMap,
) -> Result<TileOutputs> {
    let fetch_window_opt = window_opt.as_fetch_window_spec();
    let max_fragment_reach_bp = configured_max_fragment_reach_bp(opt, length_axis.as_ref()) as u64;
    // One BAM reader per tile
    let (mut reader, _tid_check, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;
    tile.ensure_matches_bam_tid(_tid_check)?;

    // Counters
    let mut counter = LengthsCounters::default();

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

    // Adapt reference-coordinate fetch bounds to the present windows using the same bounded
    // reference reach used for tile halos.
    let Some(fetch_span) = fetch_span_for_tile(
        tile,
        tile_window_span,
        windows_chr,
        &fetch_window_opt,
        chrom_len,
        max_fragment_reach_bp,
        BedFetchPolicy::CandidateWindowExtent,
    )?
    else {
        // Skip tiles with no relevant windows
        return Ok(TileOutputs {
            counters: counter,
            chr: tile.chr.clone(),
            global_counts: None,
            grouped_counts: None,
            partial_path: None,
            cross_path: None,
        });
    };
    let (fetch_from, fetch_to) = fetch_span.try_to_i64()?.as_tuple();

    reader
        .fetch((tile.tid, fetch_from, fetch_to))
        .context(format!("fetch {} {}-{}", &tile.chr, fetch_from, fetch_to))?;

    let left_assignment_reach_bp = if matches!(opt.clip_mode, ClipMode::Adjust) {
        opt.max_soft_clips as u64
    } else {
        0
    };
    let leftmost_reachable_start =
        (tile.core_start() as u64).saturating_sub(left_assignment_reach_bp);

    // Preallocate per-tile window counters
    // Keep indices aligned with global scan order so downstream merging works without remapping
    // Use Option to skip BED windows that cannot be hit by any fragment starting in this tile
    // Track the first and last index covered to translate back to global coordinates
    let (counts_start_idx, counts_end_idx_exclusive, mut counts_by_idx): (
        usize,
        usize,
        Vec<Option<TileCounts>>,
    ) = match window_opt {
        // Global mode has exactly one window covering the chromosome
        DistributionWindowSpec::Global => (
            0,
            1,
            vec![Some(TileCounts {
                counts: template.zeroed_like(),
                contained: false,
            })],
        ),

        // Fixed-size mode: allocate only the bins that a fragment from this tile can reach
        DistributionWindowSpec::Size(window_bp) => {
            // Total bins on the chromosome
            let chrom_bin_count = chrom_len.div_ceil(*window_bp) as usize;
            // Leftmost bin whose start is at or before the core start
            // (may begin before the core when cores are not aligned)
            let min_bin_idx = (leftmost_reachable_start / *window_bp) as usize;
            // Furthest coordinate a fragment starting in this tile can reach
            let max_reachable_end = (tile.core_end() as u64)
                .saturating_add(max_fragment_reach_bp)
                .min(chrom_len);
            // One past the last bin that could overlap that reach
            let max_bin_idx_exclusive = if max_reachable_end == 0 {
                min_bin_idx
            } else {
                (((max_reachable_end - 1) / *window_bp) + 1) as usize
            }
            .min(chrom_bin_count);

            let span_len = max_bin_idx_exclusive.saturating_sub(min_bin_idx);
            let mut counts = Vec::with_capacity(span_len);
            for idx in min_bin_idx..max_bin_idx_exclusive {
                let start = idx as u64 * *window_bp;
                let end = (start + *window_bp).min(chrom_len);
                // Contained means the bin sits fully inside the tile core
                let contained = start >= tile.core_start() as u64 && end <= tile.core_end() as u64;
                counts.push(Some(TileCounts {
                    counts: template.zeroed_like(),
                    contained,
                }));
            }
            (min_bin_idx, max_bin_idx_exclusive, counts)
        }

        // BED mode: reuse the precomputed span and skip only windows that still sit fully outside
        // the clip-adjusted left reach
        DistributionWindowSpec::Bed(_) => {
            let span = tile_window_span.context(
                "BED length counting requires a cached tile window span after fetch-span selection",
            )?;
            let wchr = windows_chr.context("BED length counting requires loaded windows")?;
            let span_len = span.last_idx_exclusive.saturating_sub(span.first_idx);
            let mut counts = Vec::with_capacity(span_len);
            for idx in span.first_idx..span.last_idx_exclusive {
                let window = wchr[idx];
                let win_start = window.start();
                let win_end = window.end();
                // Windows fully to the left of the furthest clip-adjusted start cannot be hit
                // because every counted fragment starts inside the core and can only extend left
                // by the configured soft-clip reach
                if win_end <= leftmost_reachable_start {
                    counts.push(None);
                    continue;
                }
                // Contained flags windows fully inside the core
                let contained =
                    win_start >= tile.core_start() as u64 && win_end <= tile.core_end() as u64;
                counts.push(Some(TileCounts {
                    counts: template.zeroed_like(),
                    contained,
                }));
            }
            (span.first_idx, span.last_idx_exclusive, counts)
        }
        DistributionWindowSpec::GroupedBed(_) => (0, 0, Vec::new()),
    };
    let mut counts_by_group: FxHashMap<u64, LengthCounts> = FxHashMap::default();

    // Fraction of a fragment that must overlap with a window to assign to that window
    let min_overlap_fraction: f64 = match opt.window_assignment.assign_by {
        WindowAssigner::Any | WindowAssigner::CountOverlap => {
            1. / (max_fragment_reach_bp as f64 + 1.0)
        } // +1 to avoid rounding error issues
        WindowAssigner::All | WindowAssigner::Midpoint => {
            1.0 - (1. / (max_fragment_reach_bp as f64 + 1.0))
        } // 1.0 but just below to avoid rounding errors
        WindowAssigner::Proportion(p) => p,
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
        .map(|b| IndexedInterval::from_interval(b.interval, 0_u64))
        .collect();

    // Function for filtering fragments after pairing
    let fragment_filter = {
        let length_axis = Arc::clone(length_axis);
        let indel_mode = opt.indel_mode;
        let clip_mode = opt.clip_mode;
        let max_soft_clips = opt.max_soft_clips as u32;
        let max_deletion_bases = opt.max_deletion_bases as u32;
        move |fragment: &FragmentWithIndelCounts| {
            if matches!(clip_mode, ClipMode::Adjust)
                && !fragment.soft_clips_within_limit(max_soft_clips)
            {
                return false;
            }
            if matches!(indel_mode, IndelMode::Adjust)
                && !fragment.deletion_bases_within_limit(max_deletion_bases)
            {
                return false;
            }
            if matches!(clip_mode, ClipMode::Skip) && fragment.has_soft_clipping() {
                return false;
            }
            length_axis.contains(fragment.adjusted_len(indel_mode, clip_mode))
        }
    };

    // Create fragment iterator with per-tile filtering and optional GC tag handling
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
    let mut iter = fragments_with_indel_counts_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        move |rec| include_read_fn(rec),
        opt.indel_mode,
        lengths_needs_cigar_inspection(opt),
        fragment_filter,
        unpaired,
    )
    .with_local_counters();

    let get_gc_weight = {
        let gc_corrector = gc_corrector_opt.as_ref();
        let gc_prefixes = gc_prefixes_opt.as_ref();
        let fetch_start = tile.fetch_start();
        move |fragment: &FragmentWithIndelCounts| -> Result<Option<f64>> {
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

        let aligned_fragment_interval = fragment.interval.try_to_u64()?;

        // Determine blacklist status
        let in_blacklist = is_blacklisted(
            blacklist_intervals,
            opt.blacklist_strategy,
            aligned_fragment_interval,
            max_fragment_reach_bp,
            &mut bl_ptr,
        );
        if in_blacklist {
            counter.blacklisted_fragments += 1;
            continue;
        }

        // Calculate fragment length and window-assignment interval
        // GC, blacklist, and scaling still use the aligned reference span
        let fragment_length = fragment.adjusted_len(opt.indel_mode, opt.clip_mode);
        let assignment_interval = fragment.assignment_interval_with_clip_mode(opt.clip_mode)?;

        // Find candidate count-windows from the interval implied by window assignment.
        // Midpoint mode queries only the 1 bp midpoint, while the other modes query the
        // full assignment interval.
        let window_selection_interval = match opt.window_assignment.assign_by {
            WindowAssigner::Midpoint => {
                let midpoint = midpoint_random_even_for_fragment(
                    &tile.chr,
                    assignment_interval.start() as u32,
                    fragment_length,
                );
                Interval::new(midpoint.into(), (midpoint + 1).into())?
            }
            WindowAssigner::Any
            | WindowAssigner::All
            | WindowAssigner::Proportion(_)
            | WindowAssigner::CountOverlap => assignment_interval,
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
            window_selection_interval,
            min_overlap_fraction,
            max_fragment_reach_bp,
        )?;
        let overlapping_windows = if let Some(overlaps) = overlapping_windows {
            overlaps
        } else {
            continue;
        };

        // Get GC correction weight
        let gc_weight_opt = get_gc_weight(&fragment)?;
        let gc_weight = match (gc_weight_opt, correct_gc) {
            (Some(w), true) => w,
            (None, true) => {
                // Tried but failed to make a GC correction weight for the current fragment
                // Fall back to no correction or skip
                counter.gc_failed_fragments += 1;
                if opt.gc.neutralize_invalid_gc {
                    1.0
                } else {
                    continue;
                }
            }
            (None, false) => 1.0, // No correction
            (Some(_), false) => bail!("unexpected GC weight when GC correction is disabled"),
        };

        counter.base.counted_fragments += 1;

        if !scaling_chr.is_empty() {
            // Find overlapping scaling-bins
            let overlapping_scaling_bins = find_overlapping_windows(
                chrom_len,
                &mut sf_ptr,
                Some(&scaling_with_bin_idx),
                None,
                aligned_fragment_interval, // Full aligned fragment
                1. / (max_fragment_reach_bp as f64 + 1.0), // Any overlap
                max_fragment_reach_bp,
            )
            .with_context(|| format!("finding overlapping scaling bins on chr {}", tile.chr))?
            .context("no overlapping scaling bins found")?; // Should always find >= 1 bin

            // Extract the indices of the overlapping bins
            let overlapping_scaling_bin_indices: Vec<usize> = overlapping_scaling_bins
                .windows
                .iter()
                .map(|w| w.idx)
                .collect();

            // Calculate the weight per overlapping count-window
            // NOTE: `compute_per_window_scaling_over_fragment_for_selected_windows`
            // always returns an overlap fraction of 1.0 (count full fragment)!
            let overlap_weights = match opt.window_assignment.assign_by {
                WindowAssigner::CountOverlap => {
                    if matches!(opt.clip_mode, ClipMode::Adjust) {
                        let scaling_overlaps =
                            build_reference_based_scaling_overlaps_for_assignment_overlaps(
                                &overlapping_windows,
                                aligned_fragment_interval,
                            )?;
                        compute_per_window_scaling_over_overlap(
                            &overlapping_windows,
                            Some(&scaling_overlaps),
                            &overlapping_scaling_bin_indices,
                            scaling_chr,
                        )?
                    } else {
                        compute_per_window_scaling_over_overlap(
                            &overlapping_windows,
                            None,
                            &overlapping_scaling_bin_indices,
                            scaling_chr,
                        )?
                    }
                }
                _ => compute_per_window_scaling_over_fragment_for_selected_windows(
                    aligned_fragment_interval,
                    &overlapping_windows,
                    &overlapping_scaling_bin_indices,
                    scaling_chr,
                )?,
            };

            // Count up the weight per overlapping count-window
            for window_scaling in overlap_weights {
                let overlapped_window_idx = window_scaling.window_idx;
                let count_weight = window_scaling.overlap_fraction_to_count
                    * window_scaling.scaling_weight
                    * gc_weight;
                if !should_store_length_weight(count_weight)? {
                    continue;
                }
                match window_opt {
                    DistributionWindowSpec::GroupedBed(_) => {
                        let windows_chr = windows_chr
                            .context("grouped BED length counting requires loaded windows")?;
                        let group_idx = windows_chr
                            .get(overlapped_window_idx)
                            .with_context(|| {
                                format!(
                                    "missing grouped window {} in chromosome-local window slice for {}",
                                    overlapped_window_idx, tile.chr
                                )
                            })?
                            .idx();
                        counts_by_group
                            .entry(group_idx)
                            .or_insert_with(|| template.zeroed_like())
                            .incr_weighted(fragment_length as usize, count_weight)?;
                    }
                    _ => {
                        let vec_idx = overlapped_window_idx - counts_start_idx;
                        if vec_idx >= counts_by_idx.len() {
                            bail!(
                                "Overlapping window idx {} outside [{}..{}) on {}",
                                overlapped_window_idx,
                                counts_start_idx,
                                counts_end_idx_exclusive,
                                tile.chr
                            );
                        }
                        if let Some(entry) = counts_by_idx[vec_idx].as_mut() {
                            entry
                                .counts
                                .incr_weighted(fragment_length as usize, count_weight)?;
                        }
                    }
                }
            }
        } else {
            // When no scaling, increment counter by 1.0 or by the overlap fraction
            for overlapped_window in overlapping_windows.windows {
                let count_weight = match opt.window_assignment.assign_by {
                    WindowAssigner::CountOverlap => overlapped_window.overlap_fraction,
                    _ => 1.0f64,
                } * gc_weight;
                if !should_store_length_weight(count_weight)? {
                    continue;
                }
                match window_opt {
                    DistributionWindowSpec::GroupedBed(_) => {
                        let windows_chr = windows_chr
                            .context("grouped BED length counting requires loaded windows")?;
                        let group_idx = windows_chr
                            .get(overlapped_window.idx)
                            .with_context(|| {
                                format!(
                                    "missing grouped window {} in chromosome-local window slice for {}",
                                    overlapped_window.idx, tile.chr
                                )
                            })?
                            .idx();
                        counts_by_group
                            .entry(group_idx)
                            .or_insert_with(|| template.zeroed_like())
                            .incr_weighted(fragment_length as usize, count_weight)?;
                    }
                    _ => {
                        let vec_idx = overlapped_window.idx - counts_start_idx;
                        if vec_idx >= counts_by_idx.len() {
                            bail!(
                                "Overlapping window idx {} outside [{}..{}) on {}",
                                overlapped_window.idx,
                                counts_start_idx,
                                counts_end_idx_exclusive,
                                tile.chr
                            );
                        }
                        if let Some(entry) = counts_by_idx[vec_idx].as_mut() {
                            entry
                                .counts
                                .incr_weighted(fragment_length as usize, count_weight)?;
                        }
                    }
                }
            }
        }
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    // Prepare outputs
    let mut window_idxs_chr: Vec<u64> = Vec::with_capacity(counts_by_idx.len());
    let mut counts: Vec<LengthCounts> = Vec::with_capacity(counts_by_idx.len());
    let mut contained_flags: Vec<bool> = Vec::with_capacity(counts_by_idx.len());
    let mut crossing_window_idxs_chr: Vec<u64> = Vec::new();
    for (offset, tile_counts_opt) in counts_by_idx.into_iter().enumerate() {
        if let Some(tile_counts) = tile_counts_opt {
            let idx = (counts_start_idx + offset) as u64;
            window_idxs_chr.push(idx);
            counts.push(tile_counts.counts);
            contained_flags.push(tile_counts.contained);
            // Aligned tile/window boundaries do not make non-contained bins tile-exclusive.
            // A fragment can start near the right edge of this core and still contribute to a
            // downstream bin that is fully contained in the next tile. Every non-contained row
            // must therefore stay in the cross-index so the reducer knows to merge it.
            if !tile_counts.contained {
                crossing_window_idxs_chr.push(idx);
            }
        }
    }

    if matches!(window_opt, DistributionWindowSpec::Global) {
        debug_assert_eq!(counts.len(), 1);
        let chr_counts = counts
            .into_iter()
            .next()
            .unwrap_or_else(|| template.zeroed_like());
        return Ok(TileOutputs {
            counters: counter,
            chr: tile.chr.clone(),
            global_counts: Some((tile.chr.clone(), chr_counts)),
            grouped_counts: None,
            partial_path: None,
            cross_path: None,
        });
    }

    if matches!(window_opt, DistributionWindowSpec::GroupedBed(_)) {
        return Ok(TileOutputs {
            counters: counter,
            chr: tile.chr.clone(),
            global_counts: None,
            grouped_counts: Some(counts_by_group),
            partial_path: None,
            cross_path: None,
        });
    }

    let partial_path = write_partials_npz(
        temp_dir,
        partials_prefix,
        &tile.chr,
        tile.index,
        temp_chrom_name_map,
        &window_idxs_chr,
        &contained_flags,
        &counts,
    )?;
    let cross_path = write_cross_npy(
        temp_dir,
        cross_prefix,
        &tile.chr,
        tile.index,
        temp_chrom_name_map,
        &crossing_window_idxs_chr,
    )?;

    Ok(TileOutputs {
        counters: counter,
        chr: tile.chr.clone(),
        global_counts: None,
        grouped_counts: None,
        partial_path,
        cross_path,
    })
}

#[inline]
fn lengths_needs_cigar_inspection(opt: &LengthsConfig) -> bool {
    !matches!(opt.indel_mode, IndelMode::Ignore) || !matches!(opt.clip_mode, ClipMode::Aligned)
}

#[cfg(test)]
mod tests {
    include!("lengths_tests.rs");
}
