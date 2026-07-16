use crate::{
    command_run::{CommandRunResult, RunOptions, status_info},
    commands::{
        cli_common::*,
        counters::GCCounters,
        gc_bias::{
            CORRECTION_CLAMP_RANGE,
            binning::{CollapseAggregation, bin_greedily_by_mass, collapse_counts_by_bins},
            config::GCConfig,
            counting::{
                GCCounts, GCPrefixes, apply_gc_percent_width_correction, build_gc_prefixes,
            },
            cross_tile_parts::{CrossingPart, stream_crossing_files, write_crossing_parts},
            interpolation::fill_unsupported_bins_with_polynomial,
            load_reference_bias::{ReferenceGCData, ReferenceGCMetadata, load_reference_gc_data},
            outliers::{OutlierRule, OutlierStats, apply_outliers_to_matrix},
            package::GCCorrectionPackage,
            smoothing::smoothe_counts_gaussian,
            support_masking::{
                build_extreme_bins_support_mask, set_masked_entries_to_value, stats_by_support_mask,
            },
            windows::{
                WindowState, advance_fixed_size_streaming_buffers, compute_window_stats,
                overlap_length, prepare_tile_windows, set_window_acgt_in_observed_interval,
            },
        },
        run_statistics::{
            DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions,
            print_fragment_run_statistics,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::apply_blacklist_mask_to_seq,
        constants::{GC_CORRECTION_SCHEMA_VERSION, MIN_ACGT_BASES_FOR_GC_FRACTION},
        fragment::minimal_fragment::Fragment,
        fragment_iterators::fragments_from_bam,
        interval::{IndexedInterval, Interval},
        io::{FinalOutputFiles, dot_join},
        midpoint::midpoint_random_even_for_fragment,
        overlaps::find_overlapping_windows,
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq_in_range,
        reference::{ContigFootprintEntry, twobit_contig_footprint},
        thread_pool::init_global_pool,
        tiled_run::{
            TempDirGuard, Tile, TileWindowSpan, build_tiles, precompute_tile_window_spans,
        },
        windowing::compute_window_offsets,
        windowing::ensure_plain_bed_windows_not_empty,
    },
};
use anyhow::{Context, Result, anyhow, bail, ensure};
use fxhash::FxHashSet;
use ndarray::{Array1, Array2, ArrayBase, Axis, Data, DataMut, Ix2, Zip};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{
    fs::create_dir_all,
    path::{Path, PathBuf},
    sync::Arc,
    time::Instant,
};
use tracing::{info, warn};

const COMMAND_TARGET: &str = "gc-bias";

/// Result from `gc-bias`.
///
/// The command writes a GC correction package. The result records the package path, final output
/// files, and counters from the observed-fragment run.
#[derive(Debug)]
pub struct GCBiasRunResult {
    /// Fragment counters collected while building the observed GC distribution.
    pub counters: GCCounters,
    /// Final GC correction package path.
    pub correction_package_path: PathBuf,
    /// Final output files produced by the command.
    pub output_files: Vec<PathBuf>,
}

impl CommandRunResult for GCBiasRunResult {
    type Counters = GCCounters;

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&Path> {
        Some(self.correction_package_path.as_path())
    }
}

/// Get the GC count for a fragment after trimming ends.
///
/// Returns `None` when the trimmed interval is empty, outside the loaded sequence,
/// or lacks sufficient ACGT support.
pub(crate) fn get_fragment_gc(
    fragment_interval: Interval<u64>,
    sequence_interval: Interval<u64>,
    end_offset: u32,
    gc_prefixes: &GCPrefixes,
    min_acgt_fraction: f32,
) -> Result<Option<usize>> {
    let gc_window = match fragment_interval.contract(end_offset as u64) {
        Some(interval) => interval,
        None => return Ok(None),
    };

    if !sequence_interval.contains_interval(gc_window) {
        return Ok(None);
    }
    let gc_window_local = gc_window
        .shift_left(sequence_interval.start())?
        .try_to_usize()?;
    let acgt = gc_prefixes.acgt_count(gc_window_local)?;
    if acgt < MIN_ACGT_BASES_FOR_GC_FRACTION
        || (acgt as f32 / gc_window.len() as f32) < min_acgt_fraction
    {
        return Ok(None);
    }

    let gc = gc_prefixes.gc_count(gc_window_local)?;
    ensure!(
        gc <= gc_window.len() as u32,
        "GC count exceeded interval length: {} > {}",
        gc,
        gc_window.len()
    );
    Ok(Some(gc as usize))
}

/// Build the interval used for window assignment for one fragment.
///
/// GC bias counts use the same checked interval semantics for fixed-size and
/// non-streaming windows. Keeping the assignment span as an `Interval` avoids
/// open-coding start/end math in each path.
fn fragment_assignment_interval(
    chromosome: &str,
    fragment: &Fragment,
    assign_by: WindowAssigner,
) -> Result<Interval<u64>> {
    match assign_by {
        WindowAssigner::Midpoint => {
            let midpoint =
                midpoint_random_even_for_fragment(chromosome, fragment.start(), fragment.len());
            Ok(Interval::new(midpoint.into(), (midpoint + 1).into())?)
        }
        WindowAssigner::Any
        | WindowAssigner::All
        | WindowAssigner::Proportion(_)
        | WindowAssigner::CountOverlap => Ok(fragment.interval.try_to_u64()?),
    }
}

pub(crate) fn finalize_window_buffer(
    buf: &mut WindowState,
    gc_prefixes: &GCPrefixes,
    sequence_interval: Interval<u64>,
    tile_core_interval: Interval<u64>,
    windows_aligned_to_tiles: bool,
    apply_window_scaling: bool,
    opt: &GCConfig,
    avg_window_span: f64,
    // Tile-level aggregate accumulator. This is not the current window being finalized. It
    // collects the accepted per-window contributions for the whole tile and later carries the
    // optional crossing-file path back to the reducer
    tile_output: &mut WindowState,
    crossing_parts: &mut Vec<CrossingPart>,
    crossing_window_index_offset: u64,
) -> Result<()> {
    let crosses_tile = !tile_core_interval.contains_interval(buf.interval);
    // For a window that crosses tile boundaries, this is the part of the window owned by the
    // current tile. Support must be counted from the tile core, not from the full fetched span,
    // otherwise neighbouring tiles' fetch halos would double count the same bases.
    let tile_owned_interval = if crosses_tile {
        buf.interval.intersection(tile_core_interval)
    } else {
        None
    };
    let should_spill_crossing = apply_window_scaling
        && !windows_aligned_to_tiles
        && crosses_tile
        && (buf.has_counts || tile_owned_interval.is_some());

    if !buf.has_counts && !should_spill_crossing {
        buf.counts.clear();
        return Ok(());
    }

    if !apply_window_scaling {
        tile_output.counts.merge_from(&buf.counts)?;
        tile_output.weight += 1;
    } else if should_spill_crossing {
        // Fixed-size streaming can finalize a synthetic "next" window that lies wholly outside
        // the current tile core. Fragments starting in this tile may still contribute counts to
        // that window, but this tile owns none of the window's support span. In that case we must
        // spill the counts with zero observed support here and let the owning tile contribute the
        // support bases later. Empty placeholders with neither counts nor owned support are
        // filtered out by the `should_spill_crossing` guard above.
        if let Some(tile_owned_interval) = tile_owned_interval {
            set_window_acgt_in_observed_interval(
                buf,
                gc_prefixes,
                tile_owned_interval,
                sequence_interval,
            )?;
        } else {
            buf.counts.num_acgt_out_of = (0, 0);
        }
        let crossing_window_idx = buf
            .idx
            .checked_add(crossing_window_index_offset)
            .context("crossing window index overflowed")?;
        crossing_parts.push(CrossingPart {
            idx: usize::try_from(crossing_window_idx)
                .context("crossing window index exceeded usize range")?,
            counts: buf.counts.clone(),
        });
    } else {
        set_window_acgt_in_observed_interval(buf, gc_prefixes, buf.interval, sequence_interval)?;
        if process_window_in_place(&mut buf.counts, opt, avg_window_span)? {
            tile_output.counts.merge_from(&buf.counts)?;
            tile_output.weight += 1;
        }
    }
    buf.counts.clear();
    buf.has_counts = false;
    Ok(())
}

struct ReduceState {
    scaled_sum: GCCounts,
    scaled_weight: usize,
    crossing_files: Vec<PathBuf>,
    counters: GCCounters,
}

impl ReduceState {
    fn from_scaled_sum(scaled_sum: GCCounts) -> Self {
        Self {
            scaled_sum,
            scaled_weight: 0,
            crossing_files: Vec::new(),
            counters: GCCounters::default(),
        }
    }

    fn merge_scaled(&mut self, other: &GCCounts) -> Result<()> {
        self.scaled_sum.merge_from(other)?;
        Ok(())
    }

    fn add_weight(&mut self, w: usize) {
        self.scaled_weight += w;
    }

    fn merge_counters(&mut self, other: GCCounters) {
        self.counters += other;
    }

    fn push_crossing_file(&mut self, path: Option<PathBuf>) {
        if let Some(p) = path {
            self.crossing_files.push(p);
        }
    }

    fn merge(mut self, mut other: Self) -> Result<Self> {
        self.merge_scaled(&other.scaled_sum)?;
        self.scaled_weight += other.scaled_weight;
        self.merge_counters(other.counters);
        self.crossing_files.append(&mut other.crossing_files);

        Ok(self)
    }
}

/// Run the `gc-bias` command.
///
/// This command combines an observed fragment GC distribution with a reference GC package and
/// writes a correction package for later commands. It applies configured fragment filters,
/// windowing, support masks, outlier handling, smoothing, and interpolation.
///
/// Reporting is controlled by `options`. `report_statistics` prints the final summary,
/// `show_progress` controls progress bars, and `log_statuses` controls status messages.
///
/// Parameters
/// ----------
/// - `opt`:
///   Fully resolved configuration for the `gc-bias` command.
/// - `options`:
///   Reporting controls for statistics, progress bars, and status logs.
///
/// Returns
/// -------
/// - `Ok(GCBiasRunResult)`:
///   Counters and output paths for the completed run.
///
/// Errors
/// ------
/// Returns an error when the reference package is incompatible, the configuration is invalid, an
/// input cannot be read, or the correction package cannot be written.
pub fn run_gc_bias(opt: &GCConfig, options: RunOptions) -> Result<GCBiasRunResult> {
    let start_time = Instant::now();
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let prefix = opt.output_prefix.trim();
    validate_output_prefix(prefix)?;
    let window_opt = opt.windows.resolve_windows();
    if options.log_equivalent_cli {
        let command = crate::ToCliCommand::to_cli_string(opt)?;
        let message = crate::command_run::equivalent_cli_log_message(&command);
        info!(target: COMMAND_TARGET, "{message}");
    }
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;

    status_info!(options, target: COMMAND_TARGET, "Loading reference GC bias");
    let ReferenceGCData {
        counts: reference_counts,
        unobservables_support_mask: reference_unobservables_support_mask,
        outliers_support_mask: mut reference_outliers_support_mask,
        gc_percent_widths: reference_gc_percent_widths,
        metadata: reference_metadata,
    } = load_reference_gc_data(&opt.ref_gc_file)?;
    validate_reference_chromosomes(&reference_metadata.chromosomes, &chromosomes)?;
    validate_fixed_size_window_for_reference(&window_opt, &reference_metadata)?;
    validate_reference_contig_footprint(
        &reference_metadata,
        &twobit_contig_footprint(&opt.ref_genome.ref_2bit)?,
    )?;
    let avg_norm_ref_counts = mean_scale_per_length_array(
        &reference_counts,
        0.,
        Some(&reference_outliers_support_mask),
    );
    drop(reference_counts);

    // Create output directory
    create_dir_all(&opt.ioc.output_dir).context("Cannot create output_dir")?;
    let final_temp_dir_guard = TempDirGuard::new(&opt.ioc.output_dir, "gc_bias_final")
        .context("create final output temp dir")?;
    let mut final_outputs = FinalOutputFiles::new(final_temp_dir_guard.path())?;
    let mut intermediate_saver = IntermediateFileSaver::new(
        opt.save_intermediates,
        final_outputs.temp_dir().to_path_buf(),
        opt.ioc.output_dir.clone(),
        prefix.to_string(),
        options.log_statuses,
    );

    // Load blacklist intervals if provided
    let blacklist_map = load_blacklist_map(
        opt.blacklist.as_ref(),
        1,
        0,
        &chromosomes,
        opt.ioc.n_threads > 1,
    )?;

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            status_info!(options, target: COMMAND_TARGET, "Loading window coordinates");
            let windows = load_windows_from_bed(
                bed,
                Some(chromosomes.as_slice()),
                None,
                None,
                opt.ioc.n_threads > 1,
            )?;
            ensure_plain_bed_windows_not_empty(&windows)?;
            Some(windows)
        }
        _ => None,
    };

    let window_stats =
        compute_window_stats(&window_opt, windows_map.as_ref(), &contigs, &chromosomes)?;
    let avg_window_span = window_stats.avg_span;
    let total_windows = window_stats.total_windows;
    let (_, chromosome_window_offsets) =
        compute_window_offsets(&window_opt, &chromosomes, &contigs, windows_map.as_ref())?;

    // Build tiles (core plus halo = max fragment length) to bound memory per worker
    let halo_bp = reference_metadata.max_fragment_length as u32;
    let align_bp = match &window_opt {
        WindowSpec::Size(bp) => Some(*bp),
        _ => None,
    };
    let (tiles, windows_aligned_to_tiles) =
        build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, align_bp)?;
    let progress = ProgressFactory::with_enabled(options.show_progress);
    let pb = Arc::new(progress.default_bar(tiles.len() as u64));

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(
        &tiles,
        |chr| {
            windows_lookup
                .and_then(|m| m.get(chr).map(|w| w.as_slice()))
                .unwrap_or(&[])
        },
        0,
        reference_metadata.max_fragment_length as u64,
    ));
    let tile_window_spans_for_threads = tile_window_spans.clone();

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    status_info!(options, target: COMMAND_TARGET, "Counting per tile");

    pb.set_position(0);

    let zero_counts = GCCounts::new(
        reference_metadata.min_fragment_length,
        reference_metadata.max_fragment_length,
        reference_metadata.end_offset as usize,
        (0, 0),
    )?;
    let zero_reduce_sum = zero_counts.zeroed_like()?;

    // Build temporary directory for cross-tile window partials
    let mut temp_dir_guard = TempDirGuard::new(&opt.ioc.output_dir, "gc_bias_cross")
        .context("create per-run temp dir")?;
    let temp_dir = temp_dir_guard.path().to_path_buf();

    let mut reduce_state: ReduceState = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<ReduceState> {
            // Crossing files share one temp directory across chromosomes, so their filenames need
            // the run-global tile index, not `Tile.index`, which is chromosome-local.
            let crossing_file_tile_idx =
                u32::try_from(tile_idx).context("global tile index exceeded u32 range")?;
            let chr = tile.chr.as_str();
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let windows_chr: Option<&[IndexedInterval<u64>]> = windows_map
                .as_ref()
                .and_then(|m| m.get(chr).map(|v| v.as_slice()));
            let blacklist_chr: &[Interval<u64>] =
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]);
            let crossing_window_index_offset_chr = match &window_opt {
                // Fixed-size window ids start at zero on each chromosome. BED windows already
                // carry global original ids, and global mode uses the single shared window id 0.
                WindowSpec::Size(_) => *chromosome_window_offsets
                    .get(chr)
                    .with_context(|| format!("missing fixed-size window offset for {chr}"))?,
                _ => 0,
            };

            let (tile_counts, counter) = process_tile(
                tile,
                tile_span.as_ref(),
                windows_aligned_to_tiles,
                opt,
                &reference_metadata,
                windows_chr,
                &window_opt,
                avg_window_span,
                &zero_counts,
                &temp_dir,
                blacklist_chr,
                crossing_file_tile_idx,
                crossing_window_index_offset_chr,
            )?;

            let mut state = ReduceState::from_scaled_sum(zero_reduce_sum.clone());
            state.merge_scaled(&tile_counts.counts)?;
            state.add_weight(tile_counts.weight);
            state.merge_counters(counter);
            state.push_crossing_file(tile_counts.crossing_file);

            pb.inc(1);
            Ok(state)
        })
        .try_fold(
            || ReduceState::from_scaled_sum(zero_reduce_sum.clone()),
            |mut acc, state| -> Result<ReduceState> {
                acc = acc.merge(state?)?;
                Ok(acc)
            },
        )
        .try_reduce(
            || ReduceState::from_scaled_sum(zero_reduce_sum.clone()),
            |a, b| a.merge(b),
        )?;

    if options.show_progress {
        pb.finish_with_message("| Finished counting");
    } else {
        pb.finish_and_clear();
    }

    // Release tile-level inputs before global aggregation
    drop(tile_window_spans_for_threads);
    drop(tile_window_spans);
    drop(tiles);
    drop(windows_map);
    drop(blacklist_map);

    status_info!(options, target: COMMAND_TARGET, "Processing counts");
    if !windows_aligned_to_tiles && !matches!(window_opt, WindowSpec::Global) {
        let (cross_sum, cross_weight) = stream_crossing_files(
            reduce_state.crossing_files.clone(),
            &zero_counts,
            opt,
            avg_window_span,
        )?;
        reduce_state.scaled_sum.merge_from(&cross_sum)?;
        reduce_state.scaled_weight += cross_weight;
    }

    if let Err(err) = temp_dir_guard.remove() {
        warn!(
            target: COMMAND_TARGET,
            "warning: failed to remove temp dir {}: {}",
            temp_dir.display(),
            err
        );
    }

    let counted_windows = if matches!(window_opt, WindowSpec::Global) {
        usize::from(reduce_state.scaled_weight > 0)
    } else {
        reduce_state.scaled_weight
    };

    let scaling_weight = if matches!(window_opt, WindowSpec::Global) && counted_windows > 0 {
        1usize
    } else {
        reduce_state.scaled_weight
    };

    let global_counter = reduce_state.counters;

    ensure!(
        scaling_weight > 0,
        "No usable GC bias windows produced counts. Check window settings and blacklist coverage."
    );

    let mut avg_gc_counts = reduce_state.scaled_sum.clone();
    avg_gc_counts.scale_counts(1.0 / scaling_weight as f64)?;
    drop(reduce_state);

    // Smoothe GC counts on counts-level for each fragment length
    if !reference_metadata.skip_smoothing {
        status_info!(options, target: COMMAND_TARGET, "Smoothing cfDNA GC counts");
        avg_gc_counts.smooth_length_rows_in_place(
            reference_metadata.smoothing_sigma,
            reference_metadata.smoothing_radius,
        )?;
    }

    // Convert GC counts to grid with lengths x GC percentage bins
    status_info!(
        options,
        target: COMMAND_TARGET,
        "Converting cfDNA GC counts to GC percentage bins"
    );
    let mut avg_gc_pct_counts = avg_gc_counts.to_gc_percent_grid(0, 100)?;
    drop(avg_gc_counts);

    // Debias GC% bins for uneven rounding widths before any smoothing/interpolation
    apply_gc_percent_width_correction(&mut avg_gc_pct_counts, &reference_gc_percent_widths)?;

    // Sanity check that the unobservable cells are still zero-valued
    // here. If the logic is correct, this should always be true.
    let stats_by_support_status =
        stats_by_support_mask(&avg_gc_pct_counts, &reference_unobservables_support_mask);
    ensure!(
        stats_by_support_status.sum_for_unsupported == 0.0,
        "Unsupported bins in the count matrix had non-zero coverage. Report please."
    );

    intermediate_saver.save_file(
        &mut final_outputs,
        &avg_gc_pct_counts,
        "avg_cfdna_counts",
        "average cfDNA counts",
    )?;

    status_info!(options, target: COMMAND_TARGET, "Normalizing average counts");
    // Normalize GC counts array by its mean (just to remove the weighting scaling)
    // Ignores unsupported elements when calculating the mean
    let mut norm_gc_counts =
        mean_scale_array(&avg_gc_pct_counts, Some(&reference_outliers_support_mask))
            .ok_or_else(|| anyhow!("masked mean scaling had no supported elements"))?;

    intermediate_saver.save_file(
        &mut final_outputs,
        &norm_gc_counts,
        "normalized_avg_cfdna_counts",
        "normalized average cfDNA counts",
    )?;

    // Interpolate counts for the unsupported cells
    if !reference_metadata.skip_interpolation {
        status_info!(
            options,
            target: COMMAND_TARGET,
            "Interpolating counts for unsupported cells (very low reference counts)"
        );
        for (row_idx, mut length_row) in norm_gc_counts.outer_iter_mut().enumerate() {
            // Rows are contiguous so we can safely borrow a mutable slice for interpolation
            let row_slice = length_row
                .as_slice_mut()
                .context("GC histogram rows should be contiguous")?;
            let mut mask_row = reference_outliers_support_mask.row_mut(row_idx);
            let mask_slice = mask_row
                .as_slice_mut()
                .context("support mask rows should be contiguous")?;
            fill_unsupported_bins_with_polynomial(
                row_slice, mask_slice, 2, 3, 3,
                // Update mask when cells become supported
                // (NOTE: not currently used downstream but worth having for future checks)
                true,
            )?;
        }

        intermediate_saver.save_file(
            &mut final_outputs,
            &norm_gc_counts,
            "interpolated_cfdna_counts",
            "interpolated cfDNA counts",
        )?;
    }

    // Smoothe the *normalized* counts
    let do_smoothing = false;
    let smoothed_gc_counts = if do_smoothing {
        status_info!(
            options,
            target: COMMAND_TARGET,
            "Smoothing counts with 2D Gaussian kernel"
        );
        // 5-element kernel (-2...+2)
        let radius: usize = 2;
        // Standard deviation (quite sharp so not too smoothed)
        let sigma = 0.5;
        smoothe_counts_gaussian(&norm_gc_counts, sigma, radius)
    } else {
        norm_gc_counts
    };

    if do_smoothing {
        intermediate_saver.save_file(
            &mut final_outputs,
            &smoothed_gc_counts,
            "smoothed_cfdna_counts",
            "smoothed cfDNA counts",
        )?;
    }

    // Get greedy bins for lengths and GC
    // Maps "length -> length bin" and "gc -> gc bin"
    let length_bins = bin_greedily_by_mass(
        &smoothed_gc_counts,
        0,
        opt.min_length_bin_mass as f64,
        opt.min_length_bin_width,
    )?;
    let gc_bins = bin_greedily_by_mass(&smoothed_gc_counts, 1, opt.min_gc_bin_mass as f64, 1)?;

    // Collapse GC and length bins
    let (binned_ref_counts, binned_gc_counts) = {
        let ref_gc_binned = collapse_counts_by_bins(
            &avg_norm_ref_counts,
            1,
            &gc_bins,
            CollapseAggregation::Sum,
            None,
        )?;

        let cfdna_gc_binned = collapse_counts_by_bins(
            &smoothed_gc_counts,
            1,
            &gc_bins,
            CollapseAggregation::Sum,
            None,
        )?;

        let ref_length_binned = collapse_counts_by_bins(
            &ref_gc_binned,
            0,
            &length_bins,
            CollapseAggregation::Mean,
            None,
        )?;

        let cfdna_length_binned = collapse_counts_by_bins(
            &cfdna_gc_binned,
            0,
            &length_bins,
            CollapseAggregation::Mean,
            None,
        )?;

        (ref_length_binned, cfdna_length_binned)
    };

    intermediate_saver.save_file(
        &mut final_outputs,
        &binned_ref_counts,
        "binned_ref_counts",
        "binned reference counts",
    )?;

    intermediate_saver.save_file(
        &mut final_outputs,
        &binned_gc_counts,
        "binned_cfdna_counts",
        "binned cfDNA counts",
    )?;

    // Mask extreme GC bins and the shortest lengths to avoid unstable corrections
    let correction_support_mask = build_extreme_bins_support_mask(
        binned_gc_counts.dim(),
        opt.num_extreme_gc_bins as usize,
        opt.num_short_length_bins as usize,
    );
    let mut norm_gc_counts =
        mean_scale_per_length_array(&binned_gc_counts, 0., Some(&correction_support_mask));
    let mut norm_ref_counts =
        mean_scale_per_length_array(&binned_ref_counts, 0., Some(&correction_support_mask));

    // Set extreme GC bins to 1.0 in both arrays to avoid zero-division etc.
    status_info!(
        options,
        target: COMMAND_TARGET,
        "Setting counts to 1.0 for up to 2x{} extreme GC bins and {} shortest-length bins",
        opt.num_extreme_gc_bins,
        opt.num_short_length_bins
    );
    if correction_support_mask.iter().any(|&supported| !supported) {
        set_masked_entries_to_value(&mut norm_gc_counts, &correction_support_mask, 1.0);
        set_masked_entries_to_value(&mut norm_ref_counts, &correction_support_mask, 1.0);
    }

    intermediate_saver.save_file(
        &mut final_outputs,
        &norm_gc_counts,
        "normalized_binned_cfdna_counts",
        "normalized binned cfDNA counts",
    )?;
    intermediate_saver.save_file(
        &mut final_outputs,
        &norm_ref_counts,
        "normalized_binned_ref_counts",
        "normalized binned reference counts",
    )?;

    // Resolve outlier handling configuration
    let (outlier_rule, outlier_action, outlier_scope) = opt.outlier_settings()?;
    let mut outlier_stats = OutlierStats::default();

    // TODO: Update this pipeline list
    // Calculate correction matrix
    // 1) Divide cfDNA counts by reference counts
    // 2) Normalize each fragment length to mean=1.0
    // 3) Interpolate masked bins so extremes follow neighbouring corrections
    // 4) Clamp any leftover extreme corrections
    let correction_matrix = {
        let raw_correction_matrix = &norm_gc_counts / &norm_ref_counts;

        // Normalize correction matrix per fragment length to be centered around 1.0
        // Still ignores extreme GC bins in the mean-calculations
        let mut norm_correction_matrix =
            mean_scale_per_length_array(&raw_correction_matrix, 0., Some(&correction_support_mask));

        if correction_support_mask.iter().any(|&supported| !supported) {
            status_info!(
                options,
                target: COMMAND_TARGET,
                "Interpolating corrections for up to 2x{} extreme GC bins and {} shortest-length bins",
                opt.num_extreme_gc_bins,
                opt.num_short_length_bins
            );
            interpolate_masked_corrections(&mut norm_correction_matrix, &correction_support_mask)?;
        }

        if !matches!(outlier_rule, OutlierRule::None) {
            status_info!(
                options,
                target: COMMAND_TARGET,
                "Applying outlier handling to correction matrix"
            );
            let support = Some(&correction_support_mask);
            outlier_stats = apply_outliers_to_matrix(
                &mut norm_correction_matrix,
                support,
                outlier_scope,
                outlier_rule,
                outlier_action,
            );
        }

        // Sanity clamp of corrections; tracked separately from outlier stats
        let mut hard_clamp_count = 0usize;
        norm_correction_matrix.mapv_inplace(|v| {
            if v < CORRECTION_CLAMP_RANGE.0 {
                hard_clamp_count += 1;
                CORRECTION_CLAMP_RANGE.0
            } else if v > CORRECTION_CLAMP_RANGE.1 {
                hard_clamp_count += 1;
                CORRECTION_CLAMP_RANGE.1
            } else {
                v
            }
        });
        outlier_stats.hard_clamped = hard_clamp_count;

        // Re-normalize correction matrix per fragment length to be centered around 1.0 (no mask).
        norm_correction_matrix =
            mean_scale_per_length_array(&norm_correction_matrix, 0., None::<&Array2<bool>>);

        // Make correction factors multipliers by inverting elements to 1 / x
        // Zeros remain 0s
        invert_elementwise_with_zeros_inplace(&mut norm_correction_matrix);

        Ok::<Array2<f64>, anyhow::Error>(norm_correction_matrix)
    }?;

    // Length-bin frequencies (normalized) used for length-agnostic GC correction
    let length_bin_frequencies: Array1<f64> = {
        let per_length_totals = binned_gc_counts.sum_axis(Axis(1));
        let total: f64 = per_length_totals.iter().sum();
        ensure!(total > 0.0, "Total fragment count for length bins is zero");
        per_length_totals.iter().map(|v| *v / total).collect()
    };

    // Write every final output to the temp directory before moving any of them into place
    // This keeps failed writes from leaving a mix of old and new final files

    // Save reusable correction package with metadata for downstream commands
    let correction_pkg = GCCorrectionPackage::from_components(
        GC_CORRECTION_SCHEMA_VERSION,
        &length_bins,
        &gc_bins,
        correction_matrix.clone(),
        length_bin_frequencies.clone(),
        &reference_metadata,
    )?;
    let correction_package_path = opt
        .ioc
        .output_dir
        .join(dot_join(&[prefix, "gc_bias_correction.zarr"]));
    let temp_correction_package_path = final_outputs.temp_path_for(&correction_package_path)?;
    correction_pkg.write_zarr(&temp_correction_package_path)?;
    final_outputs.record(
        temp_correction_package_path,
        correction_package_path.clone(),
    )?;

    // Plot the avg. gc-bias across lengths for quick QC
    #[cfg(feature = "plotters")]
    {
        use crate::commands::gc_bias::plotting::plot_gc_bias;

        status_info!(options, target: COMMAND_TARGET, "Plotting GC bias");

        let temp_plot_paths = plot_gc_bias(
            final_outputs.temp_dir(),
            prefix,
            &gc_bins,
            &length_bins,
            &correction_matrix,
            &length_bin_frequencies,
            &reference_metadata,
            &avg_gc_pct_counts,
            &binned_gc_counts,
        )?;
        // These plots were written to final_outputs.temp_dir() with their final filenames
        // Record each one as output_dir/<file name>, then move all outputs at the end
        final_outputs.record_temp_files_with_same_names_in(temp_plot_paths, &opt.ioc.output_dir)?;
    }

    final_outputs.move_into_place()?;

    let elapsed = start_time.elapsed();
    let mut extra_lines = vec![format!(
        "Windows processed: {} total, {} with counts",
        total_windows, counted_windows
    )];
    if !matches!(outlier_rule, OutlierRule::None) {
        extra_lines.push("Outlier handling:".to_string());
        extra_lines.push("  > Limits estimated from reference-supported bins only".to_string());
        extra_lines.push(format!(
            "  Supported cells examined: {} (winsorized: {})",
            outlier_stats.total_examined, outlier_stats.total_outliers_handled
        ));
        extra_lines.push(
            "  > 'supported' = bins the reference marks valid (used to set limits)".to_string(),
        );
        extra_lines.push(format!(
            "  Unsupported cells examined: {} (winsorized: {})",
            outlier_stats.unsupported_examined, outlier_stats.unsupported_outliers_handled
        ));
        extra_lines.push(
            "  > 'unsupported' = bins the reference masks out (winsorized after interpolation)"
                .to_string(),
        );
    }
    extra_lines.push(format!(
        "Extreme GC-bias values clamped to [{:.1},{:.1}] before final scaling: {}",
        CORRECTION_CLAMP_RANGE.0, CORRECTION_CLAMP_RANGE.1, outlier_stats.hard_clamped
    ));
    if options.report_statistics {
        print_fragment_run_statistics(
            &global_counter.base,
            elapsed,
            FragmentRunStatisticsOptions {
                include_section_header: true,
                notes: &[],
                labels: DEFAULT_FRAGMENT_STATISTICS_LABELS,
                blacklist_excluded_fragments: None,
                gc: None,
            },
            extra_lines.iter().map(String::as_str),
        );
    }
    Ok(GCBiasRunResult {
        counters: global_counter,
        correction_package_path: correction_package_path.clone(),
        output_files: vec![correction_package_path],
    })
}

fn validate_reference_chromosomes(
    reference_chromosomes: &[String],
    run_chromosomes: &[String],
) -> Result<()> {
    let reference_set: FxHashSet<&str> = reference_chromosomes.iter().map(String::as_str).collect();
    let run_set: FxHashSet<&str> = run_chromosomes.iter().map(String::as_str).collect();
    ensure!(
        reference_set == run_set,
        "Reference GC package was built for chromosomes [{}], but this run selected [{}]",
        reference_chromosomes.join(","),
        run_chromosomes.join(",")
    );
    Ok(())
}

fn validate_fixed_size_window_for_reference(
    window_opt: &WindowSpec,
    reference_metadata: &ReferenceGCMetadata,
) -> Result<()> {
    if let WindowSpec::Size(window_bp) = window_opt {
        ensure!(
            *window_bp >= reference_metadata.max_fragment_length as u64,
            "--by-size ({}) must be >= max fragment length from --ref-gc-file ({}) for gc-bias fixed-size windowing",
            window_bp,
            reference_metadata.max_fragment_length
        );
    }
    Ok(())
}

fn validate_reference_contig_footprint(
    reference_metadata: &ReferenceGCMetadata,
    run_reference_contig_footprint: &[ContigFootprintEntry],
) -> Result<()> {
    ensure!(
        reference_metadata.reference_contig_footprint == run_reference_contig_footprint,
        "Reference GC package was built against a different reference contig footprint than --ref-2bit"
    );
    Ok(())
}

fn process_tile(
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    windows_aligned_to_tiles: bool,
    opt: &GCConfig,
    reference_metadata: &ReferenceGCMetadata,
    windows_opt: Option<&[IndexedInterval<u64>]>,
    window_opt: &WindowSpec,
    avg_window_span: f64,
    template: &GCCounts,
    temp_dir: &PathBuf,
    blacklist_intervals: &[Interval<u64>],
    crossing_file_tile_idx: u32,
    crossing_window_index_offset: u64,
) -> anyhow::Result<(WindowState, GCCounters)> {
    let apply_window_scaling = !matches!(window_opt, WindowSpec::Global);

    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = GCCounters::default();

    // Reuse `WindowState` as the per-tile return value. On this return path, only `counts`,
    // `weight`, and `crossing_file` are read later. The interval-like fields are not used, so
    // keep a trivial checked placeholder here
    let dummy_interval = Interval::new(0, 1)?;
    let mut tile_output = WindowState::new(0, dummy_interval, true, template)?;
    let mut crossing_parts: Vec<CrossingPart> = Vec::new();

    // Prep tile windows:
    // - Fixed-size windows get two rolling WindowState buffers so we can slide without reallocating.
    // - BED/global windows get a Vec of WindowState built from the cached span to stay aligned with global indices.
    // - If the cached span is empty (no BED overlaps), we skip the tile before touching reference sequence or BAM.
    let prepared_windows = prepare_tile_windows(
        window_opt,
        windows_opt,
        tile,
        tile_window_span,
        chrom_len,
        template,
    )?;
    if prepared_windows.skip_tile {
        return Ok((tile_output, counter));
    }
    let mut windows = prepared_windows.windows;
    let streaming_buffers = prepared_windows.streaming_buffers;
    let tile_core_interval = tile.core.try_to_u64()?;

    let seq_start = tile.fetch_start() as u64;
    let seq_end = tile.fetch_end().min(chrom_len as u32) as u64;
    let sequence_interval = Interval::new(seq_start, seq_end)?;

    let gc_prefixes = {
        // Load only the tile span (core plus halo)
        let mut seq_bytes = read_seq_in_range(
            &opt.ref_genome.ref_2bit,
            &tile.chr,
            seq_start as usize..seq_end as usize,
        )?;
        // Blacklist GC prefixes to avoid using blacklist-overlapping fragments in bias-estimation
        // NOTE: Downstream commands don't blacklist prefixes so it's possible to correct such fragments
        // using their full GC context
        apply_blacklist_mask_to_seq(&mut seq_bytes, blacklist_intervals, seq_start);
        build_gc_prefixes(&seq_bytes)
    };

    // Decide whether we run streaming (fixed-size windows) or per-window Vec (BED/global)
    let using_streaming = streaming_buffers.is_some();

    let mut tile_window_intervals: Option<Vec<IndexedInterval<u64>>> = None;
    if !using_streaming {
        if windows.is_empty() {
            return Ok((tile_output, counter));
        }
        tile_window_intervals = Some(
            windows
                .iter()
                .enumerate()
                .map(|(local_idx, window)| {
                    IndexedInterval::new(window.start(), window.end(), local_idx as u64)
                })
                .collect::<crate::Result<Vec<_>>>()?,
        );
    }

    // Use the same assign-by threshold for streaming fixed-size windows and BED/global lookup
    let min_overlap_fraction = min_overlap_fraction_for_window_assignment(
        opt.window_assignment.assign_by,
        reference_metadata.max_fragment_length as u64,
    );

    reader
        .fetch((tid, tile.fetch_start() as i64, tile.fetch_end() as i64))
        .context(format!("fetch {}", &tile.chr))?;

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let min_len = reference_metadata.min_fragment_length as u32;
        let max_len = reference_metadata.max_fragment_length as u32;
        move |f: &Fragment| {
            let len = f.len();
            len >= min_len && len <= max_len
        }
    };

    let mut iter = if opt.unpaired.reads_are_fragments {
        let min_mapq = opt.min_mapq;
        let include_read_fn = move |r: &Record| default_include_read_unpaired(r, min_mapq);
        fragments_from_bam(
            reader.records().map(|r| r.map_err(anyhow::Error::from)),
            include_read_fn,
            None,
            fragment_filter,
            true,
        )
        .with_local_counters()
    } else {
        let min_mapq = opt.min_mapq;
        let require_proper_pair = opt.require_proper_pair;
        let include_read_fn =
            move |r: &Record| default_include_read_paired_end(r, require_proper_pair, min_mapq);
        fragments_from_bam(
            reader.records().map(|r| r.map_err(anyhow::Error::from)),
            include_read_fn,
            None,
            fragment_filter,
            false,
        )
        .with_local_counters()
    };

    if let Some((window_bp, mut current, mut next)) = streaming_buffers {
        // Fixed-size windows
        // Streaming branch: maintain a sliding pair of window buffers to avoid per-window allocations
        for fragment_res in iter.by_ref() {
            let fragment = fragment_res.context("reading fragment")?;
            let fragment_length = fragment.len();
            let fragment_interval = fragment.interval.try_to_u64()?;

            // Only count fragments whose start lies inside the tile core to avoid double counting
            if fragment.start() < tile.core_start() || fragment.start() >= tile.core_end() {
                continue;
            }

            // If the fragment is past the current window, finalize and advance the buffers
            while fragment.start() as u64 >= current.end() {
                finalize_window_buffer(
                    &mut current,
                    &gc_prefixes,
                    sequence_interval,
                    tile_core_interval,
                    windows_aligned_to_tiles,
                    apply_window_scaling,
                    opt,
                    avg_window_span,
                    &mut tile_output,
                    &mut crossing_parts,
                    crossing_window_index_offset,
                )?;
                let current_next = match next.take() {
                    Some(window) => window,
                    None => break,
                };
                // Fixed-size mode only needs two live window buffers at a time.
                // The helper owns the "promote next -> current, then build a new next"
                // transition so the chromosome-end behavior stays in one place.
                (current, next) = advance_fixed_size_streaming_buffers(
                    current,
                    current_next,
                    window_bp,
                    chrom_len,
                    tile_core_interval,
                    template,
                )?;
            }

            let Some(gc_count) = get_fragment_gc(
                fragment_interval,
                sequence_interval,
                reference_metadata.end_offset as u32,
                &gc_prefixes,
                1.0,
            )?
            else {
                continue;
            };

            let assignment_interval = fragment_assignment_interval(
                &tile.chr,
                &fragment,
                opt.window_assignment.assign_by,
            )?;
            let fragment_span_length = assignment_interval.len() as f64;

            // Choose buffers to update: at most current and next for fixed-size windows
            let mut current_weight: Option<f64> = None;
            let overlap_current = overlap_length(assignment_interval, current.interval);
            if overlap_current > 0 {
                let fraction = overlap_current as f64 / fragment_span_length;
                if fraction >= min_overlap_fraction {
                    // Increment counter by 1.0 or by the overlap fraction
                    let count_weight = match opt.window_assignment.assign_by {
                        WindowAssigner::CountOverlap => fraction,
                        _ => 1.0f64,
                    };
                    current_weight = Some(count_weight);
                }
            }

            let mut next_weight: Option<f64> = None;
            if let Some(next_window) = next.as_ref() {
                let overlap_next = overlap_length(assignment_interval, next_window.interval);
                if overlap_next > 0 {
                    let fraction = overlap_next as f64 / fragment_span_length;
                    if fraction >= min_overlap_fraction {
                        // Increment counter by 1.0 or by the overlap fraction
                        let count_weight = match opt.window_assignment.assign_by {
                            WindowAssigner::CountOverlap => fraction,
                            _ => 1.0,
                        };
                        next_weight = Some(count_weight);
                    }
                }
            }

            if current_weight.is_none() && next_weight.is_none() {
                continue;
            }

            counter.base.counted_fragments += 1;

            if let Some(count_weight) = current_weight {
                current.has_counts = true;
                current
                    .counts
                    .incr_weighted(fragment_length as usize, gc_count, count_weight);
            }
            match (next_weight, next.as_mut()) {
                (Some(count_weight), Some(next_window)) => {
                    next_window.has_counts = true;
                    next_window.counts.incr_weighted(
                        fragment_length as usize,
                        gc_count,
                        count_weight,
                    );
                }
                (Some(_), None) => {
                    bail!("computed next-window GC counts without a live next window");
                }
                (None, _) => {}
            }
        }

        // Flush the current buffer and the optional next buffer so partially
        // filled windows still reach downstream reducers.
        finalize_window_buffer(
            &mut current,
            &gc_prefixes,
            sequence_interval,
            tile_core_interval,
            windows_aligned_to_tiles,
            apply_window_scaling,
            opt,
            avg_window_span,
            &mut tile_output,
            &mut crossing_parts,
            crossing_window_index_offset,
        )?;
        if let Some(next_window) = next.as_mut() {
            finalize_window_buffer(
                next_window,
                &gc_prefixes,
                sequence_interval,
                tile_core_interval,
                windows_aligned_to_tiles,
                apply_window_scaling,
                opt,
                avg_window_span,
                &mut tile_output,
                &mut crossing_parts,
                crossing_window_index_offset,
            )?;
        }
    } else {
        let mut window_ptr = 0usize; // Genomic window pointer reused by overlap finder
        {
            let tile_window_intervals = tile_window_intervals.as_ref().context(
                "non-streaming GC-bias counting requires prepared tile window intervals",
            )?;

            // Iterate fragments and count GC
            for fragment_res in iter.by_ref() {
                let fragment = fragment_res.context("reading fragment")?;
                let fragment_length = fragment.len();
                let fragment_interval = fragment.interval.try_to_u64()?;

                // Only count fragments whose start lies inside the tile core to avoid double counting
                if fragment.start() < tile.core_start() || fragment.start() >= tile.core_end() {
                    continue;
                }

                let Some(gc_count) = get_fragment_gc(
                    fragment_interval,
                    sequence_interval,
                    reference_metadata.end_offset as u32,
                    &gc_prefixes,
                    1.0,
                )?
                else {
                    continue;
                };

                // Find candidate windows from the interval implied by window assignment.
                // The helper collapses midpoint mode to a 1 bp midpoint and otherwise returns
                // the fragment interval.
                let window_selection_interval = fragment_assignment_interval(
                    &tile.chr,
                    &fragment,
                    opt.window_assignment.assign_by,
                )?;
                let overlapping_windows = find_overlapping_windows(
                    chrom_len,
                    &mut window_ptr,
                    Some(tile_window_intervals.as_slice()),
                    None,
                    window_selection_interval,
                    min_overlap_fraction,
                    reference_metadata.max_fragment_length as u64,
                )?;
                let overlapping_windows = if let Some(overlaps) = overlapping_windows {
                    overlaps
                } else {
                    continue;
                };

                counter.base.counted_fragments += 1;

                // Increment counter by 1.0 or by the overlap fraction
                for overlapped_window in overlapping_windows.windows {
                    let count_weight = match opt.window_assignment.assign_by {
                        WindowAssigner::CountOverlap => overlapped_window.overlap_fraction,
                        _ => 1.0f64,
                    };
                    if let Some(state) = windows.get_mut(overlapped_window.idx) {
                        state.has_counts = true;
                        state.counts.incr_weighted(
                            fragment_length as usize,
                            gc_count,
                            count_weight,
                        );
                    }
                }
            }
        }

        // Release interval cache before finalizing windows
        drop(tile_window_intervals);

        for mut w in windows {
            finalize_window_buffer(
                &mut w,
                &gc_prefixes,
                sequence_interval,
                tile_core_interval,
                windows_aligned_to_tiles,
                apply_window_scaling,
                opt,
                avg_window_span,
                &mut tile_output,
                &mut crossing_parts,
                crossing_window_index_offset,
            )?;
        }
    }

    // Release prefix arrays before writing crossing parts
    drop(gc_prefixes);

    counter.add_from_snapshot(iter.counters_snapshot());

    // Write crossing parts when windows are not aligned
    // Streaming and non-streaming share the same writer
    tile_output.crossing_file = if using_streaming {
        if !windows_aligned_to_tiles {
            write_crossing_parts(temp_dir, crossing_file_tile_idx, template, &crossing_parts)?
        } else {
            None
        }
    } else if !windows_aligned_to_tiles && !matches!(window_opt, WindowSpec::Global) {
        write_crossing_parts(temp_dir, crossing_file_tile_idx, template, &crossing_parts)?
    } else {
        None
    };

    Ok((tile_output, counter))
}

pub(crate) fn process_window_in_place(
    gc_counts: &mut GCCounts,
    opt: &GCConfig,
    avg_window_span: f64,
) -> Result<bool> {
    // Check window has enough valid positions
    if gc_counts.pct_acgt() < opt.min_window_acgt_pct as f64 {
        return Ok(false);
    }
    if gc_counts.sum() == 0.0 {
        return Ok(false);
    }

    // Get mean coverage to scale by (supported cells only)
    let mean_count: f64 = gc_counts.mean();

    // Weighting one or both count distributions

    let num_acgt = gc_counts.num_acgt_out_of.0;
    if num_acgt == 0 {
        // No positions observed
        return Ok(false);
    }
    // Use avg window size to "scale the scaling"
    // So the counts don't explode in size
    // While keeping the relative scale between windows consistent
    // Mean scale, then multiply by the usable window-size (scaled to lower values)
    gc_counts.scale_counts((1. / mean_count) * (num_acgt as f64 / avg_window_span))?;

    Ok(true)
}

pub(crate) fn process_window(
    mut gc_counts: GCCounts,
    opt: &GCConfig,
    avg_window_span: f64,
) -> Result<Option<GCCounts>> {
    if process_window_in_place(&mut gc_counts, opt, avg_window_span)? {
        return Ok(Some(gc_counts));
    }
    Ok(None)
}

pub(crate) fn mean_scale_per_length_array<S, M>(
    x: &ArrayBase<S, Ix2>,
    pseudo_count: f64,
    support_mask: Option<&ArrayBase<M, Ix2>>,
) -> Array2<f64>
where
    S: Data<Elem = f64>,
    M: Data<Elem = bool>,
{
    let (n_rows, n_cols) = x.dim();
    if let Some(m) = support_mask {
        assert_eq!(
            m.dim(),
            (n_rows, n_cols),
            "Mask shape {:?} must match counts shape {:?}",
            m.dim(),
            (n_rows, n_cols)
        );
    }

    let mut out = Array2::zeros((n_rows, n_cols));

    for row_idx in 0..n_rows {
        let (row_sum, valid_count) = if let Some(mask_arr) = support_mask {
            let mask_row = mask_arr.row(row_idx);
            let mut sum = 0.0;
            let mut count = 0usize;
            for (value, &is_valid) in x.row(row_idx).iter().zip(mask_row.iter()) {
                if is_valid {
                    sum += *value;
                    count += 1;
                }
            }
            (sum, count)
        } else {
            (x.row(row_idx).sum(), n_cols)
        };

        // Keep empty rows at zero instead of producing NaNs when the mean is zero.
        if valid_count == 0 || row_sum == 0.0 {
            continue;
        }

        let denom = if valid_count > 0 {
            (row_sum / valid_count as f64) + pseudo_count * valid_count as f64
        } else {
            1.0
        };

        for (col_idx, value) in x.row(row_idx).iter().enumerate() {
            let numerator = *value + pseudo_count;
            out[(row_idx, col_idx)] = numerator / denom;
        }
    }

    out
}

pub(crate) fn interpolate_masked_corrections(
    matrix: &mut Array2<f64>,
    support_mask: &Array2<bool>,
) -> Result<()> {
    let (n_rows, n_cols) = matrix.dim();
    if n_rows == 0 || n_cols == 0 {
        return Ok(());
    }
    if !support_mask.iter().any(|&supported| !supported) {
        return Ok(());
    }

    // Work on a mutable copy so we can treat newly interpolated bins as anchors
    let mut mask = support_mask.to_owned();

    for row_idx in 0..n_rows {
        let mut row = matrix.row_mut(row_idx);
        let row_slice = row
            .as_slice_mut()
            .context("GC histogram rows should be contiguous")?;
        let mut mask_row = mask.row_mut(row_idx);
        let mask_slice = mask_row
            .as_slice_mut()
            .context("support mask rows should be contiguous")?;
        fill_unsupported_bins_with_polynomial(row_slice, mask_slice, 2, 3, 3, true)?;
    }

    for col_idx in 0..n_cols {
        let mut column_values: Vec<f64> = (0..n_rows)
            .map(|row_idx| matrix[(row_idx, col_idx)])
            .collect();
        let mut column_mask: Vec<bool> = (0..n_rows)
            .map(|row_idx| mask[(row_idx, col_idx)])
            .collect();

        fill_unsupported_bins_with_polynomial(&mut column_values, &mut column_mask, 2, 3, 3, true)?;

        for row_idx in 0..n_rows {
            matrix[(row_idx, col_idx)] = column_values[row_idx];
            mask[(row_idx, col_idx)] = column_mask[row_idx];
        }
    }

    Ok(())
}

// Overall scaling
// Elements that are marked as `false` in the support mask are
// still scaled but do not contribute to the mean
pub(crate) fn mean_scale_array<S, M>(
    x: &ArrayBase<S, Ix2>,
    support_mask: Option<&ArrayBase<M, Ix2>>,
) -> Option<Array2<f64>>
where
    S: Data<Elem = f64>,
    M: Data<Elem = bool>,
{
    let (n_rows, n_cols) = x.dim();
    if let Some(m) = support_mask {
        assert_eq!(
            m.dim(),
            (n_rows, n_cols),
            "Mask shape {:?} must match counts shape {:?}",
            m.dim(),
            (n_rows, n_cols)
        );
    }

    let mut out = x.to_owned();

    let mean = if let Some(mask) = support_mask {
        let mut total_val = 0f64;
        let mut num_elements = 0u64;
        Zip::from(mask).and(x).for_each(|&use_element, &value| {
            if use_element {
                total_val += value;
                num_elements += 1;
            }
        });
        if num_elements == 0 {
            return None;
        }
        total_val / num_elements as f64
    } else {
        x.mean()?
    };

    if mean == 0.0 {
        return None;
    }
    out /= mean;
    Some(out)
}

/// Invert elements in an array (x) to `1 / x`, keeping 0s as 0.
fn invert_elementwise_with_zeros_inplace<S>(x: &mut ArrayBase<S, Ix2>)
where
    S: DataMut<Elem = f64>,
{
    x.mapv_inplace(|v| if v == 0.0 { 0.0 } else { 1.0 / v });
}

/// Writes optional GC-bias intermediate arrays through the staged final-output writer.
///
/// The saver owns the intermediate naming state and the temp/final directory roots. Each
/// `save_file` call writes to `temp_dir` and records the matching final path in `FinalOutputFiles`.
/// The caller moves all recorded outputs into place once the command has finished writing.
struct IntermediateFileSaver {
    save_intermediates: bool,
    temp_dir: PathBuf,
    final_dir: PathBuf,
    prefix: String,
    log_statuses: bool,
    previously_saved: usize,
}

impl IntermediateFileSaver {
    /// Create a saver for intermediate arrays.
    ///
    /// `temp_dir` is where files are written immediately. `final_dir` is where they will be moved
    /// by `FinalOutputFiles::move_into_place`.
    fn new(
        save_intermediates: bool,
        temp_dir: PathBuf,
        final_dir: PathBuf,
        prefix: String,
        log_statuses: bool,
    ) -> Self {
        IntermediateFileSaver {
            save_intermediates,
            temp_dir,
            final_dir,
            prefix,
            log_statuses,
            previously_saved: 0,
        }
    }

    /// Write one intermediate array to the temp directory and record where it should end up.
    ///
    /// `final_outputs` is only used to record the completed temp-to-final path pair. Path naming
    /// stays inside this saver so call sites do not need to know intermediate filenames.
    fn save_file<S>(
        &mut self,
        final_outputs: &mut FinalOutputFiles,
        x: &ArrayBase<S, Ix2>,
        file_tag: &str,
        msg_tag: &str,
    ) -> Result<()>
    where
        S: Data<Elem = f64>,
    {
        if self.save_intermediates {
            status_info!(self, target: COMMAND_TARGET, "Intermediate file: Saving {}", msg_tag);
            let file_name = format!("gc_bias.{}.{}.npy", file_tag, self.previously_saved);
            let output_name = dot_join(&[self.prefix.as_str(), file_name.as_str()]);
            let temp_path = self.temp_dir.join(output_name.as_str());
            let final_path = self.final_dir.join(output_name);
            write_npy(&temp_path, x).context(format!(
                "Failed to write intermediate file {}",
                self.previously_saved
            ))?;
            final_outputs.record(temp_path, final_path)?;
            self.previously_saved += 1;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    include!("gc_bias_tests.rs");
}
