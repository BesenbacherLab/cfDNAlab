use crate::{
    commands::{
        cli_common::*,
        counters::GCCounters,
        gc_bias::{
            CORRECTION_CLAMP_RANGE, GC_CORRECTION_SCHEMA_VERSION,
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
                WindowState, compute_window_acgt, compute_window_stats, fixed_size_window_bounds,
                overlap_length, prepare_tile_windows,
            },
        },
    },
    shared::{
        bam::create_chromosome_reader,
        bed::{Windows, load_windows_from_bed},
        blacklist::apply_blacklist_mask_to_seq,
        fragment::minimal_fragment::Fragment,
        fragment_iterator::fragments_from_bam,
        midpoint::midpoint_random_even_with_thread_rng,
        overlaps::find_overlapping_windows,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        reference::read_seq_in_range,
        thread_pool::init_global_pool,
        tiled_run::{
            Tile, TileWindowSpan, build_tiles, make_temp_dir, precompute_tile_window_spans,
        },
    },
};
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array1, Array2, ArrayBase, Axis, Data, DataMut, Ix2, Zip};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{fs::create_dir_all, path::PathBuf, sync::Arc, time::Instant};

/// Get the GC count for a fragment after trimming ends.
///
/// Returns `None` when the trimmed interval is empty, outside the loaded sequence,
/// or lacks sufficient ACGT support.
pub fn get_fragment_gc(
    fragment: &Fragment,
    seq_start: u64,
    seq_end: u64,
    end_offset: u32,
    gc_prefixes: &GCPrefixes,
    min_acgt_fraction: f32,
) -> Result<Option<usize>> {
    let gc_window_start = fragment.start.saturating_add(end_offset);
    let gc_window_end = fragment.end.saturating_sub(end_offset);
    if gc_window_end <= gc_window_start {
        return Ok(None);
    }
    if gc_window_end as u64 > seq_end || (gc_window_start as u64) < seq_start {
        return Ok(None);
    }

    let acgt = gc_prefixes.acgt[(gc_window_end - seq_start as u32) as usize]
        - gc_prefixes.acgt[(gc_window_start - seq_start as u32) as usize];
    if acgt < MIN_ACGT_BASES_FOR_GC_FRACTION
        || (acgt as f32 / (gc_window_end - gc_window_start) as f32) < min_acgt_fraction
    {
        return Ok(None);
    }

    let gc = gc_prefixes.gc[(gc_window_end - seq_start as u32) as usize]
        - gc_prefixes.gc[(gc_window_start - seq_start as u32) as usize];
    ensure!(
        gc <= (gc_window_end - gc_window_start),
        "GC count exceeded interval length: {} > {}",
        gc,
        gc_window_end - gc_window_start
    );
    Ok(Some(gc as usize))
}

pub fn finalize_window_buffer(
    buf: &mut WindowState,
    gc_prefixes: &GCPrefixes,
    seq_start: u64,
    seq_end: u64,
    tile_core_start: u64,
    tile_core_end: u64,
    windows_aligned_to_tiles: bool,
    apply_window_scaling: bool,
    opt: &GCConfig,
    avg_window_span: f64,
    out: &mut WindowState, // TODO: Rename to something meaningful
    crossing_parts: &mut Vec<CrossingPart>,
) -> Result<()> {
    if !buf.has_counts {
        buf.counts.clear();
        return Ok(());
    }

    compute_window_acgt(buf, gc_prefixes, seq_start, seq_end)?;

    if !apply_window_scaling {
        out.counts.merge_from(&buf.counts)?;
        out.weight += 1;
        buf.counts.clear();
        buf.has_counts = false;
        return Ok(());
    }

    let crosses_tile = buf.start < tile_core_start || buf.end > tile_core_end;
    if !windows_aligned_to_tiles && crosses_tile {
        crossing_parts.push(CrossingPart {
            idx: buf.idx as usize,
            counts: buf.counts.clone(),
        });
    } else if process_window_in_place(&mut buf.counts, opt, Some(avg_window_span))? {
        out.counts.merge_from(&buf.counts)?;
        out.weight += 1;
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
    fn new(template: &GCCounts) -> Self {
        Self {
            scaled_sum: template
                .zeroed_like()
                .expect("failed to create zeroed reduce state"),
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
        self.crossing_files.extend(other.crossing_files.drain(..));

        Ok(self)
    }
}

pub fn run(opt: &GCConfig) -> Result<()> {
    let start_time = Instant::now();
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, &opt.ioc.bam.as_path())?;
    let window_opt = opt.windows.resolve_windows();
    let mut intermediate_saver =
        IntermediateFileSaver::new(opt.save_intermediates, opt.ioc.output_dir.clone());

    println!("Start: Loading reference GC bias");
    let ReferenceGCData {
        counts: reference_counts,
        unobservables_support_mask: reference_unobservables_support_mask,
        outliers_support_mask: mut reference_outliers_support_mask,
        gc_percent_widths: reference_gc_percent_widths,
        metadata: reference_metadata,
    } = load_reference_gc_data(&opt.ref_gc_dir)?;
    let avg_norm_ref_counts = mean_scale_per_length_array(
        &reference_counts,
        0.,
        Some(&reference_outliers_support_mask),
    );
    drop(reference_counts);

    // Create output directory
    create_dir_all(&opt.ioc.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = load_blacklist_map(opt.blacklist.as_ref(), 1, 0, &chromosomes)?;

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            println!("Start: Loading window coordinates");
            let mut wds = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            println!("Start: Merging overlapping/touching windows");
            let mut merged: FxHashMap<String, Windows> =
                FxHashMap::with_capacity_and_hasher(wds.len(), Default::default());
            let mut next_idx = 0u64;
            for chr in &chromosomes {
                if let Some(ws) = wds.remove(chr) {
                    // Flatten in-place
                    let (flat, next) = ws.into_flattened_reindexed(next_idx);
                    next_idx = next;
                    merged.insert(chr.clone(), flat);
                }
            }
            Some(merged)
        }
        _ => None,
    };

    let window_stats =
        compute_window_stats(&window_opt, windows_map.as_ref(), &contigs, &chromosomes)?;
    let avg_window_span = window_stats.avg_span;
    let total_windows = window_stats.total_windows;

    // Build tiles (core plus halo = max fragment length) to bound memory per worker
    let halo_bp = reference_metadata.max_fragment_length as u32;
    let align_bp = match &window_opt {
        WindowSpec::Size(bp) => Some(*bp),
        _ => None,
    };
    let (tiles, windows_aligned_to_tiles) =
        build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, align_bp)?;
    let pb = Arc::new(ProgressBar::new(tiles.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(&tiles, |chr| {
        windows_lookup
            .and_then(|m| m.get(chr).map(|w| w.as_slice()))
            .unwrap_or(&[])
    }));
    let tile_window_spans_for_threads = tile_window_spans.clone();

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads as usize)?;

    println!("Start: Counting per tile");

    pb.set_position(0);

    let zero_counts = GCCounts::new(
        reference_metadata.min_fragment_length as usize,
        reference_metadata.max_fragment_length as usize,
        reference_metadata.end_offset as usize,
        (0, 0),
    )?;

    // Build temporary directory for cross-tile window partials
    let temp_dir =
        make_temp_dir(&opt.ioc.output_dir, "gc_bias_cross").context("create per-run temp dir")?;

    let mut reduce_state: ReduceState = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<ReduceState> {
            let chr = tile.chr.as_str();
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let windows_chr: Option<&[(u64, u64, u64)]> = windows_map
                .as_ref()
                .and_then(|m| m.get(chr).map(|v| v.as_slice()));
            let blacklist_chr: &[(u64, u64)] =
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]);

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
            )?;

            let mut state = ReduceState::new(&zero_counts);
            state.merge_scaled(&tile_counts.counts)?;
            state.add_weight(tile_counts.weight);
            state.merge_counters(counter);
            state.push_crossing_file(tile_counts.crossing_file);

            pb.inc(1);
            Ok(state)
        })
        .try_fold(
            || ReduceState::new(&zero_counts),
            |mut acc, state| -> Result<ReduceState> {
                acc = acc.merge(state?)?;
                Ok(acc)
            },
        )
        .try_reduce(|| ReduceState::new(&zero_counts), |a, b| a.merge(b))?;

    pb.finish_with_message("| Finished counting");

    // Release tile-level inputs before global aggregation
    drop(tile_window_spans_for_threads);
    drop(tile_window_spans);
    drop(tiles);
    drop(windows_map);
    drop(blacklist_map);

    println!("Start: Processing counts");
    if !windows_aligned_to_tiles && !matches!(window_opt, WindowSpec::Global) {
        let (cross_sum, cross_weight) = stream_crossing_files(
            reduce_state.crossing_files.clone(),
            &zero_counts,
            &opt,
            avg_window_span,
        )?;
        reduce_state.scaled_sum.merge_from(&cross_sum)?;
        reduce_state.scaled_weight += cross_weight;
    }

    let keep_temp = false;
    if !keep_temp {
        if let Err(e) = std::fs::remove_dir_all(&temp_dir) {
            eprintln!(
                "warning: failed to remove temp dir {}: {}",
                temp_dir.display(),
                e
            );
        }
    } else {
        eprintln!("kept temp tiles in {}", temp_dir.display());
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
        println!("Start: Smoothing cfDNA GC counts");
        avg_gc_counts.smooth_length_rows_in_place(
            reference_metadata.smoothing_sigma,
            reference_metadata.smoothing_radius,
        );
    }

    // Convert GC counts to grid with lengths x GC percentage bins
    println!("Start: Converting cfDNA GC counts to GC percentage bins");
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
        &avg_gc_pct_counts,
        "avg_cfdna_counts",
        "average cfDNA counts",
    )?;

    println!("Start: Normalizing average counts");
    // Normalize GC counts array by its mean (just to remove the weighting scaling)
    // Ignores unsupported elements when calculating the mean
    let mut norm_gc_counts =
        mean_scale_array(&avg_gc_pct_counts, Some(&reference_outliers_support_mask))
            .expect("failed to perform masked mean scaling");

    intermediate_saver.save_file(
        &norm_gc_counts,
        "normalized_avg_cfdna_counts",
        "normalized average cfDNA counts",
    )?;

    // Interpolate counts for the unsupported cells
    if !reference_metadata.skip_interpolation {
        println!("Start: Interpolating counts for unsupported cells (very low reference counts)");
        for (row_idx, mut length_row) in norm_gc_counts.outer_iter_mut().enumerate() {
            // Rows are contiguous so we can safely borrow a mutable slice for interpolation
            let row_slice = length_row
                .as_slice_mut()
                .expect("GC histogram rows should be contiguous");
            let mut mask_row = reference_outliers_support_mask.row_mut(row_idx);
            let mask_slice = mask_row
                .as_slice_mut()
                .expect("Support mask rows should be contiguous");
            fill_unsupported_bins_with_polynomial(
                row_slice, mask_slice, 2, 3, 3,
                // Update mask when cells become supported
                // (TODO: not currently used downstream but worth having for future checks)
                true,
            )?;
        }

        intermediate_saver.save_file(
            &norm_gc_counts,
            "interpolated_cfdna_counts",
            "interpolated cfDNA counts",
        )?;
    }

    // Smoothe the *normalized* counts
    let do_smoothing = false;
    let smoothed_gc_counts = if do_smoothing {
        println!("Start: Smoothing counts with 2D Gaussian kernel");
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
        &binned_ref_counts,
        "binned_ref_counts",
        "binned reference counts",
    )?;

    intermediate_saver.save_file(
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
    println!(
        "Start: Setting counts to 1.0 for up to 2x{} extreme GC bins and {} shortest-length bins",
        opt.num_extreme_gc_bins, opt.num_short_length_bins
    );
    if correction_support_mask.iter().any(|&supported| !supported) {
        set_masked_entries_to_value(&mut norm_gc_counts, &correction_support_mask, 1.0);
        set_masked_entries_to_value(&mut norm_ref_counts, &correction_support_mask, 1.0);
    }

    intermediate_saver.save_file(
        &norm_gc_counts,
        "normalized_binned_cfdna_counts",
        "normalized binned cfDNA counts",
    )?;
    intermediate_saver.save_file(
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
            println!(
                "Start: Interpolating corrections for up to 2x{} extreme GC bins and {} shortest-length bins",
                opt.num_extreme_gc_bins, opt.num_short_length_bins
            );
            interpolate_masked_corrections(&mut norm_correction_matrix, &correction_support_mask)?;
        }

        if !matches!(outlier_rule, OutlierRule::None) {
            println!("Start: Applying outlier handling to correction matrix");
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

    // Save reusable correction package with metadata for downstream commands
    let correction_pkg = GCCorrectionPackage::from_components(
        GC_CORRECTION_SCHEMA_VERSION,
        &length_bins,
        &gc_bins,
        correction_matrix.clone(),
        length_bin_frequencies.clone(),
        &reference_metadata,
    )?;
    correction_pkg.write_npz(&opt.ioc.output_dir.join("gc_bias_correction.npz"))?;

    // Plot the avg. gc-bias across lengths for quick QC
    #[cfg(feature = "plotters")]
    {
        use crate::{
            commands::gc_bias::binning::compute_bin_edges,
            shared::plotters::{
                heatmap::{HeatmapFormat, write_heatmap},
                lineplot::write_line_plot_png,
            },
        };

        println!("Start: Plotting avg. bias across lengths");
        let gc_edges = compute_bin_edges(&gc_bins, 0, 100)?;
        let x_values: Vec<f64> = gc_edges
            .windows(2)
            .map(|window| {
                let start = window[0] as f64;
                let end = window[1] as f64;
                (start + end) / 2.0
            })
            .collect();
        let gc_edges_f: Vec<f64> = gc_edges.iter().map(|v| *v as f64).collect();
        let length_edges = compute_bin_edges(
            &length_bins,
            reference_metadata.min_fragment_length as u32,
            reference_metadata.max_fragment_length as u32,
        )?;
        let length_edges_f: Vec<f64> = length_edges.iter().map(|v| *v as f64).collect();

        // Unweighted average bias

        let num_gc_bins = correction_matrix.ncols();
        let bias_matrix = correction_matrix.mapv(|cf| if cf == 0.0 { 0.0 } else { 1.0 / cf });
        let mut unweighted_bias = vec![0.0; num_gc_bins];
        let mut unweighted_counts = vec![0usize; num_gc_bins];

        for length_biases in bias_matrix.outer_iter() {
            for (gc_idx, &bias) in length_biases.iter().enumerate() {
                if bias == 0.0 {
                    continue;
                }
                unweighted_bias[gc_idx] += bias;
                unweighted_counts[gc_idx] += 1;
            }
        }

        for (bias, count) in unweighted_bias.iter_mut().zip(unweighted_counts.iter()) {
            if *count > 0 {
                *bias /= *count as f64;
            }
        }

        // Weighted average bias

        let mut weighted_bias = vec![0.0; num_gc_bins];
        let mut weight_per_gc = vec![0.0; num_gc_bins];

        for (length_biases, &length_weight) in
            bias_matrix.outer_iter().zip(length_bin_frequencies.iter())
        {
            if length_weight == 0.0 {
                continue;
            }
            for (gc_idx, &bias) in length_biases.iter().enumerate() {
                if bias == 0.0 {
                    continue;
                }
                weight_per_gc[gc_idx] += length_weight;
                weighted_bias[gc_idx] += length_weight * bias;
            }
        }

        for (bias, weight) in weighted_bias.iter_mut().zip(weight_per_gc.iter()) {
            if *weight > 0.0 {
                *bias /= *weight;
            }
        }

        // Plot the bias

        let plot_path_unweighted = opt
            .ioc
            .output_dir
            .join("avg_gc_bias_across_lengths_unweighted.png");
        write_line_plot_png(
            &plot_path_unweighted,
            "Average GC bias across fragment lengths (unweighted)",
            "GC bin (%)",
            "GC bias",
            &x_values,
            &unweighted_bias,
            1600,
            1000,
        )
        .with_context(|| format!("writing GC bias plot to {}", plot_path_unweighted.display()))?;

        let plot_path_weighted = opt
            .ioc
            .output_dir
            .join("avg_gc_bias_across_lengths_weighted.png");
        write_line_plot_png(
            &plot_path_weighted,
            "Average GC bias across fragment lengths (weighted by length frequency)",
            "GC bin (%)",
            "GC bias",
            &x_values,
            &weighted_bias,
            1600,
            1000,
        )
        .with_context(|| format!("writing GC bias plot to {}", plot_path_weighted.display()))?;

        // Heatmap sizes
        let hm_width: u32 = 1000;
        let hm_height: u32 = 700;
        let scaling_factor = (hm_height as f32 / bias_matrix.nrows() as f32)
            .max(hm_width as f32 / bias_matrix.ncols() as f32)
            .ceil() as usize;

        // Heatmap of bias across length and GC contents
        let heatmap_path = opt.ioc.output_dir.join("gc_bias_heatmap.png");
        write_heatmap(
            &heatmap_path,
            "GC bias per length and GC %",
            "GC (%)",
            "Fragment length (bp)",
            &bias_matrix,
            Some(&gc_edges_f),
            Some(&length_edges_f),
            None,
            None,
            Some(1.0),
            None,
            None,
            true,
            scaling_factor,
            hm_width,
            hm_height,
            HeatmapFormat::Png,
        )
        .with_context(|| format!("writing GC bias heatmap to {}", heatmap_path.display()))?;

        // Heatmap of bias across length bins and GC bins
        let heatmap_path = opt.ioc.output_dir.join("gc_bias_heatmap.bins.png");
        write_heatmap(
            &heatmap_path,
            "GC bias per length bin and GC bin",
            "GC bin",
            "Fragment length bin",
            &bias_matrix,
            None,
            None,
            None,
            None,
            Some(1.0),
            None,
            None,
            true,
            scaling_factor,
            hm_width,
            hm_height,
            HeatmapFormat::Png,
        )
        .with_context(|| {
            format!(
                "writing GC bias heatmap (bins) to {}",
                heatmap_path.display()
            )
        })?;
    }

    println!("");
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
        "  Fragments counted one or more times: {}",
        global_counter.base.counted_fragments
    );
    println!(
        "  Windows processed: {} total, {} with counts",
        total_windows, counted_windows
    );
    if !matches!(outlier_rule, OutlierRule::None) {
        println!("  Outlier handling:");
        println!("    > Limits estimated from reference-supported bins only");
        println!(
            "    Supported cells examined: {} (winsorized: {})",
            outlier_stats.total_examined, outlier_stats.total_outliers_handled
        );
        println!("    > 'supported' = bins the reference marks valid (used to set limits)");
        println!(
            "    Unsupported cells examined: {} (winsorized: {})",
            outlier_stats.unsupported_examined, outlier_stats.unsupported_outliers_handled
        );
        println!(
            "    > 'unsupported' = bins the reference masks out (winsorized after interpolation)"
        );
        println!("    Clamped to [0.1,10.0]: {}", outlier_stats.hard_clamped);
    }
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn process_tile(
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    windows_aligned_to_tiles: bool,
    opt: &GCConfig,
    reference_metadata: &ReferenceGCMetadata,
    windows_opt: Option<&[(u64, u64, u64)]>,
    window_opt: &WindowSpec,
    avg_window_span: f64,
    template: &GCCounts,
    temp_dir: &PathBuf,
    blacklist_intervals: &[(u64, u64)],
) -> anyhow::Result<(WindowState, GCCounters)> {
    let apply_window_scaling = !matches!(window_opt, WindowSpec::Global);

    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = GCCounters::default();

    let mut out = WindowState::new(0, 0, 0, true, template)?;
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
        return Ok((out, counter));
    }
    let mut windows = prepared_windows.windows;
    let streaming_buffers = prepared_windows.streaming_buffers;

    let seq_start = tile.fetch_start as u64;
    let seq_end = tile.fetch_end.min(chrom_len as u32) as u64;

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
        apply_blacklist_mask_to_seq(&mut seq_bytes, &blacklist_intervals, seq_start);
        build_gc_prefixes(&seq_bytes)
    };

    // Decide whether we run streaming (fixed-size windows) or per-window Vec (BED/global)
    let using_streaming = streaming_buffers.is_some();

    let mut tile_window_intervals: Option<Vec<(u64, u64, u64)>> = None;
    if !using_streaming {
        if windows.is_empty() {
            return Ok((out, counter));
        }
        tile_window_intervals = Some(
            windows
                .iter()
                .enumerate()
                .map(|(local_idx, w)| (w.start, w.end, local_idx as u64))
                .collect(),
        );
    }

    // Fraction of a fragment that must overlap with a window to assign to that window
    // Keeps overlap assignment consistent between streaming and non-streaming paths
    let min_overlap_fraction: f64 = match opt.window_assignment.assign_by {
        WindowAssigner::Any | WindowAssigner::CountOverlap => {
            1. / (reference_metadata.max_fragment_length as f64 + 1.0)
        } // +1 to avoid rounding error issues
        WindowAssigner::All | WindowAssigner::Midpoint => {
            1.0 - (1. / (reference_metadata.max_fragment_length as f64 + 1.0))
        } // 1.0 but just below to avoid rounding errors
        WindowAssigner::Proportion(p) => p,
    };

    reader
        .fetch((tid, tile.fetch_start as i64, tile.fetch_end as i64))
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

            // Only count fragments whose start lies inside the tile core to avoid double counting
            if fragment.start < tile.core_start || fragment.start >= tile.core_end {
                continue;
            }

            // If the fragment is past the current window, finalize and advance the buffers
            while fragment.start as u64 >= current.end {
                finalize_window_buffer(
                    &mut current,
                    &gc_prefixes,
                    seq_start,
                    seq_end,
                    tile.core_start as u64,
                    tile.core_end as u64,
                    windows_aligned_to_tiles,
                    apply_window_scaling,
                    opt,
                    avg_window_span,
                    &mut out,
                    &mut crossing_parts,
                )?;
                let mut recycled = current;
                current = next;

                // Prepare the next buffer using the recycled allocation
                let next_idx = current.idx + 1;
                let (next_start, next_end) =
                    fixed_size_window_bounds(next_idx, window_bp, chrom_len);
                let next_contained =
                    next_start >= tile.core_start as u64 && next_end <= tile.core_end as u64;
                recycled.reset(next_idx, next_start, next_end, next_contained, template)?;
                next = recycled;
            }

            let Some(gc_count) = get_fragment_gc(
                &fragment,
                seq_start,
                seq_end,
                reference_metadata.end_offset as u32,
                &gc_prefixes,
                1.0,
            )?
            else {
                continue;
            };

            let (interval_start, interval_end) = match opt.window_assignment.assign_by {
                WindowAssigner::Midpoint => {
                    let midpoint =
                        midpoint_random_even_with_thread_rng(fragment.start, fragment_length);
                    (midpoint, midpoint + 1)
                }
                WindowAssigner::Any
                | WindowAssigner::All
                | WindowAssigner::Proportion(_)
                | WindowAssigner::CountOverlap => (fragment.start, fragment.end),
            };

            let fragment_span_length = (interval_end - interval_start) as f64;

            // Choose buffers to update: at most current and next for fixed-size windows
            let mut current_weight: Option<f64> = None;
            let overlap_current = overlap_length(
                interval_start as u64,
                interval_end as u64,
                current.start,
                current.end,
            );
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
            let overlap_next = overlap_length(
                interval_start as u64,
                interval_end as u64,
                next.start,
                next.end,
            );
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
            if let Some(count_weight) = next_weight {
                next.has_counts = true;
                next.counts
                    .incr_weighted(fragment_length as usize, gc_count, count_weight);
            }
        }

        // Flush both buffers so partially filled windows still reach downstream reducers
        finalize_window_buffer(
            &mut current,
            &gc_prefixes,
            seq_start,
            seq_end,
            tile.core_start as u64,
            tile.core_end as u64,
            windows_aligned_to_tiles,
            apply_window_scaling,
            opt,
            avg_window_span,
            &mut out,
            &mut crossing_parts,
        )?;
        finalize_window_buffer(
            &mut next,
            &gc_prefixes,
            seq_start,
            seq_end,
            tile.core_start as u64,
            tile.core_end as u64,
            windows_aligned_to_tiles,
            apply_window_scaling,
            opt,
            avg_window_span,
            &mut out,
            &mut crossing_parts,
        )?;
    } else {
        let mut window_ptr = 0usize; // Genomic window pointer reused by overlap finder
        {
            let tile_window_intervals = tile_window_intervals
                .as_ref()
                .expect("tile window intervals missing for non-streaming windows");

            // Iterate fragments and count GC
            for fragment_res in iter.by_ref() {
                let fragment = fragment_res.context("reading fragment")?;
                let fragment_length = fragment.len();

                // Only count fragments whose start lies inside the tile core to avoid double counting
                if fragment.start < tile.core_start || fragment.start >= tile.core_end {
                    continue;
                }

                let Some(gc_count) = get_fragment_gc(
                    &fragment,
                    seq_start,
                    seq_end,
                    reference_metadata.end_offset as u32,
                    &gc_prefixes,
                    1.0,
                )?
                else {
                    continue;
                };

                // Find all overlapping windows
                let (interval_start, interval_end) = match opt.window_assignment.assign_by {
                    WindowAssigner::Midpoint => {
                        let midpoint =
                            midpoint_random_even_with_thread_rng(fragment.start, fragment_length);
                        (midpoint, midpoint + 1)
                    }
                    WindowAssigner::Any
                    | WindowAssigner::All
                    | WindowAssigner::Proportion(_)
                    | WindowAssigner::CountOverlap => (fragment.start, fragment.end),
                };
                let overlapping_windows = find_overlapping_windows(
                    chrom_len,
                    &mut window_ptr,
                    Some(tile_window_intervals.as_slice()),
                    None,
                    interval_start.into(),
                    interval_end.into(),
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
                        WindowAssigner::CountOverlap => overlapped_window.overlap_fraction as f64,
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
                seq_start,
                seq_end,
                tile.core_start as u64,
                tile.core_end as u64,
                windows_aligned_to_tiles,
                apply_window_scaling,
                opt,
                avg_window_span,
                &mut out,
                &mut crossing_parts,
            )?;
        }
    }

    // Release prefix arrays before writing crossing parts
    drop(gc_prefixes);

    counter.add_from_snapshot(iter.counters_snapshot());

    // Write crossing parts when windows are not aligned
    // Streaming and non-streaming share the same writer
    out.crossing_file = if using_streaming {
        if !windows_aligned_to_tiles {
            write_crossing_parts(temp_dir, tile.index, template, &crossing_parts)?
        } else {
            None
        }
    } else if !windows_aligned_to_tiles && !matches!(window_opt, WindowSpec::Global) {
        write_crossing_parts(temp_dir, tile.index, template, &crossing_parts)?
    } else {
        None
    };

    Ok((out, counter))
}

pub fn process_window_in_place(
    gc_counts: &mut GCCounts,
    opt: &GCConfig,
    avg_window_size: Option<f64>,
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
    let avg_window_span = avg_window_size.expect("valid-positions needs window spans");

    // Mean scale, then multiply by the usable window-size (scaled to lower values)
    gc_counts.scale_counts((1. / mean_count) * (num_acgt as f64 / avg_window_span))?;

    Ok(true)
}

pub fn process_window(
    mut gc_counts: GCCounts,
    opt: &GCConfig,
    avg_window_size: Option<f64>,
) -> Result<Option<GCCounts>> {
    if process_window_in_place(&mut gc_counts, opt, avg_window_size)? {
        return Ok(Some(gc_counts));
    }
    Ok(None)
}

pub fn mean_scale_per_length_array<S, M>(
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

pub fn interpolate_masked_corrections(
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
            .expect("GC histogram rows should be contiguous");
        let mut mask_row = mask.row_mut(row_idx);
        let mask_slice = mask_row
            .as_slice_mut()
            .expect("Support mask rows should be contiguous");
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
pub fn mean_scale_array<S, M>(
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
    return Some(out);
}

/// Invert elements in an array (x) to `1 / x`, keeping 0s as 0.
fn invert_elementwise_with_zeros_inplace<S>(x: &mut ArrayBase<S, Ix2>)
where
    S: DataMut<Elem = f64>,
{
    x.mapv_inplace(|v| if v == 0.0 { 0.0 } else { 1.0 / v });
}

pub struct IntermediateFileSaver {
    pub save_intermediates: bool,
    pub out_dir: PathBuf,
    previously_saved: usize,
}

impl IntermediateFileSaver {
    pub fn new(save_intermediates: bool, out_dir: PathBuf) -> Self {
        IntermediateFileSaver {
            save_intermediates: save_intermediates,
            out_dir: out_dir,
            previously_saved: 0,
        }
    }

    pub fn save_file<S>(
        &mut self,
        x: &ArrayBase<S, Ix2>,
        file_tag: &str,
        msg_tag: &str,
    ) -> Result<()>
    where
        S: Data<Elem = f64>,
    {
        if self.save_intermediates {
            println!("Intermediate file: Saving {}", msg_tag);
            write_npy(
                &self.out_dir.join(format!(
                    "gc_bias.{}.{}.npy",
                    file_tag, self.previously_saved
                )),
                x,
            )
            .context(format!(
                "Failed to write intermediate file {}",
                self.previously_saved
            ))?;
            self.previously_saved += 1;
        }
        Ok(())
    }
}
