use crate::{
    commands::{
        cli_common::*,
        counters::GCCounters,
        gc_bias::{
            CORRECTION_CLAMP_RANGE, GC_CORRECTION_SCHEMA_VERSION,
            binning::{CollapseAggregation, bin_greedily_by_mass, collapse_counts_by_bins},
            config::{GCConfig, WindowWeightingSchemes},
            counting::{GCCounts, apply_gc_percent_width_correction, build_gc_prefixes},
            interpolation::fill_unsupported_bins_with_polynomial,
            load_reference_bias::{ReferenceGCData, ReferenceGCMetadata, load_reference_gc_data},
            outliers::{OutlierRule, OutlierStats, apply_outliers_to_matrix},
            package::GCCorrectionPackage,
            smoothing::smoothe_counts_gaussian,
            support_masking::{
                build_extreme_bins_support_mask, set_masked_entries_to_value, stats_by_support_mask,
            },
        },
    },
    shared::{
        bam::create_chromosome_reader,
        blacklist::apply_blacklist_mask_to_seq,
        fragment::minimal_fragment::Fragment,
        fragment_iterator::fragments_from_bam,
        midpoint::midpoint_random_even_with_thread_rng,
        overlaps::find_overlapping_windows,
        read::default_include_read,
        reference::read_seq,
        scale_genome::{compute_window_scaling_over_fragment, compute_window_scaling_over_overlap},
    },
};
use anyhow::{Context, Result, anyhow, bail, ensure};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array2, ArrayBase, Axis, Data, DataMut, Ix2, Zip};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{fs::create_dir_all, path::PathBuf, sync::Arc, time::Instant};

pub fn run(opt: &GCConfig) -> Result<()> {
    let start_time = Instant::now();
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, &opt.ioc.bam.as_path())?;
    let mut intermediate_saver =
        IntermediateFileSaver::new(opt.save_intermediates, opt.ioc.output_dir.clone());

    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    println!("Start: Loading reference GC bias");
    let ReferenceGCData {
        window_spec: window_opt,
        windows_map,
        window_indices_by_chr,
        counts: reference_counts,
        unobservables_support_mask: reference_unobservables_support_mask,
        outliers_support_mask: mut reference_outliers_support_mask,
        avg_window_size,
        gc_percent_widths: reference_gc_percent_widths,
        metadata: reference_metadata,
    } = load_reference_gc_data(
        &opt.ref_gc_dir,
        Some(&chromosomes),
        100 - opt.min_window_acgt_pct,
    )?;

    if matches!(window_opt, WindowSpec::Global)
        && matches!(opt.window_weighting, WindowWeightingSchemes::ValidPositions)
    {
        bail!(
            "Window weighting scheme 'valid-positions' requires genomic windows. \
             It cannot be used when running in global mode"
        );
    }

    // Load genomic scaling factors
    if opt.scale_genome.scaling_factors.is_some() {
        println!("Start: Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> =
        load_scaling_map(&opt.scale_genome, &chromosomes, &contigs)?;

    // Create output directory
    create_dir_all(&opt.ioc.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = load_blacklist_map(opt.blacklist.as_ref(), 1, 0, &chromosomes)?;

    // Configure global thread‐pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.ioc.n_threads as usize)
        .build_global()
        .context("building Rayon thread pool")?;

    // Prepare per-bin counts and metadata
    let mut all_bins = Vec::new();
    let mut global_counter = GCCounters::default();

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let counts_by_window: Vec<(Vec<GCCounts>, GCCounters)> = chromosomes
        .par_iter()
        .map(|chr| -> Result<_, _> {
            let out = process_chrom(
                &chr,
                opt,
                &reference_metadata,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                &window_opt,
                scaling_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    println!("Start: Processing counts");

    // Collect results (in chromosome order) back into the global vectors
    for (counts_by_bin, counter) in counts_by_window {
        all_bins.extend(counts_by_bin);
        global_counter += counter;
    }

    // Convert to single `GCCounts` for global
    // Keep wrapped in vector to simplify writer
    let mut all_bins = if matches!(window_opt, WindowSpec::Global) {
        vec![GCCounts::collapse(&all_bins)?]
    } else {
        all_bins
    };

    // Prepare each window for averaging

    // First, we make a map of windows so we can get the correct reference bias windows
    // (not guaranteed to be sorted the same)
    let avg_counts_tuple: (Array2<f64>, GCCounts) = if let Some(ws_by_chr) = window_indices_by_chr {
        let mut window_tuples: Vec<(usize, usize)> = vec![];
        let mut win_idx: usize = 0;
        for chrom in &chromosomes {
            let chrom_orig_indices = ws_by_chr.get(chrom).unwrap();
            for orig_idx in chrom_orig_indices.iter() {
                window_tuples.push((win_idx, *orig_idx as usize));
                win_idx += 1;
            }
        }

        ensure!(
            all_bins.len() == window_tuples.len(),
            "window/count mismatch"
        );

        println!("Start: Processing counts per windows");

        pb.reset();
        pb.set_length(win_idx as u64);
        pb.set_position(0);

        // Run preparation on windows in parallel
        let counts_tuples: Vec<Option<(Array2<f64>, GCCounts)>> = all_bins
            .into_par_iter()
            .enumerate()
            .map(|(count_idx, gc_counts)| -> Result<_, _> {
                let ref_idx = window_tuples[count_idx].1;
                let ref_counts_view = reference_counts.index_axis(Axis(0), ref_idx);
                let out = process_window(
                    gc_counts,
                    &ref_counts_view,
                    &reference_outliers_support_mask,
                    &opt,
                    avg_window_size,
                )?;
                pb.inc(1);
                Ok(out)
            })
            .collect::<Result<_>>()?; // short-circuits on the first Err

        pb.finish_with_message("| Finished processing");

        // Unpack to separate vectors that have the same non-None windows left
        let mut ref_counts = Vec::new();
        let mut gc_counts = Vec::new();
        for opt in counts_tuples {
            if let Some((arr, gc)) = opt {
                ref_counts.push(arr);
                gc_counts.push(gc);
            }
        }
        ensure!(
            !ref_counts.is_empty(),
            "No GC bias windows produced usable counts. \
                Check settings such as `--min-window-acgt-pct` relative many positions are blacklisted. \
                To limited fragment lengths or GC content ranges may also produce this problem."
        );

        // Get the average ref count arrays and cfDNA GC counts
        let num_windows_left = gc_counts.len();
        let avg_ref_counts = mean_of_arrays(&ref_counts)?;
        let mut avg_gc_counts = GCCounts::collapse(&gc_counts)?;
        avg_gc_counts.scale_counts(1.0 / num_windows_left as f64)?;

        (avg_ref_counts, avg_gc_counts)
    } else {
        // Global window
        let gc_counts = all_bins.remove(0);
        let ref_counts_view = reference_counts.index_axis(Axis(0), 0);
        process_window(
            gc_counts,
            &ref_counts_view,
            &reference_outliers_support_mask,
            &opt,
            avg_window_size)?.ok_or_else(|| {
            anyhow!(
                "Produced no usable GC bias counts. \
                Check settings such as `--min-window-acgt-pct` relative many positions are blacklisted. \
                To limited fragment lengths or GC content ranges may also produce this problem."
            )
        })?
    };

    let (avg_norm_ref_counts, mut avg_gc_counts) = avg_counts_tuple;

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
    println!("Start: Smoothing counts with 2D Gaussian kernel");
    let do_smoothing = false;
    let smoothed_gc_counts = if do_smoothing {
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

    // a) Collapse row-mean-scaled reference counts into the length and GC bins
    // We average the values at the collapsed indices, weighted by the occurence of the lengths in the cfDNA.
    // b) Collapse cfDNA counts into the length and GC bins.
    // We sum them to answer the question: "What is the probability of seeing this bin".
    // TODO: Since only the length dimension is normalized, these operations make a big difference. Reconsider them carefully!
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

        // Sanity clamp of corrections
        norm_correction_matrix =
            norm_correction_matrix.clamp(CORRECTION_CLAMP_RANGE.0, CORRECTION_CLAMP_RANGE.1);

        // Re-normalize correction matrix per fragment length to be centered around 1.0
        // Still ignores extreme GC bins in the mean-calculations
        let mut norm_correction_matrix = mean_scale_per_length_array(
            &norm_correction_matrix,
            0.,
            Some(&correction_support_mask),
        );

        // Make correction factors multipliers by inverting elements to 1 / x
        // Zeros remain 0s
        invert_elementwise_with_zeros_inplace(&mut norm_correction_matrix);

        Ok::<Array2<f64>, anyhow::Error>(norm_correction_matrix)
    }?;

    // Length-bin frequencies (normalized) used for length-agnostic GC correction
    let length_bin_frequencies = {
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
        length_bin_frequencies,
        &reference_metadata,
    )?;
    correction_pkg.write_npz(&opt.ioc.output_dir.join("gc_bias_correction.npz"))?;

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
    if !matches!(outlier_rule, OutlierRule::None) {
        println!(
            "  Outlier handling: examined {} ({} adjusted), unsupported examined {} ({} adjusted)",
            outlier_stats.total_examined,
            outlier_stats.total_outliers_handled,
            outlier_stats.unsupported_examined,
            outlier_stats.unsupported_outliers_handled
        );
    }
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn process_chrom(
    chr: &str,
    opt: &GCConfig,
    reference_metadata: &ReferenceGCMetadata,
    windows_opt: Option<&[(u64, u64, u64)]>,
    window_opt: &WindowSpec,
    scaling_chr: &[(u64, u64, f32)],
    // gc_bins: usize,
    blacklist_intervals: &[(u64, u64)],
) -> anyhow::Result<(Vec<GCCounts>, GCCounters)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    let mut seq_bytes = read_seq(&opt.ref_genome.ref_2bit, chr)?;
    apply_blacklist_mask_to_seq(&mut seq_bytes, &blacklist_intervals, 0);

    let gc_prefixes = build_gc_prefixes(&seq_bytes);

    // Initialize counters (default -> 0s)
    let mut counter = GCCounters::default();

    let bed_windows: Option<&[(u64, u64, u64)]> = match window_opt {
        WindowSpec::Bed(_) => match windows_opt {
            Some(slice) => {
                if slice.is_empty() {
                    return Ok((Vec::new(), counter));
                }
                Some(slice)
            }
            None => {
                bail!("Window specification is BED, but no windows provided for chromosome {chr}");
            }
        },
        WindowSpec::Global => None,
        _ => unreachable!("Only Bed and Global are used"),
    };

    let num_bins = match window_opt {
        WindowSpec::Bed(_) => bed_windows.unwrap().len(),
        WindowSpec::Global => 1,
        _ => unreachable!("Only `Bed` and `Global` are used"),
    };

    // Initialize count arrays
    let mut counts_by_bin = vec![
        GCCounts::new(
            reference_metadata.min_fragment_length,
            reference_metadata.max_fragment_length,
            reference_metadata.end_offset as usize,
            (0, 0)
        )?;
        num_bins
    ];

    // Count the number of ACGT bases and total number of bases per window
    if let Some(win_coords) = bed_windows {
        for (bin_idx, (start, end, _)) in win_coords.iter().enumerate() {
            let acgt_count = gc_prefixes.acgt[*end as usize] - gc_prefixes.acgt[*start as usize];
            counts_by_bin[bin_idx].num_acgt_out_of = (acgt_count as u64, end - start);
        }
    } else {
        let total_acgt = *gc_prefixes.acgt.last().expect("prefix has sentinel") as u64;
        counts_by_bin[0].num_acgt_out_of = (total_acgt, chrom_len);
    }

    // Fraction of a fragment that must overlap with a window to assign to that window
    let min_overlap_fraction: f64 = match opt.window_assignment.assign_by {
        WindowAssigner::Any | WindowAssigner::CountOverlap => {
            1. / (reference_metadata.max_fragment_length as f64 + 1.0)
        } // +1 to avoid rounding error issues
        WindowAssigner::All | WindowAssigner::Midpoint => {
            1.0 - (1. / (reference_metadata.max_fragment_length as f64 + 1.0))
        } // 1.0 but just below to avoid rounding errors
        WindowAssigner::Proportion(p) => p,
    };

    // Replace scaling factor with unused index (for compatibility with overlap finder)
    let scaling_with_bin_idx: Vec<(u64, u64, u64)> =
        scaling_chr.iter().map(|(s, e, _)| (*s, *e, 0u64)).collect();

    // Streaming pointers
    let mut wd_ptr = 0; // Genomic window
    let mut sf_ptr = 0; // Scaling factor bin

    // Get coordinates to fetch reads from and to
    let (fetch_from, fetch_to) = match window_opt {
        WindowSpec::Bed(_) => {
            let wn = bed_windows.expect("validated above");
            let fetch_start = wn[0].0 as i64;
            let fetch_end = wn.iter().map(|w| w.1).max().unwrap() as i64;
            (
                (fetch_start - reference_metadata.max_fragment_length as i64).max(0i64),
                (fetch_end + reference_metadata.max_fragment_length as i64).min(chrom_len as i64),
            )
        }
        _ => (0i64, chrom_len as i64),
    };

    reader
        .fetch((tid, fetch_from, fetch_to))
        .context(format!("fetch {}", chr))?;

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

    // Wrap to use opt
    let include_read_fn = {
        let opt = (*opt).clone();
        move |r: &Record| default_include_read(r, opt.require_proper_pair, opt.min_mapq)
    };

    // Create fragment iterator
    let mut iter = fragments_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        include_read_fn,
        None,
        fragment_filter,
    )
    .with_local_counters();

    // Convert variables once
    let end_offset_u32 = reference_metadata.end_offset as u32;
    let min_acgt_fraction = 1.0f32;

    // Iterate fragments and count GC
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;
        let fragment_length = fragment.len();

        // Apply fragment end offsets
        let gc_window_start = fragment.start.saturating_add(end_offset_u32);
        let gc_window_end = fragment.end.saturating_sub(end_offset_u32);
        if gc_window_end <= gc_window_start {
            continue;
        }

        // Extract GC count in the interval. We store raw GC counts per length.
        let gc_count = {
            let acgt = gc_prefixes.acgt[gc_window_end as usize]
                - gc_prefixes.acgt[gc_window_start as usize];
            if acgt < 10
                || (acgt as f32 / (gc_window_end - gc_window_start) as f32) < min_acgt_fraction
            {
                continue;
            }

            let gc =
                gc_prefixes.gc[gc_window_end as usize] - gc_prefixes.gc[gc_window_start as usize];
            ensure!(
                gc <= (gc_window_end - gc_window_start),
                "GC count exceeded interval length: {} > {}",
                gc,
                gc_window_end - gc_window_start
            );
            gc as usize
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
            &mut wd_ptr,
            bed_windows,
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

        // Find all overlapping scaling-factor bins
        // And count up the weight
        if !scaling_chr.is_empty() {
            // Find overlapping scaling-bins
            let overlapping_scaling_bins = find_overlapping_windows(
                chrom_len,
                &mut sf_ptr,
                Some(&scaling_with_bin_idx),
                None,
                fragment.start.into(), // Full fragment
                fragment.end.into(),
                1. / (reference_metadata.max_fragment_length as f64 + 1.0), // Any overlap
                reference_metadata.max_fragment_length as u64,
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
            let overlap_weights = match opt.window_assignment.assign_by {
                WindowAssigner::CountOverlap => compute_window_scaling_over_overlap(
                    &overlapping_windows,
                    &overlapping_scaling_bin_indices,
                    scaling_chr,
                )?,
                _ => compute_window_scaling_over_fragment(
                    &overlapping_windows,
                    &overlapping_scaling_bin_indices,
                    scaling_chr,
                )?,
            };

            // Count up the weight per overlapping count-window
            for (overlapped_window_idx, scaling_weight, overlap_fraction_to_count) in
                overlap_weights
            {
                counts_by_bin[overlapped_window_idx].incr_weighted(
                    fragment_length as usize,
                    gc_count,
                    overlap_fraction_to_count * scaling_weight,
                );
            }
        } else {
            // When no scaling, increment counter by 1.0 or by the overlap fraction
            // NOTE: incrementer handles min-gc-pct offsetting!
            for overlapped_window in overlapping_windows.windows {
                let count_weight = match opt.window_assignment.assign_by {
                    WindowAssigner::CountOverlap => overlapped_window.overlap_fraction as f64,
                    _ => 1.0f64,
                };
                counts_by_bin[overlapped_window.idx].incr_weighted(
                    fragment_length as usize,
                    gc_count,
                    count_weight,
                );
            }
        }
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    Ok((counts_by_bin, counter))
}

/// Scale the GC counts and reference counts
/// to ensure the averaging follows the weighting scheme.
///
/// ## Support masks
///
///  - The unobservables is a theoretical mask. We use it to check that no bins that are in-theory
///    impossible to observe (based on the combination of GC and length) get non-zero counts.
///
///  - The outliers are very rarely occuring combinations of GC and length that we don't have
///    enough support for in the reference counts. We avoid these elements when calculating
///    fragment length means (ref counts) for normalization.
///
/// ## Pipeline
///
/// * Reference counts are normalized (mean-scaled) per fragment length.
///
/// * The two arrays are weighted depending on the weighting scheme:
///
///     - `Coverage`: The normalized reference counts are multiplied by the overall mean GC count,
///       ensuring that GC counts and reference counts are weighted the same in the global matrices after averaging.
///       The GC counts already reflect the coverage.
///
///     - `ValidPositions`: The GC counts are divided by the overall mean GC count.
///       Both arrays are then multiplied by the number of ACGT bases in the window and
///       divided by the average window size (to avoid exploding the count sizes).
///
///     - `Equal`: The GC counts are divided by the overall mean GC count,
///       ensuring GC counts have the same weight for all windows.
///       The reference count normalization made those matrices equally weighted.
fn process_window<S, M>(
    mut gc_counts: GCCounts,
    ref_counts: &ArrayBase<S, Ix2>,
    ref_support_mask_outliers: &ArrayBase<M, Ix2>,
    opt: &GCConfig,
    avg_window_size: Option<f64>,
) -> Result<Option<(Array2<f64>, GCCounts)>>
where
    S: Data<Elem = f64>,
    M: Data<Elem = bool>,
{
    // Check window has enough valid positions
    if gc_counts.pct_acgt() < opt.min_window_acgt_pct as f64 {
        return Ok(None);
    }
    if gc_counts.sum() == 0.0 {
        return Ok(None);
    }

    // Get mean coverage to scale by (supported cells only)
    let mean_count: f64 = gc_counts.mean();

    // Normalize the reference bias
    // Cast to f64 once
    let ref_counts_f = ref_counts.mapv(|v| v as f64);
    // Row-wise mean-scaling
    let ref_counts_norm =
        mean_scale_per_length_array(&ref_counts_f, 0.0, Some(ref_support_mask_outliers));

    // Weighting one or both count distributions
    let weighted_ref_counts = match opt.window_weighting {
        WindowWeightingSchemes::Equal => {
            // Make sure both arrays normalized
            gc_counts.scale_counts(1. / mean_count)?;
            ref_counts_norm
        }
        WindowWeightingSchemes::Coverage => {
            // Scale normalized ref counts by average coverage
            ref_counts_norm * mean_count
        }
        WindowWeightingSchemes::ValidPositions => {
            let num_acgt = gc_counts.num_acgt_out_of.0;
            if num_acgt == 0 {
                // No positions observed
                return Ok(None);
            }
            // Use avg window size to "scale the scaling"
            // So the counts don't explode in size
            // While keeping the relative scale between windows consistent
            let avg_window_span = avg_window_size.expect("valid-positions needs window spans");

            gc_counts.scale_counts(1. / (mean_count * (num_acgt as f64 / avg_window_span)))?;
            let ref_counts_norm_scaled = ref_counts_norm * (num_acgt as f64 / avg_window_span);
            ref_counts_norm_scaled
        }
    };

    Ok(Some((weighted_ref_counts, gc_counts)))
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

fn mean_of_arrays(arrs: &[Array2<f64>]) -> anyhow::Result<Array2<f64>> {
    let mut iter = arrs.iter();
    let mut sum = iter
        .next()
        .cloned()
        .ok_or_else(|| anyhow::anyhow!("mean_of_arrays requires at least one array"))?;
    let mut count = 1usize;

    for a in iter {
        ensure!(sum.dim() == a.dim(), "array shape mismatch");
        sum += a;
        count += 1;
    }

    let scale = 1.0 / count as f64;
    sum.mapv_inplace(|v| v * scale);
    Ok(sum)
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
