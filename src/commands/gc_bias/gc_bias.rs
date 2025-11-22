use crate::{
    commands::{
        cli_common::*,
        counters::GCCounters,
        gc_bias::{
            config::{GCConfig, WindowWeightingSchemes},
            counting::{GCCounts, build_gc_prefixes, get_gc_integer_percentage_for_window},
            load_reference_bias::{ReferenceGCData, load_reference_gc_data},
            smoothing::smoothe_counts_gaussian,
        },
        reference_gc::interpolation::fill_unsupported_bins_with_polynomial,
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
use ndarray::{Array1, Array2, ArrayBase, ArrayView2, Axis, Data, Ix2, Zip};
use ndarray_npy::{NpzWriter, write_npy};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{
    fs::{File, create_dir_all},
    path::PathBuf,
    sync::Arc,
    time::Instant,
};

const CORRECTION_CLAMP_RANGE: (f64, f64) = (0.1, 10.0);
const GC_CORRECTION_SCHEMA_VERSION: u32 = 1;

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
    let all_bins = if matches!(window_opt, WindowSpec::Global) {
        vec![GCCounts::collapse(&all_bins)?]
    } else {
        all_bins
    };

    // Prepare each window for averaging

    // First, we make a map of windows so we can get the correct reference bias windows
    // (not guaranteed to be sorted the same)
    let avg_counts_tuple = if let Some(ws_by_chr) = window_indices_by_chr {
        let mut window_tuples: Vec<(usize, usize)> = vec![];
        let mut win_idx: usize = 0;
        for chrom in &chromosomes {
            let chrom_orig_indices = ws_by_chr.get(chrom).unwrap();
            for orig_idx in chrom_orig_indices.iter() {
                window_tuples.push((win_idx, *orig_idx as usize));
                win_idx += 1;
            }
        }

        println!("Start: Processing counts per windows");

        pb.reset();
        pb.set_length(win_idx as u64);
        pb.set_position(0);

        // Run preparation on windows in parallel
        let avg_counts_tuples: Vec<Option<(Array2<f64>, Array2<f64>)>> = window_tuples
            .par_iter()
            .map(|(count_idx, ref_idx)| -> Result<_, _> {
                let gc_counts = &all_bins[*count_idx];
                let ref_counts_view = reference_counts.index_axis(Axis(0), *ref_idx);
                let out = process_window(
                    gc_counts,
                    &ref_counts_view,
                    &reference_unobservables_support_mask,
                    &reference_outliers_support_mask,
                    &opt,
                    avg_window_size,
                )?;
                pb.inc(1);
                Ok(out)
            })
            .collect::<Result<_>>()?; // short-circuits on the first Err

        pb.finish_with_message("| Finished processing");

        mean_of_arrays(&avg_counts_tuples).ok_or_else(|| {
            anyhow!(
                "No GC bias windows produced usable counts. \
                Check settings such as `--min-window-acgt-pct` relative many positions are blacklisted. \
                To limited fragment lengths or GC content ranges may also produce this problem."
            )
        })?
    } else {
        // Global window
        let gc_counts = &all_bins[0];
        let ref_counts_view = reference_counts.index_axis(Axis(0), 0);
        process_window(gc_counts, &ref_counts_view,&reference_unobservables_support_mask,
                    &reference_outliers_support_mask, &opt, avg_window_size)?.ok_or_else(|| {
            anyhow!(
                "Produced no usable GC bias counts. \
                Check settings such as `--min-window-acgt-pct` relative many positions are blacklisted. \
                To limited fragment lengths or GC content ranges may also produce this problem."
            )
        })?
    };

    let (avg_gc_counts, avg_norm_ref_counts) = avg_counts_tuple;

    println!("Start: Normalizing average counts");
    // Normalize GC counts array by its mean (just to remove the weighting scaling)
    // Ignores unsupported elements when calculating the mean
    let mut norm_gc_counts =
        mean_scale_array(&avg_gc_counts, Some(&reference_outliers_support_mask))
            .expect("failed to perform masked mean scaling");

    intermediate_saver.save_file(
        &norm_gc_counts,
        "normalized_avg_cfdna_counts",
        "normalized average cfDNA counts",
    )?;

    // Interpolate counts for the unsupported cells
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

    // Smoothe the *normalized* counts
    println!("Start: Smoothing counts with 2D Gaussian kernel");
    let smoothed_gc_counts = {
        // 5-element kernel (-2...+2)
        let radius: usize = 2;
        // Standard deviation (quite sharp so not too smoothed)
        let sigma = 0.5;
        smoothe_counts_gaussian(&norm_gc_counts, sigma, radius)
    };

    intermediate_saver.save_file(
        &smoothed_gc_counts,
        "smoothed_cfdna_counts",
        "smoothed cfDNA counts",
    )?;

    // Get greedy bins for lengths and GC
    // Maps "length -> length bin" and "gc -> gc bin"
    let length_bins = bin_greedily_by_mass(&smoothed_gc_counts, 0, opt.min_length_bin_mass as f64)?;
    let gc_bins = bin_greedily_by_mass(&smoothed_gc_counts, 1, opt.min_gc_bin_mass as f64)?;

    // Collapse row-mean-scaled reference counts into the length and GC bins
    // We *average* the values at the collapsed indices. For length bin collapsing,
    // we weight the average by the occurence of the lengths in the cfDNA
    let binned_ref_counts = {
        let length_binned = collapse_counts_by_bins(
            &avg_norm_ref_counts,
            0,
            &length_bins,
            CollapseAggregation::Mean,
            // Weight the average by how common the lengths are in the cfDNA
            Some(smoothed_gc_counts.view()),
        )?;
        collapse_counts_by_bins(&length_binned, 1, &gc_bins, CollapseAggregation::Mean, None)?
    };

    intermediate_saver.save_file(
        &binned_ref_counts,
        "binned_ref_counts",
        "binned reference counts",
    )?;

    // Collapse GC counts into the length and GC bins
    // We sum the values at the collapsed indices
    let binned_gc_counts = {
        let length_binned = collapse_counts_by_bins(
            &smoothed_gc_counts,
            0,
            &length_bins,
            CollapseAggregation::Sum,
            None,
        )?;
        collapse_counts_by_bins(&length_binned, 1, &gc_bins, CollapseAggregation::Sum, None)?
    };

    intermediate_saver.save_file(
        &binned_gc_counts,
        "binned_cfdna_counts",
        "binned cfDNA counts",
    )?;

    // Mask extreme GC bins to avoid unstable corrections
    let correction_support_mask =
        build_extreme_gc_support_mask(binned_gc_counts.dim(), opt.num_extreme_gc_bins as usize);
    let mut norm_gc_counts =
        mean_scale_per_length_array(&binned_gc_counts, 0., Some(&correction_support_mask));
    let mut norm_ref_counts =
        mean_scale_per_length_array(&binned_ref_counts, 0., Some(&correction_support_mask));

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

    // Set extreme GC bins to 1.0 in both arrays to avoid zero-division etc.
    println!("Start: Setting extreme GC bins to 1.0");
    if correction_support_mask.iter().any(|&supported| !supported) {
        set_masked_entries_to_value(&mut norm_gc_counts, &correction_support_mask, 1.0);
        set_masked_entries_to_value(&mut norm_ref_counts, &correction_support_mask, 1.0);
    }

    // Calculate correction matrix
    // 1) Divide cfDNA counts by reference counts
    // 2) Normalize each fragment length to mean=1.0
    // 3) Set final correction factors for extreme GC bins to 1.0 (no corrections)
    // 4) Clamp any leftover extreme corrections
    let correction_matrix = {
        let raw_correction_matrix = &norm_gc_counts / &norm_ref_counts;

        // Normalize correction matrix per fragment length to be centered around 1.0
        // Still ignore extreme GC bins in the mean-calculations
        let mut norm_correction_matrix =
            mean_scale_per_length_array(&raw_correction_matrix, 0., Some(&correction_support_mask));

        // Refill extreme GC bins to 1.0
        if correction_support_mask.iter().any(|&supported| !supported) {
            set_masked_entries_to_value(&mut norm_correction_matrix, &correction_support_mask, 1.0);
        }

        // Sanity clamp of corrections
        norm_correction_matrix.clamp(CORRECTION_CLAMP_RANGE.0, CORRECTION_CLAMP_RANGE.1)
    };

    // Save reusable correction package with metadata for downstream commands
    let correction_pkg = GCCorrectionPackage::from_components(
        GC_CORRECTION_SCHEMA_VERSION,
        &length_bins,
        &gc_bins,
        correction_matrix.clone(),
        opt,
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
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn process_chrom(
    chr: &str,
    opt: &GCConfig,
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
            0,
            100,
            opt.fragment_lengths.min_fragment_length as usize,
            opt.fragment_lengths.max_fragment_length as usize,
            (0, 0)
        );
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
            1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0)
        } // +1 to avoid rounding error issues
        WindowAssigner::All | WindowAssigner::Midpoint => {
            1.0 - (1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0))
        } // 1.0 but just below to avoid rounding errors
        WindowAssigner::Proportion(p) => p,
    };

    // Replace scaling factor with unused index (for compatibility with overlap finder)
    let scaling_with_bin_idx: Vec<(u64, u64, u64)> =
        scaling_chr.iter().map(|(s, e, _)| (*s, *e, 0u64)).collect();

    // Streaming pointers and single fetch for this chr
    let mut wd_ptr = 0; // Genomic window
    let mut sf_ptr = 0; // Scaling factor bin

    // Get coordinates to fetch reads from and to
    let (fetch_from, fetch_to) = match window_opt {
        WindowSpec::Bed(_) => {
            let wn = bed_windows.expect("validated above");
            let fetch_start = wn[0].0 as i64;
            let fetch_end = wn.iter().map(|w| w.1).max().unwrap() as i64;
            (
                (fetch_start - opt.fragment_lengths.max_fragment_length as i64).max(0i64),
                (fetch_end + opt.fragment_lengths.max_fragment_length as i64).min(chrom_len as i64),
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
        let min_len = opt.fragment_lengths.min_fragment_length;
        let max_len = opt.fragment_lengths.max_fragment_length;
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
        fragment_filter,
    )
    .with_local_counters();

    // Convert variables once
    let end_offset_u32 = opt.end_offset as u32;
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

        // Extract GC fraction in the interval
        let gc_bin = get_gc_integer_percentage_for_window(
            &gc_prefixes,
            gc_window_start as usize,
            gc_window_end as usize,
            min_acgt_fraction,
            1, // Fraction is 100% so this is ignored
        );

        // Unpack GC fraction (or continue)
        let gc_bin = match gc_bin {
            Some(v) => v,
            None => continue,
        };

        ensure!(
            (0..=100).contains(&gc_bin),
            "GC fraction out of [0,100]: {}",
            gc_bin
        );

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
            opt.fragment_lengths.max_fragment_length.into(),
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
                    gc_bin,
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
                    gc_bin,
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
    gc_counts: &GCCounts,
    ref_counts: &ArrayBase<S, Ix2>,
    ref_support_mask_unobservables: &ArrayBase<M, Ix2>,
    ref_support_mask_outliers: &ArrayBase<M, Ix2>,
    opt: &GCConfig,
    avg_window_size: Option<f64>,
) -> Result<Option<(Array2<f64>, Array2<f64>)>>
where
    S: Data<Elem = f64>,
    M: Data<Elem = bool>,
{
    // Check window has enough valid positions
    if gc_counts.pct_acgt() < opt.min_window_acgt_pct as f64 {
        return Ok(None);
    }

    // Get counts as 2d array with shape: (n_lengths, n_gc_bins)
    let counts_mat = gc_counts.to_array2();

    // Find total counts
    let total_count: f64 = counts_mat.sum();

    if total_count == 0.0 {
        return Ok(None);
    }

    // Get counts and value-sums for supported and unsupported cells
    let stats_by_support_status =
        stats_by_support_mask(&counts_mat, &ref_support_mask_unobservables);

    // Our assumption is that unsupported bins are impossible to see
    // when all fragment positions have valid ACGT bases (in the reference genome)
    // If this is wrong, we need to reconsider these assumptions!
    if stats_by_support_status.sum_for_unsupported > 0.0 {
        bail!("Unsupported bins in the count matrix had non-zero coverage. Report please.");
    }

    // Get mean coverage to scale by (supported cells only)
    let mean_count: f64 =
        stats_by_support_status.sum_for_supported / stats_by_support_status.n_supported as f64;

    // Normalize the reference bias
    // Cast to f64 once
    let ref_counts_f = ref_counts.mapv(|v| v as f64);
    // Row-wise mean-scaling
    let ref_counts_norm =
        mean_scale_per_length_array(&ref_counts_f, 0.0, Some(ref_support_mask_outliers));

    // Weighting one or both count distributions
    let (weighted_counts, weighted_ref_counts) = match opt.window_weighting {
        WindowWeightingSchemes::Equal => {
            // Return both arrays normalized
            let counts_mat_norm = counts_mat / mean_count;
            (counts_mat_norm, ref_counts_norm)
        }
        WindowWeightingSchemes::Coverage => {
            // Return GC counts as are and scale normalized ref counts by average coverage
            (counts_mat, ref_counts_norm * mean_count)
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
            let counts_mat_norm_scaled =
                counts_mat / mean_count * (num_acgt as f64 / avg_window_span);
            let ref_counts_norm_scaled = ref_counts_norm * (num_acgt as f64 / avg_window_span);
            (counts_mat_norm_scaled, ref_counts_norm_scaled)
        }
    };

    Ok(Some((weighted_counts, weighted_ref_counts)))
}

fn mean_scale_per_length_array<S, M>(
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

pub struct StatsBySupportMask {
    pub sum_for_supported: f64,
    pub sum_for_unsupported: f64,
    pub n_supported: u64,
    pub n_unsupported: u64,
}

/// Get count and value-sums for all supported/unsupported bins.
pub fn stats_by_support_mask<S, M>(
    matrix: &ArrayBase<S, Ix2>,
    support_mask: &ArrayBase<M, Ix2>,
) -> StatsBySupportMask
where
    S: Data<Elem = f64>,
    M: Data<Elem = bool>,
{
    assert_eq!(
        matrix.dim(),
        support_mask.dim(),
        "Mask shape {:?} must match matrix shape {:?}",
        support_mask.dim(),
        matrix.dim()
    );

    let mut total_supported = 0.0;
    let mut total_unsupported = 0.0;
    let mut count_supported = 0;
    let mut count_unsupported = 0;

    Zip::from(matrix)
        .and(support_mask)
        .for_each(|value, &is_supported| {
            if is_supported {
                total_supported += *value;
                count_supported += 1;
            } else {
                total_unsupported += *value;
                count_unsupported += 1;
            }
        });

    StatsBySupportMask {
        sum_for_supported: total_supported,
        sum_for_unsupported: total_unsupported,
        n_supported: count_supported,
        n_unsupported: count_unsupported,
    }
}

// Overall scaling
// Elements that are marked as `false` in the support mask are
// still scaled but do not contribute to the mean
fn mean_scale_array<S, M>(
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

/// Average a list of optional `(matrix_a, matrix_b)` pairs independently.
///
/// Each `Some((A, B))` contributes to two separate accumulators: one for every `A`
/// (the first matrix) and one for every `B` (the second matrix). Entries that are
/// `None` are skipped entirely. All present matrices within a given slot must share
/// the same shape; otherwise the function panics in debug builds. Returns `None` if
/// no `Some` elements exist; otherwise returns the pair of mean matrices.
fn mean_of_arrays(
    arrays: &[Option<(Array2<f64>, Array2<f64>)>],
) -> Option<(Array2<f64>, Array2<f64>)> {
    let mut iter = arrays.iter().filter_map(|opt| opt.as_ref());

    let first = iter.next()?;
    let mut sum_a = first.0.clone();
    let mut sum_b = first.1.clone();
    let mut count = 1usize;

    for (arr_a, arr_b) in iter {
        debug_assert_eq!(
            sum_a.dim(),
            arr_a.dim(),
            "All first-array components must share shape"
        );
        debug_assert_eq!(
            sum_b.dim(),
            arr_b.dim(),
            "All second-array components must share shape"
        );

        Zip::from(&mut sum_a).and(arr_a).for_each(|s, &v| *s += v);
        Zip::from(&mut sum_b).and(arr_b).for_each(|s, &v| *s += v);

        count += 1;
    }

    let factor = count as f64;
    sum_a /= factor;
    sum_b /= factor;

    Some((sum_a, sum_b))
}

fn set_masked_entries_to_value(matrix: &mut Array2<f64>, mask: &Array2<bool>, fill_value: f64) {
    Zip::from(matrix).and(mask).for_each(|value, &is_valid| {
        if !is_valid {
            *value = fill_value;
        }
    });
}

pub fn build_extreme_gc_support_mask(
    shape: (usize, usize),
    extreme_bins_per_side: usize,
) -> Array2<bool> {
    let (num_length_bins, num_gc_bins) = shape;
    let bins_to_mask = extreme_bins_per_side.min(num_gc_bins);
    let column_is_supported: Vec<bool> = (0..num_gc_bins)
        .map(|col_idx| {
            if bins_to_mask == 0 {
                true
            } else {
                let mask_left = col_idx < bins_to_mask;
                let mask_right = col_idx >= num_gc_bins.saturating_sub(bins_to_mask);
                !(mask_left || mask_right)
            }
        })
        .collect();
    Array2::from_shape_fn((num_length_bins, num_gc_bins), |(_, col_idx)| {
        column_is_supported[col_idx]
    })
}

#[derive(Debug, Clone)]
pub struct BinnedAxis {
    pub index_to_bin: FxHashMap<usize, usize>,
    pub bin_to_indices: FxHashMap<usize, Vec<usize>>,
    pub num_bins: usize,
}

pub enum CollapseAggregation {
    Sum,
    Mean,
}

pub fn bin_greedily_by_mass<S>(
    counts: &ArrayBase<S, Ix2>,
    axis: usize,
    min_mass_pct: f64,
) -> Result<BinnedAxis>
where
    S: Data<Elem = f64>,
{
    ensure!(axis < 2, "axis must be 0 or 1");
    ensure!(
        (0.0..=100.0).contains(&min_mass_pct),
        "min_mass_pct must be within 0..=100"
    );

    // Sum along the other axis to get per-index mass
    let masses = match axis {
        0 => counts.sum_axis(Axis(1)),
        _ => counts.sum_axis(Axis(0)),
    };

    let total_mass: f64 = masses.iter().sum();
    if total_mass == 0.0 {
        return Ok(BinnedAxis {
            index_to_bin: FxHashMap::default(),
            bin_to_indices: FxHashMap::default(),
            num_bins: 0,
        });
    }

    let min_mass = total_mass * (min_mass_pct / 100.0);
    let mut bins: Vec<Vec<usize>> = Vec::new();
    let mut running_mass = 0.0;
    let mut current_bin_indices: Vec<usize> = Vec::new();

    for (idx, &mass) in masses.iter().enumerate() {
        running_mass += mass;
        current_bin_indices.push(idx);

        if running_mass >= min_mass {
            bins.push(current_bin_indices.clone());
            current_bin_indices.clear();
            running_mass = 0.0;
        }
    }

    if !current_bin_indices.is_empty() {
        if bins.is_empty() {
            bins.push(current_bin_indices);
        } else {
            bins.last_mut().unwrap().extend(current_bin_indices);
        }
    }

    let mut index_to_bin = FxHashMap::default();
    let mut bin_to_indices = FxHashMap::default();

    for (bin_idx, indices) in bins.iter().enumerate() {
        bin_to_indices.insert(bin_idx, indices.clone());
        for &idx in indices {
            index_to_bin.insert(idx, bin_idx);
        }
    }

    let num_bins = bin_to_indices.len();

    Ok(BinnedAxis {
        index_to_bin,
        bin_to_indices,
        num_bins,
    })
}

pub fn collapse_counts_by_bins<S>(
    counts: &ArrayBase<S, Ix2>,
    axis: usize,
    bins: &BinnedAxis,
    agg: CollapseAggregation,
    mass_counts: Option<ArrayView2<'_, f64>>,
) -> Result<Array2<f64>>
where
    S: Data<Elem = f64>,
{
    ensure!(axis < 2, "axis must be 0 or 1");
    if let Some(mass) = mass_counts.as_ref() {
        ensure!(
            mass.dim() == counts.dim(),
            "mass_counts must have same shape as counts"
        );
        if matches!(agg, CollapseAggregation::Sum) {
            bail!("mass_counts provided for Sum aggregation; weighted sums are unsupported");
        }
    }

    let (n_rows, n_cols) = counts.dim();
    match axis {
        0 => {
            let weights = mass_counts.as_ref().map(|m| m.sum_axis(Axis(1)));
            let mut out = Array2::<f64>::zeros((bins.num_bins, n_cols));
            for bin_idx in 0..bins.num_bins {
                if let Some(indices) = bins.bin_to_indices.get(&bin_idx) {
                    let mut denom = 0.0;
                    let mut count = 0usize;
                    for &row_idx in indices {
                        let source = counts.row(row_idx);
                        let mut dest = out.row_mut(bin_idx);
                        match agg {
                            CollapseAggregation::Sum => {
                                dest += &source;
                            }
                            CollapseAggregation::Mean => {
                                if let Some(ref weights_vec) = weights {
                                    let weight = weights_vec[row_idx];
                                    denom += weight;
                                    dest.scaled_add(weight, &source);
                                } else {
                                    dest += &source;
                                    count += 1;
                                }
                            }
                        }
                    }
                    if matches!(agg, CollapseAggregation::Mean) {
                        let mut dest = out.row_mut(bin_idx);
                        if weights.is_some() {
                            if denom > 0.0 {
                                dest /= denom;
                            } else if !indices.is_empty() {
                                dest /= indices.len() as f64;
                            }
                        } else if count > 0 {
                            dest /= count as f64;
                        }
                    }
                }
            }
            Ok(out)
        }
        _ => {
            let weights = mass_counts.as_ref().map(|m| m.sum_axis(Axis(0)));
            let mut out = Array2::<f64>::zeros((n_rows, bins.num_bins));
            for bin_idx in 0..bins.num_bins {
                if let Some(indices) = bins.bin_to_indices.get(&bin_idx) {
                    let mut denom = 0.0;
                    let mut count = 0usize;
                    for &col_idx in indices {
                        let source = counts.column(col_idx);
                        let mut dest = out.column_mut(bin_idx);
                        match agg {
                            CollapseAggregation::Sum => {
                                dest += &source;
                            }
                            CollapseAggregation::Mean => {
                                if let Some(ref weights_vec) = weights {
                                    let weight = weights_vec[col_idx];
                                    denom += weight;
                                    dest.scaled_add(weight, &source);
                                } else {
                                    dest += &source;
                                    count += 1;
                                }
                            }
                        }
                    }
                    if matches!(agg, CollapseAggregation::Mean) {
                        let mut dest = out.column_mut(bin_idx);
                        if let Some(_) = weights {
                            if denom > 0.0 {
                                dest /= denom;
                            } else if !indices.is_empty() {
                                dest /= indices.len() as f64;
                            }
                        } else if count > 0 {
                            dest /= count as f64;
                        }
                    }
                }
            }
            Ok(out)
        }
    }
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

pub struct GCCorrectionPackage {
    pub version: u32,
    pub end_offset: u32,
    pub length_edges: Vec<u32>,
    pub gc_edges: Vec<u32>,
    pub correction_matrix: Array2<f64>,
}

impl GCCorrectionPackage {
    pub fn from_components(
        version: u32,
        length_bins: &BinnedAxis,
        gc_bins: &BinnedAxis,
        correction_matrix: Array2<f64>,
        opt: &GCConfig,
    ) -> Result<Self> {
        let length_edges = compute_bin_edges(
            length_bins,
            opt.fragment_lengths.min_fragment_length as u32,
            opt.fragment_lengths.max_fragment_length as u32,
        )?;
        let gc_edges = compute_bin_edges(gc_bins, 0, 100)?;
        Ok(Self {
            version,
            end_offset: opt.end_offset as u32,
            length_edges,
            gc_edges,
            correction_matrix,
        })
    }

    pub fn write_npz<P: AsRef<std::path::Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut npz = NpzWriter::new(file);
        npz.add_array("correction_matrix", &self.correction_matrix)?;
        npz.add_array("length_edges", &Array1::from(self.length_edges.clone()))?;
        npz.add_array("gc_edges", &Array1::from(self.gc_edges.clone()))?;
        npz.add_array("version", &Array1::from(vec![self.version]))?;
        npz.add_array("end_offset", &Array1::from(vec![self.end_offset]))?;
        npz.finish()?;
        Ok(())
    }
}

fn compute_bin_edges(bins: &BinnedAxis, start_value: u32, max_value: u32) -> Result<Vec<u32>> {
    ensure!(
        bins.num_bins > 0,
        "Bin definition must contain at least one bin"
    );
    let mut edges = Vec::with_capacity(bins.num_bins + 1);
    for bin_idx in 0..bins.num_bins {
        let indices = bins
            .bin_to_indices
            .get(&bin_idx)
            .context("Missing indices for bin")?;
        let min_idx = indices
            .iter()
            .min()
            .copied()
            .context("Bin indices cannot be empty")?;
        edges.push(start_value + min_idx as u32);
    }
    edges.push(max_value);
    Ok(edges)
}
