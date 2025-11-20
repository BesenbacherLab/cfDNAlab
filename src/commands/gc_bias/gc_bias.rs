use crate::{
    commands::{
        cli_common::*,
        counters::GCCounters,
        gc_bias::{
            config::{GCConfig, WindowWeightingSchemes},
            counting::{GCCounts, build_gc_prefixes, get_gc_fraction_in_window},
            load_reference_bias::{ReferenceGcData, load_reference_gc_data},
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
use anyhow::{Context, Result, bail, ensure};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array2, ArrayBase, Axis, Data, Ix2, Zip};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{fs::create_dir_all, sync::Arc, time::Instant};

const CORRECTION_CLAMP_RANGE: (f64, f64) = (0.1, 10.0);

pub fn run(opt: &GCConfig) -> Result<()> {
    let start_time = Instant::now();
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, &opt.ioc.bam.as_path())?;
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    println!("Start: Loading reference GC bias");
    let ReferenceGcData {
        window_spec: window_opt,
        windows_map,
        window_indices_by_chr,
        counts: reference_counts,
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

    let results: Vec<(Vec<GCCounts>, GCCounters)> = chromosomes
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
    for (counts_by_bin, counter) in results {
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
    let avg_counts_tuple_opt = if let Some(ws_by_chr) = window_indices_by_chr {
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
                let out = process_window(gc_counts, &ref_counts_view, &opt, avg_window_size)?;
                pb.inc(1);
                Ok(out)
            })
            .collect::<Result<_>>()?; // short-circuits on the first Err

        pb.finish_with_message("| Finished processing");

        mean_of_arrays(&avg_counts_tuples)
    } else {
        // Global window
        let gc_counts = &all_bins[0];
        let ref_counts_view = reference_counts.index_axis(Axis(0), 0);
        process_window(gc_counts, &ref_counts_view, &opt, avg_window_size)?
    };

    let (avg_gc_counts, avg_norm_ref_counts) =
        avg_counts_tuple_opt.expect("avg count matrices should have been produced");

    // TODO: Smoothe the counts (use a pseudo count of 1.0)
    let smoothed_gc_counts = avg_gc_counts.clone();

    // TODO: Bin lengths and GCs by mass
    let binned_gc_counts = smoothed_gc_counts.clone();
    let binned_ref_counts = avg_norm_ref_counts.clone();

    // TODO: Outlier detection -> Smoothing

    let norm_gc_counts = mean_scale_per_length_array(&binned_gc_counts, 0.);
    let norm_ref_counts = mean_scale_per_length_array(&binned_ref_counts, 0.);

    // Calculate correction matrix
    let correction_matrix = norm_gc_counts / norm_ref_counts;

    // Sanity clamp of corrections
    let correction_matrix =
        correction_matrix.clamp(CORRECTION_CLAMP_RANGE.0, CORRECTION_CLAMP_RANGE.1);

    // TODO: Save info needed to apply correction

    // // Write final counts to output_dir
    write_npy(
        &opt.ioc.output_dir.join("gc_bias_correction.npy"),
        &correction_matrix,
    )
    .context("Write final fail")?;

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
            opt.gc_min_pct as usize,
            opt.gc_max_pct as usize,
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
        let min_len = opt
            .fragment_lengths
            .min_fragment_length
            // Minimum length to live up to --min-acgt-count and --end-offset
            .max(opt.min_fragment_acgt_count as u32 + 2 * opt.end_offset as u32);
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
    let min_fragment_acgt_count_u32 = opt.min_fragment_acgt_count as u32;
    let min_acgt_fraction = opt.min_fragment_acgt_pct as f32 / 100f32;

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
        let gc = get_gc_fraction_in_window(
            &gc_prefixes,
            gc_window_start as usize,
            gc_window_end as usize,
            min_acgt_fraction,
            min_fragment_acgt_count_u32,
        );

        // Unpack GC fraction (or continue)
        let gc = match gc {
            Some(v) => v,
            None => continue,
        };

        ensure!(gc.is_finite(), "GC non-finite: {}", gc);
        ensure!(
            (0.0..=1.0).contains(&gc),
            "GC fraction out of [0,1]: {}",
            gc
        );

        // Get GC in [0,100] as usize
        let gc_bin = (gc * 100.0).round() as usize;

        if gc_bin < opt.gc_min_pct as usize || gc_bin > opt.gc_max_pct as usize {
            continue;
        }

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
/// 1) Reference counts are normalized (mean-scaled) per fragment length.
///
/// 2) The two arrays are weighted depending on the weighting scheme:
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
fn process_window<S>(
    gc_counts: &GCCounts,
    ref_counts: &ArrayBase<S, Ix2>,
    opt: &GCConfig,
    avg_window_size: Option<f64>,
) -> Result<Option<(Array2<f64>, Array2<f64>)>>
where
    S: Data<Elem = u64>,
{
    // Check window has enough valid positions
    if gc_counts.pct_acgt() < opt.min_window_acgt_pct as f64 {
        return Ok(None);
    }

    // Get counts as 2d array with shape: (n_lengths, n_gc_bins)
    let counts_mat = gc_counts.to_array2();
    let (n_rows, n_cols) = counts_mat.dim();

    let total_count: f64 = counts_mat.sum();
    let mean_count: f64 = total_count / (n_rows * n_cols) as f64;

    if total_count == 0.0 {
        return Ok(None);
    }

    // Normalize the reference bias
    // Cast to f64 once
    let ref_counts_f = ref_counts.mapv(|v| v as f64);
    // Row-wise mean-scaling to avoid
    let ref_counts_norm = mean_scale_per_length_array(&ref_counts_f, 1.0);

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

fn mean_scale_per_length_array<S>(x: &ArrayBase<S, Ix2>, pseudo_count: f64) -> Array2<f64>
where
    S: Data<Elem = f64>,
{
    let (_, n_cols) = x.dim();
    let n_cols_f = n_cols as f64;

    // Row means: shape (n_rows,)
    let row_means = x.mean_axis(Axis(1)).expect("mean_axis on non-empty axis");

    // Denominator per row: row_mean + pseudo_count * num_cols
    // shape (n_rows,)
    let denom = &row_means + pseudo_count * n_cols_f;

    // Broadcast to (n_rows, 1) so it can divide the full matrix
    let denom_2d = denom.insert_axis(Axis(1));

    // Numerator: value + pseudo_count (elementwise)
    let num = x + pseudo_count;

    // Elementwise division with broadcasting
    num / &denom_2d
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
