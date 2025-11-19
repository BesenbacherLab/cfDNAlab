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
        100 - opt.min_window_positions_pct,
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
    let global_correction_opt = if let Some(ws_by_chr) = window_indices_by_chr {
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
        let window_corrections: Vec<Option<Array2<f64>>> = window_tuples
            .par_iter()
            .map(|(count_idx, ref_idx)| -> Result<_> {
                let gc_counts = &all_bins[*count_idx];
                let ref_counts_view = reference_counts.index_axis(Axis(0), *ref_idx);
                let out = process_window(gc_counts, &ref_counts_view, &opt, avg_window_size)?;
                pb.inc(1);
                Ok(out)
            })
            .collect::<Result<_>>()?; // short-circuits on the first Err

        pb.finish_with_message("| Finished processing");

        mean_of_arrays(&window_corrections)
    } else {
        // Global window
        let gc_counts = &all_bins[0];
        let ref_counts_view = reference_counts.index_axis(Axis(0), 0);
        process_window(gc_counts, &ref_counts_view, &opt, avg_window_size)?
    };

    let global_correction =
        global_correction_opt.expect("a correction matrix should have been produced");

    // TODO: Outlier detection -> Smoothing

    // TODO: Save info needed to apply correction

    // // Write final counts to output_dir
    write_npy(
        &opt.ioc.output_dir.join("gc_bias_correction.npy"),
        &global_correction,
    )
    .context("Write final fail")?;

    // TODO: Create actual correction factor here!

    // println!("Start: Writing mismatch rates to disk");
    // write_mismatch_rate_tsvs(&prepared_counts, &mismatch_rates_by_k, &opt.ioc.output_dir)?;

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
            .max(opt.min_acgt_count as u32 + 2 * opt.end_offset as u32);
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
    let min_acgt_count_u32 = opt.min_acgt_count as u32;
    let min_acgt_fraction = opt.min_acgt_pct as f32 / 100f32;

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
            min_acgt_count_u32,
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

fn process_window<S>(
    gc_counts: &GCCounts,
    ref_bias: &ArrayBase<S, Ix2>,
    opt: &GCConfig,
    avg_window_size: Option<f64>,
) -> Result<Option<Array2<f64>>>
where
    S: Data<Elem = u64>,
{
    let pseudo_count = 1.0;
    let min_bin_count = opt.min_bin_count as f64;

    let counts_mat = gc_counts.to_array2(); // shape: (n_lengths, n_gc_bins)
    let (n_rows, n_cols) = counts_mat.dim();

    let mut too_sparse_mask = Array2::<bool>::default((n_rows, n_cols));

    let mut total: f64 = 0.0;
    let mut row_totals = vec![0.0; n_rows];

    for (row_idx, row) in counts_mat.outer_iter().enumerate() {
        let mut row_sum = 0.0;
        for (col_idx, &val) in row.iter().enumerate() {
            let below = val < min_bin_count;
            too_sparse_mask[(row_idx, col_idx)] = below;
            row_sum += val;
            total += val;
        }
        row_totals[row_idx] = row_sum;
    }

    if total < opt.min_window_count as f64 {
        return Ok(None);
    }

    // TODO: Should we exclude masked counts from mean here as well?
    // Normalize each fragment length (row) by mean-scaling
    let mut mean_scaled = counts_mat.clone();
    for (row_idx, mut row) in mean_scaled.outer_iter_mut().enumerate() {
        let row_mean = row_totals[row_idx] / n_cols as f64;
        let denom = row_mean + pseudo_count * n_cols as f64;
        row.map_inplace(|v| *v = (*v + pseudo_count) / denom);
    }

    // Normalize the reference bias

    // Cast to f64 once
    let ref_bias_f = ref_bias.mapv(|v| v as f64);
    let ref_bias_norm = mean_scale_per_length_array(&ref_bias_f, pseudo_count);

    // Calculate the GC bias in this window!
    let bias = mean_scaled / ref_bias_norm;

    // Row-normalize using only non-masked elements, then set masked to 1.0
    let masked_scaled_bias = scale_bias_with_mask(&bias, &too_sparse_mask);

    let weighted_bias = match opt.window_weighting {
        WindowWeightingSchemes::Equal => masked_scaled_bias,
        WindowWeightingSchemes::Coverage => {
            let (n_rows, n_cols) = masked_scaled_bias.dim();
            let mean_coverage = total / (n_rows * n_cols) as f64;
            masked_scaled_bias * mean_coverage
        }
        WindowWeightingSchemes::ValidPositions => {
            let num_acgt = gc_counts.num_acgt_out_of.0;
            if num_acgt == 0 {
                // No positions observed
                return Ok(None);
            }
            let avg_window_span = avg_window_size.expect("valid-positions needs window spans");
            masked_scaled_bias * (num_acgt as f64 / avg_window_span)
        }
    };

    // TODO: Per window: Normalizations, Reference div, Scaling, weight-scheme scaling

    Ok(Some(weighted_bias))
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

/// Scale by row means calculated from non-masked elements
/// and then set masked elements to 1.0
fn scale_bias_with_mask<S>(bias: &ArrayBase<S, Ix2>, mask: &Array2<bool>) -> Array2<f64>
where
    S: Data<Elem = f64>,
{
    let (n_rows, n_cols) = bias.dim();
    debug_assert_eq!(mask.dim(), (n_rows, n_cols));

    let mut scaled = bias.to_owned();

    for row_idx in 0..n_rows {
        let bias_row = bias.row(row_idx);
        let mask_row = mask.row(row_idx);

        let mut sum = 0.0;
        let mut count = 0usize;

        for col_idx in 0..n_cols {
            if !mask_row[col_idx] {
                sum += bias_row[col_idx];
                count += 1;
            }
        }

        let row_mean = if count > 0 { sum / count as f64 } else { 1.0 };

        let mut scaled_row = scaled.row_mut(row_idx);
        for col_idx in 0..n_cols {
            scaled_row[col_idx] /= row_mean;
            if mask_row[col_idx] {
                scaled_row[col_idx] = 1.0;
            }
        }
    }

    scaled
}

/// Average the arrays and skip `None`s
fn mean_of_arrays(arrays: &[Option<Array2<f64>>]) -> Option<Array2<f64>> {
    // Find first present array to initialize accumulator
    let mut iter = arrays.iter().filter_map(|opt| opt.as_ref());

    let first = iter.next()?;
    let mut sum = first.clone();
    let mut count: usize = 1;

    // Accumulate remaining Some(Array2)
    for arr in iter {
        debug_assert_eq!(sum.dim(), arr.dim(), "All arrays must have same shape");

        Zip::from(&mut sum).and(arr).for_each(|s, &v| {
            *s += v;
        });

        count += 1;
    }

    // Average
    let factor = count as f64;
    sum /= factor;

    Some(sum)
}
