use crate::{
    commands::{
        cli_common::*,
        gc_bias::{
            counting::{
                GCCounts, apply_gc_percent_width_correction, build_gc_prefixes,
                count_reference_gc_and_length_by_window, gc_percent_widths, stack_gc_counts,
            },
            interpolation::fill_unsupported_bins_with_polynomial,
            support_masking::{
                build_theoretical_support_mask, create_support_mask_threshold_per_mb,
            },
        },
        ref_gc_counts::config::RefGCCountsConfig,
    },
    shared::{
        bed::load_windows_from_bed,
        blacklist::{apply_blacklist_mask_to_seq, compute_blacklist_overlap},
        interval::{IndexedInterval, Interval},
        io::create_text_writer,
        reference::{read_seq, twobit_contig_lengths},
        sampling::sample_starts_per_chrom,
    },
};
use anyhow::{Context, Result, ensure};
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array1, Array2, Array3};
use ndarray_npy::NpzWriter;
use rand::{SeedableRng, rngs::StdRng};
use rayon::prelude::*;
use std::{fs::create_dir_all, io::Write, sync::Arc, time::Instant};

pub fn run(opt: &RefGCCountsConfig) -> Result<()> {
    let start_time = Instant::now();
    let chromosomes = opt.chromosomes.resolve_chromosomes(None)?;
    let window_opt = opt.windows.resolve_windows();
    opt.check_smoothing_settings()?;
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    let min_effective_len = opt
        .fragment_lengths
        .min_fragment_length
        .saturating_sub(2 * u32::from(opt.end_offset));
    ensure!(
        min_effective_len >= 10,
        "Requires at least 10 bases for GC calculation. --min-fragment-length ({}) - 2x --end-offset ({}) is < 10. Please adjust --min-fragment-length.",
        opt.fragment_lengths.min_fragment_length,
        opt.end_offset
    );

    // Create output directory
    create_dir_all(&opt.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = load_blacklist_map(opt.blacklist.as_ref(), 1, 0, &chromosomes)?;

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            println!("Start: Loading window coordinates");
            Some(load_windows_from_bed(
                bed,
                Some(chromosomes.as_slice()),
                None,
                None,
            )?)
        }
        _ => None,
    };

    // Precompute GC% bin widths (gc_count -> percent) per fragment length
    let gc_percent_widths = Arc::new(gc_percent_widths(
        opt.fragment_lengths.min_fragment_length as usize,
        opt.fragment_lengths.max_fragment_length as usize,
        opt.end_offset as usize,
    ));

    let starts_per_chrom = {
        let mut rng1 = if let Some(seed) = opt.seed {
            StdRng::seed_from_u64(seed)
        } else {
            let mut thread_rng = rand::rng();
            StdRng::from_rng(&mut thread_rng)
        };
        sample_starts_per_chrom(
            &mut rng1,
            &twobit_contig_lengths(opt.ref_genome.ref_2bit.clone(), &chromosomes)?,
            opt.n_positions,
            opt.fragment_lengths.max_fragment_length as usize,
        )?
    };

    // Configure global thread‐pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.n_threads as usize)
        .build_global()
        .context("building Rayon thread pool")?;

    // Prepare per-bin counts and metadata
    let mut all_bins = Vec::new();
    let mut bin_info = Vec::new();
    let mut total_covered_acgt_positions = 0u64;

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let results: Vec<(
        Vec<GCCounts>,
        Option<Vec<(String, u64, u64, u64, f64)>>,
        u64,
    )> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _, _)> {
            let out = process_chrom(
                &chr,
                opt,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                &window_opt,
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                &starts_per_chrom.get(chr).unwrap_or(&vec![]),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    println!("Start: Processing counts");

    // Collect results (in chromosome order) back into the global vectors
    for (counts_by_bin, bin_vec, num_chrom_acgt) in results {
        all_bins.extend(counts_by_bin);
        total_covered_acgt_positions += num_chrom_acgt;
        if !matches!(window_opt, WindowSpec::Global) {
            bin_info.extend(bin_vec.unwrap());
        }
    }

    // Combine counts, convert to Array2 and interpolate zero-counts
    // Nesting to clean up GCCounts version when no longer needed
    let (mut all_bin_arrays, outlier_support_mask) = {
        // Convert to single `GCCounts` for global
        // Keep wrapped in vector to simplify writer
        let mut all_bins = if matches!(window_opt, WindowSpec::Global) {
            vec![GCCounts::collapse(&all_bins)?]
        } else {
            all_bins
        };

        // Convert counts to Array2
        println!("Start: Converting counts to arrays");
        pb.reset();
        pb.set_length(all_bins.len() as u64);
        pb.set_position(0);
        let mut all_bin_grids: Vec<Array2<f64>> = all_bins
            .par_iter_mut()
            .map(|gc_counts: &mut GCCounts| -> Result<_> {
                if !opt.skip_smoothing {
                    gc_counts
                        .smooth_length_rows_in_place(opt.smoothing_sigma, opt.smoothing_radius);
                }
                let mut out = gc_counts.to_gc_percent_grid(0, 100)?;
                apply_gc_percent_width_correction(&mut out, &gc_percent_widths)?;
                pb.inc(1);
                Ok(out)
            })
            .collect::<Result<_>>()?; // short-circuits on the first Err
        pb.finish_with_message("| Finished array conversion");

        // Create mask of supported count bins BEFORE interpolation
        // Elements that are seen less than 3 times per 1Mb are considered unsupported.
        // These include the theoretically unobservable combinations of fragment lengths and GC percentage bins.
        let mut outlier_support_mask = create_support_mask_threshold_per_mb(
            all_bin_grids.as_slice(),
            total_covered_acgt_positions,
            3.0,
        )
        .expect("support mask should be created");

        if !opt.skip_interpolation {
            println!("Start: Interpolating missing counts");
            pb.reset();
            pb.set_length(all_bins.len() as u64);
            pb.set_position(0);

            for array in all_bin_grids.iter_mut() {
                debug_assert_eq!(
                    array.dim(),
                    outlier_support_mask.dim(),
                    "Support mask and histograms must match shape"
                );
                for (row_idx, mut length_row) in array.outer_iter_mut().enumerate() {
                    // Rows are contiguous so we can safely borrow a mutable slice for interpolation
                    let row_slice = length_row
                        .as_slice_mut()
                        .expect("GC histogram rows should be contiguous");
                    let mut mask_row = outlier_support_mask.row_mut(row_idx);
                    let mask_slice = mask_row
                        .as_slice_mut()
                        .expect("Support mask rows should be contiguous");
                    fill_unsupported_bins_with_polynomial(
                        row_slice, mask_slice, 2, 3, 3,
                        // Do not update mask, we need the raw mask
                        false,
                    )?;
                }
                pb.inc(1);
            }
        }

        pb.finish_with_message("| Finished interpolating missing counts");
        (all_bin_grids, outlier_support_mask)
    };

    // Sort by original index (when given a bed file)
    if matches!(window_opt, WindowSpec::Bed(_)) {
        println!("Start: Reordering counts by original window index in BED file");

        // Zip into a single Vec to allow sorting together
        let mut paired: Vec<_> = bin_info
            .into_iter()
            .zip(all_bin_arrays.into_iter())
            .collect(); // (bin info, final counts)

        // Sort primarily by original window index
        paired.sort_unstable_by_key(|(info, _)| info.3);

        // Unzip back out if you need separate Vecs again
        (bin_info, all_bin_arrays) = paired.into_iter().unzip();
    }

    // Write support mask to output_dir
    let unobservable_support_mask = build_theoretical_support_mask(
        opt.fragment_lengths.min_fragment_length as usize,
        opt.fragment_lengths.max_fragment_length as usize,
        0,
        all_bin_arrays
            .first()
            .map(|arr| arr.dim().1 - 1)
            .expect("at least one GC histogram should exist"),
        opt.end_offset as usize,
    );

    debug_assert_eq!(
        outlier_support_mask.dim(),
        unobservable_support_mask.dim(),
        "Outlier support mask shape {:?} must match unobservable support mask shape {:?}",
        outlier_support_mask.dim(),
        unobservable_support_mask.dim()
    );

    if let Some(first_hist) = all_bin_arrays.first() {
        debug_assert_eq!(
            unobservable_support_mask.dim(),
            first_hist.dim(),
            "Support mask shape {:?} must match histogram shape {:?}",
            unobservable_support_mask.dim(),
            first_hist.dim()
        );
    }

    // Write reference GC package (counts, masks, widths, metadata)
    write_reference_gc_package(
        &opt.output_dir.join("ref_gc_package.npz"),
        &stack_gc_counts(&all_bin_arrays),
        &unobservable_support_mask,
        &outlier_support_mask,
        &*gc_percent_widths,
        opt.fragment_lengths.min_fragment_length as usize,
        opt.fragment_lengths.max_fragment_length as usize,
        opt.end_offset,
        opt.skip_interpolation,
        opt.smoothing_radius,
        opt.smoothing_sigma,
        opt.skip_smoothing,
    )
    .context("Writing reference GC package failed")?;

    // Write bins BED file
    if !matches!(window_opt, WindowSpec::Global) {
        println!("Start: Writing window coordinates to disk");
        let bins_path = opt.output_dir.join("ref_gc_bins.bed");
        let mut bed_writer = create_text_writer(&bins_path).context("Create bed fail")?;
        for (chr, start, end, _, overlap_perc) in &bin_info {
            writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
                .context("Write bed line fail")?;
        }
        bed_writer.finish().context("Finalize bins.bed writer")?;
    }

    // Print execution time
    let elapsed = start_time.elapsed();
    println!(
        "Windows covered {} total ACGT bases (may have duplicate positions if windows overlap)",
        total_covered_acgt_positions
    );
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn write_reference_gc_package(
    path: &std::path::Path,
    counts: &Array3<f64>,
    support_unobservables: &Array2<bool>,
    support_outliers: &Array2<bool>,
    gc_percent_widths: &Array2<u16>,
    length_min: usize,
    length_max: usize,
    end_offset: u8,
    skip_interpolation: bool,
    smoothing_radius: u8,
    smoothing_sigma: f64,
    skip_smoothing: bool,
) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut npz = NpzWriter::new(file);
    npz.add_array("counts", counts)?;
    npz.add_array("support_mask_unobservables", support_unobservables)?;
    npz.add_array("support_mask_outliers", support_outliers)?;
    npz.add_array("gc_percent_widths", gc_percent_widths)?;
    npz.add_array(
        "length_range",
        &Array1::from(vec![length_min as u32, length_max as u32]),
    )?;
    npz.add_array("end_offset", &Array1::from(vec![end_offset as u32]))?;
    npz.add_array(
        "skip_interpolation",
        &Array1::from(vec![skip_interpolation]),
    )?;
    npz.add_array(
        "smoothing_radius",
        &Array1::from(vec![smoothing_radius as u32]),
    )?;
    npz.add_array("smoothing_sigma", &Array1::from(vec![smoothing_sigma]))?;
    npz.add_array("skip_smoothing", &Array1::from(vec![skip_smoothing]))?;
    npz.finish()?;
    Ok(())
}

fn process_chrom(
    chr: &str,
    opt: &RefGCCountsConfig,
    windows: Option<&[IndexedInterval<u64>]>,
    window_opt: &WindowSpec,
    // gc_bins: usize,
    blacklist_intervals: &[Interval<u64>],
    start_positions: &[usize],
) -> anyhow::Result<(
    Vec<GCCounts>,
    Option<Vec<(String, u64, u64, u64, f64)>>,
    u64,
)> {
    let mut seq_bytes = read_seq(&opt.ref_genome.ref_2bit, chr)?;
    apply_blacklist_mask_to_seq(&mut seq_bytes, blacklist_intervals, 0);
    let chrom_len = seq_bytes.len() as u64;

    let gc_prefixes = build_gc_prefixes(&seq_bytes);

    // Delete seq_bytes from memory
    drop(seq_bytes);

    // Calculate window coordinates for all windowing options
    let windows: Vec<IndexedInterval<u64>> = match window_opt {
        WindowSpec::Bed(_) => windows.unwrap().to_owned(),
        WindowSpec::Size(sz) => {
            let num_windows = ((chrom_len + sz - 1) / sz) as u64;
            (0..num_windows)
                .map(|window_index| {
                    IndexedInterval::new(
                        window_index * sz,
                        (sz + window_index * sz).min(chrom_len),
                        window_index,
                    )
                })
                .collect::<std::result::Result<Vec<_>, _>>()?
        }
        WindowSpec::Global => vec![IndexedInterval::new(0, chrom_len, 0u64)?],
    };
    let num_bins = windows.len();

    // Initialize count arrays
    let mut counts_by_bin = vec![
        GCCounts::new(
            opt.fragment_lengths.min_fragment_length as usize,
            opt.fragment_lengths.max_fragment_length as usize,
            opt.end_offset as usize,
            (0, 0) // Not used in this command
        )?;
        num_bins
    ];

    count_reference_gc_and_length_by_window(
        &mut counts_by_bin,
        &gc_prefixes,
        (
            opt.fragment_lengths.min_fragment_length as u64,
            opt.fragment_lengths.max_fragment_length as u64 + 1, // make exclusive
        ),
        windows.as_slice(),
        start_positions,
        chrom_len,
        // We only count those where all bases are proper, so non-supported combinations
        // of GC and fragment lengths are truly 0 (enabling clean interpolation)
        1.0,
        MIN_ACGT_BASES_FOR_GC_FRACTION,
        opt.end_offset as usize,
    );

    // Calculate total number of ACGT positions covered
    // NOTE: This does not deduplicate positions, so overlaps
    // count per occurence - but this is okay for normalization
    let total_acgt_in_chrom = {
        let mut total_acgt = 0u64;
        for window in &windows {
            let clamped_end = window.end().min(chrom_len);
            let clamped_start = window.start().min(chrom_len);
            if clamped_end <= clamped_start {
                continue;
            }
            let acgt =
                gc_prefixes.acgt[clamped_end as usize] - gc_prefixes.acgt[clamped_start as usize];
            total_acgt += acgt as u64;
        }
        total_acgt
    };

    let bin_info = if let Some(size) = opt.windows.by_size {
        // Build bin information for chromosome
        // chrom,start,end,total_count
        let mut bl_ptr = 0;
        let mut bin_info = Vec::with_capacity(num_bins);
        for window in &windows {
            let overlap_perc = compute_blacklist_overlap(
                blacklist_intervals,
                Interval::new(start, end)?,
                0u64,
                &mut bl_ptr,
            );
            bin_info.push((
                chr.to_string(),
                window.start(),
                window.end(),
                // Note: index is a placeholder that is removed later
                window.idx(),
                overlap_perc,
            ));
        }
        Some(bin_info)
    } else if opt.windows.by_bed.is_some() {
        // build bin_info from the exact BED windows
        let mut bl_ptr = 0;
        let mut bin_info = Vec::with_capacity(num_bins);
        for window in &windows {
            let overlap_perc = compute_blacklist_overlap(
                blacklist_intervals,
                Interval::new(start, end)?,
                0u64,
                &mut bl_ptr,
            );
            bin_info.push((
                chr.to_string(),
                window.start(),
                window.end(),
                window.idx(),
                overlap_perc,
            ));
        }
        Some(bin_info)
    } else {
        None
    };

    Ok((counts_by_bin, bin_info, total_acgt_in_chrom))
}
