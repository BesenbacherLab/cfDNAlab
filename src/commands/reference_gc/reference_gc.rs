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
        reference_gc::config::RefGCConfig,
    },
    shared::{
        bed::load_windows_from_bed,
        blacklist::{apply_blacklist_mask_to_seq, compute_blacklist_overlap},
        io::create_text_writer,
        reference::{read_seq, twobit_contig_lengths},
        sampling::sample_starts_per_chrom,
    },
};
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use ndarray::Array2;
use ndarray_npy::write_npy;
use rand::{SeedableRng, rngs::StdRng};
use rayon::prelude::*;
use std::{fs::create_dir_all, io::Write, sync::Arc, time::Instant};

pub fn run(opt: &RefGCConfig) -> Result<()> {
    let start_time = Instant::now();
    let chromosomes = opt.chromosomes.resolve_chromosomes(None)?;
    let window_opt = opt.windows.resolve_windows();
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
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
        let all_bins = if matches!(window_opt, WindowSpec::Global) {
            vec![GCCounts::collapse(&all_bins)?]
        } else {
            all_bins
        };

        // Convert counts to Array2
        println!("Start: Converting counts to arrays");
        pb.reset();
        pb.set_length(all_bins.len() as u64);
        pb.set_position(0);
        let mut all_bin_arrays: Vec<Array2<f64>> = all_bins
            .par_iter()
            .map(|gc_counts| -> Result<_> {
                let mut out = gc_counts.to_array2();
                apply_gc_percent_width_correction(&mut out, &gc_percent_widths)?;
                pb.inc(1);
                Ok(out)
            })
            .collect::<Result<_>>()?; // short-circuits on the first Err
        pb.finish_with_message("| Finished array conversion");

        // Create mask of supported count bins BEFORE interpolation
        // Elements that are zero in all windows are considered impossible
        // for the combination of GC and fragment lengths when
        // all positions are valid ACGT bases
        // NOTE: These are found empirically here but could be calculated theoretically?
        let mut outlier_support_mask = create_support_mask_threshold_per_mb(
            all_bin_arrays.as_slice(),
            total_covered_acgt_positions,
            2.0,
        )
        .expect("support mask should be created");

        if !opt.skip_interpolation {
            println!("Start: Interpolating missing counts");
            pb.reset();
            pb.set_length(all_bins.len() as u64);
            pb.set_position(0);

            for array in all_bin_arrays.iter_mut() {
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
        (all_bin_arrays, outlier_support_mask)
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

    write_npy(
        &opt.output_dir.join("ref_support_mask.outliers.npy"),
        &outlier_support_mask,
    )
    .context("Writing outlier support mask failed")?;

    write_npy(
        &opt.output_dir.join("ref_support_mask.unobservables.npy"),
        &unobservable_support_mask,
    )
    .context("Writing unobservables support mask failed")?;

    write_npy(
        &opt.output_dir.join("ref_gc_percent_widths.npy"),
        &*gc_percent_widths,
    )
    .context("Writing GC percent widths failed")?;

    // Write final counts to output_dir
    write_npy(
        &opt.output_dir.join("ref_gc_counts.npy"),
        &stack_gc_counts(&all_bin_arrays),
    )
    .context("Write final fail")?;

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

fn process_chrom(
    chr: &str,
    opt: &RefGCConfig,
    windows: Option<&[(u64, u64, u64)]>,
    window_opt: &WindowSpec,
    // gc_bins: usize,
    blacklist_intervals: &[(u64, u64)],
    start_positions: &[usize],
) -> anyhow::Result<(
    Vec<GCCounts>,
    Option<Vec<(String, u64, u64, u64, f64)>>,
    u64,
)> {
    let mut seq_bytes = read_seq(&opt.ref_genome.ref_2bit, chr)?;
    apply_blacklist_mask_to_seq(&mut seq_bytes, &blacklist_intervals, 0);
    let chrom_len = seq_bytes.len() as u64;

    let gc_prefixes = build_gc_prefixes(&seq_bytes);

    // Delete seq_bytes from memory
    drop(seq_bytes);

    // Calculate window coordinates for all windowing options
    let windows: Vec<(u64, u64, u64)> = match window_opt {
        WindowSpec::Bed(_) => windows.unwrap().to_owned(),
        WindowSpec::Size(sz) => {
            let num_windows = ((chrom_len + sz - 1) / sz) as u64;
            (0..num_windows)
                .map(|s| ((s * sz) as u64, (sz + s * sz) as u64, s as u64))
                .collect()
        }
        WindowSpec::Global => vec![(0, chrom_len as u64, 0u64)],
    };
    let num_bins = windows.len();

    // Initialize count arrays
    let mut counts_by_bin = vec![
        GCCounts::new(
            0usize,
            100usize,
            opt.fragment_lengths.min_fragment_length as usize,
            opt.fragment_lengths.max_fragment_length as usize,
            (0, 0) // Not used in this command
        );
        num_bins
    ];

    count_reference_gc_and_length_by_window(
        &mut counts_by_bin,
        &gc_prefixes,
        (
            opt.fragment_lengths.min_fragment_length as u64,
            opt.fragment_lengths.max_fragment_length as u64 + 1, // make exclusive
        ),
        &windows,
        start_positions,
        chrom_len,
        // We only count those where all bases are proper, so non-supported combinations
        // of GC and fragment lengths are truly 0 (enabling clean interpolation)
        1.0,
        1u32,
    );

    // Calculate total number of ACGT positions covered
    // NOTE: This does not deduplicate positions, so overlaps
    // count per occurence - but this is okay for normalization
    let total_acgt_in_chrom = {
        let mut total_acgt = 0u64;
        for (start, end, _) in windows.iter() {
            let clamped_end = *end.min(&chrom_len);
            let clamped_start = *start.min(&chrom_len);
            if clamped_end <= clamped_start {
                continue;
            }
            let acgt =
                gc_prefixes.acgt[clamped_end as usize] - gc_prefixes.acgt[clamped_start as usize];
            total_acgt += acgt as u64;
        }
        total_acgt
    };

    // TODO: To use this for filtering, we need to know the N-positions as well (we need percentage usable positions)

    let bin_info = if let Some(size) = opt.windows.by_size {
        // Build bin information for chromosome
        // chrom,start,end,total_count
        let mut bl_ptr = 0;
        let mut bin_info = Vec::with_capacity(num_bins);
        for b in 0..num_bins {
            let start = b as u64 * size;
            let end = (start + size).min(chrom_len);
            let overlap_perc =
                compute_blacklist_overlap(blacklist_intervals, start, end, 0u64, &mut bl_ptr);
            // Note: b (index) is a placeholder that is removed later
            bin_info.push((chr.to_string(), start, end, b as u64, overlap_perc));
        }
        Some(bin_info)
    } else if opt.windows.by_bed.is_some() {
        // build bin_info from the exact BED windows
        let mut bl_ptr = 0;
        let mut bin_info = Vec::with_capacity(num_bins);
        for (_b, (wstart, wend, original_win_idx)) in windows.iter().cloned().enumerate() {
            let overlap_perc =
                compute_blacklist_overlap(blacklist_intervals, wstart, wend, 0u64, &mut bl_ptr);
            bin_info.push((
                chr.to_string(),
                wstart,
                wend,
                original_win_idx as u64,
                overlap_perc,
            ));
        }
        Some(bin_info)
    } else {
        None
    };

    Ok((counts_by_bin, bin_info, total_acgt_in_chrom))
}
