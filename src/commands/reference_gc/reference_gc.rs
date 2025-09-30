use crate::{
    commands::gc_bias::counting::{
        GCCounts, build_gc_prefixes, count_reference_gc_and_length_by_window, stack_gc_counts,
    },
    commands::{cli_common::*, reference_gc::config::RefGCConfig},
    shared::{
        bed::load_windows_from_bed,
        blacklist::{apply_blacklist_mask_to_seq, compute_blacklist_overlap, load_blacklists},
        reference::{read_seq, twobit_contig_lengths},
        sampling::sample_starts_per_chrom,
    },
};
use anyhow::{Context, Result};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use std::{
    fs::{File, create_dir_all},
    io::{BufWriter, Write},
    sync::Arc,
    time::Instant,
};

pub fn run(opt: RefGCConfig) -> Result<()> {
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
    let blacklist_map = if let Some(beds) = &opt.blacklist {
        println!("Start: Loading blacklists");
        load_blacklists(beds, 1, &chromosomes)?
    } else {
        FxHashMap::default()
    };

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            println!("Start: Loading window coordinates");
            Some(load_windows_from_bed(bed, &chromosomes, None)?)
        }
        _ => None,
    };

    let starts_per_chrom = {
        let mut rng1 = rand::rng();
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

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let results: Vec<(Vec<GCCounts>, Option<Vec<(String, u64, u64, u64, f64)>>)> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _)> {
            let out = process_chrom(
                &chr,
                &opt,
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
    for (counts_by_bin, bin_vec) in results {
        all_bins.extend(counts_by_bin);
        if !matches!(window_opt, WindowSpec::Global) {
            bin_info.extend(bin_vec.unwrap());
        }
    }

    // Convert to single `GCCounts` for global
    // Keep wrapped in vector to simplify writer
    let mut all_bins = if matches!(window_opt, WindowSpec::Global) {
        vec![GCCounts::collapse(&all_bins)?]
    } else {
        all_bins
    };

    // Sort by original index (when given a bed file)
    if matches!(window_opt, WindowSpec::Bed(_)) {
        println!("Start: Reordering counts by original window index in BED file");

        // Zip into a single Vec to allow sorting together
        let mut paired: Vec<_> = bin_info.into_iter().zip(all_bins.into_iter()).collect(); // (BinInfo, DecodedCounts)

        // Sort primarily by original window index
        paired.sort_unstable_by_key(|(info, _)| info.3);

        // Unzip back out if you need separate Vecs again
        (bin_info, all_bins) = paired.into_iter().unzip();
    }

    // Write final counts to output_dir
    write_npy(
        &opt.output_dir.join("all_ref_gc_counts.npy"),
        &stack_gc_counts(&all_bins),
    )
    .context("Write final fail")?;

    // Write bins BED file
    if !matches!(window_opt, WindowSpec::Global) {
        println!("Start: Writing window coordinates to disk");
        let mut bed_writer = BufWriter::new(
            File::create(&opt.output_dir.join("bins.bed")).context("Create bed fail")?,
        );
        for (chr, start, end, _, overlap_perc) in &bin_info {
            writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
                .context("Write bed line fail")?;
        }
    }

    // Print execution time
    let elapsed = start_time.elapsed();
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
) -> anyhow::Result<(Vec<GCCounts>, Option<Vec<(String, u64, u64, u64, f64)>>)> {
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
            opt.fragment_lengths.max_fragment_length as usize
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
        opt.min_acgt_pct as f32 / 100f32,
        opt.min_acgt_count as u32,
    );

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

    Ok((counts_by_bin, bin_info))
}
