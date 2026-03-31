use crate::{
    commands::{
        cli_common::{ensure_output_dir, load_blacklist_map, resolve_chromosomes_and_contigs},
        counters::CoverageWeightsCounters,
        coverage_weights::config::CoverageWeightsConfig,
        coverage_weights::striding::{
            StrideBin, fill_triangular_overlap, normalize_avg_overlap_by_global_mean,
        },
    },
    shared::{
        bam::create_chromosome_reader,
        coverage::Coverage,
        fragment::minimal_fragment::Fragment,
        fragment_iterators::fragments_from_bam,
        interval::Interval,
        io::dot_join,
        progress::ProgressFactory,
        read::{default_include_read_paired_end, default_include_read_unpaired},
        thread_pool::init_global_pool,
    },
};
use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{
    fs::File,
    io::{BufWriter, Write},
    sync::Arc,
    time::Instant,
};

/// Calculates weights for genomic smoothing using large bins and a stride.
///
/// Technical details:
/// - Resolves chromosomes, prepares output directories, and loads optional blacklists before
///   scanning each chromosome in parallel.
/// - Converts fragments into coverage profiles, smooths them with a triangular kernel, and writes
///   the resulting statistics to a TSV file.
/// - Tracks iterator counters so the printed summary reflects accepted fragments, blacklist hits,
///   and other bookkeeping numbers.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `coverage-weights` command.
///
/// Returns
/// -------
/// - `Ok(())`:
///     Scaling factors were written successfully.
///
/// Errors
/// ------
/// - Returns an error if the BAM cannot be read, blacklist files are invalid, or the output file
///   cannot be created.
pub fn run(opt: &CoverageWeightsConfig) -> Result<()> {
    let start_time = Instant::now();
    if opt.unpaired.reads_are_fragments && opt.require_proper_pair {
        bail!("--require-proper-pair cannot be used with --reads-are-fragments");
    }
    let (chromosomes, _contigs) =
        resolve_chromosomes_and_contigs(&opt.chromosomes, opt.ioc.bam.as_path())?;
    opt.check_bin_sizes()?;
    let progress = ProgressFactory::new();
    let pb = Arc::new(progress.default_bar(chromosomes.len() as u64));

    // Create output directory
    ensure_output_dir(&opt.ioc.output_dir)?;

    // Load blacklist intervals if provided
    if opt.blacklist.is_some() {
        println!("Start: Loading blacklists");
    }
    let blacklist_map = load_blacklist_map(
        opt.blacklist.as_ref(),
        opt.blacklist_min_size,
        0,
        &chromosomes,
    )?;

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads)?;

    // Prepare output containers
    let mut bins_by_chr =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());
    let mut global_counter = CoverageWeightsCounters::default();

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let results: Vec<(String, Vec<StrideBin>, CoverageWeightsCounters)> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _, _)> {
            let out = process_chrom(
                chr,
                opt,
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    // Collect results (in chromosome order) back into the global vectors
    for (chr, stride_bins, counter) in results {
        bins_by_chr.insert(chr, stride_bins);
        global_counter += counter;
    }

    // Normalize by global mean and invert to scaling factors (keeping 0s intact)
    let global_avg_overlap_coverage =
        normalize_avg_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

    println!(
        "Calculated the global average overlapping position-coverage: {}",
        global_avg_overlap_coverage
    );

    // Write stride-bin coordinates and scaling factors as TSV to output_dir

    println!("Start: Writing stride-bin coordinates and scaling factors to disk");
    let file_name = dot_join(&[opt.output_prefix.as_str(), "scaling_factors.tsv"]);
    let mut tsv_writer = BufWriter::new(
        File::create(opt.ioc.output_dir.join(file_name)).context("Creating tsv failed")?,
    );
    writeln!(
        tsv_writer,
        "chromosome\tstart\tend\tavg_pos_cov\tavg_overlapping_pos_cov\tscaling_factor"
    )
    .context("Write TSV header fail")?;
    for chr in chromosomes {
        let bins = bins_by_chr
            .get(&chr)
            .with_context(|| format!("missing bins for chromosome: {}", chr))?;

        for bin in bins.iter() {
            writeln!(
                tsv_writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                chr,
                bin.start(),
                bin.end(),
                bin.avg_coverage,
                bin.avg_overlap_coverage,
                bin.scaling_factor
            )
            .context("Write TSV row fail")?;
        }
    }

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
    // if opt.gc.bin_by_gc {
    //     println!("GC-excluded reads: {}", global_counter.base.gc_excl);
    // }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.base.counted_fragments
    );
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn process_chrom(
    chr: &str,
    opt: &CoverageWeightsConfig,
    blacklist_intervals: &[Interval<u64>],
) -> anyhow::Result<(String, Vec<StrideBin>, CoverageWeightsCounters)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = CoverageWeightsCounters::default();

    let mut bins: Vec<StrideBin> = {
        let mut v = Vec::new();
        let mut pos = 0u32;
        while pos < chrom_len as u32 {
            v.push(StrideBin {
                interval: Interval::new(pos, pos.saturating_add(opt.stride).min(chrom_len as u32))?,
                avg_coverage: 0.0,
                avg_overlap_coverage: 0.0,
                scaling_factor: 0.0,
            });
            pos = pos.saturating_add(opt.stride);
        }
        v
    };

    reader
        .fetch((tid, 0, chrom_len))
        .context(format!("fetch {}", chr))?;

    // Initialize coverage counter
    let mut cp = Coverage::new(chrom_len as u32);

    // Function for filtering fragments after pairing
    // Note: We need to own the data in the fn (not just pass `opt` that could disappear)
    let fragment_filter = {
        let lengths = opt.fragment_lengths.clone();
        move |f: &Fragment| lengths.contains(f.len())
    };

    // Create fragment iterator
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

    // Iterate fragments and add fragment to coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;

        counter.base.counted_fragments += 1;

        // Add to coverage prefix
        cp.add_fragment(fragment)?;
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    // Add blacklist
    if !blacklist_intervals.is_empty() {
        cp.set_blacklist_mask(blacklist_intervals)?;
    }

    // Get ready to extract average coverage per stride-bin
    cp.finalize_coverage(true);
    cp.build_indexes(true)?;

    // Decide once whether to exclude blacklisted bases
    let exclude_blacklisted = cp.blacklist_mask().is_some();

    for bin in bins.iter_mut() {
        // Calculate total coverage in bin
        bin.avg_coverage = cp.avg_coverage(bin.interval, exclude_blacklisted)?;
    }

    // Update the avg_overlap_coverage per bin
    fill_triangular_overlap(&mut bins, opt.bin_size, opt.stride);

    Ok((chr.to_string(), bins, counter))
}
