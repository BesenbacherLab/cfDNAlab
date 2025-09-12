use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{
    collections::HashMap,
    fs::{File, create_dir_all},
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Arc,
    time::Instant,
};

use crate::{
    cli_common::{ChromosomeArgs, FragmentLengthArgs, IOCArgs},
    counters::NormalizeGenomeCounters,
    utils::{
        bam::create_chromosome_reader,
        blacklist::{BlacklistStrategy, load_blacklists},
        coverage::CoveragePrefix,
        fragment::{MinimalReadInfo, collect_fragment},
        normalize_genome::{
            StrideBin, fill_triangular_overlap, normalize_avg_overlap_by_global_mean,
        },
    },
};

/// Extract fragment coverage in large genomic bins ("megabins") with a rolling
/// window and calculate normalizing scaling factors for smoothing the
/// genome.
///
/// Smoothing is performed as a triangular moving average, where we calculate
/// a weighted average of coverages from all bins overlapping a stride.
///
/// ## Triangular weighting scheme visualization
///
/// Assuming a bin-size of 6 and stride size of 2 (normally defaults to 5Mb and 0.5Mb respectively):
///
/// --------------------------
/// Stride bins (fixed along genome, each have an average coverage value):
///
/// `[A] [B] [C] [D] [E] [F] [G] ...`
///
/// Overlapping megabins (`MB*`) (each covers 3 stride-bins).
///
/// `W_D`, the number of overlapping megabins, is the (unnormalized) weight of each stride-bin
/// in the weighted-average coverage for stride-bin `D`:
///
/// ```
/// 
/// MB1: [A][B][C]
/// 
/// MB2:    [B][C][D]
/// 
/// MB3:       [C][D][E]
/// 
/// MB4:          [D][E][F]
/// 
/// MB5:             [E][F][G]
/// 
/// W_D: [0][1][2][3][2][1][0]
/// 
/// ```
///
/// At chromosome edges, the weights are truncated (e.g., `W_D: [2][3][2][1][0]`).
///
/// The weights are normalized by their sum (after potential truncation at edges).
#[cfg_attr(feature = "cli", derive(clap::Args))]
pub struct NormalizeGenomeConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    ioc: IOCArgs,

    /// Size (bp) of large genomic bins to calculate coverage in [integer]
    ///
    /// Larger values lead to a more smooth coverage across the genome.
    ///
    /// NOTE: The normalizing scaling factors are calculated per stride-sized overlap
    /// of these bins. Technically, we only count the coverage per stride-sized bin
    /// and then calculate the overlap with a triangular weighting scheme.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "5000000", value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub bin_size: u32,

    /// Size (bp) of stride [integer]
    ///
    /// NOTE: `--bin_size` must be divisible by `stride`. I.e., `bin_size % stride` == 0`.
    ///
    /// A normalizing scaling factor is calculated per stride as the weighted average coverage of the overlapping large-scale bins.
    ///
    /// Smaller values lead to a higher precision in the downstream normalization
    /// but also requires saving a larger BED file in the end (one line per stride-bin)
    /// and take longer to compute.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "500000", value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub stride: u32,

    #[cfg_attr(feature = "cli", clap(flatten))]
    chromosomes: ChromosomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    fragment_lengths: FragmentLengthArgs,

    /// Minimum mapping quality to include [integer]
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads [flag]
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions [path]
    #[cfg_attr(
        feature = "cli", clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum size of blacklist intervals to load (bp) [integer]
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-min-size",
            default_value = "1",
            help_heading = "Filtering"
        )
    )]
    pub blacklist_min_size: u64,

    /// The fragment positions that should overlap blacklisted regions for it to be excluded [string]
    ///
    /// Possible values:
    ///     "any", "all", "midpoint", or "proportion=<threshold>" [string]
    ///
    /// Example of proportion: `--blacklist-strategy proportion=0.2` (no space around `=`)
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            alias = "bl-strategy",
            default_value = "any",
            ignore_case = true,
            help_heading = "Filtering"
        )
    )]
    pub blacklist_strategy: BlacklistStrategy,
}

impl NormalizeGenomeConfig {
    pub fn check_bin_sizes(&self) -> anyhow::Result<()> {
        let stride = self.stride.clone();
        let bin_size = self.bin_size.clone();

        if stride > bin_size {
            bail!(
                "stride ({}) cannot be higher than bin_size ({})",
                stride,
                bin_size,
            );
        }
        if bin_size % stride != 0 {
            bail!(
                "bin_size ({}) must be divisible by stride ({})",
                bin_size,
                stride
            );
        }

        Ok(())
    }
}

/// Whether to include the read or continue
fn include_read(rec: &Record, opt: &NormalizeGenomeConfig) -> bool {
    !(rec.is_unmapped()
        || rec.is_mate_unmapped()
        || rec.tid() != rec.mtid()
        || rec.is_secondary()
        || rec.is_supplementary()
        || rec.is_duplicate()
        || rec.is_quality_check_failed()
        || (opt.require_proper_pair && !rec.is_proper_pair())
        || rec.mapq() < opt.min_mapq) as bool
}

pub fn run(opt: NormalizeGenomeConfig) -> Result<()> {
    let start_time = Instant::now();
    let chromosomes = opt
        .chromosomes
        .resolve_chromosomes(Some(&opt.ioc.bam.as_path()))?;
    opt.check_bin_sizes()?;
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    // Create output directory
    create_dir_all(&opt.ioc.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = if let Some(beds) = &opt.blacklist {
        println!("Start: Loading blacklists");
        load_blacklists(beds, opt.blacklist_min_size, &chromosomes)?
    } else {
        HashMap::new()
    };

    // Cconfigure global thread‐pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.ioc.n_threads as usize)
        .build_global()
        .context("building Rayon thread pool")?;

    // Prepare per-bin counts and metadata
    let mut bins_by_chr =
        FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());
    let mut global_counter = NormalizeGenomeCounters::default();

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let results: Vec<(String, Vec<StrideBin>, NormalizeGenomeCounters)> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _, _)> {
            let out = process_chrom(
                &chr,
                &opt,
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

    let global_avg_overlap_coverage = normalize_avg_overlap_by_global_mean(&mut bins_by_chr, true)?;

    println!(
        "Calculated the global average overlapping position-coverage: {}",
        global_avg_overlap_coverage
    );

    // Write window coordinates as BED file to output_dir
    // Write bins BED file

    println!("Start: Writing stride-bin coordinates and scaling factors to disk");
    let mut bed_writer = BufWriter::new(
        File::create(&opt.ioc.output_dir.join("coverage_scaling_factors.tsv"))
            .context("Creating tsv failed")?,
    );
    writeln!(
        bed_writer,
        "{}\t{}\t{}\t{}\t{}\t{}",
        "chromosome", "start", "end", "avg_pos_cov", "avg_overlapping_pos_cov", "scaling_factor"
    )
    .context("Write bed line fail")?;
    for chr in chromosomes {
        let bins = bins_by_chr
            .get(&chr)
            .with_context(|| format!("missing bins for chromosome: {}", chr))?;

        for bin in bins.into_iter() {
            writeln!(
                bed_writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                chr,
                bin.start,
                bin.end,
                bin.avg_coverage,
                bin.avg_overlap_coverage,
                bin.scaling_factor
            )
            .context("Write tsv line fail")?;
        }
    }

    // Print summary statistics and execution time
    let elapsed = start_time.elapsed();
    println!("  Total reads: {}", global_counter.total_reads);
    println!(
        "  Initially accepted reads: {} ({:.2}%, forward: {}, reverse: {})",
        global_counter.accepted_forward + global_counter.accepted_reverse,
        (global_counter.accepted_forward + global_counter.accepted_reverse) as f64
            / global_counter.total_reads as f64
            * 100.0,
        global_counter.accepted_forward,
        global_counter.accepted_reverse
    );
    println!(
        "Out-of-length-range-excluded fragments: {}",
        global_counter.illegal_length_fragments
    );
    // if opt.gc.bin_by_gc {
    //     println!("GC-excluded reads: {}", global_counter.gc_excl);
    // }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.counted_fragments
    );
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn process_chrom(
    chr: &str,
    opt: &NormalizeGenomeConfig,
    blacklist_intervals: &[(u64, u64)],
) -> anyhow::Result<(String, Vec<StrideBin>, NormalizeGenomeCounters)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    // Initialize counters (default -> 0s)
    let mut counter = NormalizeGenomeCounters::default();

    let mut bins: Vec<StrideBin> = {
        let mut v = Vec::new();
        let mut pos = 0u32;
        while pos < chrom_len as u32 {
            v.push(StrideBin {
                start: pos,
                end: pos.saturating_add(opt.stride).min(chrom_len as u32),
                avg_coverage: 0.0,
                avg_overlap_coverage: 0.0,
                scaling_factor: 0.0,
            });
            pos = pos.saturating_add(opt.stride);
        }
        v
    };

    // Stash for keeping reads until mates arrive
    let mut stash = FxHashMap::<Vec<u8>, MinimalReadInfo>::default();

    reader
        .fetch((tid, 0, chrom_len))
        .context(format!("fetch {}", chr))?;

    // Initialize coverage counter
    let mut cp = CoveragePrefix::initialize_coverage_prefix(chrom_len as u32);

    // Loop over records and count
    for res in reader.records() {
        let rec = res.context("reading bam record")?;
        counter.total_reads += 1;

        if rec.tid() != tid as i32 || !include_read(&rec, opt) {
            continue;
        }

        match rec.is_reverse() {
            true => counter.accepted_reverse += 1,
            false => counter.accepted_forward += 1,
        }

        if let Some(mate) = stash.remove(rec.qname()) {
            // Combine reads to a fragment
            let fragment = if let Some(f) = collect_fragment(&MinimalReadInfo::from(&rec), &mate) {
                f
            } else {
                continue;
            };

            counter.collected_fragments += 1;

            // Check length is within allowed range
            let fragment_length = fragment.len();
            if fragment_length < opt.fragment_lengths.min_fragment_length
                || fragment_length > opt.fragment_lengths.max_fragment_length
            {
                counter.illegal_length_fragments += 1;
                continue;
            }

            counter.counted_fragments += 1;

            // Add to coverage prefix
            cp.add_fragment_to_prefix(fragment)?;
        } else {
            // Stash read if new qname
            stash.insert(rec.qname().to_vec(), MinimalReadInfo::from(&rec));
        }
    }

    // Add blacklist
    if !blacklist_intervals.is_empty() {
        cp.initialize_blacklist_prefix();
        cp.add_blacklist_many_to_prefix(blacklist_intervals)?;
        cp.finalize_blacklist_prefix();
    }

    // Get ready to extract average coverage per stride-bin
    cp.finalize_coverage();
    cp.build_query_index()?;

    // Decide once whether to exclude blacklisted bases
    let exclude_blacklisted = cp.blacklist_mask().is_some();

    for bin in bins.iter_mut() {
        // Calculate total coverage in bin
        bin.avg_coverage = cp.avg_coverage(bin.start, bin.end, exclude_blacklisted)?;
    }

    // Update the avg_overlap_coverage per bin
    fill_triangular_overlap(&mut bins, opt.bin_size, opt.stride);

    Ok((chr.to_string(), bins, counter))
}
