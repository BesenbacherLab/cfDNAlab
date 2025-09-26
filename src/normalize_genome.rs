use anyhow::{Context, Result, bail};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{
    fs::File,
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
        blacklist::BlacklistStrategy,
        command::{ensure_output_dir, load_blacklist_map, resolve_chromosomes_and_contigs},
        coverage::coverage_prefix::Coverage,
        fragment::minimal_fragment::Fragment,
        fragment_iterator::fragments_from_bam,
        normalize_genome::{
            StrideBin, fill_triangular_overlap, normalize_avg_overlap_by_global_mean,
        },
        read::default_include_read,
        thread_pool::init_global_pool,
    },
};

// TODO: Improve docstring - hard to understand

/// Extract fragment coverage in large genomic bins ("megabins") with a rolling
/// window and calculate normalizing scaling factors for smoothing the genome.
///
/// Outputs scaling factors per stride to allow other methods to apply the normalization (by weighting fragment counts).
///
/// The scaling factors are *inverted*, so normalization becomes multiplication.
/// Zero-valued coverages lead to zero-valued scaling factors. Non-zero factors have `mean == 1.0`.
///
/// ## Coverage
///
/// The full fragment span `[forward.pos, reverse.end)` is counted without consideration of deletions and gaps.
/// This is fine for genome-scale normalization that reduces relative changes in coverage across the genome.
///
/// ## Smoothing
///
/// Smoothing is performed as a triangular moving average, calculating
/// a weighted average of coverages from all bins overlapping a stride.  
///
/// ### Example
///
/// Assuming a bin-size of 6 and stride size of 2 (normally defaults to 5Mb and 0.5Mb respectively).
///
/// **Stride bins** (fixed along genome, each with an average coverage):
///
/// `[A] [B] [C] [D] [E] [F] [G] ...`
///
/// **Overlapping megabins** (`MB*`) (each covers 3 stride-bins). **`W_D`**, the number of overlapping megabins,
/// is the (unnormalized) weight of each stride-bin in the weighted-average coverage for stride-bin `D`:
///
/// ```text
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
///
/// ## Always-on exclusion criteria
///
/// The following criteria always exclude a read:
///
/// The read or mate read is unmapped.
/// The read is mapped to a different `tid` than the mate.
/// The read is secondary, supplementary or duplicate.
/// The read failed quality check.
/// The paired reads are not inwardly directed (we require: `start(forward) <= start(reverse)`).
#[cfg_attr(feature = "cli", derive(clap::Args))]
#[derive(Clone)]
pub struct NormalizeGenomeConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    ioc: IOCArgs,

    /// Size (bp) of large genomic bins to calculate coverage in [integer]
    ///
    /// Larger values lead to a more smooth coverage across the genome.
    ///
    /// **NOTE**: The normalizing scaling factors are calculated per stride-sized overlap
    /// of these bins. Technically, we only count the coverage per stride-sized bin
    /// and then calculate the overlap with a triangular weighting scheme.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "5000000", value_parser = clap::value_parser!(u32).range(1..), help_heading="Filtering"))]
    pub bin_size: u32,

    /// Size (bp) of stride [integer]
    ///
    /// **NOTE**: `--bin_size` must be divisible by `stride`. I.e., `bin_size % stride` == 0`.
    ///
    /// A normalizing scaling factor is calculated per stride as the (inverse) weighted average coverage of the overlapping large-scale bins.
    ///
    /// Smaller values lead to a higher precision in the downstream normalization
    /// but also require saving a larger BED file in the end (one line per stride-bin)
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
    pub fn new(ioc: IOCArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            bin_size: 5_000_000,
            stride: 500_000,
            chromosomes,
            fragment_lengths: FragmentLengthArgs {
                min_fragment_length: 20,
                max_fragment_length: 1000,
            },
            min_mapq: 30,
            require_proper_pair: false,
            blacklist: None,
            blacklist_min_size: 1,
            blacklist_strategy: BlacklistStrategy::default(),
        }
    }

    pub fn set_bin_size(&mut self, bin_size: u32) {
        self.bin_size = bin_size;
    }

    pub fn set_stride(&mut self, stride: u32) {
        self.stride = stride;
    }

    pub fn fragment_lengths_mut(&mut self) -> &mut FragmentLengthArgs {
        &mut self.fragment_lengths
    }

    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    pub fn set_require_proper_pair(&mut self, require: bool) {
        self.require_proper_pair = require;
    }

    pub fn set_blacklist(&mut self, blacklist: Option<Vec<PathBuf>>) {
        self.blacklist = blacklist;
    }

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

/// Execute the genome-normalisation pipeline and emit stride-level scaling factors.
///
/// Implementation details:
/// - Resolves chromosomes, prepares output directories, and loads optional blacklists before
///   scanning each chromosome in parallel.
/// - Converts fragments into coverage profiles, smooths them with a triangular kernel, and writes
///   the resulting statistics to a TSV file.
/// - Tracks iterator counters so the summary in the terminal reflects acceptance, blacklist hits,
///   and fragment counts.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `normalize-genome` command.
///
/// Returns:
/// - `Ok(())` when the scaling-factor TSV is written successfully.
///
/// Errors:
/// - Returns an error if the BAM cannot be read, blacklist files are invalid, or the output file
///   cannot be created.
pub fn run(opt: NormalizeGenomeConfig) -> Result<()> {
    let start_time = Instant::now();
    let (chromosomes, _contigs) = resolve_chromosomes_and_contigs(&opt.chromosomes, &opt.ioc)?;
    opt.check_bin_sizes()?;
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    // Create output directory
    ensure_output_dir(&opt.ioc.output_dir)?;

    // Load blacklist intervals if provided
    if opt.blacklist.is_some() {
        println!("Start: Loading blacklists");
    }
    let blacklist_map =
        load_blacklist_map(opt.blacklist.as_ref(), opt.blacklist_min_size, &chromosomes)?;

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads as usize)?;

    // Prepare output containers
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

    // Normalize by global mean and invert to scaling factors (keeping 0s intact)
    let global_avg_overlap_coverage =
        normalize_avg_overlap_by_global_mean(&mut bins_by_chr, true, true)?;

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
        bin.avg_coverage = cp.avg_coverage(bin.start, bin.end, exclude_blacklisted)?;
    }

    // Update the avg_overlap_coverage per bin
    fill_triangular_overlap(&mut bins, opt.bin_size, opt.stride);

    Ok((chr.to_string(), bins, counter))
}
