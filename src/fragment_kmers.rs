use crate::{
    cli_common::{
        ChromosomeArgs, FragmentLengthArgs, IOCArgs, Ref2BitRequiredArgs, ScaleGenomeArgs,
        WindowSpec, WindowsArgs,
    },
    counters::FragmentKmersCounters,
    utils::{
        bam::create_chromosome_reader,
        bed::load_windows_from_bed,
        blacklist::{
            BlacklistStrategy, apply_blacklist_mask_to_seq, apply_mask::BLACKLIST_BYTE,
            compute_blacklist_overlap, is_blacklisted,
        },
        command::{
            ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        coverage::scale_genome::apply_scaling_to_coverage_in_place,
        fragment::segment_kmer_fragment::FragmentWithKmerSegments,
        fragment_iterator::fragments_with_kmer_segments_from_bam,
        indel_mode::IndelMode,
        kmers::{
            kmer_codec::{
                Kmer, KmerCodes, KmerSpec, build_kmer_specs, build_left_aligned_codes_per_k,
            },
            process_counts::{
                DecodedCounts, merge_decoded_counts, prepare_decoded_counts,
                split_and_decode_counts,
            },
            write::write_decoded_counts_matrix,
        },
        overlaps::find_overlapping_windows,
        read::default_include_read,
        reference::read_seq,
        thread_pool::init_global_pool,
    },
};
use anyhow::{Context, Result};
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

/// Count kmers within the fragments in a BAM-file.
///
/// Whereas the `cfdna ends` tool extracts end-motifs, this tool extracts all kmers
/// in a sliding window across the fragment.
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
pub struct FragmentKmersConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    ioc: IOCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    pub ref_genome: Ref2BitRequiredArgs,

    /// Prefix for output files (e.g., a sample name) `[string]`
    ///
    /// E.g., specify to enable writing to the same output directory from multiple calls to this software.
    ///
    /// Examples produce files like:
    ///   `<prefix>.k3_counts.npy`,
    ///   `<prefix>.k3_motifs.txt`,
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            short = 'x',
            default_value = "fragment_kmers",
            help_heading = "Core"
        )
    )]
    pub output_prefix: String,

    /// List of K-mer sizes [integer].
    ///
    /// When counting for many kmer-sizes (>8), consider splitting
    /// into multiple runs to reduce memory consumption at a time.
    ///
    /// Example: `--kmer-sizes 3 5 11`
    #[cfg_attr(
        feature = "cli",
        clap(short = 'k', long, num_args = 1.., value_parser = clap::value_parser!(u8).range(1..28), required=true, help_heading="Core"))]
    pub kmer_sizes: Vec<u8>,

    /// Number of bases to exclude from each end of fragments `[integer]`
    ///
    /// This allows not counting end-motifs, to focus only on the center kmers.
    /// For pure end-motif counting, use `cfdna ends` instead.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "0", value_parser = clap::value_parser!(u32).range(0..), help_heading="Core"))]
    pub end_offset: u32,

    /// How to handle insertions and deletions in fragments `[string]`
    ///
    /// Deletions: Both 'D' and 'N' in the cigar string are considered deletions.
    ///
    /// Possible values:
    ///
    /// - `"ignore"`:
    ///   Ignore whether indels are present or not.
    ///   Kmers are extracted for the full/offset fragment span from the reference genome.
    ///
    /// - `"adjust"`:
    ///   Adjust the counts by excluding kmers overlapping positions with observed insertions and deletions in the
    ///   observed bases (we cannot adjust in mate-gaps).
    ///   Outside the mate-overlap, all indels and deletions are adjusted for.
    ///   **Overlap**: In the mate-overlap, both reads must agree on the position-level.
    ///   Only overlap-positions were both reads have the indel are excluded.
    ///   **NOTE**: Blacklist exclusion and calculation of scaling weights (--scaling-factors)
    ///   use the full reference span.
    ///
    /// - `"skip"`:
    ///   Skip fragments with any insertion or deletion present.
    #[cfg_attr(
        feature = "cli",
        clap(
            long,
            default_value = "ignore",
            ignore_case = true,
            help_heading = "Core"
        )
    )]
    pub indel_mode: IndelMode,

    /// Ignore inter-mate gap `[flag]`
    ///
    /// Disable counting in the gap between reads (i.e., `[forward.end, reverse.start)`)
    /// when the two reads do not overlap.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub ignore_gap: bool,

    /// Collapse each kmer with its reverse-complement. [flag]
    ///
    /// Odd-sized kmers are collapsed such that the middle base is `A` or `C`.
    /// Even-sized kmers are collapsed to the lexicographically lowest motif.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    canonical: bool,

    /// Save counts as sparse-array. [flag]
    ///
    /// For large kmer-sizes, we cannot save dense arrays with all motifs
    /// unless we have a LOT of RAM and storage space. Enable this
    /// flag to save as a COO sparse array that can be opened in
    /// python via `scipy.sparse.load_npz()`.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Core"))]
    pub save_sparse: bool,

    #[cfg_attr(feature = "cli", clap(flatten))]
    windows: WindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    chromosomes: ChromosomeArgs,

    // TODO: Add that we use the scaling weight for the first kmer-position
    // And that sf=0 for any kmer base guarantees the kmer is excluded
    #[cfg_attr(feature = "cli", clap(flatten))]
    scale_genome: ScaleGenomeArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    fragment_lengths: FragmentLengthArgs,

    /// Minimum mapping quality to include `[integer]`
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "30", value_parser = clap::value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads `[flag]`
    ///
    /// This is NOT recommended by default as it trims the tails of the length distribution.
    #[cfg_attr(feature = "cli", clap(long, help_heading = "Filtering"))]
    pub require_proper_pair: bool,

    /// Optional BED file(s) with blacklisted regions `[path]`
    ///
    /// Two levels of filtering are performed. First, all blacklisted regions are assigned
    /// the N-"base" to exclude kmers that include the positions. Then, depending on the `--blacklist-strategy`,
    /// fragments overlapping blacklisted regions with some fraction are excluded.
    #[cfg_attr(
        feature = "cli", clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading = "Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum size of blacklist intervals to load (bp) `[integer]`
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

    /// The fragment positions that should overlap blacklisted regions for it to be excluded `[string]`
    ///
    /// Possible values:
    ///     "any", "all", "midpoint", or "proportion=<threshold>"
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
    // #[cfg_attr(feature = "cli", clap(flatten))]
    // gc: GCArgs,

    // #[cfg_attr(feature = "cli", clap(flatten))]
    // two_bit: TwoBitArgs,
}

impl FragmentKmersConfig {
    pub fn new(ioc: IOCArgs, ref_genome: Ref2BitRequiredArgs, chromosomes: ChromosomeArgs) -> Self {
        Self {
            ioc,
            ref_genome,
            output_prefix: "fragment_kmers".to_string(),
            kmer_sizes: vec![3u8],
            end_offset: 0,
            indel_mode: IndelMode::Ignore,
            ignore_gap: false,
            canonical: false,
            save_sparse: false,
            windows: WindowsArgs::default(),
            chromosomes,
            scale_genome: ScaleGenomeArgs::default(),
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

    pub fn set_output_prefix(&mut self, output_prefix: String) {
        self.output_prefix = output_prefix;
    }

    pub fn set_kmer_sizes(&mut self, kmer_sizes: Vec<u8>) {
        self.kmer_sizes = kmer_sizes;
    }

    pub fn set_end_offset(&mut self, end_offset: u32) {
        self.end_offset = end_offset;
    }

    pub fn set_ignore_gap(&mut self, ignore_gap: bool) {
        self.ignore_gap = ignore_gap;
    }

    pub fn set_canonical(&mut self, canonical: bool) {
        self.canonical = canonical;
    }

    pub fn set_save_sparse(&mut self, save_sparse: bool) {
        self.save_sparse = save_sparse;
    }

    pub fn set_indel_mode(&mut self, indel_mode: IndelMode) {
        self.indel_mode = indel_mode;
    }

    pub fn set_windows(&mut self, windows: WindowsArgs) {
        self.windows = windows;
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
}

/// Execute the fragment kmers counting pipeline end-to-end.
///
/// Implementation details:
/// - Resolves chromosomes, prepares optional windows/blacklists/scaling data, and then processes
///   each chromosome in parallel tiles using Rayon.
/// - Streams fragments through per-window accumulators, enumerating the requested k-mers inside
///   every counted window and writing dense (or optional sparse) count matrices plus motif lists.
/// - Applies fragment-length, blacklist, indel, scaling, and strand handling policies consistently
///   across threads.
///
/// Parameters:
/// - `opt`: Fully resolved configuration for the `fragment-kmers` command.
///
/// Returns:
/// - `Ok(())` when the counts and accompanying metadata files are written successfully.
///
/// Errors:
/// - Propagates IO and parsing errors when reading inputs or writing results, aborting the run on
///   the first failure.
pub fn run(opt: FragmentKmersConfig) -> Result<()> {
    let start_time = Instant::now();
    let (chromosomes, contigs) = resolve_chromosomes_and_contigs(&opt.chromosomes, &opt.ioc)?;
    let window_opt = opt.windows.resolve_windows();
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

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            println!("Start: Loading window coordinates");
            Some(load_windows_from_bed(bed, &chromosomes, None)?)
        }
        _ => None,
    };

    let kmer_specs: FxHashMap<u8, KmerSpec> = build_kmer_specs(&opt.kmer_sizes)?;

    // Load genomic scaling factors
    if opt.scale_genome.scaling_factors.is_some() {
        println!("Start: Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> =
        load_scaling_map(&opt.scale_genome, &chromosomes, &contigs)?;

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads as usize)?;

    // Prepare per-bin counts and metadata
    let mut all_bins = Vec::new();
    let mut bin_info = Vec::new();
    let mut global_counter = FragmentKmersCounters::default();

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let results: Vec<(
        Vec<FxHashMap<Kmer, f64>>,
        Option<Vec<(String, u64, u64, u64, f64)>>,
        FragmentKmersCounters,
    )> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _, _)> {
            let out = process_chrom(
                &chr,
                &opt,
                &kmer_specs,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                &window_opt,
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
                scaling_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    // Collect results (in chromosome order) back into the global vectors
    for (counts_by_bin, bin_vec, counter) in results {
        let counts_decoded: Vec<DecodedCounts> = counts_by_bin
            .iter()
            .map(|c| split_and_decode_counts(c, &kmer_specs))
            .collect();
        all_bins.extend(counts_decoded);
        if !matches!(window_opt, WindowSpec::Global) {
            bin_info.extend(bin_vec.unwrap());
        }
        global_counter += counter;
    }

    // Convert to single map for global
    // Keep wrapped in vector to simplify writer
    let all_bins = if matches!(window_opt, WindowSpec::Global) {
        vec![merge_decoded_counts(all_bins)]
    } else {
        all_bins
    };

    // Prepare counts to get correct motifs (collapsed, N-filtered, etc.)
    let (mut prepared_counts, motifs_by_k) =
        prepare_decoded_counts(&all_bins, opt.canonical, &kmer_specs);

    // Sort by original index (when given a bed file)
    if matches!(window_opt, WindowSpec::Bed(_)) {
        println!("Start: Reordering counts by original window index in BED file");

        // Zip into a single Vec to allow sorting together
        let mut paired: Vec<_> = bin_info
            .into_iter()
            .zip(prepared_counts.into_iter())
            .collect(); // (BinInfo, DecodedCounts)

        // Sort primarily by original window index
        paired.sort_unstable_by_key(|(info, _)| info.3);

        // Unzip back out if you need separate Vecs again
        (bin_info, prepared_counts) = paired.into_iter().unzip();
    }

    // Write final counts to output_dir
    println!("Start: Writing counts to disk");
    write_decoded_counts_matrix(
        &prepared_counts,
        &kmer_specs,
        &motifs_by_k,
        &opt.ioc.output_dir,
        &opt.output_prefix,
        opt.save_sparse,
    )?;

    // Write window coordinates as BED file to output_dir
    // Write bins BED file
    if !matches!(window_opt, WindowSpec::Global) {
        println!("Start: Writing window coordinates to disk");
        let mut bed_writer = BufWriter::new(
            File::create(&opt.ioc.output_dir.join("bins.bed")).context("Create bed fail")?,
        );
        for (chr, start, end, _, overlap_perc) in &bin_info {
            writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
                .context("Write bed line fail")?;
        }
    }

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
        "  Blacklist-excluded fragments: {}",
        global_counter.blacklisted_fragments
    );
    // if opt.gc.bin_by_gc {
    //     println!("GC-excluded reads: {}", global_counter.base.gc_excl);
    // }
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
    opt: &FragmentKmersConfig,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    windows: Option<&[(u64, u64, u64)]>,
    window_opt: &WindowSpec,
    blacklist_intervals: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
) -> anyhow::Result<(
    Vec<FxHashMap<Kmer, f64>>,
    Option<Vec<(String, u64, u64, u64, f64)>>,
    FragmentKmersCounters,
)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    let mut seq_bytes = read_seq(&opt.ref_genome.ref_2bit, chr)?;
    apply_blacklist_mask_to_seq(&mut seq_bytes, &blacklist_intervals);

    // Scaled weights to count up
    let positional_scaling_weights = if !scaling_chr.is_empty() {
        let mut scaling_weights = vec![1.0; seq_bytes.len()];
        apply_scaling_to_coverage_in_place(&mut scaling_weights, 0, scaling_chr);
        // "Blacklist" positions with scaling factors of 0, so they don't get counted
        for (base, weight) in seq_bytes.iter_mut().zip(&scaling_weights) {
            if *weight == 0.0 {
                *base = BLACKLIST_BYTE;
            }
        }
        Some(scaling_weights)
    } else {
        None
    };

    // Prepare left-aligned kmer-codes for each kmer-size
    let positional_codes_by_k: FxHashMap<u8, KmerCodes> =
        build_left_aligned_codes_per_k(&seq_bytes, kmer_specs);

    // Initialize counters (default -> 0s)
    let mut counter = FragmentKmersCounters::default();

    let num_bins = match window_opt {
        WindowSpec::Bed(_) => windows.unwrap().len(),
        WindowSpec::Size(s) => ((chrom_len + s - 1) / s) as usize,
        WindowSpec::Global => 1,
    };

    // Initialize count arrays
    let mut counts_by_bin = vec![FxHashMap::<Kmer, f64>::default(); num_bins];

    // Streaming pointers and single fetch for this chr
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = 0; // Genomic window

    // Get coordinates to fetch reads from and to
    let (fetch_from, fetch_to) = match window_opt {
        WindowSpec::Bed(_) => {
            let wn = windows.unwrap();
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
        let lengths = opt.fragment_lengths.clone();
        move |f: &FragmentWithKmerSegments| lengths.contains(f.len())
    };

    // Wrap to use opt
    let include_read_fn = {
        let opt = (*opt).clone();
        move |r: &Record| default_include_read(r, opt.require_proper_pair, opt.min_mapq)
    };

    // Create fragment iterator
    let mut iter = fragments_with_kmer_segments_from_bam(
        reader.records().map(|r| r.map_err(anyhow::Error::from)),
        include_read_fn,
        opt.indel_mode,
        !opt.ignore_gap,
        opt.end_offset,
        fragment_filter,
    )
    .with_local_counters();

    // Iterate fragments and add coverage
    for fragment_res in iter.by_ref() {
        let fragment = fragment_res.context("reading fragment")?;

        // Determine blacklist status
        let in_blacklist = is_blacklisted(
            blacklist_intervals,
            opt.blacklist_strategy.clone(),
            fragment.start.into(),
            fragment.end.into(),
            opt.fragment_lengths.max_fragment_length as u64,
            &mut bl_ptr,
        );
        if in_blacklist {
            counter.blacklisted_fragments += 1;
            continue;
        }

        // Find all overlapping count-windows
        let overlapping_windows = find_overlapping_windows(
            chrom_len,
            &mut wd_ptr,
            windows,
            opt.windows.by_size,
            (fragment.start + opt.end_offset).into(), // Should only get fragments where this is okay
            (fragment.end - opt.end_offset).into(),
            1. / (opt.fragment_lengths.max_fragment_length as f64 + 1.0), // Any overlap
            (opt.fragment_lengths.max_fragment_length + opt.end_offset).into(),
        )?;
        let overlapping_windows = if let Some(overlaps) = overlapping_windows {
            overlaps
        } else {
            continue;
        };

        counter.base.counted_fragments += 1;

        for overlapped_window in overlapping_windows.windows {
            let idx = overlapped_window.idx;
            let counts = &mut counts_by_bin[idx];
            count_kmers_in_segments(
                &fragment,
                &positional_codes_by_k,
                kmer_specs,
                counts,
                positional_scaling_weights.as_ref(),
            );
        }
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    let bin_info = if let Some(size) = opt.windows.by_size {
        // Build bin information for chromosome
        // chrom,start,end,blacklist_overlap
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
        let windows = windows.unwrap();
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

    Ok((counts_by_bin, bin_info, counter))
}

fn count_kmers_in_segments(
    fragment: &FragmentWithKmerSegments,
    positional_codes_by_k: &FxHashMap<u8, KmerCodes>,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    counts: &mut FxHashMap<Kmer, f64>,
    weights: Option<&Vec<f32>>,
) {
    for (&k, _) in kmer_specs {
        let codes = positional_codes_by_k
            .get(&k)
            .expect("missing positional codes for requested k");
        let k_span = k as u32;

        for (seg_start, seg_end) in &fragment.segments {
            let Some(last_start) = seg_end.checked_sub(k_span) else {
                continue;
            };

            for idx in *seg_start..=last_start {
                let idx_usize = idx as usize;
                let w = weights.map_or(1.0, |w| unsafe { *w.get_unchecked(idx_usize) });

                *counts
                    .entry(Kmer {
                        k,
                        code: codes.get(idx_usize),
                    })
                    .or_insert(0.) += w as f64;
            }
        }
    }
}
