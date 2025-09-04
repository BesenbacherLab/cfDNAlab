use crate::{
    cli_common::*,
    utils::{bed::load_windows_from_bed, blacklist::load_blacklists},
};
use anyhow::{Context, Result};
use ndarray::{Array3, s};
use ndarray_npy::write_npy;
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use std::{
    collections::HashMap,
    fs::{File, create_dir_all},
    io::{BufWriter, Write},
    path::Path,
    path::PathBuf,
    time::Instant,
};
use tempfile::Builder;

#[cfg_attr(feature = "cli", derive(clap::Args))]
#[cfg_attr(
    feature = "cli",
    clap(
        group = clap::ArgGroup::new("min_acgt")
            .args(&["min_acgt_pct", "min_acgt_count"])
            .multiple(true)))]
struct GCConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    ioc: IOCArgs,

    /// 2bit reference file [path]
    /// E.g., "hg38.2bit"
    #[cfg_attr(
        feature = "cli",
        clap(
            short = 'r',
            long,
            clap::value_parser,
            required = true,
            help_heading = "Core"
        )
    )]
    pub ref_2bit: PathBuf,

    #[cfg_attr(feature = "cli", clap(flatten))]
    windows: WindowsArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    chromosomes: ChromosomeArgs,

    /// Optional BED files of blacklisted regions [path]
    #[cfg_attr(
        feature = "cli",
        clap(short = 'b', long, value_parser, num_args = 1.., action = clap::ArgAction::Append, help_heading="Filtering"))]
    pub blacklist: Option<Vec<PathBuf>>,

    /// Minimum mapping quality to include [integer]
    #[cfg_attr(
        feature = "cli",
        clap(long, alias = "mq", default_value = "60", value_parser = value_parser!(u8).range(0..), help_heading="Filtering"))]
    pub min_mapq: u8,

    /// Only count properly paired reads [flag]
    #[cfg_attr(feature = "cli", clap(long))]
    pub require_proper_pair: bool,

    /// Minimum fragment length to include (default: 20) [integer]
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20", value_parser = value_parser!(u16).range(1..), help_heading="Filtering"))]
    pub min_fragment_length: u16,

    /// Maximum fragment length to include (default: 600) [integer]
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "600", value_parser = value_parser!(u16).range(1..), help_heading="Filtering"))]
    pub max_fragment_length: u16,

    /// Minimum GC % to consider [integer]
    ///
    /// Fragments with lower GC % are ignored.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "0", 
             value_parser = value_parser!(u8).range(0..100), help_heading="Filtering"))]
    pub gc_min_pct: u8,

    /// Maximum GC % to consider [integer]
    ///
    /// Fragments with higher GC % are ignored.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "100", 
             value_parser = value_parser!(u8).range(0..101), help_heading="Filtering"))]
    pub gc_max_pct: u8,

    /// Minimum **percentage** of ACGT bases in a fragment after blacklist exclusion [integer]
    ///
    /// Fragments where a lower percentage of bases are ACGT (not blacklisted or 'N') are ignored.
    ///
    /// When both `min_acgt_*` arguments are specified, the highest threshold (most remaining bases) must be met.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "0", group = "min_acgt", 
             value_parser = value_parser!(u8).range(0..101), help_heading="Minimum ACGT (select 0-2 args)"))]
    pub min_acgt_pct: u8,

    /// Minimum **count** of ACGT bases in a fragment after blacklist exclusion [integer]
    ///
    /// Fragments where fewer bases are ACGT (not blacklisted or 'N') are ignored.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "0", group = "min_acgt", 
             value_parser = value_parser!(u8).range(0..), help_heading="Minimum ACGT (select 0-2 args)"))]
    pub min_acgt_count: u8,
}

/// Whether to include the read or continue
fn include_read(rec: &Record, opt: &GCConfig) -> bool {
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

fn run(opt: GCConfig) -> Result<()> {
    let start_time = Instant::now();
    let chromosomes = opt.chromosomes.resolve_chromosomes()?;
    let window_opt = opt.windows.resolve_windows()?;
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    // Create output directory
    create_dir_all(&opt.io.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = if let Some(beds) = &opt.blacklist {
        println!("Start: Loading blacklists");
        load_blacklists(beds, 1, &chromosomes)?
    } else {
        HashMap::new()
    };

    let windows_map = match window_opt {
        WindowSpec::Bed => {
            println!("Start: Loading window coordinates");
            Some(load_windows_from_bed(bed, &chromosomes)?)
        }
        _ => None,
    };

    // Configure global thread‐pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.ioc.n_threads as usize)
        .build_global()
        .context("building Rayon thread pool")?;

    // Prepare per-bin counts and metadata
    let mut all_bins = Vec::new();
    let mut bin_info = Vec::new();
    let mut global_counter = ConsensusDepthCounters::default();

    // Main loop: process each autosome
    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let results: Vec<(
        Vec<FxHashMap<Kmer, BigCount>>,
        Option<Vec<(String, u64, u64, u64, f64)>>,
        ConsensusDepthCounters,
    )> = chromosomes
        .par_iter()
        .map(|chr| -> Result<(_, _, _)> {
            let out = process_chrom(
                &chr,
                &opt,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    println!("Start: Processing counts");

    // Collect results (in chromosome order) back into the global vectors
    for (counts_by_bin, bin_vec, counter) in results {
        let counts_decoded: Vec<DecodedCounts> = counts_by_bin
            .iter()
            .map(|c| split_and_decode_counts(c, &kmer_specs))
            .collect();
        all_bins.extend(counts_decoded);
        if !opt.global {
            bin_info.extend(bin_vec.unwrap());
        }
        global_counter += counter;
    }

    // Convert to single hashmap for global
    // Keep wrapped in vector to simplify writer
    let all_bins = if opt.global {
        vec![merge_decoded_counts(all_bins)]
    } else {
        all_bins
    };

    // Prepare counts to get correct motifs (collapsed, N-filtered, etc.)
    let mut prepared_counts = prepare_decoded_counts(&all_bins, !opt.strand_aware, &kmer_specs);

    // Sort by original index (when given a bed file)
    if opt.by_bed.is_some() {
        println!("Start: Reordering counts by original window index in bed file");

        // Zip into a single Vec
        let mut paired: Vec<_> = bin_info
            .into_iter()
            .zip(prepared_counts.into_iter())
            .collect(); // (BinInfo, DecodedCounts)

        // Sort primarily by original window index
        paired.sort_unstable_by_key(|(info, _)| info.3);

        // Unzip back out if you need separate Vecs again
        (bin_info, prepared_counts) = paired.into_iter().unzip();
    }

    println!("Start: Calculating mismatch rates");

    // Calculate mismatch rates
    let mismatch_rates_by_k = compute_mismatch_rates(&prepared_counts);

    println!("Start: Writing mismatch rates to disk");

    write_mismatch_rate_tsvs(&prepared_counts, &mismatch_rates_by_k, &opt.io.output_dir)?;

    // Write bins BED file
    if !opt.global {
        println!("Start: Writing window coordinates to disk");
        let mut bed_writer = BufWriter::new(
            File::create(&opt.io.output_dir.join("bins.bed")).context("Create bed fail")?,
        );
        for (chr, start, end, _, overlap_perc) in &bin_info {
            writeln!(bed_writer, "{}\t{}\t{}\t{}", chr, start, end, overlap_perc)
                .context("Write bed line fail")?;
        }
    }

    println!("Statistics:");

    // Print summary statistics and execution time
    let elapsed = start_time.elapsed();
    println!("  Total reads: {}", global_counter.total);
    println!(
        "  Accepted reads: {} ({:.2}%) | Left reads: {} | Right mates: {}",
        global_counter.accepted,
        global_counter.accepted as f64 / global_counter.total as f64 * 100.0,
        global_counter.left,
        global_counter.right_mate
    );
    // if opt.gc.bin_by_gc {
    //     println!("GC-excluded reads: {}", global_counter.gc_excl);
    // }
    if global_counter.missing_md > 0 {
        println!("  Reads missing MD tag: {}", global_counter.missing_md);
    }
    println!(
        "  Fragments counted one or more times: {}",
        global_counter.counted
    );
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn process_chrom(
    chr: &str,
    opt: &GCConfig,
    windows: Option<&[(u64, u64, u64)]>,
    // gc_bins: usize,
    blacklist_intervals: &[(u64, u64)],
) -> anyhow::Result<(
    Vec<FxHashMap<Kmer, BigCount>>,
    Option<Vec<(String, u64, u64, u64, f64)>>,
    ConsensusDepthCounters,
)> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, chr)?;

    let mut seq_bytes = read_seq(&opt.ref_2bit, chr)?;
    apply_blacklist_mask_to_seq(&mut seq_bytes, &blacklist_intervals);

    let positional_codes_by_k: HashMap<u8, KmerCodes> = build_codes_per_k(&seq_bytes, kmer_specs);

    // Not needed yet
    // let gc_prefix = build_gc_prefix(&seq_bytes);

    // Initialize counters (default -> 0s)
    let mut counter = ConsensusDepthCounters::default();

    let num_bins = if let Some(sz) = opt.by_size {
        // by-size
        ((chrom_len + sz - 1) / sz) as usize
    } else if opt.by_bed.is_some() {
        // by-bed
        windows.unwrap().len()
    } else {
        // global
        1
    };

    let mut counts_by_bin = vec![FxHashMap::<Kmer, BigCount>::default(); num_bins];

    let mut stash: HashMap<Vec<u8>, ReadInfo> = HashMap::new();

    // Streaming pointers and single fetch for this chr
    let mut wd_ptr = 0; // genomic window

    reader
        .fetch((tid, 0, chrom_len as i64))
        .context(format!("fetch {}", chr))?;

    for res in reader.records() {
        let rec = res.context("reading bam record")?;
        counter.total += 1;
        if rec.tid() != tid as i32 || filter_read(&rec, &opt.read_filters).is_none() {
            continue;
        }
        counter.accepted += 1;

        let ascii = rec.seq().as_bytes();
        let qname = rec.qname().to_vec();
        let base_qualities = rec.qual().to_owned();

        // Check the number of mismatches (NM tag)
        // We use this to check whether to parse the MD tag
        let nm_opt = read_nm_tag(&rec);
        let nm_tag = match nm_opt {
            Some(nm) if opt.read_filters.min_nm <= nm && nm <= opt.read_filters.max_nm => nm,
            _ => continue,
        };

        // TODO: cfDNA-PRO talks about these reverse fragments where
        // perhaps the right read starts before the left read somehow?
        // Does that affect whether "insert_size() > 0" is always the
        // most-left read and/or do we skip this subset of fragments?

        if rec.insert_size() > 0 {
            // Left-most read
            let s = rec.pos() as u64;
            let e = s + rec.seq_len() as u64;
            counter.left += 1;

            let md = if nm_tag > 0 {
                // Only try to parse MD when nm_tag indicates there could be mismatches
                match read_md_tag(&rec) {
                    Some(s) => s,
                    None => {
                        counter.missing_md += 1;
                        continue;
                    }
                }
            } else {
                // No mismatches by NM → we can use an empty MD string (or skip entirely later)
                String::new()
            };

            // Left-most read → stash
            stash.insert(
                qname,
                ReadInfo {
                    seq: ascii,
                    pos: s,
                    end: e,
                    base_qualities,
                    md_tag: md,
                },
            );
        } else if let Some(left) = stash.remove(&qname) {
            counter.right_mate += 1;

            // .pos() is always most left so sequence must start there
            let right_pos = rec.pos() as u64;
            let right_end = right_pos + rec.seq_len() as u64;

            let overlaps_opt = find_overlaps(
                chrom_len,
                &mut wd_ptr,
                windows,
                opt.by_size,
                left.pos,
                left.end,
                right_pos,
                right_end,
                opt.min_overlap as usize,
                opt.read_filters.max_fragment_length as u64,
            );

            let overlaps = if let Some(overlaps) = overlaps_opt {
                overlaps
            } else {
                continue;
            };

            // Get MD tag
            let md = if nm_tag > 0 {
                // Only try to parse MD when nm_tag indicates there could be mismatches
                match read_md_tag(&rec) {
                    Some(s) => s,
                    None => {
                        counter.missing_md += 1;
                        continue; // diverges, so this arm never needs to return a String
                    }
                }
            } else {
                // No mismatches by NM → we can use an empty MD string (or skip entirely later)
                String::new()
            };

            let mismatch_coordinates = if !md.is_empty() && !left.md_tag.is_empty() {
                // Parse both MD tags (after read filtering)
                let (mismatch_starts, mismatch_ends) = parse_md_tag(&md, rec.pos() as u32);
                let (left_mismatch_starts, left_mismatch_ends) =
                    parse_md_tag(&left.md_tag, left.pos as u32);

                let mismatch_coordinates: Vec<(u32, u32)> = intersect_mismatch_runs(
                    &left_mismatch_starts,
                    &left_mismatch_ends,
                    &mismatch_starts,
                    &mismatch_ends,
                )
                .into_iter()
                .filter(|&(s, e)| e - s <= 2)
                .collect();
                Some(mismatch_coordinates)
            } else {
                None
            };

            counter.counted += 1;

            // TODO: Increase documentation of this - using overlaps.overlap_start, etc.

            // TODO: If the windows overlap many times, we are doing a lot of redundant
            // processing when counting per window separately
            // Probably faster to pass the overlapping windows to the counter
            // and count in each of them after we know what to count for a position
            for overlapped_window in overlaps.windows.iter() {
                let bmap = &mut counts_by_bin[overlapped_window.idx.clone()];

                // Calculate the positions in the reads' vectors to use
                // for the current window
                let (left_offset, right_offset, count_start, count_size) =
                    calculate_overlap_coordinates(
                        left.pos,
                        right_pos,
                        overlaps.overlap_start,
                        overlaps.overlap_end,
                        overlapped_window.win_start,
                        overlapped_window.win_end,
                    );

                // Check consensus and return trinucleotide key per position
                count_overlap(
                    &left.seq[left_offset..left_offset + count_size],
                    &left.base_qualities[left_offset..left_offset + count_size],
                    &ascii[right_offset..right_offset + count_size],
                    &base_qualities[right_offset..right_offset + count_size],
                    count_start,
                    kmer_specs,
                    &positional_codes_by_k,
                    mismatch_coordinates.as_deref(),
                    bmap,
                    opt.min_base_quality,
                );
            }
        }
    }

    let bin_info = if let Some(size) = opt.by_size {
        // Build bin information for chromosome
        // chrom,start,end,total_count
        let mut bl_ptr = 0;
        let mut bin_info = Vec::with_capacity(num_bins);
        for b in 0..num_bins {
            let start = b as u64 * size;
            let end = (start + size).min(chrom_len);
            let overlap_perc =
                compute_blacklist_overlap(blacklist_intervals, start, end, &mut bl_ptr);
            // Note: b (index) is a placeholder that is removed later
            bin_info.push((chr.to_string(), start, end, b as u64, overlap_perc));
        }
        Some(bin_info)
    } else if opt.by_bed.is_some() {
        // build bin_info from the exact BED windows
        let mut bl_ptr = 0;
        let windows = windows.unwrap();
        let mut bin_info = Vec::with_capacity(num_bins);
        for (_b, (wstart, wend, original_win_idx)) in windows.iter().cloned().enumerate() {
            let overlap_perc =
                compute_blacklist_overlap(blacklist_intervals, wstart, wend, &mut bl_ptr);
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
