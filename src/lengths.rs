use anyhow::{Context, Result};
use fragsizer::cli::Count;
use fragsizer::cli::counters::FragsizeExtractionCounters;
use fragsizer::cli::io::create_chromosome_reader;
use fragsizer::cli::io::read_seq;
use fragsizer::cli::opts::*;
use fragsizer::extractors::bed::load_windows;
use fragsizer::extractors::blacklist::*;
use fragsizer::extractors::gc::*;
use ndarray::{Array3, s};
use ndarray_npy::{read_npy, write_npy};
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
        group = clap::ArgGroup::new("windows")
            .required(true)
            .args(&["by_size", "by_bed"])
    )
)]
struct LengthsConfig {
    #[cfg_attr(feature = "cli", clap(flatten))]
    io: IOArgs,

    /// Window definition: a fixed window size [integer]
    #[cfg_attr(
        feature = "cli",
        clap(long = "by-size", alias = "by", value_parser, group = "windows")
    )]
    by_size: Option<u64>,

    /// Window definition: a BED file of windows [path]
    #[cfg_attr(
        feature = "cli",
        clap(long = "by-bed", value_parser, group = "windows")
    )]
    by_bed: Option<PathBuf>,

    // /// Window definition: a single genome-wide window [flag]
    // #[cfg_attr(feature = "cli", clap(long = "global", group = "windows"))]
    // global: bool,
    #[cfg_attr(feature = "cli", clap(flatten))]
    lengths: LengthArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    filter: FilterArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    blacklisting: BlacklistArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    gc: GCArgs,

    #[cfg_attr(feature = "cli", clap(flatten))]
    two_bit: TwoBitArgs,
}


/// Whether to include the read or continue
fn include_read(rec: &Record, opt: &LengthsConfig) -> bool {
    if rec.is_unmapped()
        || rec.is_mate_unmapped()
        || rec.tid() != rec.mtid()
        || rec.is_secondary()
        || rec.is_supplementary()
        || rec.is_duplicate()
        || rec.is_quality_check_failed()
        || (opt.filter.require_proper_pair && !rec.is_proper_pair())
        || rec.mapq() < opt.filter.min_mapq
    {
        return false;
    }
    true
}

pub fn run(opt: LengthsConfig) -> Result<()> {
    let start_time = Instant::now();
    // Create output directory
    create_dir_all(&opt.io.output_dir).context("Cannot create output_dir")?;

    // Load blacklist intervals if provided
    let blacklist_map = if let Some(beds) = &opt.blacklisting.blacklist {
        load_blacklists(beds, opt.blacklisting.blacklist_min_size, AUTOSOMES)?
    } else {
        HashMap::new()
    };

    // Determine number of length bins and GC bins
    let num_lengths = (opt.lengths.max_length - opt.lengths.min_length + 1) as usize;
    let gc_bins = if opt.gc.bin_by_gc {
        let range = (opt.gc.gc_max_pct as i32 - opt.gc.gc_min_pct as i32).max(0) as u32;
        ((range + opt.gc.gc_bin_size_pct as u32) / opt.gc.gc_bin_size_pct as u32) as usize
    } else {
        1
    };

    let windows_map = if let Some(bed) = &opt.by_bed {
        Some(load_windows(bed, &AUTOSOMES)?)
    } else {
        None
    };

    // Cconfigure global thread‐pool size
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.io.n_threads as usize)
        .build_global()
        .context("building Rayon thread pool")?;

    // Prepare per-bin metadata
    let mut bin_info = Vec::new();
    let mut global_counter = FragsizeExtractionCounters::default();

    // Chromosome-wise tmp paths
    let tmp_paths: HashMap<&str, PathBuf> = AUTOSOMES
        .iter()
        .map(|&chr| {
            let path = tmp_dir.join(format!("{}_counts.npy", chr));
            (chr, path)
        })
        .collect();

    // Main loop: process each autosome

    let results: Vec<(
        Vec<(String, u64, u64, u64, f64)>,
        FragsizeExtractionCounters,
    )> = AUTOSOMES
        .par_iter()
        .enumerate()
        .map(|(i, &chr)| -> Result<(_, _)> {
            Ok(process_chrom(
                i,
                chr,
                &opt,
                &tmp_paths.get(chr).expect("missing tmp path for chromosome"),
                num_lengths,
                windows_map
                    .as_ref()
                    .and_then(|m| m.get(chr).map(|v| v.as_slice())),
                gc_bins,
                blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]),
            )?)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    // 3) collect results (in AUTOSOME order) back into your global vectors
    for (bin_vec, counter) in results {
        bin_info.extend(bin_vec);
        global_counter += counter;
    }

    // Stack all per-chrom arrays into one large array
    let total_bins = bin_info.len();
    let mut all_counts = Array3::<Count>::zeros((total_bins, num_lengths, gc_bins));
    let mut offset = 0;
    for chr in AUTOSOMES {
        let arr: Array3<Count> = read_npy(&tmp_paths[chr]).context("Read tmp npy fail")?;
        let bins = arr.dim().0;
        all_counts
            .slice_mut(s![offset..offset + bins, .., ..])
            .assign(&arr);
        offset += bins;
    }

    // Write final counts and BED file to output_dir
    write_npy(&opt.io.output_dir.join("all_counts.npy"), &all_counts)
        .context("Write final fail")?;
    let mut bed_writer = BufWriter::new(
        File::create(&opt.io.output_dir.join("bins.bed")).context("Create bed fail")?,
    );
    for (chr, start, end, total, overlap_perc) in &bin_info {
        writeln!(
            bed_writer,
            "{}\t{}\t{}\t{}\t{}",
            chr, start, end, total, overlap_perc
        )
        .context("Write bed line fail")?;
    }

    // Print summary statistics and execution time
    let elapsed = start_time.elapsed();
    println!("Total reads: {}", global_counter.total);
    println!(
        "Accepted reads: {} ({:.2}%)  # Only the first read in a pair is counted",
        global_counter.accepted,
        global_counter.accepted as f64 / global_counter.total as f64
    );
    println!("Blacklisted reads: {}", global_counter.blacklisted);
    if opt.gc.bin_by_gc {
        println!("GC-excluded reads: {}", global_counter.gc_excl);
    }
    println!("Counted reads: {}", global_counter.counted);
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn process_chrom(
    idx: usize,
    chr: &str,
    opt: &Cli,
    tmp_out_path: &Path,
    num_lengths: usize,
    windows: Option<&[(u64, u64)]>,
    gc_bins: usize,
    blacklist_intervals: &[(u64, u64)],
) -> anyhow::Result<(
    Vec<(String, u64, u64, u64, f64)>,
    FragsizeExtractionCounters,
)> {
    println!("Processing {}/22: {}", idx + 1, chr);
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.io.bam, chr)?;

    let seq_bytes = if opt.gc.bin_by_gc {
        Some(read_seq(opt.two_bit.ref_2bit.as_ref().unwrap(), chr)?)
    } else {
        None
    };

    let gc_prefix = seq_bytes
        .as_ref()
        .filter(|_| opt.gc.bin_by_gc)
        .map(|bytes| build_gc_prefix(bytes));

    // Initialize counters (default -> 0s)
    let mut counter = FragsizeExtractionCounters::default();

    let num_bins = if let Some(size) = opt.by_size {
        ((chrom_len + size - 1) / size) as usize
    } else {
        windows.unwrap().len()
    };

    // — per-chrom arrays and prefix-sum
    let mut arr = Array3::<Count>::zeros((num_bins, num_lengths, gc_bins));

    // Streaming pointers and single fetch for this chr
    let mut bl_ptr = 0; // blacklist interval
    let mut wd_ptr = 0; // genomic window

    reader
        .fetch((tid, 0, chrom_len as i64))
        .context(format!("fetch {}", chr))?;

    // Loop over records and count
    for res in reader.records() {
        let rec = res.context("reading bam record")?;
        counter.total += 1;
        if rec.tid() != tid as i32 {
            continue;
        }

        // Apply basic filters and get fragment coordinates
        if let Some((s, e, tlen)) = filter_read(&rec, &opt) {
            counter.accepted += 1;
            // Determine blacklist status
            let in_blacklist = is_blacklisted(
                blacklist_intervals,
                opt.blacklisting.blacklist_strategy.clone(),
                s,
                e,
                &mut bl_ptr,
            );
            if in_blacklist {
                counter.blacklisted += 1;
            }

            // Compute GC bin index via prefix-sum (O(1) range sum)
            let gc_bin = if let Some(pref) = &gc_prefix {
                let gc_cnt = pref[e as usize] - pref[s as usize];
                let pct = ((gc_cnt as f64 / tlen as f64) * 100.0).round() as i32;
                if pct < opt.gc.gc_min_pct as i32 || pct > opt.gc.gc_max_pct as i32 {
                    counter.gc_excl += 1;
                    continue;
                }
                ((pct as u32 - opt.gc.gc_min_pct as u32) / opt.gc.gc_bin_size_pct as u32) as usize
            } else {
                0
            };

            if let Some(size) = opt.by_size {
                // Assign fragment to window bin by midpoint
                let mid = s + (tlen as u64) / 2;
                let bin_idx: usize = (mid / size) as usize;
                let len_idx = (tlen - opt.lengths.min_length) as usize;
                arr[(bin_idx, len_idx, gc_bin)] += 1;
                counter.counted += 1;
            } else {
                let windows = windows.unwrap();
                let mid = s + (tlen as u64) / 2;
                let len_idx = (tlen - opt.lengths.min_length) as usize;
                // Skip any intervals that end entirely before the fragment start
                while wd_ptr < windows.len() && windows[wd_ptr].1 <= s {
                    wd_ptr += 1;
                }
                // Iterate over every interval that could cover the midpoint
                let mut bin_idx = wd_ptr;
                let mut overlaps: Vec<usize> = vec![];
                while bin_idx < windows.len() && windows[bin_idx].0 <= mid {
                    let (s, e) = windows[bin_idx];
                    if s <= mid && mid < e {
                        overlaps.push(bin_idx);
                    }
                    bin_idx += 1;
                }
                // Only if we get an overlapping window do
                for overlap_idx in overlaps.iter() {
                    arr[(*overlap_idx, len_idx, gc_bin)] += 1;
                }
                if overlaps.len() > 0 {
                    counter.counted += 1;
                }
            };
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
            let total: u64 = arr.slice(s![b, .., ..]).map(|&v| v as u64).sum();
            let overlap_perc =
                compute_blacklist_overlap(blacklist_intervals, start, end, &mut bl_ptr);
            bin_info.push((chr.to_string(), start, end, total, overlap_perc));
        }
        bin_info
    } else {
        // build bin_info from the exact BED windows
        let mut bl_ptr = 0;
        let windows = windows.unwrap();
        let mut bin_info = Vec::with_capacity(num_bins);
        for (b, (wstart, wend)) in windows.iter().cloned().enumerate() {
            let total: u64 = arr.slice(s![b, .., ..]).map(|&v| v as u64).sum();
            let overlap_perc =
                compute_blacklist_overlap(blacklist_intervals, wstart, wend, &mut bl_ptr);
            bin_info.push((chr.to_string(), wstart, wend, total, overlap_perc));
        }
        bin_info
    };

    // — write out your .npy
    write_npy(&tmp_out_path, &arr).context("write npy")?;

    Ok((bin_info, counter))
}
