use crate::{
    cli_common::{
        ChromosomeArgs, FragmentLengthArgs, IOCArgs, Ref2BitRequiredArgs, ScaleGenomeArgs,
        WindowSpec, WindowsArgs,
    },
    counters::FragmentKmersCounters,
    utils::{
        bam::{Contigs, create_chromosome_reader},
        bed::{Windows, load_windows_from_bed},
        blacklist::{
            BlacklistStrategy, apply_blacklist_mask_to_seq, apply_mask::BLACKLIST_BYTE,
            compute_blacklist_overlap, is_blacklisted,
        },
        command::{
            ensure_output_dir, load_blacklist_map, load_scaling_map,
            resolve_chromosomes_and_contigs,
        },
        coverage::{
            scale_genome::apply_scaling_to_coverage_in_place,
            tiled_run::{
                Tile, TileWindowSpan, build_tiles, clamp_fetch_to_window_span, make_temp_dir,
                precompute_tile_window_spans, tile_window_min_max,
            },
        },
        fragment::segment_kmer_fragment::FragmentWithKmerSegments,
        fragment_iterator::fragments_with_kmer_segments_from_bam,
        indel_mode::IndelMode,
        kmers::{
            kmer_codec::{
                Kmer, KmerCodes, KmerSpec, build_kmer_specs, build_left_aligned_codes_per_k,
            },
            process_counts::{DecodedCounts, prepare_decoded_counts, split_and_decode_counts},
            write::write_decoded_counts_matrix,
        },
        overlaps::find_overlapping_windows,
        read::default_include_read,
        reference::read_seq_in_range,
        thread_pool::init_global_pool,
    },
};
use anyhow::{Context, Result, bail};
use bincode::{
    config::standard,
    serde::{decode_from_std_read, encode_into_std_write},
};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rust_htslib::bam::{Read, Record};
use serde::{Deserialize, Serialize};
use std::{
    convert::TryInto,
    fs::{self, File},
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
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

    /// Size of tiles to parallelize over `[integer]`
    ///
    /// Chromosomes are processed in tiles of this size to reduce memory usage.
    #[cfg_attr(
        feature = "cli",
        clap(long, default_value = "20000000", value_parser = clap::value_parser!(u32).range(1000000..), help_heading="Core"))]
    pub tile_size: u32,

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
            tile_size: 20_000_000,
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

    pub fn set_tile_size(&mut self, tile_size: u32) {
        self.tile_size = tile_size;
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

    pub fn set_scale_genome(&mut self, scale: ScaleGenomeArgs) {
        self.scale_genome = scale;
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

#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Serialize, Deserialize)]
pub struct TileKmerCountEntry {
    pub k: u8,
    pub code: u64,
    pub value: f64,
}

#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Serialize, Deserialize)]
pub struct TileWindowCounts {
    pub original_idx: u64,
    pub entries: Vec<TileKmerCountEntry>,
}

/// Persist per-tile k-mer counts so they can be merged after parallel tile processing.
fn serialize_tile_counts(path: &Path, payload: &[TileWindowCounts]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("creating tile counts file: {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(512 * 1024, file);
    encode_into_std_write(payload, &mut writer, standard())
        .with_context(|| format!("serialising tile counts to {}", path.display()))?;
    writer.flush().with_context(|| {
        format!(
            "flushing tile counts file after serialisation: {}",
            path.display()
        )
    })
}

/// Load counts created by [`serialize_tile_counts`] during the reduction phase.
fn deserialize_tile_counts(path: &Path) -> Result<Vec<TileWindowCounts>> {
    let file = File::open(path)
        .with_context(|| format!("opening tile counts file: {}", path.display()))?;
    let mut reader = BufReader::with_capacity(512 * 1024, file);
    decode_from_std_read(&mut reader, standard())
        .with_context(|| format!("deserialising tile counts from {}", path.display()))
}

/// Per-tile bookkeeping for intermediate count files and fragment counters.
struct TileResult {
    chr: String,
    counts_path: Option<PathBuf>,
    counter: FragmentKmersCounters,
}

/// Lightweight view into the window configuration for a given tile.
///
/// Stores the context needed to convert chromosome-local indices into global window ids for a tile.
struct WindowContext<'a> {
    spec: &'a WindowSpec,
    windows: Option<&'a [(u64, u64, u64)]>,
    chr_idx_offset: u64,
}

impl<'a> WindowContext<'a> {
    #[inline]
    /// Return the per-chromosome windows slice when operating in BED mode.
    fn windows_slice(&self) -> Option<&'a [(u64, u64, u64)]> {
        self.windows
    }

    #[inline]
    /// Map the provided chromosome-local window index to the global window identifier expected
    /// downstream.
    ///
    /// Parameters
    /// -----------
    /// `chrom_window_idx`: Index of the window relative to the *chromosome*, as supplied by
    /// [`find_overlapping_windows`] (it counts from the start of the chromosome, not the tile).
    ///
    fn original_idx(&self, chrom_window_idx: usize) -> u64 {
        match self.spec {
            WindowSpec::Global => 0,
            WindowSpec::Size(window_bp) => {
                debug_assert_ne!(*window_bp, 0);
                self.chr_idx_offset
                    .checked_add(chrom_window_idx as u64)
                    .expect("window index overflow for size-based windows")
            }
            WindowSpec::Bed(_) => {
                self.windows.expect("windows slice required for BED mode")[chrom_window_idx].2
            }
        }
    }
}

/// Determine the genomic span to request from the BAM reader for a tile.
fn determine_fetch_span(
    tile: &Tile,
    window_ctx: &WindowContext,
    tile_window_span: Option<&TileWindowSpan>,
    chrom_len: u64,
) -> Option<(i64, i64)> {
    let chrom_len_u32 = chrom_len.min(u32::MAX as u64) as u32;
    match window_ctx.spec {
        WindowSpec::Global | WindowSpec::Size(_) => {
            Some((tile.fetch_start as i64, tile.fetch_end as i64))
        }
        WindowSpec::Bed(_) => {
            let windows = window_ctx.windows_slice()?;
            let (min_ws, max_we) = tile_window_min_max(windows, tile, tile_window_span)?;
            clamp_fetch_to_window_span(tile, chrom_len.min(chrom_len_u32 as u64), min_ws, max_we)
        }
    }
}

/// Compute the global window count together with per-chromosome index offsets.
///
/// Returns `(total_windows, chr_idx_offsets)` where `total_windows` is the expected length of the
/// window output vector and `chr_idx_offsets` maps each chromosome to the index of its first window
/// in that global vector.
///
/// Mode specifics:
/// - `Global`: there is exactly one window covering the whole genome; all chromosomes map to offset
///   `0`.
/// - `Size(window_bp)`: windows are created from the chromosome start in contiguous bins, so the
///   offset for chromosome *k* is the cumulative number of bins emitted for chromosomes `< k`.
/// - `Bed`: offsets remain `0` because each BED entry already carries its own globally unique
///   `original_idx`; consumers must use that `original_idx` when addressing global arrays.
fn compute_window_offsets(
    window_opt: &WindowSpec,
    chromosomes: &[String],
    contigs: &Contigs,
    windows_map: Option<&FxHashMap<String, Windows>>,
) -> Result<(u64, FxHashMap<String, u64>)> {
    let mut offsets: FxHashMap<String, u64> = FxHashMap::default();

    match window_opt {
        WindowSpec::Global => {
            for chr in chromosomes {
                offsets.insert(chr.clone(), 0);
            }
            Ok((1, offsets))
        }
        WindowSpec::Size(size) => {
            let mut running_window_idx_offset = 0u64; // number of windows seen so far
            for chr in chromosomes {
                offsets.insert(chr.clone(), running_window_idx_offset);
                let &(_, len_u32) = contigs
                    .contigs
                    .get(chr)
                    .with_context(|| format!("missing contig length for '{}'", chr))?;
                let len = len_u32 as u64;
                let bins = if len == 0 {
                    0
                } else {
                    (len + *size - 1) / *size
                };
                running_window_idx_offset = running_window_idx_offset.saturating_add(bins);
            }
            Ok((running_window_idx_offset, offsets))
        }
        WindowSpec::Bed(_) => {
            let win_map =
                windows_map.with_context(|| "window map required for --by-bed mode".to_string())?;
            let mut total = 0u64;
            for chr in chromosomes {
                // BED entries already encode their global `original_idx`, so the reducer should use
                // that value instead of a chromosome offset.
                offsets.insert(chr.clone(), 0);
                if let Some(windows) = win_map.get(chr) {
                    total = total.saturating_add(windows.len() as u64);
                }
            }
            Ok((total, offsets))
        }
    }
}

/// Build per-window metadata (coordinates, blacklist overlap, etc.) for downstream consumers.
///
/// When running in BED mode the `original_idx` embedded in the loaded windows is preserved so the
/// caller must continue using that identifier when addressing global vectors.
fn build_bin_info(
    window_opt: &WindowSpec,
    chromosomes: &[String],
    contigs: &Contigs,
    windows_map: Option<&FxHashMap<String, Windows>>,
    blacklist_map: &FxHashMap<String, Vec<(u64, u64)>>,
    chr_offsets: &FxHashMap<String, u64>,
) -> Result<Vec<(String, u64, u64, u64, f64)>> {
    let mut out = Vec::new();

    match window_opt {
        WindowSpec::Global => Ok(out),
        WindowSpec::Size(size) => {
            for chr in chromosomes {
                let &(_, len_u32) = contigs
                    .contigs
                    .get(chr)
                    .with_context(|| format!("missing contig length for '{}'", chr))?;
                let len = len_u32 as u64;
                let mut start = 0u64;
                let mut local_idx = 0u64;
                let mut bl_ptr = 0usize;
                let blacklist_intervals =
                    blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]);
                let chr_window_idx_offset = *chr_offsets.get(chr).unwrap_or(&0);

                while start < len {
                    let end = (start + *size).min(len);
                    let overlap =
                        compute_blacklist_overlap(blacklist_intervals, start, end, 0, &mut bl_ptr);
                    out.push((
                        chr.clone(),
                        start,
                        end,
                        chr_window_idx_offset + local_idx,
                        overlap,
                    ));
                    start += *size;
                    local_idx += 1;
                }
            }
            Ok(out)
        }
        WindowSpec::Bed(_) => {
            let win_map =
                windows_map.with_context(|| "window map required for --by-bed mode".to_string())?;
            for chr in chromosomes {
                let windows = win_map.get(chr).map(|w| w.as_slice()).unwrap_or(&[]);
                let mut bl_ptr = 0usize;
                let blacklist_intervals =
                    blacklist_map.get(chr).map(|v| v.as_slice()).unwrap_or(&[]);
                for &(start, end, original_idx) in windows {
                    let overlap =
                        compute_blacklist_overlap(blacklist_intervals, start, end, 0, &mut bl_ptr);
                    out.push((chr.clone(), start, end, original_idx, overlap));
                }
            }
            out.sort_unstable_by_key(|entry| entry.3);
            Ok(out)
        }
    }
}

/// Reduce per-tile count payloads into a dense vector aligned with the global window order.
#[cfg_attr(not(test), doc(hidden))]
pub fn merge_tile_counts<I>(
    payloads: I,
    total_windows: usize,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
) -> Result<Vec<DecodedCounts>>
where
    I: IntoIterator<Item = Vec<TileWindowCounts>>,
{
    let mut aggregated_counts: FxHashMap<u64, FxHashMap<Kmer, f64>> = FxHashMap::default();

    for payload in payloads {
        for window_counts in payload {
            let entry = aggregated_counts
                .entry(window_counts.original_idx)
                .or_insert_with(FxHashMap::default);
            for count in window_counts.entries {
                let kmer = Kmer {
                    k: count.k,
                    code: count.code,
                };
                *entry.entry(kmer).or_insert(0.0) += count.value;
            }
        }
    }

    let empty_counts: FxHashMap<Kmer, f64> = FxHashMap::default();
    let mut all_bins: Vec<DecodedCounts> = Vec::with_capacity(total_windows);
    for idx in 0..total_windows {
        if let Some(counts) = aggregated_counts.remove(&(idx as u64)) {
            all_bins.push(split_and_decode_counts(&counts, kmer_specs));
        } else {
            all_bins.push(split_and_decode_counts(&empty_counts, kmer_specs));
        }
    }

    if !aggregated_counts.is_empty() {
        bail!(
            "Received counts for unexpected window indices: {:?}",
            aggregated_counts.keys().collect::<Vec<&u64>>()
        );
    }

    Ok(all_bins)
}

fn reduce_chromosome_tile_results(tile_results: Vec<TileResult>) -> Result<Vec<TileWindowCounts>> {
    let mut aggregated: FxHashMap<u64, FxHashMap<Kmer, f64>> = FxHashMap::default();

    for tile_result in tile_results {
        let Some(path) = tile_result.counts_path else {
            continue;
        };
        let tile_payload = deserialize_tile_counts(&path)?;
        let _ = fs::remove_file(&path);

        for window_counts in tile_payload {
            let entry = aggregated
                .entry(window_counts.original_idx)
                .or_insert_with(FxHashMap::default);
            for count in window_counts.entries {
                let kmer = Kmer {
                    k: count.k,
                    code: count.code,
                };
                *entry.entry(kmer).or_insert(0.0) += count.value;
            }
        }
    }

    let mut merged: Vec<TileWindowCounts> = aggregated
        .into_iter()
        .map(|(original_idx, counts_map)| {
            let mut entries: Vec<TileKmerCountEntry> = Vec::with_capacity(counts_map.len());
            for (kmer, value) in counts_map {
                entries.push(TileKmerCountEntry {
                    k: kmer.k,
                    code: kmer.code,
                    value,
                });
            }
            TileWindowCounts {
                original_idx,
                entries,
            }
        })
        .collect();

    merged.sort_unstable_by_key(|window| window.original_idx);
    Ok(merged)
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
    let prefix = opt.output_prefix.trim();

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

    // Build temporary directory
    let temp_dir = make_temp_dir(&opt.ioc.output_dir, prefix).context("create per-run temp dir")?;

    // Window size when --by-size (otherwise None)
    let by_size_bp: Option<u64> = match &window_opt {
        WindowSpec::Size(bp) => Some(*bp as u64),
        _ => None,
    };

    // Build tiles
    let halo_bp: u32 = opt.fragment_lengths.max_fragment_length; // Safe halo for pairing/segments
    let (tiles, _tile_and_window_boundaries_align) =
        build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, by_size_bp)?;

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(&tiles, |chr| {
        windows_lookup
            .and_then(|m| m.get(chr))
            .map(|w| w.as_slice())
            .unwrap_or(&[])
    }));

    // Compute per-chromosome window offsets and overall window count. In BED mode these offsets are
    // zero because windows already carry their global `original_idx` values.
    let (total_windows, chr_offsets_map) =
        compute_window_offsets(&window_opt, &chromosomes, &contigs, windows_map.as_ref())?;
    let chr_offsets = Arc::new(chr_offsets_map);

    let total_tiles = tiles.len();
    let temp_dir = Arc::new(temp_dir);

    // Create progress bar
    let pb = Arc::new(ProgressBar::new(total_tiles as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    // Configure global thread‐pool size
    init_global_pool(opt.ioc.n_threads as usize)?;

    println!("Start: Counting per chromosome");

    pb.set_position(0);

    let tile_window_spans_for_threads = tile_window_spans.clone();

    let tile_results: Vec<TileResult> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<TileResult> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let counts_path = temp_dir.join(format!(
                "{prefix}.{chr}.{idx}.counts.bin",
                prefix = prefix,
                chr = tile.chr.as_str(),
                idx = tile.index
            ));

            let window_ctx = WindowContext {
                spec: &window_opt,
                windows: windows_map
                    .as_ref()
                    .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice())),
                chr_idx_offset: *chr_offsets.get(&tile.chr).unwrap_or(&0),
            };

            let blacklist_chr = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            let out = process_tile(
                &opt,
                tile,
                &kmer_specs,
                &window_ctx,
                tile_span.as_ref(),
                blacklist_chr,
                scaling_chr,
                counts_path.as_path(),
            )?;
            pb.inc(1);
            Ok(out)
        })
        .collect::<Result<_>>()?; // short-circuits on the first Err

    pb.finish_with_message("| Finished counting");

    println!("Start: Reducing per-tile counts");

    let mut global_counter = FragmentKmersCounters::default();
    let mut tile_results_by_chr: FxHashMap<String, Vec<TileResult>> = FxHashMap::default();

    for tile_result in tile_results {
        global_counter += tile_result.counter;
        tile_results_by_chr
            .entry(tile_result.chr.clone())
            .or_default()
            .push(tile_result);
    }

    let mut payloads: Vec<Vec<TileWindowCounts>> = Vec::with_capacity(tile_results_by_chr.len());
    for chr in &chromosomes {
        if let Some(chr_tile_results) = tile_results_by_chr.remove(chr) {
            payloads.push(reduce_chromosome_tile_results(chr_tile_results)?);
        }
    }
    if !tile_results_by_chr.is_empty() {
        let unexpected_chr = tile_results_by_chr.keys().next().unwrap();
        bail!(
            "tile results produced for unexpected chromosome '{}'",
            unexpected_chr
        );
    }

    let total_windows_usize: usize = total_windows
        .try_into()
        .context("number of windows exceeds addressable size")?;

    let all_bins = merge_tile_counts(payloads, total_windows_usize, &kmer_specs)?;

    // Prepare counts to get correct motifs (collapsed, N-filtered, etc.)
    let (prepared_counts, motifs_by_k) =
        prepare_decoded_counts(&all_bins, opt.canonical, &kmer_specs);

    // Build bin metadata when windowed
    let bin_info = build_bin_info(
        &window_opt,
        &chromosomes,
        &contigs,
        windows_map.as_ref(),
        &blacklist_map,
        chr_offsets.as_ref(),
    )?;

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

/// Process a single tile: stream fragments, accumulate per-window counts, and persist results.
fn process_tile(
    opt: &FragmentKmersConfig,
    tile: &Tile,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    window_ctx: &WindowContext,
    tile_window_span: Option<&TileWindowSpan>,
    blacklist_intervals: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    counts_path: &Path,
) -> anyhow::Result<TileResult> {
    // Open a fresh BAM reader for this thread
    let (mut reader, tid, chrom_len) = create_chromosome_reader(&opt.ioc.bam, &tile.chr)?;

    let fetch_span = determine_fetch_span(tile, window_ctx, tile_window_span, chrom_len);
    let Some((fetch_from, fetch_to)) = fetch_span else {
        return Ok(TileResult {
            chr: tile.chr.clone(),
            counts_path: None,
            counter: FragmentKmersCounters::default(),
        });
    };

    let mut seq_bytes = read_seq_in_range(
        &opt.ref_genome.ref_2bit,
        &tile.chr,
        (tile.core_start as usize)..(tile.core_end as usize),
    )?;

    apply_blacklist_mask_to_seq(&mut seq_bytes, &blacklist_intervals, tile.core_start as u64);

    // Scaled weights to count up
    let positional_scaling_weights = if !scaling_chr.is_empty() {
        let mut scaling_weights = vec![1.0; seq_bytes.len()];
        apply_scaling_to_coverage_in_place(
            &mut scaling_weights,
            tile.core_start as u32,
            scaling_chr,
        );
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

    // Sparse map keyed by original window index -> kmer counts
    let mut counts_by_window: FxHashMap<u64, FxHashMap<Kmer, f64>> = FxHashMap::default();

    // Streaming pointers and single fetch for this chr
    let mut bl_ptr = 0; // Blacklist interval
    let mut wd_ptr = tile_window_span
        .and_then(|span| (!span.is_empty()).then_some(span.first_idx))
        .unwrap_or(0);

    reader
        .fetch((tid, fetch_from, fetch_to))
        .context(format!("fetch {}", &tile.chr))?;

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

    // Initialize counters (default -> 0s)
    let mut counter = FragmentKmersCounters::default();

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
            window_ctx.windows_slice(),
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
            let original_idx = window_ctx.original_idx(overlapped_window.idx);
            let counts = counts_by_window
                .entry(original_idx)
                .or_insert_with(FxHashMap::default);
            count_kmers_in_segments_clipped(
                &fragment,
                &positional_codes_by_k,
                kmer_specs,
                counts,
                positional_scaling_weights.as_deref(),
                tile.core_start,
                tile.core_end,
            );
        }
    }

    // Get counters from iterator
    counter.add_from_snapshot(iter.counters_snapshot());

    let mut payload: Vec<TileWindowCounts> = counts_by_window
        .into_iter()
        .filter_map(|(original_idx, hm)| {
            if hm.is_empty() {
                return None;
            }
            let mut entries: Vec<TileKmerCountEntry> = Vec::with_capacity(hm.len());
            for (kmer, value) in hm {
                entries.push(TileKmerCountEntry {
                    k: kmer.k,
                    code: kmer.code,
                    value,
                });
            }
            Some(TileWindowCounts {
                original_idx,
                entries,
            })
        })
        .collect();
    payload.sort_unstable_by_key(|w| w.original_idx);

    serialize_tile_counts(counts_path, &payload)?;

    Ok(TileResult {
        chr: tile.chr.clone(),
        counts_path: Some(counts_path.to_path_buf()),
        counter,
    })
}

/// Count kmers within the fragment’s usable segments, respecting tile core boundaries.
fn count_kmers_in_segments_clipped(
    fragment: &FragmentWithKmerSegments,
    positional_codes_by_k: &FxHashMap<u8, KmerCodes>,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
    counts: &mut FxHashMap<Kmer, f64>,
    weights: Option<&[f32]>,
    tile_core_start: u32,
    tile_core_end: u32,
) {
    for (&k, _) in kmer_specs {
        let codes = positional_codes_by_k
            .get(&k)
            .expect("missing positional codes for requested k");
        let k_span = k as u32;

        for &(seg_start, seg_end) in &fragment.segments {
            let seg_start = seg_start.max(tile_core_start);
            let seg_end = seg_end.min(tile_core_end);
            if seg_start >= seg_end {
                continue;
            }

            let Some(last_start) = seg_end.checked_sub(k_span) else {
                continue;
            };
            if last_start < seg_start {
                continue;
            }

            for idx_abs in seg_start..=last_start {
                let idx_local = (idx_abs - tile_core_start) as usize;
                let w = weights.map_or(1.0, |weights| unsafe { *weights.get_unchecked(idx_local) });

                *counts
                    .entry(Kmer {
                        k,
                        code: codes.get(idx_local),
                    })
                    .or_insert(0.) += w as f64;
            }
        }
    }
}
