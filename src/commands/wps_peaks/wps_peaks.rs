//! Runner for WPS peak calling from BAM file.
//!
//! The intended logic is specified in the `peak_calling_logic.md` document.

use crate::commands::cli_common::{
    WindowSpec, ensure_output_dir, load_blacklist_map, load_scaling_map,
    resolve_chromosomes_and_contigs,
};
use crate::commands::counters::FCoverageCounters;
use crate::commands::gc_bias::correct::{GCCorrector, load_gc_corrector};
use crate::commands::wps::wps::wps_for_tile;
use crate::commands::wps_peaks::call_peaks::{PeakCall, call_peaks};
use crate::commands::wps_peaks::config::WPSPeaksConfig;
use crate::commands::wps_peaks::normalize_wps::{normalize_wps, smoothe_wps};
use crate::commands::wps_peaks::window_peak_results::PeaksWindowAction;
use crate::shared::bam::Contigs;
use crate::shared::bed::load_windows_from_bed;
use crate::shared::thread_pool::init_global_pool;
use crate::shared::tiled_run::{
    Tile, TileMode, TileWindowSpan, build_tiles, make_temp_dir, precompute_tile_window_spans,
};
use crate::shared::writers::open_zstd_auto_writer;
use anyhow::{Context, Result, anyhow, ensure};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::BTreeMap;
use std::fs::{File, remove_dir_all};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

const EXTRA_PEAK_HALO_BP: u32 = 450;

/// Execute the Snyder-style peak calling pipeline on top of windowed protection scores.
///
/// This command shares most of the setup with `cfdna wps`: we resolve chromosomes,
/// prepare blacklist and scaling lookups, then iterate tiles to compute WPS values.
/// The positional WPS stay in memory per tile so we can smooth, normalize, and call peaks
/// without writing intermediate files.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `wps-peaks` command.
///
/// Returns
/// -------
/// - `Result<()>`:
///     Indicates whether peak calling finished and all outputs were written successfully.
pub fn run(opt: &WPSPeaksConfig) -> Result<()> {
    let start_time = Instant::now();
    let (chromosomes, contigs) = resolve_chromosomes_and_contigs(
        &opt.shared_args.chromosomes,
        &opt.shared_args.ioc.bam.as_path(),
    )?;
    let prefix = opt.shared_args.output_prefix.trim();
    let window_opt = opt.shared_args.windows.resolve_windows();
    let windowed = matches!(window_opt, WindowSpec::Bed(_) | WindowSpec::Size(_));

    if windowed {
        ensure!(
            opt.per_window.is_some(),
            "when using --by-bed/--by-size, please also specify --per-window"
        );
    }

    ensure!(
        opt.shared_args.min_fragment_length >= opt.shared_args.window_size,
        "min-fragment-length ({}) must be >= window-size ({})",
        opt.shared_args.min_fragment_length,
        opt.shared_args.window_size
    );
    ensure!(
        opt.shared_args.window_size <= opt.shared_args.max_fragment_length,
        "window-size ({}) must be <= max-fragment-length ({})",
        opt.shared_args.window_size,
        opt.shared_args.max_fragment_length
    );

    // Create output directory if needed
    ensure_output_dir(&opt.shared_args.ioc.output_dir)?;

    if opt.shared_args.blacklist.is_some() {
        println!("Start: Loading blacklists");
    }
    // Dilate blacklists so fragments that could reach them do not affect the WPS baseline
    let blacklist_halo =
        (opt.shared_args.max_fragment_length + (opt.shared_args.window_size + 1) / 2) as u64;
    let blacklist_map = load_blacklist_map(
        opt.shared_args.blacklist.as_ref(),
        1,
        blacklist_halo,
        &chromosomes,
    )?;

    let fixed_window_bp = if let WindowSpec::Size(bp) = window_opt {
        Some(bp as u64)
    } else {
        None
    };

    // Load windows from BED file
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            println!("Start: Loading window coordinates");
            let wds = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            if matches!(
                opt.per_window,
                Some(PeaksWindowAction::OnlyIncludeThesePositionsUnique)
            ) {
                // Merge in-place to avoid double memory-usage
                println!("Start: Merging overlapping/touching windows");
                // Take ownership so we can remove entries by chromosome
                let mut wds_owned: FxHashMap<String, crate::shared::bed::Windows> = wds;
                let mut out: FxHashMap<String, crate::shared::bed::Windows> =
                    FxHashMap::with_capacity_and_hasher(wds_owned.len(), Default::default());
                let mut next_idx: u64 = 0;

                // Use the user-provided `chromosomes` order to assign indices deterministically
                for chr in &chromosomes {
                    if let Some(ws) = wds_owned.remove(chr) {
                        // Flatten in-place
                        let (flat, next) = ws.into_flattened_reindexed(next_idx);
                        next_idx = next;
                        out.insert(chr.clone(), flat);
                    }
                }
                Some(out)
            } else {
                Some(wds)
            }
        }
        _ => None,
    };

    // Load genomic scaling factors
    if opt.shared_args.scale_genome.scaling_factors.is_some() {
        println!("Start: Loading scaling factors");
    }
    let scaling_map: FxHashMap<String, Vec<(u64, u64, f32)>> =
        load_scaling_map(&opt.shared_args.scale_genome, &chromosomes, &contigs)?;

    // Load GC correction package if specified
    if opt.shared_args.gc.gc_file.is_some() {
        println!("Start: Loading GC correction matrix");
    }
    let gc_corrector = load_gc_corrector(
        opt.shared_args.gc.gc_file.as_ref(),
        opt.shared_args.min_fragment_length,
        opt.shared_args.max_fragment_length,
    )?;

    // Halo large enough for normalization and to keep peaks that straddle tile edges
    let normalize_halo = opt.normalize_bp / 2;
    let extra_halo = normalize_halo.saturating_add(EXTRA_PEAK_HALO_BP);
    let base_halo = opt
        .shared_args
        .max_fragment_length
        .saturating_add(opt.shared_args.window_size);
    let halo_bp = base_halo.saturating_add(extra_halo);

    let (tiles, tile_and_window_boundaries_align) = build_tiles(
        &chromosomes,
        &contigs,
        opt.shared_args.tile_size,
        halo_bp,
        match window_opt {
            WindowSpec::Size(bp) => Some(bp as u64),
            _ => None,
        },
    )?;

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(&tiles, |chr| {
        windows_lookup
            .and_then(|m| m.get(chr).map(|w| w.as_slice()))
            .unwrap_or(&[])
    }));

    let total_tiles = tiles.len();
    let pb = Arc::new(ProgressBar::new(total_tiles as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("       {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    println!("Start: Calling peaks per tile");

    // Configure global thread‐pool size
    init_global_pool(opt.shared_args.ioc.n_threads as usize)?;

    let temp_dir = make_temp_dir(&opt.shared_args.ioc.output_dir, prefix)
        .context("create per-run temp dir for peaks")?;
    let tile_window_spans_for_threads = tile_window_spans.clone();
    let stats_mode = matches!(opt.per_window, Some(PeaksWindowAction::Stats));

    let tile_results: Vec<TileResult> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| -> Result<TileResult> {
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let windows_chr: Option<&[(u64, u64, u64)]> = windows_map
                .as_ref()
                .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice()));
            let blacklist_chr = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let scaling_chr = scaling_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);

            let (counter, peaks) = peaks_for_tile(
                opt,
                tile,
                tile_span.as_ref(),
                windows_chr,
                blacklist_chr,
                scaling_chr,
                gc_corrector.clone(),
                extra_halo,
                opt.min_peak_height,
            )?;

            let stats_contributions = if stats_mode {
                if let Some(bin_size) = fixed_window_bp {
                    let chrom_len = contigs
                        .contigs
                        .get(&tile.chr)
                        .map(|(_, len)| *len as u64)
                        .ok_or_else(|| anyhow!("missing contig info for {}", tile.chr))?;
                    let windows = build_fixed_size_windows_for_tile(
                        bin_size,
                        chrom_len,
                        tile.core_start as u64,
                        tile.core_end as u64,
                    );
                    Some(compute_window_stats_contributions(
                        windows.as_slice(),
                        &peaks,
                    ))
                } else if let (Some(span), Some(windows_chr_slice)) =
                    (tile_span.as_ref(), windows_chr)
                {
                    let window_slice = &windows_chr_slice[span.first_idx..span.last_idx_exclusive];
                    Some(compute_window_stats_contributions(window_slice, &peaks))
                } else {
                    Some(Vec::new())
                }
            } else {
                None
            };

            let path = tile_peaks_path(&temp_dir, prefix, tile, tile_idx);
            if peaks.is_empty() {
                File::create(&path).context("create empty tile peak file")?;
            } else {
                // TODO: We should probably use the existing writers for this?
                let mut writer =
                    BufWriter::new(File::create(&path).context("create tile peak file")?);
                for peak in &peaks {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}",
                        peak.chromosome,
                        peak.start,
                        peak.end,
                        peak.peak_position,
                        format_float(peak.height, opt.shared_args.decimals as usize)
                    )?;
                }
                writer.flush()?;
            }

            pb.inc(1);
            Ok(TileResult {
                counter,
                peak_file_path: path,
                stats: stats_contributions,
            })
        })
        .collect::<Result<_>>()?;

    pb.finish_with_message("| Finished peak calling");

    let mut total_counter = FCoverageCounters::default();
    match opt.per_window {
        None => {
            let mut writer = GlobalWriter::new(
                opt.shared_args
                    .ioc
                    .output_dir
                    .join(format!("{prefix}.wps.peaks.tsv.zst")),
                opt.shared_args.ioc.n_threads as u32,
            )?;
            for result in &tile_results {
                total_counter += result.counter;
                writer.write_tile_file(result.peak_file_path.as_path())?;
            }
            writer.finish()?;
        }
        Some(action) => {
            let window_source = match &window_opt {
                WindowSpec::Bed(_) => {
                    let windows_src = windows_map
                        .as_ref()
                        .ok_or_else(|| anyhow!("window map required for --by-bed outputs"))?;
                    WindowSource::Bed(
                        windows_src
                            .iter()
                            .map(|(chr, ws)| (chr.clone(), ws.as_slice().to_vec()))
                            .collect(),
                    )
                }
                WindowSpec::Size(bp) => {
                    if tile_and_window_boundaries_align {
                        WindowSource::FixedSizeAligned(Arc::new(FixedSizeAlignedWindows::new(
                            *bp as u64, &contigs,
                        )))
                    } else {
                        WindowSource::FixedSizeBuffered(FixedSizeWindows::new(*bp as u64, &contigs))
                    }
                }
                WindowSpec::Global => {
                    anyhow::bail!("per-window outputs require --by-bed or --by-size");
                }
            };
            let mut writer = WindowOutputWriter::new(
                &opt.shared_args.ioc.output_dir,
                prefix,
                opt.shared_args.decimals as usize,
                opt.shared_args.ioc.n_threads as u32,
                action,
                window_source,
            )?;
            for (tile, result) in tiles.iter().zip(tile_results.iter()) {
                total_counter += result.counter;
                writer.process_tile(
                    tile,
                    result.peak_file_path.as_path(),
                    result.stats.as_deref(),
                )?;
            }
            writer.finish()?;
        }
    }

    if let Err(e) = remove_dir_all(&temp_dir) {
        eprintln!(
            "warning: failed to remove temp dir {}: {}",
            temp_dir.display(),
            e
        );
    }

    println!();
    println!("Statistics");
    println!("----------");
    println!(
        "  Note: A few reads/fragments may be counted twice in the statistics (only) around the parallelization tiles."
    );
    let elapsed = start_time.elapsed();
    println!("  Total reads: {}", total_counter.base.total_reads);
    println!(
        "  Initially accepted reads: {} ({:.2}%, forward: {}, reverse: {})",
        total_counter.base.accepted_forward + total_counter.base.accepted_reverse,
        (total_counter.base.accepted_forward + total_counter.base.accepted_reverse) as f64
            / total_counter.base.total_reads as f64
            * 100.0,
        total_counter.base.accepted_forward,
        total_counter.base.accepted_reverse
    );
    if opt.shared_args.gc.gc_file.is_some() {
        println!(
            "  GC correction failures (fragment counted with weight 1.0): {}",
            total_counter.gc_failed_fragments
        );
    }
    println!(
        "  Fragments counted one or more times: {}",
        total_counter.base.counted_fragments
    );
    println!("----------");
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

pub struct GlobalWriter {
    path: PathBuf,
    writer: BufWriter<Box<dyn Write>>,
}

impl GlobalWriter {
    fn new(path: PathBuf, threads: u32) -> Result<Self> {
        let mut writer = open_zstd_auto_writer(&path, 3, Some(threads))?;
        writeln!(writer, "chromosome\tstart\tend\tpeak_position\theight")?;
        Ok(Self { path, writer })
    }

    fn write_tile_file(&mut self, path: &Path) -> Result<()> {
        let mut reader = File::open(path).context("open tile peak file for copy")?;
        io::copy(&mut reader, &mut self.writer).context("copy tile peaks into output")?;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.writer.flush()?;
        println!("Saved peaks to: {}", self.path.display());
        Ok(())
    }
}

pub struct TileResult {
    pub counter: FCoverageCounters,
    pub peak_file_path: PathBuf,
    pub stats: Option<Vec<WindowStatsContribution>>,
}

// TODO: Document all arguments
pub struct WindowStatsContribution {
    pub window_idx: u64,
    pub count: u32,
    pub first_peak: Option<u64>,
    pub last_peak: Option<u64>,
    pub first_segment: Option<u64>,
    pub last_segment: Option<u64>,
    pub distance_sum: f64,
    pub distance_histogram: BTreeMap<u32, u32>,
}

pub fn tile_peaks_path(root: &Path, prefix: &str, tile: &Tile, tile_idx: usize) -> PathBuf {
    let safe_chr = tile.chr.replace('/', "_");
    root.join(format!(
        "{prefix}.tile.{safe_chr}.{tile_idx:05}.peaks.tmp",
        prefix = prefix,
        safe_chr = safe_chr,
        tile_idx = tile_idx
    ))
}

pub fn stream_tile_peaks<F>(path: &Path, mut handler: F) -> Result<()>
where
    F: FnMut(PeakCall) -> Result<()>,
{
    let file = File::open(path).context("open tile peak file")?;
    let reader = BufReader::new(file);
    for (line_idx, line_res) in reader.lines().enumerate() {
        let line = line_res.with_context(|| format!("read peak line {}", line_idx + 1))?;
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let chr = fields
            .next()
            .ok_or_else(|| anyhow!("missing chromosome in peak line"))?;
        let start: u64 = fields
            .next()
            .ok_or_else(|| anyhow!("missing start in peak line"))?
            .parse()
            .context("parse peak start")?;
        let end: u64 = fields
            .next()
            .ok_or_else(|| anyhow!("missing end in peak line"))?
            .parse()
            .context("parse peak end")?;
        let peak_position: u64 = fields
            .next()
            .ok_or_else(|| anyhow!("missing peak position"))?
            .parse()
            .context("parse peak position")?;
        let height: f32 = fields
            .next()
            .ok_or_else(|| anyhow!("missing peak height"))?
            .parse()
            .context("parse peak height")?;
        handler(PeakCall {
            chromosome: chr.to_string(),
            start,
            end,
            peak_position,
            height,
            segment_id: 0,
        })?;
    }
    Ok(())
}

enum WindowOutputMode {
    Unique,
    Indexed,
    Stats,
}

enum WindowSource {
    Bed(FxHashMap<String, Vec<(u64, u64, u64)>>),
    FixedSizeBuffered(FixedSizeWindows),
    FixedSizeAligned(Arc<FixedSizeAlignedWindows>),
}

struct FixedSizeWindows {
    size: u64,
    chrom_lengths: FxHashMap<String, u64>,
    progress: FxHashMap<String, FixedChromProgress>,
}

struct FixedSizeAlignedWindows {
    size: u64,
    chrom_lengths: FxHashMap<String, u64>,
}

#[derive(Default)]
struct FixedChromProgress {
    next_start: u64,
    next_idx: u64,
}

impl FixedSizeWindows {
    fn new(size: u64, contigs: &Contigs) -> Self {
        let mut chrom_lengths = FxHashMap::default();
        for (chr, (_, len)) in contigs.contigs.iter() {
            chrom_lengths.insert(chr.clone(), *len as u64);
        }
        Self {
            size,
            chrom_lengths,
            progress: FxHashMap::default(),
        }
    }

    fn ensure_progress(&mut self, chr: &str) {
        self.progress
            .entry(chr.to_string())
            .or_insert_with(FixedChromProgress::default);
    }

    fn add_windows_for_tile(
        &mut self,
        chr: &str,
        accumulator: &mut WindowAccumulator,
        tile_start: u64,
        tile_end: u64,
    ) -> Result<()> {
        let chrom_len = *self
            .chrom_lengths
            .get(chr)
            .ok_or_else(|| anyhow!("missing contig length for {}", chr))?;
        self.ensure_progress(chr);
        let state = self.progress.get_mut(chr).expect("progress initialized");

        while state.next_start < tile_end && state.next_start < chrom_len {
            let window_start = state.next_start;
            let window_end = (window_start + self.size).min(chrom_len);
            state.next_start = window_end;
            let idx = state.next_idx;
            state.next_idx += 1;

            if window_end <= tile_start {
                continue;
            }

            accumulator.push_window((window_start, window_end, idx));
        }

        Ok(())
    }
}

impl FixedSizeAlignedWindows {
    fn new(size: u64, contigs: &Contigs) -> Self {
        let mut chrom_lengths = FxHashMap::default();
        for (chr, (_, len)) in contigs.contigs.iter() {
            chrom_lengths.insert(chr.clone(), *len as u64);
        }
        Self {
            size,
            chrom_lengths,
        }
    }

    fn chrom_len(&self, chr: &str) -> Result<u64> {
        self.chrom_lengths
            .get(chr)
            .copied()
            .ok_or_else(|| anyhow!("missing contig length for {}", chr))
    }
}

impl From<PeaksWindowAction> for WindowOutputMode {
    fn from(action: PeaksWindowAction) -> Self {
        match action {
            PeaksWindowAction::OnlyIncludeThesePositionsUnique => Self::Unique,
            PeaksWindowAction::OnlyIncludeThesePositionsIndexed => Self::Indexed,
            PeaksWindowAction::Stats => Self::Stats,
        }
    }
}

// TODO: document this with pedagogical explanations
pub struct WindowOutputWriter {
    path: PathBuf,
    writer: BufWriter<Box<dyn Write>>,
    accumulator: WindowAccumulator,
    window_source: WindowSource,
    current_chr: Option<String>,
    next_idx: usize,
    mode: WindowOutputMode,
    decimals: usize,
}

impl WindowOutputWriter {
    fn new(
        output_dir: &Path,
        prefix: &str,
        decimals: usize,
        threads: u32,
        action: PeaksWindowAction,
        window_source: WindowSource,
    ) -> Result<Self> {
        let mode = WindowOutputMode::from(action);
        let path = output_dir.join(match mode {
            WindowOutputMode::Unique => format!("{prefix}.wps.peaks.unique.tsv.zst"),
            WindowOutputMode::Indexed => format!("{prefix}.wps.peaks.indexed.tsv.zst"),
            WindowOutputMode::Stats => format!("{prefix}.wps.peaks.stats.tsv.zst"),
        });
        let mut writer = open_zstd_auto_writer(&path, 3, Some(threads))?;
        match mode {
            WindowOutputMode::Unique => {
                writeln!(writer, "chromosome\tstart\tend\theight")?;
            }
            WindowOutputMode::Indexed => {
                writeln!(
                    writer,
                    "chromosome\tstart\tend\tpeak_position\theight\twindow_index"
                )?;
            }
            // TODO: We need to also add standard deviation of distances
            WindowOutputMode::Stats => {
                writeln!(
                    writer,
                    "chromosome\tstart\tend\twindow_index\tcount\tavg_distance\tmedian_distance"
                )?;
            }
        }

        Ok(Self {
            path,
            writer,
            accumulator: WindowAccumulator::new(action, decimals),
            window_source,
            current_chr: None,
            next_idx: 0,
            mode,
            decimals,
        })
    }

    fn process_tile(
        &mut self,
        tile: &Tile,
        path: &Path,
        contributions: Option<&[WindowStatsContribution]>,
    ) -> Result<()> {
        let aligned_source = if let WindowSource::FixedSizeAligned(aligned) = &self.window_source {
            Some(Arc::clone(aligned))
        } else {
            None
        };
        if let Some(aligned) = aligned_source {
            self.process_aligned_fixed_size_tile(tile, path, contributions, aligned.as_ref())?;
            return Ok(());
        }

        self.ensure_chromosome(tile)?;

        match &mut self.window_source {
            WindowSource::Bed(map) => {
                let windows_chr = map.get(&tile.chr).map(|v| v.as_slice()).unwrap_or(&[]);
                self.accumulator.add_windows_for_tile(
                    windows_chr,
                    &mut self.next_idx,
                    tile.core_start as u64,
                    tile.core_end as u64,
                );
            }
            WindowSource::FixedSizeBuffered(fixed) => {
                fixed.add_windows_for_tile(
                    &tile.chr,
                    &mut self.accumulator,
                    tile.core_start as u64,
                    tile.core_end as u64,
                )?;
            }
            WindowSource::FixedSizeAligned(_) => unreachable!("aligned windows handled earlier"),
        }

        match self.mode {
            WindowOutputMode::Unique | WindowOutputMode::Indexed => {
                stream_tile_peaks(path, |peak| {
                    self.accumulator.push_peak(&peak);
                    Ok(())
                })?;
            }
            WindowOutputMode::Stats => {
                if let Some(contribs) = contributions {
                    for contribution in contribs {
                        self.accumulator.apply_stats_contribution(contribution)?;
                    }
                } else {
                    stream_tile_peaks(path, |peak| {
                        self.accumulator.push_peak(&peak);
                        Ok(())
                    })?;
                }
            }
        }

        self.accumulator
            .flush_completed_windows(tile.core_end as u64, &mut self.writer)?;
        Ok(())
    }

    fn process_aligned_fixed_size_tile(
        &mut self,
        tile: &Tile,
        path: &Path,
        contributions: Option<&[WindowStatsContribution]>,
        fixed: &FixedSizeAlignedWindows,
    ) -> Result<()> {
        if self.current_chr.as_deref() != Some(tile.chr.as_str()) {
            self.current_chr = Some(tile.chr.clone());
        }

        match self.mode {
            WindowOutputMode::Unique => {
                self.write_aligned_unique(path, tile.chr.as_str())?;
            }
            WindowOutputMode::Indexed => {
                self.write_aligned_indexed(path, fixed.size)?;
            }
            WindowOutputMode::Stats => {
                let chrom_len = fixed.chrom_len(tile.chr.as_str())?;
                let windows = build_fixed_size_windows_for_tile(
                    fixed.size,
                    chrom_len,
                    tile.core_start as u64,
                    tile.core_end as u64,
                );
                self.write_aligned_stats(
                    tile.chr.as_str(),
                    &windows,
                    contributions.unwrap_or(&[]),
                )?;
            }
        }
        Ok(())
    }

    fn write_aligned_unique(&mut self, path: &Path, chr: &str) -> Result<()> {
        let best_by_position = Self::collect_aligned_unique_peaks(path)?;
        for (pos, height) in best_by_position {
            writeln!(
                self.writer,
                "{}\t{}\t{}\t{}",
                chr,
                pos,
                pos + 1,
                format_float(height, self.decimals)
            )?;
        }
        Ok(())
    }

    fn write_aligned_indexed(&mut self, path: &Path, bin_size: u64) -> Result<()> {
        stream_tile_peaks(path, |peak| {
            let idx = peak.peak_position / bin_size;
            writeln!(
                self.writer,
                "{}\t{}\t{}\t{}\t{}\t{}",
                peak.chromosome,
                peak.start,
                peak.end,
                peak.peak_position,
                format_float(peak.height, self.decimals),
                idx
            )?;
            Ok(())
        })
    }

    fn write_aligned_stats(
        &mut self,
        chr: &str,
        windows: &[(u64, u64, u64)],
        contributions: &[WindowStatsContribution],
    ) -> Result<()> {
        let mut contrib_lookup: FxHashMap<u64, &WindowStatsContribution> =
            FxHashMap::with_capacity_and_hasher(contributions.len(), Default::default());
        for contribution in contributions {
            contrib_lookup.insert(contribution.window_idx, contribution);
        }

        for &(start, end, idx) in windows {
            if let Some(contribution) = contrib_lookup.get(&idx) {
                let (avg, median) = stats_distance_summary(
                    contribution.distance_sum,
                    &contribution.distance_histogram,
                );
                writeln!(
                    self.writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    chr,
                    start,
                    end,
                    idx,
                    contribution.count,
                    format_float(avg, self.decimals),
                    format_float(median, self.decimals)
                )?;
            } else {
                writeln!(
                    self.writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    chr,
                    start,
                    end,
                    idx,
                    0,
                    format_float(f32::NAN, self.decimals),
                    format_float(f32::NAN, self.decimals)
                )?;
            }
        }
        Ok(())
    }

    pub fn collect_aligned_unique_peaks(path: &Path) -> Result<BTreeMap<u64, f32>> {
        let mut best_by_position = BTreeMap::<u64, f32>::new();
        stream_tile_peaks(path, |peak| {
            let entry = best_by_position
                .entry(peak.peak_position)
                .or_insert(peak.height);
            if peak.height > *entry {
                *entry = peak.height;
            }
            Ok(())
        })?;
        Ok(best_by_position)
    }

    fn finish(&mut self) -> Result<()> {
        self.accumulator.flush_all(&mut self.writer)?;
        self.writer.flush()?;
        match self.mode {
            WindowOutputMode::Stats => {
                println!("Saved window stats to: {}", self.path.display());
            }
            _ => println!("Saved peaks to: {}", self.path.display()),
        }
        Ok(())
    }

    fn ensure_chromosome(&mut self, tile: &Tile) -> Result<()> {
        if self.current_chr.as_deref() != Some(tile.chr.as_str()) {
            if self.current_chr.is_some() {
                self.accumulator.flush_all(&mut self.writer)?;
            }
            self.current_chr = Some(tile.chr.clone());
            self.accumulator.reset_for_chromosome(tile.chr.clone());
            self.next_idx = 0;
            match &mut self.window_source {
                WindowSource::FixedSizeBuffered(fixed) => fixed.ensure_progress(tile.chr.as_str()),
                _ => {}
            };
        }
        Ok(())
    }
}

pub fn peaks_for_tile(
    opt: &WPSPeaksConfig,
    tile: &Tile,
    tile_span: Option<&TileWindowSpan>,
    windows: Option<&[(u64, u64, u64)]>,
    blacklist_chr: &[(u64, u64)],
    scaling_chr: &[(u64, u64, f32)],
    gc_corrector_opt: Option<GCCorrector>,
    extra_halo: u32,
    min_peak_height: f32,
) -> Result<(FCoverageCounters, Vec<PeakCall>)> {
    let dummy_path = PathBuf::new();
    let tile_mode = TileMode::Positional {
        windows,
        out_path: dummy_path,
        indexed: false,
    };

    let (counter, wps_opt, mask_opt) = wps_for_tile(
        &opt.shared_args,
        &None, // Ignored
        false, // Ignored
        tile,
        tile_span,
        blacklist_chr,
        scaling_chr,
        gc_corrector_opt,
        tile_mode,
        opt.shared_args.decimals as i32, // Ignored
        extra_halo,
        true,
    )?;

    let wps_values = match wps_opt {
        Some(v) => v,
        None => return Ok((counter, Vec::new())),
    };
    // TODO: Handle edges in mask being masked out? Is masked to the core or just to the parts that cannot have WPS scores?
    let mut mask = mask_opt.unwrap_or_else(|| vec![0; wps_values.len()]);
    if mask.len() != wps_values.len() {
        mask.resize(wps_values.len(), 0);
    }

    let half_window = opt.shared_args.window_size / 2;
    let left_span = half_window.saturating_add(extra_halo);
    let dilated_start = tile.core_start.saturating_sub(left_span) as u64;
    let initial_segment_marker = last_mask_end_before(blacklist_chr, dilated_start);

    let processing_opts = PeakSignalProcessingOptions {
        smoothing: !opt.no_smoothing,
        normalization_bp: Some(opt.normalize_bp as usize),
        min_unmasked: opt.min_unmasked as usize,
        min_peak_height,
        initial_segment_marker,
    };
    let mut peaks = Vec::new();
    let peaks_all = peaks_from_wps_values(
        &tile.chr,
        dilated_start,
        &wps_values,
        Some(&mask),
        &processing_opts,
    );

    for mut peak in peaks_all {
        if peak.peak_position >= tile.core_start as u64 && peak.peak_position < tile.core_end as u64
        {
            peak.start = peak.start.max(tile.core_start as u64);
            peak.end = peak.end.min(tile.core_end as u64);
            peaks.push(peak);
        }
    }
    peaks.sort_by_key(|p| p.start);
    Ok((counter, peaks))
}

pub fn compute_window_stats_contributions(
    windows: &[(u64, u64, u64)],
    peaks: &[PeakCall],
) -> Vec<WindowStatsContribution> {
    if windows.is_empty() || peaks.is_empty() {
        return Vec::new();
    }

    let mut contributions = Vec::new();
    let mut start_idx = 0usize;
    let mut end_idx;

    // Windows arrive sorted by start (they may overlap), so `start_idx` can increase monotonically
    for &(start, end, idx) in windows {
        while start_idx < peaks.len() && peaks[start_idx].peak_position < start {
            start_idx += 1;
        }
        end_idx = start_idx.clone();
        while end_idx < peaks.len() && peaks[end_idx].peak_position < end {
            end_idx += 1;
        }
        if start_idx == end_idx {
            continue;
        }
        let slice = &peaks[start_idx..end_idx];
        // BTreeMap keeps distances sorted deterministically so medians/formatting stay stable
        let mut histogram = BTreeMap::new();
        let mut distance_sum = 0.0;
        let mut prev: Option<(u64, u64)> = None;
        for peak in slice {
            if let Some((previous, previous_segment)) = prev {
                if previous_segment == peak.segment_id {
                    let distance = (peak.peak_position - previous) as u32;
                    distance_sum += distance as f64;
                    let _ = histogram
                        .entry(distance)
                        .and_modify(|freq| *freq += 1)
                        .or_insert(1);
                }
            }
            prev = Some((peak.peak_position, peak.segment_id));
        }
        contributions.push(WindowStatsContribution {
            window_idx: idx,
            count: slice.len() as u32,
            first_peak: slice.first().map(|p| p.peak_position),
            last_peak: slice.last().map(|p| p.peak_position),
            first_segment: slice.first().map(|p| p.segment_id),
            last_segment: slice.last().map(|p| p.segment_id),
            distance_sum,
            distance_histogram: histogram,
        });
    }

    contributions
}

fn build_fixed_size_windows_for_tile(
    bin_size: u64,
    chrom_len: u64,
    tile_start: u64,
    tile_end: u64,
) -> Vec<(u64, u64, u64)> {
    if bin_size == 0 || tile_start >= chrom_len {
        return Vec::new();
    }
    // Windows are aligned to multiples of bin_size (0, bin_size, 2*bin_size, ...).
    // We start from the alignment that covers tile_start so we never skip windows
    // that begin before the tile but extend into it; those will already have been
    // processed by previous tiles.
    let mut start = (tile_start / bin_size) * bin_size;
    let mut windows = Vec::new();
    while start < tile_end && start < chrom_len {
        let window_start = start;
        let end = (start + bin_size).min(chrom_len);
        let idx = window_start / bin_size;
        windows.push((window_start, end, idx));
        start = start.saturating_add(bin_size);
    }
    windows
}

fn last_mask_end_before(intervals: &[(u64, u64)], position: u64) -> u64 {
    intervals
        .iter()
        .rev()
        .find_map(|&(_, end)| {
            if end <= position {
                Some(end.saturating_sub(1))
            } else {
                None
            }
        })
        .unwrap_or(0)
}

/// Tracks windows currently spanned by the active tile and accumulates per-window outputs.
pub struct WindowAccumulator {
    kind: WindowAccumulatorKind,
    decimals: usize,
    current_chr: String,
    active: Vec<WindowState>,
    // Index of the first window that may still contain upcoming peaks
    scan_start: usize,
}

enum WindowAccumulatorKind {
    Unique,
    Indexed,
    Stats,
}

pub struct WindowState {
    entry: (u64, u64, u64),
    data: WindowStateData,
}

/// Configuration for running the smoothing/normalization + peak calling pipeline on an
/// in-memory WPS signal.
#[derive(Debug, Clone)]
pub struct PeakSignalProcessingOptions {
    /// Whether Savitzky-Golay smoothing should be applied before normalization.
    pub smoothing: bool,
    /// Size of the rolling-median window used for baseline subtraction. When `None`, the
    /// incoming values are treated as residuals.
    pub normalization_bp: Option<usize>,
    /// Minimum usable bases required inside the normalization window.
    pub min_unmasked: usize,
    /// Minimum residual height required to keep a peak.
    pub min_peak_height: f32,
    /// Segment marker inherited from upstream tiles (usually the end of the last masked region).
    pub initial_segment_marker: u64,
}

/// Process an in-memory WPS signal (with optional mask) and return Snyder-style peaks.
///
/// This helper mirrors the per-tile pipeline executed by `peaks_for_tile`, but skips all BAM IO.
/// It enables unit tests to exercise the smoothing + normalization + peak-calling logic using
/// synthetic fixtures.
pub fn peaks_from_wps_values(
    chromosome: &str,
    start_offset: u64,
    wps_values: &[f32],
    mask: Option<&[u8]>,
    options: &PeakSignalProcessingOptions,
) -> Vec<PeakCall> {
    if let Some(mask_slice) = mask {
        assert_eq!(
            mask_slice.len(),
            wps_values.len(),
            "mask length must match WPS series length"
        );
    }
    let mask_cow: Cow<[u8]> = match mask {
        Some(existing) => Cow::Borrowed(existing),
        None => Cow::Owned(vec![0u8; wps_values.len()]),
    };
    let mask_slice: &[u8] = mask_cow.as_ref();
    let smoothed = if options.smoothing {
        smoothe_wps(wps_values, Some(mask_slice))
    } else {
        wps_values.to_vec()
    };
    let residuals = if let Some(window_bp) = options.normalization_bp {
        normalize_wps(
            &smoothed,
            wps_values,
            Some(mask_slice),
            window_bp,
            1,
            options.min_unmasked,
        )
    } else {
        smoothed
    };
    call_peaks(
        chromosome,
        start_offset,
        &residuals,
        mask_slice,
        options.min_peak_height,
        options.initial_segment_marker,
    )
}

enum WindowStateData {
    Unique(BTreeMap<u64, f32>),
    Indexed(Vec<PeakCall>),
    Stats {
        count: u32,
        first_peak: Option<u64>,
        last_peak: Option<u64>,
        first_segment: Option<u64>,
        last_segment: Option<u64>,
        distance_sum: f64,
        distance_histogram: BTreeMap<u32, u32>,
    },
}

impl WindowAccumulator {
    pub fn new(action: PeaksWindowAction, decimals: usize) -> Self {
        let kind = match action {
            PeaksWindowAction::OnlyIncludeThesePositionsUnique => WindowAccumulatorKind::Unique,
            PeaksWindowAction::OnlyIncludeThesePositionsIndexed => WindowAccumulatorKind::Indexed,
            PeaksWindowAction::Stats => WindowAccumulatorKind::Stats,
        };
        Self {
            kind,
            decimals,
            current_chr: String::new(),
            active: Vec::new(),
            scan_start: 0,
        }
    }

    pub fn reset_for_chromosome(&mut self, chr: String) {
        self.current_chr = chr;
        self.active.clear();
        self.scan_start = 0;
    }

    /// Register a new window with the accumulator (used by both BED and fixed-size paths).
    ///
    /// Stats windows keep enough state (first/last peaks plus histogram) to merge contributions
    /// across tiles without replaying peaks on tile boundaries.
    fn push_window(&mut self, entry: (u64, u64, u64)) {
        self.active.push(WindowState {
            entry,
            data: match self.kind {
                WindowAccumulatorKind::Unique => WindowStateData::Unique(BTreeMap::new()),
                WindowAccumulatorKind::Indexed => WindowStateData::Indexed(Vec::new()),
                WindowAccumulatorKind::Stats => WindowStateData::Stats {
                    count: 0,
                    first_peak: None,
                    first_segment: None,
                    last_peak: None,
                    last_segment: None,
                    distance_sum: 0.0,
                    distance_histogram: BTreeMap::new(),
                },
            },
        });
    }

    pub fn add_windows_for_tile(
        &mut self,
        windows: &[(u64, u64, u64)],
        next_idx: &mut usize,
        tile_start: u64,
        tile_end: u64,
    ) {
        // Windows come sorted by start, so we append in order and skip ones that end before the tile
        while *next_idx < windows.len() && windows[*next_idx].0 < tile_end {
            let entry = windows[*next_idx];
            if entry.1 <= tile_start {
                *next_idx += 1;
                continue;
            }
            self.push_window(entry);
            *next_idx += 1;
        }
    }

    /// Feed a peak directly into the active windows.
    ///
    /// This is used whenever we do not have a pre-computed contribution (e.g., unique/indexed
    /// outputs or stats for the final partial window). For stats we store the first peak seen so
    /// that subsequent tiles can add one cross-tile distance when merging contributions.
    pub fn push_peak(&mut self, peak: &PeakCall) {
        while self.scan_start < self.active.len()
            && self.active[self.scan_start].entry.1 <= peak.peak_position
        {
            self.scan_start += 1;
        }

        let mut window_idx = self.scan_start;
        while window_idx < self.active.len() {
            let state = &mut self.active[window_idx];
            let (start, end, _) = state.entry;
            if start > peak.peak_position {
                break;
            }
            if peak.peak_position >= start && peak.peak_position < end {
                match &mut state.data {
                    WindowStateData::Unique(map) => {
                        map.entry(peak.peak_position)
                            .and_modify(|h| {
                                if peak.height > *h {
                                    *h = peak.height;
                                }
                            })
                            .or_insert(peak.height);
                    }
                    WindowStateData::Indexed(list) => {
                        list.push(peak.clone());
                    }
                    WindowStateData::Stats {
                        count,
                        first_peak,
                        last_peak,
                        first_segment,
                        last_segment,
                        distance_sum,
                        distance_histogram,
                    } => {
                        if first_peak.is_none() {
                            *first_peak = Some(peak.peak_position);
                            *first_segment = Some(peak.segment_id);
                        }
                        if let (Some(last), Some(last_seg)) = (*last_peak, *last_segment) {
                            if last_seg == peak.segment_id {
                                let distance = (peak.peak_position - last) as u32;
                                *distance_sum += distance as f64;
                                let _ = distance_histogram
                                    .entry(distance)
                                    .and_modify(|freq| *freq += 1)
                                    .or_insert(1);
                            }
                        }
                        *count += 1;
                        *last_peak = Some(peak.peak_position);
                        *last_segment = Some(peak.segment_id);
                    }
                }
            }
            window_idx += 1;
        }
    }

    pub fn flush_completed_windows<W: Write>(
        &mut self,
        tile_end: u64,
        writer: &mut W,
    ) -> Result<()> {
        let mut remove_count = 0;
        while remove_count < self.active.len() && self.active[remove_count].entry.1 <= tile_end {
            remove_count += 1;
        }
        if remove_count == 0 {
            return Ok(());
        }
        let drained: Vec<WindowState> = self.active.drain(0..remove_count).collect();
        for state in drained {
            self.write_window(writer, state)?;
        }
        if self.scan_start < remove_count {
            self.scan_start = 0;
        } else {
            self.scan_start -= remove_count;
        }
        Ok(())
    }

    pub fn flush_all<W: Write>(&mut self, writer: &mut W) -> Result<()> {
        let drained: Vec<WindowState> = self.active.drain(..).collect();
        for state in drained {
            self.write_window(writer, state)?;
        }
        self.scan_start = 0;
        Ok(())
    }

    /// Merge statistics computed per tile into the live window state.
    ///
    /// Each contribution already contains all intra-tile distances. We only need to stitch the gap
    /// between the previous tile's last peak and this contribution's first peak (if both exist),
    /// then accumulate the histograms. This produces the same result as streaming every peak in
    /// order but avoids replaying them when the per-tile contribution is available.
    pub fn apply_stats_contribution(
        &mut self,
        contribution: &WindowStatsContribution,
    ) -> Result<()> {
        if contribution.count == 0 {
            return Ok(());
        }
        let state = self
            .active
            .iter_mut()
            .find(|st| st.entry.2 == contribution.window_idx)
            .ok_or_else(|| {
                anyhow!(
                    "stats contribution references inactive window {}",
                    contribution.window_idx
                )
            })?;
        match &mut state.data {
            WindowStateData::Stats {
                count,
                first_peak,
                first_segment,
                last_peak,
                last_segment,
                distance_sum,
                distance_histogram,
            } => {
                if let (
                    Some(prev_last),
                    Some(prev_segment),
                    Some(first),
                    Some(first_segment_contrib),
                ) = (
                    *last_peak,
                    *last_segment,
                    contribution.first_peak,
                    contribution.first_segment,
                ) {
                    if prev_segment == first_segment_contrib && first >= prev_last {
                        let distance = first - prev_last;
                        *distance_sum += distance as f64;
                        let _ = distance_histogram
                            .entry(distance as u32)
                            .and_modify(|freq| *freq += 1)
                            .or_insert(1);
                    }
                }
                if first_peak.is_none() {
                    *first_peak = contribution.first_peak;
                    *first_segment = contribution.first_segment;
                }
                *count += contribution.count;
                *distance_sum += contribution.distance_sum;
                for (distance, freq) in &contribution.distance_histogram {
                    let _ = distance_histogram
                        .entry(*distance)
                        .and_modify(|existing| *existing += *freq)
                        .or_insert(*freq);
                }
                if contribution.last_peak.is_some() {
                    *last_peak = contribution.last_peak;
                    *last_segment = contribution.last_segment;
                }
                Ok(())
            }
            _ => Err(anyhow!("stats contribution applied to non-stats window")),
        }
    }

    pub fn write_window<W: Write>(&self, writer: &mut W, state: WindowState) -> Result<()> {
        let (start, end, idx) = state.entry;
        match state.data {
            WindowStateData::Unique(map) => {
                for (pos, height) in map {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}",
                        self.current_chr,
                        pos,
                        pos + 1,
                        format_float(height, self.decimals)
                    )?;
                }
            }
            WindowStateData::Indexed(list) => {
                for peak in list {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        peak.chromosome,
                        peak.start,
                        peak.end,
                        peak.peak_position,
                        format_float(peak.height, self.decimals),
                        idx
                    )?;
                }
            }
            WindowStateData::Stats {
                count,
                distance_sum,
                ref distance_histogram,
                ..
            } => {
                let (avg, median) = stats_distance_summary(distance_sum, distance_histogram);
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    self.current_chr,
                    start,
                    end,
                    idx,
                    count,
                    format_float(avg, self.decimals),
                    format_float(median, self.decimals)
                )?;
            }
        }
        Ok(())
    }
}

fn format_float(value: f32, decimals: usize) -> String {
    if value.is_nan() {
        "NaN".to_string()
    } else {
        format!("{:.*}", decimals, value)
    }
}

/// Compute the median distance recorded in a histogram.
///
/// The histogram is expected to store monotonically increasing distance keys.
/// We traverse the bins once, tracking cumulative rank until we reach the lower and upper median
/// targets (needed for even-sized samples), then average those values. Empty histograms produce
/// `NaN`, signaling the caller that no distances were recorded.
///
/// Parameters
/// ----------
/// - `hist`:
///     BTreeMap keyed by distance (in bp) with frequencies as values.
///
/// Returns
/// -------
/// - `f32`:
///     Median distance in base pairs, or `NaN` if the histogram has no entries.
pub fn histogram_median(hist: &BTreeMap<u32, u32>) -> f32 {
    let total_count: u32 = hist.values().copied().sum();
    if total_count == 0 {
        return f32::NAN;
    }

    // Even counts require the average of the middle pair, so capture both ranks.
    let lower_target_rank = (total_count + 1) / 2;
    let upper_target_rank = (total_count + 2) / 2;

    let mut cumulative_rank = 0u32;
    let mut lower_median_value: Option<u32> = None;
    let mut upper_median_value: Option<u32> = None;

    for (distance_bp, frequency) in hist {
        cumulative_rank += *frequency;
        if cumulative_rank >= lower_target_rank && lower_median_value.is_none() {
            lower_median_value = Some(*distance_bp);
        }
        if cumulative_rank >= upper_target_rank {
            upper_median_value = Some(*distance_bp);
            break;
        }
    }

    match (lower_median_value, upper_median_value) {
        (Some(lower), Some(upper)) => (lower as f32 + upper as f32) * 0.5,
        (Some(value), None) | (None, Some(value)) => value as f32,
        _ => f32::NAN,
    }
}

pub fn stats_distance_summary(distance_sum: f64, histogram: &BTreeMap<u32, u32>) -> (f32, f32) {
    let total_distances: u32 = histogram.values().copied().sum();
    if total_distances == 0 {
        (f32::NAN, f32::NAN)
    } else {
        (
            (distance_sum / total_distances as f64) as f32,
            histogram_median(histogram),
        )
    }
}
