//! Runner for WPS peak calling from BAM file.
//!
//! The intended logic is specified in the `peak_calling_logic.md` document.

use crate::commands::cli_common::{
    ensure_output_dir, load_blacklist_map, load_scaling_map, resolve_chromosomes_and_contigs,
    WindowSpec,
};
use crate::commands::counters::FCoverageCounters;
use crate::commands::wps::wps::wps_for_tile;
use crate::commands::wps_peaks::call_peaks::{call_peaks, PeakCall};
use crate::commands::wps_peaks::config::WPSPeaksConfig;
use crate::commands::wps_peaks::normalize_wps::{normalize_wps, smoothe_wps};
use crate::commands::wps_peaks::window_peak_results::PeaksWindowAction;
use crate::shared::bam::Contigs;
use crate::shared::bed::load_windows_from_bed;
use crate::shared::thread_pool::init_global_pool;
use crate::shared::tiled_run::{
    build_tiles, make_temp_dir, precompute_tile_window_spans, Tile, TileMode, TileWindowSpan,
};
use crate::shared::writers::open_zstd_auto_writer;
use anyhow::{anyhow, ensure, Context, Result};
use fxhash::FxHashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::fs::{remove_dir_all, File};
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
    let (chromosomes, contigs) =
        resolve_chromosomes_and_contigs(&opt.shared_args.chromosomes, &opt.shared_args.ioc)?;
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
            let wds = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None)?;
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
        // TODO: If per-window is None, there's no windows, and then tile_and_window_boundaries_align is irrelevant?
        None => {
            let mut writer = GlobalWriter::new(
                opt.shared_args
                    .ioc
                    .output_dir
                    .join(format!("{prefix}.wps.peaks.tsv.zst")),
                opt.shared_args.decimals as usize,
                opt.shared_args.ioc.n_threads as u32,
            )?;
            if tile_and_window_boundaries_align {
                for result in &tile_results {
                    total_counter += result.counter;
                    writer.write_tile_file(result.peak_file_path.as_path())?;
                }
            } else {
                for result in &tile_results {
                    total_counter += result.counter;
                    stream_tile_peaks(result.peak_file_path.as_path(), |peak| {
                        writer.write_peak(&peak)
                    })?;
                }
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
                    WindowSource::FixedSize(FixedSizeWindows::new(*bp as u64, &contigs))
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
    decimals: usize,
}

impl GlobalWriter {
    fn new(path: PathBuf, decimals: usize, threads: u32) -> Result<Self> {
        let mut writer = open_zstd_auto_writer(&path, 3, Some(threads))?;
        writeln!(writer, "chromosome\tstart\tend\tpeak_position\theight")?;
        Ok(Self {
            path,
            writer,
            decimals,
        })
    }

    fn write_peak(&mut self, peak: &PeakCall) -> Result<()> {
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}",
            peak.chromosome,
            peak.start,
            peak.end,
            peak.peak_position,
            format_float(peak.height, self.decimals)
        )?;
        Ok(())
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
    FixedSize(FixedSizeWindows),
}

struct FixedSizeWindows {
    size: u64,
    chrom_lengths: FxHashMap<String, u64>,
    progress: FxHashMap<String, FixedChromProgress>,
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
        })
    }

    fn process_tile(
        &mut self,
        tile: &Tile,
        path: &Path,
        contributions: Option<&[WindowStatsContribution]>,
    ) -> Result<()> {
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
            WindowSource::FixedSize(fixed) => {
                fixed.add_windows_for_tile(
                    &tile.chr,
                    &mut self.accumulator,
                    tile.core_start as u64,
                    tile.core_end as u64,
                )?;
            }
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
            if let WindowSource::FixedSize(fixed) = &mut self.window_source {
                fixed.ensure_progress(tile.chr.as_str());
            }
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

    let smoothed = if opt.no_smoothing {
        wps_values.clone()
    } else {
        smoothe_wps(&wps_values, Some(&mask))
    };

    let normalized = normalize_wps(
        &smoothed,
        &wps_values,
        Some(&mask), // TODO: will never be None due to the dilation?
        opt.normalize_bp as usize,
        1,
        opt.min_unmasked as usize,
    );

    let half_window = opt.shared_args.window_size / 2;
    let left_span = half_window.saturating_add(extra_halo);
    let dilated_start = tile.core_start.saturating_sub(left_span) as u64;

    let mut peaks = Vec::new();
    let peaks_all = call_peaks(
        &tile.chr,
        dilated_start,
        &normalized,
        &mask,
        min_peak_height,
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

    // TODO: How do we handle distances between peaks of neighbouring tiles? They should be considered in stats as well

    // TODO: We need to handle by-size windows

    windows
        .iter()
        .filter_map(|&(start, end, idx)| {
            // TODO: Unacceptable not to use index pointers for this
            let start_idx = peaks.partition_point(|peak| peak.peak_position < start);
            let end_idx = peaks.partition_point(|peak| peak.peak_position < end);
            let slice = &peaks[start_idx..end_idx];
            if slice.is_empty() {
                return None;
            }
            // TODO: Why btreemap instead of fxhashmap? No hashing?
            let mut histogram = BTreeMap::new();
            let mut distance_sum = 0.0;
            let mut prev: Option<u64> = None;
            for peak in slice {
                if let Some(previous) = prev {
                    let distance = (peak.peak_position - previous) as u32;
                    distance_sum += distance as f64;
                    let _ = histogram
                        .entry(distance)
                        .and_modify(|freq| *freq += 1)
                        .or_insert(1);
                }
                prev = Some(peak.peak_position);
            }
            Some(WindowStatsContribution {
                window_idx: idx,
                count: slice.len() as u32,
                first_peak: slice.first().map(|p| p.peak_position),
                last_peak: slice.last().map(|p| p.peak_position),
                distance_sum,
                distance_histogram: histogram,
            })
        })
        .collect()
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

enum WindowStateData {
    Unique(BTreeMap<u64, f32>),
    Indexed(Vec<PeakCall>),
    Stats {
        count: u32,
        first_peak: Option<u64>,
        last_peak: Option<u64>,
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
                    last_peak: None,
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
                        distance_sum,
                        distance_histogram,
                    } => {
                        if first_peak.is_none() {
                            *first_peak = Some(peak.peak_position);
                        }
                        if let Some(last) = *last_peak {
                            let distance = (peak.peak_position - last) as u32;
                            *distance_sum += distance as f64;
                            let _ = distance_histogram
                                .entry(distance)
                                .and_modify(|freq| *freq += 1)
                                .or_insert(1);
                        }
                        *count += 1;
                        *last_peak = Some(peak.peak_position);
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
                last_peak,
                distance_sum,
                distance_histogram,
            } => {
                if let (Some(prev_last), Some(first)) = (*last_peak, contribution.first_peak) {
                    if first >= prev_last {
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
                }
                *count += contribution.count;
                *distance_sum += contribution.distance_sum;
                for (distance, freq) in &contribution.distance_histogram {
                    let _ = distance_histogram
                        .entry(*distance)
                        .and_modify(|existing| *existing += *freq)
                        .or_insert(*freq);
                }
                if let Some(last) = contribution.last_peak.or(contribution.first_peak) {
                    *last_peak = Some(last);
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
                let total_distances: u32 = distance_histogram.values().copied().sum();
                let avg = if total_distances == 0 {
                    f32::NAN
                } else {
                    (distance_sum / total_distances as f64) as f32
                };
                let median = histogram_median(distance_histogram);
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

// TODO: Make this function readable with proper variable names and comments
pub fn histogram_median(hist: &BTreeMap<u32, u32>) -> f32 {
    let total: u32 = hist.values().copied().sum();
    if total == 0 {
        return f32::NAN;
    }
    let target1 = (total + 1) / 2;
    let target2 = (total + 2) / 2;
    let mut cumulative = 0u32;
    // TODO: m1/m2 are terrible variable names. Fix
    let mut m1: Option<u32> = None;
    let mut m2: Option<u32> = None;
    for (distance, freq) in hist {
        cumulative += *freq;
        if cumulative >= target1 && m1.is_none() {
            m1 = Some(*distance);
        }
        if cumulative >= target2 {
            m2 = Some(*distance);
            break;
        }
    }
    match (m1, m2) {
        (Some(d1), Some(d2)) => (d1 as f32 + d2 as f32) * 0.5,
        (Some(d), None) | (None, Some(d)) => d as f32,
        _ => f32::NAN,
    }
}
