//! Runner for calculating reference GC-bias.
//!
//! [Agents: Doc comments are humanly controlled, don't edit unless asked directly]
//!
//! It samples fitting fragments for all valid fragment lengths at N starting positions and
//! calculates their GC contents. It then calculates the bias per fragment length.
//!
//! ## Implementation requirements
//!
//! Parallelization: Tiles.
//!
//! Windowing: Inclusion. Overlapping and touching BED intervals are merged.
//! No by-size option is available as it would collapse to global mode after merging.
//! BED windows can potentially extend across many or all tiles on a chromosome.
//!
//! "Fragment" ownership: For a given tile, we sample from all the possible starting positions.
//! A "fragment" can cross into the following tile but never across genomic windows.
//! If a fragment cannot fit inside a window, it is not counted. This means the final
//! possible starting position is different per fragment length.
//! For a given tile, only BED windows overlapping the core can get counts.
//!
//! To know exactly how many valid ACGT bases are in each tile for normalization,
//! we calculate these for the bases included in the tile (non-blacklisted bases and
//! within the windows / global). Overlapping windows reaching outside the tile are
//! thus clipped to the core (for this ACGT base count only).

use crate::{
    commands::{
        cli_common::*,
        gc_bias::{
            counting::{
                GCCounts, apply_gc_percent_width_correction, build_gc_prefixes,
                count_reference_gc_and_length_by_window, gc_percent_widths,
            },
            interpolation::fill_unsupported_bins_with_polynomial,
            support_masking::{
                build_theoretical_support_mask, create_support_mask_threshold_per_mb,
            },
        },
        ref_gc_bias::config::RefGCBiasConfig,
    },
    shared::{
        bam::Contigs,
        bed::{Windows, load_windows_from_bed},
        blacklist::apply_blacklist_mask_to_seq,
        constants::GC_CORRECTION_SCHEMA_VERSION,
        interval::{IndexedInterval, Interval},
        io::dot_join,
        progress::ProgressFactory,
        reference::{
            ContigFootprintEntry, read_seq_in_range, twobit_contig_footprint, twobit_contig_lengths,
        },
        sampling::{sample_starts_in_core, sampling_density},
        thread_pool::init_global_pool,
        tiled_run::{
            Tile, TileWindowSpan, build_tiles, overlapping_windows_for_tile,
            precompute_tile_window_spans,
        },
        windowing::ensure_plain_bed_windows_not_empty,
    },
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use ndarray::{Array1, Array2};
use ndarray_npy::NpzWriter;
use rand::{Rng, SeedableRng, rngs::StdRng};
use rayon::prelude::*;
use std::{sync::Arc, time::Instant};
use tracing::info;

const COMMAND_TARGET: &str = "ref-gc-bias";

pub fn run(opt: &RefGCBiasConfig) -> Result<()> {
    let start_time = Instant::now();
    opt.fragment_lengths.validate()?;
    ensure!(
        opt.n_positions > 0,
        "--n-positions must be greater than zero"
    );
    let prefix = opt.output_prefix.trim();
    let chromosomes = opt
        .chromosomes
        .resolve_chromosomes(Some(ContigSource::ref_2bit(&opt.ref_genome.ref_2bit)))?;
    let window_opt = opt.windows.resolve_windows();
    opt.check_smoothing_settings()?;

    let min_effective_len = opt
        .fragment_lengths
        .min_fragment_length
        .saturating_sub(2 * u32::from(opt.end_offset));
    ensure!(
        min_effective_len >= 10,
        "Requires at least 10 bases for GC calculation. --min-fragment-length ({}) - 2x --end-offset ({}) is < 10. Please adjust --min-fragment-length.",
        opt.fragment_lengths.min_fragment_length,
        opt.end_offset
    );

    // Create output directory
    ensure_output_dir(&opt.output_dir)?;

    // Load blacklist intervals if provided
    if opt.blacklist.is_some() {
        info!(target: COMMAND_TARGET, "Loading blacklists");
    }
    let blacklist_map = load_blacklist_map(opt.blacklist.as_ref(), 1, 0, &chromosomes)?;

    // Load windows from BED file and merge overlapping/touching intervals (unique positions)
    let windows_map = match &window_opt {
        WindowSpec::Bed(bed) => {
            info!(target: COMMAND_TARGET, "Loading window coordinates");
            let mut wds = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            ensure_plain_bed_windows_not_empty(&wds)?;
            info!(
                target: COMMAND_TARGET,
                "Merging overlapping/touching windows"
            );
            let mut merged: FxHashMap<String, Windows> =
                FxHashMap::with_capacity_and_hasher(wds.len(), Default::default());
            let mut next_idx = 0u64;
            for chr in &chromosomes {
                if let Some(ws) = wds.remove(chr) {
                    // Flatten in-place
                    let (flat, next) = ws.into_flattened_reindexed(next_idx);
                    next_idx = next;
                    merged.insert(chr.clone(), flat);
                }
            }
            Some(merged)
        }
        _ => None,
    };

    // Build chromosome lengths and contigs for tiling without opening BAMs
    let chrom_lengths = twobit_contig_lengths(opt.ref_genome.ref_2bit.clone(), &chromosomes)?;
    let contigs = {
        let mut map: FxHashMap<String, (i32, u32)> =
            FxHashMap::with_capacity_and_hasher(chromosomes.len(), Default::default());
        // Synthetic `tid` index to satisfy `Contigs` even though we do not read BAM headers here
        for (idx, chr) in chromosomes.iter().enumerate() {
            let len = *chrom_lengths
                .get(chr)
                .ok_or_else(|| anyhow::anyhow!("missing chromosome length for {}", chr))?;
            map.insert(chr.clone(), (idx as i32, len as u32));
        }
        Contigs { contigs: map }
    };

    // Precompute GC% bin widths (gc_count -> percent) per fragment length
    let gc_percent_widths = Arc::new(gc_percent_widths(
        opt.fragment_lengths.min_fragment_length as usize,
        opt.fragment_lengths.max_fragment_length as usize,
        opt.end_offset as usize,
    ));
    let total_valid_start_positions: u64 = chrom_lengths
        .values()
        .map(|&chrom_len| {
            if chrom_len as u64 >= opt.fragment_lengths.max_fragment_length as u64 {
                chrom_len as u64 - opt.fragment_lengths.max_fragment_length as u64 + 1
            } else {
                0
            }
        })
        .sum();
    ensure!(
        total_valid_start_positions > 0,
        "No selected chromosome is long enough for --max-fragment-length ({})",
        opt.fragment_lengths.max_fragment_length
    );
    let start_position_sampling_density = sampling_density(
        &chrom_lengths,
        opt.fragment_lengths.max_fragment_length as u64,
        opt.n_positions,
    );
    ensure!(
        start_position_sampling_density <= 1.0,
        "Sampling density {:.4} exceeds 1.0. Reduce --n-positions or increase reference span.",
        start_position_sampling_density
    );

    let mut seed_rng = if let Some(seed) = opt.seed {
        StdRng::seed_from_u64(seed)
    } else {
        let mut thread_rng = rand::rng();
        StdRng::from_rng(&mut thread_rng)
    };

    // Build tiles (core plus padding = max fragment length) to bound memory per worker
    // NOTE: We technically only need padding to the right, but the shared tile-builder
    // applies it symmetrically
    let halo_bp: u32 = opt.fragment_lengths.max_fragment_length;
    let (tiles, _) = build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, None)?;
    // Derive per-tile seeds to keep sampling deterministic without storing all start positions
    let tile_seeds: Vec<u64> = (0..tiles.len()).map(|_| seed_rng.random()).collect();
    let progress = ProgressFactory::new();
    let pb = Arc::new(progress.default_bar(tiles.len() as u64));

    let windows_lookup = windows_map.as_ref();
    let tile_window_spans = Arc::new(precompute_tile_window_spans(
        &tiles,
        |chr| {
            windows_lookup
                .and_then(|m| m.get(chr).map(|w| w.as_slice()))
                .unwrap_or(&[])
        },
        0,
        0, // Windows can extend outside tiles but must overlap
    ));

    // Configure global thread‐pool size
    init_global_pool(opt.n_threads)?;

    let tile_window_spans_for_threads = tile_window_spans.clone();

    info!(target: COMMAND_TARGET, "Counting in tiles");

    // Identity accumulator used by the reducer when no tiles have produced output yet
    let zero_counts = GCCounts::new(
        opt.fragment_lengths.min_fragment_length as usize,
        opt.fragment_lengths.max_fragment_length as usize,
        opt.end_offset as usize,
        (0, 0),
    )?;

    pb.set_position(0);

    let (total_counts, total_covered_acgt_positions) = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| {
            let chr = tile.chr.as_str();
            let tile_span = tile_window_spans_for_threads[tile_idx];
            let windows_chr: Option<&[IndexedInterval<u64>]> = windows_map
                .as_ref()
                .and_then(|m| m.get(&tile.chr).map(|v| v.as_slice()));
            let blacklist_chr: &[Interval<u64>] = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let chr_len = *chrom_lengths
                .get(chr)
                .ok_or_else(|| anyhow::anyhow!("missing chromosome length for {}", chr))?;
            let mut tile_rng = StdRng::seed_from_u64(tile_seeds[tile_idx]);
            let starts = sample_starts_in_core(
                &mut tile_rng,
                tile.core_start() as u64,
                tile.core_end() as u64,
                chr_len as u64,
                opt.fragment_lengths.max_fragment_length as u64,
                start_position_sampling_density,
            );
            // Count tile-local windows using the shared counter logic on window slices
            let res = process_tile(
                tile,
                tile_span.as_ref(),
                chr_len as u64,
                windows_chr,
                starts.as_slice(),
                blacklist_chr,
                opt,
            );
            pb.inc(1);
            res
        })
        .try_reduce(
            || (zero_counts.clone(), 0u64),
            |(mut acc_counts, acc_acgt), (tile_counts, tile_acgt)| {
                acc_counts.merge_from(&tile_counts)?;
                Ok((acc_counts, acc_acgt + tile_acgt))
            },
        )?;

    pb.finish_with_message("| Finished counting");

    // Release tile-level inputs before global aggregation
    drop(tile_window_spans_for_threads);
    drop(tile_window_spans);
    drop(tile_seeds);
    drop(tiles);
    drop(windows_map);
    drop(blacklist_map);

    info!(target: COMMAND_TARGET, "Processing counts");

    let used_start_positions =
        total_counts.sum_for_length(opt.fragment_lengths.min_fragment_length as usize)?;
    ensure!(
        used_start_positions > 0.0,
        "No sampled start positions contributed usable reference GC counts"
    );
    ensure!(
        total_covered_acgt_positions > 0,
        "Selected reference windows contained no usable ACGT bases"
    );

    // Convert counts to Array2 and interpolate zero-counts (single global grid)
    let mut global_counts = total_counts;
    if !opt.skip_smoothing {
        global_counts.smooth_length_rows_in_place(opt.smoothing_sigma, opt.smoothing_radius);
    }
    let mut global_grid = global_counts.to_gc_percent_grid(0, 100)?;
    apply_gc_percent_width_correction(&mut global_grid, &gc_percent_widths)?;

    // Create mask of supported count bins BEFORE interpolation
    // Elements seen less than N times per 1Mb are considered unsupported.
    // These include the theoretically unobservable combinations of fragment lengths and GC percentage bins.
    let threshold_per_mb = 1 + opt.n_positions / 100000000;
    let mut outlier_support_mask = create_support_mask_threshold_per_mb(
        std::slice::from_ref(&global_grid),
        total_covered_acgt_positions,
        threshold_per_mb as f64,
    )
    .expect("support mask should be created");

    if !opt.skip_interpolation {
        info!(target: COMMAND_TARGET, "Interpolating missing counts");
        debug_assert_eq!(
            global_grid.dim(),
            outlier_support_mask.dim(),
            "Support mask and histograms must match shape"
        );
        for (row_idx, mut length_row) in global_grid.outer_iter_mut().enumerate() {
            let row_slice = length_row
                .as_slice_mut()
                .expect("GC histogram rows should be contiguous");
            let mut mask_row = outlier_support_mask.row_mut(row_idx);
            let mask_slice = mask_row
                .as_slice_mut()
                .expect("Support mask rows should be contiguous");
            fill_unsupported_bins_with_polynomial(row_slice, mask_slice, 2, 3, 3, false)?;
        }
    }

    let unobservable_support_mask = build_theoretical_support_mask(
        opt.fragment_lengths.min_fragment_length as usize,
        opt.fragment_lengths.max_fragment_length as usize,
        0,
        global_grid.dim().1 - 1,
        opt.end_offset as usize,
    );

    debug_assert_eq!(
        outlier_support_mask.dim(),
        unobservable_support_mask.dim(),
        "Outlier support mask shape {:?} must match unobservable support mask shape {:?}",
        outlier_support_mask.dim(),
        unobservable_support_mask.dim()
    );

    debug_assert_eq!(
        unobservable_support_mask.dim(),
        global_grid.dim(),
        "Support mask shape {:?} must match histogram shape {:?}",
        unobservable_support_mask.dim(),
        global_grid.dim()
    );

    write_reference_gc_package(
        &opt.output_dir
            .join(dot_join(&[prefix, "ref_gc_package.npz"])),
        &global_grid,
        &unobservable_support_mask,
        &outlier_support_mask,
        &gc_percent_widths,
        opt.fragment_lengths.min_fragment_length as usize,
        opt.fragment_lengths.max_fragment_length as usize,
        opt.end_offset,
        opt.skip_interpolation,
        opt.smoothing_radius,
        opt.smoothing_sigma,
        opt.skip_smoothing,
        &chromosomes,
        &twobit_contig_footprint(&opt.ref_genome.ref_2bit)?,
    )
    .context("Writing reference GC package failed")?;

    let elapsed = start_time.elapsed();
    info!(
        target: COMMAND_TARGET,
        "Windows covered {} total ACGT bases",
        total_covered_acgt_positions
    );
    info!(
        target: COMMAND_TARGET,
        "Used {:.0} start positions at length {}",
        used_start_positions, opt.fragment_lengths.min_fragment_length
    );
    info!(target: COMMAND_TARGET, "Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn write_reference_gc_package(
    path: &std::path::Path,
    counts: &Array2<f64>,
    support_unobservables: &Array2<bool>,
    support_outliers: &Array2<bool>,
    gc_percent_widths: &Array2<u16>,
    length_min: usize,
    length_max: usize,
    end_offset: u8,
    skip_interpolation: bool,
    smoothing_radius: u8,
    smoothing_sigma: f64,
    skip_smoothing: bool,
    chromosomes: &[String],
    reference_contig_footprint: &[ContigFootprintEntry],
) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut npz = NpzWriter::new(file);
    npz.add_array("counts", counts)?;
    npz.add_array("support_mask_unobservables", support_unobservables)?;
    npz.add_array("support_mask_outliers", support_outliers)?;
    npz.add_array("gc_percent_widths", gc_percent_widths)?;
    npz.add_array("version", &Array1::from(vec![GC_CORRECTION_SCHEMA_VERSION]))?;
    npz.add_array(
        "length_range",
        &Array1::from(vec![length_min as u32, length_max as u32]),
    )?;
    npz.add_array("end_offset", &Array1::from(vec![end_offset as u32]))?;
    npz.add_array(
        "skip_interpolation",
        &Array1::from(vec![skip_interpolation]),
    )?;
    npz.add_array(
        "smoothing_radius",
        &Array1::from(vec![smoothing_radius as u32]),
    )?;
    npz.add_array("smoothing_sigma", &Array1::from(vec![smoothing_sigma]))?;
    npz.add_array("skip_smoothing", &Array1::from(vec![skip_smoothing]))?;
    let chromosomes_json = serde_json::to_vec(chromosomes)?;
    npz.add_array("chromosomes_json", &Array1::from(chromosomes_json))?;
    npz.add_array(
        "reference_contig_footprint_json",
        &Array1::from(serde_json::to_vec(reference_contig_footprint)?),
    )?;
    npz.finish()?;
    Ok(())
}

fn process_tile(
    tile: &Tile,
    tile_window_span: Option<&TileWindowSpan>,
    chrom_len: u64,
    windows: Option<&[IndexedInterval<u64>]>,
    start_positions: &[usize],
    blacklist_intervals: &[Interval<u64>],
    opt: &RefGCBiasConfig,
) -> Result<(GCCounts, u64)> {
    let core_start = tile.core_start() as u64;
    let core_end = tile.core_end() as u64;
    if core_start >= core_end || core_start >= chrom_len {
        let empty = GCCounts::new(
            opt.fragment_lengths.min_fragment_length as usize,
            opt.fragment_lengths.max_fragment_length as usize,
            opt.end_offset as usize,
            (0, 0),
        )?;
        return Ok((empty, 0));
    }

    // Absolute genomic sequence loaded for this tile via 2bit
    let seq_start = tile.fetch_start().min(tile.core_start()) as u64;
    let seq_end = tile.fetch_end().min(chrom_len as u32) as u64;

    // Load only the tile span (core plus padding) so starts in the core have full context
    let mut seq_bytes = read_seq_in_range(
        &opt.ref_genome.ref_2bit,
        &tile.chr,
        seq_start as usize..seq_end as usize,
    )?;
    apply_blacklist_mask_to_seq(&mut seq_bytes, blacklist_intervals, seq_start);
    let gc_prefixes = build_gc_prefixes(&seq_bytes);

    // Delete seq_bytes from memory
    drop(seq_bytes);

    let core_start_usize = core_start as usize;
    let core_end_usize = core_end as usize;

    // Build/redefine windows to have tile-local coordinates.
    // These windows are clipped to start in the core but may extend into the right halo.
    // We keep starts inside the core (so starts are unique per tile) while letting fragment ends
    // reach into the fetched halo, which carries the needed sequence context
    // Absolute -> tile-local coordinates: We offset coordinates by sequence start
    // so they are relative to the loaded reference sequence
    let mut tile_windows: Vec<IndexedInterval<u64>> = Vec::new();
    if let Some(win_chr) = windows {
        let iter = overlapping_windows_for_tile(win_chr, tile, tile_window_span);
        for window in iter {
            let start_abs = window.start().max(core_start).max(seq_start);
            let end_abs = window.end().min(seq_end);
            if end_abs <= start_abs {
                continue;
            }
            tile_windows.push(IndexedInterval::new(
                start_abs - seq_start,
                end_abs - seq_start,
                // Preserve the original window index so downstream counts map back
                // to the same BED window identity
                window.idx(),
            )?);
        }
    } else {
        // Global mode: one tile-local window starting from the tile core start
        // and ending at the end of the loaded reference sequence
        tile_windows.push(IndexedInterval::new(
            // Starts at the core start in the loaded sequence (offset by the left halo)
            core_start.saturating_sub(seq_start),
            seq_end.saturating_sub(seq_start),
            0,
        )?);
    }
    if tile_windows.is_empty() {
        let empty = GCCounts::new(
            opt.fragment_lengths.min_fragment_length as usize,
            opt.fragment_lengths.max_fragment_length as usize,
            opt.end_offset as usize,
            (0, 0),
        )?;
        return Ok((empty, 0));
    }

    // Filter starts to the core and shift to tile-local coordinates
    let core_start_idx = start_positions.partition_point(|&s| s < core_start_usize);
    let core_end_idx = start_positions.partition_point(|&s| s < core_end_usize);
    let tile_starts: Vec<usize> = start_positions[core_start_idx..core_end_idx]
        .iter()
        .map(|s| s.saturating_sub(seq_start as usize))
        .collect();

    // Allocate per-tile window accumulators (global sum comes from merging them)
    let mut counts_by_bin = vec![
        GCCounts::new(
            opt.fragment_lengths.min_fragment_length as usize,
            opt.fragment_lengths.max_fragment_length as usize,
            opt.end_offset as usize,
            (0, 0),
        )?;
        tile_windows.len()
    ];

    // Count the reference GC by length and windows
    // The tile_windows, tile_starts, and passed chrom_len (seq_end - seq_start)
    // all use tile-local coordinates / size
    count_reference_gc_and_length_by_window(
        &mut counts_by_bin,
        &gc_prefixes,
        (
            opt.fragment_lengths.min_fragment_length as u64,
            opt.fragment_lengths.max_fragment_length as u64 + 1,
        ),
        tile_windows.as_slice(),
        tile_starts.as_slice(),
        seq_end - seq_start,
        1.0,
        1u32,
        opt.end_offset as usize,
    )?;

    // Release per-tile start positions after counting
    drop(tile_starts);

    // Compute ACGT coverage only within the core so bases are not double-counted across tiles
    let mut total_acgt_in_core = 0u64;
    let core_start_local = core_start - seq_start;
    let core_end_local = core_end - seq_start;
    for window in &tile_windows {
        let clipped_start = window.start().max(core_start_local);
        let clipped_end = window.end().min(core_end_local);
        if clipped_end <= clipped_start {
            continue;
        }
        let acgt =
            gc_prefixes.acgt[clipped_end as usize] - gc_prefixes.acgt[clipped_start as usize];
        total_acgt_in_core += acgt as u64;
    }

    // Release per-tile buffers before merging counts
    drop(gc_prefixes);
    drop(tile_windows);

    // Sum tile-local windows into a single accumulator
    let merged = counts_by_bin.into_iter().try_fold(
        GCCounts::new(
            opt.fragment_lengths.min_fragment_length as usize,
            opt.fragment_lengths.max_fragment_length as usize,
            opt.end_offset as usize,
            (0, 0),
        )?,
        |mut acc, c| -> Result<GCCounts> {
            acc.merge_from(&c)?;
            Ok(acc)
        },
    )?;

    Ok((merged, total_acgt_in_core))
}

#[cfg(test)]
mod tests {
    include!("ref_gc_bias_tests.rs");
}
