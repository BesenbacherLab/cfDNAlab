//! Runner for calculating reference k-mer backgrounds.
//!
//! Counts reference k-mers in a tiled pass, with per-window scaling factors
//! enabling downstream count reconstruction.
//!
//! K-mers are only counted in the left->right direction.
//! For reverse-orientation k-mers, you can reverse-complement
//! the k-mers downstream.

// TODO: Perhaps make it an option to add the revcomp counts to the output?
// It's literally just a loop with:
//   counts[revcomp(motif)] += counts[motif]  # (though with a tmp dict in-between)

use crate::{
    bam::Contigs,
    command_run::{CommandRunResult, RunOptions, status_info},
    commands::{
        cli_common::*,
        ref_kmers::{
            aggregation::{CollectedRefKmerFrequencies, collect_ref_kmer_frequencies},
            config::RefKmersConfig,
            counting::{
                Enc, KmerCountsByWindow, SelectedKmerCountsByWindow, count_kmers_by_window,
                prepare_ref_kmer_window_source,
            },
            tiling::{
                TileResult, build_selected_tile_count_records, build_tile_count_records,
                serialize_selected_tile_counts, serialize_tile_counts,
            },
            zarr::{
                RefKmerRowMetadata, RefKmerWindowRowMode, RefKmerZarrPackage,
                ensure_all_ref_kmers_output_size, ensure_dense_ref_kmer_output_size,
                grouped_ref_kmer_row_metadata, ref_kmer_axis_is_complete, ref_kmer_zarr_path,
                write_ref_kmer_zarr,
            },
        },
    },
    shared::{
        bed::{load_grouped_windows_from_bed, load_windows_from_bed},
        blacklist::apply_blacklist_mask_to_seq,
        interval::{IndexedInterval, Interval},
        io::{FinalOutputFiles, dot_join},
        kmers::{
            kmer_codec::build_optional_kmer_spec,
            motifs_file::{
                SelectedMotifHalfSpec, SelectedMotifLookup, parse_selected_ref_kmers_file,
            },
        },
        overlaps::{
            DEFAULT_BROAD_WINDOW_MIN_BP, TileBedWindowView, build_bed_windows_by_chr,
            precompute_tile_bed_window_spans,
        },
        progress::ProgressFactory,
        reference::{read_seq_in_range, twobit_contig_footprint, twobit_contig_lengths},
        temp_chrom_names::TempChromNameMap,
        thread_pool::init_global_pool,
        tiled_run::{TempDirGuard, Tile, build_tiles},
        windowing::{
            DistributionWindowContext, build_bin_info, compute_window_offsets,
            ensure_plain_bed_windows_not_empty,
        },
    },
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use rayon::prelude::*;
use std::{
    path::{Path, PathBuf},
    sync::Arc,
    time::Instant,
};
use tracing::info;

const COMMAND_TARGET: &str = "ref-kmers";

/// Result from `ref-kmers`.
///
/// The command writes a Zarr package with frequencies, row scaling factors, and metadata needed
/// for downstream reference-bias correction. The result records the package path and final output
/// file list.
#[derive(Debug)]
pub struct RefKmersRunResult {
    /// Empty counter placeholder for the shared command result interface.
    pub counters: (),
    /// Final reference k-mer frequency package path.
    pub ref_kmer_counts_path: PathBuf,
    /// Final output files produced by the command.
    pub output_files: Vec<PathBuf>,
}

impl CommandRunResult for RefKmersRunResult {
    type Counters = ();

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&Path> {
        Some(self.ref_kmer_counts_path.as_path())
    }
}

/// Run the `ref-kmers` command.
///
/// This command counts reference k-mers from eligible positions of the reference genome,
/// optionally restricted to a motifs file. It writes row-wise frequencies and scaling factors so
/// downstream code can reconstruct counts when needed.
///
/// Reporting is controlled by `options`. `show_progress` controls progress bars and
/// `log_statuses` controls status messages. This command does not print a statistics summary.
///
/// Parameters
/// ----------
/// - `opt`:
///   Fully resolved configuration for the `ref-kmers` command.
/// - `options`:
///   Reporting controls for progress bars and status logs.
///
/// Returns
/// -------
/// - `Ok(RefKmersRunResult)`:
///   Output path information for the completed run.
///
/// Errors
/// ------
/// Returns an error when the configuration is invalid, the reference or BED input cannot be read,
/// or the output cannot be written.
pub fn run_ref_kmers(opt: &RefKmersConfig, options: RunOptions) -> Result<RefKmersRunResult> {
    let start_time = Instant::now();
    opt.validate()?;
    let prefix = opt.output_prefix.trim();
    let window_opt = opt.windows.resolve_windows();
    let fetch_window_opt = window_opt.as_fetch_window_spec();

    let selected_motifs = match opt.motifs_file.as_deref() {
        None => None,
        Some(motifs_file) => {
            ensure!(
                !opt.canonical,
                "`--motifs-file` cannot be combined with `--canonical`. Use the motifs file group column to define collapsed targets."
            );
            Some(parse_selected_ref_kmers_file(
                motifs_file,
                opt.kmer_size as usize,
            )?)
        }
    };
    if let Some(selected_motifs) = selected_motifs.as_ref() {
        status_info!(
            options,
            target: COMMAND_TARGET,
            "Loaded {} selected reference k-mer target(s)",
            selected_motifs.labels.len()
        );
    }

    if options.log_equivalent_cli {
        let command = crate::ToCliCommand::to_cli_string(opt)?;
        let message = crate::command_run::equivalent_cli_log_message(&command);
        info!(target: COMMAND_TARGET, "{message}");
    }

    let chromosomes = opt
        .chromosomes
        .resolve_chromosomes(Some(ContigSource::ref_2bit(&opt.ref_genome.ref_2bit)))?;

    // Create output directory
    ensure_output_dir(&opt.output_dir)?;

    // Load blacklist intervals if provided
    if opt.blacklist.is_some() {
        status_info!(options, target: COMMAND_TARGET, "Loading blacklists");
    }
    let blacklist_map = load_blacklist_map(opt.blacklist.as_ref(), 1, 0, &chromosomes)?;

    // Load windows from BED file
    let windows_map = match &window_opt {
        DistributionWindowSpec::Bed(bed) => {
            status_info!(options, target: COMMAND_TARGET, "Loading window coordinates");
            let windows = load_windows_from_bed(bed, Some(chromosomes.as_slice()), None, None)?;
            ensure_plain_bed_windows_not_empty(&windows)?;
            Some(windows)
        }
        _ => None,
    };
    let (grouped_windows_map, group_idx_to_name) = match &window_opt {
        DistributionWindowSpec::GroupedBed(bed) => {
            status_info!(options, target: COMMAND_TARGET, "Loading grouped window coordinates");
            let (windows_map, group_idx_to_name, _strand_detection) =
                load_grouped_windows_from_bed(
                    bed,
                    Some(chromosomes.as_slice()),
                    false,
                    None,
                    None,
                )?;
            ensure!(
                !group_idx_to_name.is_empty(),
                "grouped BED file did not contain any valid windows on the selected chromosomes"
            );
            (Some(windows_map), Some(group_idx_to_name))
        }
        _ => (None, None),
    };
    let indexed_windows_map: Option<FxHashMap<String, Vec<IndexedInterval<u64>>>> =
        if let Some(windows_map) = windows_map.as_ref() {
            Some(
                windows_map
                    .iter()
                    .map(|(chr, windows)| (chr.clone(), windows.as_slice().to_vec()))
                    .collect(),
            )
        } else {
            grouped_windows_map.as_ref().map(|windows_map| {
                windows_map
                    .iter()
                    .map(|(chr, windows)| (chr.clone(), windows.windows_as_slice().to_vec()))
                    .collect()
            })
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

    let halo_bp = opt.kmer_size as u32;
    let (tiles, _) = build_tiles(&chromosomes, &contigs, opt.tile_size, halo_bp, None)?;
    let progress = ProgressFactory::with_enabled(options.show_progress);
    let pb = Arc::new(progress.default_bar(tiles.len() as u64));

    let bed_windows_by_chr = indexed_windows_map
        .as_ref()
        .map(|windows| build_bed_windows_by_chr(windows, DEFAULT_BROAD_WINDOW_MIN_BP));
    let tile_bed_window_spans = Arc::new(bed_windows_by_chr.as_ref().map(|bed_windows| {
        precompute_tile_bed_window_spans(&tiles, bed_windows, 0, halo_bp as u64)
    }));

    // Configure global thread‐pool size
    init_global_pool(opt.n_threads)?;

    let tile_bed_window_spans_for_threads = tile_bed_window_spans.clone();

    // DistributionWindowSpec preserves grouped BED row identity for output, while WindowSpec is the
    // plain coordinate view used for fetch narrowing and fixed-size offset calculation.
    // Fixed-size windows need a per-chromosome row offset because overlap lookup returns
    // chromosome-local bin indices. BED-like windows already carry their final output row ids.
    let (total_windows, chr_offsets_map): (u64, FxHashMap<String, u64>) = match &window_opt {
        DistributionWindowSpec::GroupedBed(_) => (
            group_idx_to_name
                .as_ref()
                .context("group_idx_to_name missing for grouped BED mode")?
                .len() as u64,
            chromosomes.iter().map(|chr| (chr.clone(), 0_u64)).collect(),
        ),
        _ => compute_window_offsets(
            &fetch_window_opt,
            &chromosomes,
            &contigs,
            windows_map.as_ref(),
        )?,
    };
    let total_windows_usize: usize = total_windows
        .try_into()
        .context("reference k-mer output row count does not fit in usize")?;
    let chr_offsets = Arc::new(chr_offsets_map);
    let chr_offsets_for_threads = chr_offsets.clone();

    // Final sparse output decodes observed radix-5 motif keys into strings. Motifs-file output
    // uses public labels from the file, so these are intentionally `None` when `--motifs-file` is
    // used.
    let kmer_decode_spec = if selected_motifs.is_none() {
        build_optional_kmer_spec(opt.kmer_size as usize, "kmer")?
    } else {
        None
    };

    // Tile counting needs an encoder for the enabled k-mer space. Without a motifs file this is
    // the full radix-5 codec. With a motifs file, k-mers up to the radix-5 limit still use full
    // radix-5 codes, while larger k-mers use byte-backed selected subspaces.
    let kmer_counting_spec = match selected_motifs.as_ref() {
        Some(lookup) => lookup.inside_spec.clone(),
        None => kmer_decode_spec
            .as_ref()
            .cloned()
            .map(SelectedMotifHalfSpec::from_radix5),
    };

    if opt.all_motifs {
        match selected_motifs.as_ref() {
            Some(lookup) => {
                ensure_dense_ref_kmer_output_size(total_windows_usize, lookup.labels.len())?
            }
            None => {
                let kmer_decode_spec = kmer_decode_spec
                    .as_ref()
                    .context("missing k-mer decode spec for all-motifs reference k-mer output")?;
                ensure_all_ref_kmers_output_size(
                    kmer_decode_spec.k,
                    opt.canonical,
                    total_windows_usize,
                )?;
            }
        }
    }

    // Create droppable temporary directory
    let temp_dir_guard =
        TempDirGuard::new(&opt.output_dir, prefix).context("create per-run temp dir")?;
    let temp_dir = temp_dir_guard.path();
    let mut final_outputs = FinalOutputFiles::new(temp_dir)?;
    let temp_chrom_name_map = TempChromNameMap::from_contigs(&chromosomes)?;

    let counts_prefix = &dot_join(&[prefix, "counts"]);

    status_info!(options, target: COMMAND_TARGET, "Counting per tile");

    pb.set_position(0);

    let tile_results: Vec<TileResult> = tiles
        .par_iter()
        .enumerate()
        .map(|(tile_idx, tile)| {
            let chr = tile.chr.as_str();
            let blacklist_chr: &[Interval<u64>] = blacklist_map
                .get(&tile.chr)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
            let chr_len = *chrom_lengths
                .get(chr)
                .ok_or_else(|| anyhow::anyhow!("missing chromosome length for {}", chr))?;
            let tile_bed_window_view = match bed_windows_by_chr
                .as_ref()
                .and_then(|windows_by_chromosome| windows_by_chromosome.get(&tile.chr))
            {
                Some(chromosome_windows) => {
                    let spans = tile_bed_window_spans_for_threads
                        .as_ref()
                        .as_ref()
                        .context("BED reference k-mer counting requires tile BED window spans")?;
                    Some(TileBedWindowView {
                        chromosome_windows,
                        spans: &spans[tile_idx],
                    })
                }
                None => None,
            };

            // Count k-mers that start in the tile core. The k-mer span may extend into the halo.
            let tile_result = process_tile(
                opt,
                tile,
                tile_bed_window_view,
                chr_len as u64,
                *chr_offsets_for_threads.get(&tile.chr).unwrap_or(&0),
                &window_opt,
                blacklist_chr,
                temp_dir,
                counts_prefix,
                &temp_chrom_name_map,
                kmer_counting_spec.as_ref(),
                selected_motifs.as_ref(),
            );
            pb.inc(1);
            tile_result
        })
        .collect::<Result<Vec<_>>>()? // Short-circuits on the first Err
        .into_iter()
        .flatten()
        .collect();

    if options.show_progress {
        pb.finish_with_message("| Finished counting");
    } else {
        pb.finish_and_clear();
    }

    // Release tile-level views before global aggregation
    drop(tile_bed_window_spans_for_threads);
    drop(tile_bed_window_spans);
    drop(tiles);
    drop(bed_windows_by_chr);
    drop(indexed_windows_map);

    status_info!(options, target: COMMAND_TARGET, "Processing counts");
    let reference_contig_footprint = twobit_contig_footprint(&opt.ref_genome.ref_2bit)?;

    let bin_info = if matches!(&window_opt, DistributionWindowSpec::GroupedBed(_)) {
        Vec::new()
    } else {
        build_bin_info(
            &fetch_window_opt,
            &chromosomes,
            &contigs,
            windows_map.as_ref(),
            &blacklist_map,
            chr_offsets.as_ref(),
        )?
    };
    let row_metadata = match &window_opt {
        DistributionWindowSpec::Global => RefKmerRowMetadata::Global,
        DistributionWindowSpec::GroupedBed(_) => {
            let group_idx_to_name = group_idx_to_name
                .as_ref()
                .context("group_idx_to_name missing when writing grouped outputs")?;
            let grouped_windows_map = grouped_windows_map
                .as_ref()
                .context("grouped windows missing when writing grouped outputs")?;
            RefKmerRowMetadata::Groups(grouped_ref_kmer_row_metadata(
                group_idx_to_name,
                &chromosomes,
                grouped_windows_map,
                &blacklist_map,
            )?)
        }
        DistributionWindowSpec::Size(_) => RefKmerRowMetadata::Windows {
            bin_info: &bin_info,
            row_mode: RefKmerWindowRowMode::Size,
        },
        DistributionWindowSpec::Bed(_) => RefKmerRowMetadata::Windows {
            bin_info: &bin_info,
            row_mode: RefKmerWindowRowMode::Bed,
        },
    };

    let CollectedRefKmerFrequencies {
        frequency_bins,
        motif_order,
        column_kind,
    } = collect_ref_kmer_frequencies(
        &tile_results,
        selected_motifs.as_ref(),
        kmer_decode_spec.as_ref(),
        total_windows_usize,
        opt.canonical,
        opt.all_motifs,
    )?;
    let non_empty_rows = frequency_bins
        .frequency_bins
        .iter()
        .filter(|row| !row.is_empty())
        .count();
    status_info!(
        options,
        target: COMMAND_TARGET,
        "Prepared reference k-mer frequencies for {} non-empty row(s) across {} output row(s) and {} motif column(s)",
        non_empty_rows,
        total_windows,
        motif_order.len()
    );
    // Dense storage is forced by `--all-motifs`. It is also used when an observed or selected
    // motif axis is already complete, even if the user did not request `--all-motifs`
    let write_dense_output = opt.all_motifs
        || ref_kmer_axis_is_complete(opt.kmer_size, opt.canonical, column_kind, &motif_order);

    // Write every final output to the temp directory before moving any of them into place
    // This keeps failed writes from appearing completed
    let ref_kmer_counts_path = ref_kmer_zarr_path(&opt.output_dir, prefix);
    let temp_ref_kmer_counts_path = final_outputs.temp_path_for(&ref_kmer_counts_path)?;
    write_ref_kmer_zarr(
        &temp_ref_kmer_counts_path,
        RefKmerZarrPackage {
            frequency_bins: &frequency_bins.frequency_bins,
            row_scaling_factors: &frequency_bins.row_scaling_factors,
            motif_labels: &motif_order,
            column_kind,
            row_metadata,
            write_dense_output,
            kmer_size: opt.kmer_size,
            canonical: opt.canonical,
            // This is the user option stored in metadata, not the dense-storage decision above
            all_motifs: opt.all_motifs,
            assign_by: opt.assign_by,
            reference_contig_footprint: &reference_contig_footprint,
        },
    )?;

    final_outputs.record(temp_ref_kmer_counts_path, ref_kmer_counts_path.clone())?;
    final_outputs.move_into_place()?;

    let elapsed = start_time.elapsed();
    status_info!(options, target: COMMAND_TARGET, "Elapsed time: {:.2?}", elapsed);
    Ok(RefKmersRunResult {
        counters: (),
        ref_kmer_counts_path: ref_kmer_counts_path.clone(),
        output_files: vec![ref_kmer_counts_path],
    })
}

fn process_tile(
    opt: &RefKmersConfig,
    tile: &Tile,
    tile_bed_window_view: Option<TileBedWindowView<'_>>,
    chrom_len: u64,
    chr_window_idx_offset: u64,
    window_opt: &DistributionWindowSpec,
    blacklist_intervals: &[Interval<u64>],
    temp_dir: &Path,
    counts_prefix: &str,
    temp_chrom_name_map: &TempChromNameMap,
    kmers_spec: Option<&SelectedMotifHalfSpec>,
    selected_motifs: Option<&SelectedMotifLookup>,
) -> Result<Option<TileResult>> {
    let core_start = tile.core_start() as u64;
    let core_end = tile.core_end() as u64;
    if core_start >= core_end || core_start >= chrom_len {
        return Ok(None);
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

    // Path for this tile's temporary outputs
    let chr_token = temp_chrom_name_map.token_for(tile.chr.as_str())?;
    let counts_path = temp_dir.join(format!(
        "{prefix}.{chr}.{idx}.ref_kmer_counts.bin",
        prefix = counts_prefix,
        chr = chr_token,
        idx = tile.index
    ));

    let k = opt.kmer_size;
    let spec = kmers_spec.context("missing k-mer counting spec")?;
    let positional_codes = spec.build_left_aligned_codes(&seq_bytes);
    let enc = Enc {
        k,
        codes: &positional_codes,
        none: spec.missing_reference_code(),
        n: spec.masked_reference_code(),
    };

    let mut counts_by_window = KmerCountsByWindow::default();
    let mut selected_counts_by_window = SelectedKmerCountsByWindow::default();
    let owned_start = core_start.saturating_sub(seq_start);
    let owned_end = core_end.min(seq_end).saturating_sub(seq_start);
    let window_context = DistributionWindowContext {
        spec: window_opt,
        chr_idx_offset: chr_window_idx_offset,
    };
    let owned_starts = owned_start..owned_end;
    let window_source = prepare_ref_kmer_window_source(
        &window_context,
        tile_bed_window_view,
        owned_starts.clone(),
        seq_start,
        chrom_len,
        k as u64,
        opt.assign_by,
    )?;
    if let Some(window_source) = window_source {
        count_kmers_by_window(
            &mut counts_by_window,
            &mut selected_counts_by_window,
            &enc,
            &window_context,
            window_source,
            owned_starts,
            seq_start,
            chrom_len,
            opt.assign_by,
            selected_motifs,
        )?;
    }

    if selected_motifs.is_some() {
        let count_records = build_selected_tile_count_records(selected_counts_by_window);
        serialize_selected_tile_counts(&counts_path, &count_records)?;
    } else {
        let count_records = build_tile_count_records(counts_by_window);
        serialize_tile_counts(&counts_path, &count_records)?;
    }

    Ok(Some(TileResult { counts_path }))
}

#[cfg(test)]
mod tests {
    include!("ref_kmers_tests.rs");
}
