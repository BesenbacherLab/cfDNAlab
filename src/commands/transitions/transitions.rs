use crate::{
    command_run::{CommandRunResult, RunOptions},
    commands::{
        cli_common::{ensure_output_dir, validate_output_prefix},
        counters::FragmentKmersCounters,
        fragment_kmers::{config::*, fragment_kmers},
        run_statistics::{
            DEFAULT_FRAGMENT_STATISTICS_LABELS, FragmentRunStatisticsOptions, GCStatisticsSummary,
            TILE_DOUBLE_COUNT_NOTE, print_fragment_run_statistics,
        },
        transitions::config::TransitionsConfig,
    },
    shared::{
        io::{FinalOutputFiles, dot_join},
        tiled_run::TempDirGuard,
    },
};
use anyhow::{Context, Result, ensure};
use fxhash::FxHashMap;
use ndarray::{Array3, Axis, s};
use ndarray_npy::{read_npy, write_npy};
use std::collections::hash_map::Entry;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Normalise positional k-mer counts into transition probabilities.
///
/// Converts the (possibly scaled and/or corrected) counts produced by `fragment-kmers` into
/// per-prefix transition frequencies so every slice along the motif axis sums to one for a
/// fixed `(window, position, prefix)` triple.
///
/// Parameters
/// ----------
/// - `counts`:
///     Dense `(windows, positions, motifs)` cube of absolute counts.
///
/// - `order`:
///     Transition depth, equal to `k - 1` when the motifs have length `k`.
///
/// - `motifs`:
///     Motif list that matches the fastest-changing axis of `counts`.
///
/// Returns
/// -------
/// - `freqs`:
///     Dense cube with identical shape where each motif entry stores its conditional probability.
pub(crate) fn compute_transition_frequencies(
    counts: &Array3<f64>,
    order: u8,
    motifs: &[String],
) -> Result<Array3<f64>> {
    ensure!(
        !motifs.is_empty(),
        "transition frequency calculation requires at least one motif"
    );
    ensure!(
        counts.len_of(Axis(2)) == motifs.len(),
        "motifs axis must match counts axis ({} != {})",
        motifs.len(),
        counts.len_of(Axis(2))
    );

    let prefix_len = order as usize;
    ensure!(
        motifs.iter().all(|motif| motif.len() >= prefix_len),
        "motif shorter than requested order detected for order {}",
        order
    );

    // Build prefix lookup for fast bucket accumulation
    let mut prefix_map: FxHashMap<Vec<u8>, usize> =
        FxHashMap::with_capacity_and_hasher(motifs.len(), Default::default());
    let mut motif_prefix_ids: Vec<usize> = Vec::with_capacity(motifs.len());
    let mut next_prefix_id = 0usize;

    for motif in motifs {
        let key = motif.as_bytes()[..prefix_len].to_vec();
        let prefix_id = match prefix_map.entry(key) {
            Entry::Occupied(existing) => *existing.get(),
            Entry::Vacant(vacant) => {
                let id = next_prefix_id;
                next_prefix_id += 1;
                vacant.insert(id);
                id
            }
        };
        motif_prefix_ids.push(prefix_id);
    }

    let mut freqs = Array3::<f64>::zeros(counts.raw_dim());
    let prefix_slots = prefix_map.len().max(1);
    let mut totals = vec![0.0f64; prefix_slots];

    let windows = counts.len_of(Axis(0));
    let positions = counts.len_of(Axis(1));

    for w_idx in 0..windows {
        for p_idx in 0..positions {
            // Reset per-prefix totals before scanning this position
            totals.iter_mut().for_each(|total| *total = 0.0);
            let counts_slice = counts.slice(s![w_idx, p_idx, ..]);
            let mut freq_slice = freqs.slice_mut(s![w_idx, p_idx, ..]);

            // Accumulate absolute counts per prefix
            for (motif_idx, &value) in counts_slice.iter().enumerate() {
                totals[motif_prefix_ids[motif_idx]] += value;
            }

            // Convert counts into conditional probabilities
            for (motif_idx, slot) in freq_slice.iter_mut().enumerate() {
                let total = totals[motif_prefix_ids[motif_idx]];
                if total > 0.0 {
                    *slot = counts_slice[motif_idx] / total;
                } else {
                    *slot = 0.0;
                }
            }
        }
    }

    Ok(freqs)
}

/// Result from `transitions`.
///
/// The command writes transition-frequency arrays derived from fragment k-mer counts. The result
/// records all final output files and the counters from the underlying k-mer producer.
#[derive(Debug)]
pub struct TransitionsRunResult {
    /// Fragment counters from the internal `fragment-kmers` run.
    pub counters: FragmentKmersCounters,
    /// Final output files produced by the command.
    pub output_files: Vec<PathBuf>,
}

impl CommandRunResult for TransitionsRunResult {
    type Counters = FragmentKmersCounters;

    fn counters(&self) -> &Self::Counters {
        &self.counters
    }

    fn output_files(&self) -> &[PathBuf] {
        &self.output_files
    }

    fn primary_output(&self) -> Option<&Path> {
        self.output_files.first().map(PathBuf::as_path)
    }
}

/// Run the `transitions` command.
///
/// This command runs `fragment-kmers` as an internal producer, then normalizes positional k-mer
/// counts into transition probabilities. The output arrays keep the same window and position axes
/// as the source k-mer counts.
///
/// Reporting is controlled by `options.report_statistics`, which prints the final summary. The
/// internal `fragment-kmers` producer is run without its own reporting side effects.
///
/// Parameters
/// ----------
/// - `opt`:
///     Fully resolved configuration for the `transitions` command.
/// - `options`:
///     Reporting control for the transition summary.
///
/// Returns
/// -------
/// - `Ok(TransitionsRunResult)`:
///     Counters and output paths for the completed run.
///
/// Errors
/// ------
/// Returns an error when the internal k-mer run fails, a count array cannot be read, or transition
/// outputs cannot be written.
pub fn run_transitions(
    opt: &TransitionsConfig,
    options: RunOptions,
) -> Result<TransitionsRunResult> {
    let start_time = Instant::now();
    opt.shared_args.fragment_lengths.validate()?;

    let prefix = opt.shared_args.output_prefix.trim();
    validate_output_prefix(prefix)?;

    if options.log_equivalent_cli {
        let command = crate::ToCliCommand::to_cli_string(opt)?;
        let message = crate::command_run::equivalent_cli_log_message(&command);
        tracing::info!(target: "transitions", "{message}");
    }

    // Create output directory
    ensure_output_dir(&opt.shared_args.ioc.output_dir)?;

    // Build temporary directory
    let mut temp_root_guard = TempDirGuard::new(&opt.shared_args.ioc.output_dir, prefix)
        .context("create per-run temp dir")?;
    let temp_root = temp_root_guard.path().to_path_buf();
    let mut final_outputs = FinalOutputFiles::new(&temp_root)?;

    let mut tmp_ioc = opt.shared_args.ioc.clone();
    tmp_ioc.output_dir = temp_root.join("kmer_counts");

    let kmer_sizes = opt.orders.iter().map(|o| o + 1).collect();

    let mut fk_shared_args = opt.shared_args.clone();
    fk_shared_args.ioc = tmp_ioc;

    let mut fk_cfg = FragmentKmersConfig {
        shared_args: fk_shared_args,
        kmer_sizes,
        canonical: false,
        positional_counts: true,
        save_sparse: false, // TODO: Might be necessary later?
    };
    fk_cfg.set_output_prefix(prefix.to_string());

    let fragment_kmers_result = fragment_kmers::run_fragment_kmers(
        &fk_cfg,
        RunOptions {
            log_equivalent_cli: options.log_equivalent_cli,
            ..RunOptions::new_quiet()
        },
    )?;
    let global_counter = fragment_kmers_result.counters;

    let counts_dir = fk_cfg.shared_args.ioc.output_dir.as_path();
    let final_dir = &opt.shared_args.ioc.output_dir;
    let prefix = opt.shared_args.output_prefix.trim();
    let groups = ["left", "right", "mid"];
    let mut output_files = Vec::new();

    // Process each requested transition order
    for &order in &opt.orders {
        let k = order + 1;
        let mut wrote_group = false;
        for group in groups {
            let counts_path =
                counts_dir.join(dot_join(&[prefix, &format!("k{k}_{group}_counts.npy")]));
            if !counts_path.exists() {
                continue;
            }

            // Load positional cube for this orientation
            let counts: Array3<f64> = read_npy(&counts_path)
                .with_context(|| format!("loading {}", counts_path.display()))?;

            let motifs_path =
                counts_dir.join(dot_join(&[prefix, &format!("k{k}_{group}_motifs.txt")]));
            let motifs_raw = fs::read_to_string(&motifs_path)
                .with_context(|| format!("reading {}", motifs_path.display()))?;
            let motifs: Vec<String> = motifs_raw.lines().map(|line| line.to_owned()).collect();

            // Normalise counts into conditional frequencies
            let freqs = compute_transition_frequencies(&counts, order, &motifs)?;

            // Write final outputs to the temp folder first
            // They move into output_dir after all requested transition files have been written
            let output_path =
                final_dir.join(dot_join(&[prefix, &format!("k{k}_{group}_freqs.npy")]));
            let temp_output_path = final_outputs.temp_path_for(&output_path)?;
            write_npy(&temp_output_path, &freqs)
                .with_context(|| format!("writing {}", temp_output_path.display()))?;
            final_outputs.record(temp_output_path, output_path.clone())?;
            output_files.push(output_path);

            // Copy motif metadata so downstream tools share identical ordering
            let motif_dest =
                final_dir.join(dot_join(&[prefix, &format!("k{k}_{group}_motifs.txt")]));
            let temp_motif_dest = final_outputs.temp_path_for(&motif_dest)?;
            fs::copy(&motifs_path, &temp_motif_dest)
                .with_context(|| format!("copying to {}", temp_motif_dest.display()))?;
            final_outputs.record(temp_motif_dest, motif_dest.clone())?;
            output_files.push(motif_dest);

            wrote_group = true;
        }
        ensure!(
            wrote_group,
            "no positional counts found for k={} when computing order {} transitions",
            k,
            order
        );
    }

    // Position files describe the left/right/mid coordinate grids and are not order-specific
    // Copy each existing group file once after all transition arrays have been recorded
    for group in groups {
        let positions_file_name = dot_join(&[prefix, &format!("{group}_positions.txt")]);
        let positions_src = counts_dir.join(&positions_file_name);
        if positions_src.exists() {
            let positions_dest = final_dir.join(&positions_file_name);
            let temp_positions_dest = final_outputs.temp_path_for(&positions_dest)?;
            fs::copy(&positions_src, &temp_positions_dest)
                .with_context(|| format!("copying to {}", temp_positions_dest.display()))?;
            final_outputs.record(temp_positions_dest, positions_dest.clone())?;
            output_files.push(positions_dest);
        }
    }

    final_outputs.move_into_place()?;

    // Remove temporary staging directory once final outputs are written
    temp_root_guard
        .remove()
        .context("remove transitions temp directory")?;

    let elapsed = start_time.elapsed();
    if options.report_statistics {
        print_fragment_run_statistics(
            &global_counter.base,
            elapsed,
            FragmentRunStatisticsOptions {
                include_section_header: true,
                notes: &[TILE_DOUBLE_COUNT_NOTE],
                labels: DEFAULT_FRAGMENT_STATISTICS_LABELS,
                blacklist_excluded_fragments: Some(global_counter.blacklisted_fragments),
                gc: (opt.shared_args.gc.gc_file.is_some() || opt.shared_args.gc.gc_tag.is_some())
                    .then_some(GCStatisticsSummary {
                        neutralize_invalid_gc: opt.shared_args.gc.neutralize_invalid_gc,
                        failed_fragments: global_counter.gc_failed_fragments,
                        missing_tags: opt
                            .shared_args
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(global_counter.gc_missing_tags),
                        out_of_range_tags: opt
                            .shared_args
                            .gc
                            .gc_tag
                            .is_some()
                            .then_some(global_counter.gc_out_of_range_tags),
                    }),
            },
            std::iter::empty::<&str>(),
        );
    }
    Ok(TransitionsRunResult {
        counters: global_counter,
        output_files,
    })
}

#[cfg(test)]
mod tests {
    include!("transitions_tests.rs");
}
