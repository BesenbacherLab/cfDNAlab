use crate::commands::prepare_windows::{
    config::PrepareConfig,
    mergers::merge_windows,
    postprocess::{
        deduplicate_identical, enforce_min_distance_within_group, partition_safe_and_tail,
    },
    prepare_windows::FinalWindow,
    writers::{ChromTempWriter, ensure_temp_writer_for_chrom, write_windows},
};
use anyhow::Result;
use fxhash::FxHashMap;
use std::path::Path;

/// Process a chunk for a chromosome, writing the safe prefix to disk.
///
/// Combines the previous tail with the current batch, runs post-processing
/// (dedupe, spacing, merging), writes the safe portion to the chromosome temp
/// file, and retains a boundary tail for the next chunk.
///
/// Sorting keeps behavior deterministic, and every emitted window with a group
/// label increments the `global_group_counts` map used later for
/// `min_per_group` filtering.
///
/// # Parameters
/// - `chrom`: chromosome identifier being processed.
/// - `carryover_tail`: tail windows kept from the prior chunk.
/// - `current_batch`: windows accumulated for this chunk.
/// - `temp_writers`: map of chromosome temp writers.
/// - `temp_dir`: base directory for temp files.
/// - `global_group_counts`: running group counts.
/// - `cfg`: user configuration controlling post-processing.
///
/// # Returns
/// `Ok(())` on success or an error if writing fails.
pub fn process_and_write_chunk(
    chrom: &str,
    carryover_tail: &mut Vec<FinalWindow>,
    current_batch: &mut Vec<FinalWindow>,
    temp_writers: &mut FxHashMap<String, ChromTempWriter>,
    temp_dir: &Path,
    global_group_counts: &mut FxHashMap<String, u64>,
    cfg: &PrepareConfig,
) -> Result<()> {
    let mut windows: Vec<FinalWindow> =
        Vec::with_capacity(carryover_tail.len() + current_batch.len());
    windows.append(carryover_tail);
    windows.append(current_batch);

    windows.sort_unstable_by(|a, b| {
        a.group
            .cmp(&b.group)
            .then(a.chrom.cmp(&b.chrom))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    let windows = deduplicate_identical(windows, cfg.deduplicate, cfg.score_col.is_some());

    let windows = enforce_min_distance_within_group(
        windows,
        cfg.min_distance_within_group,
        cfg.distance_ties,
        cfg.score_col.is_some(),
    );

    let windows = merge_windows(windows, cfg.merge_scope, cfg.merge_gap, cfg.merge_label);

    let (safe_prefix, tail) = partition_safe_and_tail(
        windows,
        cfg.min_distance_within_group,
        cfg.merge_scope,
        cfg.merge_gap,
    );

    let writer = ensure_temp_writer_for_chrom(chrom, temp_dir, temp_writers)?;

    write_windows(writer.writer(), &safe_prefix, cfg.separator)?;
    for w in &safe_prefix {
        if !w.group.is_empty() {
            *global_group_counts.entry(w.group.clone()).or_insert(0) += 1;
        }
    }

    *carryover_tail = tail;
    Ok(())
}

/// Flush remaining windows when finishing a chromosome stream.
///
/// Invokes [`process_and_write_chunk`] one last time and then writes any
/// residual tail because no future chunk can modify it.
///
/// After processing, the tail is emitted directly to the chromosome writer and
/// group counts are updated for each window.
///
/// # Parameters
/// - `chrom`: chromosome identifier.
/// - `carryover_tail`: tail windows kept from prior chunks.
/// - `current_batch`: windows accumulated for the final chunk.
/// - `temp_writers`: map of chromosome writers.
/// - `temp_dir`: temporary directory containing writer files.
/// - `global_group_counts`: running group counts.
/// - `cfg`: configuration controlling post-processing.
///
/// # Returns
/// `Ok(())` on success or an error if writing fails.
pub fn flush_chromosome(
    chrom: &str,
    carryover_tail: &mut Vec<FinalWindow>,
    current_batch: &mut Vec<FinalWindow>,
    temp_writers: &mut FxHashMap<String, ChromTempWriter>,
    temp_dir: &Path,
    global_group_counts: &mut FxHashMap<String, u64>,
    cfg: &PrepareConfig,
) -> Result<()> {
    if carryover_tail.is_empty() && current_batch.is_empty() {
        return Ok(());
    }
    process_and_write_chunk(
        chrom,
        carryover_tail,
        current_batch,
        temp_writers,
        temp_dir,
        global_group_counts,
        cfg,
    )?;

    if !carryover_tail.is_empty() {
        let writer = ensure_temp_writer_for_chrom(chrom, temp_dir, temp_writers)?;
        write_windows(writer.writer(), carryover_tail, cfg.separator)?;
        for w in carryover_tail.drain(..) {
            if !w.group.is_empty() {
                *global_group_counts.entry(w.group.clone()).or_insert(0) += 1;
            }
        }
    }
    Ok(())
}
