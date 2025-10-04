use crate::commands::prepare_windows::{
    config::{MergeScope, PrepareConfig},
    mergers::merge_windows,
    order::{WindowSortOrder, sort_windows_in_place},
    postprocess::{
        deduplicate_identical, enforce_min_distance_within_group, partition_safe_and_tail,
    },
    prepare_windows::{BlacklistCursor, FinalWindow},
    writers::{ChromTempWriter, ensure_temp_writer_for_chrom, write_windows},
};
use crate::shared::blacklist::{BlacklistStrategy, is_blacklisted};
use anyhow::Result;
use fxhash::FxHashMap;
use std::path::Path;

// TODO: Rename "safe prefix" to something understandable
/// Process a chunk for a chromosome, writing the safe prefix to disk.
///
/// Combines the previous tail with the current batch, runs post-processing
/// (deduplication, spacing, merging), writes the safe portion to the chromosome temp
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
/// - `blacklist_cursor`: streaming cursor for blacklist intervals on this chromosome.
///   **NOTE**: Should be reset and ready for a second stream-through.
/// - `blacklist_look_back`: halo applied when checking blacklist overlap.
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
    blacklist_cursor: Option<&mut BlacklistCursor>,
    blacklist_look_back: u64,
    cfg: &PrepareConfig,
) -> Result<()> {
    let mut windows: Vec<FinalWindow> =
        Vec::with_capacity(carryover_tail.len() + current_batch.len());
    windows.append(carryover_tail);
    windows.append(current_batch);

    // Sort by either group-first or group-last
    let mut current_sort = if cfg.min_distance_within_group.is_some()
        || !matches!(cfg.merge_scope, MergeScope::Within)
    {
        // Order by `(group, chrom, start, end)`
        sort_windows_in_place(&mut windows, WindowSortOrder::GroupChromStartEnd);
        WindowSortOrder::GroupChromStartEnd
    } else {
        // Order by `(chrom, start, end, group)`
        sort_windows_in_place(&mut windows, WindowSortOrder::ChromStartEndGroup);
        WindowSortOrder::ChromStartEndGroup
    };

    // Remove duplicates
    // Works with both group-first and group-last order
    let mut windows = deduplicate_identical(windows, cfg.deduplicate, cfg.score_col.is_some());

    // TODO: Ensure we properly document that this happens before merging? Should it??
    // Spacing operates on group-first ordering
    // Has early return when min_distance_within_group is `None`
    windows = enforce_min_distance_within_group(
        windows,
        cfg.min_distance_within_group,
        cfg.distance_ties,
        cfg.score_col.is_some(),
    );

    // Prepare ordering for merge stage. Within-group merging reuses the existing
    // grouping order, whereas across-group merging requires genomic sorting
    let windows = match (cfg.merge_scope, current_sort) {
        // Already sorted properly for merging or downstream
        (MergeScope::None, WindowSortOrder::ChromStartEndGroup)
        | (MergeScope::Across, WindowSortOrder::ChromStartEndGroup)
        | (MergeScope::Within, WindowSortOrder::GroupChromStartEnd) => windows,
        // Sort to group-last
        (MergeScope::Across, WindowSortOrder::GroupChromStartEnd)
        | (MergeScope::None, WindowSortOrder::GroupChromStartEnd) => {
            sort_windows_in_place(&mut windows, WindowSortOrder::ChromStartEndGroup);
            current_sort = WindowSortOrder::ChromStartEndGroup;
            windows
        }
        // Sort to group-first
        (MergeScope::Within, WindowSortOrder::ChromStartEndGroup) => {
            sort_windows_in_place(&mut windows, WindowSortOrder::GroupChromStartEnd);
            current_sort = WindowSortOrder::GroupChromStartEnd;
            windows
        }
    };

    let mut merged_windows = merge_windows(
        windows,
        cfg.merge_scope,
        cfg.merge_gap,
        cfg.merge_label,
        true,
    );

    // Sort group-last if not already the case (merging should not affect order)
    if !matches!(current_sort, WindowSortOrder::ChromStartEndGroup) {
        sort_windows_in_place(&mut merged_windows, WindowSortOrder::ChromStartEndGroup);
        // current_sort = WindowSortOrder::ChromStartEndGroup;
    }

    let windows = filter_blacklisted_post_merge(
        merged_windows,
        blacklist_cursor,
        cfg.blacklist_strategy,
        blacklist_look_back,
    );

    let (safe_prefix, tail) = partition_safe_and_tail(
        windows,
        cfg.min_distance_within_group,
        cfg.merge_scope,
        cfg.merge_gap,
    );

    let writer = ensure_temp_writer_for_chrom(chrom, temp_dir, temp_writers)?;

    write_windows(writer.writer(), &safe_prefix, cfg.sep)?;
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
/// - `blacklist_cursor`: streaming cursor for blacklist intervals on this chromosome.
/// - `blacklist_look_back`: halo applied when checking blacklist overlap.
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
    blacklist_cursor: Option<&mut BlacklistCursor>,
    blacklist_look_back: u64,
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
        blacklist_cursor,
        blacklist_look_back,
        cfg,
    )?;

    if !carryover_tail.is_empty() {
        let writer = ensure_temp_writer_for_chrom(chrom, temp_dir, temp_writers)?;
        sort_windows_in_place(carryover_tail, WindowSortOrder::ChromStartEndGroup);
        write_windows(writer.writer(), carryover_tail, cfg.sep)?;
        for w in carryover_tail.drain(..) {
            if !w.group.is_empty() {
                *global_group_counts.entry(w.group.clone()).or_insert(0) += 1;
            }
        }
    }
    Ok(())
}

// TODO: Add docstring
/// windows should be sorted by chrom,start,end
fn filter_blacklisted_post_merge(
    windows: Vec<FinalWindow>,
    blacklist_cursor: Option<&mut BlacklistCursor>,
    strategy: BlacklistStrategy,
    look_back: u64,
) -> Vec<FinalWindow> {
    match blacklist_cursor {
        Some(cursor) if !cursor.intervals.is_empty() => {
            let intervals = cursor.intervals.as_slice();
            let mut retained: Vec<FinalWindow> = Vec::with_capacity(windows.len());
            for entry in windows.into_iter() {
                if entry.merged {
                    if is_blacklisted(
                        intervals,
                        strategy,
                        entry.start as u64,
                        entry.end as u64,
                        look_back,
                        &mut cursor.post_cursor,
                    ) {
                        continue;
                    }
                }
                retained.push(entry);
            }
            retained
        }
        _ => windows,
    }
}
