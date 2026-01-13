use crate::commands::prepare_windows::{
    config::{CoordinateSet, DedupKeep, DistSign, MergeScope, NearTiePolicy, PrepareConfig},
    intermediate::write_intermediate_windows,
    labels::{
        AtomicLabelPart, LabelKey, LabelSchema, LabelTuple, NO_NEAR_BIN_LABEL, NO_NEAR_LABEL,
        build_tuple_compositions, normalize_label_tuples, render_label_for_key,
    },
    mergers::merge_windows,
    near_file::{NearHit, NearIndex, NearSide, NearTie, NearestResult, nearest_edge_distance},
    order::{WindowSortOrder, sort_windows_in_place},
    parsers::DistanceBins,
    postprocess::{
        apply_cluster_labels, deduplicate_identical, enforce_min_distance_within_group,
        partition_safe_and_tail,
    },
    prepare_windows::{BlacklistCursor, Window},
    resizers::apply_size_transform,
    writers::{ChromTempWriter, ensure_temp_writer_for_chrom},
};
use crate::shared::blacklist::{BlacklistStrategy, is_blacklisted};
use anyhow::Result;
use fxhash::FxHashMap;
use std::path::Path;

/// Per-hit annotation values used to expand label tuples
#[derive(Clone)]
struct NearAnnotation {
    near_side: String,
    near_name: Option<String>,
    bin: Option<String>,
}

// Tab is disallowed in labels, so it is safe as a stable separator for sort keys
const OUTPUT_SORT_SEPARATOR: char = '\t';

impl NearAnnotation {
    // Keep this as a small helper to avoid repeating lookup and filtering logic
    fn from_hit(
        hit: &NearHit,
        near_index: &NearIndex,
        distance_bins: Option<&DistanceBins>,
        distance_sign: DistSign,
        distance_max: Option<u32>,
    ) -> Option<Self> {
        if let Some(max_abs) = distance_max {
            if hit.distance.unsigned_abs() > max_abs {
                return None;
            }
        }

        let distance_for_bin = match distance_sign {
            DistSign::Absolute => hit.distance.abs(),
            DistSign::Signed => hit.distance,
        };

        let bin = distance_bins.and_then(|bins| bins.match_label(distance_for_bin));
        let near_side = match hit.side {
            NearSide::Upstream => "-",
            NearSide::Downstream => "+",
            NearSide::Overlap => "=",
        }
        .to_string();
        let near_name = hit
            .group_id
            .map(|id| near_index.group_id_to_name[id as usize].clone());

        Some(Self {
            near_side,
            near_name,
            bin: bin.map(str::to_string),
        })
    }

    fn push_from_hit(
        annotations: &mut Vec<Self>,
        hit: &NearHit,
        near_index: &NearIndex,
        distance_bins: Option<&DistanceBins>,
        distance_sign: DistSign,
        distance_max: Option<u32>,
    ) {
        if let Some(annotation) =
            Self::from_hit(hit, near_index, distance_bins, distance_sign, distance_max)
        {
            annotations.push(annotation);
        }
    }
}

/// Process a chunk for a chromosome, writing the processed region to disk.
///
/// Combines the previous tail with the current batch, runs post-processing
/// (merging, deduplication, minimum-distance filtering), writes the processed region to the
/// chromosome temp file, and retains a boundary tail for the next chunk.
///
/// The processed region contains windows that cannot be affected by any future chunk,
/// so they can be written immediately. Later min-per or exclusion filters can still
/// remove them in the final pass.
///
/// Sorting keeps behavior deterministic.
///
/// # Parameters
/// - `chrom`: chromosome identifier being processed.
/// - `carryover_tail`: tail windows kept from the prior chunk.
/// - `current_batch`: windows accumulated for this chunk.
/// - `temp_writers`: map of chromosome temp writers.
/// - `temp_dir`: base directory for temp files.
/// - `blacklist_cursor`: streaming cursor for blacklist intervals on this chromosome.
///   **NOTE**: Should be reset and ready for a second stream-through.
/// - `blacklist_look_back`: halo applied when checking blacklist overlap.
/// - `chrom_size`: chromosome size for resize and flank logic.
/// - `cfg`: user configuration controlling post-processing.
/// - `near_index`: optional near interval index used for distance labeling.
/// - `distance_bins`: optional distance bin rules.
/// - `label_schema`: resolved label compositions.
/// - `merge_key`: label key used for merging.
/// - `out_labels`: label keys that define output ordering and columns.
///
/// # Returns
/// `Ok(())` on success or an error if writing fails.
pub fn process_and_write_chunk(
    chrom: &str,
    carryover_tail: &mut Vec<Window>,
    current_batch: &mut Vec<Window>,
    temp_writers: &mut FxHashMap<String, ChromTempWriter>,
    temp_dir: &Path,
    blacklist_cursor: Option<&mut BlacklistCursor>,
    blacklist_look_back: u64,
    chrom_size: Option<u32>,
    cfg: &PrepareConfig,
    near_index: &mut Option<NearIndex>,
    distance_bins: Option<&DistanceBins>,
    label_schema: &LabelSchema,
    merge_key: &LabelKey,
    out_labels: &[LabelKey],
) -> Result<()> {
    let mut windows: Vec<Window> = Vec::with_capacity(carryover_tail.len() + current_batch.len());
    windows.append(carryover_tail);
    windows.append(current_batch);

    let has_size_transform = cfg.resize.is_some() || cfg.flank.is_some();

    // Coordinate sets to use per transformation
    let dedup_coord_set = if has_size_transform {
        CoordinateSet::Resized
    } else {
        CoordinateSet::Original
    };
    let merge_coord_set = cfg.merge_on;
    let distance_coord_set = cfg.distance_from;
    let cluster_coord_set = cfg.cluster_on;
    let output_coord_set = CoordinateSet::Resized;

    let input_key = LabelKey::Atomic(AtomicLabelPart::Input);
    // Tracks whether window.group_key currently reflects the input label
    let mut is_input_keyed = false;

    // Merging requires both a scope and a gap threshold
    let merge_within_enabled =
        cfg.merge_gap.is_some() && matches!(cfg.merge_scope, MergeScope::Within);
    let merge_across_enabled =
        cfg.merge_gap.is_some() && matches!(cfg.merge_scope, MergeScope::Across);
    let dedup_enabled = !matches!(cfg.deduplicate, DedupKeep::None);

    // Track current ordering and coordinate set to avoid unnecessary resorting
    let mut current_order: Option<WindowSortOrder> = None;
    let mut current_coord_set: Option<CoordinateSet> = None;

    let mut windows = if dedup_enabled {
        // Deduplication keys on the input label and accepts either group-first or chrom-first ordering
        update_group_keys(&mut windows, &input_key, label_schema);
        is_input_keyed = true;

        let prefer_group_first = if merge_within_enabled {
            merge_key == &input_key && merge_coord_set == dedup_coord_set
        } else {
            cfg.min_distance_within_group.is_some() && !cfg.cluster_before_min_distance
        };
        let dedup_order = if prefer_group_first {
            WindowSortOrder::GroupChromStartEnd
        } else {
            WindowSortOrder::ChromStartEndGroup
        };
        ensure_sorted(
            &mut windows,
            dedup_order,
            dedup_coord_set,
            &mut current_order,
            &mut current_coord_set,
        );
        deduplicate_identical(
            windows,
            cfg.deduplicate,
            cfg.score_col.is_some(),
            dedup_coord_set,
        )
    } else {
        windows
    };

    if merge_within_enabled {
        let merge_key_is_input = merge_key == &input_key;
        // Merge grouping can be different from input
        if !merge_key_is_input {
            update_group_keys(&mut windows, merge_key, label_schema);
            current_order = None;
            is_input_keyed = false;
        }

        // Prepare ordering for the within-group merge pass
        ensure_sorted(
            &mut windows,
            WindowSortOrder::GroupChromStartEnd,
            merge_coord_set,
            &mut current_order,
            &mut current_coord_set,
        );

        windows = merge_windows(
            windows,
            MergeScope::Within,
            cfg.merge_gap,
            cfg.merge_label,
            cfg.merge_on,
            true,
        );

        if matches!(cfg.merge_on, CoordinateSet::Original) && has_size_transform {
            let mut resized_after_merge: Vec<Window> = Vec::with_capacity(windows.len());
            for mut window in windows {
                // NOTE: Falls back to original coordinates when no resizing is specified
                if let Some((resized_start, resized_end)) = apply_size_transform(
                    window.original_start,
                    window.original_end,
                    chrom_size,
                    cfg,
                ) {
                    window.resized_start = resized_start;
                    window.resized_end = resized_end;
                    resized_after_merge.push(window);
                }
            }
            windows = resized_after_merge;
        }
    }

    let cluster_before_min_distance = cfg.cluster_before_min_distance;
    if cluster_before_min_distance {
        if let Some(min_overlaps) = cfg.cluster_min_overlaps {
            // Cluster labels are based on overlap depth across groups within a chromosome
            ensure_sorted(
                &mut windows,
                WindowSortOrder::ChromStartEnd,
                cluster_coord_set,
                &mut current_order,
                &mut current_coord_set,
            );
            apply_cluster_labels(&mut windows, min_overlaps, cluster_coord_set);
        }
    }

    if cfg.min_distance_within_group.is_some() {
        // Ensure input-based grouping before minimum-distance filtering
        if !is_input_keyed {
            update_group_keys(&mut windows, &input_key, label_schema);
            current_order = None;
            //is_input_keyed = true;
        }
        // Minimum-distance filtering uses input groups and group-first ordering
        ensure_sorted(
            &mut windows,
            WindowSortOrder::GroupChromStartEnd,
            distance_coord_set,
            &mut current_order,
            &mut current_coord_set,
        );

        windows = enforce_min_distance_within_group(
            windows,
            cfg.min_distance_within_group,
            cfg.distance_policy,
            cfg.score_col.is_some(),
            distance_coord_set,
        );
    }

    if !cluster_before_min_distance {
        if let Some(min_overlaps) = cfg.cluster_min_overlaps {
            // Cluster labels are based on overlap depth across groups within a chromosome
            ensure_sorted(
                &mut windows,
                WindowSortOrder::ChromStartEnd,
                cluster_coord_set,
                &mut current_order,
                &mut current_coord_set,
            );
            apply_cluster_labels(&mut windows, min_overlaps, cluster_coord_set);
        }
    }

    if merge_across_enabled {
        // Across-group merges happen after minimum-distance filtering and clustering so cluster labels
        // reflect overlap density between original groups
        ensure_sorted(
            &mut windows,
            WindowSortOrder::ChromStartEndGroup,
            merge_coord_set,
            &mut current_order,
            &mut current_coord_set,
        );

        windows = merge_windows(
            windows,
            MergeScope::Across,
            cfg.merge_gap,
            cfg.merge_label,
            cfg.merge_on,
            true,
        );

        if matches!(cfg.merge_on, CoordinateSet::Original) && has_size_transform {
            let mut resized_after_merge: Vec<Window> = Vec::with_capacity(windows.len());
            for mut window in windows {
                // NOTE: Falls back to original coordinates when no resizing is specified
                if let Some((resized_start, resized_end)) = apply_size_transform(
                    window.original_start,
                    window.original_end,
                    chrom_size,
                    cfg,
                ) {
                    window.resized_start = resized_start;
                    window.resized_end = resized_end;
                    resized_after_merge.push(window);
                }
            }
            windows = resized_after_merge;
        }
    }

    // Blacklist checks expect chrom-first ordering on resized coordinates
    ensure_sorted(
        &mut windows,
        WindowSortOrder::ChromStartEnd,
        output_coord_set,
        &mut current_order,
        &mut current_coord_set,
    );

    let mut windows = filter_blacklisted_post_merge(
        windows,
        blacklist_cursor,
        cfg.blacklist_strategy,
        blacklist_look_back,
        output_coord_set,
    );

    // Tail detection for the next chunk (and min-distance safety) needs chrom-first ordering on the distance coordinate set
    ensure_sorted(
        &mut windows,
        WindowSortOrder::ChromStartEnd,
        distance_coord_set,
        &mut current_order,
        &mut current_coord_set,
    );

    let merge_group_keys = if merge_within_enabled {
        // Tail detection for within-group merging needs merge-key values, not output labels
        Some(build_group_keys_for_label(
            &windows,
            label_schema,
            merge_key,
        ))
    } else {
        None
    };

    // Debug: surface the largest window in this chunk before partitioning
    if let Some(max_window) = windows.iter().max_by_key(|w| w.length_for(CoordinateSet::Original)) {
        eprintln!(
            "Debug: largest window before partition {}:{}-{} len={} group={}",
            max_window.chrom,
            max_window.start_for(CoordinateSet::Original),
            max_window.end_for(CoordinateSet::Original),
            max_window.length_for(CoordinateSet::Original),
            max_window.group_key
        );
    }

    // Split into a processed region that cannot change and a tail that might still merge with the next chunk
    let (mut safe_prefix, tail) = partition_safe_and_tail(
        windows,
        cfg.min_distance_within_group,
        cfg.merge_scope,
        cfg.merge_gap,
        cfg.merge_on,
        distance_coord_set,
        cfg.cluster_min_overlaps,
        cluster_coord_set,
        merge_group_keys.as_deref(),
    );

    println!(
        "Safe length: {} | Tail length: {}",
        safe_prefix.len(),
        tail.len()
    );

    // Only annotate the safe prefix now. Tail will be annotated on final flush
    // Note: Uses the existing ChromStartEnd ordering
    safe_prefix = apply_near_annotations(
        safe_prefix,
        near_index,
        cfg,
        distance_bins,
        distance_coord_set,
    );

    // Output ordering uses the label columns requested by the user, not group_key
    sort_windows_by_output_labels(&mut safe_prefix, output_coord_set, out_labels, label_schema);

    let writer = ensure_temp_writer_for_chrom(chrom, temp_dir, temp_writers)?;

    write_intermediate_windows(writer.writer(), &safe_prefix, cfg.sep)?;
    *carryover_tail = tail;
    Ok(())
}

/// Flush remaining windows when finishing a chromosome stream.
///
/// Invokes [`process_and_write_chunk`] one last time and then writes any
/// residual tail because no future chunk can modify it.
///
/// After processing, the tail is send directly to the chromosome writer.
///
/// # Parameters
/// - `chrom`: chromosome identifier.
/// - `carryover_tail`: tail windows kept from prior chunks.
/// - `current_batch`: windows accumulated for the final chunk.
/// - `temp_writers`: map of chromosome writers.
/// - `temp_dir`: temporary directory containing writer files.
/// - `blacklist_cursor`: streaming cursor for blacklist intervals on this chromosome.
/// - `blacklist_look_back`: halo applied when checking blacklist overlap.
/// - `chrom_size`: chromosome size for resize and flank logic.
/// - `cfg`: configuration controlling post-processing.
/// - `near_index`: optional near interval index used for distance labeling.
/// - `distance_bins`: optional distance bin rules.
/// - `label_schema`: resolved label compositions.
/// - `merge_key`: label key used for merging.
/// - `out_labels`: label keys that define output ordering and columns.
///
/// # Returns
/// `Ok(())` on success or an error if writing fails.
pub fn flush_chromosome(
    chrom: &str,
    carryover_tail: &mut Vec<Window>,
    current_batch: &mut Vec<Window>,
    temp_writers: &mut FxHashMap<String, ChromTempWriter>,
    temp_dir: &Path,
    blacklist_cursor: Option<&mut BlacklistCursor>,
    blacklist_look_back: u64,
    chrom_size: Option<u32>,
    cfg: &PrepareConfig,
    near_index: &mut Option<NearIndex>,
    distance_bins: Option<&DistanceBins>,
    label_schema: &LabelSchema,
    merge_key: &LabelKey,
    out_labels: &[LabelKey],
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
        blacklist_cursor,
        blacklist_look_back,
        chrom_size,
        cfg,
        near_index,
        distance_bins,
        label_schema,
        merge_key,
        out_labels,
    )?;

    if !carryover_tail.is_empty() {
        let mut tail = std::mem::take(carryover_tail);
        tail = apply_near_annotations(tail, near_index, cfg, distance_bins, cfg.distance_from);
        sort_windows_by_output_labels(&mut tail, CoordinateSet::Resized, out_labels, label_schema);

        let writer = ensure_temp_writer_for_chrom(chrom, temp_dir, temp_writers)?;
        write_intermediate_windows(writer.writer(), &tail, cfg.sep)?;
    }
    Ok(())
}

fn update_group_keys(windows: &mut [Window], merge_key: &LabelKey, label_schema: &LabelSchema) {
    let needs_compositions = matches!(merge_key, LabelKey::Composition(_));
    for window in windows {
        if window.label_tuples.is_empty() {
            window.group_key.clear();
            continue;
        }
        let tuple_compositions = if needs_compositions {
            build_tuple_compositions(&window.label_tuples, label_schema)
        } else {
            Vec::new()
        };
        window.group_key = render_label_for_key(
            &window.label_tuples,
            &tuple_compositions,
            merge_key,
            label_schema,
        );
    }
}

fn build_group_keys_for_label(
    windows: &[Window],
    label_schema: &LabelSchema,
    key: &LabelKey,
) -> Vec<String> {
    // Precompute merge-group labels because output ordering does not track merge-key ordering
    let needs_compositions = matches!(key, LabelKey::Composition(_));
    let mut keys = Vec::with_capacity(windows.len());
    for window in windows {
        let tuple_compositions = if needs_compositions {
            build_tuple_compositions(&window.label_tuples, label_schema)
        } else {
            Vec::new()
        };
        let label =
            render_label_for_key(&window.label_tuples, &tuple_compositions, key, label_schema);
        keys.push(label);
    }
    keys
}

fn sort_windows_by_output_labels(
    windows: &mut [Window],
    coord_set: CoordinateSet,
    out_labels: &[LabelKey],
    label_schema: &LabelSchema,
) {
    // Output order is chrom, coordinates, then label columns in the user-specified order
    if windows.len() <= 1 {
        return;
    }

    let needs_compositions = out_labels
        .iter()
        .any(|key| matches!(key, LabelKey::Composition(_)));

    windows.sort_by_cached_key(|window| {
        // Cache rendered label values to avoid recomputing during comparisons
        let tuple_compositions = if needs_compositions {
            build_tuple_compositions(&window.label_tuples, label_schema)
        } else {
            Vec::new()
        };
        let mut label_key = String::new();
        for (idx, key) in out_labels.iter().enumerate() {
            if idx > 0 {
                label_key.push(OUTPUT_SORT_SEPARATOR);
            }
            let label =
                render_label_for_key(&window.label_tuples, &tuple_compositions, key, label_schema);
            label_key.push_str(&label);
        }
        (
            window.chrom.clone(),
            window.start_for(coord_set),
            window.end_for(coord_set),
            label_key,
        )
    });
}

/// Ensure windows are sorted in the required order and coordinate space.
#[inline]
fn ensure_sorted(
    windows: &mut [Window],
    desired_order: WindowSortOrder,
    desired_coord_set: CoordinateSet,
    current_order: &mut Option<WindowSortOrder>,
    current_coord_set: &mut Option<CoordinateSet>,
) {
    let order_matches = match (*current_order, desired_order) {
        (Some(current), desired) if current == desired => true,
        // Group-tiebreak ordering still satisfies pure chrom/start/end ordering
        (Some(WindowSortOrder::ChromStartEndGroup), WindowSortOrder::ChromStartEnd) => true,
        _ => false,
    };
    let coord_matches = current_coord_set.map_or(false, |coord| coord == desired_coord_set);
    if !order_matches || !coord_matches {
        sort_windows_in_place(windows, desired_order, desired_coord_set);
        *current_order = Some(desired_order);
        *current_coord_set = Some(desired_coord_set);
    }
}

pub fn apply_near_annotations(
    windows: Vec<Window>,
    near_index: &mut Option<NearIndex>,
    cfg: &PrepareConfig,
    distance_bins: Option<&DistanceBins>,
    coord_set: CoordinateSet,
) -> Vec<Window> {
    let Some(near_idx) = near_index.as_mut() else {
        return windows;
    };

    // Rewind the near cursor for this chunk's chromosome using existing cursor state.
    if let Some(first_window) = windows.first() {
        let min_start = first_window.start_for(coord_set);
        let chunk_chrom = first_window.chrom.clone();
        if let Some(chrom) = near_idx.per_chrom.get_mut(chunk_chrom.as_ref()) {
            if !chrom.intervals.is_empty() {
                let len = chrom.intervals.len();
                let mut cursor = chrom.cursor.min(len.saturating_sub(1));
                if cursor != 0 {
                    while cursor > 0 && chrom.intervals[cursor].end > min_start {
                        cursor -= 1;
                    }
                }
                chrom.cursor = cursor;
            }
        }
    }

    let is_signed_mode = matches!(cfg.distance_sign, DistSign::Signed);

    let mut retained: Vec<Window> = Vec::with_capacity(windows.len());

    for mut window in windows {
        let chrom = window.chrom.as_ref();
        let near_chrom = near_idx.per_chrom.get_mut(chrom);
        let no_near_intervals = near_chrom
            .as_ref()
            .map(|chrom| chrom.intervals.is_empty())
            .unwrap_or(true);
        if no_near_intervals {
            let first_warning = near_idx.warned_no_near.is_empty();
            let should_warn = near_idx.warned_no_near.insert(chrom.to_string());
            if should_warn {
                if first_warning {
                    let include_name = !cfg.near_group_cols.is_empty();
                    if cfg.distance_max.is_some() {
                        eprintln!(
                            "Warning: Chromosome '{}' has no near intervals. Windows on this chromosome will be dropped due to --distance-max.",
                            chrom
                        );
                    } else if distance_bins.is_some() {
                        if include_name {
                            eprintln!(
                                "Warning: Chromosome '{}' has no near intervals. Windows will keep near-side/near-name as '{}' and bin as '{}'.",
                                chrom, NO_NEAR_LABEL, NO_NEAR_BIN_LABEL
                            );
                        } else {
                            eprintln!(
                                "Warning: Chromosome '{}' has no near intervals. Windows will keep near-side as '{}' and bin as '{}'.",
                                chrom, NO_NEAR_LABEL, NO_NEAR_BIN_LABEL
                            );
                        }
                    } else if include_name {
                        eprintln!(
                            "Warning: Chromosome '{}' has no near intervals. Windows will keep near-side/near-name as '{}'.",
                            chrom, NO_NEAR_LABEL
                        );
                    } else {
                        eprintln!(
                            "Warning: Chromosome '{}' has no near intervals. Windows will keep near-side as '{}'.",
                            chrom, NO_NEAR_LABEL
                        );
                    }
                } else {
                    eprintln!("Warning: Chromosome '{}' has no near intervals.", chrom);
                }
            }

            if cfg.distance_max.is_some() {
                continue;
            }

            for tuple in &mut window.label_tuples {
                tuple.near_side = Some(NO_NEAR_LABEL.to_string());
                if !cfg.near_group_cols.is_empty() {
                    tuple.near_name = Some(NO_NEAR_LABEL.to_string());
                }
                if distance_bins.is_some() {
                    tuple.bin = Some(NO_NEAR_BIN_LABEL.to_string());
                }
            }

            retained.push(window);
            continue;
        }
        let near_chrom = near_chrom.expect("near_chrom should exist");

        let window_start = window.start_for(coord_set);
        let window_end = window.end_for(coord_set);

        let Some(nearest_result) = nearest_edge_distance(
            window_start,
            window_end,
            near_chrom,
            &cfg.near_edge,
            &cfg.near_direction,
            is_signed_mode,
        ) else {
            retained.push(window);
            continue;
        };

        let mut annotations: Vec<NearAnnotation> = Vec::new();

        match nearest_result {
            NearestResult::Single(hit) => {
                NearAnnotation::push_from_hit(
                    &mut annotations,
                    &hit,
                    near_idx,
                    distance_bins,
                    cfg.distance_sign,
                    cfg.distance_max,
                );
            }
            NearestResult::Tie(NearTie {
                upstream,
                downstream,
            }) => {
                if matches!(cfg.near_ties, NearTiePolicy::Drop) {
                    continue;
                }
                // Keep both sides when ties are allowed, so downstream filters can decide
                if let Some(hit) = upstream.as_ref() {
                    NearAnnotation::push_from_hit(
                        &mut annotations,
                        hit,
                        near_idx,
                        distance_bins,
                        cfg.distance_sign,
                        cfg.distance_max,
                    );
                }
                if let Some(hit) = downstream.as_ref() {
                    NearAnnotation::push_from_hit(
                        &mut annotations,
                        hit,
                        near_idx,
                        distance_bins,
                        cfg.distance_sign,
                        cfg.distance_max,
                    );
                }
            }
        }

        if annotations.is_empty() {
            continue;
        }

        if annotations.len() == 1 {
            // Apply a single near label to every tuple without duplicating the window
            let annotation = annotations.pop().unwrap();
            for tuple in &mut window.label_tuples {
                tuple.near_side = Some(annotation.near_side.clone());
                tuple.near_name = annotation.near_name.clone();
                tuple.bin = annotation.bin.clone();
            }
        } else {
            // Expand tuples so each label tuple keeps one near annotation
            let mut expanded: Vec<LabelTuple> =
                Vec::with_capacity(window.label_tuples.len() * annotations.len());
            for tuple in window.label_tuples.iter() {
                for annotation in &annotations {
                    let mut updated = tuple.clone();
                    updated.near_side = Some(annotation.near_side.clone());
                    updated.near_name = annotation.near_name.clone();
                    updated.bin = annotation.bin.clone();
                    expanded.push(updated);
                }
            }
            window.label_tuples = expanded;
        }

        normalize_label_tuples(&mut window.label_tuples);
        retained.push(window);
    }
    retained
}

// TODO: Add docstring
/// Windows should be sorted by the output coordinate set
fn filter_blacklisted_post_merge(
    windows: Vec<Window>,
    blacklist_cursor: Option<&mut BlacklistCursor>,
    strategy: BlacklistStrategy,
    look_back: u64,
    coord_set: CoordinateSet,
) -> Vec<Window> {
    match blacklist_cursor {
        Some(cursor) if !cursor.intervals.is_empty() => {
            let intervals = cursor.intervals.as_slice();
            let mut retained: Vec<Window> = Vec::with_capacity(windows.len());
            for entry in windows.into_iter() {
                if entry.merged {
                    if is_blacklisted(
                        intervals,
                        strategy,
                        entry.start_for(coord_set) as u64,
                        entry.end_for(coord_set) as u64,
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
