// Streaming preparation pipeline for BED-like genomic windows.
//
// This module implements a memory-bounded, chromosome-streaming pipeline that:
// 1) Validates and loads a `near` interval set (non-overlapping, unique edges).
// 2) Loads and coalesces blacklist intervals (with halo).
// 3) Streams the primary input by chromosome in chunks, applying early filters,
//    nearest-distance binning (with `-/+/=` prefixes that reflect direction), spacing,
//    merging, and deduplication.
// 4) Writes per-chromosome temporary files and finally concatenates them,
//    enforcing `min_per_group` in a final pass.
//
// The implementation favors determinism, clear rules, and low memory usage. It
// assumes input is sorted by (chrom, start). If it is not, you should either
// sort upstream or add an explicit sort prepass.

use crate::commands::prepare_windows::chunk::{flush_chromosome, process_and_write_chunk};
use crate::commands::prepare_windows::config::*;
use crate::commands::prepare_windows::near_file::{
    NearHit, NearIndex, NearSide, NearTie, NearestResult, load_near_index, nearest_edge_distance,
};
use crate::commands::prepare_windows::parsers::{
    DistanceBins, parse_distance_bins, parse_record_line, parse_score_filter,
    resolve_column_indices,
};
use crate::commands::prepare_windows::resizers::apply_size_transform;
use crate::commands::prepare_windows::writers::{
    ChromTempWriter, concatenate_temps_enforcing_min_per_group, finalize_temp_writers,
};
use crate::shared::bed::{detect_header, line_looks_like_header};
use crate::shared::blacklist::{is_blacklisted, load_blacklists};
use crate::shared::io::open_text_reader;
use crate::shared::reference::load_chrom_sizes;
use crate::shared::tiled_run::make_temp_dir;
use anyhow::{Context, Result, bail};
use fxhash::{FxHashMap, FxHashSet};
use std::collections::hash_map::Entry;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;
use std::{env, fs, mem};

/// Final window representation used throughout the pipeline.
///
/// Coordinates are 0-based, half-open `[start, end)`.
#[derive(Debug, Clone)]
pub struct FinalWindow {
    pub chrom: Arc<str>,
    pub start: u32,
    pub end: u32,
    pub merged: bool,
    pub group: String,      // Empty means "no group label"
    pub score: Option<f32>, // Present only if parsed and requested
}

impl FinalWindow {
    #[inline]
    pub fn length(&self) -> u32 {
        self.end - self.start
    }
}

/// Streaming cursor for blacklist interval sweeps.
#[derive(Debug, Default)]
pub struct BlacklistCursor {
    pub intervals: Vec<(u64, u64)>,
    pub pre_cursor: usize,  // Early filtering
    pub post_cursor: usize, // Post-merge filtering
}

/// Run the prepare pipeline using the provided configuration.
///
/// This function orchestrates near and blacklist loading, streams the main input,
/// performs early filtering and annotation, enforces spacing and merging with
/// chunked writes, and finally concatenates chromosome-temporaries applying
/// the `min_per_group` filter.
///
/// Parameters
/// ----------
/// - cfg:
///     Command-line configuration.
///
/// Returns
/// -------
/// - ok:
///     Success or error.
pub fn run(cfg: &PrepareConfig) -> Result<()> {
    // Compile distance bins (if any)
    let distance_bins = if let Some(specs) = &cfg.distance_bins {
        Some(parse_distance_bins(specs)?)
    } else {
        None
    };

    // Compile score filter (if any)
    let score_filter = if let Some(expr) = &cfg.score_filter {
        Some(parse_score_filter(expr)?)
    } else {
        None
    };

    // How to handle missing scores
    let drop_missing_scores =
        matches!(cfg.score_missing, MissingScore::Drop) && cfg.score_filter.is_some();

    // Load near index (validated)
    let mut near_index = if let Some(path) = &cfg.near {
        let has_header_final = match cfg.near_header {
            HeaderMode::Present => true,
            HeaderMode::Absent => false,
            HeaderMode::Auto => detect_header(path, cfg.sep).unwrap_or(false),
        };
        let group_col_present = true; // If your near has optional group, you can make this configurable
        Some(load_near_index(
            path,
            cfg.sep,
            has_header_final,
            group_col_present,
        )?)
    } else {
        None
    };

    // Load blacklist intervals (optional)
    let mut blacklist_cursors: FxHashMap<String, BlacklistCursor> = FxHashMap::default();
    if let Some(paths) = &cfg.blacklist {
        let loaded = load_blacklists(paths.as_slice(), 1, cfg.blacklist_halo as u64, None)?;
        for (chrom, intervals) in loaded.into_iter() {
            blacklist_cursors.insert(
                chrom,
                BlacklistCursor {
                    intervals,
                    pre_cursor: 0,
                    post_cursor: 0,
                },
            );
        }
    }
    let blacklist_look_back: u64 = 0;

    // Open input and initial reader
    let input_reader: Box<dyn BufRead> = if cfg.input.as_os_str() == "-" {
        Box::new(BufReader::with_capacity(1 << 20, std::io::stdin()))
    } else {
        open_text_reader(&cfg.input)?
    };

    // TODO: Add info to config about this temporary directory, so user knows to not call (with "-") in a folder with no storage or write permissions!
    // Prepare per-run temporary directory and chromosome writers
    let base_output_dir = match &cfg.output {
        path if path.as_os_str() != "-" => path
            .parent()
            .map(PathBuf::from)
            .unwrap_or_else(|| env::current_dir().expect("current_dir")),
        _ => env::temp_dir(),
    };
    let temp_dir =
        make_temp_dir(&base_output_dir, "prepare_windows").context("create per-run temp dir")?;
    let mut temp_writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();

    // Group counts across all chromosomes (after spacing/merge)
    let mut global_group_counts: FxHashMap<String, u64> = FxHashMap::default();

    // State for streaming by chromosome with chunking
    let chunk_size: usize = 5_000_000; // You can expose as a flag if desired
    let mut current_chrom: String = String::new();
    let mut processed_chromosomes: FxHashSet<String> = FxHashSet::default();
    let mut last_start_for_current: Option<u32> = None;
    let mut carryover_tail: Vec<FinalWindow> = Vec::new();
    let mut current_batch: Vec<FinalWindow> = Vec::with_capacity(chunk_size + 1024);
    let mut current_chrom_size: Option<u32> = None;
    let mut chrom_intern: FxHashMap<String, Arc<str>> = FxHashMap::default();

    // Optional header handling for input: skip if present
    let mut reader = input_reader;
    let mut line_buffer = String::new();
    let mut pending_line: Option<String> = None;

    match cfg.header {
        HeaderMode::Present => {
            line_buffer.clear();
            reader.read_line(&mut line_buffer)?; // Discard header
            line_buffer.clear();
        }
        HeaderMode::Auto => loop {
            // TODO: Use detect_header() to the degree possible
            line_buffer.clear();
            let bytes = reader.read_line(&mut line_buffer)?;
            if bytes == 0 {
                break;
            }
            let candidate = line_buffer.trim_end_matches(['\n', '\r']);
            if candidate.is_empty() {
                continue;
            }
            if line_looks_like_header(candidate, cfg.sep) {
                line_buffer.clear();
                break;
            } else {
                pending_line = Some(candidate.to_string());
                line_buffer.clear();
                break;
            }
        },
        HeaderMode::Absent => {}
    }

    // Map from chromosome to size if trimming/dropping is enabled
    let chrom_sizes_map: FxHashMap<String, u32> =
        if matches!(cfg.oob, OobPolicy::Trim | OobPolicy::Drop) && cfg.chrom_sizes.is_some() {
            load_chrom_sizes(cfg.chrom_sizes.as_ref().unwrap())?
        } else {
            FxHashMap::default()
        };

    let column_indices =
        resolve_column_indices(&cfg.cols, &cfg.group_cols, cfg.score_col.as_deref())?;

    // Stream input records
    loop {
        // TODO: Explain this if-else statement
        if let Some(mut pending) = pending_line.take() {
            mem::swap(&mut line_buffer, &mut pending);
        } else {
            line_buffer.clear();
            let bytes = reader.read_line(&mut line_buffer)?;
            if bytes == 0 {
                break;
            }
        }

        let line = line_buffer.trim_end_matches(['\n', '\r']);
        if line.is_empty() {
            continue;
        }
        if line.trim_start().starts_with('#') {
            continue;
        }

        let (chrom, start, end, input_group, score) =
            parse_record_line(line, cfg.sep, &column_indices)?;

        if current_chrom.is_empty() {
            current_chrom = chrom.clone();
            current_chrom_size = chrom_sizes_map.get(&current_chrom).copied();
            last_start_for_current = None;
        }

        if chrom != current_chrom {
            if !current_chrom.is_empty() {
                processed_chromosomes.insert(current_chrom.clone());
            }
            if processed_chromosomes.contains(&chrom) {
                bail!(
                    "Input is not sorted: chromosome '{}' appears after it was already processed",
                    chrom
                );
            }

            // Flush remaining chunk for previous chromosome
            flush_chromosome(
                &current_chrom,
                &mut carryover_tail,
                &mut current_batch,
                &mut temp_writers,
                &temp_dir,
                &mut global_group_counts,
                blacklist_cursors.get_mut(&current_chrom),
                blacklist_look_back,
                cfg,
            )?;

            // Move to new chromosome
            current_chrom.clear();
            current_chrom = chrom.clone();
            current_chrom_size = chrom_sizes_map.get(&current_chrom).copied();
            last_start_for_current = None;
        }

        if let Some(last) = last_start_for_current {
            if chrom == current_chrom && start < last {
                bail!(
                    "Input is not sorted: chromosome '{}' has start {} before previous {}",
                    chrom,
                    start,
                    last
                );
            }
        }
        last_start_for_current = Some(start);

        // Early score filtering
        if let Some(sf) = &score_filter {
            match score {
                Some(sv) => {
                    if !sf.eval(sv) {
                        continue;
                    }
                }
                None => {
                    if drop_missing_scores {
                        continue;
                    }
                }
            }
        }

        // Transform to final coordinates
        let chrom_size = current_chrom_size;
        let Some((final_start, final_end)) = apply_size_transform(start, end, chrom_size, cfg)
        else {
            continue;
        };

        // Blacklist pre-check on pre-merge full-size coordinates
        if let Some(cursor) = blacklist_cursors.get_mut(&chrom) {
            if !cursor.intervals.is_empty()
                && is_blacklisted(
                    cursor.intervals.as_slice(),
                    cfg.blacklist_strategy,
                    final_start as u64,
                    final_end as u64,
                    blacklist_look_back,
                    &mut cursor.pre_cursor,
                )
            {
                continue;
            }
        }

        // Find nearest intervals and bin by distance and update the group name
        let composed_group = if let Some(composed_group) = add_near_group_annotation(
            &mut near_index,
            &chrom,
            final_start,
            final_end,
            cfg,
            distance_bins.as_ref(),
            &input_group,
        ) {
            composed_group
        } else {
            continue;
        };

        // Intern chromosome identifiers so every window for the same chromosome
        // shares a single heap allocation; the Arc keeps those copies alive
        // while allowing cheap cloning across batches.
        let chrom_arc = match chrom_intern.entry(chrom.clone()) {
            Entry::Occupied(entry) => entry.get().clone(),
            Entry::Vacant(entry) => {
                let arc: Arc<str> = chrom.into_boxed_str().into();
                entry.insert(arc.clone());
                arc
            }
        };

        // Accumulate into batch
        current_batch.push(FinalWindow {
            chrom: chrom_arc,
            start: final_start,
            end: final_end,
            merged: false,
            group: composed_group,
            score,
        });

        // If batch exceeds chunk size, process and write safe prefix
        if current_batch.len() >= chunk_size {
            process_and_write_chunk(
                &current_chrom,
                &mut carryover_tail,
                &mut current_batch,
                &mut temp_writers,
                &temp_dir,
                &mut global_group_counts,
                blacklist_cursors.get_mut(&current_chrom),
                blacklist_look_back,
                cfg,
            )?;
        }
    }

    // Flush final chromosome
    if !current_chrom.is_empty() {
        flush_chromosome(
            &current_chrom,
            &mut carryover_tail,
            &mut current_batch,
            &mut temp_writers,
            &temp_dir,
            &mut global_group_counts,
            blacklist_cursors.get_mut(&current_chrom),
            blacklist_look_back,
            cfg,
        )?;
    }

    // Final pass: concatenate temps, enforcing min_per_group
    let temp_entries = finalize_temp_writers(&mut temp_writers)?;
    concatenate_temps_enforcing_min_per_group(cfg, &temp_entries, &global_group_counts)?;
    fs::remove_dir_all(&temp_dir).ok();

    Ok(())
}

/// Compose the output group label from optional parts, omitting empty segments.
///
/// The order is `{input_group?}.{near_group?}.{bin_label?}`.
///
/// Parameters
/// ----------
/// - input_group:
///     Group string from `--group-cols` concatenation (may be empty).
/// - near_group:
///     Optional near group name.
/// - bin_label:
///     Optional distance bin label.
///
/// Returns
/// -------
/// - group:
///     Composed group label (possibly empty string).
fn compose_group_label(
    input_group: &str,
    near_group: Option<&str>,
    bin_label: Option<&str>,
) -> String {
    let mut parts: Vec<&str> = Vec::with_capacity(3);
    if !input_group.is_empty() {
        parts.push(input_group);
    }
    if let Some(ng) = near_group {
        if !ng.is_empty() {
            parts.push(ng);
        }
    }
    if let Some(label) = bin_label {
        if !label.is_empty() {
            parts.push(label);
        }
    }
    parts.join(".")
}

/// Build the near-annotation label for a window using the cursor-based nearest lookup.
/// Returns the composed group label (input group + optional near group + optional bin),
/// or `None` when the window should be dropped by policy (e.g., ties=drop, out of range).
fn add_near_group_annotation(
    near_index: &mut Option<NearIndex>,
    chrom: &str,
    final_start: u32,
    final_end: u32,
    cfg: &PrepareConfig,
    distance_bins: Option<&DistanceBins>,
    input_group: &str,
) -> Option<String> {
    // If no near index (or no data for this chrom), keep the input group as-is
    let Some(near_idx) = near_index.as_mut() else {
        return Some(input_group.to_owned());
    };
    let Some(near_chrom) = near_idx.per_chrom.get_mut(chrom) else {
        return Some(input_group.to_owned());
    };
    if near_chrom.intervals.is_empty() {
        return Some(input_group.to_owned());
    }

    // Helpers
    let is_signed_mode = matches!(cfg.distance_sign, DistSign::Signed);

    let within_max_distance = |distance_bp: i32| -> bool {
        if let Some(max_abs) = cfg.distance_max {
            return distance_bp.unsigned_abs() <= max_abs;
        }
        true
    };

    let normalize_for_binning = |distance_bp: i32| -> i32 {
        if matches!(cfg.distance_sign, DistSign::Absolute) {
            distance_bp.abs()
        } else {
            distance_bp
        }
    };

    // Compute nearest using the in-struct cursor
    let Some(nearest_result) = nearest_edge_distance(
        final_start,
        final_end,
        near_chrom,
        &cfg.near_edge,
        &cfg.near_direction,
        is_signed_mode,
    ) else {
        // No nearby intervals within the considered directions -> keep input group only
        return Some(input_group.to_owned());
    };

    // Compose a side-prefixed label for a single hit.
    let make_side_label = |hit: &NearHit| -> String {
        let side_prefix = match hit.side {
            NearSide::Upstream => "-",
            NearSide::Downstream => "+",
            NearSide::Overlap => "=",
        };
        match hit.group_id {
            Some(id) => {
                let name = near_idx.group_id_to_name[id as usize].as_str();
                format!("{side_prefix}{name}")
            }
            None => side_prefix.to_string(),
        }
    };

    match nearest_result {
        NearestResult::Single(mut hit) => {
            // Distance filtering
            if !within_max_distance(hit.distance) {
                return None;
            }
            // Normalize for binning/output according to cfg.distance_sign
            hit.distance = normalize_for_binning(hit.distance);

            // Optional bin label
            let bin_label = distance_bins.and_then(|bins| bins.match_label(hit.distance));

            // Label for the near group (+-=)
            let near_label = make_side_label(&hit);

            // Compose final group label
            Some(compose_group_label(
                input_group,
                Some(near_label.as_str()),
                bin_label,
            ))
        }

        NearestResult::Tie(NearTie {
            mut upstream,
            mut downstream,
        }) => {
            // Caller policy: drop on ties if requested
            if matches!(cfg.near_ties, NearTiePolicy::Drop) {
                return None;
            }

            // Apply distance threshold and normalization to both (when present)
            let mut kept_labels: Vec<String> = Vec::with_capacity(2);

            if let Some(ref mut upstream_hit) = upstream {
                if within_max_distance(upstream_hit.distance) {
                    upstream_hit.distance = normalize_for_binning(upstream_hit.distance);
                    kept_labels.push(make_side_label(upstream_hit));
                }
            }
            if let Some(ref mut downstream_hit) = downstream {
                if within_max_distance(downstream_hit.distance) {
                    downstream_hit.distance = normalize_for_binning(downstream_hit.distance);
                    kept_labels.push(make_side_label(downstream_hit));
                }
            }

            if kept_labels.is_empty() {
                return None;
            }

            // Distances in a tie have the same absolute value by construction; pick one for binning
            let bin_distance = upstream
                .as_ref()
                .map(|h| h.distance)
                .or_else(|| downstream.as_ref().map(|h| h.distance))
                .unwrap();

            let bin_label = distance_bins.and_then(|bins| bins.match_label(bin_distance));

            let near_combo = kept_labels.join("/");

            Some(compose_group_label(
                input_group,
                Some(near_combo.as_str()),
                bin_label,
            ))
        }
    }
}
