// Streaming preparation pipeline for BED-like genomic windows.
//
// This module implements a memory-bounded, chromosome-streaming pipeline that:
// 1) Validates and loads a `near` interval set (non-overlapping, unique edges).
// 2) Loads and coalesces blacklist intervals (with halo).
// 3) Streams the primary input by chromosome in chunks, applying early filters,
//    nearest-distance binning, spacing, merging, and deduplication.
// 4) Writes per-chromosome temporary files and finally concatenates them,
//    enforcing `min_per_group` in a final pass.
//
// The implementation favors determinism, clear rules, and low memory usage. It
// assumes input is sorted by (chrom, start). If it is not, you should either
// sort upstream or add an explicit sort prepass.

use crate::commands::prepare_windows::chunk::{flush_chromosome, process_and_write_chunk};
use crate::commands::prepare_windows::config::*;
use crate::commands::prepare_windows::near_file::{load_near_index, nearest_edge_distance};
use crate::commands::prepare_windows::parsers::{
    parse_distance_bins, parse_record_line, parse_score_filter, resolve_column_indices,
};
use crate::commands::prepare_windows::resizers::apply_size_transform;
use crate::commands::prepare_windows::writers::{
    ChromTempWriter, concatenate_temps_enforcing_min_per_group, finalize_temp_writers,
};
use crate::shared::bed::{detect_header, line_looks_like_header};
use crate::shared::blacklist::{is_blacklisted, load_blacklists};
use crate::shared::reference::load_chrom_sizes;
use crate::shared::tiled_run::make_temp_dir;
use anyhow::{Context, Result, bail};
use fxhash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::{env, fs, mem};

/// Final window representation used throughout the pipeline.
///
/// Coordinates are 0-based, half-open `[start, end)`.
#[derive(Debug, Clone)]
pub struct FinalWindow {
    pub chrom: String,
    pub start: u32,
    pub end: u32,
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
struct BlacklistCursor {
    intervals: Vec<(u64, u64)>,
    cursor: usize,
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

    // Load near index (validated)
    let near_index = if let Some(path) = &cfg.near {
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
                    cursor: 0,
                },
            );
        }
    }
    let blacklist_look_back: u64 = 0;

    // Open input and initial reader
    let input_reader: Box<dyn BufRead> = if let Some(path) = &cfg.input {
        if path.as_os_str() == "-" {
            Box::new(BufReader::with_capacity(1 << 20, std::io::stdin()))
        } else {
            Box::new(BufReader::with_capacity(
                1 << 20,
                File::open(path).with_context(|| format!("Opening input {:?}", path))?,
            ))
        }
    } else {
        Box::new(BufReader::with_capacity(1 << 20, std::io::stdin()))
    };

    // Prepare per-run temporary directory and chromosome writers
    let base_output_dir = match &cfg.output {
        Some(path) if path.as_os_str() != "-" => path
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

    // Optional header handling for input: skip if present
    let mut reader = input_reader;
    let mut line_buffer = String::new();
    let mut pending_line: Option<String> = None;

    match cfg.header {
        HeaderMode::Present => {
            line_buffer.clear();
            reader.read_line(&mut line_buffer)?; // discard header
            line_buffer.clear();
        }
        HeaderMode::Auto => loop {
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

    // Pre-resolve score filter missing behavior
    let drop_missing_scores =
        matches!(cfg.score_missing, MissingScore::Drop) && cfg.score_filter.is_some();

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
                cfg,
            )?;

            // Move to new chromosome
            current_chrom.clear();
            current_chrom = chrom.clone();
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
        let chrom_size = chrom_sizes_map.get(&chrom).copied();
        let Some((final_start, final_end)) = apply_size_transform(start, end, chrom_size, cfg)
        else {
            continue;
        };

        // Blacklist pre-check on final coords
        if let Some(cursor) = blacklist_cursors.get_mut(&chrom) {
            if !cursor.intervals.is_empty()
                && is_blacklisted(
                    cursor.intervals.as_slice(),
                    cfg.blacklist_strategy,
                    final_start as u64,
                    final_end as u64,
                    blacklist_look_back,
                    &mut cursor.cursor,
                )
            {
                continue;
            }
        }

        // Nearest distance and binning
        let (composed_group, maybe_skip) = if let Some(near_idx) = &near_index {
            if let Some(near_chrom) = near_idx.per_chrom.get(&chrom) {
                if near_chrom.intervals.is_empty() {
                    (input_group.clone(), false)
                } else {
                    // Compute distance (0 if overlapping)
                    let signed = matches!(cfg.distance_sign, DistSign::Signed);
                    if let Some((mut dist, maybe_group_id)) = nearest_edge_distance(
                        final_start,
                        final_end,
                        near_chrom,
                        &cfg.near_edge,
                        signed,
                    ) {
                        // Direction filter
                        if matches!(cfg.near_direction, NearDirection::Upstream) && dist > 0 {
                            continue;
                        }
                        if matches!(cfg.near_direction, NearDirection::Downstream) && dist < 0 {
                            continue;
                        }

                        // Absolute distance if requested
                        if matches!(cfg.distance_sign, DistSign::Absolute) {
                            dist = dist.abs();
                        }

                        // Max distance
                        if let Some(maxd) = cfg.distance_max {
                            if dist.unsigned_abs() > maxd {
                                continue;
                            }
                        }

                        // Distance bin label
                        let bin_label = if let Some(bins) = &distance_bins {
                            bins.match_label(dist)
                        } else {
                            None
                        };

                        // Compose group label:
                        // Motivation for omitting near group when equidistant to two intervals:
                        // we validated near intervals to avoid overlaps and duplicate edges,
                        // but exact equidistance across a gap can still occur due to geometry.
                        // Attaching a specific near group in that rare case would imply a false
                        // association. We therefore keep only the input group + bin label.
                        let near_group_name = maybe_group_id
                            .map(|id| near_idx.group_id_to_name[id as usize].as_str());
                        let group_label = if near_group_name.is_some() {
                            compose_group_label(&input_group, near_group_name, bin_label)
                        } else {
                            compose_group_label(&input_group, None, bin_label)
                        };
                        (group_label, false)
                    } else {
                        (input_group.clone(), false)
                    }
                }
            } else {
                (input_group.clone(), false)
            }
        } else {
            // No near file; group may be from input or empty
            (input_group.clone(), false)
        };

        if maybe_skip {
            continue;
        }

        // Accumulate into batch
        current_batch.push(FinalWindow {
            chrom,
            start: final_start,
            end: final_end,
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
