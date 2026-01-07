//! Streaming preparation pipeline for BED-like genomic windows.
//!
//! The intended labeling and filtering logic is specified in the `label_and_filter_logic.md` document.
//!
//!  This module implements a memory-bounded, chromosome-streaming pipeline that:
//!  - Validates and loads a `near` interval set (half-open, duplicate edges handled by `--near-duplicates`).
//!  - Loads and combines blacklist intervals (with an optional halo).
//!  - Streams the BAM file by chromosome in chunks, applying early filters,
//!    nearest-distance binning (with `-/+/=` prefixes that reflect direction),
//!    minimum-distance filtering, merging, and deduplication.
//!  - Writes per-chromosome temporary files and concatenates them in a final pass.
//!
//!  The implementation favors determinism, clear rules, and low memory usage. It
//!  assumes input is sorted by (chrom, start).

use crate::commands::prepare_windows::chunk::{flush_chromosome, process_and_write_chunk};
use crate::commands::prepare_windows::config::*;
use crate::commands::prepare_windows::filters::{
    ExcludeRule, filter_and_write_output, parse_exclude_rules, parse_min_per_rules,
};
use crate::commands::prepare_windows::labels::{LabelSchema, LabelTuple};
use crate::commands::prepare_windows::near_file::load_near_index;
use crate::commands::prepare_windows::parsers::{
    parse_distance_bins, parse_record_line, parse_score_filter, resolve_column_indices,
};
use crate::commands::prepare_windows::resizers::apply_size_transform;
use crate::commands::prepare_windows::writers::{ChromTempWriter, finalize_temp_writers};
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
use std::time::Instant;
use std::{env, fs, mem};

/// Final window representation used throughout the pipeline.
///
/// Coordinates are 0-based, half-open `[start, end)`.
#[derive(Debug, Clone)]
pub struct FinalWindow {
    pub chrom: Arc<str>,
    pub original_start: u32,
    pub original_end: u32,
    pub resized_start: u32,
    pub resized_end: u32,
    pub merged: bool,
    pub label_tuples: Vec<LabelTuple>, // Atomic label tuples for this window
    pub group_key: String,             // Grouping key for merge and minimum-distance filtering
    pub score: Option<f32>,            // Present only if parsed and requested
}

impl FinalWindow {
    #[inline]
    pub fn start_for(&self, coord_set: CoordinateSet) -> u32 {
        match coord_set {
            CoordinateSet::Original => self.original_start,
            CoordinateSet::Resized => self.resized_start,
        }
    }

    #[inline]
    pub fn end_for(&self, coord_set: CoordinateSet) -> u32 {
        match coord_set {
            CoordinateSet::Original => self.original_end,
            CoordinateSet::Resized => self.resized_end,
        }
    }

    #[inline]
    pub fn length_for(&self, coord_set: CoordinateSet) -> u32 {
        self.end_for(coord_set) - self.start_for(coord_set)
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
pub fn run(cfg: &PrepareConfig) -> Result<()> {
    let start_time = Instant::now();

    println!("Preparing BED-like file");

    // TODO: Print pipeline that will be applied

    // TODO: Validate IO paths early (and other args)

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

    // Resolve label schema and key references
    let label_schema = LabelSchema::new(&cfg.compose)?;
    let merge_key = label_schema.resolve_key(&cfg.merge_key)?;
    let out_labels = label_schema.resolve_keys(&cfg.out_labels)?;

    let min_per_rules = parse_min_per_rules(&cfg.min_per, &label_schema)?;
    let exclude_rules: Vec<ExcludeRule> = parse_exclude_rules(&cfg.exclude_labels, &label_schema)?;

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

        let strand_col_present = true; // TODO: configure as optional
        let group_col_present = true; // TODO: configure as optional
        Some(load_near_index(
            path,
            cfg.sep,
            has_header_final,
            strand_col_present,
            group_col_present,
            matches!(cfg.near_edge, NearEdge::Upstream | NearEdge::Downstream),
            cfg.near_duplicates,
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

    // State for streaming by chromosome with chunking
    // Chunk size in windows to limit memory while keeping sequential IO fast
    let chunk_size: usize = 5_000_000;
    // Current chromosome name for grouping and change detection
    let mut current_chrom: String = String::new();
    // Tracks chromosomes already seen to guard against out-of-order input
    let mut processed_chromosomes: FxHashSet<String> = FxHashSet::default();
    // Previous window start for the current chromosome for ordering checks
    let mut last_start_for_current: Option<u32> = None;
    // Tail windows that must carry into the next chunk for merging
    let mut carryover_tail: Vec<FinalWindow> = Vec::new();
    // Batch of windows collected for processing before flush
    let mut current_batch: Vec<FinalWindow> = Vec::with_capacity(chunk_size + 1024);
    // Known size for the current chromosome for trimming and bounds checks
    let mut current_chrom_size: Option<u32> = None;
    // Interned chromosome strings to avoid repeated allocations
    let mut chromosome_intern_pool: FxHashMap<String, Arc<str>> = FxHashMap::default();

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
                blacklist_cursors.get_mut(&current_chrom),
                blacklist_look_back,
                current_chrom_size,
                cfg,
                &mut near_index,
                distance_bins.as_ref(),
                &label_schema,
                &merge_key,
                &out_labels,
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

        // Resized coordinates are computed per record even when merging uses originals
        // Transform to resized coordinates
        let chrom_size = current_chrom_size;
        let Some((resized_start, resized_end)) = apply_size_transform(start, end, chrom_size, cfg)
        else {
            continue;
        };

        // TODO: Get the cursor along with the chromosome sizes to avoid 1B get_mut calls (hashing)
        // Blacklist pre-check on resized coordinates
        if let Some(cursor) = blacklist_cursors.get_mut(&chrom) {
            if !cursor.intervals.is_empty()
                && is_blacklisted(
                    cursor.intervals.as_slice(),
                    cfg.blacklist_strategy,
                    resized_start as u64,
                    resized_end as u64,
                    blacklist_look_back,
                    &mut cursor.pre_cursor,
                )
            {
                continue;
            }
        }

        let label_tuples = vec![LabelTuple::new(input_group.clone())];

        // Intern chromosome identifiers so every window for the same chromosome
        // shares a single heap allocation. The Arc keeps those copies alive
        // while allowing cheap cloning across batches.
        let chrom_arc = match chromosome_intern_pool.entry(chrom.clone()) {
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
            original_start: start,
            original_end: end,
            resized_start,
            resized_end,
            merged: false,
            label_tuples,
            group_key: input_group,
            score,
        });

        // If batch exceeds chunk size, process and write the processed region
        if current_batch.len() >= chunk_size {
            process_and_write_chunk(
                &current_chrom,
                &mut carryover_tail,
                &mut current_batch,
                &mut temp_writers,
                &temp_dir,
                blacklist_cursors.get_mut(&current_chrom),
                blacklist_look_back,
                current_chrom_size,
                cfg,
                &mut near_index,
                distance_bins.as_ref(),
                &label_schema,
                &merge_key,
                &out_labels,
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
            blacklist_cursors.get_mut(&current_chrom),
            blacklist_look_back,
            current_chrom_size,
            cfg,
            &mut near_index,
            distance_bins.as_ref(),
            &label_schema,
            &merge_key,
            &out_labels,
        )?;
    }

    // Final pass: apply filtering and write output
    let temp_entries = finalize_temp_writers(&mut temp_writers)?;
    filter_and_write_output(
        cfg,
        &temp_entries,
        &label_schema,
        &out_labels,
        &min_per_rules,
        &exclude_rules,
        &base_output_dir,
    )?;
    if temp_dir.exists() {
        fs::remove_dir_all(&temp_dir).context("Removing temp dir")?;
    }

    let elapsed = start_time.elapsed();
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}
