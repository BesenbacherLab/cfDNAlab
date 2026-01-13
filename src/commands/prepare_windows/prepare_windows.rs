//! Streaming preparation pipeline for BED-like genomic windows.
//!
//! The intended labeling and filtering logic is specified in the `label_and_filter_logic.md` document.
//!
//!  This module implements a memory-bounded, chromosome-streaming pipeline that:
//!  - Validates and loads a `near` interval set (half-open, duplicate edges handled by `--near-duplicates`).
//!  - Loads and combines blacklist intervals (with an optional halo).
//!  - Streams the BED file by chromosome in chunks, applying early filters,
//!    resize or flank transforms, blacklist checks, deduplication, merging,
//!    clustering, minimum-distance filtering, and near annotation with distance bins.
//!  - Writes per-chromosome temporary files and concatenates them in a final pass.
//!
//!  The implementation favors determinism, clear rules, and low memory usage. It
//!  assumes input is sorted by (chrom, start).

use crate::commands::prepare_windows::chunk::{flush_chromosome, process_and_write_chunk};
use crate::commands::prepare_windows::config::*;
use crate::commands::prepare_windows::filters::{
    ExcludeRule, filter_and_write_output, parse_exclude_rules, parse_min_per_rules,
    validate_available_keys, validate_compositions_available,
};
use crate::commands::prepare_windows::labels::{AtomicLabelPart, LabelSchema, LabelTuple};
use crate::commands::prepare_windows::near_file::load_near_index;
use crate::commands::prepare_windows::parsers::{
    parse_distance_bins, parse_record_line, parse_score_filter, parse_single_index,
    resolve_column_indices,
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
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::hash_map::Entry;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{Duration, Instant};
use std::{env, fs, mem};

/// Window representation used throughout the prepare pipeline.
///
/// Coordinates are 0-based, half-open `[start, end)`.
#[derive(Debug, Clone)]
pub struct Window {
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

impl Window {
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

// Removes the temp directory on drop
struct TempDirGuard {
    path: PathBuf,
    removed: bool,
}

impl TempDirGuard {
    fn new(base_dir: &Path, prefix: &str) -> Result<Self> {
        let path = make_temp_dir(base_dir, prefix).context("create per-run temp dir")?;
        Ok(Self {
            path,
            removed: false,
        })
    }

    fn path(&self) -> &Path {
        &self.path
    }

    fn remove(&mut self) -> Result<()> {
        if self.removed {
            return Ok(());
        }
        if self.path.exists() {
            fs::remove_dir_all(&self.path).context("Removing temp dir")?;
        }
        self.removed = true;
        Ok(())
    }
}

impl Drop for TempDirGuard {
    fn drop(&mut self) {
        let _ = self.remove();
    }
}

/// Run the prepare pipeline using the provided configuration.
pub fn run(cfg: &PrepareConfig) -> Result<()> {
    let start_time = Instant::now();

    println!("Preparing windows");

    // Prepare the temp directory early so we fail fast on missing write permissions
    println!("Start: Creating temporary directory");
    let base_output_dir = match &cfg.output {
        path if path.as_os_str() != "-" => path
            .parent()
            .map(PathBuf::from)
            .unwrap_or_else(|| env::current_dir().expect("current_dir")),
        _ => env::temp_dir(),
    };
    let mut temp_dir_guard = TempDirGuard::new(&base_output_dir, "prepare_windows")?;
    // Keep per-chrom temp files in a subdirectory so they can be removed once filtered
    let stream_temp_dir = make_temp_dir(temp_dir_guard.path(), "prepare_windows_stream")
        .context("create stream temp dir")?;

    // TODO: Validate IO paths early (and other args)

    if cfg.distance_bins.is_some() && cfg.near.is_none() {
        bail!("--distance-bins requires --near");
    }

    let use_input_chrom_order = cfg.chromosomes.chromosomes_file.is_none()
        && matches!(
            cfg.chromosomes.chromosomes.as_deref(),
            Some([single]) if single.eq_ignore_ascii_case("all")
        );
    let chromosomes = if use_input_chrom_order {
        Vec::new()
    } else {
        cfg.chromosomes.resolve_chromosomes(None)?
    };
    let allowed_chromosomes = if use_input_chrom_order {
        None
    } else {
        Some(chromosomes.iter().cloned().collect::<FxHashSet<String>>())
    };
    let chromosomes_for_blacklist = if use_input_chrom_order {
        None
    } else {
        Some(chromosomes.as_slice())
    };

    // Compile distance bins (if any)
    let distance_bins = if let Some(specs) = &cfg.distance_bins {
        println!("Start: Parsing distance bins");
        Some(parse_distance_bins(specs)?)
    } else {
        None
    };

    // Compile score filter (if any)
    let score_filter = if let Some(expr) = &cfg.score_filter {
        println!("Start: Parsing score filter");
        Some(parse_score_filter(expr)?)
    } else {
        None
    };

    // Resolve label schema and key references
    println!("Start: Resolving label schema");
    let label_schema = LabelSchema::new(&cfg.compose)?;
    let available_parts = available_atomic_parts(cfg);
    validate_compositions_available(&label_schema, &available_parts)?;
    let merge_key = label_schema.resolve_key(&cfg.merge_key)?;
    validate_available_keys(
        std::slice::from_ref(&merge_key),
        &label_schema,
        &available_parts,
        "merge-key",
    )?;
    let out_labels = label_schema.resolve_keys(&cfg.out_labels)?;
    validate_available_keys(&out_labels, &label_schema, &available_parts, "out-labels")?;

    if !cfg.min_per.is_empty() {
        println!("Start: Parsing min-per rules");
    }
    let min_per_rules = parse_min_per_rules(&cfg.min_per, &label_schema, &available_parts)?;

    if !cfg.exclude_labels.is_empty() {
        println!("Start: Parsing exclude rules");
    }
    let exclude_rules: Vec<ExcludeRule> =
        parse_exclude_rules(&cfg.exclude_labels, &label_schema, &available_parts)?;

    // How to handle missing scores
    let drop_missing_scores =
        matches!(cfg.score_missing, MissingScore::Drop) && cfg.score_filter.is_some();

    // Load near index (validated)
    let near_strand_col = cfg
        .near_strand_col
        .as_deref()
        .map(parse_single_index)
        .transpose()?;
    let mut near_group_cols: Vec<usize> = Vec::new();
    for col in &cfg.near_group_cols {
        near_group_cols.push(parse_single_index(col)?);
    }
    let near_group_cols = if near_group_cols.is_empty() {
        None
    } else {
        Some(near_group_cols)
    };

    let mut near_index = if let Some(path) = &cfg.near {
        println!("Start: Loading near file");
        let has_header_final = match cfg.near_header {
            HeaderMode::Present => true,
            HeaderMode::Absent => false,
            HeaderMode::Auto => detect_header(path, cfg.sep).unwrap_or(false),
        };

        let consider_strand_in_dups =
            matches!(cfg.near_edge, NearEdge::Upstream | NearEdge::Downstream)
                && near_strand_col.is_some();
        Some(load_near_index(
            path,
            cfg.sep,
            has_header_final,
            near_strand_col,
            near_group_cols.as_deref(),
            consider_strand_in_dups,
            cfg.near_duplicates,
        )?)
    } else {
        None
    };

    // Load blacklist intervals (optional)
    let mut blacklist_cursors: FxHashMap<String, BlacklistCursor> = FxHashMap::default();
    if let Some(paths) = &cfg.blacklist {
        println!("Start: Loading blacklist intervals");
        let loaded = load_blacklists(
            paths.as_slice(),
            1,
            cfg.blacklist_halo as u64,
            chromosomes_for_blacklist,
        )?;
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
    println!("Start: Opening input reader");
    let input_reader: Box<dyn BufRead> = if cfg.input.as_os_str() == "-" {
        Box::new(BufReader::with_capacity(1 << 20, std::io::stdin()))
    } else {
        open_text_reader(&cfg.input)?
    };

    // TODO: Add info to config about this temporary directory, so user knows to not call (with "-") in a folder with no storage or write permissions!
    let mut temp_writers: FxHashMap<String, ChromTempWriter> = FxHashMap::default();

    // State for streaming by chromosome with chunking
    // Chunk size in windows to limit memory while keeping sequential IO fast
    let chunk_size: usize = 5_000_000;
    // Current chromosome name for grouping and change detection
    let mut current_chrom: String = String::new();
    // Chromosome order as it appears in the input stream
    let mut chrom_order: Vec<String> = Vec::new();
    // Tracks chromosomes already seen to guard against out-of-order input
    let mut processed_chromosomes: FxHashSet<String> = FxHashSet::default();
    // Previous window start for the current chromosome for ordering checks
    let mut prev_start_for_current: Option<u32> = None;
    // Tail windows that must carry into the next chunk for merging
    let mut carryover_tail: Vec<Window> = Vec::new();
    // Batch of windows collected for processing before flush
    let mut current_batch: Vec<Window> = Vec::with_capacity(chunk_size + 1024);
    // Known size for the current chromosome for trimming and bounds checks
    let mut current_chrom_size: Option<u32> = None;
    // Interned chromosome strings to avoid repeated allocations
    let mut chromosome_intern_pool: FxHashMap<String, Arc<str>> = FxHashMap::default();

    // Optional header handling for input: skip if present
    let mut reader = input_reader;
    let mut line_buffer = String::new();
    let mut pending_line: Option<String> = None;

    let has_known_chroms = !chromosomes.is_empty();
    let pb = if has_known_chroms {
        let bar = Arc::new(ProgressBar::new(chromosomes.len() as u64));
        bar.set_style(
            ProgressStyle::default_bar()
                .template("Chromosomes {bar:40} {pos}/{len} [{elapsed_precise}] {msg}")
                .unwrap(),
        );
        bar.set_position(0);
        bar
    } else {
        let spinner = Arc::new(ProgressBar::new_spinner());
        spinner.set_style(
            ProgressStyle::default_spinner()
                .template("Chromosomes {spinner} {msg} [{elapsed_precise}]")
                .unwrap(),
        );
        spinner.enable_steady_tick(Duration::from_millis(100));
        spinner.set_message("0 processed");
        spinner
    };

    let mut processed_chrom_count: u64 = 0;

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
    let has_size_transform = cfg.resize.is_some() || cfg.flank.is_some();
    let chrom_sizes_map: FxHashMap<String, u32> = if has_size_transform
        && matches!(cfg.oob, OobPolicy::Trim | OobPolicy::Drop)
        && cfg.chrom_sizes.is_some()
    {
        println!("Start: Loading chromosome sizes");
        load_chrom_sizes(cfg.chrom_sizes.as_ref().unwrap())?
    } else {
        FxHashMap::default()
    };

    println!("Start: Resolving input columns");
    let column_indices =
        resolve_column_indices(&cfg.cols, &cfg.group_cols, cfg.score_col.as_deref())?;

    // Stream input records
    println!("Start: Streaming input records");
    loop {
        // Header auto-detection may have already read the first data line, so consume it here
        // Swap to avoid losing the buffered line and to reuse the existing allocation
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

        if let Some(allowed) = allowed_chromosomes.as_ref() {
            if !allowed.contains(&chrom) {
                continue;
            }
        }

        if current_chrom.is_empty() {
            current_chrom = chrom.clone();
            current_chrom_size = chrom_sizes_map.get(&current_chrom).copied();
            prev_start_for_current = None;
            if use_input_chrom_order {
                chrom_order.push(current_chrom.clone());
            }
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
                &stream_temp_dir,
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
            prev_start_for_current = None;
            if use_input_chrom_order {
                chrom_order.push(current_chrom.clone());
            }

            processed_chrom_count += 1;
            if has_known_chroms {
                pb.inc(1);
                pb.set_message(format!("Last {}", current_chrom));
            } else {
                pb.set_message(format!(
                    "{} processed (last {})",
                    processed_chrom_count, current_chrom
                ));
            }
        }

        if let Some(prev_start) = prev_start_for_current {
            if chrom == current_chrom && start < prev_start {
                bail!(
                    "Input is not sorted: chromosome '{}' has start {} before previous {}",
                    chrom,
                    start,
                    prev_start
                );
            }
        }
        prev_start_for_current = Some(start);

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
        current_batch.push(Window {
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
                &stream_temp_dir,
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
            &stream_temp_dir,
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
        processed_chrom_count += 1;
        if has_known_chroms {
            pb.inc(1);
            pb.set_message(format!("Last {}", current_chrom));
        } else {
            pb.set_message(format!(
                "{} processed (last {})",
                processed_chrom_count, current_chrom
            ));
        }
    }

    if has_known_chroms {
        pb.finish_with_message(format!(
            "{} processed of {}",
            processed_chrom_count,
            chromosomes.len()
        ));
    } else {
        pb.finish_with_message(format!(
            "{} processed (input order)",
            processed_chrom_count
        ));
    }

    // Final pass: apply filtering and write output
    println!("Start: Finalizing output");
    let temp_entries = finalize_temp_writers(&mut temp_writers)?;
    let output_chromosomes = if use_input_chrom_order {
        chrom_order
    } else {
        chromosomes
    };

    filter_and_write_output(
        cfg,
        &temp_entries,
        &label_schema,
        &out_labels,
        &min_per_rules,
        &exclude_rules,
        temp_dir_guard.path(),
        &output_chromosomes,
    )?;
    println!("Start: Removing temporary directory");
    temp_dir_guard.remove()?;

    let elapsed = start_time.elapsed();
    println!("Elapsed time: {:.2?}", elapsed);
    Ok(())
}

fn available_atomic_parts(cfg: &PrepareConfig) -> FxHashSet<AtomicLabelPart> {
    let mut parts: FxHashSet<AtomicLabelPart> = FxHashSet::default();
    parts.insert(AtomicLabelPart::Input);
    if cfg.near.is_some() {
        parts.insert(AtomicLabelPart::NearSide);
        if !cfg.near_group_cols.is_empty() {
            parts.insert(AtomicLabelPart::NearName);
        }
    }
    if cfg.distance_bins.is_some() {
        parts.insert(AtomicLabelPart::Bin);
    }
    if cfg.cluster_min_overlaps.is_some() {
        parts.insert(AtomicLabelPart::Cluster);
    }
    parts
}
