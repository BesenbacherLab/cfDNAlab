use crate::commands::gc_bias::config::GCConfig;
use crate::commands::gc_bias::counting::GCCounts;
use crate::commands::gc_bias::gc_bias::process_window;
use crate::shared::tiled_run::parse_tile_index;
use anyhow::{Result, ensure};
use fxhash::{FxHashMap, FxHashSet};
use ndarray::{Array1, Array2};
use ndarray_npy::{NpzReader, NpzWriter};
use std::collections::hash_map::Entry;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

#[derive(Clone)]
pub struct CrossingPart {
    pub idx: usize,
    pub counts: GCCounts,
}

pub fn write_crossing_parts(
    temp_dir: &PathBuf,
    tile_idx: u32,
    template: &GCCounts,
    parts: &[CrossingPart],
) -> Result<Option<PathBuf>> {
    if parts.is_empty() {
        return Ok(None);
    }
    let path = temp_dir.join(format!("cross.{}.npz", tile_idx));
    let file = File::create(&path)?;
    let mut npz = NpzWriter::new(file);

    let counts_len = template.borrow_raw_counts().len();
    let mut idxs: Vec<u64> = Vec::with_capacity(parts.len());
    let mut acgt0: Vec<u64> = Vec::with_capacity(parts.len());
    let mut acgt1: Vec<u64> = Vec::with_capacity(parts.len());
    let mut counts_arr = Array2::zeros((parts.len(), counts_len));
    for (row_idx, part) in parts.iter().enumerate() {
        idxs.push(part.idx as u64);
        acgt0.push(part.counts.num_acgt_out_of.0);
        acgt1.push(part.counts.num_acgt_out_of.1);
        let counts = part.counts.borrow_raw_counts();
        ensure!(
            counts.len() == counts_len,
            "Crossing part counts length {} did not match expected {}",
            counts.len(),
            counts_len
        );
        counts_arr
            .row_mut(row_idx)
            .assign(&ndarray::ArrayView1::from(counts.as_slice()));
    }

    npz.add_array("idx", &Array1::from(idxs))?;
    npz.add_array("acgt0", &Array1::from(acgt0))?;
    npz.add_array("acgt1", &Array1::from(acgt1))?;
    npz.add_array("counts", &counts_arr)?;

    npz.finish()?;
    Ok(Some(path))
}

/// Merge per-tile crossing window fragments into fully scaled windows.
///
/// Each NPZ file contains partial counts for windows that extended outside a tile.
/// Files are sorted by tile index so we can stream them once: when a window index
/// disappears from the current file, all its parts have been seen and the merged
/// counts can be scaled, accumulated and released to reduce memory.
/// Every NPZ row corresponds to a single window within that tile. Duplicates are
/// merged across files, not within a file.
///
/// Parameters
/// ----------
/// - `files`:
///     NPZ archives written by `write_crossing_parts`, one per tile that had spillover.
/// - `template`:
///     Shape/source for constructing zeroed counts and validating incoming buffers.
/// - `opt`:
///     GC bias configuration used when scaling a finalized window.
/// - `avg_window_span`:
///     Average genomic span of a window, used to normalize counts.
///
/// Returns
/// -------
/// - `(GCCounts, usize)`:
///     Sum of scaled windows and the number of windows contributing to that sum.
pub fn stream_crossing_files(
    mut files: Vec<PathBuf>,
    template: &GCCounts,
    opt: &GCConfig,
    avg_window_span: f64,
) -> Result<(GCCounts, usize)> {
    if files.is_empty() {
        return Ok((template.zeroed_like()?, 0));
    }
    // Sort files so tile order matches genomic order for streaming flush logic
    files.sort_by_key(|p| {
        parse_tile_index(p.file_name().and_then(|s| s.to_str()).unwrap_or("")).unwrap_or(u32::MAX)
    });

    // Accumulates scaled windows and the total window count
    let mut total_sum = template.zeroed_like()?;
    let mut total_weight = 0usize;

    // Tracks partial windows keyed by window index while streaming across tiles
    let mut window_parts: FxHashMap<usize, GCCounts> = FxHashMap::default();

    for path in files {
        // Marks which window indices appeared in this tile's file
        let mut seen_in_file: FxHashSet<usize> = FxHashSet::default();
        let file = File::open(&path)?;
        let mut npz = NpzReader::new(BufReader::new(file))?;
        let idxs: Array1<u64> = npz.by_name("idx")?;
        let acgt0: Array1<u64> = npz.by_name("acgt0")?;
        let acgt1: Array1<u64> = npz.by_name("acgt1")?;
        let counts: Array2<f64> = npz.by_name("counts")?;
        ensure!(
            counts.dim().1 == template.borrow_raw_counts().len(),
            "Counts matrix width mismatch in {:?}",
            path
        );
        ensure!(
            idxs.len() == counts.nrows() && acgt0.len() == idxs.len() && acgt1.len() == idxs.len(),
            "Crossing file {:?} had inconsistent vector lengths",
            path
        );

        window_parts.reserve(idxs.len());

        for (((idx, &ac0), &ac1), row) in idxs
            .iter()
            .zip(acgt0.iter())
            .zip(acgt1.iter())
            .zip(counts.outer_iter())
        {
            let row = row.to_vec();
            let counts = GCCounts::from_parts(
                row,
                template.length_min,
                template.length_max,
                template.end_offset(),
                (ac0, ac1),
            )?;

            let idx = *idx as usize;
            seen_in_file.insert(idx);
            // Merge parts for the same window index across tiles
            match window_parts.entry(idx) {
                Entry::Vacant(slot) => {
                    slot.insert(counts);
                }
                Entry::Occupied(mut slot) => {
                    slot.get_mut().merge_from(&counts)?;
                }
            }
        }

        if !window_parts.is_empty() {
            // Anything absent from this file is finished and can be flushed
            let mut finished_idxs: Vec<usize> = Vec::new();
            finished_idxs.reserve(window_parts.len());
            for idx in window_parts.keys() {
                if !seen_in_file.contains(idx) {
                    finished_idxs.push(*idx);
                }
            }
            for idx in finished_idxs {
                if let Some(counts) = window_parts.remove(&idx) {
                    // Scale and accumulate completed window before dropping it
                    if let Some(scaled) = process_window(counts, opt, Some(avg_window_span))? {
                        total_sum.merge_from(&scaled)?;
                        total_weight += 1;
                    }
                }
            }
        }

        let _ = std::fs::remove_file(&path);
    }

    for counts in window_parts.into_values() {
        // Flush any windows that persisted through every tile
        if let Some(scaled) = process_window(counts, opt, Some(avg_window_span))? {
            total_sum.merge_from(&scaled)?;
            total_weight += 1;
        }
    }

    Ok((total_sum, total_weight))
}
