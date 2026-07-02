use crate::{
    commands::ref_kmers::counting::{KmerCounts, KmerCountsByWindow, SelectedKmerCountsByWindow},
    shared::kmers::kmer_codec::{Kmer, KmerOrientation},
};
use anyhow::{Context, Result};
use bincode::{
    config::standard,
    serde::{decode_from_std_read, encode_into_std_write},
};
use serde::{Deserialize, Serialize};
use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

/// Per-tile bookkeeping for intermediate reference k-mer count files.
///
/// Parallel tile workers write sparse count records to disk and return this
/// lightweight handle so the outer reducer can merge the files after tile processing.
#[derive(Debug)]
pub(crate) struct TileResult {
    pub(crate) counts_path: PathBuf,
}

/// Serialized tile entry for one reference k-mer.
///
/// This is the compact on-disk form used between tile counting and final reduction. The orientation
/// is currently always forward for `ref-kmers`, but keeping it in the record matches the shared
/// k-mer key and avoids inventing another temporary key shape.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct TileKmerCountEntry {
    /// K-mer size
    pub(crate) k: u8,
    /// Encoded k-mer code
    pub(crate) code: u64,
    /// K-mer orientation in the shared key
    pub(crate) orientation: KmerOrientation,
    /// Count accumulated for this key in one output row
    pub(crate) value: f64,
}

impl From<(Kmer, f64)> for TileKmerCountEntry {
    fn from((kmer, value): (Kmer, f64)) -> Self {
        Self {
            k: kmer.k,
            code: kmer.code,
            orientation: kmer.orientation,
            value,
        }
    }
}

impl From<&TileKmerCountEntry> for Kmer {
    fn from(entry: &TileKmerCountEntry) -> Self {
        Self {
            k: entry.k,
            code: entry.code,
            orientation: entry.orientation,
        }
    }
}

/// Serialized sparse counts for one output row in one tile.
///
/// `original_idx` is the final output row id. For fixed-size windows this already includes the
/// chromosome row offset. For BED-like windows it is the original BED row id or grouped BED group id.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct TileWindowKmerCounts {
    /// Global output-row index for this tile-local count row
    pub(crate) original_idx: u64,
    /// Non-empty reference k-mer counts for the row
    pub(crate) entries: Vec<TileKmerCountEntry>,
}

/// Serialized tile entry for one selected motifs-file target.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct TileSelectedKmerCountEntry {
    /// Original motifs-file target index
    pub(crate) target_idx: u32,
    /// Count accumulated for this target in one output row
    pub(crate) value: f64,
}

/// Serialized selected-target counts for one output row in one tile.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct TileWindowSelectedKmerCounts {
    /// Global output-row index for this tile-local count row
    pub(crate) original_idx: u64,
    /// Non-empty selected target counts for the row
    pub(crate) entries: Vec<TileSelectedKmerCountEntry>,
}

/// Persist per-tile reference k-mer counts so they can be reduced later.
pub(crate) fn serialize_tile_counts(
    path: &Path,
    count_records: &[TileWindowKmerCounts],
) -> Result<()> {
    let file = File::create(path).with_context(|| {
        format!(
            "creating reference k-mer tile counts file: {}",
            path.display()
        )
    })?;
    let mut writer = BufWriter::with_capacity(512 * 1024, file);
    encode_into_std_write(count_records, &mut writer, standard()).with_context(|| {
        format!(
            "serialising reference k-mer tile counts to {}",
            path.display()
        )
    })?;
    writer.flush().with_context(|| {
        format!(
            "flushing reference k-mer tile counts file after serialisation: {}",
            path.display()
        )
    })
}

/// Load counts created by [`serialize_tile_counts`] during reduction.
pub(crate) fn deserialize_tile_counts(path: &Path) -> Result<Vec<TileWindowKmerCounts>> {
    let file = File::open(path).with_context(|| {
        format!(
            "opening reference k-mer tile counts file: {}",
            path.display()
        )
    })?;
    let mut reader = BufReader::with_capacity(512 * 1024, file);
    decode_from_std_read(&mut reader, standard()).with_context(|| {
        format!(
            "deserialising reference k-mer tile counts from {}",
            path.display()
        )
    })
}

/// Persist per-tile selected reference k-mer counts so they can be reduced later.
pub(crate) fn serialize_selected_tile_counts(
    path: &Path,
    count_records: &[TileWindowSelectedKmerCounts],
) -> Result<()> {
    let file = File::create(path).with_context(|| {
        format!(
            "creating selected reference k-mer tile counts file: {}",
            path.display()
        )
    })?;
    let mut writer = BufWriter::with_capacity(512 * 1024, file);
    encode_into_std_write(count_records, &mut writer, standard()).with_context(|| {
        format!(
            "serialising selected reference k-mer tile counts to {}",
            path.display()
        )
    })?;
    writer.flush().with_context(|| {
        format!(
            "flushing selected reference k-mer tile counts file after serialisation: {}",
            path.display()
        )
    })
}

/// Load selected counts created by [`serialize_selected_tile_counts`] during reduction.
pub(crate) fn deserialize_selected_tile_counts(
    path: &Path,
) -> Result<Vec<TileWindowSelectedKmerCounts>> {
    let file = File::open(path).with_context(|| {
        format!(
            "opening selected reference k-mer tile counts file: {}",
            path.display()
        )
    })?;
    let mut reader = BufReader::with_capacity(512 * 1024, file);
    decode_from_std_read(&mut reader, standard()).with_context(|| {
        format!(
            "deserialising selected reference k-mer tile counts from {}",
            path.display()
        )
    })
}

/// Convert sparse per-window k-mer counts into stable serialized tile records.
pub(crate) fn build_tile_count_records(
    counts_by_window: KmerCountsByWindow,
) -> Vec<TileWindowKmerCounts> {
    let mut count_records: Vec<TileWindowKmerCounts> = counts_by_window
        .into_iter()
        .filter_map(|(original_idx, counts)| {
            let mut entries: Vec<TileKmerCountEntry> = counts
                .counts
                .into_iter()
                .map(TileKmerCountEntry::from)
                .collect();
            if entries.is_empty() {
                return None;
            }

            entries.sort_unstable_by_key(|entry| {
                (
                    entry.k,
                    entry.code,
                    kmer_orientation_sort_key(entry.orientation),
                )
            });
            Some(TileWindowKmerCounts {
                original_idx,
                entries,
            })
        })
        .collect();

    count_records.sort_unstable_by_key(|window_counts| window_counts.original_idx);
    count_records
}

/// Convert sparse per-window selected-target counts into stable serialized tile records.
pub(crate) fn build_selected_tile_count_records(
    counts_by_window: SelectedKmerCountsByWindow,
) -> Vec<TileWindowSelectedKmerCounts> {
    let mut count_records: Vec<TileWindowSelectedKmerCounts> = counts_by_window
        .into_iter()
        .filter_map(|(original_idx, counts)| {
            let mut entries: Vec<TileSelectedKmerCountEntry> = counts
                .counts
                .into_iter()
                .map(|(target_idx, value)| TileSelectedKmerCountEntry { target_idx, value })
                .collect();
            if entries.is_empty() {
                return None;
            }

            entries.sort_unstable_by_key(|entry| entry.target_idx);
            Some(TileWindowSelectedKmerCounts {
                original_idx,
                entries,
            })
        })
        .collect();

    count_records.sort_unstable_by_key(|window_counts| window_counts.original_idx);
    count_records
}

/// Merge one tile's serialized k-mer counts into reduced sparse counts.
pub(crate) fn merge_tile_count_records(
    merged: &mut KmerCountsByWindow,
    tile_count_records: Vec<TileWindowKmerCounts>,
) -> Result<()> {
    for window_counts in tile_count_records {
        for entry in window_counts.entries {
            if KmerCounts::should_store_weight(entry.value)? {
                merged
                    .entry(window_counts.original_idx)
                    .or_default()
                    .incr_weighted(Kmer::from(&entry), entry.value);
            }
        }
    }
    Ok(())
}

/// Merge one tile's serialized selected-target counts into reduced sparse counts.
pub(crate) fn merge_selected_tile_count_records(
    merged: &mut SelectedKmerCountsByWindow,
    tile_count_records: Vec<TileWindowSelectedKmerCounts>,
) -> Result<()> {
    for window_counts in tile_count_records {
        for entry in window_counts.entries {
            if KmerCounts::should_store_weight(entry.value)? {
                merged
                    .entry(window_counts.original_idx)
                    .or_default()
                    .incr_weighted(entry.target_idx, entry.value);
            }
        }
    }
    Ok(())
}

fn kmer_orientation_sort_key(orientation: KmerOrientation) -> u8 {
    match orientation {
        KmerOrientation::Forward => 0,
        KmerOrientation::Reverse => 1,
    }
}

#[cfg(test)]
mod tests {
    include!("tiling_tests.rs");
}
