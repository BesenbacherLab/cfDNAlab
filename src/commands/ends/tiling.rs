use crate::{
    commands::{
        counters::EndsCounters,
        ends::counting::{
            EndCountsByWindow, EndMotifCounts, SelectedEndCountsByWindow, TileEndMotifCountEntry,
            TileWindowEndCounts,
        },
    },
    shared::kmers::motifs_file::EncodedMotifKey,
};
use anyhow::{Context, Result};
use bincode::{
    config::standard,
    serde::{decode_from_std_read, encode_into_std_write},
};
use fxhash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

/// Per-tile bookkeeping for intermediate sparse motif count files and fragment counters.
///
/// Each parallel tile writes its sparse counts to disk and returns one of these
/// structs so the outer reducer can later merge the files and combine the
/// command statistics.
pub(crate) struct TileResult {
    pub(crate) counts_path: PathBuf,
    pub(crate) counter: EndsCounters,
}

/// Persist per-tile end-motif counts so they can be merged after parallel tile processing.
///
/// Parameters
/// ----------
/// - `path`:
///   Destination for the serialized tile count records
/// - `count_records`:
///   Sparse per-window counts for one tile
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after the tile count records have been flushed to disk
pub(crate) fn serialize_tile_counts(
    path: &Path,
    count_records: &[TileWindowEndCounts],
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("creating tile counts file: {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(512 * 1024, file);
    encode_into_std_write(count_records, &mut writer, standard())
        .with_context(|| format!("serialising tile counts to {}", path.display()))?;
    writer.flush().with_context(|| {
        format!(
            "flushing tile counts file after serialisation: {}",
            path.display()
        )
    })
}

/// Load counts created by [`serialize_tile_counts`] during the reduction phase.
///
/// Parameters
/// ----------
/// - `path`:
///   Location of the serialized tile count records
///
/// Returns
/// -------
/// - `Result<Vec<TileWindowEndCounts>>`:
///   Decoded sparse tile count records
pub(crate) fn deserialize_tile_counts(path: &Path) -> Result<Vec<TileWindowEndCounts>> {
    let file = File::open(path)
        .with_context(|| format!("opening tile counts file: {}", path.display()))?;
    let mut reader = BufReader::with_capacity(512 * 1024, file);
    decode_from_std_read(&mut reader, standard())
        .with_context(|| format!("deserialising tile counts from {}", path.display()))
}

/// Serialized tile entry for one selected motif-file target.
///
/// This is the compact on-disk form used only for `--motifs-file` runs. `target_idx` is the
/// original target index assigned by the parser, not the final compact output column.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct TileSelectedEndMotifCountEntry {
    /// Original motifs-file target index
    pub target_idx: u32,
    /// Weighted count accumulated for this target in one output row
    pub value: f64,
}

/// Serialized selected-target counts for one output window in one tile.
///
/// Tile workers write these structs to temporary files so parallel counting does not hold all tile
/// results in memory. `entries` are sorted by target index before serialization for deterministic
/// count files and easier test fixtures.
#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct TileWindowSelectedEndCounts {
    /// Global output-row index for this tile-local count row
    pub original_idx: u64,
    /// Non-empty selected target counts for the row
    pub entries: Vec<TileSelectedEndMotifCountEntry>,
}

/// Persist per-tile selected end-motif counts so they can be reduced later.
///
/// This is the selected-motif equivalent of [`serialize_tile_counts`]. It uses a separate record
/// shape because selected counts are already mapped to numeric target indices and do not need to
/// preserve encoded inside and outside motif halves.
///
/// Parameters
/// ----------
/// - `path`:
///   Destination for the serialized selected-count records
/// - `count_records`:
///   Sparse selected-target counts for one tile
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after the selected-count records have been flushed to disk
pub(crate) fn serialize_selected_tile_counts(
    path: &Path,
    count_records: &[TileWindowSelectedEndCounts],
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("creating selected tile counts file: {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(512 * 1024, file);
    encode_into_std_write(count_records, &mut writer, standard())
        .with_context(|| format!("serialising selected tile counts to {}", path.display()))?;
    writer.flush().with_context(|| {
        format!(
            "flushing selected tile counts file after serialisation: {}",
            path.display()
        )
    })
}

/// Load selected counts created by [`serialize_selected_tile_counts`].
///
/// Parameters
/// ----------
/// - `path`:
///   Location of the serialized selected-count records
///
/// Returns
/// -------
/// - `Result<Vec<TileWindowSelectedEndCounts>>`:
///   Decoded sparse selected-count records
pub(crate) fn deserialize_selected_tile_counts(
    path: &Path,
) -> Result<Vec<TileWindowSelectedEndCounts>> {
    let file = File::open(path)
        .with_context(|| format!("opening selected tile counts file: {}", path.display()))?;
    let mut reader = BufReader::with_capacity(512 * 1024, file);
    decode_from_std_read(&mut reader, standard())
        .with_context(|| format!("deserialising selected tile counts from {}", path.display()))
}

/// Convert sparse per-window count maps into stable serialized count records.
///
/// Parameters
/// ----------
/// - `counts_by_window`:
///   Sparse per-window motif counts accumulated for one tile
///
/// Returns
/// -------
/// - `Vec<TileWindowEndCounts>`:
///   Stable, sorted count records ready for serialization
pub(crate) fn build_tile_count_records(
    counts_by_window: FxHashMap<u64, EndMotifCounts>,
) -> Vec<TileWindowEndCounts> {
    let mut count_records: Vec<TileWindowEndCounts> = counts_by_window
        .into_iter()
        .filter_map(|(original_idx, counts)| {
            let mut entries: Vec<TileEndMotifCountEntry> = counts
                .counts
                .into_iter()
                .map(TileEndMotifCountEntry::from)
                .collect();
            if entries.is_empty() {
                return None;
            }

            entries.sort_unstable_by_key(|entry| {
                (
                    entry.inside_code,
                    entry.outside_code,
                    entry.reverse_on_decode,
                )
            });

            Some(TileWindowEndCounts {
                original_idx,
                entries,
            })
        })
        .collect();

    count_records.sort_unstable_by_key(|window_counts| window_counts.original_idx);
    count_records
}

/// Convert selected per-window count maps into stable serialized count records.
///
/// Empty rows are omitted because the final writer already knows the full row
/// count. Rows and target entries are sorted to keep temporary files deterministic across hash-map
/// iteration orders.
///
/// Parameters
/// ----------
/// - `counts_by_window`:
///   Sparse selected-target counts accumulated for one tile
///
/// Returns
/// -------
/// - `Vec<TileWindowSelectedEndCounts>`:
///   Stable, sorted selected-count records ready for serialization
pub(crate) fn build_selected_tile_count_records(
    counts_by_window: SelectedEndCountsByWindow,
) -> Vec<TileWindowSelectedEndCounts> {
    let mut count_records: Vec<TileWindowSelectedEndCounts> = counts_by_window
        .into_iter()
        .filter_map(|(original_idx, counts)| {
            let mut entries: Vec<TileSelectedEndMotifCountEntry> = counts
                .into_iter()
                .map(|(target_idx, value)| TileSelectedEndMotifCountEntry { target_idx, value })
                .collect();
            if entries.is_empty() {
                // Empty rows carry no information because output shape is tracked separately
                return None;
            }

            // Stable order is useful for deterministic temporary count files and tests
            entries.sort_unstable_by_key(|entry| entry.target_idx);

            Some(TileWindowSelectedEndCounts {
                original_idx,
                entries,
            })
        })
        .collect();

    // Stable row order also makes the final reduction independent of hash-map iteration order
    count_records.sort_unstable_by_key(|window_counts| window_counts.original_idx);
    count_records
}

/// Merge one tile's serialized count records into the reduced sparse counts.
///
/// Tile count records already carry global window ids, so merging is just a sparse
/// sum over `(window, motif-key)` entries.
///
/// Parameters
/// ----------
/// - `merged`:
///   Reduced sparse counts updated in place
/// - `tile_count_records`:
///   One tile's serialized sparse counts
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after all counts have been merged
pub(crate) fn merge_tile_count_records(
    merged: &mut EndCountsByWindow,
    tile_count_records: Vec<TileWindowEndCounts>,
) -> Result<()> {
    for window_counts in tile_count_records {
        for entry in window_counts.entries {
            if EndMotifCounts::should_store_weight(entry.value)? {
                merged
                    .entry(window_counts.original_idx)
                    .or_default()
                    .incr_weighted(EncodedMotifKey::from(&entry), entry.value);
            }
        }
    }

    Ok(())
}

/// Merge one tile's serialized selected-count records into reduced selected counts.
///
/// The selected-count records already use global row ids, so reduction is a sparse sum over
/// `(row, target_idx)`. Weight validation is repeated here so corrupt or future temporary count
/// files cannot silently create invalid counts.
///
/// Parameters
/// ----------
/// - `merged`:
///   Reduced sparse selected-target counts updated in place
/// - `tile_count_records`:
///   One tile's serialized selected-target counts
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after all valid weights have been merged
pub(crate) fn merge_selected_tile_count_records(
    merged: &mut SelectedEndCountsByWindow,
    tile_count_records: Vec<TileWindowSelectedEndCounts>,
) -> Result<()> {
    for window_counts in tile_count_records {
        for entry in window_counts.entries {
            // Keep the same sparse-count rules as full motif counting
            if EndMotifCounts::should_store_weight(entry.value)? {
                merged
                    .entry(window_counts.original_idx)
                    .or_default()
                    .entry(entry.target_idx)
                    .and_modify(|value| *value += entry.value)
                    .or_insert(entry.value);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    include!("tiling_tests.rs");
}
