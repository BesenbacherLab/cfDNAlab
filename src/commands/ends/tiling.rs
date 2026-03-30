use crate::commands::{
    counters::EndsCounters,
    ends::counting::{
        EncodedEndMotifKey, EndCountsByWindow, EndMotifCounts, TileEndMotifCountEntry,
        TileWindowEndCounts,
    },
};
use anyhow::{Context, Result};
use bincode::{
    config::standard,
    serde::{decode_from_std_read, encode_into_std_write},
};
use fxhash::FxHashMap;
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
pub struct TileResult {
    pub chr: String,
    pub counts_path: PathBuf,
    pub counter: EndsCounters,
}

/// Persist per-tile end-motif counts so they can be merged after parallel tile processing.
///
/// Parameters
/// ----------
/// - `path`:
///   Destination for the serialized tile payload
/// - `payload`:
///   Sparse per-window counts for one tile
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after the tile payload has been flushed to disk
pub fn serialize_tile_counts(path: &Path, payload: &[TileWindowEndCounts]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("creating tile counts file: {}", path.display()))?;
    let mut writer = BufWriter::with_capacity(512 * 1024, file);
    encode_into_std_write(payload, &mut writer, standard())
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
///   Location of the serialized tile payload
///
/// Returns
/// -------
/// - `Result<Vec<TileWindowEndCounts>>`:
///   The decoded sparse tile payload
pub fn deserialize_tile_counts(path: &Path) -> Result<Vec<TileWindowEndCounts>> {
    let file = File::open(path)
        .with_context(|| format!("opening tile counts file: {}", path.display()))?;
    let mut reader = BufReader::with_capacity(512 * 1024, file);
    decode_from_std_read(&mut reader, standard())
        .with_context(|| format!("deserialising tile counts from {}", path.display()))
}

/// Convert sparse per-window count maps into a stable serialized payload.
///
/// Parameters
/// ----------
/// - `counts_by_window`:
///   Sparse per-window motif counts accumulated for one tile
///
/// Returns
/// -------
/// - `Vec<TileWindowEndCounts>`:
///   Stable, sorted payload ready for serialization
pub fn build_tile_payload(
    counts_by_window: FxHashMap<u64, EndMotifCounts>,
) -> Vec<TileWindowEndCounts> {
    let mut payload: Vec<TileWindowEndCounts> = counts_by_window
        .into_iter()
        .map(|(original_idx, counts)| {
            let mut entries: Vec<TileEndMotifCountEntry> = counts
                .counts
                .into_iter()
                .map(TileEndMotifCountEntry::from)
                .collect();

            entries.sort_unstable_by_key(|entry| {
                (
                    entry.inside_code,
                    entry.outside_code,
                    entry.reverse_on_decode,
                )
            });

            TileWindowEndCounts {
                original_idx,
                entries,
            }
        })
        .collect();

    payload.sort_unstable_by_key(|window_counts| window_counts.original_idx);
    payload
}

/// Merge one serialized tile payload into the reduced sparse counts.
///
/// Tile payloads already carry global window ids, so merging is just a sparse
/// sum over `(window, motif-key)` entries.
///
/// Parameters
/// ----------
/// - `merged`:
///   Reduced sparse counts updated in place
/// - `tile_payload`:
///   One tile's serialized sparse counts
///
/// Returns
/// -------
/// - `Result<()>`:
///   `Ok(())` after all counts have been merged
pub fn merge_tile_payload(
    merged: &mut EndCountsByWindow,
    tile_payload: Vec<TileWindowEndCounts>,
) -> Result<()> {
    for window_counts in tile_payload {
        let dst = merged.entry(window_counts.original_idx).or_default();
        for entry in window_counts.entries {
            dst.incr_weighted(EncodedEndMotifKey::from(&entry), entry.value);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    include!("tiling_tests.rs");
}
