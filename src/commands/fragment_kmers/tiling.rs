use crate::{
    commands::counters::FragmentKmersCounters,
    shared::{
        kmers::{
            kmer_codec::{Kmer, KmerOrientation, KmerSpec},
            process_counts::{DecodedCounts, split_and_decode_counts},
        },
        positioning::PositionGroup,
    },
};
use anyhow::{Context, Result, bail};
use bincode::{
    config::standard,
    serde::{decode_from_std_read, encode_into_std_write},
};
use fxhash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::{
    fs::{self, File},
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct TileKmerCountEntry {
    pub(crate) k: u8,
    pub(crate) code: u64,
    pub(crate) position: Option<i32>,
    pub(crate) group: PositionGroup,
    pub(crate) value: f64,
}

#[cfg_attr(not(test), doc(hidden))]
#[derive(Debug, Serialize, Deserialize)]
pub(crate) struct TileWindowCounts {
    pub(crate) original_idx: u64,
    pub(crate) entries: Vec<TileKmerCountEntry>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct CountKey {
    pub(crate) k: u8,
    pub(crate) code: u64,
    pub(crate) position: Option<i32>,
    pub(crate) group: PositionGroup,
}

impl CountKey {
    #[inline]
    pub(crate) fn as_kmer(self) -> Kmer {
        Kmer {
            k: self.k,
            code: self.code,
            orientation: self.orientation(),
        }
    }

    // Get orientation based on `group`
    pub(crate) fn orientation(&self) -> KmerOrientation {
        KmerOrientation::from_position_group(self.group)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct PositionDescriptor {
    pub(crate) group: PositionGroup,
    pub(crate) offset: i32,
}

impl From<&TileKmerCountEntry> for CountKey {
    fn from(entry: &TileKmerCountEntry) -> Self {
        Self {
            k: entry.k,
            code: entry.code,
            position: entry.position,
            group: entry.group,
        }
    }
}

impl From<(CountKey, f64)> for TileKmerCountEntry {
    fn from((key, value): (CountKey, f64)) -> Self {
        Self {
            k: key.k,
            code: key.code,
            position: key.position,
            group: key.group,
            value,
        }
    }
}

/// Persist per-tile k-mer counts so they can be merged after parallel tile processing.
pub(crate) fn serialize_tile_counts(path: &Path, payload: &[TileWindowCounts]) -> Result<()> {
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
pub(crate) fn deserialize_tile_counts(path: &Path) -> Result<Vec<TileWindowCounts>> {
    let file = File::open(path)
        .with_context(|| format!("opening tile counts file: {}", path.display()))?;
    let mut reader = BufReader::with_capacity(512 * 1024, file);
    decode_from_std_read(&mut reader, standard())
        .with_context(|| format!("deserialising tile counts from {}", path.display()))
}

/// Per-tile bookkeeping for intermediate count files and fragment counters.
pub(crate) struct TileResult {
    pub(crate) chr: String,
    pub(crate) counts_path: Option<PathBuf>,
    pub(crate) counter: FragmentKmersCounters,
}

/// Reduce per-tile count payloads into a dense vector aligned with the global window order.
#[cfg_attr(not(test), doc(hidden))]
pub(crate) fn merge_tile_counts<I>(
    payloads: I,
    total_windows: usize,
    kmer_specs: &FxHashMap<u8, KmerSpec>,
) -> Result<Vec<DecodedCounts>>
where
    I: IntoIterator<Item = Vec<TileWindowCounts>>,
{
    let mut aggregated_counts: FxHashMap<u64, FxHashMap<CountKey, f64>> = FxHashMap::default();

    for payload in payloads {
        for window_counts in payload {
            let entry = aggregated_counts
                .entry(window_counts.original_idx)
                .or_default();
            for count in window_counts.entries {
                let key = CountKey::from(&count);
                debug_assert!(
                    key.position.is_none(),
                    "merge_tile_counts received positional entry in non-positional mode"
                );
                *entry.entry(key).or_insert(0.0) += count.value;
            }
        }
    }

    let empty_counts: FxHashMap<Kmer, f64> = FxHashMap::default();
    let mut all_bins: Vec<DecodedCounts> = Vec::with_capacity(total_windows);
    for idx in 0..total_windows {
        if let Some(counts) = aggregated_counts.remove(&(idx as u64)) {
            let mut plain_counts: FxHashMap<Kmer, f64> = FxHashMap::default();
            plain_counts.reserve(counts.len());
            for (key, value) in counts {
                debug_assert!(
                    key.position.is_none(),
                    "merge_tile_counts received positional entry in non-positional mode"
                );
                let kmer = key.as_kmer();
                *plain_counts.entry(kmer).or_insert(0.0) += value;
            }
            all_bins.push(split_and_decode_counts(&plain_counts, kmer_specs));
        } else {
            all_bins.push(split_and_decode_counts(&empty_counts, kmer_specs));
        }
    }

    if !aggregated_counts.is_empty() {
        bail!(
            "Received counts for unexpected window indices: {:?}",
            aggregated_counts.keys().collect::<Vec<&u64>>()
        );
    }

    Ok(all_bins)
}

#[cfg_attr(not(test), doc(hidden))]
pub(crate) fn merge_tile_counts_positional<I>(
    payloads: I,
    total_windows: usize,
) -> Result<Vec<FxHashMap<PositionDescriptor, FxHashMap<Kmer, f64>>>>
where
    I: IntoIterator<Item = Vec<TileWindowCounts>>,
{
    let mut aggregated_counts: FxHashMap<u64, FxHashMap<CountKey, f64>> = FxHashMap::default();

    for payload in payloads {
        for window_counts in payload {
            let entry = aggregated_counts
                .entry(window_counts.original_idx)
                .or_default();
            for count in window_counts.entries {
                let key = CountKey::from(&count);
                *entry.entry(key).or_insert(0.0) += count.value;
            }
        }
    }

    let mut all_bins: Vec<FxHashMap<PositionDescriptor, FxHashMap<Kmer, f64>>> =
        Vec::with_capacity(total_windows);
    for idx in 0..total_windows {
        if let Some(counts) = aggregated_counts.remove(&(idx as u64)) {
            let mut by_position: FxHashMap<PositionDescriptor, FxHashMap<Kmer, f64>> =
                FxHashMap::default();
            for (key, value) in counts {
                let group = key.group;
                let offset = match key.position {
                    Some(offset) => offset,
                    _ => bail!(
                        "Positional merge encountered entry without position for window {}",
                        idx
                    ),
                };
                let descriptor = PositionDescriptor { group, offset };
                let kmer = key.as_kmer();
                let entry = by_position.entry(descriptor).or_default();
                *entry.entry(kmer).or_insert(0.0) += value;
            }
            all_bins.push(by_position);
        } else {
            all_bins.push(FxHashMap::default());
        }
    }

    if !aggregated_counts.is_empty() {
        bail!(
            "Received counts for unexpected window indices: {:?}",
            aggregated_counts.keys().collect::<Vec<&u64>>()
        );
    }

    Ok(all_bins)
}

pub(crate) fn reduce_chromosome_tile_results(
    tile_results: Vec<TileResult>,
) -> Result<Vec<TileWindowCounts>> {
    let mut aggregated: FxHashMap<u64, FxHashMap<CountKey, f64>> = FxHashMap::default();

    for tile_result in tile_results {
        let Some(path) = tile_result.counts_path else {
            continue;
        };
        let tile_payload = deserialize_tile_counts(&path)?;
        let _ = fs::remove_file(&path);

        for window_counts in tile_payload {
            let entry = aggregated.entry(window_counts.original_idx).or_default();
            for count in window_counts.entries {
                let key = CountKey::from(&count);
                *entry.entry(key).or_insert(0.0) += count.value;
            }
        }
    }

    let mut merged: Vec<TileWindowCounts> = aggregated
        .into_iter()
        .map(|(original_idx, counts_map)| {
            let mut entries: Vec<TileKmerCountEntry> = Vec::with_capacity(counts_map.len());
            for (kmer, value) in counts_map {
                entries.push(TileKmerCountEntry::from((kmer, value)));
            }
            TileWindowCounts {
                original_idx,
                entries,
            }
        })
        .collect();

    merged.sort_unstable_by_key(|window| window.original_idx);
    Ok(merged)
}

#[cfg(test)]
mod tests {
    include!("tiling_tests.rs");
}
